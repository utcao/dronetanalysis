#!/usr/bin/env python3
"""
Stage 7b: Collect Permutation Test Results and Compute Empirical P-values.

Reads per-gene HDF5 permutation null distributions (produced by Stage 7),
computes two-sided empirical p-values, and writes two output tables:

  permutation_pvals.tsv         — per-gene empirical p-values + null summaries
  rewiring_hubs_with_pvals.tsv  — rewiring_hubs.tsv augmented with pval_* columns

The empirical p-value formula mirrors utils_permutation_net.R::permutate_p_two():
    p = max(sum(|null| >= |obs|) / n_perm,  1 / n_perm)

Pipeline position
-----------------
Stage 7   07_permutation_test.py    →  permutation_null/{gi}_{gene_id}.h5
Stage 7b  THIS SCRIPT               →  permutation_pvals.tsv
                                        rewiring_hubs_with_pvals.tsv

Usage
-----
python 08_collect_permutation_pvals.py \\
    --null-dir results/permutation_null \\
    --observed-tsv results/rewiring_hubs.tsv \\
    --out-pvals-tsv results/permutation_pvals.tsv \\
    --out-augmented-tsv results/rewiring_hubs_with_pvals.tsv \\
    --alpha 0.05
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# Metrics tested (must match 07_permutation_test.py base_metrics)
BASE_METRICS = [
    "n_diff_edges",
    "focus_deg_low",
    "focus_deg_high",
    "mean_abs_delta",
    "max_abs_delta",
    "L1_n_nodes",
    "L1_rewire",
    "L1_frac_rewire",
]
L2_METRICS = ["L2L1_deg", "L2L1_rewire"]


def permutation_pval(obs: float, null: np.ndarray) -> float:
    """
    Two-sided empirical p-value.

    p = max(sum(|null| >= |obs|) / n_perm,  1 / n_perm)

    Mirrors utils_permutation_net.R::permutate_p_two().
    Floor at 1/n_perm so p can never be exactly 0.
    """
    n = len(null)
    if n == 0:
        return float("nan")
    count = int(np.sum(np.abs(null) >= np.abs(obs)))
    return max(count / n, 1.0 / n)


def load_gene_null(h5_path: Path) -> dict | None:
    """
    Load null distribution and observed values from a per-gene HDF5 file.

    Returns None if the file does not exist.
    Returns a dict with keys: gene_id, n_permutations, skipped, observed, null.
    """
    if not h5_path.exists():
        return None

    with h5py.File(h5_path, "r") as h5:
        meta = h5["meta"]
        gene_id = meta.attrs["gene_id"]
        if isinstance(gene_id, bytes):
            gene_id = gene_id.decode()
        gene_index   = int(meta.attrs["gene_index"])
        n_perm       = int(meta.attrs["n_permutations"])
        skipped      = bool(meta.attrs.get("skipped", False))
        metrics_str  = str(meta.attrs.get("metrics_stored", ",".join(BASE_METRICS)))
        stored_metrics = [m.strip() for m in metrics_str.split(",")]

        observed = {}
        null_dists = {}

        if not skipped and n_perm > 0:
            obs_grp = h5["observed"]
            for m in stored_metrics:
                observed[m] = float(obs_grp.attrs.get(m, float("nan")))

            null_grp = h5["null"]
            for m in stored_metrics:
                if m in null_grp:
                    null_dists[m] = null_grp[m][:]
                else:
                    null_dists[m] = np.array([], dtype=np.float32)

    return {
        "gene_id":       gene_id,
        "gene_index":    gene_index,
        "n_permutations": n_perm,
        "skipped":       skipped,
        "stored_metrics": stored_metrics,
        "observed":      observed,
        "null":          null_dists,
    }


def compute_pvals_for_gene(gene_data: dict, alpha: float = 0.05) -> dict:
    """
    Compute p-values and null summaries for all metrics of one gene.

    Returns a flat dict suitable for one row of the output TSV.
    """
    row = {
        "gene_id":       gene_data["gene_id"],
        "gene_idx":      gene_data["gene_index"],
        "n_permutations": gene_data["n_permutations"],
        "skipped":       gene_data["skipped"],
    }

    all_metrics = gene_data["stored_metrics"]

    if gene_data["skipped"] or gene_data["n_permutations"] == 0:
        # Fill NAs
        for m in BASE_METRICS + L2_METRICS:
            row[f"obs_{m}"]           = float("nan")
            row[f"pval_{m}"]          = float("nan")
            row[f"null_mean_{m}"]     = float("nan")
            row[f"null_std_{m}"]      = float("nan")
            row[f"null_all_zero_{m}"] = False
        row["n_metrics_sig"] = 0
        row["any_sig"]       = False
        return row

    n_sig = 0
    for m in all_metrics:
        obs  = gene_data["observed"].get(m, float("nan"))
        null = gene_data["null"].get(m, np.array([]))

        row[f"obs_{m}"] = obs

        if len(null) == 0 or np.isnan(obs):
            row[f"pval_{m}"]          = float("nan")
            row[f"null_mean_{m}"]     = float("nan")
            row[f"null_std_{m}"]      = float("nan")
            row[f"null_all_zero_{m}"] = False
        else:
            pval = permutation_pval(obs, null)
            row[f"pval_{m}"]          = pval
            row[f"null_mean_{m}"]     = float(np.mean(null))
            row[f"null_std_{m}"]      = float(np.std(null))
            row[f"null_all_zero_{m}"] = bool(np.all(null == 0))
            if pval < alpha:
                n_sig += 1

    row["n_metrics_sig"] = n_sig
    row["any_sig"]       = n_sig > 0
    return row


def collect_permutation_pvals(
    null_dir: Path,
    observed_tsv: Path,
    out_pvals_tsv: Path,
    out_augmented_tsv: Path | None,
    alpha: float = 0.05,
) -> None:
    """
    Main collection function.

    Scans null_dir for per-gene HDF5 files, computes p-values, writes output.

    Parameters
    ----------
    null_dir : directory containing permutation_null/{gi}_{gene_id}.h5
    observed_tsv : path to rewiring_hubs.tsv (Stage 3b output)
    out_pvals_tsv : output permutation_pvals.tsv
    out_augmented_tsv : optional output rewiring_hubs_with_pvals.tsv
    alpha : significance threshold for n_metrics_sig / any_sig columns
    """
    null_dir = Path(null_dir)
    out_pvals_tsv = Path(out_pvals_tsv)
    out_pvals_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Discover all per-gene HDF5 files
    h5_files = sorted(null_dir.glob("????_*.h5"))
    # Exclude files that end with _null_dist.tsv.gz pattern (shouldn't happen, but safe)
    h5_files = [f for f in h5_files if f.suffix == ".h5"]

    if not h5_files:
        print(f"No permutation null HDF5 files found in {null_dir}")
        return

    print(f"Found {len(h5_files)} permutation null files.")

    # Process each gene
    rows = []
    for h5_path in h5_files:
        gene_data = load_gene_null(h5_path)
        if gene_data is None:
            print(f"  WARNING: Could not load {h5_path}")
            continue
        row = compute_pvals_for_gene(gene_data, alpha=alpha)
        rows.append(row)
        status = "SKIPPED" if row["skipped"] else f"{row['n_metrics_sig']} sig"
        print(f"  {gene_data['gene_id']:20s}  n_perm={gene_data['n_permutations']:4d}  "
              f"obs_n_diff={row.get('obs_n_diff_edges', 'NA')}  [{status}]")

    if not rows:
        print("No valid permutation results found.")
        return

    pvals_df = pd.DataFrame(rows)

    # Reorder columns: gene info first, then per-metric groups
    id_cols = ["gene_id", "gene_idx", "n_permutations", "skipped"]
    metric_cols = []
    all_m_names = BASE_METRICS + [m for m in L2_METRICS if f"obs_{L2_METRICS[0]}" in pvals_df.columns]
    for m in all_m_names:
        for prefix in ["obs_", "pval_", "null_mean_", "null_std_", "null_all_zero_"]:
            col = f"{prefix}{m}"
            if col in pvals_df.columns:
                metric_cols.append(col)
    summary_cols = ["n_metrics_sig", "any_sig"]
    ordered_cols = id_cols + metric_cols + summary_cols
    ordered_cols = [c for c in ordered_cols if c in pvals_df.columns]
    pvals_df = pvals_df[ordered_cols]

    pvals_df.to_csv(out_pvals_tsv, sep="\t", index=False, float_format="%.6g")
    print(f"\nWrote {len(pvals_df)} genes to {out_pvals_tsv}")

    # Summary statistics
    sig_genes = pvals_df[pvals_df["any_sig"] == True]
    print(f"Genes with >= 1 significant metric (p < {alpha}): {len(sig_genes)}")

    # ---- Augmented rewiring_hubs TSV ----
    if out_augmented_tsv is not None and observed_tsv is not None:
        out_augmented_tsv = Path(out_augmented_tsv)
        out_augmented_tsv.parent.mkdir(parents=True, exist_ok=True)

        if not Path(observed_tsv).exists():
            print(f"WARNING: {observed_tsv} not found — skipping augmented TSV.")
            return

        hubs_df = pd.read_csv(observed_tsv, sep="\t")

        # Determine join key: gene_name or gene_id column
        join_col_hubs = None
        for candidate in ["gene_name", "gene_id", "gene"]:
            if candidate in hubs_df.columns:
                join_col_hubs = candidate
                break

        if join_col_hubs is None:
            print(f"WARNING: Cannot find gene ID column in {observed_tsv}. Columns: {list(hubs_df.columns)}")
            return

        # Keep only pval_* columns from pvals_df to merge
        pval_cols = [c for c in pvals_df.columns if c.startswith("pval_") or c in ["n_metrics_sig", "any_sig"]]
        merge_df = pvals_df[["gene_id"] + pval_cols].rename(columns={"gene_id": join_col_hubs})

        augmented = hubs_df.merge(merge_df, on=join_col_hubs, how="left")
        augmented.to_csv(out_augmented_tsv, sep="\t", index=False, float_format="%.6g")
        print(f"Wrote augmented table to {out_augmented_tsv} "
              f"({len(augmented)} rows, added {len(pval_cols)} pval columns)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Stage 7b: Collect permutation test p-values across all focus genes."
    )
    p.add_argument("--null-dir",         required=True,
                   help="Directory containing per-gene permutation null HDF5 files")
    p.add_argument("--observed-tsv",     default=None,
                   help="Path to rewiring_hubs.tsv (Stage 3b output; used for augmented TSV)")
    p.add_argument("--out-pvals-tsv",    required=True,
                   help="Output permutation_pvals.tsv path")
    p.add_argument("--out-augmented-tsv", default=None,
                   help="Optional output rewiring_hubs_with_pvals.tsv path")
    p.add_argument("--alpha",            type=float, default=0.05,
                   help="Significance threshold for n_metrics_sig / any_sig columns (default: 0.05)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    collect_permutation_pvals(
        null_dir=Path(args.null_dir),
        observed_tsv=Path(args.observed_tsv) if args.observed_tsv else None,
        out_pvals_tsv=Path(args.out_pvals_tsv),
        out_augmented_tsv=Path(args.out_augmented_tsv) if args.out_augmented_tsv else None,
        alpha=args.alpha,
    )


if __name__ == "__main__":
    main()
