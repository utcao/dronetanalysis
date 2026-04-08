#!/usr/bin/env python3
"""
Stage 7: Permutation Test for Differential Co-expression Metrics.

For a given focus gene, this script tests whether the observed differential
co-expression metrics (rewiring, edge counts, effect sizes) exceed what is
expected under random sample assignment — i.e., whether the observed
high/low differences are driven by the true expression gradient rather than
chance.

Method
------
For each permutation replicate:
  1. Draw k_low + k_high samples uniformly at random from ALL samples
     (not restricted to the extreme-pool used in the real run).
  2. Split randomly into pseudo-low (k_low) and pseudo-high (k_high) groups.
  3. Rank-normalize ALL N genes within each pseudo-group (correct for Spearman).
  4. Extract only the focus gene's row of correlations (O(N × k) not O(N² × k)).
  5. Apply Fisher's Z-test and BH FDR correction on the N-1 focus-gene pairs.
  6. Compute the same metrics as Stage 3b: n_diff_edges, focus_deg_low/high,
     mean/max_abs_delta, L1_n_nodes, L1_rewire, L1_frac_rewire.
  7. Optionally compute L1×L1 sub-matrix for L2 expansion metrics.

The null distribution is stored alongside the observed values loaded from the
existing Stage 3 network file.

FDR Note
--------
The real pipeline applies BH FDR across all N(N-1)/2 ≈ 38M tests (full matrix).
Here we apply BH on N-1 ≈ 8762 tests (focus-gene row only). With fdr_alpha=0.7
used in practice, both procedures are near-saturated and the difference is
negligible. For stricter alphas this asymmetry could modestly inflate null
n_diff_edges — in that case the test is conservative (harder to reach
significance), which is a safe direction. The asymmetry is intentional and
documented here.

Pipeline position
-----------------
Stage 3   03_reconstruct_diff_network.py     →  networks/{gi}_{gene_id}.h5
Stage 7   THIS SCRIPT (per focus gene)       →  permutation_null/{gi}_{gene_id}.h5
                                                 permutation_null/{gi}_{gene_id}_null_dist.tsv.gz
Stage 7b  08_collect_permutation_pvals.py    →  permutation_pvals.tsv
                                                 rewiring_hubs_with_pvals.tsv

Usage
-----
python 07_permutation_test.py \\
    --expr-h5 results/expression.h5 \\
    --indices-h5 results/bootstrap_indices.h5 \\
    --network-h5 results/networks/0042_FBgn0001234.h5 \\
    --out-h5 results/permutation_null/0042_FBgn0001234.h5 \\
    --out-tsv results/permutation_null/0042_FBgn0001234_null_dist.tsv.gz \\
    --gene-index 42 \\
    --n-permutations 100 \\
    --fdr-alpha 0.7 \\
    --corr-threshold 0.0001 \\
    --seed 9999
"""

from __future__ import annotations

import argparse
import gzip
import io
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import rankdata, t as t_dist, norm
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------------------
# Statistical utilities (copied from 02a_calc_base_correlations.py)
# ---------------------------------------------------------------------------

def rank_zscore_rows(x: np.ndarray, ddof: int = 1) -> np.ndarray:
    """Rank each row, then z-score with the given ddof."""
    r = rankdata(x, axis=1, method="average").astype(np.float32)
    mean = r.mean(axis=1, keepdims=True)
    r_centered = r - mean
    denom = r.std(axis=1, keepdims=True, ddof=ddof)
    denom = np.where(denom < np.finfo(np.float32).eps, np.nan, denom)
    return (r_centered / denom).astype(np.float32)


def corr_to_pvals(r: np.ndarray, n_samples: int) -> np.ndarray:
    """Convert correlation values to p-values using t-distribution."""
    df = n_samples - 2
    r = np.clip(r, -1.0, 1.0)
    eps = np.finfo(np.float32).eps
    denom = np.maximum(1.0 - r * r, eps)
    t_stat = r * np.sqrt(df / denom)
    p = 2.0 * t_dist.sf(np.abs(t_stat), df=df)
    return p.astype(np.float32)


def fisher_z_test(
    r1: np.ndarray, r2: np.ndarray, n1: int, n2: int
) -> tuple[np.ndarray, np.ndarray]:
    """Fisher's Z-transform test for comparing two correlations."""
    r1_clip = np.clip(r1, -0.9999, 0.9999)
    r2_clip = np.clip(r2, -0.9999, 0.9999)
    z1 = 0.5 * np.log((1 + r1_clip) / (1 - r1_clip))
    z2 = 0.5 * np.log((1 + r2_clip) / (1 - r2_clip))
    se_diff = np.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    z_score = (z2 - z1) / se_diff
    pval = 2 * (1 - norm.cdf(np.abs(z_score)))
    return z_score.astype(np.float32), pval.astype(np.float32)


# ---------------------------------------------------------------------------
# Qualitative classification (inlined from 03_reconstruct_diff_network.py)
# ---------------------------------------------------------------------------

QUAL_DISAPPEAR   = 1
QUAL_NEW         = 2
QUAL_SIGN_CHANGE = 3
QUAL_STRENGTHEN  = 4
QUAL_WEAKEN      = 5


def classify_qualitative_counts(
    r_low: np.ndarray,
    r_high: np.ndarray,
    sig_low: np.ndarray,
    sig_high: np.ndarray,
    corr_threshold: float = 0.0001,
) -> dict:
    """
    Classify qualitative changes and return counts only (no labels needed here).

    Parameters
    ----------
    r_low, r_high : (n_edges,) correlation arrays for significant edges
    sig_low, sig_high : (n_edges,) boolean significance masks
    corr_threshold : minimum |r| to consider edge "present"

    Returns
    -------
    dict: n_disappear, n_new, n_sign_change, n_strengthen, n_weaken
    """
    present_low  = sig_low  & (np.abs(r_low)  >= corr_threshold)
    present_high = sig_high & (np.abs(r_high) >= corr_threshold)

    sign_low  = np.sign(r_low)
    sign_high = np.sign(r_high)

    mask_disappear  = present_low & ~present_high
    mask_new        = ~present_low & present_high
    mask_both       = present_low & present_high
    mask_sign_chg   = mask_both & (sign_low != sign_high)
    abs_low, abs_high = np.abs(r_low), np.abs(r_high)
    mask_strengthen = mask_both & ~mask_sign_chg & (abs_high > abs_low)
    mask_weaken     = mask_both & ~mask_sign_chg & (abs_high < abs_low)

    return {
        "n_disappear":   int(mask_disappear.sum()),
        "n_new":         int(mask_new.sum()),
        "n_sign_change": int(mask_sign_chg.sum()),
        "n_strengthen":  int(mask_strengthen.sum()),
        "n_weaken":      int(mask_weaken.sum()),
    }


# ---------------------------------------------------------------------------
# Per-permutation computation
# ---------------------------------------------------------------------------

def compute_focus_row_correlations(
    expr: np.ndarray,
    sample_indices: np.ndarray,
    focus_gene: int,
) -> np.ndarray:
    """
    Compute Spearman correlations between focus_gene and ALL other genes
    using only the given sample subset.

    All N genes are rank-normalized within the subset (required for correct
    Spearman). Only the focus gene's row is extracted — O(N × k), not O(N² × k).

    Parameters
    ----------
    expr : (n_genes, n_samples) full expression matrix
    sample_indices : (k,) indices into the sample axis
    focus_gene : row index of the focus gene

    Returns
    -------
    (n_genes,) float32 Spearman r; position focus_gene = 1.0 (self-corr)
    """
    k = len(sample_indices)
    z = rank_zscore_rows(np.ascontiguousarray(expr[:, sample_indices]))  # (N, k)
    row_corr = (z[focus_gene] @ z.T) / float(k - 1)    # (N,)
    row_corr = np.clip(row_corr, -1.0, 1.0)
    return row_corr.astype(np.float32)


def run_single_permutation(
    expr: np.ndarray,
    k_low: int,
    k_high: int,
    focus_gene: int,
    fdr_alpha: float,
    corr_threshold: float,
    rng: np.random.Generator,
    compute_l2: bool = False,
) -> dict:
    """
    Run one permutation replicate and return focus-gene metrics.

    Parameters
    ----------
    expr : (n_genes, n_samples)
    k_low, k_high : group sizes (match the real run)
    focus_gene : index of the focus gene
    fdr_alpha : BH FDR threshold
    corr_threshold : minimum |r| to consider an edge "present"
    rng : seeded random generator
    compute_l2 : whether to compute L1×L1 sub-matrix for L2 metrics

    Returns
    -------
    dict of scalar metrics for this replicate
    """
    n_samples = expr.shape[1]

    # --- 1. Random group assignment from ALL samples ---
    perm_indices = rng.choice(n_samples, size=k_low + k_high, replace=False)
    pseudo_low  = perm_indices[:k_low]
    pseudo_high = perm_indices[k_low:]

    # --- 2. Focus-gene row correlations (all N genes ranked, one row extracted) ---
    r_low  = compute_focus_row_correlations(expr, pseudo_low,  focus_gene)
    r_high = compute_focus_row_correlations(expr, pseudo_high, focus_gene)

    # Remove self-correlation
    mask = np.ones(len(r_low), dtype=bool)
    mask[focus_gene] = False
    r_low_partner  = r_low[mask]
    r_high_partner = r_high[mask]

    # --- 3. Significance tests on N-1 pairs ---
    pval_low  = corr_to_pvals(r_low_partner,  k_low)
    pval_high = corr_to_pvals(r_high_partner, k_high)
    _, pval_diff = fisher_z_test(r_low_partner, r_high_partner, k_low, k_high)

    _, qval_low,  _, _ = multipletests(pval_low,  alpha=fdr_alpha, method="fdr_bh")
    _, qval_high, _, _ = multipletests(pval_high, alpha=fdr_alpha, method="fdr_bh")
    _, qval_diff, _, _ = multipletests(pval_diff, alpha=fdr_alpha, method="fdr_bh")

    sig_low_p  = (qval_low  < fdr_alpha).astype(bool)
    sig_high_p = (qval_high < fdr_alpha).astype(bool)
    sig_diff_p = (qval_diff < fdr_alpha).astype(bool)
    sig_edges  = (sig_low_p | sig_high_p) & sig_diff_p

    n_diff_edges = int(sig_edges.sum())

    # --- 4. Focus-gene degree in low/high condition networks ---
    # Degree = how many partners are "present" (sig AND |r| >= threshold)
    focus_deg_low  = int((sig_low_p  & (np.abs(r_low_partner)  >= corr_threshold)).sum())
    focus_deg_high = int((sig_high_p & (np.abs(r_high_partner) >= corr_threshold)).sum())

    # --- 5. Effect size metrics over significant edges ---
    if n_diff_edges > 0:
        delta_sig = r_high_partner[sig_edges] - r_low_partner[sig_edges]
        mean_abs_delta = float(np.mean(np.abs(delta_sig)))
        max_abs_delta  = float(np.max(np.abs(delta_sig)))
    else:
        mean_abs_delta = 0.0
        max_abs_delta  = 0.0

    # --- 6. Qualitative classification for L1 rewiring ---
    if n_diff_edges > 0:
        qual_counts = classify_qualitative_counts(
            r_low_partner[sig_edges],
            r_high_partner[sig_edges],
            sig_low_p[sig_edges],
            sig_high_p[sig_edges],
            corr_threshold=corr_threshold,
        )
        L1_rewire = qual_counts["n_disappear"] + qual_counts["n_new"] + qual_counts["n_sign_change"]
    else:
        L1_rewire = 0

    L1_frac_rewire = float(L1_rewire / n_diff_edges) if n_diff_edges > 0 else 0.0

    result = {
        "n_diff_edges":    n_diff_edges,
        "focus_deg_low":   focus_deg_low,
        "focus_deg_high":  focus_deg_high,
        "mean_abs_delta":  mean_abs_delta,
        "max_abs_delta":   max_abs_delta,
        "L1_n_nodes":      n_diff_edges,    # L1 node count = direct partner count
        "L1_rewire":       L1_rewire,
        "L1_frac_rewire":  L1_frac_rewire,
    }

    # --- 7. Optional: L1×L1 sub-matrix for L2 expansion metrics ---
    if compute_l2 and n_diff_edges > 1:
        # partner_indices are positions in the masked (N-1) array
        partner_positions = np.where(sig_edges)[0]

        # Map back to original gene indices (skip focus_gene)
        all_gene_indices = np.delete(np.arange(expr.shape[0]), focus_gene)
        l1_gene_indices = all_gene_indices[partner_positions]

        # Sub-expression matrices for L1 genes only
        z_low_l1  = rank_zscore_rows(np.ascontiguousarray(expr[l1_gene_indices][:, pseudo_low]))
        z_high_l1 = rank_zscore_rows(np.ascontiguousarray(expr[l1_gene_indices][:, pseudo_high]))
        nl1 = len(l1_gene_indices)

        r_l1_low  = (z_low_l1  @ z_low_l1.T)  / float(k_low  - 1)
        r_l1_high = (z_high_l1 @ z_high_l1.T) / float(k_high - 1)

        # Upper triangle of L1×L1 block
        triu_r, triu_c = np.triu_indices(nl1, k=1)
        r_l1l1_low  = r_l1_low[triu_r, triu_c]
        r_l1l1_high = r_l1_high[triu_r, triu_c]

        pval_l1_low  = corr_to_pvals(r_l1l1_low,  k_low)
        pval_l1_high = corr_to_pvals(r_l1l1_high, k_high)
        _, pval_l1_diff = fisher_z_test(r_l1l1_low, r_l1l1_high, k_low, k_high)

        if len(pval_l1_low) > 1:
            _, qval_l1_low,  _, _ = multipletests(pval_l1_low,  alpha=fdr_alpha, method="fdr_bh")
            _, qval_l1_high, _, _ = multipletests(pval_l1_high, alpha=fdr_alpha, method="fdr_bh")
            _, qval_l1_diff, _, _ = multipletests(pval_l1_diff, alpha=fdr_alpha, method="fdr_bh")
            sig_l1_low  = qval_l1_low  < fdr_alpha
            sig_l1_high = qval_l1_high < fdr_alpha
            sig_l1_diff = (sig_l1_low | sig_l1_high) & (qval_l1_diff < fdr_alpha)
            n_l1l1_diff = int(sig_l1_diff.sum())
        else:
            n_l1l1_diff = 0
            sig_l1_diff = np.zeros(len(pval_l1_low), dtype=bool)

        # L2 nodes: genes connected to any L1 partner but not the focus gene or other L1s
        # (approximation in permutation: count unique partners of L1 genes in L1×L1 sig edges)
        # Use L2L1_deg = L2_n_nodes / L1_n_nodes as ratio metric
        L2_n_nodes_approx = max(n_l1l1_diff, 0)   # edges in L1 block ≈ connectivity measure
        L2L1_deg    = float(L2_n_nodes_approx / n_diff_edges) if n_diff_edges > 0 else 0.0
        L2L1_rewire = 0.0  # not computed in this approximation

        result["L2L1_deg"]    = L2L1_deg
        result["L2L1_rewire"] = L2L1_rewire

    return result


# ---------------------------------------------------------------------------
# Main permutation loop
# ---------------------------------------------------------------------------

def load_observed_values(network_h5: Path) -> dict:
    """
    Load observed focus-gene metrics from an existing Stage 3 network file.

    Parameters
    ----------
    network_h5 : path to networks/{gi}_{gene_id}.h5

    Returns
    -------
    dict of observed metric values, or {'skipped': True} if no significant edges
    """
    with h5py.File(network_h5, "r") as h5:
        n_sig = int(h5["meta"].attrs["n_significant"])
        k_low  = int(h5["meta"].attrs["k_low"])
        k_high = int(h5["meta"].attrs["k_high"])
        focus_gene_index = int(h5["meta"].attrs.get("focus_gene_index",
                                                     h5["meta"].attrs.get("gene_index_used", -1)))

        if n_sig == 0:
            return {
                "skipped": True,
                "k_low": k_low,
                "k_high": k_high,
                "focus_gene_index": focus_gene_index,
            }

        # Focus gene metrics stored as attrs under focus_gene/metrics/
        metrics = {}
        if "focus_gene/metrics" in h5:
            for key, val in h5["focus_gene/metrics"].attrs.items():
                metrics[key] = val

        # Effect size from edges
        if "edges/delta_base" in h5:
            delta = h5["edges/delta_base"][:]
            obs_mean_abs_delta = float(np.mean(np.abs(delta))) if len(delta) > 0 else 0.0
            obs_max_abs_delta  = float(np.max(np.abs(delta)))  if len(delta) > 0 else 0.0
        else:
            obs_mean_abs_delta = 0.0
            obs_max_abs_delta  = 0.0

        return {
            "skipped": False,
            "k_low":  k_low,
            "k_high": k_high,
            "focus_gene_index": focus_gene_index,
            "n_diff_edges":   n_sig,
            "focus_deg_low":  int(metrics.get("focus_deg_low",  0)),
            "focus_deg_high": int(metrics.get("focus_deg_high", 0)),
            "mean_abs_delta": obs_mean_abs_delta,
            "max_abs_delta":  obs_max_abs_delta,
            "L1_n_nodes":     int(metrics.get("L1_n_nodes",     0)),
            "L1_rewire":      int(metrics.get("L1_rewire",      0)),
            "L1_frac_rewire": float(metrics.get("L1_frac_rewire", 0.0)),
            "L2L1_deg":       float(metrics.get("L2L1_deg",     0.0)),
            "L2L1_rewire":    float(metrics.get("L2L1_rewire",  0.0)),
        }


def run_permutation_test(
    expr: np.ndarray,
    observed: dict,
    out_h5_path: Path,
    out_tsv_path: Path,
    gene_id: str,
    gene_index: int,
    n_permutations: int = 100,
    fdr_alpha: float = 0.7,
    corr_threshold: float = 0.0001,
    seed: int = 9999,
    compute_l2: bool = False,
) -> None:
    """
    Run the permutation test for one focus gene and write outputs.

    Parameters
    ----------
    expr : (n_genes, n_samples) full expression matrix
    observed : dict returned by load_observed_values()
    out_h5_path : output HDF5 path
    out_tsv_path : output TSV (gzipped) path for R plotting
    gene_id : gene name string (e.g. "FBgn0001234")
    gene_index : index of this gene in the expression matrix
    n_permutations : number of permutation replicates
    fdr_alpha : BH FDR threshold (should match the real pipeline run)
    corr_threshold : minimum |r| to consider edge "present"
    seed : base seed; actual seed for this gene = seed + gene_index
    compute_l2 : whether to compute L1×L1 sub-matrix for L2 metrics
    """
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)

    k_low  = observed["k_low"]
    k_high = observed["k_high"]
    focus_gene = gene_index

    # Determine which metrics to store
    base_metrics = [
        "n_diff_edges", "focus_deg_low", "focus_deg_high",
        "mean_abs_delta", "max_abs_delta",
        "L1_n_nodes", "L1_rewire", "L1_frac_rewire",
    ]
    l2_metrics = ["L2L1_deg", "L2L1_rewire"] if compute_l2 else []
    all_metrics = base_metrics + l2_metrics

    if observed.get("skipped", False):
        print(f"  Gene {gene_id} has 0 significant edges in Stage 3. Writing skipped HDF5.")
        with h5py.File(out_h5_path, "w") as h5:
            meta = h5.create_group("meta")
            meta.attrs["gene_index"]     = gene_index
            meta.attrs["gene_id"]        = gene_id
            meta.attrs["n_permutations"] = 0
            meta.attrs["skipped"]        = True
            meta.attrs["k_low"]          = k_low
            meta.attrs["k_high"]         = k_high
        # Write empty TSV
        with gzip.open(out_tsv_path, "wt") as f:
            f.write("\t".join(all_metrics) + "\n")
        return

    rng = np.random.default_rng(seed + gene_index)

    print(f"  Running {n_permutations} permutations for {gene_id} "
          f"(k_low={k_low}, k_high={k_high}, n_genes={expr.shape[0]}, "
          f"n_samples={expr.shape[1]})...")

    # Collect null distributions
    null: dict[str, list] = {m: [] for m in all_metrics}

    for p in range(n_permutations):
        res = run_single_permutation(
            expr=expr,
            k_low=k_low,
            k_high=k_high,
            focus_gene=focus_gene,
            fdr_alpha=fdr_alpha,
            corr_threshold=corr_threshold,
            rng=rng,
            compute_l2=compute_l2,
        )
        for m in base_metrics:
            null[m].append(res[m])
        for m in l2_metrics:
            null[m].append(res.get(m, 0.0))

        if (p + 1) % 10 == 0:
            print(f"    {p + 1}/{n_permutations} done", flush=True)

    # ---- Write HDF5 ----
    with h5py.File(out_h5_path, "w") as h5:
        # Metadata
        meta = h5.create_group("meta")
        meta.attrs["gene_index"]     = gene_index
        meta.attrs["gene_id"]        = gene_id
        meta.attrs["n_permutations"] = n_permutations
        meta.attrs["n_genes"]        = expr.shape[0]
        meta.attrs["n_samples"]      = expr.shape[1]
        meta.attrs["k_low"]          = k_low
        meta.attrs["k_high"]         = k_high
        meta.attrs["fdr_alpha"]      = fdr_alpha
        meta.attrs["corr_threshold"] = corr_threshold
        meta.attrs["seed"]           = seed
        meta.attrs["compute_l2"]     = compute_l2
        meta.attrs["skipped"]        = False
        meta.attrs["metrics_stored"] = ",".join(all_metrics)

        # Observed values
        obs_grp = h5.create_group("observed")
        for m in all_metrics:
            obs_grp.attrs[m] = observed.get(m, 0.0)

        # Null distributions
        null_grp = h5.create_group("null")
        chunk = min(n_permutations, 500)
        for m in all_metrics:
            arr = np.array(null[m], dtype=np.float32)
            null_grp.create_dataset(
                m, data=arr,
                chunks=(chunk,),
                compression="gzip", compression_opts=6,
            )

    # ---- Write TSV (for R plotting) ----
    rows = []
    header = ["gene_id", "perm_idx"] + all_metrics
    for p in range(n_permutations):
        row = [gene_id, str(p)] + [str(null[m][p]) for m in all_metrics]
        rows.append("\t".join(row))

    with gzip.open(out_tsv_path, "wt") as f:
        f.write("\t".join(header) + "\n")
        f.write("\n".join(rows) + "\n")

    print(f"  Done. Null distributions written to {out_h5_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Stage 7: Permutation test for differential co-expression metrics."
    )
    p.add_argument("--expr-h5",     required=True,  help="Path to expression.h5")
    p.add_argument("--indices-h5",  required=True,  help="Path to bootstrap_indices.h5")
    p.add_argument("--network-h5",  required=True,  help="Path to networks/{gi}_{gene_id}.h5 (Stage 3 output)")
    p.add_argument("--out-h5",      required=True,  help="Output HDF5 path")
    p.add_argument("--out-tsv",     required=True,  help="Output null-dist TSV (.tsv.gz) path for R plotting")
    p.add_argument("--gene-index",  required=True, type=int, help="Index of this gene in expression matrix")
    p.add_argument("--n-permutations", type=int, default=100, help="Number of permutation replicates (default: 100)")
    p.add_argument("--fdr-alpha",   type=float, default=0.7,    help="BH FDR threshold (default: 0.7)")
    p.add_argument("--corr-threshold", type=float, default=0.0001, help="Minimum |r| to consider edge present (default: 0.0001)")
    p.add_argument("--seed",        type=int, default=9999,  help="Base random seed; actual seed = seed + gene_index")
    p.add_argument("--compute-l2",  action="store_true",     help="Compute L1×L1 sub-matrix for L2 expansion metrics (slower)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Load expression matrix
    print(f"Loading expression matrix from {args.expr_h5}...")
    with h5py.File(args.expr_h5, "r") as h5:
        expr = h5["expr"][:]  # (n_genes, n_samples)
        if "gene_names" in h5:
            gene_names = [
                x.decode() if isinstance(x, bytes) else x
                for x in h5["gene_names"][:]
            ]
        else:
            gene_names = [f"gene_{i}" for i in range(expr.shape[0])]
    print(f"  Expression matrix: {expr.shape[0]:,} genes × {expr.shape[1]:,} samples")

    gene_id = gene_names[args.gene_index]
    print(f"  Focus gene: [{args.gene_index}] {gene_id}")

    # Load observed values from Stage 3 output
    print(f"Loading observed values from {args.network_h5}...")
    observed = load_observed_values(Path(args.network_h5))
    observed["focus_gene_index"] = args.gene_index  # ensure correct index

    # Run permutation test
    run_permutation_test(
        expr=expr,
        observed=observed,
        out_h5_path=Path(args.out_h5),
        out_tsv_path=Path(args.out_tsv),
        gene_id=gene_id,
        gene_index=args.gene_index,
        n_permutations=args.n_permutations,
        fdr_alpha=args.fdr_alpha,
        corr_threshold=args.corr_threshold,
        seed=args.seed,
        compute_l2=args.compute_l2,
    )


if __name__ == "__main__":
    main()
