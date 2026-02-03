#!/usr/bin/env python3
"""
Per-gene SGE array task: read bootstrap indices → subset expression → compute
Spearman upper-triangle correlations for every replicate.

Pipeline position
-----------------
Stage 1  01subset/01subset_bootstrap.py   →  bootstrap_indices.h5   (one job)
Stage 2  THIS SCRIPT                      →  corr/gene_XXXXX.h5     (one job per gene)

I/O per gene
------------
Input   bootstrap_indices.h5          (shared, read-only)
            indices/low          (n_genes, k_low)
            indices/high         (n_genes, k_high)
            indices/low_boot     (n_genes, n_bootstrap, k_resample_low)
            indices/high_boot    (n_genes, n_bootstrap, k_resample_high)

Output  corr/gene_XXXXX.h5      (one file per gene)
            meta/                    params + gene label
            low/base/corr_triu       base-partition Spearman triu  (n_tests,)
            low/boot/corr_triu       (n_bootstrap, n_tests)
            high/base/corr_triu      base-partition Spearman triu  (n_tests,)
            high/boot/corr_triu      (n_bootstrap, n_tests)

Design choices
--------------
* rank_zscore_rows is copied locally so that the digit-prefixed
  01cal_corr_edge.py is never imported (importlib dance, heavy statsmodels).
* Bootstrap subsets are small (k_resample ~ 0.8 * 0.2 * N_samples), so a single
  matmul without block-tiling stays comfortably in L2/L3 cache.
* corr_triu is extracted via np.triu_indices, which produces the same row-major
  upper-triangle order as TriuIndex in 01cal_corr_edge.py.
* --toy uses the identical 5 × 50 RNG seed as 01subset_bootstrap.py --toy so
  that the two stages are compatible for end-to-end smoke tests.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import rankdata


# ---------------------------------------------------------------------------
# Local rank-zscore (mirrors 01cal_corr_edge.py:rank_zscore_rows exactly)
# ---------------------------------------------------------------------------
def rank_zscore_rows(x: np.ndarray, ddof: int = 1) -> np.ndarray:
    """Rank each row, then z-score with the given ddof."""
    r = rankdata(x, axis=1, method="average").astype(np.float32)
    mean = r.mean(axis=1, keepdims=True)
    r_centered = r - mean
    denom = r.std(axis=1, keepdims=True, ddof=ddof)
    denom = np.where(denom < np.finfo(np.float32).eps, np.nan, denom)
    return (r_centered / denom).astype(np.float32)


# ---------------------------------------------------------------------------
# Core: rank → zscore → matmul → extract upper triangle
# ---------------------------------------------------------------------------
def compute_corr_triu(x_subset: np.ndarray, triu_rows: np.ndarray, triu_cols: np.ndarray) -> np.ndarray:
    """Full Spearman triu for a small (n_genes_subset, n_samples) matrix.

    Returns
    -------
    corr_triu : float32 (n_tests,) in row-major upper-triangle order. It means values are stored row by row and only for i < j,
           e.g. starting from the first gene pair (0,1), then (0,2), ..., (0,N-1), then (1,2), ..., (1,N-1), and so on.
    """
    z = rank_zscore_rows(x_subset)
    n_samples = x_subset.shape[1]
    c_full = (z @ z.T) / float(n_samples - 1)
    return c_full[triu_rows, triu_cols].astype(np.float32) # use indices to get pairs


# ---------------------------------------------------------------------------
# Per-gene orchestration
# ---------------------------------------------------------------------------
def run_gene(
    gene_id: int,
    expr: np.ndarray,  # full expression matrix (n_genes_total, n_samples)
    indices_h5_path: Path,
    out_dir: Path,
) -> None:
    """Read indices for *gene_id*, subset, compute, write one HDF5."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"gene_{gene_id:05d}.h5"

    # --- read indices (one chunk per dataset thanks to chunks=(1, …)) ---
    with h5py.File(indices_h5_path, "r") as idx_h5:
        low_base = idx_h5["indices/low"][gene_id]          # (k_low,)
        high_base = idx_h5["indices/high"][gene_id]        # (k_high,)
        low_boot = idx_h5["indices/low_boot"][gene_id]     # (n_bootstrap, k_resample_low)
        high_boot = idx_h5["indices/high_boot"][gene_id]   # (n_bootstrap, k_resample_high)

        # copy bootstrap params from metadata
        meta_src = idx_h5["meta"]
        meta_attrs = {k: meta_src.attrs[k] for k in meta_src.attrs}

    n_bootstrap = low_boot.shape[0]

    # --- column-subset helper ------------------------------------------------
    # expr rows are genes in the *full* dataset; the subset selects *samples*
    # (columns).  All genes participate in each subset's co-expression.
    n_genes_total = expr.shape[0]

    def subset_expr(sample_indices: np.ndarray) -> np.ndarray:
        """
        Return expr[:, sample_indices] as contiguous float32.
                ensures an array is stored in contiguous memory
        """
        return np.ascontiguousarray(expr[:, sample_indices])

    # --- pre-compute triu index arrays for the full gene set ------------------
    # np.triu_indices returns the row and column indices of the upper triangle of a square matrix,
    #                 excluding the diagonal.
    triu_rows, triu_cols = np.triu_indices(n_genes_total, k=1)
    n_tests = len(triu_rows)

    # --- write output HDF5 ---------------------------------------------------
    with h5py.File(out_path, "w") as h5:
        # metadata
        meta = h5.create_group("meta")
        meta.attrs["gene_id"] = gene_id
        for k, v in meta_attrs.items():
            meta.attrs[k] = v

        # helper to create one group with base + boot datasets
        def write_partition(
            group_name: str,  # "low" or "high"
            base_indices: np.ndarray,
            boot_indices: np.ndarray,
        ) -> None:
            grp = h5.create_group(group_name)

            # --- base correlation (single subset) ---
            base_corr = compute_corr_triu(subset_expr(base_indices), triu_rows, triu_cols)
            grp_base = grp.create_group("base")
            grp_base.create_dataset(
                "corr_triu",
                data=base_corr,
                chunks=(n_tests,),
                compression="gzip",
                compression_opts=4,
            )

            # --- bootstrap replicates ---
            grp_boot = grp.create_group("boot")
            # pre-allocate; chunked along replicate axis so each replicate is one chunk read
            ds = grp_boot.create_dataset(
                "corr_triu",
                shape=(n_bootstrap, n_tests),
                dtype=np.float32,
                chunks=(1, n_tests),
                compression="gzip",
                compression_opts=4,
            )
            for rep in range(n_bootstrap):
                ds[rep] = compute_corr_triu(
                    subset_expr(boot_indices[rep]), triu_rows, triu_cols
                )

        write_partition("low", low_base, low_boot)
        write_partition("high", high_base, high_boot)

    print(f"  gene {gene_id:05d} → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Per-gene bootstrap correlation (SGE array task)."
    )
    parser.add_argument(
        "--gene-id",
        type=int,
        default=None,
        help="Gene index (0-based).  Defaults to $SGE_TASK_ID if unset.",
    )
    parser.add_argument(
        "--indices-h5",
        type=str,
        default=None,
        help="Path to bootstrap_indices.h5 from Stage 1.",
    )
    parser.add_argument(
        "--expr-tsv",
        type=str,
        default=None,
        help="Input expression TSV (genes as rows, samples as columns).",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="results/corr",
        help="Output directory for per-gene HDF5 files (default: results/corr).",
    )
    parser.add_argument(
        "--toy",
        action="store_true",
        help="Use the same 5×50 toy matrix as 01subset_bootstrap.py --toy.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # --- resolve gene_id (CLI flag > SGE_TASK_ID) ---
    gene_id = args.gene_id
    if gene_id is None:
        sge_id = os.environ.get("SGE_TASK_ID")
        if sge_id is None:
            raise SystemExit(
                "Provide --gene-id or run inside an SGE job array (SGE_TASK_ID)."
            )
        gene_id = int(sge_id)

    # --- load expression matrix ---
    if args.toy:
        # MUST match 01subset_bootstrap.py --toy exactly: seed=1, shape (5, 50)
        expr = np.random.default_rng(1).normal(size=(5, 50)).astype(np.float32)
        indices_h5_path = Path(args.indices_h5) if args.indices_h5 else Path("results/bootstrap_indices.h5")
    else:
        if args.expr_tsv is None or args.indices_h5 is None:
            raise SystemExit("Provide --expr-tsv and --indices-h5, or use --toy.")
        import pandas as pd

        df = pd.read_csv(args.expr_tsv, sep="\t", index_col=0)
        expr = df.to_numpy(dtype=np.float32)
        indices_h5_path = Path(args.indices_h5)

    run_gene(
        gene_id=gene_id,
        expr=expr,
        indices_h5_path=indices_h5_path,
        out_dir=Path(args.out_dir),
    )


if __name__ == "__main__":
    main()
