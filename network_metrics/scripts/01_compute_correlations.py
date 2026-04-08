#!/usr/bin/env python3
"""
Stage 1: Full-matrix Spearman correlations with single BH-FDR correction.

Computes all-vs-all Spearman correlations on the FULL expression matrix
(all samples), applies ONE global BH-FDR correction, and stores only the
significant edges in sparse HDF5 format.

Key differences from the differential pipeline (02a_calc_base_correlations.py):
  1. Uses ALL samples — no low/high subpopulation partitioning.
  2. Applies ONE global FDR correction (not 3 independent corrections).
  3. No bootstrap indices needed.
  4. Single output: (gene_i, gene_j, corr, qval) for significant edges only.

Pipeline position
-----------------
Stage 0  preprocess (00convert_expr_to_hdf5.py)    →  expression.h5
Stage 1  THIS SCRIPT                               →  corr_significant.h5
Stage 2  02_build_network.py                       →  network_edges.h5

Output HDF5 layout
------------------
    meta/
        n_genes, n_samples, n_tests, n_significant, fdr_alpha, corr_threshold
    edges/
        gene_i   (n_sig,) int32   — row index (always < gene_j, upper-triangle)
        gene_j   (n_sig,) int32   — col index (always > gene_i)
        corr     (n_sig,) float32 — Spearman r
        qval     (n_sig,) float32 — BH-FDR q-value
    gene_names   (n_genes,) str  — stored as dataset (not attribute) to avoid 64 KB HDF5 limit
    sample_names (n_samples,) str

Memory budget (18K genes × 100 samples)
-----------------------------------------
  z @ z.T intermediate  : 18K×18K×4B ≈ 1.2 GB  (freed after triu extraction)
  corr_triu             : 162M×4B = 648 MB
  pval_triu             : 162M×4B = 648 MB
  qval_triu             : 162M×4B = 648 MB
  Peak                  : ~3.8 GB  →  request mem_mb=8000 in Snakefile

Low-RAM option
--------------
  Use --chunk-size (e.g. 1000) to process the full matrix in row blocks.
  This avoids holding the 1.2 GB z@z.T matrix simultaneously, reducing
  peak RAM to ~2 GB, at the cost of ~2× runtime.
  Note: BH-FDR must be applied globally, so all 162M p-values still
  accumulate before correction regardless of chunking.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import rankdata, t as t_dist
from statsmodels.stats.multitest import multipletests


# =============================================================================
# Spearman computation (same algorithm as differential pipeline)
# =============================================================================

def _rank_zscore_rows(x: np.ndarray) -> np.ndarray:
    """Rank each row (average ties), then z-score (ddof=1)."""
    r = rankdata(x, axis=1, method="average").astype(np.float32)
    r -= r.mean(axis=1, keepdims=True)
    denom = r.std(axis=1, keepdims=True, ddof=1)
    denom = np.where(denom < np.finfo(np.float32).eps, np.nan, denom)
    return (r / denom).astype(np.float32)


def compute_spearman_triu(expr: np.ndarray) -> np.ndarray:
    """
    All-vs-all Spearman via rank→zscore→BLAS DGEMM.

    Parameters
    ----------
    expr : (n_genes, n_samples) float32

    Returns
    -------
    corr_triu : (n_genes*(n_genes-1)//2,) float32  — upper-triangle flat array
    """
    z = _rank_zscore_rows(expr)
    n_genes = expr.shape[0]
    n_samples = expr.shape[1]
    c_full = (z @ z.T) / float(n_samples - 1)
    ti, tj = np.triu_indices(n_genes, k=1)
    return c_full[ti, tj].astype(np.float32)


def compute_spearman_triu_chunked(expr: np.ndarray, chunk_size: int) -> np.ndarray:
    """
    Chunked variant: processes C rows at a time to reduce peak RAM.

    Avoids allocating the full n_genes×n_genes matrix simultaneously.
    Still requires accumulating all n_tests p-values for global FDR.
    """
    z = _rank_zscore_rows(expr)
    n_genes, n_samples = expr.shape
    n_tests = n_genes * (n_genes - 1) // 2
    corr_triu = np.empty(n_tests, dtype=np.float32)
    ptr = 0
    for start in range(0, n_genes, chunk_size):
        end = min(start + chunk_size, n_genes)
        block = (z[start:end] @ z.T) / float(n_samples - 1)  # (chunk, n_genes)
        for local_i, gi in enumerate(range(start, end)):
            n_entries = n_genes - gi - 1
            if n_entries <= 0:
                continue
            corr_triu[ptr: ptr + n_entries] = block[local_i, gi + 1:]
            ptr += n_entries
    return corr_triu


def corr_to_pvals(r: np.ndarray, n_samples: int) -> np.ndarray:
    """t-distribution p-values from Spearman r. df = n_samples - 2."""
    df = n_samples - 2
    r_clip = np.clip(r, -1.0, 1.0)
    eps = np.finfo(np.float32).eps
    denom = np.maximum(1.0 - r_clip * r_clip, eps)
    t_stat = r_clip * np.sqrt(df / denom)
    p = 2.0 * t_dist.sf(np.abs(t_stat), df=df)
    return p.astype(np.float32)


def flat_idx_to_pair(k: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert flat upper-triangle index to (row i, col j) gene pairs.

    Uses quadratic formula — no n×n matrix required.
    k = i*n - i*(i+1)//2 + j - i - 1  where i < j.
    """
    k = k.astype(np.float64)
    i = (n - 2 - np.floor(np.sqrt(-8.0 * k + 4.0 * n * (n - 1) - 7) / 2.0 - 0.5)).astype(np.int32)
    j = (k + i + 1 - n * (n - 1) // 2 + (n - i) * (n - i - 1) // 2).astype(np.int32)
    return i, j


# =============================================================================
# Main computation
# =============================================================================

def compute_full_correlations(
    expr: np.ndarray,
    gene_names: list[str] | None,
    sample_names: list[str] | None,
    out_h5_path: Path,
    fdr_alpha: float = 0.05,
    corr_threshold: float = 0.0001,
    chunk_size: int | None = None,
    compression_level: int = 6,
) -> None:
    """
    Compute all-vs-all Spearman + single BH-FDR; save significant edges.

    Parameters
    ----------
    expr              : (n_genes, n_samples) expression matrix
    gene_names        : gene name list (stored as HDF5 string dataset)
    sample_names      : sample name list (stored as HDF5 string dataset)
    out_h5_path       : output path
    fdr_alpha         : BH-FDR threshold
    corr_threshold    : minimum |r| to retain (post-FDR pre-filter)
    chunk_size        : row chunk size for low-RAM mode; None = single pass
    compression_level : gzip level 1-9
    """
    n_genes, n_samples = expr.shape
    n_tests = n_genes * (n_genes - 1) // 2
    print(f"Full-matrix Spearman: {n_genes:,} genes × {n_samples} samples → {n_tests:,} tests")
    print(f"  Float32 per array: {n_tests * 4 / 1e9:.2f} GB")

    # --- Step 1: Spearman correlations ---
    if chunk_size:
        print(f"  Chunked mode (chunk_size={chunk_size}) to reduce peak RAM...")
        corr_triu = compute_spearman_triu_chunked(expr, chunk_size)
    else:
        print("  Single-pass (z @ z.T)...")
        corr_triu = compute_spearman_triu(expr)

    # --- Step 2: p-values ---
    print("  Computing p-values (t-dist, df=n-2)...")
    pval_triu = corr_to_pvals(corr_triu, n_samples)

    # --- Step 3: Single global BH-FDR ---
    print(f"  BH-FDR correction (alpha={fdr_alpha}) on {n_tests:,} p-values...")
    _, qval_triu, _, _ = multipletests(pval_triu, alpha=fdr_alpha, method="fdr_bh")
    qval_triu = qval_triu.astype(np.float32)
    del pval_triu

    # --- Step 4: Identify significant edges ---
    sig_mask = (qval_triu < fdr_alpha) & (np.abs(corr_triu) >= corr_threshold)
    flat_sig = np.where(sig_mask)[0].astype(np.int64)
    n_sig = len(flat_sig)
    del sig_mask
    print(f"\n  Significant edges: {n_sig:,} ({100.0 * n_sig / n_tests:.3f}% of {n_tests:,} tests)")

    # --- Step 5: Convert flat indices to (gene_i, gene_j) pairs ---
    print("  Converting flat indices to (gene_i, gene_j) pairs...")
    gene_i, gene_j = flat_idx_to_pair(flat_sig, n_genes)
    corr_sig = corr_triu[flat_sig]
    qval_sig = qval_triu[flat_sig]
    del corr_triu, qval_triu, flat_sig

    # --- Step 6: Write output HDF5 ---
    print(f"\n  Writing {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)
    chunk_sz = max(1, min(1_000_000, n_sig))
    kw = dict(compression="gzip", compression_opts=compression_level)

    with h5py.File(out_h5_path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["n_genes"] = n_genes
        meta.attrs["n_samples"] = n_samples
        meta.attrs["n_tests"] = n_tests
        meta.attrs["n_significant"] = n_sig
        meta.attrs["fdr_alpha"] = fdr_alpha
        meta.attrs["corr_threshold"] = corr_threshold

        grp = h5.create_group("edges")
        if n_sig > 0:
            grp.create_dataset("gene_i", data=gene_i,    chunks=(chunk_sz,), **kw)
            grp.create_dataset("gene_j", data=gene_j,    chunks=(chunk_sz,), **kw)
            grp.create_dataset("corr",   data=corr_sig,  chunks=(chunk_sz,), **kw)
            grp.create_dataset("qval",   data=qval_sig,  chunks=(chunk_sz,), **kw)
        else:
            for name, dt in [("gene_i", np.int32), ("gene_j", np.int32),
                              ("corr", np.float32), ("qval", np.float32)]:
                grp.create_dataset(name, data=np.array([], dtype=dt), **kw)

        # Store names as datasets (not attributes) to avoid 64 KB HDF5 attribute limit
        if gene_names is not None:
            h5.create_dataset("gene_names",   data=gene_names,   dtype=h5py.string_dtype())
        if sample_names is not None:
            h5.create_dataset("sample_names", data=sample_names, dtype=h5py.string_dtype())

    print("Done.")


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="[Network Metrics Stage 1] Full-matrix Spearman + single BH-FDR."
    )
    p.add_argument("--expr-h5", type=str, default=None,
                   help="Path to expression.h5 (required unless --toy).")
    p.add_argument("--out-h5", type=str, required=True,
                   help="Output path for corr_significant.h5.")
    p.add_argument("--fdr-alpha", type=float, default=0.05,
                   help="BH-FDR significance threshold (default: 0.05).")
    p.add_argument("--corr-threshold", type=float, default=0.0001,
                   help="Minimum |r| to retain (default: 0.0001).")
    p.add_argument("--chunk-size", type=int, default=None,
                   help="Row chunk size for low-RAM mode (default: None = single pass).")
    p.add_argument("--compression-level", type=int, default=6,
                   help="Gzip compression level 1-9 (default: 6).")
    p.add_argument("--toy", action="store_true",
                   help="Generate toy data (20 genes, 80 samples) for testing.")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.toy:
        rng = np.random.default_rng(42)
        expr = rng.normal(size=(20, 80)).astype(np.float32)
        gene_names = [f"gene_{i:03d}" for i in range(20)]
        sample_names = [f"sample_{j:03d}" for j in range(80)]
        print("[TOY] 20 genes × 80 samples")
    else:
        if args.expr_h5 is None:
            raise SystemExit("ERROR: --expr-h5 is required unless --toy is used.")
        print(f"Loading {args.expr_h5}...")
        with h5py.File(args.expr_h5, "r") as f:
            expr = f["expr"][:]
            gene_names = (
                [x.decode() if isinstance(x, bytes) else x for x in f["gene_names"][:]]
                if "gene_names" in f else None
            )
            sample_names = (
                [x.decode() if isinstance(x, bytes) else x for x in f["sample_names"][:]]
                if "sample_names" in f else None
            )
        print(f"  {expr.shape[0]:,} genes × {expr.shape[1]} samples")

    compute_full_correlations(
        expr=expr,
        gene_names=gene_names,
        sample_names=sample_names,
        out_h5_path=Path(args.out_h5),
        fdr_alpha=args.fdr_alpha,
        corr_threshold=args.corr_threshold,
        chunk_size=args.chunk_size,
        compression_level=args.compression_level,
    )


if __name__ == "__main__":
    main()
