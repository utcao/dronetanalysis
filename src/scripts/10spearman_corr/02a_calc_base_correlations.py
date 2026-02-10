#!/usr/bin/env python3
"""
Stage 2a: Compute BASE correlations for low/high subpopulations.

This script computes correlations for ALL gene pairs using base (non-bootstrap)
sample indices, then performs three significance tests:
1. Individual significance for low correlations (t-test → p-value → q-value)
2. Individual significance for high correlations (t-test → p-value → q-value)
3. Differential significance: r_low vs r_high (Fisher's Z test → p-value → q-value)

Pipeline position
-----------------
Stage 1   01subset/01get_extreme_pop_bootstrap.py  →  bootstrap_indices.h5
Stage 2a  THIS SCRIPT                              →  base_correlations.h5
Stage 2b  02b_bootstrap_significant_edges.py       →  bootstrap_significant.h5
Stage 3   03analyze_differential_coexpression.py   →  differential_edges.h5

Output HDF5 layout
------------------
    meta/
        n_genes, n_samples, n_tests, k_low, k_high
    low/
        corr_triu      (n_tests,) - Spearman correlations
        pval_triu      (n_tests,) - p-values from t-test
        qval_triu      (n_tests,) - FDR-corrected q-values
    high/
        corr_triu      (n_tests,)
        pval_triu      (n_tests,)
        qval_triu      (n_tests,)
    diff/
        fisher_z       (n_tests,) - Fisher's Z statistic
        pval_triu      (n_tests,) - p-values
        qval_triu      (n_tests,) - FDR-corrected q-values
    significant/
        sig_low        (n_tests,) bool - low_qval < threshold
        sig_high       (n_tests,) bool - high_qval < threshold
        sig_individual (n_tests,) bool - sig_low | sig_high
        sig_differential (n_tests,) bool - diff_qval < threshold
        sig_edges      (n_tests,) bool - sig_individual & sig_differential
        indices        (n_sig_edges,) int - flat indices of sig_edges

Efficiency
----------
Single job, O(N² × M) compute for two subsets.
Much faster than computing bootstrap for all edges.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import rankdata, t as t_dist, norm
from statsmodels.stats.multitest import multipletests


def rank_zscore_rows(x: np.ndarray, ddof: int = 1) -> np.ndarray:
    """Rank each row, then z-score with the given ddof."""
    r = rankdata(x, axis=1, method="average").astype(np.float32)
    mean = r.mean(axis=1, keepdims=True)
    r_centered = r - mean
    denom = r.std(axis=1, keepdims=True, ddof=ddof)
    denom = np.where(denom < np.finfo(np.float32).eps, np.nan, denom)
    return (r_centered / denom).astype(np.float32)


def compute_spearman_triu(expr_subset: np.ndarray) -> np.ndarray:
    """
    Compute Spearman correlation upper triangle for expression subset.

    Parameters
    ----------
    expr_subset : (n_genes, n_samples_subset)

    Returns
    -------
    corr_triu : (n_tests,) float32
    """
    z = rank_zscore_rows(expr_subset)
    n_samples = expr_subset.shape[1]
    n_genes = expr_subset.shape[0]

    # Full correlation matrix
    c_full = (z @ z.T) / float(n_samples - 1)

    # Extract upper triangle
    triu_rows, triu_cols = np.triu_indices(n_genes, k=1)
    return c_full[triu_rows, triu_cols].astype(np.float32)


def corr_to_pvals(r: np.ndarray, n_samples: int) -> np.ndarray:
    """
    Convert correlation values to p-values using t-distribution.

    t = r * sqrt((n-2) / (1-r²)), df = n-2
    """
    df = n_samples - 2
    r = np.clip(r, -1.0, 1.0)

    eps = np.finfo(np.float32).eps
    denom = np.maximum(1.0 - r * r, eps)
    t_stat = r * np.sqrt(df / denom)

    p = 2.0 * t_dist.sf(np.abs(t_stat), df=df)
    return p.astype(np.float32)


def fisher_z_test(r1: np.ndarray, r2: np.ndarray, n1: int, n2: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Fisher's Z-transform test for comparing two correlations.

    Parameters
    ----------
    r1, r2 : correlation arrays
    n1, n2 : sample sizes

    Returns
    -------
    z_score : Fisher's Z statistic
    pval : two-tailed p-value
    """
    # Clip to avoid log(0) or log(inf)
    r1_clip = np.clip(r1, -0.9999, 0.9999)
    r2_clip = np.clip(r2, -0.9999, 0.9999)

    # Fisher's Z transform: z = 0.5 * ln((1+r)/(1-r))
    z1 = 0.5 * np.log((1 + r1_clip) / (1 - r1_clip))
    z2 = 0.5 * np.log((1 + r2_clip) / (1 - r2_clip))

    # Standard error of difference
    se_diff = np.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))

    # Z-score for difference
    z_score = (z2 - z1) / se_diff

    # Two-tailed p-value
    pval = 2 * (1 - norm.cdf(np.abs(z_score)))

    return z_score.astype(np.float32), pval.astype(np.float32)


def compute_base_correlations(
    expr: np.ndarray,
    indices_h5_path: Path,
    out_h5_path: Path,
    fdr_alpha: float = 0.05,
    gene_index: int = 0,
    gene_names: list = None,
) -> None:
    """
    Compute base correlations and significance tests.

    Parameters
    ----------
    expr : (n_genes, n_samples) expression matrix
    indices_h5_path : path to bootstrap_indices.h5
    out_h5_path : output path
    fdr_alpha : FDR threshold for significance
    gene_index : which gene's indices to use (default: 0, uses first gene)
    gene_names : list of gene name strings (propagated from expression.h5)
    """
    n_genes, n_samples_total = expr.shape
    n_tests = n_genes * (n_genes - 1) // 2

    print(f"Computing base correlations for {n_genes:,} genes ({n_tests:,} edges)...")

    # Load base indices for the specified gene
    with h5py.File(indices_h5_path, "r") as h5:
        low_indices = h5["indices/low"][gene_index]    # (k_low,)
        high_indices = h5["indices/high"][gene_index]  # (k_high,)

        # Get parameters
        low_frac = h5["meta"].attrs["low_frac"]
        high_frac = h5["meta"].attrs["high_frac"]

    k_low = len(low_indices)
    k_high = len(high_indices)

    print(f"  Low subset: {k_low} samples ({low_frac*100:.0f}%)")
    print(f"  High subset: {k_high} samples ({high_frac*100:.0f}%)")

    # Subset expression matrices
    expr_low = np.ascontiguousarray(expr[:, low_indices])
    expr_high = np.ascontiguousarray(expr[:, high_indices])

    # Compute correlations
    print("  Computing low correlations...")
    corr_low = compute_spearman_triu(expr_low)

    print("  Computing high correlations...")
    corr_high = compute_spearman_triu(expr_high)

    # Compute p-values for individual correlations
    print("  Computing p-values...")
    pval_low = corr_to_pvals(corr_low, k_low)
    pval_high = corr_to_pvals(corr_high, k_high)

    # Fisher's Z test for differential correlation
    print("  Computing Fisher's Z test...")
    fisher_z, pval_diff = fisher_z_test(corr_low, corr_high, k_low, k_high)

    # FDR correction
    print("  Applying FDR correction...")
    _, qval_low, _, _ = multipletests(pval_low, alpha=fdr_alpha, method='fdr_bh')
    _, qval_high, _, _ = multipletests(pval_high, alpha=fdr_alpha, method='fdr_bh')
    _, qval_diff, _, _ = multipletests(pval_diff, alpha=fdr_alpha, method='fdr_bh')

    qval_low = qval_low.astype(np.float32)
    qval_high = qval_high.astype(np.float32)
    qval_diff = qval_diff.astype(np.float32)

    # Compute significance masks
    sig_low = qval_low < fdr_alpha
    sig_high = qval_high < fdr_alpha
    sig_individual = sig_low | sig_high
    sig_differential = qval_diff < fdr_alpha
    sig_edges = sig_individual & sig_differential

    # Get indices of significant edges
    sig_indices = np.where(sig_edges)[0].astype(np.int32)

    print(f"\n  Significance summary (FDR < {fdr_alpha}):")
    print(f"    sig_low:          {sig_low.sum():,} ({100*sig_low.sum()/n_tests:.2f}%)")
    print(f"    sig_high:         {sig_high.sum():,} ({100*sig_high.sum()/n_tests:.2f}%)")
    print(f"    sig_individual:   {sig_individual.sum():,} ({100*sig_individual.sum()/n_tests:.2f}%)")
    print(f"    sig_differential: {sig_differential.sum():,} ({100*sig_differential.sum()/n_tests:.2f}%)")
    print(f"    sig_edges (both): {sig_edges.sum():,} ({100*sig_edges.sum()/n_tests:.2f}%)")

    # Save results
    print(f"\n  Saving to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)

    chunk_size = min(1_000_000, n_tests)

    with h5py.File(out_h5_path, "w") as h5:
        # Metadata
        meta = h5.create_group("meta")
        meta.attrs["n_genes"] = n_genes
        meta.attrs["n_samples_total"] = n_samples_total
        meta.attrs["n_tests"] = n_tests
        meta.attrs["k_low"] = k_low
        meta.attrs["k_high"] = k_high
        meta.attrs["low_frac"] = low_frac
        meta.attrs["high_frac"] = high_frac
        meta.attrs["fdr_alpha"] = fdr_alpha
        meta.attrs["gene_index_used"] = gene_index

        # Low correlations
        grp_low = h5.create_group("low")
        grp_low.create_dataset("corr_triu", data=corr_low, chunks=(chunk_size,),
                               compression="gzip", compression_opts=4)
        grp_low.create_dataset("pval_triu", data=pval_low, chunks=(chunk_size,),
                               compression="gzip", compression_opts=4)
        grp_low.create_dataset("qval_triu", data=qval_low, chunks=(chunk_size,),
                               compression="gzip", compression_opts=4)

        # High correlations
        grp_high = h5.create_group("high")
        grp_high.create_dataset("corr_triu", data=corr_high, chunks=(chunk_size,),
                                compression="gzip", compression_opts=4)
        grp_high.create_dataset("pval_triu", data=pval_high, chunks=(chunk_size,),
                                compression="gzip", compression_opts=4)
        grp_high.create_dataset("qval_triu", data=qval_high, chunks=(chunk_size,),
                                compression="gzip", compression_opts=4)

        # Differential test
        grp_diff = h5.create_group("diff")
        grp_diff.create_dataset("fisher_z", data=fisher_z, chunks=(chunk_size,),
                                compression="gzip", compression_opts=4)
        grp_diff.create_dataset("pval_triu", data=pval_diff, chunks=(chunk_size,),
                                compression="gzip", compression_opts=4)
        grp_diff.create_dataset("qval_triu", data=qval_diff, chunks=(chunk_size,),
                                compression="gzip", compression_opts=4)

        # Delta (for convenience)
        delta = corr_high - corr_low
        grp_diff.create_dataset("delta_triu", data=delta, chunks=(chunk_size,),
                                compression="gzip", compression_opts=4)

        # Significance masks
        grp_sig = h5.create_group("significant")
        grp_sig.create_dataset("sig_low", data=sig_low, compression="gzip")
        grp_sig.create_dataset("sig_high", data=sig_high, compression="gzip")
        grp_sig.create_dataset("sig_individual", data=sig_individual, compression="gzip")
        grp_sig.create_dataset("sig_differential", data=sig_differential, compression="gzip")
        grp_sig.create_dataset("sig_edges", data=sig_edges, compression="gzip")
        grp_sig.create_dataset("indices", data=sig_indices, compression="gzip")

        grp_sig.attrs["n_sig_low"] = int(sig_low.sum())
        grp_sig.attrs["n_sig_high"] = int(sig_high.sum())
        grp_sig.attrs["n_sig_individual"] = int(sig_individual.sum())
        grp_sig.attrs["n_sig_differential"] = int(sig_differential.sum())
        grp_sig.attrs["n_sig_edges"] = int(sig_edges.sum())

        # Gene names (propagated from expression.h5)
        if gene_names is not None:
            h5.create_dataset("gene_names", data=gene_names, dtype=h5py.string_dtype())

    print("Done!")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 2a: Compute base correlations and significance tests."
    )
    parser.add_argument(
        "--expr-h5",
        type=str,
        default=None,
        help="Path to expression HDF5 file (required unless --toy is used).",
    )
    parser.add_argument(
        "--indices-h5",
        type=str,
        required=True,
        help="Path to bootstrap_indices.h5 from Stage 1.",
    )
    parser.add_argument(
        "--out-h5",
        type=str,
        default=None,
        help="Output HDF5 path (single-file mode).",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=None,
        help="Output directory for per-gene mode. Filename: {index:04d}_{gene_id}.h5",
    )
    parser.add_argument(
        "--fdr-alpha",
        type=float,
        default=0.05,
        help="FDR threshold for significance (default: 0.05).",
    )
    parser.add_argument(
        "--gene-index",
        type=int,
        default=0,
        help="Which gene's indices to use for subsetting (default: 0).",
    )
    parser.add_argument(
        "--toy",
        action="store_true",
        help="Use toy data for testing.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.toy:
        # Toy data: 5 genes, 50 samples (matches Stage 1 --toy)
        expr = np.random.default_rng(1).normal(size=(5, 50)).astype(np.float32)
        gene_names = None
        indices_h5_path = Path(args.indices_h5)
    else:
        if args.expr_h5 is None:
            raise SystemExit("ERROR: --expr-h5 is required when not using --toy.")
        # Load expression from HDF5
        print(f"Loading expression from {args.expr_h5}...")
        with h5py.File(args.expr_h5, "r") as f:
            expr = f["expr"][:]
            if "gene_names" in f:
                gene_names = [x.decode() if isinstance(x, bytes) else x for x in f["gene_names"][:]]
            else:
                gene_names = None
        print(f"  Shape: {expr.shape[0]} genes × {expr.shape[1]} samples")
        if gene_names:
            print(f"  Gene names: {gene_names[0]} ... {gene_names[-1]}")
        indices_h5_path = Path(args.indices_h5)

    # Resolve output path
    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        gene_id = gene_names[args.gene_index] if gene_names else f"gene_{args.gene_index}"
        out_h5_path = out_dir / f"{args.gene_index:04d}_{gene_id}.h5"
    elif args.out_h5:
        out_h5_path = Path(args.out_h5)
    else:
        raise SystemExit("ERROR: provide --out-h5 or --out-dir.")

    compute_base_correlations(
        expr=expr,
        indices_h5_path=indices_h5_path,
        out_h5_path=out_h5_path,
        fdr_alpha=args.fdr_alpha,
        gene_index=args.gene_index,
        gene_names=gene_names,
    )


if __name__ == "__main__":
    main()
