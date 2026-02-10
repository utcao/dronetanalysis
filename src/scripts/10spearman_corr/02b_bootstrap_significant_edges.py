#!/usr/bin/env python3
"""
Stage 2b: Bootstrap ONLY significant edges.

Instead of computing bootstrap correlations for ALL ~40M edges, this script
computes bootstrap only for edges identified as significant in Stage 2a.
This provides ~100× speedup when ~1% of edges are significant.

Pipeline position
-----------------
Stage 1   01subset/01get_extreme_pop_bootstrap.py  →  bootstrap_indices.h5
Stage 2a  02a_calc_base_correlations.py            →  base_correlations.h5
Stage 2b  THIS SCRIPT                              →  bootstrap_significant.h5
Stage 3   03analyze_differential_coexpression.py   →  differential_edges.h5

Input
-----
    expression.h5           - expression matrix
    bootstrap_indices.h5    - bootstrap sample indices from Stage 1
    base_correlations.h5    - from Stage 2a, contains significant edge indices

Output HDF5 layout
------------------
    meta/
        n_sig_edges, n_bootstrap, edge_selection_mode
    edges/
        indices         (n_sig_edges,) - flat indices in triu array
        gene_i          (n_sig_edges,) - gene i for each edge; (gene_i, gene_j) pairs
        gene_j          (n_sig_edges,) - gene j for each edge
    base/
        delta           (n_sig_edges,) - original delta (high - low)
        r_low           (n_sig_edges,) - original low correlation
        r_high          (n_sig_edges,) - original high correlation
    boot/
        delta           (n_sig_edges, n_bootstrap) - bootstrap delta
        delta_mean      (n_sig_edges,) - mean of bootstrap delta
        delta_std       (n_sig_edges,) - std of bootstrap delta
        ci_low          (n_sig_edges,) - 95% CI lower bound
        ci_high         (n_sig_edges,) - 95% CI upper bound
        bias            (n_sig_edges,) - delta_mean - delta_base (bootstrap bias)
    pval/
        bootstrap_pval  (n_sig_edges,) - bootstrap p-value for delta ≠ 0

Design choices
--------------
* Only computes correlations for significant edges, not all edges
* Uses selective indexing to extract relevant pairs from correlation matrix
* Stores both base delta and bootstrap delta for comparison
* Computes bootstrap bias: difference between original delta and mean bootstrap delta
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import rankdata


def rank_zscore_rows(x: np.ndarray, ddof: int = 1) -> np.ndarray:
    """Rank each row, then z-score with the given ddof."""
    r = rankdata(x, axis=1, method="average").astype(np.float32)
    mean = r.mean(axis=1, keepdims=True)
    r_centered = r - mean
    denom = r.std(axis=1, keepdims=True, ddof=ddof)
    denom = np.where(denom < np.finfo(np.float32).eps, np.nan, denom)
    return (r_centered / denom).astype(np.float32)


def flat_idx_to_pair(k: np.ndarray, n_genes: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert flat upper-triangle indices to (gene_i, gene_j) pairs.

    For upper triangle (excluding diagonal), row-major order:
        k = i * (2*N - i - 1) // 2 + (j - i - 1)

    Inverse formula (approximate, then verify):
        i ≈ N - 1 - floor(sqrt(2 * (K - k) - 0.25) - 0.5)
        where K = N*(N-1)//2
    """
    K = n_genes * (n_genes - 1) // 2
    # Use quadratic formula to find i
    # k = i*N - i*(i+1)/2 - 1 + (j - i)
    # Rearranging: i*(2N - i - 1)/2 <= k
    # Approximate: i ≈ (2N - 1 - sqrt((2N-1)² - 8k)) / 2

    k = np.asarray(k)
    # More robust formula
    i = n_genes - 2 - np.floor(np.sqrt(-8 * k + 4 * n_genes * (n_genes - 1) - 7) / 2.0 - 0.5).astype(np.int32)
    j = k + i + 1 - n_genes * (n_genes - 1) // 2 + (n_genes - i) * ((n_genes - i) - 1) // 2
    j = j.astype(np.int32)

    return i, j


def compute_correlations_for_edges(
    expr_subset: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> np.ndarray:
    """
    Compute Spearman correlations for specific edges only.

    Parameters
    ----------
    expr_subset : (n_genes, n_samples_subset)
    edge_i, edge_j : arrays of gene indices for edges

    Returns
    -------
    corr : (n_edges,) correlations for specified edges
    """
    z = rank_zscore_rows(expr_subset)
    n_samples = expr_subset.shape[1]

    # Compute correlations only for specified pairs
    # corr[e] = sum(z[i] * z[j]) / (n-1)
    corr = np.sum(z[edge_i] * z[edge_j], axis=1) / float(n_samples - 1)

    return corr.astype(np.float32)


def bootstrap_significant_edges(
    expr: np.ndarray,
    indices_h5_path: Path,
    base_h5_path: Path,
    out_h5_path: Path,
    gene_index: int = 0,
    edge_selection: str = "sig_edges",
    ci_alpha: float = 0.05,
) -> None:
    """
    Compute bootstrap correlations for significant edges only.

    Parameters
    ----------
    expr : (n_genes, n_samples) expression matrix
    indices_h5_path : path to bootstrap_indices.h5
    base_h5_path : path to base_correlations.h5 from Stage 2a
    out_h5_path : output path
    gene_index : which gene's bootstrap indices to use
    edge_selection : which edges to use: "sig_edges" or "sig_differential"
    ci_alpha : confidence interval alpha (default: 0.05 for 95% CI)
    """
    n_genes = expr.shape[0]

    # Load significant edge indices from Stage 2a
    print(f"Loading significant edges from {base_h5_path}...")
    with h5py.File(base_h5_path, "r") as h5:
        if edge_selection == "sig_differential":
            sig_mask = h5["significant/sig_differential"][:]
        else:  # default: sig_edges
            sig_mask = h5["significant/sig_edges"][:]

        # Get base correlations for significant edges
        sig_indices = np.where(sig_mask)[0].astype(np.int32)
        corr_low_base = h5["low/corr_triu"][:][sig_indices]
        corr_high_base = h5["high/corr_triu"][:][sig_indices]
        delta_base = corr_high_base - corr_low_base

        # Propagate gene names
        if "gene_names" in h5:
            gene_names = [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
        else:
            gene_names = None

    n_sig_edges = len(sig_indices)
    print(f"  Found {n_sig_edges:,} significant edges ({edge_selection})")

    if n_sig_edges == 0:
        print("  No significant edges to bootstrap — writing empty output.")
        out_h5_path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(out_h5_path, "w") as h5:
            meta = h5.create_group("meta")
            meta.attrs["n_sig_edges"] = 0
            meta.attrs["n_bootstrap"] = 0
            meta.attrs["edge_selection_mode"] = edge_selection
            meta.attrs["ci_alpha"] = ci_alpha
            meta.attrs["gene_index_used"] = gene_index
            h5.create_dataset("edges/indices", data=np.array([], dtype=np.int32))
            h5.create_dataset("edges/gene_i", data=np.array([], dtype=np.int32))
            h5.create_dataset("edges/gene_j", data=np.array([], dtype=np.int32))
            h5.create_dataset("base/delta", data=np.array([], dtype=np.float32))
            h5.create_dataset("base/r_low", data=np.array([], dtype=np.float32))
            h5.create_dataset("base/r_high", data=np.array([], dtype=np.float32))
            if gene_names is not None:
                h5.create_dataset("gene_names", data=gene_names, dtype=h5py.string_dtype())
        print(f"  Wrote {out_h5_path}")
        return

    # Convert flat indices to gene pairs
    edge_i, edge_j = flat_idx_to_pair(sig_indices, n_genes)

    # Load bootstrap indices
    print(f"Loading bootstrap indices from {indices_h5_path}...")
    with h5py.File(indices_h5_path, "r") as h5:
        low_boot_indices = h5["indices/low_boot"][gene_index]    # (n_bootstrap, k_resample)
        high_boot_indices = h5["indices/high_boot"][gene_index]  # (n_bootstrap, k_resample)
        n_bootstrap = low_boot_indices.shape[0]

    print(f"  {n_bootstrap} bootstrap replicates")

    # Compute bootstrap correlations for significant edges only
    print(f"Computing bootstrap correlations for {n_sig_edges:,} edges...")
    delta_boot = np.zeros((n_sig_edges, n_bootstrap), dtype=np.float32)

    for b in range(n_bootstrap):
        if (b + 1) % 10 == 0 or b == 0:
            print(f"  Bootstrap {b+1}/{n_bootstrap}...")

        # Subset expression for this bootstrap
        expr_low = np.ascontiguousarray(expr[:, low_boot_indices[b]])
        expr_high = np.ascontiguousarray(expr[:, high_boot_indices[b]])

        # Compute correlations for significant edges only
        corr_low = compute_correlations_for_edges(expr_low, edge_i, edge_j)
        corr_high = compute_correlations_for_edges(expr_high, edge_i, edge_j)

        delta_boot[:, b] = corr_high - corr_low

    # Compute statistics
    print("Computing bootstrap statistics...")
    delta_mean = np.mean(delta_boot, axis=1)
    delta_std = np.std(delta_boot, axis=1, ddof=1)

    # Confidence intervals (percentile method)
    ci_low = np.percentile(delta_boot, 100 * ci_alpha / 2, axis=1).astype(np.float32)
    ci_high = np.percentile(delta_boot, 100 * (1 - ci_alpha / 2), axis=1).astype(np.float32)

    # Bootstrap bias: difference between original delta and mean bootstrap delta
    bias = delta_mean - delta_base

    # Bootstrap p-value for delta ≠ 0
    p_pos = np.mean(delta_boot > 0, axis=1)
    p_neg = np.mean(delta_boot < 0, axis=1)
    pval = 2 * np.minimum(p_pos, p_neg)
    pval = np.maximum(pval, 1.0 / n_bootstrap).astype(np.float32)

    # Summary statistics
    print(f"\n  Bootstrap summary:")
    print(f"    Mean |Δr|: {np.mean(np.abs(delta_mean)):.4f}")
    print(f"    Mean bias: {np.mean(bias):.4f}")
    print(f"    Mean CI width: {np.mean(ci_high - ci_low):.4f}")
    print(f"    Edges with CI excluding 0: {np.sum((ci_low > 0) | (ci_high < 0)):,}")

    # Save results
    print(f"\nSaving to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out_h5_path, "w") as h5:
        # Metadata
        meta = h5.create_group("meta")
        meta.attrs["n_genes"] = n_genes
        meta.attrs["n_sig_edges"] = n_sig_edges
        meta.attrs["n_bootstrap"] = n_bootstrap
        meta.attrs["edge_selection_mode"] = edge_selection
        meta.attrs["ci_alpha"] = ci_alpha
        meta.attrs["gene_index_used"] = gene_index

        # Edge indices
        edges = h5.create_group("edges")
        edges.create_dataset("indices", data=sig_indices, compression="gzip")
        edges.create_dataset("gene_i", data=edge_i, compression="gzip")
        edges.create_dataset("gene_j", data=edge_j, compression="gzip")

        # Base (original) values
        base = h5.create_group("base")
        base.create_dataset("delta", data=delta_base, compression="gzip", compression_opts=4)
        base.create_dataset("r_low", data=corr_low_base, compression="gzip", compression_opts=4)
        base.create_dataset("r_high", data=corr_high_base, compression="gzip", compression_opts=4)

        # Bootstrap results
        boot = h5.create_group("boot")
        boot.create_dataset("delta", data=delta_boot, chunks=(n_sig_edges, 1),
                           compression="gzip", compression_opts=4)
        boot.create_dataset("delta_mean", data=delta_mean.astype(np.float32),
                           compression="gzip", compression_opts=4)
        boot.create_dataset("delta_std", data=delta_std.astype(np.float32),
                           compression="gzip", compression_opts=4)
        boot.create_dataset("ci_low", data=ci_low, compression="gzip", compression_opts=4)
        boot.create_dataset("ci_high", data=ci_high, compression="gzip", compression_opts=4)
        boot.create_dataset("bias", data=bias.astype(np.float32),
                           compression="gzip", compression_opts=4)

        # P-values
        pval_grp = h5.create_group("pval")
        pval_grp.create_dataset("bootstrap_pval", data=pval, compression="gzip", compression_opts=4)

        # Gene names (propagated from base_correlations.h5)
        if gene_names is not None:
            h5.create_dataset("gene_names", data=gene_names, dtype=h5py.string_dtype())

    print("Done!")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 2b: Bootstrap significant edges only."
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
        "--base-h5",
        type=str,
        default=None,
        help="Path to base_correlations.h5 from Stage 2a (single-file mode).",
    )
    parser.add_argument(
        "--base-dir",
        type=str,
        default=None,
        help="Directory of per-gene base_correlations from Stage 2a. Finds file by gene-index prefix.",
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
        "--gene-index",
        type=int,
        default=0,
        help="Which gene's bootstrap indices to use (default: 0).",
    )
    parser.add_argument(
        "--edge-selection",
        type=str,
        choices=["sig_edges", "sig_differential"],
        default="sig_edges",
        help="Which edges to bootstrap: 'sig_edges' (individual & differential) or 'sig_differential' only.",
    )
    parser.add_argument(
        "--ci-alpha",
        type=float,
        default=0.05,
        help="Confidence interval alpha (default: 0.05 for 95%% CI).",
    )
    parser.add_argument(
        "--toy",
        action="store_true",
        help="Use toy data for testing.",
    )
    return parser.parse_args()


def _find_gene_file(directory: Path, gene_index: int) -> Path:
    """Find per-gene h5 file by index prefix (e.g. '0003_*.h5')."""
    prefix = f"{gene_index:04d}_"
    matches = list(directory.glob(f"{prefix}*.h5"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly 1 file matching {prefix}*.h5 in {directory}, found {len(matches)}"
        )
    return matches[0]


def main() -> None:
    args = parse_args()

    if args.toy:
        expr = np.random.default_rng(1).normal(size=(5, 50)).astype(np.float32)
        gene_names = None
    else:
        if args.expr_h5 is None:
            raise SystemExit("ERROR: --expr-h5 is required when not using --toy.")
        print(f"Loading expression from {args.expr_h5}...")
        with h5py.File(args.expr_h5, "r") as f:
            expr = f["expr"][:]
            if "gene_names" in f:
                gene_names = [x.decode() if isinstance(x, bytes) else x for x in f["gene_names"][:]]
            else:
                gene_names = None
        print(f"  Shape: {expr.shape[0]} genes × {expr.shape[1]} samples")

    # Resolve base input path
    if args.base_dir:
        base_h5_path = _find_gene_file(Path(args.base_dir), args.gene_index)
    elif args.base_h5:
        base_h5_path = Path(args.base_h5)
    else:
        raise SystemExit("ERROR: provide --base-h5 or --base-dir.")

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

    bootstrap_significant_edges(
        expr=expr,
        indices_h5_path=Path(args.indices_h5),
        base_h5_path=base_h5_path,
        out_h5_path=out_h5_path,
        gene_index=args.gene_index,
        edge_selection=args.edge_selection,
        ci_alpha=args.ci_alpha,
    )


if __name__ == "__main__":
    main()
