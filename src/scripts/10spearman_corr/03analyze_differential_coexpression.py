#!/usr/bin/env python3
"""
Stage 3: Differential Co-expression Analysis

Detect edges with significantly different correlation between low/high subpopulations
using a two-stage bootstrap-based approach.

Input (from Stage 2):
    corr/gene_XXXXX.h5
        low/base/corr_triu    (n_tests,)
        low/boot/corr_triu    (n_bootstrap, n_tests)
        high/base/corr_triu   (n_tests,)
        high/boot/corr_triu   (n_bootstrap, n_tests)

Output:
    differential_edges.h5
        edges/idx              (n_significant, 2) - gene pair indices
        edges/delta_mean       (n_significant,) - mean difference
        edges/ci_low           (n_significant,) - 95% CI lower
        edges/ci_high          (n_significant,) - 95% CI upper  
        edges/pval             (n_significant,) - bootstrap p-value
        edges/qval             (n_significant,) - FDR-corrected q-value
        edges/r_low_base       (n_significant,) - base corr in low
        edges/r_high_base      (n_significant,) - base corr in high

Strategy:
    1. Load all per-gene HDF5 files (parallelizable)
    2. Screen with Fisher's Z test (fast, all edges)
    3. Test candidates with bootstrap difference (robust)
    4. FDR correction (Benjamini-Hochberg)
    5. Filter by effect size and significance
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import h5py
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


def load_correlation_data(corr_dir: Path, n_genes: int) -> Tuple[np.ndarray, ...]:
    """
    Load correlation data from all per-gene HDF5 files.
    
    Returns
    -------
    low_base, high_base : (n_tests,) float32
        Base correlations for each edge
    low_boot, high_boot : (n_tests, n_bootstrap) float32
        Bootstrap correlations for each edge
    """
    # Determine n_tests from first file
    first_file = corr_dir / "gene_00000.h5"
    with h5py.File(first_file, "r") as h5:
        n_tests = h5["low/base/corr_triu"].shape[0]
        n_bootstrap = h5["low/boot/corr_triu"].shape[0]
    
    print(f"Loading {n_genes} genes, {n_tests:,} edges, {n_bootstrap} bootstraps...")
    
    # Pre-allocate arrays
    low_base = np.zeros(n_tests, dtype=np.float32)
    high_base = np.zeros(n_tests, dtype=np.float32)
    low_boot = np.zeros((n_tests, n_bootstrap), dtype=np.float32)
    high_boot = np.zeros((n_tests, n_bootstrap), dtype=np.float32)
    
    # Load from first gene file (all genes should have same edges)
    # In your current setup, each gene file has correlations for ALL gene pairs
    # So we just need to load from one representative file
    with h5py.File(first_file, "r") as h5:
        low_base[:] = h5["low/base/corr_triu"][:]
        high_base[:] = h5["high/base/corr_triu"][:]
        low_boot[:, :] = h5["low/boot/corr_triu"][:].T  # Transpose to (n_tests, n_bootstrap)
        high_boot[:, :] = h5["high/boot/corr_triu"][:].T
    
    print(f"  Loaded {n_tests:,} edges")
    return low_base, high_base, low_boot, high_boot


def fisher_z_screening(
    low_base: np.ndarray,
    high_base: np.ndarray,
    n_low: int,
    n_high: int,
    p_threshold: float = 0.001
) -> np.ndarray:
    """
    Stage 1: Fast screening with Fisher's Z-transform test.
    
    Returns
    -------
    candidate_idx : array of int
        Indices of edges passing screening threshold
    """
    print(f"\n[Stage 1: Fisher's Z Screening]")
    
    # Fisher's Z transform (handle edge cases near ±1)
    r_low_clip = np.clip(low_base, -0.9999, 0.9999)
    r_high_clip = np.clip(high_base, -0.9999, 0.9999)
    
    z_low = 0.5 * np.log((1 + r_low_clip) / (1 - r_low_clip))
    z_high = 0.5 * np.log((1 + r_high_clip) / (1 - r_high_clip))
    
    # Standard error of difference
    se_diff = np.sqrt(1.0 / (n_low - 3) + 1.0 / (n_high - 3))
    
    # Z-score and p-value
    z_score = (z_high - z_low) / se_diff
    pvals = 2 * (1 - stats.norm.cdf(np.abs(z_score)))
    
    # Select candidates
    candidate_idx = np.where(pvals < p_threshold)[0]
    
    print(f"  Tested {len(pvals):,} edges")
    print(f"  Candidates (p < {p_threshold}): {len(candidate_idx):,} ({100*len(candidate_idx)/len(pvals):.2f}%)")
    
    return candidate_idx


def bootstrap_difference_test(
    low_boot: np.ndarray,
    high_boot: np.ndarray,
    alpha: float = 0.05
) -> Tuple[np.ndarray, ...]:
    """
    Stage 2: Bootstrap-based difference test for candidate edges.
    
    Parameters
    ----------
    low_boot, high_boot : (n_candidates, n_bootstrap) array
        Bootstrap correlations
    
    Returns
    -------
    delta_mean : (n_candidates,) - mean difference
    ci_low, ci_high : (n_candidates,) - confidence interval
    pvals : (n_candidates,) - bootstrap p-values
    """
    print(f"\n[Stage 2: Bootstrap Difference Test]")
    
    # Compute difference for each bootstrap replicate
    delta_boot = high_boot - low_boot  # (n_candidates, n_bootstrap)
    
    # Mean difference
    delta_mean = np.mean(delta_boot, axis=1)
    
    # 95% confidence intervals (percentile method)
    ci_low = np.percentile(delta_boot, 100 * alpha / 2, axis=1)
    ci_high = np.percentile(delta_boot, 100 * (1 - alpha / 2), axis=1)
    
    # Two-tailed p-value: proportion of bootstrap samples crossing 0
    p_pos = np.mean(delta_boot > 0, axis=1)
    p_neg = np.mean(delta_boot < 0, axis=1)
    pvals = 2 * np.minimum(p_pos, p_neg)
    
    # Clip p-values to minimum resolution (1 / n_bootstrap)
    n_bootstrap = delta_boot.shape[1]
    pvals = np.maximum(pvals, 1.0 / n_bootstrap)
    
    print(f"  Tested {len(delta_mean):,} candidates")
    print(f"  Mean |Δr|: {np.mean(np.abs(delta_mean)):.3f}")
    print(f"  Median p-value: {np.median(pvals):.4f}")
    
    return delta_mean, ci_low, ci_high, pvals


def analyze_differential_coexpression(
    corr_dir: Path,
    n_genes: int,
    n_low: int,
    n_high: int,
    fdr_threshold: float = 0.05,
    effect_size_threshold: float = 0.1,
    screening_p: float = 0.001,
) -> dict:
    """
    Two-stage differential co-expression analysis.
    
    Parameters
    ----------
    corr_dir : Path
        Directory with per-gene correlation HDF5 files
    n_genes : int
        Total number of genes
    n_low, n_high : int
        Sample sizes in low/high subpopulations
    fdr_threshold : float
        FDR threshold for significance (default: 0.05)
    effect_size_threshold : float
        Minimum |Δr| for biological relevance (default: 0.1)
    screening_p : float
        P-value threshold for Stage 1 screening (default: 0.001)
    
    Returns
    -------
    results : dict
        Significant differential edges with statistics
    """
    # Load data
    low_base, high_base, low_boot, high_boot = load_correlation_data(corr_dir, n_genes)
    n_edges = len(low_base)
    
    # Stage 1: Fisher's Z screening
    candidate_idx = fisher_z_screening(
        low_base, high_base, n_low, n_high, p_threshold=screening_p
    )
    
    if len(candidate_idx) == 0:
        print("\nNo candidates passed screening. Try increasing screening_p.")
        return None
    
    # Stage 2: Bootstrap testing on candidates
    delta_mean, ci_low, ci_high, pvals = bootstrap_difference_test(
        low_boot[candidate_idx, :],
        high_boot[candidate_idx, :]
    )
    
    # FDR correction
    print(f"\n[FDR Correction]")
    reject, qvals, _, _ = multipletests(
        pvals, alpha=fdr_threshold, method='fdr_bh'
    )
    
    n_sig_fdr = np.sum(reject)
    print(f"  Significant at FDR < {fdr_threshold}: {n_sig_fdr:,} / {len(pvals):,}")
    
    # Filter by effect size
    effect_filter = np.abs(delta_mean) >= effect_size_threshold
    final_filter = reject & effect_filter
    
    n_final = np.sum(final_filter)
    print(f"  After effect size filter (|Δr| ≥ {effect_size_threshold}): {n_final:,}")
    
    if n_final == 0:
        print("\nNo edges passed both FDR and effect size thresholds.")
        return None
    
    # Extract significant edges
    sig_candidate_idx = candidate_idx[final_filter]
    
    # Convert flat indices to (gene_i, gene_j) pairs
    # Assuming edges are stored in row-major upper-triangle order
    edge_pairs = np.zeros((n_final, 2), dtype=np.int32)
    for i, edge_idx in enumerate(sig_candidate_idx):
        # Inverse of np.triu_indices: convert flat index to (row, col)
        # This is approximate; you may need to adjust based on your indexing
        gene_i = int((-1 + np.sqrt(1 + 8 * edge_idx)) / 2)
        gene_j = edge_idx - gene_i * (2 * n_genes - gene_i - 1) // 2
        edge_pairs[i] = [gene_i, gene_j]
    
    results = {
        'edge_pairs': edge_pairs,
        'delta_mean': delta_mean[final_filter],
        'ci_low': ci_low[final_filter],
        'ci_high': ci_high[final_filter],
        'pval': pvals[final_filter],
        'qval': qvals[final_filter],
        'r_low_base': low_base[sig_candidate_idx],
        'r_high_base': high_base[sig_candidate_idx],
        'n_edges_total': n_edges,
        'n_candidates': len(candidate_idx),
        'n_significant': n_final,
        'fdr_threshold': fdr_threshold,
        'effect_size_threshold': effect_size_threshold
    }
    
    # Sort by effect size (absolute difference)
    sort_idx = np.argsort(np.abs(results['delta_mean']))[::-1]
    for key in ['edge_pairs', 'delta_mean', 'ci_low', 'ci_high', 'pval', 'qval', 'r_low_base', 'r_high_base']:
        if isinstance(results[key], np.ndarray):
            results[key] = results[key][sort_idx]
    
    return results


def save_results(results: dict, out_path: Path) -> None:
    """Save differential co-expression results to HDF5."""
    print(f"\n[Saving Results]")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    with h5py.File(out_path, "w") as h5:
        # Metadata
        meta = h5.create_group("meta")
        meta.attrs['n_edges_total'] = results['n_edges_total']
        meta.attrs['n_candidates'] = results['n_candidates']
        meta.attrs['n_significant'] = results['n_significant']
        meta.attrs['fdr_threshold'] = results['fdr_threshold']
        meta.attrs['effect_size_threshold'] = results['effect_size_threshold']
        
        # Edge data
        edges = h5.create_group("edges")
        for key in ['edge_pairs', 'delta_mean', 'ci_low', 'ci_high', 'pval', 'qval', 'r_low_base', 'r_high_base']:
            edges.create_dataset(key, data=results[key], compression='gzip', compression_opts=4)
    
    print(f"  Saved {results['n_significant']:,} edges to {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Differential co-expression analysis using bootstrap."
    )
    parser.add_argument(
        "--corr-dir",
        type=str,
        required=True,
        help="Directory with per-gene correlation HDF5 files",
    )
    parser.add_argument(
        "--n-genes",
        type=int,
        required=True,
        help="Total number of genes",
    )
    parser.add_argument(
        "--n-low",
        type=int,
        required=True,
        help="Number of samples in low expression subpopulation",
    )
    parser.add_argument(
        "--n-high",
        type=int,
        required=True,
        help="Number of samples in high expression subpopulation",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="results/differential_edges.h5",
        help="Output HDF5 file path",
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=0.05,
        help="FDR threshold (default: 0.05)",
    )
    parser.add_argument(
        "--min-effect",
        type=float,
        default=0.1,
        help="Minimum effect size |Δr| (default: 0.1)",
    )
    parser.add_argument(
        "--screening-p",
        type=float,
        default=0.001,
        help="P-value threshold for Fisher's Z screening (default: 0.001)",
    )
    
    args = parser.parse_args()
    
    results = analyze_differential_coexpression(
        corr_dir=Path(args.corr_dir),
        n_genes=args.n_genes,
        n_low=args.n_low,
        n_high=args.n_high,
        fdr_threshold=args.fdr,
        effect_size_threshold=args.min_effect,
        screening_p=args.screening_p,
    )
    
    if results is not None:
        save_results(results, Path(args.out))
        
        # Print summary
        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)
        print(f"Total edges tested: {results['n_edges_total']:,}")
        print(f"Candidates after screening: {results['n_candidates']:,}")
        print(f"Significant (FDR < {results['fdr_threshold']}, |Δr| ≥ {results['effect_size_threshold']}): {results['n_significant']:,}")
        print(f"\nTop 5 edges by effect size:")
        for i in range(min(5, results['n_significant'])):
            print(f"  {i+1}. Genes ({results['edge_pairs'][i, 0]}, {results['edge_pairs'][i, 1]}): "
                  f"Δr = {results['delta_mean'][i]:.3f}, "
                  f"q = {results['qval'][i]:.4f}, "
                  f"[{results['ci_low'][i]:.3f}, {results['ci_high'][i]:.3f}]")


if __name__ == "__main__":
    main()
