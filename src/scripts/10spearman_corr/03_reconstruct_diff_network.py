#!/usr/bin/env python3
"""
todo lists:
- to clarify that both degree/connectivity for low/high pop and diff net are needed
- use parameter to control when calculate per gene matrix (much less unc=necessary), by default not;
- analyze efficiency of current scripts

Next steps:
- collect information of focus genes across all hdf5 files
- visualization
- calculate network on whole network with edged filtered by fisher-z test and obtain all network metrics
    - combine with rewiring metrics table


Stage 3: Reconstruct Differential Co-expression Network with Topological Analysis.

Combines results from Stage 2a (base correlations) and Stage 2b (bootstrap)
to produce final differential network with:
- Quantitative changes (delta = r_high - r_low)
- Qualitative changes (disappear, new, sign_change)
- Global topological features (connectivity, degree distribution)
- Focus gene neighborhood analysis (1st and 2nd layer partners)
- Rewiring hub metrics for all genes

Pipeline position
-----------------
Stage 1   01subset/01get_extreme_pop_bootstrap.py  →  bootstrap_indices.h5
Stage 2a  02a_calc_base_correlations.py            →  base_correlations.h5
Stage 2b  02b_bootstrap_significant_edges.py       →  bootstrap_significant.h5
Stage 3   THIS SCRIPT                              →  differential_network.h5

Qualitative Changes (with low as reference)
-------------------------------------------
- disappear:   significant in low (|r_low| > threshold) but weak in high
- new:         weak in low but significant in high (|r_high| > threshold)
- sign_change: sign(r_low) != sign(r_high) with both being significant
- strengthen:  same sign, |r_high| > |r_low|
- weaken:      same sign, |r_high| < |r_low|

Output HDF5 layout
------------------
    meta/
        n_genes, n_significant, edge_selection, focus_gene_index
    edges/
        gene_i, gene_j, delta_base, r_low, r_high, ...
        qualitative_score  (n_sig,) int: 0=unchanged, 1=disappear, 2=new, 3=sign_change
        qualitative_label  (n_sig,) str: category name
    topology/
        global/
            n_nodes, n_edges, density, avg_degree, ...
        per_gene/
            degree          (n_genes,) - degree in diff network
            degree_low      (n_genes,) - degree in low network
            degree_high     (n_genes,) - degree in high network
            n_disappear     (n_genes,) - edges that disappeared
            n_new           (n_genes,) - new edges
            n_sign_change   (n_genes,) - sign-changed edges
            rewiring_score  (n_genes,) - total qualitative changes
    focus_gene/
        gene_index
        direct_partners      (n_direct,) - 1st layer neighbors
        indirect_partners    (n_indirect,) - 2nd layer neighbors
        direct_stats         - summary for 1st layer
        two_layer_stats      - summary for 1st + 2nd layers
"""

from __future__ import annotations

import argparse
from pathlib import Path
from collections import defaultdict

import h5py
import numpy as np
from scipy import sparse


# =============================================================================
# Qualitative Change Classification
# =============================================================================

QUAL_UNCHANGED = 0
QUAL_DISAPPEAR = 1
QUAL_NEW = 2
QUAL_SIGN_CHANGE = 3
QUAL_STRENGTHEN = 4
QUAL_WEAKEN = 5

QUAL_LABELS = {
    QUAL_UNCHANGED: "unchanged",
    QUAL_DISAPPEAR: "disappear",
    QUAL_NEW: "new",
    QUAL_SIGN_CHANGE: "sign_change",
    QUAL_STRENGTHEN: "strengthen",
    QUAL_WEAKEN: "weaken",
}


def classify_qualitative_change(
    r_low: np.ndarray,
    r_high: np.ndarray,
    sig_low: np.ndarray,
    sig_high: np.ndarray,
    corr_threshold: float = 0.001,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Classify qualitative changes for each edge (low as reference).

    Parameters
    ----------
    r_low, r_high : correlation values
    sig_low, sig_high : boolean significance masks
    corr_threshold : minimum |r| to be considered "present"

    Returns
    -------
    qual_score : int array with category codes
    qual_label : string array with category names
    """
    n_edges = len(r_low)
    qual_score = np.zeros(n_edges, dtype=np.int8)

    # Determine if edge is "present" (significant or above threshold)
    present_low = sig_low | (np.abs(r_low) >= corr_threshold)
    present_high = sig_high | (np.abs(r_high) >= corr_threshold)

    # Sign of correlations
    sign_low = np.sign(r_low)
    sign_high = np.sign(r_high)

    # Classify each edge
    for i in range(n_edges):
        if present_low[i] and not present_high[i]:
            # Edge disappeared
            qual_score[i] = QUAL_DISAPPEAR
        elif not present_low[i] and present_high[i]:
            # New edge
            qual_score[i] = QUAL_NEW
        elif present_low[i] and present_high[i]:
            if sign_low[i] != sign_high[i] and sign_low[i] != 0 and sign_high[i] != 0:
                # Sign changed
                qual_score[i] = QUAL_SIGN_CHANGE
            elif np.abs(r_high[i]) > np.abs(r_low[i]):
                # Strengthened
                qual_score[i] = QUAL_STRENGTHEN
            elif np.abs(r_high[i]) < np.abs(r_low[i]):
                # Weakened
                qual_score[i] = QUAL_WEAKEN
            else:
                qual_score[i] = QUAL_UNCHANGED
        else:
            # Neither present
            qual_score[i] = QUAL_UNCHANGED

    # Create label array
    qual_label = np.array([QUAL_LABELS[s] for s in qual_score])

    return qual_score, qual_label


# =============================================================================
# Topological Features
# =============================================================================

def compute_global_topology(
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    values: np.ndarray,
    n_genes: int,
) -> dict:
    """
    Compute global topological features of the differential network.
    """
    n_edges = len(gene_i)

    # Build adjacency list
    adj = defaultdict(set)
    for i, j in zip(gene_i, gene_j):
        adj[i].add(j)
        adj[j].add(i)

    # Nodes with at least one edge
    nodes_with_edges = set(gene_i) | set(gene_j)
    n_nodes = len(nodes_with_edges)

    # Degree distribution
    degrees = np.array([len(adj[g]) for g in range(n_genes)])

    # Global metrics
    density = 2 * n_edges / (n_genes * (n_genes - 1)) if n_genes > 1 else 0
    avg_degree = 2 * n_edges / n_genes if n_genes > 0 else 0

    # Degree statistics (only for nodes with edges)
    active_degrees = degrees[degrees > 0]

    return {
        "n_nodes_total": n_genes,
        "n_nodes_active": n_nodes,
        "n_edges": n_edges,
        "density": density,
        "avg_degree": avg_degree,
        "max_degree": int(degrees.max()) if len(degrees) > 0 else 0,
        "median_degree": float(np.median(active_degrees)) if len(active_degrees) > 0 else 0,
        "degrees": degrees,
    }


def compute_per_gene_metrics(
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    qual_score: np.ndarray,
    r_low: np.ndarray,
    r_high: np.ndarray,
    n_genes: int,
) -> dict:
    """
    Compute per-gene topological and rewiring metrics.
    """
    # Initialize arrays
    degree = np.zeros(n_genes, dtype=np.int32)
    degree_low = np.zeros(n_genes, dtype=np.int32)
    degree_high = np.zeros(n_genes, dtype=np.int32)
    n_disappear = np.zeros(n_genes, dtype=np.int32)
    n_new = np.zeros(n_genes, dtype=np.int32)
    n_sign_change = np.zeros(n_genes, dtype=np.int32)
    n_strengthen = np.zeros(n_genes, dtype=np.int32)
    n_weaken = np.zeros(n_genes, dtype=np.int32)
    sum_delta = np.zeros(n_genes, dtype=np.float32)
    sum_abs_delta = np.zeros(n_genes, dtype=np.float32)

    # Threshold for "present"
    corr_threshold = 0.1

    for idx, (i, j) in enumerate(zip(gene_i, gene_j)):
        # Both genes get this edge
        for g in [i, j]:
            degree[g] += 1

            # Count by presence in low/high
            if np.abs(r_low[idx]) >= corr_threshold:
                degree_low[g] += 1
            if np.abs(r_high[idx]) >= corr_threshold:
                degree_high[g] += 1

            # Count by qualitative category
            if qual_score[idx] == QUAL_DISAPPEAR:
                n_disappear[g] += 1
            elif qual_score[idx] == QUAL_NEW:
                n_new[g] += 1
            elif qual_score[idx] == QUAL_SIGN_CHANGE:
                n_sign_change[g] += 1
            elif qual_score[idx] == QUAL_STRENGTHEN:
                n_strengthen[g] += 1
            elif qual_score[idx] == QUAL_WEAKEN:
                n_weaken[g] += 1

            # Delta metrics
            delta = r_high[idx] - r_low[idx]
            sum_delta[g] += delta
            sum_abs_delta[g] += np.abs(delta)

    # Rewiring score: total qualitative changes (excluding strengthen/weaken)
    rewiring_score = n_disappear + n_new + n_sign_change

    return {
        "degree": degree,
        "degree_low": degree_low,
        "degree_high": degree_high,
        "n_disappear": n_disappear,
        "n_new": n_new,
        "n_sign_change": n_sign_change,
        "n_strengthen": n_strengthen,
        "n_weaken": n_weaken,
        "rewiring_score": rewiring_score,
        "sum_delta": sum_delta,
        "sum_abs_delta": sum_abs_delta,
    }


# =============================================================================
# Focus Gene Analysis
# =============================================================================

def analyze_focus_gene(
    focus_gene: int,
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    qual_score: np.ndarray,
    r_low: np.ndarray,
    r_high: np.ndarray,
    delta: np.ndarray,
    n_genes: int,
) -> dict:
    """
    Analyze neighborhood of the focus gene (1st and 2nd layer).

    Returns
    -------
    results : dict with direct partners, indirect partners, and statistics
    """
    # Build adjacency list
    adj = defaultdict(set)
    edge_map = {}  # (i, j) -> edge index
    for idx, (i, j) in enumerate(zip(gene_i, gene_j)):
        adj[i].add(j)
        adj[j].add(i)
        edge_map[(min(i, j), max(i, j))] = idx

    # 1st layer: direct partners
    direct_partners = sorted(adj[focus_gene])

    # 2nd layer: partners of partners (excluding focus gene and direct partners)
    indirect_partners = set()
    for partner in direct_partners:
        for second in adj[partner]:
            if second != focus_gene and second not in direct_partners:
                indirect_partners.add(second)
    indirect_partners = sorted(indirect_partners)

    # Get edges involving focus gene (1st layer edges)
    direct_edges = []
    for partner in direct_partners:
        key = (min(focus_gene, partner), max(focus_gene, partner))
        if key in edge_map:
            direct_edges.append(edge_map[key])

    direct_edges = np.array(direct_edges, dtype=np.int32)

    # Get edges in 2-layer neighborhood (1st + edges between 1st layer)
    two_layer_edges = set(direct_edges)
    for i, p1 in enumerate(direct_partners):
        for p2 in direct_partners[i+1:]:
            key = (min(p1, p2), max(p1, p2))
            if key in edge_map:
                two_layer_edges.add(edge_map[key])

    two_layer_edges = np.array(sorted(two_layer_edges), dtype=np.int32)

    # Statistics for direct edges (1st layer)
    def compute_edge_stats(edge_indices):
        if len(edge_indices) == 0:
            return {
                "n_edges": 0,
                "n_disappear": 0,
                "n_new": 0,
                "n_sign_change": 0,
                "n_strengthen": 0,
                "n_weaken": 0,
                "mean_delta": 0.0,
                "mean_abs_delta": 0.0,
                "mean_r_low": 0.0,
                "mean_r_high": 0.0,
            }

        qs = qual_score[edge_indices]
        return {
            "n_edges": len(edge_indices),
            "n_disappear": int(np.sum(qs == QUAL_DISAPPEAR)),
            "n_new": int(np.sum(qs == QUAL_NEW)),
            "n_sign_change": int(np.sum(qs == QUAL_SIGN_CHANGE)),
            "n_strengthen": int(np.sum(qs == QUAL_STRENGTHEN)),
            "n_weaken": int(np.sum(qs == QUAL_WEAKEN)),
            "mean_delta": float(np.mean(delta[edge_indices])),
            "mean_abs_delta": float(np.mean(np.abs(delta[edge_indices]))),
            "mean_r_low": float(np.mean(r_low[edge_indices])),
            "mean_r_high": float(np.mean(r_high[edge_indices])),
        }

    direct_stats = compute_edge_stats(direct_edges)
    two_layer_stats = compute_edge_stats(two_layer_edges)

    # Per-partner breakdown for direct partners
    partner_details = []
    for partner in direct_partners:
        key = (min(focus_gene, partner), max(focus_gene, partner))
        if key in edge_map:
            idx = edge_map[key]
            partner_details.append({
                "partner": partner,
                "r_low": float(r_low[idx]),
                "r_high": float(r_high[idx]),
                "delta": float(delta[idx]),
                "qual_score": int(qual_score[idx]),
                "qual_label": QUAL_LABELS[qual_score[idx]],
            })

    return {
        "focus_gene": focus_gene,
        "n_direct_partners": len(direct_partners),
        "n_indirect_partners": len(indirect_partners),
        "direct_partners": np.array(direct_partners, dtype=np.int32),
        "indirect_partners": np.array(indirect_partners, dtype=np.int32),
        "direct_stats": direct_stats,
        "two_layer_stats": two_layer_stats,
        "partner_details": partner_details,
    }


# =============================================================================
# Data Loading and Filtering
# =============================================================================

def load_and_process(
    base_h5_path: Path,
    boot_h5_path: Path,
    edge_selection: str = "sig_edges",
    min_effect: float = 0.0,
    require_ci_exclude_zero: bool = True,
    corr_threshold: float = 0.1,
) -> dict:
    """
    Load edges, compute qualitative changes, and apply filters.
    """
    print(f"Loading base correlations from {base_h5_path}...")
    with h5py.File(base_h5_path, "r") as h5:
        n_genes = h5["meta"].attrs["n_genes"]
        n_tests = h5["meta"].attrs["n_tests"]
        fdr_alpha = h5["meta"].attrs["fdr_alpha"]
        k_low = h5["meta"].attrs["k_low"]
        k_high = h5["meta"].attrs["k_high"]

        # Load significance masks
        if edge_selection == "sig_differential":
            base_sig_mask = h5["significant/sig_differential"][:]
        else:
            base_sig_mask = h5["significant/sig_edges"][:]

        # Load full arrays for qualitative analysis
        sig_low_full = h5["significant/sig_low"][:]
        sig_high_full = h5["significant/sig_high"][:]

        # Load correlations for significant edges
        sig_indices = np.where(base_sig_mask)[0]
        corr_low = h5["low/corr_triu"][:][sig_indices]
        corr_high = h5["high/corr_triu"][:][sig_indices]
        pval_diff = h5["diff/pval_triu"][:][sig_indices]
        qval_diff = h5["diff/qval_triu"][:][sig_indices]
        delta_base = h5["diff/delta_triu"][:][sig_indices]

        # Significance for selected edges
        sig_low = sig_low_full[sig_indices]
        sig_high = sig_high_full[sig_indices]

    print(f"  {len(sig_indices):,} edges from base selection ({edge_selection})")

    print(f"Loading bootstrap results from {boot_h5_path}...")
    with h5py.File(boot_h5_path, "r") as h5:
        boot_indices = h5["edges/indices"][:]

        # Verify indices match
        if not np.array_equal(sig_indices, boot_indices):
            print("  Warning: Edge indices mismatch, finding common subset...")
            common_mask = np.isin(sig_indices, boot_indices)
            sig_indices = sig_indices[common_mask]
            corr_low = corr_low[common_mask]
            corr_high = corr_high[common_mask]
            pval_diff = pval_diff[common_mask]
            qval_diff = qval_diff[common_mask]
            delta_base = delta_base[common_mask]
            sig_low = sig_low[common_mask]
            sig_high = sig_high[common_mask]

        gene_i = h5["edges/gene_i"][:]
        gene_j = h5["edges/gene_j"][:]
        delta_boot_mean = h5["boot/delta_mean"][:]
        delta_boot_std = h5["boot/delta_std"][:]
        ci_low = h5["boot/ci_low"][:]
        ci_high = h5["boot/ci_high"][:]
        bias = h5["boot/bias"][:]
        pval_boot = h5["pval/bootstrap_pval"][:]
        ci_alpha = h5["meta"].attrs["ci_alpha"]
        gene_index_used = h5["meta"].attrs.get("gene_index_used", 0)

    print(f"  {len(gene_i):,} edges with bootstrap data")

    # Compute qualitative changes
    print("Computing qualitative changes...")
    qual_score, qual_label = classify_qualitative_change(
        corr_low, corr_high, sig_low, sig_high, corr_threshold
    )

    # Print qualitative summary
    for code, label in QUAL_LABELS.items():
        count = np.sum(qual_score == code)
        if count > 0:
            print(f"    {label}: {count:,} ({100*count/len(qual_score):.1f}%)")

    # Apply filters
    filter_mask = np.ones(len(gene_i), dtype=bool)

    if min_effect > 0:
        effect_mask = np.abs(delta_base) >= min_effect
        filter_mask &= effect_mask
        print(f"  After effect size filter (|Δr| ≥ {min_effect}): {filter_mask.sum():,}")

    if require_ci_exclude_zero:
        ci_mask = (ci_low > 0) | (ci_high < 0)
        filter_mask &= ci_mask
        print(f"  After CI filter (excludes 0): {filter_mask.sum():,}")

    # Apply filter to all arrays
    results = {
        "n_genes": n_genes,
        "n_tests": n_tests,
        "n_significant": int(filter_mask.sum()),
        "edge_selection": edge_selection,
        "min_effect": min_effect,
        "ci_alpha": ci_alpha,
        "fdr_alpha": fdr_alpha,
        "k_low": k_low,
        "k_high": k_high,
        "corr_threshold": corr_threshold,
        "require_ci_exclude_zero": require_ci_exclude_zero,
        "gene_index_used": gene_index_used,
        # Filtered edge data
        "gene_i": gene_i[filter_mask],
        "gene_j": gene_j[filter_mask],
        "delta_base": delta_base[filter_mask],
        "delta_boot_mean": delta_boot_mean[filter_mask],
        "delta_boot_std": delta_boot_std[filter_mask],
        "ci_low": ci_low[filter_mask],
        "ci_high": ci_high[filter_mask],
        "bias": bias[filter_mask],
        "r_low": corr_low[filter_mask],
        "r_high": corr_high[filter_mask],
        "pval_diff": pval_diff[filter_mask],
        "qval_diff": qval_diff[filter_mask],
        "pval_boot": pval_boot[filter_mask],
        "qual_score": qual_score[filter_mask],
        "qual_label": qual_label[filter_mask],
        "sig_low": sig_low[filter_mask],
        "sig_high": sig_high[filter_mask],
    }

    return results


# =============================================================================
# Save Results
# =============================================================================

def save_results(
    results: dict,
    global_topo: dict,
    per_gene: dict,
    focus_analysis: dict,
    out_h5_path: Path,
    out_tsv_path: Path = None,
    gene_ids: list = None,
) -> None:
    """Save all results to HDF5."""
    print(f"\nSaving results to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out_h5_path, "w") as h5:
        # Metadata
        meta = h5.create_group("meta")
        for key in ["n_genes", "n_tests", "n_significant", "edge_selection",
                    "min_effect", "ci_alpha", "fdr_alpha", "k_low", "k_high",
                    "corr_threshold", "require_ci_exclude_zero", "gene_index_used"]:
            meta.attrs[key] = results[key]

        # Edge data
        edges = h5.create_group("edges")
        for key in ["gene_i", "gene_j", "delta_base", "delta_boot_mean", "delta_boot_std",
                    "ci_low", "ci_high", "bias", "r_low", "r_high",
                    "pval_diff", "qval_diff", "pval_boot", "qual_score"]:
            edges.create_dataset(key, data=results[key], compression="gzip", compression_opts=4)

        # Qualitative labels (string)
        dt = h5py.string_dtype(encoding="utf-8")
        edges.create_dataset("qual_label", data=results["qual_label"], dtype=dt)

        # Qualitative summary
        qual_summary = edges.create_group("qual_summary")
        for code, label in QUAL_LABELS.items():
            qual_summary.attrs[f"n_{label}"] = int(np.sum(results["qual_score"] == code))

        # Global topology
        topo = h5.create_group("topology")
        global_grp = topo.create_group("global")
        for key in ["n_nodes_total", "n_nodes_active", "n_edges", "density",
                    "avg_degree", "max_degree", "median_degree"]:
            global_grp.attrs[key] = global_topo[key]

        # Per-gene metrics
        per_gene_grp = topo.create_group("per_gene")
        for key in ["degree", "degree_low", "degree_high", "n_disappear", "n_new",
                    "n_sign_change", "n_strengthen", "n_weaken", "rewiring_score",
                    "sum_delta", "sum_abs_delta"]:
            per_gene_grp.create_dataset(key, data=per_gene[key], compression="gzip")

        # Focus gene analysis
        if focus_analysis is not None:
            focus = h5.create_group("focus_gene")
            focus.attrs["gene_index"] = focus_analysis["focus_gene"]
            focus.attrs["n_direct_partners"] = focus_analysis["n_direct_partners"]
            focus.attrs["n_indirect_partners"] = focus_analysis["n_indirect_partners"]

            focus.create_dataset("direct_partners", data=focus_analysis["direct_partners"])
            focus.create_dataset("indirect_partners", data=focus_analysis["indirect_partners"])

            # Direct stats
            direct_grp = focus.create_group("direct_stats")
            for key, val in focus_analysis["direct_stats"].items():
                direct_grp.attrs[key] = val

            # Two-layer stats
            two_layer_grp = focus.create_group("two_layer_stats")
            for key, val in focus_analysis["two_layer_stats"].items():
                two_layer_grp.attrs[key] = val

        # Sparse matrix for reconstruction
        if results["n_significant"] > 0:
            mat = sparse.csr_matrix(
                (results["delta_base"], (results["gene_i"], results["gene_j"])),
                shape=(results["n_genes"], results["n_genes"])
            )
            # Make symmetric
            mat = mat + mat.T

            matrices = h5.create_group("matrices")
            matrices.create_dataset("delta_data", data=mat.data, compression="gzip")
            matrices.create_dataset("delta_indices", data=mat.indices, compression="gzip")
            matrices.create_dataset("delta_indptr", data=mat.indptr, compression="gzip")
            matrices.attrs["shape"] = mat.shape
            matrices.attrs["format"] = "csr"

    print(f"  Saved {results['n_significant']:,} edges")

    # Save rewiring hub table (TSV)
    if out_tsv_path is not None:
        save_rewiring_hub_table(per_gene, out_tsv_path, gene_ids)


def save_rewiring_hub_table(per_gene: dict, out_path: Path, gene_ids: list = None) -> None:
    """Save per-gene rewiring metrics as TSV."""
    print(f"Saving rewiring hub table to {out_path}...")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_genes = len(per_gene["degree"])

    with open(out_path, "w") as f:
        # Header
        f.write("gene_idx\tgene_id\tdegree\tdegree_low\tdegree_high\t"
                "n_disappear\tn_new\tn_sign_change\tn_strengthen\tn_weaken\t"
                "rewiring_score\tsum_delta\tsum_abs_delta\n")

        # Sort by rewiring score (descending)
        sort_idx = np.argsort(per_gene["rewiring_score"])[::-1]

        for i in sort_idx:
            # Skip genes with no edges
            if per_gene["degree"][i] == 0:
                continue

            gene_id = gene_ids[i] if gene_ids else f"gene_{i}"
            f.write(f"{i}\t{gene_id}\t"
                   f"{per_gene['degree'][i]}\t"
                   f"{per_gene['degree_low'][i]}\t"
                   f"{per_gene['degree_high'][i]}\t"
                   f"{per_gene['n_disappear'][i]}\t"
                   f"{per_gene['n_new'][i]}\t"
                   f"{per_gene['n_sign_change'][i]}\t"
                   f"{per_gene['n_strengthen'][i]}\t"
                   f"{per_gene['n_weaken'][i]}\t"
                   f"{per_gene['rewiring_score'][i]}\t"
                   f"{per_gene['sum_delta'][i]:.4f}\t"
                   f"{per_gene['sum_abs_delta'][i]:.4f}\n")

    n_with_edges = np.sum(per_gene["degree"] > 0)
    print(f"  Saved {n_with_edges:,} genes with edges")


# =============================================================================
# Summary Printing
# =============================================================================

def print_summary(results: dict, global_topo: dict, per_gene: dict, focus_analysis: dict) -> None:
    """Print comprehensive summary."""
    print("\n" + "=" * 70)
    print("DIFFERENTIAL CO-EXPRESSION NETWORK ANALYSIS SUMMARY")
    print("=" * 70)

    print(f"\n[NETWORK OVERVIEW]")
    print(f"  Total genes: {results['n_genes']:,}")
    print(f"  Total possible edges: {results['n_tests']:,}")
    print(f"  Significant differential edges: {results['n_significant']:,} "
          f"({100 * results['n_significant'] / results['n_tests']:.4f}%)")

    print(f"\n[GLOBAL TOPOLOGY]")
    print(f"  Active nodes (with edges): {global_topo['n_nodes_active']:,}")
    print(f"  Network density: {global_topo['density']:.6f}")
    print(f"  Average degree: {global_topo['avg_degree']:.2f}")
    print(f"  Max degree: {global_topo['max_degree']}")
    print(f"  Median degree: {global_topo['median_degree']:.1f}")

    print(f"\n[QUALITATIVE CHANGES] (low → high)")
    for code, label in QUAL_LABELS.items():
        count = np.sum(results["qual_score"] == code)
        if count > 0:
            pct = 100 * count / results["n_significant"]
            print(f"  {label:12s}: {count:>8,} ({pct:5.1f}%)")

    print(f"\n[QUANTITATIVE CHANGES]")
    print(f"  Mean Δr: {np.mean(results['delta_base']):+.4f}")
    print(f"  Mean |Δr|: {np.mean(np.abs(results['delta_base'])):.4f}")
    print(f"  Max |Δr|: {np.max(np.abs(results['delta_base'])):.4f}")
    print(f"  Mean bootstrap bias: {np.mean(results['bias']):.4f}")

    print(f"\n[TOP REWIRING HUBS]")
    top_idx = np.argsort(per_gene["rewiring_score"])[::-1][:10]
    print(f"  {'Gene':>8} {'Degree':>8} {'Disappear':>10} {'New':>6} {'SignChg':>8} {'RewireScore':>12}")
    for i in top_idx:
        if per_gene["degree"][i] > 0:
            print(f"  {i:>8} {per_gene['degree'][i]:>8} "
                  f"{per_gene['n_disappear'][i]:>10} "
                  f"{per_gene['n_new'][i]:>6} "
                  f"{per_gene['n_sign_change'][i]:>8} "
                  f"{per_gene['rewiring_score'][i]:>12}")

    if focus_analysis is not None:
        fg = focus_analysis["focus_gene"]
        print(f"\n[FOCUS GENE ANALYSIS] (gene {fg})")
        print(f"  Direct partners (1st layer): {focus_analysis['n_direct_partners']}")
        print(f"  Indirect partners (2nd layer): {focus_analysis['n_indirect_partners']}")

        ds = focus_analysis["direct_stats"]
        print(f"\n  1st Layer Stats:")
        print(f"    Edges: {ds['n_edges']}")
        print(f"    Disappear: {ds['n_disappear']}, New: {ds['n_new']}, Sign change: {ds['n_sign_change']}")
        print(f"    Mean Δr: {ds['mean_delta']:+.4f}, Mean |Δr|: {ds['mean_abs_delta']:.4f}")
        print(f"    Mean r_low: {ds['mean_r_low']:.4f}, Mean r_high: {ds['mean_r_high']:.4f}")

        ts = focus_analysis["two_layer_stats"]
        print(f"\n  1st + 2nd Layer Stats:")
        print(f"    Edges: {ts['n_edges']}")
        print(f"    Disappear: {ts['n_disappear']}, New: {ts['n_new']}, Sign change: {ts['n_sign_change']}")

        if focus_analysis["partner_details"]:
            print(f"\n  Direct Partner Details:")
            print(f"    {'Partner':>8} {'r_low':>8} {'r_high':>8} {'Δr':>8} {'Category':>12}")
            for pd in focus_analysis["partner_details"][:10]:  # Show top 10
                print(f"    {pd['partner']:>8} {pd['r_low']:>8.3f} {pd['r_high']:>8.3f} "
                      f"{pd['delta']:>+8.3f} {pd['qual_label']:>12}")


# =============================================================================
# Main
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 3: Reconstruct differential network with topology analysis."
    )
    parser.add_argument(
        "--base-h5", type=str, required=True,
        help="Path to base_correlations.h5 from Stage 2a.",
    )
    parser.add_argument(
        "--boot-h5", type=str, required=True,
        help="Path to bootstrap_significant.h5 from Stage 2b.",
    )
    parser.add_argument(
        "--out-h5", type=str, default="results/differential_network.h5",
        help="Output HDF5 path.",
    )
    parser.add_argument(
        "--out-tsv", type=str, default=None,
        help="Output TSV path for rewiring hub table.",
    )
    parser.add_argument(
        "--edge-selection", type=str, choices=["sig_edges", "sig_differential"],
        default="sig_edges",
        help="Edge selection mode.",
    )
    parser.add_argument(
        "--min-effect", type=float, default=0.0,
        help="Minimum |Δr| threshold.",
    )
    parser.add_argument(
        "--corr-threshold", type=float, default=0.1,
        help="Correlation threshold for qualitative 'present' (default: 0.1).",
    )
    parser.add_argument(
        "--no-ci-filter", action="store_true",
        help="Don't require bootstrap CI to exclude 0.",
    )
    parser.add_argument(
        "--focus-gene", type=int, default=None,
        help="Focus gene index for neighborhood analysis (default: gene used for subsetting).",
    )
    parser.add_argument(
        "--gene-ids", type=str, default=None,
        help="Path to file with gene IDs (one per line) for labeling.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Load and process data
    results = load_and_process(
        base_h5_path=Path(args.base_h5),
        boot_h5_path=Path(args.boot_h5),
        edge_selection=args.edge_selection,
        min_effect=args.min_effect,
        require_ci_exclude_zero=not args.no_ci_filter,
        corr_threshold=args.corr_threshold,
    )

    if results["n_significant"] == 0:
        print("\nNo significant edges after filtering. Exiting.")
        return

    # Compute global topology
    print("\nComputing global topology...")
    global_topo = compute_global_topology(
        results["gene_i"], results["gene_j"],
        results["delta_base"], results["n_genes"]
    )

    # Compute per-gene metrics
    print("Computing per-gene metrics...")
    per_gene = compute_per_gene_metrics(
        results["gene_i"], results["gene_j"],
        results["qual_score"], results["r_low"], results["r_high"],
        results["n_genes"]
    )

    # Focus gene analysis
    focus_gene = args.focus_gene if args.focus_gene is not None else results["gene_index_used"]
    print(f"Analyzing focus gene {focus_gene}...")
    focus_analysis = analyze_focus_gene(
        focus_gene,
        results["gene_i"], results["gene_j"],
        results["qual_score"], results["r_low"], results["r_high"],
        results["delta_base"], results["n_genes"]
    )

    # Load gene IDs if provided
    gene_ids = None
    if args.gene_ids:
        with open(args.gene_ids) as f:
            gene_ids = [line.strip() for line in f]

    # Save results
    save_results(
        results, global_topo, per_gene, focus_analysis,
        out_h5_path=Path(args.out_h5),
        out_tsv_path=Path(args.out_tsv) if args.out_tsv else None,
        gene_ids=gene_ids,
    )

    # Print summary
    print_summary(results, global_topo, per_gene, focus_analysis)


if __name__ == "__main__":
    main()
