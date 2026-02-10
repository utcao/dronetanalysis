#!/usr/bin/env python3
"""
Stage 3: Reconstruct Differential Co-expression Network with Topological Analysis.

Combines results from Stage 2a (base correlations) and Stage 2b (bootstrap)
to produce final differential network with:
- Quantitative changes (delta = r_high - r_low)
- Qualitative changes (disappear, new, sign_change)
- Global topological features for LOW, HIGH, and DIFF networks
- Focus gene neighborhood analysis (1st and 2nd layer partners)
- Optional: per-gene metrics (controlled by --calc-per-gene-metrics)

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
        global_low/    - Low network: n_nodes, n_edges, density, avg_degree, ...
        global_high/   - High network: n_nodes, n_edges, density, avg_degree, ...
        global_diff/   - Differential network: n_nodes, n_edges, density, avg_degree, ...
        per_gene/      - (Optional, --calc-per-gene-metrics)
            degree, degree_low, degree_high, n_disappear, n_new, n_sign_change, rewiring_score
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
    Vectorized implementation for efficiency.

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

    # Determine if edge is "present" (significant and above threshold)
    present_low = sig_low & (np.abs(r_low) >= corr_threshold)
    present_high = sig_high & (np.abs(r_high) >= corr_threshold)

    # Sign of correlations
    sign_low = np.sign(r_low)
    sign_high = np.sign(r_high)

    # Vectorized classification
    # Case 1: Disappear (present in low, absent in high)
    mask_disappear = present_low & ~present_high
    qual_score[mask_disappear] = QUAL_DISAPPEAR

    # Case 2: New (absent in low, present in high)
    mask_new = ~present_low & present_high
    qual_score[mask_new] = QUAL_NEW

    # Case 3-5: Both present
    mask_both = present_low & present_high

    # Sign change (both present, signs differ, neither zero)
    mask_sign_change = mask_both & (sign_low != sign_high) & (sign_low != 0) & (sign_high != 0)
    qual_score[mask_sign_change] = QUAL_SIGN_CHANGE

    # Strengthen (same sign, |r_high| > |r_low|)
    abs_low, abs_high = np.abs(r_low), np.abs(r_high)
    mask_strengthen = mask_both & ~mask_sign_change & (abs_high > abs_low)
    qual_score[mask_strengthen] = QUAL_STRENGTHEN

    # Weaken (same sign, |r_high| < |r_low|)
    mask_weaken = mask_both & ~mask_sign_change & (abs_high < abs_low)
    qual_score[mask_weaken] = QUAL_WEAKEN

    # Remaining are QUAL_UNCHANGED (already 0)

    # Create label array (vectorized lookup)
    label_arr = np.array([QUAL_LABELS[i] for i in range(6)])
    qual_label = label_arr[qual_score]

    return qual_score, qual_label


# =============================================================================
# Topological Features
# =============================================================================

def compute_global_topology(
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    values: np.ndarray,
    n_genes: int,
    corr_threshold: float = 0.0001,
) -> dict:
    """
    Compute global topological features of a network.

    Parameters
    ----------
    gene_i, gene_j : edge indices
    values : correlation values for each edge
    n_genes : total number of genes
    corr_threshold : minimum |r| for edge to be considered "present"

    Returns
    -------
    dict with n_nodes_active, n_edges, density, avg_degree, max_degree, median_degree, degrees
    """
    # Filter edges by threshold
    mask = np.abs(values) >= corr_threshold
    gi, gj = gene_i[mask], gene_j[mask]
    n_edges = len(gi)

    if n_edges == 0:
        return {
            "n_nodes_total": n_genes,
            "n_nodes_active": 0,
            "n_edges": 0,
            "density": 0.0,
            "avg_degree": 0.0,
            "max_degree": 0,
            "median_degree": 0.0,
            "degrees": np.zeros(n_genes, dtype=np.int32),
        }

    # Build adjacency list
    adj = defaultdict(set)
    for i, j in zip(gi, gj):
        adj[i].add(j)
        adj[j].add(i)

    # Nodes with at least one edge
    nodes_with_edges = set(gi) | set(gj)
    n_nodes = len(nodes_with_edges)

    # Degree distribution
    degrees = np.array([len(adj[g]) for g in range(n_genes)], dtype=np.int32)

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


def compute_all_global_topologies(
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    r_low: np.ndarray,
    r_high: np.ndarray,
    delta: np.ndarray,
    n_genes: int,
    corr_threshold: float = 0.0001,
) -> dict:
    """
    Compute global topology for low, high, and differential networks.

    Returns
    -------
    dict with 'low', 'high', 'diff' keys, each containing topology stats
    """
    return {
        "low": compute_global_topology(gene_i, gene_j, r_low, n_genes, corr_threshold),
        "high": compute_global_topology(gene_i, gene_j, r_high, n_genes, corr_threshold),
        "diff": compute_global_topology(gene_i, gene_j, delta, n_genes, corr_threshold=0.0),  # all sig edges for diff
    }


def compute_per_gene_metrics(
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    qual_score: np.ndarray,
    r_low: np.ndarray,
    r_high: np.ndarray,
    n_genes: int,
    corr_threshold: float = 0.0001, # Threshold for "present"
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
                "mean_delta": 0.0,
                "mean_abs_delta": 0.0,
                "mean_r_low": 0.0,
                "mean_r_high": 0.0,
                "n_strengthen": 0,
                "n_weaken": 0,
            }

        qs = qual_score[edge_indices]
        return {
            "n_edges": len(edge_indices),
            "n_disappear": int(np.sum(qs == QUAL_DISAPPEAR)),
            "n_new": int(np.sum(qs == QUAL_NEW)),
            "n_sign_change": int(np.sum(qs == QUAL_SIGN_CHANGE)),
            "mean_delta": float(np.mean(delta[edge_indices])),
            "mean_abs_delta": float(np.mean(np.abs(delta[edge_indices]))),
            "mean_r_low": float(np.mean(r_low[edge_indices])),
            "mean_r_high": float(np.mean(r_high[edge_indices])),
            "n_strengthen": int(np.sum(qs == QUAL_STRENGTHEN)),
            "n_weaken": int(np.sum(qs == QUAL_WEAKEN)),
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
    corr_threshold: float = 0.0001,
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

        # Propagate gene names
        if "gene_names" in h5:
            gene_names = [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
        else:
            gene_names = None

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

    n_sig = len(sig_indices)
    print(f"  {n_sig:,} edges from base selection ({edge_selection}): \
        obtain the indecies of edges in 1D array")

    print(f"Loading bootstrap results from {boot_h5_path}...")
    with h5py.File(boot_h5_path, "r") as h5:
        boot_n_sig = int(h5["meta"].attrs["n_sig_edges"])
        boot_indices = h5["edges/indices"][:]

        if boot_n_sig == 0 or n_sig == 0:
            # No significant edges — return empty results
            print("  No significant edges found. Writing empty network.")
            gene_i = np.array([], dtype=np.int32)
            gene_j = np.array([], dtype=np.int32)
            corr_low = np.array([], dtype=np.float32)
            corr_high = np.array([], dtype=np.float32)
            delta_base = np.array([], dtype=np.float32)
            pval_diff = np.array([], dtype=np.float32)
            qval_diff = np.array([], dtype=np.float32)
            sig_low = np.array([], dtype=bool)
            sig_high = np.array([], dtype=bool)
            delta_boot_mean = np.array([], dtype=np.float32)
            delta_boot_std = np.array([], dtype=np.float32)
            ci_low = np.array([], dtype=np.float32)
            ci_high = np.array([], dtype=np.float32)
            bias = np.array([], dtype=np.float32)
            pval_boot = np.array([], dtype=np.float32)
            ci_alpha = h5["meta"].attrs["ci_alpha"]
            # Error bug: no gene_index_used
            gene_index_used = h5["meta"].attrs.get("gene_index_used", 0)
            sig_indices = np.array([], dtype=np.int32)
        else:
            # Verify indices match
            if not np.array_equal(sig_indices, boot_indices):
                print("  Warning: Edge indices mismatch, finding common subset...")
                raise ValueError(
                    f"Critical Error: Edge indices mismatch! "
                    f"Significant indices ({len(sig_indices)}) do not match "
                    f"Bootstrap indices ({len(boot_indices)}). Pipeline stopped."
                )
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
        "gene_names": gene_names,
    }

    return results


# =============================================================================
# Save Results
# =============================================================================

def save_results(
    results: dict,
    all_topo: dict,
    per_gene: dict | None,
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

        # Qualitative labels (string) - convert numpy U-strings to bytes for h5py
        dt = h5py.string_dtype(encoding="utf-8")
        qual_labels = np.array(results["qual_label"], dtype="S")
        edges.create_dataset("qual_label", data=qual_labels, dtype=dt)

        # Gene names (propagated from base_correlations.h5)
        if gene_ids is not None:
            h5.create_dataset("gene_names", data=gene_ids, dtype=h5py.string_dtype())

        # Qualitative summary
        qual_summary = edges.create_group("qual_summary")
        for code, label in QUAL_LABELS.items():
            qual_summary.attrs[f"n_{label}"] = int(np.sum(results["qual_score"] == code))

        # Global topology for low, high, diff networks
        topo = h5.create_group("topology")
        for net_name in ["low", "high", "diff"]:
            net_topo = all_topo[net_name]
            grp = topo.create_group(f"global_{net_name}")
            for key in ["n_nodes_total", "n_nodes_active", "n_edges", "density",
                        "avg_degree", "max_degree", "median_degree"]:
                grp.attrs[key] = net_topo[key]
            # Also save degree array for each network
            grp.create_dataset("degrees", data=net_topo["degrees"], compression="gzip")

        # Per-gene metrics (optional)
        if per_gene is not None:
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

    # Save rewiring hub table (TSV) - only if per-gene metrics computed
    if out_tsv_path is not None and per_gene is not None:
        save_rewiring_hub_table(per_gene, out_tsv_path, gene_ids)
    elif out_tsv_path is not None:
        print(f"  Skipping TSV output (per-gene metrics not computed)")


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
# Per-Gene Collection (directory mode)
# =============================================================================

def _find_gene_file(directory: Path, gene_index: int) -> Path | None:
    """Find per-gene h5 file by index prefix. Returns None if not found."""
    prefix = f"{gene_index:04d}_"
    matches = list(directory.glob(f"{prefix}*.h5"))
    return matches[0] if len(matches) == 1 else None


def collect_per_gene_networks(
    base_dir: Path,
    boot_dir: Path,
    out_h5_path: Path,
    out_tsv_path: Path | None = None,
    edge_selection: str = "sig_edges",
    min_effect: float = 0.0,
    require_ci_exclude_zero: bool = True,
    corr_threshold: float = 0.0001,
) -> None:
    """
    Collect per-gene network results from directories of per-gene h5 files.

    For each reference gene g, reads base_correlations and bootstrap_significant
    files, applies filters, computes qualitative classification and topology,
    and writes a summary.
    """
    # Discover per-gene files
    base_files = sorted(base_dir.glob("*_*.h5"))
    boot_files = sorted(boot_dir.glob("*_*.h5"))

    if not base_files:
        raise FileNotFoundError(f"No per-gene h5 files found in {base_dir}")

    # Extract gene indices and names from filenames
    gene_indices = []
    gene_file_names = []
    for f in base_files:
        parts = f.stem.split("_", 1)  # "0003_FBgn0025702" → ["0003", "FBgn0025702"]
        gene_indices.append(int(parts[0]))
        gene_file_names.append(parts[1] if len(parts) > 1 else f"gene_{parts[0]}")

    n_ref_genes = len(gene_indices)
    print(f"Found {n_ref_genes} per-gene base_correlations files")
    print(f"Found {len(boot_files)} per-gene bootstrap_significant files")

    # Read n_genes from first file (same in the rest)
    with h5py.File(base_files[0], "r") as h5:
        n_genes = h5["meta"].attrs["n_genes"]
        # Read propagated gene names if available
        if "gene_names" in h5:
            all_gene_names = [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
        else:
            all_gene_names = None

    print(f"Network size: {n_genes} genes, {n_genes * (n_genes - 1) // 2} possible edges")

    # Initialize per-reference-gene summary arrays
    summary = {
        "gene_index": np.array(gene_indices, dtype=np.int32),
        "n_sig_edges": np.zeros(n_ref_genes, dtype=np.int32),
        "n_after_filter": np.zeros(n_ref_genes, dtype=np.int32),
        "n_disappear": np.zeros(n_ref_genes, dtype=np.int32),
        "n_new": np.zeros(n_ref_genes, dtype=np.int32),
        "n_sign_change": np.zeros(n_ref_genes, dtype=np.int32),
        "n_strengthen": np.zeros(n_ref_genes, dtype=np.int32),
        "n_weaken": np.zeros(n_ref_genes, dtype=np.int32),
        "n_edges_low": np.zeros(n_ref_genes, dtype=np.int32),
        "n_edges_high": np.zeros(n_ref_genes, dtype=np.int32),
        "n_edges_diff": np.zeros(n_ref_genes, dtype=np.int32),
        "mean_abs_delta": np.zeros(n_ref_genes, dtype=np.float32),
        "max_abs_delta": np.zeros(n_ref_genes, dtype=np.float32),
    }

    # Process each reference gene
    for ref_idx, (g_idx, g_name) in enumerate(zip(gene_indices, gene_file_names)):
        base_path = base_files[ref_idx]
        boot_path = _find_gene_file(boot_dir, g_idx)

        if boot_path is None:
            print(f"  [{ref_idx+1}/{n_ref_genes}] Gene {g_idx} ({g_name}): no bootstrap file, skipping")
            continue

        if (ref_idx + 1) % max(1, n_ref_genes // 20) == 0 or ref_idx == 0:
            print(f"  [{ref_idx+1}/{n_ref_genes}] Gene {g_idx} ({g_name})...")

        # Use existing load_and_process iterating for each .h5/gene
        results = load_and_process(
            base_h5_path=base_path,
            boot_h5_path=boot_path,
            edge_selection=edge_selection,
            min_effect=min_effect,
            require_ci_exclude_zero=require_ci_exclude_zero,
            corr_threshold=corr_threshold,
        )

        n_sig = results["n_significant"]
        summary["n_sig_edges"][ref_idx] = n_sig

        if n_sig == 0:
            continue

        summary["n_after_filter"][ref_idx] = n_sig

        # Qualitative counts
        qs = results["qual_score"]
        summary["n_disappear"][ref_idx] = int(np.sum(qs == QUAL_DISAPPEAR))
        summary["n_new"][ref_idx] = int(np.sum(qs == QUAL_NEW))
        summary["n_sign_change"][ref_idx] = int(np.sum(qs == QUAL_SIGN_CHANGE))
        summary["n_strengthen"][ref_idx] = int(np.sum(qs == QUAL_STRENGTHEN))
        summary["n_weaken"][ref_idx] = int(np.sum(qs == QUAL_WEAKEN))

        # Topology counts
        topo = compute_all_global_topologies(
            results["gene_i"], results["gene_j"],
            results["r_low"], results["r_high"], results["delta_base"],
            n_genes, corr_threshold,
        )
        summary["n_edges_low"][ref_idx] = topo["low"]["n_edges"]
        summary["n_edges_high"][ref_idx] = topo["high"]["n_edges"]
        summary["n_edges_diff"][ref_idx] = topo["diff"]["n_edges"]

        # Delta stats
        abs_delta = np.abs(results["delta_base"])
        summary["mean_abs_delta"][ref_idx] = float(np.mean(abs_delta))
        summary["max_abs_delta"][ref_idx] = float(np.max(abs_delta))

    # Write summary HDF5
    print(f"\nSaving summary to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out_h5_path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["n_genes"] = n_genes
        meta.attrs["n_ref_genes"] = n_ref_genes
        meta.attrs["edge_selection"] = edge_selection
        meta.attrs["min_effect"] = min_effect
        meta.attrs["corr_threshold"] = corr_threshold
        meta.attrs["require_ci_exclude_zero"] = require_ci_exclude_zero

        if all_gene_names is not None:
            h5.create_dataset("gene_names", data=all_gene_names, dtype=h5py.string_dtype())

        per_gene = h5.create_group("per_gene")
        for key, arr in summary.items():
            per_gene.create_dataset(key, data=arr, compression="gzip")

        # Store reference gene names
        per_gene.create_dataset(
            "gene_name", data=gene_file_names, dtype=h5py.string_dtype()
        )

    # Print summary
    n_with_edges = np.sum(summary["n_sig_edges"] > 0)
    print(f"\n{'='*60}")
    print(f"PER-GENE NETWORK COLLECTION SUMMARY")
    print(f"{'='*60}")
    print(f"  Reference genes processed: {n_ref_genes}")
    print(f"  Reference genes with edges: {n_with_edges}")
    print(f"  Total sig edges across all genes: {summary['n_sig_edges'].sum():,}")
    if n_with_edges > 0:
        active = summary["n_sig_edges"] > 0
        print(f"  Mean sig edges per active gene: {summary['n_sig_edges'][active].mean():.1f}")
        print(f"  Max sig edges: {summary['n_sig_edges'].max():,}")
        print(f"  Mean |Δr| (active genes): {summary['mean_abs_delta'][active].mean():.4f}")

    # Top rewiring genes
    if n_with_edges > 0:
        rewiring = summary["n_disappear"] + summary["n_new"] + summary["n_sign_change"]
        top_idx = np.argsort(rewiring)[::-1][:min(20, n_ref_genes)]
        print(f"\n  Top rewirinwsdg reference genes:")
        print(f"  {'Idx':>6} {'Gene':>20} {'SigEdges':>10} {'Disappear':>10} {'New':>6} {'SignChg':>8}")
        for i in top_idx:
            if summary["n_sig_edges"][i] > 0:
                print(f"  {summary['gene_index'][i]:>6} {gene_file_names[i]:>20} "
                      f"{summary['n_sig_edges'][i]:>10} {summary['n_disappear'][i]:>10} "
                      f"{summary['n_new'][i]:>6} {summary['n_sign_change'][i]:>8}")

    # Write TSV
    if out_tsv_path:
        print(f"\nSaving TSV to {out_tsv_path}...")
        out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_tsv_path, "w") as f:
            f.write("gene_idx\tgene_id\tn_sig_edges\tn_disappear\tn_new\tn_sign_change\t"
                    "n_strengthen\tn_weaken\tn_edges_low\tn_edges_high\tn_edges_diff\t"
                    "mean_abs_delta\tmax_abs_delta\n")
            sort_idx = np.argsort(summary["n_sig_edges"])[::-1]
            for i in sort_idx:
                if summary["n_sig_edges"][i] == 0:
                    continue
                f.write(f"{summary['gene_index'][i]}\t{gene_file_names[i]}\t"
                       f"{summary['n_sig_edges'][i]}\t{summary['n_disappear'][i]}\t"
                       f"{summary['n_new'][i]}\t{summary['n_sign_change'][i]}\t"
                       f"{summary['n_strengthen'][i]}\t{summary['n_weaken'][i]}\t"
                       f"{summary['n_edges_low'][i]}\t{summary['n_edges_high'][i]}\t"
                       f"{summary['n_edges_diff'][i]}\t"
                       f"{summary['mean_abs_delta'][i]:.4f}\t"
                       f"{summary['max_abs_delta'][i]:.4f}\n")
        print(f"  Saved {n_with_edges} genes to TSV")


# =============================================================================
# Summary Printing
# =============================================================================

def print_summary(results: dict, all_topo: dict, per_gene: dict | None, focus_analysis: dict) -> None:
    """Print comprehensive summary."""
    print("\n" + "=" * 70)
    print("DIFFERENTIAL CO-EXPRESSION NETWORK ANALYSIS SUMMARY")
    print("=" * 70)

    print(f"\n[NETWORK OVERVIEW]")
    print(f"  Total genes: {results['n_genes']:,}")
    print(f"  Total possible edges: {results['n_tests']:,}")
    print(f"  Significant differential edges: {results['n_significant']:,} "
          f"({100 * results['n_significant'] / results['n_tests']:.4f}%)")

    print(f"\n[GLOBAL TOPOLOGY - LOW NETWORK]")
    topo_low = all_topo["low"]
    print(f"  Active nodes: {topo_low['n_nodes_active']:,}, Edges: {topo_low['n_edges']:,}")
    print(f"  Density: {topo_low['density']:.6f}, Avg degree: {topo_low['avg_degree']:.2f}")
    print(f"  Max degree: {topo_low['max_degree']}, Median degree: {topo_low['median_degree']:.1f}")

    print(f"\n[GLOBAL TOPOLOGY - HIGH NETWORK]")
    topo_high = all_topo["high"]
    print(f"  Active nodes: {topo_high['n_nodes_active']:,}, Edges: {topo_high['n_edges']:,}")
    print(f"  Density: {topo_high['density']:.6f}, Avg degree: {topo_high['avg_degree']:.2f}")
    print(f"  Max degree: {topo_high['max_degree']}, Median degree: {topo_high['median_degree']:.1f}")

    print(f"\n[GLOBAL TOPOLOGY - DIFFERENTIAL NETWORK]")
    topo_diff = all_topo["diff"]
    print(f"  Active nodes: {topo_diff['n_nodes_active']:,}, Edges: {topo_diff['n_edges']:,}")
    print(f"  Density: {topo_diff['density']:.6f}, Avg degree: {topo_diff['avg_degree']:.2f}")
    print(f"  Max degree: {topo_diff['max_degree']}, Median degree: {topo_diff['median_degree']:.1f}")

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

    if per_gene is not None:
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
    else:
        print(f"\n[PER-GENE METRICS]")
        print(f"  Skipped (use --calc-per-gene-metrics to enable)")

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
    # Single-file mode
    parser.add_argument(
        "--base-h5", type=str, default=None,
        help="Path to base_correlations.h5 from Stage 2a (single-gene mode).",
    )
    parser.add_argument(
        "--boot-h5", type=str, default=None,
        help="Path to bootstrap_significant.h5 from Stage 2b (single-gene mode).",
    )
    # Per-gene directory mode
    parser.add_argument(
        "--base-dir", type=str, default=None,
        help="Directory of per-gene base_correlations h5 files (per-gene mode).",
    )
    parser.add_argument(
        "--boot-dir", type=str, default=None,
        help="Directory of per-gene bootstrap_significant h5 files (per-gene mode).",
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
        default="sig_differential",
        help="Edge selection mode.",
    )
    parser.add_argument(
        "--min-effect", type=float, default=0.0001,
        help="Minimum |Δr| threshold.",
    )
    parser.add_argument(
        "--corr-threshold", type=float, default=0.0001,
        help="Correlation threshold for qualitative 'present' (default: 0.0001).",
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
    parser.add_argument(
        "--calc-per-gene-metrics", action="store_true",
        help="Calculate per-gene metrics for ALL genes (expensive, default: False).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Dispatch: per-gene directory mode vs single-file mode
    if args.base_dir and args.boot_dir:
        collect_per_gene_networks(
            base_dir=Path(args.base_dir),
            boot_dir=Path(args.boot_dir),
            out_h5_path=Path(args.out_h5),
            out_tsv_path=Path(args.out_tsv) if args.out_tsv else None,
            edge_selection=args.edge_selection,
            min_effect=args.min_effect,
            require_ci_exclude_zero=not args.no_ci_filter,
            corr_threshold=args.corr_threshold,
        )
        return

    if not args.base_h5 or not args.boot_h5:
        raise SystemExit("ERROR: provide --base-h5/--boot-h5 or --base-dir/--boot-dir.")

    # Single-file mode (original behavior)
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
        print("\nNo significant edges after filtering. Writing empty network.")
        out_h5_path = Path(args.out_h5)
        out_h5_path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(out_h5_path, "w") as h5:
            meta = h5.create_group("meta")
            for key in ["n_genes", "n_tests", "n_significant", "edge_selection",
                        "min_effect", "ci_alpha", "fdr_alpha", "k_low", "k_high",
                        "corr_threshold", "require_ci_exclude_zero", "gene_index_used"]:
                meta.attrs[key] = results[key]
            h5.create_group("edges")
            h5.create_group("topology")
            if results.get("gene_names") is not None:
                h5.create_dataset("gene_names", data=results["gene_names"], dtype=h5py.string_dtype())
        print(f"  Saved empty network to {out_h5_path}")
        return

    # Compute global topology for all three networks (low, high, diff)
    print("\nComputing global topology for low/high/diff networks...")
    all_topo = compute_all_global_topologies(
        results["gene_i"], results["gene_j"],
        results["r_low"], results["r_high"], results["delta_base"],
        results["n_genes"], args.corr_threshold
    )

    # Compute per-gene metrics (optional - expensive for large networks)
    per_gene = None
    if args.calc_per_gene_metrics:
        print("Computing per-gene metrics for all genes...")
        per_gene = compute_per_gene_metrics(
            results["gene_i"], results["gene_j"],
            results["qual_score"], results["r_low"], results["r_high"],
            results["n_genes"]
        )
    else:
        print("Skipping per-gene metrics (use --calc-per-gene-metrics to enable)")

    # Focus gene analysis (always computed)
    focus_gene = args.focus_gene if args.focus_gene is not None else results["gene_index_used"]
    print(f"Analyzing focus gene {focus_gene}...")
    focus_analysis = analyze_focus_gene(
        focus_gene,
        results["gene_i"], results["gene_j"],
        results["qual_score"], results["r_low"], results["r_high"],
        results["delta_base"], results["n_genes"]
    )

    # Gene IDs: prefer propagated names from HDF5 chain, fall back to --gene-ids file
    gene_ids = results.get("gene_names")
    if gene_ids is None and args.gene_ids:
        with open(args.gene_ids) as f:
            gene_ids = [line.strip() for line in f]

    # Save results
    save_results(
        results, all_topo, per_gene, focus_analysis,
        out_h5_path=Path(args.out_h5),
        out_tsv_path=Path(args.out_tsv) if args.out_tsv else None,
        gene_ids=gene_ids,
    )

    # Print summary
    print_summary(results, all_topo, per_gene, focus_analysis)


if __name__ == "__main__":
    main()
