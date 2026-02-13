#!/usr/bin/env python3
"""
Stage 3: Reconstruct Differential Co-expression Network with Topological Analysis.

Combines results from Stage 2a (base correlations) and Stage 2b (bootstrap)
to produce final differential network with:
- Quantitative changes (delta = r_high - r_low)
- Qualitative changes (disappear, new, sign_change)
- Global topological features for LOW, HIGH, and DIFF networks
- Focus gene neighborhood analysis (1st and 2nd layer partners)

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
from scipy import stats as scipy_stats


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

def _edge_key(i: int, j: int) -> tuple:
    """Return canonical edge key (smaller index first)."""
    return (min(i, j), max(i, j))


def compute_scale_free_metrics(degrees: np.ndarray) -> dict:
    """
    Test if network follows scale-free topology.
    Biological networks typically have gamma = 2-3.
    """
    active_degrees = degrees[degrees > 0]
    if len(active_degrees) < 10:
        return {"power_law_exponent": float("nan"), "power_law_r_squared": float("nan"), "is_scale_free": False}

    unique_degrees, counts = np.unique(active_degrees, return_counts=True)
    mask = counts >= 2
    if mask.sum() < 5:
        return {"power_law_exponent": float("nan"), "power_law_r_squared": float("nan"), "is_scale_free": False}

    log_k = np.log10(unique_degrees[mask].astype(float))
    log_p = np.log10(counts[mask].astype(float))

    slope, intercept, r_value, _, _ = scipy_stats.linregress(log_k, log_p)
    gamma = -slope

    is_scale_free = (2.0 <= gamma <= 3.5) and (r_value**2 > 0.8)

    return {
        "power_law_exponent": float(gamma),
        "power_law_r_squared": float(r_value**2),
        "is_scale_free": bool(is_scale_free),
    }


def compute_clustering_coefficient(adj: dict, n_genes: int) -> float:
    """
    Compute global clustering coefficient (average of local clustering).
    High clustering = functional modules/pathways.
    """
    local_clustering = []
    for node in range(n_genes):
        neighbors = list(adj.get(node, set()))
        k = len(neighbors)
        if k < 2:
            continue
        triangles = sum(
            1 for i_idx, n1 in enumerate(neighbors)
            for n2 in neighbors[i_idx + 1:]
            if n2 in adj.get(n1, set())
        )
        possible = k * (k - 1) / 2
        local_clustering.append(triangles / possible)

    return float(np.mean(local_clustering)) if local_clustering else 0.0


def compute_connected_components(adj: dict, n_genes: int) -> dict:
    """Find disconnected subnetworks using BFS."""
    from collections import deque

    visited = set()
    components = []

    for start in range(n_genes):
        if start in visited or len(adj.get(start, set())) == 0:
            continue
        component = set()
        queue = deque([start])
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            component.add(node)
            queue.extend(n for n in adj.get(node, set()) if n not in visited)
        components.append(len(component))

    components = sorted(components, reverse=True)
    return {
        "n_components": len(components),
        "largest_component_size": components[0] if components else 0,
        "largest_component_fraction": float(components[0] / n_genes) if components else 0.0,
    }


def compute_assortativity(gene_i: np.ndarray, gene_j: np.ndarray, degrees: np.ndarray) -> float:
    """
    Pearson correlation of degrees at edge endpoints.
    r > 0: hubs connect to hubs (assortative)
    r < 0: hubs connect to periphery (disassortative)
    """
    if len(gene_i) < 2:
        return 0.0
    deg_i = degrees[gene_i].astype(float)
    deg_j = degrees[gene_j].astype(float)
    if np.std(deg_i) == 0 or np.std(deg_j) == 0:
        return 0.0
    r, _ = scipy_stats.pearsonr(deg_i, deg_j)
    return float(r)


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

    _empty = {
        "n_nodes_total": n_genes,
        "n_nodes_active": 0,
        "fraction_active": 0.0,
        "n_edges": 0,
        "density": 0.0,
        "avg_degree": 0.0,
        "max_degree": 0,
        "min_degree": 0,
        "median_degree": 0.0,
        "std_degree": 0.0,
        "degree_p25": 0.0, "degree_p50": 0.0, "degree_p75": 0.0,
        "degree_p90": 0.0, "degree_p95": 0.0, "degree_p99": 0.0,
        "power_law_exponent": float("nan"),
        "power_law_r_squared": float("nan"),
        "is_scale_free": False,
        "global_clustering": 0.0,
        "n_components": 0,
        "largest_component_size": 0,
        "largest_component_fraction": 0.0,
        "assortativity": 0.0,
        "degrees": np.zeros(n_genes, dtype=np.int32),
    }
    if n_edges == 0:
        return _empty

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

    # Basic metrics
    density = 2 * n_edges / (n_genes * (n_genes - 1)) if n_genes > 1 else 0
    avg_degree = 2 * n_edges / n_genes if n_genes > 0 else 0

    active_degrees = degrees[degrees > 0]
    max_degree = int(degrees.max()) if len(degrees) > 0 else 0
    min_degree = int(active_degrees.min()) if len(active_degrees) > 0 else 0
    median_degree = float(np.median(active_degrees)) if len(active_degrees) > 0 else 0
    std_degree = float(active_degrees.std()) if len(active_degrees) > 0 else 0

    # Degree percentiles
    percentiles = {}
    if len(active_degrees) > 0:
        for p in [25, 50, 75, 90, 95, 99]:
            percentiles[f"degree_p{p}"] = float(np.percentile(active_degrees, p))
    else:
        for p in [25, 50, 75, 90, 95, 99]:
            percentiles[f"degree_p{p}"] = 0.0

    # Scale-free properties
    scale_free = compute_scale_free_metrics(degrees)

    # Clustering coefficient
    global_clustering = compute_clustering_coefficient(adj, n_genes)

    # Connected components
    components = compute_connected_components(adj, n_genes)

    # Assortativity
    assortativity = compute_assortativity(gi, gj, degrees)

    return {
        "n_nodes_total": n_genes,
        "n_nodes_active": n_nodes,
        "fraction_active": float(n_nodes / n_genes) if n_genes > 0 else 0.0,
        "n_edges": n_edges,
        "density": density,
        "avg_degree": avg_degree,
        "max_degree": max_degree,
        "min_degree": min_degree,
        "median_degree": median_degree,
        "std_degree": std_degree,
        **percentiles,
        "power_law_exponent": scale_free["power_law_exponent"],
        "power_law_r_squared": scale_free["power_law_r_squared"],
        "is_scale_free": scale_free["is_scale_free"],
        "global_clustering": global_clustering,
        "n_components": components["n_components"],
        "largest_component_size": components["largest_component_size"],
        "largest_component_fraction": components["largest_component_fraction"],
        "assortativity": assortativity,
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
    sum_r_low = np.zeros(n_genes, dtype=np.float32)
    sum_r_high = np.zeros(n_genes, dtype=np.float32)

    for idx, (i, j) in enumerate(zip(gene_i, gene_j)):
        # Both genes get this edge
        for g in [i, j]:
            degree[g] += 1

            # Count by presence in low/high
            if np.abs(r_low[idx]) >= corr_threshold:
                degree_low[g] += 1
                sum_r_low[g] += np.abs(r_low[idx])
            if np.abs(r_high[idx]) >= corr_threshold:
                degree_high[g] += 1
                sum_r_high[g] += np.abs(r_high[idx])

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

    # Mean metrics (avoid divide-by-zero)
    mean_r_low = np.where(degree_low > 0, sum_r_low / degree_low, 0.0).astype(np.float32)
    mean_r_high = np.where(degree_high > 0, sum_r_high / degree_high, 0.0).astype(np.float32)
    mean_abs_delta = np.where(degree > 0, sum_abs_delta / degree, 0.0).astype(np.float32)

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
        "sum_r_low": sum_r_low,
        "sum_r_high": sum_r_high,
        "mean_r_low": mean_r_low,
        "mean_r_high": mean_r_high,
        "mean_abs_delta": mean_abs_delta,
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
    sig_low: np.ndarray | None = None,
    sig_high: np.ndarray | None = None,
    corr_threshold: float = 0.0001,
) -> dict:
    """
    Analyze neighborhood of the focus gene (1st and 2nd layer)
    across all three networks: low, high, and differential.

    Returns
    -------
    results : dict with direct partners, indirect partners, statistics,
              low/high network stats, and L1/L2 two-layer metrics with ratios
    """
    # Build THREE adjacency lists + edge map
    adj_diff = defaultdict(set)  # differential network (all sig edges)
    adj_low = defaultdict(set)   # low network (|r_low| >= threshold)
    adj_high = defaultdict(set)  # high network (|r_high| >= threshold)
    edge_map = {}  # (i, j) -> edge index
    r_low_map = {}  # edge_key -> r_low value
    r_high_map = {}  # edge_key -> r_high value
    delta_map = {}  # edge_key -> delta value

    for idx, (i, j) in enumerate(zip(gene_i, gene_j)):
        adj_diff[i].add(j)
        adj_diff[j].add(i)
        key = _edge_key(i, j)
        edge_map[key] = idx
        r_low_map[key] = r_low[idx]
        r_high_map[key] = r_high[idx]
        delta_map[key] = delta[idx]

        # "Present" = significant AND |r| >= threshold (consistent with classify_qualitative_change)
        is_sig_low = sig_low[idx] if sig_low is not None else True
        is_sig_high = sig_high[idx] if sig_high is not None else True
        if is_sig_low and abs(r_low[idx]) >= corr_threshold:
            adj_low[i].add(j)
            adj_low[j].add(i)
        if is_sig_high and abs(r_high[idx]) >= corr_threshold:
            adj_high[i].add(j)
            adj_high[j].add(i)

    # 1st layer: direct partners (in diff network)
    direct_partners = sorted(adj_diff[focus_gene])

    # 2nd layer: partners of partners (excluding focus gene and direct partners)
    direct_set = set(direct_partners)
    indirect_partners = set()
    for partner in direct_partners:
        for second in adj_diff[partner]:
            if second != focus_gene and second not in direct_set:
                indirect_partners.add(second)
    indirect_partners = sorted(indirect_partners)

    # Get edges involving focus gene (1st layer edges)
    direct_edges = []
    for partner in direct_partners:
        key = _edge_key(focus_gene, partner)
        if key in edge_map:
            direct_edges.append(edge_map[key])
    direct_edges = np.array(direct_edges, dtype=np.int32)

    # Get edges in 2-layer neighborhood (1st + edges between 1st layer)
    two_layer_edges = set(direct_edges.tolist())
    for i_idx, p1 in enumerate(direct_partners):
        for p2 in direct_partners[i_idx + 1:]:
            key = _edge_key(p1, p2)
            if key in edge_map:
                two_layer_edges.add(edge_map[key])
    two_layer_edges = np.array(sorted(two_layer_edges), dtype=np.int32)

    # Get FULL two-layer edge set (focus→L1 + L1↔L1 + L1→L2)
    # This extends two_layer_edges with edges from L1 partners to L2-only partners
    indirect_set = set(indirect_partners)
    full_two_layer_edges = set(two_layer_edges.tolist())
    for p1 in direct_partners:
        for p2 in indirect_set:
            key = _edge_key(p1, p2)
            if key in edge_map:
                full_two_layer_edges.add(edge_map[key])
    full_two_layer_edges = np.array(sorted(full_two_layer_edges), dtype=np.int32)

    # Statistics for edge subsets (L1, L2, full two-layer)
    def compute_edge_stats(edge_indices):
        if len(edge_indices) == 0:
            return {
                "n_edges": 0,
                "n_disappear": 0, "n_new": 0, "n_sign_change": 0,
                "n_strengthen": 0, "n_weaken": 0,
                "mean_abs_delta": 0.0,
                "sum_abs_delta": 0.0, "str_low": 0.0, "str_high": 0.0,
            }
        qs = qual_score[edge_indices]
        rl = np.abs(r_low[edge_indices])
        rh = np.abs(r_high[edge_indices])

        # Qualitative counts from qual_score
        n_disappear = int(np.sum(qs == QUAL_DISAPPEAR))
        n_new = int(np.sum(qs == QUAL_NEW))
        n_sign_change = int(np.sum(qs == QUAL_SIGN_CHANGE))

        # Strengthen/weaken INDEPENDENT of sign_change:
        # For all "both present" edges, compare |r_high| vs |r_low|
        mask_both = ((qs == QUAL_SIGN_CHANGE) | (qs == QUAL_STRENGTHEN) | (qs == QUAL_WEAKEN))
        n_strengthen = int(np.sum(mask_both & (rh > rl)))
        n_weaken = int(np.sum(mask_both & (rh < rl)))

        # str_low / str_high: sum|r| only for "present" edges
        if sig_low is not None:
            present_low = sig_low[edge_indices] & (rl >= corr_threshold)
            present_high = sig_high[edge_indices] & (rh >= corr_threshold)
        else:
            present_low = rl >= corr_threshold
            present_high = rh >= corr_threshold
        str_low = float(np.sum(rl[present_low]))
        str_high = float(np.sum(rh[present_high]))

        return {
            "n_edges": len(edge_indices),
            "n_disappear": n_disappear,
            "n_new": n_new,
            "n_sign_change": n_sign_change,
            "n_strengthen": n_strengthen,
            "n_weaken": n_weaken,
            "mean_abs_delta": float(np.mean(np.abs(delta[edge_indices]))),
            "sum_abs_delta": float(np.sum(np.abs(delta[edge_indices]))),
            "str_low": str_low,
            "str_high": str_high,
        }

    direct_stats = compute_edge_stats(direct_edges)
    two_layer_stats = compute_edge_stats(two_layer_edges)
    full_two_layer_stats = compute_edge_stats(full_two_layer_edges)

    # --- L1/L2 degree in diff network ---
    L1_deg_diff = len(adj_diff.get(focus_gene, set()))
    # L2 = L1 + partners-of-L1-partners (combined neighborhood size)
    L1_set = set(adj_diff.get(focus_gene, set()))
    L2_set = set()
    for p in L1_set:
        L2_set.update(adj_diff.get(p, set()))
    L2_set.discard(focus_gene)
    L2_deg_diff = len(L2_set)  # all unique nodes in L1+L2

    # --- Rewiring scores (disappear + new + sign_change) ---
    L1_rewire = (direct_stats["n_disappear"] + direct_stats["n_new"]
                 + direct_stats["n_sign_change"])
    L2_rewire = (full_two_layer_stats["n_disappear"] + full_two_layer_stats["n_new"]
                 + full_two_layer_stats["n_sign_change"])

    # --- Strength metrics ---
    L1_str_low = direct_stats["str_low"]
    L1_str_high = direct_stats["str_high"]
    L1_str_diff = direct_stats["sum_abs_delta"]
    L1_mean_abs_dr = direct_stats["mean_abs_delta"]
    L2_str_low = full_two_layer_stats["str_low"]
    L2_str_high = full_two_layer_stats["str_high"]
    L2_str_diff = full_two_layer_stats["sum_abs_delta"]

    # --- L2/L1 expansion ratios (diff network) ---
    L2L1_deg = L2_deg_diff / L1_deg_diff if L1_deg_diff > 0 else 0.0
    L2L1_rewire = L2_rewire / L1_rewire if L1_rewire > 0 else 0.0
    L2L1_str = L2_str_diff / L1_str_diff if L1_str_diff > 0 else 0.0

    # --- High/Low condition strength ratios ---
    HL_str_L1 = L1_str_high / L1_str_low if L1_str_low > 0 else 0.0
    HL_str_L2 = L2_str_high / L2_str_low if L2_str_low > 0 else 0.0

    # Per-partner breakdown for direct partners
    partner_details = []
    for partner in direct_partners:
        key = _edge_key(focus_gene, partner)
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
        "full_two_layer_stats": full_two_layer_stats,
        "partner_details": partner_details,
        # --- Flat metrics for TSV collection ---
        "L1_deg_diff": L1_deg_diff,
        "L1_n_disappear": direct_stats["n_disappear"],
        "L1_n_new": direct_stats["n_new"],
        "L1_n_sign_chg": direct_stats["n_sign_change"],
        "L1_n_strengthen": direct_stats["n_strengthen"],
        "L1_n_weaken": direct_stats["n_weaken"],
        "L1_rewire": L1_rewire,
        "L1_str_low": float(L1_str_low),
        "L1_str_high": float(L1_str_high),
        "L1_str_diff": float(L1_str_diff),
        "L1_mean_abs_dr": float(L1_mean_abs_dr),
        "L2_deg_diff": L2_deg_diff,
        "L2_rewire": L2_rewire,
        "L2_str_low": float(L2_str_low),
        "L2_str_high": float(L2_str_high),
        "L2_str_diff": float(L2_str_diff),
        "L2L1_deg": float(L2L1_deg),
        "L2L1_rewire": float(L2L1_rewire),
        "L2L1_str": float(L2L1_str),
        "HL_str_L1": float(HL_str_L1),
        "HL_str_L2": float(HL_str_L2),
    }


# =============================================================================
# Data Loading and Filtering
# =============================================================================

def load_and_process(
    base_h5_path: Path,
    boot_h5_path: Path,
    edge_selection: str = "sig_differential",
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
            # Read gene_index_used with warning if missing
            if "meta" in h5 and "gene_index_used" in h5["meta"].attrs:
                gene_index_used = h5["meta"].attrs["gene_index_used"]
            else:
                print("  WARNING: gene_index_used not found in bootstrap file!")
                print("  Defaulting to gene_index=0. This may indicate a Stage 2b issue.")
                gene_index_used = 0
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
            # Read gene_index_used with warning if missing
            if "meta" in h5 and "gene_index_used" in h5["meta"].attrs:
                gene_index_used = h5["meta"].attrs["gene_index_used"]
            else:
                print("  WARNING: gene_index_used not found in bootstrap file!")
                print("  Defaulting to gene_index=0. This may indicate a Stage 2b issue.")
                gene_index_used = 0

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
    focus_analysis: dict,
    out_h5_path: Path,
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
        global_topo_attrs = [
            "n_nodes_total", "n_nodes_active", "fraction_active",
            "n_edges", "density",
            "avg_degree", "max_degree", "min_degree", "median_degree", "std_degree",
            "degree_p25", "degree_p50", "degree_p75", "degree_p90", "degree_p95", "degree_p99",
            "power_law_exponent", "power_law_r_squared", "is_scale_free",
            "global_clustering",
            "n_components", "largest_component_size", "largest_component_fraction",
            "assortativity",
        ]
        for net_name in ["low", "high", "diff"]:
            net_topo = all_topo[net_name]
            grp = topo.create_group(f"global_{net_name}")
            for key in global_topo_attrs:
                val = net_topo.get(key)
                if val is not None:
                    # Handle NaN for HDF5 attrs
                    if isinstance(val, float) and np.isnan(val):
                        grp.attrs[key] = float("nan")
                    else:
                        grp.attrs[key] = val
            # Also save degree array for each network
            grp.create_dataset("degrees", data=net_topo["degrees"], compression="gzip")

        # Focus gene analysis
        if focus_analysis is not None:
            focus = h5.create_group("focus_gene")
            focus.attrs["gene_index"] = focus_analysis["focus_gene"]
            focus.attrs["n_direct_partners"] = focus_analysis["n_direct_partners"]
            focus.attrs["n_indirect_partners"] = focus_analysis["n_indirect_partners"]

            focus.create_dataset("direct_partners", data=focus_analysis["direct_partners"])
            focus.create_dataset("indirect_partners", data=focus_analysis["indirect_partners"])

            # Direct stats (L1)
            direct_grp = focus.create_group("direct_stats")
            for key, val in focus_analysis["direct_stats"].items():
                direct_grp.attrs[key] = val

            # Two-layer stats (L1+L1↔L1)
            two_layer_grp = focus.create_group("two_layer_stats")
            for key, val in focus_analysis["two_layer_stats"].items():
                two_layer_grp.attrs[key] = val

            # Flat metrics (L1/L2/ratios)
            flat_keys = [k for k in focus_analysis if k.startswith(("L1_", "L2_", "L2L1_", "HL_"))]
            metrics_grp = focus.create_group("metrics")
            for key in flat_keys:
                metrics_grp.attrs[key] = focus_analysis[key]

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
    out_focus_tsv_path: Path | None = None,
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
        "n_sig_total": np.zeros(n_ref_genes, dtype=np.int32),
        "n_after_filter": np.zeros(n_ref_genes, dtype=np.int32),
        # Whole-network qualitative counts
        "n_disappear": np.zeros(n_ref_genes, dtype=np.int32),
        "n_new": np.zeros(n_ref_genes, dtype=np.int32),
        "n_sign_change": np.zeros(n_ref_genes, dtype=np.int32),
        "n_strengthen": np.zeros(n_ref_genes, dtype=np.int32),
        "n_weaken": np.zeros(n_ref_genes, dtype=np.int32),
        "mean_abs_delta": np.zeros(n_ref_genes, dtype=np.float32),
        "max_abs_delta": np.zeros(n_ref_genes, dtype=np.float32),
        # Focus gene L1 metrics
        "L1_deg_diff": np.zeros(n_ref_genes, dtype=np.int32),
        "L1_n_disappear": np.zeros(n_ref_genes, dtype=np.int32),
        "L1_n_new": np.zeros(n_ref_genes, dtype=np.int32),
        "L1_n_sign_chg": np.zeros(n_ref_genes, dtype=np.int32),
        "L1_n_strengthen": np.zeros(n_ref_genes, dtype=np.int32),
        "L1_n_weaken": np.zeros(n_ref_genes, dtype=np.int32),
        "L1_rewire": np.zeros(n_ref_genes, dtype=np.int32),
        "L1_str_low": np.zeros(n_ref_genes, dtype=np.float32),
        "L1_str_high": np.zeros(n_ref_genes, dtype=np.float32),
        "L1_str_diff": np.zeros(n_ref_genes, dtype=np.float32),
        "L1_mean_abs_dr": np.zeros(n_ref_genes, dtype=np.float32),
        # Focus gene L2 metrics
        "L2_deg_diff": np.zeros(n_ref_genes, dtype=np.int32),
        "L2_rewire": np.zeros(n_ref_genes, dtype=np.int32),
        "L2_str_low": np.zeros(n_ref_genes, dtype=np.float32),
        "L2_str_high": np.zeros(n_ref_genes, dtype=np.float32),
        "L2_str_diff": np.zeros(n_ref_genes, dtype=np.float32),
        # L2/L1 expansion ratios (diff network)
        "L2L1_deg": np.zeros(n_ref_genes, dtype=np.float32),
        "L2L1_rewire": np.zeros(n_ref_genes, dtype=np.float32),
        "L2L1_str": np.zeros(n_ref_genes, dtype=np.float32),
        # High/Low condition strength ratios
        "HL_str_L1": np.zeros(n_ref_genes, dtype=np.float32),
        "HL_str_L2": np.zeros(n_ref_genes, dtype=np.float32),
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
        try:
            results = load_and_process(
                base_h5_path=base_path,
                boot_h5_path=boot_path,
                edge_selection=edge_selection,
                min_effect=min_effect,
                require_ci_exclude_zero=require_ci_exclude_zero,
                corr_threshold=corr_threshold,
            )
        except (ValueError, KeyError) as e:
            print(f"  [{ref_idx+1}/{n_ref_genes}] Gene {g_idx} ({g_name}): ERROR - {e}, skipping")
            continue

        n_sig = results["n_significant"]
        summary["n_sig_total"][ref_idx] = n_sig

        if n_sig == 0:
            continue

        summary["n_after_filter"][ref_idx] = n_sig

        # Whole-network qualitative counts
        qs = results["qual_score"]
        summary["n_disappear"][ref_idx] = int(np.sum(qs == QUAL_DISAPPEAR))
        summary["n_new"][ref_idx] = int(np.sum(qs == QUAL_NEW))
        summary["n_sign_change"][ref_idx] = int(np.sum(qs == QUAL_SIGN_CHANGE))
        summary["n_strengthen"][ref_idx] = int(np.sum(qs == QUAL_STRENGTHEN))
        summary["n_weaken"][ref_idx] = int(np.sum(qs == QUAL_WEAKEN))

        # Delta stats
        abs_delta = np.abs(results["delta_base"])
        summary["mean_abs_delta"][ref_idx] = float(np.mean(abs_delta))
        summary["max_abs_delta"][ref_idx] = float(np.max(abs_delta))

        # Focus gene analysis (reference gene = focus gene)
        fa = analyze_focus_gene(
            focus_gene=g_idx,
            gene_i=results["gene_i"],
            gene_j=results["gene_j"],
            qual_score=results["qual_score"],
            r_low=results["r_low"],
            r_high=results["r_high"],
            delta=results["delta_base"],
            n_genes=n_genes,
            sig_low=results["sig_low"],
            sig_high=results["sig_high"],
            corr_threshold=corr_threshold,
        )

        # Copy flat metrics from analyze_focus_gene into summary arrays
        for key in ["L1_deg_diff", "L1_n_disappear", "L1_n_new", "L1_n_sign_chg",
                     "L1_n_strengthen", "L1_n_weaken", "L1_rewire",
                     "L1_str_low", "L1_str_high", "L1_str_diff", "L1_mean_abs_dr",
                     "L2_deg_diff", "L2_rewire", "L2_str_low", "L2_str_high", "L2_str_diff",
                     "L2L1_deg", "L2L1_rewire", "L2L1_str",
                     "HL_str_L1", "HL_str_L2"]:
            summary[key][ref_idx] = fa[key]

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
    n_with_edges = np.sum(summary["n_sig_total"] > 0)
    print(f"\n{'='*60}")
    print(f"PER-GENE NETWORK COLLECTION SUMMARY")
    print(f"{'='*60}")
    print(f"  Reference genes processed: {n_ref_genes}")
    print(f"  Reference genes with edges: {n_with_edges}")
    print(f"  Total sig edges across all genes: {summary['n_sig_total'].sum():,}")
    if n_with_edges > 0:
        active = summary["n_sig_total"] > 0
        print(f"  Mean sig edges per active gene: {summary['n_sig_total'][active].mean():.1f}")
        print(f"  Max sig edges: {summary['n_sig_total'].max():,}")
        print(f"  Mean |Δr| (active genes): {summary['mean_abs_delta'][active].mean():.4f}")

    # Top genes by L2L1_deg
    if n_with_edges > 0:
        top_idx = np.argsort(summary["L2L1_deg"])[::-1][:min(20, n_ref_genes)]
        print(f"\n  Top genes by L2/L1 degree ratio (diff network):")
        print(f"  {'Idx':>6} {'Gene':>20} {'SigTotal':>10} {'L1_deg':>7} {'L2_deg':>7} "
              f"{'L2L1':>7} {'L1_rew':>7} {'L2_rew':>7}")
        for i in top_idx:
            if summary["n_sig_total"][i] > 0:
                print(f"  {summary['gene_index'][i]:>6} {gene_file_names[i]:>20} "
                      f"{summary['n_sig_total'][i]:>10} {summary['L1_deg_diff'][i]:>7} "
                      f"{summary['L2_deg_diff'][i]:>7} {summary['L2L1_deg'][i]:>7.2f} "
                      f"{summary['L1_rewire'][i]:>7} {summary['L2_rewire'][i]:>7}")

    # Write focus gene TSV (sorted by L2L1_deg descending)
    if out_focus_tsv_path:
        print(f"\nSaving focus gene TSV to {out_focus_tsv_path}...")
        out_focus_tsv_path.parent.mkdir(parents=True, exist_ok=True)

        # TSV columns (order matches revised plan)
        focus_tsv_cols = [
            "gene_idx", "gene_id", "n_sig_total",
            "L1_deg_diff",
            "L1_n_disappear", "L1_n_new", "L1_n_sign_chg",
            "L1_n_strengthen", "L1_n_weaken", "L1_rewire",
            "L1_str_low", "L1_str_high", "L1_str_diff", "L1_mean_abs_dr",
            "L2_deg_diff", "L2_rewire",
            "L2_str_low", "L2_str_high", "L2_str_diff",
            "L2L1_deg", "L2L1_rewire", "L2L1_str",
            "HL_str_L1", "HL_str_L2",
        ]

        # Sort by L2L1_deg descending
        sort_idx = np.argsort(summary["L2L1_deg"])[::-1]

        # Integer columns (no decimal formatting)
        int_cols = {"gene_idx", "n_sig_total", "L1_deg_diff",
                    "L1_n_disappear", "L1_n_new", "L1_n_sign_chg",
                    "L1_n_strengthen", "L1_n_weaken", "L1_rewire",
                    "L2_deg_diff", "L2_rewire"}

        with open(out_focus_tsv_path, "w") as f:
            f.write("\t".join(focus_tsv_cols) + "\n")
            for i in sort_idx:
                if summary["n_sig_total"][i] == 0:
                    continue
                vals = []
                for col in focus_tsv_cols:
                    if col == "gene_idx":
                        vals.append(str(summary["gene_index"][i]))
                    elif col == "gene_id":
                        vals.append(gene_file_names[i])
                    elif col in int_cols:
                        vals.append(str(summary[col][i]))
                    else:
                        vals.append(f"{summary[col][i]:.4f}")
                f.write("\t".join(vals) + "\n")
        print(f"  Saved {n_with_edges} genes to focus gene TSV")


# =============================================================================
# Summary Printing
# =============================================================================

def _fmt_val(val) -> str:
    """Format a metric value for display."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return "N/A"
    elif isinstance(val, bool):
        return "Yes" if val else "No"
    elif isinstance(val, float):
        if abs(val) < 0.01 and val != 0:
            return f"{val:.2e}"
        return f"{val:.4f}"
    elif isinstance(val, (int, np.integer)):
        return f"{val:,}"
    return str(val)


def _print_network_topology(label: str, t: dict) -> None:
    """Print comprehensive topology for one network."""
    print(f"\n[GLOBAL TOPOLOGY - {label}]")
    print(f"  Nodes: {t['n_nodes_active']:,} / {t['n_nodes_total']:,} "
          f"({100*t.get('fraction_active', 0):.1f}% active)")
    print(f"  Edges: {t['n_edges']:,}, Density: {t['density']:.6f}")
    print(f"  Degree: avg={t['avg_degree']:.2f}, median={t['median_degree']:.1f}, "
          f"max={t['max_degree']}, std={t.get('std_degree', 0):.2f}")
    clustering = t.get("global_clustering", 0)
    print(f"  Clustering: {clustering:.4f}", end="")
    if clustering > 0.3:
        print(" (highly modular)")
    elif clustering > 0.1:
        print(" (moderate)")
    else:
        print(" (sparse)")
    n_comp = t.get("n_components", 0)
    lc_frac = t.get("largest_component_fraction", 0)
    print(f"  Components: {n_comp}, largest={100*lc_frac:.1f}% of genes")
    gamma = t.get("power_law_exponent", float("nan"))
    r2 = t.get("power_law_r_squared", float("nan"))
    sf = t.get("is_scale_free", False)
    print(f"  Scale-free: {_fmt_val(sf)} (gamma={_fmt_val(gamma)}, R2={_fmt_val(r2)})")
    assort = t.get("assortativity", 0)
    print(f"  Assortativity: {assort:+.4f}", end="")
    if assort > 0.3:
        print(" (hubs connect to hubs)")
    elif assort < -0.3:
        print(" (hubs connect to periphery)")
    else:
        print(" (neutral)")


def print_summary(results: dict, all_topo: dict, focus_analysis: dict) -> None:
    """Print comprehensive summary."""
    print("\n" + "=" * 70)
    print("DIFFERENTIAL CO-EXPRESSION NETWORK ANALYSIS SUMMARY")
    print("=" * 70)

    print(f"\n[NETWORK OVERVIEW]")
    print(f"  Total genes: {results['n_genes']:,}")
    print(f"  Total possible edges: {results['n_tests']:,}")
    print(f"  Significant differential edges: {results['n_significant']:,} "
          f"({100 * results['n_significant'] / results['n_tests']:.4f}%)")

    _print_network_topology("LOW NETWORK", all_topo["low"])
    _print_network_topology("HIGH NETWORK", all_topo["high"])
    _print_network_topology("DIFFERENTIAL NETWORK", all_topo["diff"])

    print(f"\n[QUALITATIVE CHANGES] (low -> high)")
    for code, label in QUAL_LABELS.items():
        count = np.sum(results["qual_score"] == code)
        if count > 0:
            pct = 100 * count / results["n_significant"]
            print(f"  {label:12s}: {count:>8,} ({pct:5.1f}%)")

    print(f"\n[QUANTITATIVE CHANGES]")
    print(f"  Mean dr: {np.mean(results['delta_base']):+.4f}")
    print(f"  Mean |dr|: {np.mean(np.abs(results['delta_base'])):.4f}")
    print(f"  Max |dr|: {np.max(np.abs(results['delta_base'])):.4f}")
    print(f"  Mean bootstrap bias: {np.mean(results['bias']):.4f}")

    if focus_analysis is not None:
        fg = focus_analysis["focus_gene"]
        print(f"\n[FOCUS GENE ANALYSIS] (gene {fg})")
        print(f"  Direct partners (1st layer): {focus_analysis['n_direct_partners']}")
        print(f"  Indirect partners (2nd layer): {focus_analysis['n_indirect_partners']}")

        ds = focus_analysis["direct_stats"]
        fa = focus_analysis
        print(f"\n  L1 (Direct) Metrics:")
        print(f"    Degree (diff): {fa['L1_deg_diff']}")
        print(f"    Disappear: {ds['n_disappear']}, New: {ds['n_new']}, "
              f"Sign change: {ds['n_sign_change']}")
        print(f"    Strengthen: {ds['n_strengthen']}, Weaken: {ds['n_weaken']}")
        print(f"    Rewiring: {fa['L1_rewire']}")
        print(f"    str_low: {fa['L1_str_low']:.4f}, str_high: {fa['L1_str_high']:.4f}, "
              f"str_diff: {fa['L1_str_diff']:.4f}")
        print(f"    Mean |Δr|: {fa['L1_mean_abs_dr']:.4f}")

        print(f"\n  L2 (Full Two-Layer) Metrics:")
        print(f"    Degree (diff): {fa['L2_deg_diff']}")
        print(f"    Rewiring: {fa['L2_rewire']}")
        print(f"    str_low: {fa['L2_str_low']:.4f}, str_high: {fa['L2_str_high']:.4f}, "
              f"str_diff: {fa['L2_str_diff']:.4f}")

        print(f"\n  Ratios:")
        print(f"    L2/L1 deg: {fa['L2L1_deg']:.2f}, rewire: {fa['L2L1_rewire']:.2f}, "
              f"str: {fa['L2L1_str']:.2f}")
        print(f"    H/L str L1: {fa['HL_str_L1']:.4f}, L2: {fa['HL_str_L2']:.4f}")

        if focus_analysis["partner_details"]:
            print(f"\n  Direct Partner Details:")
            print(f"    {'Partner':>8} {'r_low':>8} {'r_high':>8} {'dr':>8} {'Category':>12}")
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
        "--out-focus-tsv", type=str, default=None,
        help="Output TSV path for focus gene metrics table (per-gene mode).",
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Dispatch: per-gene directory mode vs single-file mode
    if args.base_dir and args.boot_dir:
        collect_per_gene_networks(
            base_dir=Path(args.base_dir),
            boot_dir=Path(args.boot_dir),
            out_h5_path=Path(args.out_h5),
            out_focus_tsv_path=Path(args.out_focus_tsv) if args.out_focus_tsv else None,
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

    # Focus gene analysis (always computed)
    focus_gene = args.focus_gene if args.focus_gene is not None else results["gene_index_used"]
    print(f"Analyzing focus gene {focus_gene}...")
    focus_analysis = analyze_focus_gene(
        focus_gene,
        results["gene_i"], results["gene_j"],
        results["qual_score"], results["r_low"], results["r_high"],
        results["delta_base"], results["n_genes"],
        sig_low=results["sig_low"],
        sig_high=results["sig_high"],
        corr_threshold=args.corr_threshold,
    )

    # Gene IDs: prefer propagated names from HDF5 chain, fall back to --gene-ids file
    gene_ids = results.get("gene_names")
    if gene_ids is None and args.gene_ids:
        with open(args.gene_ids) as f:
            gene_ids = [line.strip() for line in f]

    # Save results
    save_results(
        results, all_topo, focus_analysis,
        out_h5_path=Path(args.out_h5),
        gene_ids=gene_ids,
    )

    # Print summary
    print_summary(results, all_topo, focus_analysis)


if __name__ == "__main__":
    main()
