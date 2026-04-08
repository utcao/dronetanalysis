#!/usr/bin/env python3
"""
Stage 3: Per-gene L1/L2 neighbourhood metrics (single co-expression network).

Runs in single-gene mode, one Snakemake job per gene (rule per_gene_metrics).
Uses the CSR adjacency in network_edges.h5 to load only O(degree + L1_degrees)
bytes per job — no full edge list read required.

Differences from analyze_focus_gene() in the differential pipeline:
  - Single adjacency only (not three: low/high/diff).
  - No qualitative change categories (disappear/new/sign_change/strengthen/weaken).
  - No rewiring metrics — those require two conditions to compare.
  - L2L1_conn uses |r| sums (not delta-r sums).
  - Adds positive_fraction (fraction of L1 edges with r > 0).
  - local_clustering_coefficient = L1_clique_density (triangles through focus gene).

Pipeline position
-----------------
Stage 2  02_build_network.py          →  network_edges.h5
Stage 3  THIS SCRIPT (per gene)       →  per_gene/{gi}_{gene_id}.h5
Stage 3b 03b_collect_metrics.py       →  network_metrics_summary.h5 + metrics.tsv

Output HDF5 layout (per_gene/{gi}_{gene_id}.h5)
-------------------------------------------------
    meta/
        gene_index, gene_name, n_genes, n_edges
    metrics/     ← HDF5 attributes, collected by Stage 3b
        (all METRIC_KEYS listed below)
    partners/
        direct   (L1_n_nodes,) int32
        indirect (L2_n_nodes,) int32

Metric keys (METRIC_KEYS)
--------------------------
    degree, weighted_degree
    mean_abs_corr, max_abs_corr, std_abs_corr
    positive_edge_count, negative_edge_count, positive_fraction
    L1_n_nodes, L1_conn_sum, L1_conn_mean
    L1_clique_density, n_l1_to_l1_edges
    local_clustering_coefficient
    L2_n_nodes, L2_n_edges
    L2_conn_sum, L2_conn_mean, n_l1_to_l2_edges
    L2L1_deg, L2L1_conn
    full_n_edges, full_conn_sum
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


# =============================================================================
# Exported constants (used by Stage 3b)
# =============================================================================

METRIC_KEYS: list[str] = [
    "degree", "weighted_degree",
    "mean_abs_corr", "max_abs_corr", "std_abs_corr",
    "positive_edge_count", "negative_edge_count", "positive_fraction",
    "L1_n_nodes", "L1_conn_sum", "L1_conn_mean",
    "L1_clique_density", "n_l1_to_l1_edges",
    "local_clustering_coefficient",
    "L2_n_nodes", "L2_n_edges",
    "L2_conn_sum", "L2_conn_mean", "n_l1_to_l2_edges",
    "L2L1_deg", "L2L1_conn",
    "full_n_edges", "full_conn_sum",
]

_INT_KEYS: frozenset = frozenset([
    "degree", "positive_edge_count", "negative_edge_count",
    "L1_n_nodes", "n_l1_to_l1_edges",
    "L2_n_nodes", "L2_n_edges", "n_l1_to_l2_edges",
    "full_n_edges",
])


# =============================================================================
# CSR neighbour loader
# =============================================================================

def _load_neighbors(
    indptr: np.ndarray,
    h5: h5py.File,
    gene: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Slice (neighbor_indices, corr_values) for one gene from the open HDF5.
    Reads only the relevant HDF5 slice — no full array load.
    """
    s, e = int(indptr[gene]), int(indptr[gene + 1])
    if s == e:
        return np.array([], dtype=np.int32), np.array([], dtype=np.float32)
    return (h5["adjacency/indices"][s:e].astype(np.int32),
            h5["adjacency/data"][s:e].astype(np.float32))


# =============================================================================
# L1/L2 analysis
# =============================================================================

def analyze_gene(focus: int, indptr: np.ndarray, h5_net: h5py.File) -> dict:
    """
    Compute all per-gene metrics for one focus gene.

    Parameters
    ----------
    focus   : gene index (0-based)
    indptr  : (n_genes+1,) int64 — preloaded full array (small, ~144 KB for 18K genes)
    h5_net  : open h5py.File for network_edges.h5

    Returns
    -------
    Flat dict of all metric values plus private _direct_partners / _indirect_partners.
    """
    # --- L1: direct partners ---
    l1_idx, l1_corr = _load_neighbors(indptr, h5_net, focus)
    L1_n = len(l1_idx)

    if L1_n == 0:
        return _empty_metrics()

    l1_set = set(l1_idx.tolist())
    abs_l1 = np.abs(l1_corr)

    # Basic connectivity
    weighted_degree  = float(abs_l1.sum())
    mean_abs_corr    = float(abs_l1.mean())
    max_abs_corr     = float(abs_l1.max())
    std_abs_corr     = float(abs_l1.std()) if L1_n > 1 else 0.0
    pos_count        = int((l1_corr > 0).sum())
    neg_count        = int((l1_corr < 0).sum())
    positive_fraction = float(pos_count / L1_n)

    # --- Load neighbour lists for all L1 partners (needed for L1-L1 and L1-L2) ---
    l1_nb:   dict[int, list[int]]   = {}
    l1_corr_abs: dict[int, dict[int, float]] = {}

    for partner in l1_idx.tolist():
        nb_idx, nb_corr = _load_neighbors(indptr, h5_net, partner)
        l1_nb[partner] = nb_idx.tolist()
        l1_corr_abs[partner] = {
            int(nb): float(abs(rv))
            for nb, rv in zip(nb_idx.tolist(), nb_corr.tolist())
        }

    # --- L1-L1 edges (both endpoints in L1, neither is the focus gene) ---
    l1l1_pairs: dict[tuple, float] = {}
    l1_list = l1_idx.tolist()
    for a_idx, pa in enumerate(l1_list):
        for pb in l1_list[a_idx + 1:]:
            if pb in l1_nb.get(pa, []):
                key = (min(pa, pb), max(pa, pb))
                if key not in l1l1_pairs:
                    l1l1_pairs[key] = l1_corr_abs[pa].get(pb, 0.0)

    n_l1l1 = len(l1l1_pairs)
    max_l1_pairs = L1_n * (L1_n - 1) / 2
    L1_clique_density = float(n_l1l1 / max_l1_pairs) if max_l1_pairs > 0 else 0.0

    # --- L2 nodes: neighbours of L1 that are not focus and not L1 ---
    l2_set: set[int] = set()
    for partner in l1_list:
        for nb in l1_nb.get(partner, []):
            if nb != focus and nb not in l1_set:
                l2_set.add(nb)
    L2_n = len(l2_set)

    # --- L1-L2 edges ---
    l1l2_pairs: dict[tuple, float] = {}
    for partner in l1_list:
        for nb, rv in l1_corr_abs[partner].items():
            if nb in l2_set:
                key = (min(partner, nb), max(partner, nb))
                if key not in l1l2_pairs:
                    l1l2_pairs[key] = rv

    n_l1l2 = len(l1l2_pairs)

    # --- L2 connectivity (outer ring = L1-L1 + L1-L2) ---
    L2_n_edges = n_l1l1 + n_l1l2
    l2_corr_all = np.array(
        list(l1l1_pairs.values()) + list(l1l2_pairs.values()),
        dtype=np.float32,
    )
    L2_conn_sum  = float(l2_corr_all.sum())  if len(l2_corr_all) > 0 else 0.0
    L2_conn_mean = float(l2_corr_all.mean()) if len(l2_corr_all) > 0 else 0.0

    # --- Topological ratios ---
    L2L1_deg  = float(L2_n / L1_n)
    L2L1_conn = float(L2_conn_sum / weighted_degree) if weighted_degree > 0 else 0.0

    # --- Full 2-layer neighbourhood ---
    full_n_edges = L1_n + L2_n_edges
    full_conn_sum = weighted_degree + L2_conn_sum

    return {
        "degree": L1_n,
        "weighted_degree": weighted_degree,
        "mean_abs_corr": mean_abs_corr,
        "max_abs_corr": max_abs_corr,
        "std_abs_corr": std_abs_corr,
        "positive_edge_count": pos_count,
        "negative_edge_count": neg_count,
        "positive_fraction": positive_fraction,
        "L1_n_nodes": L1_n,
        "L1_conn_sum": weighted_degree,
        "L1_conn_mean": mean_abs_corr,
        "L1_clique_density": L1_clique_density,
        "n_l1_to_l1_edges": n_l1l1,
        "local_clustering_coefficient": L1_clique_density,
        "L2_n_nodes": L2_n,
        "L2_n_edges": L2_n_edges,
        "L2_conn_sum": L2_conn_sum,
        "L2_conn_mean": L2_conn_mean,
        "n_l1_to_l2_edges": n_l1l2,
        "L2L1_deg": L2L1_deg,
        "L2L1_conn": L2L1_conn,
        "full_n_edges": full_n_edges,
        "full_conn_sum": full_conn_sum,
        "_direct_partners":   np.array(sorted(l1_set), dtype=np.int32),
        "_indirect_partners": np.array(sorted(l2_set), dtype=np.int32),
    }


def _empty_metrics() -> dict:
    return {
        k: (0 if k in _INT_KEYS else 0.0) for k in METRIC_KEYS
    } | {
        "_direct_partners":   np.array([], dtype=np.int32),
        "_indirect_partners": np.array([], dtype=np.int32),
    }


# =============================================================================
# I/O
# =============================================================================

def run_per_gene(
    network_h5_path: Path,
    gene_index: int,
    out_h5_path: Path,
) -> None:
    with h5py.File(network_h5_path, "r") as h5:
        n_genes = int(h5["meta"].attrs["n_genes"])
        n_edges = int(h5["meta"].attrs["n_edges"])
        # indptr is small (~144 KB for 18K genes) — load fully
        indptr = h5["adjacency/indptr"][:]
        gene_names = (
            [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
            if "gene_names" in h5 else None
        )
        gene_name = gene_names[gene_index] if gene_names else f"gene_{gene_index}"
        metrics = analyze_gene(gene_index, indptr, h5)

    out_h5_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_h5_path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["gene_index"] = gene_index
        meta.attrs["gene_name"]  = gene_name
        meta.attrs["n_genes"]    = n_genes
        meta.attrs["n_edges"]    = n_edges

        met = h5.create_group("metrics")
        for key in METRIC_KEYS:
            met.attrs[key] = metrics[key]

        grp_p = h5.create_group("partners")
        grp_p.create_dataset("direct",   data=metrics["_direct_partners"],
                             compression="gzip", compression_opts=6)
        grp_p.create_dataset("indirect", data=metrics["_indirect_partners"],
                             compression="gzip", compression_opts=6)


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="[Network Metrics Stage 3] Per-gene L1/L2 neighbourhood metrics."
    )
    p.add_argument("--network-h5",  type=str, required=True,
                   help="Path to network_edges.h5 from Stage 2.")
    p.add_argument("--gene-index",  type=int, required=True,
                   help="Zero-based gene index.")
    p.add_argument("--out-h5",      type=str, required=True,
                   help="Output per-gene HDF5 path.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_per_gene(
        network_h5_path=Path(args.network_h5),
        gene_index=args.gene_index,
        out_h5_path=Path(args.out_h5),
    )


if __name__ == "__main__":
    main()
