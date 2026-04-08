#!/usr/bin/env python3
"""
Stage 2: Build network adjacency index and compute global topology metrics.

Reads significant edges from corr_significant.h5 (Stage 1) and produces
network_edges.h5 which contains:
  1. Global topology metrics (degree distribution, scale-free exponent,
     clustering coefficient, connected components, assortativity, etc.)
  2. CSR-format adjacency (indptr / indices / data) enabling per-gene Stage 3
     jobs to load only O(degree) bytes rather than the full edge list.

Pipeline position
-----------------
Stage 1  01_compute_correlations.py   →  corr_significant.h5
Stage 2  THIS SCRIPT                  →  network_edges.h5
Stage 3  03_per_gene_metrics.py       →  per_gene/{gi}_{gene_id}.h5

CSR adjacency rationale
------------------------
For 18K genes and ~5M significant edges, the full edge list is ~60 MB.
Each Stage 3 per-gene job needs only its gene's neighbours.
CSR layout stores both directions (gene_i→gene_j and gene_j→gene_i), sorted
within each row, so a single HDF5 slice read gives O(degree) values per gene.

Output HDF5 layout
------------------
    meta/
        n_genes, n_edges, n_genes_active, density, fdr_alpha, corr_threshold
    edges/
        gene_i, gene_j, corr, qval  (n_edges,) — copied from Stage 1 output
    gene_names  (n_genes,) str dataset
    adjacency/
        indptr   (n_genes+1,) int64  — indptr[g]:indptr[g+1] gives slice for gene g
        indices  (2*n_edges,) int32  — neighbour indices, sorted within each row
        data     (2*n_edges,) float32 — correlation values aligned to indices
    global_topology/   (HDF5 attributes + degrees dataset)
        n_genes, n_genes_active, fraction_active
        n_edges, density
        avg_degree, median_degree, max_degree, min_degree, std_degree
        degree_p25/50/75/90/95/99
        power_law_exponent, power_law_r_squared, is_scale_free
        global_clustering
        n_components, largest_component_size, largest_component_fraction
        assortativity
        mean_abs_corr, median_abs_corr, std_abs_corr
        positive_edge_fraction, negative_edge_fraction
        degrees  (n_genes,) int32  — per-gene degree array
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from pathlib import Path

import h5py
import numpy as np
from scipy import stats as scipy_stats


# =============================================================================
# Topology helpers
# =============================================================================

def _scale_free_metrics(degrees: np.ndarray) -> dict:
    """Power-law fit on degree distribution. γ = 2-3 → scale-free network."""
    active = degrees[degrees > 0]
    nan_result = {"power_law_exponent": float("nan"),
                  "power_law_r_squared": float("nan"),
                  "is_scale_free": False}
    if len(active) < 10:
        return nan_result
    unique, counts = np.unique(active, return_counts=True)
    mask = counts >= 2
    if mask.sum() < 5:
        return nan_result
    slope, _, r, _, _ = scipy_stats.linregress(
        np.log10(unique[mask].astype(float)),
        np.log10(counts[mask].astype(float))
    )
    gamma = -slope
    return {
        "power_law_exponent": float(gamma),
        "power_law_r_squared": float(r ** 2),
        "is_scale_free": bool((2.0 <= gamma <= 3.5) and (r ** 2 > 0.8)),
    }


def _clustering_coefficient(adj: dict[int, set]) -> float:
    """Global clustering = mean of local clustering over nodes with degree ≥ 2."""
    vals = []
    for node, neighbors in adj.items():
        k = len(neighbors)
        if k < 2:
            continue
        nb_list = list(neighbors)
        tri = sum(
            1 for a, n1 in enumerate(nb_list)
            for n2 in nb_list[a + 1:]
            if n2 in adj.get(n1, set())
        )
        vals.append(tri / (k * (k - 1) / 2))
    return float(np.mean(vals)) if vals else 0.0


def _connected_components(adj: dict[int, set], n_genes: int) -> dict:
    """Number of connected subgraphs + largest component size via BFS."""
    visited: set[int] = set()
    sizes: list[int] = []
    for start in range(n_genes):
        if start in visited or not adj.get(start):
            continue
        cc: set[int] = set()
        q = deque([start])
        while q:
            node = q.popleft()
            if node in visited:
                continue
            visited.add(node)
            cc.add(node)
            q.extend(nb for nb in adj.get(node, set()) if nb not in visited)
        sizes.append(len(cc))
    sizes = sorted(sizes, reverse=True)
    return {
        "n_components": len(sizes),
        "largest_component_size": sizes[0] if sizes else 0,
        "largest_component_fraction": float(sizes[0] / n_genes) if sizes else 0.0,
    }


def compute_global_topology(
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    corr: np.ndarray,
    n_genes: int,
) -> dict:
    """Compute all global topology metrics for the co-expression network."""
    n_edges = len(gene_i)
    _zero = {
        "n_genes": n_genes, "n_genes_active": 0, "fraction_active": 0.0,
        "n_edges": 0, "density": 0.0,
        "avg_degree": 0.0, "median_degree": 0.0, "max_degree": 0,
        "min_degree": 0, "std_degree": 0.0,
        **{f"degree_p{p}": 0.0 for p in [25, 50, 75, 90, 95, 99]},
        "power_law_exponent": float("nan"), "power_law_r_squared": float("nan"),
        "is_scale_free": False,
        "global_clustering": 0.0,
        "n_components": 0, "largest_component_size": 0,
        "largest_component_fraction": 0.0,
        "assortativity": 0.0,
        "mean_abs_corr": 0.0, "median_abs_corr": 0.0, "std_abs_corr": 0.0,
        "positive_edge_fraction": 0.0, "negative_edge_fraction": 0.0,
        "degrees": np.zeros(n_genes, dtype=np.int32),
    }
    if n_edges == 0:
        return _zero

    # Degree distribution (vectorised)
    degrees = np.zeros(n_genes, dtype=np.int32)
    np.add.at(degrees, gene_i, 1)
    np.add.at(degrees, gene_j, 1)
    active_degrees = degrees[degrees > 0]
    n_active = len(active_degrees)

    # Build adjacency for clustering + components
    adj: dict[int, set] = defaultdict(set)
    for ii, jj in zip(gene_i.tolist(), gene_j.tolist()):
        adj[ii].add(jj)
        adj[jj].add(ii)

    density = 2 * n_edges / (n_genes * (n_genes - 1)) if n_genes > 1 else 0.0
    percentiles = (
        {f"degree_p{p}": float(np.percentile(active_degrees, p)) for p in [25, 50, 75, 90, 95, 99]}
        if n_active > 0 else {f"degree_p{p}": 0.0 for p in [25, 50, 75, 90, 95, 99]}
    )

    sf = _scale_free_metrics(degrees)
    cc = _connected_components(adj, n_genes)

    deg_i = degrees[gene_i].astype(float)
    deg_j = degrees[gene_j].astype(float)
    assortativity = 0.0
    if len(deg_i) >= 2 and deg_i.std() > 0 and deg_j.std() > 0:
        assortativity, _ = scipy_stats.pearsonr(deg_i, deg_j)

    abs_c = np.abs(corr)
    return {
        "n_genes": n_genes,
        "n_genes_active": n_active,
        "fraction_active": float(n_active / n_genes) if n_genes > 0 else 0.0,
        "n_edges": n_edges,
        "density": float(density),
        "avg_degree": float(2 * n_edges / n_genes) if n_genes > 0 else 0.0,
        "median_degree": float(np.median(active_degrees)) if n_active > 0 else 0.0,
        "max_degree": int(degrees.max()),
        "min_degree": int(active_degrees.min()) if n_active > 0 else 0,
        "std_degree": float(active_degrees.std()) if n_active > 0 else 0.0,
        **percentiles, **sf, **cc,
        "global_clustering": float(_clustering_coefficient(adj)),
        "assortativity": float(assortativity),
        "mean_abs_corr": float(abs_c.mean()),
        "median_abs_corr": float(np.median(abs_c)),
        "std_abs_corr": float(abs_c.std()),
        "positive_edge_fraction": float((corr > 0).sum() / n_edges),
        "negative_edge_fraction": float((corr < 0).sum() / n_edges),
        "degrees": degrees,
    }


def build_csr_adjacency(
    gene_i: np.ndarray,
    gene_j: np.ndarray,
    corr: np.ndarray,
    n_genes: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build CSR-format undirected adjacency from edge list.

    Both directions are stored so each gene g has all its neighbours in
    indices[indptr[g]:indptr[g+1]] — sorted for efficient lookup.

    Returns
    -------
    indptr  : (n_genes+1,) int64
    indices : (2*n_edges,) int32
    data    : (2*n_edges,) float32
    """
    n_edges = len(gene_i)
    # Duplicate edges for both directions
    src = np.concatenate([gene_i, gene_j]).astype(np.int32)
    dst = np.concatenate([gene_j, gene_i]).astype(np.int32)
    dat = np.concatenate([corr, corr]).astype(np.float32)

    # Sort by source gene
    order = np.argsort(src, kind="stable")
    src, dst, dat = src[order], dst[order], dat[order]

    # Build indptr
    indptr = np.zeros(n_genes + 1, dtype=np.int64)
    for g in src:
        indptr[g + 1] += 1
    np.cumsum(indptr, out=indptr)

    # Sort within each row by destination index
    indices = np.empty(2 * n_edges, dtype=np.int32)
    data = np.empty(2 * n_edges, dtype=np.float32)
    for g in range(n_genes):
        s, e = int(indptr[g]), int(indptr[g + 1])
        if s == e:
            continue
        row_ord = np.argsort(dst[s:e], kind="stable")
        indices[s:e] = dst[s:e][row_ord]
        data[s:e] = dat[s:e][row_ord]

    return indptr, indices, data


# =============================================================================
# Main entry
# =============================================================================

def build_network(
    corr_h5_path: Path,
    out_h5_path: Path,
    compression_level: int = 6,
) -> None:
    print(f"Loading significant edges from {corr_h5_path}...")
    with h5py.File(corr_h5_path, "r") as h5:
        n_genes   = int(h5["meta"].attrs["n_genes"])
        n_sig     = int(h5["meta"].attrs["n_significant"])
        fdr_alpha = float(h5["meta"].attrs["fdr_alpha"])
        corr_thr  = float(h5["meta"].attrs["corr_threshold"])
        gene_names = (
            [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
            if "gene_names" in h5 else None
        )
        if n_sig > 0:
            gene_i = h5["edges/gene_i"][:]
            gene_j = h5["edges/gene_j"][:]
            corr   = h5["edges/corr"][:]
            qval   = h5["edges/qval"][:]
        else:
            gene_i = np.array([], dtype=np.int32)
            gene_j = np.array([], dtype=np.int32)
            corr   = np.array([], dtype=np.float32)
            qval   = np.array([], dtype=np.float32)

    print(f"  {n_genes:,} genes, {n_sig:,} significant edges")

    # --- Global topology ---
    print("Computing global topology...")
    topo = compute_global_topology(gene_i, gene_j, corr, n_genes)
    print(f"  Active genes     : {topo['n_genes_active']:,} ({100*topo['fraction_active']:.1f}%)")
    print(f"  Density          : {topo['density']:.2e}")
    print(f"  Avg / max degree : {topo['avg_degree']:.1f} / {topo['max_degree']}")
    print(f"  Global clustering: {topo['global_clustering']:.4f}")
    print(f"  Components       : {topo['n_components']}, "
          f"largest {topo['largest_component_fraction']:.3f} of network")
    if not np.isnan(topo["power_law_exponent"]):
        print(f"  Scale-free       : γ={topo['power_law_exponent']:.2f}, "
              f"R²={topo['power_law_r_squared']:.3f}, "
              f"is_scale_free={topo['is_scale_free']}")
    print(f"  Assortativity    : {topo['assortativity']:.4f}")
    print(f"  Positive edges   : {topo['positive_edge_fraction']:.1%}")

    # --- CSR adjacency ---
    if n_sig > 0:
        print("Building CSR adjacency...")
        indptr, indices, csr_data = build_csr_adjacency(gene_i, gene_j, corr, n_genes)
        print(f"  CSR sizes: indptr {indptr.nbytes/1e6:.1f} MB, "
              f"indices {indices.nbytes/1e6:.1f} MB, data {csr_data.nbytes/1e6:.1f} MB")
    else:
        indptr   = np.zeros(n_genes + 1, dtype=np.int64)
        indices  = np.array([], dtype=np.int32)
        csr_data = np.array([], dtype=np.float32)

    # --- Write output ---
    print(f"\nSaving to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)
    kw = dict(compression="gzip", compression_opts=compression_level)

    with h5py.File(out_h5_path, "w") as h5:
        # Meta
        meta = h5.create_group("meta")
        meta.attrs["n_genes"]        = n_genes
        meta.attrs["n_edges"]        = n_sig
        meta.attrs["n_genes_active"] = topo["n_genes_active"]
        meta.attrs["density"]        = topo["density"]
        meta.attrs["fdr_alpha"]      = fdr_alpha
        meta.attrs["corr_threshold"] = corr_thr

        # Gene names
        if gene_names is not None:
            h5.create_dataset("gene_names", data=gene_names, dtype=h5py.string_dtype())

        # Raw edges (for downstream use)
        grp_e = h5.create_group("edges")
        chunk_e = max(1, min(1_000_000, n_sig))
        if n_sig > 0:
            grp_e.create_dataset("gene_i", data=gene_i, chunks=(chunk_e,), **kw)
            grp_e.create_dataset("gene_j", data=gene_j, chunks=(chunk_e,), **kw)
            grp_e.create_dataset("corr",   data=corr,   chunks=(chunk_e,), **kw)
            grp_e.create_dataset("qval",   data=qval,   chunks=(chunk_e,), **kw)
        else:
            for name, dt in [("gene_i", np.int32), ("gene_j", np.int32),
                              ("corr", np.float32), ("qval", np.float32)]:
                grp_e.create_dataset(name, data=np.array([], dtype=dt), **kw)

        # CSR adjacency (key for Stage 3 efficiency)
        grp_adj = h5.create_group("adjacency")
        grp_adj.create_dataset("indptr", data=indptr, compression="gzip", compression_opts=4)
        if len(indices) > 0:
            chunk_a = max(1, min(1_000_000, len(indices)))
            grp_adj.create_dataset("indices", data=indices,  chunks=(chunk_a,), **kw)
            grp_adj.create_dataset("data",    data=csr_data, chunks=(chunk_a,), **kw)
        else:
            grp_adj.create_dataset("indices", data=np.array([], dtype=np.int32),   **kw)
            grp_adj.create_dataset("data",    data=np.array([], dtype=np.float32), **kw)

        # Global topology scalars as HDF5 attributes
        grp_t = h5.create_group("global_topology")
        for key, val in topo.items():
            if key == "degrees":
                grp_t.create_dataset("degrees", data=val,
                                     compression="gzip", compression_opts=compression_level)
            elif isinstance(val, (bool, np.bool_)):
                grp_t.attrs[key] = int(val)
            elif isinstance(val, (float, int, np.integer, np.floating)):
                grp_t.attrs[key] = val

    print("Done.")


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="[Network Metrics Stage 2] Build network + global topology + CSR adjacency."
    )
    p.add_argument("--corr-h5", type=str, required=True,
                   help="Path to corr_significant.h5 from Stage 1.")
    p.add_argument("--out-h5", type=str, required=True,
                   help="Output path for network_edges.h5.")
    p.add_argument("--compression-level", type=int, default=6,
                   help="Gzip compression level 1-9 (default: 6).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    build_network(
        corr_h5_path=Path(args.corr_h5),
        out_h5_path=Path(args.out_h5),
        compression_level=args.compression_level,
    )


if __name__ == "__main__":
    main()
