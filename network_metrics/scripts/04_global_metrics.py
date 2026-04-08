#!/usr/bin/env python3
"""
Stage 4 (optional): Eigenvector centrality, PageRank, and approximate betweenness.

This stage is disabled by default (skip_global_metrics: true in config).
Enable it for a more complete per-gene centrality profile.

Two-tier betweenness strategy
------------------------------
Tier 1 (always computed): eigenvector_centrality via scipy.sparse power iteration.
  - O(k × E) where k ≈ 50-100 iterations, E = n_edges.
  - For 5M edges: ~500M operations ≈ minutes.
  - Best for scale-free networks: captures hub-of-hubs structure.

Tier 2 (requires --compute-betweenness): approximate Brandes via `networkit`.
  - Requires `networkit` package (C++ extension, pip install networkit).
  - EstimateBetweenness(G, nSamples=500) → error O(1/√nSamples).
  - For 18K nodes, 5M edges, nSamples=500: ~10-30 min.
  - Falls back to None if networkit not installed.

Why not exact betweenness?
  Exact Brandes is O(V × E). For 18K genes and 5M edges: 18K × 5M = 9×10¹⁰ ops.
  Even at 10⁹ ops/sec, this is ~90,000 seconds. Not feasible in Python.

Pipeline position
-----------------
Stage 3b  03b_collect_metrics.py     →  network_metrics_summary.h5
Stage 4   THIS SCRIPT (optional)     →  global_topology.h5

Output HDF5 layout (global_topology.h5)
-----------------------------------------
    meta/    n_genes, n_edges, computed_betweenness
    per_gene/
        eigenvector_centrality  (n_genes,) float32
        pagerank                (n_genes,) float32
        approx_betweenness      (n_genes,) float32    [only if --compute-betweenness]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp


# =============================================================================
# Eigenvector centrality
# =============================================================================

def compute_eigenvector_centrality(
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
    n_genes: int,
    max_iter: int = 100,
    tol: float = 1e-6,
) -> np.ndarray:
    """
    Power iteration for eigenvector centrality on the weighted adjacency.

    Converges in O(max_iter × n_edges). Each iteration is a sparse matrix-vector
    multiply (scipy CSR), which is cache-friendly and well-optimised.

    Parameters
    ----------
    indptr, indices, data : CSR adjacency from network_edges.h5
    n_genes               : number of genes
    max_iter              : maximum power iterations (default: 100)
    tol                   : convergence tolerance

    Returns
    -------
    centrality : (n_genes,) float32, normalised to unit L2 norm
    """
    if len(indices) == 0:
        return np.zeros(n_genes, dtype=np.float32)

    A = sp.csr_matrix(
        (data.astype(np.float64), indices.astype(np.int32), indptr.astype(np.int32)),
        shape=(n_genes, n_genes),
    )

    # Use absolute values for centrality (negative correlations count too)
    A = A.copy()
    A.data = np.abs(A.data)

    x = np.ones(n_genes, dtype=np.float64)
    x /= np.linalg.norm(x)

    for iteration in range(max_iter):
        x_new = A @ x
        norm = np.linalg.norm(x_new)
        if norm == 0:
            break
        x_new /= norm
        if np.linalg.norm(x_new - x) < tol:
            print(f"  Eigenvector centrality converged after {iteration + 1} iterations.")
            break
        x = x_new
    else:
        print(f"  Warning: eigenvector centrality did not converge in {max_iter} iterations.")

    return x_new.astype(np.float32)


# =============================================================================
# PageRank
# =============================================================================

def compute_pagerank(
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
    n_genes: int,
    damping: float = 0.85,
    max_iter: int = 100,
    tol: float = 1e-6,
) -> np.ndarray:
    """
    Power-iteration PageRank on the (undirected) co-expression network.

    Uses |r| as edge weights. For undirected networks, each gene is both
    a link source and a link target, so PageRank measures influence equally
    in both directions.
    """
    if len(indices) == 0:
        return np.ones(n_genes, dtype=np.float32) / n_genes

    # Build row-normalised transition matrix
    A = sp.csr_matrix(
        (np.abs(data).astype(np.float64), indices.astype(np.int32), indptr.astype(np.int32)),
        shape=(n_genes, n_genes),
    )
    # Row-normalise: divide each row by its sum (so rows sum to 1)
    row_sums = np.asarray(A.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0  # avoid division by zero for isolated nodes
    D_inv = sp.diags(1.0 / row_sums)
    M = D_inv @ A  # row-stochastic matrix

    pr = np.ones(n_genes, dtype=np.float64) / n_genes
    teleport = (1 - damping) / n_genes

    for iteration in range(max_iter):
        pr_new = damping * (M.T @ pr) + teleport
        if np.linalg.norm(pr_new - pr, 1) < tol:
            print(f"  PageRank converged after {iteration + 1} iterations.")
            break
        pr = pr_new
    else:
        print(f"  Warning: PageRank did not converge in {max_iter} iterations.")

    return pr_new.astype(np.float32)


# =============================================================================
# Approximate betweenness (networkit)
# =============================================================================

def compute_approx_betweenness(
    indptr: np.ndarray,
    indices: np.ndarray,
    data: np.ndarray,
    n_genes: int,
    n_samples: int = 500,
) -> np.ndarray | None:
    """
    Approximate betweenness centrality via Brandes sampling (networkit).

    Requires `networkit` (C++ extension). Install with:  pip install networkit

    Parameters
    ----------
    n_samples : number of pivot nodes for Brandes sampling (default: 500)
                Error scales as O(1/√n_samples). Use 1000+ for publication quality.

    Returns
    -------
    betweenness : (n_genes,) float32 normalised scores, or None if networkit unavailable
    """
    try:
        import networkit as nk
    except ImportError:
        print("  networkit not installed — skipping approximate betweenness.")
        print("  Install with: pip install networkit")
        return None

    if len(indices) == 0:
        return np.zeros(n_genes, dtype=np.float32)

    print(f"  Building networkit graph ({n_genes} nodes, {len(indices)//2} edges)...")
    G = nk.Graph(n_genes, weighted=True, directed=False)

    # Add edges (use indptr to iterate upper triangle only)
    for g in range(n_genes):
        s, e = int(indptr[g]), int(indptr[g + 1])
        for idx in range(s, e):
            nb = int(indices[idx])
            if nb > g:  # upper triangle only to avoid duplicates
                G.addEdge(g, nb, abs(float(data[idx])))

    print(f"  Running EstimateBetweenness (nSamples={n_samples})...")
    est = nk.centrality.EstimateBetweenness(G, nSamples=n_samples, normalized=True)
    est.run()
    scores = np.array(est.scores(), dtype=np.float32)
    print("  Approximate betweenness done.")
    return scores


# =============================================================================
# Main entry
# =============================================================================

def compute_global_metrics(
    network_h5_path: Path,
    out_h5_path: Path,
    compute_betweenness: bool = False,
    betweenness_n_samples: int = 500,
    eig_max_iter: int = 100,
    eig_tol: float = 1e-6,
) -> None:
    print(f"Loading CSR adjacency from {network_h5_path}...")
    with h5py.File(network_h5_path, "r") as h5:
        n_genes = int(h5["meta"].attrs["n_genes"])
        n_edges = int(h5["meta"].attrs["n_edges"])
        indptr  = h5["adjacency/indptr"][:]
        indices = h5["adjacency/indices"][:]
        data    = h5["adjacency/data"][:]

    print(f"  {n_genes:,} genes, {n_edges:,} edges  (CSR: {indptr.nbytes + indices.nbytes + data.nbytes:.0f} B)")

    # --- Eigenvector centrality ---
    print("\nComputing eigenvector centrality (power iteration)...")
    eig = compute_eigenvector_centrality(indptr, indices, data, n_genes,
                                         max_iter=eig_max_iter, tol=eig_tol)

    # --- PageRank ---
    print("Computing PageRank (power iteration)...")
    pr = compute_pagerank(indptr, indices, data, n_genes)

    # --- Approximate betweenness (optional) ---
    betw = None
    if compute_betweenness:
        print(f"\nComputing approximate betweenness (n_samples={betweenness_n_samples})...")
        betw = compute_approx_betweenness(indptr, indices, data, n_genes,
                                           n_samples=betweenness_n_samples)

    # --- Write output ---
    print(f"\nSaving to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)
    kw = dict(compression="gzip", compression_opts=6)

    with h5py.File(out_h5_path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["n_genes"]              = n_genes
        meta.attrs["n_edges"]              = n_edges
        meta.attrs["computed_betweenness"] = int(betw is not None)
        if betw is not None:
            meta.attrs["betweenness_n_samples"] = betweenness_n_samples

        pg = h5.create_group("per_gene")
        pg.create_dataset("eigenvector_centrality", data=eig,  **kw)
        pg.create_dataset("pagerank",               data=pr,   **kw)
        if betw is not None:
            pg.create_dataset("approx_betweenness", data=betw, **kw)

    print("Done.")
    print(f"  Top 5 by eigenvector centrality: {np.argsort(eig)[::-1][:5].tolist()}")
    print(f"  Top 5 by pagerank: {np.argsort(pr)[::-1][:5].tolist()}")


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="[Network Metrics Stage 4] Eigenvector centrality + optional betweenness."
    )
    p.add_argument("--network-h5",    type=str, required=True,
                   help="Path to network_edges.h5 from Stage 2.")
    p.add_argument("--out-h5",        type=str, required=True,
                   help="Output path for global_topology.h5.")
    p.add_argument("--compute-betweenness", action="store_true",
                   help="Compute approximate betweenness (requires networkit).")
    p.add_argument("--betweenness-n-samples", type=int, default=500,
                   help="Pivot nodes for Brandes sampling (default: 500).")
    p.add_argument("--eig-max-iter",  type=int, default=100,
                   help="Max power iterations for eigenvector centrality (default: 100).")
    p.add_argument("--eig-tol",       type=float, default=1e-6,
                   help="Convergence tolerance for eigenvector centrality (default: 1e-6).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    compute_global_metrics(
        network_h5_path=Path(args.network_h5),
        out_h5_path=Path(args.out_h5),
        compute_betweenness=args.compute_betweenness,
        betweenness_n_samples=args.betweenness_n_samples,
        eig_max_iter=args.eig_max_iter,
        eig_tol=args.eig_tol,
    )


if __name__ == "__main__":
    main()
