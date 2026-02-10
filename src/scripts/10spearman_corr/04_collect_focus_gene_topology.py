#!/usr/bin/env python3
"""
Stage 4: Collect Focus Gene Topology Across All Gene-Wise HDF5 Files.

Aggregates topology information for focus genes across all per-gene networks,
where each network is defined by a different reference gene's low/high split.

Pipeline position
-----------------
Stage 1   01subset/01get_extreme_pop_bootstrap.py  →  bootstrap_indices.h5
Stage 2a  02a_calc_base_correlations.py            →  base_correlations.h5
Stage 2b  02b_bootstrap_significant_edges.py       →  bootstrap_significant.h5
Stage 3   03_reconstruct_diff_network.py           →  differential_network_gene_*.h5
Stage 4   THIS SCRIPT                              →  focus_gene_topology.h5

Output Structure
----------------
For each focus gene × reference gene combination:
- degree_low[focus, ref]:  Focus gene's degree in reference gene's low network
- degree_high[focus, ref]: Focus gene's degree in reference gene's high network
- degree_diff[focus, ref]: Focus gene's rewiring in reference gene's context

Usage
-----
python 04_collect_focus_gene_topology.py \\
    --network-dir results/networks \\
    --focus-genes 0,1,2,3,4  OR  top:50 \\
    --n-genes 20000 \\
    --out-h5 results/focus_gene_topology.h5 \\
    --n-jobs 8
"""

from __future__ import annotations

import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import re

import h5py
import numpy as np


def extract_focus_gene_stats_from_file(
    h5_path: Path,
    focus_genes: list[int],
) -> dict:
    """
    Extract topology stats for focus genes from a single differential network file.

    Parameters
    ----------
    h5_path : Path to differential_network_gene_*.h5
    focus_genes : List of focus gene indices to extract

    Returns
    -------
    dict with:
        ref_gene: int - reference gene index (from filename or metadata)
        focus_stats: dict[focus_gene] -> {degree_low, degree_high, degree_diff, ...}
    """
    # Extract reference gene index from filename (e.g., differential_network_gene_123.h5)
    match = re.search(r"gene[_-]?(\d+)", h5_path.stem)
    if match:
        ref_gene = int(match.group(1))
    else:
        ref_gene = 0

    focus_stats = {}

    try:
        with h5py.File(h5_path, "r") as h5:
            # Check if file has topology data
            if "topology" not in h5:
                return {"ref_gene": ref_gene, "focus_stats": focus_stats}

            # Get degree arrays from global topology
            topo = h5["topology"]

            # Degree arrays for each network type
            degrees_low = topo["global_low/degrees"][:] if "global_low/degrees" in topo else None
            degrees_high = topo["global_high/degrees"][:] if "global_high/degrees" in topo else None
            degrees_diff = topo["global_diff/degrees"][:] if "global_diff/degrees" in topo else None

            # Also try to get focus_gene stats if available
            has_focus = "focus_gene" in h5
            if has_focus:
                focus_grp = h5["focus_gene"]
                file_focus_gene = focus_grp.attrs.get("gene_index", -1)
            else:
                file_focus_gene = -1

            # Extract stats for each focus gene
            for fg in focus_genes:
                stats = {
                    "degree_low": int(degrees_low[fg]) if degrees_low is not None else 0,
                    "degree_high": int(degrees_high[fg]) if degrees_high is not None else 0,
                    "degree_diff": int(degrees_diff[fg]) if degrees_diff is not None else 0,
                }

                # If this focus gene matches the file's focus gene, add detailed stats
                if fg == file_focus_gene and has_focus:
                    ds = focus_grp["direct_stats"]
                    stats.update({
                        "n_direct_partners": focus_grp.attrs.get("n_direct_partners", 0),
                        "n_indirect_partners": focus_grp.attrs.get("n_indirect_partners", 0),
                        "direct_n_disappear": ds.attrs.get("n_disappear", 0),
                        "direct_n_new": ds.attrs.get("n_new", 0),
                        "direct_n_sign_change": ds.attrs.get("n_sign_change", 0),
                        "direct_mean_delta": ds.attrs.get("mean_delta", 0.0),
                        "direct_mean_abs_delta": ds.attrs.get("mean_abs_delta", 0.0),
                    })

                focus_stats[fg] = stats

    except Exception as e:
        print(f"  Warning: Failed to read {h5_path}: {e}")
        return {"ref_gene": ref_gene, "focus_stats": {}}

    return {"ref_gene": ref_gene, "focus_stats": focus_stats}


def collect_focus_gene_topology(
    network_dir: Path,
    focus_genes: list[int],
    n_genes: int,
    n_jobs: int = 1,
    file_pattern: str = "differential_network_gene_*.h5",
) -> dict:
    """
    Collect topology for focus genes across all per-gene network files.

    Parameters
    ----------
    network_dir : Directory containing differential_network_gene_*.h5 files
    focus_genes : List of focus gene indices
    n_genes : Total number of genes (for reference gene axis)
    n_jobs : Number of parallel jobs
    file_pattern : Glob pattern for network files

    Returns
    -------
    dict with:
        focus_genes: array of focus gene indices
        ref_genes: array of reference gene indices (found in files)
        degree_low: (n_focus, n_ref) matrix
        degree_high: (n_focus, n_ref) matrix
        degree_diff: (n_focus, n_ref) matrix
    """
    # Find all network files
    h5_files = sorted(network_dir.glob(file_pattern))
    print(f"Found {len(h5_files)} network files in {network_dir}")

    if len(h5_files) == 0:
        raise ValueError(f"No files matching {file_pattern} in {network_dir}")

    n_focus = len(focus_genes)
    n_ref = len(h5_files)

    # Initialize result matrices (focus_gene x ref_gene)
    degree_low = np.zeros((n_focus, n_ref), dtype=np.int32)
    degree_high = np.zeros((n_focus, n_ref), dtype=np.int32)
    degree_diff = np.zeros((n_focus, n_ref), dtype=np.int32)
    ref_genes = np.zeros(n_ref, dtype=np.int32)

    # Process files
    if n_jobs > 1:
        print(f"Processing with {n_jobs} parallel workers...")
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = {
                executor.submit(extract_focus_gene_stats_from_file, f, focus_genes): i
                for i, f in enumerate(h5_files)
            }
            for future in as_completed(futures):
                file_idx = futures[future]
                result = future.result()
                ref_genes[file_idx] = result["ref_gene"]
                for fg_idx, fg in enumerate(focus_genes):
                    if fg in result["focus_stats"]:
                        stats = result["focus_stats"][fg]
                        degree_low[fg_idx, file_idx] = stats["degree_low"]
                        degree_high[fg_idx, file_idx] = stats["degree_high"]
                        degree_diff[fg_idx, file_idx] = stats["degree_diff"]
    else:
        print("Processing sequentially...")
        for file_idx, h5_file in enumerate(h5_files):
            if file_idx % 100 == 0:
                print(f"  Processing file {file_idx + 1}/{n_ref}...")
            result = extract_focus_gene_stats_from_file(h5_file, focus_genes)
            ref_genes[file_idx] = result["ref_gene"]
            for fg_idx, fg in enumerate(focus_genes):
                if fg in result["focus_stats"]:
                    stats = result["focus_stats"][fg]
                    degree_low[fg_idx, file_idx] = stats["degree_low"]
                    degree_high[fg_idx, file_idx] = stats["degree_high"]
                    degree_diff[fg_idx, file_idx] = stats["degree_diff"]

    return {
        "focus_genes": np.array(focus_genes, dtype=np.int32),
        "ref_genes": ref_genes,
        "degree_low": degree_low,
        "degree_high": degree_high,
        "degree_diff": degree_diff,
    }


def save_results(results: dict, out_h5_path: Path, gene_ids: list = None) -> None:
    """Save results to HDF5."""
    print(f"\nSaving results to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out_h5_path, "w") as h5:
        # Metadata
        h5.attrs["n_focus_genes"] = len(results["focus_genes"])
        h5.attrs["n_ref_genes"] = len(results["ref_genes"])

        # Gene indices
        h5.create_dataset("focus_genes", data=results["focus_genes"])
        h5.create_dataset("ref_genes", data=results["ref_genes"])

        # Degree matrices
        h5.create_dataset("degree_low", data=results["degree_low"], compression="gzip")
        h5.create_dataset("degree_high", data=results["degree_high"], compression="gzip")
        h5.create_dataset("degree_diff", data=results["degree_diff"], compression="gzip")

        # Gene IDs if provided
        if gene_ids:
            dt = h5py.string_dtype(encoding="utf-8")
            focus_ids = [gene_ids[i] for i in results["focus_genes"]]
            h5.create_dataset("focus_gene_ids", data=focus_ids, dtype=dt)

    print(f"  Saved {len(results['focus_genes'])} focus genes × {len(results['ref_genes'])} reference genes")


def print_summary(results: dict, gene_ids: list = None) -> None:
    """Print summary statistics."""
    print("\n" + "=" * 70)
    print("FOCUS GENE TOPOLOGY COLLECTION SUMMARY")
    print("=" * 70)

    focus_genes = results["focus_genes"]
    degree_low = results["degree_low"]
    degree_high = results["degree_high"]
    degree_diff = results["degree_diff"]

    print(f"\n[OVERVIEW]")
    print(f"  Focus genes: {len(focus_genes)}")
    print(f"  Reference genes (contexts): {len(results['ref_genes'])}")

    print(f"\n[DEGREE STATISTICS ACROSS CONTEXTS]")
    print(f"  {'FocusGene':>12} {'LowDeg':>10} {'HighDeg':>10} {'DiffDeg':>10} {'Variance':>10}")

    for i, fg in enumerate(focus_genes[:20]):  # Show top 20
        gene_id = gene_ids[fg] if gene_ids else f"gene_{fg}"
        low_mean = np.mean(degree_low[i])
        high_mean = np.mean(degree_high[i])
        diff_mean = np.mean(degree_diff[i])
        diff_var = np.var(degree_diff[i])
        print(f"  {gene_id:>12} {low_mean:>10.1f} {high_mean:>10.1f} {diff_mean:>10.1f} {diff_var:>10.1f}")

    # Identify context-dependent genes (high variance)
    diff_var = np.var(degree_diff, axis=1)
    top_var_idx = np.argsort(diff_var)[::-1][:5]

    print(f"\n[TOP CONTEXT-DEPENDENT GENES] (highest variance in degree_diff)")
    for idx in top_var_idx:
        fg = focus_genes[idx]
        gene_id = gene_ids[fg] if gene_ids else f"gene_{fg}"
        print(f"  {gene_id}: variance = {diff_var[idx]:.1f}, "
              f"range = [{degree_diff[idx].min()}, {degree_diff[idx].max()}]")


def parse_focus_genes(spec: str, n_genes: int, rewiring_scores: np.ndarray = None) -> list[int]:
    """
    Parse focus gene specification.

    Formats:
    - "0,1,2,3,4" - explicit list
    - "top:50" - top N by rewiring score
    - "range:0:100" - range of indices
    """
    if spec.startswith("top:"):
        n = int(spec[4:])
        if rewiring_scores is not None:
            return list(np.argsort(rewiring_scores)[::-1][:n])
        else:
            # Default to first N genes
            return list(range(min(n, n_genes)))
    elif spec.startswith("range:"):
        parts = spec[6:].split(":")
        start, end = int(parts[0]), int(parts[1])
        return list(range(start, min(end, n_genes)))
    else:
        return [int(x) for x in spec.split(",")]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 4: Collect focus gene topology across all gene-wise networks."
    )
    parser.add_argument(
        "--network-dir", type=str, required=True,
        help="Directory containing differential_network_gene_*.h5 files.",
    )
    parser.add_argument(
        "--focus-genes", type=str, required=True,
        help="Focus genes: comma-separated indices, 'top:N', or 'range:start:end'.",
    )
    parser.add_argument(
        "--n-genes", type=int, required=True,
        help="Total number of genes in the dataset.",
    )
    parser.add_argument(
        "--out-h5", type=str, default="results/focus_gene_topology.h5",
        help="Output HDF5 path.",
    )
    parser.add_argument(
        "--gene-ids", type=str, default=None,
        help="Path to file with gene IDs (one per line).",
    )
    parser.add_argument(
        "--n-jobs", type=int, default=1,
        help="Number of parallel jobs (default: 1).",
    )
    parser.add_argument(
        "--file-pattern", type=str, default="differential_network_gene_*.h5",
        help="Glob pattern for network files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Load gene IDs if provided
    gene_ids = None
    if args.gene_ids:
        with open(args.gene_ids) as f:
            gene_ids = [line.strip() for line in f]

    # Parse focus genes
    focus_genes = parse_focus_genes(args.focus_genes, args.n_genes)
    print(f"Focus genes: {len(focus_genes)} genes")
    print(f"  First 10: {focus_genes[:10]}")

    # Collect topology
    results = collect_focus_gene_topology(
        network_dir=Path(args.network_dir),
        focus_genes=focus_genes,
        n_genes=args.n_genes,
        n_jobs=args.n_jobs,
        file_pattern=args.file_pattern,
    )

    # Save results
    save_results(results, Path(args.out_h5), gene_ids)

    # Print summary
    print_summary(results, gene_ids)


if __name__ == "__main__":
    main()
