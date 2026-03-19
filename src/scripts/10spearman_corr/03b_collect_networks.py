#!/usr/bin/env python3
"""
Stage 3b: Collect Per-Gene Differential Networks into a Summary.

Aggregates per-gene differential_network.h5 files (from Stage 3,
reconstruct_single Snakemake rule) into a single summary HDF5 and an
optional TSV of rewiring hub metrics.

Reads either:
  --networks-dir  Directory of per-gene differential_network.h5 files
                  written by Stage 3 (preferred, avoids re-computation).
  --base-dir /    Fallback: raw per-gene correlation + bootstrap directories
  --boot-dir      (equivalent to old Stage 3 per-gene batch mode).

Pipeline position
-----------------
Stage 3   03_reconstruct_diff_network.py        →  networks/{gi}_{gene_id}.h5  (×N genes)
Stage 3b  THIS SCRIPT                           →  differential_network_summary.h5
                                                    rewiring_hubs.tsv
Stage 4   04_collect_focus_gene_topology.py     →  focus_gene_topology.h5

Output HDF5 layout (summary.h5)
--------------------------------
    meta/
        n_genes, n_ref_genes, edge_selection, min_effect,
        corr_threshold, require_ci_exclude_zero
    gene_names    (n_genes,) str
    per_gene/
        gene_index          (n_ref,) int32
        gene_name           (n_ref,) str
        n_sig_edges_diff    (n_ref,) int32   # # significant differential edges
        n_disappear, n_new, n_sign_change, n_strengthen, n_weaken  (n_ref,) int32
        mean_abs_delta      (n_ref,) float32
        max_abs_delta       (n_ref,) float32
        focus_deg_low       (n_ref,) int32
        focus_deg_high      (n_ref,) int32
        L1_n_nodes, L1_rewire, L1_frac_rewire, L1_clique_density, L1_conn_low, ...   (n_ref,) float32 / int32
        L2_n_nodes, L2_n_edges, L2_rewire, L2_conn_low, ...   (n_ref,) float32 / int32
        L2L1_deg, L2L1_rewire, L2L1_conn           (n_ref,) float32
        HL_conn_L1, HL_conn_L2                     (n_ref,) float32
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
from functools import partial
from pathlib import Path
from collections import defaultdict

import h5py
import numpy as np

# Import shared computation from Stage 3 script
from importlib.util import spec_from_file_location, module_from_spec as _mfs

def _import_stage3():
    """Lazily import from 03_reconstruct_diff_network.py in the same directory."""
    script_dir = Path(__file__).resolve().parent
    spec = spec_from_file_location(
        "stage3", script_dir / "03_reconstruct_diff_network.py"
    )
    mod = _mfs(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Mode A: Read from per-gene differential_network.h5 files (Stage 3 output)
# ---------------------------------------------------------------------------

_SUMMARY_KEYS = [
    "n_sig_edges_diff",
    "n_disappear", "n_new", "n_sign_change", "n_strengthen", "n_weaken",
    "mean_abs_delta", "max_abs_delta",
    "focus_deg_low", "focus_deg_high",
    "L1_n_nodes", "L1_n_disappear", "L1_n_new", "L1_n_sign_chg",
    "L1_n_strengthen", "L1_n_weaken", "L1_rewire",
    "L1_frac_disappear", "L1_frac_new", "L1_frac_rewire",
    "L1_conn_low", "L1_conn_high", "L1_conn_diff",
    "L1_conn_mean_low", "L1_conn_mean_high", "L1_mean_abs_dr", "L1_mean_delta",
    "L1_clique_density",
    "L2_n_nodes", "L2_n_edges",
    "L2_n_disappear", "L2_n_new", "L2_n_sign_chg",
    "L2_n_strengthen", "L2_n_weaken", "L2_rewire",
    "L2_conn_low", "L2_conn_high", "L2_conn_diff",
    "L2_conn_mean_low", "L2_conn_mean_high", "L2_mean_abs_dr", "L2_mean_delta",
    "full_n_edges", "full_rewire",
    "full_conn_low", "full_conn_high", "full_conn_diff",
    "full_mean_abs_dr", "full_mean_delta",
    "n_l1_to_l1_edges", "n_l1_to_l2_edges",
    "L2L1_deg", "L2L1_rewire", "L2L1_conn",
    "HL_conn_L1", "HL_conn_L2",
]

_INT_KEYS = {
    "n_sig_edges_diff",
    "n_disappear", "n_new", "n_sign_change", "n_strengthen", "n_weaken",
    "focus_deg_low", "focus_deg_high",
    "L1_n_nodes", "L1_n_disappear", "L1_n_new", "L1_n_sign_chg",
    "L1_n_strengthen", "L1_n_weaken", "L1_rewire",
    "L2_n_nodes", "L2_n_edges",
    "L2_n_disappear", "L2_n_new", "L2_n_sign_chg",
    "L2_n_strengthen", "L2_n_weaken", "L2_rewire",
    "full_n_edges", "full_rewire",
    "n_l1_to_l1_edges", "n_l1_to_l2_edges",
}


def _read_gene_summary_from_network_file(h5_path: Path) -> dict | None:
    """
    Extract per-gene summary metrics from a single Stage 3 output HDF5.
    Returns a flat dict of metrics, or None if the file has no significant edges.
    """
    try:
        with h5py.File(h5_path, "r") as h5:
            n_sig = int(h5["meta"].attrs.get("n_significant", 0))
            gene_index = int(h5["meta"].attrs.get("gene_index_used", 0))

            if n_sig == 0:
                return {"gene_index": gene_index, "n_sig_edges_diff": 0}

            result = {
                "gene_index": gene_index,
                "n_sig_edges_diff": n_sig,
            }

            # Qualitative counts from edges/qual_summary
            if "edges/qual_summary" in h5:
                qs = h5["edges/qual_summary"]
                for label in ["disappear", "new", "sign_change", "strengthen", "weaken"]:
                    result[f"n_{label}"] = int(qs.attrs.get(f"n_{label}", 0))
            else:
                for label in ["disappear", "new", "sign_change", "strengthen", "weaken"]:
                    result[f"n_{label}"] = 0

            # Delta stats from edges
            delta = h5["edges/delta_base"][:]
            result["mean_abs_delta"] = float(np.mean(np.abs(delta)))
            result["max_abs_delta"] = float(np.max(np.abs(delta)))

            # Focus gene degree per condition
            result["focus_deg_low"] = 0
            result["focus_deg_high"] = 0
            if "focus_gene" in h5:
                fg = h5["focus_gene"]
                if "metrics" in fg:
                    met = fg["metrics"]
                    for key in _SUMMARY_KEYS:
                        if key in met.attrs:
                            result[key] = met.attrs[key]
                # focus_deg_low/high are stored directly in focus_gene group
                # (not in metrics/); they are computed per-condition
                if "direct_stats" in fg:
                    pass  # already loaded via metrics/

                # focus degree per condition comes from analyze_focus_gene flat keys
                # stored under focus_gene/metrics/ by save_results
                if "focus_deg_low" not in result:
                    result["focus_deg_low"] = 0
                if "focus_deg_high" not in result:
                    result["focus_deg_high"] = 0

            return result

    except Exception as e:
        print(f"  Warning: failed to read {h5_path.name}: {e}")
        return None


def collect_from_network_dir(
    networks_dir: Path,
    out_h5_path: Path,
    out_focus_tsv_path: Path | None = None,
    annotate: bool = False,
) -> None:
    """
    Collect summary from per-gene differential_network.h5 files in networks_dir.
    These are the files written by Stage 3 reconstruct_single rule.
    """
    h5_files = sorted(networks_dir.glob("*.h5"))
    if not h5_files:
        raise FileNotFoundError(f"No .h5 files found in {networks_dir}")

    print(f"Found {len(h5_files)} network files in {networks_dir}")

    # Read n_genes and gene_names from first file
    with h5py.File(h5_files[0], "r") as h5:
        n_genes = int(h5["meta"].attrs["n_genes"])
        if "gene_names" in h5:
            all_gene_names = [
                x.decode() if isinstance(x, bytes) else x
                for x in h5["gene_names"][:]
            ]
        else:
            all_gene_names = None

    # Extract gene index and name from filenames: {gi}_{gene_id}.h5
    gene_indices = []
    gene_file_names = []
    for f in h5_files:
        parts = f.stem.split("_", 1)
        gene_indices.append(int(parts[0]))
        gene_file_names.append(parts[1] if len(parts) > 1 else f"gene_{parts[0]}")

    n_ref = len(h5_files)

    # Initialize summary arrays
    summary = {k: np.zeros(n_ref, dtype=np.int32 if k in _INT_KEYS else np.float32)
               for k in _SUMMARY_KEYS}
    summary["gene_index"] = np.array(gene_indices, dtype=np.int32)

    n_processed = 0
    n_with_edges = 0
    for ref_idx, h5_path in enumerate(h5_files):
        res = _read_gene_summary_from_network_file(h5_path)
        if res is None:
            print(f"  {h5_path.name}: skipped (read error)")
            continue

        n_processed += 1
        n_sig = res.get("n_sig_edges_diff", 0)
        summary["n_sig_edges_diff"][ref_idx] = n_sig

        if n_sig == 0:
            continue

        n_with_edges += 1
        for key in _SUMMARY_KEYS:
            if key in res:
                summary[key][ref_idx] = res[key]

    print(f"  Processed {n_processed} genes ({n_with_edges} with significant edges)")

    _write_summary_h5(
        out_h5_path, summary, gene_file_names, all_gene_names,
        n_genes=n_genes,
        edge_selection="from_network_files",
        min_effect=0.0,
        corr_threshold=0.0,
        require_ci_exclude_zero=True,
    )

    _print_summary(summary, gene_file_names, n_ref)

    if out_focus_tsv_path:
        _write_focus_tsv(out_focus_tsv_path, summary, gene_file_names)
        if annotate:
            _annotate_tsv(out_focus_tsv_path)


# ---------------------------------------------------------------------------
# Mode B: Read from raw base/boot directories (fallback, uses Stage 3 logic)
# ---------------------------------------------------------------------------

def _find_gene_file(directory: Path, gene_index: int) -> Path | None:
    """Find per-gene h5 file by index prefix. Returns None if not found."""
    prefix = f"{gene_index:04d}_"
    matches = list(directory.glob(f"{prefix}*.h5"))
    return matches[0] if len(matches) == 1 else None


def _worker(
    ref_idx: int,
    g_idx: int,
    g_name: str,
    base_path: Path,
    boot_dir: Path,
    n_genes: int,
    edge_selection: str,
    min_effect: float,
    require_ci_exclude_zero: bool,
    corr_threshold: float,
) -> tuple[int, dict | None]:
    """Worker for parallel per-gene processing from raw base/boot files."""
    s3 = _import_stage3()
    try:
        boot_path = _find_gene_file(boot_dir, g_idx)
        if boot_path is None:
            return (ref_idx, None)

        results = s3.load_and_process(
            base_h5_path=base_path,
            boot_h5_path=boot_path,
            edge_selection=edge_selection,
            min_effect=min_effect,
            require_ci_exclude_zero=require_ci_exclude_zero,
            corr_threshold=corr_threshold,
        )

        n_sig = results["n_significant"]
        if n_sig == 0:
            return (ref_idx, {"n_sig_edges_diff": 0})

        qs = results["qual_score"]
        qual_counts = {
            "n_disappear": int(np.sum(qs == s3.QUAL_DISAPPEAR)),
            "n_new": int(np.sum(qs == s3.QUAL_NEW)),
            "n_sign_change": int(np.sum(qs == s3.QUAL_SIGN_CHANGE)),
            "n_strengthen": int(np.sum(qs == s3.QUAL_STRENGTHEN)),
            "n_weaken": int(np.sum(qs == s3.QUAL_WEAKEN)),
        }

        abs_delta = np.abs(results["delta_base"])
        delta_stats = {
            "mean_abs_delta": float(np.mean(abs_delta)),
            "max_abs_delta": float(np.max(abs_delta)),
        }

        fa = s3.analyze_focus_gene(
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

        focus_metrics = {k: fa[k] for k in [
            "focus_deg_low", "focus_deg_high",
            "L1_n_nodes", "L1_n_disappear", "L1_n_new", "L1_n_sign_chg",
            "L1_n_strengthen", "L1_n_weaken", "L1_rewire",
            "L1_frac_disappear", "L1_frac_new", "L1_frac_rewire",
            "L1_conn_low", "L1_conn_high", "L1_conn_diff",
            "L1_conn_mean_low", "L1_conn_mean_high", "L1_mean_abs_dr", "L1_mean_delta",
            "L1_clique_density",
            "L2_n_nodes", "L2_n_edges",
            "L2_n_disappear", "L2_n_new", "L2_n_sign_chg",
            "L2_n_strengthen", "L2_n_weaken", "L2_rewire",
            "L2_conn_low", "L2_conn_high", "L2_conn_diff",
            "L2_conn_mean_low", "L2_conn_mean_high", "L2_mean_abs_dr", "L2_mean_delta",
            "full_n_edges", "full_rewire",
            "full_conn_low", "full_conn_high", "full_conn_diff",
            "full_mean_abs_dr", "full_mean_delta",
            "n_l1_to_l1_edges", "n_l1_to_l2_edges",
            "L2L1_deg", "L2L1_rewire", "L2L1_conn",
            "HL_conn_L1", "HL_conn_L2",
        ]}

        return (ref_idx, {
            "n_sig_edges_diff": n_sig,
            **qual_counts,
            **delta_stats,
            **focus_metrics,
        })

    except (ValueError, KeyError):
        return (ref_idx, None)


def collect_from_raw_dirs(
    base_dir: Path,
    boot_dir: Path,
    out_h5_path: Path,
    out_focus_tsv_path: Path | None = None,
    edge_selection: str = "sig_edges",
    min_effect: float = 0.0,
    require_ci_exclude_zero: bool = True,
    corr_threshold: float = 0.0001,
    annotate: bool = False,
) -> None:
    """
    Collect per-gene networks from raw base_correlations and bootstrap_significant
    directories, computing topology on-the-fly (fallback mode, no Stage 3 outputs needed).
    """
    base_files = sorted(base_dir.glob("*_*.h5"))
    boot_files = sorted(boot_dir.glob("*_*.h5"))

    if not base_files:
        raise FileNotFoundError(f"No per-gene h5 files found in {base_dir}")

    gene_indices = []
    gene_file_names = []
    for f in base_files:
        parts = f.stem.split("_", 1)
        gene_indices.append(int(parts[0]))
        gene_file_names.append(parts[1] if len(parts) > 1 else f"gene_{parts[0]}")

    n_ref = len(gene_indices)
    print(f"Found {n_ref} per-gene base_correlations files")
    print(f"Found {len(boot_files)} per-gene bootstrap_significant files")

    with h5py.File(base_files[0], "r") as h5:
        n_genes = int(h5["meta"].attrs["n_genes"])
        if "gene_names" in h5:
            all_gene_names = [
                x.decode() if isinstance(x, bytes) else x
                for x in h5["gene_names"][:]
            ]
        else:
            all_gene_names = None

    print(f"Network size: {n_genes} genes")

    summary = {k: np.zeros(n_ref, dtype=np.int32 if k in _INT_KEYS else np.float32)
               for k in _SUMMARY_KEYS}
    summary["gene_index"] = np.array(gene_indices, dtype=np.int32)

    n_cpus = max(1, mp.cpu_count() - 2)
    n_workers = min(n_cpus, 16, n_ref)
    print(f"  Using {n_workers} parallel workers for {n_ref} genes...")

    worker_args = [
        (ref_idx, gene_indices[ref_idx], gene_file_names[ref_idx],
         base_files[ref_idx], boot_dir, n_genes,
         edge_selection, min_effect, require_ci_exclude_zero, corr_threshold)
        for ref_idx in range(n_ref)
    ]

    if n_workers > 1:
        with mp.Pool(processes=n_workers) as pool:
            results_list = pool.starmap(_worker, worker_args)
    else:
        results_list = [_worker(*args) for args in worker_args]

    n_processed = 0
    n_with_edges = 0
    for ref_idx, gene_result in results_list:
        if gene_result is None:
            print(f"  Gene {gene_indices[ref_idx]} ({gene_file_names[ref_idx]}): skipped")
            continue
        n_processed += 1
        n_sig = gene_result.get("n_sig_edges_diff", 0)
        summary["n_sig_edges_diff"][ref_idx] = n_sig
        if n_sig == 0:
            continue
        n_with_edges += 1
        for key in _SUMMARY_KEYS:
            if key in gene_result:
                summary[key][ref_idx] = gene_result[key]

    print(f"  Processed {n_processed} genes ({n_with_edges} with significant edges)")

    _write_summary_h5(
        out_h5_path, summary, gene_file_names, all_gene_names,
        n_genes=n_genes,
        edge_selection=edge_selection,
        min_effect=min_effect,
        corr_threshold=corr_threshold,
        require_ci_exclude_zero=require_ci_exclude_zero,
    )

    _print_summary(summary, gene_file_names, n_ref)

    if out_focus_tsv_path:
        _write_focus_tsv(out_focus_tsv_path, summary, gene_file_names)
        if annotate:
            _annotate_tsv(out_focus_tsv_path)


# ---------------------------------------------------------------------------
# Shared output helpers
# ---------------------------------------------------------------------------

def _write_summary_h5(
    out_h5_path: Path,
    summary: dict,
    gene_file_names: list,
    all_gene_names: list | None,
    n_genes: int,
    edge_selection: str,
    min_effect: float,
    corr_threshold: float,
    require_ci_exclude_zero: bool,
) -> None:
    print(f"\nSaving summary to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)
    n_ref = len(gene_file_names)

    with h5py.File(out_h5_path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["n_genes"] = n_genes
        meta.attrs["n_ref_genes"] = n_ref
        meta.attrs["edge_selection"] = edge_selection
        meta.attrs["min_effect"] = min_effect
        meta.attrs["corr_threshold"] = corr_threshold
        meta.attrs["require_ci_exclude_zero"] = require_ci_exclude_zero

        if all_gene_names is not None:
            h5.create_dataset("gene_names", data=all_gene_names, dtype=h5py.string_dtype())

        per_gene = h5.create_group("per_gene")
        for key, arr in summary.items():
            per_gene.create_dataset(key, data=arr, compression="gzip")
        per_gene.create_dataset("gene_name", data=gene_file_names, dtype=h5py.string_dtype())

    print(f"  Saved {n_ref} reference genes to {out_h5_path}")


def _print_summary(summary: dict, gene_file_names: list, n_ref: int) -> None:
    n_with_edges = int(np.sum(summary["n_sig_edges_diff"] > 0))
    print(f"\n{'='*60}")
    print("PER-GENE NETWORK COLLECTION SUMMARY")
    print(f"{'='*60}")
    print(f"  Reference genes processed: {n_ref}")
    print(f"  Reference genes with edges: {n_with_edges}")
    print(f"  Total sig diff edges across all genes: {summary['n_sig_edges_diff'].sum():,}")
    if n_with_edges > 0:
        active = summary["n_sig_edges_diff"] > 0
        print(f"  Mean sig diff edges per active gene: {summary['n_sig_edges_diff'][active].mean():.1f}")
        print(f"  Max sig diff edges: {summary['n_sig_edges_diff'].max():,}")
        print(f"  Mean |Δr| (active genes): {summary['mean_abs_delta'][active].mean():.4f}")

    if n_with_edges > 0 and "L2L1_deg" in summary:
        top_idx = np.argsort(summary["L2L1_deg"])[::-1][:min(20, n_ref)]
        print(f"\n  Top genes by L2/L1 node ratio (diff network):")
        print(f"  {'Idx':>6} {'Gene':>20} {'DiffEdges':>10} {'L1_nodes':>8} "
              f"{'L2_nodes':>8} {'L2L1':>7} {'L1_rew':>7} {'L2_rew':>7}")
        for i in top_idx:
            if summary["n_sig_edges_diff"][i] > 0:
                print(f"  {summary['gene_index'][i]:>6} {gene_file_names[i]:>20} "
                      f"{summary['n_sig_edges_diff'][i]:>10} {summary['L1_n_nodes'][i]:>8} "
                      f"{summary['L2_n_nodes'][i]:>8} {summary['L2L1_deg'][i]:>7.2f} "
                      f"{summary['L1_rewire'][i]:>7} {summary['L2_rewire'][i]:>7}")


_FOCUS_TSV_COLS = [
    "gene_idx", "gene_id", "n_sig_edges_diff",
    # Focus degree per condition (among differential edges only)
    "focus_deg_low", "focus_deg_high",
    # L1: direct partners of focus gene in diff network
    "L1_n_nodes",
    "L1_n_disappear", "L1_n_new", "L1_n_sign_chg",
    "L1_n_strengthen", "L1_n_weaken", "L1_rewire",
    "L1_frac_disappear", "L1_frac_new", "L1_frac_rewire",
    "L1_conn_low", "L1_conn_high", "L1_conn_diff",
    "L1_conn_mean_low", "L1_conn_mean_high",
    "L1_mean_abs_dr", "L1_mean_delta",
    "L1_clique_density",
    # L2: pure L2 nodes; outer-layer edges (L1↔L1 + L1→L2)
    "L2_n_nodes", "L2_n_edges",
    "L2_n_disappear", "L2_n_new", "L2_n_sign_chg",
    "L2_n_strengthen", "L2_n_weaken", "L2_rewire",
    "L2_conn_low", "L2_conn_high", "L2_conn_diff",
    "L2_conn_mean_low", "L2_conn_mean_high",
    "L2_mean_abs_dr", "L2_mean_delta",
    # Full 2-layer neighbourhood (direct + L1↔L1 + L1→L2 combined)
    "full_n_edges", "full_rewire",
    "full_conn_low", "full_conn_high", "full_conn_diff",
    "full_mean_abs_dr", "full_mean_delta",
    # Edge sublayer breakdown
    "n_l1_to_l1_edges", "n_l1_to_l2_edges",
    # Ratios
    "L2L1_deg", "L2L1_rewire", "L2L1_conn",
    "HL_conn_L1", "HL_conn_L2",
]


def _write_focus_tsv(
    out_path: Path,
    summary: dict,
    gene_file_names: list,
) -> None:
    print(f"\nSaving focus gene TSV to {out_path}...")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sort_idx = np.argsort(summary["L2L1_deg"])[::-1]

    with open(out_path, "w") as f:
        f.write("\t".join(_FOCUS_TSV_COLS) + "\n")
        for i in sort_idx:
            if summary["n_sig_edges_diff"][i] == 0:
                continue
            vals = []
            for col in _FOCUS_TSV_COLS:
                if col == "gene_idx":
                    vals.append(str(summary["gene_index"][i]))
                elif col == "gene_id":
                    vals.append(gene_file_names[i])
                elif col in _INT_KEYS or col in {"n_sig_edges_diff", "gene_idx"}:
                    vals.append(str(int(summary[col][i])))
                else:
                    vals.append(f"{summary[col][i]:.4f}")
            f.write("\t".join(vals) + "\n")

    n_written = int(np.sum(summary["n_sig_edges_diff"] > 0))
    print(f"  Saved {n_written} genes to focus gene TSV")


def _annotate_tsv(tsv_path: Path, gene_id_col: str = "gene_id") -> None:
    """Annotate TSV with gene symbols using 06_annotate_rewiring_table.R."""
    import subprocess

    script_dir = Path(__file__).resolve().parent
    r_script = script_dir / "06_annotate_rewiring_table.R"

    if not r_script.exists():
        print(f"  WARNING: Annotation script not found at {r_script}")
        return

    out_path = tsv_path.with_name(tsv_path.stem + "_annotated.tsv")
    cmd = [
        "Rscript", str(r_script),
        "--input-tsv", str(tsv_path),
        "--output-tsv", str(out_path),
        "--gene-id-col", gene_id_col,
    ]
    print(f"\nAnnotating {tsv_path.name} with gene names...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print(f"  Annotated TSV saved to {out_path}")
        else:
            print(f"  WARNING: Annotation failed:\n{result.stderr}")
    except FileNotFoundError:
        print("  WARNING: Rscript not found. Skipping annotation.")
    except subprocess.TimeoutExpired:
        print("  WARNING: Annotation timed out after 120s.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 3b: Collect per-gene differential networks into a summary."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--networks-dir", type=str, default=None,
        help=(
            "Directory of per-gene differential_network.h5 files from Stage 3 "
            "(reconstruct_single Snakemake rule).  Preferred over --base-dir/--boot-dir."
        ),
    )
    mode.add_argument(
        "--base-dir", type=str, default=None,
        help="Fallback: directory of per-gene base_correlations h5 files.",
    )
    parser.add_argument(
        "--boot-dir", type=str, default=None,
        help="Required when --base-dir is used: bootstrap_significant directory.",
    )
    parser.add_argument(
        "--out-h5", type=str, default="results/differential_network_summary.h5",
        help="Output summary HDF5 path.",
    )
    parser.add_argument(
        "--out-focus-tsv", type=str, default=None,
        help="Output TSV path for rewiring hubs table.",
    )
    parser.add_argument(
        "--edge-selection", type=str, choices=["sig_edges", "sig_differential"],
        default="sig_differential",
        help="Edge selection mode (used only with --base-dir/--boot-dir).",
    )
    parser.add_argument(
        "--min-effect", type=float, default=0.0,
        help="Minimum |Δr| threshold (used only with --base-dir/--boot-dir).",
    )
    parser.add_argument(
        "--corr-threshold", type=float, default=0.0001,
        help="Correlation threshold for qualitative 'present'.",
    )
    parser.add_argument(
        "--no-ci-filter", action="store_true",
        help="Don't require bootstrap CI to exclude 0 (--base-dir/--boot-dir mode).",
    )
    parser.add_argument(
        "--annotate", action="store_true",
        help="Annotate the output TSV using 06_annotate_rewiring_table.R (requires R).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    out_h5 = Path(args.out_h5)
    out_tsv = Path(args.out_focus_tsv) if args.out_focus_tsv else None

    if args.networks_dir:
        collect_from_network_dir(
            networks_dir=Path(args.networks_dir),
            out_h5_path=out_h5,
            out_focus_tsv_path=out_tsv,
            annotate=args.annotate,
        )
    else:
        if not args.boot_dir:
            raise SystemExit("ERROR: --boot-dir is required when using --base-dir.")
        collect_from_raw_dirs(
            base_dir=Path(args.base_dir),
            boot_dir=Path(args.boot_dir),
            out_h5_path=out_h5,
            out_focus_tsv_path=out_tsv,
            edge_selection=args.edge_selection,
            min_effect=args.min_effect,
            require_ci_exclude_zero=not args.no_ci_filter,
            corr_threshold=args.corr_threshold,
            annotate=args.annotate,
        )


if __name__ == "__main__":
    main()
