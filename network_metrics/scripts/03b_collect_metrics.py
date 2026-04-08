#!/usr/bin/env python3
"""
Stage 3b: Aggregate per-gene metric HDF5 files into a summary table.

Reads all per_gene/{gi}_{gene_id}.h5 files produced by Stage 3 and writes:
  - network_metrics_summary.h5  (structured HDF5 with per_gene/ group)
  - network_metrics.tsv         (flat TSV sorted by degree descending)

Mirrors 03b_collect_networks.py in the differential pipeline, adapted for
the single-network metric set.

Pipeline position
-----------------
Stage 3   03_per_gene_metrics.py  →  per_gene/{gi}_{gene_id}.h5  (×n_genes)
Stage 3b  THIS SCRIPT             →  network_metrics_summary.h5
                                      network_metrics.tsv
Stage 4   04_global_metrics.py    →  global_topology.h5  (optional)

Output HDF5 layout
------------------
    meta/
        n_genes, n_ref_genes, n_edges
    gene_names  (n_genes,) str dataset
    per_gene/
        gene_index  (n_ref,) int32
        gene_name   (n_ref,) str
        <all METRIC_KEYS>  (n_ref,) int32 or float32
"""

from __future__ import annotations

import argparse
from importlib.util import spec_from_file_location, module_from_spec as _mfs
from pathlib import Path

import h5py
import numpy as np


# =============================================================================
# Import METRIC_KEYS / _INT_KEYS from Stage 3 script
# =============================================================================

def _import_stage3(script_dir: Path):
    spec = spec_from_file_location("nm_stage3", script_dir / "03_per_gene_metrics.py")
    mod = _mfs(spec)
    spec.loader.exec_module(mod)
    return mod


# =============================================================================
# Read one per-gene file
# =============================================================================

def _read_gene_file(h5_path: Path, metric_keys: list[str]) -> dict | None:
    try:
        with h5py.File(h5_path, "r") as h5:
            gene_index = int(h5["meta"].attrs["gene_index"])
            result = {"gene_index": gene_index}
            met = h5.get("metrics")
            for key in metric_keys:
                result[key] = met.attrs[key] if (met is not None and key in met.attrs) else 0
            return result
    except Exception as e:
        print(f"  Warning: skipping {h5_path.name}: {e}")
        return None


# =============================================================================
# Collect
# =============================================================================

def collect_from_per_gene_dir(
    per_gene_dir: Path,
    network_h5_path: Path,
    out_h5_path: Path,
    out_tsv_path: Path | None = None,
    script_dir: Path | None = None,
) -> None:
    if script_dir is None:
        script_dir = Path(__file__).resolve().parent

    s3 = _import_stage3(script_dir)
    metric_keys = s3.METRIC_KEYS
    int_keys    = s3._INT_KEYS

    h5_files = sorted(per_gene_dir.glob("*.h5"))
    if not h5_files:
        raise FileNotFoundError(f"No .h5 files found in {per_gene_dir}")
    print(f"Found {len(h5_files)} per-gene files")

    # Global metadata from network_edges.h5
    with h5py.File(network_h5_path, "r") as h5:
        n_genes = int(h5["meta"].attrs["n_genes"])
        n_edges = int(h5["meta"].attrs["n_edges"])
        all_gene_names = (
            [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
            if "gene_names" in h5 else None
        )

    # Parse {gi:04d}_{gene_id} filenames
    gene_indices:    list[int] = []
    gene_file_names: list[str] = []
    for f in h5_files:
        parts = f.stem.split("_", 1)
        gene_indices.append(int(parts[0]))
        gene_file_names.append(parts[1] if len(parts) > 1 else f"gene_{parts[0]}")

    n_ref = len(h5_files)
    summary: dict[str, np.ndarray] = {
        k: np.zeros(n_ref, dtype=np.int32 if k in int_keys else np.float32)
        for k in metric_keys
    }
    summary["gene_index"] = np.array(gene_indices, dtype=np.int32)

    n_ok = 0
    for ref_idx, h5_path in enumerate(h5_files):
        res = _read_gene_file(h5_path, metric_keys)
        if res is None:
            continue
        n_ok += 1
        for key in metric_keys:
            if key in res:
                summary[key][ref_idx] = res[key]

    print(f"  Processed {n_ok}/{n_ref} files successfully.")

    # --- Write HDF5 ---
    print(f"Saving summary HDF5 to {out_h5_path}...")
    out_h5_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_h5_path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["n_genes"]     = n_genes
        meta.attrs["n_ref_genes"] = n_ref
        meta.attrs["n_edges"]     = n_edges
        if all_gene_names is not None:
            h5.create_dataset("gene_names", data=all_gene_names, dtype=h5py.string_dtype())
        pg = h5.create_group("per_gene")
        pg.create_dataset("gene_index", data=summary["gene_index"], compression="gzip")
        pg.create_dataset("gene_name",  data=gene_file_names, dtype=h5py.string_dtype())
        for key in metric_keys:
            dtype = np.int32 if key in int_keys else np.float32
            pg.create_dataset(key, data=summary[key].astype(dtype), compression="gzip")
    print(f"  Saved {n_ref} genes.")

    # --- Print summary ---
    n_active = int((summary["degree"] > 0).sum())
    print(f"\n{'='*60}")
    print("NETWORK METRICS COLLECTION SUMMARY")
    print(f"{'='*60}")
    print(f"  Reference genes   : {n_ref}")
    print(f"  Genes with degree>0: {n_active}")
    if n_active > 0:
        active = summary["degree"] > 0
        print(f"  Mean degree (active)    : {summary['degree'][active].mean():.1f}")
        print(f"  Max degree              : {int(summary['degree'].max())}")
        print(f"  Mean L2L1_deg (active)  : {summary['L2L1_deg'][active].mean():.3f}")
        print(f"  Mean clique density     : {summary['L1_clique_density'][active].mean():.4f}")
        top = np.argsort(summary["L2L1_deg"])[::-1][:20]
        print(f"\n  Top genes by L2L1_deg:")
        print(f"  {'Idx':>6} {'Gene':>22} {'Degree':>7} "
              f"{'L1':>5} {'L2':>5} {'L2L1':>7} {'Clique':>7} {'L2L1c':>7}")
        for i in top:
            if summary["degree"][i] > 0:
                print(f"  {int(summary['gene_index'][i]):>6} {gene_file_names[i]:>22} "
                      f"{int(summary['degree'][i]):>7} "
                      f"{int(summary['L1_n_nodes'][i]):>5} "
                      f"{int(summary['L2_n_nodes'][i]):>5} "
                      f"{summary['L2L1_deg'][i]:>7.3f} "
                      f"{summary['L1_clique_density'][i]:>7.4f} "
                      f"{summary['L2L1_conn'][i]:>7.3f}")

    # --- Write TSV ---
    if out_tsv_path:
        print(f"\nSaving TSV to {out_tsv_path}...")
        out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
        sort_idx = np.argsort(summary["degree"])[::-1]
        cols = ["gene_idx", "gene_id"] + metric_keys
        with open(out_tsv_path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for i in sort_idx:
                row = [str(int(summary["gene_index"][i])), gene_file_names[i]]
                for key in metric_keys:
                    v = summary[key][i]
                    row.append(str(int(v)) if key in int_keys else f"{float(v):.6g}")
                f.write("\t".join(row) + "\n")
        print(f"  Saved {n_ref} rows.")


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="[Network Metrics Stage 3b] Collect per-gene HDF5s into summary."
    )
    p.add_argument("--per-gene-dir", type=str, required=True,
                   help="Directory of per_gene/{gi}_{gene_id}.h5 files.")
    p.add_argument("--network-h5",   type=str, required=True,
                   help="network_edges.h5 from Stage 2 (for n_genes, gene_names).")
    p.add_argument("--out-h5",       type=str, required=True,
                   help="Output summary HDF5 path.")
    p.add_argument("--out-tsv",      type=str, default=None,
                   help="Output TSV path (optional).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    collect_from_per_gene_dir(
        per_gene_dir=Path(args.per_gene_dir),
        network_h5_path=Path(args.network_h5),
        out_h5_path=Path(args.out_h5),
        out_tsv_path=Path(args.out_tsv) if args.out_tsv else None,
    )


if __name__ == "__main__":
    main()
