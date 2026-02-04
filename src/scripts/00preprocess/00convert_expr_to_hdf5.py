#!/usr/bin/env python3
"""
Preprocessing: Convert expression TSV to HDF5 for fast array-job access.

Usage
-----
    python 00convert_expr_to_hdf5.py --expr-tsv data/expression.tsv --out-h5 data/expression.h5

Benefits
--------
* Parse TSV once instead of N_GENES times
* Memory-mapped access (~instant load vs. seconds of pandas parsing)
* Compression reduces disk usage while maintaining fast access

Output schema
-------------
expression.h5
    /expr              (n_genes, n_samples) float32 dataset
    /expr.attrs        metadata: gene_names, sample_names, shape info
"""

import argparse
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def convert_tsv_to_hdf5(
    tsv_path: Path,
    h5_path: Path,
    compression: str = "gzip",
    compression_opts: int = 4,
) -> None:
    """
    Convert expression TSV to HDF5 with metadata preservation.
    
    Parameters
    ----------
    tsv_path : Path
        Input TSV file (genes as rows, samples as columns)
    h5_path : Path
        Output HDF5 file path
    compression : str
        HDF5 compression algorithm (gzip, lzf, or None)
    compression_opts : int
        Compression level for gzip (1-9, lower=faster, higher=smaller)
    """
    print(f"Reading TSV: {tsv_path}")
    t0 = time.time()
    
    df = pd.read_csv(tsv_path, sep="\t", index_col=0)
    expr = df.to_numpy(dtype=np.float32)
    
    t_read = time.time() - t0
    print(f"  ✓ Loaded {expr.shape[0]:,} genes × {expr.shape[1]:,} samples in {t_read:.2f}s")
    print(f"  Memory footprint: {expr.nbytes / 1024**2:.1f} MB")
    
    # Write HDF5
    print(f"\nWriting HDF5: {h5_path}")
    t0 = time.time()
    
    with h5py.File(h5_path, "w") as f:
        # Main dataset with optional compression
        ds = f.create_dataset(
            "expr",
            data=expr,
            dtype=np.float32,
            compression=compression if compression != "none" else None,
            compression_opts=compression_opts if compression == "gzip" else None,
            chunks=True,  # Auto-chunking for good performance
        )
        
        # Store metadata as attributes
        ds.attrs["n_genes"] = expr.shape[0]
        ds.attrs["n_samples"] = expr.shape[1]
        ds.attrs["gene_names"] = df.index.astype(str).tolist()
        ds.attrs["sample_names"] = df.columns.astype(str).tolist()
        
        # Optional: store gene/sample names as separate datasets for easier access
        # (attributes have size limits, though usually fine for typical datasets)
        f.create_dataset("gene_names", data=df.index.astype(str).tolist(), dtype=h5py.string_dtype())
        f.create_dataset("sample_names", data=df.columns.astype(str).tolist(), dtype=h5py.string_dtype())
    
    t_write = time.time() - t0
    h5_size = h5_path.stat().st_size / 1024**2
    
    print(f"  ✓ Written in {t_write:.2f}s")
    print(f"  File size: {h5_size:.1f} MB")
    if compression:
        compression_ratio = expr.nbytes / h5_path.stat().st_size
        print(f"  Compression ratio: {compression_ratio:.2f}x")
    
    # Quick validation
    print(f"\nValidating...")
    with h5py.File(h5_path, "r") as f:
        assert f["expr"].shape == expr.shape
        assert np.allclose(f["expr"][0, :10], expr[0, :10], rtol=1e-6)
    print(f"  ✓ Validation passed")
    
    print(f"\n{'='*60}")
    print(f"SUCCESS: {h5_path}")
    print(f"Speedup estimate: {t_read:.1f}s → ~0.1s per array job")
    print(f"Total time saved for N jobs: {t_read:.1f}s × N array tasks")
    print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert expression TSV to HDF5 for fast parallel access."
    )
    parser.add_argument(
        "--expr-tsv",
        type=str,
        required=True,
        help="Input expression TSV (genes as rows, samples as columns)",
    )
    parser.add_argument(
        "--out-h5",
        type=str,
        required=True,
        help="Output HDF5 file path",
    )
    parser.add_argument(
        "--compression",
        type=str,
        default="gzip",
        choices=["gzip", "lzf", "none"],
        help="Compression algorithm (default: gzip). Use 'none' for fastest access.",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=4,
        choices=range(1, 10),
        help="GZIP compression level 1-9 (default: 4). Lower=faster, higher=smaller.",
    )
    
    args = parser.parse_args()
    
    tsv_path = Path(args.expr_tsv)
    h5_path = Path(args.out_h5)
    
    if not tsv_path.exists():
        raise FileNotFoundError(f"TSV not found: {tsv_path}")
    
    # Create output directory if needed
    h5_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert
    convert_tsv_to_hdf5(
        tsv_path=tsv_path,
        h5_path=h5_path,
        compression=args.compression,
        compression_opts=args.compression_level if args.compression == "gzip" else None,
    )


if __name__ == "__main__":
    main()
