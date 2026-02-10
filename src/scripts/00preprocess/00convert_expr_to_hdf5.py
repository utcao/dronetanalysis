#!/usr/bin/env python3
"""
Preprocessing: Convert expression TSV to HDF5 for fast array-job access.

Usage
-----
    python 00convert_expr_to_hdf5.py --expr-tsv data/test/voomdataCtrl_test.txt --out-h5 data/expression.h5

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


def generate_toy_data(h5_path: Path, n_genes: int = 10, n_samples: int = 60, seed: int = 42) -> None:
    """
    Generate toy expression data with known differential co-expression structure.

    Design (samples split by gene 0 expression, low/high 20%):
      - Gene 0: gradient gene (defines the low/high split)
      - Genes 1-2: correlated (r~0.8) in LOW, uncorrelated in HIGH  -> "disappear"
      - Genes 3-4: uncorrelated in LOW, correlated (r~0.8) in HIGH  -> "new"
      - Genes 5-6: positive in LOW, negative in HIGH                -> "sign_change"
      - Genes 7-8: weakly correlated in both                        -> "unchanged"
      - Gene 9: noise
    """
    rng = np.random.default_rng(seed)
    expr = np.zeros((n_genes, n_samples), dtype=np.float32)

    # Gene 0: gradient (sorted so low=first 20%, high=last 20%)
    expr[0] = np.linspace(-3, 3, n_samples) + rng.normal(0, 0.1, n_samples)

    # Determine low/high sample boundaries (20% each)
    k_low = int(n_samples * 0.2)   # 12 samples
    k_high = int(n_samples * 0.2)  # 12 samples
    low_idx = np.arange(k_low)
    high_idx = np.arange(n_samples - k_high, n_samples)
    mid_idx = np.arange(k_low, n_samples - k_high)

    # Helper: generate correlated pair for a subset
    def corr_pair(idx, r=0.8):
        n = len(idx)
        z = rng.normal(size=n)
        a = z
        b = r * z + np.sqrt(1 - r**2) * rng.normal(size=n)
        return a.astype(np.float32), b.astype(np.float32)

    # Genes 1-2: correlated in LOW, uncorrelated in HIGH
    expr[1] = rng.normal(size=n_samples)
    expr[2] = rng.normal(size=n_samples)
    a, b = corr_pair(low_idx, r=0.85)
    expr[1, low_idx] = a
    expr[2, low_idx] = b

    # Genes 3-4: uncorrelated in LOW, correlated in HIGH
    expr[3] = rng.normal(size=n_samples)
    expr[4] = rng.normal(size=n_samples)
    a, b = corr_pair(high_idx, r=0.85)
    expr[3, high_idx] = a
    expr[4, high_idx] = b

    # Genes 5-6: positive in LOW (r~0.7), negative in HIGH (r~-0.7)
    expr[5] = rng.normal(size=n_samples)
    expr[6] = rng.normal(size=n_samples)
    a, b = corr_pair(low_idx, r=0.7)
    expr[5, low_idx] = a
    expr[6, low_idx] = b
    a, b = corr_pair(high_idx, r=-0.7)
    expr[5, high_idx] = a
    expr[6, high_idx] = b

    # Genes 7-8: weakly correlated in both (r~0.3)
    a_all, b_all = corr_pair(np.arange(n_samples), r=0.3)
    expr[7] = a_all
    expr[8] = b_all

    # Gene 9: pure noise
    expr[9] = rng.normal(size=n_samples).astype(np.float32)

    # Fill remaining genes if n_genes > 10
    for i in range(10, n_genes):
        expr[i] = rng.normal(size=n_samples).astype(np.float32)

    gene_names = [f"toy_gene_{i}" for i in range(n_genes)]
    sample_names = [f"sample_{i}" for i in range(n_samples)]

    # Write HDF5
    h5_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(h5_path, "w") as f:
        ds = f.create_dataset("expr", data=expr, dtype=np.float32)
        ds.attrs["n_genes"] = n_genes
        ds.attrs["n_samples"] = n_samples
        ds.attrs["gene_names"] = gene_names
        ds.attrs["sample_names"] = sample_names
        f.create_dataset("gene_names", data=gene_names, dtype=h5py.string_dtype())
        f.create_dataset("sample_names", data=sample_names, dtype=h5py.string_dtype())

    print(f"Generated toy data: {n_genes} genes x {n_samples} samples")
    print(f"  Expected edges: (1,2)=disappear, (3,4)=new, (5,6)=sign_change, (7,8)=unchanged")
    print(f"  Saved to {h5_path}")

    # Also write TSV for inspection
    tsv_path = h5_path.with_suffix(".tsv")
    import pandas as pd
    df = pd.DataFrame(expr, index=gene_names, columns=sample_names)
    df.index.name = "gene_id"
    df.to_csv(tsv_path, sep="\t")
    print(f"  Also saved TSV: {tsv_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert expression TSV to HDF5 for fast parallel access."
    )
    parser.add_argument(
        "--expr-tsv",
        type=str,
        default=None,
        help="Input expression TSV (genes as rows, samples as columns).",
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
    parser.add_argument(
        "--toy",
        action="store_true",
        help="Generate toy data with known differential edges (10 genes, 60 samples).",
    )

    args = parser.parse_args()
    h5_path = Path(args.out_h5)

    if args.toy:
        generate_toy_data(h5_path)
    elif args.expr_tsv:
        tsv_path = Path(args.expr_tsv)
        if not tsv_path.exists():
            raise FileNotFoundError(f"TSV not found: {tsv_path}")
        h5_path.parent.mkdir(parents=True, exist_ok=True)
        convert_tsv_to_hdf5(
            tsv_path=tsv_path,
            h5_path=h5_path,
            compression=args.compression,
            compression_opts=args.compression_level if args.compression == "gzip" else None,
        )
    else:
        raise SystemExit("ERROR: provide --expr-tsv or use --toy.")


if __name__ == "__main__":
    main()
