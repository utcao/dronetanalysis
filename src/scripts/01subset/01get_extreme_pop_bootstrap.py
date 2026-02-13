#!/usr/bin/env python3
"""
Bootstrap index generation for low/high expression subpopulations.

For every gene in the expression matrix:
1. Partitions samples into the lowest low_frac and highest high_frac expression
   groups using argpartition (O(n), no full sort).
2. Draws n_bootstrap resamples (with replacement, bootstrap_frac size) from
   each group — fully vectorised across all genes.

Output HDF5 layout (all datasets int32):
    indices/low            (n_genes, k_low)
    indices/high           (n_genes, k_high)
    indices/low_boot       (n_genes, n_bootstrap, k_resample_low)
    indices/high_boot      (n_genes, n_bootstrap, k_resample_high)

Each dataset is chunked along the gene axis (chunk size 1) so that a single
SGE job can read one gene's indices with one HDF5 chunk read.

**UPDATED**: Can now read from HDF5 expression files for faster loading.

"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


@dataclass(frozen=True)
class GradientParams:
    low_frac: float  # e.g. 0.2 or 0.3
    high_frac: float  # usually same as low_frac
    n_bootstrap: int  # e.g. 50
    bootstrap_frac: float  # e.g. 0.8
    seed: int = 0


def _estimate_memory_gb(n_genes: int, n_bootstrap: int, k_resample: int) -> float:
    """
    Estimate peak memory usage in GB for bootstrap index generation.

    Accounts for:
    - rand_low, rand_high: 2 × (n_genes, n_bootstrap, k_resample) × 4 bytes (int32)
    - all_low_boot, all_high_boot: 2 × (n_genes, n_bootstrap, k_resample) × 4 bytes
    - Plus overhead for base indices and expression matrix
    """
    # Random indices (int32)
    rand_arrays = 2 * n_genes * n_bootstrap * k_resample * 4
    # Bootstrap indices (int32)
    boot_arrays = 2 * n_genes * n_bootstrap * k_resample * 4
    # Total in GB with 20% overhead
    total_bytes = (rand_arrays + boot_arrays) * 1.2
    return total_bytes / (1024**3)


def build_gradient_indices_h5(
    expr: np.ndarray,  # shape (n_genes, n_samples)
    out_h5: Path,
    params: GradientParams,
    batch_size: int = None,
) -> None:
    """
    Build low/high and bootstrap sample indices for every gene and save into one HDF5.

    Parameters
    ----------
    expr : (n_genes, n_samples) expression matrix
    out_h5 : output HDF5 path
    params : bootstrap parameters
    batch_size : genes per batch (default: None = vectorized, all at once)
                 Use batch processing for large datasets or limited memory.
                 Recommended: 100-500 genes per batch.

    Output datasets (all int32, chunked along gene axis):
        indices/low            (n_genes, k_low)
        indices/high           (n_genes, k_high)
        indices/low_boot       (n_genes, n_bootstrap, k_resample_low)
        indices/high_boot      (n_genes, n_bootstrap, k_resample_high)
    """
    if expr.ndim != 2:
        raise ValueError("expr must be 2D (n_genes, n_samples)")

    n_genes, n_samples = expr.shape
    k_low = int(np.floor(n_samples * params.low_frac))
    k_high = int(np.floor(n_samples * params.high_frac))
    k_resample_low = int(np.floor(k_low * params.bootstrap_frac))
    k_resample_high = int(np.floor(k_high * params.bootstrap_frac))

    if k_low <= 0 or k_high <= 0:
        raise ValueError("low/high fraction too small for number of samples")
    if k_resample_low <= 0 or k_resample_high <= 0:
        raise ValueError("bootstrap_frac too small for base index set")

    out_h5.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(params.seed)

    # Determine processing mode
    if batch_size is None:
        # --- VECTORIZED MODE: Process all genes at once (faster, more memory) ---
        print(f"Vectorized mode: processing all {n_genes} genes at once...")
        print(f"  Estimated peak memory: ~{_estimate_memory_gb(n_genes, params.n_bootstrap, k_resample_low):.1f} GB")

        # vectorised partition: O(n_samples) per gene, all genes at once
        all_low = np.argpartition(expr, k_low, axis=1)[:, :k_low].astype(np.int32)
        all_high = np.argpartition(expr, n_samples - k_high, axis=1)[:, -k_high:].astype(np.int32)

        # vectorised bootstrap: random positions within each partition → gather
        gene_idx = np.arange(n_genes)[:, None, None]  # (n_genes, 1, 1)

        # picks random positions within each gene's partition
        rand_low = rng.integers(0, k_low, size=(n_genes, params.n_bootstrap, k_resample_low), dtype=np.int32)
        all_low_boot = all_low[gene_idx, rand_low]

        rand_high = rng.integers(0, k_high, size=(n_genes, params.n_bootstrap, k_resample_high), dtype=np.int32)
        all_high_boot = all_high[gene_idx, rand_high]

        # write: one dataset per category, chunked by gene for SGE access
        with h5py.File(out_h5, "w") as h5:
            meta = h5.create_group("meta")
            meta.attrs["low_frac"] = params.low_frac
            meta.attrs["high_frac"] = params.high_frac
            meta.attrs["n_bootstrap"] = params.n_bootstrap
            meta.attrs["bootstrap_frac"] = params.bootstrap_frac
            meta.attrs["seed"] = params.seed
            meta.attrs["processing_mode"] = "vectorized"

            indices = h5.create_group("indices")
            for name, arr in [
                ("low", all_low),
                ("high", all_high),
                ("low_boot", all_low_boot),
                ("high_boot", all_high_boot),
            ]:
                indices.create_dataset(
                    name,
                    data=arr,
                    chunks=(1,) + arr.shape[1:],
                    compression="gzip",
                    compression_opts=4,
                    shuffle=True,
                )

        print("Done!")

    else:
        # --- BATCH MODE: Process genes in batches (memory-efficient) ---
        n_batches = (n_genes + batch_size - 1) // batch_size
        print(f"Batch mode: processing {n_genes} genes in {n_batches} batches of {batch_size}...")
        print(f"  Estimated peak memory: ~{_estimate_memory_gb(batch_size, params.n_bootstrap, k_resample_low):.1f} GB per batch")

        # Create HDF5 file and datasets first
        with h5py.File(out_h5, "w") as h5:
            meta = h5.create_group("meta")
            meta.attrs["low_frac"] = params.low_frac
            meta.attrs["high_frac"] = params.high_frac
            meta.attrs["n_bootstrap"] = params.n_bootstrap
            meta.attrs["bootstrap_frac"] = params.bootstrap_frac
            meta.attrs["seed"] = params.seed
            meta.attrs["processing_mode"] = f"batch_{batch_size}"

            indices = h5.create_group("indices")

            # Create empty datasets with final shape
            ds_low = indices.create_dataset(
                "low",
                shape=(n_genes, k_low),
                dtype=np.int32,
                chunks=(1, k_low),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            ds_high = indices.create_dataset(
                "high",
                shape=(n_genes, k_high),
                dtype=np.int32,
                chunks=(1, k_high),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            ds_low_boot = indices.create_dataset(
                "low_boot",
                shape=(n_genes, params.n_bootstrap, k_resample_low),
                dtype=np.int32,
                chunks=(1, params.n_bootstrap, k_resample_low),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            ds_high_boot = indices.create_dataset(
                "high_boot",
                shape=(n_genes, params.n_bootstrap, k_resample_high),
                dtype=np.int32,
                chunks=(1, params.n_bootstrap, k_resample_high),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )

            # Process genes in batches
            for batch_idx in range(n_batches):
                start_gene = batch_idx * batch_size
                end_gene = min(start_gene + batch_size, n_genes)
                batch_n_genes = end_gene - start_gene

                if batch_idx % max(1, n_batches // 10) == 0 or batch_idx == n_batches - 1:
                    print(f"  Batch {batch_idx+1}/{n_batches}: genes {start_gene}-{end_gene-1}")

                # Get expression for this batch
                expr_batch = expr[start_gene:end_gene, :]

                # Vectorised partition for this batch
                batch_low = np.argpartition(expr_batch, k_low, axis=1)[:, :k_low].astype(np.int32)
                batch_high = np.argpartition(expr_batch, n_samples - k_high, axis=1)[:, -k_high:].astype(np.int32)

                # Write base indices
                ds_low[start_gene:end_gene, :] = batch_low
                ds_high[start_gene:end_gene, :] = batch_high

                # Generate bootstrap indices for this batch
                gene_idx_batch = np.arange(batch_n_genes)[:, None, None]

                rand_low = rng.integers(0, k_low, size=(batch_n_genes, params.n_bootstrap, k_resample_low), dtype=np.int32)
                batch_low_boot = batch_low[gene_idx_batch, rand_low]
                ds_low_boot[start_gene:end_gene, :, :] = batch_low_boot

                rand_high = rng.integers(0, k_high, size=(batch_n_genes, params.n_bootstrap, k_resample_high), dtype=np.int32)
                batch_high_boot = batch_high[gene_idx_batch, rand_high]
                ds_high_boot[start_gene:end_gene, :, :] = batch_high_boot

        print("Done!")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate low/high bootstrap sample indices -> HDF5."
    )
    parser.add_argument(
        "--in-tsv",
        type=str,
        default=None,
        help="Input TSV (genes as rows, samples as columns).",
    )
    parser.add_argument(
        "--in-h5",
        type=str,
        default=None,
        help="Input HDF5 expression file (alternative to --in-tsv). Use /expr dataset.",
    )
    parser.add_argument(
        "--out-h5", type=str, required=True, help="Output HDF5 path for indices."
    )
    parser.add_argument(
        "--low-frac",
        type=float,
        default=0.2,
        help="Fraction of samples in low-expression group (default: 0.2).",
    )
    parser.add_argument(
        "--high-frac",
        type=float,
        default=0.2,
        help="Fraction of samples in high-expression group (default: 0.2).",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=50,
        help="Number of bootstrap replicates (default: 50).",
    )
    parser.add_argument(
        "--bootstrap-frac",
        type=float,
        default=0.8,
        help="Fraction of partition kept per bootstrap resample (default: 0.8).",
    )
    parser.add_argument("--seed", type=int, default=0, help="RNG seed (default: 0).")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="Process genes in batches (for limited memory). Default: None (vectorized mode, faster but more memory). Recommended: 100-500 for large datasets.",
    )
    parser.add_argument(
        "--toy",
        action="store_true",
        help="Run on toy data (5 genes, 50 samples) instead of input file.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.toy:
        expr = np.random.default_rng(1).normal(size=(5, 50)).astype(np.float32)
    else:
        # *** UPDATED: Support both TSV and HDF5 input ***
        if args.in_h5:
            print(f"Loading expression from HDF5: {args.in_h5}")
            with h5py.File(args.in_h5, "r") as f:
                expr = f["expr"][:]
            print(f"  Loaded {expr.shape[0]} genes × {expr.shape[1]} samples")
        elif args.in_tsv:
            print(f"Loading expression from TSV: {args.in_tsv}")
            import pandas as pd
            df = pd.read_csv(args.in_tsv, sep="\t", index_col=0)
            expr = df.to_numpy(dtype=np.float32)
            print(f"  Loaded {expr.shape[0]} genes × {expr.shape[1]} samples")
        else:
            raise SystemExit("Provide --in-tsv, --in-h5, or use --toy.")

    params = GradientParams(
        low_frac=args.low_frac,
        high_frac=args.high_frac,
        n_bootstrap=args.n_bootstrap,
        bootstrap_frac=args.bootstrap_frac,
        seed=args.seed,
    )
    build_gradient_indices_h5(expr, Path(args.out_h5), params, batch_size=args.batch_size)
    print(f"Saved: {args.out_h5}")


if __name__ == "__main__":
    main()
