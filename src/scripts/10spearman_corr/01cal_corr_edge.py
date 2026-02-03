#!/usr/bin/env python3
"""
Blockwise spearman co-expression pipeline with HDF5 storage. (rank->zscore->matmul), following Black style guide.

For test:

- python src/test/test_coexpr_data.py --out-h5 results/test_toy.hdf5 --toy --validate-toy
- python src/test/test_coexpr_data.py --in-tsv dataset/test/voomdataCtrl_test.txt --out-h5 resul
ts/test_voomct.hdf5

Computes all pairwise Spearman correlations for a gene-expression matrix (genes x samples) and stores results in HDF5.

The core trick is the identity Spearman(x,y) == Pearson(rank(x), rank(y)), which lets it convert the problem to:

1. Rank each gene row, then z-score the ranks (ddof=1)
2. Matrix-multiply to get the full correlation: C = Z @ Z.T / (M-1)
3. Compute p-values via the t-distribution, then apply BH-FDR correction
4. Store only the upper triangle as a flat 1D array in HDF5 (row-major order), avoiding the cost of an N×N dense matrix

The blockwise loop tiles the Z @ Z.T multiply so that large gene counts stay within memory.

Storage:
- Store only the upper triangle (excluding diagonal) as a 1D vector of length K=N*(N-1)/2
  in row-major-by-left-gene order:
    (0,1)(0,2)...(0,N-1)(1,2)...(N-2,N-1)

for index:

1. gene index: gene_i, gene_j are row/column numbers (0..N-1).
2. pair index: k is the position in the 1D arrays (corr_triu[k], pval_triu[k]).

In the loop, we always:

1. choose a gene row (gene_i)
2. choose a range of partner genes (partner_start_gene .. partner_end_gene-1)
3. map that gene range to a contiguous k slice [k_start:k_end]

"""

from __future__ import annotations  # Tuple: just use  tuple[…]  inline

import argparse
import os
import time
from dataclasses import dataclass
from collections.abc import Sequence  # replaces Iterable (see below)

import h5py
import numpy as np
from scipy.stats import rankdata, t
from statsmodels.stats.multitest import multipletests


@dataclass
class Progress:
    """Simple progress tracker with ETA for blockwise loops."""

    total_steps: int
    label: str = "progress"
    report_every: int = 1

    def __post_init__(self) -> None:
        self.start_time = time.time()
        self.last_report_step = 0

    def update(self, step: int) -> None:
        """Report progress at the given step (1-based or 0-based is fine if consistent)."""
        if (
            step - self.last_report_step < self.report_every
            and step != self.total_steps
        ):
            return

        elapsed = time.time() - self.start_time
        frac = step / self.total_steps if self.total_steps else 1.0
        frac = min(max(frac, 0.0), 1.0)

        eta = (elapsed * (1.0 - frac) / frac) if frac > 0 else float("inf")
        rate = step / elapsed if elapsed > 0 else 0.0

        print(
            f"[{self.label}] {step}/{self.total_steps} "
            f"({frac * 100:.1f}%) | elapsed {elapsed:.1f}s | ETA {eta:.1f}s | {rate:.2f} steps/s"
        )
        self.last_report_step = step


@dataclass(frozen=True)
class TriuIndex:
    """Index mapping for upper-triangle (excluding diagonal) stored as 1D.

    Order:
    (0,1)(0,2)...(0,N-1)(1,2)...(N-2,N-1)

    Mapping (i < j):
        row_start(i) = i*(N-1) - i*(i-1)//2
        k(i,j) = row_start(i) + (j - i - 1)
    """

    n_genes: int

    def k(self, i: int, j: int) -> int:
        """Map gene pair (i, j) to 1D index k. Requires i < j."""
        if not (0 <= i < j < self.n_genes):
            raise ValueError(
                f"Require 0 <= i < j < N. Got i={i}, j={j}, N={self.n_genes}"
            )
        row_start = i * (self.n_genes - 1) - (i * (i - 1)) // 2
        return row_start + (j - i - 1)

    def row_start(self, i: int) -> int:
        """Start index in triu vector for row i (pairs (i, i+1..N-1))."""
        if not (0 <= i < self.n_genes):
            raise ValueError(f"i out of range: {i}")
        return i * (self.n_genes - 1) - (i * (i - 1)) // 2

    @property
    def n_tests(self) -> int:
        """Total number of unique pairs (upper triangle excluding diagonal)."""
        n = self.n_genes
        return (n * (n - 1)) // 2


def rank_zscore_rows(x: np.ndarray, ddof: int = 1) -> np.ndarray:
    """Rank-transform each row and z-score it (ddof=1 is recommended).

    Args:
        x: Array shape (N_genes, M_samples).
        ddof: Degrees of freedom for std. Use 1 to align with correlation scaling.

    Returns:
        Z: float32 array shape (N_genes, M_samples), row-wise ranked and z-scored.
    """
    # rankdata supports axis; method="average" handles ties.
    r = rankdata(x, axis=1, method="average").astype(np.float32)

    # Row-wise mean and std
    mean = r.mean(axis=1, keepdims=True)
    r_centered = r - mean
    # Sample std with ddof=1: std^2 = sum((x-mean)^2)/(n-1)
    # Equivalent, shorter, and uses NumPy's internal optimised path:
    denom = r.std(axis=1, keepdims=True, ddof=ddof)

    # Avoid division by zero if a row has zero variance.
    # A constant-variance row can produce a value very close to—but not exactly—zero
    # due to floating-point arithmetic. Use a threshold:
    denom = np.where(denom < np.finfo(np.float32).eps, np.nan, denom)

    z = r_centered / denom
    return z.astype(np.float32)


def corr_to_pvals_from_t(r: np.ndarray, n_samples: int) -> np.ndarray:
    """Approximate two-sided p-values for correlation values using t distribution.

    For Pearson correlation:
      t = r * sqrt((n-2)/(1-r^2)), df = n-2

    Args:
        r: correlation array (any shape).
        n_samples: number of samples (M).

    Returns:
        p-values array same shape as r (float32).
    """
    df = n_samples - 2
    r = np.clip(r, -1.0, 1.0)

    # Handle r close to +/-1 safely.
    eps = np.finfo(np.float32).eps
    denom = np.maximum(1.0 - r * r, eps)
    t_stat = r * np.sqrt(df / denom)

    p = 2.0 * t.sf(np.abs(t_stat), df=df)
    return p.astype(np.float32)


def create_h5_datasets(
    h5: h5py.File,
    n_genes: int,
    n_samples: int,
    gene_ids: Sequence[str],  # or simply list[str]
    sample_ids: Sequence[str],
    chunk_len_1d: int = 1_000_000,
    store_z: bool = False,
    skip_pval: bool = False,
) -> None:
    """Create HDF5 groups/datasets with attributes."""
    idx = TriuIndex(n_genes=n_genes)

    h5.attrs["tool"] = "spearman_h5"
    h5.attrs["method"] = "spearman"
    h5.attrs["triangle"] = "upper_no_diag"
    h5.attrs["indexing"] = "row_major_by_left_gene"
    h5.attrs["n_genes"] = n_genes
    h5.attrs["n_samples"] = n_samples
    h5.attrs["n_tests"] = idx.n_tests
    h5.attrs["k_formula"] = (
        "row_start(i)=i*(N-1)-i*(i-1)//2; k(i,j)=row_start(i)+(j-i-1), for 0<=i<j<N"
    )

    meta = h5.require_group("meta")
    meta.create_dataset(
        "gene_ids", data=np.array(list(gene_ids), dtype=h5py.string_dtype())
    )
    meta.create_dataset(
        "sample_ids", data=np.array(list(sample_ids), dtype=h5py.string_dtype())
    )

    results = h5.require_group("results")
    results.create_dataset(
        "corr_triu",
        shape=(idx.n_tests,),
        dtype=np.float32,
        chunks=(min(chunk_len_1d, idx.n_tests),),
        compression="gzip",
        shuffle=True,
    )
    if not skip_pval:
        results.create_dataset(
            "pval_triu",
            shape=(idx.n_tests,),
            dtype=np.float32,
            chunks=(min(chunk_len_1d, idx.n_tests),),
            compression="gzip",
            shuffle=True,
        )

    if store_z:
        dev = h5.require_group("dev")
        # Chunk by gene blocks; keep full sample axis for contiguous math/debug.
        dev.create_dataset(
            "Z_rank_zscore",
            shape=(n_genes, n_samples),
            dtype=np.float32,
            chunks=(min(256, n_genes), n_samples),
            compression="gzip",
            shuffle=True,
        )
        dev.attrs["rank_method"] = "average"
        dev.attrs["zscore_ddof"] = 1


def write_block_to_triu_1d(
    triu_index: TriuIndex,
    corr_triu: h5py.Dataset,
    pval_triu: h5py.Dataset,
    c_block: np.ndarray,
    p_block: np.ndarray,
    row_gene_start: int,
    col_gene_start: int,
) -> None:
    """Write a correlation/pvalue block into 1D upper-triangle datasets.

    Args:
        triu_index: TriuIndex mapping for N genes.
        corr_triu: HDF5 1D dataset for correlations.
        pval_triu: HDF5 1D dataset for p-values.
        c_block: (Br, Bc) correlation block for rows row_gene_start.. and cols col_gene_start..
        p_block: (Br, Bc) p-values block aligned with c_block.
        row_gene_start: global gene index of c_block row 0.
        col_gene_start: global gene index of c_block col 0.
    """
    br, bc = c_block.shape

    for local_r in range(br):
        gene_i = row_gene_start + local_r

        partner_start_gene = max(col_gene_start, gene_i + 1)
        partner_end_gene = col_gene_start + bc  # exclusive

        if partner_start_gene >= partner_end_gene:
            continue

        local_c_start = partner_start_gene - col_gene_start
        local_c_end = partner_end_gene - col_gene_start

        k_start = triu_index.k(gene_i, partner_start_gene)
        k_end = triu_index.k(gene_i, partner_end_gene - 1) + 1

        corr_triu[k_start:k_end] = c_block[local_r, local_c_start:local_c_end]
        if pval_triu is not None:
            pval_triu[k_start:k_end] = p_block[local_r, local_c_start:local_c_end]


def compute_and_store(
    x: np.ndarray,
    gene_ids: Sequence[str],
    sample_ids: Sequence[str],
    out_h5: str,
    block_size: int,
    store_z: bool,
    fdr_alpha: float,
    skip_pval: bool = False,
) -> None:
    """Compute Spearman correlation/p-values in blocks and store to HDF5."""
    n_genes, n_samples = x.shape
    triu_index = TriuIndex(n_genes=n_genes)

    z = rank_zscore_rows(x, ddof=1)

    os.makedirs(os.path.dirname(out_h5) or ".", exist_ok=True)
    with h5py.File(out_h5, "w") as h5:
        create_h5_datasets(
            h5,
            n_genes=n_genes,
            n_samples=n_samples,
            gene_ids=gene_ids,
            sample_ids=sample_ids,
            store_z=store_z,
            skip_pval=skip_pval,
        )

        if store_z:
            h5["dev/Z_rank_zscore"][:, :] = z

        corr_triu = h5["results/corr_triu"]
        pval_triu = h5["results/pval_triu"] if not skip_pval else None

        # display progress
        n_row_blocks = (n_genes + block_size - 1) // block_size
        progress = Progress(
            total_steps=n_row_blocks, label="corr+pval blocks", report_every=1
        )
        # Blockwise upper triangle by blocks
        for row_block_idx, row_start in enumerate(
            range(0, n_genes, block_size), start=1
        ):
            progress.update(row_block_idx)

            row_end = min(n_genes, row_start + block_size)
            z_rows = z[row_start:row_end, :]

            for col_start in range(row_start, n_genes, block_size):
                col_end = min(n_genes, col_start + block_size)
                z_cols = z[col_start:col_end, :]

                # Correlation block
                c_block = (z_rows @ z_cols.T) / float(n_samples - 1)
                c_block = c_block.astype(np.float32)

                # P-values block (skip when running bootstrap replicates)
                p_block = None
                if not skip_pval:
                    p_block = corr_to_pvals_from_t(c_block, n_samples=n_samples)

                write_block_to_triu_1d(
                    triu_index=triu_index,
                    corr_triu=corr_triu,
                    pval_triu=pval_triu,
                    c_block=c_block,
                    p_block=p_block,
                    row_gene_start=row_start,
                    col_gene_start=col_start,
                )

        # BH-FDR on full p-value vector (skipped in bootstrap-replicate mode)
        if not skip_pval:
            pvals = pval_triu[:].astype(np.float64)
            reject, qvals, _, _ = multipletests(pvals, alpha=fdr_alpha, method="fdr_bh")

            h5["results"].create_dataset(
                "qval_triu",
                data=qvals.astype(np.float32),
                chunks=(min(1_000_000, triu_index.n_tests),),
                compression="gzip",
                shuffle=True,
            )
            h5["results"].create_dataset(
                "reject_triu",
                data=reject.astype(np.uint8),
                chunks=(min(1_000_000, triu_index.n_tests),),
                compression="gzip",
                shuffle=True,
            )

        # Summary
        summary = h5.require_group("summary")
        corr_arr = corr_triu[:]  # single disk read
        summary.create_dataset("corr_min", data=np.float32(np.nanmin(corr_arr)))
        summary.create_dataset("corr_max", data=np.float32(np.nanmax(corr_arr)))
        summary.create_dataset(
            "mean_abs_corr", data=np.float32(np.nanmean(np.abs(corr_arr)))
        )
        if not skip_pval:
            summary.attrs["fdr_alpha"] = fdr_alpha
            summary.create_dataset("n_reject", data=np.int64(reject.sum()))


def triu_to_full_matrix(
    triu_1d: np.ndarray, n_genes: int, diag_value: float = 1.0
) -> np.ndarray:
    """Reconstruct full symmetric matrix from upper-triangle 1D storage."""
    idx = TriuIndex(n_genes=n_genes)
    mat = np.full((n_genes, n_genes), np.nan, dtype=np.float32)
    np.fill_diagonal(mat, np.float32(diag_value))

    for i in range(n_genes - 1):
        k0 = idx.k(i, i + 1)
        k1 = idx.k(i, n_genes - 1) + 1
        row_seg = triu_1d[k0:k1]
        mat[i, i + 1 :] = row_seg
        mat[i + 1 :, i] = row_seg

    return mat


def load_tsv_expression(path: str) -> tuple[np.ndarray, list[str], list[str]]:
    """Load TSV with genes as rows and samples as columns; first column is gene id."""
    import pandas as pd  # local import; optional dependency

    df = pd.read_csv(path, sep="\t", index_col=0)
    x = df.to_numpy(dtype=np.float32)
    gene_ids = df.index.astype(str).tolist()
    sample_ids = df.columns.astype(str).tolist()
    return x, gene_ids, sample_ids


def make_toy_data(
    n_genes: int = 10, n_samples: int = 10, seed: int = 0
) -> tuple[np.ndarray, list[str], list[str]]:
    """Generate toy expression data with some ties."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n_genes, n_samples)).astype(np.float32)

    # Introduce ties for realism (Spearman tie handling)
    if n_genes >= 3 and n_samples >= 5:
        x[0, :3] = 1.0
        x[1, 2:5] = -2.0

    gene_ids = [f"gene_{i}" for i in range(n_genes)]
    sample_ids = [f"sample_{j}" for j in range(n_samples)]
    return x, gene_ids, sample_ids


def validate_toy(out_h5: str, x: np.ndarray) -> None:
    """Validate toy output against SciPy spearmanr full matrix."""
    from scipy.stats import spearmanr

    n_genes, _ = x.shape

    with h5py.File(out_h5, "r") as h5:
        corr_triu = h5["results/corr_triu"][:]
        corr_full = triu_to_full_matrix(corr_triu, n_genes=n_genes, diag_value=1.0)

    # SciPy spearmanr: compute correlation of rows (genes)
    # spearmanr returns correlation matrix when given a 2D array; set axis=1 to treat rows as variables.
    sci_corr, _ = spearmanr(x, axis=1)

    # sci_corr is (n_genes, n_genes)
    diff = np.nanmax(np.abs(corr_full - sci_corr.astype(np.float32)))
    print(f"[validate] max|difference| vs scipy spearmanr: {diff:.6g}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Spearman coexpression to HDF5 (upper triangle 1D)."
    )
    parser.add_argument(
        "--in-tsv",
        type=str,
        default=None,
        help="Input TSV (genes as rows, samples as columns).",
    )
    parser.add_argument("--out-h5", type=str, required=True, help="Output HDF5 path.")
    parser.add_argument(
        "--block-size",
        type=int,
        default=1024,
        help="Gene block size for matmul (default: 1024).",
    )
    parser.add_argument(
        "--store-z",
        action="store_true",
        help="Store /dev/Z_rank_zscore (development/debug only).",
    )
    parser.add_argument(
        "--fdr-alpha", type=float, default=0.05, help="BH FDR alpha (default: 0.05)."
    )
    parser.add_argument(
        "--skip-pval",
        action="store_true",
        help="Skip p-value and FDR (use for bootstrap replicates).",
    )

    parser.add_argument(
        "--toy",
        action="store_true",
        help="Run on toy data (10x10) instead of input file.",
    )
    parser.add_argument("--toy-seed", type=int, default=0, help="Seed for toy data.")
    parser.add_argument(
        "--validate-toy",
        action="store_true",
        help="After running toy, compare vs scipy.spearmanr().",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.toy:
        x, gene_ids, sample_ids = make_toy_data(
            n_genes=10, n_samples=10, seed=args.toy_seed
        )
    else:
        if args.in_tsv is None:
            raise SystemExit("Provide --in-tsv or use --toy.")
        x, gene_ids, sample_ids = load_tsv_expression(args.in_tsv)

    compute_and_store(
        x=x,
        gene_ids=gene_ids,
        sample_ids=sample_ids,
        out_h5=args.out_h5,
        block_size=args.block_size,
        store_z=args.store_z,
        fdr_alpha=args.fdr_alpha,
        skip_pval=args.skip_pval,
    )

    print(f"Saved: {args.out_h5}")

    if args.toy and args.validate_toy:
        validate_toy(args.out_h5, x)


if __name__ == "__main__":
    main()
