Here is a summary of every change and the terminal examples.

# Changes made

## 01get_extreme_pop_bootstrap.py

Added `import argparse`, a `parse_args()`, and a `main()` that replaces the old `_toy()`. The argument set mirrors the fields of GradientParams exactly — nothing extra, nothing missing:

|flag|	maps to|	default|
|---|---|---|
|in-tsv|	input expression file|	(none)|
|out-h5|	output indices file|	required|
|low-frac|	GradientParams.low_frac|	0.2|
|high-frac|	GradientParams.high_frac|	0.2|
|n-bootstrap|	GradientParams.n_bootstrap|	50|
|bootstrap-frac|	GradientParams.bootstrap_frac|	0.8|
|seed|	GradientParams.seed|	0|
|toy|	use 5×50 random matrix|	—|

Key Features Implemented

|Script|	Purpose|
|---|---|
|02a_calc_base_correlations.py|	Compute base correlations + 3 significance tests|
|02b_bootstrap_significant_edges.py|	Bootstrap ONLY significant edges|
|03_reconstruct_diff_network.py|	Final filtering + network reconstruction|

Stage 2a outputs significance masks:

- sig_low          = low_qval < 0.05
- sig_high         = high_qval < 0.05
- sig_individual   = sig_low | sig_high           # at least one significant
- sig_differential = diff_qval < 0.05             # significant difference
- sig_edges        = sig_individual & sig_differential  # both criteria

Stage 2b computes:

Bootstrap delta for each significant edge
Confidence intervals (CI)
Bootstrap bias: bias = delta_boot_mean - delta_base
Stage 3 supports:

- `--edge-selection` sig_edges (default) or `--edge-selection` sig_differential
- Effect size filter: --min-effect 0.1
- CI filter: require CI excludes 0 (default, disable with `--no-ci-filter`)
- Sparse matrix reconstruction for network visualization
- Optional TSV output

Efficiency Gain

|Metric	Old Pipeline|	New Pipeline|
|---|---|
|Bootstrap computations|	40M × 50 × 2|	~400K × 50 × 2|
|Speedup|	1×|	~100×|

## 01calc_corr_edge.py

Added `--skip-pval`. When set, three things are removed from each job:

1. `corr_to_pvals_from_t` — the `scipy.stats.t.sf call` that is the single slowest operation in the script (not efficiently vectorised internally)
2. The entire `multipletests` FDR pass over N×(N-1)/2 p-values
3. The `pval_triu`, `qval_triu`, `reject_triu` datasets — the output HDF5 contains only corr_triu plus the correlation summary stats

The output difference is visible in the test above: full mode writes 10 datasets, `--skip-pval` writes 6. Black also reformatted all the pre-existing line-length overflows in this file.

## Terminal examples

### Stage 1 — generate bootstrap indices (single job, runs once):

```bash
# From real data
python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
    --in-tsv dataset/voomdataCtrl.txt \
    --out-h5 results/bootstrap_indices.h5 \
    --low-frac 0.2 --high-frac 0.2 \
    --n-bootstrap 50 --bootstrap-frac 0.8 \
    --seed 42

# Quick smoke test (no input file needed)
python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
    --toy --out-h5 results/test_indices.h5
Stage 2 — per-replicate correlations (one SGE array task per gene, --skip-pval on):
```

### Stage 2 — per-replicate correlations (one SGE array task per gene, --skip-pval on):

```bash
# Single replicate, manually (what each SGE task will do)
python src/scripts/10spearman_corr/01calc_corr_edge.py \
    --in-tsv subset_gene0_boot3_low.tsv \
    --out-h5 results/corr/gene0_boot3_low.h5 \
    --skip-pval \
    --block-size 1024

# Full pipeline toy test (pval included, for validation)
python src/scripts/10spearman_corr/01calc_corr_edge.py \
    --toy --out-h5 results/test_edge.h5 --validate-toy

# Same toy, but in bootstrap-replicate mode (no pval, faster)
python src/scripts/10spearman_corr/01calc_corr_edge.py \
    --toy --out-h5 results/test_edge_nopval.h5 --skip-pval
SGE job array template (Stage 2 depends on Stage 1):
```

## SGE job array template (Stage 2 depends on Stage 1):

### Stage 1 — single job

`-hold_jid 0`:

1. `-hold_jid` means: "Wait until the job(s) with these IDs finish before starting."
2. `0` means: "Wait for no job — i.e., no dependency."
3. is equivalent to no dependency at all.

wrap="...":

1. This tells SGE to run the command inside a shell (like bash), and the `wrap` argument is the entire command string to execute.
2. The `wrap` is useful when your command is long or contains special characters (like `--`, `>`, `|`, etc.).

```bash
qsub -N bootstrap_indices -l mem=8G \
    -hold_jid 0 \  # no dependency for stage 1
    wrap="python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
          --in-tsv dataset/voomdataCtrl.txt \
          --out-h5 results/bootstrap_indices.h5 \
          --n-bootstrap 50 --seed 42"
```


### Stage 2 — job array, one task per gene

```bash
# Each task reads its gene's indices from the shared HDF5,
# subsets the expression, writes its own output HDF5.
qsub -N corr_array -t 1-<N_GENES> \
    -hold_jid bootstrap_indices \
    -l mem=16G \
    -b y "python src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py \
          --gene-id \$((SGE_TASK_ID - 1)) \
          --indices-h5 results/bootstrap_indices.h5 \
          --expr-tsv dataset/voomdataCtrl.txt \
          --out-dir results/corr"
```

SGE task IDs are 1-based; `$((SGE_TASK_ID - 1))` converts back to the 0-based gene index used in the HDF5. The full two-stage submission is now handled by `src/SGE_scripts/run_bootstrap_pipeline.sh`.