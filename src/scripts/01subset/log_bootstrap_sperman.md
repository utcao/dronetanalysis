Here is a summary of every change and the terminal examples.

# Changes made

## 01subset_bootstrap.py

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


## 01cal_corr_edge.py

Added `--skip-pval`. When set, three things are removed from each job:

1. `corr_to_pvals_from_t` — the `scipy.stats.t.sf call` that is the single slowest operation in the script (not efficiently vectorised internally)
2. The entire `multipletests` FDR pass over N×(N-1)/2 p-values
3. The `pval_triu`, `qval_triu`, `reject_triu` datasets — the output HDF5 contains only corr_triu plus the correlation summary stats

The output difference is visible in the test above: full mode writes 10 datasets, `--skip-pval` writes 6. Black also reformatted all the pre-existing line-length overflows in this file.

## Terminal examples

### Stage 1 — generate bootstrap indices (single job, runs once):

```bash
# From real data
python src/scripts/01subset/01subset_bootstrap.py \
    --in-tsv dataset/voomdataCtrl.txt \
    --out-h5 results/bootstrap_indices.h5 \
    --low-frac 0.2 --high-frac 0.2 \
    --n-bootstrap 50 --bootstrap-frac 0.8 \
    --seed 42

# Quick smoke test (no input file needed)
python src/scripts/01subset/01subset_bootstrap.py \
    --toy --out-h5 results/test_indices.h5
Stage 2 — per-replicate correlations (one SGE array task per gene, --skip-pval on):
```

### Stage 2 — per-replicate correlations (one SGE array task per gene, --skip-pval on):

```bash
# Single replicate, manually (what each SGE task will do)
python src/scripts/10spearman_corr/01cal_corr_edge.py \
    --in-tsv subset_gene0_boot3_low.tsv \
    --out-h5 results/corr/gene0_boot3_low.h5 \
    --skip-pval \
    --block-size 1024

# Full pipeline toy test (pval included, for validation)
python src/scripts/10spearman_corr/01cal_corr_edge.py \
    --toy --out-h5 results/test_edge.h5 --validate-toy

# Same toy, but in bootstrap-replicate mode (no pval, faster)
python src/scripts/10spearman_corr/01cal_corr_edge.py \
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
    wrap="python src/scripts/01subset/01subset_bootstrap.py \
          --in-tsv dataset/voomdataCtrl.txt \
          --out-h5 results/bootstrap_indices.h5 \
          --n-bootstrap 50 --seed 42"
```


### Stage 2 — job array, one task per gene

```bash
# Each task reads its gene's indices from the shared HDF5,
# subsets the expression, writes its own output HDF5.
qsub -N corr_array -t 0-<N_GENES> \
    -hold_jid bootstrap_indices \
    -l mem=16G \
    wrap="python src/scripts/10spearman_corr/run_one_gene.py \
          --gene-id \$SGE_TASK_ID \
          --indices-h5 results/bootstrap_indices.h5 \
          --expr-tsv dataset/voomdataCtrl.txt \
          --out-dir results/corr \
          --skip-pval"
```

The `run_one_gene.py` wrapper (not yet written) would: read indices/low_boot[gene_id] and indices/high_boot[gene_id] from the shared HDF5, loop over replicates, column-subset the expression, and call compute_and_store(..., skip_pval=True) for each. That is the natural next piece once you are ready to wire the two stages together.