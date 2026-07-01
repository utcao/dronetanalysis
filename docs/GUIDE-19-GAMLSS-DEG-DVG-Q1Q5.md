# GAMLSS DEG/DVG Q1 vs Q5 Analysis Guide

## Overview

This guide explains how to identify differentially expressed genes (DEG) and differentially
variable genes (DVG) by stratifying samples based on the expression level of a **focal gene**.

For a given focal gene, samples within a condition (CT or HS) are ranked by the focal gene's
CPM expression. The bottom 20% form group **Q1** and the top 20% form group **Q5**; the
middle 60% are discarded. GAMLSS models with a negative binomial distribution (NBI family)
then test whether each of the other ~8,786 genes in the count matrix differs in mean
expression (DEG) or expression variability (DVG) between Q1 and Q5.

This approach asks: *what is co-regulated with, or buffered against, the expression of the
focal gene?*

---

## Prerequisites

### Software environments

| Environment | Purpose | Path |
|---|---|---|
| `base` (miniforge3) | Snakemake 9 + cluster-generic executor plugin | `/tmp/global2/caoyt/miniforge3/bin/snakemake` |
| `ganlss` conda env | R + GAMLSS libraries (Rscript) | `/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript` |

Install the executor plugin in the base environment if not already present:
```bash
/tmp/global2/caoyt/miniforge3/bin/pip install snakemake-executor-plugin-cluster-generic
```

### R packages (ganlss env)

`gamlss`, `gamlss.dist`, `edgeR`, `data.table`, `BiocParallel`, `glue`, `stringr`,
`purrr`, `dplyr`, `foreach`, `futile.logger`, `argparse`

### Input files

| File | Description |
|---|---|
| `data/count/RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt` | Raw count matrix (genes × samples); row 1 = sample IDs only (no `gene_id` column header) |
| `data/count/Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt` | Sample metadata; `id` column matches count file column names |
| `data/count/TMM_Voom_sv{1..4}_10.txt` | Surrogate variable files; positional rows align to count file columns |
| `src/gamlss_q1q5/focal_genes.csv` | CSV with `gene_id` column (FlyBase IDs) + annotation columns (symbol, description, etc.) |

### Working directory

All commands below assume `code/dronetanalysis/` as the working directory.

---

## Architecture

The analysis is split across four files in `src/gamlss_q1q5/`:

```
tools_gamlss.R                    fit_one(), fit_sum() model-fitting primitives
gamlss_dvgdeg_core.R              run_gamlss_dvgdeg(): TMM norm, model fitting, LR tests, CSV output
preprocess_quintile.R             prepare_quintile_groups(): condition filter → quintile split → SV merge
run_gamlss_quintile_analysis.R    driver: parameters, argparse CLI, calls preprocess + core
```

The Snakemake workflow `compute_degdvg_gamlss_q1q5.snakemake` submits one SGE job per
focal gene × condition pair by calling the driver with named CLI flags.

---

## Step-by-Step Instructions

### Step 1: Choose focal genes

Edit `src/gamlss_q1q5/focal_genes.csv`. The required column is `gene_id`; additional
annotation columns (symbol, description, etc.) are carried for reference but not used
by the workflow:

```csv
gene_id,symbol,description
FBgn0027560,ken,ken and barbie
FBgn0039752,CG12345,CG12345
```

Adding or removing rows changes which jobs the Snakemake workflow creates; the Snakefile
itself does not need editing. The `gene_id` values become the `{gene}` wildcard.

### Step 2: Test a single gene interactively

Before a full cluster run, verify the pipeline on one gene with test mode (caps the
number of *genes fitted*, not samples):

```bash
/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript \
    src/gamlss_q1q5/run_gamlss_quintile_analysis.R \
    --focal-gene FBgn0027560 \
    --condition-val 1 \
    --n-test-genes 20
```

Expected output: a CSV in `../../results/gamlss/` with 20 rows. Check that Q1/Q5 sample
counts are non-zero (printed to stderr) and that the CSV has all expected columns.

### Step 3: Run a full single-gene analysis

```bash
/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript \
    src/gamlss_q1q5/run_gamlss_quintile_analysis.R \
    --focal-gene FBgn0027560 \
    --condition-val 1 \
    --full-run
```

CLI flags:

| Flag | Values | Default |
|---|---|---|
| `--focal-gene` | FlyBase ID | `FBgn0027560` |
| `--condition-val` | `1` = CT, `2` = HS | `1` |
| `--full-run` | (flag, no value) | off — test mode (1 gene) |
| `--n-test-genes N` | integer | `1` — only used when `--full-run` is absent |

### Step 4: Dry-run the Snakemake workflow

Verify that the correct number of jobs are planned (2 × N for N focal genes):

```bash
/tmp/global2/caoyt/miniforge3/bin/snakemake \
    -s src/gamlss_q1q5/compute_degdvg_gamlss_q1q5.snakemake \
    --profile src/gamlss_q1q5/profiles/sge \
    -n -p
```

Each listed job should show the correct `--focal-gene`, `--condition-val`, and `--full-run`
flags in the shell command printed by `-p`.

### Step 5: Submit to the cluster

**Option A — Snakemake directly (requires an interactive session for the full run duration):**

```bash
/tmp/global2/caoyt/miniforge3/bin/snakemake \
    -s src/gamlss_q1q5/compute_degdvg_gamlss_q1q5.snakemake \
    --profile src/gamlss_q1q5/profiles/sge
```

**Option B — qsub master job (recommended; survives connection drops):**

```bash
qsub src/gamlss_q1q5/submit_gamlss_q1q5.sh
```

The master job (`gamlss_master`) runs Snakemake, which submits one child job per
focal gene × condition. Master job walltime is set to 560 h — it must exceed the
longest expected child job (48 h) plus any queue wait time.

### Step 6: Monitor progress

```bash
qstat -u $USER | grep gamlss
```

Log files per child job: `../../logs/gamlss_0608_2026/{gene}_{cond}.log`

---

## Output Files

Each focal gene × condition pair produces two CSV files in `results/gamlss_0608_2026/`:

| File | Content |
|---|---|
| `{gene}_{cond}_Q1vsQ5_gamlss_dvgdeg.csv` | One row per tested gene; all mu/sigma/BCV estimates and p/q values |
| `{gene}_{cond}_Q1vsQ5_gamlss_dvgdeg_sum.csv` | Count summary: DEG / DVG / Both / NS across all tested genes |

### Column reference — main CSV

All value columns are suffixed with `_Q5` (the non-reference treatment level).

| Column | Description |
|---|---|
| `gene_id` | FlyBase ID of the tested gene |
| `cpm_count_ctrl_Q5` | Expected CPM in Q1 (reference group) |
| `cpm_count_trt_Q5` | Expected CPM in Q5 group |
| `mu_ctrl_Q5` / `mu_trt_Q5` | GAMLSS mu coefficients (log scale) |
| `logFC_mu_Q5` | log2 fold-change in mu (Q5 vs Q1) |
| `logFC_cpm_Q5` | log2 fold-change in CPM (Q5 vs Q1) |
| `BCV_ctrl_Q5` / `BCV_trt_Q5` | Biological coefficient of variation in Q1 / Q5 |
| `logFC_BCV_Q5` | log2 fold-change in BCV (Q5 vs Q1) |
| `diff_BCV_Q5` | Arithmetic difference BCV_trt − BCV_ctrl |
| `sigma_ctrl_Q5` / `sigma_trt_Q5` | GAMLSS sigma coefficients (log scale) |
| `p_mu_Q5` / `q_mu_Q5` | LR test p-value / BH-adjusted q-value for mean difference |
| `p_sigma_Q5` / `q_sigma_Q5` | LR test p-value / BH-adjusted q-value for variability difference |
| `p_mu_single_Q5` / `q_mu_single_Q5` | Single-model t-test p/q for mean (from `m_full` summary) |
| `p_sigma_single_Q5` / `q_sigma_single_Q5` | Single-model t-test p/q for variability |
| `class_single_Q5` | Classification using t-test q-values: DEG / DVG / Both / NS |
| `class_modelcomp_Q5` | Classification using LR test q-values: DEG / DVG / Both / NS (preferred) |

---

## Statistical Design

### Model family

GAMLSS with the NBI (negative binomial type I) family:

- **mu**: mean expression (log-link)
- **sigma**: overdispersion / BCV² (log-link; NBI parameterisation: σ = log(dispersion))

BCV is recovered as `sqrt(exp(sigma_coefficient))`.

### Three-model LR test design

| Model | mu formula | sigma formula | Purpose |
|---|---|---|---|
| `m_full` | `~expr_quintile + <batch covariates> + sv1-sv4` | `~expr_quintile` | Baseline full model |
| `m_sigma_int` | `~expr_quintile + <batch covariates> + sv1-sv4` | `~1` | Null for sigma — DVG test |
| `m_mu_nontrt` | `~<batch covariates> + sv1-sv4` | `~expr_quintile` | Null for mu — DEG test |

- **DEG test**: `LR(m_full, m_mu_nontrt)` — removing `expr_quintile` from mu
- **DVG test**: `LR(m_full, m_sigma_int)` — collapsing sigma to an intercept

LR statistic = `2 × (logLik(m_full) − logLik(m_reduced))`, distributed χ² with
df = difference in fitted degrees of freedom.

### Normalization

TMM normalization (edgeR) is applied to raw counts before GAMLSS fitting. Per-sample
offsets `log(lib.size × norm.factor)` are passed as a fixed `offset()` term in the
GAMLSS mu formula, adjusting for library size while preserving absolute scale.

### Contrast coding

`expr_quintile` uses **standard treatment contrasts** (`~expr_quintile`, not `~0+`):
Q1 is the reference level (intercept); the `expr_quintileQ5` coefficient captures the
Q5 effect relative to Q1. Using `~0+` would produce collinearity and GAMLSS may silently
drop a coefficient, yielding missing or nonsensical estimates.

Sum-to-zero contrasts (`contr.sum`) are applied to the categorical batch covariates
(`RNAlibBatch`, `RNAseqBatch`, `egglayBatch`, `platingBatch`, `well`) for better
identifiability. Surrogate variables `sv1`–`sv4` are continuous and are not affected
by contrast coding.

### Classification

Genes are classified at BH-adjusted q < 0.05:

| Class | q_mu | q_sigma |
|---|---|---|
| DEG | < 0.05 | ≥ 0.05 |
| DVG | ≥ 0.05 | < 0.05 |
| Both | < 0.05 | < 0.05 |
| NS | ≥ 0.05 | ≥ 0.05 |

Two independent classifications are reported per gene: `class_modelcomp` (preferred;
uses LR test q-values) and `class_single` (uses single-model t-test q-values from
the `m_full` coefficient table; useful for comparison or when LR tests fail to converge).

---

## Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| `focal_gene '...' not found in count matrix` | Gene ID typo or wrong file | Verify ID against `rownames(count_mat)` |
| `BiocParallel errors: Error in str2lang unexpected '='` | `mu_sum_covars` set to a formula string instead of a character vector | Use `c("RNAlibBatch", "RNAseqBatch", ...)` not a single string |
| `ERROR contrasts apply only to factors` | Continuous SV listed in `mu_sum_covars`, or integer-coded batch column not coerced to factor | Check `mu_sum_covars` contains only batch variable names, not `sv1`–`sv4`; the core script applies unconditional `as.factor()` to all `mu_sum_covars` |
| `WorkflowError: Resource 'runtime' could not be parsed as minutes` | Used Snakemake-reserved `runtime` resource with an `HH:MM:SS` string | Rename to `walltime` (already done; see `compute_degdvg_gamlss_q1q5.snakemake`) |
| `MissingRuleException: No rule to produce --cluster=qsub ...` | Profile uses old `cluster:` key, removed in Snakemake 9 | Switch to `executor: cluster-generic` + `cluster-generic-submit-cmd:` in the profile YAML |
| `No module named snakemake_executor_plugin_cluster_generic` | Plugin missing from base env | `pip install snakemake-executor-plugin-cluster-generic` in the base conda env |
| One or both quintile groups are empty | Highly skewed focal gene expression | Inspect CPM distribution; consider reducing `frac_low`/`frac_high` cutoffs in the driver |
| Models converge for almost no genes | Wrong formula or missing covariates | Check that all columns in `formula_mu` exist in `meta_dt` before calling `run_gamlss_dvgdeg()` |
| `Error in { : task 1 failed - "subscript out of bounds"` (CT jobs) | Coefficient index derived from one gene's model applied to another gene with a different model-matrix shape | Already fixed via per-gene `grep()` lookup; see [FIX-04-GAMLSS-Q1Q5-CT-HS-Errors.md](FIX-04-GAMLSS-Q1Q5-CT-HS-Errors.md) |
| `Unable to allocate 8 MB * 1 thread buffers ... Cannot allocate memory` in `fwrite` (HS jobs, larger n) | Three full-size GAMLSS model lists held concurrently in memory during fitting | See ranked mitigations in [OPTIMIZATION-05-GAMLSS-Q1Q5-Memory-Mitigation.md](OPTIMIZATION-05-GAMLSS-Q1Q5-Memory-Mitigation.md) |

---

## Related Reading

- [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) — Snakemake workflow setup and SGE cluster submission (dronetanalysis bootstrap pipeline)
- [OPTIMIZATION-04-Snakemake-Efficiency.md](OPTIMIZATION-04-Snakemake-Efficiency.md) — Modular Snakemake patterns, cluster-generic executor plugin, multi-file rule organization
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — Statistical methodology reference for differential testing and multiple comparisons
- [FIX-03-HPC-SGE-Pipeline.md](FIX-03-HPC-SGE-Pipeline.md) — SGE-specific fixes (NFS locking, `latency-wait`, memory allocation)
- [FIX-04-GAMLSS-Q1Q5-CT-HS-Errors.md](FIX-04-GAMLSS-Q1Q5-CT-HS-Errors.md) — Root cause analysis for the CT `subscript out of bounds` and HS `fwrite` OOM errors observed in the 2026-06-23 run
- [OPTIMIZATION-05-GAMLSS-Q1Q5-Memory-Mitigation.md](OPTIMIZATION-05-GAMLSS-Q1Q5-Memory-Mitigation.md) — Ranked memory mitigation options for the HS OOM failures

---

**Last Updated:** 2026-06-09
**Status:** ✅ Active — workflow validated with 8-gene dry run (16 jobs); FBgn0027560 CT test run confirmed
