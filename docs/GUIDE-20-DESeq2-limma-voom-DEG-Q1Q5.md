# DESeq2 and limma-voom DEG Q1 vs Q5 Analysis Guide

## Overview

This guide explains how to identify differentially expressed genes (DEG) by stratifying
samples based on the expression level of a **focal gene**, using either DESeq2 or
limma-voom as the statistical engine.

For a given focal gene, samples within a condition (CT or HS) are ranked by the focal
gene's CPM expression. The bottom 20% form group **Q1** and the top 20% form group
**Q5**; the middle 60% are discarded. The test asks: *what is co-expressed with, or
buffered against, the focal gene?*

DESeq2 uses a negative-binomial GLM with empirical Bayes dispersion shrinkage.
limma-voom converts counts to precision-weighted log-CPM and fits a linear model.
Both methods test for mean expression differences only (DEG); neither detects DVG.
Use the GAMLSS workflow (GUIDE-19) when DVG detection is needed.

---

## Prerequisites

### Software environments

| Environment | Purpose | Path |
|---|---|---|
| `base` (miniforge3) | Snakemake 9 + cluster-generic executor plugin | `/tmp/global2/caoyt/miniforge3/bin/snakemake` |
| `ganlss` conda env | R + DESeq2/limma packages (Rscript) | `/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript` |

Install the executor plugin in the base environment if not already present:
```bash
/tmp/global2/caoyt/miniforge3/bin/pip install snakemake-executor-plugin-cluster-generic
```

### R packages (ganlss env)

**DESeq2 workflow:** `DESeq2`, `edgeR`, `BiocParallel`, `data.table`, `dplyr`, `glue`,
`stringr`, `purrr`, `foreach`, `argparse`

**limma-voom workflow:** `limma`, `edgeR`, `data.table`, `dplyr`, `glue`, `stringr`,
`purrr`, `foreach`, `argparse`

### Input files

| File | Description |
|---|---|
| `data/count/RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt` | Raw count matrix (genes × samples); row 1 = sample IDs only (no `gene_id` column header) |
| `data/count/Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt` | Sample metadata; `id` column matches count file column names |
| `data/count/TMM_Voom_sv{1..4}_10.txt` | Surrogate variable files; positional rows align to count file columns |
| `src/q1q5_degdvg/shared/focal_genes.csv` | CSV with `gene_id` column (FlyBase IDs) + annotation columns (symbol, description, etc.) |

### Working directory

All commands below assume `code/dronetanalysis/` as the working directory.

---

## Architecture

All three Q1/Q5 workflows (GAMLSS, DESeq2, limma-voom) share a common preprocessor
and gene list. Per-method code lives in separate subfolders:

```
src/q1q5_degdvg/
  shared/
    preprocess_quintile.R       condition filter → quintile split → SV merge
    focal_genes.csv             gene_id (required) + annotation columns

  q1q5_deseq2/
    deseq2_deg_core.R           run_deseq2_deg(): DDS construction, DESeq(), lfcShrink
    run_deseq2_quintile_analysis.R  driver: parameters, argparse CLI
    compute_deg_deseq2_q1q5.snakemake
    submit_deseq2_q1q5.sh
    profiles/sge/config.yaml

  q1q5_limma_voom/
    limma_voom_deg_core.R       run_limma_voom_deg(): TMM, voom, lmFit, eBayes, topTable
    run_limma_voom_quintile_analysis.R  driver: parameters, argparse CLI
    compute_deg_limma_voom_q1q5.snakemake
    submit_limma_voom_q1q5.sh
    profiles/sge/config.yaml
```

The shared preprocessor `prepare_quintile_groups()` returns a **raw integer count
matrix** and an **annotated metadata table** — it never normalises. Each core function
applies its own normalisation (size factors for DESeq2; TMM + voom weights for limma).

---

## Step-by-Step Instructions

### Step 1: Choose focal genes

Edit `src/q1q5_degdvg/shared/focal_genes.csv`. The required column is `gene_id`;
additional annotation columns (symbol, description, etc.) are carried for reference
but not used by the workflow:

```csv
gene_id,symbol,description
FBgn0027560,ken,ken and barbie
FBgn0016694,foxo,forkhead box, sub-group O
FBgn0001235,Ilp2,Insulin-like peptide 2
```

Both Snakefiles read only the `gene_id` column — no other changes are needed.

---

### Step 2: Test run (interactive, single gene)

Run one gene interactively to check that paths and packages are correct.

**DESeq2 — test mode (1 gene, CT condition):**
```bash
RSCRIPT=/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript

$RSCRIPT src/q1q5_degdvg/q1q5_deseq2/run_deseq2_quintile_analysis.R \
  --focal-gene FBgn0027560 \
  --condition-val 1 \
  --n-test-genes 10
```

**limma-voom — test mode (1 gene, CT condition):**
```bash
$RSCRIPT src/q1q5_degdvg/q1q5_limma_voom/run_limma_voom_quintile_analysis.R \
  --focal-gene FBgn0027560 \
  --condition-val 1 \
  --n-test-genes 10
```

Expected output in `../../results/deseq2/` and `../../results/limma_voom/`:
- `FBgn0027560_CT_Q1vsQ5_deseq2_deg.csv` / `FBgn0027560_CT_Q1vsQ5_limma_voom_deg.csv`
- `FBgn0027560_CT_Q1vsQ5_deseq2_deg_sum.csv` / `FBgn0027560_CT_Q1vsQ5_limma_voom_deg_sum.csv`

---

### Step 3: Full interactive run (all genes, one focal gene)

```bash
$RSCRIPT src/q1q5_degdvg/q1q5_deseq2/run_deseq2_quintile_analysis.R \
  --focal-gene FBgn0027560 \
  --condition-val 1 \
  --full-run

$RSCRIPT src/q1q5_degdvg/q1q5_limma_voom/run_limma_voom_quintile_analysis.R \
  --focal-gene FBgn0027560 \
  --condition-val 1 \
  --full-run
```

Typical wall times for a full run (~8,786 genes):
- DESeq2: ~2–3 h with 4 MulticoreParam workers
- limma-voom: ~10–20 min (lmFit is vectorized; 1 worker is sufficient)

---

### Step 4: Snakemake dry run

Verify the workflow DAG (should show 16 child jobs + 1 `all` = 17 total):

```bash
SNAKEMAKE=/tmp/global2/caoyt/miniforge3/bin/snakemake

$SNAKEMAKE -s src/q1q5_degdvg/q1q5_deseq2/compute_deg_deseq2_q1q5.snakemake -n -p
$SNAKEMAKE -s src/q1q5_degdvg/q1q5_limma_voom/compute_deg_limma_voom_q1q5.snakemake -n -p
```

---

### Step 5: Cluster submission

Submit via the qsub master script. The master job runs Snakemake; child jobs run
the R analysis per focal gene × condition pair.

```bash
# from code/dronetanalysis/
qsub src/q1q5_degdvg/q1q5_deseq2/submit_deseq2_q1q5.sh
qsub src/q1q5_degdvg/q1q5_limma_voom/submit_limma_voom_q1q5.sh
```

Cluster resources per child job:

| Parameter | DESeq2 | limma-voom |
|---|---|---|
| Threads / slots | 4 | 1 |
| Memory | 20 GB | 10 GB |
| Walltime | 12 h | 2 h |
| Master walltime | 72 h | 24 h |
| SGE account (child) | `child_deseq2_q1q5` | `child_limma_voom_q1q5` |

---

## Output Column Reference

### DESeq2 result table (`*_deseq2_deg.csv`)

Each row is one gene. Numeric columns are suffixed with the treatment level (e.g. `_Q5`).

| Column | Type | Description |
|---|---|---|
| `gene_id` | chr | FlyBase gene ID |
| `baseMean_Q5` | num | Mean normalised count across all samples |
| `log2FC_Q5` | num | Log2 fold-change Q5 vs Q1 (MLE; not shrunk) |
| `log2FC_shrunken_Q5` | num | Log2FC after `lfcShrink(type="normal")`; NA if `--no-lfc-shrink` |
| `lfcSE_Q5` | num | Standard error of the log2FC estimate |
| `stat_Q5` | num | Wald statistic (or LRT statistic when `--test LRT`) |
| `pvalue_Q5` | num | Raw p-value |
| `padj_Q5` | num | BH-adjusted p-value |
| `class_Q5` | chr | `up_deg` / `down_deg` / `NS` |

Classification rule: `padj < 0.05` and `log2FC > 0` → `up_deg`; `log2FC ≤ 0` → `down_deg`.

### limma-voom result table (`*_limma_voom_deg.csv`)

| Column | Type | Description |
|---|---|---|
| `gene_id` | chr | FlyBase gene ID |
| `AveExpr_Q5` | num | Average log2-CPM across all samples |
| `logFC_Q5` | num | Log2 fold-change Q5 vs Q1 (from precision-weighted linear model) |
| `t_Q5` | num | Moderated t-statistic |
| `pvalue_Q5` | num | Raw p-value |
| `adj_pvalue_Q5` | num | BH-adjusted p-value |
| `B_Q5` | num | Log-odds that the gene is differentially expressed (B-statistic) |
| `class_Q5` | chr | `up_deg` / `down_deg` / `NS` |

Classification rule: `adj_pvalue < 0.05` and `logFC > 0` → `up_deg`; `logFC ≤ 0` → `down_deg`.

### Summary table (`*_deg_sum.csv`)

One row per comparison (e.g. `Q5`), columns for each class label (`up_deg`, `down_deg`,
`NS`) plus `total`.

---

## Statistical Design

Both workflows fit the same covariate structure as the GAMLSS workflow:

**Full model formula:**
```
~ expr_quintile + RNAlibBatch + RNAseqBatch + egglayBatch + platingBatch + well
  + sv1 + sv2 + sv3 + sv4
```

**Contrast of interest:** Q5 vs Q1 (reference) for `expr_quintile`.

**Batch covariates** (`RNAlibBatch`, `RNAseqBatch`, `egglayBatch`, `platingBatch`,
`well`) use **sum-to-zero (effects) coding** (`contr.sum`). This makes the intercept
represent the grand mean rather than one specific batch level, improving interpretability
of the `expr_quintile` coefficient.

**Surrogate variables** (`sv1`–`sv4`) are continuous and require no contrast recoding.

**Key implementation differences between the two methods:**

| Aspect | DESeq2 | limma-voom |
|---|---|---|
| Count model | Negative binomial GLM | Precision-weighted normal (log-CPM) |
| Normalisation | Internal size factors (median-of-ratios) | TMM factors + voom precision weights |
| Dispersion | Per-gene shrinkage toward trend (empirical Bayes) | Not modelled explicitly |
| LFC shrinkage | `lfcShrink(type="normal")` → `log2FC_shrunken` | Not applicable |
| Contrast timing | `contrasts()` on colData **before** `DESeqDataSetFromMatrix()` | `contrasts()` before `model.matrix()` |
| Coefficient name | `expr_quintile_Q5_vs_Q1` | `expr_quintileQ5` (no separator) |
| Parallelism | `MulticoreParam(4)` via `parallel=TRUE` in `DESeq()` | None needed (vectorized `lmFit`) |

For concept-level detail on how these methods compare to GAMLSS and to each other,
see [REFERENCE-03-DEG-DVG-Identification.md](REFERENCE-03-DEG-DVG-Identification.md).

---

## Directional Class Values

Both methods annotate each gene with a directional class label:

| Class | Meaning |
|---|---|
| `up_deg` | Significant DEG; higher mean in Q5 (logFC > 0) |
| `down_deg` | Significant DEG; lower mean in Q5 (logFC ≤ 0) |
| `NS` | Not significant at the chosen FDR threshold |

Unlike the GAMLSS workflow, there are no DVG classes (`up_dvg`, `down_dvg`, `Both`)
because DESeq2 and limma-voom do not model the dispersion (sigma) component.

---

## Troubleshooting

**`the model matrix is not full rank`** (DESeq2)
After Q1/Q5 filtering, some batch levels (especially `well`) can appear exclusively
in Q1 or exclusively in Q5, making the design matrix rank-deficient. DESeq2's
`checkFullRank()` hard-fails in this case. The core function pre-builds the design
matrix with `model.matrix()` and drops aliased columns via QR decomposition before
calling `DESeqDataSetFromMatrix()`, mirroring what `lm()` and limma do automatically.
Dropped columns are printed as a message in the log. If the `expr_quintile` treatment
coefficient itself is dropped, the job stops — the Q1/Q5 grouping is perfectly confounded
with batch structure for that gene, and DEG cannot be estimated.

**`Error: Expected coefficient(s) not found in resultsNames(dds)`**
Coefficient names follow `model.matrix()` format (`expr_quintileQ5`, no `_vs_`
separator). Verify that `trt_col_name` and `label_high` are consistent with the
column names produced by `model.matrix()`.

**`Error: Expected coefficient(s) not found in design matrix`** (limma-voom)
Same `model.matrix()` naming convention: `expr_quintileQ5` (no separator). The
core function checks this explicitly and prints available column names. If the factor
level name changed, update `trt_col_name` or `label_high`.

**Output files written to wrong directory**
The Snakemake rules pass `--out-dir {params.out_dir}` to each driver script.
`params.out_dir` is set to `DESEQ2_DIR` / `LIMMA_DIR` / `GAMLSS_DIR` at the top of
each Snakefile — that variable is the single source of truth. The default directory
in the PARAMETERS block of each driver script is only used for interactive runs.

**DESeq2 `mode(count_mat) must be integer`**
`prepare_quintile_groups()` returns raw integer counts, but subsetting with R can
promote integers to doubles. The core function guards against this with
`mode(count_mat) <- "integer"`.

**`lfcShrink` fails: package not found**
The `ganlss` conda environment does not have `apeglm` or `ashr`. The core function
always uses `type="normal"` (the MAP normal-prior shrinkage available in the base
`DESeq2` package). Do not change this to `apeglm` without first installing it.

**Snakemake dry run shows 0 jobs for a method that was already run**
Output CSVs already exist. Use `--forceall` to force re-execution, or delete
the relevant result files.

---

## Related Reading

- [GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md](GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md) — GAMLSS workflow (DEG + DVG)
- [REFERENCE-03-DEG-DVG-Identification.md](REFERENCE-03-DEG-DVG-Identification.md) — Concept guide: what DEG/DVG are, Q1/Q5 design, method comparison
- [OPTIMIZATION-04-Snakemake-Efficiency.md](OPTIMIZATION-04-Snakemake-Efficiency.md) — Snakemake patterns used here (executor plugin, walltime resource, master qsub)
- [FIX-03-HPC-SGE-Pipeline.md](FIX-03-HPC-SGE-Pipeline.md) — SGE-specific fixes: NFS locking, `latency-wait`, memory

---

**Last Updated:** 2026-06-09
**Status:** ✅ Active
