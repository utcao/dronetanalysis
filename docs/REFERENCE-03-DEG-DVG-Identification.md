# DEG and DVG Identification: Concepts and Method Comparison

## What Is a DEG?

A **differentially expressed gene (DEG)** shows a statistically significant difference
in **mean expression** between two groups. In the Q1/Q5 design the two groups are
samples with low focal-gene expression (Q1) vs. samples with high focal-gene expression
(Q5). A DEG in this context is a gene whose transcriptional level co-varies with —
or is buffered against — the focal gene.

- `up_deg`: significantly higher mean in Q5 (logFC > 0; the gene is upregulated when the focal gene is highly expressed)
- `down_deg`: significantly lower mean in Q5 (logFC ≤ 0)

## What Is a DVG?

A **differentially variable gene (DVG)** shows a statistically significant difference
in **expression variability** between two groups, with no requirement that the means differ.
Variability is quantified as the biological coefficient of variation (BCV):

```
BCV = sqrt(exp(sigma))
```

where `sigma` is the NBI dispersion parameter. A DVG identifies genes whose
expression becomes more or less tightly controlled as the focal gene's activity changes —
a signal of altered buffering capacity.

DVG detection requires a model with a separate equation for the dispersion (variance)
component. Only the **GAMLSS** workflow provides this; DESeq2 and limma-voom test the
mean only and cannot detect pure DVGs.

- `up_dvg`: significantly higher BCV in Q5 (more variable when focal gene is highly expressed)
- `down_dvg`: significantly lower BCV in Q5 (more tightly buffered)

---

## Why the Q1/Q5 Design?

The Q1/Q5 stratification uses a focal gene's expression level as a **proxy for the
activity state** of a modifier gene or pathway. Samples with high focal-gene expression
(Q5) are presumed to have high modifier activity; Q1 samples have low activity. This
design:

1. **Avoids direct genetic manipulation** — it works on observational expression data
   from natural variants (DGRP inbred lines).
2. **Captures continuous variation** — the focal gene's expression varies across inbred
   lines and across individuals within lines; Q1/Q5 extracts the extremes of this
   natural gradient.
3. **Drops the middle 60%** — removing samples near the median improves statistical
   power by maximising the contrast between groups.
4. **Preserves the full transcriptome as the response** — the focal gene defines the
   groups; all other ~8,786 genes are tested for differential expression or variability.

The implicit hypothesis is: *genes that change systematically between individuals with
low vs. high focal-gene expression are candidates for regulation downstream of — or in
parallel to — the focal gene.*

---

## The Three Models: GAMLSS, DESeq2, limma-voom

### GAMLSS / NBI

**What it models:** Two separate linear predictors — one for the mean (mu, log link)
and one for the overdispersion (sigma, log link) — fitted jointly to raw counts under
the negative binomial (NBI) family.

**How DEG/DVG are tested:** Three models are fitted per gene:
- `m_full`: full model (covariates + treatment in mu; treatment in sigma)
- `m_sigma_int`: intercept-only sigma (removes treatment from sigma)
- `m_mu_nontrt`: non-treatment mu (removes treatment from mu)

Two LR tests:
- DVG: `2 × (logLik(m_full) − logLik(m_sigma_int))`, df = Δ parameters
- DEG: `2 × (logLik(m_full) − logLik(m_mu_nontrt))`, df = Δ parameters

**Normalisation:** TMM offsets passed as a fixed term (`o(offset)`) in the log-link mu formula.

**Strengths:** Only method that can detect pure DVG; allows asymmetric sigma–mean coupling.

**Weaknesses:** Slowest (~48 h/gene on cluster with 4 cores); convergence can fail
for genes with sparse counts or unusual distributions; no built-in empirical Bayes
for dispersion.

---

### DESeq2 / Negative Binomial GLM

**What it models:** A single GLM for the count mean under the negative binomial family
with a log link. The dispersion is estimated per-gene and then shrunk toward a
genome-wide trend using an empirical Bayes procedure.

**How DEG is tested:** By default, a Wald test on the `expr_quintile_Q5_vs_Q1`
coefficient. Alternatively, a likelihood ratio test against a reduced model
(`formula_reduced`) can be used via `--test LRT`.

**Normalisation:** Internally estimated size factors (median-of-ratios method); these
replace TMM and do not require pre-computing offsets.

**LFC shrinkage:** `lfcShrink(type="normal")` produces a MAP estimate that shrinks
large LFCs from genes with high uncertainty toward zero. This is stored in
`log2FC_shrunken`; the unshrunk MLE is in `log2FC`.

**Contrast timing:** `contrasts()` must be stamped on colData factor columns **before**
`DESeqDataSetFromMatrix()` because DESeq2 copies colData at construction. Unlike GAMLSS
which receives `contrasts.arg` into `model.matrix()` at fit time, a late contrasts call
has no effect on a DESeq2 object.

**Strengths:** Empirical Bayes dispersion shrinkage improves stability for low-count
genes; widely adopted and well-validated; built-in size factors.

**Weaknesses:** Slower than limma-voom because each gene's dispersion is estimated
iteratively; benefits most from moderate-to-large sample sizes; does not test variability.

---

### limma-voom / Precision-Weighted Linear Model

**What it models:** After TMM normalisation, raw counts are converted to log2-CPM
values. The mean-variance relationship of log-CPM values is estimated and used to
assign precision weights to each observation. These weights are passed to a standard
weighted least-squares linear model.

**How DEG is tested:** A moderated t-test via `eBayes()` on the `expr_quintileQ5`
coefficient (note: `model.matrix()` with treatment coding produces `expr_quintileQ5`
with no separator, unlike DESeq2's `expr_quintile_Q5_vs_Q1`). The moderated F and
t statistics borrow information across genes for variance estimation.

**Normalisation:** TMM factors computed via `edgeR::calcNormFactors()`, then passed
to `voom()` which also estimates and applies precision weights per-observation.

**Coefficient naming:** `model.matrix()` concatenates factor name and level directly:
```
~expr_quintile → columns: (Intercept)  expr_quintileQ5
```
The core function verifies the coefficient name explicitly and stops with a diagnostic
message if it is missing.

**Strengths:** Fastest of the three (~minutes per gene on a single core); fully
vectorized `lmFit` operates on all genes simultaneously; mature, well-tested; B-statistic
(log-odds) provides additional ranking beyond p-value.

**Weaknesses:** Log-normal approximation to count data; less powerful than DESeq2 for
very low counts; no dispersion (variability) testing.

---

## Normalisation Strategies Compared

| Method | Strategy | Notes |
|---|---|---|
| GAMLSS | TMM offsets as `o(log(lib.size × norm.factor))` in log-link mu | Added as fixed term in model formula |
| DESeq2 | Median-of-ratios size factors (internal to `DESeq()`) | Applied before GLM fitting; no user action needed |
| limma-voom | TMM factors → `calcNormFactors()` → passed to `voom()` | voom also computes per-observation precision weights |

All three use library-size correction. DESeq2 uses its own size factor algorithm
(robust to outliers via the geometric mean); GAMLSS and limma-voom use edgeR's TMM
(trims extreme fold-changes before computing normalisation factors).

---

## Contrast Coding for Batch Covariates

All three models include the same batch covariates:
`RNAlibBatch`, `RNAseqBatch`, `egglayBatch`, `platingBatch`, `well`

These use **sum-to-zero (effects) coding** (`contr.sum`). Under treatment coding
(the R default), the intercept represents one arbitrary baseline batch level —
the `expr_quintile` coefficient then estimates the Q1/Q5 contrast at that specific
batch, not at the grand mean. Sum-to-zero coding makes the intercept the grand mean
across all batch levels, so the `expr_quintile` coefficient is interpretable as
the average Q1/Q5 effect across the whole experiment.

**Why this matters in practice:** If one batch had unusually high or low expression,
treatment coding would absorb that into the intercept and push it into the treatment
effect estimate. Effects coding distributes the batch contribution symmetrically and
keeps the treatment coefficient clean.

Surrogate variables (`sv1`–`sv4`) are **continuous** and require no contrast recoding.

**Timing of contrast application:**

| Method | When to apply | How |
|---|---|---|
| GAMLSS | At fit time | Passed as `contrasts.arg` to `model.matrix()` inside `fit_one()` |
| DESeq2 | **Before** `DESeqDataSetFromMatrix()` | `contrasts(meta_df[[col]]) <- contr.sum(...)` |
| limma-voom | Before `model.matrix()` | `contrasts(meta_df[[col]]) <- contr.sum(...)` |

DESeq2 is the strictest: it copies `colData` at construction, so contrasts set
after construction are silently ignored.

---

## Directional Annotation

All three methods annotate each gene with a directional class label stored in a
`class` column (suffixed with the treatment level, e.g. `class_Q5`).

### DESeq2 and limma-voom

```
up_deg   — padj < fdr_threshold AND logFC > 0
down_deg — padj < fdr_threshold AND logFC ≤ 0
NS       — padj ≥ fdr_threshold (or NA)
```

### GAMLSS

```
up_deg   — q_mu < fdr AND q_sigma ≥ fdr AND logFC_mu > 0
down_deg — q_mu < fdr AND q_sigma ≥ fdr AND logFC_mu ≤ 0
up_dvg   — q_mu ≥ fdr AND q_sigma < fdr AND logFC_BCV > 0
down_dvg — q_mu ≥ fdr AND q_sigma < fdr AND logFC_BCV ≤ 0
Both     — q_mu < fdr AND q_sigma < fdr   (direction: read logFC_mu and logFC_BCV columns)
NS       — neither q_mu nor q_sigma < fdr
```

`Both` genes are significant for both mean and variability change. The `logFC_mu` and
`logFC_BCV` columns contain the directional information; `Both` is kept undirected in
the class label because the biological interpretation requires examining both dimensions.

---

## When to Prefer Each Method

| Scenario | Recommended |
|---|---|
| Need DVG (variability) detection | GAMLSS |
| Sparse counts, small sample sizes | DESeq2 (best dispersion shrinkage) |
| Large dataset, speed is priority | limma-voom |
| Want LFC shrinkage for ranking | DESeq2 (`log2FC_shrunken`) |
| Cross-validating with two methods | DESeq2 + limma-voom |
| Full biological characterisation | All three |

For the Q1/Q5 stratified design used here (typically ~80–160 samples per group after
quintile splitting from ~395–790 condition-filtered samples), all three methods have
adequate power. DESeq2 and limma-voom typically produce concordant DEG lists;
discordant genes are worth examining manually.

---

## References

- Love MI, Huber W, Anders S (2014). *Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2.* Genome Biology 15:550.
  doi:10.1186/s13059-014-0550-8

- Ritchie ME, Phipson B, Wu D, et al. (2015). *limma powers differential expression
  analyses for RNA-sequencing and microarray studies.* Nucleic Acids Research 43:e47.
  doi:10.1093/nar/gkv007

- Rigby RA, Stasinopoulos DM (2005). *Generalized additive models for location, scale
  and shape.* Applied Statistics 54:507–554. doi:10.1111/j.1467-9876.2005.00510.x

- Robinson MD, McCarthy DJ, Smyth GK (2010). *edgeR: a Bioconductor package for
  differential expression analysis of digital gene expression data.*
  Bioinformatics 26:139–140. doi:10.1093/bioinformatics/btp616

---

## Related Reading

- [GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md](GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md) — GAMLSS workflow guide
- [GUIDE-20-DESeq2-limma-voom-DEG-Q1Q5.md](GUIDE-20-DESeq2-limma-voom-DEG-Q1Q5.md) — DESeq2 and limma-voom workflow guide
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — Statistical methods overview for the broader project

---

**Last Updated:** 2026-06-09
**Status:** ✅ Active
