# Guide: Sample Quintile Overlap — Expression-Gradient Partitioning

**Category:** Guide — Per-sample characterisation of systematic expression ranking across the transcriptome

---

## Purpose

This analysis asks: **does each sample occupy the same expression tier (low / mid / high) consistently across genes?**

For every gene in the expression matrix, samples are ranked by expression value and assigned to one of five equal-width quintiles:

| Quintile | Expression tier |
|---|---|
| Q1 | Bottom 20% — lowest-expressing samples for this gene |
| Q2 | 20–40% |
| Q3 | 40–60% (middle) |
| Q4 | 60–80% |
| Q5 | Top 20% — highest-expressing samples for this gene |

The assignment is repeated across all genes and the counts accumulated per sample.  Under the null hypothesis of no systematic bias, each sample lands in each quintile with equal probability ≈ 0.20.  Samples that consistently land in Q1 or Q5 across many genes represent individuals with a globally depressed or elevated transcriptome relative to the population.

---

## Uniformity Test

For each sample, a **chi-square goodness-of-fit test** (df = 4) tests whether the observed Q1–Q5 counts differ from the uniform expectation:

```
H0: freq(Q1) = freq(Q2) = freq(Q3) = freq(Q4) = freq(Q5) = 0.20
chi_sq_i = sum_q (count_iq - n_genes/5)^2 / (n_genes/5)
```

BH correction is applied across all samples.  Two supplementary measures quantify the *degree* of non-uniformity:

| Measure | Formula | Range | Interpretation |
|---|---|---|---|
| `entropy_norm` | −∑ f·log₂(f) / log₂(5) | 0–1 | 1 = perfectly uniform; lower = more skewed |
| `max_deviation` | max \|f_q − 0.20\| | 0–0.80 | Worst-quintile deviation from expected |

---

## Script

`code/dronetanalysis/src/scripts/15analysis/plot_sample_quintile_overlap.R`

**Required packages:** `data.table`, `ggplot2`, `argparse`, `pheatmap` — all available in the `dronetanalysis` conda environment.

### Arguments

| Argument | Default | Description |
|---|---|---|
| `--expr-file` | required | VOOM expression matrix (tab-separated, genes × samples) |
| `--output-dir` | required | Destination directory (created if absent) |
| `--condition-label` | `"Condition"` | Label for plot titles |
| `--seed` | `42` | Random seed for tie-breaking in rank assignment |

### Usage

**Control condition:**
```bash
conda run -n dronetanalysis \
Rscript code/dronetanalysis/src/scripts/15analysis/plot_sample_quintile_overlap.R \
  --expr-file    data/processed/VOOM/voomdataCtrl.txt \
  --output-dir   results/quintile_overlap_ct \
  --condition-label "Control"
```

**HS condition:**
```bash
conda run -n dronetanalysis \
Rscript code/dronetanalysis/src/scripts/15analysis/plot_sample_quintile_overlap.R \
  --expr-file    data/processed/VOOM/voomdataHS.txt \
  --output-dir   results/quintile_overlap_hs \
  --condition-label "HS"
```

---

## Output Files

All outputs land in `--output-dir`.

### Tables

#### `sample_quintile_counts.tsv`

One row per sample, sorted by `Q1_freq` descending.  Use this to identify samples that are systematically low- or high-expressing across the transcriptome.

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `Q1`–`Q5` | Raw count of genes for which this sample fell in each quintile |
| `Q1_freq`–`Q5_freq` | Frequency = count / n_genes (expected ≈ 0.20) |
| `chisq_stat` | Chi-square test statistic (df = 4) |
| `chisq_p` | Unadjusted p-value |
| `chisq_padj` | BH-adjusted p-value across all samples |
| `entropy_norm` | Normalised Shannon entropy (1 = uniform, lower = more skewed) |
| `max_deviation` | Maximum \|freq_q − 0.20\| across Q1–Q5 |

**Typical usage — find Q1-biased samples:**
```r
dt <- data.table::fread("results/quintile_overlap_ct/sample_quintile_counts.tsv")
dt[Q1_freq > 0.30 & chisq_padj < 0.05]   # samples > 30% Q1, significantly non-uniform
```

#### `sample_noise_metrics.tsv`

Per-sample transcriptomic noise computed across all genes.  One row per sample, sorted by `mad_log2` descending.

| Column | Scale | Description |
|---|---|---|
| `sample_id` | — | Sample identifier |
| `n_genes` | — | Number of non-missing genes used |
| `mean_log2` | log2 | Mean expression across genes, then apply log2 |
| `median_log2` | log2 | Median expression across genes, then apply log2 |
| `sd_log2` | log2 | Standard deviation across genes, then apply log2 |
| `mad_log2` | log2 | Median absolute deviation — robust noise measure, then apply log2 |
| `iqr_log2` | log2 | Interquartile range across genes, then apply log2 |
| `cv_log2` | log2 | Coefficient of variation (sd / mean) on log2 scale, then apply log2 |


**Why two CV scales?**  
`cv_log2` is straightforward and depth-comparable.  `cv2_linear` uses the linear (CPM-like) scale and is the standard transcriptomic noise measure for count data; it is also consistent with the per-gene CV² computed in `compute_mad_variability_ranks.R`, enabling direct cross-table comparisons.

**Typical usage:**
```r
noise <- data.table::fread("results/quintile_overlap_ct/sample_noise_metrics.tsv")
noise[mad_log2 > quantile(mad_log2, 0.95)]   # top 5% noisiest samples
```

---

#### `gene_q1q5_assignments.rds`

An R named list mapping each gene to the samples assigned to its Q1 and Q5.  Contains the **top 5 000 genes** ranked by `gene_max_dev` (see [Gene Filtering](#gene-filtering) below).  Use for O(1) lookup by gene ID.

```r
assignments <- readRDS("results/quintile_overlap_ct/gene_q1q5_assignments.rds")

# Which samples are in Q1 for this gene?
assignments[["FBgn0000008"]]$Q1

# Which samples are in Q5?
assignments[["FBgn0000008"]]$Q5

# How many of the top-5000 genes place sample X in Q5?
q5_count <- sum(sapply(assignments, function(g) "112_G11" %in% g$Q5))
```

Storage: xz-compressed RDS, 1.8 MB for 5 000 genes × 938 samples (vs. 9.3 MB for the equivalent long-format TSV.gz with all genes).

#### `gene_q1q5_wide.tsv.gz`

Human-readable companion to the RDS.  Same top-5 000 gene filter.  One row per gene; Q1/Q5 sample IDs are pipe-separated strings.  Useful for inspection in a spreadsheet or cross-language use.

| Column | Description |
|---|---|
| `gene_id` | FlyBase gene ID |
| `gene_max_dev` | Gene-level extremity score (see [Gene Filtering](#gene-filtering)) |
| `Q1_n_samples` | Number of samples in Q1 for this gene |
| `Q1_samples` | Pipe-separated Q1 sample IDs |
| `Q5_n_samples` | Number of samples in Q5 for this gene |
| `Q5_samples` | Pipe-separated Q5 sample IDs |

### Plots

#### `quintile_distribution.pdf`

Five-panel density plot (one panel per quintile) of per-sample membership frequencies.  The dashed line marks the 0.20 expected frequency.  The subtitle reports how many samples are significantly non-uniform (BH < 0.05).

A broad or right-shifted Q1 distribution indicates that many samples fall below the median expression for most genes.  A narrow distribution centred on 0.20 indicates a well-balanced cohort with no strong individual biases.

#### `quintile_heatmap.pdf`

Samples (rows, unlabelled) × Q1–Q5 (columns).  Rows are clustered by complete-linkage hierarchical clustering on Euclidean distance of the frequency profiles.  Columns retain the natural Q1→Q5 order.

Colour scale is diverging and centred at 0.20:
- **Blue** — sample is under-represented in this quintile relative to expectation
- **White** — ~uniform (freq ≈ 0.20)
- **Red** — sample is over-represented in this quintile

Clusters of red-Q1 / blue-Q5 rows identify samples with a globally suppressed transcriptome.  Clusters of blue-Q1 / red-Q5 identify globally elevated samples.  Scattered samples with no cluster affinity (near-white rows) show no systematic quintile bias.

---

## Gene Filtering

The `gene_q1q5_assignments.rds` and `gene_q1q5_wide.tsv.gz` files store only the **top 5 000 genes** — those whose Q1/Q5 groupings are most informative about the global sample-expression pattern.

**Why not score genes the same way as samples?**  
Every gene has exactly ~n_samples/5 members per quintile by construction, so the gene's own marginal Q1–Q5 distribution is always uniform and gives no information.

**`gene_max_dev` definition:**

```
gene_q1_bias = mean(Q1_freq of the gene's Q1 samples) - 0.20
gene_q5_bias = mean(Q5_freq of the gene's Q5 samples) - 0.20
gene_max_dev = max(gene_q1_bias, gene_q5_bias)
```

`Q1_freq` and `Q5_freq` here are the *sample-level* frequencies from `sample_quintile_counts.tsv` — how often each sample sits in Q1 (or Q5) across all genes.

A gene with high `gene_max_dev` draws its Q1 group disproportionately from samples that are *globally* Q1-dominant (e.g. a sample with Q1_freq = 0.36 is in Q1 for 36% of all genes; if it also appears in this gene's Q1, that reinforces the pattern).  These genes most faithfully reflect the transcriptome-wide expression-tier structure, making them the most useful for downstream analyses.

**Vectorised computation:**

```r
gene_q1_bias = (quintile_mat == 1) %*% Q1_freq_vec / rowSums(quintile_mat == 1) - 0.20
```

No per-gene loop — the matrix product gives all 8 763 gene scores at once.

---

## Implementation Notes

### `ceiling(rank * 5 / n)` — rank-to-quintile formula

This is the standard rank-to-bin mapping used by `dplyr::ntile()`.  For ranks 1..n assigned to K bins:

```
quintile = ceiling(rank * K / n)
```

It produces bin sizes of `floor(n/K)` or `ceil(n/K)`, differing by at most 1 — the smallest possible imbalance for any equal-width partition.

Why not `cut(rank, quantile-breaks)`?  Both approaches are equivalent in result, but `cut()` allocates a factor and computes quantiles separately; the formula avoids both, making it faster in the `apply(expr_mat, 1, ...)` loop over thousands of genes.

### `pmin(..., 5L)` — vectorised cap

`pmin` takes the element-wise minimum of two vectors.  It is needed because `ceiling(n * 5 / n)` can evaluate to 6 under floating-point rounding (e.g. when `n * 5 / n` yields `5.000...001`).  `pmin` is preferred over `ifelse(q > 5, 5L, q)` because it avoids allocating an intermediate logical vector.

### `vapply` over `lapply` / `sapply`

```r
count_mat <- vapply(1:5, function(q) colSums(...), numeric(n_samples))
```

- **`lapply`** returns a list — needs `do.call(cbind, ...)` to become a matrix.
- **`sapply`** infers its return type from the first call; if the function accidentally returns the wrong length it silently produces a list instead of a matrix.
- **`vapply`** requires an explicit `FUN.VALUE` template (`numeric(n_samples)`).  R pre-allocates the result matrix to that shape and raises an error immediately if any column has the wrong type or length — making type bugs impossible to miss.

---

## Interpreting Results

**927/938 samples non-uniform (CT, chi-sq BH < 0.05):**  
With ~8 800 genes, the test is highly powered.  Near-universal significance means most samples deviate at least slightly from uniform, but this does not imply biological relevance.  Focus on samples with `Q1_freq > 0.30` or `entropy_norm < 0.95` to identify practically meaningful outliers.

**Q1 and Q5 show wider spread than Q2–Q4:**  
The sd of Q1/Q5 frequencies (~0.055/0.064) exceeds that of middle quintiles (~0.034–0.046).  This is expected: expression extremes are more variable across the cohort because a single high-leverage sample or batch effect can dominate the tails of each gene's distribution.

**Heatmap clusters:**  
Row clusters correspond to samples that share the same global expression tier.  A "Q1-high, Q5-low" cluster may represent samples with low overall expression (possibly biological, e.g. a distinct genotype, or technical, e.g. lower library depth).  Cross-reference with metadata (line, replicate, batch) to distinguish.

**Using `entropy_norm` vs. chi-sq p-value:**  
`chisq_padj` answers *is this sample non-uniform?* — it will be significant for almost all samples with large n_genes.  `entropy_norm` answers *how non-uniform?* and is comparable across datasets of different sizes.  Use `entropy_norm < 0.95` as a practical threshold for "meaningfully skewed."

---

## Related Scripts

| Script | Purpose |
|---|---|
| `compute_mad_transcriptomic_variability.R` | Gene-level MAD: does a focus gene's expression level co-segregate with global transcriptomic variability? |
| `compute_mad_permute_transcriptomic_noise.R` | Permutation validation of gene-level MAD results |
| `plot_itv_sample_comparison.R` | Per-sample ITV scores for LOW vs HIGH sub-groups of a specific focus gene |
