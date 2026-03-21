# Expression Variability Analysis Guide

## Overview

This guide explains two standalone R scripts that quantify how transcriptomic
variability differs between the **LOW** and **HIGH** sample groups defined by a
focus gene's expression rank. These analyses complement the differential
co-expression network results by asking a related but distinct question:

> *When a focus gene is highly expressed, does the rest of the transcriptome
> become more or less variable?*

Two complementary views are provided:

| Script | Unit of observation | Variability measure |
|--------|--------------------|--------------------|
| `plot_gene_mad_variability.R` | One dot per gene | MAD across group samples |
| `plot_sample_variability.R`  | One dot per sample | ITV — deviation from population mean |

Both scripts are self-contained CLI tools. They derive the LOW/HIGH sample
groups directly from the expression file (no pipeline HDF5 outputs required)
and produce a PDF boxplot with a Wilcoxon rank-sum test annotation.

---

## Prerequisites

- R ≥ 4.0 with the following packages installed:
  ```r
  install.packages(c("data.table", "dplyr", "ggplot2", "argparse"))
  ```
- A VOOM-normalized expression matrix (tab-separated, genes × samples). The
  full matrices are at `data/processed/VOOM/voomdataCtrl.txt` and
  `voomdataHS.txt`; test subsets are in `dataset/test/`.
- The FlyBase ID of the focus gene driving the LOW/HIGH split.

---

## Concepts

### LOW and HIGH sample groups

For focus gene **g**, the `n_samples × low_frac` samples with the *lowest*
expression of **g** form the **LOW** group, and the `n_samples × high_frac`
samples with the *highest* expression form the **HIGH** group. By default
`low_frac = high_frac = 0.2` (bottom and top 20%).

This matches the definition used in Stage 1 of the bootstrap pipeline
(`01get_extreme_pop_bootstrap.py`), ensuring consistency between network
analysis and variability analysis.

### Gene-level MAD (Script 1)

For each gene *other than* the focus gene, **MAD** (Median Absolute Deviation)
is computed across the samples in each group:

```
MAD_low[gene]  = mad( expr[gene, low_samples]  )
MAD_high[gene] = mad( expr[gene, high_samples] )
```

Each dot in the boxplot represents one gene. The comparison asks: are gene
expression values more or less spread across samples when the focus gene is
highly expressed?

### Sample-level ITV (Script 2)

For each sample, an **Individual Transcriptomic Variability (ITV)** score
measures how far that sample deviates from the population average across all
genes. Because the expression matrix is on a log scale, subtracting the
gene-wise mean is equivalent to log(A / geometric_mean):

```
mean_expr[gene]  = mean( expr[gene, ALL samples] )   # population baseline
deviation[gene]  = |expr[gene, sample] - mean_expr[gene]|
ITV(sample)      = summary_fn( deviation )            # median / mean / sum
```

Each dot represents one sample. The comparison asks: do samples in the HIGH
group show more or less individual deviation from the population norm?

---

## Step-by-Step Instructions

### 1. Script 1 — Gene-level MAD variability

**Basic usage:**
```bash
Rscript src/scripts/15analysis/plot_gene_mad_variability.R \
  --expr-file data/processed/VOOM/voomdataCtrl.txt \
  --focus-gene FBgn0002563 \
  --output-dir results/variability \
  --condition-label "Control" \
  --gene-symbol "l(2)gl" \
  --save-csv
```

**Key arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--expr-file` | *(required)* | VOOM expression matrix (TSV) |
| `--focus-gene` | *(required)* | FlyBase ID driving the LOW/HIGH split |
| `--output-dir` | *(required)* | Output directory for PDF and CSV |
| `--low-frac` | 0.2 | Fraction of samples in LOW group |
| `--high-frac` | 0.2 | Fraction of samples in HIGH group |
| `--condition-label` | "Condition" | Label shown in plot title |
| `--gene-symbol` | *none* | Human-readable gene name for plot label |
| `--gene-subset` | *none* | Path to text file of gene IDs (one per line) to restrict MAD computation |
| `--partner-type` | "all" | `all` / `direct` / `indirect` — see [Partner Restriction](#partner-restriction-future) |
| `--network-file` | *none* | Required when `--partner-type` ≠ all (future) |
| `--save-csv` | *flag* | Also write MAD values to CSV |

**Output files:**
- `{output_dir}/{focus_gene}_mad_variability.pdf`
- `{output_dir}/{focus_gene}_mad_variability.csv` (if `--save-csv`)

---

### 2. Script 2 — Sample-level ITV variability

**Basic usage:**
```bash
Rscript src/scripts/15analysis/plot_sample_variability.R \
  --expr-file data/processed/VOOM/voomdataCtrl.txt \
  --focus-gene FBgn0002563 \
  --output-dir results/variability \
  --condition-label "Control" \
  --summary-metric median \
  --save-csv
```

**Key arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--expr-file` | *(required)* | VOOM expression matrix (TSV) |
| `--focus-gene` | *(required)* | FlyBase ID driving the LOW/HIGH split |
| `--output-dir` | *(required)* | Output directory for PDF and CSV |
| `--low-frac` | 0.2 | Fraction of samples in LOW group |
| `--high-frac` | 0.2 | Fraction of samples in HIGH group |
| `--condition-label` | "Condition" | Label shown in plot title |
| `--gene-symbol` | *none* | Human-readable gene name for plot label |
| `--summary-metric` | "median" | How to collapse per-gene deviations: `median` / `mean` / `sum` |
| `--gene-subset` | *none* | Path to text file of gene IDs to restrict ITV computation |
| `--save-csv` | *flag* | Also write ITV values to CSV |

**Output files** (metric name embedded so parallel runs do not overwrite):
- `{output_dir}/{focus_gene}_sample_variability_median.pdf`
- `{output_dir}/{focus_gene}_sample_variability_median.csv` (if `--save-csv`)

---

### 3. Running all focus genes

The 15 focus genes are listed under `gene_subset:` in
`config/ct_voom_snakemake.yaml` and `config/hs_voom_snakemake.yaml`. To run
both scripts for every focus gene in a condition, use a shell loop:

```bash
# Extract focus genes from YAML (requires yq)
GENES=$(yq '.gene_subset[]' config/ct_voom_snakemake.yaml)

for gene in $GENES; do
  Rscript src/scripts/15analysis/plot_gene_mad_variability.R \
    --expr-file data/processed/VOOM/voomdataCtrl.txt \
    --focus-gene "$gene" \
    --output-dir results/variability_ct \
    --condition-label "Control" \
    --save-csv

  Rscript src/scripts/15analysis/plot_sample_variability.R \
    --expr-file data/processed/VOOM/voomdataCtrl.txt \
    --focus-gene "$gene" \
    --output-dir results/variability_ct \
    --condition-label "Control" \
    --summary-metric median \
    --save-csv
done
```

For the heat-shock condition, replace `ct` with `hs` throughout.

---

### 4. Restricting to a gene subset

Both scripts accept `--gene-subset` pointing to a plain-text file with one
gene ID per line. This lets you focus MAD or ITV calculations on a biologically
relevant set (e.g. chaperones, a WGCNA module, known interaction partners):

```bash
# Compute MAD only for chaperone genes
Rscript src/scripts/15analysis/plot_gene_mad_variability.R \
  --expr-file data/processed/VOOM/voomdataCtrl.txt \
  --focus-gene FBgn0002563 \
  --output-dir results/variability \
  --gene-subset dataset/flybase/FlyBase_IDs_CHPN_CHAPERONES_AND_CO-CHAPERONES.txt \
  --condition-label "Control"
```

The LOW/HIGH sample split is always derived from the full expression matrix
regardless of `--gene-subset`.

---

## Plot Elements

Both plots share the same visual style:

| Element | Meaning |
|---------|---------|
| Blue violin + box (LOW) | Distribution for the low-expression sample group |
| Magenta violin + box (HIGH) | Distribution for the high-expression sample group |
| Dots | Individual observations (genes or samples) |
| Subtitle | Wilcoxon rank-sum p-value (`p < 0.001 ***`, etc.) |
| Star + bracket | Significance annotation above the boxplots |
| Caption | Group sizes and (for Script 2) the summary metric used |

**Dot interpretation differs between scripts:**
- Script 1: each dot = one gene's MAD value (expect thousands of dots, shown small)
- Script 2: each dot = one sample's ITV score (expect ~k_low + k_high dots, shown larger)

---

## Partner Restriction (Future)

Script 1 (`plot_gene_mad_variability.R`) accepts `--partner-type direct` or
`--partner-type indirect`. This will restrict the MAD gene pool to **degree-1**
or **degree-2** neighbors of the focus gene in the differential co-expression
network, allowing the variability question to be asked specifically within the
focus gene's interaction neighborhood.

This feature is currently a **stub** — it requires a `--network-file` argument
pointing to a partner gene list or adjacency file derived from the bootstrap
pipeline output. The interface is designed and reserved; implementation will
follow once differential network outputs are finalized.

```bash
# Future usage (not yet active):
Rscript src/scripts/15analysis/plot_gene_mad_variability.R \
  --expr-file data/processed/VOOM/voomdataCtrl.txt \
  --focus-gene FBgn0002563 \
  --output-dir results/variability \
  --partner-type direct \
  --network-file results_ct_voom/networks/0001_FBgn0002563.h5
```

---

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|---------|
| `Focus gene not found` | FlyBase ID missing from expression matrix | Confirm the gene ID is in the first column; check for trailing whitespace |
| `LOW group has < 2 samples` | Matrix too small or `--low-frac` too small | Use test data with ≥ 10 samples, or increase `--low-frac` |
| `Warning: partner-type 'direct' not yet implemented` | Partner restriction stub triggered | Use `--partner-type all` or omit the argument for now |
| `argparse` error: `could not find function "ArgumentParser"` | `argparse` R package not installed | `install.packages("argparse")` |
| PDF produced but violins are very flat | Too few samples or low variance genes | Expected with small test datasets; use full expression matrices for production |

---

## Related Reading

- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) — Full pipeline from expression data to differential networks; Stage 1 defines the same LOW/HIGH groups used here
- [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) — Config reference for `low_frac`, `high_frac`, and `gene_subset` parameters
- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) — Network topology metrics that complement these variability results
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — Statistical background for the bootstrap pipeline and significance testing

---

**Last Updated:** 2026-03-21
**Scripts:** `src/scripts/15analysis/plot_gene_mad_variability.R`, `src/scripts/15analysis/plot_sample_variability.R`
**Status:** ✅ Both scripts implemented; partner-type extension is a planned stub
