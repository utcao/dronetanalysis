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
  All packages are already present in the `dronet` conda environment.

- **Python on PATH** (required by R's `argparse` package): when running
  inside a conda environment, pass the Python path explicitly:
  ```bash
  # Linux / HPC
  export PYTHON=$(which python)

  # Windows + conda (dronet env)
  conda run -n dronet env PYTHON="$CONDA_PREFIX/python.exe" Rscript ...
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

### 3. Running via Snakemake (Stage 6 — recommended for gene_subset configs)

When `gene_subset` is set in a config file, the Snakemake pipeline can run both
variability scripts automatically as **Stage 6**, one job per focus gene in parallel:

```yaml
# In your config YAML (e.g. config/ct_voom_snakemake.yaml):
skip_variability: false        # enable Stage 6
condition_label: "Control (VOOM)"
variability_save_csv: false    # set true to also write CSV files
```

```bash
# Run Stage 6 alongside the full pipeline:
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    -j 15

# Or run Stage 6 only (skipping network stages 2a–3b):
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --forcerun gene_mad_variability gene_sample_variability \
    -j 15
```

Outputs are written to `{out_dir}/variability/`. Snakemake handles parallelism
and skips genes whose output PDFs already exist.

> **Requirements:** `expr_tsv` must be set in the config (the R scripts read
> the expression TSV directly, not the pipeline HDF5). `gene_subset` must be
> non-empty, otherwise Stage 6 is silently skipped.

---

### 4. Batch analysis: both conditions + Excel summary

For a complete cross-condition analysis (CT + HS, all three metrics, Excel output)
use `run_variability_batch.py` with `config/variability_batch.yaml`.

#### Standalone usage

```bash
# Dry run — print all 152 Rscript commands without executing:
python code/dronetanalysis/src/scripts/15analysis/run_variability_batch.py \
    --config config/variability_batch.yaml --dry-run

# Full run (4 parallel workers):
python code/dronetanalysis/src/scripts/15analysis/run_variability_batch.py \
    --config config/variability_batch.yaml --jobs 4

# Skip genes that already have a PDF (resume partial run):
python code/dronetanalysis/src/scripts/15analysis/run_variability_batch.py \
    --config config/variability_batch.yaml --skip-existing --jobs 4

# Re-generate Excel only (all PDFs and CSVs already exist):
python code/dronetanalysis/src/scripts/15analysis/run_variability_batch.py \
    --config config/variability_batch.yaml --excel-only
```

#### Job counts

| Type | Jobs |
|------|------|
| Gene-level MAD | 2 conditions × 19 genes = 38 |
| Sample-level ITV | 2 conditions × 19 genes × 3 metrics = 114 |
| **Total** | **152** |

#### Output layout

```
results/variability_combined/
├── CT/
│   ├── FBgn0039562_mad_variability.pdf / .csv
│   ├── FBgn0039562_sample_variability_median.pdf / .csv
│   ├── FBgn0039562_sample_variability_mean.pdf / .csv
│   ├── FBgn0039562_sample_variability_sum.pdf / .csv
│   └── ... (same for all 19 genes)
├── HS/
│   └── ... (same structure)
└── variability_summary.xlsx
```

#### Excel summary (`variability_summary.xlsx`)

**Sheet `gene_mad`** — one row per (focus_gene × condition):

| Column | Description |
|--------|-------------|
| `focus_gene`, `symbol`, `condition` | Gene and condition identifiers |
| `n_genes` | Number of genes in MAD comparison |
| `mean_mad_low` / `mean_mad_high` | Mean MAD across LOW / HIGH groups |
| `median_mad_low` / `median_mad_high` | Median MAD |
| `wilcoxon_p` | Two-sided Wilcoxon rank-sum p-value |

**Sheet `sample_itv`** — one row per (focus_gene × condition × metric):

| Column | Description |
|--------|-------------|
| `focus_gene`, `symbol`, `condition`, `metric` | Identifiers |
| `n_low`, `n_high` | Sample counts per group |
| `mean_itv_low` / `mean_itv_high` | Mean ITV score per group |
| `median_itv_low` / `median_itv_high` | Median ITV score per group |
| `wilcoxon_p` | Two-sided Wilcoxon rank-sum p-value |

#### Optional Snakemake trigger (Stage 7)

To call the batch script from within `Snakefile_bootstrap`, add to any condition config:

```yaml
variability_batch_config: config/variability_batch.yaml
variability_batch_outdir: results/variability_combined   # optional
variability_batch_jobs:   4                              # optional
```

The `variability_batch` rule (Stage 7) will then run after Stage 6 completes.

---

### 5. Running all focus genes manually

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

### 6. Restricting to a gene subset

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
| `Couldn't find a sufficient Python binary` (from argparse) | R's `argparse` package calls Python internally; Python not on PATH | Set the `PYTHON` environment variable: `export PYTHON=$(which python)` before running, or `env PYTHON=$(which python) Rscript ...`; in conda: `conda run -n dronet env PYTHON=... Rscript ...` |
| PDF produced but violins are very flat | Too few samples or low variance genes | Expected with small test datasets; use full expression matrices for production |

---

## Related Reading

- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) — Full pipeline from expression data to differential networks; Stage 1 defines the same LOW/HIGH groups used here
- [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) — Config reference for `low_frac`, `high_frac`, and `gene_subset` parameters
- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) — Network topology metrics that complement these variability results
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — Statistical background for the bootstrap pipeline and significance testing

---

**Last Updated:** 2026-03-23
**Scripts:** `src/scripts/15analysis/plot_gene_mad_variability.R`, `src/scripts/15analysis/plot_sample_variability.R`, `src/scripts/15analysis/run_variability_batch.py`
**Status:** ✅ Both scripts implemented; integrated as Stage 6 in `Snakefile_bootstrap`; batch script with Excel output available; partner-type extension is a planned stub
