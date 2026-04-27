# Expression Variability Analysis Guide

## Overview

This guide explains two standalone R scripts that quantify how transcriptomic
variability differs between the **LOW** and **HIGH** sample groups defined by a
focus gene's expression rank. These analyses complement the differential
co-expression network results by asking a related but distinct question:

> *When a focus gene is highly expressed, does the rest of the transcriptome
> become more or less variable?*

Three complementary tools are provided:

| Script | Scope | Unit of observation | Variability measure | Output |
|--------|-------|--------------------|--------------------|--------|
| `plot_mad_gene_comparison.R` | Single focus gene | One dot per gene | MAD across group samples | PDF (+ optional CSV) |
| `plot_itv_sample_comparison.R`  | Single focus gene | One dot per sample | ITV — deviation from population mean | PDF (+ optional CSV) |
| `compute_mad_transcriptomic_variability.R` | **All genes at once** | One row per gene | MAD + Wilcoxon + BH-FDR | xlsx table |

The first two scripts are self-contained CLI tools. They derive the LOW/HIGH
sample groups directly from the expression file (no pipeline HDF5 outputs
required) and produce a PDF boxplot with a Wilcoxon rank-sum test annotation.
The third script runs the MAD analysis genome-wide and consolidates all
statistics into a single ranked xlsx table.

---

## Prerequisites

- R ≥ 4.0 with the following packages installed:
  ```r
  install.packages(c("data.table", "dplyr", "ggplot2", "argparse", "openxlsx"))
  ```
  All packages are already present in the `dronet` conda environment.
  `openxlsx` requires no Java dependency.

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
Rscript src/scripts/15analysis/plot_mad_gene_comparison.R \
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

### 2. Script 2 — Genome-wide batch MAD summary (all genes → xlsx)

`compute_mad_transcriptomic_variability.R` runs the same LOW/HIGH MAD analysis as
Script 1, but iterates over **every gene** in the expression matrix. The
expression matrix is loaded once (not once per gene), making this efficient
even for genome-scale data. A gene_id → SYMBOL mapping file is joined to add
readable gene names.

**Usage:**
```bash
Rscript src/scripts/15analysis/compute_mad_transcriptomic_variability.R \
  --expr-file    data/processed/VOOM/voomdataCtrl.txt \
  --mapping-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --output-file  results/variability/all_genes_mad_summary_ct.xlsx \
  --condition-label "Control"
```

**Key arguments:**

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--expr-file` | yes | — | VOOM expression matrix (TSV, genes × samples) |
| `--mapping-file` | yes | — | TSV/CSV with `gene_id` and `SYMBOL` columns |
| `--output-file` | yes | — | Output `.xlsx` path |
| `--low-frac` | no | 0.2 | Fraction of samples in LOW group |
| `--high-frac` | no | 0.2 | Fraction of samples in HIGH group |
| `--condition-label` | no | "Condition" | Added as a constant column in the table |

**Output xlsx columns** (sheet `MAD_variability`, sorted by `wilcoxon_p_adj` ascending):

| Column | Description |
|--------|-------------|
| `gene_id` | FlyBase ID |
| `SYMBOL` | Gene symbol from mapping table (`NA` if absent) |
| `wilcoxon_p` | Raw Wilcoxon rank-sum p-value |
| `wilcoxon_stars` | `***` / `**` / `*` / `ns` based on raw p |
| `wilcoxon_p_adj` | BH-FDR adjusted p-value |
| `wilcoxon_adj_stars` | Stars based on `wilcoxon_p_adj` |
| `mean_mad_low` | Mean MAD across LOW-group samples |
| `mean_mad_high` | Mean MAD across HIGH-group samples |
| `delta_mean_mad` | `mean_mad_high − mean_mad_low` (positive = more variable when HIGH) |
| `median_mad_low` | Median MAD (LOW) |
| `median_mad_high` | Median MAD (HIGH) |
| `n_genes_compared` | Number of other genes in the MAD pool |
| `n_low_samples` | Number of samples in LOW group |
| `n_high_samples` | Number of samples in HIGH group |
| `condition_label` | Value of `--condition-label` |

**Console summary** printed after the xlsx is saved:

```
=== Summary ===
Total genes tested:            N
Significant (FDR < 0.05):      X  (XX.X%)
  ↑ Increased variability (HIGH > LOW):  A  (AA.A% of sig)
  ↓ Decreased variability (LOW > HIGH):  B  (BB.B% of sig)
Results saved to: <output-file>
```

**Spot-check:** pick any gene from the xlsx and run Script 1 on that focus gene.
The `wilcoxon_p`, `mean_mad_low`, and `mean_mad_high` values should match
exactly (BH-FDR correction changes `wilcoxon_p_adj` only).

**Requires:** `openxlsx` R package (no Java dependency; already present in the
`dronet` conda environment). All other packages (`data.table`, `argparse`) are
likewise available.

---

### 3. Script 3 — Sample-level ITV variability

**Basic usage:**
```bash
Rscript src/scripts/15analysis/plot_itv_sample_comparison.R \
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

### 4. Running via Snakemake (Stage 6 — recommended for gene_subset configs)

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

### 5. Batch analysis: both conditions + Excel summary (subset of focus genes)

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

### 6. Running all focus genes manually

The 15 focus genes are listed under `gene_subset:` in
`config/ct_voom_snakemake.yaml` and `config/hs_voom_snakemake.yaml`. To run
both scripts for every focus gene in a condition, use a shell loop:

```bash
# Extract focus genes from YAML (requires yq)
GENES=$(yq '.gene_subset[]' config/ct_voom_snakemake.yaml)

for gene in $GENES; do
  Rscript src/scripts/15analysis/plot_mad_gene_comparison.R \
    --expr-file data/processed/VOOM/voomdataCtrl.txt \
    --focus-gene "$gene" \
    --output-dir results/variability_ct \
    --condition-label "Control" \
    --save-csv

  Rscript src/scripts/15analysis/plot_itv_sample_comparison.R \
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

### 7. Restricting to a gene subset

Both scripts accept `--gene-subset` pointing to a plain-text file with one
gene ID per line. This lets you focus MAD or ITV calculations on a biologically
relevant set (e.g. chaperones, a WGCNA module, known interaction partners):

```bash
# Compute MAD only for chaperone genes
Rscript src/scripts/15analysis/plot_mad_gene_comparison.R \
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
- Script 3: each dot = one sample's ITV score (expect ~k_low + k_high dots, shown larger)

---

## Correlating MAD Results with SNP Data

The expression-based MAD and ITV scores produced by these scripts can be directly correlated with genetic variation from the SNP VCF file. Both data sources share the same sample naming convention (`{lineID}_{wellID}`, e.g. `106_A10`), making the join straightforward.

### Key SNP dataset facts (for correlation planning)

See [`docs/dataset_snp_structure.md`](../../../../docs/dataset_snp_structure.md) for full details. Summary:

| Property | Value |
|---|---|
| File | `data/snp/Dmel_head_hs_ct_Miss80_MAF5_LD8_HWE_1975ind.vcf` |
| SNPs | 413,348 biallelic sites; IDs in `{CHROM}_{POS}` format |
| Samples | 1,975 (same `{lineID}_{wellID}` naming as the expression matrix) |
| Genotype encoding | 0/0 → 0, 0/1 → 1, 1/1 → 2, ./. → missing |
| Chromosomes | 2L, 2R, 3L, 3R, 4, and **23** (= X in PLINK encoding) |

### Gene-level MAD ↔ SNP variation

The `gene_id` in the expression MAD results (FlyBase IDs) can be linked to genomic positions via the dm6 annotation, then matched to SNPs by genomic overlap:

```bash
# Step 1: convert FlyBase gene IDs to coordinates
# (requires dm6 GTF or a pre-built gene BED)
awk '$3=="gene"' dm6.gtf | awk '{print $1"\t"$4"\t"$5"\t"$NF}' > dm6_genes.bed

# Step 2: link SNPs to genes
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' \
    data/snp/Dmel_head_hs_ct_Miss80_MAF5_LD8_HWE_1975ind.vcf > snps.bed
bedtools intersect -a snps.bed -b dm6_genes.bed -wa -wb > snps_with_genes.bed

# Step 3: compute per-gene SNP diversity (mean dosage variance across samples)
# → join with per-gene mean_mad_high/mean_mad_low from the xlsx output
```

Once SNP dosage variance and expression MAD are both aggregated at gene level, a Pearson or Spearman correlation can test whether genes with high allelic diversity also show high transcriptomic variability.

> **Important:** Chromosome `23` in the VCF = X chromosome. When intersecting with dm6 annotation (which uses `chrX` or `X`), recode the VCF chromosome label before `bedtools` or filter it separately.

### Sample-level ITV ↔ individual SNP diversity

The per-sample ITV score from `plot_itv_sample_comparison.R` can be correlated with a per-sample measure of genetic heterozygosity:

```python
import allel, numpy as np

callset = allel.read_vcf(
    'data/snp/Dmel_head_hs_ct_Miss80_MAF5_LD8_HWE_1975ind.vcf',
    fields=['calldata/GT', 'samples']
)
gt = allel.GenotypeArray(callset['calldata/GT'])   # (n_snps, n_samples, 2)
# Count heterozygous sites per sample
het_count = gt.count_het(axis=0)                  # shape: (n_samples,)
samples = callset['samples']
# Join with ITV CSV output by sample name → correlate het_count with ITV
```

> **Caution:** multiple samples share the same `lineID` (inbred line replicates). Samples from the same line should be genetically near-identical — differences in ITV between replicates of the same line reflect **technical noise**, not biological variation. For meaningful genetic–transcriptomic correlations, aggregate by `lineID` (mean ITV per line) before correlating with SNP-based diversity, which should also be aggregated per line.

---

## Partner Restriction (Future)

Script 1 (`plot_mad_gene_comparison.R`) accepts `--partner-type direct` or
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
Rscript src/scripts/15analysis/plot_mad_gene_comparison.R \
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
- [docs/dataset_snp_structure.md](../../../../docs/dataset_snp_structure.md) — Full SNP VCF structure reference, sample naming, genotype encoding, and MAD-relevant fields
- [docs/guide_snp_operations.md](../../../../docs/guide_snp_operations.md) — SNP file operations: subsetting, numeric conversion, sequence reconstruction
- [GUIDE-11-PCA-Gene-Metrics.md](GUIDE-11-PCA-Gene-Metrics.md) — PCA combining condition-specific expression stats (mean, median, MAD, CV² for CT and HS) with L2L1 network metrics; uses outputs of both `compute_mad_transcriptomic_variability.R` and `compute_full_mad_cv2_ranks.R`

---

**Last Updated:** 2026-04-17
**Scripts:** `src/scripts/15analysis/plot_mad_gene_comparison.R`, `src/scripts/15analysis/compute_mad_transcriptomic_variability.R`, `src/scripts/15analysis/plot_itv_sample_comparison.R`, `src/scripts/15analysis/run_variability_batch.py`
**Status:** ✅ All scripts implemented; integrated as Stage 6 in `Snakefile_bootstrap`; genome-wide batch summary (xlsx + BH-FDR) added 2026-04-10; partner-type extension is a planned stub
