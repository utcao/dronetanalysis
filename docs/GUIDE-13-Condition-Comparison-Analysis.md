# Condition Comparison Analysis Guide

## Overview

`analyze_conditions.R` is the most comprehensive post-pipeline analysis script
in Stage 15. It loads gene-level network metrics and expression statistics for
two conditions (Control = CT, treatment = High Sugar or HS) and produces a full
suite of comparative visualizations and summary tables.

> **Note:** Previously named `compare_conditions.R`. See
> [REFERENCE-02-Script-Naming-Convention.md](REFERENCE-02-Script-Naming-Convention.md)
> for the rename rationale.

---

## Prerequisites

- R ≥ 4.0 with packages: `data.table`, `dplyr`, `ggplot2`, `argparse`,
  `ggrepel`, `pheatmap`, `factoextra`, `FactoMineR`, `openxlsx`
- Gene metrics CSV files from Stage 13 (`src/scripts/13gene_metrics/`) for both
  CT and HS conditions
- (Optional) expression statistics files and module assignment CSVs

---

## Inputs

| Argument | Description |
|---|---|
| `--ct-metrics` | Gene metrics CSV for control condition (CT) |
| `--hs-metrics` | Gene metrics CSV for treatment condition (HS or other) |
| `--output-dir` | Directory for all output files |
| `--focus-module` | (Optional) module ID to highlight in plots |
| `--ct-label` | Display label for CT condition (default: `CT`) |
| `--hs-label` | Display label for treatment condition (default: `HS`) |

---

## Outputs

All files are written to `--output-dir`.

### Figures (PDFs)

| File | Content |
|---|---|
| `pca_network_metrics.pdf` | PCA on network metrics only |
| `pca_combined.pdf` | PCA combining network metrics + expression stats |
| `delta_barplot.pdf` | Mean change in each metric from CT to HS |
| `violin_expression_metrics.pdf` | Violin/box plots of expression metric distributions with Wilcoxon tests |
| `density_overlay.pdf` | Overlaid density plots for each metric, CT vs HS |
| `centrality_vs_connectivity.pdf` | Scatterplot of centrality vs connectivity per gene |
| `ct_vs_hs_scatter.pdf` | Per-metric CT vs HS scatterplot with diagonal reference |
| `expression_vs_network.pdf` | Expression level vs network metric scatterplots |
| `heatmap_correlation.pdf` | Correlation heatmap across all metrics |
| `heatmap_delta.pdf` | Heatmap of per-gene ΔMetric (HS − CT) |
| `heatmap_module_summary.pdf` | Module-level mean metric summary |

### Tables (CSV / xlsx)

| File | Content |
|---|---|
| `top_increased_connectivity.csv` | Genes with largest increase in connectivity CT → HS |
| `top_decreased_connectivity.csv` | Genes with largest decrease in connectivity CT → HS |
| `top_expression_change.csv` | Genes with largest expression change |
| `condition_summary.xlsx` | Full per-gene delta table and module summaries |

---

## Usage

```bash
Rscript src/scripts/15analysis/analyze_conditions.R \
    --ct-metrics results/ct/gene_metrics.csv \
    --hs-metrics results/hs/gene_metrics.csv \
    --output-dir results/condition_comparison/ \
    --ct-label "Control" \
    --hs-label "High Sugar"
```

With optional module focus:

```bash
Rscript src/scripts/15analysis/analyze_conditions.R \
    --ct-metrics results/ct/gene_metrics.csv \
    --hs-metrics results/hs/gene_metrics.csv \
    --output-dir results/condition_comparison/ \
    --focus-module "M3"
```

---

## Interpreting the Outputs

### Delta barplot

Shows mean ΔMetric (HS − CT) across all genes for each metric. A positive bar
means the average gene increases in that metric from CT to HS. The error bars
are SE across genes.

### PCA

The combined PCA (network + expression) uses log1p transformation on
right-skewed count metrics then z-score scaling. Points are coloured by module
if `--focus-module` is set. Features (arrows) are coloured by category
(CT network, HS network, expression).

### Connectivity vs centrality scatterplot

Useful for identifying outliers: genes with high centrality but low connectivity
have few but critical edges. These may be regulatory hubs that rewire rather
than simply gain connections.

---

## Troubleshooting

| Issue | Fix |
|---|---|
| `Error: metric columns not found` | Check that the CT and HS metric files share the same column names |
| PCA drops too many genes | Genes with NA in any metric are excluded; pre-filter or impute |
| Empty violin plots | Verify that `--focus-module` module ID exists in both condition files |

---

## Related Reading

- [GUIDE-12-15Analysis-Overview.md](GUIDE-12-15Analysis-Overview.md) — Full script inventory
- [REFERENCE-02-Script-Naming-Convention.md](REFERENCE-02-Script-Naming-Convention.md) — Naming convention
- [GUIDE-14-Module-Analysis.md](GUIDE-14-Module-Analysis.md) — Module-level enrichment and overlap
- [GUIDE-11-PCA-Gene-Metrics.md](GUIDE-11-PCA-Gene-Metrics.md) — PCA on L2L1 + variability metrics

---

**Last Updated:** 2026-04-27
**Status:** ✅ Active
