# 15Analysis Scripts Overview

## Overview

The `src/scripts/15analysis/` folder contains the final-stage standalone
analysis and visualization scripts for the dronet project. Unlike the numbered
pipeline stages (00–14), these scripts are run **manually and selectively**
after the main pipeline completes. They consume the pipeline's outputs (gene
metrics TSVs, expression matrices, network results) and produce figures and
summary tables.

Scripts follow the `<verb>_<metric/method>_<description>` naming convention —
see [REFERENCE-02-Script-Naming-Convention.md](REFERENCE-02-Script-Naming-Convention.md)
for the full rationale and rename history.

---

## Script Inventory

### Condition Comparison

| Script | Primary output | Guide |
|---|---|---|
| `analyze_conditions.R` | 10+ PDFs (PCA, violin, heatmaps, scatterplots), multiple CSV summary tables | [GUIDE-13](GUIDE-13-Condition-Comparison-Analysis.md) |

Compares two conditions (CT vs High Sugar) across all gene-level network metrics
and expression statistics. The single most comprehensive post-pipeline analysis
script.

---

### MAD / Expression Variability

| Script | Primary output | Guide |
|---|---|---|
| `compute_mad_variability_ranks.R` | `mad_ranks.xlsx` — per-gene MAD, CV², ranks for CT, HS, and full matrix | [GUIDE-11](GUIDE-11-PCA-Gene-Metrics.md) |
| `compute_mad_transcriptomic_variability.R` | `mad_variability.xlsx` — one row per gene: Wilcoxon stats, mean MAD (LOW/HIGH), BH-FDR | [GUIDE-05](GUIDE-05-Expression-Variability-Analysis.md) |
| `plot_mad_gene_comparison.R` | PDF boxplot — MAD distribution across LOW vs HIGH groups for a single focus gene | [GUIDE-05](GUIDE-05-Expression-Variability-Analysis.md) |
| `plot_itv_sample_comparison.R` | PDF boxplot — ITV per sample across LOW vs HIGH groups for a single focus gene | [GUIDE-05](GUIDE-05-Expression-Variability-Analysis.md) |
| `run_variability_batch.py` | Batch runner — calls `plot_mad_gene_comparison.R` across many focus genes | [GUIDE-05](GUIDE-05-Expression-Variability-Analysis.md) |

---

### Hub / Pathway Enrichment

| Script | Primary output | Guide |
|---|---|---|
| `analyze_hub_enrichment.R` | `.xlsx` with GO/KEGG/GSEA results, PDF enrichment dotplots | [GUIDE-07](GUIDE-07-Pathway-Enrichment.md) |

Selects top/bottom N hub genes ranked by any numeric column (default:
`L2L1_rewire`) and performs GO, KEGG, and optional GSEA enrichment for
*Drosophila*.

---

### PCA

| Script | Primary output | Guide |
|---|---|---|
| `plot_pca_l2l1_variability.R` | PDF static biplots + HTML interactive biplots for CT, HS, and merged PCAs | [GUIDE-11](GUIDE-11-PCA-Gene-Metrics.md) |

Runs PCA on L2L1 network rewiring metrics combined with condition-specific
expression variability statistics (MAD, CV²).

---

### Permutation Test Plots

| Script | Primary output | Guide |
|---|---|---|
| `plot_permutation_null_dist.R` | Multi-page PDF — one histogram per gene showing null distribution + observed value | [GUIDE-08](GUIDE-08-Permutation-Test.md) |

---

### Expression Violin Plot

| Script | Primary output | Guide |
|---|---|---|
| `plot_expression_violin.py` | PDF/SVG/PNG violin+boxplot for individual gene expression across samples | [GUIDE-06](GUIDE-06-Expression-Violin-Plot.md) |

---

### Module Analysis (unchanged filenames)

| Script | Primary output | Guide |
|---|---|---|
| `module_enrichment.R` | Per-module GO/KEGG enrichment PDFs and CSV tables | [GUIDE-14](GUIDE-14-Module-Analysis.md) |
| `module_overlap_analysis.R` | Overlap heatmaps, preservation barplot, statistics CSV | [GUIDE-14](GUIDE-14-Module-Analysis.md) |
| `module_preservation_netrep.R` | Module preservation statistics and plots (NetRep) | [GUIDE-14](GUIDE-14-Module-Analysis.md) |
| `visualise_network.R` | Network graph PDFs/PNGs (igraph, multiple layout algorithms) | [GUIDE-14](GUIDE-14-Module-Analysis.md) |

---

### `figures/` — Ad-Hoc Post-Pipeline Visualizations

Scripts in `figures/` are not part of any main workflow. They are run manually
to generate exploratory or publication-ready figures from existing result tables.

| Script | Primary output | Guide |
|---|---|---|
| `figures/plot_mad_rank_shift.R` | PDF dumbbell plot — genome-wide MAD rank with selected genes highlighted | [GUIDE-15](GUIDE-15-MAD-Rank-Shift-Plot.md) |

---

## Typical Workflow

After the main pipeline completes and gene metrics TSVs are available:

```
1. compute_mad_variability_ranks.R   → mad_ranks.xlsx (genome-wide MAD/CV² ranks)
2. plot_pca_l2l1_variability.R            → PCA biplots on network + variability features
3. analyze_conditions.R              → Cross-condition comparison (CT vs HS)
4. analyze_hub_enrichment.R          → GO/KEGG enrichment on top rewiring hubs
5. compute_mad_transcriptomic_variability.R  → Genome-wide MAD batch summary
6. module_enrichment.R               → Per-module functional enrichment
7. module_overlap_analysis.R         → Cross-condition module overlap
8. figures/plot_mad_rank_shift.R     → Rank shift visualization for selected genes
```

Scripts 1–7 are typically run in order. Scripts in `figures/` are run on demand
as needed for exploration or manuscript figures.

---

## Related Reading

- [REFERENCE-02-Script-Naming-Convention.md](REFERENCE-02-Script-Naming-Convention.md) — Naming convention and rename history
- [GUIDE-05-Expression-Variability-Analysis.md](GUIDE-05-Expression-Variability-Analysis.md) — MAD and ITV analysis scripts
- [GUIDE-11-PCA-Gene-Metrics.md](GUIDE-11-PCA-Gene-Metrics.md) — PCA and MAD variability ranks
- [GUIDE-13-Condition-Comparison-Analysis.md](GUIDE-13-Condition-Comparison-Analysis.md) — Condition comparison
- [GUIDE-14-Module-Analysis.md](GUIDE-14-Module-Analysis.md) — Module enrichment, overlap, preservation, network vis
- [GUIDE-15-MAD-Rank-Shift-Plot.md](GUIDE-15-MAD-Rank-Shift-Plot.md) — Downstream rank shift figure

---

**Last Updated:** 2026-04-27
**Status:** ✅ Active
