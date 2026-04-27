# Module Analysis Guide

## Overview

Four scripts in `src/scripts/15analysis/` handle module-level analysis. Three
perform statistical and functional analysis of WGCNA modules; one visualizes
the gene co-expression network directly. All four retained their original
filenames (see naming convention note below).

| Script | Purpose |
|---|---|
| `module_enrichment.R` | GO / KEGG functional enrichment per module |
| `module_overlap_analysis.R` | Cross-condition module overlap (Jaccard, Fisher's exact) |
| `module_preservation_netrep.R` | Module preservation via NetRep permutation test |
| `visualise_network.R` | igraph visualization of the TOM-derived co-expression network |

> **Naming note:** These four scripts begin with `module_` or contain
> "network", and were excluded from the 2026-04-27 rename batch to preserve
> continuity with existing downstream references. They do not follow the
> `<verb>_<metric>_<description>` convention but are otherwise unchanged.

---

## Prerequisites

All scripts require R ≥ 4.0. Additional package requirements per script:

| Script | Key packages |
|---|---|
| `module_enrichment.R` | `clusterProfiler`, `org.Dm.eg.db`, `DOSE`, `ggplot2`, `openxlsx` |
| `module_overlap_analysis.R` | `data.table`, `ggplot2`, `pheatmap` |
| `module_preservation_netrep.R` | `NetRep`, `data.table`, `ggplot2` |
| `visualise_network.R` | `igraph`, `data.table`, `ggplot2` |

---

## `module_enrichment.R`

### Purpose

For each WGCNA module, performs GO (BP, MF, CC) and KEGG pathway enrichment
using `clusterProfiler`. Produces per-module text reports, dotplot PDFs, and a
combined summary table.

### Inputs

| Argument | Description |
|---|---|
| `--gene-metrics` | Gene metrics CSV with a `module` column |
| `--output-dir` | Directory for all output files |
| `--p-cutoff` | Initial BH-adjusted p-value cutoff (default: 0.05); relaxed automatically if no results |
| `--simplify` | Flag to reduce redundant GO terms |

### Outputs

```
output-dir/
├── enrichment/
│   ├── module_M1_GO_BP.txt       # per-module text summary
│   ├── module_M1_enrichment.pdf  # dotplot
│   └── ...
├── module_enrichment_combined.csv
└── module_summary_barchart.pdf
```

### Usage

```bash
Rscript src/scripts/15analysis/module_enrichment.R \
    --gene-metrics results/ct/gene_metrics_with_modules.csv \
    --output-dir results/enrichment/ \
    --simplify
```

---

## `module_overlap_analysis.R`

### Purpose

Quantifies how well WGCNA modules are preserved between two conditions by
computing gene-membership overlap statistics.

### Statistics computed

- **Overlap count** — number of shared genes between every CT-module × HS-module pair
- **Jaccard similarity** — |A ∩ B| / |A ∪ B|
- **Overlap coefficient** — |A ∩ B| / min(|A|, |B|)
- **Fisher's exact test** — significance of overlap given background gene set
- **Bonferroni correction** — applied across all module pairs

### Inputs

| Argument | Description |
|---|---|
| `--ct-modules` | CSV with gene_id and module columns for CT |
| `--hs-modules` | CSV with gene_id and module columns for HS |
| `--output-dir` | Directory for all output files |

### Outputs

```
output-dir/
├── overlap_count_heatmap.pdf
├── jaccard_heatmap.pdf
├── overlap_coefficient_heatmap.pdf
├── significance_heatmap.pdf
├── preservation_barplot.pdf
├── overlap_statistics.csv
└── module_best_matches.csv
```

Preservation is categorized: **strong** (Jaccard > 0.5), **moderate** (0.3–0.5),
**weak** (0.1–0.3), **none** (< 0.1).

### Usage

```bash
Rscript src/scripts/15analysis/module_overlap_analysis.R \
    --ct-modules results/ct/module_assignments.csv \
    --hs-modules results/hs/module_assignments.csv \
    --output-dir results/module_overlap/
```

---

## `module_preservation_netrep.R`

### Purpose

Tests module preservation using the NetRep algorithm — a permutation-based
approach that is 10–100× faster and more memory-efficient than WGCNA's
`modulePreservation()` function.

Control-condition module assignments are used as reference; the HS expression
matrix is used to test whether the same gene sets show similar network
structure.

### Inputs

| Argument | Description |
|---|---|
| `--ct-expr` | CT VOOM expression matrix (genes × samples) |
| `--hs-expr` | HS VOOM expression matrix |
| `--ct-modules` | Module assignment CSV (gene_id, module) for CT |
| `--output-dir` | Directory for output |
| `--n-perm` | Number of permutations (default: 1000) |
| `--n-threads` | Parallel threads (default: 4) |
| `--corr-type` | Correlation type: `pearson` or `bicor` (default: `bicor`) |

### Usage

```bash
Rscript src/scripts/15analysis/module_preservation_netrep.R \
    --ct-expr data/processed/VOOM/voomdataCtrl.txt \
    --hs-expr data/processed/VOOM/voomdataHS.txt \
    --ct-modules results/ct/module_assignments.csv \
    --output-dir results/preservation/ \
    --n-perm 1000 \
    --n-threads 8
```

---

## `visualise_network.R`

### Purpose

Loads a Topological Overlap Matrix (TOM) and visualizes the resulting
co-expression network using igraph. Nodes are coloured by module membership.

### Inputs

| Argument | Description |
|---|---|
| `--tom` | TOM matrix CSV (genes × genes) |
| `--gene-metrics` | Gene metrics CSV (for module colors) |
| `--module-assignments` | Module assignment CSV |
| `--threshold` | Minimum TOM edge weight to include (default: 0.1) |
| `--layout` | Graph layout: `fr` (Fruchterman-Reingold), `kk` (Kamada-Kawai), `drl` (DrL) |
| `--hsp-genes` | (Optional) text file of HSP/cochaperone gene IDs to highlight |
| `--output-dir` | Directory for output |

### Usage

```bash
Rscript src/scripts/15analysis/visualise_network.R \
    --tom results/ct/tom_matrix.csv \
    --gene-metrics results/ct/gene_metrics.csv \
    --module-assignments results/ct/module_assignments.csv \
    --threshold 0.15 \
    --layout fr \
    --output-dir results/network_vis/
```

---

## Troubleshooting

| Issue | Fix |
|---|---|
| `module_enrichment.R` finds no results | Relaxed p-value cutoffs are tried automatically; check that gene IDs are valid FlyBase IDs |
| `org.Dm.eg.db` not installed | `BiocManager::install("org.Dm.eg.db")` |
| NetRep crashes with memory error | Reduce `--n-perm` or `--n-threads`; ensure sufficient RAM (≥16 GB for large networks) |
| Network plot too dense | Increase `--threshold` to reduce edge count |

---

## Related Reading

- [GUIDE-12-15Analysis-Overview.md](GUIDE-12-15Analysis-Overview.md) — Full script inventory
- [REFERENCE-02-Script-Naming-Convention.md](REFERENCE-02-Script-Naming-Convention.md) — Naming convention (why these files are unchanged)
- [GUIDE-13-Condition-Comparison-Analysis.md](GUIDE-13-Condition-Comparison-Analysis.md) — Cross-condition comparison
- [GUIDE-07-Pathway-Enrichment.md](GUIDE-07-Pathway-Enrichment.md) — Hub-level (not module-level) enrichment

---

**Last Updated:** 2026-04-27
**Status:** ✅ Active
