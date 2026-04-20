# Pathway Enrichment Analysis Guide

## Overview

`src/scripts/15analysis/pathway_enrichment_hubs.R` is a standalone R CLI tool that performs
functional enrichment analysis on rewiring hub genes from the differential co-expression pipeline.

Key features:

- Selects **top N, bottom N, or both** genes from the hub metrics table ranked by any numeric column
- Primary sorting columns are the **L2L1 ratio metrics** (`L2L1_rewire`, `L2L1_conn`, `L2L1_deg`)
  which capture how much a gene's extended neighbourhood rewires relative to its direct partners
- Runs **GO enrichment** (BP, MF, CC, or all three) and **KEGG pathway enrichment** via clusterProfiler
- Optional **GSEA** (`gseGO` + `gseKEGG`) using the full ranked gene list
- Applies **cascading p-value cutoffs** so you always get output even for sparsely annotated gene sets
- Accepts an **external background gene list** (`--universe-file`) for a biologically correct universe
- Outputs an Excel workbook (one sheet per ontology), a PDF of dot/bar plots, and a TSV of the
  selected genes for traceability

---

## Prerequisites

R packages (all present in the `dronetanalysis` conda environment):

```
clusterProfiler   enrichplot   org.Dm.eg.db
data.table        dplyr        ggplot2
stringr           argparse     openxlsx
```

Input: the annotated rewiring hub TSV produced by Stage 6 of the pipeline
(`results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv`). The file must contain:
- A gene identifier column (default `gene_id` — FlyBase FBgn IDs, e.g. `FBgn0002563`)
- An `ENTREZID` column (populated by Stage 6 annotation; used directly to skip remapping)
- A numeric sort column (e.g. `L2L1_rewire`)

---

## L2L1 Ratio Metrics — Choosing a Sort Column

These three columns measure how much a gene's **outer-layer (L2) network** differs from its
**direct-partner (L1) network**, capturing whether the gene acts as a broad network organiser:

| Column | Formula | Interpretation |
|--------|---------|----------------|
| `L2L1_rewire` | `full_rewire / L1_rewire` | Ratio of total rewiring to direct-partner rewiring. High → gene's extended neighbourhood rewires disproportionately |
| `L2L1_conn` | `full_conn_diff / L1_conn_diff` | Same but for connectivity change (sum of \|Δr\|). Captures strength of rewiring, not just counts |
| `L2L1_deg` | `L2_n_nodes / L1_n_nodes` | Ratio of outer-layer partner count to direct partner count. Measures topological reach |

> **Recommended default:** `L2L1_rewire` — most directly interpretable as differential hub behaviour.

Other useful sort columns: `L1_rewire`, `n_sig_edges_diff`, `full_rewire`.

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input-file FILE` | Path to annotated rewiring hub TSV |

### Gene selection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gene-col COL` | `gene_id` | Column with gene identifiers. FBgn IDs or ENTREZID accepted |
| `--sort-col COL` | `L2L1_rewire` | Column to rank genes by (descending). Rows with NA are excluded |
| `--n-genes N` | `500` | Number of genes to select |
| `--direction` | `top` | `top` — highest N values; `bottom` — lowest N; `both` — union of top N and bottom N |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output-dir DIR` | `results/enrichment` | Output directory (created if missing) |
| `--output-prefix STR` | `hubs_enrich` | Prefix for all output file names |

### Enrichment

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ontology` | `BP` | GO ontology subset: `BP`, `MF`, `CC`, or `ALL` (runs all three) |
| `--pvalue-cutoff` | `0.05` | Initial enrichment p-value cutoff |
| `--qvalue-cutoff` | `0.2` | FDR q-value cutoff (BH correction) |
| `--simplify-go` | off | If set, removes redundant GO terms via semantic similarity (cutoff 0.7) |
| `--run-gsea` | off | If set, runs `gseGO` + `gseKEGG` on the full background ranked by `--sort-col` |
| `--min-genesetsize` | `10` | Minimum gene set size for any analysis |
| `--max-genesetsize` | `500` | Maximum gene set size for any analysis |
| `--show-category` | `20` | Number of top terms to display in plots |

### Background universe

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--universe-file FILE` | `NULL` | Path to a background gene list (one FBgn per line, or TSV/CSV with a `gene_id` column). If omitted, the universe falls back to all genes in the input hub file |

> **Recommended:** always supply `--universe-file dataset/background_gene_list.tsv` (all expressed genes
> from the voom matrix). Using only the hub file as background under-represents the true universe and
> can inflate enrichment significance.

---

## Output Files

All outputs land in `--output-dir`. File names follow the pattern:

```
{prefix}_{direction}{n}_{sort_col}.{ext}
```

| File | Description |
|------|-------------|
| `*_selected_genes.tsv` | The selected gene rows from the input file with all original columns |
| `*.xlsx` | Excel workbook — one sheet per GO ontology (`GO_BP`, `GO_MF`, `GO_CC`) plus `KEGG` |
| `*_plots.pdf` | Dot plot(s) for GO, bar plot for KEGG (one page per analysis) |
| `*_gsea.xlsx` | GSEA results (sheets `GSEA_GO_*`, `GSEA_KEGG`) — only if `--run-gsea` |
| `*_gsea_plots.pdf` | Ridge plots + GSEA enrichment plots for top terms — only if `--run-gsea` |

### Excel sheet columns

Each sheet contains the standard clusterProfiler result columns:

| Column | Description |
|--------|-------------|
| `ID` | Term or pathway identifier (e.g. GO:0006952, dme04151) |
| `Description` | Term/pathway name (organism suffix stripped for KEGG) |
| `GeneRatio` | k/n — fraction of query genes annotated to this term |
| `BgRatio` | M/N — fraction of background genes annotated to this term |
| `pvalue` | Raw enrichment p-value (hypergeometric test) |
| `p.adjust` | BH-FDR adjusted p-value |
| `qvalue` | q-value (Storey's method) |
| `geneID` | `/`-separated gene symbols (or IDs) in the term |
| `Count` | Number of query genes in this term |

---

## Cascading P-value Cutoffs

GO and KEGG enrichment are attempted at three levels in sequence; the first level that returns
any results is used:

1. `p < --pvalue-cutoff` and `q < --qvalue-cutoff` (strict)
2. `p < 0.20` and `q < 1.0` (relaxed)
3. `p < 1.00` and `q < 1.0` (no cutoff — report all overlapping terms)

> This is especially important for Drosophila KEGG, which has limited pathway coverage. Without
> cascading cutoffs, KEGG often returns zero results even when genuine pathway overlaps exist.

---

## Generating the Background Gene List

The voom expression matrix (`data/processed/VOOM/voomdataCtrl.txt`) contains all expressed genes;
CT and HS share the same gene set. Run this once to create the universe file:

```r
library(data.table)

expr <- fread("data/processed/VOOM/voomdataCtrl.txt", header = TRUE, select = 1L)
setnames(expr, 1, "gene_id")

fwrite(expr[, .(gene_id)],
       "dataset/background_gene_list.tsv",
       sep = "\t")
cat("Saved", nrow(expr), "background genes\n")
```

Save the result to `dataset/background_gene_list.tsv` and pass it via `--universe-file` in every
enrichment run (see examples below).

---

## Step-by-Step Examples

### 1. Top 500 genes by L2L1_rewire — BP enrichment

```bash
Rscript src/scripts/15analysis/pathway_enrichment_hubs.R \
  --input-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --sort-col L2L1_rewire \
  --n-genes 500 \
  --direction top \
  --output-dir results/enrichment \
  --output-prefix ct_hubs \
  --ontology BP \
  --universe-file dataset/background_gene_list.tsv
```

Outputs:
```
results/enrichment/ct_hubs_top500_L2L1_rewire_selected_genes.tsv
results/enrichment/ct_hubs_top500_L2L1_rewire.xlsx          (sheets: GO_BP, KEGG)
results/enrichment/ct_hubs_top500_L2L1_rewire_plots.pdf
```

---

### 2. All GO ontologies + simplified terms

```bash
Rscript src/scripts/15analysis/pathway_enrichment_hubs.R \
  --input-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --sort-col L2L1_rewire \
  --n-genes 500 \
  --direction top \
  --output-dir results/enrichment \
  --output-prefix ct_hubs_all_go \
  --ontology ALL \
  --simplify-go \
  --universe-file dataset/background_gene_list.tsv
```

Sheets in Excel: `GO_BP`, `GO_MF`, `GO_CC`, `KEGG`.
Redundant GO terms are collapsed using semantic similarity (GOSemSim, cutoff 0.7).

---

### 3. Bottom 200 genes by L2L1_conn — connectivity-driven selection

```bash
Rscript src/scripts/15analysis/pathway_enrichment_hubs.R \
  --input-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --sort-col L2L1_conn \
  --n-genes 200 \
  --direction bottom \
  --output-dir results/enrichment \
  --output-prefix ct_hubs_low_conn \
  --ontology BP \
  --universe-file dataset/background_gene_list.tsv
```

---

### 4. Both ends (union of top + bottom 300)

Useful to compare whether the functional composition of high-rewiring vs low-rewiring extremes differs.
Inspect the `_selected_genes.tsv` to see which end each gene came from.

```bash
Rscript src/scripts/15analysis/pathway_enrichment_hubs.R \
  --input-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --sort-col L2L1_rewire \
  --n-genes 300 \
  --direction both \
  --output-dir results/enrichment \
  --output-prefix ct_hubs_extremes \
  --ontology BP \
  --universe-file dataset/background_gene_list.tsv
```

---

### 5. Full run with GSEA

GSEA uses all background genes ranked by `--sort-col` (not just the selected top/bottom subset).

```bash
Rscript src/scripts/15analysis/pathway_enrichment_hubs.R \
  --input-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --sort-col L2L1_rewire \
  --n-genes 500 \
  --direction top \
  --output-dir results/enrichment \
  --output-prefix ct_hubs_gsea \
  --ontology BP \
  --run-gsea \
  --universe-file dataset/background_gene_list.tsv
```

Additional outputs:
```
results/enrichment/ct_hubs_gsea_top500_L2L1_rewire_gsea.xlsx
results/enrichment/ct_hubs_gsea_top500_L2L1_rewire_gsea_plots.pdf
```

---

## Statistical Methods

### Over-Representation Analysis (GO and KEGG)

Uses a **hypergeometric test** (Fisher's exact test for 2×2 tables) as implemented in
`clusterProfiler::enrichGO()` and `clusterProfiler::enrichKEGG()`.

- **Query set:** selected top/bottom N genes with valid Entrez IDs
- **Background universe:** if `--universe-file` is supplied, all genes in that file (recommended:
  all expressed genes from the voom matrix); otherwise falls back to all genes in the input hub file.
  Using the full expression universe is more accurate than restricting to hub-table genes.
- **Multiple testing:** BH-FDR correction within each ontology separately
- **GO simplification:** if `--simplify-go`, redundant terms are removed using GOSemSim semantic
  similarity (`clusterProfiler::simplify()`, Wang method, cutoff 0.7, keeping the term with
  lowest adjusted p-value per cluster)

### Gene Set Enrichment Analysis (GSEA)

Uses a **Kolmogorov-Smirnov-like running-sum statistic** as implemented in
`clusterProfiler::gseGO()` and `clusterProfiler::gseKEGG()`.

- **Ranked list:** all background genes sorted descending by `--sort-col` value
- **Interpretation:** positive enrichment score → gene set members cluster near the top of the
  ranked list (high `sort-col` value); negative → cluster near the bottom
- **Multiple testing:** BH-FDR within each analysis

---

## Gene ID Mapping

The script detects the gene ID format automatically:

1. **ENTREZID column present** (standard for Stage 6 output): used directly — no remapping needed
2. **FBgn IDs in `--gene-col`**: mapped to ENTREZID via `AnnotationDbi::mapIds(org.Dm.eg.db, keytype = "FLYBASE")`
3. **Numeric IDs in `--gene-col`**: treated as ENTREZID directly

One-to-many mappings keep the first match. Unmapped genes are silently excluded from both query and universe.

---

## Troubleshooting

| Issue | Cause | Fix |
|-------|-------|-----|
| `Columns not found: L2L1_rewire` | Sort column absent from input file | Check `colnames(fread(file))` — use `L2L1_rewire`, `L2L1_conn`, or `L2L1_deg` |
| `No query genes could be mapped` | All ENTREZID values are NA | Run Stage 6 annotation (`06_annotate_rewiring_table.R`) to populate ENTREZID |
| `Universe file not found` | Wrong path to `--universe-file` | Generate it first (see "Generating the Background Gene List" above) |
| Zero KEGG results (even relaxed) | Limited Drosophila KEGG annotation | Expected — only ~80 pathways for *D. melanogaster*; try `--ontology ALL` for GO instead |
| `simplify()` fails | Python not found in environment | Set `Sys.setenv(PYTHON = ...)` at the top of the script, or skip `--simplify-go` |
| Very large GSEA output | Many significant terms | Reduce `--max-genesetsize` or increase `--pvalue-cutoff` |
| Excel write fails on Windows with `Permission denied` | Output file open in Excel | Close the xlsx file in Excel, then rerun |
| `Error in pdf(...): cannot open file` | Output directory doesn't exist | Script creates it automatically — check write permissions |

---

## Related Reading

- [GUIDE-04-Qualitative-Change-Metrics.md](GUIDE-04-Qualitative-Change-Metrics.md) — full
  definitions of the L1/L2 rewiring metrics and the L2L1 ratio columns used for gene selection
- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) — topological context for
  degree, connectivity and rewiring measures
- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) — how the hub TSV is generated
  (Stage 3b → Stage 6 annotation)
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — bootstrap and
  differential testing methods used upstream of enrichment

---

**Last Updated:** 2026-04-20
**Script:** `src/scripts/15analysis/pathway_enrichment_hubs.R`
**Status:** ✅ Implemented; run against `rewiring_hubs_ct_anno_0408_2026.tsv`
