# Guide: PCA on L2L1 Network Metrics + Expression Variability

**Category:** Guide — Dimensionality reduction and gene characterisation across network rewiring and expression variability axes

---

## Purpose

This analysis combines two orthogonal views of each gene's behaviour:

1. **Network rewiring structure** — how the gene's local hub-layer topology (L2L1 ratios) changes between low- and high-expression sub-groups, in each condition (CT and HS).
2. **Expression variability** — how variable the gene is across all samples (global MAD and CV²) and within condition-specific sub-groups (delta MAD between HIGH and LOW expression groups).

Two independent PCAs are run, one per condition (CT and HS), so that condition-specific structure is not masked by between-condition variance. A curated set of biologically relevant genes is annotated on all plots.

---

## Feature Matrix (8 columns per PCA)

| Feature | Source | Transformation |
|---|---|---|
| `L2L1_deg` | CT or HS annotation TSV | log1p → z-score |
| `L2L1_rewire` | CT or HS annotation TSV | log1p → z-score |
| `L2L1_conn` | CT or HS annotation TSV | log1p → z-score |
| `HL_conn_L1` | CT or HS annotation TSV | log1p → z-score |
| `HL_conn_L2` | CT or HS annotation TSV | log1p → z-score |
| `delta_mean_mad` | CT or HS variability summary xlsx | z-score only |
| `mad_full` | Full-matrix variability xlsx | z-score only |
| `cv2_full` | Full-matrix variability xlsx | z-score only |

`mad_full` and `cv2_full` come from the combined CT+HS VOOM matrix and are shared between both PCAs as global context.

---

## Preprocessing Rationale

### Why log1p on L2L1 / HL columns?

Network degree and rewiring counts are right-skewed: a small number of hub genes have values orders of magnitude larger than average. Without transformation, the first PC is dominated by these extreme hubs rather than capturing meaningful biological gradients. `log1p` compresses the tail while preserving rank order and keeping zeros at zero.

### Why z-score all features?

The 8 features operate on completely different scales (e.g., `L2L1_rewire` can reach thousands after log1p, `delta_mean_mad` is typically ±0.3). Z-scoring (subtract mean, divide by SD) gives each feature equal weight in the PCA covariance matrix. Without it, high-variance features dominate all PCs regardless of biological relevance.

### Why not use ranks for PCA?

Rank transformation imposes a uniform marginal distribution, which destroys information about *how extreme* a hub is relative to others. A gene at rank 1 could be 100× more extreme than rank 2, or only marginally different — ranks cannot tell you. The log1p + z-score pipeline achieves the same outlier compression but retains this magnitude information in PC scores. Ranks are better suited for non-parametric tests, not PCA.

### Why not mix ranks with raw metrics?

Mixing ranks (uniform distribution) with z-scored raw values (approximately normal) creates a feature matrix with inconsistent marginal distributions. PCA implicitly assumes all features are on a comparable scale; mixing breaks this assumption and makes PC loadings harder to interpret.

---

## Missing Data

**Complete cases only per PCA independently.** A gene is included in the CT PCA if all 8 CT features are non-NA, and likewise for HS. The two gene sets may differ slightly if some genes are annotated in only one condition's network result. The number of genes retained per PCA is reported in the console output.

---

## Script

**Path:** `code/dronetanalysis/src/scripts/15analysis/pca_l2l1_variability.R`

**Required R packages:** `data.table`, `openxlsx`, `argparse`, `ggplot2`, `ggrepel`, `factoextra`, `pheatmap`

Install missing packages (once):
```r
install.packages(c("openxlsx", "argparse", "factoextra", "pheatmap"))
```

### Input Files

| Argument | File |
|---|---|
| `--ct-file` | `run_voomct/results_ct_voom/rewiring_hubs_ct_anno_0408_2026.tsv` |
| `--hs-file` | `run_voomhs/results_hs_voom/rewiring_hubs_hs_anno_0413_2026.tsv` |
| `--ct-var-file` | `results/variability/voomct_all_genes_mad_summary.xlsx` |
| `--hs-var-file` | `results/variability/voomhs_all_genes_mad_summary.xlsx` |
| `--full-var-file` | `results/variability/full_mad_cv2_ranks.xlsx` |
| `--gene-list` | Plain text file, one gene SYMBOL per line |
| `--output-dir` | Destination directory (created if absent) |

The annotation TSVs use `ENSEMBL` (FBgn IDs) as the gene key; the script joins all sources on the `gene_id` column. CT and HS annotations are loaded independently — their gene SYMBOL coverage differs, which is why a union mapping approach was used when building the full variability file (see `compute_full_mad_cv2_ranks.R`).

### Gene List Format

Plain text, one SYMBOL per line, no header. Matching is **case-insensitive**.

```
Hsp83
Trap1
Ago1
Ago2
dicer1
dicer2
```

Genes from the list not found in the dataset are reported in the console but do not cause errors.

### Run Command

```bash
Rscript code/dronetanalysis/src/scripts/15analysis/pca_l2l1_variability.R \
  --ct-file       run_voomct/results_ct_voom/rewiring_hubs_ct_anno_0408_2026.tsv \
  --hs-file       run_voomhs/results_hs_voom/rewiring_hubs_hs_anno_0413_2026.tsv \
  --ct-var-file   results/variability/voomct_all_genes_mad_summary.xlsx \
  --hs-var-file   results/variability/voomhs_all_genes_mad_summary.xlsx \
  --full-var-file results/variability/full_mad_cv2_ranks.xlsx \
  --gene-list     gene_list.txt \
  --output-dir    results/pca_gene_metrics
```

> The script refuses to overwrite existing output files. Remove or rename the output directory before re-running.

---

## Output Files

All outputs land in `--output-dir` with `ct_` or `hs_` prefix.

### Tables (xlsx)

| File | Content |
|---|---|
| `{cond}_feature_matrix.xlsx` | z-scored feature matrix used as PCA input (genes × 8 features) |
| `{cond}_pca_eigenvalues.xlsx` | Per-PC: std dev, variance, % variance, cumulative %, broken-stick threshold |
| `{cond}_pca_loadings.xlsx` | Variable loadings for all 8 PCs |
| `{cond}_pca_scores.xlsx` | Per-gene PC scores + SYMBOL + `is_annotated` flag |

### Plots (PDF)

| File | Description |
|---|---|
| `{cond}_pca_scree.pdf` | Bar chart of % variance per PC. Red line = cumulative variance. Green dashed = broken-stick threshold (PCs above the line are considered significant). |
| `{cond}_pca_biplot_PC1_PC2.pdf` | Biplot: all genes as grey points, feature loading arrows in blue, annotated genes in red with `ggrepel` labels. |
| `{cond}_pca_biplot_PC2_PC3.pdf` | Same for PC2 vs PC3. |
| `{cond}_pca_correlation_circle.pdf` | Variable correlation circle coloured by cos² (quality of representation on the two displayed PCs). |
| `{cond}_pca_loadings_heatmap.pdf` | Heatmap of loadings (features × first 5 PCs) with values printed in each cell. Blue = negative loading, red = positive loading. |
| `{cond}_pca_annotated_scores.pdf` | PC1 vs PC2 scores: background genes in grey, annotated genes highlighted in red with repelled labels. Cleaner view than the biplot for inspecting gene positions. |

---

## Interpreting Results

### Scree plot
- Bars above the green broken-stick line represent PCs that explain more variance than expected by chance.
- If PC1 explains > 70% of variance, inspect the loadings heatmap to check if a single feature dominates — this may indicate a confound rather than biological signal.

### Biplot
- Features (arrows) pointing in the same direction are positively correlated.
- Genes near the tip of an arrow score high on that feature.
- Genes near the centre are average on all features.
- Annotated genes (red) far from the centre are outliers on one or more feature axes.

### Correlation circle
- Variables close to the unit circle are well-represented on these two PCs.
- Variables near the centre are not well captured by the displayed PCs — check subsequent PCs.
- Variables pointing in opposite directions are negatively correlated.

### Loadings heatmap
- Shows which features drive each PC.
- A PC dominated by a single feature usually reflects a technical axis; biologically interesting PCs have mixed contributions from multiple features.

---

## Related Scripts

| Script | Purpose |
|---|---|
| `compute_full_mad_cv2_ranks.R` | Computes per-gene MAD and CV² on the combined CT+HS VOOM matrix |
| `summarize_all_genes_mad_variability.R` | Batch Wilcoxon MAD analysis (LOW vs HIGH sub-groups) for all genes |
| `pca_l2l1_variability.R` | This analysis |
