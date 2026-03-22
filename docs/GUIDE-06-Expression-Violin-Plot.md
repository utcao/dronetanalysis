# Expression Violin/Boxplot Guide

## Overview

`src/scripts/plot_expr_violin.py` is a standalone Python CLI tool that produces
publication-quality **violin + boxplot** figures for gene expression data.

Key features:

- Works with any number of expression matrices (one violin per group)
- Selects **top X% and bottom X%** individuals per group by expression rank
- Overlays all dots by default; optionally **highlights** top/bottom subsets in colour
- Runs pairwise **Wilcoxon rank-sum** tests (= Mann-Whitney U) with **BH-FDR** correction
- Exports **PDF + SVG + PNG** and a companion **`_stats.csv`** into `results/expr_violin/`

---

## Prerequisites

Python packages (all present in the `dronet` conda environment):

```
numpy  pandas  matplotlib  seaborn
```

> `scipy` and `statsmodels` are **not required** — the Mann-Whitney U test and
> BH correction are implemented in pure NumPy/math to avoid compiled-binary
> import failures on Windows.

Input CSV format: **rows = genes** (FlyBase IDs), **columns = individuals**,
values = logCPM (or any continuous expression measure).

```
dataset/raw/logCPM_Ctrl_Dros.csv   # 8 763 genes × 938 individuals
dataset/raw/logCPM_HS_Dros.csv     # 8 763 genes × 1 037 individuals
```

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--files FILE [FILE …]` | One or more input CSV files |
| `--groups NAME [NAME …]` | Group label for each file (same order) |

### Individual selection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--individual_sel` | `1.0` | `1` → use all individuals; `0 < x < 1` → select top x% **and** bottom x% per file, ranked by expression of `--gene` (or mean across all genes) |

### X-axis behaviour

| Condition | X-axis categories |
|-----------|-------------------|
| Multiple files | Group names (e.g. CT, HS) |
| Single file + `individual_sel < 1` | `Bottom X%` and `Top X%` |
| Single file + `individual_sel = 1` | Single violin (no comparison) |

### Dot highlighting

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--highlight_sel` | `0` | `0` or `1` → no highlighting; `0 < x < 1` → within the selected individuals, highlight top x% (red) and bottom x% (blue) |

Highlighting is applied **after** individual selection.

### Gene target

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gene GENEID` | *(none)* | FlyBase ID of the gene to plot. If omitted, per-individual mean across all genes is used for both ranking and the y-axis |
| `--gene_name NAME` | *(none)* | Human-readable gene symbol (e.g. `Hsp70Ab`). When provided, `{GENEID}_{NAME}` is automatically appended to every output file name, making files self-describing without editing `--output_prefix` |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output_prefix` | `expr_violin` | File name stem (no extension) |
| `--output_dir` | `results/expr_violin/` | Directory for all outputs |
| `--dpi` | `150` | PNG resolution |

### Plot appearance

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--figsize W H` | auto | Figure size in inches |
| `--violin_alpha` | `0.55` | Violin fill transparency |
| `--dot_alpha` | `0.35` | Non-highlighted dot transparency |
| `--dot_size` | `8` | Non-highlighted dot size (pt²) |
| `--no_violin` | off | Show box plot only |
| `--no_box` | off | Show violin only |

---

## Dot colour scheme

| Dot type | Colour |
|----------|--------|
| All-selected (when `individual_sel = 1`), not highlighted | light gray |
| Top-selected, not highlighted | salmon (#f4a582) |
| Bottom-selected, not highlighted | steel blue (#92c5de) |
| Top-highlighted | red |
| Bottom-highlighted | royal blue |

---

## Step-by-Step Examples

### 1. All individuals — two-group comparison

```bash
python src/scripts/plot_expr_violin.py \
    --files dataset/raw/logCPM_Ctrl_Dros.csv dataset/raw/logCPM_HS_Dros.csv \
    --groups CT HS \
    --individual_sel 1 \
    --output_prefix CT_vs_HS_all
```

Output: `results/expr_violin/CT_vs_HS_all.{pdf,svg,png}` + `_stats.csv`

---

### 2. Top/bottom 10% selection — two-group comparison, with highlighting

Selects the 10% highest- and 10% lowest-expressing individuals from each file,
then highlights the extreme 5% within that selected pool.

```bash
python src/scripts/plot_expr_violin.py \
    --files dataset/raw/logCPM_Ctrl_Dros.csv dataset/raw/logCPM_HS_Dros.csv \
    --groups CT HS \
    --individual_sel 0.1 \
    --highlight_sel 0.05 \
    --output_prefix CT_vs_HS_top10pct
```

CT violin: 94 top + 94 bottom = 188 individuals
HS violin: 104 top + 104 bottom = 208 individuals

---

### 3. Single gene, readable name in output files

```bash
python src/scripts/plot_expr_violin.py \
    --files dataset/raw/logCPM_Ctrl_Dros.csv dataset/raw/logCPM_HS_Dros.csv \
    --groups CT HS \
    --individual_sel 1 \
    --highlight_sel 0.2 \
    --gene FBgn0039562 \
    --gene_name Hsp70Ab \
    --output_prefix CT_vs_HS
```

Output files: `results/expr_violin/CT_vs_HS_FBgn0039562_Hsp70Ab.{pdf,svg,png}`
Title and y-axis label: `FBgn0039562 (Hsp70Ab)`

---

### 4. Single file — top vs bottom on x-axis

```bash
python src/scripts/plot_expr_violin.py \
    --files dataset/raw/logCPM_Ctrl_Dros.csv \
    --groups CT \
    --individual_sel 0.1 \
    --highlight_sel 0.3 \
    --gene FBgn0039562 \
    --gene_name Hsp70Ab \
    --output_prefix CT_topvsbottom
```

X-axis: `Bottom 10%` (n=94) vs `Top 10%` (n=94)
Output: `results/expr_violin/CT_topvsbottom_FBgn0039562_Hsp70Ab.{pdf,svg,png}`

---

## Output Files

All outputs land in `results/expr_violin/` (or `--output_dir`):

| File | Description |
|------|-------------|
| `{prefix}.pdf` | Vector PDF (publication quality) |
| `{prefix}.svg` | Scalable SVG (editable in Inkscape/Illustrator) |
| `{prefix}.png` | Raster PNG (default 150 dpi) |
| `{prefix}_stats.csv` | Pairwise test results (raw p, adj p BH-FDR, group n) |

File name structure when `--gene_name` is provided:

```
{output_prefix}_{gene_id}_{gene_name}.pdf
```

---

## Statistical Method

Pairwise **Wilcoxon rank-sum test** (two-sided), also known as the Mann-Whitney U
test — the two names refer to the same test. For independent groups of
different individuals, this is the appropriate non-parametric comparison.

> Note: this is *not* the Wilcoxon signed-rank test, which requires paired samples.

Multiple comparison correction: **Benjamini-Hochberg FDR** (BH). Exact adjusted
p-values are displayed on the figure and saved in `_stats.csv`.

Both tests are implemented in pure NumPy/math (no scipy dependency) using:

- `_rankdata()` — average-rank with tie handling
- `mannwhitneyu()` — normal approximation of U statistic
- `bh_correction()` — step-up BH procedure

---

## Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `ModuleNotFoundError: seaborn` | seaborn not installed | `pip install seaborn` |
| `ImportError: cannot import _spropack` | scipy DLL broken (Windows base env) | Script avoids scipy entirely — ensure you are not importing it elsewhere |
| `ERROR: gene 'X' not found` | Gene ID absent from index | Check the exact row label with `head -2 file.csv` |
| `ERROR: --files and --groups must have the same number` | Mismatched argument counts | Provide one `--groups` name per `--files` entry |
| Bimodal violins when `individual_sel < 1` and multiple groups | Expected — top+bottom selected individuals pooled per violin | Intentional: use dot colours (salmon/blue) to distinguish subsets visually |

---

## Related Reading

- [GUIDE-05-Expression-Variability-Analysis.md](GUIDE-05-Expression-Variability-Analysis.md) — complementary R scripts for gene-level MAD and sample-level ITV variability across LOW/HIGH groups
- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) — full pipeline context for how expression data is generated
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) — deeper discussion of statistical methods used across the project

---

**Last Updated:** 2026-03-22
**Script:** `src/scripts/plot_expr_violin.py`
**Status:** ✅ Tested on CT (n=938) and HS (n=1037) expression matrices
