# MAD Rank Shift Plot Guide

## Overview

`figures/plot_mad_rank_shift.R` is a downstream visualization script that shows
where a selected group of genes sits within the genome-wide MAD rank (from a
reference condition, CT) and how their position shifts in a treatment condition
(TR, e.g. HS or high-sugar).

This is a **dumbbell plot**:

- **Y-axis** — selected gene labels, ordered by the selection metric
  (descending; the gene with the largest value is at the top)
- **X-axis** — MAD rank from 1 (least variable) to N (most variable)
- **● solid dot** — each gene's CT MAD rank
- **○ hollow dot** — the same gene's TR MAD rank
- **Segment** — connects CT rank to TR rank, showing direction and magnitude of the shift
- **Grey rug** — all background genes at their CT MAD rank (density reference)

The script lives in `figures/` because it is an ad-hoc post-pipeline
visualization tool, not part of the main analysis workflow.

---

## Prerequisites

```r
install.packages(c("data.table", "ggplot2", "argparse", "openxlsx"))
```

Input file: the xlsx output of `compute_mad_variability_ranks.R`, which
contains the columns `mad_ct`, `mad_hs`, `mad_hs_minus_ct`, `gene_id`, and
`SYMBOL`. Any other table with numeric MAD columns works, provided you specify
the correct column names.

---

## Arguments

| Argument | Default | Description |
|---|---|---|
| `-i` / `--input` | required | Input table: xlsx, CSV, or TSV |
| `--ct-mad-col` | `mad_ct` | Column with reference-condition MAD values |
| `--tr-mad-col` | `mad_hs` | Column with treatment MAD values (any condition) |
| `--highlight-col` | required | Column used to rank and select genes |
| `--top-n` | `100` | Number of top genes to highlight |
| `--bottom-n` | `0` | Number of bottom genes to highlight in a contrasting colour |
| `--gene-col` | `gene_id` | Column for gene identifiers |
| `--label-col` | `SYMBOL` | Column for gene labels on Y-axis |
| `--ct-label` | `CT` | Display label for the reference condition |
| `--tr-label` | `TR` | Display label for the treatment condition |
| `--title` | (auto) | Optional plot title |
| `-o` / `--output` | required | Output PDF path |

### `--highlight-col` choices

| Use case | `--highlight-col` value |
|---|---|
| Top MAD increase CT → HS | `mad_hs_minus_ct` |
| Top MAD in CT | `mad_ct` |
| Top rewiring hubs | Column from hub metrics TSV (join first) |

---

## Usage Examples

### Top 100 genes with largest MAD increase in HS

```bash
Rscript src/scripts/15analysis/figures/plot_mad_rank_shift.R \
    --input results/variability/mad_ranks.xlsx \
    --ct-mad-col mad_ct \
    --tr-mad-col mad_hs \
    --highlight-col mad_hs_minus_ct \
    --top-n 100 \
    --tr-label HS \
    --output results/figures/mad_rank_shift_top100_delta_mad.pdf
```

### Top and bottom 50 genes (bidirectional shift)

```bash
Rscript src/scripts/15analysis/figures/plot_mad_rank_shift.R \
    --input results/variability/mad_ranks.xlsx \
    --highlight-col mad_hs_minus_ct \
    --top-n 50 --bottom-n 50 \
    --tr-label HS \
    --output results/figures/mad_rank_shift_top_bottom_50.pdf
```

The top group (increased MAD) is shown in purple-red; the bottom group
(decreased MAD) in blue. Background genes remain grey.

### Using a different treatment condition

Change `--tr-mad-col` and `--tr-label` to point to any other condition's MAD
column and its display name:

```bash
Rscript src/scripts/15analysis/figures/plot_mad_rank_shift.R \
    --input results/variability/mad_ranks.xlsx \
    --ct-mad-col mad_ct \
    --tr-mad-col mad_hsugar \
    --highlight-col mad_hsugar_minus_ct \
    --tr-label "High Sugar" \
    --output results/figures/mad_rank_shift_hsugar.pdf
```

### Visualizing top rewiring hubs

The rewiring hub TSV does not contain MAD columns directly. Join it with the
`compute_mad_variability_ranks.R` xlsx first (in R or using `data.table`
merge), then pass the combined table as `--input` with `--highlight-col` set to
the rewiring metric (e.g. `L2L1_rewire`).

---

## Interpreting the Plot

### Segment direction

| Solid dot (CT) | Hollow dot (TR) | Interpretation |
|---|---|---|
| Left of segment | Right | Gene becomes *more* variable in TR |
| Right of segment | Left | Gene becomes *less* variable in TR |
| Overlapping | | Little change in variability rank |

### Y-axis ordering

Genes are sorted by `--highlight-col` descending, so the gene at the top of the
plot has the highest value of the selection metric. For `mad_hs_minus_ct`, the
top gene gained the most variability; for a rewiring score, the top gene was the
most rewired hub.

### Grey rug

The rug at the bottom shows the density of all background genes across the CT
MAD rank axis. A gap in the rug suggests few genes have intermediate variability
(which would be unusual). The rug is plotted at CT ranks only.

---

## Output Details

- **Format:** PDF (vector, print-quality)
- **Dimensions:** automatically scaled — `height = max(8, min(28, n_selected × 0.13 + 3))` inches; width = 9 inches. For 100 genes this gives approximately 16 × 9 inches.
- **Resolution:** unlimited (vector PDF)

---

## Troubleshooting

| Issue | Fix |
|---|---|
| `Missing columns in input table` | Run `names(readxl::read_xlsx(...))` to list available columns; pass correct `--ct-mad-col` / `--tr-mad-col` |
| Gene labels overlap | Reduce `--top-n` or the plot will scale height automatically; alternatively filter to fewer genes |
| `--highlight-col` has many NAs | Script automatically excludes genes with NA in highlight/MAD columns before selecting top N |
| Label column is all NA | Pass `--label-col gene_id` to use gene IDs instead of symbols |

---

## Related Reading

- [GUIDE-12-15Analysis-Overview.md](GUIDE-12-15Analysis-Overview.md) — Full script inventory and when to use `figures/`
- [GUIDE-11-PCA-Gene-Metrics.md](GUIDE-11-PCA-Gene-Metrics.md) — `compute_mad_variability_ranks.R` (produces the input xlsx)
- [GUIDE-05-Expression-Variability-Analysis.md](GUIDE-05-Expression-Variability-Analysis.md) — Single-gene MAD and ITV comparison plots
- [REFERENCE-02-Script-Naming-Convention.md](REFERENCE-02-Script-Naming-Convention.md) — Script naming convention

---

**Last Updated:** 2026-04-27
**Status:** ✅ Active
