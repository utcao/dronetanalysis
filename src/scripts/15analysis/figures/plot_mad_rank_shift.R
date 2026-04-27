#!/usr/bin/env Rscript
# ==============================================================================
# MAD Rank Shift Plot: Reference vs Treatment Condition
#
# Dumbbell plot showing where a selected gene group (e.g., top rewiring genes
# or top MAD-change genes) sits within the genome-wide MAD rank from a
# reference condition (CT) and how their rank shifts in a treatment (TR).
#
# Visualization:
#   Y-axis  : selected gene labels, ordered by highlight column (descending)
#   X-axis  : MAD rank  (1 = least variable → N = most variable)
#   ●  solid : each gene's CT MAD rank
#   ○ hollow : same gene's TR MAD rank
#   segment  : connects CT rank to TR rank for each gene
#   grey rug : all background genes at their CT MAD rank (bottom of plot)
#
# If --bottom-n > 0, top and bottom groups are shown in contrasting colours.
#
# Usage:
#   Rscript figures/plot_mad_rank_shift.R \
#     --input  results/variability/mad_ranks.xlsx \
#     --ct-mad-col mad_ct \
#     --tr-mad-col mad_hs \
#     --highlight-col mad_hs_minus_ct \
#     --top-n 100 \
#     --tr-label HS \
#     --output results/figures/mad_rank_shift_top100.pdf
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
  library(openxlsx)
})

# ----- 1. CLI arguments -----
parser <- ArgumentParser(
  description = "MAD rank shift dumbbell plot: reference condition vs treatment condition"
)
parser$add_argument("-i", "--input", required = TRUE,
    help = "Input table (xlsx / CSV / TSV) from compute_mad_variability_ranks.R")
parser$add_argument("--ct-mad-col", default = "mad_ct",
    help = "Column with reference-condition MAD values (default: mad_ct)")
parser$add_argument("--tr-mad-col", default = "mad_hs",
    help = "Column with treatment MAD values (default: mad_hs; works for any condition)")
parser$add_argument("--highlight-col", required = TRUE,
    help = "Column used to rank/select genes (e.g., mad_hs_minus_ct, rewiring_score)")
parser$add_argument("--top-n", type = "integer", default = 100,
    help = "Number of top genes to highlight (default: 100)")
parser$add_argument("--bottom-n", type = "integer", default = 0,
    help = "Number of bottom genes to also highlight in a contrasting colour (default: 0)")
parser$add_argument("--gene-col", default = "gene_id",
    help = "Column for gene identifiers (default: gene_id)")
parser$add_argument("--label-col", default = "SYMBOL",
    help = "Column for gene labels shown on Y-axis (default: SYMBOL)")
parser$add_argument("--ct-label", default = "CT",
    help = "Display label for the reference condition (default: CT)")
parser$add_argument("--tr-label", default = "TR",
    help = "Display label for the treatment condition (default: TR)")
parser$add_argument("--title", default = "",
    help = "Optional plot title (auto-generated if omitted)")
parser$add_argument("-o", "--output", required = TRUE,
    help = "Output PDF path")
args <- parser$parse_args()

cat("=== MAD Rank Shift Plot ===\n")
cat("Input:          ", args$input,          "\n")
cat("CT MAD column:  ", args$ct_mad_col,     "\n")
cat("TR MAD column:  ", args$tr_mad_col,     "\n")
cat("Highlight col:  ", args$highlight_col,  "\n")
cat("Top N:          ", args$top_n,          "\n")
cat("Bottom N:       ", args$bottom_n,       "\n")
cat("CT label:       ", args$ct_label,       "\n")
cat("TR label:       ", args$tr_label,       "\n")
cat("Output:         ", args$output,         "\n\n")

# ----- 2. Load input table -----
ext <- tolower(tools::file_ext(args$input))
if (!file.exists(args$input))
  stop("Input file not found: ", args$input)

cat("Loading input table...\n")
if (ext == "xlsx") {
  dt <- as.data.table(read.xlsx(args$input))
} else {
  sep <- if (ext == "csv") "," else "\t"
  dt  <- fread(args$input, sep = sep, data.table = TRUE)
}
cat("  Rows:", nrow(dt), "  Cols:", ncol(dt), "\n")

# ----- 3. Validate columns -----
required_cols <- c(args$gene_col, args$ct_mad_col, args$tr_mad_col, args$highlight_col)
missing_cols  <- setdiff(required_cols, colnames(dt))
if (length(missing_cols) > 0)
  stop("Missing columns in input table: ", paste(missing_cols, collapse = ", "))

# Non-fatal: fall back to gene_col if label_col is absent
if (!args$label_col %in% colnames(dt)) {
  warning("Label column '", args$label_col,
          "' not found — falling back to '", args$gene_col, "'")
  args$label_col <- args$gene_col
}

# ----- 4. Compute condition-specific MAD ranks -----
# Ascending rank: rank 1 = lowest MAD (least variable), rank N = highest MAD (most variable)
n_genes <- nrow(dt)
dt[, ct_mad_rank := rank(get(args$ct_mad_col), ties.method = "average", na.last = "keep")]
dt[, tr_mad_rank := rank(get(args$tr_mad_col), ties.method = "average", na.last = "keep")]
cat("  Genes with CT rank: ", sum(!is.na(dt$ct_mad_rank)), "\n")
cat("  Genes with TR rank: ", sum(!is.na(dt$tr_mad_rank)), "\n\n")

# ----- 5. Select genes -----
# Genes must have non-NA values in both MAD columns and the highlight column
dt_valid <- dt[!is.na(ct_mad_rank) & !is.na(tr_mad_rank) & !is.na(get(args$highlight_col))]

dt_top <- dt_valid[order(-get(args$highlight_col))][seq_len(min(args$top_n, nrow(dt_valid)))]
dt_top[, grp := "top"]

dt_bot <- data.table()
if (args$bottom_n > 0) {
  dt_bot <- dt_valid[order(get(args$highlight_col))][seq_len(min(args$bottom_n, nrow(dt_valid)))]
  # Remove any overlap with top group
  dt_bot <- dt_bot[!get(args$gene_col) %in% dt_top[[args$gene_col]]]
  dt_bot[, grp := "bottom"]
}

selected   <- rbind(dt_top, dt_bot)
background <- dt[!get(args$gene_col) %in% selected[[args$gene_col]]]

n_top <- sum(selected$grp == "top")
n_bot <- if (args$bottom_n > 0) sum(selected$grp == "bottom") else 0
cat("Selected genes:  top =", n_top, "  bottom =", n_bot, "\n")
cat("Background genes:", nrow(background), "\n\n")

# ----- 6. Prepare gene labels for Y-axis -----
selected[, gene_label := as.character(get(args$label_col))]
# Fall back to gene_id for missing symbols
selected[is.na(gene_label) | gene_label == "" | gene_label == "NA",
         gene_label := as.character(get(args$gene_col))]
# Deduplicate: if two genes share a symbol, append the gene_id
dup_labels <- selected$gene_label[duplicated(selected$gene_label)]
if (length(dup_labels) > 0) {
  selected[gene_label %in% dup_labels,
           gene_label := paste0(gene_label, " (", get(args$gene_col), ")")]
}

# Y-axis factor: order ascending by highlight_col so the gene with the
# highest value appears at the TOP of the plot (ggplot plots last level highest)
y_levels <- selected[order(get(args$highlight_col), na.last = TRUE), gene_label]
selected[, gene_label := factor(gene_label, levels = y_levels)]

# ----- 7. Colour palette -----
col_top    <- "#A23B72"   # purple-red for top genes
col_bot    <- "#2E86AB"   # blue for bottom genes
col_bg     <- "grey65"

colour_vals <- if (n_bot > 0) {
  c(top = col_top, bottom = col_bot)
} else {
  c(top = col_top)
}
colour_labs <- if (n_bot > 0) {
  c(top    = paste0("Top ",    n_top, " (high ", args$highlight_col, ")"),
    bottom = paste0("Bottom ", n_bot, " (low ",  args$highlight_col, ")"))
} else {
  c(top = paste0("Top ", n_top, " by ", args$highlight_col))
}

# ----- 8. Axis labels and title -----
x_label <- paste0(
  "MAD rank    (1 = least variable → ", n_genes, " = most variable)\n",
  "● = ", args$ct_label, " rank    ○ = ", args$tr_label, " rank"
)

plot_title <- if (nchar(trimws(args$title)) > 0) {
  args$title
} else {
  paste0("MAD rank shift: ", args$ct_label, " → ", args$tr_label)
}

plot_subtitle <- paste0(
  "Genes selected by: ", args$highlight_col,
  "    |    background n = ", n_genes
)

# ----- 9. Build plot -----
cat("Building plot...\n")

p <- ggplot() +
  # Background rug: all genes at CT MAD rank
  geom_rug(
    data    = background,
    mapping = aes(x = ct_mad_rank),
    colour  = col_bg,
    alpha   = 0.25,
    sides   = "b",
    length  = unit(0.025, "npc"),
    linewidth = 0.3
  ) +
  # Segment: CT rank → TR rank
  geom_segment(
    data    = selected,
    mapping = aes(x     = ct_mad_rank,
                  xend  = tr_mad_rank,
                  y     = gene_label,
                  yend  = gene_label,
                  colour = grp),
    alpha     = 0.55,
    linewidth = 0.7
  ) +
  # Solid dot: reference condition (CT) rank
  geom_point(
    data    = selected,
    mapping = aes(x = ct_mad_rank, y = gene_label, colour = grp),
    shape   = 16,
    size    = 2.5
  ) +
  # Hollow dot: treatment condition (TR) rank
  geom_point(
    data    = selected,
    mapping = aes(x = tr_mad_rank, y = gene_label, colour = grp),
    shape   = 21,
    size    = 2.5,
    fill    = "white",
    stroke  = 0.9
  ) +
  scale_x_continuous(
    name   = x_label,
    limits = c(0, n_genes + 1),
    expand = expansion(mult = 0.01)
  ) +
  scale_y_discrete(name = NULL) +
  scale_colour_manual(values = colour_vals, labels = colour_labs,
                      name = "Selected genes") +
  labs(title = plot_title, subtitle = plot_subtitle) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(hjust = 0, face = "bold", size = 12),
    plot.subtitle    = element_text(hjust = 0, size = 9, colour = "grey40"),
    axis.title.x     = element_text(size = 9, colour = "grey30"),
    axis.text.y      = element_text(size = 7),
    axis.text.x      = element_text(size = 8),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9, face = "bold")
  )

# ----- 10. Save output -----
n_selected <- nrow(selected)
pdf_height <- max(8, min(28, n_selected * 0.13 + 3))
pdf_width  <- 9

out_dir <- dirname(args$output)
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

ggsave(args$output, p, width = pdf_width, height = pdf_height, device = "pdf")
cat("Saved: ", args$output, "\n")
cat("  PDF size: ", pdf_width, "x", round(pdf_height, 1), "inches\n")
cat("\nDone.\n")
