#!/usr/bin/env Rscript
# ==============================================================================
# Distribution Plot with Node Highlighting
#
# Plots the density distribution of a chosen numeric column from an input
# table.  Specific entries can be highlighted by exact match against the ID
# column or (when --ann-col is supplied) the annotation column; they appear
# as coloured points on the density curve and are labelled with ggrepel.
#
# Highlight matching rule:
#   An entry is highlighted when its ID column value OR its annotation column
#   value exactly matches any element of --highlight (comma-separated).
#   Labels shown are the annotation column value when --ann-col is supplied,
#   otherwise the ID column value.
#
# Usage:
#   Rscript figures/plot_distribution_highlight.R \
#     --input     run_voomct/results_ct_voom/rewiring_hubs_ct_anno_0408_2026.tsv \
#     --data-col  L2L1_rewire \
#     --id-col    gene_idx \
#     --ann-col   SYMBOL \
#     --highlight "Tps1,Hsp83" \
#     --output    results/figures/dist_L2L1_rewire_tps1.pdf
#
#   # Without annotation column — ID values used as labels:
#   Rscript figures/plot_distribution_highlight.R \
#     --input     data/metrics.tsv \
#     --data-col  score \
#     --highlight "gene_42,gene_7" \
#     --output    results/figures/score_dist.pdf
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(argparse)
  library(openxlsx)
})

# ----- 1. CLI arguments -----
parser <- ArgumentParser(
  description = "Density distribution plot with highlighted entry annotation"
)
parser$add_argument("-i", "--input", required = TRUE,
    help = "Input table (xlsx / CSV / TSV)")
parser$add_argument("--data-col", required = TRUE,
    help = "Numeric column whose distribution is plotted")
parser$add_argument("--id-col", default = NULL,
    help = "ID column for matching highlights (default: first column)")
parser$add_argument("--ann-col", default = NULL,
    help = paste("Optional annotation column; used as display label and also",
                 "matched against --highlight values"))
parser$add_argument("--highlight", default = NULL,
    help = "Comma-separated IDs or annotation values to highlight (exact match)")
parser$add_argument("--title", default = "",
    help = "Optional plot title (auto-generated from --data-col if omitted)")
parser$add_argument("--tag", default = NULL,
    help = "Optional free-text tag appended to the plot title and inserted before the file extension in the output path")
parser$add_argument("-o", "--output", required = TRUE,
    help = "Output PDF path")
args <- parser$parse_args()

cat("=== Distribution Plot with Node Highlighting ===\n")
cat("Input:     ", args$input,    "\n")
cat("Data col:  ", args$data_col, "\n")
cat("ID col:    ", if (is.null(args$id_col)) "(first column)" else args$id_col, "\n")
cat("Ann col:   ", if (is.null(args$ann_col)) "(none)" else args$ann_col, "\n")
cat("Highlight: ", if (is.null(args$highlight)) "(none)" else args$highlight, "\n")
cat("Tag:       ", if (is.null(args$tag)) "(none)" else args$tag, "\n")
cat("Output:    ", args$output,   "\n\n")

# ----- 2. Load input table -----
if (!file.exists(args$input))
  stop("Input file not found: ", args$input)

ext <- tolower(tools::file_ext(args$input))
cat("Loading input table...\n")
if (ext == "xlsx") {
  dt <- as.data.table(read.xlsx(args$input))
} else {
  sep <- if (ext == "csv") "," else "\t"
  dt  <- fread(args$input, sep = sep, data.table = TRUE)
}
cat("  Rows:", nrow(dt), "  Cols:", ncol(dt), "\n")

# ----- 3. Resolve column names -----
id_col <- if (!is.null(args$id_col)) args$id_col else colnames(dt)[1]

required_cols <- c(id_col, args$data_col)
missing_cols  <- setdiff(required_cols, colnames(dt))
if (length(missing_cols) > 0)
  stop("Missing columns in input: ", paste(missing_cols, collapse = ", "))

if (!is.null(args$ann_col) && !args$ann_col %in% colnames(dt))
  stop("Annotation column not found: ", args$ann_col)

# ----- 4. Extract vectors -----
vals    <- as.numeric(dt[[args$data_col]])
ids     <- as.character(dt[[id_col]])
ann_vec <- if (!is.null(args$ann_col)) as.character(dt[[args$ann_col]]) else NULL

n_total <- length(vals)
n_valid <- sum(!is.na(vals))
cat("  Non-NA values:", n_valid, "/", n_total, "\n")

if (n_valid < 2)
  stop("Too few non-NA values in '", args$data_col, "' to compute a density.")

# ----- 5. Resolve highlight index -----
# --- annotation: computed once for all rows ---
is_highlight <- logical(n_total)
if (!is.null(args$highlight)) {
  highlight_targets <- trimws(strsplit(args$highlight, ",")[[1]])
  if (!is.null(ann_vec)) {
    is_highlight <- is_highlight | (ann_vec %in% highlight_targets)
  }
  is_highlight <- is_highlight | (ids %in% highlight_targets)

  cat("  Highlight targets:", paste(highlight_targets, collapse = ", "), "\n")
  cat("  Matched entries:  ", sum(is_highlight), "\n")

  ann_match <- if (!is.null(ann_vec)) highlight_targets %in% ann_vec else rep(FALSE, length(highlight_targets))
  id_match  <- highlight_targets %in% ids
  unmatched <- highlight_targets[!id_match & !ann_match]
  if (length(unmatched) > 0)
    warning("Highlight targets not found in data: ",
            paste(unmatched, collapse = ", "))
}

# Display labels: annotation column if available, otherwise ID
if (!is.null(ann_vec)) {
  display_labels <- ifelse(
    is.na(ann_vec) | ann_vec == "" | ann_vec == "NA",
    ids, ann_vec
  )
} else {
  display_labels <- ids
}

# ----- 6. Compute density and interpolate highlight positions -----
vals_complete  <- vals[!is.na(vals)]
dens           <- density(vals_complete)
ecdf_fn        <- ecdf(vals_complete)

# --- per-highlight: positions on the density curve for geom_point ---
is_hl_complete <- is_highlight & !is.na(vals)
hl_vals        <- vals[is_hl_complete]
hl_labels      <- display_labels[is_hl_complete]
hl_dens        <- approx(dens$x, dens$y, xout = hl_vals)$y
hl_pct         <- (1 - ecdf_fn(hl_vals)) * 100

hl_labels_pct  <- paste0(hl_labels, "\n", sprintf("%.3f%%", hl_pct))

# ----- 7. Assemble plot data frames -----
plot_df <- data.frame(
  val = vals_complete,
  stringsAsFactors = FALSE
)

plot_df_hl <- data.frame(
  val  = hl_vals,
  yval = hl_dens,
  lbl  = hl_labels_pct,
  stringsAsFactors = FALSE
)

# Calculate the 95th percentile (top 5% threshold)
top_5_threshold <- quantile(vals_complete, 0.95)

# ----- 8. Plot title and output path -----
plot_title <- if (nchar(trimws(args$title)) > 0) args$title else args$data_col
if (!is.null(args$tag))
  plot_title <- paste0(plot_title, " — ", args$tag)

if (!is.null(args$tag)) {
  tag_safe <- gsub("[^A-Za-z0-9._-]", "_", args$tag)
  out_path <- paste0(
    tools::file_path_sans_ext(args$output), "_", tag_safe,
    ".", tools::file_ext(args$output)
  )
} else {
  out_path <- args$output
}

# ----- 9. Build plot -----
cat("Building plot...\n")

p <- ggplot(plot_df, aes(x = val)) +
  geom_density(fill = "grey85", color = "grey50", alpha = 0.7, linewidth = 0.6) +
  # Add the vertical dashed line for the top 5%
  geom_vline(xintercept = top_5_threshold, 
             color = "darkred", 
             linetype = "dashed", 
             linewidth = 0.8) +
  # The "Top 5%" text label
  geom_text(aes(x = top_5_threshold, 
                y = max(dens$y) * 0.9, # Places it at 90% of the plot's maximum height
                label = "Top 5%"),
            color = "darkred", 
            hjust = -0.1,            # Nudges text slightly to the right of the line
            vjust = 0,               # Centers vertically on the anchor point
            size = 5) +
  geom_point(data = plot_df_hl,
             aes(x = val, y = yval),
             color = "#D7191C", size = 2.5) +
  geom_text_repel(data = plot_df_hl,
                  aes(x = val*1.2, y = yval, label = lbl),
                  color = "#D7191C",
                  size = 5, max.overlaps = 30) +
  labs(title = plot_title,
       x = args$data_col,
       y = "Density") +
  theme_classic(base_size = 20)

# ----- 10. Save output -----
out_dir <- dirname(args$output)
if (nchar(out_dir) > 0 && !dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

ggsave(args$output, p, width = 9, height = 6, device = "pdf")
cat("Saved:", args$output, "\n\nDone.\n")
