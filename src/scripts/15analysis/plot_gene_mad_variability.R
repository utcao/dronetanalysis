#!/usr/bin/env Rscript
# ==============================================================================
# Gene-Level MAD Variability: LOW vs HIGH Expression Groups
#
# For a given focus gene, samples are split into LOW (bottom k%) and HIGH
# (top k%) groups by that gene's expression rank. For every other gene in the
# matrix (or an optional partner-gene subset), MAD (Median Absolute Deviation)
# is computed across each group's samples. The two distributions are compared
# via a boxplot and Wilcoxon rank-sum test.
#
# One dot per gene; n_dots ≈ n_genes (or n_partner_genes).
#
# Usage:
#   Rscript plot_gene_mad_variability.R \
#     --expr-file data/processed/VOOM/voomdataCtrl.txt \
#     --focus-gene FBgn0002563 \
#     --output-dir results/variability \
#     --condition-label "Control" \
#     --gene-symbol "l(2)gl" \
#     --save-csv
#
# Future extension (--partner-type direct/indirect):
#   Provide --network-file with a two-column TSV (gene1, gene2) of edges, or a
#   one-column text file of partner gene IDs. Currently a stub.
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(argparse)
})

# ----- 2. Helper functions -----
format_pval <- function(p) {
  if (p < 0.001) return(sprintf("p = %s ***", formatC(p, format = "e", digits = 2)))
  if (p < 0.01)  return(sprintf("p = %.3f **",  p))
  if (p < 0.05)  return(sprintf("p = %.3f *",   p))
  return(sprintf("p = %.3f (ns)", p))
}

get_sig_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# ----- 3. Command-line arguments -----
parser <- ArgumentParser(description = "Gene-level MAD variability: LOW vs HIGH expression groups")
parser$add_argument("--expr-file",       required = TRUE,
                    help = "VOOM expression matrix (tab-separated, genes x samples)")
parser$add_argument("--focus-gene",      required = TRUE,
                    help = "FlyBase ID of the focus gene driving the LOW/HIGH split")
parser$add_argument("--output-dir",      required = TRUE,
                    help = "Directory for output PDF and optional CSV")
parser$add_argument("--low-frac",        type = "double",  default = 0.2,
                    help = "Fraction of samples in LOW group (default: 0.2)")
parser$add_argument("--high-frac",       type = "double",  default = 0.2,
                    help = "Fraction of samples in HIGH group (default: 0.2)")
parser$add_argument("--condition-label", default = "Condition",
                    help = "Human-readable condition label for plot title")
parser$add_argument("--gene-symbol",     default = NULL,
                    help = "Optional gene symbol for plot label (e.g. 'Hsp83')")
parser$add_argument("--gene-subset",     default = NULL,
                    help = "Optional path to a text file of gene IDs (one per line) to restrict MAD computation")
parser$add_argument("--partner-type",    default = "all",
                    choices = c("all", "direct", "indirect"),
                    help = "Gene set to include: all genes, direct (degree-1), or indirect (degree-2) partners (default: all)")
parser$add_argument("--network-file",    default = NULL,
                    help = "Path to partner gene list or edge file — required when --partner-type != all (stub)")
parser$add_argument("--save-csv",        action = "store_true",
                    help = "Also save MAD values as a CSV file")
args <- parser$parse_args()

cat("=== Gene-Level MAD Variability Analysis ===\n")
cat("Expression file:  ", args$expr_file,       "\n")
cat("Focus gene:       ", args$focus_gene,       "\n")
cat("Output directory: ", args$output_dir,       "\n")
cat("LOW fraction:     ", args$low_frac,         "\n")
cat("HIGH fraction:    ", args$high_frac,        "\n")
cat("Condition label:  ", args$condition_label,  "\n")
if (!is.null(args$gene_symbol))
  cat("Gene symbol:      ", args$gene_symbol, "\n")
if (!is.null(args$gene_subset))
  cat("Gene subset file: ", args$gene_subset, "\n")
cat("Partner type:     ", args$partner_type,     "\n\n")

# ----- 4. Validate inputs -----
if (!file.exists(args$expr_file))
  stop("Expression file not found: ", args$expr_file)
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = TRUE)
if (args$low_frac  <= 0 || args$low_frac  >= 1)
  stop("--low-frac must be between 0 and 1 (exclusive)")
if (args$high_frac <= 0 || args$high_frac >= 1)
  stop("--high-frac must be between 0 and 1 (exclusive)")

# ----- 5. Load expression matrix -----
cat("Loading expression matrix...\n")
expr_dt <- fread(args$expr_file, data.table = TRUE)

# Auto-detect gene ID column (first column; accept fly_id, gene_id, or any name)
gene_col <- colnames(expr_dt)[1]
cat("  Gene ID column detected: '", gene_col, "'\n", sep = "")
cat("  Dimensions:", nrow(expr_dt), "genes x", ncol(expr_dt) - 1, "samples\n")

if (!args$focus_gene %in% expr_dt[[gene_col]])
  stop("Focus gene not found in expression matrix: ", args$focus_gene)

# ----- 6. Optional: restrict to gene subset -----
if (!is.null(args$gene_subset)) {
  subset_ids <- readLines(args$gene_subset)
  subset_ids <- trimws(subset_ids[nzchar(subset_ids)])
  # Keep focus gene for group-splitting regardless of subset
  keep_ids   <- unique(c(subset_ids, args$focus_gene))
  expr_dt    <- expr_dt[get(gene_col) %in% keep_ids]
  cat("  After gene-subset filter:", nrow(expr_dt), "genes retained\n")
}

# ----- 7. Partner-type stub -----
# "direct" and "indirect" are reserved for future network-based filtering.
# When activated, this block reads a partner gene list from --network-file
# and further restricts the MAD gene pool to degree-1 or degree-2 neighbors.
if (args$partner_type != "all") {
  if (is.null(args$network_file))
    stop("--network-file is required when --partner-type is '", args$partner_type, "'")
  cat("  WARNING: --partner-type '", args$partner_type,
      "' is not yet implemented. Falling back to 'all'.\n", sep = "")
  args$partner_type <- "all"
}

# ----- 8. Identify LOW and HIGH sample groups -----
# Use sort-based rank selection (consistent with Python np.argpartition in Stage 1):
# take the k_low samples with the lowest and k_high samples with the highest
# expression of the focus gene. No tied-value ambiguity from quantile thresholds.
sample_cols <- setdiff(colnames(expr_dt), gene_col)
n_samples   <- length(sample_cols)

k_low  <- floor(n_samples * args$low_frac)
k_high <- floor(n_samples * args$high_frac)

if (k_low  < 2) stop("LOW group has < 2 samples; increase --low-frac or use a larger expression matrix.")
if (k_high < 2) stop("HIGH group has < 2 samples; increase --high-frac or use a larger expression matrix.")

focus_vec <- as.numeric(expr_dt[get(gene_col) == args$focus_gene, ..sample_cols])
names(focus_vec) <- sample_cols

sorted_asc  <- names(sort(focus_vec))
sorted_desc <- names(sort(focus_vec, decreasing = TRUE))

low_samples  <- sorted_asc[1:k_low]
high_samples <- sorted_desc[1:k_high]

cat("  Total samples:", n_samples, "\n")
cat("  LOW group  (bottom", k_low, "samples): expression range [",
    round(min(focus_vec[low_samples]),  3), ",",
    round(max(focus_vec[low_samples]),  3), "]\n")
cat("  HIGH group (top   ", k_high, "samples): expression range [",
    round(min(focus_vec[high_samples]), 3), ",",
    round(max(focus_vec[high_samples]), 3), "]\n\n")

# ----- 9. Compute MAD per gene for each group -----
cat("Computing MAD...\n")

# Exclude focus gene from the MAD pool (comparing it to itself is uninformative)
other_ids  <- setdiff(expr_dt[[gene_col]], args$focus_gene)
expr_other <- as.matrix(expr_dt[get(gene_col) %in% other_ids, ..sample_cols])
rownames(expr_other) <- expr_dt[get(gene_col) %in% other_ids, get(gene_col)]

mad_low  <- apply(expr_other[, low_samples,  drop = FALSE], 1, mad, na.rm = TRUE)
mad_high <- apply(expr_other[, high_samples, drop = FALSE], 1, mad, na.rm = TRUE)

n_genes <- length(mad_low)
cat("  MAD computed for", n_genes, "genes\n")
cat("  MAD LOW  — median:", round(median(mad_low),  4),
    " | range: [", round(min(mad_low),  4), ",", round(max(mad_low),  4), "]\n")
cat("  MAD HIGH — median:", round(median(mad_high), 4),
    " | range: [", round(min(mad_high), 4), ",", round(max(mad_high), 4), "]\n\n")

# ----- 10. Wilcoxon rank-sum test -----
wt <- wilcox.test(mad_low, mad_high, paired = FALSE)
cat("  Wilcoxon rank-sum test:\n")
cat("  ", format_pval(wt$p.value), "\n\n")

# ----- 11. Assemble long-format data frame -----
plot_df <- data.frame(
  gene      = c(names(mad_low), names(mad_high)),
  mad_value = c(mad_low, mad_high),
  group     = factor(
    c(rep("LOW", length(mad_low)), rep("HIGH", length(mad_high))),
    levels = c("LOW", "HIGH")
  ),
  stringsAsFactors = FALSE
)

# ----- 12. Build plot -----
cat("Building plot...\n")

gene_label <- if (!is.null(args$gene_symbol)) {
  paste0(args$gene_symbol, " (", args$focus_gene, ")")
} else {
  args$focus_gene
}

mean_low  <- mean(mad_low,  na.rm = TRUE)
mean_high <- mean(mad_high, na.rm = TRUE)

title_str    <- paste0("Gene-level MAD Variability\n",
                       "Focus: ", gene_label, "  |  Condition: ", args$condition_label)
subtitle_str <- paste0(format_pval(wt$p.value),
                       "  |  mean LOW = ", round(mean_low, 4),
                       ",  mean HIGH = ", round(mean_high, 4))

# Append partner type to subtitle if not "all"
if (args$partner_type != "all")
  subtitle_str <- paste0(subtitle_str, "  [", args$partner_type, " partners only]")

colors_groups <- c("LOW" = "#2E86AB", "HIGH" = "#A23B72")

# Star annotation y-position: just above the violin max
y_max  <- max(plot_df$mad_value, na.rm = TRUE)
y_star <- y_max * 1.08

p <- ggplot(plot_df, aes(x = group, y = mad_value, fill = group)) +
  geom_violin(alpha = 0.55, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, alpha = 0.85, outlier.alpha = 0.25, linewidth = 0.4) +
  geom_jitter(width = 0.07, size = 0.45, alpha = 0.12, color = "black") +
  annotate("text", x = 1.5, y = y_star,
           label = get_sig_stars(wt$p.value), size = 6, vjust = 0) +
  annotate("segment", x = 1, xend = 2, y = y_max * 1.03, yend = y_max * 1.03,
           linewidth = 0.35, color = "grey40") +
  scale_fill_manual(values = colors_groups) +
  labs(
    title    = title_str,
    subtitle = subtitle_str,
    x        = paste0("Sample group (split by ", gene_label, " expression)"),
    y        = "MAD across samples (VOOM)",
    caption  = paste0("n_genes = ", n_genes,
                      "  |  n_LOW = ", k_low, " samples",
                      "  |  n_HIGH = ", k_high, " samples")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey30"),
    plot.caption    = element_text(hjust = 0.5, size = 8,  color = "grey50"),
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# ----- 13. Save outputs -----
partner_suffix <- if (args$partner_type == "all") "" else paste0("_", args$partner_type)
out_stem <- file.path(args$output_dir,
                      paste0(args$focus_gene, "_mad_variability", partner_suffix))

pdf_path <- paste0(out_stem, ".pdf")
ggsave(pdf_path, p, width = 6, height = 7)
cat("  Saved plot: ", pdf_path, "\n")

if (args$save_csv) {
  csv_path <- paste0(out_stem, ".csv")
  write.csv(plot_df, csv_path, row.names = FALSE)
  cat("  Saved data: ", csv_path, "\n")
}

cat("\nDone.\n")
