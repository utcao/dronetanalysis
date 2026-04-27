#!/usr/bin/env Rscript
# ==============================================================================
# Sample-Level Individual Transcriptomic Variability (ITV): LOW vs HIGH Groups
#
# For a given focus gene, samples are split into LOW (bottom k%) and HIGH
# (top k%) groups by that gene's expression rank. For each sample, an ITV score
# is computed as the summary (median/mean/sum) of absolute deviations from the
# gene-wise population mean across all genes (log-scale subtraction ≡ log(A/mean)).
# The ITV distributions of both groups are compared via a boxplot and Wilcoxon
# rank-sum test.
#
# One dot per sample; n_dots = n_LOW + n_HIGH samples.
#
# Usage:
#   Rscript plot_sample_variability.R \
#     --expr-file data/processed/VOOM/voomdataCtrl.txt \
#     --focus-gene FBgn0002563 \
#     --output-dir results/variability \
#     --condition-label "Control" \
#     --summary-metric median \
#     --save-csv
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
parser <- ArgumentParser(description = "Sample-level ITV: LOW vs HIGH expression groups")
parser$add_argument("--expr-file",       required = TRUE,
                    help = "VOOM expression matrix (tab-separated, genes x samples)")
parser$add_argument("--focus-gene",      required = TRUE,
                    help = "FlyBase ID of the focus gene driving the LOW/HIGH split")
parser$add_argument("--output-dir",      required = TRUE,
                    help = "Directory for output PDF and optional CSV")
parser$add_argument("--low-frac",        type = "double", default = 0.2,
                    help = "Fraction of samples in LOW group (default: 0.2)")
parser$add_argument("--high-frac",       type = "double", default = 0.2,
                    help = "Fraction of samples in HIGH group (default: 0.2)")
parser$add_argument("--condition-label", default = "Condition",
                    help = "Human-readable condition label for plot title")
parser$add_argument("--gene-symbol",     default = NULL,
                    help = "Optional gene symbol for plot label (e.g. 'Hsp83')")
parser$add_argument("--summary-metric",  default = "median",
                    choices = c("median", "mean", "sum"),
                    help = "How to collapse per-gene deviations into one ITV value per sample (default: median)")
parser$add_argument("--gene-subset",     default = NULL,
                    help = "Optional path to a text file of gene IDs (one per line) to restrict ITV computation")
parser$add_argument("--save-csv",        action = "store_true",
                    help = "Also save ITV values as a CSV file")
args <- parser$parse_args()

cat("=== Sample-Level Individual Transcriptomic Variability (ITV) ===\n")
cat("Expression file:  ", args$expr_file,       "\n")
cat("Focus gene:       ", args$focus_gene,       "\n")
cat("Output directory: ", args$output_dir,       "\n")
cat("LOW fraction:     ", args$low_frac,         "\n")
cat("HIGH fraction:    ", args$high_frac,        "\n")
cat("Condition label:  ", args$condition_label,  "\n")
cat("Summary metric:   ", args$summary_metric,   "\n")
if (!is.null(args$gene_symbol))
  cat("Gene symbol:      ", args$gene_symbol, "\n")
if (!is.null(args$gene_subset))
  cat("Gene subset file: ", args$gene_subset, "\n")
cat("\n")

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
expr_dt  <- fread(args$expr_file, data.table = TRUE)

# Auto-detect gene ID column (first column; accept fly_id, gene_id, or any name)
gene_col <- colnames(expr_dt)[1]
cat("  Gene ID column detected: '", gene_col, "'\n", sep = "")
cat("  Dimensions:", nrow(expr_dt), "genes x", ncol(expr_dt) - 1, "samples\n")

if (!args$focus_gene %in% expr_dt[[gene_col]])
  stop("Focus gene not found in expression matrix: ", args$focus_gene)

# All sample column names (used for group splitting AND reference mean)
all_sample_cols <- setdiff(colnames(expr_dt), gene_col)
n_samples       <- length(all_sample_cols)

# ----- 6. Optional: restrict gene pool for ITV computation -----
# The focus gene still drives sample group assignment using the full matrix.
# The gene subset only affects which genes contribute to the ITV calculation.
if (!is.null(args$gene_subset)) {
  subset_ids  <- readLines(args$gene_subset)
  subset_ids  <- trimws(subset_ids[nzchar(subset_ids)])
  itv_gene_ids <- intersect(subset_ids, expr_dt[[gene_col]])
  if (length(itv_gene_ids) == 0)
    stop("No genes from --gene-subset were found in the expression matrix.")
  cat("  Gene subset for ITV:", length(itv_gene_ids), "genes retained\n")
} else {
  itv_gene_ids <- expr_dt[[gene_col]]
}

# ----- 7. Identify LOW and HIGH sample groups -----
# Sort-based rank selection (consistent with Python np.argpartition in Stage 1):
# bottom k_low and top k_high samples by focus gene expression.
k_low  <- floor(n_samples * args$low_frac)
k_high <- floor(n_samples * args$high_frac)

if (k_low  < 2) stop("LOW group has < 2 samples; increase --low-frac or use a larger expression matrix.")
if (k_high < 2) stop("HIGH group has < 2 samples; increase --high-frac or use a larger expression matrix.")

focus_vec   <- as.numeric(expr_dt[get(gene_col) == args$focus_gene, ..all_sample_cols])
names(focus_vec) <- all_sample_cols

low_samples  <- names(sort(focus_vec))[1:k_low]
high_samples <- names(sort(focus_vec, decreasing = TRUE))[1:k_high]

cat("  Total samples:", n_samples, "\n")
cat("  LOW group  (bottom", k_low, "samples): expression range [",
    round(min(focus_vec[low_samples]),  3), ",",
    round(max(focus_vec[low_samples]),  3), "]\n")
cat("  HIGH group (top   ", k_high, "samples): expression range [",
    round(min(focus_vec[high_samples]), 3), ",",
    round(max(focus_vec[high_samples]), 3), "]\n\n")

# ----- 8. Build ITV gene matrix and reference mean -----
cat("Computing reference mean and ITV scores...\n")

itv_mat <- as.matrix(
  expr_dt[get(gene_col) %in% itv_gene_ids, ..all_sample_cols]
)
rownames(itv_mat) <- expr_dt[get(gene_col) %in% itv_gene_ids, get(gene_col)]

n_itv_genes <- nrow(itv_mat)
cat("  ITV gene pool:", n_itv_genes, "genes\n")

# Gene-wise mean across ALL samples (global population baseline)
# In log scale: subtracting the mean ≡ log(expr / geometric_mean)
mean_expr <- rowMeans(itv_mat, na.rm = TRUE)

# ----- 9. Compute ITV per sample -----
# ITV(s) = summary_fn( |expr[genes, s] - mean[genes]| )
summary_fn <- switch(args$summary_metric,
  "median" = function(x) median(x, na.rm = TRUE),
  "mean"   = function(x) mean(x,   na.rm = TRUE),
  "sum"    = function(x) sum(x,    na.rm = TRUE)
)

compute_itv <- function(sample_id) {
  deviation <- abs(itv_mat[, sample_id] - mean_expr)
  summary_fn(deviation)
}

itv_low  <- sapply(low_samples,  compute_itv)
itv_high <- sapply(high_samples, compute_itv)

cat("  ITV LOW  (n =", length(itv_low),  "):  median =", round(median(itv_low),  4),
    " | range: [", round(min(itv_low),  4), ",", round(max(itv_low),  4), "]\n")
cat("  ITV HIGH (n =", length(itv_high), "):  median =", round(median(itv_high), 4),
    " | range: [", round(min(itv_high), 4), ",", round(max(itv_high), 4), "]\n\n")

# ----- 10. Wilcoxon rank-sum test -----
wt <- wilcox.test(itv_low, itv_high, paired = FALSE)
cat("  Wilcoxon rank-sum test:\n")
cat("  ", format_pval(wt$p.value), "\n\n")

# ----- 11. Assemble long-format data frame -----
plot_df <- data.frame(
  sample    = c(names(itv_low), names(itv_high)),
  itv       = c(itv_low, itv_high),
  group     = factor(
    c(rep("LOW", length(itv_low)), rep("HIGH", length(itv_high))),
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

# Y-axis label reflects the summary metric
metric_y_labels <- c(
  "median" = "ITV: median |log-expr \u2212 population mean|",
  "mean"   = "ITV: mean |log-expr \u2212 population mean|",
  "sum"    = "ITV: sum |log-expr \u2212 population mean|"
)
y_label <- metric_y_labels[[args$summary_metric]]

title_str    <- paste0("Individual Transcriptomic Variability\n",
                       "Focus: ", gene_label, "  |  Condition: ", args$condition_label)
subtitle_str <- format_pval(wt$p.value)

colors_groups <- c("LOW" = "#2E86AB", "HIGH" = "#A23B72")

# Star annotation y-position: just above the violin max
y_max  <- max(plot_df$itv, na.rm = TRUE)
y_star <- y_max * 1.08

p <- ggplot(plot_df, aes(x = group, y = itv, fill = group)) +
  geom_violin(alpha = 0.55, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, alpha = 0.85, outlier.alpha = 0.4, linewidth = 0.4) +
  # Larger, more visible jitter: n_dots = k_low + k_high (~40-80 samples typical)
  geom_jitter(width = 0.10, size = 1.6, alpha = 0.55, color = "black") +
  annotate("text", x = 1.5, y = y_star,
           label = format_pval(wt$p.value), size = 6, vjust = 0) +
  annotate("segment", x = 1, xend = 2, y = y_max * 1.03, yend = y_max * 1.03,
           linewidth = 0.35, color = "grey40") +
  scale_fill_manual(values = colors_groups) +
  labs(
    title    = title_str,
    subtitle = subtitle_str,
    x        = paste0("Sample group (split by ", gene_label, " expression)"),
    y        = y_label,
    caption  = paste0("n_LOW = ", k_low, " samples  |  n_HIGH = ", k_high, " samples",
                      "  |  ", n_itv_genes, " genes  |  metric: ", args$summary_metric)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "grey30"),
    plot.caption     = element_text(hjust = 0.5, size = 8,  color = "grey50"),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# ----- 13. Save outputs -----
out_stem <- file.path(args$output_dir,
  paste0(args$focus_gene, "_sample_variability_", args$summary_metric))

pdf_path <- paste0(out_stem, ".pdf")
ggsave(pdf_path, p, width = 6, height = 7)
cat("  Saved plot: ", pdf_path, "\n")

if (args$save_csv) {
  csv_path <- paste0(out_stem, ".csv")
  write.csv(plot_df, csv_path, row.names = FALSE)
  cat("  Saved data: ", csv_path, "\n")
}

cat("\nDone.\n")
