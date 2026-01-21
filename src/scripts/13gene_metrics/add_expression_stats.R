#!/usr/bin/env Rscript
# ==============================================================================
# Add expression statistics (mean and variance) to gene metrics
#
# Script by Gabriel Thornes
#
# Last Updated: 11/11/2025
#
# This script:
#   1. Loads processed expression data (VOOM or VST)
#   2. Calculates mean and variance of expression for each gene
#   3. Merges these statistics with existing gene metrics CSV
#   4. Outputs enriched gene metrics file
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(yaml)
  library(argparse)
})

# ----- 2. Command-line arguments -----
parser <- ArgumentParser(description = 'Add expression statistics to gene metrics')
parser$add_argument('--gene-metrics-file', help = 'Path to gene metrics CSV file', 
                   default = 'results/network_features/gene_metrics/HS_adjacency_gene_metrics.csv')
parser$add_argument('--expression-files', nargs = '+', help = 'Path(s) to expression data files (space-separated)', 
                   default = 'dataset/processed/VOOM/voomdataHS.txt')
parser$add_argument('--output-file', help = 'Path to write enriched gene metrics CSV', 
                   default = 'results/network_features/gene_metrics/HS_adjacency_gene_metrics_with_expression.csv')
parser$add_argument('--expression-type', help = "Type of expression data ('VOOM' or 'VST')", 
                   default = 'VOOM')
args <- parser$parse_args()

cat("=== Adding Expression Statistics to Gene Metrics ===\n")
cat("Gene metrics file:", args$gene_metrics_file, "\n")
cat("Expression files:", paste(args$expression_files, collapse = ", "), "\n")
cat("Output file:", args$output_file, "\n\n")

# ----- 3. Load gene metrics -----
cat("Loading gene metrics...\n")
gene_metrics <- read.csv(args$gene_metrics_file, stringsAsFactors = FALSE)
cat("  Loaded metrics for", nrow(gene_metrics), "genes\n")

# ----- 4. Load and combine expression data -----
cat("\nLoading expression data...\n")

expr_list <- lapply(args$expression_files, function(f) {
  if (!file.exists(f)) {
    stop("Expression file not found: ", f)
  }
  cat("  Reading:", f, "\n")
  
  # Read expression file (assumes gene names in first column/row)
  expr_data <- read.table(f, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  cat("    Dimensions:", nrow(expr_data), "x", ncol(expr_data), "\n")
  return(expr_data)
})

# ----- 5. Combine expression data (if multiple files) -----
if (length(expr_list) == 1) {
  combined_expr <- expr_list[[1]]
} else {
  cat("\nCombining multiple expression files...\n")
  # Combine by taking all samples from all files
  combined_expr <- do.call(cbind, expr_list)
  cat("  Combined dimensions:", nrow(combined_expr), "x", ncol(combined_expr), "\n")
}

# ----- 6. Calculate expression statistics -----
cat("\nCalculating expression statistics...\n")

expr_stats <- data.frame(
  gene = rownames(combined_expr),
  mean_expression = rowMeans(combined_expr, na.rm = TRUE),
  variance_expression = apply(combined_expr, 1, var, na.rm = TRUE),
  sd_expression = apply(combined_expr, 1, sd, na.rm = TRUE),
  mad_expression = apply(combined_expr, 1, mad, na.rm = TRUE),
  median_expression = apply(combined_expr, 1, median, na.rm = TRUE),
  min_expression = apply(combined_expr, 1, min, na.rm = TRUE),
  max_expression = apply(combined_expr, 1, max, na.rm = TRUE),
  stringsAsFactors = FALSE
)

cat("  Expression statistics computed for", nrow(expr_stats), "genes\n")
cat("  Mean expression range:", round(range(expr_stats$mean_expression, na.rm = TRUE), 4), "\n")
cat("  Variance expression range:", round(range(expr_stats$variance_expression, na.rm = TRUE), 4), "\n")
cat("  MAD expression range:", round(range(expr_stats$mad_expression, na.rm = TRUE), 4), "\n")

# ----- 7. Merge with gene metrics -----
cat("\nMerging expression statistics with gene metrics...\n")

# Check if gene column names match
cat("  Gene metrics file has", nrow(gene_metrics), "genes\n")
cat("  Expression file has", nrow(expr_stats), "genes\n")

# Merge by gene name
enriched_metrics <- merge(gene_metrics, expr_stats, by = "gene", all.x = TRUE)

# Check for unmatched genes
n_unmatched <- sum(is.na(enriched_metrics$mean_expression))
if (n_unmatched > 0) {
  cat("  WARNING:", n_unmatched, "genes in metrics file not found in expression data (NA values)\n")
}

cat("  Final merged data:", nrow(enriched_metrics), "genes x", ncol(enriched_metrics), "columns\n")

# ----- 8. Save enriched metrics -----
cat("\nSaving enriched gene metrics...\n")
write.csv(enriched_metrics, args$output_file, row.names = FALSE)
cat("  Saved to:", args$output_file, "\n")

# ----- 9. Print summary -----
cat("\n=== Summary ===\n")
cat("Columns in output file:\n")
cat(paste(colnames(enriched_metrics), collapse = ", "), "\n\n")

cat("Expression statistics summary:\n")
cat("  Mean expression: mean =", round(mean(enriched_metrics$mean_expression, na.rm = TRUE), 4), 
    ", range = [", round(min(enriched_metrics$mean_expression, na.rm = TRUE), 4), ", ",
    round(max(enriched_metrics$mean_expression, na.rm = TRUE), 4), "]\n")
cat("  Variance expression: mean =", round(mean(enriched_metrics$variance_expression, na.rm = TRUE), 4), 
    ", range = [", round(min(enriched_metrics$variance_expression, na.rm = TRUE), 4), ", ",
    round(max(enriched_metrics$variance_expression, na.rm = TRUE), 4), "]\n")
cat("  MAD expression: mean =", round(mean(enriched_metrics$mad_expression, na.rm = TRUE), 4), 
    ", range = [", round(min(enriched_metrics$mad_expression, na.rm = TRUE), 4), ", ",
    round(max(enriched_metrics$mad_expression, na.rm = TRUE), 4), "]\n")

cat("\nMerge complete.\n")
