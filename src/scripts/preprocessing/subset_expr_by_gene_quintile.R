#!/usr/bin/env Rscript
# ==============================================================================
# Subset Expression Matrix by Gene Expression Quintiles
#
# Script by Gabriel Thornes
#
# Last Updated: 21/11/2025
#
# This script:
#   1. Takes a target gene and full expression matrix
#   2. Bins samples into quintiles based on target gene expression
#   3. Extracts top quintile (highest expression) and bottom quintile (lowest)
#   4. Saves subsetted expression matrices for downstream network analysis
#
# Input:
#   - Full expression matrix (genes × samples)
#   - Target gene ID
#
# Output:
#   - Top quintile expression matrix (samples × genes)
#   - Bottom quintile expression matrix (samples × genes)
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages and utilities -----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(argparse)
})

# Source utility functions
source("src/utils/utils_extract_extremes.R")

# ----- 2. Command-line arguments -----
parser <- ArgumentParser(description = 'Subset expression matrix by gene quintiles')
parser$add_argument('--expr-matrix', help = 'Full expression matrix (genes x samples)', 
                   required = TRUE)
parser$add_argument('--target-gene', help = 'Gene ID to bin samples by', 
                   required = TRUE)
parser$add_argument('--output-5th', help = 'Output file for top quintile expression', 
                   required = TRUE)
parser$add_argument('--output-1st', help = 'Output file for bottom quintile expression', 
                   required = TRUE)
parser$add_argument('--n-bins', type = 'integer', help = 'Number of expression bins (default 5 for quintiles)', 
                   default = 5)

args <- parser$parse_args()

cat("=== Gene Expression Quintile Subsetting ===\n")
cat("Expression matrix:", args$expr_matrix, "\n")
cat("Target gene:", args$target_gene, "\n")
cat("Number of bins:", args$n_bins, "\n")
cat("Output (top quintile):", args$output_top, "\n")
cat("Output (bottom quintile):", args$output_bottom, "\n\n")

# ----- 3. Load expression data -----
cat("Loading expression data...\n")

expr_data <- fread(args$expr_matrix, data.table = TRUE)

# Check if first column needs renaming
if (colnames(expr_data)[1] %in% c("V1", "")) {
  setnames(expr_data, 1, "fly_id")
}

if (!"fly_id" %in% colnames(expr_data)) {
  stop("Expression matrix must have 'fly_id' as the first column (gene IDs)")
}

cat("  Loaded:", nrow(expr_data), "genes x", ncol(expr_data) - 1, "samples\n")

# Check if target gene exists
if (!args$target_gene %in% expr_data$fly_id) {
  stop(paste0("Target gene '", args$target_gene, "' not found in expression matrix"))
}

cat("  Target gene found in expression matrix\n\n")

# ----- 4. Bin samples by target gene expression -----
cat("Binning samples into", args$n_bins, "quantiles based on", args$target_gene, "expression...\n")

# Use utility function to bin expression
binned_expr <- bin_expr(expr_data, n = args$n_bins)

cat("  Samples binned successfully\n")

# Check distribution of samples across bins
bin_counts <- binned_expr[fly_id == args$target_gene, .N, by = tile][order(tile)]
cat("\n  Sample distribution across bins for", args$target_gene, ":\n")
print(bin_counts)

# ----- 5. Extract quintile expression matrices -----
cat("\nExtracting expression matrices for extreme quintiles...\n")

# Filter binned_expr to only the target gene (this determines which samples to keep)
target_gene_bins <- binned_expr[fly_id == args$target_gene]

# Bottom quintile (tile = 1, lowest expression)
# Get sample IDs for the bottom quintile based on target gene expression
bottom_samples <- target_gene_bins[tile == 1, as.character(sample)]
cat("  Bottom quintile: selecting", length(bottom_samples), "samples with lowest", 
    args$target_gene, "expression\n")

# Extract those samples from full expression matrix (all genes)
bottom_quintile_mat <- expr_data[, .SD, .SDcols = c("fly_id", bottom_samples)] %>%
  column_to_rownames("fly_id") %>%
  t()

cat("  Bottom quintile matrix:", nrow(bottom_quintile_mat), "samples x", 
    ncol(bottom_quintile_mat), "genes\n")

# Top quintile (tile = n_bins, highest expression)
top_samples <- target_gene_bins[tile == args$n_bins, as.character(sample)]
cat("  Top quintile: selecting", length(top_samples), "samples with highest", 
    args$target_gene, "expression\n")

top_quintile_mat <- expr_data[, .SD, .SDcols = c("fly_id", top_samples)] %>%
  column_to_rownames("fly_id") %>%
  t()

cat("  Top quintile matrix:", nrow(top_quintile_mat), "samples x", 
    ncol(top_quintile_mat), "genes\n")

# ----- 6. Validate and report -----
cat("\nValidating subsets...\n")

# Check target gene expression in each subset
target_gene_bottom <- bottom_quintile_mat[, args$target_gene]
target_gene_top <- top_quintile_mat[, args$target_gene]

cat("  Target gene expression summary:\n")
cat("    Bottom quintile - Mean:", round(mean(target_gene_bottom), 3), 
    "Range: [", round(min(target_gene_bottom), 3), ",", round(max(target_gene_bottom), 3), "]\n")
cat("    Top quintile    - Mean:", round(mean(target_gene_top), 3), 
    "Range: [", round(min(target_gene_top), 3), ",", round(max(target_gene_top), 3), "]\n")

if (mean(target_gene_top) <= mean(target_gene_bottom)) {
  warning("Top quintile mean is not greater than bottom quintile - check data")
}

# ----- 7. Save output files -----
cat("\nSaving output files...\n")

# Create output directories if they don't exist
dir.create(dirname(args$output_bottom), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$output_top), recursive = TRUE, showWarnings = FALSE)

# Convert matrices back to data.table with gene IDs
# Note: matrices are samples × genes, need to transpose for standard format (genes × samples)
bottom_dt <- as.data.table(t(bottom_quintile_mat), keep.rownames = "fly_id")
top_dt <- as.data.table(t(top_quintile_mat), keep.rownames = "fly_id")

# Save as CSV
fwrite(bottom_dt, args$output_bottom)
cat("  Saved bottom quintile:", args$output_bottom, "\n")

fwrite(top_dt, args$output_top)
cat("  Saved top quintile:", args$output_top, "\n")

# ----- 8. Summary report -----
cat("\n=== Subsetting Complete ===\n\n")
cat("Summary:\n")
cat("  Target gene:", args$target_gene, "\n")
cat("  Total samples in full matrix:", ncol(expr_data) - 1, "\n")
cat("  Samples per quintile: ~", round((ncol(expr_data) - 1) / args$n_bins), "\n")
cat("  Bottom quintile samples:", nrow(bottom_quintile_mat), "\n")
cat("  Top quintile samples:", nrow(top_quintile_mat), "\n")
cat("  Genes retained:", ncol(bottom_quintile_mat), "\n\n")

cat("Next steps:\n")
cat("  1. Use these matrices for correlation network construction\n")
cat("  2. Compare networks between high/low expression contexts\n")
cat("  3. Identify differential co-expression patterns\n\n")

cat("Output files:\n")
cat("  ", args$output_bottom, " (genes x samples, bottom quintile)\n")
cat("  ", args$output_top, " (genes x samples, top quintile)\n")
