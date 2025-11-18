#!/usr/bin/env Rscript
# ==============================================================================
# Extract Top and Bottom Expression Quintiles by Diet
#
# Script by Gabriel Thornes
#
# Last Updated: 18/11/2025
#
# This script:
#   1. Loads expression data for Control and High Sugar diets
#   2. Bins expression values into quintiles for each gene
#   3. Extracts samples in the top quintile (highest expression)
#   4. Extracts samples in the bottom quintile (lowest expression)
#   5. Saves expression matrices for each quintile and diet combination
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(argparse)
})

# ----- 2. Source utility functions -----
source("src/utils/utils_extract_extremes.R")

# ----- 3. Command-line arguments -----
parser <- ArgumentParser(description = 'Extract top and bottom expression quintiles for each diet')
parser$add_argument('--control-file', help = 'Path to Control diet expression data file',
                   default = 'dataset/processed/VOOM/voomdataCtrl.txt')
parser$add_argument('--highsugar-file', help = 'Path to High Sugar diet expression data file',
                   default = 'dataset/processed/VOOM/voomdataHS.txt')
parser$add_argument('--output-dir', help = 'Directory to save quintile expression matrices', 
                   default = 'dataset/subset/expression_quintiles')
parser$add_argument('--n-bins', help = 'Number of bins to divide expression into (default: 5 for quintiles)', 
                   default = 5, type = 'integer')
args <- parser$parse_args()

output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Expression Quintile Extraction ===\n")
cat("Control file:", args$control_file, "\n")
cat("High Sugar file:", args$highsugar_file, "\n")
cat("Number of bins:", args$n_bins, "\n")
cat("Output directory:", output_dir, "\n\n")

# ----- 4. Load Control diet expression data -----
cat("Loading Control diet expression data...\n")
ctrl_expr <- fread(args$control_file, header = TRUE)

# Handle row names - first column should be gene IDs
if (colnames(ctrl_expr)[1] %in% c("V1", "")) {
  setnames(ctrl_expr, old = colnames(ctrl_expr)[1], new = "fly_id")
} else if (colnames(ctrl_expr)[1] != "fly_id") {
  setnames(ctrl_expr, old = colnames(ctrl_expr)[1], new = "fly_id")
}

cat("  Control samples:", ncol(ctrl_expr) - 1, "\n")
cat("  Control genes:", nrow(ctrl_expr), "\n\n")

# ----- 5. Load High Sugar diet expression data -----
cat("Loading High Sugar diet expression data...\n")
hs_expr <- fread(args$highsugar_file, header = TRUE)

# Handle row names - first column should be gene IDs
if (colnames(hs_expr)[1] %in% c("V1", "")) {
  setnames(hs_expr, old = colnames(hs_expr)[1], new = "fly_id")
} else if (colnames(hs_expr)[1] != "fly_id") {
  setnames(hs_expr, old = colnames(hs_expr)[1], new = "fly_id")
}

cat("  High Sugar samples:", ncol(hs_expr) - 1, "\n")
cat("  High Sugar genes:", nrow(hs_expr), "\n\n")

# ----- 6. Process Control diet -----
cat("=== Processing Control Diet ===\n")

# Bin expression into quintiles
cat("Binning Control expression into", args$n_bins, "bins...\n")
ctrl_binned <- bin_expr(ctrl_expr, n = args$n_bins)

cat("  Binned data dimensions:", nrow(ctrl_binned), "rows\n")
cat("  Samples per bin (approximate):\n")
bin_counts <- ctrl_binned[, .N, by = tile][order(tile)]
print(bin_counts)

# Extract bottom quintile (tile = 1)
cat("\nExtracting bottom quintile (tile = 1)...\n")
ctrl_bottom <- extract_biase_expr(ctrl_binned, ctrl_expr, nth_tile = 1)
cat("  Bottom quintile dimensions:", nrow(ctrl_bottom), "samples ×", ncol(ctrl_bottom), "genes\n")

# Extract top quintile (tile = n_bins)
cat("Extracting top quintile (tile =", args$n_bins, ")...\n")
ctrl_top <- extract_biase_expr(ctrl_binned, ctrl_expr, nth_tile = args$n_bins)
cat("  Top quintile dimensions:", nrow(ctrl_top), "samples ×", ncol(ctrl_top), "genes\n")

# Save Control quintiles
ctrl_bottom_file <- file.path(output_dir, "control_bottom_quintile.csv")
ctrl_top_file <- file.path(output_dir, "control_top_quintile.csv")

# Convert to data.table with sample names as first column
ctrl_bottom_dt <- data.table(sample = rownames(ctrl_bottom), as.data.frame(ctrl_bottom))
ctrl_top_dt <- data.table(sample = rownames(ctrl_top), as.data.frame(ctrl_top))

fwrite(ctrl_bottom_dt, ctrl_bottom_file)
fwrite(ctrl_top_dt, ctrl_top_file)

cat("  Saved:", ctrl_bottom_file, "\n")
cat("  Saved:", ctrl_top_file, "\n\n")

# ----- 7. Process High Sugar diet -----
cat("=== Processing High Sugar Diet ===\n")

# Bin expression into quintiles
cat("Binning High Sugar expression into", args$n_bins, "bins...\n")
hs_binned <- bin_expr(hs_expr, n = args$n_bins)

cat("  Binned data dimensions:", nrow(hs_binned), "rows\n")
cat("  Samples per bin (approximate):\n")
bin_counts_hs <- hs_binned[, .N, by = tile][order(tile)]
print(bin_counts_hs)

# Extract bottom quintile (tile = 1)
cat("\nExtracting bottom quintile (tile = 1)...\n")
hs_bottom <- extract_biase_expr(hs_binned, hs_expr, nth_tile = 1)
cat("  Bottom quintile dimensions:", nrow(hs_bottom), "samples ×", ncol(hs_bottom), "genes\n")

# Extract top quintile (tile = n_bins)
cat("Extracting top quintile (tile =", args$n_bins, ")...\n")
hs_top <- extract_biase_expr(hs_binned, hs_expr, nth_tile = args$n_bins)
cat("  Top quintile dimensions:", nrow(hs_top), "samples ×", ncol(hs_top), "genes\n")

# Save High Sugar quintiles
hs_bottom_file <- file.path(output_dir, "highsugar_bottom_quintile.csv")
hs_top_file <- file.path(output_dir, "highsugar_top_quintile.csv")

# Convert to data.table with sample names as first column
hs_bottom_dt <- data.table(sample = rownames(hs_bottom), as.data.frame(hs_bottom))
hs_top_dt <- data.table(sample = rownames(hs_top), as.data.frame(hs_top))

fwrite(hs_bottom_dt, hs_bottom_file)
fwrite(hs_top_dt, hs_top_file)

cat("  Saved:", hs_bottom_file, "\n")
cat("  Saved:", hs_top_file, "\n\n")

# ----- 8. Generate summary statistics -----
cat("=== Summary Statistics ===\n\n")

summary_stats <- data.frame(
  Diet = c("Control", "Control", "High Sugar", "High Sugar"),
  Quintile = c("Bottom", "Top", "Bottom", "Top"),
  n_samples = c(nrow(ctrl_bottom), nrow(ctrl_top), nrow(hs_bottom), nrow(hs_top)),
  n_genes = c(ncol(ctrl_bottom), ncol(ctrl_top), ncol(hs_bottom), ncol(hs_top)),
  mean_expression = c(
    mean(as.matrix(ctrl_bottom), na.rm = TRUE),
    mean(as.matrix(ctrl_top), na.rm = TRUE),
    mean(as.matrix(hs_bottom), na.rm = TRUE),
    mean(as.matrix(hs_top), na.rm = TRUE)
  ),
  median_expression = c(
    median(as.matrix(ctrl_bottom), na.rm = TRUE),
    median(as.matrix(ctrl_top), na.rm = TRUE),
    median(as.matrix(hs_bottom), na.rm = TRUE),
    median(as.matrix(hs_top), na.rm = TRUE)
  )
)

print(summary_stats)

# Save summary
summary_file <- file.path(output_dir, "quintile_summary.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat("\nSummary saved to:", summary_file, "\n")

cat("\n=== Extraction Complete ===\n")
cat("Output files:\n")
cat("  - control_bottom_quintile.csv (lowest 20% expression per gene)\n")
cat("  - control_top_quintile.csv (highest 20% expression per gene)\n")
cat("  - highsugar_bottom_quintile.csv (lowest 20% expression per gene)\n")
cat("  - highsugar_top_quintile.csv (highest 20% expression per gene)\n")
cat("  - quintile_summary.csv (summary statistics)\n")
