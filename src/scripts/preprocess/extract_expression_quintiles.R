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
ctrl_expr <- fread(args$control_file)

# fread detects missing column name and assigns "V1" to first column (gene IDs)
setnames(ctrl_expr, "V1", "fly_id")

cat("  Control samples:", ncol(ctrl_expr) - 1, "\n")
cat("  Control genes:", nrow(ctrl_expr), "\n")
cat("  Sample ID examples:", paste(colnames(ctrl_expr)[2:min(6, ncol(ctrl_expr))], collapse = ", "), "\n\n")

# ----- 5. Load High Sugar diet expression data -----
cat("Loading High Sugar diet expression data...\n")
hs_expr <- fread(args$highsugar_file)

# fread detects missing column name and assigns "V1" to first column (gene IDs)
setnames(hs_expr, "V1", "fly_id")

cat("  High Sugar samples:", ncol(hs_expr) - 1, "\n")
cat("  High Sugar genes:", nrow(hs_expr), "\n")
cat("  Sample ID examples:", paste(colnames(hs_expr)[2:min(6, ncol(hs_expr))], collapse = ", "), "\n\n")

# ----- 6. Process Control diet -----
cat("=== Processing Control Diet ===\n")

# Bin expression into quintiles (per gene)
cat("Binning Control expression into", args$n_bins, "bins per gene...\n")
ctrl_binned <- bin_expr(ctrl_expr, n = args$n_bins)

# For each gene, get the samples in bottom and top quintiles
cat("Extracting bottom quintile samples for each gene...\n")
ctrl_bottom_long <- ctrl_binned[tile == 1, .(gene = fly_id, sample = as.character(sample), value = logCPM)]

cat("Extracting top quintile samples for each gene...\n")
ctrl_top_long <- ctrl_binned[tile == args$n_bins, .(gene = fly_id, sample = as.character(sample), value = logCPM)]

# Convert to wide format: rows = genes, columns = samples (with expression values)
cat("Converting to wide format...\n")
ctrl_bottom_wide <- dcast(ctrl_bottom_long, gene ~ sample, value.var = "value")
ctrl_top_wide <- dcast(ctrl_top_long, gene ~ sample, value.var = "value")

cat("  Bottom quintile:", nrow(ctrl_bottom_wide), "genes ×", ncol(ctrl_bottom_wide)-1, "unique samples\n")
cat("  Top quintile:", nrow(ctrl_top_wide), "genes ×", ncol(ctrl_top_wide)-1, "unique samples\n")

# Save Control quintiles
ctrl_bottom_file <- file.path(output_dir, "control_bottom_quintile.csv")
ctrl_top_file <- file.path(output_dir, "control_top_quintile.csv")

fwrite(ctrl_bottom_wide, ctrl_bottom_file)
fwrite(ctrl_top_wide, ctrl_top_file)

cat("  Saved:", ctrl_bottom_file, "\n")
cat("  Saved:", ctrl_top_file, "\n\n")

# ----- 7. Process High Sugar diet -----
cat("=== Processing High Sugar Diet ===\n")

# Bin expression into quintiles (per gene)
cat("Binning High Sugar expression into", args$n_bins, "bins per gene...\n")
hs_binned <- bin_expr(hs_expr, n = args$n_bins)

# For each gene, get the samples in bottom and top quintiles
cat("Extracting bottom quintile samples for each gene...\n")
hs_bottom_long <- hs_binned[tile == 1, .(gene = fly_id, sample = as.character(sample), value = logCPM)]

cat("Extracting top quintile samples for each gene...\n")
hs_top_long <- hs_binned[tile == args$n_bins, .(gene = fly_id, sample = as.character(sample), value = logCPM)]

# Convert to wide format: rows = genes, columns = samples (with expression values)
cat("Converting to wide format...\n")
hs_bottom_wide <- dcast(hs_bottom_long, gene ~ sample, value.var = "value")
hs_top_wide <- dcast(hs_top_long, gene ~ sample, value.var = "value")

cat("  Bottom quintile:", nrow(hs_bottom_wide), "genes ×", ncol(hs_bottom_wide)-1, "unique samples\n")
cat("  Top quintile:", nrow(hs_top_wide), "genes ×", ncol(hs_top_wide)-1, "unique samples\n")

# Save High Sugar quintiles
hs_bottom_file <- file.path(output_dir, "highsugar_bottom_quintile.csv")
hs_top_file <- file.path(output_dir, "highsugar_top_quintile.csv")

fwrite(hs_bottom_wide, hs_bottom_file)
fwrite(hs_top_wide, hs_top_file)

cat("  Saved:", hs_bottom_file, "\n")
cat("  Saved:", hs_top_file, "\n\n")

# ----- 8. Generate summary statistics -----
cat("=== Summary Statistics ===\n\n")

summary_stats <- data.frame(
  Diet = c("Control", "Control", "High Sugar", "High Sugar"),
  Quintile = c("Bottom", "Top", "Bottom", "Top"),
  n_genes = c(nrow(ctrl_bottom_wide), nrow(ctrl_top_wide), 
              nrow(hs_bottom_wide), nrow(hs_top_wide)),
  n_samples_per_gene = c(
    ncol(ctrl_bottom_wide) - 1,
    ncol(ctrl_top_wide) - 1,
    ncol(hs_bottom_wide) - 1,
    ncol(hs_top_wide) - 1
  ),
  mean_expression = c(
    mean(as.matrix(ctrl_bottom_wide[, -1]), na.rm = TRUE),
    mean(as.matrix(ctrl_top_wide[, -1]), na.rm = TRUE),
    mean(as.matrix(hs_bottom_wide[, -1]), na.rm = TRUE),
    mean(as.matrix(hs_top_wide[, -1]), na.rm = TRUE)
  )
)

print(summary_stats)

# Save summary
summary_file <- file.path(output_dir, "quintile_summary.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat("\nSummary saved to:", summary_file, "\n")

cat("\n=== Extraction Complete ===\n")
cat("Output files:\n")
cat("  - control_bottom_quintile.csv (lowest 20% samples per gene)\n")
cat("  - control_top_quintile.csv (highest 20% samples per gene)\n")
cat("  - highsugar_bottom_quintile.csv (lowest 20% samples per gene)\n")
cat("  - highsugar_top_quintile.csv (highest 20% samples per gene)\n")
cat("  - quintile_summary.csv (summary statistics)\n")
cat("\nNote: Each gene has its own set of top/bottom samples.\n")
cat("      Matrix format: rows = genes, columns = samples\n")
cat("      NA values indicate samples not in that gene's quintile.\n")
