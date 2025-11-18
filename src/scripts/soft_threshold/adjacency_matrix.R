# !/usr/bin/env Rscript
# ==============================================================================
# Create adjacency matrices and apply soft thresholding
#
# Script by Gabriel Thornes
#
# Last Updated: 18/11/2025
#
# This script:
#   1. Takes correlation matrix as input
#   2. Creates binary matrix to preserve negative correlations
#   3. Calculates adjacency matrices
#   4. Automatically selects soft-thresholding power based on criteria
#   5. Applies soft threshold and outputs final adjacency matrix
#   6. Outputs adjacency matrices, plots and fit indices as CSV files
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
library(argparse)
library(data.table)
library(dplyr)
library(tidyr)
library(WGCNA)
library(yaml)
})

# ----- 2. Parse command line arguments -----
parser <- ArgumentParser(description = 'Create adjacency matrices and apply soft threshold')
parser$add_argument('--input', type="character", default = "results/spearman_correlation/spearman_matrices/HS_spearman_correlation_matrix.csv",
                   help='Input correlation matrix file path')
parser$add_argument('--output-dir', type="character", default = "results/network_features",
                   help='Output directory for adjacency matrices and plots')
parser$add_argument('--power-range', type="character", default="1:20",
                   help='Power range for soft thresholding (e.g. "1:20")')
parser$add_argument('--network-type', type="character", default="unsigned",
                   help='Network type (unsigned or signed)')
parser$add_argument('--r2-threshold', type="double", default=0.9,
                   help='R^2 threshold for scale-free topology fit (default: 0.9)')
parser$add_argument('--mean-k-threshold', type="double", default=100,
                   help='Maximum mean connectivity threshold (default: 100)')

args <- parser$parse_args()

# Parse power range
power_range <- eval(parse(text = args$power_range))

# ----- 3. Set Input and Output Paths -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")

matrix_file <- args$input
output_dir <- ifelse(is.null(args$output_dir), dirname(dirname(args$input)), args$output_dir)

cat("Input file:", matrix_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Power range:", paste(range(power_range), collapse = "-"), "\n")
cat("R^2 threshold:", args$r2_threshold, "\n")
cat("Mean connectivity threshold:", args$mean_k_threshold, "\n\n")

# ----- 4. Load correlation matrix -----
corr_matrix <- fread(matrix_file)
cat("File read:", matrix_file,"\n")
df.corr.m <- as.matrix(corr_matrix[,-1, with=FALSE])
rownames(df.corr.m) <- corr_matrix[[1]]

# ----- 5. Create binary matrix to preserve negative correlations -----
binary_output_dir <- file.path(output_dir, "binary_matrix")
create_directories(binary_output_dir)
binary_output_file <- file.path(binary_output_dir, "HS_binary_signed_matrix.csv")
binary_matrix <- create_correlation_sign_matrix(df.corr.m, output_file = binary_output_file)
cat("Binary sign matrix created and saved to:", binary_output_file, "\n\n")

# ----- 6. Calculate adjacency matrix -----
adjacency <- abs(df.corr.m)

# Save unsigned adjacency matrix
unsigned_dir <- file.path(output_dir, "soft_threshold")
create_directories(unsigned_dir)
unsigned_output_file <- file.path(unsigned_dir, "HS_unsigned_adjacency_matrix.csv")
write.csv(adjacency, unsigned_output_file, row.names = TRUE)
cat("Unsigned adjacency matrix saved to:", unsigned_output_file, "\n\n")

# ----- 7. Analyze network topology and select soft threshold power -----
cat("=== Analyzing Network Topology ===\n")
sft <- pickSoftThreshold(adjacency, powerVector = power_range,
                         networkType = args$network_type, verbose = 5)

# Save fit indices
fit_indices_file <- file.path(unsigned_dir, "HS_soft_threshold_fit_indices.csv")
write.csv(sft$fitIndices, file = fit_indices_file, row.names = FALSE)
cat("\nFit indices saved to:", fit_indices_file, "\n")

# Generate diagnostic plots
plot_file <- file.path(unsigned_dir, "HS_soft_thresholding_diagnostics.pdf")
sft_plot(sft, plot_file, power_range)
cat("Diagnostic plots saved to:", plot_file, "\n\n")

# ----- 8. Automatically select power based on criteria -----
cat("=== Selecting Optimal Soft Threshold Power ===\n")

# Find powers that meet both criteria
fit_df <- sft$fitIndices
valid_powers <- fit_df %>%
  filter(SFT.R.sq >= args$r2_threshold, 
         mean.k. <= args$mean_k_threshold) %>%
  arrange(Power)

if (nrow(valid_powers) > 0) {
  # Choose the lowest power that meets criteria (most conservative)
  selected_power <- valid_powers$Power[1]
  cat("Selected power:", selected_power, "\n")
  cat("  - R^2 =", valid_powers$SFT.R.sq[1], "(threshold:", args$r2_threshold, ")\n")
  cat("  - Mean connectivity =", valid_powers$mean.k.[1], "(threshold:", args$mean_k_threshold, ")\n\n")
} else {
  # No power meets both criteria - select based on relaxed criteria
  cat("WARNING: No power meets both criteria. Selecting best available option...\n")
  
  # Prioritize R^2 threshold
  r2_met <- fit_df %>% filter(SFT.R.sq >= args$r2_threshold)
  
  if (nrow(r2_met) > 0) {
    # Choose power with lowest mean connectivity among those meeting R^2
    selected_power <- r2_met %>% arrange(mean.k.) %>% pull(Power) %>% .[1]
    cat("Selected power:", selected_power, "(based on R^2 threshold)\n")
  } else {
    # Choose power with highest R^2
    selected_power <- fit_df %>% arrange(desc(SFT.R.sq)) %>% pull(Power) %>% .[1]
    cat("Selected power:", selected_power, "(highest R^2 available)\n")
  }
  
  selected_stats <- fit_df %>% filter(Power == selected_power)
  cat("  - R^2 =", selected_stats$SFT.R.sq, "\n")
  cat("  - Mean connectivity =", selected_stats$mean.k., "\n\n")
}

# Save selected power
power_file <- file.path(unsigned_dir, "HS_selected_soft_power.txt")
writeLines(as.character(selected_power), power_file)
cat("Selected power saved to:", power_file, "\n\n")

# ----- 9. Apply soft threshold -----
cat("=== Applying Soft Threshold ===\n")
soft_adj <- adjacency^selected_power

soft_output_file <- file.path(unsigned_dir, "HS_soft_thresholded_adjacency_matrix.csv")
write.csv(soft_adj, soft_output_file, row.names = TRUE)
cat("Soft-thresholded adjacency matrix saved to:", soft_output_file, "\n\n")

# ----- 10. Summary -----
cat("=== Summary ===\n")
cat("Input correlation matrix:", matrix_file, "\n")
cat("Network type:", args$network_type, "\n")
cat("Selected soft threshold power:", selected_power, "\n")
cat("Output files:\n")
cat("  - Binary sign matrix:", binary_output_file, "\n")
cat("  - Unsigned adjacency:", unsigned_output_file, "\n")
cat("  - Fit indices:", fit_indices_file, "\n")
cat("  - Diagnostic plots:", plot_file, "\n")
cat("  - Selected power:", power_file, "\n")
cat("  - Soft-thresholded adjacency:", soft_output_file, "\n")
cat("\nPipeline complete!\n")