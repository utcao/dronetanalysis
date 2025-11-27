#!/usr/bin/env Rscript
# ==============================================================================
# Calculate Differential TOM Matrix and Identify Modules
#
# Script by Gabriel Thornes
#
# Last Updated: 24/11/2025
#
# This script:
#   1. Loads two TOM matrices (e.g., from top and bottom quintile networks)
#   2. Calculates the absolute difference between them
#   3. Saves the differential TOM matrix
#   4. Identifies modules in the differential network
#   5. Calculates hub genes for differential modules
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(argparse)
})

# ----- 2. Source utility functions -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")

# ----- 3. Command-line arguments -----
parser <- ArgumentParser(description = 'Calculate differential TOM and identify modules')
parser$add_argument('--tom_1st', help = 'Path to first TOM matrix (bottom quintile)',
                   required = TRUE)
parser$add_argument('--tom_5th', help = 'Path to second TOM matrix (top quintile)',
                   required = TRUE)
parser$add_argument('--output_dir', help = 'Directory to save differential TOM results',
                   required = TRUE)
parser$add_argument('--connection_threshold', help = 'Connection threshold for network metrics',
                   default = 0.01, type = 'double')
parser$add_argument('--min_module_size', help = 'Minimum module size for WGCNA',
                   default = 20, type = 'integer')
parser$add_argument('--merge_threshold', help = 'Module merge threshold',
                   default = 0.15, type = 'double')
parser$add_argument('--deep_split', help = 'Deep split parameter for module detection (0-4)',
                   default = 3, type = 'integer')
parser$add_argument('--top_n_hubs', help = 'Number of top hub genes per module',
                   default = 5, type = 'integer')
args <- parser$parse_args()

output_dir <- args$output_dir
module_output_dir <- file.path(output_dir, "modules")

# Extract gene name from output directory (assumes dir ends with gene name)
gene_name <- basename(output_dir)

cat("=== Differential TOM Analysis ===\n")
cat("Target gene:", gene_name, "\n")
cat("TOM matrix 1 (bottom quintile):", args$tom_1st, "\n")
cat("TOM matrix 2 (top quintile):", args$tom_5th, "\n")
cat("Output directory:", output_dir, "\n")
cat("Connection threshold:", args$connection_threshold, "\n")
cat("Min module size:", args$min_module_size, "\n")
cat("Merge threshold:", args$merge_threshold, "\n")
cat("Deep split:", args$deep_split, "\n")
cat("Top hubs per module:", args$top_n_hubs, "\n\n")

# Create output directories
create_directories(output_dir)
create_directories(module_output_dir)

# ----- 4. Load TOM matrices -----
cat("Loading TOM matrices...\n")

# Load first TOM matrix (bottom quintile)
tom_1st <- fread(args$tom_1st, data.table = FALSE)
# First column contains gene names
gene_names_1st <- as.character(tom_1st[[1]])
tom_1st_mat <- as.matrix(tom_1st[, -1])
rownames(tom_1st_mat) <- gene_names_1st
colnames(tom_1st_mat) <- gene_names_1st
storage.mode(tom_1st_mat) <- "double"

cat("  TOM 1st quintile:", nrow(tom_1st_mat), "x", ncol(tom_1st_mat), "genes\n")

# Load second TOM matrix (top quintile)
tom_5th <- fread(args$tom_5th, data.table = FALSE)
gene_names_5th <- as.character(tom_5th[[1]])
tom_5th_mat <- as.matrix(tom_5th[, -1])
rownames(tom_5th_mat) <- gene_names_5th
colnames(tom_5th_mat) <- gene_names_5th
storage.mode(tom_5th_mat) <- "double"

cat("  TOM 5th quintile:", nrow(tom_5th_mat), "x", ncol(tom_5th_mat), "genes\n")

# ----- 5. Validate matrices -----
cat("\nValidating matrices...\n")

# Check dimensions match
if (nrow(tom_1st_mat) != nrow(tom_5th_mat) || ncol(tom_1st_mat) != ncol(tom_5th_mat)) {
  stop("ERROR: TOM matrix dimensions do not match!\n",
       "  TOM 1st: ", nrow(tom_1st_mat), "x", ncol(tom_1st_mat), "\n",
       "  TOM 5th: ", nrow(tom_5th_mat), "x", ncol(tom_5th_mat))
}

# Check gene names match
if (!all(gene_names_1st == gene_names_5th)) {
  cat("  Warning: Gene names do not match exactly. Using intersection...\n")
  common_genes <- intersect(gene_names_1st, gene_names_5th)
  cat("  Common genes:", length(common_genes), "\n")
  
  tom_1st_mat <- tom_1st_mat[common_genes, common_genes]
  tom_5th_mat <- tom_5th_mat[common_genes, common_genes]
} else {
  cat("  Gene names match: ", length(gene_names_1st), "genes\n")
}

# Check for NA/Inf values
if (any(is.na(tom_1st_mat)) || any(is.na(tom_5th_mat))) {
  cat("  Warning: NA values detected. Replacing with 0...\n")
  tom_1st_mat[is.na(tom_1st_mat)] <- 0
  tom_5th_mat[is.na(tom_5th_mat)] <- 0
}

if (any(is.infinite(tom_1st_mat)) || any(is.infinite(tom_5th_mat))) {
  cat("  Warning: Infinite values detected. Replacing with 0...\n")
  tom_1st_mat[is.infinite(tom_1st_mat)] <- 0
  tom_5th_mat[is.infinite(tom_5th_mat)] <- 0
}

cat("  Validation complete\n\n")

# ----- 6. Calculate differential TOM -----
cat("Calculating differential TOM (absolute difference)...\n")

# Calculate absolute difference
diff_tom <- abs(tom_5th_mat - tom_1st_mat)

# Ensure symmetry (should already be symmetric, but ensure it)
diff_tom <- (diff_tom + t(diff_tom)) / 2

# Ensure diagonal is zero (genes don't differ from themselves)
diag(diff_tom) <- 0

cat("  Differential TOM calculated\n")
cat("  Matrix dimensions:", dim(diff_tom), "\n")
cat("  Value range:", round(range(diff_tom), 4), "\n")
cat("  Mean difference:", round(mean(diff_tom[upper.tri(diff_tom)]), 4), "\n")
cat("  Median difference:", round(median(diff_tom[upper.tri(diff_tom)]), 4), "\n\n")

# ----- 7. Save differential TOM -----
cat("Saving differential TOM matrix...\n")

diff_tom_file <- file.path(output_dir, paste0(gene_name, "_diff_tom_matrix.csv"))
diff_tom_df <- as.data.frame(diff_tom)
diff_tom_df <- cbind(gene = rownames(diff_tom), diff_tom_df)
write.csv(diff_tom_df, file = diff_tom_file, row.names = FALSE)

cat("  Saved:", diff_tom_file, "\n\n")

# ----- 9. Calculate network metrics -----
cat("\nCalculating network metrics for differential TOM...\n")
diff_tom_metrics <- calculate_network_metrics(diff_tom,
                                              "Differential TOM",
                                              connection_threshold = args$connection_threshold)

# Create summary data frame
summary_stats <- data.frame(
  Matrix_Type = "Differential_TOM",
  Mean_Abs_Difference = diff_tom_metrics$mean_abs_corr,
  Median_Abs_Difference = diff_tom_metrics$median_abs_corr,
  Max_Difference = diff_tom_metrics$max_corr,
  Min_Difference = diff_tom_metrics$min_corr,
  Mean_Connectivity = mean(diff_tom_metrics$weighted_connectivity, na.rm = TRUE),
  Mean_Degree = mean(diff_tom_metrics$degree, na.rm = TRUE)
)

# Save summary statistics
summary_file <- file.path(output_dir, paste0(gene_name, "_diff_tom_network_metrics.csv"))
write.csv(summary_stats, file = summary_file, row.names = FALSE)

cat("\n=== Differential TOM Network Statistics ===\n")
print(summary_stats)
cat("\nSummary saved to:", summary_file, "\n")

# ----- 10. Identify modules in differential TOM -----
cat("\n=== Module Detection in Differential TOM ===\n")
cat("Genes with large TOM differences will cluster into differential modules\n\n")

diff_modules <- create_modules(
  adjacency_matrix = diff_tom,
  min_module_size = args$min_module_size,
  merge_threshold = args$merge_threshold,
  deep_split = args$deep_split,
  output_dir = module_output_dir,
  save_plots = TRUE
)

# ----- 11. Identify hub genes in differential modules -----
cat("\n=== Identifying Hub Genes in Differential Modules ===\n")
cat("Hub genes = genes with largest TOM changes across conditions\n\n")

diff_hubs <- identify_hubs(
  module_results = diff_modules,
  top_n_hubs = args$top_n_hubs,
  output_dir = module_output_dir,
  save_files = TRUE
)

# ----- 12. Save additional differential information -----
cat("\n=== Saving Additional Results ===\n")

# Save gene-level differential metrics
gene_diff_metrics <- data.frame(
  gene = rownames(diff_tom),
  mean_tom_1st = rowMeans(tom_1st_mat),
  mean_tom_5th = rowMeans(tom_5th_mat),
  mean_diff_tom = rowMeans(diff_tom),
  max_diff_tom = apply(diff_tom, 1, max),
  module = diff_modules$final_colours,
  stringsAsFactors = FALSE
)

# Add connectivity from differential network
gene_diff_metrics$diff_connectivity <- diff_tom_metrics$weighted_connectivity
gene_diff_metrics$diff_degree <- diff_tom_metrics$degree

# Sort by mean differential TOM (genes with largest average changes)
gene_diff_metrics <- gene_diff_metrics[order(-gene_diff_metrics$mean_diff_tom), ]

gene_metrics_file <- file.path(output_dir, paste0(gene_name, "_gene_differential_metrics.csv"))
write.csv(gene_diff_metrics, file = gene_metrics_file, row.names = FALSE)

cat("  Gene-level metrics saved:", gene_metrics_file, "\n")

# ----- 13. Summary report -----
cat("\n=== Differential TOM Analysis Complete ===\n\n")

cat("Summary:\n")
cat("  Total genes analyzed:", nrow(diff_tom), "\n")
cat("  Differential modules identified:", length(unique(diff_modules$final_colours[diff_modules$final_colours != "grey"])), "\n")
cat("  Genes in grey module (no differential pattern):", sum(diff_modules$final_colours == "grey"), "\n")
cat("  Mean TOM difference:", round(mean(diff_tom[upper.tri(diff_tom)]), 4), "\n\n")

cat("Files generated:\n")
cat("  1.", paste0(gene_name, "_differential_tom_matrix.csv"), "- Full differential TOM matrix\n")
cat("  2.", paste0(gene_name, "_diff_tom_network_metrics.csv"), "- Network-level statistics\n")
cat("  3.", paste0(gene_name, "_gene_differential_metrics.csv"), "- Gene-level differential metrics\n")
cat("  4. modules/ - Module detection results and hub genes\n\n")

cat("Interpretation:\n")
cat("  - Differential TOM shows which genes change their network position\n")
cat("  - High values = gene's network neighborhood changed substantially\n")
cat("  - Modules = groups of genes with coordinated TOM changes\n")
cat("  - Hub genes = genes with largest differential connectivity (key rewiring)\n\n")

cat("Results saved to:", output_dir, "\n")