# !/usr/bin/env Rscript
# ==============================================================================
# Calculate network features with soft thresholded adjacency matrix
#
# Script by Gabriel Thornes
#
# Last Updated: 21/12/2025
#
# This script::
#   1. Takes wide correlation matrix following soft thresholding as input
#   2. Calculates network-level metrics
#   3. Outputs network metrics and gene metrics as summary statistics CSV file
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(WGCNA)
  library(argparse)
})

# ----- 2. Set up utility functions -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")

# ----- 3. Command-line arguments -----
parser <- ArgumentParser(description = 'Calculate network features and identify modules from adjacency matrix')
parser$add_argument('--adjacency_file', help = 'Path to soft-thresholded adjacency matrix CSV file',
                   default = 'results/network_features/soft_threshold/_st_adjacency_matrix.csv')
parser$add_argument('--output_dir', help = 'Directory to save network metrics and module results', 
                   default = 'results/network_features/features_calc/_adjacency')
parser$add_argument('--prefix', help = 'Prefix name of results')
parser$add_argument('--connection_threshold', help = 'Connection threshold for degree calculation', 
                   default = 0.05, type = 'double')
parser$add_argument('--min_module_size', help = 'Minimum module size for WGCNA', 
                   default = 20, type = 'integer')
parser$add_argument('--merge_threshold', help = 'Module merge threshold (height to merge similar modules)', 
                   default = 0.15, type = 'double')
parser$add_argument('--deep_split', help = 'Deep split parameter for module detection (0-4)', 
                   default = 3, type = 'integer')
parser$add_argument('--top_n_hubs', help = 'Number of top hub genes to identify per module', 
                   default = 5, type = 'integer')
args <- parser$parse_args()

prefix <- args$prefix
output_dir <- args$output_dir

parent_dir <- dirname(output_dir)
if (str_detect(parent_dir, "modules$")) {
  tom_basedir <- basename(output_dir)
  tom_dir <- file.path(parent_dir, tom_basedir)
  dir_tom_file <- output_dir
  module_output_dir <- file.path(tom_dir, prefix)

}else{
  stop("Output directory must be within a 'modules' folder.")
}

cat("=== Network Metrics Analysis ===\n")
cat("Adjacency file:", args$adjacency_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Connection threshold:", args$connection_threshold, "\n")
cat("Min module size:", args$min_module_size, "\n")
cat("Merge threshold:", args$merge_threshold, "\n")
cat("Deep split:", args$deep_split, "\n")
cat("Top hubs per module:", args$top_n_hubs, "\n\n")

# Read adjacency matrix
adjacency <- read.csv(args$adjacency_file)
cat("Adjacency matrix loaded:", nrow(adjacency), "x", ncol(adjacency), "\n\n")

# Create output directory if it doesn't exist
walk(c(output_dir, module_output_dir), create_directories)

# ----- 4. Calculate network-level summary statistics -----

# Convert adjacency to proper matrix format
adjacency_matrix <- as.matrix(adjacency[,-1])
rownames(adjacency_matrix) <- adjacency[[1]]
colnames(adjacency_matrix) <- adjacency[[1]]

cat("Matrix conversion completed:\n")
cat("- Adjacency matrix dimensions:", dim(adjacency_matrix), "\n\n")

cat("Generating threshold analysis plot...\n")
plot_threshold_analysis(adjacency_matrix, 
                        output_file = file.path(output_dir,  paste0(prefix,"_adjacency_threshold_analysis.pdf")), 
                        threshold_range = seq(0.02, 0.4, by = 0.02),
                        matrix_name = "Adjacency Matrix")

# Calculate network metrics
cat("\nCalculating network metrics...\n")
adjacency_metrics <- calculate_network_metrics(adjacency_matrix, 
                                               "Soft-Thresholded Adjacency", 
                                               connection_threshold = args$connection_threshold)

# Create summary data frame for export
summary_stats <- data.frame(
    Matrix_Type = "Adjacency",
    Mean_Abs_Correlation = adjacency_metrics$mean_abs_corr,
    Median_Abs_Correlation = adjacency_metrics$median_abs_corr,
    Max_Correlation = adjacency_metrics$max_corr,
    Min_Correlation = adjacency_metrics$min_corr,
    Mean_Connectivity = mean(adjacency_metrics$weighted_connectivity, na.rm = TRUE),
    Mean_Degree = mean(adjacency_metrics$degree, na.rm = TRUE)
)

# Save summary statistics
summary_file <- file.path(output_dir,  paste0(prefix,"_network_metrics_summary.csv"))
write.csv(summary_stats, file = summary_file, row.names = FALSE)

cat("\n=== Network Summary Statistics ===\n")
print(summary_stats)
cat("\nSummary statistics saved to:", summary_file, "\n")

# ----- 5. Calculate modules and hub genes via WGCNA -----

cat("\n=== WGCNA Module Detection ===\n")
adjacency_modules <- create_modules(
    adjacency_matrix = adjacency_matrix,
    min_module_size = args$min_module_size,
    merge_threshold = args$merge_threshold,
    deep_split = args$deep_split,
    output_dir = module_output_dir,
    save_plots = TRUE
)

# Identify hub genes
cat("\n=== Identifying Hub Genes ===\n")
adjacency_hubs <- identify_hubs(
    module_results = adjacency_modules,
    top_n_hubs = args$top_n_hubs,
    output_dir = module_output_dir,
    save_files = TRUE
)

# Save TOM matrices (from the module results)
tom_file <- file.path(tom_dir,  paste0(prefix,"_tom.csv"))
tomd_file <- file.path(tom_dir,  paste0(prefix,"_tomd.csv"))

write.csv(adjacency_modules$tom_matrix, file = tom_file, row.names = TRUE)
write.csv(adjacency_modules$tomd_matrix, file = tomd_file, row.names = TRUE)

cat("\nTOM matrices saved:\n")
cat("- TOM matrix:", tom_file, "\n")
cat("- TOMD matrix:", tomd_file, "\n")

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", output_dir, "\n")
