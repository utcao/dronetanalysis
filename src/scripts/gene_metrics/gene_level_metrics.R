#!/usr/bin/env Rscript
# ==============================================================================
# Calculate comprehensive gene-level network metrics
#
# Script by Gabriel Thornes
#
# Last Updated: 11/11/2025
#
# This script:
#   1. Loads correlation and adjacency matrices
#   2. Calculates basic network metrics for each gene (degree, connectivity, centrality)
#   3. Integrates module-based metrics from WGCNA results
#   4. Creates comprehensive gene metrics matrix
#   5. Outputs results as CSV files
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
library(data.table)
library(dplyr)
library(tidyr)
library(WGCNA)
library(igraph)
library(yaml)
library(argparse)
})

# ----- 2. Source utilities -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")
source("src/utils/utils_gene_metrics.R")

# ----- 3. Set paths -----
config <- yaml::read_yaml("config/config.yaml")

# ----- Command-line arguments (override config values) -----
parser <- ArgumentParser(description = 'Calculate gene-level network metrics')
parser$add_argument('--adjacency-file', help = 'Path to soft-thresholded adjacency CSV (wide format).', default = config$network_feature_files$soft_threshold_files)
parser$add_argument('--output-dir', help = 'Directory to write gene metrics output.', default = file.path(config$output_dirs$network_features_dir, 'gene_metrics'))
parser$add_argument('--matrix-type', help = "Which matrix to process: 'both', 'adjacency', or 'spearman'", default = 'adjacency')
parser$add_argument('--threshold', type = 'double', help = 'Connection threshold for binary metrics (passed to calculate_gene_level_metrics)', default = 0.01)
args <- parser$parse_args()

adjacency_file <- args$adjacency_file
output_dir <- args$output_dir
matrix_type <- tolower(args$matrix_type)
threshold <- args$threshold
adjacency_results_dir <- file.path("results/network_features/features_calc/HS_adjacency")

create_directories(output_dir)

# ----- 4. Load matrices -----
cat("Loading matrices...\n")
adjacency <- read.csv(adjacency_file)

# Convert to proper matrix format
adjacency_matrix <- as.matrix(adjacency[,-1])
rownames(adjacency_matrix) <- adjacency[[1]]
colnames(adjacency_matrix) <- adjacency[[1]]

cat("Matrices loaded successfully:\n")
cat("- Genes in analysis:", nrow(adjacency_matrix), "\n")

# ----- 5. Calculate basic gene metrics for each matrix -----

cat("\n=== Calculating Basic Network Metrics ===\n")

adjacency_gene_metrics <- calculate_gene_level_metrics(
    matrix_data = adjacency_matrix,
    matrix_name = "Adjacency", 
    threshold   = threshold,
    matrix_type = "adjacency"
)

# ----- 6. Load module results and add module-based metrics -----

cat("\n=== Adding Module-Based Metrics ===\n")

# Try to load existing module results
adjacency_connectivity_file <- file.path(adjacency_results_dir, "gene_connectivity.csv")

if (file.exists(adjacency_connectivity_file)) {
    cat("Loading adjacency network module results...\n")
    adjacency_module_results <- read.csv(adjacency_connectivity_file)
    
    # Merge with basic metrics
    adjacency_gene_metrics <- merge(adjacency_gene_metrics,
                                   adjacency_module_results[, c("gene", "kWithin", "kOut", "kDiff", "module")],
                                   by = "gene", all.x = TRUE)
} else {
    cat("Adjacency module results not found. Run network_metrics script first.\n")
}

# ----- 7. Save comprehensive gene metrics -----

cat("\n=== Saving Results ===\n")

write.csv(adjacency_gene_metrics, 
          file.path(output_dir, "HS_adjacency_gene_metrics.csv"),
          row.names = FALSE)
    
save(adjacency_gene_metrics, file = file.path(output_dir, "HS_gene_metrics.RData"))

cat("Gene-level metrics saved to:", output_dir, "\n")
cat("Files created:\n")
cat("- HS_adjacency_gene_metrics.csv\n")
cat("- HS_gene_metrics.RData\n")

# ----- 8. Generate summary statistics -----
cat("\n=== Summary Statistics ===\n")

summarize_gene_metrics(adjacency_gene_metrics, "Adjacency Network")

cat("\nAnalysis complete.\n")
