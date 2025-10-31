#!/usr/bin/env Rscript
# ==============================================================================
# Calculate comprehensive gene-level network metrics
#
# Script by Gabriel Thornes
#
# Last Updated: 31/10/2025
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
})

# ----- 2. Source utilities -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")
source("src/utils/utils_gene_metrics.R")

# ----- 3. Set paths -----
config <- yaml::read_yaml("config/config.yaml")
spearman_input_file <- file.path(config$spearman_correlation_files$permutation_files, "sig_matrix_wide.csv")
adjacency_file <- config$network_feature_files$soft_threshold_files

# Module results directory
adjacency_results_dir <- file.path(config$output_dirs$network_features_dir, "features_calc/adjacency")
output_dir <- file.path(config$output_dirs$network_features_dir, "gene_metrics")

create_directories(output_dir)

# ----- 4. Load matrices -----
cat("Loading matrices...\n")
spearman <- read.csv(spearman_input_file)
adjacency <- read.csv(adjacency_file)

# Convert to proper matrix format
spearman_matrix <- as.matrix(spearman[,-1])
rownames(spearman_matrix) <- spearman[[1]]
colnames(spearman_matrix) <- spearman[[1]]

adjacency_matrix <- as.matrix(adjacency[,-1])
rownames(adjacency_matrix) <- adjacency[[1]]
colnames(adjacency_matrix) <- adjacency[[1]]

cat("Matrices loaded successfully:\n")
cat("- Genes in analysis:", nrow(spearman_matrix), "\n")

# ----- 5. Calculate basic gene metrics for each matrix -----

cat("\n=== Calculating Basic Network Metrics ===\n")

# Calculate for both matrices
spearman_gene_metrics <- calculate_gene_level_metrics(
    matrix_data = spearman_matrix,
    matrix_name = "Spearman",
    threshold   = 0.2,
    matrix_type = "correlation"
)

adjacency_gene_metrics <- calculate_gene_level_metrics(
    matrix_data = adjacency_matrix,
    matrix_name = "Adjacency", 
    threshold   = 0.2,
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

write.csv(spearman_gene_metrics, 
          file.path(output_dir, "spearman_gene_metrics.csv"), 
          row.names = FALSE)

write.csv(adjacency_gene_metrics, 
          file.path(output_dir, "adjacency_gene_metrics.csv"),
          row.names = FALSE)

# Create combined summary
combined_metrics <- list(
    spearman = spearman_gene_metrics,
    adjacency = adjacency_gene_metrics
)

save(combined_metrics, file = file.path(output_dir, "all_gene_metrics.RData"))

cat("Gene-level metrics saved to:", output_dir, "\n")
cat("Files created:\n")
cat("- spearman_gene_metrics.csv\n")
cat("- adjacency_gene_metrics.csv\n")
cat("- all_gene_metrics.RData\n")

# ----- 8. Generate summary statistics -----
cat("\n=== Summary Statistics ===\n")

summarize_gene_metrics(spearman_gene_metrics, "Spearman Correlation Network")
summarize_gene_metrics(adjacency_gene_metrics, "Adjacency Network")

cat("\nAnalysis complete.\n")
