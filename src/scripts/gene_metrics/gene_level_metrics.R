#!/usr/bin/env Rscript
# ==============================================================================
# Calculate comprehensive gene-level network metrics
#
# Script by Gabriel Thornes
#
# Last Updated: 27/10/2025
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
library(data.table)
library(dplyr)
library(tidyr)
library(WGCNA)
library(igraph)
library(yaml)

# ----- 2. Source utilities -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")
source("src/utils/utils_gene_metrics.R")

# ----- 3. Set paths -----
config <- yaml::read_yaml("config/config.yaml")
spearman_input_file <- file.path(config$spearman_correlation_files$permutation_files, "sig_matrix_wide.csv")
signed_adj_file <- config$network_feature_files$soft_threshold_files[1]
unsigned_adj_file <- config$network_feature_files$soft_threshold_files[2]

# Module results directories
signed_results_dir <- file.path(config$output_dirs$network_features_dir, "features_calc/signed")
unsigned_results_dir <- file.path(config$output_dirs$network_features_dir, "features_calc/unsigned")
output_dir <- file.path(config$output_dirs$network_features_dir, "gene_metrics")

create_directories(output_dir)

# ----- 4. Load matrices -----
cat("Loading matrices...\n")
spearman <- read.csv(spearman_input_file)
signed_adj <- read.csv(signed_adj_file)
unsigned_adj <- read.csv(unsigned_adj_file)

signed_adj<- signed_adj[1:1000,1:1001]
unsigned_adj<- unsigned_adj[1:1000,1:1001]

# Convert to proper matrix format
spearman_matrix <- as.matrix(spearman[,-1])
rownames(spearman_matrix) <- spearman[[1]]
colnames(spearman_matrix) <- spearman[[1]]

signed_matrix <- as.matrix(signed_adj[,-1])  
rownames(signed_matrix) <- signed_adj[[1]]
colnames(signed_matrix) <- signed_adj[[1]]

unsigned_matrix <- as.matrix(unsigned_adj[,-1])
rownames(unsigned_matrix) <- unsigned_adj[[1]]
colnames(unsigned_matrix) <- unsigned_adj[[1]]

cat("Matrices loaded successfully:\n")
cat("- Genes in analysis:", nrow(spearman_matrix), "\n")

# ----- 5. Calculate basic gene metrics for each matrix -----

cat("\n=== Calculating Basic Network Metrics ===\n")

# Calculate for all three matrices
spearman_gene_metrics <- calculate_gene_level_metrics(
    matrix_data = spearman_matrix,
    matrix_name = "Spearman_Correlation",
    threshold   = 0.5,
    matrix_type = "correlation"
)

signed_gene_metrics <- calculate_gene_level_metrics(
    matrix_data = signed_matrix, 
    matrix_name = "Signed_Adjacency",
    threshold   = 0.1,
    matrix_type = "adjacency"
)

unsigned_gene_metrics <- calculate_gene_level_metrics(
    matrix_data = unsigned_matrix,
    matrix_name = "Unsigned_Adjacency", 
    threshold   = 0.2,
    matrix_type = "adjacency"
)

# ----- 6. Load module results and add module-based metrics -----

cat("\n=== Adding Module-Based Metrics ===\n")

# Try to load existing module results
signed_connectivity_file <- file.path(signed_results_dir, "signed_gene_connectivity.csv")
unsigned_connectivity_file <- file.path(unsigned_results_dir, "unsigned_gene_connectivity.csv")

if (file.exists(signed_connectivity_file)) {
    cat("Loading signed network module results...\n")
    signed_module_results <- read.csv(signed_connectivity_file)
    
    # Merge with basic metrics
    signed_gene_metrics <- merge(signed_gene_metrics, 
                                signed_module_results[, c("gene", "kWithin", "kOut", "kDiff", "module")], 
                                by = "gene", all.x = TRUE)
} else {
    cat("Signed module results not found. Run network_metrics script first.\n")
}

if (file.exists(unsigned_connectivity_file)) {
    cat("Loading unsigned network module results...\n") 
    unsigned_module_results <- read.csv(unsigned_connectivity_file)
    
    # Merge with basic metrics
    unsigned_gene_metrics <- merge(unsigned_gene_metrics,
                                  unsigned_module_results[, c("gene", "kWithin", "kOut", "kDiff", "module")],
                                  by = "gene", all.x = TRUE)
} else {
    cat("Unsigned module results not found. Run network_metrics script first.\n")
}

# ----- 7. Save comprehensive gene metrics -----

cat("\n=== Saving Results ===\n")

write.csv(spearman_gene_metrics, 
          file.path(output_dir, "spearman_gene_metrics.csv"), 
          row.names = FALSE)

write.csv(signed_gene_metrics,
          file.path(output_dir, "signed_gene_metrics.csv"),
          row.names = FALSE)

write.csv(unsigned_gene_metrics, 
          file.path(output_dir, "unsigned_gene_metrics.csv"),
          row.names = FALSE)

# Create combined summary
combined_metrics <- list(
    spearman = spearman_gene_metrics,
    signed = signed_gene_metrics, 
    unsigned = unsigned_gene_metrics
)

save(combined_metrics, file = file.path(output_dir, "all_gene_metrics.RData"))

cat("Gene-level metrics saved to:", output_dir, "\n")
cat("Files created:\n")
cat("- spearman_gene_metrics.csv\n")
cat("- signed_gene_metrics.csv\n") 
cat("- unsigned_gene_metrics.csv\n")
cat("- all_gene_metrics.RData\n")

# ----- 8. Generate summary statistics -----
cat("\n=== Summary Statistics ===\n")

summarize_gene_metrics(spearman_gene_metrics, "Spearman Correlation Network")
summarize_gene_metrics(signed_gene_metrics, "Signed Adjacency Network") 
summarize_gene_metrics(unsigned_gene_metrics, "Unsigned Adjacency Network")

cat("\nAnalysis complete.\n")
