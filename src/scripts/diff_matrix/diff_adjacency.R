#!/usr/bin/env Rscript
# ==============================================================================
# Calculate Differential Adjacency Matrix
#
# Script by Gabriel Thornes
#
# Last Updated: 26/11/2025
#
# This script:
#   1. Loads two adjacency matrices
#   2. Calculates the absolute difference between them
#   3. Saves the differential adjacency matrix
# ==============================================================================

rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
})

# Command-line arguments
parser <- ArgumentParser(description = 'Calculate differential adjacency matrix')
parser$add_argument('--adj-1st', help = 'Path to first adjacency matrix',
                   required = TRUE)
parser$add_argument('--adj-5th', help = 'Path to second adjacency matrix',
                   required = TRUE)
parser$add_argument('--output', help = 'Path to save differential adjacency matrix',
                   required = TRUE)
parser$add_argument('--threshold', type = 'double', help = 'Hard threshold to apply before differencing',
                   default = 0.1)
args <- parser$parse_args()

# Load first adjacency matrix
adj_1st <- fread(args$adj_1st, data.table = FALSE)
gene_names <- as.character(adj_1st[[1]])
adj_1st_mat <- as.matrix(adj_1st[, -1])
rownames(adj_1st_mat) <- gene_names
colnames(adj_1st_mat) <- gene_names

# Load second adjacency matrix
adj_5th <- fread(args$adj_5th, data.table = FALSE)
adj_5th_mat <- as.matrix(adj_5th[, -1])
rownames(adj_5th_mat) <- gene_names
colnames(adj_5th_mat) <- gene_names

# Apply hard thresholding to both matrices
hard_threshold <- args$threshold
adj_1st_mat[adj_1st_mat < hard_threshold] <- 0
adj_5th_mat[adj_5th_mat < hard_threshold] <- 0

# Calculate absolute difference
diff_adj <- abs(adj_5th_mat - adj_1st_mat)

# Save result
diff_adj_df <- as.data.frame(diff_adj)
diff_adj_df <- cbind(gene = rownames(diff_adj), diff_adj_df)
write.csv(diff_adj_df, file = args$output, row.names = FALSE)

cat("Differential adjacency matrix saved to:", args$output, "\n")