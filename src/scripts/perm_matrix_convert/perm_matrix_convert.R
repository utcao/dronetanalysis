# !/usr/bin/env Rscript
# ==============================================================================
# Convert long p values of gene pairs matrix to wide matrix
#
# Script by Gabriel Thornes
#
# Last Updated: 04/11/2025
#
# This script::
#   1. Takes p values of gene pairs matrix datasets as input
#   2. Filters by FDR < 0.05
#   3. Outputs p x p matrix
# ==============================================================================

rm(list = ls())

source("src/utils/utils_io.R")
source("src/utils/utils_permutation_net.R")

# Load required packages
library(data.table)
library(dplyr)
library(tidyr)
library(argparse)

# ----- 1. Add simple command line argument parsing -----
parser <- ArgumentParser(description = 'Convert gene pairs to correlation matrix')
parser$add_argument('--input', type="character", 
                   default="results/spearman_correlation/perm/permutation_test/sig_edges_coexpr_net.csv",
                   help='Input file path')
parser$add_argument('--output', type="character",
                   default="results/spearman_correlation/perm/permutation_test/sig_matrix_wide.csv", 
                   help='Output file path')
parser$add_argument('--gene-pairs-col', type="character", default="gene_pairs",
                   help='Gene pairs column name')
parser$add_argument('--value-col', type="character", default="rho",
                   help='Value column name')

args <- parser$parse_args()

# ----- 2. Use arguments instead of hardcoded paths -----
pval_tab_file <- args$input
out_file <- args$output
gene_pairs_col <- args$gene_pairs_col
value_col <- args$value_col

cat("Input file:", pval_tab_file, "\n")
cat("Output file:", out_file, "\n")

pval_tab <- fread(pval_tab_file)

cat("Original data dimensions:", dim(pval_tab), "\n")

# Create matrix from original significant edges (p-value based)
square_correlation_matrix <- long_to_square_dcast(pval_tab, 
                                                  gene_pairs_col = gene_pairs_col, 
                                                  value_col = value_col)
square_correlation_matrix[square_correlation_matrix == 0] <- NA  # Set zero correlations to NA

# Save matrix
write.csv(square_correlation_matrix, out_file, row.names = TRUE)

cat("Square matrix saved to:", out_file, "\n")
cat("Matrix dimensions:", dim(square_correlation_matrix), "\n")
