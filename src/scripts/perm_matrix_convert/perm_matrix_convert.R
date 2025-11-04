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

# ----- 2. Set Input and Output Paths from config.yaml -----

pval_tab_dir <- "results/spearman_correlation/perm/permutation_test"
pval_tab_file <- file.path(pval_tab_dir, "sig_edges_coexpr_net.csv")
out_file <- file.path(pval_tab_dir, "sig_matrix_wide.csv")

pval_tab <- fread(pval_tab_file)

cat("Original data dimensions:", dim(pval_tab), "\n")

# Create matrix from original significant edges (p-value based)
square_correlation_matrix <- long_to_square_dcast(pval_tab, gene_pairs_col = "gene_pairs", value_col = "rho")
square_correlation_matrix[square_correlation_matrix == 0] <- NA  # Set zero correlations to NA

# Save matrix
write.csv(square_correlation_matrix, out_file, row.names = TRUE)

cat("Square matrix saved to:", out_file, "\n")
cat("Matrix dimensions:", dim(square_correlation_matrix), "\n")
