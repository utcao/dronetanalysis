# !/usr/bin/env Rscript
# ==============================================================================
# Convert long p values of gene pairs matrix to wide matrix
#
# Script by Gabriel Thornes
#
# Last Updated: 15/10/2025
#
# This script::
#   1. Takes p values of gene pairs matrix datasets as input
#   2. Outputs p x p matrix
# ==============================================================================

rm(list = ls())

source("src/utils/utils_io.R")
source("src/utils/utils_permutation_net.R")

# Load required packages
library(data.table)
library(dplyr)
library(tidyr)

# ----- 2. Set Input and Output Paths from config.yaml -----

pval_tab_dir <- "results/spearman_correlation/permutation_test/"
pval_tab_file <- file.path(pval_tab_dir, "sig_edges_coexpr_net.csv")
out_dir <- "results/spearman_correlation/permutation_test/"
out_file <- file.path(out_dir, "sig_matrix_wide.csv")

pval_tab <- fread(pval_tab_file)

square_correlation_matrix <- long_to_square_dcast(pval_tab, gene_pairs_col = "gene_pairs", value_col = "rho")

# Set zeros to NA
square_correlation_matrix[square_correlation_matrix == 0] <- NA

# Save the matrix
write.csv(square_correlation_matrix, out_file, row.names = TRUE)

cat("Square matrix saved to:", out_file, "\n")
cat("Matrix dimensions:", dim(square_correlation_matrix), "\n")
