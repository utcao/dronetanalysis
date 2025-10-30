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

pval_tab_dir <- "results/spearman_correlation/permutation_test/"
pval_tab_file <- file.path(pval_tab_dir, "sig_edges_coexpr_net.csv")
out_dir <- "results/spearman_correlation/permutation_test/"
out_file <- file.path(out_dir, "sig_matrix_wide.csv")
out_file_fdr <- file.path(out_dir, "sig_matrix_wide_fdr005.csv")

pval_tab <- fread(pval_tab_file)

cat("Original data dimensions:", dim(pval_tab), "\n")
cat("FDR range:", range(pval_tab$q_twotail_fdr, na.rm = TRUE), "\n")

# Filter by FDR < 0.05
pval_tab_fdr <- pval_tab[q_twotail_fdr < 0.05]

cat("After FDR < 0.05 filtering:", dim(pval_tab_fdr), "\n")
cat("Percentage of edges retained:", round(100 * nrow(pval_tab_fdr) / nrow(pval_tab), 2), "%\n")

# Create matrix from original significant edges (p-value based)
square_correlation_matrix <- long_to_square_dcast(pval_tab, gene_pairs_col = "gene_pairs", value_col = "rho")

# Create matrix from FDR-filtered edges
square_correlation_matrix_fdr <- long_to_square_dcast(pval_tab_fdr, gene_pairs_col = "gene_pairs", value_col = "rho")

# Save both matrices
write.csv(square_correlation_matrix, out_file, row.names = TRUE)
write.csv(square_correlation_matrix_fdr, out_file_fdr, row.names = TRUE)

cat("Square matrix saved to:", out_file, "\n")
cat("Matrix dimensions:", dim(square_correlation_matrix), "\n")

cat("FDR-filtered matrix saved to:", out_file_fdr, "\n")
cat("FDR-filtered matrix dimensions:", dim(square_correlation_matrix_fdr), "\n")

# Summary statistics
non_na_original <- sum(!is.na(square_correlation_matrix))
non_na_fdr <- sum(!is.na(square_correlation_matrix_fdr))

cat("\n=== Filtering Summary ===\n")
cat("Original p-value filtered edges:", non_na_original, "\n")
cat("FDR < 0.05 filtered edges:", non_na_fdr, "\n")
cat("Reduction due to FDR filtering:", round(100 * (1 - non_na_fdr/non_na_original), 2), "%\n")

if (non_na_fdr > 0) {
    fdr_range <- range(square_correlation_matrix_fdr, na.rm = TRUE)
    cat("FDR-filtered correlation range:", round(fdr_range, 3), "\n")
}
