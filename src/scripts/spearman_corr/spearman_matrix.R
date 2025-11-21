#!/usr/bin/env Rscript
# ==============================================================================
# Spearman Correlation Matrix Calculation
#
# Script by Gabriel Thornes
#
# Last Updated: 21/11/2025
#
# This script:
#   1. Takes expression matrix as input (genes × samples)
#   2. Calculates Spearman rank correlations
#   3. Saves correlation matrix output
# ==============================================================================

#################################
##### Packages and Setup ########
#################################
suppressPackageStartupMessages({
    library(argparse)
    library(WGCNA)
    library(data.table)
})

# ----- Parse command line arguments -----
parser <- ArgumentParser(description = 'Calculate Spearman correlation matrix')
parser$add_argument('--input', type="character", required=TRUE,
                   help='Input expression matrix file (genes x samples)')
parser$add_argument('--output', type="character", required=TRUE,
                   help='Output file path for Spearman correlation matrix')
parser$add_argument('--method', type="character", default="spearman",
                   help='Correlation method (default: spearman)')
parser$add_argument('--use', type="character", default="pairwise.complete.obs",
                   help='Handling of missing values (default: pairwise.complete.obs)')

args <- parser$parse_args()

# Options for the analysis
options(stringsAsFactors = FALSE)

# Create output directory if it doesn't exist
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)

cat("=== Spearman Correlation Matrix Calculation ===\n")
cat("Input file:", args$input, "\n")
cat("Output file:", args$output, "\n")
cat("Method:", args$method, "\n")
cat("Missing value handling:", args$use, "\n\n")

###############################
##### Data Loading ############
###############################

cat("Loading expression data...\n")

# Load data using fread for better handling
data <- fread(args$input, data.table = FALSE)

# Check if first column needs renaming
if (colnames(data)[1] %in% c("V1", "")) {
  colnames(data)[1] <- "gene_id"
}

# Set first column as row names
rownames(data) <- data[[1]]
data <- data[, -1]

cat("  Loaded:", nrow(data), "genes x", ncol(data), "samples\n")

# Transpose data for correlation calculation (samples as rows, genes as columns)
datExpr <- as.data.frame(t(data))

cat("  Transposed to:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

# Check for genes and samples with too many missing values
cat("\nChecking data quality...\n")
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  cat("  Removing genes/samples with too many missing values\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("  After filtering:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")
} else {
  cat("  All samples and genes pass quality check\n")
}

#################################
##### Spearman Correlation ######
#################################

cat("\n=== CALCULATING SPEARMAN CORRELATION ===\n")

# Calculate Spearman correlation matrix
cat("Computing correlation matrix (", args$method, ")...\n")
spearman_cor <- cor(datExpr, method = args$method, use = args$use)

cat("  Correlation matrix dimensions:", nrow(spearman_cor), "x", ncol(spearman_cor), "\n")

# Basic correlation statistics
cat("\nCorrelation statistics:\n")
upper_tri_vals <- spearman_cor[upper.tri(spearman_cor)]
cat("  Mean absolute correlation:", round(mean(abs(upper_tri_vals)), 4), "\n")
cat("  Median absolute correlation:", round(median(abs(upper_tri_vals)), 4), "\n")
cat("  Min correlation:", round(min(upper_tri_vals), 4), "\n")
cat("  Max correlation:", round(max(upper_tri_vals), 4), "\n")

# Save the correlation matrix
cat("\nSaving correlation matrix...\n")
write.csv(spearman_cor, file = args$output, row.names = TRUE)
cat("  Saved to:", args$output, "\n")

cat("\n=== CORRELATION CALCULATION COMPLETE ===\n")
