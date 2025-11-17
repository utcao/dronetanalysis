# !/usr/bin/env Rscript
# ==============================================================================
# Spearman-based WGCNA for drosophila transcriptomic data
#
# Script by Gabriel Thornes
#
# Last Updated: 11/11/2025
#
# This script::
#   1. Takes subset datasets as input
#   2. Calculates Spearman rank correlations
#   4. Saves matrix outputs
# ==============================================================================

#################################
##### Packages and Setup ########
#################################
suppressPackageStartupMessages({
    library(argparse)
    library(WGCNA)
    library(yaml)
})
# ----- Parse command line arguments -----
parser <- ArgumentParser(description = 'Convert gene pairs to correlation matrix')
parser$add_argument('--input', type="character", 
                   default="/tmp/global2/gthornes/dronetanalysis/dataset/processed/VOOM/voomdataHS.txt",
                   help='Input file path (e.g. VOOM or VST processed data)')
parser$add_argument('--output', type="character",
                   default="/tmp/global2/gthornes/dronetanalysis/results/spearman_correlation/spearman_matrices", 
                   help='Output file path for Spearman correlation matrices')
parser$add_argument('--method', type="character", default="spearman",
                   help='Correlation method')
parser$add_argument('--use', type="character", default="complete.obs",
                   help='Value column name')

args <- parser$parse_args()

data_file <- args$input
out_file <- args$output
method <- args$method
use <- args$use

source("src/utils/utils_io.R")

# Load configuration
config <- yaml::read_yaml("config/config.yaml")

# Options for the analysis
options(stringsAsFactors = FALSE)

# Create output directories if they don't exist  
dir.create(file.path(config$output_dirs$spearman_dir, "spearman_matrices"), recursive = TRUE, showWarnings = FALSE)

###############################
##### Data Loading ############
###############################

# Load data
data <- read.table(data_file, header=TRUE, sep="\t", row.names=1)

# Transpose data for WGCNA format (samples as rows)
datExpr <- as.data.frame(t(data))

# Store original gene names before any filtering
original_gene_names <- colnames(datExpr)

# Check for genes and samples with too many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  # Remove the offending genes and samples
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  # Update the gene names list to match filtered data
  original_gene_names <- original_gene_names[gsg$goodGenes]
}

cat("Working with", nrow(datExpr), "samples and", ncol(datExpr), "genes\n")

#################################
##### Spearman Correlation ######
#################################

cat("\n=== SPEARMAN CORRELATION CALCULATION ===\n")

# Calculate Spearman correlation matrix
cat("Calculating Spearman correlation matrix...\n")
spearman_cor <- cor(datExpr, method = method, use = use)

# Save the correlation matrix
spearman_file <- file.path(out_file, "/HS_spearman_correlation_matrix.csv")
write.csv(spearman_cor, file = spearman_file)
subset_spearman_file <- file.path(out_file, "/HS_spearman_correlation_matrix_subset.csv")
write.csv(spearman_cor[1:500, 1:500], file = subset_spearman_file) # choose subset matrix size
cat("Spearman correlation matrix saved to:", spearman_file, "\n")
cat("Subset of Spearman correlation matrix (500x500) saved to:", subset_spearman_file, "\n")

# Basic correlation statistics
cat("Correlation statistics:\n")
cat("Mean absolute correlation:", mean(abs(spearman_cor[upper.tri(spearman_cor)])), "\n")
cat("Median absolute correlation:", median(abs(spearman_cor[upper.tri(spearman_cor)])), "\n")
