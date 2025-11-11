# !/usr/bin/env Rscript
# ==============================================================================
# Create adjacency matrices from correlation matrix
#
# Script by Gabriel Thornes
#
# Last Updated: 04/11/2025
#
# This script::
#   1. Takes spearman correlation matrix as input
#   2. Creates binary matrix to preserve negative correlations
#   3. Calculates adjacency matrices
#   4. Generates plots to decide soft-thresholding power
#   5. Outputs adjacency matrices, plots and fit indices as CSV files in appropriate folder
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
library(argparse)
library(data.table)
library(dplyr)
library(tidyr)
library(WGCNA)
library(yaml)
})

# ----- Parse command line arguments -----
parser <- ArgumentParser(description = 'Create adjacency matrices from correlation matrix')
parser$add_argument('--input', type="character", 
                   default="results/spearman_correlation/spearman_matrices/HS_spearman_correlation_matrix.csv",
                   help='Input correlation matrix file path')
parser$add_argument('--output-dir', type="character",
                   default="results/network_features", 
                   help='Output directory for adjacency matrices and plots')
parser$add_argument('--power-range', type="character", default="1:20",
                   help='Power range for soft thresholding (e.g. "1:20")')
parser$add_argument('--network-type', type="character", default="unsigned",
                   help='Network type (unsigned or signed)')

args <- parser$parse_args()

# Parse power range
power_range <- eval(parse(text = args$power_range))

# ----- 2. Set Input and Output Paths -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")

# Use command line arguments if provided, otherwise fall back to config
if (!is.null(args$input) && args$input != "results/spearman_correlation/spearman_matrices/HS_spearman_correlation_matrix.csv") {
    matrix_file <- args$input
    output_dir <- args$output_dir
} else {
    # Fall back to config file
    config <- yaml::read_yaml("config/config.yaml")
    matrix_dir <- config$output_dirs$spearman_dir
    matrix_file <- file.path(matrix_dir, "spearman_matrices/HS_spearman_correlation_matrix.csv")
    output_dir <- config$output_dirs$network_features_dir
}

cat("Input file:", matrix_file, "\n")
cat("Output directory:", output_dir, "\n")

binary_output_file <- file.path(output_dir, "binary_matrix/HS_binary_signed_matrix.csv")
unsigned_output_file <- file.path(output_dir, "soft_threshold/HS_unsigned_adjacency_matrix.csv")

corr_matrix <- fread(matrix_file)
cat("File read:", matrix_file,"\n")
df.corr.m <- as.matrix(corr_matrix[,-1, with=FALSE])
    rownames(df.corr.m) <- corr_matrix[[1]]

# ----- 3. Create binary matrix to preserve negative correlations -----
create_directories(file.path(output_dir, "binary_matrix/"))
binary_matrix <- create_correlation_sign_matrix(df.corr.m, output_file = binary_output_file)
cat("Binary sign matrix created.\n")

# ----- 4. Calculate adjacency matrices and generate soft thresholding plots -----
adjacency <- abs(df.corr.m) # adjacency matrix

# Save the adjacency matrices
write.csv(adjacency, unsigned_output_file, row.names = TRUE)
cat("Adjacency matrix saved to:", unsigned_output_file, "\n")

# ----- 5. Generate soft thresholding plots to decide on appropriate power -----

# Generate soft-thresholded matrices and plots

# Call the network topology analysis function for unsigned network
sft_unsigned <- pickSoftThreshold(adjacency, powerVector = power_range,
                                  networkType = args$network_type, verbose = 5)

# Generate plots to decide soft-thresholding power
u_output_plot <- file.path(output_dir, "soft_threshold/HS_soft_thresholding.pdf")
sft_plot(sft_unsigned, u_output_plot, power_range)

# Write sft to file for analysis
unsigned_output_file <- file.path(output_dir, "soft_threshold/HS_soft_threshold.csv")
write.csv(sft_unsigned$fitIndices, file = unsigned_output_file, row.names = FALSE)

############################################################################################
#### Please inspect plots and output file to decide on appropriate soft threshold power ####
############################################################################################