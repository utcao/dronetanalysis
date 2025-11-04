# !/usr/bin/env Rscript
# ==============================================================================
# Apply soft threshold to adjacency matrices
#
# Script by Gabriel Thornes
#
# Last Updated: 30/10/2025
#
# This script::
#   1. Takes adjacency matrices as input
#   2. Applies soft-thresholding power based on plot analysis
#   3. Outputs soft threshold-transformed matrices as CSV files in appropriate folder
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
library(argparse)
library(data.table)
library(dplyr)
library(tidyr)
library(WGCNA)
library(yaml)

# ----- Parse command line arguments -----
parser <- ArgumentParser(description = 'Apply soft threshold to adjacency matrices')
parser$add_argument('--input', type="character", 
                   default="results/network_features/soft_threshold/unsigned/unsigned_adjacency_matrix.csv",
                   help='Input adjacency matrix file path')
parser$add_argument('--output-dir', type="character",
                   default="results/network_features", 
                   help='Output directory for soft-thresholded matrices')
parser$add_argument('--soft-power', type="integer", default=4,
                   help='Soft-thresholding power to apply')

args <- parser$parse_args()

# ----- 2. Set Input and Output Paths -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")

# Use command line arguments if provided, otherwise fall back to config
if (!is.null(args$input) && args$input != "results/network_features/soft_threshold/unsigned/unsigned_adjacency_matrix.csv") {
    soft_threshold_file <- args$input
    output_dir <- args$output_dir
} else {
    # Fall back to config file
    config <- yaml::read_yaml("config/config.yaml")
    soft_threshold_dir <- config$output_dirs$network_features_dir
    soft_threshold_file <- file.path(soft_threshold_dir, "soft_threshold/unsigned/unsigned_adjacency_matrix.csv")
    output_dir <- config$output_dirs$network_features_dir
}

cat("Input file:", soft_threshold_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Soft power:", args$soft_power, "\n")

# ----- 3. Load adjacency matrices -----

sft_unsigned <- read.csv(soft_threshold_file, row.names = 1)
cat("File read:", soft_threshold_file,"\n")

# ----- 4. Apply soft-thresholding -----

## Apply soft-thresholding power from command line argument ##

soft_power <- args$soft_power
soft_adj <- sft_unsigned^soft_power
output_file <- file.path(output_dir, "soft_threshold/unsigned/unsigned_soft_thresholded_adjacency_matrix.csv")
write.csv(soft_adj, output_file, row.names = TRUE)
cat("Soft-thresholded adjacency matrix saved to:", output_file, "\n")
