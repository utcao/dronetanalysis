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
library(data.table)
library(dplyr)
library(tidyr)
library(WGCNA)
library(yaml)

# ----- 2. Set Input and Output Paths from config.yaml -----
source("src/utils/utils_io.R")
source("src/utils/utils_network_feats.R")
config <- yaml::read_yaml("config/config.yaml")
soft_threshold_dir <- config$output_dirs$network_features_dir
soft_threshold_file <- file.path(soft_threshold_dir, "soft_threshold/unsigned/unsigned_adjacency_matrix.csv")
output_dir <- config$output_dirs$network_features_dir

# ----- 3. Load adjacency matrices -----

sft_unsigned <- read.csv(soft_threshold_file, row.names = 1)
cat("File read:", soft_threshold_file,"\n")

# ----- 4. Apply soft-thresholding -----

## Apply soft-thresholding power based on plot analysis ##

soft_power <- 4  ### Example power, adjust based on plot analysis ###
soft_adj <- sft_unsigned^soft_power
output_file <- file.path(output_dir, "soft_threshold/unsigned/unsigned_soft_thresholded_adjacency_matrix.csv")
write.csv(soft_adj, output_file, row.names = TRUE)
cat("Unsigned soft-thresholded adjacency matrix saved to:", output_file, "\n")
