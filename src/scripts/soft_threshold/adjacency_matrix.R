# !/usr/bin/env Rscript
# ==============================================================================
# Create adjacency matrices from correlation matrix
#
# Script by Gabriel Thornes
#
# Last Updated: 27/10/2025
#
# This script::
#   1. Takes spearman correlation matrix as input
#   2. Calculates adjacency matrices
#   3. Generates plots to decide soft-thresholding power
#   4. Outputs adjacency matrices, plots and fit indices as CSV files in appropriate folder
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
matrix_dir <- config$output_dirs$spearman_dir
matrix_file <- file.path(matrix_dir, "spearman_correlation_matrix.csv")
output_dir <- config$output_dirs$network_features_dir
signed_output_file <- file.path(output_dir, "soft_threshold/signed/signed_adjacency_matrix.csv")
unsigned_output_file <- file.path(output_dir, "soft_threshold/unsigned/unsigned_adjacency_matrix.csv")

corr_matrix <- fread(matrix_file)
cat("File read:", matrix_file,"\n")
df.corr.m <- as.matrix(corr_matrix[,-1, with=FALSE])
    rownames(df.corr.m) <- corr_matrix[[1]]

# ----- 3. Calculate adjacency matrices and generate soft thresholding plots -----

# Option 1: Signed network (preserves directionality)
# Convert correlations to signed adjacency (keeps positive/negative relationships)
# Values range from -1 to +1, then transformed to 0 to 1 as: (1 + cor) / 2
signed_adjacency <- (1 + df.corr.m) / 2

# Option 2: Traditional unsigned network (magnitude only)
unsigned_adjacency <- abs(df.corr.m) # adjacency matrix

# Save the adjacency matrices
write.csv(signed_adjacency, signed_output_file, row.names = TRUE)
write.csv(unsigned_adjacency, unsigned_output_file, row.names = TRUE)
cat("Signed adjacency matrix saved to:", signed_output_file, "\n")
cat("Unsigned adjacency matrix saved to:", unsigned_output_file, "\n")

# -----4. Generate soft thresholding plots to decide on appropriate power -----

# Generate soft-thresholded matrices and plots

##############################
## Option 1: Signed network ##
##############################

# Call the network topology analysis function for signed network
sft_signed <- pickSoftThreshold(signed_adjacency, powerVector = c(1:30), 
                                networkType = "signed", verbose = 5)

# Generate plots to decide soft-thresholding power
s_output_plot <- file.path(output_dir, "soft_threshold/signed/soft_thresholding_signed.pdf")
sft_plot(sft_signed, s_output_plot, c(1:30))

# Write sft to file for analysis
signed_output_file <- file.path(output_dir, "soft_threshold/signed/signed_soft_threshold.csv")
write.csv(sft_signed$fitIndices, file = signed_output_file, row.names = FALSE)

################################
## Option 2: Unsigned network ##
################################

# Call the network topology analysis function for unsigned network
sft_unsigned <- pickSoftThreshold(unsigned_adjacency, powerVector = c(1:20), # lower power range for unsigned
                                  networkType = "unsigned", verbose = 5)

# Generate plots to decide soft-thresholding power
u_output_plot <- file.path(output_dir, "soft_threshold/unsigned/soft_thresholding_unsigned.pdf")
sft_plot(sft_unsigned, u_output_plot, c(1:20))

# Write sft to file for analysis
unsigned_output_file <- file.path(output_dir, "soft_threshold/unsigned/unsigned_soft_threshold.csv")
write.csv(sft_unsigned$fitIndices, file = unsigned_output_file, row.names = FALSE)

############################################################################################
#### Please inspect plots and output file to decide on appropriate soft threshold power ####
############################################################################################