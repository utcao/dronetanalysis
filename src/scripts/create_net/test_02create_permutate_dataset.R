# !/usr/bin/env Rscript
# ==============================================================================
# Trim unsignificant edge in permutation test before transcriptomic network analysis
#
# Script by Yu-tao Cao
#
# Last Updated: 08/10/2025
#
# This script::
#   1. create permutation datasets
#   2. trim connection accroding to the wilcox test
#   3. write out result table
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
    library(purrr)
    library(data.table)
    library(glue)
    library(stringr)
    library(ggplot2)
    library(tibble)  # for rownames_to_column / column_to_rownames
    library(futile.logger)
})

# ----- 1. Source Configuration and Utility Scripts -----
library(yaml)  # to read the YAML config
config <- yaml::read_yaml("config/config.yaml")

source("src/scripts/create_net/wrap_02create_permutate_dataset.R")
source("src/utils/utils_io.R")  # e.g., create_directories()
source("src/utils/utils_permutation_net.R")

# ----- 2. Extract Key Settings from config.yaml -----
expcor_tab_dir <- "results/spearman_correlation/"

expcor_tab_file <- file.path(expcor_tab_dir, "spearman_correlation_matrix.csv")

# ----- 3. Load data -----
coexp_expr_tab <- fread(expcor_tab_file)[1:3000, 1:3001]

# ----- 4. create permutation datasets by shuffle gene_id -----
sig_edge_tab_file <- "sig_edges_coexpr_net.csv"

sig_coexp_pairs <- permutation_test_stat_tab(coexp_expr_tab[1:300, 1:301],
                        sig_edge_tab_file,
                        expcor_tab_dir,
                        obs_col = "rho",
                        permut_cols_pattern = "seed",
                        permutation_num = 20,
                        alpha = 0.5)