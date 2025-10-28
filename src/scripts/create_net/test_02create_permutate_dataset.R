#!/usr/bin/env Rscript
# ==============================================================================
# Co-expression Network Filtering via Permutation Testing
#
# Author: Yu-tao Cao
# Last Updated: 23/10/2025
#
# Description:
#   This script performs statistical filtering of co-expression networks using
#   permutation testing to identify biologically significant gene-gene interactions.
#   By comparing observed correlation patterns against null distributions generated
#   through random gene label shuffling, it removes spurious correlations that
#   may arise by chance, enhancing the biological relevance of the resulting network.
#
# Key Functionality:
#   1. Generates permutation datasets by shuffling gene identifiers
#   2. Performs statistical testing (Wilcoxon test) to identify significant edges
#   3. Filters co-expression network based on permutation p-values
#   4. Outputs validated network edges and visualization plots
#   5. Provides quality control through distribution visualization
#
# Input Requirements:
#   - Spearman correlation matrix (genes x genes)
#   - Configuration file (config.yaml) with analysis parameters
#
# Output Products:
#   - Filtered co-expression network edges (CSV format)
#   - Permutation test results and statistics
#   - Distribution plots for selected gene pairs
#
# Dependencies: purrr, data.table, glue, stringr, ggplot2, tibble, futile.logger, yaml
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
expcor_tab_dir <- "results/spearman_correlation"

expcor_tab_file <- file.path(expcor_tab_dir, "spearman_correlation_matrix.csv")

# ----- 3. Load data -----
coexp_expr_tab <- fread(expcor_tab_file)[1:1000,1:1001]

# ----- 4. create permutation datasets by shuffle gene_id -----
sig_edge_tab_file <- "sig_edges_coexpr_net.csv"

setnames(coexp_expr_tab, "V1", "gene_id")
sig_coexp_pairs <- permutation_test_stat(
                    coexp_expr_tab = coexp_expr_tab,
                    sig_edge_tab_file = sig_edge_tab_file,
                    tab_output_dir = expcor_tab_dir,
                    obs_col = "rho",
                    pair_id_col = "gene_pairs",
                    permu_cols_pattern = "seed",
                    permu_n = 300,
                    alpha = 0.05, keep_permu = TRUE)

id_col <- "gene_pairs"
permu_cols_pat <- "seed"
obs_col <- "rho"
plot_output_dir <- file.path("results/spearman_correlation/permutation_test/plots")

# plot
sig_coexp_pairs_ltab <- sig_coexp_pairs[, .SD, .SDcols = patterns(glue("{id_col}|{permu_cols_pat}|{obs_col}"))] |>
                melt(id.vars = id_col, value.name = obs_col, variable.name = "datasets")

coexpr_pairs <- sig_coexp_pairs_ltab[, unique(get(id_col))]
choose_pairs_n <- 3  # Number of pairs to plot
choose_pairs <- sample(coexpr_pairs, choose_pairs_n)

map(choose_pairs, ~ plot_permu_distr(
            coexp_pairs_ltab = sig_coexp_pairs_ltab,
            id_col = "gene_pairs",
            chose_pair = .x,
            val_col = "rho", var_col = "datasets",
            plot_output_dir = plot_output_dir, prefix = 1))

# # Check if observed correlations look reasonable
cat("Summary of observed correlation coefficients (rho):\n")
print(summary(sig_coexp_pairs$rho))
