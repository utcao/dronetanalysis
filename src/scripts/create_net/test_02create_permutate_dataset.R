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

# If you have utility scripts (optional):
source("src/utils/utils_io.R")  # e.g., create_directories()
source("src/utils/utils_permutation_net.R")

# ----- 2. Extract Key Settings from config.yaml -----
expcor_tab_dir <- "results/spearman_correlation/"
permutate_response_dir <- file.path(expcor_tab_dir, "shuffle")
permutate_res_tab_dir <- file.path(expcor_tab_dir, "permutation_test")

create_directories(c(permutate_response_dir, permutate_res_tab_dir))

expcor_tab_file <- file.path(expcor_tab_dir, "spearman_correlation_matrix.csv")
sig_edge_tab_file <- file.path(permutate_res_tab_dir,
                                "sig_edges_coexpr_net.csv")

# ----- 3. Load data -----
coexp_expr_tab <- fread(expcor_tab_file)[1:3000, 1:3001]

# ----- 4. create permutation datasets by shuffle gene_id -----
# or consider split shuffle and permutation test into separate scripts
permutation_test_stat_tab <- function(coexp_expr_tab, sig_edge_tab_file,
                                    permutation_num = 50){
    
    permute_coexp_list <- permutate_col(coexp_expr_tab, "gene_id", permutation_num,
                permutate_response_dir, prefix="permute_corexp_gene_id")

    permutate_seeds <- map(permute_coexp_list, str_extract, pattern = "seed[0-9]+")
    permute_coexp_tab_list <- map2(permute_coexp_list, permutate_seeds,
            function(x, y) get_lower_pair_coexp(x, y))

    permute_coexp_tab <- permute_coexp_tab_list |>
        purrr::reduce(merge, by = "gene_pairs")

    # test example
    # permute_coexp_tab_all <- permute_coexp_tab_list[1:3] |>
    #     purrr::reduce(merge, by = "gene_pairs", all = TRUE)

    # permutation test
    for( i in permute_coexp_tab[, .I]){
        set(permute_coexp_tab, i, "p_twotail",
                permutate_p_two(permute_coexp_tab[i, unlist(.SD),.SDcols=patterns("seed50")],
                    permute_coexp_tab[i, unlist(.SD),.SDcols=patterns("seed[1-4][0-9]*")]) )
        if (i %% 10000 == 0) cat("Processed", i, "genes\n")
    }

    sig_coexp_pairs <- permute_coexp_tab[p_twotail < 0.05, list(gene_pairs, p_twotail)]
    fwrite(tab_output_dir, sig_edge_tab_file)
}

permutation_test_stat_tab(coexp_expr_tab, permutation_test_stat_tab,
                        permutation_num = 50)