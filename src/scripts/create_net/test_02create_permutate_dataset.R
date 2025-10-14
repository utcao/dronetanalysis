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

expcor_tab_file <- file.path(expcor_tab_dir, "spearman_correlation_matrix.csv")

# ----- 3. Load data -----
coexp_expr_tab <- fread(expcor_tab_file)[1:3000, 1:3001]

# ----- 4. create permutation datasets by shuffle gene_id -----
# or consider split shuffle and permutation test into separate scripts
permutation_test_stat_tab <- function(coexp_expr_tab,
                                    sig_edge_tab_file, tab_output_dir,
                                    obs_col_pattern,
                                    permutate_col = "gene_id",
                                    permut_cols_pattern = "seed",
                                    permutation_num = 50, alpha = 0.05){
    stopifnot({
        is.data.table(coexp_expr_tab)
        permutate_col %in% colnames(coexp_expr_tab)
    })
    permutate_dir <- file.path(tab_output_dir, "shuffle")
    permutate_res_tab_dir <- file.path(tab_output_dir, "permutation_test")

    create_directories(c(permutate_dir, permutate_res_tab_dir))
    
    permute_coexp_list <- permutate_col(coexp_expr_tab, permutate_col, permutation_num,
                permutate_dir, prefix="permute_corexp_gene_id")

    permutate_seeds <- map(permute_coexp_list, str_extract, pattern = "seed[0-9]+")
    permute_coexp_tab_list <- map2(permute_coexp_list, permutate_seeds,
            function(x, y) get_lower_pair_coexp(x, y))

    permute_coexp_tab <- permute_coexp_tab_list |>
        purrr::reduce(merge, by = "gene_pairs")

    # permutation test between obeservation and permutation
    for( i in permute_coexp_tab[, .I]){
        set(permute_coexp_tab, i, "p_twotail",
                permutate_p_two(
                    permute_coexp_tab[i, unlist(.SD),
                                        .SDcols=patterns(obs_col_pattern)],
                    permute_coexp_tab[i, unlist(.SD),
                                        .SDcols=patterns(permut_cols_pattern)])
                )
        if (i %% 10000 == 0) cat("Processed", i, "genes\n")
    }
    # keeps observation columns and result
    permut_cols <- permute_coexp_tab[, colnames(.SD),
                                        .SDcols=patterns(permut_cols_pattern)]
    flog.info("The permutation comes from:")
    print(permut_cols)
    
    keep_cols <- setdiff(colnames(permute_coexp_tab), permut_cols)
    flog.info("The obersavation value will be kept in the output:")
    print(keep_cols)

    sig_coexp_pairs <- permute_coexp_tab[p_twotail <= alpha, .SD, .SDcols = keep_cols]
    sig_coexp_pairs[, q_twotail_fdr:=p.adjust(p_twotail, method = "fdr")]
    fwrite(sig_coexp_pairs, file.path(permutate_res_tab_dir, sig_edge_tab_file))
    return(sig_coexp_pairs[])
}

sig_edge_tab_file <- "sig_edges_coexpr_net.csv"

sig_coexp_pairs <- permutation_test_stat_tab(coexp_expr_tab[1:30, 1:31],
                        sig_edge_tab_file,
                        expcor_tab_dir,
                        obs_col_pattern = "seed1",
                        permut_cols_pattern = "(seed1[0-9])|(seed[2-9])",
                        permutation_num = 20,
                        alpha = 1)