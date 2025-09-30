#!/usr/bin/env Rscript
# ==============================================================================
# subset dataset for drosophila transcriptomic network analysis
#
# This script performs:
#   1. format strings in data table
#   2. subset dataset
# ==============================================================================
suppressPackageStartupMessages({
    library(purrr)
    library(data.table)
    library(glue)
    library(stringr)
    library(ggplot2)
    library(tibble)  # for rownames_to_column / column_to_rownames
})

source("src/utils/utils_io.R")
# -----  Extract Key Settings from config.yaml -----
library(yaml)  # to read the YAML config
config <- yaml::read_yaml("config/config.yaml")
# Output prefixes
rawdata_dir <- config$project_dirs$rawdata_dir
processed_data_dir <- config$project_dirs$processed_data_dir
test_data_dir <- config$project_dirs$test_data_dir

create_directories(c(rawdata_dir, processed_data_dir, test_data_dir))

# ----- subset dataset -----
ct_file <- file.path(rawdata_dir, "logCPM_Ctrl_Dros.csv")
hs_file <- file.path(rawdata_dir, "logCPM_HS_Dros.csv")
subset_ct_file <- file.path(processed_data_dir, "logCPM_ct_dros_subset.csv")
subset_hs_file <- file.path(processed_data_dir, "logCPM_hs_dros_subset.csv")
test_ct_file <- file.path(test_data_dir, "logCPM_ct_dros_subset.csv")
test_hs_file <- file.path(test_data_dir, "logCPM_hs_dros_subset.csv")
# data load
dro_tab_list <- purrr::map(c(ct_file, hs_file), function(f){
    tab <- fread(f)
    setnames(tab, "V1", "gene_id")
})
# random samples
set.seed(1234)
subset_size <- 100
dro_subtab_list <- purrr::map(dro_tab_list, function(tab){
    samples <- colnames(tab)[-1]
    sub_samples <- sample(samples, subset_size)
    tab[, .SD, .SDcols = c("gene_id", sub_samples)]
})
# write subsets for both conditions
purrr::walk2(dro_subtab_list,
            c(subset_ct_file, subset_hs_file),
            ~fwrite(.x, .y))

# create test dataset: 30 genes * 30 samples
test_size <- 30
dro_test_tab_list <- purrr::map(dro_tab_list, function(tab){
    samples <- colnames(tab)[-1]
    genes <- tab$gene_id
    sub_samples <- sample(samples, test_size)
    sub_genes <- sample(genes, test_size)
    tab[gene_id %in% sub_genes, .SD, .SDcols = c("gene_id", sub_samples)]
})

# write test dataset for both conditions
purrr::walk2(dro_test_tab_list,
            c(test_ct_file, test_hs_file),
            ~fwrite(.x, .y))
