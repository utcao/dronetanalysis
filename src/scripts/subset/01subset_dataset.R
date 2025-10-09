# !/usr/bin/env Rscript
# ==============================================================================
# Subset dataset for drosophila transcriptomic network analysis
#
# Script by Gabriel Thornes and Yu-tao Cao
#
# Last Updated: 08/10/2025
#
# This script::
#   1. Takes VOOM and VST processed matrices as input
#   2. Formats strings in data table
#   3. Subsets dataset for easier testing and debugging
#   4. Creates test datasets with smaller dimensions
# ==============================================================================


#################################
##### Packages and Setup ########
#################################

rm(list = ls())

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
subset_data_dir <- config$project_dirs$subset_data_dir
test_data_dir <- config$project_dirs$test_data_dir

# Create directories if they do not exist
create_directories(c(rawdata_dir, processed_data_dir, test_data_dir))

#####################
##### Datasets ######
#####################

# VST processed files
vst_files <- config$processed_files$vst_processed_files

# VOOM processed files
voom_files <- config$processed_files$voom_processed_files

# Output files for subsets
subset_vst_ctrl_file <- file.path(subset_data_dir, "VSTdataCtrl_subset.txt")
subset_vst_hs_file <- file.path(subset_data_dir, "VSTdataHS_subset.txt")
subset_voom_ctrl_file <- file.path(subset_data_dir, "voomdataCtrl_subset.txt")
subset_voom_hs_file <- file.path(subset_data_dir, "voomdataHS_subset.txt")

# Test files (smaller subsets)
test_vst_ctrl_file <- file.path(test_data_dir, "VSTdataCtrl_test.txt")
test_vst_hs_file <- file.path(test_data_dir, "VSTdataHS_test.txt")
test_voom_ctrl_file <- file.path(test_data_dir, "voomdataCtrl_test.txt")
test_voom_hs_file <- file.path(test_data_dir, "voomdataHS_test.txt")

# data load - VST data
vst_data_list <- purrr::map(c(vst_files), function(f){
    cat("Loading:", f, "\n")
    tab <- fread(f, header = TRUE)
    # First column gene names/IDs
    if(!"gene_id" %in% colnames(tab)) {
        setnames(tab, 1, "gene_id")  # rename first column to gene_id
    }
    return(tab)
})

# data load - VOOM data
voom_data_list <- purrr::map(c(voom_files), function(f){
    cat("Loading:", f, "\n")
    tab <- fread(f, header = TRUE)
    # First column gene names/IDs
    if(!"gene_id" %in% colnames(tab)) {
        setnames(tab, 1, "gene_id")  # rename first column to gene_id
    }
    return(tab)
})


###########################
##### Create subsets ######
###########################


# random samples
set.seed(1234)
subset_size <- 100

# Create subsets for VST data
vst_subset_list <- purrr::map(vst_data_list, function(tab){
    samples <- colnames(tab)[-1]  # exclude gene_id column
    sub_samples <- sample(samples, min(subset_size, length(samples)))
    tab[, .SD, .SDcols = c("gene_id", sub_samples)]
})

# Create subsets for VOOM data
voom_subset_list <- purrr::map(voom_data_list, function(tab){
    samples <- colnames(tab)[-1]  # exclude gene_id column
    sub_samples <- sample(samples, min(subset_size, length(samples)))
    tab[, .SD, .SDcols = c("gene_id", sub_samples)]
})

# write subsets for VST data
purrr::walk2(vst_subset_list,
            c(subset_vst_ctrl_file, subset_vst_hs_file),
            ~fwrite(.x, .y, sep = "\t"))

# write subsets for VOOM data
purrr::walk2(voom_subset_list,
            c(subset_voom_ctrl_file, subset_voom_hs_file),
            ~fwrite(.x, .y, sep = "\t"))


#################################
##### Create test datasets ######
#################################

# set test matrix size = 30*30
test_size <- 30

# Create test datasets for VST data
vst_test_list <- purrr::map(vst_data_list, function(tab){
    samples <- colnames(tab)[-1]
    genes <- tab$gene_id
    sub_samples <- sample(samples, min(test_size, length(samples)))
    sub_genes <- sample(genes, min(test_size, length(genes)))
    tab[gene_id %in% sub_genes, .SD, .SDcols = c("gene_id", sub_samples)]
})

# Create test datasets for VOOM data
voom_test_list <- purrr::map(voom_data_list, function(tab){
    samples <- colnames(tab)[-1]
    genes <- tab$gene_id
    sub_samples <- sample(samples, min(test_size, length(samples)))
    sub_genes <- sample(genes, min(test_size, length(genes)))
    tab[gene_id %in% sub_genes, .SD, .SDcols = c("gene_id", sub_samples)]
})

# write test datasets for VST data
purrr::walk2(vst_test_list,
            c(test_vst_ctrl_file, test_vst_hs_file),
            ~fwrite(.x, .y, sep = "\t"))

# write test datasets for VOOM data
purrr::walk2(voom_test_list,
            c(test_voom_ctrl_file, test_voom_hs_file),
            ~fwrite(.x, .y, sep = "\t"))

cat("Subset creation complete\n")