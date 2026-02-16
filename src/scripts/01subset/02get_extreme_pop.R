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
# ==============================================================================

rm(list = ls())

#######################################################
#  ||  Enter subset and test dataset sizes below  ||  #
#  \/    enter NULL to select all samples/genes   \/  #
#######################################################

subset_sample_size <- NULL
subset_gene_size   <- NULL
test_sample_size   <- 939
test_gene_size     <- 500

# Gene IDs that MUST be included in subsets (NULL = no requirement).
# These genes are included first, then random selection fills the rest.
# Genes not found in the dataset are reported as warnings.
# Can be a character vector, e.g.: c("FBgn0000014", "FBgn0000015")
# or read from a file, e.g.: readLines("config/required_genes.txt")
# Dicer 1 & 2:  FBgn0039016, FBgn0034246
# Chaperon
# Hsp83, Trap1, Gp93: FBgn0026761, FBgn0001233, FBgn0039562
# DnaJ-1: FBgn0263106 
# Kinase_ids <- readLines("dataset/flybase/FlyBase_IDs_KIN_KINASES.txt")
hsp_ids <- readLines("dataset/flybase/FlyBase_IDs_HSP_HEAT_SHOCK_PROTEINS.txt")
required_gene_ids  <- c("FBgn0039016", "FBgn0034246", hsp_ids)
#################################
##### Packages and Setup ########
#################################


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
vst_files <- config$processed_files$vst_files

# VOOM processed files
voom_files <- config$processed_files$voom_files

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
    
    # First, read the header line to get sample IDs
    header_line <- readLines(f, n = 1)
    sample_ids <- unlist(strsplit(header_line, "\t"))
    
    # Now read the data without treating first line as header
    tab <- fread(f, header = FALSE, skip = 1)
    
    # Set proper column names: gene_id + sample_ids
    setnames(tab, c("gene_id", sample_ids))
    
    cat("Loaded data with", nrow(tab), "genes and", ncol(tab)-1, "samples\n")
    cat("Sample IDs:", paste(head(sample_ids, 5), collapse = ", "), "...\n")
    
    return(tab)
})

# data load - VOOM data
voom_data_list <- purrr::map(c(voom_files), function(f){
    cat("Loading:", f, "\n")
    
    # First, read the header line to get sample IDs
    header_line <- readLines(f, n = 1)
    sample_ids <- unlist(strsplit(header_line, "\t"))
    
    # Now read the data without treating first line as header
    tab <- fread(f, header = FALSE, skip = 1)
    
    # Set proper column names: gene_id + sample_ids
    setnames(tab, c("gene_id", sample_ids))
    
    cat("Loaded data with", nrow(tab), "genes and", ncol(tab)-1, "samples\n")
    cat("Sample IDs:", paste(head(sample_ids, 5), collapse = ", "), "...\n")
    
    return(tab)
})


###########################
##### Create subsets ######
###########################

# Define subset configurations
subset_configs <- list(
    list(name = "subset", n_genes = subset_gene_size, n_samples = subset_sample_size),
    list(name = "test", n_genes = test_gene_size, n_samples = test_sample_size)
)

# Create subsets for VST data
cat("Creating VST subsets...\n")
vst_subsets <- create_multiple_subsets(vst_data_list, subset_configs, seed = 1234,
                                      required_gene_ids = required_gene_ids)

# Create subsets for VOOM data
cat("Creating VOOM subsets...\n")
voom_subsets <- create_multiple_subsets(voom_data_list, subset_configs, seed = 1234,
                                       required_gene_ids = required_gene_ids)

# Write VST subsets
cat("Writing VST subsets...\n")
fwrite(vst_subsets$subset[[1]], subset_vst_ctrl_file, sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(vst_subsets$subset[[2]], subset_vst_hs_file, sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(vst_subsets$test[[1]], test_vst_ctrl_file, sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(vst_subsets$test[[2]], test_vst_hs_file, sep = "\t", col.names = TRUE, row.names = FALSE)

# Write VOOM subsets
cat("Writing VOOM subsets...\n")
fwrite(voom_subsets$subset[[1]], subset_voom_ctrl_file, sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(voom_subsets$subset[[2]], subset_voom_hs_file, sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(voom_subsets$test[[1]], test_voom_ctrl_file, sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(voom_subsets$test[[2]], test_voom_hs_file, sep = "\t", col.names = TRUE, row.names = FALSE)

cat("Subset creation complete\n")