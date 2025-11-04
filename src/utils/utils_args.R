# utils_args.R - Shared argument parsing utilities
suppressPackageStartupMessages(library(argparse))

#' Create base parser with common arguments
#' @param description Script description
#' @return ArgumentParser object with common args added
create_base_parser <- function(description = "Bioinformatics analysis script") {
    parser <- ArgumentParser(description = description)
    
    # Common arguments across all scripts
    parser$add_argument('--config', type="character", 
                       default="config/config.yaml",
                       help='Configuration file path [default: config/config.yaml]')
    parser$add_argument('--output-dir', type="character", 
                       default="results",
                       help='Base output directory [default: results]')
    parser$add_argument('--verbose', action="store_true",
                       help='Enable verbose output')
    parser$add_argument('--dry-run', action="store_true",
                       help='Show what would be done without executing')
    
    return(parser)
}

#' Add HPC-specific arguments
#' @param parser ArgumentParser object
#' @return Modified parser with HPC args
add_hpc_args <- function(parser) {
    parser$add_argument('--task-id', type="integer",
                       default=as.integer(Sys.getenv("SGE_TASK_ID", 
                                                   Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))),
                       help='Task ID for array job [default: 1 or from SGE_TASK_ID/SLURM_ARRAY_TASK_ID]')
    parser$add_argument('--total-tasks', type="integer", default=1,
                       help='Total number of array tasks [default: 1]')
    return(parser)
}

#' Add permutation-specific arguments
#' @param parser ArgumentParser object
#' @return Modified parser with permutation args
add_permutation_args <- function(parser) {
    parser$add_argument('--permutations', type="integer", default=1000,
                       help='Number of permutations [default: 1000]')
    parser$add_argument('--alpha', type="double", default=0.05,
                       help='Significance threshold [default: 0.05]')
    parser$add_argument('--seed', type="integer", default=42,
                       help='Random seed for reproducibility [default: 42]')
    parser$add_argument('--keep-permutations', action="store_true",
                       help='Keep all permutation columns in output')
    return(parser)
}

#' Add network analysis arguments
#' @param parser ArgumentParser object
#' @return Modified parser with network args
add_network_args <- function(parser) {
    parser$add_argument('--min-module-size', type="integer", default=30,
                       help='Minimum module size for WGCNA [default: 30]')
    parser$add_argument('--merge-threshold', type="double", default=0.25,
                       help='Module merge threshold [default: 0.25]')
    parser$add_argument('--adjacency-threshold', type="double", default=0.1,
                       help='Adjacency threshold for network construction [default: 0.1]')
    parser$add_argument('--deep-split', type="integer", default=2,
                       help='Deep split parameter for module detection [default: 2]')
    return(parser)
}

#' Add file input/output arguments
#' @param parser ArgumentParser object
#' @return Modified parser with I/O args
add_io_args <- function(parser) {
    parser$add_argument('--input', type="character", required=TRUE,
                       help='Input file path [REQUIRED]')
    parser$add_argument('--correlation-matrix', type="character",
                       help='Correlation matrix file')
    parser$add_argument('--adjacency-matrix', type="character", 
                       help='Adjacency matrix file')
    parser$add_argument('--gene-subset', type="integer",
                       help='Number of genes to subset for testing (e.g., 1000)')
    return(parser)
}

#' Add matrix conversion arguments
#' @param parser ArgumentParser object
#' @return Modified parser with matrix conversion args
add_matrix_args <- function(parser) {
    parser$add_argument('--gene-pairs-col', type="character", default="gene_pairs",
                       help='Column name for gene pairs [default: gene_pairs]')
    parser$add_argument('--value-col', type="character", default="rho",
                       help='Column name for correlation values [default: rho]')
    parser$add_argument('--separator', type="character", default="_",
                       help='Separator for gene pairs [default: _]')
    return(parser)
}
