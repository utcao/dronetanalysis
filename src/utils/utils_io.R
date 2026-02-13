## utils_io.R
library(data.table)
library(purrr)
library(glue)

#################################
##### IO Utility Functions ######
#################################

#' Create directories if they do not exist
#' 
#' @param dir_list A character vector (or list) of directory paths
#' @return None (side effect: creates directories on disk)
create_directories <- function(dir_list) {
  purrr::walk(dir_list, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))
}

#################################
##### Subsetting Functions ######
#################################

#' Subset dataset with specified gene and sample numbers
#'
#' @param data_table A data.table with gene_id as first column
#' @param n_genes Number of genes to select (default: all genes)
#' @param n_samples Number of samples to select (default: all samples)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param required_gene_ids Character vector of gene IDs that must be included
#'   in the subset (default: NULL). Genes not found in the dataset are reported.
#'   These genes are included first, then random selection fills the remainder.
#' @return Subsetted data.table
subset_dataset <- function(data_table, n_genes = NULL, n_samples = NULL,
                           seed = NULL, required_gene_ids = NULL) {

    if (!is.null(seed)) set.seed(seed)

    # Get available genes and samples
    genes <- data_table$gene_id
    samples <- colnames(data_table)[-1]  # exclude gene_id column

    # Handle required gene IDs
    required_found <- character(0)
    if (!is.null(required_gene_ids)) {
        missing <- setdiff(required_gene_ids, genes)
        if (length(missing) > 0) {
            cat("  WARNING: required gene IDs not found in dataset (",
                length(missing), "/", length(required_gene_ids), "):\n")
            cat("    ", paste(head(missing, 20), collapse = ", "))
            if (length(missing) > 20) cat(" ... and", length(missing) - 20, "more")
            cat("\n")
        }
        required_found <- intersect(required_gene_ids, genes)
        if (length(required_found) > 0) {
            cat("  Including", length(required_found), "required gene(s)\n")
        }
    }

    # Determine number of genes to select
    if (is.null(n_genes)) {
        n_genes <- length(genes)
    } else {
        n_genes <- min(n_genes, length(genes))
    }

    # Determine number of samples to select
    if (is.null(n_samples)) {
        n_samples <- length(samples)
    } else {
        n_samples <- min(n_samples, length(samples))
    }

    # Select genes: required first, then random from the rest
    if (length(required_found) > 0) {
        n_required <- length(required_found)
        if (n_required >= n_genes) {
            # More required genes than requested total — use all required
            selected_genes <- required_found
            cat("  Note: required genes (", n_required,
                ") >= requested n_genes (", n_genes,
                "), using all required genes\n")
        } else {
            # Fill remaining slots with random selection from non-required genes
            remaining_pool <- setdiff(genes, required_found)
            n_random <- n_genes - n_required
            random_genes <- sample(remaining_pool, min(n_random, length(remaining_pool)))
            selected_genes <- c(required_found, random_genes)
        }
    } else {
        selected_genes <- sample(genes, n_genes)
    }

    selected_samples <- sample(samples, n_samples)

    # Create subset
    subset_data <- data_table[gene_id %in% selected_genes,
                             .SD, .SDcols = c("gene_id", selected_samples)]

    cat("Created subset with", nrow(subset_data), "genes and",
        ncol(subset_data) - 1, "samples\n")

    return(subset_data)
}

#' Create multiple subsets with different dimensions
#'
#' @param data_list List of data.tables to subset
#' @param subset_configs List of configurations (each with n_genes, n_samples,
#'   name, and optionally required_gene_ids)
#' @param seed Random seed for reproducibility
#' @param required_gene_ids Character vector of gene IDs that must be included
#'   in ALL subsets (default: NULL). Can also be set per-config via
#'   config$required_gene_ids, which takes precedence over this global default.
#' @return Named list of subsetted data.tables
create_multiple_subsets <- function(data_list, subset_configs, seed = 1234,
                                   required_gene_ids = NULL) {

    result_list <- list()

    for (config in subset_configs) {
        config_name <- config$name
        cat("Creating subset:", config_name, "\n")

        # Per-config required_gene_ids takes precedence over global
        config_required <- config$required_gene_ids
        if (is.null(config_required)) {
            config_required <- required_gene_ids
        }

        subsets <- purrr::map(data_list, function(data_table) {
            subset_dataset(
                data_table = data_table,
                n_genes = config$n_genes,
                n_samples = config$n_samples,
                seed = seed,
                required_gene_ids = config_required
            )
        })

        result_list[[config_name]] <- subsets
    }

    return(result_list)
}