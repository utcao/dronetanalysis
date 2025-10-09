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
#' @return Subsetted data.table
subset_dataset <- function(data_table, n_genes = NULL, n_samples = NULL, seed = NULL) {
    
    if (!is.null(seed)) set.seed(seed)
    
    # Get available genes and samples
    genes <- data_table$gene_id
    samples <- colnames(data_table)[-1]  # exclude gene_id column
    
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
    
    # Randomly select genes and samples
    selected_genes <- sample(genes, n_genes)
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
#' @param subset_configs List of configurations (each with n_genes, n_samples, name)
#' @param seed Random seed for reproducibility
#' @return Named list of subsetted data.tables
create_multiple_subsets <- function(data_list, subset_configs, seed = 1234) {
    
    result_list <- list()
    
    for (config in subset_configs) {
        config_name <- config$name
        cat("Creating subset:", config_name, "\n")
        
        subsets <- purrr::map(data_list, function(data_table) {
            subset_dataset(
                data_table = data_table,
                n_genes = config$n_genes,
                n_samples = config$n_samples,
                seed = seed
            )
        })
        
        result_list[[config_name]] <- subsets
    }
    
    return(result_list)
}

#################################
##### WGCNA Functions ###########
#################################

# Hub gene identification
hubGenes <- function(datExpr, moduleColors, module) {
  # Select genes in the specified module
  moduleGenes <- (moduleColors == module)
  if (sum(moduleGenes) == 0) {
    return(NULL)
  }
  moduleExpr <- datExpr[, moduleGenes]
  # Calculate module eigengene
  ME <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes[, paste0("ME", module) ]
  # Calculate module membership (correlation with module eigengene)
  MM <- cor(moduleExpr, ME, use = "p")
  # Identify the gene with the highest module membership
  hubGene <- rownames(moduleExpr)[which.max(MM)]
  return(hubGene)
}