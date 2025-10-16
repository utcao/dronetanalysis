suppressPackageStartupMessages({
    library(purrr)
    library(data.table)
    library(glue)
    library(stringr)
    library(ggplot2)
    library(tibble)  # for rownames_to_column / column_to_rownames
    library(futile.logger)
})

#' Calculate Two-Sided Permutation Test P-value
#'
#' This function performs a two-sided permutation test by comparing the observed
#' test statistic to a null distribution generated through permutations.
#'
#' @param obs Numeric. The observed test statistic.
#' @param permu Numeric vector. The permutation statistics (null distribution).
#'
#' @return Numeric. The two-sided permutation p-value.
#'
#' @details
#' The two-sided p-value is calculated as:
#' p = (number of permutation statistics with absolute values >= absolute observed value) / total permutations
#'
#' This tests the null hypothesis that the observed effect is zero against the
#' alternative that it is different from zero (in either direction).
#'
#' @examples
#' \dontrun{
#' # Example: Correlation test
#' observed_cor <- 0.75
#' permuted_cors <- c(0.12, -0.08, 0.05, -0.15, 0.22, -0.07, 0.18, -0.03)
#' p_value <- permutate_p_two(observed_cor, permuted_cors)
#' }
#'
#' @export
permutate_p_two <- function(obs, permu){
    stopifnot({
        !is.numeric(obs)
        length(obs) != 1
        all(!is.numeric(permu))
        length(permu) > 0
    })
    (sum(abs(permu) >= abs(obs)))/(length(permu))
}

#' Permutate a Column in a Data table
#'
#' This function takes a data table and a column name, and returns a new data table
#' in which the values in the specified column are randomly shuffled while all other
#' columns remain unchanged. This is useful for creating permutation tests or randomization
#' experiments.
#'
#' @param dat_tab A data table The input table whose column will be permuted.
#' @param shuffle_column A string or symbol. The name of the column to shuffle.
#'
#' @return a name list of files. Meanwhile, it will write out data table with the same shape of the input
#'         but with the specified column permuted randomly.
#'
#' @examples
#' df <- data.table(a = 1:5, b = letters[1:5])
#' permuted_df <- permutate_col(df, "a")
permutate_col <- function(dat_tab, shuffle_col = "sample_name", permu_n = 30,
                        tab_output_dir, prefix = NULL, tag = NULL){
    prefix <- ifelse(is.null(prefix), "", glue("{prefix}_"))
    tag <- ifelse(is.null(tag), "", glue("_{tag}"))
    seeds <- seq(1, permu_n)
    map(seeds, function(seed){
        permutation_file <- file.path(tab_output_dir,
                        glue("{prefix}seed{seed}{tag}.csv"))
        set.seed(seed)
        dat_tab[ , (shuffle_col):=(sample(get(shuffle_col), nrow(dat_tab)))]
        fwrite(dat_tab, permutation_file)
        permutation_file
    })
}

#' Extract Upper Triangle Coexpression Pairs from Symmetric Matrix
#' 
#' This function reads a symmetric coexpression matrix, extracts the upper triangle 
#' (excluding diagonal), and returns a melted data table of gene pairs with their 
#' coexpression values.
#'
#' @param coexp_file Path to the coexpression file. Expected to be a symmetric matrix
#'   with gene IDs as row names and columns.
#' @param value_name_melt Character string specifying the name for the coexpression 
#'   value column in the output. Default is "coexpr".
#' @param pair_id_col Character string specifying the name for the gene pair ID column.
#'   Default is "gene_pairs".
#'
#' @return A data.table with two columns:
#'   \item{coexpression_value}{Column named according to \code{value_name_melt} 
#'         containing coexpression values}
#'   \item{gene_pair_ids}{Column named according to \code{pair_id_col} containing
#'         concatenated gene IDs in the format "gene1_gene2"}
#'
#' @details The function performs the following steps:
#'   \enumerate{
#'     \item Reads the symmetric matrix using \code{\link[data.table]{fread}}
#'     \item Converts to matrix format with gene IDs as row names
#'     \item Sets lower triangle and diagonal to NA
#'     \item Melts the matrix to long format keeping only non-NA values
#'     \item Creates unique gene pair identifiers
#'     \item Removes redundant columns
#'   }
#'   Note: The input matrix must be symmetric for meaningful results.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' coexp_pairs <- get_lower_pair_coexp("coexpression_matrix.csv")
#' 
#' # Custom column names
#' coexp_pairs <- get_lower_pair_coexp("coexpression_matrix.csv", 
#'                                    value_name_melt = "correlation",
#'                                    pair_id_col = "gene_pairs")
#' }
#'
#' @importFrom data.table fread as.data.table
#' @importFrom tibble column_to_rownames
#' @importFrom reshape2 melt
#' @export
get_lower_pair_coexp <- function(coexp_file, value_name_melt = "coexpr",
                                pair_id_col = "gene_pairs"){
    flog.info("The input matrix must be symmetric")
    coexp_reduncols <- c("gene_id", "gene_pair")
    if(is.character(coexp_file)){
        coexp_tab <- fread(coexp_file) |>
                tibble::column_to_rownames("gene_id") |>
                as.matrix()
    }else{
        coexp_tab <- coexp_file |>
                tibble::column_to_rownames("gene_id") |>
                as.matrix()
    }
    colnames(coexp_tab) <- row.names(coexp_tab)
    coexp_tab[lower.tri(coexp_tab, diag = TRUE)] <- NA
    coexp_upper_l <- coexp_tab |>
                as.data.table(keep.rownames = "gene_id") |>
                melt(id.vars = "gene_id", variable.name = "gene_pair",
                    value.name = value_name_melt, na.rm = TRUE)
    coexp_upper_l[, (pair_id_col):= paste0(gene_id, "_", gene_pair)]
    coexp_upper_l[, (coexp_reduncols):=NULL]
    coexp_upper_l
}

#' Convert long matrix of significant gene pairs to wide matrix
#' 
#' This function converts a long-format data table of gene pairs and their associated values
#' into a wide-format square matrix. The resulting matrix has genes as both row and column
#' names, with the values filled in from the input data table. The matrix is symmetric
#' with self-correlations (diagonal) set to 1.
#' @param long_dt A data.table in long format containing gene pairs and their values.
#' @param gene_pairs_col The name of the column containing gene pair identifiers.
#' @param value_col The name of the column containing the values to fill in the matrix.
#' @param separator The separator used in the gene pair identifiers (default is "_").
#' @return A wide-format square matrix with genes as both row and column names.
#' \item{gene1}{Row gene names}
#' \item{gene2}{Column gene names}
#' \item{value}{Values from the input data table}
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Splits the gene pair identifiers into two separate columns.
#'   \item Creates both directions for symmetry (gene1-gene2 and gene2-gene1).
#'   \item Adds diagonal entries (self-correlations = 1).
#'   \item Casts the combined data table to wide format using \code{\link[reshape2]{dcast}}.
#'   \item Converts the resulting data table to a matrix with proper row names.
#' }
#' @examples
#' \dontrun{
#' # Example usage
#' long_dt <- data.table(gene_pairs = c("geneA_geneB", "geneA_geneC", "geneB_geneC"),
#'                       rho = c(0.8, 0.5, 0.6))
#' wide_matrix <- long_to_square_dcast(long_dt, gene_pairs_col = "gene_pairs", value_col = "rho")
#' }
#' @importFrom data.table tstrsplit dcast
#' @export
long_to_square_dcast <- function(long_dt, gene_pairs_col = "gene_pairs", 
                                value_col = "rho", separator = "_") {
    
    # Split gene pairs and duplicate for symmetry
    long_dt[, c("gene1", "gene2") := tstrsplit(get(gene_pairs_col), separator, fixed = TRUE)]
    
    # Create both directions for symmetry
    dt1 <- long_dt[, .(gene1, gene2, value = get(value_col))]
    dt2 <- long_dt[, .(gene1 = gene2, gene2 = gene1, value = get(value_col))]
    combined_dt <- rbind(dt1, dt2)
    
    # Add diagonal (self-correlations = 1)
    all_genes <- unique(c(combined_dt$gene1, combined_dt$gene2))
    diagonal_dt <- data.table(gene1 = all_genes, gene2 = all_genes, value = 1)
    combined_dt <- rbind(combined_dt, diagonal_dt)
    
    # Cast to wide format
    square_matrix <- dcast(combined_dt, gene1 ~ gene2, value.var = "value", fill = NA)
    
    # Convert to matrix with proper row names
    gene_names <- square_matrix$gene1
    square_matrix[, gene1 := NULL]
    square_matrix <- as.matrix(square_matrix)
    rownames(square_matrix) <- gene_names
    
    return(square_matrix)
}