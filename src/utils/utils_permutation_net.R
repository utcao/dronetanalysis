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
        !is.numeric(obs) || length(obs) != 1
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
    coexp_tab <- fread(coexp_file) |>
                tibble::column_to_rownames("gene_id") |>
                as.matrix()

    coexp_tab[lower.tri(coexp_tab, diag = TRUE)] <- NA
    coexp_upper_l <- coexp_tab |>
                as.data.table(keep.rownames = "gene_id") |>
                melt(id.vars = "gene_id", variable.name = "gene_pair",
                    value.name = value_name_melt, na.rm = TRUE)
    coexp_upper_l[, (pair_id_col):= paste0(gene_id, "_", gene_pair)]
    coexp_upper_l[, (coexp_reduncols):=NULL]
    coexp_upper_l
}