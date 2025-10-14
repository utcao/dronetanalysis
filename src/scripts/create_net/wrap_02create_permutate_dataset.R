# or consider split shuffle and permutation test into separate scripts

#' Perform Permutation Test for Coexpression Network Edge Significance
#'
#' This function conducts a comprehensive permutation test to identify statistically 
#' significant edges in a coexpression network. It compares observed coexpression 
#' values against a null distribution generated through column permutation, 
#' performs hypothesis testing, and returns significant edges after multiple 
#' testing correction.
#'
#' @param coexp_expr_tab A data.table containing gene expression data with genes 
#'   as rows and samples as columns. Must contain the column specified in 
#'   `permutate_col`.
#' @param sig_edge_tab_file Character string specifying the output filename for 
#'   significant edges results (without path).
#' @param tab_output_dir Character string specifying the base directory for 
#'   output files.
#' @param obs_col Character string or regular expression pattern to 
#'   identify the column(s) containing observed coexpression values.
#' @param permutate_col Character string specifying the column name to permute. 
#'   Default is "gene_id".
#' @param permut_cols_pattern Character string or regular expression pattern to 
#'   identify columns containing permutation results. Default is "seed".
#' @param permutation_num Integer specifying the number of permutations to 
#'   perform. Default is 50.
#' @param alpha Numeric significance level for hypothesis testing. Default is 0.05.
#'
#' @return A data.table containing significant coexpression pairs with the 
#'   following columns:
#'   \item{gene_pairs}{Character identifier for each gene pair}
#'   \item{obs_col}{One or more columns with observed coexpression values}
#'   \item{p_twotail}{Two-sided permutation p-values}
#'   \item{q_twotail_fdr}{FDR-adjusted q-values}
#'   Also writes the results to a file in the permutation_test subdirectory.
#'
#' @details
#' The function performs the following workflow:
#' \enumerate{
#'   \item \strong{Validation}: Checks that input is a data.table and contains 
#'         the specified permutation column
#'   \item \strong{Directory Setup}: Creates "shuffle" and "permutation_test" 
#'         subdirectories for intermediate and final results
#'   \item \strong{Permutation}: Generates multiple permuted datasets by shuffling 
#'         the specified column using `permutate_col()`
#'   \item \strong{Coexpression Calculation}: Computes coexpression for each 
#'         permuted dataset using `get_lower_pair_coexp()`
#'   \item \strong{Merging}: Combines all permutation results with observed data
#'   \item \strong{Hypothesis Testing}: Performs two-sided permutation tests for 
#'         each gene pair using `permutate_p_two()`
#'   \item \strong{Multiple Testing Correction}: Applies FDR correction using 
#'         Benjamini-Hochberg method
#'   \item \strong{Output}: Returns significant results and writes to file
#' }
#'
#' The permutation test creates a null distribution by breaking the relationship 
#' between genes through column shuffling, then compares observed coexpression 
#' values against this distribution to assess statistical significance.
#'
#' @note
#' \itemize{
#'   \item Requires helper functions: `create_directories()`, `permutate_col()`, 
#'         `get_lower_pair_coexp()`, and `permutate_p_two()`
#'   \item Progress is reported every 10,000 gene pairs processed
#'   \item Intermediate permutation files are saved in "shuffle" subdirectory
#'   \item Only columns matching `obs_col` are kept in final output
#'   \item Uses strict inequality (p_twotail <= alpha) for significance testing
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' sig_edges <- permutation_test_stat_tab(
#'   coexp_expr_tab = expression_data,
#'   sig_edge_tab_file = "significant_edges.csv",
#'   tab_output_dir = "/path/to/results",
#'   obs_col = "seed20"
#' )
#'
#' # Customized permutation test
#' sig_edges <- permutation_test_stat_tab(
#'   coexp_expr_tab = expression_data,
#'   sig_edge_tab_file = "network_edges_fdr05.csv",
#'   tab_output_dir = "/analysis/coexpression",
#'   obs_col = "observed_cor",
#'   permutate_col = "gene_symbol",
#'   permut_cols_pattern = "perm_[0-9]+",
#'   permutation_num = 100,
#'   alpha = 0.01
#' )
#' }
#'
#' @seealso
#' \code{\link{get_lower_pair_coexp}} for coexpression matrix calculation
#' \code{\link{permutate_p_two}} for permutation p-value computation
#' \code{\link[stats]{p.adjust}} for FDR correction
#'
#' @importFrom data.table .SD .I :=
#' @importFrom purrr map map2 reduce
#' @importFrom stringr str_extract
#' @importFrom stats p.adjust
#' @importFrom futile.logger flog.info
#' @export
permutation_test_stat_tab <- function(coexp_expr_tab,
                                    sig_edge_tab_file, tab_output_dir,
                                    obs_col = "rho",
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

    coexp_pairs <- get_lower_pair_coexp(coexp_expr_tab, obs_col)

    permute_coexp_tab <- c(list(coexp_pairs), permute_coexp_tab_list) |>
        purrr::reduce(merge, by = "gene_pairs")

    # permutation test between obeservation and permutation
    for( i in permute_coexp_tab[, .I]){
        set(permute_coexp_tab, i, "p_twotail",
                permutate_p_two(
                    permute_coexp_tab[i, unlist(.SD),
                                        .SDcols=patterns(obs_col)],
                    permute_coexp_tab[i, unlist(.SD),
                                        .SDcols=patterns(permut_cols_pattern)])
                )
        if (i %% 10000 == 0) cat("Processed", i, "gene pairs\n")
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