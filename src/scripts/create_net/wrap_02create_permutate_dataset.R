# or consider split shuffle and permutation test into separate scripts

#' Perform Permutation Test for Co-expression Network Analysis
#'
#' This function conducts a permutation-based significance test for co-expression networks
#' by shuffling gene labels and comparing observed co-expression values against a null
#' distribution generated from multiple permutations.
#'
#' @param coexp_expr_tab A data.table containing co-expression values with gene pairs
#'   and their correlation coefficients. Must contain columns specified by `id_col` and
#'   `pair_id_col`.
#' @param sig_edge_tab_file Character string specifying the output filename for
#'   significant edges after permutation testing.
#' @param tab_output_dir Character string specifying the directory where output files
#'   will be saved.
#' @param obs_col Character string specifying the column name containing observed
#'   co-expression values (e.g., correlation coefficients). Default is "rho".
#' @param id_col Character string specifying the column name containing gene identifiers
#'   to be shuffled during permutation. Default is "gene_id".
#' @param pair_id_col Character string specifying the column name containing unique
#'   gene pair identifiers. Default is "gene_pairs".
#' @param permu_cols_pattern Character string pattern used to identify permutation
#'   columns in the results. Default is "seed".
#' @param permu_n Integer specifying the number of permutations to perform. Default is 50.
#' @param alpha Numeric value specifying the significance threshold for p-values.
#'   Default is 0.05.
#' @param keep_permu Logical indicating whether to keep all permutation columns in the
#'   output. If FALSE, only observed values and test results are kept. Default is FALSE.
#'
#' @return A data.table containing significant co-expression edges that passed the
#'   permutation test, with the following columns:
#'   \itemize{
#'     \item \code{gene_pairs}: Unique identifier for each gene pair
#'     \item \code{rho}: Observed co-expression value (or specified \code{obs_col})
#'     \item \code{p_twotail}: Two-tailed empirical p-value from permutation test
#'     \item \code{q_twotail_fdr}: FDR-adjusted q-values for multiple testing correction
#'   }
#'   If \code{keep_permu = TRUE}, additional columns with permutation values are included.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input data structure and creates output directories
#'   \item Generates permuted datasets by shuffling gene labels
#'   \item Computes co-expression pairs for both observed and permuted data
#'   \item Performs two-tailed permutation tests for each gene pair
#'   \item Applies FDR correction for multiple testing
#'   \item Writes significant edges to file and returns results
#' }
#'
#' The permutation test assesses whether observed co-expression values are significantly
#' different from what would be expected by chance, providing a robust non-parametric
#' approach for identifying biologically relevant co-expression relationships.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' sig_edges <- permutation_test_stat(
#'   coexp_expr_tab = my_coexpression_data,
#'   sig_edge_tab_file = "significant_edges.csv",
#'   tab_output_dir = "results/coexpression",
#'   obs_col = "correlation",
#'   permu_n = 100,
#'   alpha = 0.01
#' )
#' }
#'
#' @seealso
#' \code{\link{permutate_col}}, \code{\link{get_lower_pair_coexp}}, 
#' \code{\link{permutate_p_two}}
#'
#' @export
permutation_test_stat <- function(coexp_expr_tab,
                                    sig_edge_tab_file, tab_output_dir,
                                    obs_col = "rho",
                                    id_col = "gene_id",
                                    pair_id_col = "gene_pairs",
                                    permu_cols_pattern = "seed",
                                    permu_n = 50, alpha = 0.05,
                                    keep_permu = FALSE){
    stopifnot({
        is.data.table(coexp_expr_tab)
        id_col %in% colnames(coexp_expr_tab)
    })
    col_feats <- setdiff(colnames(coexp_expr_tab), id_col)
    row_feats <- coexp_expr_tab[, get(id_col)]
    stopifnot({
        all(col_feats == row_feats)
    })

    permutate_dir <- file.path(tab_output_dir, "shuffle")
    permutate_res_tab_dir <- file.path(tab_output_dir, "permutation_test")

    create_directories(c(permutate_dir, permutate_res_tab_dir))
    
    permute_coexp_list <- permutate_col(
                            dat_tab = coexp_expr_tab,
                            shuffle_col = id_col,
                            permu_n = permu_n,
                            tab_output_dir = permutate_dir,
                            prefix="permute_corexp_gene_id")

    permutate_seeds <- map(permute_coexp_list, str_extract, pattern = "seed[0-9]+")
    permute_coexp_tab_list <- map2(permute_coexp_list, permutate_seeds,
                                 ~ get_lower_pair_coexp(
                                    coexp_file = .x, value_name_melt = .y,
                                    id_col = id_col, pair_id_col = pair_id_col))

    coexp_pairs <- get_lower_pair_coexp(coexp_expr_tab, value_name_melt = obs_col,
                                    id_col = id_col, pair_id_col = pair_id_col)

    permute_coexp_tab <- c(list(coexp_pairs), permute_coexp_tab_list) |>
        purrr::reduce(merge, by = "gene_pairs")

    # permutation test between obeservation and permutation
    for( i in permute_coexp_tab[, .I]){
        set(permute_coexp_tab, i, "p_twotail",
                permutate_p_two(
                    permute_coexp_tab[i, unlist(.SD),
                                        .SDcols=patterns(obs_col)],
                    permute_coexp_tab[i, unlist(.SD),
                                        .SDcols=patterns(permu_cols_pattern)])
                )
        if (i %% 10000 == 0) cat("Processed", i, "gene pairs\n")
    }

    # keeps observation columns and result
    permut_cols <- permute_coexp_tab[, colnames(.SD),
                                        .SDcols=patterns(permu_cols_pattern)]
    flog.info(glue("The permutation comes from {length(permut_cols)} datasets:"))
    print(permut_cols)
    if(keep_permu){
        keep_cols <- colnames(permute_coexp_tab)
    }else{
        keep_cols <- setdiff(colnames(permute_coexp_tab), permut_cols)
        flog.info("The columns kept in the output are:")
        print(keep_cols)
    }

    sig_coexp_pairs <- permute_coexp_tab[p_twotail <= alpha, .SD, .SDcols = keep_cols]
    sig_coexp_pairs[, q_twotail_fdr:=p.adjust(p_twotail, method = "fdr")]
    fwrite(sig_coexp_pairs, file.path(permutate_res_tab_dir, sig_edge_tab_file))

    raw_edges_n <- nrow(permute_coexp_tab)
    valid_edges_n <- nrow(sig_coexp_pairs)
    vaild_perc <- signif(valid_edges_n/raw_edges_n, 3)*100
    flog.info(glue(" {vaild_perc}% ( {valid_edges_n} out of {raw_edges_n}) edges survived from {permu_n} times of permutation test"))

    return(sig_coexp_pairs[])
}

#' Plot Distribution of Co-expression Values for Selected Gene Pairs
#'
#' This function visualizes the distribution of co-expression values (e.g., correlation coefficients)
#' for specified gene pairs, highlighting the observed value against the background distribution
#' from permutations or multiple datasets.
#'
#' @param coexp_pairs_ltab A data.table containing co-expression values for multiple
#'   gene pairs across different datasets or permutations. Must contain columns specified
#'   by `id_col`, `val_col`, and `var_col`.
#' @param id_col Character string specifying the column name containing gene pair
#'   identifiers. Default is typically "gene_pairs".
#' @param choose_paris_n Integer specifying the number of random gene pairs to select
#'   when `chose_pair` is NULL. Default is 1.
#' @param chose_pair Character string or vector specifying specific gene pair(s) to plot.
#'   If NULL, random pairs are selected based on `choose_paris_n`. Default is NULL.
#' @param val_col Character string specifying the column name containing co-expression
#'   values to plot (e.g., correlation coefficients). Default is "rho".
#' @param var_col Character string specifying the column name that distinguishes between
#'   different datasets or permutation types. Used to identify the observed value.
#'   Default is "datasets".
#' @param plot_output_dir Character string specifying the directory where the plot
#'   will be saved. If NULL, plot is not saved to file. Default is NULL.
#' @param prefix Character string specifying a prefix for the output filename.
#'   Default is NULL.
#'
#' @return A ggplot object displaying:
#'   \itemize{
#'     \item Density distribution of co-expression values for the selected gene pair
#'     \item Vertical dashed red line indicating the observed value
#'     \item Text annotation labeling the observed value
#'   }
#'   The plot can be further customized using standard ggplot2 functions.
#'
#' @details
#' This function is particularly useful for:
#' \itemize{
#'   \item Visualizing the distribution of permuted co-expression values for specific gene pairs
#'   \item Comparing observed co-expression values against null distributions
#'   \item Quality control of permutation tests by examining individual pairs
#'   \identifying outliers or unusual patterns in co-expression relationships
#' }
#'
#' When `chose_pair` is NULL, the function randomly selects gene pairs using a fixed
#' seed (set.seed(1)) for reproducibility. The observed value is identified as the
#' record where `var_col` equals `val_col`.
#'
#' @examples
#' \dontrun{
#' # Plot a random gene pair
#' p1 <- plot_permu_distr(
#'   coexp_pairs_ltab = permutation_results,
#'   id_col = "gene_pairs",
#'   val_col = "correlation",
#'   var_col = "dataset_type"
#' )
#'
#' # Plot specific gene pairs
#' p2 <- plot_permu_distr(
#'   coexp_pairs_ltab = permutation_results,
#'   id_col = "gene_pairs", 
#'   chose_pair = c("geneA_geneB", "geneC_geneD"),
#'   val_col = "rho",
#'   plot_output_dir = "results/plots",
#'   prefix = "validation_"
#' )
#'
#' # Customize the returned ggplot
#' p1 + theme_minimal() + labs(title = "Custom Title")
#' }
#'
#' @seealso
#' \code{\link{permutation_test_stat}}, \code{\link[ggplot2]{geom_density}},
#' \code{\link[ggplot2]{geom_vline}}
#'
#' @export
plot_permu_distr <- function(coexp_pairs_ltab, id_col,
                            choose_paris_n = 1,
                            chose_pair = NULL,
                            val_col = "rho",
                            var_col = "datasets",
                            plot_output_dir = NULL, prefix = NULL){
    if(is.null(chose_pair)){
        set.seed(1)
        coexpr_pairs <- coexp_pairs_ltab[, unique(get(id_col))]
        chose_pair <- sample(coexpr_pairs, choose_paris_n)
    }

    coexpr_pair_ltab <- coexp_pairs_ltab[get(id_col) == chose_pair]
    highlight_value <- coexp_pairs_ltab[get(id_col) == chose_pair &
                                            get(var_col) == val_col, get(val_col)]

    p <- ggplot(coexpr_pair_ltab, aes(x = get(val_col))) +
            geom_density(fill = "skyblue", alpha = 0.5) +
            # geom_rug(alpha = 0.2) +
            labs(title = glue("{val_col} distribution of {chose_pair}"),
                x = val_col, y = "") +
                geom_vline(xintercept = highlight_value, color = "red",
                            linetype = "dashed", linewidth = 1) +
                annotate("text", x = highlight_value, y = 0.025,
                                label = val_col, color = "red",
                                angle = 90, vjust = -0.5,
                                size = 6) +
            theme_bw()
    if(!is.null(plot_output_dir)){
        ggsave(file.path(plot_output_dir, glue("{prefix}distribution_{val_col}_{chose_pair}.png")), p)
    }
    return(p)
}