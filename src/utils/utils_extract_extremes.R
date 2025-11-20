suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})
#' Bin expression values by gene across samples
#'
#' @param dat data.table with first column "fly_id" (gene IDs)
#'            and remaining columns expression values per sample
#' @param n integer, number of bins to divide expression values into (default = 5)
#'
#' @return A long-format data.table with columns: fly_id, sample, logCPM, and tile
#' @export
bin_expr <- function(dat, n = 5) {
  # Input checks
  stopifnot(is.data.table(dat))
  if (!"fly_id" %in% names(dat))
    stop("dat must have a column named 'fly_id'")
  if (n <= 1) stop("'n' must be greater than 1")
  # Melt into long format
  dat_l <- melt(
    dat,
    id.vars = "fly_id",
    variable.name = "sample",
    value.name = "logCPM"
  )
  # Bin expression within each gene (fly_id)
  dat_l[, tile := dplyr::ntile(logCPM, n), by = "fly_id"]
  setcolorder(dat_l, c("fly_id", "sample", "logCPM", "tile"))
  dat_l[]
}
#' Extract expression matrix for a specific expression bin
#'
#' @param binned_expr data.table from \code{bin_expr()} (must include fly_id, sample, tile)
#' @param full_expr full expression data.table with "fly_id" as the first column
#' @param nth_tile integer, the tile number to extract (e.g. 1 = lowest, n = highest)
#'
#' @return A numeric matrix of expression values (samples × genes)
#' @export
extract_biase_expr <- function(binned_expr, full_expr, nth_tile) {
  req_cols <- c("fly_id", "tile", "sample")
  if (!all(req_cols %in% colnames(binned_expr)))
    stop("binned_expr must contain columns: fly_id, tile, sample")
  if (!"fly_id" %in% names(full_expr))
    stop("full_expr must contain a column 'fly_id'")
  # Get samples belonging to the target tile
  target_samples <- binned_expr[tile == nth_tile, sample] |> as.character() |> unique()
  # Subset the full expression table
  exp_mat <- full_expr[, .SD, .SDcols = c("fly_id", target_samples)] |>
    column_to_rownames("fly_id") |>
    t()
  return(exp_mat)
}