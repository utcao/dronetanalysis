# ==============================================================================
# limma_voom_deg_core.R
#
# Defines run_limma_voom_deg(): a self-contained function that takes a
# pre-aligned count matrix + annotated metadata and returns a per-gene
# DEG classification table using limma-voom.
#
# INPUTS: standardised — no file I/O, no config.yaml, no preprocessing logic.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glue)
  library(stringr)
  library(purrr)
  library(edgeR)
  library(limma)
  library(foreach)
})

# ------------------------------------------------------------------------------
# run_limma_voom_deg
#
# Arguments:
#   count_mat       matrix; rownames = gene IDs, colnames = sample IDs (integers)
#   meta_dt         data.table; one row per sample — must contain sample_id_col,
#                   trt_col, and every covariate in formula_full
#   trt_col         column name in meta_dt holding treatment labels
#   ref_group       reference level value (must be present in meta_dt[[trt_col]])
#   formula_full    formula string for model.matrix, e.g.
#                   "~ expr_quintile + RNAlibBatch + ... + sv4"
#   mu_sum_covars   character vector of batch covariates to apply sum-to-zero coding;
#                   applied to model.matrix via contrasts() before model.matrix() call
#   sample_id_col   column in meta_dt whose values match colnames(count_mat)
#   fdr_threshold   BH-adjusted p-value cutoff for DEG classification
#   n_workers       kept for API symmetry; lmFit is vectorized and does not use it
#   out_prefix      if non-empty, write result CSVs to out_dir
#   out_dir         output directory (used only when out_prefix is given)
#
# Returns:
#   data.table with one row per gene and columns for AveExpr, logFC, t,
#   pvalue, adj_pvalue, B, and class (up_deg / down_deg / NS).
#   Columns suffixed with the treatment level.
# ------------------------------------------------------------------------------
run_limma_voom_deg <- function(
  count_mat,
  meta_dt,
  trt_col,
  ref_group,
  formula_full,
  mu_sum_covars = NULL,
  sample_id_col = "id",
  fdr_threshold = 0.05,
  n_workers     = 1,
  out_prefix    = "",
  out_dir       = "."
) {
  # --- parse all variables from formula ---
  parse_formula_vars <- function(f) {
    vars <- unlist(str_split(f, "~|\\+|\\s+"))
    vars[!(vars %in% c("", "0", "1"))]
  }
  all_model_vars <- unique(parse_formula_vars(formula_full))

  # --- prepare metadata: factors, sum-to-zero, treatment reference ---
  # contrasts() must be set on factor columns BEFORE model.matrix() because
  # model.matrix() reads the contrast attribute at call time.
  meta_df <- as.data.frame(copy(meta_dt))
  rownames(meta_df) <- meta_df[[sample_id_col]]

  for (col in all_model_vars) {
    if (col %in% names(meta_df) && !is.numeric(meta_df[[col]])) {
      meta_df[[col]] <- as.factor(meta_df[[col]])
    }
  }
  if (!is.null(mu_sum_covars)) {
    for (col in mu_sum_covars) {
      if (col %in% names(meta_df)) meta_df[[col]] <- as.factor(meta_df[[col]])
    }
  }

  meta_df[[trt_col]] <- relevel(factor(meta_df[[trt_col]]), ref = ref_group)

  if (!is.null(mu_sum_covars)) {
    for (col in mu_sum_covars) {
      if (col %in% names(meta_df))
        contrasts(meta_df[[col]]) <- contr.sum(nlevels(meta_df[[col]]))
    }
  }

  design <- model.matrix(as.formula(formula_full), data = meta_df)

  # --- align count_mat ---
  count_mat <- count_mat[, rownames(meta_df), drop = FALSE]

  # --- verify treatment coefficient name ---
  trt_groups <- setdiff(levels(meta_df[[trt_col]]), ref_group)
  # model.matrix() with treatment coding produces e.g. "expr_quintileQ5" (no separator)
  result_coef <- paste0(trt_col, trt_groups)
  missing <- setdiff(result_coef, colnames(design))
  if (length(missing) > 0) {
    stop(glue(
      "Expected coefficient(s) not found in design matrix:\n",
      "  missing: {paste(missing, collapse=', ')}\n",
      "  available: {paste(colnames(design), collapse=', ')}"
    ))
  }

  # --- TMM normalisation + voom ---
  dge      <- DGEList(counts = count_mat)
  dge      <- calcNormFactors(dge, method = "TMM")
  voom_obj <- voom(dge, design = design, plot = FALSE)

  # --- lmFit → eBayes ---
  fit <- lmFit(voom_obj, design)
  fit <- eBayes(fit)

  # --- per-treatment topTable ---
  classify_deg <- function(q, logfc) {
    dplyr::case_when(
      !is.na(q) & q < fdr_threshold & logfc > 0  ~ "up_deg",
      !is.na(q) & q < fdr_threshold & logfc <= 0 ~ "down_deg",
      TRUE                                         ~ "NS"
    )
  }

  result_list <- foreach(trt = trt_groups, coef = result_coef) %do% {
    tt <- topTable(fit, coef = coef, number = Inf,
                   adjust.method = "BH", sort.by = "none")

    result_dt <- data.table(
      gene_id    = rownames(tt),
      AveExpr    = tt$AveExpr,
      logFC      = tt$logFC,
      t          = tt$t,
      pvalue     = tt$P.Value,
      adj_pvalue = tt$adj.P.Val,
      B          = tt$B,
      class      = classify_deg(tt$adj.P.Val, tt$logFC)
    )
    trt_result_cols <- setdiff(names(result_dt), "gene_id")
    setnames(result_dt, trt_result_cols, paste0(trt_result_cols, "_", trt))
    result_dt
  }

  # --- merge all treatment comparisons ---
  deg_tab <- purrr::reduce(result_list,
    function(x, y) merge(x, y, by = "gene_id", all = TRUE))

  # --- optional output ---
  if (nchar(out_prefix) > 0) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    fwrite(deg_tab,
      file = file.path(out_dir, glue("{out_prefix}_limma_voom_deg.csv")))

    class_cols <- grep("^class_", names(deg_tab), value = TRUE)
    sum_list   <- lapply(class_cols, function(col) {
      as.data.frame(t(unclass(table(deg_tab[[col]]))), check.names = FALSE)
    })
    sum_dt <- data.table(
      comparison = gsub("^class_", "", class_cols),
      rbindlist(sum_list, fill = TRUE)
    )
    sum_dt[, total := rowSums(.SD, na.rm = TRUE), .SDcols = is.numeric]
    fwrite(sum_dt,
      file = file.path(out_dir, glue("{out_prefix}_limma_voom_deg_sum.csv")))
  }

  invisible(deg_tab)
}
