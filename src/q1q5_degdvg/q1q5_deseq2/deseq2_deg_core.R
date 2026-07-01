# ==============================================================================
# deseq2_deg_core.R
#
# Defines run_deseq2_deg(): a self-contained function that takes a
# pre-aligned count matrix + annotated metadata and returns a per-gene
# DEG classification table using DESeq2.
#
# INPUTS: standardised — no file I/O, no config.yaml, no preprocessing logic.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glue)
  library(stringr)
  library(purrr)
  library(DESeq2)
  library(BiocParallel)
  library(foreach)
})

# ------------------------------------------------------------------------------
# run_deseq2_deg
#
# Arguments:
#   count_mat       matrix; rownames = gene IDs, colnames = sample IDs (integers)
#   meta_dt         data.table; one row per sample — must contain sample_id_col,
#                   trt_col, and every covariate in formula_full
#   trt_col         column name in meta_dt holding treatment labels
#   ref_group       reference level value (must be present in meta_dt[[trt_col]])
#   formula_full    formula string for DESeq2 design, e.g.
#                   "~ expr_quintile + RNAlibBatch + ... + sv4"
#   formula_reduced reduced formula for LRT; NULL when test = "Wald"
#   mu_sum_covars   character vector of batch covariates to apply sum-to-zero coding;
#                   must be applied BEFORE DESeqDataSetFromMatrix() construction
#   sample_id_col   column in meta_dt whose values match colnames(count_mat)
#   fdr_threshold   BH-adjusted p-value cutoff for DEG classification
#   lfc_shrink      logical; apply lfcShrink(type="normal") to shrunken LFC column
#   test            "Wald" (default) or "LRT"
#   n_workers       workers for BiocParallel (MulticoreParam when >1)
#   out_prefix      if non-empty, write result CSVs to out_dir
#   out_dir         output directory (used only when out_prefix is given)
#
# Returns:
#   data.table with one row per gene and columns for baseMean, log2FC,
#   log2FC_shrunken (NA when lfc_shrink=FALSE), lfcSE, stat, pvalue, padj,
#   and class (up_deg / down_deg / NS). Columns suffixed with the treatment level.
# ------------------------------------------------------------------------------
run_deseq2_deg <- function(
  count_mat,
  meta_dt,
  trt_col,
  ref_group,
  formula_full,
  formula_reduced = NULL,
  mu_sum_covars   = NULL,
  sample_id_col   = "id",
  fdr_threshold   = 0.05,
  lfc_shrink      = TRUE,
  test            = "Wald",
  n_workers       = 1,
  out_prefix      = "",
  out_dir         = "."
) {
  # --- parse all variables from formula ---
  parse_formula_vars <- function(f) {
    vars <- unlist(str_split(f, "~|\\+|\\s+"))
    vars[!(vars %in% c("", "0", "1"))]
  }
  all_model_vars <- unique(parse_formula_vars(formula_full))

  # --- prepare colData: factors, sum-to-zero, treatment reference ---
  # Contrasts MUST be stamped on factor columns BEFORE DESeqDataSetFromMatrix().
  # DESeq2 copies colData at construction — late contrasts() calls have no effect.
  meta_df <- as.data.frame(copy(meta_dt))
  rownames(meta_df) <- meta_df[[sample_id_col]]

  for (col in all_model_vars) {
    if (col %in% names(meta_df) && !is.numeric(meta_df[[col]])) {
      meta_df[[col]] <- as.factor(meta_df[[col]])
    }
  }
  # integer-coded batch columns pass is.numeric(); coerce mu_sum_covars explicitly
  if (!is.null(mu_sum_covars)) {
    for (col in mu_sum_covars) {
      if (col %in% names(meta_df)) meta_df[[col]] <- as.factor(meta_df[[col]])
    }
  }

  # treatment reference level
  meta_df[[trt_col]] <- relevel(factor(meta_df[[trt_col]]), ref = ref_group)

  # sum-to-zero contrast coding for batch covariates
  if (!is.null(mu_sum_covars)) {
    for (col in mu_sum_covars) {
      if (col %in% names(meta_df))
        contrasts(meta_df[[col]]) <- contr.sum(nlevels(meta_df[[col]]))
    }
  }

  # --- align count_mat and ensure integer mode ---
  count_mat <- count_mat[, rownames(meta_df), drop = FALSE]
  mode(count_mat) <- "integer"

  # --- build full-rank design matrix ---
  # DESeq2's checkFullRank() hard-fails if any batch level appears exclusively in
  # Q1 or Q5 (a common occurrence after quintile filtering). Build the design
  # matrix with model.matrix() first and drop linearly dependent columns via QR
  # decomposition — the same thing model.matrix()/lm() does silently for limma.
  # When columns are dropped, DESeq2 uses the matrix directly and coefficient names
  # come from colnames(design_mat), e.g. "expr_quintileQ5" (no "_vs_" separator).
  design_mat <- model.matrix(as.formula(formula_full), data = meta_df)
  qr_decomp  <- qr(design_mat)
  if (qr_decomp$rank < ncol(design_mat)) {
    kept_idx   <- qr_decomp$pivot[seq_len(qr_decomp$rank)]
    dropped    <- colnames(design_mat)[-kept_idx]
    message("Dropped rank-deficient column(s) from design matrix: ",
            paste(dropped, collapse = ", "))
    design_mat <- design_mat[, kept_idx, drop = FALSE]
  }

  # Verify treatment coefficient survived rank reduction
  trt_groups  <- setdiff(levels(meta_df[[trt_col]]), ref_group)
  result_coef <- paste0(trt_col, trt_groups)   # e.g. "expr_quintileQ5"
  missing_trt <- setdiff(result_coef, colnames(design_mat))
  if (length(missing_trt) > 0) {
    stop(glue(
      "Treatment coefficient dropped during rank reduction: ",
      "{paste(missing_trt, collapse=', ')}. ",
      "Q1/Q5 groups are perfectly confounded with batch structure — ",
      "cannot test for DEG with this focal gene × condition."
    ))
  }

  # --- construct DESeqDataSet ---
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = meta_df,
    design    = design_mat
  )

  # --- run DESeq2 ---
  bpp <- if (n_workers > 1) MulticoreParam(n_workers) else SerialParam()
  if (test == "Wald") {
    dds <- DESeq(dds, test = "Wald", parallel = (n_workers > 1), BPPARAM = bpp,
                 quiet = TRUE)
  } else {
    if (is.null(formula_reduced))
      stop("formula_reduced required when test = 'LRT'")
    # For LRT with a design matrix, build the reduced matrix the same way
    design_red <- model.matrix(as.formula(formula_reduced), data = meta_df)
    qr_red     <- qr(design_red)
    if (qr_red$rank < ncol(design_red))
      design_red <- design_red[, qr_red$pivot[seq_len(qr_red$rank)], drop = FALSE]
    dds <- DESeq(dds, test = "LRT", reduced = design_red,
                 parallel = (n_workers > 1), BPPARAM = bpp, quiet = TRUE)
  }

  # --- per-treatment extraction ---
  available_coefs <- resultsNames(dds)
  missing <- setdiff(result_coef, available_coefs)
  if (length(missing) > 0) {
    stop(glue(
      "Expected coefficient(s) not found in resultsNames(dds):\n",
      "  missing: {paste(missing, collapse=', ')}\n",
      "  available: {paste(available_coefs, collapse=', ')}"
    ))
  }

  classify_deg <- function(q, logfc) {
    dplyr::case_when(
      !is.na(q) & q < fdr_threshold & logfc > 0  ~ "up_deg",
      !is.na(q) & q < fdr_threshold & logfc <= 0 ~ "down_deg",
      TRUE                                         ~ "NS"
    )
  }

  result_list <- foreach(trt = trt_groups, coef = result_coef) %do% {
    res_raw <- results(dds, name = coef, alpha = fdr_threshold)

    log2FC_shrunken <- if (lfc_shrink) {
      # type="normal" only — apeglm/ashr are not available in the ganlss conda env
      suppressMessages(
        lfcShrink(dds, coef = coef, type = "normal", quiet = TRUE)$log2FoldChange
      )
    } else {
      rep(NA_real_, nrow(res_raw))
    }

    result_dt <- data.table(
      gene_id         = rownames(res_raw),
      baseMean        = as.numeric(res_raw$baseMean),
      log2FC          = as.numeric(res_raw$log2FoldChange),
      log2FC_shrunken = log2FC_shrunken,
      lfcSE           = as.numeric(res_raw$lfcSE),
      stat            = as.numeric(res_raw$stat),
      pvalue          = as.numeric(res_raw$pvalue),
      padj            = as.numeric(res_raw$padj),
      class           = classify_deg(as.numeric(res_raw$padj),
                                     as.numeric(res_raw$log2FoldChange))
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
      file = file.path(out_dir, glue("{out_prefix}_deseq2_deg.csv")))

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
      file = file.path(out_dir, glue("{out_prefix}_deseq2_deg_sum.csv")))
  }

  invisible(deg_tab)
}
