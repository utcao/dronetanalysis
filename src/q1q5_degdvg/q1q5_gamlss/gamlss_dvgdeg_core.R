# ==============================================================================
# gamlss_dvgdeg_core.R
#
# Defines run_gamlss_dvgdeg(): a self-contained function that takes a
# pre-aligned count matrix + annotated metadata and returns a per-gene
# DEG/DVG classification table.
#
# DEPENDENCIES: tools_gamlss.R must be sourced before calling this function.
# INPUTS: standardised — no file I/O, no config.yaml, no preprocessing logic.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glue)
  library(stringr)
  library(purrr)
  library(edgeR)
  library(gamlss)
  library(gamlss.dist)
  library(BiocParallel)
  library(foreach)
  library(futile.logger)
})

# ------------------------------------------------------------------------------
# run_gamlss_dvgdeg
#
# Arguments:
#   count_mat      matrix; rownames = gene IDs, colnames = sample IDs
#   meta_dt        data.table; one row per sample — must contain sample_id_col,
#                  trt_col, and every covariate mentioned in the three formulas
#   trt_col        column name in meta_dt holding treatment labels
#   ref_group      reference level value (must be present in meta_dt[[trt_col]])
#   formula_mu     formula string for the mean; must include trt_col, e.g.
#                  "~0 + expr_quintile + RNAseqBatch + sv1"
#   formula_sigma  formula string for overdispersion; typically "~0 + trt_col"
#   formula_nontrt formula_mu without the trt_col term (used for LR test on mu)
#   mu_sum_covars  character vector of covariates to apply sum-to-zero coding;
#                  pass NULL to use default treatment contrasts
#   sample_id_col  column in meta_dt whose values match colnames(count_mat)
#   fdr_threshold  BH-adjusted p-value cutoff for DEG/DVG classification
#   n_workers      number of parallel SnowParam workers
#   out_prefix     if non-empty, write result CSVs to gamlss_dir
#   gamlss_dir     output directory (used only when out_prefix is given)
#
# Returns:
#   data.table with one row per gene and columns for mu, sigma, BCV,
#   p/q values (LR test and single-model t-test), and classification
#   (DEG / DVG / Both / NS). Column names are suffixed with the treatment level.
# ------------------------------------------------------------------------------
run_gamlss_dvgdeg <- function(
  count_mat,
  meta_dt,
  trt_col,
  ref_group,
  formula_mu,
  formula_sigma,
  formula_nontrt,
  mu_sum_covars  = NULL,
  sample_id_col  = "id",
  fdr_threshold  = 0.05,
  n_workers      = 1,
  out_prefix     = "",
  gamlss_dir     = "."
) {
  # --- parse covariates from all three formulas ---
  parse_formula_vars <- function(f) {
    vars <- unlist(str_split(f, "~|\\+|\\s+"))
    vars[!(vars %in% c("", "0", "1"))]
  }
  all_covars <- unique(c(
    parse_formula_vars(formula_mu),
    parse_formula_vars(formula_sigma),
    parse_formula_vars(formula_nontrt)
  ))

  # Coerce non-numeric covariates to factor; leave SVs (already numeric) untouched
  meta_dt <- copy(meta_dt)
  for (col in all_covars) {
    if (col %in% names(meta_dt) && !is.numeric(meta_dt[[col]])) {
      meta_dt[[col]] <- as.factor(meta_dt[[col]])
    }
    if(col %in% mu_sum_covars) {
      meta_dt[[col]] <- as.factor(meta_dt[[col]])
    }
  }
  # mu_sum_covars must be factors for contr.sum regardless of storage type:
  # integer-coded batch columns pass is.numeric() and would be missed above.
  if (!is.null(mu_sum_covars)) {
    for (col in mu_sum_covars) {
      if (col %in% names(meta_dt)) meta_dt[[col]] <- as.factor(meta_dt[[col]])
    }
  }

  # Set treatment factor levels: ref_group first so it becomes the intercept
  group_level     <- unique(as.character(meta_dt[[trt_col]]))
  meta_dt[[trt_col]] <- factor(meta_dt[[trt_col]],
    levels = c(ref_group, setdiff(group_level, ref_group)))

  # Align count_mat columns to meta_dt row order
  sample_order <- meta_dt[[sample_id_col]]
  count_mat    <- count_mat[, sample_order, drop = FALSE]

  # --- TMM normalisation ---
  # offsets = log(lib.size * norm.factor), passed as fixed term in the log-link
  dge_obj    <- DGEList(count_mat)
  dge_obj    <- calcNormFactors(dge_obj, method = "TMM")
  tmm_counts <- dge_obj$counts
  n_genes    <- nrow(tmm_counts)
  offsets    <- log(dge_obj$samples$lib.size * dge_obj$samples$norm.factors)

  # --- fit three NBI models per gene ---
  # MulticoreParam (fork-based) is far more stable than SnowParam (socket-based)
  # for long-running GAMLSS jobs on Linux. SerialParam avoids all IPC overhead
  # when only one worker is requested.
  bpp <- if (n_workers > 1) MulticoreParam(n_workers) else SerialParam()

  # suppress summary() printing — summary.gamlss prints as a side effect
  fit_sum_quiet <- function(m) {
    if (is.null(m)) return(NULL)
    result <- NULL
    utils::capture.output(result <- summary(m))
    result
  }

  # m_full: all covariates in mu, treatment in sigma
  m_full <- bplapply(seq_len(n_genes), fit_one,
    dat           = tmm_counts,
    meta          = meta_dt,
    mu_formula    = formula_mu,
    sigma_formula = formula_sigma,
    offsets       = offsets,
    mu_sum_vars   = mu_sum_covars,
    BPPARAM       = bpp)

  m_full_sum <- bplapply(m_full, fit_sum_quiet, BPPARAM = bpp)

  # m_sigma_int: intercept-only sigma — LR test isolates treatment effect on sigma
  m_sigma_int <- bplapply(seq_len(n_genes), fit_one,
    dat           = tmm_counts,
    meta          = meta_dt,
    mu_formula    = formula_mu,
    sigma_formula = "~ 1",
    offsets       = offsets,
    mu_sum_vars   = mu_sum_covars,
    BPPARAM       = bpp)

  # m_mu_nontrt: non-treatment mu — LR test isolates treatment effect on mu
  m_mu_nontrt <- bplapply(seq_len(n_genes), fit_one,
    dat           = tmm_counts,
    meta          = meta_dt,
    mu_formula    = formula_nontrt,
    sigma_formula = formula_sigma,
    offsets       = offsets,
    mu_sum_vars   = mu_sum_covars,
    BPPARAM       = bpp)

  # --- per-treatment extraction ---
  # m_full_sum rows: combined mu + sigma coefficient table from summary.gamlss
  coef_names <- rownames(m_full_sum[[1]])
  trt_groups <- setdiff(levels(meta_dt[[trt_col]]), ref_group)

  # Directional classification: direction encoded in the label.
  # logfc_mu  > 0 → higher mean in treatment (up_deg)
  # logfc_bcv > 0 → higher dispersion/BCV in treatment (up_dvg)
  # "Both" genes (both mu and sigma significant) can be further examined
  # via the logFC_mu and logFC_BCV columns in the output table.
  classify_gamlss <- function(q_mean, q_var, logfc_mu, logfc_bcv) {
    dplyr::case_when(
      q_mean < fdr_threshold & q_var >= fdr_threshold & logfc_mu > 0  ~ "up_deg",
      q_mean < fdr_threshold & q_var >= fdr_threshold & logfc_mu <= 0 ~ "down_deg",
      q_mean >= fdr_threshold & q_var < fdr_threshold & logfc_bcv > 0  ~ "up_dvg",
      q_mean >= fdr_threshold & q_var < fdr_threshold & logfc_bcv <= 0 ~ "down_dvg",
      q_mean < fdr_threshold & q_var < fdr_threshold                   ~ "Both",
      TRUE                                                              ~ "NS"
    )
  }

  result_list <- foreach(trt = trt_groups) %do% {
    # coefficient name pattern, e.g. "expr_quintileQ5"
    coef_trt <- paste0(trt_col, trt)
    coef_idx <- grep(coef_trt, coef_names)
    stopifnot(
      "expected exactly 2 coefficient matches (one mu, one sigma)" =
        length(coef_idx) == 2
    )

    # single-model t-tests from the m_full summary table
    p_mu_single    <- vapply(m_full_sum,
      function(m) ifelse(length(m), m[coef_idx[1], "Pr(>|t|)"], NA_real_), numeric(1))
    p_sigma_single <- vapply(m_full_sum,
      function(m) ifelse(length(m), m[coef_idx[2], "Pr(>|t|)"], NA_real_), numeric(1))

    # expected CPM at the mean offset
    cpm_count_ctrl <- vapply(m_full,
      function(m) ifelse(length(m),
        exp(m$mu.coefficients[1] + mean(offsets)), NA_real_), numeric(1))
    cpm_count_trt  <- vapply(m_full,
      function(m) ifelse(length(m),
        exp(m$mu.coefficients[coef_trt] + m$mu.coefficients[1] + mean(offsets)),
        NA_real_), numeric(1))

    mu_ctrl <- vapply(m_full,
      function(m) ifelse(length(m), m$mu.coefficients[1], NA_real_), numeric(1))
    mu_trt  <- vapply(m_full,
      function(m) ifelse(length(m),
        m$mu.coefficients[coef_trt] + m$mu.coefficients[1], NA_real_), numeric(1))

    logFC_mu  <- log2(mu_trt  / mu_ctrl)
    logFC_cpm <- log2(cpm_count_trt / cpm_count_ctrl)

    # LR test: treatment effect on mu (m_full vs m_mu_nontrt)
    p_mu <- map2_dbl(m_full, m_mu_nontrt,
      function(m_ful, m_red) {
        if (length(m_ful) && length(m_red)) {
          pchisq(2 * (logLik(m_ful) - logLik(m_red)),
                 df = m_ful$df.fit - m_red$df.fit, lower.tail = FALSE)[[1]]
        } else NA_real_
      })

    # BCV = sqrt(exp(sigma)); NBI parameterisation: sigma = log(dispersion)
    BCV_ctrl <- vapply(m_full,
      function(m) ifelse(length(m), sqrt(exp(m$sigma.coefficients[1])), NA_real_), numeric(1))
    BCV_trt  <- vapply(m_full,
      function(m) ifelse(length(m),
        sqrt(exp(m$sigma.coefficients[1] + m$sigma.coefficients[coef_trt])),
        NA_real_), numeric(1))

    sigma_ctrl <- vapply(m_full,
      function(m) ifelse(length(m), m$sigma.coefficients[1], NA_real_), numeric(1))
    sigma_trt  <- vapply(m_full,
      function(m) ifelse(length(m),
        m$sigma.coefficients[1] + m$sigma.coefficients[coef_trt], NA_real_), numeric(1))

    logFC_BCV <- log2(BCV_trt / BCV_ctrl)
    diff_BCV  <- BCV_trt - BCV_ctrl

    # LR test: treatment effect on sigma (m_full vs m_sigma_int)
    p_sigma <- map2_dbl(m_full, m_sigma_int,
      function(m_ful, m_red) {
        if (length(m_ful) && length(m_red)) {
          pchisq(2 * (logLik(m_ful) - logLik(m_red)),
                 df = m_ful$df.fit - m_red$df.fit, lower.tail = FALSE)[[1]]
        } else NA_real_
      })

    q_mu           <- p.adjust(p_mu,           method = "BH")
    q_sigma        <- p.adjust(p_sigma,        method = "BH")
    q_mu_single    <- p.adjust(p_mu_single,    method = "BH")
    q_sigma_single <- p.adjust(p_sigma_single, method = "BH")

    if (n_genes != length(q_mu)) {
      flog.error("q_mu length does not match n_genes")
    }

    result_dt <- data.table(
      gene_id        = rownames(tmm_counts)[seq_len(n_genes)],
      cpm_count_ctrl = cpm_count_ctrl,
      cpm_count_trt  = cpm_count_trt,
      mu_ctrl        = mu_ctrl,
      mu_trt         = mu_trt,
      logFC_mu       = logFC_mu,
      logFC_cpm      = logFC_cpm,
      BCV_ctrl       = BCV_ctrl,
      BCV_trt        = BCV_trt,
      logFC_BCV      = logFC_BCV,
      diff_BCV       = diff_BCV,
      sigma_ctrl     = sigma_ctrl,
      sigma_trt      = sigma_trt,
      p_mu           = p_mu,
      p_sigma        = p_sigma,
      q_mu           = q_mu,
      q_sigma        = q_sigma,
      p_mu_single    = p_mu_single,
      p_sigma_single = p_sigma_single,
      q_mu_single    = q_mu_single,
      q_sigma_single = q_sigma_single,
      class_single    = classify_gamlss(q_mu_single,  q_sigma_single,  logFC_mu, logFC_BCV),
      class_modelcomp = classify_gamlss(q_mu,         q_sigma,         logFC_mu, logFC_BCV)
    )

    # suffix all columns except gene_id with the treatment level name
    trt_result_cols <- setdiff(names(result_dt), "gene_id")
    setnames(result_dt, trt_result_cols, paste0(trt_result_cols, "_", trt))
    result_dt
  }

  # --- merge all treatment comparisons ---
  dvdg_tab <- purrr::reduce(result_list,
    function(x, y) merge(x, y, by = "gene_id", all = TRUE))

  # --- optional output ---
  if (nchar(out_prefix) > 0) {
    fwrite(dvdg_tab,
      file = file.path(gamlss_dir,
        glue("{out_prefix}_gamlss_dvgdeg.csv")))

    class_cols  <- grep("class", names(dvdg_tab), value = TRUE)
    sum_class   <- purrr::map_df(class_cols, function(col) {
      table(dvdg_tab[[col]]) |> unclass()
    })
    sum_class_dt <- dplyr::mutate(as.data.frame(sum_class),
        comparison = gsub("class_", "", class_cols, fixed = TRUE)) |>
      dplyr::select(comparison, dplyr::everything()) |>
      as.data.table()
    sum_class_dt[, total := rowSums(.SD, na.rm = TRUE), .SDcols = is.numeric]

    fwrite(sum_class_dt,
      file = file.path(gamlss_dir,
        glue("{out_prefix}_gamlss_dvgdeg_sum.csv")))
  }

  invisible(dvdg_tab)
}
