#!/usr/bin/env Rscript
# ==============================================================================
# Permutation Test for Gene-Level MAD Transcriptomic Noise
#
# For each focus gene G, the observed delta_mean_mad (mean(MAD_HIGH) -
# mean(MAD_LOW) across all other genes) is compared against a null distribution
# built by permuting G's expression values across samples.
#
# HOW THE LINK IS BROKEN:
#   In the original analysis, LOW/HIGH groups are defined by G's expression rank
#   across samples — the biological hypothesis being that G's level co-segregates
#   with global transcriptomic variability.  The permutation shuffles G's
#   expression values across samples before re-ranking.  Because the shuffle
#   breaks the mapping between G's values and the actual samples, the re-derived
#   LOW/HIGH split is effectively a random partition of equal size (k_low /
#   k_high).  The null distribution answers: "how large a MAD difference arises
#   under ANY random split of this size — not one driven by G's biology?"
#
#   A focus gene whose observed |delta| exceeds the null at the chosen alpha is
#   genuinely organising samples into high-variability vs. low-variability states.
#
# PERMUTATION p-value (two-sided):
#   perm_p = #{|null_delta| >= |obs_delta|} / n_perms
#   (replicates the logic of permutate_p_two() in utils_permutation_net.R)
#
# OUTPUTS:
#   - XLSX: observed stats + permutation p-values, BH-FDR-corrected
#   - Optional TSV.GZ per gene: null distribution of delta_mean_mad
#     (compatible with plot_permutation_null_dist.R)
#   - Optional permutation_pvals.tsv for plot_permutation_null_dist.R
#
# Usage:
#   Rscript compute_mad_permute_transcriptomic_noise.R \
#     --expr-file    data/processed/VOOM/voomdataCtrl.txt \
#     --mapping-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
#     --output-file  results/variability/mad_permutation_ct.xlsx \
#     --n-perms      500 \
#     --seed         42 \
#     --condition-label "Control"
#
# TIP: Use --focus-genes to restrict permutation to a subset (e.g. originally
#      significant genes) to reduce compute time.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
  library(openxlsx)
})

# ----- Helper: significance stars -----
get_sig_stars <- function(p) {
  ifelse(is.na(p), NA_character_,
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01,  "**",
        ifelse(p < 0.05,  "*", "ns"))))
}

# ----- Row-wise MAD: use matrixStats if available, else apply fallback -----
row_mads <- if (requireNamespace("matrixStats", quietly = TRUE)) {
  function(mat) matrixStats::rowMads(mat, na.rm = TRUE)
} else {
  function(mat) apply(mat, 1, mad, na.rm = TRUE)
}

# ----- Core: permutation test for one focus gene -----
run_permutation_test <- function(focus_gene, expr_mat, gene_ids,
                                 k_low, k_high, n_perms) {
  tryCatch({
    focus_idx  <- which(gene_ids == focus_gene)
    focus_vec  <- expr_mat[focus_idx, ]        # length = N samples
    other_idx  <- which(gene_ids != focus_gene)
    expr_other <- expr_mat[other_idx, , drop = FALSE]

    # ---- Observed LOW/HIGH groups (expression-rank-based) ----
    obs_low  <- order(focus_vec)[seq_len(k_low)]
    obs_high <- order(focus_vec, decreasing = TRUE)[seq_len(k_high)]

    mad_obs_low  <- row_mads(expr_other[, obs_low,  drop = FALSE])
    mad_obs_high <- row_mads(expr_other[, obs_high, drop = FALSE])

    obs_delta      <- mean(mad_obs_high, na.rm = TRUE) - mean(mad_obs_low, na.rm = TRUE)
    obs_mean_low   <- mean(mad_obs_low,  na.rm = TRUE)
    obs_mean_high  <- mean(mad_obs_high, na.rm = TRUE)

    # Parametric Wilcoxon on observed groups (for direct comparison with
    # compute_mad_transcriptomic_variability.R output)
    wt           <- wilcox.test(mad_obs_low, mad_obs_high, paired = FALSE, exact = FALSE)
    obs_wilcox_p <- wt$p.value

    # ---- Permutation null ----
    # Shuffle focus_vec across samples, re-rank → random LOW/HIGH partition.
    # seed is already set globally before calling this function, so each gene
    # gets a reproducible but independent sequence of shuffles.
    null_deltas <- vapply(seq_len(n_perms), function(p) {
      focus_perm <- sample(focus_vec)
      perm_low   <- order(focus_perm)[seq_len(k_low)]
      perm_high  <- order(focus_perm, decreasing = TRUE)[seq_len(k_high)]
      mad_l <- row_mads(expr_other[, perm_low,  drop = FALSE])
      mad_h <- row_mads(expr_other[, perm_high, drop = FALSE])
      mean(mad_h, na.rm = TRUE) - mean(mad_l, na.rm = TRUE)
    }, FUN.VALUE = numeric(1L))

    # Two-sided permutation p-value
    perm_p <- sum(abs(null_deltas) >= abs(obs_delta)) / n_perms

    list(
      focus_gene         = focus_gene,
      obs_delta_mean_mad = obs_delta,
      obs_mean_mad_low   = obs_mean_low,
      obs_mean_mad_high  = obs_mean_high,
      obs_wilcoxon_p     = obs_wilcox_p,
      perm_p             = perm_p,
      n_perms            = n_perms,
      null_deltas        = null_deltas
    )
  }, error = function(e) {
    list(
      focus_gene         = focus_gene,
      obs_delta_mean_mad = NA_real_,
      obs_mean_mad_low   = NA_real_,
      obs_mean_mad_high  = NA_real_,
      obs_wilcoxon_p     = NA_real_,
      perm_p             = NA_real_,
      n_perms            = n_perms,
      null_deltas        = NULL
    )
  })
}

# ----- CLI arguments -----
parser <- ArgumentParser(
  description = "Permutation test for gene-level MAD transcriptomic noise"
)
parser$add_argument("--expr-file",
  required = TRUE,
  help = "VOOM expression matrix (tab-separated, genes x samples)"
)
parser$add_argument("--mapping-file",
  required = TRUE,
  help = "TSV/CSV with gene_id and SYMBOL columns"
)
parser$add_argument("--output-file",
  required = TRUE,
  help = "Output .xlsx file path"
)
parser$add_argument("--null-dist-dir",
  default = NULL,
  help = paste(
    "Optional: directory to write per-gene *_null_dist.tsv.gz files.",
    "Also writes permutation_pvals.tsv for use with plot_permutation_null_dist.R"
  )
)
parser$add_argument("--focus-genes",
  default = NULL,
  help = paste(
    "Optional: comma-separated gene IDs to restrict permutation.",
    "Useful to limit compute to originally-significant genes."
  )
)
parser$add_argument("--low-frac",
  type = "double", default = 0.2,
  help = "Fraction of samples in LOW group (default: 0.2)"
)
parser$add_argument("--high-frac",
  type = "double", default = 0.2,
  help = "Fraction of samples in HIGH group (default: 0.2)"
)
parser$add_argument("--n-perms",
  type = "integer", default = 500,
  help = "Number of permutations per gene (default: 500)"
)
parser$add_argument("--seed",
  type = "integer", default = 42,
  help = "Random seed for reproducibility (default: 42)"
)
parser$add_argument("--condition-label",
  default = "Condition",
  help = "Human-readable condition label added as a column"
)
args <- parser$parse_args()

cat("=== Permutation Test: Gene-Level MAD Transcriptomic Noise ===\n")
cat("Expression file:  ", args$expr_file,       "\n")
cat("Mapping file:     ", args$mapping_file,     "\n")
cat("Output file:      ", args$output_file,      "\n")
cat("Permutations:     ", args$n_perms,          "\n")
cat("Seed:             ", args$seed,             "\n")
cat("LOW fraction:     ", args$low_frac,         "\n")
cat("HIGH fraction:    ", args$high_frac,        "\n")
cat("Condition label:  ", args$condition_label,  "\n\n")

# ----- Validate inputs -----
if (!file.exists(args$expr_file))
  stop("Expression file not found: ", args$expr_file)
if (!file.exists(args$mapping_file))
  stop("Mapping file not found: ", args$mapping_file)

out_dir <- dirname(args$output_file)
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

if (!is.null(args$null_dist_dir))
  dir.create(args$null_dist_dir, recursive = TRUE, showWarnings = FALSE)

if (args$low_frac  <= 0 || args$low_frac  >= 1)
  stop("--low-frac must be between 0 and 1 (exclusive)")
if (args$high_frac <= 0 || args$high_frac >= 1)
  stop("--high-frac must be between 0 and 1 (exclusive)")

# ----- Load expression matrix -----
cat("Loading expression matrix...\n")
expr_dt  <- fread(args$expr_file, data.table = TRUE)
gene_col <- colnames(expr_dt)[1]
cat("  Gene ID column: '", gene_col, "'\n", sep = "")
cat("  Dimensions:", nrow(expr_dt), "genes x", ncol(expr_dt) - 1, "samples\n\n")

sample_cols <- setdiff(colnames(expr_dt), gene_col)
n_samples   <- length(sample_cols)
k_low       <- floor(n_samples * args$low_frac)
k_high      <- floor(n_samples * args$high_frac)

if (k_low  < 2) stop("LOW group has < 2 samples; increase --low-frac.")
if (k_high < 2) stop("HIGH group has < 2 samples; increase --high-frac.")

cat("  Total samples:", n_samples, "\n")
cat("  LOW group:  k =", k_low,  "samples\n")
cat("  HIGH group: k =", k_high, "samples\n\n")

gene_ids <- expr_dt[[gene_col]]
expr_mat <- as.matrix(expr_dt[, ..sample_cols])
rownames(expr_mat) <- gene_ids

# ----- Load mapping table -----
cat("Loading mapping table...\n")
mapping_dt <- fread(args$mapping_file, data.table = TRUE)
if (!all(c("gene_id", "SYMBOL") %in% colnames(mapping_dt)))
  stop("Mapping file must have 'gene_id' and 'SYMBOL' columns.")
mapping_dt <- mapping_dt[, .(gene_id, SYMBOL)]
cat("  Mapping entries:", nrow(mapping_dt), "\n\n")

# ----- Determine focus genes -----
if (!is.null(args$focus_genes)) {
  focus_gene_list <- trimws(strsplit(args$focus_genes, ",")[[1]])
  focus_gene_list <- intersect(focus_gene_list, gene_ids)
  if (length(focus_gene_list) == 0)
    stop("None of the --focus-genes IDs found in expression matrix.")
  cat("Restricting permutation to", length(focus_gene_list),
      "specified focus genes.\n\n")
} else {
  focus_gene_list <- gene_ids
}

n_total <- length(focus_gene_list)
cat("Running permutation MAD analysis for", n_total, "genes",
    "(", args$n_perms, "permutations each)...\n")

# Set global seed once so shuffles are reproducible across genes
set.seed(args$seed)

# ----- Main loop -----
results_list <- vector("list", n_total)
pvals_rows   <- vector("list", n_total)   # for permutation_pvals.tsv

for (i in seq_along(focus_gene_list)) {
  if (i %% 50 == 0 || i == n_total)
    cat("  Progress:", i, "/", n_total, "\n")

  res <- run_permutation_test(
    focus_gene = focus_gene_list[i],
    expr_mat   = expr_mat,
    gene_ids   = gene_ids,
    k_low      = k_low,
    k_high     = k_high,
    n_perms    = args$n_perms
  )

  # Optionally write per-gene null distribution TSV
  if (!is.null(args$null_dist_dir) && !is.null(res$null_deltas)) {
    gene_tag <- gsub("[^A-Za-z0-9_]", "_", res$focus_gene)
    tsv_path <- file.path(
      args$null_dist_dir,
      sprintf("%04d_%s_null_dist.tsv.gz", i, gene_tag)
    )
    null_dt <- data.table(
      gene_id        = res$focus_gene,
      perm_idx       = seq_len(args$n_perms),
      delta_mean_mad = res$null_deltas
    )
    fwrite(null_dt, tsv_path, compress = "gzip")

    # Row for the pvals TSV expected by plot_permutation_null_dist.R
    pvals_rows[[i]] <- data.table(
      gene_id             = res$focus_gene,
      obs_delta_mean_mad  = res$obs_delta_mean_mad,
      pval_delta_mean_mad = res$perm_p
    )
  }

  # Drop null_deltas before accumulating to avoid memory bloat
  res$null_deltas <- NULL
  results_list[[i]] <- as.data.table(res)
}

# Write pvals TSV companion (for plot_permutation_null_dist.R)
if (!is.null(args$null_dist_dir) && any(!sapply(pvals_rows, is.null))) {
  fwrite(
    rbindlist(pvals_rows, fill = TRUE),
    file.path(args$null_dist_dir, "permutation_pvals.tsv"),
    sep = "\t"
  )
  cat("Null distribution TSVs saved to:", args$null_dist_dir, "\n\n")
}

# ----- Assemble results table -----
results_dt <- rbindlist(results_list, fill = TRUE)

# BH-FDR correction
results_dt[, obs_p_adj  := p.adjust(obs_wilcoxon_p, method = "BH")]
results_dt[, perm_p_adj := p.adjust(perm_p,          method = "BH")]

# Significance stars
results_dt[, obs_stars      := get_sig_stars(obs_wilcoxon_p)]
results_dt[, obs_adj_stars  := get_sig_stars(obs_p_adj)]
results_dt[, perm_stars     := get_sig_stars(perm_p)]
results_dt[, perm_adj_stars := get_sig_stars(perm_p_adj)]

# Validated: significant by both parametric (BH) and permutation (BH) tests
results_dt[, validated := (!is.na(obs_p_adj)  & obs_p_adj  < 0.05) &
                          (!is.na(perm_p_adj) & perm_p_adj < 0.05)]

# Join gene symbols
results_dt <- merge(
  mapping_dt,
  results_dt,
  by.x  = "gene_id",
  by.y  = "focus_gene",
  all.y = TRUE
)

results_dt[, condition_label := args$condition_label]
setorder(results_dt, perm_p_adj, na.last = TRUE)

col_order <- c(
  "gene_id", "SYMBOL",
  "obs_delta_mean_mad", "obs_mean_mad_low", "obs_mean_mad_high",
  "obs_wilcoxon_p", "obs_stars", "obs_p_adj", "obs_adj_stars",
  "perm_p", "perm_stars", "perm_p_adj", "perm_adj_stars",
  "n_perms", "validated", "condition_label"
)
col_order  <- intersect(col_order, colnames(results_dt))
results_dt <- results_dt[, ..col_order]

# ----- Write xlsx -----
wb <- createWorkbook()
addWorksheet(wb, "MAD_permutation")
writeData(wb, "MAD_permutation", as.data.frame(results_dt),
          headerStyle = createStyle(textDecoration = "bold"))
setColWidths(wb, "MAD_permutation", cols = seq_len(ncol(results_dt)), widths = "auto")
saveWorkbook(wb, args$output_file, overwrite = TRUE)

# ----- Console summary -----
n_tested    <- nrow(results_dt)
n_obs_sig   <- sum(!is.na(results_dt$obs_p_adj)  & results_dt$obs_p_adj  < 0.05, na.rm = TRUE)
n_perm_sig  <- sum(!is.na(results_dt$perm_p_adj) & results_dt$perm_p_adj < 0.05, na.rm = TRUE)
n_validated <- sum(results_dt$validated, na.rm = TRUE)
pct_obs     <- round(100 * n_obs_sig   / n_tested, 1)
pct_perm    <- round(100 * n_perm_sig  / n_tested, 1)
pct_val     <- round(100 * n_validated / n_tested, 1)

cat("\n=== Summary ===\n")
cat("Total genes tested:                   ", n_tested,   "\n")
cat("Significant — parametric BH < 0.05:   ",
    n_obs_sig,  " (", pct_obs,  "%)\n", sep = "")
cat("Significant — permutation BH < 0.05: ",
    n_perm_sig, " (", pct_perm, "%)\n", sep = "")
cat("Validated (both):                     ",
    n_validated, " (", pct_val, "%)\n", sep = "")
cat("Results saved to:", args$output_file, "\n")
