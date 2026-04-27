#!/usr/bin/env Rscript
# ==============================================================================
# Batch Gene-Level MAD Variability: ALL Genes → xlsx Summary
#
# For every gene in the expression matrix, samples are split into LOW (bottom
# k%) and HIGH (top k%) groups by that gene's expression rank. MAD (Median
# Absolute Deviation) is computed for all other genes in each group and compared
# via a Wilcoxon rank-sum test. BH-FDR correction is applied across all genes.
# Results are written to a single xlsx file (one row per focus gene).
#
# Usage:
#   Rscript summarize_all_genes_mad_variability.R \
#     --expr-file    data/processed/VOOM/voomdataCtrl.txt \
#     --mapping-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
#     --output-file  results/variability/all_genes_mad_summary.xlsx \
#     --condition-label "Control"
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
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns"))))
}

# ----- Core: compute MAD stats for one focus gene -----
compute_mad_stats <- function(focus_gene, expr_mat, gene_ids,
                              sample_cols, k_low, k_high) {
  tryCatch({
    focus_idx <- which(gene_ids == focus_gene)
    focus_vec <- expr_mat[focus_idx, ]

    sorted_asc  <- order(focus_vec)
    sorted_desc <- order(focus_vec, decreasing = TRUE)

    low_idx  <- sorted_asc[seq_len(k_low)]
    high_idx <- sorted_desc[seq_len(k_high)]

    other_idx <- which(gene_ids != focus_gene)
    expr_other <- expr_mat[other_idx, , drop = FALSE]

    mad_low  <- apply(expr_other[, low_idx,  drop = FALSE], 1, mad, na.rm = TRUE)
    mad_high <- apply(expr_other[, high_idx, drop = FALSE], 1, mad, na.rm = TRUE)

    wt <- wilcox.test(mad_low, mad_high, paired = FALSE, exact = FALSE)

    list(
      focus_gene       = focus_gene,
      wilcoxon_p       = wt$p.value,
      mean_mad_low     = mean(mad_low,    na.rm = TRUE),
      mean_mad_high    = mean(mad_high,   na.rm = TRUE),
      delta_mean_mad   = mean(mad_high,   na.rm = TRUE) - mean(mad_low, na.rm = TRUE),
      median_mad_low   = median(mad_low,  na.rm = TRUE),
      median_mad_high  = median(mad_high, na.rm = TRUE),
      n_genes_compared = length(mad_low),
      n_low_samples    = k_low,
      n_high_samples   = k_high
    )
  }, error = function(e) {
    list(
      focus_gene       = focus_gene,
      wilcoxon_p       = NA_real_,
      mean_mad_low     = NA_real_,
      mean_mad_high    = NA_real_,
      delta_mean_mad   = NA_real_,
      median_mad_low   = NA_real_,
      median_mad_high  = NA_real_,
      n_genes_compared = NA_integer_,
      n_low_samples    = k_low,
      n_high_samples   = k_high
    )
  })
}

# ----- CLI arguments -----
parser <- ArgumentParser(description = "Batch gene-level MAD variability summary")
parser$add_argument("--expr-file",       required = TRUE,
                    help = "VOOM expression matrix (tab-separated, genes x samples)")
parser$add_argument("--mapping-file",    required = TRUE,
                    help = "TSV/CSV with gene_id and SYMBOL columns")
parser$add_argument("--output-file",     required = TRUE,
                    help = "Output .xlsx file path")
parser$add_argument("--low-frac",        type = "double", default = 0.2,
                    help = "Fraction of samples in LOW group (default: 0.2)")
parser$add_argument("--high-frac",       type = "double", default = 0.2,
                    help = "Fraction of samples in HIGH group (default: 0.2)")
parser$add_argument("--condition-label", default = "Condition",
                    help = "Human-readable condition label added as a column")
args <- parser$parse_args()

cat("=== Batch Gene-Level MAD Variability Analysis ===\n")
cat("Expression file:  ", args$expr_file,       "\n")
cat("Mapping file:     ", args$mapping_file,     "\n")
cat("Output file:      ", args$output_file,      "\n")
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

# Convert to plain matrix for fast row/column indexing
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

# ----- Loop over all genes -----
all_genes <- gene_ids
n_total   <- length(all_genes)
cat("Running MAD analysis for", n_total, "genes...\n")

results <- vector("list", n_total)
for (i in seq_along(all_genes)) {
  if (i %% 100 == 0 || i == n_total)
    cat("  Progress:", i, "/", n_total, "\n")
  results[[i]] <- compute_mad_stats(
    focus_gene  = all_genes[i],
    expr_mat    = expr_mat,
    gene_ids    = gene_ids,
    sample_cols = sample_cols,
    k_low       = k_low,
    k_high      = k_high
  )
}

results_dt <- rbindlist(results, fill = TRUE)

# ----- BH-FDR correction -----
results_dt[, wilcoxon_p_adj := p.adjust(wilcoxon_p, method = "BH")]

# Stars for raw and adjusted p
results_dt[, wilcoxon_stars     := get_sig_stars(wilcoxon_p)]
results_dt[, wilcoxon_adj_stars := get_sig_stars(wilcoxon_p_adj)]

# ----- Join gene symbols (mapping_dt on left → gene_id + SYMBOL first) -----
results_dt <- merge(
  mapping_dt,
  results_dt,
  by.x  = "gene_id",
  by.y  = "focus_gene",
  all.y = TRUE   # keep genes not in mapping (SYMBOL → NA)
)

# Add condition label and sort by adjusted p
results_dt[, condition_label := args$condition_label]
setorder(results_dt, wilcoxon_p_adj, na.last = TRUE)

# Reorder columns for readability
col_order <- c(
  "gene_id", "SYMBOL",
  "wilcoxon_p", "wilcoxon_stars",
  "wilcoxon_p_adj", "wilcoxon_adj_stars",
  "mean_mad_low", "mean_mad_high", "delta_mean_mad",
  "median_mad_low", "median_mad_high",
  "n_genes_compared", "n_low_samples", "n_high_samples",
  "condition_label"
)
col_order <- intersect(col_order, colnames(results_dt))
results_dt <- results_dt[, ..col_order]

# ----- Write xlsx -----
wb <- createWorkbook()
addWorksheet(wb, "MAD_variability")
writeData(wb, "MAD_variability", as.data.frame(results_dt),
          headerStyle = createStyle(textDecoration = "bold"))
setColWidths(wb, "MAD_variability", cols = seq_len(ncol(results_dt)), widths = "auto")
saveWorkbook(wb, args$output_file, overwrite = TRUE)

# ----- Console summary -----
n_tested <- nrow(results_dt)
sig_dt   <- results_dt[!is.na(wilcoxon_p_adj) & wilcoxon_p_adj < 0.05]
n_sig    <- nrow(sig_dt)
pct_sig  <- round(100 * n_sig / n_tested, 1)

n_increased <- sum(sig_dt$delta_mean_mad > 0, na.rm = TRUE)
n_decreased <- sum(sig_dt$delta_mean_mad < 0, na.rm = TRUE)
pct_inc  <- if (n_sig > 0) round(100 * n_increased / n_sig, 1) else 0
pct_dec  <- if (n_sig > 0) round(100 * n_decreased / n_sig, 1) else 0

cat("\n=== Summary ===\n")
cat("Total genes tested:            ", n_tested, "\n")
cat("Significant (FDR < 0.05):      ", n_sig, " (", pct_sig, "%)\n", sep = "")
cat("  \u2191 Increased variability (HIGH > LOW):  ",
    n_increased, " (", pct_inc, "% of sig)\n", sep = "")
cat("  \u2193 Decreased variability (LOW > HIGH):  ",
    n_decreased, " (", pct_dec, "% of sig)\n", sep = "")
cat("Results saved to: ", args$output_file, "\n")
