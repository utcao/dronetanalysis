#!/usr/bin/env Rscript
# ==============================================================================
# Bootstrap Identification of Stable High-MAD Genes
#
# Tests whether certain genes have constitutively high MAD (median absolute
# deviation) regardless of which individuals are included.  For each bootstrap
# iteration, a random subset of samples (--sample-frac) is drawn and per-gene
# MAD is computed and ranked.  Genes that consistently land in the top MAD
# percentile across iterations are "structurally variable" — not driven by any
# particular expression gradient.
#
# WHY THIS MATTERS:
#   The original MAD analysis splits samples by a focus gene's expression rank.
#   If a core set of genes have inherently high MAD, any random LOW/HIGH split
#   (regardless of which gene drives it) will show MAD differences for those
#   genes.  This inflates the Wilcoxon test across all focus genes.
#
#   Genes flagged here (is_stable_high_mad = TRUE) can be:
#     1. Excluded from the comparison pool when re-running
#        compute_mad_transcriptomic_variability.R → compare n_sig before/after
#     2. Used in residual-MAD analysis (subtract bootstrap-expected MAD before
#        Wilcoxon) to correct for baseline variability
#     3. Cross-tabulated with 'validated' from compute_mad_permute_transcriptomic_noise.R:
#        truly focus-gene-driven effects should be validated AND non-stable
#
# KEY OUTPUT COLUMNS:
#   mean_mad_rank_norm  : mean normalized MAD rank across iterations, 0-1
#                         (1 = always highest MAD, 0 = always lowest)
#   pct_in_top<X>       : fraction of iterations where gene is in top X% MAD
#                         (X = 100 * --top-pct; default 10)
#   pct_in_top20        : same for the top 20% fixed reference
#   is_stable_high_mad  : TRUE when pct_in_top<X> >= --stability-threshold
#
# Usage:
#   Rscript identify_stable_high_mad_genes.R \
#     --expr-file    data/processed/VOOM/voomdataCtrl.txt \
#     --mapping-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
#     --output-file  results/variability/stable_high_mad_ct.xlsx \
#     --n-bootstrap  500 \
#     --sample-frac  0.4 \
#     --seed         42
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
  library(openxlsx)
})

# ----- Row-wise MAD: use matrixStats if available, else apply fallback -----
row_mads <- if (requireNamespace("matrixStats", quietly = TRUE)) {
  function(mat) matrixStats::rowMads(mat, na.rm = TRUE)
} else {
  function(mat) apply(mat, 1, mad, na.rm = TRUE)
}

# ----- CLI arguments -----
parser <- ArgumentParser(
  description = "Bootstrap identification of stable high-MAD genes"
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
parser$add_argument("--n-bootstrap",
  type = "integer", default = 500,
  help = "Number of bootstrap iterations (default: 500)"
)
parser$add_argument("--sample-frac",
  type = "double", default = 0.4,
  help = "Fraction of samples drawn per bootstrap iteration (default: 0.4)"
)
parser$add_argument("--top-pct",
  type = "double", default = 0.10,
  help = "Top MAD percentile that defines 'high MAD' per iteration (default: 0.10)"
)
parser$add_argument("--stability-threshold",
  type = "double", default = 0.90,
  help = paste(
    "Minimum fraction of iterations a gene must be in the top --top-pct",
    "to be flagged is_stable_high_mad (default: 0.90)"
  )
)
parser$add_argument("--seed",
  type = "integer", default = 42,
  help = "Random seed for reproducibility (default: 42)"
)
parser$add_argument("--condition-label",
  default = "Condition",
  help = "Condition label added as a column"
)
args <- parser$parse_args()

top_pct_label <- sprintf("pct_in_top%02d", round(100 * args$top_pct))

cat("=== Bootstrap Identification of Stable High-MAD Genes ===\n")
cat("Expression file:      ", args$expr_file,          "\n")
cat("Mapping file:         ", args$mapping_file,        "\n")
cat("Output file:          ", args$output_file,         "\n")
cat("Bootstrap iterations: ", args$n_bootstrap,         "\n")
cat("Sample fraction:      ", args$sample_frac,         "\n")
cat("Top MAD percentile:   ", args$top_pct,
    sprintf("(top %s column)\n", top_pct_label))
cat("Stability threshold:  ", args$stability_threshold, "\n")
cat("Seed:                 ", args$seed,                "\n\n")

# ----- Validate inputs -----
if (!file.exists(args$expr_file))
  stop("Expression file not found: ", args$expr_file)
if (!file.exists(args$mapping_file))
  stop("Mapping file not found: ", args$mapping_file)

out_dir <- dirname(args$output_file)
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

if (args$sample_frac       <= 0 || args$sample_frac       >= 1)
  stop("--sample-frac must be between 0 and 1 (exclusive)")
if (args$top_pct           <= 0 || args$top_pct           >= 1)
  stop("--top-pct must be between 0 and 1 (exclusive)")
if (args$stability_threshold <= 0 || args$stability_threshold > 1)
  stop("--stability-threshold must be in (0, 1]")

# ----- Load expression matrix -----
cat("Loading expression matrix...\n")
expr_dt  <- fread(args$expr_file, data.table = TRUE)
gene_col <- colnames(expr_dt)[1]
cat("  Gene ID column: '", gene_col, "'\n", sep = "")
cat("  Dimensions:", nrow(expr_dt), "genes x", ncol(expr_dt) - 1, "samples\n\n")

sample_cols <- setdiff(colnames(expr_dt), gene_col)
n_samples   <- length(sample_cols)
k_draw      <- floor(n_samples * args$sample_frac)

if (k_draw < 3)
  stop("Bootstrap sample size (", k_draw, ") < 3; increase --sample-frac.")

gene_ids <- expr_dt[[gene_col]]
n_genes  <- length(gene_ids)
expr_mat <- as.matrix(expr_dt[, ..sample_cols])
rownames(expr_mat) <- gene_ids

cat("  Total samples:              ", n_samples, "\n")
cat("  Samples per bootstrap draw: ", k_draw,    "\n\n")

# ----- Load mapping table -----
cat("Loading mapping table...\n")
mapping_dt <- fread(args$mapping_file, data.table = TRUE)
if (!all(c("gene_id", "SYMBOL") %in% colnames(mapping_dt)))
  stop("Mapping file must have 'gene_id' and 'SYMBOL' columns.")
mapping_dt <- mapping_dt[, .(gene_id, SYMBOL)]
cat("  Mapping entries:", nrow(mapping_dt), "\n\n")

# ----- Bootstrap loop -----
# For each iteration:
#   1. Draw k_draw samples without replacement
#   2. Compute per-gene MAD
#   3. Rank genes by MAD (rank 1 = highest MAD, using negation)
#   4. Record normalized rank and whether gene is in the top percentile

n_top_x  <- floor(n_genes * args$top_pct)        # count threshold for top_pct
n_top_20 <- floor(n_genes * 0.20)                 # fixed top-20% count

# Accumulators
rank_sum    <- numeric(n_genes)   # sum of normalized ranks across iterations
rank_sq_sum <- numeric(n_genes)   # sum of squared normalized ranks (for SD)
top_x_count <- integer(n_genes)   # how often gene lands in top top_pct
top_20_count <- integer(n_genes)  # how often gene lands in top 20%

cat("Running", args$n_bootstrap, "bootstrap iterations...\n")
set.seed(args$seed)

for (b in seq_len(args$n_bootstrap)) {
  if (b %% 100 == 0 || b == args$n_bootstrap)
    cat("  Iteration:", b, "/", args$n_bootstrap, "\n")

  samp_idx  <- sample(n_samples, k_draw, replace = FALSE)
  sub_mat   <- expr_mat[, samp_idx, drop = FALSE]
  gene_mads <- row_mads(sub_mat)

  # rank(-mads): gene with highest MAD gets rank 1
  gene_ranks <- rank(-gene_mads, ties.method = "average")
  norm_ranks <- gene_ranks / n_genes     # normalize to (0, 1], lower = higher MAD

  rank_sum    <- rank_sum    + norm_ranks
  rank_sq_sum <- rank_sq_sum + norm_ranks^2
  top_x_count <- top_x_count  + as.integer(gene_ranks <= n_top_x)
  top_20_count <- top_20_count + as.integer(gene_ranks <= n_top_20)
}

# ----- Summarize across iterations -----
cat("\nSummarizing bootstrap results...\n")

mean_rank  <- rank_sum / args$n_bootstrap
# Variance: E[X^2] - E[X]^2
var_rank   <- pmax(0, rank_sq_sum / args$n_bootstrap - mean_rank^2)
sd_rank    <- sqrt(var_rank)

# Invert rank so that mean_mad_rank_norm = 1 means always the highest MAD gene
mean_mad_rank_norm <- 1 - mean_rank
sd_mad_rank_norm   <- sd_rank                 # SD is invariant to inversion

pct_in_top_x  <- top_x_count  / args$n_bootstrap
pct_in_top20  <- top_20_count / args$n_bootstrap

is_stable <- pct_in_top_x >= args$stability_threshold

# ----- Build results data.table -----
results_dt <- data.table(
  gene_id            = gene_ids,
  mean_mad_rank_norm = round(mean_mad_rank_norm, 4),
  sd_mad_rank_norm   = round(sd_mad_rank_norm,   4),
  pct_in_top20       = round(pct_in_top20,       4),
  is_stable_high_mad = is_stable
)
# Dynamic column name for the user-specified top percentile
results_dt[, (top_pct_label) := round(pct_in_top_x, 4)]

# Join gene symbols
results_dt <- merge(
  mapping_dt,
  results_dt,
  by    = "gene_id",
  all.y = TRUE
)

results_dt[, condition_label := args$condition_label]
setorder(results_dt, -mean_mad_rank_norm, na.last = TRUE)

col_order <- c(
  "gene_id", "SYMBOL",
  "mean_mad_rank_norm", "sd_mad_rank_norm",
  top_pct_label, "pct_in_top20",
  "is_stable_high_mad", "condition_label"
)
col_order  <- intersect(col_order, colnames(results_dt))
results_dt <- results_dt[, ..col_order]

# ----- Write xlsx -----
wb <- createWorkbook()
addWorksheet(wb, "stable_high_MAD")
writeData(wb, "stable_high_MAD", as.data.frame(results_dt),
          headerStyle = createStyle(textDecoration = "bold"))
setColWidths(wb, "stable_high_MAD", cols = seq_len(ncol(results_dt)), widths = "auto")
saveWorkbook(wb, args$output_file, overwrite = TRUE)

# ----- Console summary -----
n_flagged   <- sum(is_stable, na.rm = TRUE)
pct_flagged <- round(100 * n_flagged / n_genes, 1)

cat("\n=== Summary ===\n")
cat("Total genes evaluated:        ", n_genes,     "\n")
cat("Stable high-MAD genes flagged:", n_flagged,
    " (", pct_flagged, "%)\n", sep = "")
cat(sprintf(
  "  (in top %d%% MAD in >= %d%% of %d iterations)\n",
  round(100 * args$top_pct),
  round(100 * args$stability_threshold),
  args$n_bootstrap
))
cat("Results saved to:", args$output_file, "\n")
cat("\nNEXT STEPS:\n")
cat("  1. Exclude flagged genes from the 'other genes' comparison pool in\n")
cat("     compute_mad_transcriptomic_variability.R and compare n_sig before/after.\n")
cat("  2. Cross-tabulate is_stable_high_mad with the 'validated' column from\n")
cat("     compute_mad_permute_transcriptomic_noise.R.\n")
cat("     Genes that are validated AND NOT stable are the strongest candidates.\n")
