#!/usr/bin/env Rscript
# ==============================================================================
# Sample Quintile Overlap: Expression-Gradient Partitioning
#
# For every gene in the expression matrix, samples are ranked by expression and
# assigned to quintiles Q1 (bottom 20%, lowest expression) through Q5 (top 20%,
# highest expression).  The script then tallies how often each sample lands in
# each quintile across all genes and visualises the result.
#
# LOGIC
#   For each gene: rank(sample_1, ..., sample_N) -> quintile label Q1-Q5.
#   Q1 = lowest-expressing 20% of samples for that gene.
#   Count accumulated over all genes -> per-sample Q1-Q5 membership count.
#   Frequency = count / n_genes  (expected ~0.20 under uniform distribution).
#
# UNIFORMITY TEST (per sample)
#   A chi-square goodness-of-fit test (df = 4, expected = n_genes / 5 per bin)
#   tests whether a sample's Q1-Q5 distribution deviates from uniform.
#   Effect-size measures:
#     cramer_v      -- Cramér's V = sqrt(chi2 / (n_genes * 4)); bounded [0,1];
#                      0 = perfectly uniform, 1 = all mass in one quintile.
#                      Comparable across datasets of different sizes.
#     entropy_norm  -- normalised Shannon entropy; 1.0 = perfectly uniform,
#                      0.0 = entirely concentrated in one quintile
#     max_deviation -- max |freq_qi - 0.20| across the five quintiles
#
# OUTPUTS  (<output-dir>/)
#   sample_quintile_counts.tsv     -- one row per sample: Q1-Q5 raw counts,
#                                     frequencies, chi-sq stat/p/padj,
#                                     cramer_v, entropy_norm, max_deviation.
#                                     Sorted by Q1_freq desc.
#   sample_noise_metrics.tsv       -- one row per sample: MAD, SD, IQR, CV,
#                                     CV2 computed across all genes.
#                                     Sorted by mad_log2 desc.
#   quintile_matrix.rds            -- full integer matrix (n_genes x n_samples),
#                                     values 1-5; rownames = gene IDs,
#                                     colnames = sample IDs.  Covers all genes
#                                     and all quintiles; use for downstream
#                                     driver-gene analysis.
#   gene_q1q5_wide.tsv.gz          -- human-readable companion: one row per
#                                     gene (top 5000), Q1_samples / Q5_samples
#                                     as pipe-separated strings; includes
#                                     gene_max_dev column.
#   quintile_distribution.pdf      -- page 1: per-quintile density of sample
#                                     frequencies; page 2: violin + boxplot
#   quintile_heatmap.pdf           -- clustered heatmap: samples (rows,
#                                     unlabelled) x Q1-Q5 (cols)
#
# GENE FILTERING (quintile_matrix.rds / gene_q1q5_wide.tsv.gz)
#   Every gene always has ~n_samples/5 members per quintile, so its own marginal
#   distribution is trivially uniform and uninformative.  Instead, genes are
#   scored by how "globally extreme" their Q1/Q5 groups are:
#     gene_max_dev = max(mean(Q1_freq of Q1 samples), mean(Q5_freq of Q5 samples)) - 0.20
#   A gene whose Q1 group is populated by samples that are globally Q1-dominant
#   has a high gene_max_dev -- it reinforces the transcriptome-wide tier pattern.
#   Only the top 5000 genes by gene_max_dev are stored.
#
# Usage:
#   Rscript plot_sample_quintile_overlap.R \
#     --expr-file    data/processed/VOOM/voomdataCtrl.txt \
#     --output-dir   results/quintile_overlap_ct \
#     --condition-label "Control"
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
  library(pheatmap)
})

# ---- CLI --------------------------------------------------------------------
parser <- ArgumentParser(description = "Sample quintile overlap analysis")
parser$add_argument("--expr-file",       required = TRUE,
                    help = "VOOM/VST expression matrix (tab-sep, genes x samples)")
parser$add_argument("--output-dir",      required = TRUE,
                    help = "Directory for all output files (created if absent)")
parser$add_argument("--condition-label", default = "Condition",
                    help = "Label used in plot titles (default: 'Condition')")
parser$add_argument("--select-gene", default = NA,
                    help = "file path to the gene list for restricting analysis (default: NA - no filter)")
parser$add_argument("--seed",            type = "integer", default = 42,
                    help = "Random seed for tie-breaking in ranks (default: 42)")
args <- parser$parse_args()

cat("=== Sample Quintile Overlap Analysis ===\n")
cat("Expression file:  ", args$expr_file,       "\n")
cat("Output directory: ", args$output_dir,       "\n")
cat("Condition label:  ", args$condition_label,  "\n\n")

if (!file.exists(args$expr_file))
  stop("Expression file not found: ", args$expr_file)

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(args$seed)

# ---- Source shared helpers ---------------------------------------------------
# find the directory of this script based on where the script was invoked from (Rscript --file=...) 
script_dir <- tryCatch(
  dirname(normalizePath(
    sub("--file=", "",
        grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]),
    mustWork = FALSE
  )),
  error = function(e) getwd()
)
source(file.path(script_dir, "utils_sample_quintile.R"))

n_steps <- 100L   # colour steps per side for heatmap diverging palette

# ---- Load expression matrix -------------------------------------------------
cat("Loading expression matrix...\n")
expr_dt  <- fread(args$expr_file, data.table = TRUE)
gene_col <- colnames(expr_dt)[1]
cat("  Gene column:  '", gene_col, "'\n", sep = "")

sample_cols <- setdiff(colnames(expr_dt), gene_col)
n_samples   <- length(sample_cols)
gene_ids    <- expr_dt[[gene_col]]
n_genes     <- length(gene_ids)

cat("  Genes:    ", n_genes,   "\n")
cat("  Samples:  ", n_samples, "\n\n")

expr_mat <- as.matrix(expr_dt[, ..sample_cols])
rownames(expr_mat) <- gene_ids

# ---- Quintile assignment ----------------------------------------------------
# For each gene (row), rank all samples and split into Q1-Q5.
# Formula: quintile = ceiling( rank * 5 / n ), capped at 5.
# Ties broken randomly so each sample receives a unique quintile label.
cat("Assigning quintiles across", n_genes, "genes x", n_samples, "samples...\n")

assign_quintile <- function(x) {
  r <- rank(x, ties.method = "random", na.last = "keep")
  n <- sum(!is.na(x))
  # ceiling(r * K / n) is the standard rank-to-bin formula (same as dplyr::ntile).
  # For ranks 1..n and K bins it produces bin sizes of floor(n/K) or ceil(n/K),
  # differing by at most 1 — the smallest possible imbalance for any equal division.
  # Equivalent to cut(rank, quantile-breaks) but faster (no factor allocation).
  # pmin caps at 5 because ceiling(n * 5 / n) can equal 6 under floating-point
  # rounding (e.g. ceiling(1.0 + eps)); vectorised, no logical-vector overhead.
  pmin(as.integer(ceiling(r * 5L / n)), 5L)
}

quintile_mat <- t(apply(expr_mat, 1L, assign_quintile))
# dimensions: n_genes x n_samples  (integer 1-5)
rownames(quintile_mat) <- gene_ids
colnames(quintile_mat) <- sample_cols
cat("  Done.\n\n")

# ---- Per-sample Q1-Q5 counts and frequencies --------------------------------
cat("Counting Q1-Q5 memberships per sample...\n")

# vapply is preferred over lapply (returns a list needing do.call/unlist) and
# over sapply (return type varies with input — silent failure risk).  The
# FUN.VALUE = numeric(n_samples) contract pre-allocates the result matrix and
# raises an error immediately if a column has the wrong length or type.
count_mat <- vapply(
  1:5,
  function(q) colSums(quintile_mat == q, na.rm = TRUE),
  numeric(n_samples)
)
colnames(count_mat) <- paste0("Q", 1:5)
rownames(count_mat) <- sample_cols

freq_mat <- count_mat / n_genes   # expected ~0.20 per quintile

# ---- Uniformity tests -------------------------------------------------------
# Chi-square goodness-of-fit: observed Q1-Q5 counts vs. uniform expectation.
# Vectorised: chi_sq_i = sum_q (obs_iq - E)^2 / E  where E = n_genes / 5
expected   <- n_genes / 5
chisq_stat <- rowSums((count_mat - expected)^2 / expected)
chisq_p    <- pchisq(chisq_stat, df = 4L, lower.tail = FALSE)
chisq_padj <- p.adjust(chisq_p, method = "BH")

# Cramér's V for one-way GOF with k=5 categories: V = sqrt(chi2 / (N*(k-1))).
# Bounded [0,1]; 0 = perfectly uniform, 1 = all mass in one quintile.
cramer_v <- sqrt(chisq_stat / (n_genes * 4L))

# Normalised Shannon entropy: H / log2(5).  1 = uniform, 0 = all in one bin.
safe_freq    <- pmax(freq_mat, 1e-15)
entropy_norm <- -rowSums(safe_freq * log2(safe_freq)) / log2(5)

# Maximum absolute deviation from the expected frequency 0.20
max_dev <- apply(abs(freq_mat - 0.2), 1L, max)

# ---- Per-sample noise metrics across genes ----------------------------------
# For each sample, treat the n_genes expression values as a distribution and
# summarise its spread.  All primary stats on the log2 scale (VOOM units);
# CV2 additionally on the linear scale for compatibility with per-gene stats
# elsewhere in this project.
cat("Computing per-sample noise metrics across genes...\n")

lin_mat    <- 2^expr_mat

n_valid     <- colSums(!is.na(lin_mat))
mean_log2   <- log2(colMeans(lin_mat, na.rm = TRUE)+1)
median_log2 <- log2(apply(lin_mat, 2L, median, na.rm = TRUE)+1)
sd_log2     <- log2(apply(lin_mat, 2L, sd,     na.rm = TRUE)+1)
mad_log2    <- log2(apply(lin_mat, 2L, mad,    na.rm = TRUE)+1)
iqr_log2    <- log2(apply(lin_mat, 2L, IQR,   na.rm = TRUE)+1)
cv_log2     <- sd_log2 - mean_log2

noise_dt <- data.table(
  sample_id   = sample_cols,
  n_genes     = n_valid,
  mean_log2   = mean_log2,
  median_log2 = median_log2,
  sd_log2     = sd_log2,
  mad_log2    = mad_log2,
  iqr_log2    = iqr_log2,
  cv_log2     = cv_log2
)
setorder(noise_dt, -mad_log2)

# ---- Output 1: sample_quintile_counts.tsv -----------------------------------
sample_dt <- as.data.table(count_mat, keep.rownames = "sample_id")

freq_dt <- as.data.table(freq_mat, keep.rownames = "sample_id")
setnames(freq_dt, paste0("Q", 1:5), paste0("Q", 1:5, "_freq"))

stats_dt <- data.table(
  sample_id     = sample_cols,
  chisq_stat    = chisq_stat,
  chisq_p       = chisq_p,
  chisq_padj    = chisq_padj,
  cramer_v      = cramer_v,
  entropy_norm  = entropy_norm,
  max_deviation = max_dev
)

sample_dt <- Reduce(function(a, b) merge(a, b, by = "sample_id"),
                    list(sample_dt, freq_dt, stats_dt))
setorder(sample_dt, -Q1_freq)

out_counts <- file.path(args$output_dir, "sample_quintile_counts.tsv")
fwrite(sample_dt, out_counts, sep = "\t")
cat("  Counts ->", out_counts, "\n")

# ---- Output 2: sample_noise_metrics.tsv -------------------------------------
out_noise <- file.path(args$output_dir, "sample_noise_metrics.tsv")
fwrite(noise_dt, out_noise, sep = "\t")
cat("  Noise  ->", out_noise, "\n")

# ---- Gene-level max_deviation and top-5000 filter ---------------------------
# Every gene has exactly ~n_samples/5 members per quintile by construction, so
# the gene's own marginal distribution is always uniform and uninformative.
# Instead, measure how "globally extreme" a gene's Q1/Q5 groups are:
#   gene_q1_bias = mean(Q1_freq of Q1-assigned samples) - 0.20
#   gene_q5_bias = mean(Q5_freq of Q5-assigned samples) - 0.20
# A gene whose Q1 group is populated by samples that are globally Q1-dominant
# (high Q1_freq overall) scores high — it reinforces the transcriptome-wide
# expression tier pattern.  gene_max_dev = max(q1_bias, q5_bias).
# Vectorised via matrix products: no per-gene loop needed.
cat("Computing gene-level max_deviation and selecting top 5000 genes...\n")

q1_indicator  <- quintile_mat == 1L   # n_genes x n_samples logical
q5_indicator  <- quintile_mat == 5L
q1_freq_vec   <- freq_mat[, "Q1"]     # per-sample Q1 frequency
q5_freq_vec   <- freq_mat[, "Q5"]

gene_q1_bias  <- as.vector(q1_indicator %*% q1_freq_vec) / rowSums(q1_indicator) - 0.2
gene_q5_bias  <- as.vector(q5_indicator %*% q5_freq_vec) / rowSums(q5_indicator) - 0.2
gene_max_dev  <- pmax(gene_q1_bias, gene_q5_bias)

n_top         <- min(5000L, n_genes)
top_idx       <- order(gene_max_dev, decreasing = TRUE)[seq_len(n_top)]
top_gene_ids  <- gene_ids[top_idx]
cat("  Retaining", n_top, "genes (of", n_genes, ") with highest gene_max_dev\n")
cat("  gene_max_dev range (top genes): [",
    round(min(gene_max_dev[top_idx]), 4), ",",
    round(max(gene_max_dev[top_idx]), 4), "]\n\n")

# ---- Output 2a: quintile_matrix.rds (full integer matrix) -------------------
# Full quintile assignment matrix: rows = all genes, cols = all samples,
# values 1-5.  Load with: qmat <- readRDS("quintile_matrix.rds")
cat("Saving full quintile matrix (RDS, all", n_genes, "genes)...\n")

out_qmat <- file.path(args$output_dir, "quintile_matrix.rds")
saveRDS(quintile_mat, out_qmat, compress = "xz")
cat("  Quintile matrix ->", out_qmat, "\n")

# ---- Output 2b: gene_q1q5_wide.tsv.gz (human-readable companion) ------------
cat("Building gene Q1/Q5 wide table (top", n_top, "genes)...\n")

top_qmat <- quintile_mat[top_idx, , drop = FALSE]
gene_q1q5_wide <- data.table(
  gene_id      = top_gene_ids,
  gene_max_dev = gene_max_dev[top_idx],
  Q1_n_samples = rowSums(top_qmat == 1L, na.rm = TRUE),
  Q1_samples   = apply(top_qmat, 1L, function(r)
                   paste(sample_cols[!is.na(r) & r == 1L], collapse = "|")),
  Q5_n_samples = rowSums(top_qmat == 5L, na.rm = TRUE),
  Q5_samples   = apply(top_qmat, 1L, function(r)
                   paste(sample_cols[!is.na(r) & r == 5L], collapse = "|"))
)

out_wide <- file.path(args$output_dir, "gene_q1q5_wide.tsv.gz")
fwrite(gene_q1q5_wide, out_wide, sep = "\t", compress = "gzip")
cat("  Gene Q1/Q5 wide ->", out_wide, "\n\n")

# ---- Plot 1: Distribution of per-sample frequencies per quintile ------------
cat("Plotting quintile frequency distributions...\n")

freq_long <- melt(
  as.data.table(freq_mat, keep.rownames = "sample_id"),
  id.vars       = "sample_id",
  variable.name = "quintile",
  value.name    = "frequency"
)
freq_long[, quintile := factor(quintile, levels = paste0("Q", 1:5))]

q_colours <- c(Q1 = "#4393C3", Q2 = "#92C5DE", Q3 = "#BABABA",
               Q4 = "#F4A582", Q5 = "#D6604D")

n_sig  <- sum(chisq_padj < 0.05)
p_dist <- ggplot(freq_long, aes(x = frequency, fill = quintile, colour = quintile)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(xintercept = 0.2, linetype = "dashed",
             colour = "grey30", linewidth = 0.65) +
  facet_wrap(~quintile, nrow = 1) +
  scale_fill_manual(values   = q_colours) +
  scale_colour_manual(values = q_colours) +
  labs(
    title    = sprintf("Sample Quintile Membership Frequencies  -  %s", args$condition_label),
    subtitle = sprintf(
      "n = %d samples, %d genes  |  dashed = expected 0.20  |  %d/%d samples non-uniform (chi-sq BH < 0.05)",
      n_samples, n_genes, n_sig, n_samples
    ),
    x = "Frequency  (proportion of genes)",
    y = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 11)
  )

p_box <- ggplot(freq_long, aes(x = quintile, y = frequency, fill = quintile)) +
  geom_violin(alpha = 0.35, linewidth = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.14, outlier.size = 0.7, outlier.alpha = 0.5,
               linewidth = 0.55, alpha = 0.85) +
  geom_hline(yintercept = 0.2, linetype = "dashed",
             colour = "grey30", linewidth = 0.65) +
  scale_fill_manual(values = q_colours) +
  labs(
    title    = sprintf("Sample Quintile Frequency Distribution  -  %s", args$condition_label),
    subtitle = sprintf(
      "n = %d samples, %d genes  |  dashed = expected 0.20  |  violin + box; dots = outliers",
      n_samples, n_genes
    ),
    x = "Quintile",
    y = "Frequency of assignment (proportion of genes)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold")
  )

dist_pdf <- file.path(args$output_dir, "quintile_distribution.pdf")
pdf(dist_pdf, width = 12, height = 5)
print(p_dist)
print(p_box)
invisible(dev.off())
cat("  Distribution ->", dist_pdf, "\n")

# ---- Plot 2: Clustered heatmap (samples x Q1-Q5) ----------------------------
# draw_quintile_freq_heatmap is defined in utils_sample_quintile.R
cat("Plotting clustered heatmap...\n")

hm_mat           <- freq_mat
rownames(hm_mat) <- sample_cols

hm_pdf <- file.path(args$output_dir, "quintile_heatmap.pdf")
draw_quintile_freq_heatmap(
  freq_mat   = hm_mat,
  cond_label = args$condition_label,
  n_steps    = n_steps,
  out_pdf    = hm_pdf
)
cat("\n")

# ---- Console summary --------------------------------------------------------
cat("=== Summary ===\n")
cat(sprintf("Genes:    %d\nSamples:  %d\n\n", n_genes, n_samples))

cat("Uniformity (chi-sq BH < 0.05):", n_sig, "/", n_samples, "samples non-uniform\n\n")

cat("Frequency statistics (expected ~0.200 per quintile):\n")
for (q in paste0("Q", 1:5)) {
  v <- freq_mat[, q]
  cat(sprintf("  %s  mean=%.4f  sd=%.5f  range=[%.4f, %.4f]\n",
              q, mean(v), sd(v), min(v), max(v)))
}

cat(sprintf(
  "\nEntropy (normalised): mean=%.4f  sd=%.5f  range=[%.4f, %.4f]\n",
  mean(entropy_norm), sd(entropy_norm), min(entropy_norm), max(entropy_norm)
))
cat(sprintf(
  "Cramér's V:           mean=%.4f  sd=%.5f  range=[%.4f, %.4f]\n",
  mean(cramer_v), sd(cramer_v), min(cramer_v), max(cramer_v)
))

eff_cor <- cor(cbind(cramer_v, entropy_norm, max_deviation = max_dev),
               method = "pearson")
cat("\nEffect-size metric correlations (Pearson):\n")
print(round(eff_cor, 4))

cat(sprintf(
  "\nNoise metrics (across %d genes):\n  MAD log2: mean=%.4f sd=%.5f range=[%.4f, %.4f]\n  CV  log2: mean=%.4f sd=%.5f range=[%.4f, %.4f]\n",
  n_genes,
  mean(mad_log2),  sd(mad_log2),  min(mad_log2),  max(mad_log2),
  mean(cv_log2),   sd(cv_log2),   min(cv_log2),   max(cv_log2)
))

cat("\nOutputs:\n",
    " ", out_counts, "\n",
    " ", out_noise,  "\n",
    " ", out_qmat,   "\n",
    " ", out_wide,   "\n",
    " ", dist_pdf,   "\n",
    " ", hm_pdf,     "\n")
