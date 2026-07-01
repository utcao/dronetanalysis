#!/usr/bin/env Rscript
# ==============================================================================
# Quintile imbalance vs. batch/technical covariates
#
# Correlates cramer_v, entropy_norm, and max_deviation from
# sample_quintile_counts.tsv against the batch and surrogate-variable covariates
# used in the GAMLSS quintile analysis.
#
# SV files (sv1–sv4) are positionally aligned to the raw count file column
# order (first line of count file = sample IDs).  CT and HS samples are mixed
# in the meta/SV files; conditioning on a condition-specific
# sample_quintile_counts.tsv (inner-join) separates them automatically.
#
# Outputs (<output-dir>/):
#   covariate_correlations_{label}.tsv
#   covariate_dotplot_{label}.pdf   — Spearman rho per covariate (Cleveland dot)
#   covariate_scatter_{label}.pdf   — scatter per covariate x metric
#
# Usage:
#   Rscript plot_umbalance_batch_effectors.R \
#     --sample-metrics  results/quintile_overlap_vst_ct/sample_quintile_counts.tsv \
#     --count-file      data/count/RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt \
#     --meta-file       data/count/Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt \
#     --sv1             data/count/TMM_Voom_sv1_10.txt \
#     --sv2             data/count/TMM_Voom_sv2_10.txt \
#     --sv3             data/count/TMM_Voom_sv3_10.txt \
#     --sv4             data/count/TMM_Voom_sv4_10.txt \
#     --output-dir      results/batch_imbalance_ct \
#     --condition-label CT
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

# ---- CLI ---------------------------------------------------------------------
parser <- ArgumentParser(description = "Correlate quintile imbalance with batch covariates")
parser$add_argument("--sample-metrics", required = TRUE,
                    help = "sample_quintile_counts.tsv (cramer_v, entropy_norm, max_deviation)")
parser$add_argument("--count-file", required = TRUE,
                    help = "Raw count file; first line = sample IDs (positional SV anchor)")
parser$add_argument("--meta-file", required = TRUE,
                    help = "Metadata TSV (columns: id, RNAlibBatch, RNAseqBatch, egglayBatch, platingBatch, well)")
parser$add_argument("--sv1", required = TRUE,
                    help = "SV1 file; rows positionally aligned to count-file columns")
parser$add_argument("--sv2", required = TRUE, help = "SV2 file")
parser$add_argument("--sv3", required = TRUE, help = "SV3 file")
parser$add_argument("--sv4", required = TRUE, help = "SV4 file")
parser$add_argument("--output-dir", required = TRUE,
                    help = "Output directory (created if absent)")
parser$add_argument("--condition-label", default = "CT",
                    help = "Label appended to output file names and plot titles (default: CT)")
parser$add_argument("--pval-cutoff", type = "double", default = 0.05,
                    help = "BH-adjusted p-value threshold for significance markers (default: 0.05)")
args <- parser$parse_args()

lbl <- args$condition_label
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Quintile Imbalance ~ Batch Covariates ===\n")
cat("Condition:", lbl, "\n")
cat("Output:   ", args$output_dir, "\n\n")

for (f in c(args$sample_metrics, args$count_file, args$meta_file,
            args$sv1, args$sv2, args$sv3, args$sv4))
  if (!file.exists(f)) stop("File not found: ", f)

# ---- Load imbalance metrics --------------------------------------------------
cat("Loading sample metrics:", args$sample_metrics, "\n")
metrics_dt     <- fread(args$sample_metrics, data.table = TRUE)
imbalance_cols <- c("cramer_v", "entropy_norm", "max_deviation")
missing_cols   <- setdiff(c("sample_id", imbalance_cols), colnames(metrics_dt))
if (length(missing_cols) > 0L)
  stop("Columns missing from --sample-metrics: ", paste(missing_cols, collapse = ", "))
cat("  Samples:", nrow(metrics_dt), "\n")

# ---- Read count file header (sample IDs in column order) ---------------------
# SV files have values positionally aligned to this column order — the SV files
# themselves have no sample ID column.
cat("Reading count file header:", args$count_file, "\n")
sample_order <- scan(args$count_file, what = "", nlines = 1L, quiet = TRUE)
cat("  Samples in count file:", length(sample_order), "\n")

# ---- Load SV files positionally ----------------------------------------------
load_sv_positional <- function(sv_path, sample_ids, sv_name) {
  sv_raw  <- read.table(sv_path, header = TRUE)
  sv_vals <- sv_raw[[1L]]
  n_exp   <- length(sample_ids)
  if (length(sv_vals) != n_exp)
    stop(sv_name, " has ", length(sv_vals), " rows; expected ", n_exp,
         " (count file columns)")
  out <- data.table(id = sample_ids, v = sv_vals)
  setnames(out, "v", sv_name)
  out
}

sv_dt <- Reduce(
  function(a, b) merge(a, b, by = "id"),
  list(
    load_sv_positional(args$sv1, sample_order, "sv1"),
    load_sv_positional(args$sv2, sample_order, "sv2"),
    load_sv_positional(args$sv3, sample_order, "sv3"),
    load_sv_positional(args$sv4, sample_order, "sv4")
  )
)

# ---- Load metadata and attach SVs --------------------------------------------
cat("Loading metadata:", args$meta_file, "\n")
meta_dt <- fread(args$meta_file, data.table = TRUE)
meta_dt <- merge(meta_dt, sv_dt, by = "id", all.x = FALSE)
setnames(meta_dt, "id", "sample_id")

# ---- Inner join: keeps only samples in the condition-specific metrics file ---
merged_dt <- merge(
  metrics_dt[, .SD, .SDcols = c("sample_id", imbalance_cols)],
  meta_dt,
  by = "sample_id"
)
cat("  Samples matched:", nrow(merged_dt), "\n\n")
if (nrow(merged_dt) == 0L)
  stop("No samples matched — verify sample IDs in metrics file vs meta file.")

# ---- Encode covariates as numeric for Spearman -------------------------------
# Batch columns are integer-coded in the meta file; coerce to be safe.
for (bc in c("RNAlibBatch", "RNAseqBatch", "egglayBatch", "platingBatch"))
  set(merged_dt, j = bc, value = as.integer(merged_dt[[bc]]))

# well (e.g. "B11"): decompose into plate row (A=1..H=8) and column (1..12)
# so each dimension has an interpretable ordinal scale for Spearman.
merged_dt[, well_row := match(substr(well, 1L, 1L), LETTERS[1:8])]
merged_dt[, well_col := as.integer(substr(well, 2L, nchar(well)))]

# Pre-computed covariate map: column name -> display label (constant across loops)
cov_map <- c(
  RNAlibBatch  = "RNAlibBatch",
  RNAseqBatch  = "RNAseqBatch",
  egglayBatch  = "egglayBatch",
  platingBatch = "platingBatch",
  well_row     = "well row (A->H)",
  well_col     = "well col (1->12)",
  sv1          = "sv1",
  sv2          = "sv2",
  sv3          = "sv3",
  sv4          = "sv4"
)
cov_cols   <- names(cov_map)
cov_labels <- unname(cov_map)

# ---- Spearman correlations ---------------------------------------------------
cat("Computing Spearman correlations...\n")
cor_rows <- vector("list", length(cov_cols) * length(imbalance_cols))
idx <- 0L
for (mc in imbalance_cols) {
  y_all <- merged_dt[[mc]]
  for (ci in seq_along(cov_cols)) {
    x_all <- merged_dt[[cov_cols[ci]]]
    keep  <- complete.cases(x_all, y_all)
    idx   <- idx + 1L
    if (sum(keep) < 5L) {
      cor_rows[[idx]] <- data.table(
        metric    = mc,
        covariate = cov_labels[ci],
        rho       = NA_real_,
        pval      = NA_real_,
        n         = sum(keep)
      )
      next
    }
    ct <- cor.test(x_all[keep], y_all[keep], method = "spearman", exact = FALSE)
    cor_rows[[idx]] <- data.table(
      metric    = mc,
      covariate = cov_labels[ci],
      rho       = unname(ct$estimate),
      pval      = ct$p.value,
      n         = sum(keep)
    )
  }
}
cor_dt <- rbindlist(cor_rows)
cor_dt[, padj := p.adjust(pval, method = "BH"), by = "metric"]
cor_dt[, significant := !is.na(padj) & padj < args$pval_cutoff]

cat(sprintf("  Significant (padj < %.2g): %d / %d covariate-metric pairs\n",
            args$pval_cutoff,
            sum(cor_dt$significant, na.rm = TRUE),
            nrow(cor_dt[!is.na(rho)])))

# ---- Write correlation table -------------------------------------------------
out_tsv <- file.path(args$output_dir,
                     sprintf("covariate_correlations_%s.tsv", lbl))
fwrite(cor_dt, out_tsv, sep = "\t")
cat("Correlations ->", out_tsv, "\n")

# ---- Dot plot ----------------------------------------------------------------
cor_dt[, cov_f    := factor(covariate, levels = rev(cov_labels))]
cor_dt[, metric_f := factor(metric,    levels = imbalance_cols)]

p_dot <- ggplot(
    cor_dt[!is.na(rho)],
    aes(x      = rho,
        y      = cov_f,
        colour = significant,
        size   = pmin(-log10(pval + 1e-300), 10))
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_point(alpha = 0.9) +
  facet_wrap(~metric_f, nrow = 1L, scales = "free_x") +
  scale_colour_manual(
    name   = NULL,
    values = c("TRUE" = "#D6604D", "FALSE" = "grey60"),
    labels = c("TRUE" = sprintf("p.adj < %.2g", args$pval_cutoff),
               "FALSE" = "n.s.")
  ) +
  scale_size_continuous(
    name  = expression(-log[10](p)),
    range = c(2, 6)
  ) +
  labs(
    title    = sprintf("Quintile imbalance ~ batch covariates  |  %s", lbl),
    subtitle = sprintf(
      "Spearman rho, BH-adjusted across covariates per metric  |  n = %d samples",
      nrow(merged_dt)),
    x = "Spearman rho",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold"),
    strip.background   = element_rect(fill = "grey92"),
    strip.text         = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

out_dot <- file.path(args$output_dir,
                     sprintf("covariate_dotplot_%s.pdf", lbl))
pdf(out_dot, width = 12, height = 5)
print(p_dot)
invisible(dev.off())
cat("Dot plot ->", out_dot, "\n")

# ---- Scatter plots: one PDF page per metric, one panel per covariate ---------
out_scatter <- file.path(args$output_dir,
                         sprintf("covariate_scatter_%s.pdf", lbl))
pdf(out_scatter, width = 16, height = 7)

for (mc in imbalance_cols) {
  y_all     <- merged_dt[[mc]]
  plot_list <- vector("list", length(cov_cols))
  ann_list  <- vector("list", length(cov_cols))

  for (ci in seq_along(cov_cols)) {
    x_all <- merged_dt[[cov_cols[ci]]]
    keep  <- complete.cases(x_all, y_all)
    plot_list[[ci]] <- data.table(
      x       = x_all[keep],
      y       = y_all[keep],
      cov_lbl = cov_labels[ci]
    )
    rr <- cor_dt[metric == mc & covariate == cov_labels[ci]]
    ann_list[[ci]] <- data.table(
      cov_lbl = cov_labels[ci],
      ann     = if (nrow(rr) > 0L && !is.na(rr$rho))
                  sprintf("rho = %.3f\np = %.2g", rr$rho, rr$pval)
                else "---"
    )
  }
  plot_dt <- rbindlist(plot_list)
  ann_dt  <- rbindlist(ann_list)
  plot_dt[, cov_f := factor(cov_lbl, levels = cov_labels)]
  ann_dt[,  cov_f := factor(cov_lbl, levels = cov_labels)]

  p_sc <- ggplot(plot_dt, aes(x = x, y = y)) +
    geom_jitter(alpha = 0.25, size = 0.6, width = 0.05,
                colour = "#2166AC") +
    geom_smooth(method = "lm", se = TRUE, colour = "#D6604D",
                linewidth = 0.7, fill = "#FDDBC7", alpha = 0.4) +
    geom_text(
      data        = ann_dt,
      mapping     = aes(x = -Inf, y = Inf, label = ann),
      hjust       = -0.1,
      vjust       = 1.2,
      size        = 2.8,
      colour      = "grey20",
      inherit.aes = FALSE
    ) +
    facet_wrap(~cov_f, ncol = 5L, scales = "free_x") +
    labs(
      title    = sprintf("%s ~ batch covariates  |  %s", mc, lbl),
      subtitle = "Spearman rho annotated; linear smoother shown",
      x        = "covariate value",
      y        = mc
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold", size = 8)
    )

  print(p_sc)
}

invisible(dev.off())
cat("Scatter plots ->", out_scatter, "\n")

cat("\n=== Done ===\n")
