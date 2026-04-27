#!/usr/bin/env Rscript
# ==============================================================================
# Permutation Null Distribution Plots
#
# For each focus gene, plot the null distribution (histogram) of each
# differential co-expression metric alongside the observed value (vertical
# red line). Produces one multi-page PDF per gene.
#
# Reads the *_null_dist.tsv.gz files produced by Stage 7 (07_permutation_test.py)
# — no R HDF5 dependency required.
#
# Usage:
#   Rscript 09_plot_permutation_null.R \
#     --null-tsv-dir results/permutation_null \
#     --pvals-tsv results/permutation_pvals.tsv \
#     --output-dir results/visualization_data/perm_null \
#     --condition-label "Control" \
#     --metrics "n_diff_edges,L1_rewire,L1_frac_rewire,mean_abs_delta"
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

# ----- Helper: significance stars -----
get_sig_stars <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# ----- Metric display names -----
METRIC_LABELS <- c(
  n_diff_edges   = "Differential edges (n)",
  focus_deg_low  = "Degree in LOW network",
  focus_deg_high = "Degree in HIGH network",
  mean_abs_delta = "Mean |Δr|",
  max_abs_delta  = "Max |Δr|",
  L1_n_nodes     = "L1 partner count",
  L1_rewire      = "L1 rewiring count",
  L1_frac_rewire = "L1 rewiring fraction",
  L2L1_deg       = "L2/L1 degree ratio",
  L2L1_rewire    = "L2/L1 rewiring ratio"
)

# ----- Argument parsing -----
parser <- ArgumentParser(description = "Plot permutation null distributions")
parser$add_argument("--null-tsv-dir", required = TRUE,
                    help = "Directory containing per-gene *_null_dist.tsv.gz files")
parser$add_argument("--pvals-tsv",    default = NULL,
                    help = "Optional: permutation_pvals.tsv for p-value annotation")
parser$add_argument("--output-dir",   required = TRUE,
                    help = "Directory for output PDF files")
parser$add_argument("--condition-label", default = "Condition",
                    help = "Human-readable condition label for plot titles")
parser$add_argument("--metrics", default = NULL,
                    help = "Comma-separated metric names to plot (default: all available)")
args <- parser$parse_args()

null_tsv_dir   <- args$null_tsv_dir
pvals_tsv_path <- args$pvals_tsv
output_dir     <- args$output_dir
condition_label <- args$condition_label
metrics_filter <- if (!is.null(args$metrics)) strsplit(args$metrics, ",")[[1]] else NULL

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----- Load p-values table if provided -----
pvals_dt <- NULL
if (!is.null(pvals_tsv_path) && file.exists(pvals_tsv_path)) {
  pvals_dt <- fread(pvals_tsv_path)
  cat(sprintf("Loaded p-values for %d genes from %s\n", nrow(pvals_dt), pvals_tsv_path))
}

# ----- Discover per-gene TSV files -----
tsv_files <- list.files(null_tsv_dir, pattern = "_null_dist\\.tsv\\.gz$", full.names = TRUE)
if (length(tsv_files) == 0) {
  stop(sprintf("No *_null_dist.tsv.gz files found in: %s", null_tsv_dir))
}
cat(sprintf("Found %d null distribution TSV files.\n", length(tsv_files)))

# ----- Process each gene -----
for (tsv_path in tsv_files) {
  gene_tag <- sub("_null_dist\\.tsv\\.gz$", "", basename(tsv_path))
  gene_id  <- sub("^[0-9]+_", "", gene_tag)  # strip leading NNNN_

  cat(sprintf("  Processing %s ...\n", gene_id))

  # Read null distributions
  null_dt <- tryCatch(
    fread(tsv_path),
    error = function(e) { warning(sprintf("  Cannot read %s: %s", tsv_path, e$message)); NULL }
  )
  if (is.null(null_dt) || nrow(null_dt) == 0) next

  # Determine metrics to plot
  avail_metrics <- setdiff(names(null_dt), c("gene_id", "perm_idx"))
  if (!is.null(metrics_filter)) {
    plot_metrics <- intersect(metrics_filter, avail_metrics)
  } else {
    plot_metrics <- avail_metrics
  }
  if (length(plot_metrics) == 0) {
    cat(sprintf("  No metrics to plot for %s. Skipping.\n", gene_id))
    next
  }

  # Look up p-values row for this gene
  pval_row <- NULL
  if (!is.null(pvals_dt)) {
    pval_row <- pvals_dt[gene_id == gene_id]
    if (nrow(pval_row) == 0) pval_row <- NULL
  }

  # Open PDF
  pdf_path <- file.path(output_dir, sprintf("%s_permutation_null.pdf", gene_tag))
  n_pages  <- length(plot_metrics)
  # Use a 2-column layout if many metrics
  ncol_lay <- if (n_pages > 4) 2L else 1L
  nrow_lay <- ceiling(n_pages / ncol_lay)
  pdf_width  <- 7 * ncol_lay
  pdf_height <- 4 * nrow_lay

  pdf(pdf_path, width = pdf_width, height = pdf_height)
  par_old <- par(no.readonly = TRUE)

  plot_list <- list()

  for (m in plot_metrics) {
    null_vals <- null_dt[[m]]
    n_perm <- length(null_vals)

    # Observed value
    obs_col <- paste0("obs_", m)
    if (!is.null(pval_row) && obs_col %in% names(pval_row)) {
      obs_val <- as.numeric(pval_row[[obs_col]][1])
    } else {
      obs_val <- NA_real_
    }

    # P-value
    pval_col <- paste0("pval_", m)
    if (!is.null(pval_row) && pval_col %in% names(pval_row)) {
      pval_val <- as.numeric(pval_row[[pval_col]][1])
    } else {
      pval_val <- NA_real_
    }

    # Null summary
    null_mean <- mean(null_vals, na.rm = TRUE)
    null_sd   <- sd(null_vals,   na.rm = TRUE)

    # Human-readable metric label
    m_label <- if (m %in% names(METRIC_LABELS)) METRIC_LABELS[[m]] else m

    # Build annotation string
    if (!is.na(pval_val)) {
      stars    <- get_sig_stars(pval_val)
      ann_text <- sprintf("p = %.4f %s\nNull: mean=%.3f, sd=%.3f\nn = %d permutations",
                          pval_val, stars, null_mean, null_sd, n_perm)
    } else {
      ann_text <- sprintf("Null: mean=%.3f, sd=%.3f\nn = %d permutations",
                          null_mean, null_sd, n_perm)
    }

    # Plot data frame
    plot_df <- data.frame(value = null_vals)

    # Determine x-range to include observed value
    x_min <- min(null_vals, na.rm = TRUE)
    x_max <- max(null_vals, na.rm = TRUE)
    if (!is.na(obs_val)) {
      x_min <- min(x_min, obs_val)
      x_max <- max(x_max, obs_val)
    }
    x_pad <- (x_max - x_min) * 0.1
    if (x_pad == 0) x_pad <- 1

    g <- ggplot(plot_df, aes(x = value)) +
      geom_histogram(fill = "grey65", colour = "white", bins = 30) +
      labs(
        title    = sprintf("%s — %s", gene_id, condition_label),
        subtitle = m_label,
        x        = m_label,
        y        = "Count (permutations)",
        caption  = ann_text
      ) +
      xlim(x_min - x_pad, x_max + x_pad) +
      theme_bw(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 9),
        plot.caption  = element_text(size = 8, hjust = 0)
      )

    # Add observed value line
    if (!is.na(obs_val)) {
      g <- g +
        geom_vline(xintercept = obs_val, colour = "firebrick", linewidth = 1, linetype = "dashed") +
        annotate("text",
                 x     = obs_val,
                 y     = Inf,
                 label = sprintf(" Observed:\n %.3f", obs_val),
                 hjust = -0.05, vjust = 1.5,
                 colour = "firebrick", size = 3)
    }

    plot_list[[m]] <- g
  }

  # Print all plots to PDF
  for (g in plot_list) {
    print(g)
  }

  dev.off()
  par(par_old)
  cat(sprintf("    Saved %s\n", pdf_path))
}

cat(sprintf("\nDone. PDFs written to %s\n", output_dir))
