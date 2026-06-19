#!/usr/bin/env Rscript
# ==============================================================================
# PCA Candidate Selection
#
# Reads pre-computed PCA scores (*_pca_scores.xlsx from plot_pca_l2l1_variability.R)
# and selects candidate genes in two steps:
#
#   Step 1 — Quadrant selection:
#     Baseline genes define thresholds on each axis (min or max of their PC scores
#     depending on the quadrant). Candidates must be strictly more extreme on both
#     axes simultaneously.
#
#   Step 2 — Direction band (optional, --direction-band):
#     Slope β is estimated from the baseline gene coordinates.  Primary candidates
#     within ±w of the line  y = β·x  become secondary (direction) candidates.
#
# Outputs:
#   {prefix}_candidates.xlsx        — all genes with selection flags
#   {prefix}_selection_steps.pdf    — 2- or 3-panel stepwise figure
#   {prefix}_selection_params.txt   — threshold / β / count summary
#
# Usage:
#   Rscript select_pca_candidates.R \
#     --scores-file  results/pca_gene_metrics/ct_pca_scores.xlsx \
#     --pc-x PC1 --pc-y PC2 \
#     --quadrant Q3 \
#     --baseline-genes "hsp83,Hsc70-4,Hsc70-5" \
#     --direction-band --width 0.5 \
#     --gene-list data/candidates/01_devcell_2025_gProfiler_hsapiens_dmelanogaster_0417.csv \
#     --output-dir results/pca_candidates \
#     --prefix ct_q3
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(argparse)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

# ----- Safe xlsx writer (never overwrites) -----
write_xlsx_safe <- function(dt_list, path) {
  if (file.exists(path))
    stop("Output file already exists: ", path, "\nRemove it before re-running.")
  wb <- createWorkbook()
  for (nm in names(dt_list)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, as.data.frame(dt_list[[nm]]),
              headerStyle = createStyle(textDecoration = "bold"))
    setColWidths(wb, nm, cols = seq_len(ncol(dt_list[[nm]])), widths = "auto")
  }
  saveWorkbook(wb, path, overwrite = FALSE)
}

# ==============================================================================
# Data loaders
# ==============================================================================

load_scores <- function(path) {
  cat("Loading scores from:", path, "\n")
  dt <- as.data.table(read.xlsx(path))
  missing <- setdiff(c("gene_id", "SYMBOL"), colnames(dt))
  if (length(missing))
    stop("Scores file missing required columns: ", paste(missing, collapse = ", "))
  pc_cols <- grep("^PC[0-9]+$", colnames(dt), value = TRUE)
  if (length(pc_cols) == 0)
    stop("No PC columns found (expected PC1, PC2, ...)")
  # Drop annotation columns written by plot_pca_l2l1_variability.R — this script
  # rebuilds them from --gene-list so carrying them through would cause confusion.
  dt[, c("is_annotated", "label") := NULL]
  cat("  Genes:", nrow(dt), "| Available PCs:", paste(pc_cols, collapse = ", "), "\n")
  dt
}

# --baseline-genes: comma-separated string OR path to one-symbol-per-line file
parse_baseline_genes <- function(arg) {
  if (file.exists(arg)) {
    syms <- trimws(readLines(arg))
    syms <- syms[nzchar(syms)]
    cat("  Baseline genes (file):", length(syms), "\n")
    return(syms)
  }
  syms <- trimws(strsplit(arg, ",")[[1]])
  syms <- syms[nzchar(syms)]
  cat("  Baseline genes (arg):", length(syms), "—",
      paste(syms, collapse = ", "), "\n")
  syms
}

# --gene-list: format detected from file content, not extension:
#   (a) gProfiler table  — contains columns initial_alias + ortholog_name
#   (b) Tabular gene list — contains a SYMBOL column (gene_id optional)
#   (c) Plain text        — one gene symbol per line (fallback)
load_gene_annotations <- function(path) {
  if (is.null(path)) return(character(0))

  # Try to read as a delimited table first
  tbl <- tryCatch(fread(path, data.table = TRUE, showProgress = FALSE),
                  error = function(e) NULL)

  # (a) gProfiler: has ortholog_name (+ initial_alias)
  if (!is.null(tbl) && "ortholog_name" %in% colnames(tbl)) {
    syms <- unique(tbl$ortholog_name[!is.na(tbl$ortholog_name) &
                                      nzchar(tbl$ortholog_name) &
                                      tbl$ortholog_name != "N/A"])
    cat("  Gene list (gProfiler — ortholog_name):", length(syms), "symbols\n")
    return(syms)
  }

  # (b) Tabular with SYMBOL column
  if (!is.null(tbl) && "SYMBOL" %in% colnames(tbl)) {
    syms <- unique(tbl$SYMBOL[!is.na(tbl$SYMBOL) & nzchar(tbl$SYMBOL)])
    cat("  Gene list (tabular — SYMBOL column):", length(syms), "symbols\n")
    return(syms)
  }

  # (c) Plain text: one symbol per line
  raw  <- trimws(readLines(path))
  raw  <- raw[nzchar(raw)]
  keep <- grepl("^[A-Za-z][A-Za-z0-9_.\\-]*$", raw) &
          nchar(raw) <= 30 &
          !grepl("[0-9]{4}", raw)
  if (any(!keep))
    cat("  Skipped non-gene lines:", paste(raw[!keep], collapse = ", "), "\n")
  cat("  Gene list (plain text):", sum(keep), "symbols\n")
  raw[keep]
}

# ==============================================================================
# Selection logic
# ==============================================================================

# Returns list(thr_x, thr_y, op_x, op_y) based on baseline genes' extreme scores
get_quadrant_thresholds <- function(dt, baseline_syms, quadrant, pc_x, pc_y) {
  sym_uc  <- toupper(dt$SYMBOL)
  base_uc <- toupper(baseline_syms)
  mask    <- sym_uc %in% base_uc

  missing <- baseline_syms[!base_uc %in% sym_uc]
  if (length(missing))
    warning("Baseline genes not found in scores file: ", paste(missing, collapse = ", "))
  if (sum(mask) == 0)
    stop("No baseline genes matched in scores file. Check SYMBOL spelling.")

  cat("  Baseline matched:", sum(mask), "/", length(baseline_syms),
      ":", paste(dt$SYMBOL[mask], collapse = ", "), "\n")

  bx <- dt[[pc_x]][mask]
  by <- dt[[pc_y]][mask]

  thr <- switch(quadrant,
    Q1 = list(thr_x = max(min(bx), 0), thr_y = max(min(by), 0), op_x = ">",  op_y = ">"),
    Q2 = list(thr_x = min(max(bx), 0), thr_y = max(min(by), 0), op_x = "<",  op_y = ">"),
    Q3 = list(thr_x = min(max(bx), 0), thr_y = min(max(by), 0), op_x = "<",  op_y = "<"),
    Q4 = list(thr_x = max(min(bx), 0), thr_y = max(min(by), 0), op_x = ">",  op_y = "<"),
    stop("--quadrant must be Q1, Q2, Q3, or Q4")
  )
  cat(sprintf("  Thresholds (%s): %s %s %.4f  |  %s %s %.4f\n",
              quadrant, pc_x, thr$op_x, thr$thr_x, pc_y, thr$op_y, thr$thr_y))
  thr
}

# Logical mask: genes strictly beyond both thresholds, excluding baseline genes
select_quadrant_candidates <- function(dt, baseline_mask, thr, pc_x, pc_y) {
  x <- dt[[pc_x]]
  y <- dt[[pc_y]]
  pass_x <- switch(thr$op_x, ">" = x > thr$thr_x, "<" = x < thr$thr_x)
  pass_y <- switch(thr$op_y, ">" = y > thr$thr_y, "<" = y < thr$thr_y)
  mask   <- pass_x & pass_y & !baseline_mask & !is.na(x) & !is.na(y)
  cat("  Primary candidates:", sum(mask), "\n")
  mask
}

# Estimates slope β from baseline gene PC coordinates
fit_direction_slope <- function(dt, baseline_mask, pc_x, pc_y, method) {
  bx <- dt[[pc_x]][baseline_mask]
  by <- dt[[pc_y]][baseline_mask]

  beta <- switch(method,
    "ols" = {
      if (sum(baseline_mask) < 2)
        stop("At least 2 baseline genes required for 'ols' slope estimation")
      coef(lm(by ~ bx))[["bx"]]
    },
    "ols-origin" = {
      sum(by * bx) / sum(bx^2)
    },
    "unit-vector" = {
      r <- sqrt(bx^2 + by^2)
      ok <- r > 0
      mean(by[ok] / r[ok]) / mean(bx[ok] / r[ok])
    },
    stop("--slope-method must be ols, ols-origin, or unit-vector")
  )
  cat(sprintf("  Direction slope β = %.4f  (method: %s)\n", beta, method))
  beta
}

# Logical mask over ALL genes: within ±width of y=β·x
band_mask_all <- function(dt, beta, width, pc_x, pc_y) {
  x    <- dt[[pc_x]]
  y    <- dt[[pc_y]]
  abs(y - beta * x) <= width
}

# Returns logical mask: TRUE = gene is lncRNA / asRNA / snRNA / hpRNA
is_lncrna_gene <- function(symbols) {
  grepl("^(lncRNA|asRNA|snRNA|hpRNA):", symbols)
}

# Returns logical mask: TRUE = gene is ribosomal protein (RpS / RpL)
is_ribo_gene <- function(symbols) {
  grepl("^Rp[SsLl]", symbols)
}

# Removes flagged gene types from candidate_mask and reports counts
apply_gene_type_filters <- function(candidate_mask, dt,
                                    exclude_lncrna, exclude_ribo) {
  lncrna_mask <- is_lncrna_gene(dt$SYMBOL)
  ribo_mask   <- is_ribo_gene(dt$SYMBOL)

  if (exclude_lncrna) {
    n_removed <- sum(candidate_mask & lncrna_mask)
    if (n_removed > 0) {
      removed <- dt$SYMBOL[candidate_mask & lncrna_mask]
      cat("  Excluded lncRNA/asRNA/snRNA/hpRNA candidates:", n_removed,
          "—", paste(removed, collapse = ", "), "\n")
    } else {
      cat("  Excluded lncRNA/asRNA/snRNA/hpRNA candidates: 0\n")
    }
    candidate_mask <- candidate_mask & !lncrna_mask
  }

  if (exclude_ribo) {
    n_removed <- sum(candidate_mask & ribo_mask)
    if (n_removed > 0) {
      removed <- dt$SYMBOL[candidate_mask & ribo_mask]
      cat("  Excluded ribosomal candidates:", n_removed,
          "—", paste(removed, collapse = ", "), "\n")
    } else {
      cat("  Excluded ribosomal candidates: 0\n")
    }
    candidate_mask <- candidate_mask & !ribo_mask
  }

  cat("  Candidates after type filters:", sum(candidate_mask), "\n")
  candidate_mask
}

# ==============================================================================
# Plot helpers
# ==============================================================================

COLORS <- list(
  bg        = "grey78",
  genelist  = "#9E9AC8",
  baseline  = "#E6811B",
  candidate = "#D6604D",
  secondary = "#41AB5D",
  axis      = "grey55",
  threshold = "#4A4A8A",
  direction = "#1B6B31"
)

# Shared axis limits with 5 % padding
make_lims <- function(dt, pc_x, pc_y) {
  xr <- range(dt[[pc_x]], na.rm = TRUE)
  yr <- range(dt[[pc_y]], na.rm = TRUE)
  xp <- diff(xr) * 0.05
  yp <- diff(yr) * 0.05
  list(xlim = c(xr[1] - xp, xr[2] + xp),
       ylim = c(yr[1] - yp, yr[2] + yp))
}

# Base scatter of all genes
base_scatter <- function(dt, genelist_mask, pc_x, pc_y) {
  bg <- as.data.frame(dt[!genelist_mask])
  gl <- as.data.frame(dt[genelist_mask])
  layers <- list(
    geom_point(data = bg, aes(x = .data[[pc_x]], y = .data[[pc_y]]),
               color = COLORS$bg, alpha = 0.40, size = 0.6),
    if (nrow(gl) > 0)
      geom_point(data = gl, aes(x = .data[[pc_x]], y = .data[[pc_y]]),
                 color = COLORS$genelist, alpha = 0.70, size = 1.2)
  )
  Filter(Negate(is.null), layers)
}

# Crosshair at origin
origin_lines <- function() {
  list(
    geom_hline(yintercept = 0, linetype = "dashed",
               color = COLORS$axis, linewidth = 0.35),
    geom_vline(xintercept = 0, linetype = "dashed",
               color = COLORS$axis, linewidth = 0.35)
  )
}

# Baseline points + repelled labels
baseline_layer <- function(dt, baseline_mask, pc_x, pc_y) {
  base_df <- as.data.frame(dt[baseline_mask])
  list(
    geom_point(data = base_df, aes(x = .data[[pc_x]], y = .data[[pc_y]]),
               color = COLORS$baseline, size = 2.8),
    geom_text_repel(data = base_df,
                    aes(x = .data[[pc_x]], y = .data[[pc_y]], label = SYMBOL),
                    color = COLORS$baseline, size = 3.2, fontface = "bold",
                    max.overlaps = 80, box.padding = 0.45,
                    segment.color = "grey48", segment.size = 0.3)
  )
}

# Candidate labels (with optional genelist highlight color)
candidate_labels <- function(cand_df, genelist_uc, pc_x, pc_y, label_n) {
  if (nrow(cand_df) == 0) return(list())
  if (!is.null(label_n) && nrow(cand_df) > label_n) {
    cand_df$dist_ <- sqrt(cand_df[[pc_x]]^2 + cand_df[[pc_y]]^2)
    cand_df       <- cand_df[order(cand_df$dist_, decreasing = TRUE)[seq_len(label_n)], ]
  }
  cand_df$lbl_col <- ifelse(toupper(cand_df$SYMBOL) %in% genelist_uc,
                            "#7B2D8B", COLORS$candidate)
  list(
    geom_text_repel(data = cand_df,
                    aes(x = .data[[pc_x]], y = .data[[pc_y]], label = SYMBOL),
                    color = cand_df$lbl_col, size = 2.8, max.overlaps = 120,
                    box.padding = 0.30,
                    segment.color = "grey55", segment.size = 0.3)
  )
}

# ------ Panel 1 ------
plot_step1 <- function(dt, baseline_mask, genelist_mask, pc_x, pc_y, lims,
                       title = "Step 1: All genes + baseline") {
  ggplot() +
    base_scatter(dt, genelist_mask, pc_x, pc_y) +
    origin_lines() +
    baseline_layer(dt, baseline_mask, pc_x, pc_y) +
    coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    labs(title = title, x = pc_x, y = pc_y) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          legend.position = "none")
}

# ------ Panel 2 ------
plot_step2 <- function(dt, baseline_mask, candidate_mask, genelist_mask,
                       thr, pc_x, pc_y, lims, label_n = NULL,
                       title = "Step 2: Quadrant selection") {

  shade_x <- switch(thr$op_x,
    ">" = c(thr$thr_x, lims$xlim[2]),
    "<" = c(lims$xlim[1], thr$thr_x))
  shade_y <- switch(thr$op_y,
    ">" = c(thr$thr_y, lims$ylim[2]),
    "<" = c(lims$ylim[1], thr$thr_y))
  shade_df <- data.frame(xmin = shade_x[1], xmax = shade_x[2],
                         ymin = shade_y[1], ymax = shade_y[2])

  cand_df  <- as.data.frame(dt[candidate_mask])
  gl_uc    <- toupper(dt$SYMBOL[genelist_mask])

  # annotation positions for threshold labels
  ann_x_pos <- if (thr$op_x == "<")
    thr$thr_x - diff(lims$xlim) * 0.02 else thr$thr_x + diff(lims$xlim) * 0.02
  ann_x_hjust <- if (thr$op_x == "<") 1 else 0
  ann_y_pos <- thr$thr_y + diff(lims$ylim) * (if (thr$op_y == "<") -0.02 else 0.02)

  ggplot() +
    geom_rect(data = shade_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = COLORS$candidate, alpha = 0.07, inherit.aes = FALSE) +
    base_scatter(dt, genelist_mask, pc_x, pc_y) +
    geom_point(data = cand_df, aes(x = .data[[pc_x]], y = .data[[pc_y]]),
               color = COLORS$candidate, size = 1.8, alpha = 0.85) +
    origin_lines() +
    geom_vline(xintercept = thr$thr_x, linetype = "dashed",
               color = COLORS$threshold, linewidth = 0.65) +
    geom_hline(yintercept = thr$thr_y, linetype = "dashed",
               color = COLORS$threshold, linewidth = 0.65) +
    annotate("text", x = ann_x_pos, y = lims$ylim[2] - diff(lims$ylim) * 0.02,
             label = sprintf("%s = %.2f", pc_x, thr$thr_x),
             size = 2.6, color = COLORS$threshold, hjust = ann_x_hjust, vjust = 1) +
    annotate("text", x = lims$xlim[2] - diff(lims$xlim) * 0.02, y = ann_y_pos,
             label = sprintf("%s = %.2f", pc_y, thr$thr_y),
             size = 2.6, color = COLORS$threshold, hjust = 1) +
    candidate_labels(cand_df, gl_uc, pc_x, pc_y, label_n) +
    baseline_layer(dt, baseline_mask, pc_x, pc_y) +
    coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    labs(title    = title,
         subtitle = sprintf("n = %d primary candidates", sum(candidate_mask, na.rm = TRUE)),
         x = pc_x, y = pc_y) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          legend.position = "none")
}

# ------ Panel 3 ------
plot_step3 <- function(dt, baseline_mask, candidate_mask, secondary_mask,
                       genelist_mask, beta, width, pc_x, pc_y, lims,
                       label_n = NULL, title = "Step 3: Direction band") {

  # Shaded band polygon clipped to axis range
  bx_seq   <- seq(lims$xlim[1], lims$xlim[2], length.out = 300)
  yc       <- beta * bx_seq
  band_df  <- data.frame(
    x = c(bx_seq, rev(bx_seq)),
    y = c(yc + width, rev(yc - width))
  )
  line_df  <- data.frame(
    x = lims$xlim,
    y = beta * lims$xlim
  )

  nonsec_mask <- candidate_mask & !secondary_mask
  cand_nosec  <- as.data.frame(dt[nonsec_mask])
  sec_df      <- as.data.frame(dt[secondary_mask])
  gl_uc       <- toupper(dt$SYMBOL[genelist_mask])

  # Label position for β annotation
  ann_x <- lims$xlim[2] - diff(lims$xlim) * 0.04
  ann_y <- beta * ann_x + width + diff(lims$ylim) * 0.04

  ggplot() +
    geom_polygon(data = band_df, aes(x = x, y = y),
                 fill = COLORS$secondary, alpha = 0.13, inherit.aes = FALSE) +
    base_scatter(dt, genelist_mask, pc_x, pc_y) +
    geom_point(data = cand_nosec, aes(x = .data[[pc_x]], y = .data[[pc_y]]),
               color = COLORS$candidate, size = 1.5, alpha = 0.35) +
    geom_point(data = sec_df, aes(x = .data[[pc_x]], y = .data[[pc_y]]),
               color = COLORS$secondary, size = 2.0, alpha = 0.90) +
    origin_lines() +
    geom_line(data = line_df, aes(x = x, y = y),
              color = COLORS$direction, linewidth = 0.75, linetype = "dashed") +
    annotate("text", x = ann_x, y = ann_y,
             label = sprintf("y = %.2f·x", beta),
             size = 2.8, color = COLORS$direction, hjust = 1) +
    candidate_labels(sec_df, gl_uc, pc_x, pc_y, label_n) +
    baseline_layer(dt, baseline_mask, pc_x, pc_y) +
    coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    labs(title    = title,
         subtitle = sprintf("n = %d direction candidates  (±%.2f band)",
                            sum(secondary_mask, na.rm = TRUE), width),
         x = pc_x, y = pc_y) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          legend.position = "none")
}

# ==============================================================================
# CLI
# ==============================================================================

parser <- ArgumentParser(
  description = paste(
    "Select PCA candidate genes by quadrant position and optional direction band.",
    "Reads *_pca_scores.xlsx produced by plot_pca_l2l1_variability.R."
  )
)
parser$add_argument("--scores-file",    required = TRUE,
  help = "PCA scores xlsx (*_pca_scores.xlsx)")
parser$add_argument("--pc-x",           default = "PC1",
  help = "PC column for x-axis [default: PC1]")
parser$add_argument("--pc-y",           default = "PC2",
  help = "PC column for y-axis [default: PC2]")
parser$add_argument("--quadrant",       default = "Q3",
  help = "Selection quadrant: Q1 (+x+y)  Q2 (-x+y)  Q3 (-x-y)  Q4 (+x-y) [default: Q3]")
parser$add_argument("--baseline-genes", required = TRUE,
  help = "Comma-separated SYMBOL list (e.g. 'hsp83,Hsc70-4') or path to one-per-line file")
parser$add_argument("--direction-band", action = "store_true", default = FALSE,
  help = "Enable step-2 direction band selection")
parser$add_argument("--width",          type = "double", default = 0.5,
  help = "Half-width of direction band in PC units [default: 0.5]")
parser$add_argument("--slope-method",   default = "ols",
  help = "Slope estimation: ols | ols-origin | unit-vector [default: ols]")
parser$add_argument("--gene-list",      default = NULL, nargs = "+",
  help = "One or more annotation files (space-separated). Each creates its own boolean column named after the first chars of the filename. All are unioned for plot highlighting.")
parser$add_argument("--prefix",         default = "pca_candidates",
  help = "Output file prefix [default: pca_candidates]")
parser$add_argument("--output-dir",     default = ".",
  help = "Output directory [default: current directory]")
parser$add_argument("--label-n",        type = "integer", default = NULL,
  help = "Cap on candidate gene labels (default: label all)")
parser$add_argument("--keep-lncrna",    action = "store_true", default = FALSE,
  help = "Keep lncRNA/asRNA/snRNA/hpRNA genes in candidates (excluded by default)")
parser$add_argument("--keep-ribo",      action = "store_true", default = FALSE,
  help = "Keep ribosomal protein genes (RpS/RpL) in candidates (excluded by default)")

args <- parser$parse_args()

cat("=== PCA Candidate Selection ===\n")
cat(sprintf("  %-20s %s\n",  "scores-file:",    args$scores_file))
cat(sprintf("  %-20s %s  vs  %s\n", "axes:", args$pc_x, args$pc_y))
cat(sprintf("  %-20s %s\n",  "quadrant:",       args$quadrant))
cat(sprintf("  %-20s %s\n",  "baseline-genes:", args$baseline_genes))
cat(sprintf("  %-20s %s\n",  "direction-band:", args$direction_band))
if (args$direction_band)
  cat(sprintf("  %-20s %.3f  (method: %s)\n", "width / method:",
              args$width, args$slope_method))
cat(sprintf("  %-20s %s\n",  "output-dir:",     args$output_dir))
cat(sprintf("  %-20s %s\n",  "prefix:",         args$prefix))
if (!is.null(args$label_n))
  cat(sprintf("  %-20s %d\n", "label-n:", args$label_n))

# ----- Validate -----
if (!file.exists(args$scores_file))
  stop("Scores file not found: ", args$scores_file)
if (!args$quadrant %in% c("Q1", "Q2", "Q3", "Q4"))
  stop("--quadrant must be Q1, Q2, Q3, or Q4")
if (!args$slope_method %in% c("ols", "ols-origin", "unit-vector"))
  stop("--slope-method must be ols, ols-origin, or unit-vector")

if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = TRUE)
prefix <- file.path(args$output_dir, args$prefix)

# ----- Load -----
cat("\nLoading data...\n")
dt            <- load_scores(args$scores_file)
baseline_syms <- parse_baseline_genes(args$baseline_genes)

# Load each annotation file → named list col_name → character vector of symbols
# Column name = first 10 chars of filename (no extension); de-duplicated if needed
annot_list <- list()
if (!is.null(args$gene_list)) {
  for (gpath in args$gene_list) {
    syms <- load_gene_annotations(gpath)
    if (length(syms) == 0) next
    col_nm <- substr(sub("\\.[^.]+$", "", basename(gpath)), 1, 10)
    col_nm <- make.names(col_nm)
    if (col_nm %in% names(annot_list))
      col_nm <- make.unique(c(names(annot_list), col_nm), sep = "_")[length(annot_list) + 1]
    annot_list[[col_nm]] <- syms
  }
}

for (pc in c(args$pc_x, args$pc_y))
  if (!pc %in% colnames(dt))
    stop(pc, " not found in scores file. Available: ",
         paste(grep("^PC[0-9]+$", colnames(dt), value = TRUE), collapse = ", "))

# ----- Masks -----
cat("\nRunning selection...\n")
sym_uc        <- toupper(dt$SYMBOL)
baseline_mask <- sym_uc %in% toupper(baseline_syms)
# Union of all annotation masks — used for plot label highlighting
genelist_mask <- if (length(annot_list) > 0) {
  Reduce("|", lapply(annot_list, function(s) sym_uc %in% toupper(s)))
} else {
  rep(FALSE, nrow(dt))
}

thr            <- get_quadrant_thresholds(dt, baseline_syms, args$quadrant,
                                          args$pc_x, args$pc_y)
candidate_mask <- select_quadrant_candidates(dt, baseline_mask, thr,
                                             args$pc_x, args$pc_y)

# ----- Gene-type filters (applied before direction band) -----
exclude_lncrna <- !args$keep_lncrna
exclude_ribo   <- !args$keep_ribo
if (exclude_lncrna || exclude_ribo)
  candidate_mask <- apply_gene_type_filters(candidate_mask, dt,
                                            exclude_lncrna, exclude_ribo)

beta           <- NA_real_
secondary_mask <- rep(FALSE, nrow(dt))
if (args$direction_band) {
  beta           <- fit_direction_slope(dt, baseline_mask, args$pc_x, args$pc_y,
                                        args$slope_method)
  in_band        <- band_mask_all(dt, beta, args$width, args$pc_x, args$pc_y)
  secondary_mask <- candidate_mask & in_band
  cat("  Secondary candidates (primary ∩ band):", sum(secondary_mask), "\n")
}

# ----- Build output table -----
out_dt <- copy(dt)
out_dt[, dist_from_origin   := sqrt(.SD[[1]]^2 + .SD[[2]]^2),
       .SDcols = c(args$pc_x, args$pc_y)]
out_dt[, is_lncrna           := is_lncrna_gene(SYMBOL)]
out_dt[, is_ribo             := is_ribo_gene(SYMBOL)]
out_dt[, is_baseline         := baseline_mask]
out_dt[, is_candidate        := candidate_mask]
out_dt[, is_annotated         := genelist_mask]
if (args$direction_band)
  out_dt[, is_direction_candidate := secondary_mask]
# One boolean column per annotation file
annot_cols <- character(0)
for (col_nm in names(annot_list)) {
  out_dt[, (col_nm) := sym_uc %in% toupper(annot_list[[col_nm]])]
  annot_cols <- c(annot_cols, col_nm)
}
# Reorder columns: identifiers + flagged PCs + flags + annotation cols + everything else
pc_cols   <- grep("^PC[0-9]+$", colnames(out_dt), value = TRUE)
flag_cols <- c("is_lncrna", "is_ribo", "is_baseline", "is_candidate",
               if (args$direction_band) "is_direction_candidate",
               annot_cols)
front     <- c("gene_id", "SYMBOL", args$pc_x, args$pc_y,
               "dist_from_origin", flag_cols)
rest      <- setdiff(colnames(out_dt), c(front, pc_cols))
col_order <- unique(c(front, setdiff(pc_cols, c(args$pc_x, args$pc_y)), rest))
out_dt    <- out_dt[, ..col_order]

# ----- Plots -----
cat("\nBuilding plots...\n")
lims <- make_lims(dt, args$pc_x, args$pc_y)

p1 <- plot_step1(dt, baseline_mask, genelist_mask,
                 args$pc_x, args$pc_y, lims,
                 title = "Step 1: All genes + baseline")

p2 <- plot_step2(dt, baseline_mask, candidate_mask, genelist_mask,
                 thr, args$pc_x, args$pc_y, lims,
                 label_n = args$label_n,
                 title = sprintf("Step 2: Quadrant %s selection", args$quadrant))

panel_list <- list(p1, p2)
n_panels   <- 2L

if (args$direction_band) {
  p3 <- plot_step3(dt, baseline_mask, candidate_mask, secondary_mask,
                   genelist_mask, beta, args$width,
                   args$pc_x, args$pc_y, lims,
                   label_n = args$label_n,
                   title = "Step 3: Direction band")
  panel_list <- c(panel_list, list(p3))
  n_panels   <- 3L
}

pdf_path <- paste0(prefix, "_selection_steps.pdf")
if (file.exists(pdf_path))
  stop("Output already exists: ", pdf_path, "\nRemove before re-running.")

combined <- plot_grid(plotlist   = panel_list,
                      labels     = LETTERS[seq_len(n_panels)],
                      label_size = 11,
                      nrow       = 1,
                      rel_widths = rep(1, n_panels))
ggsave(pdf_path, combined,
       width  = n_panels * 5.6,
       height = 5.6,
       device = "pdf")
cat("  Saved:", basename(pdf_path), "\n")

# ----- xlsx -----
cat("Writing candidates table...\n")
xlsx_path <- paste0(prefix, "_candidates.xlsx")
write_xlsx_safe(list(candidates = out_dt), xlsx_path)
cat("  Saved:", basename(xlsx_path), "\n")

# ----- Params summary -----
params_path <- paste0(prefix, "_selection_params.txt")
if (file.exists(params_path))
  stop("Output already exists: ", params_path, "\nRemove before re-running.")

lines <- c(
  "=== PCA Candidate Selection Parameters ===",
  sprintf("scores_file:      %s", args$scores_file),
  sprintf("pc_x:             %s", args$pc_x),
  sprintf("pc_y:             %s", args$pc_y),
  sprintf("quadrant:         %s", args$quadrant),
  sprintf("baseline_genes:   %s", paste(baseline_syms, collapse = ", ")),
  sprintf("baseline_found:   %s",
          paste(dt$SYMBOL[baseline_mask], collapse = ", ")),
  sprintf("threshold_%s:    %.6f", args$pc_x, thr$thr_x),
  sprintf("threshold_%s:    %.6f", args$pc_y, thr$thr_y),
  sprintf("n_primary:        %d", sum(candidate_mask)),
  sprintf("excl_lncrna:      %s", if (exclude_lncrna) "TRUE" else "FALSE (--keep-lncrna)"),
  sprintf("excl_ribo:        %s", if (exclude_ribo)   "TRUE" else "FALSE (--keep-ribo)"),
  "",
  if (args$direction_band)
    c("direction_band:   TRUE",
      sprintf("slope_method:     %s", args$slope_method),
      sprintf("beta:             %.6f", beta),
      sprintf("width:            %.4f", args$width),
      sprintf("n_direction:      %d", sum(secondary_mask)))
  else
    "direction_band:   FALSE"
)
writeLines(lines, params_path)
cat("  Saved:", basename(params_path), "\n")

cat("\nDone.\n")
