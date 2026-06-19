#!/usr/bin/env Rscript
# ==============================================================================
# DEG / DVG Method Comparison
#
# Compares differential-expression (DEG) results from three Q1-vs-Q5 analysis
# methods for the same focal gene and condition:
#   - limma-voom  (DEG only)
#   - DESeq2      (DEG only)
#   - GAMLSS      (DEG + DVG)
#
# Produces:
#   1. 3-way Venn diagrams for up-DEG and down-DEG overlap (PNG)
#   2. xlsx with 5 sheets:
#        summary_stats      -- DEG/DVG counts + % of total tested per method
#        upDEG_partitions   -- 7-partition overlap table for up-DEG
#        downDEG_partitions -- same for down-DEG
#        gene_list_upDEG    -- gene-level partition membership for up-DEG
#        gene_list_downDEG  -- same for down-DEG
#
# GAMLSS "Both" genes (significant in mean AND variance) are assigned to DEG
# based on logFC_mu direction and to DVG based on logFC_BCV direction.
#
# Supported GAMLSS label formats:
#   New (current code): up_deg / down_deg / up_dvg / down_dvg / Both / NS
#   Old (example files): DEG / DVG / Both / NS  (direction inferred from logFC)
#
# Usage:
#   Rscript compare_deg_dvg_methods.R \
#     --limma-file   results/limma_voom/FBgn0001235_CT_Q1vsQ5_limma_voom_deg.csv \
#     --deseq2-file  results/deseq2/FBgn0001235_CT_Q1vsQ5_deseq2_deg.csv \
#     --gamlss-file  results/gamlss/FBgn0001235_CT_Q1vsQ5_gamlss_dvgdeg.csv \
#     --output-xlsx  results/deg_comparison/FBgn0001235_CT_deg_comparison.xlsx \
#     --condition    CT \
#     --focal-gene   FBgn0001235
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
  library(openxlsx)
  library(RColorBrewer)
  library(VennDiagram)
})

# ---- Source shared Venn utility ---------------------------------------------
get_script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L)
    dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
  else
    getwd()
}
source(file.path(get_script_dir(), "shared", "utils_venn.R"))

# ---- CLI --------------------------------------------------------------------
parser <- ArgumentParser(description = "Compare DEG/DVG results across three methods")
parser$add_argument("--limma-file",   required = TRUE,
                    help = "limma-voom DEG CSV")
parser$add_argument("--deseq2-file",  required = TRUE,
                    help = "DESeq2 DEG CSV")
parser$add_argument("--gamlss-file",  required = TRUE,
                    help = "GAMLSS DEG+DVG CSV")
parser$add_argument("--output-xlsx",  required = TRUE,
                    help = "Output .xlsx path")
parser$add_argument("--condition",    required = TRUE,
                    help = "Condition label (CT / HS) — used in titles and filenames")
parser$add_argument("--venn-dir",     default = NULL,
                    help = "Directory for Venn PNG files (default: dirname of --output-xlsx)")
parser$add_argument("--gamlss-class", default = "modelcomp",
                    choices = c("modelcomp", "single"),
                    help = "Which GAMLSS class column to use (default: modelcomp)")
parser$add_argument("--focal-gene",   default = "",
                    help = "Focal gene label — optional; included in titles and filenames")
args <- parser$parse_args()

cat("=== DEG/DVG Method Comparison ===\n")
cat("limma file:   ", args$limma_file,  "\n")
cat("DESeq2 file:  ", args$deseq2_file, "\n")
cat("GAMLSS file:  ", args$gamlss_file, "\n")
cat("Output xlsx:  ", args$output_xlsx, "\n")
cat("Condition:    ", args$condition,   "\n")
cat("GAMLSS class: ", args$gamlss_class,"\n\n")

for (f in c(args$limma_file, args$deseq2_file, args$gamlss_file))
  if (!file.exists(f)) stop("File not found: ", f)

out_dir <- dirname(args$output_xlsx)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

venn_dir <- if (!is.null(args$venn_dir)) args$venn_dir else out_dir
dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

# Build a tag for output filenames and plot titles
gene_tag  <- if (nchar(args$focal_gene) > 0L) args$focal_gene else NULL
title_tag <- paste(c(gene_tag, args$condition), collapse = " | ")
file_tag  <- paste(c(gene_tag, args$condition), collapse = "_")

# ---- Helper: detect column suffix (_Q5 etc.) --------------------------------
detect_suffix <- function(dt, prefix) {
  col <- grep(paste0("^", prefix), names(dt), value = TRUE)[1L]
  if (is.na(col)) stop("Column with prefix '", prefix, "' not found in ", deparse(substitute(dt)))
  sub(paste0("^", prefix), "", col)
}

# ---- Helper: normalise GAMLSS class labels to directional form --------------
# Handles both current code labels (up_deg / down_dvg / …)
# and old example-file labels (DEG / DVG / Both / NS).
normalise_gamlss_class <- function(cls, logfc_mu, logfc_bcv) {
  data.table::fcase(
    cls == "up_deg",                            "up_deg",
    cls == "down_deg",                          "down_deg",
    cls == "up_dvg",                            "up_dvg",
    cls == "down_dvg",                          "down_dvg",
    cls == "Both"  & logfc_mu  >  0,            "both_up_deg",
    cls == "Both"  & logfc_mu  <= 0,            "both_down_deg",
    # old-format labels -------------------------------------------------------
    cls == "DEG"   & logfc_mu  >  0,            "up_deg",
    cls == "DEG"   & logfc_mu  <= 0,            "down_deg",
    cls == "DVG"   & logfc_bcv >  0,            "up_dvg",
    cls == "DVG"   & logfc_bcv <= 0,            "down_dvg",
    default = "NS"
  )
}

# ---- Helper: extract gene sets from one method's result data.table ----------
# Returns a named list:
#   up_deg, down_deg         — used in Venn comparison across methods
#   up_dvg, down_dvg, both   — used in summary stats (GAMLSS only)
extract_sets <- function(dt, method, gamlss_class_type = "modelcomp") {
  ids <- dt$gene_id

  if (method == "gamlss") {
    pfx           <- paste0("class_", gamlss_class_type, "_")
    sfx           <- detect_suffix(dt, pfx)
    class_col     <- paste0(pfx, sfx)
    logfc_mu_col  <- grep("^logFC_mu_",  names(dt), value = TRUE)[1L]
    logfc_bcv_col <- grep("^logFC_BCV_", names(dt), value = TRUE)[1L]

    norm <- normalise_gamlss_class(
      dt[[class_col]], dt[[logfc_mu_col]], dt[[logfc_bcv_col]]
    )

    list(
      up_deg   = ids[norm %in% c("up_deg",   "both_up_deg")],
      down_deg = ids[norm %in% c("down_deg", "both_down_deg")],
      up_dvg   = ids[norm %in% c("up_dvg",   "both_up_deg",   "both_down_deg")],
      down_dvg = ids[norm %in% c("down_dvg",  "both_up_deg",  "both_down_deg")],
      both     = ids[norm %in% c("both_up_deg", "both_down_deg")]
    )
    # Note: DVG direction for "Both" genes follows logFC_BCV:
    #   both_up_deg  could have up or down BCV — we include all Both in both up_dvg
    #   and down_dvg to capture them fully. Use the raw GAMLSS CSV for directional DVG.
  } else {
    sfx       <- detect_suffix(dt, "class_")
    class_col <- paste0("class_", sfx)
    cls       <- dt[[class_col]]
    list(
      up_deg   = ids[cls == "up_deg"],
      down_deg = ids[cls == "down_deg"],
      up_dvg   = character(0L),
      down_dvg = character(0L),
      both     = character(0L)
    )
  }
}

# ---- Helper: compute all 2^n-1 membership partitions -----------------------
compute_partitions <- function(sets_list) {
  n    <- length(sets_list)
  nms  <- names(sets_list)
  pool <- unique(unlist(sets_list))

  if (length(pool) == 0L) {
    empty_cols <- setNames(vector("list", n), nms)
    return(data.table::data.table(
      partition = character(0L), n_genes = integer(0L),
      gene_ids  = character(0L), .( do.call(data.table, empty_cols) )
    ))
  }

  # Membership matrix: genes x methods (logical)
  memb <- matrix(FALSE, nrow = length(pool), ncol = n,
                 dimnames = list(pool, nms))
  for (j in seq_len(n))
    memb[pool %in% sets_list[[j]], j] <- TRUE

  # Generate all 2^n membership patterns, drop all-FALSE row
  pat_grid  <- as.matrix(expand.grid(
    setNames(replicate(n, c(TRUE, FALSE), simplify = FALSE), nms)
  ))
  pat_grid  <- pat_grid[rowSums(pat_grid) > 0L, , drop = FALSE]

  rows <- lapply(seq_len(nrow(pat_grid)), function(i) {
    pat   <- as.logical(pat_grid[i, ])
    genes <- pool[apply(memb, 1L, function(x) all(x == pat))]
    c(
      list(
        partition = paste(nms[pat], collapse = " & "),
        n_genes   = length(genes),
        gene_ids  = paste(sort(genes), collapse = "|")
      ),
      setNames(as.list(as.integer(pat)), nms)
    )
  })

  data.table::rbindlist(lapply(rows, as.data.table), fill = TRUE)
}

# ---- Helper: write one xlsx sheet with bold header + auto width -------------
write_sheet <- function(wb, sheet_name, dt) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, as.data.frame(dt),
            headerStyle = createStyle(textDecoration = "bold"))
  setColWidths(wb, sheet_name, cols = seq_len(ncol(dt)), widths = "auto")
}

# ---- Load input files -------------------------------------------------------
cat("Loading result files...\n")
limma_dt  <- fread(args$limma_file,  data.table = TRUE)
deseq2_dt <- fread(args$deseq2_file, data.table = TRUE)
gamlss_dt <- fread(args$gamlss_file, data.table = TRUE)

n_genes_limma  <- nrow(limma_dt)
n_genes_deseq2 <- nrow(deseq2_dt)
n_genes_gamlss <- nrow(gamlss_dt)

cat("  limma:  ", n_genes_limma,  "genes\n")
cat("  DESeq2: ", n_genes_deseq2, "genes\n")
cat("  GAMLSS: ", n_genes_gamlss, "genes\n\n")

# ---- Extract gene sets ------------------------------------------------------
cat("Extracting gene sets...\n")
sets_limma  <- extract_sets(limma_dt,  "limma")
sets_deseq2 <- extract_sets(deseq2_dt, "deseq2")
sets_gamlss <- extract_sets(gamlss_dt, "gamlss", args$gamlss_class)

# ---- Summary stats ----------------------------------------------------------
pct <- function(n, total) ifelse(total > 0, round(100 * n / total, 2), NA_real_)

summary_dt <- data.table(
  method        = c("limma", "DESeq2", "GAMLSS"),
  total_tested  = c(n_genes_limma, n_genes_deseq2, n_genes_gamlss),
  up_DEG        = c(length(sets_limma$up_deg),
                    length(sets_deseq2$up_deg),
                    length(sets_gamlss$up_deg)),
  down_DEG      = c(length(sets_limma$down_deg),
                    length(sets_deseq2$down_deg),
                    length(sets_gamlss$down_deg)),
  up_DVG        = c(NA_integer_, NA_integer_, length(sets_gamlss$up_dvg)),
  down_DVG      = c(NA_integer_, NA_integer_, length(sets_gamlss$down_dvg)),
  Both          = c(NA_integer_, NA_integer_, length(sets_gamlss$both))
)

summary_dt[, DEG_total := up_DEG + down_DEG]
summary_dt[, DVG_total := fcase(
  method == "GAMLSS", up_DVG + down_DVG,
  default = NA_integer_
)]
summary_dt[, NS := total_tested - DEG_total - fcoalesce(
  fcase(method == "GAMLSS", DVG_total - Both, default = NA_integer_),
  0L
)]

summary_dt[, pct_up_DEG   := pct(up_DEG,   total_tested)]
summary_dt[, pct_down_DEG := pct(down_DEG, total_tested)]
summary_dt[, pct_up_DVG   := fcase(method == "GAMLSS", pct(up_DVG,   total_tested), default = NA_real_)]
summary_dt[, pct_down_DVG := fcase(method == "GAMLSS", pct(down_DVG, total_tested), default = NA_real_)]

setcolorder(summary_dt, c("method", "total_tested",
                           "up_DEG", "pct_up_DEG",
                           "down_DEG", "pct_down_DEG",
                           "DEG_total",
                           "up_DVG", "pct_up_DVG",
                           "down_DVG", "pct_down_DVG",
                           "DVG_total", "Both", "NS"))

cat("Summary stats:\n")
print(summary_dt[, .(method, up_DEG, down_DEG, up_DVG, down_DVG, Both)])
cat("\n")

# ---- 3-way partition tables -------------------------------------------------
cat("Computing overlap partitions...\n")

up_sets   <- list(limma  = sets_limma$up_deg,
                  DESeq2 = sets_deseq2$up_deg,
                  GAMLSS = sets_gamlss$up_deg)
down_sets <- list(limma  = sets_limma$down_deg,
                  DESeq2 = sets_deseq2$down_deg,
                  GAMLSS = sets_gamlss$down_deg)

up_part   <- compute_partitions(up_sets)
down_part <- compute_partitions(down_sets)

setorder(up_part,   -n_genes)
setorder(down_part, -n_genes)

cat("  Up-DEG:   ", nrow(up_part),   "partitions,",
    sum(up_part$n_genes), "total genes\n")
cat("  Down-DEG: ", nrow(down_part), "partitions,",
    sum(down_part$n_genes), "total genes\n\n")

# ---- Per-gene partition membership (long tables) ----------------------------
make_gene_list <- function(part_dt) {
  flag_cols  <- intersect(c("limma", "DESeq2", "GAMLSS"), names(part_dt))
  empty_args <- c(
    list(gene_id = character(0L), partition = character(0L)),
    setNames(replicate(length(flag_cols), integer(0L), simplify = FALSE), flag_cols)
  )
  empty <- do.call(data.table, empty_args)

  rows <- part_dt[n_genes > 0L]
  if (nrow(rows) == 0L) return(empty)

  dt <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    genes <- strsplit(rows$gene_ids[i], "|", fixed = TRUE)[[1L]]
    data.table(gene_id = genes, partition = rows$partition[i])
  }))

  flags <- rows[, c("partition", ..flag_cols)]
  dt    <- merge(dt, flags, by = "partition", all.x = TRUE)
  setorder(dt, partition, gene_id)
  dt
}

gene_list_up   <- make_gene_list(up_part)
gene_list_down <- make_gene_list(down_part)

# Drop gene_ids column from partition tables before writing
up_part_out   <- up_part[,   !c("gene_ids"), with = FALSE]
down_part_out <- down_part[, !c("gene_ids"), with = FALSE]

# ---- Venn diagrams ----------------------------------------------------------
cat("Drawing Venn diagrams...\n")

venn_up_path   <- file.path(venn_dir,
  sprintf("venn_upDEG_%s.png", file_tag))
venn_down_path <- file.path(venn_dir,
  sprintf("venn_downDEG_%s.png", file_tag))

draw_named_venn(
  sets_list   = up_sets,
  output_path = venn_up_path,
  title       = sprintf("Up-DEG overlap  |  %s", title_tag)
)
draw_named_venn(
  sets_list   = down_sets,
  output_path = venn_down_path,
  title       = sprintf("Down-DEG overlap  |  %s", title_tag)
)
cat("  Up-DEG Venn   ->", venn_up_path,   "\n")
cat("  Down-DEG Venn ->", venn_down_path, "\n\n")

# ---- Write xlsx -------------------------------------------------------------
cat("Writing xlsx...\n")
wb <- createWorkbook()

write_sheet(wb, "summary_stats",      summary_dt)
write_sheet(wb, "upDEG_partitions",   up_part_out)
write_sheet(wb, "downDEG_partitions", down_part_out)
write_sheet(wb, "gene_list_upDEG",    gene_list_up)
write_sheet(wb, "gene_list_downDEG",  gene_list_down)

saveWorkbook(wb, args$output_xlsx, overwrite = TRUE)
cat("  xlsx ->", args$output_xlsx, "\n\n")

# ---- Console summary --------------------------------------------------------
cat("=== Summary ===\n")
cat(sprintf("Condition: %s\n", args$condition))
if (nchar(args$focal_gene) > 0L)
  cat(sprintf("Focal gene: %s\n", args$focal_gene))
cat(sprintf("GAMLSS classification: %s\n\n", args$gamlss_class))

cat("Up-DEG total (union across methods):", sum(up_part$n_genes), "\n")
cat("Down-DEG total (union across methods):", sum(down_part$n_genes), "\n")
cat("All-three agreement (up)  :", up_part[partition == "limma & DESeq2 & GAMLSS", n_genes], "\n")
cat("All-three agreement (down):", down_part[partition == "limma & DESeq2 & GAMLSS", n_genes], "\n")
