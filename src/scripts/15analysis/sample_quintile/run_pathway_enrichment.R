#!/usr/bin/env Rscript
# ==============================================================================
# Generic Pathway Enrichment Script
#
# Runs GO (via clusterProfiler::enrichGO) and KEGG (via clusterProfiler::enricher
# with a locally cached pathway table) enrichment on any gene list provided as a
# TSV file.  Designed to accept output from any analysis in this project that
# writes a gene_id column with FlyBase IDs.
#
# KEY DESIGN CHOICES
#   - No progressive fallback: all enrichment runs at the specified --pval-cutoff
#     and --qval-cutoff so results from different calls are directly comparable.
#   - Consistent output: xlsx is always written (with "no results" message when
#     nothing is enriched); PDF is written only when ≥1 non-empty result exists;
#     a marker .txt file is also written when both GO and KEGG are empty, for
#     easy downstream detection.
#   - KEGG cache: pathway-to-gene table is downloaded once (download_KEGG) and
#     saved as an RDS; subsequent calls load from cache without network access.
#   - Background: if --background-file is omitted the full gene list in
#     --gene-file (all rows, before filter) is used as the universe, which
#     represents the whole measured transcriptome.
#
# OUTPUTS  (<output-dir>/)
#   enrichment_{label}.xlsx        GO and KEGG result sheets (always written)
#   enrichment_{label}.pdf         dotplot (GO) + barplot (KEGG); only if results
#   enrichment_{label}_empty.txt   marker file when no enrichment found
#
# Usage:
#   Rscript run_pathway_enrichment.R \
#     --gene-file      results/deviant_gene_drivers_ct/gene_drivers_Q1_10pct.tsv \
#     --filter-col     in_M2_topN \
#     --filter-val     TRUE \
#     --output-dir     results/enrichment_ct \
#     --label          Q1_10pct_M2 \
#     --kegg-cache     data/cache/kegg_dme.rds \
#     --condition-label "Ctrl Q1 M2 top-N"
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
  library(clusterProfiler)
  library(org.Dm.eg.db)
  library(openxlsx)
  library(grid)
  library(stringr)
})

# ---- CLI --------------------------------------------------------------------
parser <- ArgumentParser(description = "Generic pathway enrichment for any gene-list TSV")
parser$add_argument("--gene-file",       required = TRUE,
                    help = "TSV with a FlyBase gene-ID column (FBgn...)")
parser$add_argument("--id-col",          default = "gene_id",
                    help = "Column in --gene-file holding FlyBase IDs (default: 'gene_id')")
parser$add_argument("--filter-col",      default = NULL,
                    help = "Column to filter rows on before enrichment (optional)")
parser$add_argument("--filter-val",      default = "TRUE",
                    help = "Value to keep in --filter-col (default: 'TRUE')")
parser$add_argument("--background-file", default = NULL,
                    help = "TSV with background gene IDs (same --id-col). Omit to use all genes in --gene-file as universe")
parser$add_argument("--utils-dir",       default = NULL,
                    help = "Directory containing utils_pathways.R; auto-resolved from script location if NULL")
parser$add_argument("--kegg-cache",      default = NULL,
                    help = "Path to RDS cache of download_KEGG('dme'); created on first run if absent")
parser$add_argument("--output-dir",      required = TRUE,
                    help = "Destination directory (created if absent)")
parser$add_argument("--ontology",        default = "BP",
                    choices = c("BP", "CC", "MF", "ALL"),
                    help = "GO ontology (default: BP)")
parser$add_argument("--simplify-go",     type = "logical", default = TRUE,
                    help = "Simplify redundant GO terms; valid for BP/CC/MF (default: TRUE)")
parser$add_argument("--pval-cutoff",     type = "double", default = 0.05,
                    help = "p-value cutoff for both GO and KEGG (default: 0.05)")
parser$add_argument("--qval-cutoff",     type = "double", default = 0.2,
                    help = "q-value cutoff for both GO and KEGG (default: 0.2)")
parser$add_argument("--label",           default = "",
                    help = "String appended to output filenames (e.g. 'Q1_10pct_M2')")
parser$add_argument("--condition-label", default = "Condition",
                    help = "Label for plot titles (default: 'Condition')")
args <- parser$parse_args()

cat("=== Pathway Enrichment ===\n")
cat("Gene file:        ", args$gene_file,       "\n")
cat("Filter:           ", args$filter_col, "==", args$filter_val, "\n")
cat("Ontology:         ", args$ontology,         "\n")
cat("p / q cutoff:     ", args$pval_cutoff, "/", args$qval_cutoff, "\n")
cat("Output directory: ", args$output_dir,       "\n\n")

if (!file.exists(args$gene_file))
  stop("Gene file not found: ", args$gene_file)

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Resolve output filename stem
out_stem <- file.path(args$output_dir,
                      if (nzchar(args$label)) paste0("enrichment_", args$label)
                      else "enrichment")

# ---- Source utils_pathways.R ------------------------------------------------
resolve_utils_path <- function(utils_dir) {
  if (!is.null(utils_dir)) {
    path <- file.path(utils_dir, "utils_pathways.R")
    if (file.exists(path)) return(path)
    stop("utils_pathways.R not found in --utils-dir: ", utils_dir)
  }
  # Auto-resolve: walk up from this script's directory to src/utils/
  script_dir <- tryCatch(
    dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
      grep("--file=", commandArgs(trailingOnly = FALSE))
    ], mustWork = FALSE)),
    error = function(e) getwd()
  )
  candidates <- c(
    file.path(script_dir, "..", "..", "utils", "utils_pathways.R"),
    file.path(script_dir, "..", "utils", "utils_pathways.R"),
    file.path(getwd(), "src", "utils", "utils_pathways.R")
  )
  found <- candidates[file.exists(normalizePath(candidates, mustWork = FALSE))]
  if (length(found) == 0L)
    stop("Cannot auto-resolve utils_pathways.R. Pass --utils-dir explicitly.")
  normalizePath(found[1L])
}

utils_path <- resolve_utils_path(args$utils_dir)
cat("Sourcing utils:", utils_path, "\n")
source(utils_path)

# ---- KEGG cache helper ------------------------------------------------------
load_or_download_kegg <- function(cache_path) {
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cat("  Loading KEGG cache from:", cache_path, "\n")
    return(readRDS(cache_path))
  }
  cat("  Downloading KEGG pathway data for dme...\n")
  kegg_sets <- clusterProfiler::download_KEGG("dme")
  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(kegg_sets, cache_path)
    cat("  KEGG cache saved to:", cache_path, "\n")
  }
  kegg_sets
}

# ---- Load and filter gene list ----------------------------------------------
gene_dt_full <- fread(args$gene_file, data.table = TRUE)

if (!args$id_col %in% colnames(gene_dt_full))
  stop("Column '", args$id_col, "' not found in --gene-file. ",
       "Available columns: ", paste(colnames(gene_dt_full), collapse = ", "))

# Background = all genes in the full file (full measured transcriptome)
bg_ids <- unique(gene_dt_full[[args$id_col]])
bg_ids <- bg_ids[!is.na(bg_ids) & nzchar(bg_ids)]

# Override background from file if provided
if (!is.null(args$background_file)) {
  if (!file.exists(args$background_file))
    stop("Background file not found: ", args$background_file)
  bg_dt  <- fread(args$background_file, data.table = TRUE)
  bg_ids <- unique(bg_dt[[args$id_col]])
  bg_ids <- bg_ids[!is.na(bg_ids) & nzchar(bg_ids)]
}

# Filter query set
if (!is.null(args$filter_col)) {
  filter_col_present <- args$filter_col %in% colnames(gene_dt_full)
  if (!filter_col_present)
    stop("--filter-col '", args$filter_col, "' not found in gene file.")
  query_dt <- gene_dt_full[as.character(get(args$filter_col)) == args$filter_val]
} else {
  query_dt <- gene_dt_full
}
query_ids <- unique(query_dt[[args$id_col]])
query_ids <- query_ids[!is.na(query_ids) & nzchar(query_ids)]

cat(sprintf("Query genes: %d  |  Background genes: %d\n",
            length(query_ids), length(bg_ids)))

# ---- Gene ID conversion: FlyBase -> ENTREZID --------------------------------
cat("Converting FlyBase -> ENTREZID...\n")
query_conv <- tryCatch(
  bitr(query_ids, fromType = "FLYBASE", toType = "ENTREZID", OrgDb = org.Dm.eg.db,
       drop = TRUE),
  error = function(e) { message("  bitr error: ", e$message); data.frame() }
)
bg_conv <- tryCatch(
  bitr(bg_ids, fromType = "FLYBASE", toType = "ENTREZID", OrgDb = org.Dm.eg.db,
       drop = TRUE),
  error = function(e) { message("  bitr error (bg): ", e$message); data.frame() }
)

query_entrez <- if (nrow(query_conv) > 0L) query_conv$ENTREZID else character(0L)
bg_entrez    <- if (nrow(bg_conv)    > 0L) bg_conv$ENTREZID    else NULL

cat(sprintf("  Query: %d input -> %d ENTREZID (%.1f%%)\n",
            length(query_ids), length(query_entrez),
            100 * length(query_entrez) / max(length(query_ids), 1L)))
cat(sprintf("  Background: %d input -> %d ENTREZID\n",
            length(bg_ids), length(bg_entrez)))

# ---- KEGG ID conversion: ENTREZID -> KEGG -----------------------------------
cat("Converting ENTREZID -> KEGG IDs...\n")
kegg_sets <- load_or_download_kegg(args$kegg_cache)

query_kegg_conv <- tryCatch(
  bitr_kegg(query_entrez, fromType = "ncbi-geneid", toType = "kegg",
            organism = "dme"),
  error = function(e) { message("  bitr_kegg error: ", e$message); data.frame() }
)
bg_kegg_conv <- if (!is.null(bg_entrez)) {
  tryCatch(
    bitr_kegg(bg_entrez, fromType = "ncbi-geneid", toType = "kegg",
              organism = "dme"),
    error = function(e) { message("  bitr_kegg error (bg): ", e$message); data.frame() }
  )
} else {
  data.frame()
}

query_kegg <- if (nrow(query_kegg_conv) > 0L) query_kegg_conv$kegg else character(0L)
bg_kegg    <- if (nrow(bg_kegg_conv)    > 0L) bg_kegg_conv$kegg    else NULL

cat(sprintf("  KEGG: %d ENTREZID -> %d KEGG IDs\n",
            length(query_entrez), length(query_kegg)))

# ---- Early exit if no query genes converted ---------------------------------
if (length(query_entrez) < 2L) {
  cat("  Too few query genes converted (<2). Writing empty marker.\n")
  con <- file(paste0(out_stem, "_empty.txt"), open = "w")
  writeLines(c(as.character(Sys.time()),
               sprintf("Input: %d genes, converted: %d ENTREZID",
                       length(query_ids), length(query_entrez)),
               "No enrichment run (too few genes after ID conversion)"), con)
  close(con)
  openxlsx::write.xlsx(
    list(GO = data.frame(message = "No enrichment run: too few genes after ID conversion"),
         KEGG = data.frame(message = "No enrichment run: too few genes after ID conversion")),
    file = paste0(out_stem, ".xlsx")
  )
  quit(save = "no", status = 0L)
}

# ---- GO enrichment ----------------------------------------------------------
cat(sprintf("Running GO enrichment (ontology=%s, p<%.2g, q<%.2g)...\n",
            args$ontology, args$pval_cutoff, args$qval_cutoff))

# Use EnrichGo from utils_pathways.R but override cutoffs via the underlying call
# because EnrichGo hardcodes p<0.05 / q<0.2; for full flexibility we call
# enrichGO directly here while keeping simplify logic from utils.
run_go <- function(entrez_ids, bg, ont, simplify_flag, pval, qval) {
  result <- tryCatch(
    enrichGO(
      gene          = entrez_ids,
      universe      = bg,
      OrgDb         = org.Dm.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pvalueCutoff  = pval,
      qvalueCutoff  = qval,
      pAdjustMethod = "fdr",
      minGSSize     = 10L,
      maxGSSize     = 500L,
      readable      = TRUE,
      pool          = FALSE
    ),
    error = function(e) { message("  enrichGO error: ", e$message); NULL }
  )
  no_results <- is.null(result) || nrow(result) == 0L
  if (!no_results && simplify_flag && ont %in% c("BP", "CC", "MF")) {
    result <- tryCatch(
      clusterProfiler::simplify(result, cutoff = 0.7, by = "p.adjust",
                                 select_fun = min),
      error = function(e) result
    )
  }
  result
}

go_result <- run_go(query_entrez, bg_entrez, args$ontology,
                    args$simplify_go, args$pval_cutoff, args$qval_cutoff)
go_n      <- if (!is.null(go_result)) nrow(go_result) else 0L
cat(sprintf("  GO: %d terms enriched\n", go_n))

# ---- KEGG enrichment (using cached gene sets) --------------------------------
cat(sprintf("Running KEGG enrichment (p<%.2g, q<%.2g)...\n",
            args$pval_cutoff, args$qval_cutoff))

kegg_result <- if (length(query_kegg) >= 2L) {
  tryCatch(
    enricher(
      gene          = query_kegg,
      universe      = bg_kegg,
      TERM2GENE     = kegg_sets$KEGGPATHID2EXTID,
      TERM2NAME     = kegg_sets$KEGGPATHID2NAME,
      pvalueCutoff  = args$pval_cutoff,
      qvalueCutoff  = args$qval_cutoff,
      pAdjustMethod = "fdr",
      minGSSize     = 10L,
      maxGSSize     = 500L
    ),
    error = function(e) { message("  enricher KEGG error: ", e$message); NULL }
  )
} else {
  cat("  Too few KEGG IDs (<2) — skipping KEGG\n")
  NULL
}
kegg_n <- if (!is.null(kegg_result)) nrow(kegg_result) else 0L
cat(sprintf("  KEGG: %d pathways enriched\n", kegg_n))

# ---- Build output data frames for xlsx ----------------------------------------
go_sheet <- if (go_n > 0L) {
  as.data.frame(go_result)
} else {
  data.frame(message = sprintf(
    "No GO terms enriched (n_query=%d, n_converted=%d, ontology=%s, p<%.2g)",
    length(query_ids), length(query_entrez), args$ontology, args$pval_cutoff
  ))
}

kegg_sheet <- if (kegg_n > 0L) {
  df <- as.data.frame(kegg_result)
  # Strip organism suffix from KEGG pathway names for readability
  df$Description <- gsub(" - Drosophila melanogaster \\(fruit fly\\)", "",
                          df$Description)
  df
} else {
  data.frame(message = sprintf(
    "No KEGG pathways enriched (n_query=%d, n_kegg_ids=%d, p<%.2g)",
    length(query_ids), length(query_kegg), args$pval_cutoff
  ))
}

# ---- Write xlsx -------------------------------------------------------------
out_xlsx <- paste0(out_stem, ".xlsx")
sheet_names <- list()
sheet_names[[paste0("GO_", args$ontology)]] <- go_sheet
sheet_names[["KEGG"]] <- kegg_sheet
openxlsx::write.xlsx(sheet_names, file = out_xlsx)
cat("  xlsx ->", out_xlsx, "\n")

# ---- Write PDF (only when ≥1 non-empty result) --------------------------------
both_empty <- (go_n == 0L && kegg_n == 0L)

if (!both_empty) {
  out_pdf <- paste0(out_stem, ".pdf")
  pdf(out_pdf, width = 10, height = 7)

  if (go_n > 0L) {
    go_plot <- enrichplot::dotplot(go_result, showCategory = 20) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
      labs(title = sprintf("GO %s enrichment  —  %s  |  %s",
                           args$ontology, args$label, args$condition_label),
           subtitle = sprintf("n_query=%d  n_ENTREZID=%d  terms=%d",
                              length(query_ids), length(query_entrez), go_n))
    print(go_plot)
  }

  if (kegg_n > 0L) {
    kegg_result_strip <- kegg_result
    kegg_result_strip@result$Description <- gsub(
      " - Drosophila melanogaster \\(fruit fly\\)", "",
      kegg_result_strip@result$Description
    )
    kegg_plot <- barplot(kegg_result_strip, showCategory = 20) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
      labs(title = sprintf("KEGG enrichment  —  %s  |  %s",
                           args$label, args$condition_label),
           subtitle = sprintf("n_query=%d  n_KEGG=%d  pathways=%d",
                              length(query_ids), length(query_kegg), kegg_n))
    print(kegg_plot)
  }

  invisible(dev.off())
  cat("  pdf  ->", out_pdf, "\n")
}

# ---- Write empty marker when no enrichment found ----------------------------
if (both_empty) {
  out_empty <- paste0(out_stem, "_empty.txt")
  con <- file(out_empty, open = "w")
  writeLines(c(
    as.character(Sys.time()),
    sprintf("label: %s", args$label),
    sprintf("ontology: %s  p<%.2g  q<%.2g", args$ontology,
            args$pval_cutoff, args$qval_cutoff),
    sprintf("n_input=%d  n_ENTREZID=%d  n_KEGG=%d",
            length(query_ids), length(query_entrez), length(query_kegg)),
    "No pathways enriched in either GO or KEGG"
  ), con)
  close(con)
  cat("  empty marker ->", out_empty, "\n")
}

cat("\n=== Done ===\n")
