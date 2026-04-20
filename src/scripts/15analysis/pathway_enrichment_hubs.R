#!/usr/bin/env Rscript
# ==============================================================================
# Pathway Enrichment Analysis for Rewiring Hub Genes
#
# Selects top/bottom N genes from a rewiring metrics table by a sort column
# (e.g., L2L1_rewire, L2L1_conn, L2L1_deg), then performs GO, KEGG, and
# optional GSEA enrichment analysis (Drosophila melanogaster).
#
# Usage:
#   Rscript pathway_enrichment_hubs.R \
#     --input-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
#     --sort-col L2L1_rewire \
#     --n-genes 500 \
#     --direction top \
#     --output-dir results/enrichment \
#     --output-prefix hubs_top500 \
#     --ontology BP \
#     --run-gsea
# ==============================================================================

rm(list = ls())

# ----- 1. Load packages -----

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(argparse)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Dm.eg.db)
  library(openxlsx)
})

# ----- 2. Parse command-line arguments -----

parser <- ArgumentParser(description = "GO/KEGG/GSEA enrichment on top or bottom N hub genes (Drosophila)")

parser$add_argument("--input-file",
  help    = "Path to rewiring hub TSV file (must contain ENTREZID or gene_id FBgn column)",
  required = TRUE)
parser$add_argument("--gene-col",
  help    = "Column with gene IDs. FBgn IDs (e.g., gene_id) or numeric ENTREZID",
  default = "gene_id")
parser$add_argument("--sort-col",
  help    = "Column to sort by for gene selection (e.g., L2L1_rewire, L2L1_conn, L2L1_deg)",
  default = "L2L1_rewire")
parser$add_argument("--n-genes",
  help    = "Number of genes to select",
  type    = "integer",
  default = 500)
parser$add_argument("--direction",
  help    = "Which end to select: 'top', 'bottom', or 'both' (union of top + bottom)",
  default = "top",
  choices = c("top", "bottom", "both"))
parser$add_argument("--output-dir",
  help    = "Output directory (created if missing)",
  default = "results/enrichment")
parser$add_argument("--output-prefix",
  help    = "Prefix for all output file names",
  default = "hubs_enrich")
parser$add_argument("--ontology",
  help    = "GO ontology: 'BP', 'MF', 'CC', or 'ALL' (runs all three)",
  default = "BP",
  choices = c("BP", "MF", "CC", "ALL"))
parser$add_argument("--pvalue-cutoff",
  help    = "Initial p-value cutoff for enrichment",
  type    = "double",
  default = 0.05)
parser$add_argument("--qvalue-cutoff",
  help    = "FDR q-value cutoff for enrichment",
  type    = "double",
  default = 0.2)
parser$add_argument("--simplify-go",
  help    = "Simplify GO results using semantic similarity (removes redundant terms)",
  action  = "store_true",
  default = FALSE)
parser$add_argument("--run-gsea",
  help    = "Run GSEA (gseGO + gseKEGG) using all genes ranked by sort-col",
  action  = "store_true",
  default = FALSE)
parser$add_argument("--min-genesetsize",
  help    = "Minimum gene set size for enrichment",
  type    = "integer",
  default = 10)
parser$add_argument("--max-genesetsize",
  help    = "Maximum gene set size for enrichment",
  type    = "integer",
  default = 500)
parser$add_argument("--universe-file",
  help    = "Optional path to a background gene list file. One gene ID per line, or a
             TSV/CSV with a 'gene_id' column. If omitted, the universe defaults to all
             genes in the input hub file.",
  default = NULL)
parser$add_argument("--show-category",
  help    = "Number of top categories to show in plots",
  type    = "integer",
  default = 20)

args <- parser$parse_args()

# ----- 3. Setup output directory -----

if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}

# Build output file basename (encodes direction, n, and sort column)
sort_col_safe <- gsub("[^A-Za-z0-9_]", "_", args$sort_col)
out_base <- file.path(args$output_dir,
  paste0(args$output_prefix, "_", args$direction, args$n_genes, "_", sort_col_safe))

cat("=== Pathway Enrichment Analysis (Drosophila) ===\n")
cat("Input file  :", args$input_file, "\n")
cat("Sort column :", args$sort_col, "\n")
cat("Direction   :", args$direction, "\n")
cat("N genes     :", args$n_genes, "\n")
cat("Ontology    :", args$ontology, "\n")
cat("Output base :", out_base, "\n\n")

# ----- 4. Load input table -----

cat("Loading input file...\n")
if (!file.exists(args$input_file)) stop("Input file not found: ", args$input_file)
hub_tab <- fread(args$input_file)
cat("  Rows:", nrow(hub_tab), " | Columns:", ncol(hub_tab), "\n")

# Validate required columns
required_cols <- c(args$gene_col, args$sort_col)
missing_cols  <- setdiff(required_cols, colnames(hub_tab))
if (length(missing_cols) > 0) {
  stop("Columns not found in input: ", paste(missing_cols, collapse = ", "),
       "\nAvailable columns: ", paste(colnames(hub_tab), collapse = ", "))
}

# ----- 5. Select genes -----

cat("Selecting genes...\n")

# Drop rows with NA in sort column
n_before <- nrow(hub_tab)
hub_tab  <- hub_tab[!is.na(get(args$sort_col))]
if (nrow(hub_tab) < n_before) {
  cat("  Dropped", n_before - nrow(hub_tab), "rows with NA in", args$sort_col, "\n")
}

# Sort descending by sort column
hub_tab <- hub_tab[order(-get(args$sort_col))]
n_total <- nrow(hub_tab)
n_sel   <- min(args$n_genes, n_total)

selected <- switch(args$direction,
  top    = hub_tab[seq_len(n_sel)],
  bottom = hub_tab[seq(n_total - n_sel + 1L, n_total)],
  both   = unique(rbind(hub_tab[seq_len(n_sel)],
                        hub_tab[seq(n_total - n_sel + 1L, n_total)]))
)

cat("  Total genes in file :", n_total, "\n")
cat("  Selected genes      :", nrow(selected), "\n")
cat("  Sort col range (all): [",
    round(hub_tab[[args$sort_col]][n_total], 4), ",",
    round(hub_tab[[args$sort_col]][1],       4), "]\n")
cat("  Sort col range (sel): [",
    round(min(selected[[args$sort_col]]), 4), ",",
    round(max(selected[[args$sort_col]]), 4), "]\n\n")

# Write selected gene list for traceability
fwrite(selected, paste0(out_base, "_selected_genes.tsv"), sep = "\t")
cat("  Saved selected gene list:", paste0(out_base, "_selected_genes.tsv"), "\n\n")

# ----- 6. Map gene IDs to ENTREZID -----

cat("Mapping gene IDs to Entrez IDs...\n")

map_to_entrez <- function(ids, gene_col) {
  # If the column already contains numeric-style IDs, try using them directly
  if (gene_col == "ENTREZID" || !any(grepl("^FBgn", ids, ignore.case = TRUE))) {
    clean <- ids[!is.na(ids) & ids != "" & ids != "NA"]
    cat("  Using column directly as ENTREZID (", length(clean), "non-NA values)\n")
    return(as.character(clean))
  }
  # FBgn → ENTREZID via org.Dm.eg.db
  cat("  Mapping FBgn IDs via org.Dm.eg.db...\n")
  mapped <- AnnotationDbi::mapIds(org.Dm.eg.db,
    keys     = unique(ids[!is.na(ids)]),
    column   = "ENTREZID",
    keytype  = "FLYBASE",
    multiVals = "first")
  result <- mapped[!is.na(mapped)]
  cat("  Mapped:", length(result), "/", length(unique(ids[!is.na(ids)])), "FBgn IDs\n")
  as.character(result)
}

# If ENTREZID column is already present and populated, use it directly
has_entrezid_col <- "ENTREZID" %in% colnames(hub_tab) &&
                    !all(is.na(hub_tab$ENTREZID))

if (has_entrezid_col && args$gene_col != "ENTREZID") {
  cat("  Found ENTREZID column — using it directly\n")
  query_entrez <- as.character(selected$ENTREZID[!is.na(selected$ENTREZID)])
  hub_entrez   <- as.character(hub_tab$ENTREZID[!is.na(hub_tab$ENTREZID)])
} else {
  gene_ids_sel <- as.character(selected[[args$gene_col]])
  gene_ids_all <- as.character(hub_tab[[args$gene_col]])
  query_entrez <- map_to_entrez(gene_ids_sel, args$gene_col)
  hub_entrez   <- map_to_entrez(gene_ids_all, args$gene_col)
}

# Universe: load from file if provided, otherwise fall back to hub table genes
if (!is.null(args$universe_file)) {
  if (!file.exists(args$universe_file)) {
    stop("Universe file not found: ", args$universe_file)
  }
  cat("  Loading background universe from:", args$universe_file, "\n")
  uni_raw <- fread(args$universe_file, header = TRUE)
  if ("gene_id" %in% colnames(uni_raw)) {
    uni_ids <- as.character(uni_raw[["gene_id"]])
  } else {
    uni_ids <- as.character(uni_raw[[1]])
  }
  uni_ids <- uni_ids[!is.na(uni_ids) & uni_ids != ""]
  cat("  Universe IDs loaded:", length(uni_ids), "\n")
  universe_entrez <- map_to_entrez(uni_ids, "gene_id")
  cat("  Universe genes (Entrez):", length(universe_entrez), "\n")
} else {
  universe_entrez <- hub_entrez
}

cat("  Query genes (Entrez)    :", length(query_entrez), "\n")
cat("  Universe genes (Entrez) :", length(universe_entrez), "\n\n")

if (length(query_entrez) == 0) {
  stop("No query genes could be mapped to Entrez IDs. Check --gene-col and input file.")
}

# ----- 7. Helper: run GO enrichment with cascading cutoffs -----

run_enrichGO <- function(gene_ids, universe_ids, ont,
                         pcut, qcut, min_gs, max_gs, simplify_go) {
  cutoff_tries <- list(
    list(p = pcut, q = qcut,  label = sprintf("p<%.2g, q<%.2g", pcut, qcut))
    # list(p = 0.2,  q = 1.0,   label = "p<0.2,  q<1.0  (relaxed)"),
    # list(p = 1.0,  q = 1.0,   label = "p<1.0,  q<1.0  (no cutoff)")
  )

  result <- NULL
  for (try in cutoff_tries) {
    tryCatch({
      result <- enrichGO(
        gene          = gene_ids,
        universe      = universe_ids,
        OrgDb         = org.Dm.eg.db,
        ont           = ont,
        keyType       = "ENTREZID",
        pvalueCutoff  = try$p,
        qvalueCutoff  = try$q,
        pAdjustMethod = "fdr",
        minGSSize     = min_gs,
        maxGSSize     = max_gs,
        readable      = TRUE)
    }, error = function(e) {
      cat("    Warning:", e$message, "\n")
      result <<- NULL
    })
    if (!is.null(result) && nrow(result) > 0) break
    cat("    No", ont, "results at", try$label, "— trying next...\n")
  }

  if (is.null(result) || nrow(result) == 0) {
    cat("    No", ont, "enrichments found\n")
    return(data.frame())
  }

  cat("    Found", nrow(result), ont, "terms\n")

  if (simplify_go && ont %in% c("BP", "MF", "CC")) {
    tryCatch({
      simplified <- clusterProfiler::simplify(result, cutoff = 0.7,
                                              by = "p.adjust", select_fun = min)
      cat("    Simplified:", nrow(result), "→", nrow(simplified), "terms\n")
      result <- simplified
    }, error = function(e) {
      cat("    Warning: simplify() failed —", e$message, "\n")
    })
  }

  result
}

# ----- 8. Run GO enrichment -----

ont_list <- if (args$ontology == "ALL") c("BP", "MF", "CC") else args$ontology
go_results <- list()

for (ont in ont_list) {
  cat("Running GO enrichment (", ont, ")...\n", sep = "")
  go_results[[ont]] <- run_enrichGO(
    gene_ids    = query_entrez,
    universe_ids = universe_entrez,
    ont         = ont,
    pcut        = args$pvalue_cutoff,
    qcut        = args$qvalue_cutoff,
    min_gs      = args$min_genesetsize,
    max_gs      = args$max_genesetsize,
    simplify_go = args$simplify_go
  )
}

# ----- 9. Run KEGG enrichment -----

cat("Running KEGG enrichment...\n")

kegg_cutoffs <- list(
  list(p = args$pvalue_cutoff, q = args$qvalue_cutoff,
       label = sprintf("p<%.2g, q<%.2g", args$pvalue_cutoff, args$qvalue_cutoff)),
  list(p = 0.2, q = 1.0, label = "p<0.2, q<1.0 (relaxed)"),
  list(p = 1.0, q = 1.0, label = "p<1.0, q<1.0 (no cutoff)")
)

kegg_result <- NULL
for (try in kegg_cutoffs) {
  tryCatch({
    kegg_result <- enrichKEGG(
      gene          = query_entrez,
      universe      = universe_entrez,
      organism      = "dme",
      keyType       = "ncbi-geneid",
      pvalueCutoff  = try$p,
      qvalueCutoff  = try$q,
      pAdjustMethod = "fdr",
      minGSSize     = args$min_genesetsize,
      maxGSSize     = args$max_genesetsize)
  }, error = function(e) {
    cat("  Warning:", e$message, "\n")
    kegg_result <<- NULL
  })
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) break
  cat("  No KEGG results at", try$label, "— trying next...\n")
}

if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
  # Make readable: convert Entrez IDs in gene lists to symbols
  tryCatch({
    kegg_result <- setReadable(kegg_result, OrgDb = org.Dm.eg.db, keyType = "ENTREZID")
  }, error = function(e) {
    cat("  Warning: setReadable() failed for KEGG —", e$message, "\n")
  })
  # Strip organism suffix from pathway names
  kegg_result@result$Description <- gsub(
    " - Drosophila melanogaster \\(fruit fly\\)", "",
    kegg_result@result$Description)
  cat("  Found", nrow(kegg_result), "KEGG pathways\n")
} else {
  cat("  No KEGG enrichments found",
      "(common for Drosophila — limited KEGG pathway coverage)\n")
  kegg_result <- data.frame()
}
cat("\n")

# ----- 10. Optional GSEA -----

gsea_go_results  <- list()
gsea_kegg_result <- NULL

if (args$run_gsea) {
  cat("Running GSEA...\n")

  # Build ranked list: all background genes ordered by sort_col
  # Names must be ENTREZID; values are the sort metric
  entrez_for_gsea <- if (has_entrezid_col && args$gene_col != "ENTREZID") {
    hub_tab$ENTREZID
  } else {
    hub_tab[[args$gene_col]]
  }

  gsea_tab <- hub_tab[!is.na(entrez_for_gsea) & !is.na(get(args$sort_col))]
  if (has_entrezid_col && args$gene_col != "ENTREZID") {
    gsea_ids <- as.character(gsea_tab$ENTREZID)
  } else {
    gsea_ids <- as.character(gsea_tab[[args$gene_col]])
  }

  ranked_vec <- setNames(as.numeric(gsea_tab[[args$sort_col]]), gsea_ids)
  ranked_vec <- ranked_vec[!is.na(names(ranked_vec)) & names(ranked_vec) != ""]
  # Remove duplicated IDs, keeping first occurrence (highest value since sorted)
  ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]
  ranked_vec <- sort(ranked_vec, decreasing = TRUE)

  cat("  Ranked gene list size:", length(ranked_vec), "\n")

  for (ont in ont_list) {
    cat("  gseGO (", ont, ")...\n", sep = "")
    tryCatch({
      gsea_go_results[[ont]] <- gseGO(
        geneList      = ranked_vec,
        OrgDb         = org.Dm.eg.db,
        ont           = ont,
        keyType       = "ENTREZID",
        minGSSize     = args$min_genesetsize,
        maxGSSize     = args$max_genesetsize,
        pvalueCutoff  = args$pvalue_cutoff,
        pAdjustMethod = "fdr",
        verbose       = FALSE)
      cat("    Found", nrow(gsea_go_results[[ont]]), "terms\n")
    }, error = function(e) {
      cat("    Warning: gseGO failed —", e$message, "\n")
      gsea_go_results[[ont]] <<- data.frame()
    })
  }

  cat("  gseKEGG...\n")
  tryCatch({
    gsea_kegg_result <- gseKEGG(
      geneList      = ranked_vec,
      organism      = "dme",
      keyType       = "ncbi-geneid",
      minGSSize     = args$min_genesetsize,
      maxGSSize     = args$max_genesetsize,
      pvalueCutoff  = args$pvalue_cutoff,
      pAdjustMethod = "fdr",
      verbose       = FALSE)
    if (!is.null(gsea_kegg_result) && nrow(gsea_kegg_result) > 0) {
      tryCatch(
        gsea_kegg_result <- setReadable(gsea_kegg_result,
                                        OrgDb = org.Dm.eg.db, keyType = "ENTREZID"),
        error = function(e) invisible(NULL)
      )
      gsea_kegg_result@result$Description <- gsub(
        " - Drosophila melanogaster \\(fruit fly\\)", "",
        gsea_kegg_result@result$Description)
    }
    cat("    Found", nrow(gsea_kegg_result), "pathways\n")
  }, error = function(e) {
    cat("    Warning: gseKEGG failed —", e$message, "\n")
    gsea_kegg_result <<- data.frame()
  })
  cat("\n")
}

# ----- 11. Output: Excel workbook -----

cat("Writing Excel output...\n")

result_to_df <- function(res) {
  if (is.data.frame(res) || nrow(res) == 0) return(data.frame())
  tryCatch(as.data.frame(res), error = function(e) data.frame())
}

sheets <- list()
for (ont in ont_list) {
  sheets[[paste0("GO_", ont)]] <- result_to_df(go_results[[ont]])
}
sheets[["KEGG"]] <- result_to_df(kegg_result)

xlsx_file <- paste0(out_base, ".xlsx")
tryCatch({
  openxlsx::write.xlsx(sheets, file = xlsx_file)
  cat("  Saved:", xlsx_file, "\n")
}, error = function(e) {
  cat("  Warning: could not write Excel file —", e$message, "\n")
})

if (args$run_gsea) {
  gsea_sheets <- list()
  for (ont in ont_list) {
    gsea_sheets[[paste0("GSEA_GO_", ont)]] <- result_to_df(gsea_go_results[[ont]])
  }
  gsea_sheets[["GSEA_KEGG"]] <- result_to_df(gsea_kegg_result)
  gsea_xlsx <- paste0(out_base, "_gsea.xlsx")
  tryCatch({
    openxlsx::write.xlsx(gsea_sheets, file = gsea_xlsx)
    cat("  Saved:", gsea_xlsx, "\n")
  }, error = function(e) {
    cat("  Warning: could not write GSEA Excel file —", e$message, "\n")
  })
}
cat("\n")

# ----- 12. Output: PDF plots -----

has_result <- function(res) {
  !is.data.frame(res) && !is.null(res) && nrow(res) > 0
}

make_dotplot <- function(res, title, n_show) {
  tryCatch(
    enrichplot::dotplot(res, showCategory = n_show, title = title) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      theme(plot.title = element_text(size = 12, face = "bold")),
    error = function(e) NULL
  )
}

make_barplot <- function(res, title, n_show) {
  tryCatch(
    barplot(res, showCategory = n_show, title = title) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
      theme(plot.title = element_text(size = 12, face = "bold")),
    error = function(e) NULL
  )
}

cat("Writing PDF plots...\n")
pdf_file <- paste0(out_base, "_plots.pdf")

plots_go   <- Filter(Negate(is.null),
  lapply(ont_list, function(ont) {
    if (has_result(go_results[[ont]])) {
      make_dotplot(go_results[[ont]],
        title  = paste0("GO ", ont, " enrichment (", args$direction, " ",
                        nrow(selected), " by ", args$sort_col, ")"),
        n_show = args$show_category)
    }
  })
)

plots_kegg <- list()
if (has_result(kegg_result)) {
  p <- make_barplot(kegg_result,
    title  = paste0("KEGG enrichment (", args$direction, " ",
                    nrow(selected), " by ", args$sort_col, ")"),
    n_show = args$show_category)
  if (!is.null(p)) plots_kegg[["KEGG"]] <- p
}

all_plots <- c(plots_go, plots_kegg)

if (length(all_plots) > 0) {
  tryCatch({
    pdf(pdf_file, width = 12, height = 9)
    for (p in all_plots) suppressMessages(print(p))
    dev.off()
    cat("  Saved:", pdf_file, "\n")
  }, error = function(e) {
    dev.off()
    cat("  Warning: PDF plot error —", e$message, "\n")
  })
} else {
  cat("  No enrichment results to plot\n")
}

# GSEA plots
if (args$run_gsea) {
  gsea_pdf <- paste0(out_base, "_gsea_plots.pdf")

  gsea_plots <- list()
  for (ont in ont_list) {
    res <- gsea_go_results[[ont]]
    if (has_result(res) && nrow(res) >= 1) {
      # Ridgeplot: distribution of ranked genes in enriched sets
      tryCatch({
        rp <- enrichplot::ridgeplot(res, showCategory = min(15, nrow(res))) +
          labs(title = paste0("GSEA GO ", ont)) +
          theme(plot.title = element_text(size = 12, face = "bold"))
        gsea_plots[[paste0("ridge_", ont)]] <- rp
      }, error = function(e) invisible(NULL))
      # GSEA plot for the top term
      tryCatch({
        sp <- enrichplot::gseaplot2(res, geneSetID = 1,
          title = paste0("GSEA GO ", ont, ": ", res@result$Description[1]))
        gsea_plots[[paste0("gsea_", ont)]] <- sp
      }, error = function(e) invisible(NULL))
    }
  }
  if (has_result(gsea_kegg_result) && nrow(gsea_kegg_result) >= 1) {
    tryCatch({
      rp <- enrichplot::ridgeplot(gsea_kegg_result,
                                  showCategory = min(15, nrow(gsea_kegg_result))) +
        labs(title = "GSEA KEGG") +
        theme(plot.title = element_text(size = 12, face = "bold"))
      gsea_plots[["ridge_KEGG"]] <- rp
    }, error = function(e) invisible(NULL))
    tryCatch({
      sp <- enrichplot::gseaplot2(gsea_kegg_result, geneSetID = 1,
        title = paste0("GSEA KEGG: ", gsea_kegg_result@result$Description[1]))
      gsea_plots[["gsea_KEGG"]] <- sp
    }, error = function(e) invisible(NULL))
  }

  if (length(gsea_plots) > 0) {
    tryCatch({
      pdf(gsea_pdf, width = 12, height = 9)
      for (p in gsea_plots) suppressMessages(print(p))
      dev.off()
      cat("  Saved:", gsea_pdf, "\n")
    }, error = function(e) {
      dev.off()
      cat("  Warning: GSEA PDF error —", e$message, "\n")
    })
  } else {
    cat("  No GSEA results to plot\n")
  }
}

# ----- 13. Summary -----

cat("\n=== Summary ===\n")
for (ont in ont_list) {
  n <- if (is.data.frame(go_results[[ont]])) nrow(go_results[[ont]])
       else if (has_result(go_results[[ont]])) nrow(go_results[[ont]])
       else 0
  cat(sprintf("  GO %-3s enriched terms : %d\n", ont, n))
}
n_kegg <- if (is.data.frame(kegg_result)) nrow(kegg_result)
          else if (has_result(kegg_result)) nrow(kegg_result)
          else 0
cat(sprintf("  KEGG pathways         : %d\n", n_kegg))
if (args$run_gsea) {
  for (ont in ont_list) {
    n <- tryCatch(nrow(gsea_go_results[[ont]]), error = function(e) 0)
    cat(sprintf("  GSEA GO %-3s terms    : %d\n", ont, n))
  }
  n_gkegg <- tryCatch(nrow(gsea_kegg_result), error = function(e) 0)
  cat(sprintf("  GSEA KEGG pathways    : %d\n", n_gkegg))
}
cat("\nDone.\n")
