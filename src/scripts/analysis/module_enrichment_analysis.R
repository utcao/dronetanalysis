#!/usr/bin/env Rscript
# ==============================================================================
# Module Enrichment Analysis
#
# Script by Gabriel Thornes
#
# Last Updated: 12/11/2025
#
# This script:
#   1. Loads gene metrics and module assignments
#   2. Performs functional enrichment analysis (GO, KEGG) for each module
#   3. Identifies enriched biological processes, pathways, and functions
#   4. Compares enrichment patterns between modules
#   5. Generates comprehensive enrichment reports and visualizations
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(argparse)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Dm.eg.db)
})

# ----- 2. Command-line arguments -----
parser <- ArgumentParser(description = 'Perform enrichment analysis on WGCNA modules (Drosophila)')
parser$add_argument('--gene-metrics-file', help = 'Path to gene metrics CSV file with module assignments', 
                   default = 'results/network_features/gene_metrics/adjacency_gene_metrics_with_expression.csv')
parser$add_argument('--output-dir', help = 'Directory to save enrichment results', 
                   default = 'results/analysis/enrichment')
parser$add_argument('--pvalue-cutoff', help = 'P-value cutoff for enrichment significance', 
                   default = 0.05, type = 'double')
parser$add_argument('--qvalue-cutoff', help = 'Q-value (adjusted p-value) cutoff', 
                   default = 0.2, type = 'double')
parser$add_argument('--min-genesetsize', help = 'Minimum gene set size for enrichment', 
                   default = 2, type = 'integer')
parser$add_argument('--max-genesetsize', help = 'Maximum gene set size for enrichment', 
                   default = 500, type = 'integer')
args <- parser$parse_args()

output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Module Enrichment Analysis (Drosophila) ===\n")
cat("Gene metrics file:", args$gene_metrics_file, "\n")
cat("P-value cutoff:", args$pvalue_cutoff, "\n")
cat("Q-value cutoff:", args$qvalue_cutoff, "\n")
cat("Output directory:", output_dir, "\n\n")

# ----- 3. Set up organism database -----
cat("Setting up Drosophila melanogaster annotation database...\n")
org_db <- org.Dm.eg.db

# ----- 4. Load gene metrics -----
cat("Loading gene metrics...\n")
gene_metrics <- read.csv(args$gene_metrics_file, stringsAsFactors = FALSE)

cat("  Total genes:", nrow(gene_metrics), "\n")
cat("  Columns:", paste(colnames(gene_metrics), collapse = ", "), "\n\n")

# ----- 5. Get unique modules -----
modules <- unique(gene_metrics$module)
modules <- modules[!is.na(modules)]
modules <- sort(modules)

cat("Found", length(modules), "modules:\n")
cat("  ", paste(modules, collapse = ", "), "\n\n")

# ----- 6. Convert gene names to Entrez IDs -----
cat("Converting gene names to Entrez IDs...\n")
cat("  Sample gene names from data:\n")
cat("    ", paste(head(gene_metrics$gene, 5), collapse = ", "), "\n\n")

# Check what keytypes are available
available_keys <- keytypes(org_db)
cat("  Available key types in org.Dm.eg.db:\n")
cat("    ", paste(available_keys, collapse = ", "), "\n\n")

# Try to detect gene name format
# First, check if genes look like FlyBase IDs (e.g., FBgn0000001)
is_flybase <- grepl("^FBgn", gene_metrics$gene[1])
# Check if they look like standard symbols
is_symbol <- gene_metrics$gene[1] %in% keys(org_db, keytype = "SYMBOL")

if (is_flybase) {
  cat("  Detected gene name format: FlyBase IDs (FBgn...)\n")
  keytype_to_use <- "FLYBASE"
} else if (is_symbol) {
  cat("  Detected gene name format: Gene symbols\n")
  keytype_to_use <- "SYMBOL"
} else {
  # Try to find the best matching keytype
  cat("  Attempting to match gene names to available key types...\n")
  for (kt in available_keys) {
    n_matches <- sum(gene_metrics$gene %in% keys(org_db, keytype = kt))
    cat("    ", kt, ":", n_matches, "matches\n")
    if (n_matches > length(gene_metrics$gene) * 0.5) {
      keytype_to_use <- kt
      cat("  Using keytype:", kt, "\n\n")
      break
    }
  }
}

# Create gene name to Entrez ID mapping
gene2entrez <- mapIds(org_db,
                      keys = gene_metrics$gene,
                      column = "ENTREZID",
                      keytype = keytype_to_use,
                      multiVals = "first")

cat("  Successfully mapped:", sum(!is.na(gene2entrez)), "/", length(gene2entrez), "genes\n\n")

# ----- 7. Perform enrichment for each module -----
cat("Performing enrichment analysis for each module...\n\n")

enrichment_results <- list()
all_enrichments <- list()

for (module in modules) {
  cat("Processing module:", module, "\n")
  
  # Get genes in this module
  module_genes <- gene_metrics %>%
    filter(module == module) %>%
    pull(gene)
  
  cat("  Genes in module:", length(module_genes), "\n")
  
  # Get Entrez IDs for module genes
  module_entrez <- gene2entrez[names(gene2entrez) %in% module_genes]
  module_entrez <- module_entrez[!is.na(module_entrez)]
  
  cat("  Mapped to Entrez IDs:", length(module_entrez), "\n")
  
  if (length(module_entrez) == 0) {
    cat("  WARNING: No genes could be mapped for this module\n\n")
    next
  }
  
  cat("  Sample Entrez IDs:", paste(head(module_entrez, 3), collapse = ", "), "\n\n")
  
  # ----- GO Enrichment (Biological Process) -----
  cat("  Running GO enrichment (BP)...\n")
  tryCatch({
    go_bp <- enrichGO(gene = as.character(module_entrez),
                      OrgDb = org_db,
                      ont = "BP",
                      pvalueCutoff = args$pvalue_cutoff,
                      qvalueCutoff = args$qvalue_cutoff,
                      minGSSize = args$min_genesetsize,
                      maxGSSize = args$max_genesetsize,
                      readable = TRUE)
  }, error = function(e) {
    cat("    WARNING: GO BP enrichment failed -", e$message, "\n")
    go_bp <<- NULL
  })
  
  if (is.null(go_bp)) {
    go_bp <- data.frame()
  }
  
  # ----- GO Enrichment (Molecular Function) -----
  cat("  Running GO enrichment (MF)...\n")
  tryCatch({
    go_mf <- enrichGO(gene = as.character(module_entrez),
                      OrgDb = org_db,
                      ont = "MF",
                      pvalueCutoff = args$pvalue_cutoff,
                      qvalueCutoff = args$qvalue_cutoff,
                      minGSSize = args$min_genesetsize,
                      maxGSSize = args$max_genesetsize,
                      readable = TRUE)
  }, error = function(e) {
    cat("    WARNING: GO MF enrichment failed -", e$message, "\n")
    go_mf <<- NULL
  })
  
  if (is.null(go_mf)) {
    go_mf <- data.frame()
  }
  
  # ----- GO Enrichment (Cellular Component) -----
  cat("  Running GO enrichment (CC)...\n")
  tryCatch({
    go_cc <- enrichGO(gene = as.character(module_entrez),
                      OrgDb = org_db,
                      ont = "CC",
                      pvalueCutoff = args$pvalue_cutoff,
                      qvalueCutoff = args$qvalue_cutoff,
                      minGSSize = args$min_genesetsize,
                      maxGSSize = args$max_genesetsize,
                      readable = TRUE)
  }, error = function(e) {
    cat("    WARNING: GO CC enrichment failed -", e$message, "\n")
    go_cc <<- NULL
  })
  
  if (is.null(go_cc)) {
    go_cc <- data.frame()
  }
  
  # Store results
  enrichment_results[[module]] <- list(
    genes = module_genes,
    n_genes = length(module_genes),
    n_mapped = length(module_entrez),
    go_bp = go_bp,
    go_mf = go_mf,
    go_cc = go_cc
  )
  
  # Combine all GO results for summary
  if (nrow(go_bp) > 0) {
    go_bp_df <- as.data.frame(go_bp) %>% mutate(module = module, ontology = "BP")
    all_enrichments[[paste(module, "GO_BP", sep = "_")]] <- go_bp_df
  }
  
  if (nrow(go_mf) > 0) {
    go_mf_df <- as.data.frame(go_mf) %>% mutate(module = module, ontology = "MF")
    all_enrichments[[paste(module, "GO_MF", sep = "_")]] <- go_mf_df
  }
  
  if (nrow(go_cc) > 0) {
    go_cc_df <- as.data.frame(go_cc) %>% mutate(module = module, ontology = "CC")
    all_enrichments[[paste(module, "GO_CC", sep = "_")]] <- go_cc_df
  }
  
  cat("  Enrichments found:\n")
  cat("    GO BP:", nrow(go_bp), "\n")
  cat("    GO MF:", nrow(go_mf), "\n")
  cat("    GO CC:", nrow(go_cc), "\n")
  
  # Try to show actual p-value ranges if results exist
  if (!is.null(go_bp) && nrow(go_bp) > 0) {
    bp_df <- as.data.frame(go_bp)
    cat("    GO BP p-value range: [", format(min(bp_df$pvalue, na.rm=T), scientific=T), ", ", 
        format(max(bp_df$pvalue, na.rm=T), scientific=T), "]\n", sep="")
  }
  cat("\n")
}

# ----- 8. Generate summary tables -----
cat("Generating summary tables...\n")

# Combine all enrichments into one table
if (length(all_enrichments) > 0) {
  all_enrich_df <- do.call(rbind, all_enrichments) %>%
    arrange(module, pvalue) %>%
    dplyr::select(module, ontology, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
  
  write.csv(all_enrich_df, 
            file.path(output_dir, "all_module_enrichments.csv"), 
            row.names = FALSE)
  cat("  Saved: all_module_enrichments.csv\n")
} else {
  cat("  No enrichments found across all modules\n")
}

# ----- 9. Generate per-module enrichment reports -----
cat("Generating per-module enrichment reports...\n")

for (module in names(enrichment_results)) {
  results <- enrichment_results[[module]]
  
  # Create report file
  report_file <- file.path(output_dir, paste0("enrichment_", module, ".txt"))
  
  cat("Writing enrichment report for module:", module, "...\n", file = report_file, append = FALSE)
  
  # Write header
  cat("\n=== ENRICHMENT ANALYSIS FOR MODULE:", module, "===\n\n", 
      file = report_file, append = TRUE)
  cat("Genes in module:", results$n_genes, "\n",
      file = report_file, append = TRUE)
  cat("Genes mapped to Entrez IDs:", results$n_mapped, "\n\n",
      file = report_file, append = TRUE)
  
  # Write GO BP results
  cat("--- BIOLOGICAL PROCESS (GO BP) ---\n",
      file = report_file, append = TRUE)
  if (!is.null(results$go_bp) && is.data.frame(results$go_bp) && nrow(results$go_bp) > 0) {
    go_bp_df <- results$go_bp %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust) %>%
      arrange(pvalue)
    write.table(go_bp_df, file = report_file, append = TRUE, quote = FALSE, 
                sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    cat("No significant enrichments\n", file = report_file, append = TRUE)
  }
  
  cat("\n--- MOLECULAR FUNCTION (GO MF) ---\n",
      file = report_file, append = TRUE)
  if (!is.null(results$go_mf) && is.data.frame(results$go_mf) && nrow(results$go_mf) > 0) {
    go_mf_df <- results$go_mf %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust) %>%
      arrange(pvalue)
    write.table(go_mf_df, file = report_file, append = TRUE, quote = FALSE,
                sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    cat("No significant enrichments\n", file = report_file, append = TRUE)
  }
  
  cat("\n--- CELLULAR COMPONENT (GO CC) ---\n",
      file = report_file, append = TRUE)
  if (!is.null(results$go_cc) && is.data.frame(results$go_cc) && nrow(results$go_cc) > 0) {
    go_cc_df <- results$go_cc %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust) %>%
      arrange(pvalue)
    write.table(go_cc_df, file = report_file, append = TRUE, quote = FALSE,
                sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    cat("No significant enrichments\n", file = report_file, append = TRUE)
  }
}

# ----- 10. Generate enrichment plots -----
cat("Generating enrichment plots...\n")

pdf(file.path(output_dir, "enrichment_plots.pdf"), width = 14, height = 10)

for (module in names(enrichment_results)) {
  results <- enrichment_results[[module]]
  
  plots_to_draw <- list()
  
  # GO BP dotplot
  if (nrow(results$go_bp) > 0 && nrow(results$go_bp) <= 30) {
    tryCatch({
      p_bp <- dotplot(results$go_bp, showCategory = min(15, nrow(results$go_bp)), 
                      title = paste("Module", module, "- GO Biological Process"))
      plots_to_draw[[length(plots_to_draw) + 1]] <- p_bp
    }, error = function(e) {
      cat("    Warning: Could not create GO BP dotplot -", e$message, "\n")
    })
  }
  
  # GO MF dotplot
  if (nrow(results$go_mf) > 0 && nrow(results$go_mf) <= 30) {
    tryCatch({
      p_mf <- dotplot(results$go_mf, showCategory = min(15, nrow(results$go_mf)),
                      title = paste("Module", module, "- GO Molecular Function"))
      plots_to_draw[[length(plots_to_draw) + 1]] <- p_mf
    }, error = function(e) {
      cat("    Warning: Could not create GO MF dotplot -", e$message, "\n")
    })
  }
  
  # GO CC dotplot
  if (nrow(results$go_cc) > 0 && nrow(results$go_cc) <= 30) {
    tryCatch({
      p_cc <- dotplot(results$go_cc, showCategory = min(15, nrow(results$go_cc)),
                      title = paste("Module", module, "- GO Cellular Component"))
      plots_to_draw[[length(plots_to_draw) + 1]] <- p_cc
    }, error = function(e) {
      cat("    Warning: Could not create GO CC dotplot -", e$message, "\n")
    })
  }
  
  # Print plots (one per page or grouped)
  if (length(plots_to_draw) > 0) {
    for (p in plots_to_draw) {
      print(p)
    }
  } else if (nrow(results$go_bp) == 0 && nrow(results$go_mf) == 0 && nrow(results$go_cc) == 0) {
    # Create a simple text plot if no enrichments found
    plot.new()
    text(0.5, 0.5, paste("No enrichments found for module", module), 
         cex = 1.5, ha = "center")
  }
}

dev.off()
cat("  Saved: enrichment_plots.pdf\n")

# ----- 11. Generate summary statistics -----
cat("\n=== ENRICHMENT SUMMARY ===\n")

for (module in names(enrichment_results)) {
  results <- enrichment_results[[module]]
  cat("\nModule:", module, "\n")
  cat("  Genes:", results$n_genes, "(mapped:", results$n_mapped, ")\n")
  cat("  GO BP enrichments:", nrow(results$go_bp), "\n")
  cat("  GO MF enrichments:", nrow(results$go_mf), "\n")
  cat("  GO CC enrichments:", nrow(results$go_cc), "\n")
  
  # Top enrichment
  all_go <- rbind(
    if (nrow(results$go_bp) > 0) as.data.frame(results$go_bp) else NULL,
    if (nrow(results$go_mf) > 0) as.data.frame(results$go_mf) else NULL,
    if (nrow(results$go_cc) > 0) as.data.frame(results$go_cc) else NULL
  )
  
  if (!is.null(all_go) && nrow(all_go) > 0) {
    top_enrich <- all_go %>%
      arrange(pvalue) %>%
      slice_head(n = 3) %>%
      select(Description, pvalue, p.adjust)
    
    cat("  Top enrichments:\n")
    for (i in 1:nrow(top_enrich)) {
      cat("    -", top_enrich$Description[i], 
          "(p =", format(top_enrich$pvalue[i], digits = 3), ")\n")
    }
  }
}

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", output_dir, "\n")
cat("\nOutput files:\n")
cat("  - all_module_enrichments.csv (combined enrichments for all modules)\n")
cat("  - enrichment_[MODULE].txt (detailed reports per module)\n")
cat("  - enrichment_plots.pdf (visualization of top enrichments)\n")
