#!/usr/bin/env Rscript
# ==============================================================================
# Module Enrichment Analysis
#
# Script by Gabriel Thornes
#
# Last Updated: 16/11/2025
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

# Set Python path to conda environment Python for clusterProfiler::simplify() 
Sys.setenv(PYTHON = "/tmp/global2/gthornes/miniforge3/envs/dronetanalysis/bin/python")
options(python_cmd = "/tmp/global2/gthornes/miniforge3/envs/dronetanalysis/bin/python")

# ----- 2. Command-line arguments -----
parser <- ArgumentParser(description = 'Perform enrichment analysis on WGCNA modules (Drosophila)')
parser$add_argument('--gene-metrics-file', help = 'Path to gene metrics CSV file with module assignments',
                   default = 'results/network_features/gene_metrics/adjacency_gene_metrics_with_expression.csv')
parser$add_argument('--output-dir', help = 'Directory to save enrichment results', 
                   default = 'results/analysis/control_enrichment')
parser$add_argument('--pvalue-cutoff', help = 'P-value cutoff for enrichment significance', 
                   default = 0.05, type = 'double')
parser$add_argument('--qvalue-cutoff', help = 'Q-value (adjusted p-value) cutoff', 
                   default = 0.1, type = 'double')
parser$add_argument('--min-genesetsize', help = 'Minimum gene set size for enrichment', 
                   default = 10, type = 'integer')
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

cat("  Successfully mapped:", sum(!is.na(gene2entrez)), "/", length(gene2entrez), "genes to Entrez IDs\n")

# Create gene name to KEGG ID mapping for Drosophila
# KEGG for Drosophila uses "dme:Dmel_" prefix format
cat("\n  Creating KEGG ID mappings for Drosophila...\n")
gene2kegg <- tryCatch({
  mapIds(org_db,
         keys = gene_metrics$gene,
         column = "FLYBASE",  # FlyBase IDs work better for KEGG
         keytype = keytype_to_use,
         multiVals = "first")
}, error = function(e) {
  cat("    Warning: Could not map to FLYBASE, trying alternative approach\n")
  NULL
})

# If FLYBASE mapping worked, convert FlyBase to KEGG format
if (!is.null(gene2kegg)) {
  # Remove "FBgn" prefix and keep only the numeric part for KEGG
  # KEGG Drosophila genes are often just the Entrez ID
  gene2kegg_clean <- gene2entrez  # Use Entrez IDs for KEGG
  cat("  Successfully mapped:", sum(!is.na(gene2kegg_clean)), "/", length(gene2kegg_clean), "genes to KEGG-compatible IDs\n")
} else {
  gene2kegg_clean <- gene2entrez
}

# DIAGNOSTIC: Check a few successful and failed mappings
cat("\n  DIAGNOSTIC - Sample mappings:\n")
sample_genes <- head(gene_metrics$gene, 10)
for (g in sample_genes) {
  entrez <- gene2entrez[names(gene2entrez) == g]
  if (is.na(entrez)) {
    cat("    ", g, "-> FAILED\n")
  } else {
    cat("    ", g, "-> Entrez:", entrez, "\n")
  }
}
cat("\n")

# Create background universe (all genes analyzed in WGCNA)
cat("Creating background universe for enrichment...\n")
universe_entrez <- gene2entrez[!is.na(gene2entrez)]
universe_kegg <- gene2kegg_clean[!is.na(gene2kegg_clean)]
cat("  Universe size (Entrez):", length(universe_entrez), "genes\n")
cat("  Universe size (KEGG):", length(universe_kegg), "genes\n")
cat("  (Using analyzed genes as background instead of full genome)\n\n")

# ----- 7. Perform enrichment for each module -----
cat("Performing enrichment analysis for each module...\n\n")

enrichment_results <- list()
all_enrichments <- list()

for (current_module in modules) {
  cat("Processing module:", current_module, "\n")
  
  # Get genes in this module
  module_genes <- gene_metrics %>%
    filter(module == current_module) %>%
    pull(gene)
  
  cat("  Genes in module:", length(module_genes), "\n")
  
  # Get Entrez IDs for module genes
  module_entrez <- gene2entrez[names(gene2entrez) %in% module_genes]
  module_entrez <- module_entrez[!is.na(module_entrez)]
  
  # Get KEGG IDs for module genes
  module_kegg <- gene2kegg_clean[names(gene2kegg_clean) %in% module_genes]
  module_kegg <- module_kegg[!is.na(module_kegg)]
  
  cat("  Mapped to Entrez IDs:", length(module_entrez), "\n")
  cat("  Mapped to KEGG IDs:", length(module_kegg), "\n")
  
  if (length(module_entrez) == 0) {
    cat("  WARNING: No genes could be mapped for this module\n\n")
    next
  }
  
  cat("  Sample Entrez IDs:", paste(head(module_entrez, 3), collapse = ", "), "\n")
  cat("  Mapping rate:", round(100 * length(module_entrez) / length(module_genes), 1), "%\n\n")
  
  # ----- GO Enrichment (Biological Process) -----
  cat("  Running GO enrichment (BP)...\n")
  tryCatch({
    go_bp <- enrichGO(gene = as.character(module_entrez),
                      universe = as.character(universe_entrez),
                      OrgDb = org_db,
                      ont = "BP",
                      pvalueCutoff = args$pvalue_cutoff,
                      qvalueCutoff = args$qvalue_cutoff,
                      minGSSize = args$min_genesetsize,
                      maxGSSize = args$max_genesetsize,
                      readable = TRUE)
    if (is.null(go_bp) || nrow(go_bp) == 0) {
      cat("    No BP results - trying with relaxed cutoffs (p<0.2, q<1.0)...\n")
      go_bp <- enrichGO(gene = as.character(module_entrez),
                        universe = as.character(universe_entrez),
                        OrgDb = org_db,
                        ont = "BP",
                        pvalueCutoff = 0.2,
                        qvalueCutoff = 1.0,
                        minGSSize = args$min_genesetsize,
                        maxGSSize = args$max_genesetsize,
                        readable = TRUE)
    }
    if (is.null(go_bp) || nrow(go_bp) == 0) {
      cat("    Still no BP results - trying with NO cutoffs (p<1.0, q<1.0)...\n")
      go_bp <- enrichGO(gene = as.character(module_entrez),
                        universe = as.character(universe_entrez),
                        OrgDb = org_db,
                        ont = "BP",
                        pvalueCutoff = 1.0,
                        qvalueCutoff = 1.0,
                        minGSSize = args$min_genesetsize,
                        maxGSSize = args$max_genesetsize,
                        readable = TRUE)
    }
    # Simplify BP results to remove redundancy
    if (!is.null(go_bp) && nrow(go_bp) > 0) {
      go_bp_simplified <- clusterProfiler::simplify(go_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)
      cat("    Simplified from", nrow(go_bp), "to", nrow(go_bp_simplified), "terms\n")
    }
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
                      universe = as.character(universe_entrez),
                      OrgDb = org_db,
                      ont = "MF",
                      pvalueCutoff = args$pvalue_cutoff,
                      qvalueCutoff = args$qvalue_cutoff,
                      minGSSize = args$min_genesetsize,
                      maxGSSize = args$max_genesetsize,
                      readable = TRUE)
    if (is.null(go_mf) || nrow(go_mf) == 0) {
      cat("    No MF results - trying with relaxed cutoffs (p<0.2, q<1.0)...\n")
      go_mf <- enrichGO(gene = as.character(module_entrez),
                        universe = as.character(universe_entrez),
                        OrgDb = org_db,
                        ont = "MF",
                        pvalueCutoff = 0.2,
                        qvalueCutoff = 1.0,
                        minGSSize = args$min_genesetsize,
                        maxGSSize = args$max_genesetsize,
                        readable = TRUE)
    }
    if (is.null(go_mf) || nrow(go_mf) == 0) {
      cat("    Still no MF results - trying with NO cutoffs (p<1.0, q<1.0)...\n")
      go_mf <- enrichGO(gene = as.character(module_entrez),
                        universe = as.character(universe_entrez),
                        OrgDb = org_db,
                        ont = "MF",
                        pvalueCutoff = 1.0,
                        qvalueCutoff = 1.0,
                        minGSSize = args$min_genesetsize,
                        maxGSSize = args$max_genesetsize,
                        readable = TRUE)
    }
    # Simplify MF results
    if (!is.null(go_mf) && nrow(go_mf) > 0) {
      go_mf_simplified <- clusterProfiler::simplify(go_mf, cutoff = 0.7, by = "p.adjust", select_fun = min)
      cat("    Simplified from", nrow(go_mf), "to", nrow(go_mf_simplified), "terms\n")
    }
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
                      universe = as.character(universe_entrez),
                      OrgDb = org_db,
                      ont = "CC",
                      pvalueCutoff = args$pvalue_cutoff,
                      qvalueCutoff = args$qvalue_cutoff,
                      minGSSize = args$min_genesetsize,
                      maxGSSize = args$max_genesetsize,
                      readable = TRUE)
    if (is.null(go_cc) || nrow(go_cc) == 0) {
      cat("    No CC results - trying with relaxed cutoffs (p<0.2, q<1.0)...\n")
      go_cc <- enrichGO(gene = as.character(module_entrez),
                        universe = as.character(universe_entrez),
                        OrgDb = org_db,
                        ont = "CC",
                        pvalueCutoff = 0.2,
                        qvalueCutoff = 1.0,
                        minGSSize = args$min_genesetsize,
                        maxGSSize = args$max_genesetsize,
                        readable = TRUE)
    }
    if (is.null(go_cc) || nrow(go_cc) == 0) {
      cat("    Still no CC results - trying with NO cutoffs (p<1.0, q<1.0)...\n")
      go_cc <- enrichGO(gene = as.character(module_entrez),
                        universe = as.character(universe_entrez),
                        OrgDb = org_db,
                        ont = "CC",
                        pvalueCutoff = 1.0,
                        qvalueCutoff = 1.0,
                        minGSSize = args$min_genesetsize,
                        maxGSSize = args$max_genesetsize,
                        readable = TRUE)
    }
    # Simplify CC results
    if (!is.null(go_cc) && nrow(go_cc) > 0) {
      go_cc_simplified <- clusterProfiler::simplify(go_cc, cutoff = 0.7, by = "p.adjust", select_fun = min)
      cat("    Simplified from", nrow(go_cc), "to", nrow(go_cc_simplified), "terms\n")
    }
  }, error = function(e) {
    cat("    WARNING: GO CC enrichment failed -", e$message, "\n")
    go_cc <<- NULL
  })
  
  if (is.null(go_cc)) {
    go_cc <- data.frame()
  }
  
  # ----- KEGG Enrichment -----
  cat("  Running KEGG enrichment...\n")
  cat("    DIAGNOSTIC: Sample KEGG IDs:", paste(head(module_kegg, 5), collapse = ", "), "\n")
  cat("    Module has", length(module_kegg), "genes with KEGG IDs\n")
  
  # Skip KEGG if no genes mapped
  if (length(module_kegg) == 0) {
    cat("    WARNING: No genes mapped to KEGG IDs - skipping KEGG enrichment\n")
    kegg <- data.frame()
  } else {
    # Check if genes are in KEGG database
    tryCatch({
      # Try to download KEGG pathway info
      kegg_info <- clusterProfiler::download_KEGG(species = "dme", keggType = "KEGG", keyType = "kegg")
      if (!is.null(kegg_info) && length(kegg_info) > 0) {
        all_kegg_genes <- unique(unlist(kegg_info$KEGGPATHID2EXTID))
        cat("    Total KEGG genes available for dme:", length(all_kegg_genes), "\n")
        # Check overlap
        overlap <- sum(as.character(module_kegg) %in% as.character(all_kegg_genes))
        cat("    Module genes found in KEGG:", overlap, "/", length(module_kegg), "\n")
      }
    }, error = function(e) {
      cat("    Could not check KEGG gene overlap:", e$message, "\n")
    })
    
    tryCatch({
      kegg <- enrichKEGG(gene = as.character(module_kegg),
                         universe = as.character(universe_kegg),
                         organism = "dme",  # Drosophila melanogaster
                         keyType = "ncbi-geneid",  # Specify we're using Entrez IDs
                         pvalueCutoff = args$pvalue_cutoff,
                         qvalueCutoff = args$qvalue_cutoff,
                         minGSSize = args$min_genesetsize,
                         maxGSSize = args$max_genesetsize)
      if (is.null(kegg) || nrow(kegg) == 0) {
        cat("    No KEGG results - trying with relaxed cutoffs (p<0.2, q<1.0)...\n")
        kegg <- enrichKEGG(gene = as.character(module_kegg),
                           universe = as.character(universe_kegg),
                           organism = "dme",
                           keyType = "ncbi-geneid",
                           pvalueCutoff = 0.2,
                           qvalueCutoff = 1.0,
                           minGSSize = args$min_genesetsize,
                           maxGSSize = args$max_genesetsize)
      }
      if (is.null(kegg) || nrow(kegg) == 0) {
        cat("    Still no KEGG results - trying with NO cutoffs (p<1.0, q<1.0)...\n")
        kegg <- enrichKEGG(gene = as.character(module_kegg),
                           universe = as.character(universe_kegg),
                           organism = "dme",
                           keyType = "ncbi-geneid",
                           pvalueCutoff = 1.0,
                           qvalueCutoff = 1.0,
                           minGSSize = args$min_genesetsize,
                           maxGSSize = args$max_genesetsize)
      }
      # Convert Entrez IDs to gene symbols for readability
      if (!is.null(kegg) && nrow(kegg) > 0) {
        kegg <- setReadable(kegg, OrgDb = org_db, keyType = "ENTREZID")
        cat("    KEGG enrichment found:", nrow(kegg), "pathways\n")
      } else {
        cat("    No KEGG pathways enriched even with relaxed cutoffs\n")
        cat("    This is common for Drosophila - KEGG has limited pathway coverage\n")
      }
    }, error = function(e) {
      cat("    WARNING: KEGG enrichment failed -", e$message, "\n")
      cat("    This may indicate KEGG database access issues or limited pathway annotations\n")
      kegg <<- data.frame()
    })
  }
  
  if (is.null(kegg)) {
    kegg <- data.frame()
  }
  
  # Store results
  enrichment_results[[current_module]] <- list(
    genes = module_genes,
    n_genes = length(module_genes),
    n_mapped = length(module_entrez),
    go_bp = go_bp,
    go_mf = go_mf,
    go_cc = go_cc,
    kegg = kegg,
    go_bp_simplified = if (exists("go_bp_simplified")) go_bp_simplified else data.frame(),
    go_mf_simplified = if (exists("go_mf_simplified")) go_mf_simplified else data.frame(),
    go_cc_simplified = if (exists("go_cc_simplified")) go_cc_simplified else data.frame()
  )
  
  # Combine all GO results for summary
  if (nrow(go_bp) > 0) {
    go_bp_df <- as.data.frame(go_bp) %>% mutate(module = current_module, ontology = "BP")
    all_enrichments[[paste(current_module, "GO_BP", sep = "_")]] <- go_bp_df
  }
  
  if (nrow(go_mf) > 0) {
    go_mf_df <- as.data.frame(go_mf) %>% mutate(module = current_module, ontology = "MF")
    all_enrichments[[paste(current_module, "GO_MF", sep = "_")]] <- go_mf_df
  }
  
  if (nrow(go_cc) > 0) {
    go_cc_df <- as.data.frame(go_cc) %>% mutate(module = current_module, ontology = "CC")
    all_enrichments[[paste(current_module, "GO_CC", sep = "_")]] <- go_cc_df
  }
  
  if (nrow(kegg) > 0) {
    kegg_df <- as.data.frame(kegg) %>% mutate(module = current_module, ontology = "KEGG")
    all_enrichments[[paste(current_module, "KEGG", sep = "_")]] <- kegg_df
  }
  
  cat("  Enrichments found:\n")
  cat("    GO BP:", nrow(go_bp), "\n")
  cat("    GO MF:", nrow(go_mf), "\n")
  cat("    GO CC:", nrow(go_cc), "\n")
  cat("    KEGG:", nrow(kegg), "\n")
  
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
    dplyr::arrange(module, pvalue)
  
  # Convert geneID to character to preserve gene list
  all_enrich_df$geneID <- as.character(all_enrich_df$geneID)
  
  # Select and reorder columns
  all_enrich_df <- all_enrich_df %>%
    dplyr::select(module, ontology, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
  
  write.csv(all_enrich_df, 
            file.path(output_dir, "all_module_enrichments.csv"), 
            row.names = FALSE)
  cat("  Saved: all_module_enrichments.csv\n")
  cat("    Total enrichments:", nrow(all_enrich_df), "\n")
} else {
  cat("  No enrichments found across all modules\n")
}

# ----- 9. Generate per-module enrichment reports -----
cat("Generating per-module enrichment reports...\n")

for (current_module in names(enrichment_results)) {
  results <- enrichment_results[[current_module]]
  
  # Create report file
  report_file <- file.path(output_dir, paste0("enrichment_", current_module, ".txt"))
  
  cat("Writing enrichment report for module:", current_module, "...\n", file = report_file, append = FALSE)
  
  # Write header
  cat("\n=== ENRICHMENT ANALYSIS FOR MODULE:", current_module, "===\n\n", 
      file = report_file, append = TRUE)
  cat("Genes in module:", results$n_genes, "\n",
      file = report_file, append = TRUE)
  cat("Genes mapped to Entrez IDs:", results$n_mapped, "\n\n",
      file = report_file, append = TRUE)
  
  # Write GO BP results
  cat("--- BIOLOGICAL PROCESS (GO BP) ---\n",
      file = report_file, append = TRUE)
  if (!is.null(results$go_bp) && nrow(results$go_bp) > 0) {
    go_bp_df <- as.data.frame(results$go_bp) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      dplyr::arrange(pvalue)
    # Convert geneID to character
    go_bp_df$geneID <- as.character(go_bp_df$geneID)
    write.table(go_bp_df, file = report_file, append = TRUE, quote = FALSE, 
                sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    cat("No significant enrichments\n", file = report_file, append = TRUE)
  }
  
  cat("\n--- MOLECULAR FUNCTION (GO MF) ---\n",
      file = report_file, append = TRUE)
  if (!is.null(results$go_mf) && nrow(results$go_mf) > 0) {
    go_mf_df <- as.data.frame(results$go_mf) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      dplyr::arrange(pvalue)
    # Convert geneID to character
    go_mf_df$geneID <- as.character(go_mf_df$geneID)
    write.table(go_mf_df, file = report_file, append = TRUE, quote = FALSE,
                sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    cat("No significant enrichments\n", file = report_file, append = TRUE)
  }
  
  cat("\n--- CELLULAR COMPONENT (GO CC) ---\n",
      file = report_file, append = TRUE)
  if (!is.null(results$go_cc) && nrow(results$go_cc) > 0) {
    go_cc_df <- as.data.frame(results$go_cc) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      dplyr::arrange(pvalue)
    # Convert geneID to character
    go_cc_df$geneID <- as.character(go_cc_df$geneID)
    write.table(go_cc_df, file = report_file, append = TRUE, quote = FALSE,
                sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    cat("No significant enrichments\n", file = report_file, append = TRUE)
  }
  
  cat("\n--- KEGG PATHWAYS ---\n",
      file = report_file, append = TRUE)
  if (!is.null(results$kegg) && nrow(results$kegg) > 0) {
    kegg_df <- as.data.frame(results$kegg) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      dplyr::arrange(pvalue)
    # Convert geneID to character
    kegg_df$geneID <- as.character(kegg_df$geneID)
    write.table(kegg_df, file = report_file, append = TRUE, quote = FALSE,
                sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    cat("No significant enrichments\n", file = report_file, append = TRUE)
  }
}

# ----- 10. Generate enrichment plots -----
cat("Generating enrichment plots...\n")

pdf(file.path(output_dir, "enrichment_plots.pdf"), width = 14, height = 15)

for (current_module in names(enrichment_results)) {
  results <- enrichment_results[[current_module]]
  
  plots_to_draw <- list()
  
  # GO BP dotplot
  if (nrow(results$go_bp) > 0) {
    tryCatch({
      p_bp <- dotplot(results$go_bp, showCategory = min(10, nrow(results$go_bp)), 
                      title = "Biological Process") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))
      plots_to_draw[["BP"]] <- p_bp
    }, error = function(e) {
      cat("    Warning: Could not create GO BP dotplot -", e$message, "\n")
    })
  }
  
  # GO MF dotplot
  if (nrow(results$go_mf) > 0) {
    tryCatch({
      p_mf <- dotplot(results$go_mf, showCategory = min(10, nrow(results$go_mf)),
                      title = "Molecular Function") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))
      plots_to_draw[["MF"]] <- p_mf
    }, error = function(e) {
      cat("    Warning: Could not create GO MF dotplot -", e$message, "\n")
    })
  }
  
  # GO CC dotplot
  if (nrow(results$go_cc) > 0) {
    tryCatch({
      p_cc <- dotplot(results$go_cc, showCategory = min(10, nrow(results$go_cc)),
                      title = "Cellular Component") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))
      plots_to_draw[["CC"]] <- p_cc
    }, error = function(e) {
      cat("    Warning: Could not create GO CC dotplot -", e$message, "\n")
    })
  }
  
  # KEGG dotplot
  if (nrow(results$kegg) > 0) {
    tryCatch({
      p_kegg <- dotplot(results$kegg, showCategory = min(10, nrow(results$kegg)),
                        title = "KEGG Pathways") +
                        theme(plot.title = element_text(hjust = 0.5, size = 10))
      plots_to_draw[["KEGG"]] <- p_kegg
    }, error = function(e) {
      cat("    Warning: Could not create KEGG dotplot -", e$message, "\n")
    })
  }
  
  # Combine plots into one figure
  if (length(plots_to_draw) > 0) {
    # Add overall title
    library(gridExtra)
    library(grid)
    grid_title <- textGrob(paste("Module", current_module, "- Functional Enrichment"), 
                           gp = gpar(fontsize = 14, fontface = "bold"))
    
    # Arrange plots based on how many we have
    if (length(plots_to_draw) == 3) {
      combined_plot <- grid.arrange(grobs = plots_to_draw, ncol = 1, top = grid_title)
    } else if (length(plots_to_draw) == 2) {
      combined_plot <- grid.arrange(grobs = plots_to_draw, ncol = 1, top = grid_title)
    } else {
      combined_plot <- grid.arrange(grobs = plots_to_draw, ncol = 1, top = grid_title)
    }
    
    print(combined_plot)
  } else if (nrow(results$go_bp) == 0 && nrow(results$go_mf) == 0 && nrow(results$go_cc) == 0 && nrow(results$kegg) == 0) {
    # Create a simple text plot if no enrichments found
    plot.new()
    text(0.5, 0.5, paste("No enrichments found for module", current_module), 
         cex = 1.5, col = "grey50")
  }
}

dev.off()
cat("  Saved: enrichment_plots.pdf\n")

# ----- 10b. Generate module summary bar chart -----
cat("Generating module summary bar chart...\n")

# Prepare data for plotting
plot_data <- data.frame()
for (current_module in names(enrichment_results)) {
  results <- enrichment_results[[current_module]]
  
  # Get top enriched term (including KEGG)
  all_go <- rbind(
    if (nrow(results$go_bp) > 0) as.data.frame(results$go_bp) else NULL,
    if (nrow(results$go_mf) > 0) as.data.frame(results$go_mf) else NULL,
    if (nrow(results$go_cc) > 0) as.data.frame(results$go_cc) else NULL,
    if (nrow(results$kegg) > 0) as.data.frame(results$kegg) else NULL
  )
  
  top_term <- "No enrichment"
  if (!is.null(all_go) && nrow(all_go) > 0) {
    top_term <- all_go %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::pull(Description)
    
    # Truncate long terms
    if (nchar(top_term) > 40) {
      top_term <- paste0(substr(top_term, 1, 37), "...")
    }
  }
  
  plot_data <- rbind(plot_data, data.frame(
    module = current_module,
    n_genes = results$n_genes,
    top_enrichment = top_term,
    stringsAsFactors = FALSE
  ))
}

# Sort by module size (ascending)
plot_data <- plot_data %>%
  dplyr::arrange(n_genes) %>%
  dplyr::mutate(module_label = ifelse(module == "grey", "non-modular genes", module))

# Move grey to the beginning (bottom after coord_flip)
if ("grey" %in% plot_data$module) {
  grey_row <- plot_data %>% dplyr::filter(module == "grey")
  other_rows <- plot_data %>% dplyr::filter(module != "grey")
  plot_data <- rbind(grey_row, other_rows)
}

# Create proper WGCNA module colors (module name IS the colour)
module_colors <- setNames(plot_data$module, plot_data$module)

# Create legend labels with module and top enrichment
legend_labels <- paste0(plot_data$top_enrichment)
names(legend_labels) <- plot_data$module

# Create simple bar chart showing module sizes
p_summary <- ggplot(plot_data, aes(x = factor(module_label, levels = module_label), y = n_genes, fill = module)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = module_colors, 
                    labels = legend_labels,
                    name = "Module: Top Enrichment") +
  coord_flip() +
  labs(title = "Module Sizes",
       x = "",
       y = "Number of Genes") +
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(hjust = 0.5))

# Save plot
ggsave(file.path(output_dir, "module_summary_barplot.pdf"), 
       plot = p_summary, width = 12, height = 8)
cat("  Saved: module_summary_barplot.pdf\n")

# ----- 11. Generate summary statistics -----
cat("\n=== ENRICHMENT SUMMARY ===\n")

# Create summary data frame
  module_summary <- data.frame()

for (current_module in names(enrichment_results)) {
  results <- enrichment_results[[current_module]]
  cat("\nModule:", current_module, "\n")
  cat("  Genes:", results$n_genes, "(mapped:", results$n_mapped, ")\n")
  cat("  GO BP enrichments:", nrow(results$go_bp), "\n")
  cat("  GO MF enrichments:", nrow(results$go_mf), "\n")
  cat("  GO CC enrichments:", nrow(results$go_cc), "\n")
  cat("  KEGG enrichments:", nrow(results$kegg), "\n")  # Show simplified enrichment summary
  if (!is.null(results$go_bp_simplified) && nrow(results$go_bp_simplified) > 0) {
    cat("  Simplified BP terms:", nrow(results$go_bp_simplified), "\n")
  }
  
  # Top enrichment (GO + KEGG)
  all_go <- rbind(
    if (nrow(results$go_bp) > 0) as.data.frame(results$go_bp) else NULL,
    if (nrow(results$go_mf) > 0) as.data.frame(results$go_mf) else NULL,
    if (nrow(results$go_cc) > 0) as.data.frame(results$go_cc) else NULL,
    if (nrow(results$kegg) > 0) as.data.frame(results$kegg) else NULL
  )
  
  top_term_1 <- NA
  top_term_2 <- NA
  top_term_3 <- NA
  top_pval_1 <- NA
  top_pval_2 <- NA
  top_pval_3 <- NA
  
  if (!is.null(all_go) && nrow(all_go) > 0) {
    top_enrich <- all_go %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice_head(n = 3) %>%
      dplyr::select(Description, pvalue, p.adjust)
    
    cat("  Top enrichments:\n")
    for (i in 1:nrow(top_enrich)) {
      cat("    -", top_enrich$Description[i], 
          "(p =", format(top_enrich$pvalue[i], digits = 3), ")\n")
    }
    
    # Store top terms for summary file
    if (nrow(top_enrich) >= 1) {
      top_term_1 <- top_enrich$Description[1]
      top_pval_1 <- top_enrich$pvalue[1]
    }
    if (nrow(top_enrich) >= 2) {
      top_term_2 <- top_enrich$Description[2]
      top_pval_2 <- top_enrich$pvalue[2]
    }
    if (nrow(top_enrich) >= 3) {
      top_term_3 <- top_enrich$Description[3]
      top_pval_3 <- top_enrich$pvalue[3]
    }
  }
  
  # Simplified functional summary
  simplified_term_1 <- NA
  simplified_term_2 <- NA
  simplified_term_3 <- NA
  
  if (!is.null(results$go_bp_simplified) && nrow(results$go_bp_simplified) > 0) {
    simplified_summary <- as.data.frame(results$go_bp_simplified) %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice_head(n = 3) %>%
      dplyr::select(Description)
    
    cat("  Simplified functional summary (BP):\n")
    for (i in 1:nrow(simplified_summary)) {
      cat("    *", simplified_summary$Description[i], "\n")
    }
    
    # Store simplified terms
    if (nrow(simplified_summary) >= 1) simplified_term_1 <- simplified_summary$Description[1]
    if (nrow(simplified_summary) >= 2) simplified_term_2 <- simplified_summary$Description[2]
    if (nrow(simplified_summary) >= 3) simplified_term_3 <- simplified_summary$Description[3]
  }
  
  # Add to summary dataframe
  module_summary <- rbind(module_summary, data.frame(
    module = current_module,
    n_genes = results$n_genes,
    n_mapped = results$n_mapped,
    n_BP_enrichments = nrow(results$go_bp),
    n_MF_enrichments = nrow(results$go_mf),
    n_CC_enrichments = nrow(results$go_cc),
    n_KEGG_enrichments = nrow(results$kegg),
    n_BP_simplified = if (!is.null(results$go_bp_simplified)) nrow(results$go_bp_simplified) else 0,
    top_term_1 = top_term_1,
    top_pvalue_1 = top_pval_1,
    top_term_2 = top_term_2,
    top_pvalue_2 = top_pval_2,
    top_term_3 = top_term_3,
    top_pvalue_3 = top_pval_3,
    simplified_BP_1 = simplified_term_1,
    simplified_BP_2 = simplified_term_2,
    simplified_BP_3 = simplified_term_3,
    stringsAsFactors = FALSE
  ))
}

# Write summary to file
write.csv(module_summary, 
          file.path(output_dir, "module_enrichment_summary.csv"), 
          row.names = FALSE)
cat("\n  Saved: module_enrichment_summary.csv\n")

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", output_dir, "\n")
cat("\nOutput files:\n")
cat("  - all_module_enrichments.csv (combined enrichments for all modules)\n")
cat("  - module_enrichment_summary.csv (one-line summary per module)\n")
cat("  - enrichment_[MODULE].txt (detailed reports per module)\n")
cat("  - enrichment_plots.pdf (visualization of top enrichments)\n")
cat("  - module_summary_barplot.pdf (overview of module sizes and enrichment counts)\n")
