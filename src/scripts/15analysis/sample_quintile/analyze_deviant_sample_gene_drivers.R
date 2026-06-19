#!/usr/bin/env Rscript
# ==============================================================================
# Deviant-Sample Gene Driver Analysis
#
# For samples with strongly non-uniform quintile distributions (deviant samples,
# selected as the top X% by Cramér's V), this script identifies which genes
# drive each sample into its dominant quintile using three independent methods:
#
# METHOD 1  Quintile frequency
#   freq_q[g] = (# deviant_Q samples in quintile Q for gene g) / |deviant_Q|
#   Expected under null ≈ 0.20.  Rank all genes by freq_q; take top N.
#
# METHOD 2  Odds ratio + Fisher's exact test (BH-corrected)
#   2×2 contingency: deviant_Q vs. non-deviant × in-Q vs. not-in-Q for gene g.
#   OR[g] = (a*d)/(b*c) with 0.5 continuity correction.
#   p_fisher via fisher.test, BH-adjusted.  Top N restricted to padj < --padj-max.
#
# METHOD 3  Consensus rank aggregation (M1 + M2)
#   Combined rank = mean(rank_freq, rank_OR).
#   Additionally requires padj_fisher < --padj-max.  Top N by combined rank.
#
# For each method and quintile Q, two co-expression heatmaps are produced:
#   Heatmap A  full-sample basis (all samples in the expression matrix)
#   Heatmap B  quintile deviant block (only deviant_Q samples)
#
# Additionally, a bootstrap comparison boxplot shows whether driver genes are
# more co-expressed within the deviant block than random gene sets of equal
# size (--n-random replicates).  Observed pairwise correlations are compared to
# pooled correlations from random gene sets; visualised as a boxplot per method.
#
# OUTPUTS  (<output-dir>/)  [one set per threshold X in --top-pct, per quintile Q]
#   gene_drivers_Q{q}_{X}pct.tsv            per-gene metrics + method membership
#   coexpression_Q{q}_{X}pct.pdf            6-page: M1/M2/M3 × full/deviant heatmaps
#   coexpression_comparison_Q{q}_{X}pct.pdf 3-panel boxplot: driver vs random cors
#
# Usage:
#   Rscript analyze_deviant_sample_gene_drivers.R \
#     --sample-metrics  results/quintile_overlap_ct/sample_quintile_counts.tsv \
#     --quintile-matrix results/quintile_overlap_ct/quintile_matrix.rds \
#     --expr-file       data/processed/VST/vstdataCtrl.txt \
#     --output-dir      results/deviant_gene_drivers_ct \
#     --top-pct         5,10,30 \
#     --n-top-genes     1000 \
#     --n-random        100 \
#     --cor-method      spearman \
#     --condition-label "Ctrl"
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Dm.eg.db)
  library(openxlsx)
  library(grid)
  library(stringr)
})

# ---- CLI --------------------------------------------------------------------
parser <- ArgumentParser(description = "Deviant-sample gene driver analysis")
parser$add_argument("--sample-metrics",  required = TRUE,
                    help = "sample_quintile_counts.tsv (must contain cramer_v column)")
parser$add_argument("--quintile-matrix", required = TRUE,
                    help = "quintile_matrix.rds (n_genes x n_samples integer 1-5)")
parser$add_argument("--expr-file",       required = TRUE,
                    help = "VOOM expression matrix (tab-sep, genes x samples)")
parser$add_argument("--output-dir",      required = TRUE,
                    help = "Directory for all outputs (created if absent)")
parser$add_argument("--top-pct",         default = "10,30",
                    help = "Comma-separated percentile thresholds (default: '10,30')")
parser$add_argument("--n-top-genes",     type = "integer", default = 100L,
                    help = "Top-N genes per method for heatmaps (default: 100)")
parser$add_argument("--n-random",        type = "integer", default = 100L,
                    help = "Random gene sets for bootstrap comparison (default: 100)")
parser$add_argument("--cor-method",      default = "spearman",
                    choices = c("spearman", "pearson"),
                    help = "Correlation method (default: spearman)")
parser$add_argument("--padj-max",        type = "double", default = 0.05,
                    help = "BH-adjusted p-value cutoff for OR method (default: 0.05)")
parser$add_argument("--condition-label", default = "Condition",
                    help = "Label used in plot titles (default: 'Condition')")
parser$add_argument("--skip-enrichment", action = "store_true", default = FALSE,
                    help = "Disable pathway enrichment (default: enrichment runs)")
parser$add_argument("--utils-dir",       default = NULL,
                    help = "Directory containing utils_pathways.R (auto-resolved if absent)")
parser$add_argument("--kegg-cache",      default = NULL,
                    help = "Path to RDS KEGG cache; created on first run if absent")
args <- parser$parse_args()

cat("=== Deviant-Sample Gene Driver Analysis ===\n")
cat("Sample metrics:   ", args$sample_metrics,  "\n")
cat("Quintile matrix:  ", args$quintile_matrix,  "\n")
cat("Expression file:  ", args$expr_file,        "\n")
cat("Output directory: ", args$output_dir,       "\n\n")

for (f in c(args$sample_metrics, args$quintile_matrix, args$expr_file))
  if (!file.exists(f)) stop("File not found: ", f)

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

top_pcts    <- as.integer(trimws(strsplit(args$top_pct, ",")[[1L]]))
n_top       <- args$n_top_genes
n_rand      <- args$n_random
cor_method  <- args$cor_method
padj_thresh <- args$padj_max
n_steps     <- 100L

# ---- Load inputs ------------------------------------------------------------
cat("Loading sample metrics...\n")
sample_dt <- fread(args$sample_metrics, data.table = TRUE)
if (!"cramer_v" %in% colnames(sample_dt))
  stop("sample_metrics file must contain a 'cramer_v' column. ",
       "Re-run plot_sample_quintile_overlap.R to add it.")

cat("Loading quintile matrix...\n")
qmat        <- readRDS(args$quintile_matrix)
gene_ids    <- rownames(qmat)
sample_cols <- colnames(qmat)
n_genes     <- nrow(qmat)
cat("  Genes:", n_genes, "  Samples:", ncol(qmat), "\n")

cat("Loading expression matrix...\n")
expr_dt    <- fread(args$expr_file, data.table = TRUE)
gene_col_e <- colnames(expr_dt)[1L]
expr_samps <- setdiff(colnames(expr_dt), gene_col_e)
expr_mat   <- as.matrix(expr_dt[, ..expr_samps])
rownames(expr_mat) <- expr_dt[[gene_col_e]]
cat("  Expr genes:", nrow(expr_mat), "  Samples:", ncol(expr_mat), "\n\n")

# Align both matrices to their common sample set
common_samps <- intersect(sample_cols, expr_samps)
if (length(common_samps) < 2L)
  stop("Fewer than 2 samples shared between quintile matrix and expression matrix.")
qmat     <- qmat[, common_samps, drop = FALSE]
expr_mat <- expr_mat[, common_samps, drop = FALSE]
cat("  Common samples:", length(common_samps), "\n\n")

# ---- Sequential heatmap colour builder for absolute correlation ∈ [0, 1] ----
make_hm_colors <- function(vals) {
  bk <- seq(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE),
            length.out = 2L * n_steps + 1L)
  cl <- colorRampPalette(c("#FFFFFF", "#D6604D"))(2L * n_steps)
  list(breaks = bk, colors = cl)
}

# ---- Draw one correlation heatmap to the active PDF device ------------------
# Returns TRUE on success, FALSE if gene set is too small.
draw_cor_heatmap <- function(gene_set, sample_set, fq_named, or_named,
                             method_label, basis_label, q_label,
                             show_names, cond_label, cor_meth, expr_mat) {
  gs <- intersect(gene_set, rownames(expr_mat))
  if (length(gs) < 2L) {
    message("    [skip ", method_label, " (", basis_label,
            "): < 2 genes in expression matrix]")
    return(FALSE)
  }
  sub_mat <- expr_mat[gs, sample_set, drop = FALSE]
  abs_cor_mat <- cor(t(sub_mat), method = cor_meth) |> abs()
  diag(abs_cor_mat) <- NA_real_

  hm      <- make_hm_colors(abs_cor_mat[!is.na(abs_cor_mat)])

  # Continuous gene-level annotations; pheatmap expects 2-colour vectors for
  # continuous variables (low-to-high gradient).
  ann_df  <- data.frame(
    freq_q     = fq_named[gs],
    odds_ratio = pmin(or_named[gs], 50),   # cap extreme ORs for display
    row.names  = gs
  )
  ann_colors <- list(
    freq_q     = c("#FFFFFF", "#D6604D"),
    odds_ratio = c("#FFFFFF", "#4393C3")
  )

  pheatmap(
    mat               = abs_cor_mat,
    color             = hm$colors,
    breaks            = hm$breaks,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = show_names,
    show_colnames     = show_names,
    annotation_row    = ann_df,
    annotation_colors = ann_colors,
    border_color      = NA,
    na_col            = "grey90",
    main              = sprintf(
      "%s  |  %s  —  %s\nbasis: %s  |  n_genes=%d  n_samples=%d  method=%s",
      cond_label, q_label, method_label, basis_label,
      length(gs), length(sample_set), cor_meth
    )
  )
  TRUE
}

# ---- Bootstrap co-expression comparison (driver vs random gene sets) --------
# Returns a ggplot object.
make_bootstrap_boxplot <- function(driver_sets, sample_set, q_label,
                                   cond_label, pct_label, cor_meth,
                                   n_random, expr_mat) {
  rows <- vector("list", length(driver_sets))
  for (mi in seq_along(driver_sets)) {
    m_name <- names(driver_sets)[mi]
    gs     <- intersect(driver_sets[[mi]], rownames(expr_mat))
    if (length(gs) < 2L) next

    cor_obs <- cor(t(expr_mat[gs, sample_set, drop = FALSE]), method = cor_meth)
    obs_v   <- cor_obs[upper.tri(cor_obs)]

    pool    <- setdiff(rownames(expr_mat), gs)
    rand_v  <- vapply(seq_len(n_random), function(.) {
      rg <- sample(pool, min(length(gs), length(pool)), replace = FALSE)
      rc <- cor(t(expr_mat[rg, sample_set, drop = FALSE]), method = cor_meth)
      mean(rc[upper.tri(rc)], na.rm = TRUE)
    }, numeric(1L))

    rows[[mi]] <- rbind(
      data.table(method = m_name, group = "Driver genes", value = abs(obs_v)),
      data.table(method = m_name, group = "Random genes", value = abs(rand_v))
    )
  }
  plot_dt <- rbindlist(Filter(Negate(is.null), rows))
  if (nrow(plot_dt) == 0L) return(NULL)

  ggplot(plot_dt, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4,
                 linewidth = 0.5, alpha = 0.85) +
    facet_wrap(~method, nrow = 1) +
    scale_fill_manual(values = c("Driver genes" = "#D6604D",
                                  "Random genes" = "#92C5DE")) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    labs(
      title    = sprintf(
        "Co-expression: driver genes vs. random gene sets  —  %s  |  %s",
        q_label, cond_label
      ),
      subtitle = sprintf(
        "top-%s%%  |  deviant-block samples n=%d  |  %d random sets  |  method=%s",
        pct_label, length(sample_set), n_random, cor_meth
      ),
      x = NULL,
      y = sprintf("absolute %s correlation (pairwise)", cor_meth)
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position   = "none",
          plot.title        = element_text(face = "bold"),
          strip.background  = element_rect(fill = "grey92"),
          strip.text        = element_text(face = "bold"))
}

# ---- Dominant quintile per sample -------------------------------------------
freq_cols <- paste0("Q", 1:5, "_freq")
dom_q     <- apply(as.matrix(sample_dt[, ..freq_cols]), 1L, which.max)
sample_dt[, dominant_quintile := as.integer(dom_q)]

# ---- Enrichment pre-loop setup (computed once, shared across all iterations) -
# Sourcing, ID conversion, and KEGG cache are done here so each method-quintile
# call only passes pre-converted IDs — no repeated downloads or bitr calls.
enrich_active <- !args$skip_enrichment

if (enrich_active) {
  # Resolve utils_pathways.R
  resolve_utils <- function(utils_dir) {
    if (!is.null(utils_dir)) {
      path <- file.path(utils_dir, "utils_pathways.R")
      if (file.exists(path)) return(path)
      stop("utils_pathways.R not found in --utils-dir: ", utils_dir)
    }
    script_dir <- tryCatch(
      dirname(normalizePath(
        sub("--file=", "",
            grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]),
        mustWork = FALSE
      )),
      error = function(e) getwd()
    )
    candidates <- c(
      file.path(script_dir, "..", "..", "utils", "utils_pathways.R"),
      file.path(script_dir, "..", "utils", "utils_pathways.R")
    )
    found <- Filter(file.exists, normalizePath(candidates, mustWork = FALSE))
    if (length(found) == 0L)
      stop("Cannot auto-resolve utils_pathways.R. Pass --utils-dir.")
    found[1L]
  }
  utils_path <- resolve_utils(args$utils_dir)
  cat("Sourcing pathway utils:", utils_path, "\n")
  source(utils_path)

  # Convert all expression-matrix genes (full transcriptome) to ENTREZID once
  cat("Converting background genes (full expr matrix) -> ENTREZID...\n")
  bg_flybase   <- rownames(expr_mat)
  bg_conv      <- suppressMessages(
    bitr(bg_flybase, fromType = "FLYBASE", toType = "ENTREZID",
         OrgDb = org.Dm.eg.db, drop = TRUE)
  )
  bg_entrez_all <- bg_conv$ENTREZID
  cat(sprintf("  Background: %d genes -> %d ENTREZID\n",
              length(bg_flybase), length(bg_entrez_all)))

  # Build / load KEGG cache
  load_or_download_kegg <- function(cache_path) {
    if (!is.null(cache_path) && file.exists(cache_path)) {
      cat("  Loading KEGG cache from:", cache_path, "\n")
      return(readRDS(cache_path))
    }
    cat("  Downloading KEGG pathway data for dme...\n")
    sets <- clusterProfiler::download_KEGG("dme")
    if (!is.null(cache_path)) {
      dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(sets, cache_path)
      cat("  KEGG cache saved to:", cache_path, "\n")
    }
    sets
  }
  kegg_sets <- load_or_download_kegg(args$kegg_cache)

  # KEGG IDs for background
  bg_kegg_conv  <- suppressMessages(
    bitr_kegg(bg_entrez_all, fromType = "ncbi-geneid", toType = "kegg",
              organism = "dme")
  )
  bg_kegg_all   <- bg_kegg_conv$kegg
}

# ---- Enrichment helper (defined once, called per method gene set) -----------
# Runs GO + KEGG and writes xlsx + optional PDF for a single gene set.
# bg_entrez / bg_kegg / kegg_sets are pre-computed above.
run_enrichment_for_geneset <- function(gene_set_flybase, bg_entrez, bg_kegg,
                                       kegg_sets, out_stem, label,
                                       ontology = "BP", pval = 0.05, qval = 0.2,
                                       cond_label = "Condition") {
  # Query ID conversion
  query_conv <- suppressMessages(
    bitr(gene_set_flybase, fromType = "FLYBASE", toType = "ENTREZID",
         OrgDb = org.Dm.eg.db, drop = TRUE)
  )
  query_entrez <- if (nrow(query_conv) > 0L) query_conv$ENTREZID else character(0L)

  too_few <- length(query_entrez) < 2L
  if (too_few) {
    cat(sprintf("    [enrich %s] < 2 ENTREZID after conversion — writing empty marker\n",
                label))
  }

  # GO enrichment
  go_result <- if (!too_few) {
    tryCatch(
      enrichGO(
        gene          = query_entrez,
        universe      = bg_entrez,
        OrgDb         = org.Dm.eg.db,
        keyType       = "ENTREZID",
        ont           = ontology,
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
  } else { NULL }

  if (!is.null(go_result) && nrow(go_result) > 0L && ontology %in% c("BP","CC","MF")) {
    go_result <- tryCatch(
      clusterProfiler::simplify(go_result, cutoff = 0.7, by = "p.adjust",
                                 select_fun = min),
      error = function(e) go_result
    )
  }

  # KEGG enrichment
  query_kegg_conv <- if (!too_few) {
    tryCatch(
      bitr_kegg(query_entrez, fromType = "ncbi-geneid", toType = "kegg",
                organism = "dme"),
      error = function(e) data.frame()
    )
  } else { data.frame() }
  query_kegg <- if (nrow(query_kegg_conv) > 0L) query_kegg_conv$kegg else character(0L)

  kegg_result <- if (length(query_kegg) >= 2L) {
    tryCatch(
      enricher(
        gene          = query_kegg,
        universe      = bg_kegg,
        TERM2GENE     = kegg_sets$KEGGPATHID2EXTID,
        TERM2NAME     = kegg_sets$KEGGPATHID2NAME,
        pvalueCutoff  = pval,
        qvalueCutoff  = qval,
        pAdjustMethod = "fdr",
        minGSSize     = 10L,
        maxGSSize     = 500L
      ),
      error = function(e) { message("  enricher KEGG error: ", e$message); NULL }
    )
  } else { NULL }

  go_n   <- if (!is.null(go_result))   nrow(go_result)   else 0L
  kegg_n <- if (!is.null(kegg_result)) nrow(kegg_result) else 0L
  cat(sprintf("    [enrich %s] GO=%d  KEGG=%d\n", label, go_n, kegg_n))

  # xlsx — always written
  go_sheet <- if (go_n > 0L) {
    as.data.frame(go_result)
  } else {
    data.frame(message = sprintf(
      "No GO terms enriched (n_query=%d, n_ENTREZID=%d, ont=%s, p<%.2g)",
      length(gene_set_flybase), length(query_entrez), ontology, pval
    ))
  }
  kegg_sheet <- if (kegg_n > 0L) {
    df <- as.data.frame(kegg_result)
    df$Description <- gsub(" - Drosophila melanogaster \\(fruit fly\\)", "",
                            df$Description)
    df
  } else {
    data.frame(message = sprintf(
      "No KEGG pathways enriched (n_query=%d, n_KEGG=%d, p<%.2g)",
      length(gene_set_flybase), length(query_kegg), pval
    ))
  }
  openxlsx::write.xlsx(
    setNames(list(go_sheet, kegg_sheet),
             c(paste0("GO_", ontology), "KEGG")),
    file = paste0(out_stem, ".xlsx")
  )

  # PDF — only when ≥1 non-empty result
  if (go_n > 0L || kegg_n > 0L) {
    pdf(paste0(out_stem, ".pdf"), width = 10, height = 7)
    if (go_n > 0L) {
      print(
        enrichplot::dotplot(go_result, showCategory = 20) +
          scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
          labs(title = sprintf("GO %s  —  %s  |  %s", ontology, label, cond_label))
      )
    }
    if (kegg_n > 0L) {
      kegg_strip <- kegg_result
      kegg_strip@result$Description <- gsub(
        " - Drosophila melanogaster \\(fruit fly\\)", "",
        kegg_strip@result$Description
      )
      print(
        barplot(kegg_strip, showCategory = 20) +
          scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
          labs(title = sprintf("KEGG  —  %s  |  %s", label, cond_label))
      )
    }
    invisible(dev.off())
  }

  # Marker file when nothing found
  if (go_n == 0L && kegg_n == 0L) {
    con <- file(paste0(out_stem, "_empty.txt"), open = "w")
    writeLines(c(
      as.character(Sys.time()),
      sprintf("label: %s", label),
      sprintf("n_input=%d  n_ENTREZID=%d  n_KEGG=%d",
              length(gene_set_flybase), length(query_entrez), length(query_kegg)),
      "No pathways enriched in either GO or KEGG"
    ), con)
    close(con)
  }
  invisible(NULL)
}

# ---- Main loop: per top-pct threshold, per quintile -------------------------
for (pct in top_pcts) {
  cat(sprintf("\n===== Top %d%% deviant samples =====\n", pct))

  n_deviant  <- ceiling(nrow(sample_dt) * pct / 100)
  setorder(sample_dt, -cramer_v)
  deviant_dt <- sample_dt[seq_len(n_deviant)]
  cat(sprintf("  Selected %d deviant samples (cramer_v >= %.4f)\n",
              n_deviant, min(deviant_dt$cramer_v)))

  # Restrict all sample ID vectors to those present in the aligned qmat/expr_mat
  non_deviant_ids <- intersect(
    setdiff(sample_dt$sample_id, deviant_dt$sample_id), common_samps
  )

  for (q in 1:5) {
    deviant_q_ids <- intersect(
      deviant_dt[dominant_quintile == q, sample_id], common_samps
    )
    n_dq <- length(deviant_q_ids)

    if (n_dq < 5L) {
      cat(sprintf("  Q%d: only %d genes — skipping\n", q, n_dq))
      next
    }
    cat(sprintf("\n  --- Q%d  (%d deviant samples) ---\n", q, n_dq))

    qmat_dq <- qmat[, deviant_q_ids,    drop = FALSE]
    qmat_nd <- qmat[, non_deviant_ids,  drop = FALSE]
    q_int   <- as.integer(q)

    # ---- Method 1: quintile frequency ---------------------------------------
    cat("    Method 1: quintile frequency...\n")
    freq_q_m1 <- rowMeans(qmat_dq == q_int, na.rm = TRUE)
    rank_freq <- rank(-freq_q_m1, ties.method = "average")

    # ---- Method 3: odds ratio + Fisher's exact (BH) ------------------------
    cat("    Method 3: odds ratio + Fisher test (", n_genes, "genes)...\n")
    a_vec <- rowSums(qmat_dq == q_int, na.rm = TRUE)
    b_vec <- n_dq - a_vec
    c_vec <- rowSums(qmat_nd == q_int, na.rm = TRUE)
    d_vec <- ncol(qmat_nd) - c_vec

    or_vec <- ((a_vec + 0.5) * (d_vec + 0.5)) /
              ((b_vec + 0.5) * (c_vec + 0.5))

    p_fisher <- vapply(seq_len(n_genes), function(i) {
      tryCatch(
        fisher.test(matrix(c(a_vec[i], b_vec[i], c_vec[i], d_vec[i]), 2L),
                    alternative = "greater")$p.value,
        error = function(e) NA_real_
      )
    }, numeric(1L))
    padj_fisher <- p.adjust(p_fisher, method = "BH")
    rank_or     <- rank(-or_vec, ties.method = "average")

    # ---- Method 2: consensus rank aggregation --------------------------------
    rank_consensus <- (rank_freq + rank_or) / 2

    # Build full gene results table
    gene_dt <- data.table(
      gene_id        = gene_ids,
      freq_q         = freq_q_m1,
      rank_freq      = rank_freq,
      odds_ratio     = or_vec,
      p_fisher       = p_fisher,
      padj_fisher    = padj_fisher,
      rank_OR        = rank_or,
      rank_consensus = rank_consensus
    )

    # Top-N per method
    m1_genes <- gene_dt[order(rank_freq)][seq_len(min(n_top, .N)), gene_id]

    M2_pool  <- gene_dt[!is.na(padj_fisher) & padj_fisher < padj_thresh]
    M2_genes <- if (nrow(M2_pool) >= 2L) {
      M2_pool[order(rank_OR)][seq_len(min(n_top, nrow(M2_pool))), gene_id]
    } else {
      cat(sprintf("    M2: no genes pass padj < %.2f — using top %d by OR\n",
                  padj_thresh, n_top))
      gene_dt[order(rank_OR)][seq_len(min(n_top, .N)), gene_id]
    }

    M3_pool  <- gene_dt[is.na(padj_fisher) | padj_fisher < padj_thresh]
    M3_genes <- if (nrow(M3_pool) >= 2L) {
      M3_pool[order(rank_consensus)][seq_len(min(n_top, nrow(M3_pool))), gene_id]
    } else {
      gene_dt[order(rank_consensus)][seq_len(min(n_top, .N)), gene_id]
    }

    gene_dt[, in_M1_topN      := gene_id %in% m1_genes]
    gene_dt[, in_M2_topN      := gene_id %in% M2_genes]
    gene_dt[, in_M3_consensus := gene_id %in% M3_genes]
    setorder(gene_dt, rank_consensus)

    overlap_m1M2 <- sum(gene_dt$in_M1_topN & gene_dt$in_M2_topN)
    cat(sprintf("    M1 top-%d | M2 top-%d | overlap=%d | M3 consensus top-%d\n",
                length(m1_genes), length(M2_genes), overlap_m1M2, length(M3_genes)))

    # Write gene driver TSV
    out_gene <- file.path(args$output_dir,
                          sprintf("gene_drivers_Q%d_%dpct.tsv", q, pct))
    fwrite(gene_dt, out_gene, sep = "\t")
    cat("    Gene drivers ->", out_gene, "\n")

    # ---- Pathway enrichment for each method gene set -------------------------
    if (enrich_active) {
      enrich_dir <- file.path(args$output_dir, "enrichment")
      dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
      for (mset_enrich in list(
        list(genes = m1_genes, tag = sprintf("Q%d_%dpct_M1", q, pct)),
        list(genes = M2_genes, tag = sprintf("Q%d_%dpct_M2", q, pct)),
        list(genes = M3_genes, tag = sprintf("Q%d_%dpct_M3", q, pct))
      )) {
        out_enrich_stem <- file.path(enrich_dir,
                                     paste0("enrichment_", mset_enrich$tag))
        run_enrichment_for_geneset(
          gene_set_flybase = mset_enrich$genes,
          bg_entrez        = bg_entrez_all,
          bg_kegg          = bg_kegg_all,
          kegg_sets        = kegg_sets,
          out_stem         = out_enrich_stem,
          label            = mset_enrich$tag,
          cond_label       = args$condition_label
        )
      }
    }

    # Named vectors for heatmap annotations (gene_id → value)
    fq_named <- setNames(gene_dt$freq_q,     gene_dt$gene_id)
    or_named <- setNames(gene_dt$odds_ratio, gene_dt$gene_id)

    show_names <- (n_top <= 60L)
    q_label    <- sprintf("Q%d deviant", q)

    # # ---- Heatmaps ------------------------------------------------------------
    # out_hm <- file.path(args$output_dir,
    #                     sprintf("coexpression_Q%d_%dpct.pdf", q, pct))
    # cat("    Heatmaps ->", out_hm, "\n")

    # pdf(out_hm, width = 9, height = 8)
    # for (mset in list(
    #   list(genes = m1_genes, label = "M1: Quintile Frequency"),
    #   list(genes = M3_genes, label = "M3: Consensus Rank"),
    #   list(genes = M2_genes, label = "M2: Odds Ratio (Fisher)")
    # )) {
    #   draw_cor_heatmap(
    #     gene_set     = mset$genes, sample_set = common_samps,
    #     fq_named     = fq_named,   or_named   = or_named,
    #     method_label = mset$label,
    #     basis_label  = sprintf("all samples (n=%d)", length(common_samps)),
    #     q_label      = q_label,   show_names = show_names,
    #     cond_label   = args$condition_label, cor_meth = cor_method,
    #     expr_mat     = expr_mat
    #   )
    #   draw_cor_heatmap(
    #     gene_set     = mset$genes, sample_set = deviant_q_ids,
    #     fq_named     = fq_named,   or_named   = or_named,
    #     method_label = mset$label,
    #     basis_label  = sprintf("Q%d deviant block (n=%d)", q, n_dq),
    #     q_label      = q_label,   show_names = show_names,
    #     cond_label   = args$condition_label, cor_meth = cor_method,
    #     expr_mat     = expr_mat
    #   )
    # }
    # invisible(dev.off())

    # ---- Bootstrap comparison boxplot ----------------------------------------
    out_box <- file.path(args$output_dir,
                         sprintf("coexpression_comparison_Q%d_%dpct.pdf", q, pct))
    cat("    Bootstrap comparison ->", out_box, "\n")

    p_box <- make_bootstrap_boxplot(
      driver_sets = list(
        "M1: Frequency"   = m1_genes,
        "M3: Consensus"   = M3_genes,
        "M2: Odds Ratio"  = M2_genes
      ),
      sample_set  = deviant_q_ids,
      q_label     = q_label,
      cond_label  = args$condition_label,
      pct_label   = as.character(pct),
      cor_meth    = cor_method,
      n_random    = n_rand,
      expr_mat    = expr_mat
    )
    if (!is.null(p_box)) {
      pdf(out_box, width = 12, height = 5)
      print(p_box)
      invisible(dev.off())
    }
  }
}

cat("\n=== Done ===\n")
