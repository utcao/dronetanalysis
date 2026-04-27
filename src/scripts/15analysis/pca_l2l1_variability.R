#!/usr/bin/env Rscript
# ==============================================================================
# PCA: L2L1 Network Metrics + Expression Variability
#
# Three PCAs are produced:
#   CT  PCA  (15 features): CT L2L1 network + CT full-network metrics
#                           + CT delta_mean_mad + mean_ct + mad_ct
#   HS  PCA  (15 features): HS L2L1 network + HS full-network metrics
#                           + HS delta_mean_mad + mean_hs + mad_hs
#   Merged PCA (~30 features): union of both + mad_hs_minus_ct,
#     feature arrows coloured by group:
#     CT Network | HS Network | CT Full Net | HS Full Net |
#     CT Variability | HS Variability | CT Expression | HS Expression |
#     HS-CT Variability
#
# Preprocessing: log1p on L2L1/HL columns, then z-score all features.
# Missing values: complete cases only per PCA (sets may differ between CT and HS).
# Gene annotation: gProfiler ortholog CSV (ortholog_name col) or plain text (one SYMBOL/line).
#
# Usage:
#   Rscript pca_l2l1_variability.R \
#     --ct-file       run_voomct/results_ct_voom/rewiring_hubs_ct_anno_0408_2026.tsv \
#     --hs-file       run_voomhs/results_hs_voom/rewiring_hubs_hs_anno_0413_2026.tsv \
#     --ct-var-file   results/variability/voomct_all_genes_mad_summary.xlsx \
#     --hs-var-file   results/variability/voomhs_all_genes_mad_summary.xlsx \
#     --full-stats-file results/variability/full_mad_cv2_ranks.xlsx \
#     --ct-net-file   results/results_full_network_metrics_CT/network_metrics_summary.tsv \
#     --hs-net-file   results/results_full_network_metrics_HS/network_metrics_summary.tsv \
#     --gene-list     data/candidates/01_devcell_2025_gProfiler_hsapiens_dmelanogaster_0417.csv \
#     --output-dir    results/pca_gene_metrics
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(argparse)
  library(ggplot2)
  library(ggrepel)
  library(factoextra)
  library(pheatmap)
  library(plotly)
  library(htmlwidgets)
})

# ----- Feature-group colour palette (used in merged PCA) -----
GROUP_COLORS <- c(
  "CT Network"     = "#2166AC",   # L2L1 hub metrics from annotation TSV
  "HS Network"     = "#D6604D",
  "CT Full Net"    = "#41AB5D",   # full-network metrics from HDF5 (CT)
  "HS Full Net"    = "#FD8D3C",   # full-network metrics from HDF5 (HS)
  "CT Variability" = "#74ADD1",   # delta_mean_mad_ct (subgroup)
  "HS Variability" = "#F4A582",   # delta_mean_mad_hs (subgroup)
  "CT Expression"  = "#01665E",   # mean_ct, mad_ct, ...
  "HS Expression"  = "#8C510A",   # mean_hs, mad_hs, ...
  "HS−CT"     = "#762A83"    # mad_hs_minus_ct, cv2_hs_minus_ct
)

# ----- Broken-stick threshold -----
broken_stick <- function(n_vars) {
  sapply(seq_len(n_vars), function(k) sum(1 / seq(k, n_vars)) / n_vars)
}

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

# ----- Core PCA runner -----
# gene_map:       data.table with columns drosophila_symbol (match key) and label (display)
# feature_groups: named character vector mapping feature name → group label, or NULL
run_pca_condition <- function(feat_dt, cond_label, prefix, gene_map,
                              feature_groups = NULL) {
  cat("\n=== PCA:", cond_label, "===\n")

  meta_cols <- c("gene_id", "SYMBOL")
  feat_cols  <- setdiff(colnames(feat_dt), meta_cols)

  # Complete cases
  complete_idx <- complete.cases(feat_dt[, ..feat_cols])
  cat("  Complete cases:", sum(complete_idx),
      "(dropped", sum(!complete_idx), "with NA)\n")
  feat_dt <- feat_dt[complete_idx]

  # log1p on right-skewed count/sum columns
  log1p_cols <- grep("^(L2L1|HL_conn|degree|L2_n_edges|L1_conn_sum)", feat_cols, value = TRUE)
  for (col in log1p_cols) feat_dt[, (col) := log1p(get(col))]

  # z-score
  feat_mat        <- as.matrix(feat_dt[, ..feat_cols])
  rownames(feat_mat) <- feat_dt$gene_id
  feat_mat_scaled <- scale(feat_mat)

  low_sd <- names(which(attr(feat_mat_scaled, "scaled:scale") < 0.01))
  if (length(low_sd)) warning("Near-zero variance features: ", paste(low_sd, collapse = ", "))

  # PCA
  pca   <- prcomp(feat_mat_scaled, center = FALSE, scale. = FALSE)
  n_pcs <- ncol(pca$rotation)
  imp   <- summary(pca)$importance

  # Eigenvalue table
  eig_dt <- data.table(
    PC               = paste0("PC", seq_len(n_pcs)),
    std_dev          = pca$sdev,
    variance         = pca$sdev^2,
    pct_variance     = imp["Proportion of Variance", ] * 100,
    cum_pct_variance = imp["Cumulative Proportion", ]  * 100,
    broken_stick_pct = broken_stick(length(feat_cols))  * 100
  )
  cat("  Variance explained (PC1–", min(5, n_pcs), "): ",
      paste0(round(eig_dt$pct_variance[seq_len(min(5, n_pcs))], 1), "%",
             collapse = "  "), "\n", sep = "")

  # Loadings
  load_dt <- as.data.table(pca$rotation, keep.rownames = "feature")

  # Scores + gene annotation
  scores_dt <- as.data.table(pca$x, keep.rownames = "gene_id")
  scores_dt <- merge(feat_dt[, .(gene_id, SYMBOL)], scores_dt, by = "gene_id")
  # Match by Drosophila ortholog name (case-insensitive); display label from gene_map
  sym_upper <- toupper(scores_dt$SYMBOL)
  map_upper <- toupper(gene_map$drosophila_symbol)
  scores_dt[, is_annotated := sym_upper %in% map_upper]
  scores_dt[, label := {
    idx <- match(sym_upper, map_upper)
    ifelse(is_annotated, gene_map$label[idx], "")
  }]
  n_found <- sum(scores_dt$is_annotated)
  cat("  Annotated genes found:", n_found, "/", nrow(gene_map), "\n")
  missing_sym <- gene_map$drosophila_symbol[!map_upper %in% sym_upper]
  if (length(missing_sym))
    cat("  Not in dataset:", paste(missing_sym, collapse = ", "), "\n")

  # ------------------------------------------------------------------ Plots --

  ann_df <- as.data.frame(scores_dt[is_annotated == TRUE])
  non_df <- as.data.frame(scores_dt[is_annotated == FALSE])

  pct_var <- round(eig_dt$pct_variance, 1)

  # -- Scree plot --
  scree_df <- data.frame(
    PC  = factor(paste0("PC", seq_len(n_pcs)), levels = paste0("PC", seq_len(n_pcs))),
    pct = eig_dt$pct_variance,
    cum = eig_dt$cum_pct_variance,
    bs  = eig_dt$broken_stick_pct
  )
  p_scree <- ggplot(scree_df, aes(x = PC, y = pct)) +
    geom_col(fill = "#4A90D9", width = 0.7) +
    geom_line(aes(y = cum, group = 1), color = "#E05A44", linewidth = 1) +
    geom_point(aes(y = cum), color = "#E05A44", size = 2) +
    geom_line(aes(y = bs, group = 1), color = "#2CA02C",
              linewidth = 0.8, linetype = "dashed") +
    labs(title   = paste("Scree Plot —", cond_label),
         subtitle = "Red: cumulative variance; Green dashed: broken-stick threshold",
         x = NULL, y = "Variance explained (%)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(prefix, "_pca_scree.pdf"), p_scree, width = 8, height = 5)

  # -- Biplot helper --
  make_biplot <- function(pc_x, pc_y) {
    xi <- as.integer(sub("PC", "", pc_x))
    yi <- as.integer(sub("PC", "", pc_y))
    score_range <- max(abs(c(scores_dt[[pc_x]], scores_dt[[pc_y]])), na.rm = TRUE)
    load_scale  <- score_range * 0.8
    load_df <- data.frame(
      feature = load_dt$feature,
      x       = load_dt[[pc_x]] * load_scale,
      y       = load_dt[[pc_y]] * load_scale,
      stringsAsFactors = FALSE
    )
    # Attach group when feature_groups is provided
    use_groups <- !is.null(feature_groups)
    if (use_groups) {
      load_df$group <- factor(feature_groups[load_df$feature],
                              levels = names(GROUP_COLORS))
    }

    p <- ggplot() +
      geom_point(data = non_df,
                 aes(x = .data[[pc_x]], y = .data[[pc_y]]),
                 color = "grey70", alpha = 0.4, size = 0.6)

    if (use_groups) {
      p <- p +
        geom_segment(data = load_df,
                     aes(x = 0, y = 0, xend = x, yend = y, color = group),
                     arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                     linewidth = 0.85, show.legend = TRUE) +
        geom_text(data = load_df,
                  aes(x = x * 1.08, y = y * 1.08, label = feature, color = group),
                  size = 3.2, fontface = "bold", show.legend = FALSE) +
        scale_color_manual("Feature group", values = GROUP_COLORS,
                           drop = FALSE)
    } else {
      p <- p +
        geom_segment(data = load_df,
                     aes(x = 0, y = 0, xend = x, yend = y),
                     arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                     color = "#2C7BB6", linewidth = 0.75) +
        geom_text(data = load_df,
                  aes(x = x * 1.08, y = y * 1.08, label = feature),
                  color = "#2C7BB6", size = 3.2, fontface = "bold")
    }

    p +
      geom_point(data = ann_df,
                 aes(x = .data[[pc_x]], y = .data[[pc_y]]),
                 color = "#D7191C", size = 2.5, alpha = 0.9) +
      geom_text_repel(data = ann_df,
                      aes(x = .data[[pc_x]], y = .data[[pc_y]], label = label),
                      color = "#D7191C", size = 3.2, max.overlaps = 60,
                      box.padding = 0.3, segment.color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
      labs(title = paste("Biplot —", cond_label),
           x = paste0(pc_x, " (", pct_var[xi], "%)"),
           y = paste0(pc_y, " (", pct_var[yi], "%)")) +
      theme_bw(base_size = 12) +
      theme(legend.position = if (use_groups) "right" else "none")
  }
  ggsave(paste0(prefix, "_pca_biplot_PC1_PC2.pdf"),
         make_biplot("PC1", "PC2"), width = 8, height = 7)
  if (n_pcs >= 3)
    ggsave(paste0(prefix, "_pca_biplot_PC2_PC3.pdf"),
           make_biplot("PC2", "PC3"), width = 8, height = 7)

  # -- Interactive biplot (plotly → HTML) --
  make_interactive_biplot <- function(pc_x, pc_y) {
    xi <- as.integer(sub("PC", "", pc_x))
    yi <- as.integer(sub("PC", "", pc_y))
    score_range <- max(abs(c(scores_dt[[pc_x]], scores_dt[[pc_y]])), na.rm = TRUE)
    load_scale  <- score_range * 0.8

    load_df <- data.frame(
      feature = load_dt$feature,
      x       = load_dt[[pc_x]] * load_scale,
      y       = load_dt[[pc_y]] * load_scale,
      stringsAsFactors = FALSE
    )
    if (!is.null(feature_groups)) {
      load_df$group <- factor(feature_groups[load_df$feature],
                              levels = names(GROUP_COLORS))
      arrow_colors <- GROUP_COLORS[as.character(load_df$group)]
    } else {
      arrow_colors <- rep("#2C7BB6", nrow(load_df))
    }

    # Hover text for all genes: SYMBOL + PC scores
    scores_df <- as.data.frame(scores_dt)
    scores_df$hover <- paste0(
      "<b>", scores_df$SYMBOL, "</b><br>",
      pc_x, ": ", round(scores_df[[pc_x]], 3), "<br>",
      pc_y, ": ", round(scores_df[[pc_y]], 3)
    )

    non_df2 <- scores_df[!scores_df$is_annotated, ]
    ann_df2  <- scores_df[ scores_df$is_annotated, ]

    p <- plot_ly() |>
      # Background genes (grey, hover only)
      add_markers(data       = non_df2,
                  x          = ~.data[[pc_x]],
                  y          = ~.data[[pc_y]],
                  text       = ~hover,
                  hoverinfo  = "text",
                  marker     = list(color = "rgba(150,150,150,0.35)", size = 4),
                  name       = "All genes",
                  showlegend = TRUE) |>
      # Annotated genes (red, hover + permanent label)
      add_markers(data       = ann_df2,
                  x          = ~.data[[pc_x]],
                  y          = ~.data[[pc_y]],
                  text       = ~hover,
                  hoverinfo  = "text",
                  marker     = list(color = "#D7191C", size = 8,
                                    line = list(color = "white", width = 1)),
                  name       = "Annotated genes",
                  showlegend = TRUE) |>
      add_text(data      = ann_df2,
               x         = ~.data[[pc_x]],
               y         = ~.data[[pc_y]],
               text      = ~label,
               textfont  = list(color = "#D7191C", size = 11),
               textposition = "top right",
               hoverinfo = "skip",
               showlegend = FALSE)

    # Add one arrow per feature via add_annotations
    for (i in seq_len(nrow(load_df))) {
      p <- p |> add_annotations(
        x          = load_df$x[i],
        y          = load_df$y[i],
        ax         = 0, ay = 0,
        xref = "x", yref = "y", axref = "x", ayref = "y",
        text       = load_df$feature[i],
        font       = list(color = arrow_colors[i], size = 11),
        arrowcolor = arrow_colors[i],
        arrowwidth = 1.8,
        arrowhead  = 2,
        arrowsize  = 1,
        showarrow  = TRUE
      )
    }

    p <- p |> layout(
      title  = paste("Interactive Biplot —", cond_label,
                     paste0("(", pc_x, " vs ", pc_y, ")")),
      xaxis  = list(title      = paste0(pc_x, " (", pct_var[xi], "%)"),
                    zeroline   = TRUE, zerolinecolor = "#cccccc",
                    showgrid   = FALSE),
      yaxis  = list(title      = paste0(pc_y, " (", pct_var[yi], "%)"),
                    zeroline   = TRUE, zerolinecolor = "#cccccc",
                    showgrid   = FALSE),
      legend = list(orientation = "v", x = 1.02, y = 1),
      hovermode = "closest"
    )

    if (!is.null(feature_groups)) {
      # Add invisible dummy traces for the group colour legend
      present_groups <- unique(as.character(load_df$group))
      for (grp in present_groups) {
        p <- p |> add_markers(x = NA_real_, y = NA_real_,
                               marker = list(color = GROUP_COLORS[grp], size = 10),
                               name = grp, showlegend = TRUE)
      }
    }
    p
  }

  htmlwidgets::saveWidget(
    make_interactive_biplot("PC1", "PC2"),
    file         = paste0(prefix, "_pca_biplot_PC1_PC2_interactive.html"),
    selfcontained = TRUE
  )
  if (n_pcs >= 3)
    htmlwidgets::saveWidget(
      make_interactive_biplot("PC2", "PC3"),
      file         = paste0(prefix, "_pca_biplot_PC2_PC3_interactive.html"),
      selfcontained = TRUE
    )

  # -- Correlation circle --
  pdf(paste0(prefix, "_pca_correlation_circle.pdf"), width = 7, height = 7)
  if (!is.null(feature_groups)) {
    grp_vec <- factor(feature_groups[rownames(pca$rotation)],
                      levels = names(GROUP_COLORS))
    print(fviz_pca_var(pca,
                       col.var  = grp_vec,
                       palette  = GROUP_COLORS[levels(grp_vec)],
                       legend.title = "Feature group",
                       title    = paste("Correlation Circle —", cond_label),
                       repel    = TRUE))
  } else {
    print(fviz_pca_var(pca,
                       col.var       = "cos2",
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       title         = paste("Correlation Circle —", cond_label),
                       repel         = TRUE))
  }
  dev.off()

  # -- Loadings heatmap --
  n_show   <- min(5, n_pcs)
  load_mat <- as.matrix(load_dt[, paste0("PC", seq_len(n_show)), with = FALSE])
  rownames(load_mat) <- load_dt$feature
  annot_row <- if (!is.null(feature_groups)) {
    data.frame(Group = factor(feature_groups[rownames(load_mat)],
                              levels = names(GROUP_COLORS)),
               row.names = rownames(load_mat))
  } else NULL
  annot_colors <- if (!is.null(feature_groups))
    list(Group = GROUP_COLORS[names(GROUP_COLORS) %in% levels(annot_row$Group)])
  else NULL
  hm_width  <- if (!is.null(feature_groups)) 8.5 else 7
  hm_height <- if (!is.null(feature_groups)) 6   else 5
  pdf(paste0(prefix, "_pca_loadings_heatmap.pdf"), width = hm_width, height = hm_height)
  pheatmap(load_mat,
           main             = paste("Loadings —", cond_label),
           cluster_cols     = FALSE,
           color            = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
           display_numbers  = TRUE,
           number_format    = "%.2f",
           fontsize         = 10,
           border_color     = "grey85",
           annotation_row   = annot_row,
           annotation_colors = annot_colors)
  dev.off()

  # -- Annotated scores plot --
  p_ann <- ggplot() +
    geom_point(data = non_df, aes(x = PC1, y = PC2),
               color = "grey75", alpha = 0.45, size = 0.65) +
    geom_point(data = ann_df, aes(x = PC1, y = PC2),
               color = "#D7191C", size = 2.5) +
    geom_text_repel(data = ann_df, aes(x = PC1, y = PC2, label = label),
                    color = "#D7191C", size = 3.5, max.overlaps = 60,
                    box.padding = 0.4, segment.color = "grey55",
                    segment.linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    labs(title = paste("PCA Scores —", cond_label),
         x = paste0("PC1 (", pct_var[1], "%)"),
         y = paste0("PC2 (", pct_var[2], "%)")) +
    theme_bw(base_size = 12)
  ggsave(paste0(prefix, "_pca_annotated_scores.pdf"), p_ann, width = 8, height = 7)

  # ------------------------------------------------------------------ xlsx ---

  feat_out_dt <- cbind(feat_dt[, .(gene_id, SYMBOL)],
                       as.data.table(feat_mat_scaled))

  write_xlsx_safe(list(feature_matrix = feat_out_dt),
                  paste0(prefix, "_feature_matrix.xlsx"))
  write_xlsx_safe(list(eigenvalues = eig_dt),
                  paste0(prefix, "_pca_eigenvalues.xlsx"))
  write_xlsx_safe(list(loadings = load_dt),
                  paste0(prefix, "_pca_loadings.xlsx"))
  write_xlsx_safe(list(scores = scores_dt),
                  paste0(prefix, "_pca_scores.xlsx"))

  cat("  All outputs written with prefix:", basename(prefix), "\n")
  # feat_mat is post-log1p, pre-z-score; returned for downstream correlation plots
  feat_log1p_dt <- cbind(feat_dt[, .(gene_id, SYMBOL)], as.data.table(feat_mat))
  invisible(list(pca       = pca,
                 scores    = scores_dt,
                 eig       = eig_dt,
                 loadings  = load_dt,
                 feat_log1p = feat_log1p_dt))
}

# ==============================================================================
# Helper: load full-network metrics from HDF5
# ==============================================================================

load_network_tsv <- function(tsv_path, label) {
  target <- c("gene_id", "degree", "L2_n_edges", "L1_conn_mean", "L1_conn_sum", "mean_abs_corr")
  dt <- fread(tsv_path, data.table = TRUE)
  missing <- setdiff(target, names(dt))
  if (length(missing))
    stop(label, " TSV: required columns missing: ", paste(missing, collapse = ", "))
  dt <- dt[, .SD, .SDcols = target]
  dt[, gene_id := as.character(gene_id)]
  cat("  ", label, " network TSV: ", nrow(dt), " genes | metrics: degree, L2_n_edges, L1_conn_mean, L1_conn_sum, mean_abs_corr\n", sep = "")
  dt
}

# ==============================================================================
# Helpers: correlation plots
# ==============================================================================

sig_stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}

spearman_matrix <- function(mat) {
  n   <- ncol(mat)
  rho <- matrix(1,  n, n, dimnames = list(colnames(mat), colnames(mat)))
  pv  <- matrix(NA, n, n, dimnames = list(colnames(mat), colnames(mat)))
  diag(pv) <- 0
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      ok <- complete.cases(mat[, c(i, j)])
      if (sum(ok) < 5) next
      ct <- suppressWarnings(
        cor.test(mat[ok, i], mat[ok, j], method = "spearman", exact = FALSE))
      rho[i, j] <- rho[j, i] <- ct$estimate
      pv[i, j]  <- pv[j, i]  <- ct$p.value
    }
  }
  list(rho = rho, pval = pv)
}

make_corr_heatmap_gg <- function(feat_log1p_dt, cond_label, out_path,
                                  feature_groups = NULL) {
  meta  <- c("gene_id", "SYMBOL")
  fcols <- setdiff(colnames(feat_log1p_dt), meta)
  mat   <- as.matrix(feat_log1p_dt[, ..fcols])
  cr    <- spearman_matrix(mat)

  n     <- length(fcols)
  df    <- expand.grid(Var1 = fcols, Var2 = fcols, stringsAsFactors = FALSE)
  idx     <- cbind(df$Var1, df$Var2)
  df$rho  <- cr$rho[idx]
  df$pval <- cr$pval[idx]
  df$lbl  <- ifelse(df$Var1 == df$Var2, "1.00",
                    paste0(sprintf("%.2f", df$rho), "\n", sprintf("%.2e", df$pval)))

  # Order by feature group if provided, else keep original order
  if (!is.null(feature_groups)) {
    ord <- fcols[order(match(feature_groups[fcols], names(GROUP_COLORS)))]
  } else {
    ord <- fcols
  }
  df$Var1 <- factor(df$Var1, levels = ord)
  df$Var2 <- factor(df$Var2, levels = rev(ord))

  txt_sz  <- if (n > 18) 1.8 else if (n > 12) 2.2 else 2.8
  fig_sz  <- max(7, n * 0.48 + 2.5)

  p <- ggplot(df, aes(x = Var1, y = Var2, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = lbl), size = txt_sz, lineheight = 0.85) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                         midpoint = 0, limits = c(-1, 1), name = "rho") +
    labs(title    = paste("Feature Correlations —", cond_label),
         subtitle = "Spearman rho (top) | p-value in scientific notation (bottom)",
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y    = element_text(size = 7.5),
          legend.position = "right")
  ggsave(out_path, p, width = fig_sz, height = fig_sz * 0.92)
  invisible(p)
}

make_ct_vs_hs_scatter <- function(merged_log1p_dt, gene_map, out_dir) {
  all_cols <- colnames(merged_log1p_dt)
  ct_cols  <- grep("_ct$", all_cols, value = TRUE)
  hs_match <- sub("_ct$", "_hs", ct_cols)
  keep     <- hs_match %in% all_cols
  ct_cols  <- ct_cols[keep]
  hs_cols  <- hs_match[keep]
  base     <- sub("_ct$", "", ct_cols)

  if (!length(ct_cols)) {
    cat("  CT vs HS: no matched column pairs found.\n"); return(invisible(NULL))
  }

  sym_up  <- toupper(merged_log1p_dt$SYMBOL)
  map_up  <- toupper(gene_map$drosophila_symbol)
  is_ann  <- sym_up %in% map_up
  lbl_vec <- ifelse(is_ann, gene_map$label[match(sym_up, map_up)], "")

  out_path <- file.path(out_dir, "ct_vs_hs_per_feature.pdf")
  pdf(out_path, width = 6, height = 6)
  for (i in seq_along(base)) {
    xv <- merged_log1p_dt[[ct_cols[i]]]
    yv <- merged_log1p_dt[[hs_cols[i]]]
    ok <- !is.na(xv) & !is.na(yv)
    if (sum(ok) < 5) next
    ct_r    <- suppressWarnings(
      cor.test(xv[ok], yv[ok], method = "spearman", exact = FALSE))
    rho_lbl <- sprintf("rho = %.3f\np = %.2e", ct_r$estimate, ct_r$p.value)
    df_p    <- data.frame(x = xv, y = yv, ann = is_ann, lbl = lbl_vec)
    p <- ggplot(df_p, aes(x = x, y = y)) +
      geom_point(data = df_p[!df_p$ann, ],
                 color = "grey70", alpha = 0.35, size = 0.7) +
      geom_point(data = df_p[df_p$ann, ],
                 color = "#D7191C", size = 2.5) +
      geom_text_repel(data = df_p[df_p$ann, ],
                      aes(label = lbl), color = "#D7191C",
                      size = 3, max.overlaps = 30) +
      geom_smooth(method = "lm", formula = y ~ x,
                  color = "#2166AC", se = FALSE, linewidth = 0.8, linetype = "dashed") +
      annotate("text", x = -Inf, y = Inf, label = rho_lbl,
               hjust = -0.1, vjust = 1.3, size = 3.5, color = "#2166AC") +
      labs(title = paste("CT vs HS —", base[i]),
           x = ct_cols[i], y = hs_cols[i]) +
      theme_bw(base_size = 12)
    print(p)
  }
  dev.off()
  cat("  CT vs HS scatter:", basename(out_path), "(", length(ct_cols), "pages)\n")
}

make_within_diet_scatter <- function(feat_log1p_dt, cond_label, out_dir, gene_map) {
  meta  <- c("gene_id", "SYMBOL")
  fcols <- setdiff(colnames(feat_log1p_dt), meta)
  n     <- length(fcols)
  if (n < 2) return(invisible(NULL))

  sym_up  <- toupper(feat_log1p_dt$SYMBOL)
  map_up  <- toupper(gene_map$drosophila_symbol)
  is_ann  <- sym_up %in% map_up
  lbl_vec <- ifelse(is_ann, gene_map$label[match(sym_up, map_up)], "")

  slug     <- gsub("[^A-Za-z0-9]", "_", tolower(cond_label))
  out_path <- file.path(out_dir, paste0(slug, "_within_diet_feature_pairs.pdf"))
  pdf(out_path, width = 6, height = 6)
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      xv <- feat_log1p_dt[[fcols[i]]]
      yv <- feat_log1p_dt[[fcols[j]]]
      ok <- !is.na(xv) & !is.na(yv)
      if (sum(ok) < 5) next
      ct_r    <- suppressWarnings(
        cor.test(xv[ok], yv[ok], method = "spearman", exact = FALSE))
      rho_lbl <- sprintf("rho = %.3f\np = %.2e", ct_r$estimate, ct_r$p.value)
      df_p    <- data.frame(x = xv, y = yv, ann = is_ann, lbl = lbl_vec)
      p <- ggplot(df_p, aes(x = x, y = y)) +
        geom_point(data = df_p[!df_p$ann, ],
                   color = "grey70", alpha = 0.35, size = 0.7) +
        geom_point(data = df_p[df_p$ann, ],
                   color = "#D7191C", size = 2.5) +
        geom_text_repel(data = df_p[df_p$ann, ],
                        aes(label = lbl), color = "#D7191C",
                        size = 3, max.overlaps = 30) +
        geom_smooth(method = "lm", formula = y ~ x,
                    color = "grey40", se = FALSE,
                    linewidth = 0.7, linetype = "dashed") +
        annotate("text", x = -Inf, y = Inf, label = rho_lbl,
                 hjust = -0.1, vjust = 1.3, size = 3.5, color = "grey30") +
        labs(title = paste("Feature pair —", cond_label),
             x = fcols[i], y = fcols[j]) +
        theme_bw(base_size = 12)
      print(p)
    }
  }
  dev.off()
  n_pairs <- n * (n - 1) / 2
  cat("  Within-diet scatter:", basename(out_path), "(", n_pairs, "pages)\n")
}

# ==============================================================================
# CLI
# ==============================================================================

parser <- ArgumentParser(description = "PCA on L2L1 metrics + variability (CT and HS)")
parser$add_argument("--ct-file",       required = TRUE,
                    help = "CT annotation TSV (rewiring_hubs_ct_anno_*.tsv)")
parser$add_argument("--hs-file",       required = TRUE,
                    help = "HS annotation TSV (rewiring_hubs_hs_anno_*.tsv)")
parser$add_argument("--ct-var-file",   required = TRUE,
                    help = "CT variability summary xlsx (voomct_all_genes_mad_summary.xlsx)")
parser$add_argument("--hs-var-file",   required = TRUE,
                    help = "HS variability summary xlsx (voomhs_all_genes_mad_summary.xlsx)")
parser$add_argument("--full-stats-file", required = TRUE,
                    help = "Full-matrix gene stats xlsx (full_mad_cv2_ranks.xlsx): condition-specific mean, median, mad, cv2")
parser$add_argument("--ct-net-file",   required = TRUE,
                    help = "CT full-network metrics TSV (must contain: gene_id, degree, L2_n_edges, L1_conn_mean, L1_conn_sum, mean_abs_corr)")
parser$add_argument("--hs-net-file",   required = TRUE,
                    help = "HS full-network metrics TSV (must contain: gene_id, degree, L2_n_edges, L1_conn_mean, L1_conn_sum, mean_abs_corr)")
parser$add_argument("--gene-list",     required = TRUE,
                    help = "Plain text file with one gene SYMBOL per line to annotate")
parser$add_argument("--output-dir",    required = TRUE,
                    help = "Directory for all output files")
args <- parser$parse_args()

cat("=== PCA: L2L1 Metrics + Expression Variability ===\n")
for (nm in c("ct_file", "hs_file", "ct_var_file", "hs_var_file",
             "full_stats_file", "ct_net_file", "hs_net_file",
             "gene_list", "output_dir")) {
  cat(sprintf("  %-18s %s\n", paste0(nm, ":"), args[[nm]]))
}

# ----- Validate -----
for (f in c(args$ct_file, args$hs_file, args$ct_var_file, args$hs_var_file,
            args$full_stats_file, args$ct_net_file, args$hs_net_file,
            args$gene_list)) {
  if (!file.exists(f)) stop("File not found: ", f)
}
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = TRUE)

# ----- Load sources -----
cat("\nLoading data...\n")
net_cols <- c("L2L1_deg", "L2L1_rewire", "L2L1_conn")
#"L1_mean_delta", "L2_mean_delta", "HL_deg_L1", "HL_deg_L2","HL_conn_L1", "HL_conn_L2")

l2l1_cols <- c("gene_id", "SYMBOL", net_cols)

ct_ann   <- fread(args$ct_file,       data.table = TRUE)[, ..l2l1_cols]
hs_ann   <- fread(args$hs_file,       data.table = TRUE)[, ..l2l1_cols]
ct_var   <- as.data.table(read.xlsx(args$ct_var_file))[,  .(gene_id, delta_mean_mad)]
hs_var   <- as.data.table(read.xlsx(args$hs_var_file))[,  .(gene_id, delta_mean_mad)]
full_var <- as.data.table(read.xlsx(args$full_stats_file))[, .(
  gene_id,
  mean_ct, median_ct, mad_ct, cv2_ct,
  mean_hs, median_hs, mad_hs, cv2_hs,
  mad_hs_minus_ct, cv2_hs_minus_ct
)]
cat("  CT annotation:   ", nrow(ct_ann),   "genes\n")
cat("  HS annotation:   ", nrow(hs_ann),   "genes\n")
cat("  CT var summary:  ", nrow(ct_var),   "genes\n")
cat("  HS var summary:  ", nrow(hs_var),   "genes\n")
cat("  Full stats:      ", nrow(full_var), "genes\n")
cat("Loading full-network TSV metrics...\n")
ct_net <- load_network_tsv(args$ct_net_file, "CT")
hs_net <- load_network_tsv(args$hs_net_file, "HS")

# ----- Load gene list -----
# Accepts two formats:
#   (a) gProfiler ortholog CSV with columns: initial_alias, ortholog_name, gene_id, ...
#       The Drosophila ortholog_name is used for matching; it is also the plot label.
#   (b) Plain text file (one SYMBOL per line). Lines that look like citation keys
#       (contain a 4-digit year or are longer than 30 characters) are skipped.
if (grepl("\\.csv$", args$gene_list, ignore.case = TRUE)) {
  gp_dt     <- fread(args$gene_list, data.table = TRUE)
  req_cols  <- c("initial_alias", "ortholog_name")
  if (!all(req_cols %in% colnames(gp_dt)))
    stop("gProfiler CSV must contain columns: ", paste(req_cols, collapse = ", "))
  gene_map  <- unique(gp_dt[!is.na(ortholog_name) & nzchar(ortholog_name) &
                              ortholog_name != "N/A",
                             .(drosophila_symbol = ortholog_name,
                               label             = ortholog_name)])
  cat("  Gene list (gProfiler CSV):", nrow(gene_map), "ortholog symbols\n")
} else {
  raw <- trimws(readLines(args$gene_list))
  raw <- raw[nzchar(raw)]
  keep <- grepl("^[A-Za-z][A-Za-z0-9_.\\-]*$", raw) &
          nchar(raw) <= 30 &
          !grepl("[0-9]{4}", raw)
  skipped <- raw[!keep]
  if (length(skipped))
    cat("  Skipped non-gene lines:", paste(skipped, collapse = ", "), "\n")
  gene_map <- data.table(drosophila_symbol = raw[keep], label = raw[keep])
  cat("  Gene list (plain text):  ", nrow(gene_map), "symbols\n")
}

# ----- Build feature matrices -----
ct_expr_var <- full_var[, .(gene_id, mean_ct, mad_ct)]
hs_expr_var <- full_var[, .(gene_id, mean_hs, mad_hs)]
cross_var   <- full_var[, .(gene_id, mad_hs_minus_ct)]

merge_fn <- function(a, b) merge(a, b, by = "gene_id", all = FALSE)

ct_feat <- Reduce(merge_fn, list(ct_ann, ct_var, ct_expr_var, ct_net))
hs_feat <- Reduce(merge_fn, list(hs_ann, hs_var, hs_expr_var, hs_net))

cat("\nCT feature matrix before NA filter:", nrow(ct_feat), "genes\n")
cat("HS feature matrix before NA filter:", nrow(hs_feat), "genes\n")

# ----- Build merged feature matrix -----
# All per-condition columns get _ct / _hs suffix; shared cross-condition added last.
full_net_cols <- c("degree", "L2_n_edges", "L1_conn_mean", "L1_conn_sum", "mean_abs_corr")
all_cond_cols <- c(net_cols, full_net_cols)   # columns needing suffix renaming

ct_side <- copy(ct_feat)
setnames(ct_side, all_cond_cols,    paste0(all_cond_cols, "_ct"))
setnames(ct_side, "delta_mean_mad", "delta_mean_mad_ct")

hs_side <- copy(hs_feat)
setnames(hs_side, all_cond_cols,    paste0(all_cond_cols, "_hs"))
setnames(hs_side, "delta_mean_mad", "delta_mean_mad_hs")
hs_side[, SYMBOL := NULL]

merged_feat <- Reduce(merge_fn, list(ct_side, hs_side, cross_var))
cat("\nMerged feature matrix (CT+HS) before NA filter:", nrow(merged_feat), "genes\n")

# Feature → group mapping for the merged PCA
merged_feat_groups <- c(
  setNames(rep("CT Network",  length(net_cols)),      paste0(net_cols,      "_ct")),
  setNames(rep("HS Network",  length(net_cols)),      paste0(net_cols,      "_hs")),
  setNames(rep("CT Full Net", length(full_net_cols)), paste0(full_net_cols, "_ct")),
  setNames(rep("HS Full Net", length(full_net_cols)), paste0(full_net_cols, "_hs")),
  delta_mean_mad_ct = "CT Variability",
  delta_mean_mad_hs = "HS Variability",
  mean_ct           = "CT Expression",
  mad_ct            = "CT Expression",
  mean_hs           = "HS Expression",
  mad_hs            = "HS Expression",
  mad_hs_minus_ct   = "HS−CT"
)

# ----- Run PCAs -----
ct_prefix     <- file.path(args$output_dir, "ct")
hs_prefix     <- file.path(args$output_dir, "hs")
merged_prefix <- file.path(args$output_dir, "merged")

ct_res     <- run_pca_condition(ct_feat,     "CT (Control)",    ct_prefix,     gene_map)
hs_res     <- run_pca_condition(hs_feat,     "HS (High Sugar)", hs_prefix,     gene_map)
merged_res <- run_pca_condition(merged_feat, "CT + HS (Merged)",merged_prefix, gene_map,
                                feature_groups = merged_feat_groups)

# ----- Correlation plots (separate subfolders) -----
cat("\n=== Generating correlation plots ===\n")

dir_heatmap <- file.path(args$output_dir, "corr_heatmaps")
dir_ct_vs_hs <- file.path(args$output_dir, "corr_ct_vs_hs")
dir_within   <- file.path(args$output_dir, "corr_within_diet")
for (d in c(dir_heatmap, dir_ct_vs_hs, dir_within))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Correlation heatmaps (one per PCA, on post-log1p pre-z-score values)
cat("  Heatmaps...\n")
make_corr_heatmap_gg(ct_res$feat_log1p,  "CT (Control)",
                     file.path(dir_heatmap, "ct_corr_heatmap.pdf"))
make_corr_heatmap_gg(hs_res$feat_log1p,  "HS (High Sugar)",
                     file.path(dir_heatmap, "hs_corr_heatmap.pdf"))
make_corr_heatmap_gg(merged_res$feat_log1p, "CT + HS (Merged)",
                     file.path(dir_heatmap, "merged_corr_heatmap.pdf"),
                     feature_groups = merged_feat_groups)

# CT vs HS scatter (one page per matched feature pair)
cat("  CT vs HS scatter...\n")
make_ct_vs_hs_scatter(merged_res$feat_log1p, gene_map, dir_ct_vs_hs)

# Within-diet pairwise feature scatter
cat("  Within-diet scatter...\n")
make_within_diet_scatter(ct_res$feat_log1p, "CT (Control)",    dir_within, gene_map)
make_within_diet_scatter(hs_res$feat_log1p, "HS (High Sugar)", dir_within, gene_map)

cat("\n=== Done. All outputs in:", args$output_dir, "===\n")
