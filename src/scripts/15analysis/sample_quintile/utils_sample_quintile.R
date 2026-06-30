# ==============================================================================
# utils_sample_quintile.R
#
# Shared visualisation helpers for the sample-quintile analysis pipeline.
# Sourced by:
#   plot_sample_quintile_overlap.R    -- draw_quintile_freq_heatmap
#   analyze_deviant_sample_gene_drivers.R -- all functions
#
# All functions take only explicit parameters, 
# and no hidden global variables needs to be considered (no closure captures).
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# ---- Clustered sample x Q1-Q5 frequency heatmap -----------------------------
# Draws the diverging blue->white->red heatmap and writes to out_pdf.
# freq_mat:   n_samples x 5 matrix; colnames "Q1"..."Q5"; values in [0,1].
# cond_label: title string.
# n_steps:    number of colour steps on each side of the 0.20 midpoint.
# out_pdf:    output file path.
draw_quintile_freq_heatmap <- function(freq_mat, cond_label, n_steps, out_pdf) {
  hm_mat           <- freq_mat
  colnames(hm_mat) <- paste0("Q", 1:5)

  lo_breaks <- seq(max(0,    min(hm_mat)), 0.20, length.out = n_steps + 1L)
  hi_breaks <- seq(0.20, min(1, max(hm_mat)), length.out = n_steps + 1L)
  hm_breaks <- unique(c(lo_breaks, hi_breaks[-1L]))

  n_lo <- sum(hm_breaks <  0.20)
  n_hi <- sum(hm_breaks >  0.20)
  # pheatmap requires length(color) == length(breaks) - 1.
  # Both palettes approach white at the 0.20 midpoint so the transition is smooth.
  hm_colors <- c(
    colorRampPalette(c("#2166AC", "#FFFFFF"))(n_lo),
    colorRampPalette(c("#FFFFFF", "#D6604D"))(n_hi)
  )

  hm_height <- max(6, ceiling(nrow(hm_mat) * 0.018))

  pdf(out_pdf, width = 6, height = hm_height)
  pheatmap(
    mat           = hm_mat,
    color         = hm_colors,
    breaks        = hm_breaks,
    cluster_rows  = TRUE,
    cluster_cols  = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    fontsize_col  = 12,
    border_color  = NA,
    main          = sprintf(
      "Sample Quintile Frequency  -  %s\n(rows clustered; expected freq ~0.20 per quintile)",
      cond_label
    )
  )
  invisible(dev.off())
  cat("  Heatmap ->", out_pdf, "\n")
}

# ---- Sequential colour builder for absolute correlation in [0, 1] -----------
# Returns list(breaks, colors) for use with pheatmap.
# n_steps: number of colour intervals on each side of the range midpoint.
make_hm_colors <- function(vals, n_steps) {
  bk <- seq(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE),
            length.out = 2L * n_steps + 1L)
  cl <- colorRampPalette(c("#FFFFFF", "#D6604D"))(2L * n_steps)
  list(breaks = bk, colors = cl)
}

# ---- Draw one correlation heatmap to the active PDF device ------------------
# Returns TRUE on success, FALSE if gene set is too small.
draw_cor_heatmap <- function(gene_set, sample_set, fq_named, or_named,
                             method_label, basis_label, q_label,
                             show_names, cond_label, cor_meth, expr_mat,
                             n_steps) {
  gs <- intersect(gene_set, rownames(expr_mat))
  if (length(gs) < 2L) {
    message("    [skip ", method_label, " (", basis_label,
            "): < 2 genes in expression matrix]")
    return(FALSE)
  }
  sub_mat     <- expr_mat[gs, sample_set, drop = FALSE]
  abs_cor_mat <- cor(t(sub_mat), method = cor_meth) |> abs()
  diag(abs_cor_mat) <- NA_real_

  hm <- make_hm_colors(abs_cor_mat[!is.na(abs_cor_mat)], n_steps)

  # Continuous gene-level annotations; pheatmap expects 2-colour vectors for
  # continuous variables (low-to-high gradient).
  ann_df <- data.frame(
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
      "%s  |  %s  -  %s\nbasis: %s  |  n_genes=%d  n_samples=%d  method=%s",
      cond_label, q_label, method_label, basis_label,
      length(gs), length(sample_set), cor_meth
    )
  )
  TRUE
}

# ---- Bootstrap co-expression comparison (driver vs random gene sets) --------
# Returns a ggplot object, or NULL if no valid driver set is found.
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

    pool   <- setdiff(rownames(expr_mat), gs)
    # Collect individual pairwise abs correlations from each random replicate so
    # the random distribution matches the same unit as obs_v (pairwise values,
    # not per-replicate means — means would collapse variance via CLT).
    rand_v <- unlist(lapply(seq_len(n_random), function(.) {
      rg <- sample(pool, min(length(gs), length(pool)), replace = FALSE)
      rc <- cor(t(expr_mat[rg, sample_set, drop = FALSE]), method = cor_meth)
      abs(rc[upper.tri(rc)])
    }))

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
        "Co-expression: driver genes vs. random gene sets  -  %s  |  %s",
        q_label, cond_label
      ),
      subtitle = sprintf(
        "top-%s%%  |  deviant-block samples n=%d  |  %d random sets  |  method=%s",
        pct_label, length(sample_set), n_random, cor_meth
      ),
      x = NULL,
      y = sprintf("absolute %s correlation", cor_meth)
    ) +
    theme_bw(base_size = 16) +
    theme(
      legend.position  = "none",
      plot.title       = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold")
    )
}

# ---- Gene score distribution: freq_q (M1) focus + M2/M3 expansion ----------
# Panel M1: density of freq_q for all genes; selected genes highlighted;
#   vertical dashed line marks the selection cut-off.
# Panels M2/M3: same layout for log(1+OR) and consensus rank respectively.
# n_top OR freq_threshold (exactly one non-NULL) controls M1 selection.
make_distribution_plot <- function(gene_dt, q_label, cond_label,
                                   pct_label,
                                   n_top = NULL, freq_threshold = NULL) {
  dt <- copy(gene_dt)
  dt[, log_or := log1p(odds_ratio)]

  # Determine M1 selection and the vertical-line cut-off value
  if (!is.null(freq_threshold)) {
    m1_sel     <- dt$freq_q >= freq_threshold
    cutoff_x   <- freq_threshold
    sel_label  <- sprintf("freq_q ≥ %.2f  (%d genes)", freq_threshold, sum(m1_sel))
  } else {
    m1_sel     <- dt$in_M1_selected
    cutoff_x   <- if (any(m1_sel)) min(dt$freq_q[m1_sel]) else NA_real_
    sel_label  <- sprintf("top-%d  (freq_q cut-off ≈ %.3f)", sum(m1_sel), cutoff_x)
  }

  # Pre-loop: one pair of rows (background + selected) per panel
  panels <- list(
    list(score = "freq_q",         label = "M1: quintile frequency",
         sel = m1_sel,              group = "M1 selected", color = "#D6604D"),
    list(score = "log_or",         label = "M2: log(1 + odds ratio)",
         sel = dt$in_M2_topN,       group = "M2 selected", color = "#4393C3"),
    list(score = "rank_consensus", label = "M3: consensus rank",
         sel = dt$in_M3_consensus,  group = "M3 selected", color = "#74C476")
  )
  row_blocks <- vector("list", length(panels) * 2L)
  blk <- 0L
  for (pi in seq_along(panels)) {
    p    <- panels[[pi]]
    vals <- dt[[p$score]]
    blk               <- blk + 1L
    row_blocks[[blk]] <- data.table(x = vals,          panel_lbl = p$label,
                                    group = "All genes", is_bg = TRUE)
    blk               <- blk + 1L
    row_blocks[[blk]] <- data.table(x = vals[p$sel],   panel_lbl = p$label,
                                    group = p$group,    is_bg = FALSE)
  }
  plot_dt <- rbindlist(row_blocks)
  panel_levels <- sapply(panels, `[[`, "label")
  group_levels <- c("All genes", sapply(panels, `[[`, "group"))
  plot_dt[, panel_f := factor(panel_lbl, levels = panel_levels)]
  plot_dt[, group_f := factor(group,     levels = group_levels)]

  colour_map <- c("All genes"   = "grey60",
                  "M1 selected" = "#D6604D",
                  "M2 selected" = "#4393C3",
                  "M3 selected" = "#74C476")
  fill_map   <- c("All genes"   = "grey85",
                  "M1 selected" = "#D6604D",
                  "M2 selected" = "#4393C3",
                  "M3 selected" = "#74C476")

  bg_dt  <- plot_dt[is_bg == TRUE]
  sel_dt <- plot_dt[is_bg == FALSE]

  # Vertical cut-off line shown only in the M1 panel
  vline_dt <- data.table(
    panel_f    = factor("M1: quintile frequency", levels = panel_levels),
    xintercept = cutoff_x
  )

  ggplot(mapping = aes(x = x, colour = group_f, fill = group_f)) +
    geom_density(data      = bg_dt,
                 alpha     = 0.20,
                 linewidth = 0.4,
                 key_glyph = "blank") +
    geom_density(data      = sel_dt,
                 alpha     = 0.35,
                 linewidth = 0.9) +
    geom_vline(data        = vline_dt,
               mapping     = aes(xintercept = xintercept),
               linetype    = "dashed",
               colour      = "#D6604D",
               linewidth   = 0.7,
               inherit.aes = FALSE) +
    facet_wrap(~panel_f, nrow = 1L, scales = "free_x") +
    scale_colour_manual(name = NULL, values = colour_map) +
    scale_fill_manual(  name = NULL, values = fill_map) +
    guides(
      fill   = guide_legend(override.aes = list(alpha = 0.5)),
      colour = guide_legend(override.aes = list(linewidth = 1.2))
    ) +
    labs(
      title    = sprintf(
        "Gene score distributions  |  %s  |  %s", q_label, cond_label),
      subtitle = sprintf(
        "top-%s%% deviant samples  |  M1 selection: %s",
        pct_label, sel_label),
      x = "score value",
      y = "density"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "top"
    )
}
