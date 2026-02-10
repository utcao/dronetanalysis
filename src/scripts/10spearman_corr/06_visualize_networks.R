#!/usr/bin/env Rscript
#' Stage 6: Differential Network Visualization with igraph
#'
#' Visualizes low, high, and differential networks side-by-side,
#' highlighting rewiring genes and edge type changes.
#'
#' Usage:
#'   Rscript 06_visualize_networks.R --data-dir visualization_data --top-n 10 --output-prefix networks
#'
#' Input (from Stage 5):
#'   visualization_data/
#'   ├── edges_low.tsv
#'   ├── edges_high.tsv
#'   ├── edges_diff.tsv
#'   └── nodes.tsv
#'
#' Output:
#'   ├── networks_ego_comparison.pdf   - Side-by-side ego networks (low, high, diff)
#'   ├── networks_degree_analysis.pdf  - Degree distribution and changes
#'   ├── networks_qualitative.pdf      - Qualitative change summary
#'   └── networks_summary.txt          - Text summary

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(argparse)
})

# =============================================================================
# Color Schemes
# =============================================================================

QUAL_COLORS <- c(
  "unchanged" = "#999999",
  "disappear" = "#d73027",   # red
  "new" = "#1a9850",         # green
  "sign_change" = "#7570b3", # purple
  "strengthen" = "#fdae61",  # orange
  "weaken" = "#abd9e9"       # light blue
)

SIGN_COLORS <- c(
  "positive" = "#e31a1c",    # red
  "negative" = "#1f78b4"     # blue
)

# =============================================================================
# Data Loading
# =============================================================================

load_network_data <- function(data_dir) {
  cat("Loading network data from:", data_dir, "\n")

  edges_low <- read.delim(file.path(data_dir, "edges_low.tsv"), stringsAsFactors = FALSE)
  edges_high <- read.delim(file.path(data_dir, "edges_high.tsv"), stringsAsFactors = FALSE)
  edges_diff <- read.delim(file.path(data_dir, "edges_diff.tsv"), stringsAsFactors = FALSE)
  nodes <- read.delim(file.path(data_dir, "nodes.tsv"), stringsAsFactors = FALSE)

  cat("  Low network:", nrow(edges_low), "edges\n")
  cat("  High network:", nrow(edges_high), "edges\n")
  cat("  Differential:", nrow(edges_diff), "edges\n")
  cat("  Nodes:", nrow(nodes), "genes with edges\n")

  list(
    edges_low = edges_low,
    edges_high = edges_high,
    edges_diff = edges_diff,
    nodes = nodes
  )
}

# =============================================================================
# Network Construction
# =============================================================================

create_network <- function(edges, nodes, value_col = "r") {
  # Create graph from edge list
  g <- graph_from_data_frame(
    edges[, c("from", "to")],
    directed = FALSE,
    vertices = nodes
  )

  # Add edge attributes
  if (value_col %in% names(edges)) {
    E(g)$weight <- edges[[value_col]]
  }
  if ("abs_r" %in% names(edges)) {
    E(g)$abs_weight <- edges$abs_r
  }

  g
}

# =============================================================================
# Ego Network Extraction
# =============================================================================

extract_ego_subgraph <- function(g, ego_gene, edges_df = NULL) {
  if (!(ego_gene %in% V(g)$name)) {
    warning(paste("Gene", ego_gene, "not in network"))
    return(NULL)
  }

  # Get neighbors
  ego_idx <- which(V(g)$name == ego_gene)
  neighbor_idx <- neighbors(g, ego_idx)
  all_idx <- c(ego_idx, neighbor_idx)

  # Induced subgraph
  subg <- induced_subgraph(g, all_idx)

  # Mark ego node
  V(subg)$is_ego <- V(subg)$name == ego_gene

  subg
}

# =============================================================================
# Visualization: Side-by-Side Ego Networks
# =============================================================================

plot_ego_comparison <- function(data, gene_id, layout_seed = 42) {
  # Create networks
  g_low <- create_network(data$edges_low, data$nodes, "r")
  g_high <- create_network(data$edges_high, data$nodes, "r")
  g_diff <- create_network(data$edges_diff, data$nodes, "delta")

  # Extract ego subgraphs
  ego_low <- extract_ego_subgraph(g_low, gene_id)
  ego_high <- extract_ego_subgraph(g_high, gene_id)
  ego_diff <- extract_ego_subgraph(g_diff, gene_id)

  if (is.null(ego_diff)) {
    return(NULL)
  }

  # Common nodes across all networks for consistent layout
  all_nodes <- unique(c(V(ego_low)$name, V(ego_high)$name, V(ego_diff)$name))

  # Use diff network layout as reference
  set.seed(layout_seed)
  layout_ref <- layout_with_fr(ego_diff)
  rownames(layout_ref) <- V(ego_diff)$name

  # Plot function
  plot_single_ego <- function(g, title, edge_color_fn) {
    if (is.null(g) || vcount(g) == 0) {
      plot.new()
      title(title)
      return()
    }

    # Node aesthetics
    node_colors <- ifelse(V(g)$is_ego, "red", "lightblue")
    node_sizes <- pmin(V(g)$degree_diff / max(V(g)$degree_diff, 1) * 12 + 3, 20)

    # Edge aesthetics
    edge_colors <- edge_color_fn(g)
    edge_widths <- pmin(abs(E(g)$weight) * 3 + 0.5, 5)

    # Use consistent layout
    node_names <- V(g)$name
    layout_g <- matrix(0, nrow = length(node_names), ncol = 2)
    for (i in seq_along(node_names)) {
      if (node_names[i] %in% rownames(layout_ref)) {
        layout_g[i, ] <- layout_ref[node_names[i], ]
      } else {
        layout_g[i, ] <- runif(2, -1, 1)
      }
    }

    plot(g,
         layout = layout_g,
         vertex.size = node_sizes,
         vertex.color = node_colors,
         vertex.label.cex = 0.6,
         vertex.label.color = "black",
         edge.color = edge_colors,
         edge.width = edge_widths,
         main = title)
  }

  # Edge color functions
  color_by_sign <- function(g) {
    w <- E(g)$weight
    ifelse(w > 0, SIGN_COLORS["positive"], SIGN_COLORS["negative"])
  }

  color_by_qual <- function(g) {
    # Need to get qual_label from edges_diff
    edge_list <- as_edgelist(g)
    colors <- rep("gray50", nrow(edge_list))

    for (i in seq_len(nrow(edge_list))) {
      from <- edge_list[i, 1]
      to <- edge_list[i, 2]
      idx <- which((data$edges_diff$from == from & data$edges_diff$to == to) |
                   (data$edges_diff$from == to & data$edges_diff$to == from))
      if (length(idx) > 0) {
        qual <- data$edges_diff$qual_label[idx[1]]
        if (qual %in% names(QUAL_COLORS)) {
          colors[i] <- QUAL_COLORS[qual]
        }
      }
    }
    colors
  }

  # Create 3-panel plot
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))

  plot_single_ego(ego_low, paste0("LOW: ", gene_id), color_by_sign)
  plot_single_ego(ego_high, paste0("HIGH: ", gene_id), color_by_sign)
  plot_single_ego(ego_diff, paste0("DIFF: ", gene_id), color_by_qual)
}

# =============================================================================
# Visualization: Degree Analysis
# =============================================================================

plot_degree_analysis <- function(data) {
  nodes <- data$nodes

  # Panel 1: Degree distribution comparison
  p1 <- nodes %>%
    select(gene_id, degree_low, degree_high) %>%
    pivot_longer(cols = c(degree_low, degree_high),
                 names_to = "network",
                 values_to = "degree") %>%
    mutate(network = ifelse(network == "degree_low", "Low", "High")) %>%
    ggplot(aes(x = degree, fill = network)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = c("Low" = "#3182bd", "High" = "#e6550d")) +
    labs(title = "Degree Distribution: Low vs High",
         x = "Degree", y = "Count") +
    theme_minimal() +
    theme(legend.position = "top")

  # Panel 2: Degree change scatter
  p2 <- ggplot(nodes, aes(x = degree_low, y = degree_high, color = rewiring_score)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_viridis_c(name = "Rewiring\nScore") +
    labs(title = "Degree Change: Low vs High",
         x = "Degree (Low)", y = "Degree (High)") +
    theme_minimal()

  # Panel 3: Rewiring score distribution
  p3 <- ggplot(nodes, aes(x = rewiring_score)) +
    geom_histogram(bins = 50, fill = "#756bb1", alpha = 0.7) +
    labs(title = "Rewiring Score Distribution",
         x = "Rewiring Score", y = "Count") +
    theme_minimal()

  # Panel 4: Top rewiring genes
  top_genes <- nodes %>%
    arrange(desc(rewiring_score)) %>%
    head(15) %>%
    mutate(gene_id = factor(gene_id, levels = rev(gene_id)))

  p4 <- ggplot(top_genes, aes(x = gene_id, y = rewiring_score)) +
    geom_bar(stat = "identity", fill = "#756bb1") +
    coord_flip() +
    labs(title = "Top 15 Rewiring Genes",
         x = "", y = "Rewiring Score") +
    theme_minimal()

  plot_grid(p1, p2, p3, p4, ncol = 2)
}

# =============================================================================
# Visualization: Qualitative Changes
# =============================================================================

plot_qualitative_changes <- function(data) {
  edges_diff <- data$edges_diff

  # Panel 1: Count by category
  qual_summary <- edges_diff %>%
    group_by(qual_label) %>%
    summarize(count = n(), mean_abs_delta = mean(abs(delta)), .groups = "drop") %>%
    mutate(qual_label = factor(qual_label, levels = names(QUAL_COLORS)))

  p1 <- ggplot(qual_summary, aes(x = reorder(qual_label, -count), y = count, fill = qual_label)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = QUAL_COLORS) +
    labs(title = "Qualitative Edge Changes",
         x = "Change Type", y = "Number of Edges") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  # Panel 2: Delta distribution by category
  p2 <- edges_diff %>%
    filter(qual_label != "unchanged") %>%
    ggplot(aes(x = qual_label, y = abs(delta), fill = qual_label)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = QUAL_COLORS) +
    labs(title = "Effect Size by Change Type",
         x = "Change Type", y = "|Delta r|") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  # Panel 3: r_low vs r_high colored by category
  p3 <- edges_diff %>%
    filter(qual_label != "unchanged") %>%
    ggplot(aes(x = r_low, y = r_high, color = qual_label)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    scale_color_manual(values = QUAL_COLORS) +
    labs(title = "Correlation Change Pattern",
         x = "r (Low)", y = "r (High)") +
    theme_minimal() +
    theme(legend.position = "right")

  plot_grid(p1, p2, p3, ncol = 2, rel_widths = c(1, 1))
}

# =============================================================================
# Main
# =============================================================================

main <- function() {
  # Parse arguments
  parser <- ArgumentParser(description = "Visualize differential co-expression networks")
  parser$add_argument("--data-dir", type = "character", default = "visualization_data",
                      help = "Directory with TSV files from Stage 5")
  parser$add_argument("--top-n", type = "integer", default = 10,
                      help = "Number of top rewiring genes to visualize")
  parser$add_argument("--output-prefix", type = "character", default = "networks",
                      help = "Prefix for output files")

  args <- parser$parse_args()

  # Load data
  data <- load_network_data(args$data_dir)

  # Identify top rewiring genes
  top_genes <- data$nodes %>%
    arrange(desc(rewiring_score)) %>%
    head(args$top_n) %>%
    pull(gene_id)

  cat("\nTop", args$top_n, "rewiring genes:\n")
  print(data$nodes %>%
          filter(gene_id %in% top_genes) %>%
          select(gene_id, degree_low, degree_high, degree_diff, rewiring_score))

  # ==========================================================================
  # Plot 1: Ego network comparisons
  # ==========================================================================
  cat("\nGenerating ego network comparisons...\n")
  pdf(paste0(args$output_prefix, "_ego_comparison.pdf"), width = 15, height = 5 * ceiling(length(top_genes) / 2))

  for (gene in top_genes) {
    cat("  Plotting ego network for:", gene, "\n")
    tryCatch({
      plot_ego_comparison(data, gene)

      # Add legend
      par(mfrow = c(1, 1))
      plot.new()
      legend("center",
             title = "Edge Colors (DIFF network)",
             legend = names(QUAL_COLORS),
             fill = QUAL_COLORS,
             ncol = 2,
             bty = "n")
    }, error = function(e) {
      cat("    Warning: Could not plot", gene, "-", conditionMessage(e), "\n")
    })
  }
  dev.off()
  cat("Saved:", paste0(args$output_prefix, "_ego_comparison.pdf"), "\n")

  # ==========================================================================
  # Plot 2: Degree analysis
  # ==========================================================================
  cat("\nGenerating degree analysis...\n")
  pdf(paste0(args$output_prefix, "_degree_analysis.pdf"), width = 12, height = 10)
  print(plot_degree_analysis(data))
  dev.off()
  cat("Saved:", paste0(args$output_prefix, "_degree_analysis.pdf"), "\n")

  # ==========================================================================
  # Plot 3: Qualitative changes
  # ==========================================================================
  cat("\nGenerating qualitative changes plot...\n")
  pdf(paste0(args$output_prefix, "_qualitative.pdf"), width = 12, height = 8)
  print(plot_qualitative_changes(data))
  dev.off()
  cat("Saved:", paste0(args$output_prefix, "_qualitative.pdf"), "\n")

  # ==========================================================================
  # Summary statistics
  # ==========================================================================
  cat("\nWriting summary statistics...\n")
  sink(paste0(args$output_prefix, "_summary.txt"))
  cat("=== Differential Network Analysis Summary ===\n\n")

  cat("Network sizes:\n")
  cat("  Low network:", nrow(data$edges_low), "edges\n")
  cat("  High network:", nrow(data$edges_high), "edges\n")
  cat("  Differential:", nrow(data$edges_diff), "edges\n")
  cat("  Active nodes:", nrow(data$nodes), "\n\n")

  cat("Qualitative changes:\n")
  qual_counts <- table(data$edges_diff$qual_label)
  for (q in names(qual_counts)) {
    cat(sprintf("  %s: %d (%.1f%%)\n", q, qual_counts[q],
                100 * qual_counts[q] / sum(qual_counts)))
  }

  cat("\nDegree statistics:\n")
  cat(sprintf("  Mean degree (low): %.2f\n", mean(data$nodes$degree_low)))
  cat(sprintf("  Mean degree (high): %.2f\n", mean(data$nodes$degree_high)))
  cat(sprintf("  Mean degree (diff): %.2f\n", mean(data$nodes$degree_diff)))

  cat("\nTop 20 rewiring genes:\n")
  print(data$nodes %>%
          arrange(desc(rewiring_score)) %>%
          head(20) %>%
          select(gene_id, degree_low, degree_high, degree_diff, rewiring_score,
                 n_disappear, n_new, n_sign_change))

  sink()
  cat("Saved:", paste0(args$output_prefix, "_summary.txt"), "\n")

  cat("\n=== Visualization complete! ===\n")
}

# Run main
main()
