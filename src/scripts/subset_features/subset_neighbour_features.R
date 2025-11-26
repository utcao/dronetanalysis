#!/usr/bin/env Rscript
# ==============================================================================
# Extract Network Metrics for Gene neighbours from Differential Networks
#
# Script by Gabriel Thornes
#
# This script:
#   1. Loads differential adjacency matrix
#   2. Identifies neighbours of a focal gene at specified layer(s)
#   3. Extracts network metrics for those neighbours
#   4. Saves neighbour-specific feature tables
# ==============================================================================

rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(argparse)
})

# Command-line arguments
parser <- ArgumentParser(description = 'Extract neighbour features from differential networks')
parser$add_argument('--focal_gene', help = 'Focal gene to find neighbours for',
                   required = TRUE)
parser$add_argument('--diff_adj_matrix', help = 'Path to differential adjacency matrix',
                   required = TRUE)
parser$add_argument('--gene_metrics_diff', help = 'Path to gene metrics from differential network',
                   required = TRUE)
parser$add_argument('--output_dir', help = 'Directory to save neighbour features',
                   required = TRUE)
parser$add_argument('--neighbour_layers', help = 'neighbour layers to extract (comma-separated, e.g., "1,2")',
                   default = "1,2")
args <- parser$parse_args()

# Parse neighbour layers
neighbour_layers <- as.integer(strsplit(args$neighbour_layers, ",")[[1]])

cat("=== neighbour Feature Extraction ===\n")
cat("Focal gene:", args$focal_gene, "\n")
cat("Differential adjacency matrix:", args$diff_adj_matrix, "\n")
cat("Gene metrics (differential network):", args$gene_metrics_diff, "\n")
cat("Output directory:", args$output_dir, "\n")
cat("neighbour layers:", paste(neighbour_layers, collapse = ", "), "\n")

# Create output directory
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Load differential adjacency matrix
diff_adj <- fread(args$diff_adj_matrix, data.table = FALSE)
gene_names <- as.character(diff_adj[[1]])
diff_adj_mat <- as.matrix(diff_adj[, -1])
rownames(diff_adj_mat) <- gene_names
colnames(diff_adj_mat) <- gene_names

cat("Loaded differential adjacency matrix:", nrow(diff_adj_mat), "x", ncol(diff_adj_mat), "\n")

# Check if focal gene exists in matrix
if (!args$focal_gene %in% gene_names) {
  stop("ERROR: Focal gene '", args$focal_gene, "' not found in adjacency matrix!")
}

# Load gene metrics from differential network
metrics_diff <- fread(args$gene_metrics_diff, data.table = FALSE)

cat("Loaded gene metrics:\n")
cat("  Differential network:", nrow(metrics_diff), "genes\n\n")

# Create igraph object for neighbour finding
g <- graph_from_adjacency_matrix(diff_adj_mat, mode = "undirected", diag = FALSE)
V(g)$name <- gene_names

cat("  Vertices:", vcount(g), "\n")
cat("  Edges:", ecount(g), "\n")
cat("  Density:", round(edge_density(g), 4), "\n\n")

# Function to get neighbours at specific layer
get_neighbours_at_layer <- function(graph, focal_node, layer) {
  if (layer == 1) {
    # Direct neighbours
    return(neighbors(graph, focal_node))
  } else {
    # Get all nodes within 'layer' distance, then subtract closer layers
    all_within <- ego(graph, order = layer, nodes = focal_node)[[1]]
    if (layer > 1) {
      closer <- ego(graph, order = layer - 1, nodes = focal_node)[[1]]
      layer_neighbours <- setdiff(all_within, closer)
    } else {
      layer_neighbours <- setdiff(all_within, focal_node)
    }
    return(layer_neighbours)
  }
}

# Find focal gene vertex
focal_vertex <- which(V(g)$name == args$focal_gene)

if (length(focal_vertex) == 0) {
  stop("ERROR: Could not find focal gene vertex in graph!")
}

# Extract neighbours for each layer
results_list <- list()

for (layer in neighbour_layers) {
  cat("=== Processing Layer", layer, "neighbours ===\n")
  
  # Get neighbours at this layer
  neighbour_vertices <- get_neighbours_at_layer(g, focal_vertex, layer)
  neighbour_genes <- V(g)$name[neighbour_vertices]
  
  cat("Found", length(neighbour_genes), "neighbours at layer", layer, "\n")
  
  if (length(neighbour_genes) == 0) {
    cat("No neighbours found at layer", layer, ". Skipping.\n\n")
    next
  }
  
  # Get differential adjacency values for these neighbours
  diff_values <- diff_adj_mat[args$focal_gene, neighbour_genes]
  
  # Extract metrics for these neighbours from differential network
  metrics_subset <- metrics_diff[metrics_diff$gene %in% neighbour_genes, ]
  
  # Merge data
  neighbour_data <- data.frame(
    neighbour_gene = neighbour_genes,
    layer = layer,
    diff_adjacency_to_focal = diff_values,
    stringsAsFactors = FALSE
  )
  
  # Add metrics from differential network
  if (nrow(metrics_subset) > 0) {
    metrics_subset <- metrics_subset[match(neighbour_genes, metrics_subset$gene), ]
    neighbour_data <- cbind(neighbour_data, metrics_subset[, -1])
  }
  
  # Sort by differential adjacency (largest changes first)
  neighbour_data <- neighbour_data[order(-abs(neighbour_data$diff_adjacency_to_focal)), ]
  
  # Save results
  output_file <- file.path(args$output_dir, 
                           paste0(args$focal_gene, "_layer", layer, "_neighbours.csv"))
  write.csv(neighbour_data, output_file, row.names = FALSE)
  cat("Saved layer", layer, "neighbours to:", output_file, "\n")
  
  # Print summary
  cat("\nLayer", layer, "Summary:\n")
  cat("  Total neighbours:", nrow(neighbour_data), "\n")
  cat("  Mean diff adjacency:", round(mean(neighbour_data$diff_adjacency_to_focal, na.rm = TRUE), 4), "\n")
  cat("  Top 5 neighbours by diff adjacency:\n")
  print(neighbour_data[1:min(5, nrow(neighbour_data)), c("neighbour_gene", "diff_adjacency_to_focal")])
  cat("\n")
  
  results_list[[paste0("layer_", layer)]] <- neighbour_data
}

# Create combined summary across all layers
if (length(results_list) > 0) {
  combined_data <- do.call(rbind, results_list)
  summary_file <- file.path(args$output_dir, 
                            paste0(args$focal_gene, "_all_layers_neighbours.csv"))
  write.csv(combined_data, summary_file, row.names = FALSE)
  cat("Combined results saved to:", summary_file, "\n\n")
}

# Summary statistics
cat("=== Overall Summary ===\n")
cat("Focal gene:", args$focal_gene, "\n")
for (layer in neighbour_layers) {
  layer_key <- paste0("layer_", layer)
  if (layer_key %in% names(results_list)) {
    n_neighbours <- nrow(results_list[[layer_key]])
    cat("Layer", layer, ":", n_neighbours, "neighbours\n")
  } else {
    cat("Layer", layer, ": 0 neighbours\n")
  }
}

cat("\nneighbour feature extraction complete!\n")
