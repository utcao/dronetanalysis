#!/usr/bin/env Rscript
# ==============================================================================
# Network Visualisation using igraph
#
# Script by Gabriel Thornes
#
# Last Updated: 13/11/2025
#
# This script:
#   1. Loads TOM matrix and gene metrics
#   2. Creates igraph network object
#   3. Generates full network overview Visualisations
#   4. Exports high-quality plots
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
  library(argparse)
})

# ----- 2. Command-line arguments -----
parser <- ArgumentParser(description = 'Visualize gene co-expression network')
parser$add_argument('--tom-file', help = 'TOM matrix CSV file', 
                   default = 'results/network_features/features_calc/adjacency/tom_adjacency_matrix.csv')
parser$add_argument('--metrics-file', help = 'Gene metrics + expression CSV', 
                   default = 'results/network_features/gene_metrics/adjacency_gene_metrics_with_expression.csv')
parser$add_argument('--output-dir', help = 'Directory to save plots', 
                   default = 'results/analysis/network_vis/Control_modules/ctrl')
parser$add_argument('--threshold', type = 'double', help = 'TOM threshold for edge filtering', 
                   default = 0.2)
parser$add_argument('--layout', help = 'Layout algorithm: fr (Fruchterman-Reingold), kk (Kamada-Kawai), drl', 
                   default = 'kk')
parser$add_argument('--seed', type = 'integer', help = 'Random seed for reproducible layouts', 
                   default = 42)
args <- parser$parse_args()

output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Network Visualisation ===\n")
cat("TOM file:", args$tom_file, "\n")
cat("Metrics file:", args$metrics_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("TOM threshold:", args$threshold, "\n")
cat("Layout algorithm:", args$layout, "\n\n")

set.seed(args$seed)

# ----- 3. Load data -----
cat("Loading TOM matrix...\n")
tom_matrix <- read.csv(args$tom_file, row.names = 1, check.names = FALSE)
cat("  TOM dimensions:", nrow(tom_matrix), "x", ncol(tom_matrix), "\n")
cat("  TOM range:", round(min(tom_matrix, na.rm = TRUE), 4), "to", round(max(tom_matrix, na.rm = TRUE), 4), "\n")
cat("  TOM quantiles:\n")
print(round(quantile(as.vector(as.matrix(tom_matrix)), probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm = TRUE), 4))
cat("\n")

cat("Loading gene metrics...\n")
gene_metrics <- read.csv(args$metrics_file, stringsAsFactors = FALSE)
cat("  Genes:", nrow(gene_metrics), "\n")
cat("  Modules:", length(unique(gene_metrics$module)), "\n\n")

# ----- 4. Apply threshold to TOM matrix -----
cat("Applying threshold to TOM matrix...\n")
cat("  Total possible edges:", nrow(tom_matrix) * (nrow(tom_matrix) - 1) / 2, "\n")
cat("  Non-zero edges before threshold:", sum(tom_matrix > 0, na.rm = TRUE) / 2, "\n")

# Set diagonal to 0
diag(tom_matrix) <- 0

tom_matrix[tom_matrix < args$threshold] <- 0
cat("  Edges after threshold (", args$threshold, "):", sum(tom_matrix > 0, na.rm = TRUE) / 2, "\n")

if (sum(tom_matrix > 0) == 0) {
  cat("\n*** WARNING: No edges remain after threshold! ***\n")
  cat("*** Try a lower threshold. Suggested thresholds based on your data:\n")
  cat("***   95th percentile:", round(quantile(as.vector(as.matrix(tom_matrix)), 0.95, na.rm = TRUE), 4), "\n")
  cat("***   90th percentile:", round(quantile(as.vector(as.matrix(tom_matrix)), 0.90, na.rm = TRUE), 4), "\n")
  cat("***   75th percentile:", round(quantile(as.vector(as.matrix(tom_matrix)), 0.75, na.rm = TRUE), 4), "\n")
  stop("No edges in network. Please adjust --threshold parameter.")
}
cat("\n")

# ----- 5. Create igraph network -----
cat("Creating igraph network object...\n")
g <- graph_from_adjacency_matrix(
  as.matrix(tom_matrix),
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

cat("  Nodes:", vcount(g), "\n")
cat("  Edges:", ecount(g), "\n")
cat("  Density:", round(edge_density(g), 6), "\n\n")

# ----- 6. Add node attributes -----
cat("Adding node attributes...\n")

# Match gene metrics to network nodes
node_genes <- V(g)$name
gene_data <- gene_metrics %>%
  filter(gene %in% node_genes) %>%
  arrange(match(gene, node_genes))

# Module assignment
V(g)$module <- gene_data$module[match(V(g)$name, gene_data$gene)]

# Network metrics
V(g)$connectivity <- gene_data$weighted_connectivity[match(V(g)$name, gene_data$gene)]
V(g)$degree <- gene_data$degree[match(V(g)$name, gene_data$gene)]
V(g)$betweenness <- gene_data$betweenness_centrality[match(V(g)$name, gene_data$gene)]
V(g)$eigenvector <- gene_data$eigenvector_centrality[match(V(g)$name, gene_data$gene)]

# WGCNA module membership metrics
V(g)$kWithin <- gene_data$kWithin[match(V(g)$name, gene_data$gene)]
V(g)$kDiff <- gene_data$kDiff[match(V(g)$name, gene_data$gene)]

# Expression
V(g)$expression <- gene_data$mean_expression[match(V(g)$name, gene_data$gene)]
V(g)$variance <- gene_data$variance_expression[match(V(g)$name, gene_data$gene)]

cat("  Attributes added: module, connectivity, degree, betweenness, eigenvector, kWithin, kDiff, expression, variance\n\n")

# ----- 7. Set up visual properties -----
cat("Setting up visual properties...\n")

# Module colors (use module name as color where valid)
all_modules <- unique(V(g)$module)
module_color_map <- setNames(all_modules, all_modules)
for (mod in all_modules) {
  if (mod %in% colors() && mod != "grey") {
    module_color_map[mod] <- mod
  }
}

V(g)$color <- module_color_map[V(g)$module]
V(g)$frame.color <- NA  # No border

# Node size based on connectivity (range 2-15)
V(g)$size <- rescale(V(g)$connectivity, to = c(2, 15))

# Edge properties - make edges more visible
E(g)$width <- rescale(E(g)$weight, to = c(0.3, 3))
E(g)$color <- rgb(0.3, 0.3, 0.3, alpha = 0.5)  # Darker semi-transparent grey

cat("  Node colors: by module\n")
cat("  Node sizes: by weighted connectivity\n")
cat("  Edge widths: by TOM weight\n")
cat("  Edge transparency: 50% for clarity\n\n")

# ----- 8. Compute layout -----
cat("Computing network layout (", args$layout, ")...\n")
cat("  This may take a few minutes for large networks...\n")

if (args$layout == "fr") {
  layout_coords <- layout_with_fr(g)
} else if (args$layout == "kk") {
  layout_coords <- layout_with_kk(g)
} else if (args$layout == "drl") {
  layout_coords <- layout_with_drl(g)
} else {
  cat("  Warning: Unknown layout '", args$layout, "', using Fruchterman-Reingold\n")
  layout_coords <- layout_with_fr(g)
}

cat("  Layout computed.\n\n")

# ----- 9. Plot full network overview -----
cat("Generating full network overview...\n")

pdf(file.path(output_dir, "01_full_network_overview.pdf"), width = 16, height = 16)
par(mar = c(0, 0, 2, 0))
plot(g,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = V(g)$size,
     vertex.color = V(g)$color,
     vertex.frame.color = V(g)$frame.color,
     edge.width = E(g)$width,
     edge.color = E(g)$color,
     main = "Full Network Overview - Colored by Module")

# Add legend
legend("topright",
       legend = sort(unique(V(g)$module)),
       pch = 21,
       pt.bg = module_color_map[sort(unique(V(g)$module))],
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       title = "Module")
dev.off()
cat("  Saved: 01_full_network_overview.pdf\n")

# ----- 10. Plot network colored by connectivity -----
cat("Generating network colored by connectivity...\n")

# Create color gradient for connectivity
connectivity_colors <- colorRampPalette(c("lightblue", "yellow", "red"))(100)
V(g)$connectivity_color <- connectivity_colors[cut(V(g)$connectivity, breaks = 100)]

pdf(file.path(output_dir, "02_network_by_connectivity.pdf"), width = 16, height = 16)
par(mar = c(0, 0, 2, 0))
plot(g,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = V(g)$size,
     vertex.color = V(g)$connectivity_color,
     vertex.frame.color = NA,
     edge.width = E(g)$width,
     edge.color = E(g)$color,
     main = "Network Colored by Weighted Connectivity")

# Add color scale legend
legend_colors <- connectivity_colors[seq(1, 100, length.out = 5)]
legend_labels <- round(quantile(V(g)$connectivity, probs = seq(0, 1, 0.25)), 2)
legend("topright",
       legend = legend_labels,
       pch = 21,
       pt.bg = legend_colors,
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       title = "Connectivity")
dev.off()
cat("  Saved: 02_network_by_connectivity.pdf\n")

# ----- 11. Plot network colored by expression -----
cat("Generating network colored by expression...\n")

# Create color gradient for expression (blue-white-red)
expr_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
V(g)$expression_color <- expr_colors[cut(V(g)$expression, breaks = 100)]

pdf(file.path(output_dir, "03_network_by_expression.pdf"), width = 16, height = 16)
par(mar = c(0, 0, 2, 0))
plot(g,
     layout = layout_coords,
     vertex.label = NA,
     vertex.size = V(g)$size,
     vertex.color = V(g)$expression_color,
     vertex.frame.color = NA,
     edge.width = E(g)$width,
     edge.color = E(g)$color,
     main = "Network Colored by Mean Expression Level")

# Add color scale legend
legend_expr <- round(quantile(V(g)$expression, probs = seq(0, 1, 0.25), na.rm = TRUE), 2)
legend("topright",
       legend = legend_expr,
       pch = 21,
       pt.bg = expr_colors[seq(1, 100, length.out = 5)],
       pt.cex = 2,
       cex = 0.8,
       bty = "n",
       title = "Expression")
dev.off()
cat("  Saved: 03_network_by_expression.pdf\n")

# ----- 12. Highlight hub genes -----
cat("Generating network with hub genes highlighted...\n")

# Define hubs as top 5% by weighted connectivity (global modifiers)
# Weighted connectivity is the primary determinant of global modifier status
hub_threshold <- quantile(V(g)$connectivity, 0.95, na.rm = TRUE)
V(g)$is_hub <- V(g)$connectivity >= hub_threshold

# Hub gene names
hub_genes <- V(g)$name[V(g)$is_hub]
cat("  Hub genes (top 5% by weighted connectivity >=", round(hub_threshold, 3), "):", length(hub_genes), "\n")

# Color: hubs in red, others by module with transparency
V(g)$hub_color <- ifelse(V(g)$is_hub, "red", 
                          adjustcolor(V(g)$color, alpha.f = 0.6))

# Size: hubs larger
V(g)$hub_size <- ifelse(V(g)$is_hub, 12, V(g)$size)

pdf(file.path(output_dir, "04_network_hub_genes.pdf"), width = 16, height = 16)
par(mar = c(0, 0, 2, 0))
plot(g,
     layout = layout_coords,
     vertex.label = ifelse(V(g)$is_hub, V(g)$name, NA),
     vertex.label.cex = 0.5,
     vertex.label.color = "black",
     vertex.size = V(g)$hub_size,
     vertex.color = V(g)$hub_color,
     vertex.frame.color = ifelse(V(g)$is_hub, "black", NA),
     edge.width = E(g)$width,
     edge.color = E(g)$color,
     main = "Network with Hub Genes Highlighted (Top 5% by Weighted Connectivity)")

legend("topright",
       legend = c("Hub genes (global modifiers)", "Other genes"),
       pch = 21,
       pt.bg = c("red", "white"),
       pt.cex = 2,
       cex = 0.8,
       bty = "n")
dev.off()
cat("  Saved: 04_network_hub_genes.pdf\n")

# ----- 13. Network statistics summary -----
cat("\n=== Network Statistics ===\n")

cat("\nBasic Properties:\n")
cat("  Nodes:", vcount(g), "\n")
cat("  Edges:", ecount(g), "\n")
cat("  Density:", round(edge_density(g), 6), "\n")
cat("  Average degree:", round(mean(degree(g)), 2), "\n")
cat("  Transitivity (clustering):", round(transitivity(g), 4), "\n")

cat("\nConnectivity:\n")
cat("  Is connected:", is_connected(g), "\n")
cat("  Number of components:", components(g)$no, "\n")
cat("  Size of largest component:", max(components(g)$csize), "\n")

cat("\nDegree Distribution:\n")
degree_stats <- summary(degree(g))
print(degree_stats)

cat("\nModule Sizes:\n")
module_sizes <- table(V(g)$module)
print(sort(module_sizes, decreasing = TRUE))

# Save network object for later use
cat("\nSaving network object...\n")
saveRDS(g, file.path(output_dir, "network_object.rds"))
cat("  Saved: network_object.rds\n")

# Save layout coordinates
cat("Saving layout coordinates...\n")
layout_df <- data.frame(
  gene = V(g)$name,
  x = layout_coords[, 1],
  y = layout_coords[, 2]
)
write.csv(layout_df, file.path(output_dir, "network_layout_coordinates.csv"), row.names = FALSE)
cat("  Saved: network_layout_coordinates.csv\n")

# Save hub genes list
cat("Saving hub genes list...\n")
hub_df <- data.frame(
  gene = hub_genes,
  connectivity = V(g)$connectivity[V(g)$is_hub],
  kWithin = V(g)$kWithin[V(g)$is_hub],
  kDiff = V(g)$kDiff[V(g)$is_hub],
  eigenvector = V(g)$eigenvector[V(g)$is_hub],
  module = V(g)$module[V(g)$is_hub],
  expression = V(g)$expression[V(g)$is_hub]
) %>%
  arrange(desc(connectivity))
write.csv(hub_df, file.path(output_dir, "hub_genes.csv"), row.names = FALSE)
cat("  Saved: hub_genes.csv\n")

cat("\n=== Visualisation Complete ===\n")
cat("All plots saved to:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("  - 01_full_network_overview.pdf (colored by module)\n")
cat("  - 02_network_by_connectivity.pdf (colored by connectivity)\n")
cat("  - 03_network_by_expression.pdf (colored by expression)\n")
cat("  - 04_network_hub_genes.pdf (hub genes highlighted)\n")
cat("  - network_object.rds (saved igraph object)\n")
cat("  - network_layout_coordinates.csv (x,y positions for reproducibility)\n")
cat("  - hub_genes.csv (list of hub genes)\n")
