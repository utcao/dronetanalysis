#!/usr/bin/env Rscript
# ==============================================================================
# Module Preservation Analysis using NetRep
#
# Script by Gabriel Thornes
#
# Last Updated: 20/11/2025
#
# This script:
#   1. Loads expression data for both conditions
#   2. Uses Control module assignments as reference
#   3. Calculates module preservation using NetRep (much faster than WGCNA)
#   4. Generates preservation plots and summary tables
#
# NetRep advantages over WGCNA modulePreservation:
#   - 10-100x faster (uses permutation-based null)
#   - More memory efficient
#   - Better for large networks
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(NetRep)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(argparse)
})

# ----- 2. Command-line arguments -----
parser <- ArgumentParser(description = 'Module preservation analysis using NetRep')
parser$add_argument('--ctrl-expr', help = 'Control expression data (genes x samples)', 
                   default = 'dataset/processed/VOOM/voomdataCtrl.txt')
parser$add_argument('--treat-expr', help = 'Treatment expression data (genes x samples)', 
                   default = 'dataset/processed/VOOM/voomdataHS.txt')
parser$add_argument('--ctrl-modules', help = 'Control module assignments CSV', 
                   default = 'results/network_features/features_calc/adjacency/modules/gene_connectivity.csv')
parser$add_argument('--output-dir', help = 'Directory to save results', 
                   default = 'results/analysis/module_overlap/module_preservation_netrep')
parser$add_argument('--nPermutations', type = 'integer', help = 'Number of permutations for significance testing', 
                   default = 10000)
parser$add_argument('--nThreads', type = 'integer', help = 'Number of parallel threads', 
                   default = 4)
parser$add_argument('--corType', help = 'Correlation type: pearson or bicor', 
                   default = 'bicor')
args <- parser$parse_args()

output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Module Preservation Analysis (NetRep) ===\n")
cat("Control expression:", args$ctrl_expr, "\n")
cat("Treatment expression:", args$treat_expr, "\n")
cat("Control modules:", args$ctrl_modules, "\n")
cat("Output directory:", output_dir, "\n")
cat("Permutations:", args$nPermutations, "\n")
cat("Threads:", args$nThreads, "\n")
cat("Correlation type:", args$corType, "\n\n")

# ----- 3. Load data -----
cat("Loading expression data...\n")

# Control expression - first column is gene names without header
ctrl_expr <- fread(args$ctrl_expr, data.table = FALSE)
# fread adds V1 as default name when first column has no header
if (colnames(ctrl_expr)[1] %in% c("V1", "")) {
  colnames(ctrl_expr)[1] <- "gene"
}
ctrl_genes <- as.character(ctrl_expr[[1]])  # Ensure character vector
# Extract numeric columns only
ctrl_expr_numeric <- ctrl_expr[, -1, drop = FALSE]
# Convert to numeric matrix explicitly
ctrl_expr_mat <- matrix(
  as.numeric(as.matrix(ctrl_expr_numeric)),
  nrow = nrow(ctrl_expr_numeric),
  ncol = ncol(ctrl_expr_numeric)
)
rownames(ctrl_expr_mat) <- ctrl_genes
colnames(ctrl_expr_mat) <- colnames(ctrl_expr_numeric)
ctrl_expr_mat <- t(ctrl_expr_mat)  # NetRep expects samples x genes
storage.mode(ctrl_expr_mat) <- "double"  # Ensure numeric storage
cat("  Control:", nrow(ctrl_expr_mat), "samples x", ncol(ctrl_expr_mat), "genes\n")

# Treatment expression - first column is gene names without header
treat_expr <- fread(args$treat_expr, data.table = FALSE)
if (colnames(treat_expr)[1] %in% c("V1", "")) {
  colnames(treat_expr)[1] <- "gene"
}
treat_genes <- as.character(treat_expr[[1]])  # Ensure character vector
# Extract numeric columns only
treat_expr_numeric <- treat_expr[, -1, drop = FALSE]
# Convert to numeric matrix explicitly
treat_expr_mat <- matrix(
  as.numeric(as.matrix(treat_expr_numeric)),
  nrow = nrow(treat_expr_numeric),
  ncol = ncol(treat_expr_numeric)
)
rownames(treat_expr_mat) <- treat_genes
colnames(treat_expr_mat) <- colnames(treat_expr_numeric)
treat_expr_mat <- t(treat_expr_mat)  # NetRep expects samples x genes
storage.mode(treat_expr_mat) <- "double"  # Ensure numeric storage
cat("  Treatment:", nrow(treat_expr_mat), "samples x", ncol(treat_expr_mat), "genes\n")

# Ensure same genes in both datasets
common_genes <- intersect(colnames(ctrl_expr_mat), colnames(treat_expr_mat))
cat("  Common genes:", length(common_genes), "\n")

ctrl_expr_mat <- ctrl_expr_mat[, common_genes]
treat_expr_mat <- treat_expr_mat[, common_genes]

# Load control modules
cat("\nLoading control module assignments...\n")
ctrl_modules <- read.csv(args$ctrl_modules, stringsAsFactors = FALSE)
cat("  Loaded:", nrow(ctrl_modules), "genes with module assignments\n")

# Match modules to expression data
module_assignments <- ctrl_modules$module[match(common_genes, ctrl_modules$gene)]
names(module_assignments) <- common_genes

# Remove genes without module assignment
genes_with_modules <- !is.na(module_assignments)
ctrl_expr_mat <- ctrl_expr_mat[, genes_with_modules]
treat_expr_mat <- treat_expr_mat[, genes_with_modules]
module_assignments <- module_assignments[genes_with_modules]

cat("  Genes with module assignments:", length(module_assignments), "\n")
cat("  Modules:", length(unique(module_assignments)), "\n")

# Module size distribution
module_sizes <- table(module_assignments)
cat("\nModule sizes:\n")
print(sort(module_sizes, decreasing = TRUE))

# ----- 4. Calculate correlation matrices -----
cat("\nCalculating correlation matrices...\n")
cat("  This may take a few minutes...\n")

if (args$corType == "bicor") {
  # Biweight midcorrelation (more robust to outliers)
  cat("  Using biweight midcorrelation...\n")
  ctrl_cor <- WGCNA::bicor(ctrl_expr_mat, use = "pairwise.complete.obs")
  treat_cor <- WGCNA::bicor(treat_expr_mat, use = "pairwise.complete.obs")
} else {
  # Pearson correlation
  cat("  Using Pearson correlation...\n")
  ctrl_cor <- cor(ctrl_expr_mat, use = "pairwise.complete.obs")
  treat_cor <- cor(treat_expr_mat, use = "pairwise.complete.obs")
}

# Ensure matrices are proper numeric matrices (not data frames)
ctrl_cor <- matrix(as.numeric(ctrl_cor), nrow = nrow(ctrl_cor), ncol = ncol(ctrl_cor))
rownames(ctrl_cor) <- colnames(ctrl_expr_mat)
colnames(ctrl_cor) <- colnames(ctrl_expr_mat)
storage.mode(ctrl_cor) <- "double"

treat_cor <- matrix(as.numeric(treat_cor), nrow = nrow(treat_cor), ncol = ncol(treat_cor))
rownames(treat_cor) <- colnames(treat_expr_mat)
colnames(treat_cor) <- colnames(treat_expr_mat)
storage.mode(treat_cor) <- "double"

# Ensure diagonal is exactly 1 (NetRep requirement)
diag(ctrl_cor) <- 1
diag(treat_cor) <- 1

cat("  Correlation matrices calculated\n")
cat("  Control correlation matrix:", nrow(ctrl_cor), "x", ncol(ctrl_cor), "\n")
cat("  Treatment correlation matrix:", nrow(treat_cor), "x", ncol(treat_cor), "\n\n")

# ----- 5. Create adjacency matrices (networks) -----
cat("Creating adjacency matrices from correlations...\n")
# NetRep expects both correlation AND network (adjacency) matrices
# For weighted networks, we can use the absolute correlations as weights
# Or use soft-thresholding as in WGCNA

# Simple approach: use absolute correlation as adjacency
ctrl_adj <- abs(ctrl_cor)
treat_adj <- abs(treat_cor)

# Ensure proper matrix format
storage.mode(ctrl_adj) <- "double"
storage.mode(treat_adj) <- "double"

cat("  Adjacency matrices created\n\n")

# ----- 5. Run NetRep module preservation analysis -----
cat("Running NetRep module preservation analysis...\n")
cat("  Expected runtime: 5-15 minutes with", args$nPermutations, "permutations\n\n")

# Convert module assignments to list format for NetRep
module_list <- split(names(module_assignments), module_assignments)

# Remove grey module if present (unassigned genes)
if ("grey" %in% names(module_list)) {
  cat("  Removing grey module (unassigned genes) from analysis\n")
  module_list <- module_list[names(module_list) != "grey"]
}

cat("  Testing", length(module_list), "modules\n\n")

# Run module preservation test
# discovery = Control (reference), test = Treatment
# NetRep requires: network (adjacency), data (expression), and correlation matrices
preservation_result <- modulePreservation(
  network = list(
    discovery = ctrl_adj,
    test = treat_adj
  ),
  data = list(
    discovery = ctrl_expr_mat,
    test = treat_expr_mat
  ),
  correlation = list(
    discovery = ctrl_cor,
    test = treat_cor
  ),
  moduleAssignments = list(discovery = module_assignments),
  modules = names(module_list),
  nPerm = args$nPermutations,
  nThreads = args$nThreads,
  verbose = TRUE
)

cat("\nModule preservation analysis complete!\n\n")

# ----- 6. Extract and save results -----
cat("Extracting preservation statistics...\n")

# Check structure of result
cat("  Result structure:\n")
print(str(preservation_result, max.level = 2))

# NetRep returns a list with different structure depending on simplify parameter
# Extract the actual data - it should be in preservation_result[[1]] or similar
if (is.list(preservation_result) && !is.null(names(preservation_result))) {
  # If result has named elements
  pres_data <- preservation_result
} else if (is.list(preservation_result) && length(preservation_result) == 1) {
  # If result is a list with one element
  pres_data <- preservation_result[[1]]
} else {
  pres_data <- preservation_result
}

cat("\n  Extracted data structure:\n")
print(str(pres_data, max.level = 2))

# Get preservation statistics for each module
# NetRep returns preservation and p.values as matrices or data frames
preservation_stats <- data.frame(
  module = names(module_list),
  moduleSize = sapply(module_list, length),
  stringsAsFactors = FALSE
)

# Add preservation statistics (check if they exist)
if (!is.null(pres_data$preservation)) {
  if (is.matrix(pres_data$preservation) || is.data.frame(pres_data$preservation)) {
    preservation_stats$avg_weight_preservation <- pres_data$preservation[, "avg.weight"]
    preservation_stats$avg_cor_preservation <- pres_data$preservation[, "avg.cor"]
    preservation_stats$avg_contrib_preservation <- pres_data$preservation[, "avg.contrib"]
    preservation_stats$coherence_preservation <- pres_data$preservation[, "coherence"]
  } else {
    # Try accessing as list
    preservation_stats$avg_weight_preservation <- pres_data$preservation$avg.weight
    preservation_stats$avg_cor_preservation <- pres_data$preservation$avg.cor
    preservation_stats$avg_contrib_preservation <- pres_data$preservation$avg.contrib
    preservation_stats$coherence_preservation <- pres_data$preservation$coherence
  }
}

# Add p-values (check if they exist)
if (!is.null(pres_data$p.values)) {
  if (is.matrix(pres_data$p.values) || is.data.frame(pres_data$p.values)) {
    preservation_stats$avg_weight_pval <- pres_data$p.values[, "avg.weight"]
    preservation_stats$avg_cor_pval <- pres_data$p.values[, "avg.cor"]
    preservation_stats$avg_contrib_pval <- pres_data$p.values[, "avg.contrib"]
    preservation_stats$coherence_pval <- pres_data$p.values[, "coherence"]
  } else {
    # Try accessing as list
    preservation_stats$avg_weight_pval <- pres_data$p.values$avg.weight
    preservation_stats$avg_cor_pval <- pres_data$p.values$avg.cor
    preservation_stats$avg_contrib_pval <- pres_data$p.values$avg.contrib
    preservation_stats$coherence_pval <- pres_data$p.values$coherence
  }
}

preservation_stats <- preservation_stats %>%
  arrange(avg_weight_pval)

# Add FDR correction
preservation_stats$avg_weight_fdr <- p.adjust(preservation_stats$avg_weight_pval, method = "fdr")
preservation_stats$coherence_fdr <- p.adjust(preservation_stats$coherence_pval, method = "fdr")

# Add interpretation based on p-values
preservation_stats <- preservation_stats %>%
  mutate(
    preservation_level = case_when(
      avg_weight_pval < 0.001 & coherence_pval < 0.001 ~ "Strong preservation",
      avg_weight_pval < 0.05 & coherence_pval < 0.05 ~ "Moderate preservation",
      avg_weight_pval < 0.1 | coherence_pval < 0.1 ~ "Weak preservation",
      TRUE ~ "No preservation"
    )
  )

cat("\nPreservation Summary:\n")
print(preservation_stats)

# Save results
write.csv(preservation_stats, 
          file.path(output_dir, "netrep_preservation_statistics.csv"), 
          row.names = FALSE)
cat("\nSaved: netrep_preservation_statistics.csv\n")

# Save full preservation object
saveRDS(preservation_result, file.path(output_dir, "netrep_preservation_full_results.rds"))
cat("Saved: netrep_preservation_full_results.rds\n")

# ----- 7. Generate plots -----
cat("\nGenerating preservation plots...\n")

# 7a. Preservation statistics plot (avg weight vs coherence)
pdf(file.path(output_dir, "01_preservation_scatter.pdf"), width = 10, height = 8)

plot_data <- preservation_stats

par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
plot(plot_data$avg_weight_preservation, plot_data$coherence_preservation,
     pch = 21, bg = plot_data$module,
     cex = 2,
     xlab = "Average Edge Weight Preservation",
     ylab = "Module Coherence Preservation",
     main = "Module Preservation: Control → High Sugar\n(NetRep Analysis)",
     cex.lab = 1.2, cex.axis = 1.1)

# Add text labels for modules
text(plot_data$avg_weight_preservation, plot_data$coherence_preservation,
     labels = plot_data$module, pos = 3, cex = 0.7)

# Add grid
grid()

dev.off()
cat("  Saved: 01_preservation_scatter.pdf\n")

# 7b. P-value plot (-log10)
pdf(file.path(output_dir, "02_preservation_pvalues.pdf"), width = 10, height = 8)

plot_data_pval <- plot_data %>%
  mutate(
    log10p_weight = -log10(avg_weight_pval + 1e-300),
    log10p_coherence = -log10(coherence_pval + 1e-300)
  )

par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
plot(plot_data_pval$log10p_weight, plot_data_pval$log10p_coherence,
     pch = 21, bg = plot_data_pval$module,
     cex = 2,
     xlab = "-log10(p-value) Weight Preservation",
     ylab = "-log10(p-value) Coherence Preservation",
     main = "Module Preservation Significance",
     cex.lab = 1.2, cex.axis = 1.1)

# Add significance thresholds
abline(h = -log10(0.05), v = -log10(0.05), lty = 2, col = "blue", lwd = 2)
abline(h = -log10(0.001), v = -log10(0.001), lty = 2, col = "darkgreen", lwd = 2)

text(plot_data_pval$log10p_weight, plot_data_pval$log10p_coherence,
     labels = plot_data_pval$module, pos = 3, cex = 0.7)

legend("bottomright",
       legend = c("p < 0.001", "p < 0.05"),
       lty = 2, col = c("darkgreen", "blue"),
       lwd = 2, cex = 0.9, bty = "n")

dev.off()
cat("  Saved: 02_preservation_pvalues.pdf\n")

# 7c. Bar plot of preservation levels
pdf(file.path(output_dir, "03_preservation_levels_barplot.pdf"), width = 12, height = 6)

preservation_summary <- preservation_stats %>%
  mutate(module = reorder(module, avg_weight_preservation))

p <- ggplot(preservation_summary, 
            aes(x = module, y = avg_weight_preservation, fill = preservation_level)) +
  geom_col() +
  scale_fill_manual(values = c("Strong preservation" = "#2E7D32",
                               "Moderate preservation" = "#FFA726",
                               "Weak preservation" = "#FFE082",
                               "No preservation" = "#D32F2F")) +
  labs(title = "Module Preservation: Control → High Sugar\n(Average Edge Weight)",
       x = "Module", 
       y = "Average Edge Weight Preservation",
       fill = "Preservation Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right")

print(p)
dev.off()
cat("  Saved: 03_preservation_levels_barplot.pdf\n")

# 7d. Heatmap of preservation statistics
pdf(file.path(output_dir, "04_preservation_heatmap.pdf"), width = 10, height = 8)

heatmap_data <- preservation_stats %>%
  select(module, avg_weight_preservation, avg_cor_preservation, 
         avg_contrib_preservation, coherence_preservation) %>%
  column_to_rownames("module") %>%
  as.matrix()

# Scale by column (z-score)
heatmap_data_scaled <- scale(heatmap_data)

colnames(heatmap_data_scaled) <- c("Avg Weight", "Avg Correlation", 
                                    "Avg Contribution", "Coherence")

pheatmap(heatmap_data_scaled,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 8,
         main = "Module Preservation Statistics (Z-scored)\nControl → High Sugar",
         angle_col = 45)

dev.off()
cat("  Saved: 04_preservation_heatmap.pdf\n")

# 7e. Module size vs preservation
pdf(file.path(output_dir, "05_size_vs_preservation.pdf"), width = 10, height = 8)

par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
plot(plot_data$moduleSize, plot_data$avg_weight_preservation,
     pch = 21, bg = plot_data$module,
     cex = 2,
     xlab = "Module Size (number of genes)",
     ylab = "Average Edge Weight Preservation",
     main = "Module Size vs Preservation",
     cex.lab = 1.2, cex.axis = 1.1)

text(plot_data$moduleSize, plot_data$avg_weight_preservation,
     labels = plot_data$module, pos = 3, cex = 0.7)

grid()

dev.off()
cat("  Saved: 05_size_vs_preservation.pdf\n")

# ----- 8. Summary statistics -----
cat("\n=== Module Preservation Summary ===\n\n")

preservation_counts <- preservation_stats %>%
  count(preservation_level)

cat("Preservation levels:\n")
print(preservation_counts)

cat("\nModules with strong preservation (p < 0.001):\n")
strong_pres <- preservation_stats %>% 
  filter(preservation_level == "Strong preservation") %>%
  select(module, moduleSize, avg_weight_preservation, avg_weight_pval)
if (nrow(strong_pres) > 0) {
  print(strong_pres)
} else {
  cat("  None\n")
}

cat("\nModules with no preservation:\n")
no_pres <- preservation_stats %>% 
  filter(preservation_level == "No preservation") %>%
  select(module, moduleSize, avg_weight_preservation, avg_weight_pval)
if (nrow(no_pres) > 0) {
  print(no_pres)
} else {
  cat("  None\n")
}

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", output_dir, "\n\n")

cat("Interpretation:\n")
cat("  NetRep tests module preservation using permutation-based null distributions\n")
cat("  - avg.weight: Average edge weight preservation (higher = better)\n")
cat("  - coherence: Within-module connectivity preservation (higher = better)\n")
cat("  - p < 0.001: Strong evidence of preservation\n")
cat("  - p < 0.05: Moderate evidence of preservation\n")
cat("  - p > 0.05: No significant preservation (module reorganized under High Sugar)\n\n")

cat("NetRep advantages:\n")
cat("  - 10-100x faster than WGCNA modulePreservation\n")
cat("  - More memory efficient\n")
cat("  - Suitable for large networks and many permutations\n\n")

cat("Files generated:\n")
cat("  - netrep_preservation_statistics.csv (main results table with p-values)\n")
cat("  - netrep_preservation_full_results.rds (complete NetRep output)\n")
cat("  - 01_preservation_scatter.pdf (preservation statistics scatter)\n")
cat("  - 02_preservation_pvalues.pdf (significance plot)\n")
cat("  - 03_preservation_levels_barplot.pdf (preservation levels by module)\n")
cat("  - 04_preservation_heatmap.pdf (all statistics heatmap)\n")
cat("  - 05_size_vs_preservation.pdf (module size relationship)\n")
