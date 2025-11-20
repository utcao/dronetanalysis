#!/usr/bin/env Rscript
# ==============================================================================
# Module Overlap Analysis between Control and High Sugar Conditions
#
# Script by Gabriel Thornes
#
# Last Updated: 20/11/2025
#
# This script:
#   1. Loads module assignments for both conditions
#   2. Calculates overlap between modules using Fisher's exact test
#   3. Computes Jaccard similarity and overlap coefficients
#   4. Generates heatmaps and summary statistics
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(argparse)
})

# ----- 2. Command-line arguments -----
parser <- ArgumentParser(description = 'Module overlap analysis between conditions')
parser$add_argument('--ctrl-modules', help = 'Control module assignments CSV', 
                   default = 'results/network_features/features_calc/adjacency/modules/gene_connectivity.csv')
parser$add_argument('--treat-modules', help = 'Treatment module assignments CSV', 
                   default = 'results/network_features/features_calc/HS_adjacency/modules/gene_connectivity.csv')
parser$add_argument('--output-dir', help = 'Directory to save results', 
                   default = 'results/analysis/module_overlap')
parser$add_argument('--p-threshold', type = 'double', help = 'P-value threshold for significance', 
                   default = 0.05)
args <- parser$parse_args()

output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Module Overlap Analysis ===\n")
cat("Control modules:", args$ctrl_modules, "\n")
cat("Treatment modules:", args$treat_modules, "\n")
cat("Output directory:", output_dir, "\n")
cat("P-value threshold:", args$p_threshold, "\n\n")

# ----- 3. Load data -----
cat("Loading module assignments...\n")

ctrl_modules <- read.csv(args$ctrl_modules, stringsAsFactors = FALSE)
treat_modules <- read.csv(args$treat_modules, stringsAsFactors = FALSE)

cat("  Control:", nrow(ctrl_modules), "genes in", length(unique(ctrl_modules$module)), "modules\n")
cat("  Treatment:", nrow(treat_modules), "genes in", length(unique(treat_modules$module)), "modules\n")

# Find common genes
common_genes <- intersect(ctrl_modules$gene, treat_modules$gene)
cat("  Common genes:", length(common_genes), "\n\n")

# Filter to common genes only
ctrl_data <- ctrl_modules %>% 
  filter(gene %in% common_genes) %>%
  select(gene, module) %>%
  rename(module_ctrl = module)

treat_data <- treat_modules %>% 
  filter(gene %in% common_genes) %>%
  select(gene, module) %>%
  rename(module_treat = module)

# Merge
combined <- merge(ctrl_data, treat_data, by = "gene")

cat("Module assignments loaded for", nrow(combined), "common genes\n\n")

# Get unique modules
ctrl_mods <- sort(unique(combined$module_ctrl))
treat_mods <- sort(unique(combined$module_treat))

cat("Control modules:", length(ctrl_mods), "-", paste(ctrl_mods, collapse = ", "), "\n")
cat("Treatment modules:", length(treat_mods), "-", paste(treat_mods, collapse = ", "), "\n\n")

# ----- 4. Calculate overlap statistics -----
cat("Calculating module overlap statistics...\n")

# Create contingency table for each pair
overlap_results <- list()

for (ctrl_mod in ctrl_mods) {
  for (treat_mod in treat_mods) {
    # Genes in each module
    ctrl_genes <- combined$gene[combined$module_ctrl == ctrl_mod]
    treat_genes <- combined$gene[combined$module_treat == treat_mod]
    
    # Overlap
    overlap_genes <- intersect(ctrl_genes, treat_genes)
    overlap_count <- length(overlap_genes)
    
    # Contingency table
    # |              | In treat_mod | Not in treat_mod |
    # | In ctrl_mod  |      a       |        b         |
    # | Not ctrl_mod |      c       |        d         |
    
    a <- overlap_count
    b <- length(ctrl_genes) - a
    c <- length(treat_genes) - a
    d <- nrow(combined) - a - b - c
    
    # Fisher's exact test
    contingency <- matrix(c(a, b, c, d), nrow = 2)
    fisher_test <- fisher.test(contingency)
    
    # Jaccard similarity: intersection / union
    jaccard <- overlap_count / length(union(ctrl_genes, treat_genes))
    
    # Overlap coefficient: intersection / min(|A|, |B|)
    overlap_coef <- overlap_count / min(length(ctrl_genes), length(treat_genes))
    
    # Proportion of ctrl module that overlaps
    prop_ctrl_overlap <- overlap_count / length(ctrl_genes)
    
    # Proportion of treat module that overlaps
    prop_treat_overlap <- overlap_count / length(treat_genes)
    
    overlap_results[[paste0(ctrl_mod, "_", treat_mod)]] <- data.frame(
      ctrl_module = ctrl_mod,
      treat_module = treat_mod,
      ctrl_size = length(ctrl_genes),
      treat_size = length(treat_genes),
      overlap_count = overlap_count,
      jaccard_similarity = jaccard,
      overlap_coefficient = overlap_coef,
      prop_ctrl_in_treat = prop_ctrl_overlap,
      prop_treat_in_ctrl = prop_treat_overlap,
      fisher_pval = fisher_test$p.value,
      fisher_odds_ratio = fisher_test$estimate,
      significant = fisher_test$p.value < args$p_threshold
    )
  }
}

overlap_df <- bind_rows(overlap_results)

# Adjust p-values for multiple testing (Bonferroni)
overlap_df$fisher_pval_adj <- p.adjust(overlap_df$fisher_pval, method = "bonferroni")
overlap_df$significant_adj <- overlap_df$fisher_pval_adj < args$p_threshold

cat("  Overlap analysis complete for", nrow(overlap_df), "module pairs\n\n")

# ----- 5. Save results -----
cat("Saving overlap statistics...\n")

# Sort by significance and overlap
overlap_df_sorted <- overlap_df %>%
  arrange(fisher_pval, desc(overlap_count))

write.csv(overlap_df_sorted, 
          file.path(output_dir, "module_overlap_statistics.csv"), 
          row.names = FALSE)
cat("  Saved: module_overlap_statistics.csv\n")

# Save significant overlaps only
sig_overlaps <- overlap_df_sorted %>%
  filter(significant_adj == TRUE, overlap_count > 0)

if (nrow(sig_overlaps) > 0) {
  write.csv(sig_overlaps, 
            file.path(output_dir, "significant_module_overlaps.csv"), 
            row.names = FALSE)
  cat("  Saved: significant_module_overlaps.csv (", nrow(sig_overlaps), "significant overlaps)\n")
} else {
  cat("  No significant overlaps found after correction\n")
}

# ----- 6. Generate heatmaps -----
cat("\nGenerating heatmaps...\n")

# 6a. Overlap count heatmap
overlap_matrix <- overlap_df %>%
  select(ctrl_module, treat_module, overlap_count) %>%
  pivot_wider(names_from = treat_module, values_from = overlap_count, values_fill = 0) %>%
  column_to_rownames("ctrl_module") %>%
  as.matrix()

pdf(file.path(output_dir, "01_overlap_count_heatmap.pdf"), width = 10, height = 10)
pheatmap(overlap_matrix,
         color = colorRampPalette(c("white", "#FFF7BC", "#FEC44F", "#D95F0E", "#993404"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.0f",
         number_color = "black",
         fontsize_number = 8,
         main = "Module Overlap: Gene Count\n(Control modules × Treatment modules)",
         angle_col = 45)
dev.off()
cat("  Saved: 01_overlap_count_heatmap.pdf\n")

# 6b. Jaccard similarity heatmap
jaccard_matrix <- overlap_df %>%
  select(ctrl_module, treat_module, jaccard_similarity) %>%
  pivot_wider(names_from = treat_module, values_from = jaccard_similarity, values_fill = 0) %>%
  column_to_rownames("ctrl_module") %>%
  as.matrix()

pdf(file.path(output_dir, "02_jaccard_similarity_heatmap.pdf"), width = 10, height = 10)
pheatmap(jaccard_matrix,
         color = colorRampPalette(c("white", "#C6DBEF", "#6BAED6", "#2171B5", "#08519C"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 8,
         main = "Module Overlap: Jaccard Similarity\n(intersection / union)",
         angle_col = 45,
         breaks = seq(0, 1, length.out = 101))
dev.off()
cat("  Saved: 02_jaccard_similarity_heatmap.pdf\n")

# 6c. Overlap coefficient heatmap
overlap_coef_matrix <- overlap_df %>%
  select(ctrl_module, treat_module, overlap_coefficient) %>%
  pivot_wider(names_from = treat_module, values_from = overlap_coefficient, values_fill = 0) %>%
  column_to_rownames("ctrl_module") %>%
  as.matrix()

pdf(file.path(output_dir, "03_overlap_coefficient_heatmap.pdf"), width = 10, height = 10)
pheatmap(overlap_coef_matrix,
         color = colorRampPalette(c("white", "#E5F5E0", "#A1D99B", "#31A354", "#006D2C"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 8,
         main = "Module Overlap: Overlap Coefficient\n(intersection / min(size))",
         angle_col = 45,
         breaks = seq(0, 1, length.out = 101))
dev.off()
cat("  Saved: 03_overlap_coefficient_heatmap.pdf\n")

# 6d. -log10(p-value) heatmap for significance
pval_matrix <- overlap_df %>%
  mutate(log10p = pmin(-log10(fisher_pval_adj + 1e-300), 300)) %>%  # Cap at 300, avoid -log10(0)
  select(ctrl_module, treat_module, log10p) %>%
  pivot_wider(names_from = treat_module, values_from = log10p, values_fill = 0) %>%
  column_to_rownames("ctrl_module") %>%
  as.matrix()

pdf(file.path(output_dir, "04_significance_heatmap.pdf"), width = 10, height = 10)
pheatmap(pval_matrix,
         color = colorRampPalette(c("white", "#FCBBA1", "#FC9272", "#FB6A4A", "#CB181D", "#67000D"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.1f",
         number_color = "black",
         fontsize_number = 8,
         main = "Module Overlap: Statistical Significance\n(-log10 adjusted p-value, Fisher's exact test)",
         angle_col = 45)
dev.off()
cat("  Saved: 04_significance_heatmap.pdf\n")

# ----- 7. Module preservation summary -----
cat("\nGenerating module preservation summary...\n")

# For each control module, find best match in treatment
best_matches <- overlap_df %>%
  group_by(ctrl_module) %>%
  arrange(desc(jaccard_similarity)) %>%
  slice(1) %>%
  ungroup() %>%
  select(ctrl_module, treat_module, ctrl_size, overlap_count, 
         jaccard_similarity, overlap_coefficient, fisher_pval_adj) %>%
  mutate(preservation_score = overlap_coefficient) %>%
  arrange(desc(preservation_score))

write.csv(best_matches, 
          file.path(output_dir, "best_module_matches.csv"), 
          row.names = FALSE)
cat("  Saved: best_module_matches.csv\n")

# Preservation categories
preservation_summary <- best_matches %>%
  mutate(
    preservation_level = case_when(
      overlap_coefficient >= 0.7 ~ "Strong preservation",
      overlap_coefficient >= 0.4 ~ "Moderate preservation",
      overlap_coefficient >= 0.2 ~ "Weak preservation",
      TRUE ~ "No preservation"
    )
  )

# Bar plot
pdf(file.path(output_dir, "05_module_preservation_barplot.pdf"), width = 12, height = 6)

p <- ggplot(preservation_summary, 
            aes(x = reorder(ctrl_module, overlap_coefficient), 
                y = overlap_coefficient, 
                fill = preservation_level)) +
  geom_col() +
  geom_hline(yintercept = c(0.2, 0.4, 0.7), linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Strong preservation" = "#2E7D32",
                               "Moderate preservation" = "#FFA726",
                               "Weak preservation" = "#FFE082",
                               "No preservation" = "#D32F2F")) +
  labs(title = "Module Preservation: Control → High Sugar\n(Best match overlap coefficient)",
       x = "Control Module", 
       y = "Overlap Coefficient (best match in Treatment)",
       fill = "Preservation Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)
dev.off()
cat("  Saved: 05_module_preservation_barplot.pdf\n")

# ----- 8. Summary statistics -----
cat("\n=== Module Overlap Summary ===\n\n")

cat("Preservation levels (by best match overlap coefficient):\n")
pres_counts <- table(preservation_summary$preservation_level)
print(pres_counts)

cat("\nModules with strong preservation (overlap coef >= 0.7):\n")
strong_pres <- preservation_summary %>% 
  filter(preservation_level == "Strong preservation") %>%
  select(ctrl_module, treat_module, overlap_coefficient, jaccard_similarity)
if (nrow(strong_pres) > 0) {
  print(strong_pres)
} else {
  cat("  None\n")
}

cat("\nModules with no/weak preservation (overlap coef < 0.4):\n")
no_pres <- preservation_summary %>% 
  filter(overlap_coefficient < 0.4) %>%
  select(ctrl_module, treat_module, overlap_coefficient, jaccard_similarity)
if (nrow(no_pres) > 0) {
  print(no_pres)
} else {
  cat("  None\n")
}

cat("\nStatistically significant overlaps (adj p <", args$p_threshold, "):\n")
cat("  Total significant pairs:", sum(overlap_df$significant_adj), "\n")

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", output_dir, "\n\n")

cat("Interpretation:\n")
cat("  - Jaccard similarity: Overlap relative to union (0-1, higher = more similar)\n")
cat("  - Overlap coefficient: Overlap relative to smaller module (0-1, measures containment)\n")
cat("  - Fisher's exact test: Statistical significance of overlap (corrected for multiple testing)\n")
cat("  - Overlap coef >= 0.7: Strong preservation (most genes maintain module membership)\n")
cat("  - Overlap coef 0.4-0.7: Moderate preservation (partial rewiring)\n")
cat("  - Overlap coef < 0.4: Weak/no preservation (extensive reorganization)\n\n")

cat("Files generated:\n")
cat("  - module_overlap_statistics.csv (all pairwise comparisons)\n")
cat("  - significant_module_overlaps.csv (statistically significant overlaps only)\n")
cat("  - best_module_matches.csv (best treatment module match for each control module)\n")
cat("  - 01_overlap_count_heatmap.pdf\n")
cat("  - 02_jaccard_similarity_heatmap.pdf\n")
cat("  - 03_overlap_coefficient_heatmap.pdf\n")
cat("  - 04_significance_heatmap.pdf\n")
cat("  - 05_module_preservation_barplot.pdf\n")
