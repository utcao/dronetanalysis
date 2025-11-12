#!/usr/bin/env Rscript
# ==============================================================================
# Compare gene metrics between control and treatment conditions
#
# Script by Gabriel Thornes
#
# Last Updated: 12/11/2025
#
# This script:
#   1. Loads gene metrics + expression stats for two conditions
#   2. Generates exploratory visualizations (PCA, scatterplots, histograms, heatmaps)
#   3. Identifies genes with significant differences in network properties between conditions
#   4. Exports plots and summary tables
# ==============================================================================

rm(list = ls())

# ----- 1. Load required packages -----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
  library(corrplot)
  library(scales)
  library(argparse)
})

# ----- 2. Command-line arguments -----
parser <- ArgumentParser(description = 'Compare gene metrics between two conditions')
parser$add_argument('--ctrl-file', help = 'Control condition gene metrics + expression CSV', 
                   default = 'results/network_features/gene_metrics/adjacency_gene_metrics_with_expression.csv')
parser$add_argument('--treat-file', help = 'Treatment condition gene metrics + expression CSV', 
                   default = 'results/network_features/gene_metrics/HS_adjacency_gene_metrics_with_expression.csv')
parser$add_argument('--ctrl-label', help = 'Label for control condition', default = 'Control')
parser$add_argument('--treat-label', help = 'Label for treatment condition', default = 'High Sugar')
parser$add_argument('--output-dir', help = 'Directory to save plots and tables', 
                   default = 'results/analysis')
parser$add_argument('--focus-module', help = 'Optional: focus analysis on specific module (e.g., turquoise)', 
                   default = NULL)
args <- parser$parse_args()

output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Comparing Gene Metrics Between Conditions ===\n")
cat("Control file:", args$ctrl_file, "\n")
cat("Treatment file:", args$treat_file, "\n")
cat("Control label:", args$ctrl_label, "\n")
cat("Treatment label:", args$treat_label, "\n")
cat("Output directory:", output_dir, "\n\n")
if (!is.null(args$focus_module)) {
  cat("Focusing analysis on module:", args$focus_module, "\n\n")
}

# ----- 3. Load data -----
cat("Loading data...\n")
ctrl_metrics <- read.csv(args$ctrl_file, stringsAsFactors = FALSE)
treat_metrics <- read.csv(args$treat_file, stringsAsFactors = FALSE)

cat("Control:", nrow(ctrl_metrics), "genes,", ncol(ctrl_metrics), "columns\n")
cat("Treatment:", nrow(treat_metrics), "genes,", ncol(treat_metrics), "columns\n")

# Merge both conditions into one table for comparison
combined <- merge(
  ctrl_metrics %>% select(gene, weighted_connectivity, degree, betweenness_centrality, 
                          clustering_coefficient, eigenvector_centrality,
                          mean_expression, variance_expression, module) %>%
    rename_with(~paste0(., "_ctrl"), -gene),
  treat_metrics %>% select(gene, weighted_connectivity, degree, betweenness_centrality,
                           clustering_coefficient, eigenvector_centrality,
                           mean_expression, variance_expression, module) %>%
    rename_with(~paste0(., "_treat"), -gene),
  by = "gene", all = TRUE
)

cat("  Combined:", nrow(combined), "genes matched between conditions\n\n")

# ----- 3.5. Create module colour mapping -----
cat("Creating module colour mapping...\n")

# Get all unique modules from both conditions
all_modules <- unique(c(combined$module_ctrl, combined$module_treat))
all_modules <- all_modules[!is.na(all_modules)]

# Use WGCNA's labels2colors to get standard WGCNA colours for module names
# Convert module names (colours) to numeric labels, then back to ensure consistency
module_color_map <- setNames(
  WGCNA::labels2colors(as.numeric(as.factor(all_modules))),
  all_modules
)

# For named modules (like "turquoise", "blue", etc), use the name itself as the colour
# This ensures "turquoise" module shows in turquoise colour
for (mod in all_modules) {
  if (mod %in% colors()) {
    # If module name is a valid R color, use it
    module_color_map[mod] <- mod
  }
}

cat("  Module colours mapped:", length(module_color_map), "modules\n")
if (!is.null(args$focus_module)) {
  cat("Filtering to module:", args$focus_module, "\n")
  combined <- combined %>% 
    filter(module_ctrl == args$focus_module | module_treat == args$focus_module)
  cat("  Genes in module:", nrow(combined), "\n\n")
}

# ----- 5. PCA of network metrics -----
cat("Generating PCA plot...\n")

# Select numeric columns for PCA (network metrics only, not expression)
pca_vars <- c("weighted_connectivity", "degree", "betweenness_centrality", 
              "clustering_coefficient", "eigenvector_centrality")

# Extract PCA data and ensure consistent column names
pca_ctrl <- combined[, paste0(pca_vars, "_ctrl")] %>% na.omit()
colnames(pca_ctrl) <- pca_vars
pca_ctrl$condition <- args$ctrl_label

pca_treat <- combined[, paste0(pca_vars, "_treat")] %>% na.omit()
colnames(pca_treat) <- pca_vars
pca_treat$condition <- args$treat_label

# Combined PCA - now both have the same column structure
combined_pca_data <- rbind(pca_ctrl, pca_treat)

pca_result <- prcomp(combined_pca_data[, pca_vars], scale = TRUE, center = TRUE)
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  condition = combined_pca_data$condition
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, fill = condition)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_density2d(alpha = 0.2, show.legend = FALSE) +
  labs(title = "PCA of Network Metrics",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "01_pca_network_metrics.pdf"), p_pca, width = 8, height = 6)
cat("  Saved: 01_pca_network_metrics.pdf\n")

# ----- 6. Histogram comparisons: Network metrics -----
cat("Generating histogram comparisons...\n")

# Define custom colors for conditions (use setNames to create named vector dynamically)
colors_conditions <- setNames(c("#2E86AB", "#A23B72"), c(args$ctrl_label, args$treat_label))  # Blue and magenta

hist_plots <- lapply(pca_vars, function(var) {
  var_ctrl <- paste0(var, "_ctrl")
  var_treat <- paste0(var, "_treat")
  
  hist_data <- bind_rows(
    combined %>% select(all_of(var_ctrl)) %>% rename(value = !!var_ctrl) %>% mutate(condition = args$ctrl_label),
    combined %>% select(all_of(var_treat)) %>% rename(value = !!var_treat) %>% mutate(condition = args$treat_label)
  ) %>% filter(!is.na(value))
  
  ggplot(hist_data, aes(x = value, fill = condition)) +
    geom_histogram(bins = 50, position = "dodge", alpha = 0.75) +
    scale_fill_manual(values = colors_conditions) +
    facet_wrap(~condition, scales = "free_y") +
    labs(title = gsub("_", " ", var), x = "Value", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
})

p_hists <- do.call(gridExtra::grid.arrange, c(hist_plots, ncol = 2))
ggsave(file.path(output_dir, "02_histogram_network_metrics.pdf"), p_hists, width = 12, height = 10)
cat("  Saved: 02_histogram_network_metrics.pdf\n")

# ----- 7. Expression histograms -----
cat("Generating expression statistics comparisons...\n")

expr_plots <- list(
  ggplot(bind_rows(
    combined %>% select(mean_expression_ctrl) %>% rename(value = mean_expression_ctrl) %>% mutate(condition = args$ctrl_label, metric = "Mean Expression"),
    combined %>% select(mean_expression_treat) %>% rename(value = mean_expression_treat) %>% mutate(condition = args$treat_label, metric = "Mean Expression")
  ) %>% filter(!is.na(value)), aes(x = value, fill = condition)) +
    geom_histogram(bins = 50, position = "dodge", alpha = 0.75) +
    scale_fill_manual(values = colors_conditions) +
    facet_wrap(~condition, scales = "free_y") +
    labs(title = "Mean Expression", x = "Expression Level", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none"),
  
  ggplot(bind_rows(
    combined %>% select(variance_expression_ctrl) %>% rename(value = variance_expression_ctrl) %>% mutate(condition = args$ctrl_label, metric = "Variance"),
    combined %>% select(variance_expression_treat) %>% rename(value = variance_expression_treat) %>% mutate(condition = args$treat_label, metric = "Variance")
  ) %>% filter(!is.na(value)), aes(x = value, fill = condition)) +
    geom_histogram(bins = 50, position = "dodge", alpha = 0.75) +
    scale_fill_manual(values = colors_conditions) +
    facet_wrap(~condition, scales = "free_y") +
    labs(title = "Expression Variance", x = "Variance", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
)

p_expr_hists <- do.call(gridExtra::grid.arrange, c(expr_plots, ncol = 2))
ggsave(file.path(output_dir, "03_histogram_expression_stats.pdf"), p_expr_hists, width = 10, height = 5)
cat("  Saved: 03_histogram_expression_stats.pdf\n")

# ----- 7b. Overlaid density plots for expression statistics -----
cat("Generating overlaid density plots...\n")

expr_overlay_plots <- list(
  ggplot(bind_rows(
    combined %>% select(mean_expression_ctrl) %>% rename(value = mean_expression_ctrl) %>% mutate(condition = args$ctrl_label),
    combined %>% select(mean_expression_treat) %>% rename(value = mean_expression_treat) %>% mutate(condition = args$treat_label)
  ) %>% filter(!is.na(value)), aes(x = value, fill = condition, color = condition)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    scale_fill_manual(values = colors_conditions) +
    scale_color_manual(values = colors_conditions) +
    labs(title = "Mean Expression Distribution (Overlaid)", x = "Expression Level", y = "Density", fill = "Condition", color = "Condition") +
    theme_minimal() +
    theme(legend.position = "right"),
  
  ggplot(bind_rows(
    combined %>% select(variance_expression_ctrl) %>% rename(value = variance_expression_ctrl) %>% mutate(condition = args$ctrl_label),
    combined %>% select(variance_expression_treat) %>% rename(value = variance_expression_treat) %>% mutate(condition = args$treat_label)
  ) %>% filter(!is.na(value)), aes(x = value, fill = condition, color = condition)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    scale_fill_manual(values = colors_conditions) +
    scale_color_manual(values = colors_conditions) +
    labs(title = "Expression Variance Distribution (Overlaid)", x = "Variance", y = "Density", fill = "Condition", color = "Condition") +
    theme_minimal() +
    theme(legend.position = "right")
)

p_expr_overlay <- do.call(gridExtra::grid.arrange, c(expr_overlay_plots, ncol = 2))
ggsave(file.path(output_dir, "03b_overlaid_expression_distributions.pdf"), p_expr_overlay, width = 12, height = 5)
cat("  Saved: 03b_overlaid_expression_distributions.pdf\n")

# ----- 8. Scatterplots: Ctrl vs Treat network metrics -----
cat("Generating scatterplots (Ctrl vs Treat)...\n")

scatter_plots <- lapply(pca_vars, function(var) {
  var_ctrl <- paste0(var, "_ctrl")
  var_treat <- paste0(var, "_treat")
  
  scatter_data <- combined %>% 
    select(all_of(c(var_ctrl, var_treat, "module_ctrl"))) %>%
    filter(!is.na(!!sym(var_ctrl)) & !is.na(!!sym(var_treat)))
  
  ggplot(scatter_data, aes_string(x = var_ctrl, y = var_treat, color = "module_ctrl")) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", linewidth = 1) +
    scale_color_manual(values = module_color_map, name = "Module (Ctrl)") +
    labs(title = gsub("_", " ", var),
         x = paste0(args$ctrl_label, " - ", gsub("_", " ", var)),
         y = paste0(args$treat_label, " - ", gsub("_", " ", var))) +
    theme_minimal() +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 10))
})

p_scatters <- do.call(gridExtra::grid.arrange, c(scatter_plots, ncol = 2))
ggsave(file.path(output_dir, "04_scatterplots_ctrl_vs_treat.pdf"), p_scatters, width = 14, height = 12)
cat("  Saved: 04_scatterplots_ctrl_vs_treat.pdf\n")

# ----- 9. Expression vs Network metrics scatterplots -----
cat("Generating expression vs network metrics plots...\n")

expr_net_plots <- list(
  ggplot(combined %>% filter(!is.na(weighted_connectivity_ctrl) & !is.na(mean_expression_ctrl)), 
         aes(x = mean_expression_ctrl, y = weighted_connectivity_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Control: Expression vs Weighted Connectivity",
         x = "Mean Expression", y = "Weighted Connectivity") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(weighted_connectivity_treat) & !is.na(mean_expression_treat)), 
         aes(x = mean_expression_treat, y = weighted_connectivity_treat, color = module_treat)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Treatment: Expression vs Weighted Connectivity",
         x = "Mean Expression", y = "Weighted Connectivity") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(variance_expression_ctrl) & !is.na(degree_ctrl)), 
         aes(x = variance_expression_ctrl, y = degree_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Control: Expression Variance vs Degree",
         x = "Expression Variance", y = "Degree") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(variance_expression_treat) & !is.na(degree_treat)), 
         aes(x = variance_expression_treat, y = degree_treat, color = module_treat)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Treatment: Expression Variance vs Degree",
         x = "Expression Variance", y = "Degree") +
    theme_minimal(),

  ggplot(combined %>% filter(!is.na(mean_expression_ctrl) & !is.na(eigenvector_centrality_ctrl)), 
         aes(x = mean_expression_ctrl, y = eigenvector_centrality_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Control: Expression Variance vs Eigenvector Centrality",
         x = "Mean Expression", y = "Eigenvector Centrality") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(mean_expression_treat) & !is.na(eigenvector_centrality_treat)), 
         aes(x = mean_expression_treat, y = eigenvector_centrality_treat, color = module_treat)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Treatment: Expression Variance vs Eigenvector Centrality",
         x = "Mean Expression", y = "Eigenvector Centrality") +
    theme_minimal(),

    ggplot(combined %>% filter(!is.na(variance_expression_ctrl) & !is.na(eigenvector_centrality_ctrl)), 
         aes(x = variance_expression_ctrl, y = eigenvector_centrality_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Control: Expression Variance vs Eigenvector Centrality",
         x = "Expression Variance", y = "Eigenvector Centrality") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(variance_expression_treat) & !is.na(eigenvector_centrality_treat)), 
         aes(x = variance_expression_treat, y = eigenvector_centrality_treat, color = module_treat)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Treatment: Expression Variance vs Eigenvector Centrality",
         x = "Expression Variance", y = "Eigenvector Centrality") +
    theme_minimal()

)

# Save expression vs network metrics plots in multiple pages (6 plots total, 2 per page)
pdf(file.path(output_dir, "05_expression_vs_network_metrics.pdf"), width = 14, height = 6)
for (i in seq(1, length(expr_net_plots), by = 2)) {
  plots_to_arrange <- expr_net_plots[i:min(i+1, length(expr_net_plots))]
  p_expr_net_page <- do.call(gridExtra::grid.arrange, c(plots_to_arrange, ncol = 2))
  print(p_expr_net_page)
}
dev.off()
cat("  Saved: 05_expression_vs_network_metrics.pdf\n")

# ----- 10. Identify differentially important genes -----
cat("Identifying genes with condition-specific patterns...\n")

# Compute differences
diff_metrics <- combined %>%
  mutate(
    delta_connectivity = weighted_connectivity_treat - weighted_connectivity_ctrl,
    delta_degree = degree_treat - degree_ctrl,
    delta_betweenness = betweenness_centrality_treat - betweenness_centrality_ctrl,
    delta_expression = mean_expression_treat - mean_expression_ctrl,
    delta_variance = variance_expression_treat - variance_expression_ctrl,
    # Since expression values are log2-transformed by VOOM, fold change is 2^(delta_expression)
    fold_change_expression = 2^(delta_expression)
  ) %>%
  select(gene, module_ctrl, module_treat, starts_with("delta_"), fold_change_expression)

# Top genes with increased connectivity in treatment
top_increased_conn <- diff_metrics %>% 
  filter(!is.na(delta_connectivity)) %>%
  arrange(desc(delta_connectivity)) %>%
  slice_head(n = 20)

# Top genes with decreased connectivity in treatment
top_decreased_conn <- diff_metrics %>%
  filter(!is.na(delta_connectivity)) %>%
  arrange(delta_connectivity) %>%
  slice_head(n = 20)

# Top genes with fold change in expression
top_expr_ratio_up <- diff_metrics %>%
  filter(!is.na(fold_change_expression)) %>%
  arrange(desc(fold_change_expression)) %>%
  slice_head(n = 20)

# Top genes with expression changes
top_expr_up <- diff_metrics %>%
  filter(!is.na(delta_expression)) %>%
  arrange(desc(delta_expression)) %>%
  slice_head(n = 20)

# Top genes with expression variance changes
top_expr_var_up <- diff_metrics %>%
  filter(!is.na(delta_variance)) %>%
  arrange(desc(delta_variance)) %>%
  slice_head(n=20)

# Write summary tables
write.csv(top_increased_conn, file.path(output_dir, "top_increased_connectivity.csv"), row.names = FALSE)
write.csv(top_decreased_conn, file.path(output_dir, "top_decreased_connectivity.csv"), row.names = FALSE)
write.csv(top_expr_ratio_up, file.path(output_dir, "top_expression_fold_changes.csv"), row.names = FALSE)
write.csv(top_expr_var_up, file.path(output_dir, "top_variance_expression_changes.csv"), row.names = FALSE)
write.csv(top_expr_up, file.path(output_dir, "top_expression_changes.csv"), row.names = FALSE)

cat("  Saved summary tables:\n")
cat("    - top_increased_connectivity.csv\n")
cat("    - top_decreased_connectivity.csv\n")
cat("    - top_expression_changes.csv\n\n")

# ----- 11. Volcano-like plots for network metrics changes -----
cat("Generating condition-change plots...\n")

volcano_plots <- list(
  ggplot(diff_metrics %>% filter(!is.na(delta_connectivity) & !is.na(delta_expression)), 
         aes(x = delta_connectivity, y = delta_expression, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Ctrl Module") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "Connectivity vs Expression Changes (Ctrl -> Treat)",
         x = "Change in Weighted Connectivity",
         y = "Change in Mean Expression") +
    theme_minimal(),
  
  ggplot(diff_metrics %>% filter(!is.na(delta_degree) & !is.na(delta_betweenness)), 
         aes(x = delta_degree, y = delta_betweenness, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Ctrl Module") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "Degree vs Betweenness Changes (Ctrl -> Treat)",
         x = "Change in Degree",
         y = "Change in Betweenness Centrality") +
    theme_minimal(),

  ggplot(diff_metrics %>% filter(!is.na(delta_expression) & !is.na(delta_variance)), 
         aes(x = delta_expression, y = delta_variance, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Ctrl Module") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "Expression vs Variance Changes (Ctrl -> Treat)",
         x = "Change in Mean Expression",
         y = "Change in Expression Variance") +
    theme_minimal()

)

p_volcano <- do.call(gridExtra::grid.arrange, c(volcano_plots, ncol = 3))
ggsave(file.path(output_dir, "06_condition_changes.pdf"), p_volcano, width = 12, height = 20)
cat("  Saved: 06_condition_changes.pdf\n")

# ----- 12. Summary statistics -----
cat("\n=== Summary Statistics ===\n")

cat("\nNetwork Metrics Comparison:\n")
for (var in pca_vars) {
  var_ctrl <- paste0(var, "_ctrl")
  var_treat <- paste0(var, "_treat")
  
  ctrl_vals <- combined[[var_ctrl]]
  treat_vals <- combined[[var_treat]]
  
  cat(sprintf("%-30s: Ctrl mean=%.4f, Treat mean=%.4f, Δ=%.4f\n",
              var,
              mean(ctrl_vals, na.rm = TRUE),
              mean(treat_vals, na.rm = TRUE),
              mean(treat_vals - ctrl_vals, na.rm = TRUE)))
}

cat("\nExpression Statistics Comparison:\n")
expr_vars <- c("mean_expression", "variance_expression")
for (var in expr_vars) {
  var_ctrl <- paste0(var, "_ctrl")
  var_treat <- paste0(var, "_treat")
  
  ctrl_vals <- combined[[var_ctrl]]
  treat_vals <- combined[[var_treat]]
  
  cat(sprintf("%-30s: Ctrl mean=%.4f, Treat mean=%.4f, Δ=%.4f\n",
              var,
              mean(ctrl_vals, na.rm = TRUE),
              mean(treat_vals, na.rm = TRUE),
              mean(treat_vals - ctrl_vals, na.rm = TRUE)))
}

# Save combined differences table
write.csv(diff_metrics, file.path(output_dir, "gene_metric_differences.csv"), row.names = FALSE)
cat("\nSaved gene metric differences table: gene_metric_differences.csv\n")

# Save summary statistics to a CSV
summary_stats <- data.frame(
  Metric = c(pca_vars, expr_vars),
  Ctrl_Mean = sapply(c(pca_vars, expr_vars), function(var) mean(combined[[paste0(var, "_ctrl")]], na.rm = TRUE)),
  Treat_Mean = sapply(c(pca_vars, expr_vars), function(var) mean(combined[[paste0(var, "_treat")]], na.rm = TRUE)),
  Delta_Mean = sapply(c(pca_vars, expr_vars), function(var) mean(combined[[paste0(var, "_treat")]] - combined[[paste0(var, "_ctrl")]], na.rm = TRUE))
)

write.csv(summary_stats, file.path(output_dir, "Ctrl_HS_summary_statistics.csv"), row.names = FALSE)
cat("Saved summary statistics: Ctrl_HS_summary_statistics.csv\n")

cat("\n=== Analysis Complete ===\n")
cat("All plots saved to:", output_dir, "\n")
