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
  library(tibble)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
  library(corrplot)
  library(scales)
  library(argparse)
})

# Helper function to format p-values for display
format_pval <- function(p) {
  if (p < 0.001) return("p < 0.001***")
  if (p < 0.01) return(sprintf("p = %.3f**", p))
  if (p < 0.05) return(sprintf("p = %.3f*", p))
  return(sprintf("p = %.3f", p))
}

# Helper function to get significance stars
get_sig_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

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

cat("=== Comparing Gene Metrics Between Conditions ===")
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
                          mean_expression, variance_expression, mad_expression, module) %>%
    rename_with(~paste0(., "_ctrl"), -gene),
  treat_metrics %>% select(gene, weighted_connectivity, degree, betweenness_centrality,
                           clustering_coefficient, eigenvector_centrality,
                           mean_expression, variance_expression, mad_expression, module) %>%
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

# ----- 5a2. PCA with network metrics + expression -----
cat("Generating combined PCA plot (network + expression)...\n")

# Select all variables including expression
pca_all_vars <- c("weighted_connectivity", "degree", "betweenness_centrality", 
                  "clustering_coefficient", "eigenvector_centrality",
                  "mean_expression", "variance_expression")

# Extract PCA data for combined analysis
pca_all_ctrl <- combined[, paste0(pca_all_vars, "_ctrl")] %>% na.omit()
colnames(pca_all_ctrl) <- pca_all_vars
pca_all_ctrl$condition <- args$ctrl_label

pca_all_treat <- combined[, paste0(pca_all_vars, "_treat")] %>% na.omit()
colnames(pca_all_treat) <- pca_all_vars
pca_all_treat$condition <- args$treat_label

# Combined PCA
combined_pca_all_data <- rbind(pca_all_ctrl, pca_all_treat)

pca_all_result <- prcomp(combined_pca_all_data[, pca_all_vars], scale = TRUE, center = TRUE)
pca_all_df <- data.frame(
  PC1 = pca_all_result$x[, 1],
  PC2 = pca_all_result$x[, 2],
  condition = combined_pca_all_data$condition
)

p_pca_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = condition, fill = condition)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_density2d(alpha = 0.2, show.legend = FALSE) +
  labs(title = "PCA of Network Metrics + Expression",
       x = paste0("PC1 (", round(summary(pca_all_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_all_result)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "01a2_pca_network_and_expression.pdf"), p_pca_all, width = 8, height = 6)
cat("  Saved: 01a2_pca_network_and_expression.pdf\n")

# ----- 5a. Delta barplot showing mean changes across all metrics -----
cat("Generating delta barplot for metric changes...\n")

# Calculate mean deltas for all metrics
delta_summary <- data.frame(
  Metric = c(pca_vars, "mean_expression", "variance_expression", "mad_expression"),
  Delta_Mean = c(
    sapply(pca_vars, function(var) {
      mean(combined[[paste0(var, "_treat")]] - combined[[paste0(var, "_ctrl")]], na.rm = TRUE)
    }),
    mean(combined$mean_expression_treat - combined$mean_expression_ctrl, na.rm = TRUE),
    mean(combined$variance_expression_treat - combined$variance_expression_ctrl, na.rm = TRUE),
    mean(combined$mad_expression_treat - combined$mad_expression_ctrl, na.rm = TRUE)
  )
) %>%
  mutate(
    Metric_Label = gsub("_", " ", Metric),
    Metric_Label = tools::toTitleCase(Metric_Label),
    Direction = ifelse(Delta_Mean > 0, "Increased", "Decreased"),
    Metric_Type = ifelse(Metric %in% pca_vars, "Network", "Expression")
  )

# Create barplot with colors indicating direction
p_delta <- ggplot(delta_summary, aes(x = reorder(Metric_Label, Delta_Mean), y = Delta_Mean, fill = Direction)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Increased" = "#A23B72", "Decreased" = "#2E86AB")) +
  coord_flip() +
  facet_wrap(~Metric_Type, scales = "free_y", ncol = 1) +
  labs(title = paste0("Mean Change in Metrics: ", args$ctrl_label, " → ", args$treat_label),
       x = "", y = "Mean Delta (Treatment - Control)",
       fill = "Direction") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(output_dir, "01a_metric_deltas_barplot.pdf"), p_delta, width = 10, height = 8)
cat("  Saved: 01a_metric_deltas_barplot.pdf\n")

# ----- 5b. Violin/Box plots for expression metrics comparison -----
cat("Generating expression metrics distribution comparison plots...\n")

# Define custom colors for conditions (needed for violin plots)
colors_conditions <- setNames(c("#2E86AB", "#A23B72"), c(args$ctrl_label, args$treat_label))  # Blue and magenta

# Prepare expression data in long format for violin plots
expr_metrics_long <- bind_rows(
  combined %>% 
    select(gene, mean_expression_ctrl, variance_expression_ctrl, mad_expression_ctrl) %>%
    rename(mean_expression = mean_expression_ctrl, variance_expression = variance_expression_ctrl, mad_expression = mad_expression_ctrl) %>%
    pivot_longer(cols = c(mean_expression, variance_expression, mad_expression), names_to = "metric", values_to = "value") %>%
    mutate(condition = args$ctrl_label),
  combined %>%
    select(gene, mean_expression_treat, variance_expression_treat, mad_expression_treat) %>%
    rename(mean_expression = mean_expression_treat, variance_expression = variance_expression_treat, mad_expression = mad_expression_treat) %>%
    pivot_longer(cols = c(mean_expression, variance_expression, mad_expression), names_to = "metric", values_to = "value") %>%
    mutate(condition = args$treat_label)
) %>% filter(!is.na(value))

# Statistical tests for expression metrics
mean_expr_test <- wilcox.test(
  combined$mean_expression_ctrl, 
  combined$mean_expression_treat, 
  paired = FALSE
)
var_expr_test <- wilcox.test(
  combined$variance_expression_ctrl, 
  combined$variance_expression_treat, 
  paired = FALSE
)
mad_expr_test <- wilcox.test(
  combined$mad_expression_ctrl, 
  combined$mad_expression_treat, 
  paired = FALSE
)

cat("  Mean expression test:", format_pval(mean_expr_test$p.value), "\n")
cat("  Variance expression test:", format_pval(var_expr_test$p.value), "\n")
cat("  MAD expression test:", format_pval(mad_expr_test$p.value), "\n")

# Create violin plots with overlaid boxplots and significance indicators
expr_violin_plots <- list(
  # Mean expression
  ggplot(expr_metrics_long %>% filter(metric == "mean_expression"), 
         aes(x = condition, y = value, fill = condition)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_manual(values = colors_conditions) +
    labs(title = "Mean Expression Distribution",
         subtitle = format_pval(mean_expr_test$p.value),
         x = "", y = "Mean Expression (log2 CPM)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11)),
  
  # Expression variance
  ggplot(expr_metrics_long %>% filter(metric == "variance_expression"), 
         aes(x = condition, y = value, fill = condition)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_manual(values = colors_conditions) +
    scale_y_log10(labels = scales::label_scientific()) +
    labs(title = "Expression Variance Distribution",
         subtitle = format_pval(var_expr_test$p.value),
         x = "", y = "Expression Variance (log10 scale)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11)),

   # Expression MAD
  ggplot(expr_metrics_long %>% filter(metric == "mad_expression"), 
         aes(x = condition, y = value, fill = condition)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_manual(values = colors_conditions) +
    labs(title = "Expression MAD Distribution",
         subtitle = format_pval(mad_expr_test$p.value),
         x = "", y = "MAD (Median Absolute Deviation)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11)),

  # Expression log MAD
  ggplot(expr_metrics_long %>% filter(metric == "mad_expression"), 
         aes(x = condition, y = value, fill = condition)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_manual(values = colors_conditions) +
    scale_y_log10(labels = scales::label_scientific()) +
    labs(title = "Expression logMAD Distribution",
         subtitle = format_pval(mad_expr_test$p.value),
         x = "", y = "MAD (Median Absolute Deviation, log10 scale)") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11))
)

# Save the plot directly without capturing (grid.arrange plots directly to device)
pdf(file.path(output_dir, "01b_expression_metrics_distributions.pdf"), width = 18, height = 12)
gridExtra::grid.arrange(grobs = expr_violin_plots, ncol = 2, 
                        top = grid::textGrob("Expression Metrics Distribution Comparison", 
                                            gp = grid::gpar(fontsize = 16, fontface = "bold")))
dev.off()
cat("  Saved: 01b_expression_metrics_distributions.pdf\n")

# ----- 6. Histogram comparisons: Network metrics -----
cat("Generating histogram comparisons...\n")

# Define custom colors for conditions (use setNames to create named vector dynamically)
colors_conditions <- setNames(c("#2E86AB", "#A23B72"), c(args$ctrl_label, args$treat_label))  # Blue and magenta

# Perform statistical tests for network metrics
network_tests <- lapply(pca_vars, function(var) {
  var_ctrl <- paste0(var, "_ctrl")
  var_treat <- paste0(var, "_treat")
  
  ctrl_vals <- combined[[var_ctrl]][!is.na(combined[[var_ctrl]])]
  treat_vals <- combined[[var_treat]][!is.na(combined[[var_treat]])]
  
  test_result <- wilcox.test(ctrl_vals, treat_vals, paired = FALSE)
  list(var = var, p.value = test_result$p.value)
})

# Bonferroni correction for multiple testing
n_tests <- length(network_tests)
for (i in seq_along(network_tests)) {
  network_tests[[i]]$p.adjusted <- network_tests[[i]]$p.value * n_tests
  network_tests[[i]]$p.adjusted <- min(network_tests[[i]]$p.adjusted, 1)  # Cap at 1
  cat("  ", network_tests[[i]]$var, ":", format_pval(network_tests[[i]]$p.adjusted), "(Bonferroni-corrected)\n")
}

hist_plots <- lapply(seq_along(pca_vars), function(i) {
  var <- pca_vars[i]
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
    labs(title = gsub("_", " ", var),
         subtitle = format_pval(network_tests[[i]]$p.adjusted),
         x = "Value", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.subtitle = element_text(hjust = 0.5, size = 9))
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
    labs(title = "Mean Expression Distribution (Overlaid)",
         subtitle = format_pval(mean_expr_test$p.value),
         x = "Expression Level", y = "Density", fill = "Condition", color = "Condition") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.subtitle = element_text(hjust = 0.5, size = 10)),
  
  ggplot(bind_rows(
    combined %>% select(variance_expression_ctrl) %>% rename(value = variance_expression_ctrl) %>% mutate(condition = args$ctrl_label),
    combined %>% select(variance_expression_treat) %>% rename(value = variance_expression_treat) %>% mutate(condition = args$treat_label)
  ) %>% filter(!is.na(value)), aes(x = value, fill = condition, color = condition)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    scale_fill_manual(values = colors_conditions) +
    scale_color_manual(values = colors_conditions) +
    labs(title = "Expression Variance Distribution (Overlaid)",
         subtitle = format_pval(var_expr_test$p.value),
         x = "Variance", y = "Density", fill = "Condition", color = "Condition") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.subtitle = element_text(hjust = 0.5, size = 10))
)

p_expr_overlay <- do.call(gridExtra::grid.arrange, c(expr_overlay_plots, ncol = 2))
ggsave(file.path(output_dir, "03b_overlaid_expression_distributions.pdf"), p_expr_overlay, width = 12, height = 5)
cat("  Saved: 03b_overlaid_expression_distributions.pdf\n")

# ----- 7c. Centrality vs Connectivity plots -----
cat("Generating centrality vs connectivity plots...\n")

centrality_conn_plots <- list(
  # Control: Eigenvector centrality vs Weighted connectivity
  ggplot(combined %>% filter(!is.na(eigenvector_centrality_ctrl) & !is.na(weighted_connectivity_ctrl)), 
         aes(x = weighted_connectivity_ctrl, y = eigenvector_centrality_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = paste0(args$ctrl_label, ": Eigenvector Centrality vs Weighted Connectivity"),
         x = "Weighted Connectivity", y = "Eigenvector Centrality") +
    theme_minimal() +
    theme(legend.position = "right"),
  
  # Treatment: Eigenvector centrality vs Weighted connectivity (colored by Control module)
  ggplot(combined %>% filter(!is.na(eigenvector_centrality_treat) & !is.na(weighted_connectivity_treat)), 
         aes(x = weighted_connectivity_treat, y = eigenvector_centrality_treat, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = paste0(args$treat_label, ": Eigenvector Centrality vs Weighted Connectivity"),
         subtitle = "(colored by Control module assignment)",
         x = "Weighted Connectivity", y = "Eigenvector Centrality") +
    theme_minimal() +
    theme(legend.position = "right"),
  
  # Control: Betweenness centrality vs Weighted connectivity
  ggplot(combined %>% filter(!is.na(betweenness_centrality_ctrl) & !is.na(weighted_connectivity_ctrl)), 
         aes(x = weighted_connectivity_ctrl, y = betweenness_centrality_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = paste0(args$ctrl_label, ": Betweenness Centrality vs Weighted Connectivity"),
         x = "Weighted Connectivity", y = "Betweenness Centrality") +
    theme_minimal() +
    theme(legend.position = "right"),
  
  # Treatment: Betweenness centrality vs Weighted connectivity (colored by Control module)
  ggplot(combined %>% filter(!is.na(betweenness_centrality_treat) & !is.na(weighted_connectivity_treat)), 
         aes(x = weighted_connectivity_treat, y = betweenness_centrality_treat, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = paste0(args$treat_label, ": Betweenness Centrality vs Weighted Connectivity"),
         subtitle = "(colored by Control module assignment)",
         x = "Weighted Connectivity", y = "Betweenness Centrality") +
    theme_minimal() +
    theme(legend.position = "right")
)

p_centrality_conn <- do.call(gridExtra::grid.arrange, c(centrality_conn_plots, ncol = 2))
ggsave(file.path(output_dir, "03c_centrality_vs_connectivity.pdf"), p_centrality_conn, width = 14, height = 10)
cat("  Saved: 03c_centrality_vs_connectivity.pdf\n")

# ----- 8. Scatterplots: Ctrl vs Treat network metrics -----
cat("Generating scatterplots (Ctrl vs Treat)...\n")

# Control module coloring (reference modules from Control diet)
scatter_plots <- lapply(pca_vars, function(var) {
  var_ctrl <- paste0(var, "_ctrl")
  var_treat <- paste0(var, "_treat")
  
  scatter_data <- combined %>% 
    select(all_of(c(var_ctrl, var_treat, "module_ctrl"))) %>%
    filter(!is.na(!!sym(var_ctrl)) & !is.na(!!sym(var_treat)))
  
  ggplot(scatter_data, aes_string(x = var_ctrl, y = var_treat, color = "module_ctrl")) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", linewidth = 1) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
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

# VERSION 1: Control module coloring for both conditions (reference module tracking)
cat("  Creating Control-module reference plots...\n")
expr_net_plots_ctrl <- list(
  ggplot(combined %>% filter(!is.na(weighted_connectivity_ctrl) & !is.na(mean_expression_ctrl)), 
         aes(x = mean_expression_ctrl, y = weighted_connectivity_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = "Control: Expression vs Weighted Connectivity",
         x = "Mean Expression", y = "Weighted Connectivity") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(weighted_connectivity_treat) & !is.na(mean_expression_treat)), 
         aes(x = mean_expression_treat, y = weighted_connectivity_treat, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = "High Sugar: Expression vs Weighted Connectivity\n(colored by Control module assignment)",
         x = "Mean Expression", y = "Weighted Connectivity") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(variance_expression_ctrl) & !is.na(degree_ctrl)), 
         aes(x = variance_expression_ctrl, y = degree_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = "Control: Expression Variance vs Degree",
         x = "Expression Variance", y = "Degree") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(variance_expression_treat) & !is.na(degree_treat)), 
         aes(x = variance_expression_treat, y = degree_treat, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = "High Sugar: Expression Variance vs Degree\n(colored by Control module assignment)",
         x = "Expression Variance", y = "Degree") +
    theme_minimal(),

  ggplot(combined %>% filter(!is.na(mean_expression_ctrl) & !is.na(eigenvector_centrality_ctrl)), 
         aes(x = mean_expression_ctrl, y = eigenvector_centrality_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = "Control: Mean Expression vs Eigenvector Centrality",
         x = "Mean Expression", y = "Eigenvector Centrality") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(mean_expression_treat) & !is.na(eigenvector_centrality_treat)), 
         aes(x = mean_expression_treat, y = eigenvector_centrality_treat, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = "High Sugar: Mean Expression vs Eigenvector Centrality\n(colored by Control module assignment)",
         x = "Mean Expression", y = "Eigenvector Centrality") +
    theme_minimal(),

  ggplot(combined %>% filter(!is.na(variance_expression_ctrl) & !is.na(eigenvector_centrality_ctrl)), 
         aes(x = variance_expression_ctrl, y = eigenvector_centrality_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = "Control: Expression Variance vs Eigenvector Centrality",
         x = "Expression Variance", y = "Eigenvector Centrality") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(variance_expression_treat) & !is.na(eigenvector_centrality_treat)), 
         aes(x = variance_expression_treat, y = eigenvector_centrality_treat, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    labs(title = "High Sugar: Expression Variance vs Eigenvector Centrality\n(colored by Control module assignment)",
         x = "Expression Variance", y = "Eigenvector Centrality") +
    theme_minimal()
)

# VERSION 2: Native module assignments for side-by-side comparison (shows interesting patterns)
cat("  Creating native module assignment plots for side-by-side comparison...\n")
expr_net_plots_native <- list(
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
    labs(title = "High Sugar: Expression vs Weighted Connectivity",
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
    labs(title = "High Sugar: Expression Variance vs Degree",
         x = "Expression Variance", y = "Degree") +
    theme_minimal(),

  ggplot(combined %>% filter(!is.na(mean_expression_ctrl) & !is.na(eigenvector_centrality_ctrl)), 
         aes(x = mean_expression_ctrl, y = eigenvector_centrality_ctrl, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "Control: Mean Expression vs Eigenvector Centrality",
         x = "Mean Expression", y = "Eigenvector Centrality") +
    theme_minimal(),
  
  ggplot(combined %>% filter(!is.na(mean_expression_treat) & !is.na(eigenvector_centrality_treat)), 
         aes(x = mean_expression_treat, y = eigenvector_centrality_treat, color = module_treat)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Module") +
    labs(title = "High Sugar: Mean Expression vs Eigenvector Centrality",
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
    labs(title = "High Sugar: Expression Variance vs Eigenvector Centrality",
         x = "Expression Variance", y = "Eigenvector Centrality") +
    theme_minimal()
)

# Save Control-module reference plots (both conditions colored by Control modules)
pdf(file.path(output_dir, "05a_expression_vs_network_metrics_CTRL_reference.pdf"), width = 14, height = 6)
for (i in seq(1, length(expr_net_plots_ctrl), by = 2)) {
  plots_to_arrange <- expr_net_plots_ctrl[i:min(i+1, length(expr_net_plots_ctrl))]
  p_expr_net_page <- do.call(gridExtra::grid.arrange, c(plots_to_arrange, ncol = 2))
  print(p_expr_net_page)
}
dev.off()
cat("  Saved: 05a_expression_vs_network_metrics_CTRL_reference.pdf\n")

# Save native module assignment plots (side-by-side comparison)

pdf(file.path(output_dir, "05b_expression_vs_network_metrics_native_comparison.pdf"), width = 14, height = 6)
for (i in seq(1, length(expr_net_plots_native), by = 2)) {
  plots_to_arrange <- expr_net_plots_native[i:min(i+1, length(expr_net_plots_native))]
  p_expr_net_page <- do.call(gridExtra::grid.arrange, c(plots_to_arrange, ncol = 2))
  print(p_expr_net_page)
}
dev.off()
cat("  Saved: 05b_expression_vs_network_metrics_native_comparison.pdf\n")

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
top_expr_fold_change <- diff_metrics %>%
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
write.csv(top_expr_fold_change, file.path(output_dir, "top_expression_fold_changes.csv"), row.names = FALSE)
write.csv(top_expr_var_up, file.path(output_dir, "top_variance_expression_changes.csv"), row.names = FALSE)
write.csv(top_expr_up, file.path(output_dir, "top_expression_changes.csv"), row.names = FALSE)

cat("  Saved summary tables:\n")
cat("    - top_increased_connectivity.csv\n")
cat("    - top_decreased_connectivity.csv\n")
cat("    - top_expression_changes.csv\n\n")

# ----- 11. Volcano-like plots for network metrics changes -----
cat("Generating condition-change plots...\n")

# Control module coloring (reference modules)
volcano_plots <- list(
  ggplot(diff_metrics %>% filter(!is.na(delta_connectivity) & !is.na(delta_expression)), 
         aes(x = delta_connectivity, y = delta_expression, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "Connectivity vs Expression Changes (Ctrl -> HS)",
         x = "Change in Weighted Connectivity",
         y = "Change in Mean Expression") +
    theme_minimal(),
  
  ggplot(diff_metrics %>% filter(!is.na(delta_degree) & !is.na(delta_betweenness)), 
         aes(x = delta_degree, y = delta_betweenness, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "Degree vs Betweenness Changes (Ctrl -> HS)",
         x = "Change in Degree",
         y = "Change in Betweenness Centrality") +
    theme_minimal(),

  ggplot(diff_metrics %>% filter(!is.na(delta_expression) & !is.na(delta_variance)), 
         aes(x = delta_expression, y = delta_variance, color = module_ctrl)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = module_color_map, name = "Control Module") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "Expression vs Variance Changes (Ctrl -> HS)",
         x = "Change in Mean Expression",
         y = "Change in Expression Variance") +
    theme_minimal()
)

p_volcano <- do.call(gridExtra::grid.arrange, c(volcano_plots, ncol = 1))
ggsave(file.path(output_dir, "06_condition_changes.pdf"), p_volcano, width = 12, height = 20)
cat("  Saved: 06_condition_changes.pdf\n")

# ----- 12. Heatmaps -----
cat("Generating heatmaps...\n")

# 12a. Correlation heatmap - within each condition
cat("  Creating within-condition correlation heatmaps...\n")
cor_ctrl <- combined %>% 
  select(all_of(paste0(pca_vars, "_ctrl"))) %>%
  na.omit() %>%
  cor()
colnames(cor_ctrl) <- rownames(cor_ctrl) <- gsub("_ctrl", "", colnames(cor_ctrl))

cor_treat <- combined %>% 
  select(all_of(paste0(pca_vars, "_treat"))) %>%
  na.omit() %>%
  cor()
colnames(cor_treat) <- rownames(cor_treat) <- gsub("_treat", "", colnames(cor_treat))

pdf(file.path(output_dir, "07a_correlation_heatmap_within_conditions.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))
corrplot(cor_ctrl, method = "color", type = "upper", 
         title = paste0(args$ctrl_label, " - Network Metrics Correlation"),
         mar = c(0, 0, 2, 0), tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         addCoef.col = "black", number.cex = 0.7)
corrplot(cor_treat, method = "color", type = "upper", 
         title = paste0(args$treat_label, " - Network Metrics Correlation"),
         mar = c(0, 0, 2, 0), tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         addCoef.col = "black", number.cex = 0.7)
dev.off()
cat("    Saved: 07a_correlation_heatmap_within_conditions.pdf\n")

# 12b. Delta heatmap - changes for top genes
cat("  Creating delta heatmap for top differentially affected genes...\n")
delta_vars <- c("delta_connectivity", "delta_degree", "delta_betweenness", 
                "delta_expression", "delta_variance")

# Select top 50 genes by absolute change in connectivity
top_delta_genes <- diff_metrics %>%
  filter(!is.na(delta_connectivity)) %>%
  mutate(abs_delta_conn = abs(delta_connectivity)) %>%
  arrange(desc(abs_delta_conn)) %>%
  slice_head(n = 50)

# Create matrix for heatmap
delta_matrix <- top_delta_genes %>%
  select(gene, all_of(delta_vars)) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Scale by column (z-score) for better visualization
delta_matrix_scaled <- scale(delta_matrix)

# Create annotation for modules
annotation_row <- top_delta_genes %>%
  select(gene, module_ctrl, module_treat) %>%
  column_to_rownames("gene")

# Create color mapping for modules - only include modules that exist in the data
unique_modules_ctrl <- unique(annotation_row$module_ctrl)
unique_modules_treat <- unique(annotation_row$module_treat)

module_colors_ctrl <- module_color_map[unique_modules_ctrl]
module_colors_ctrl <- module_colors_ctrl[!is.na(module_colors_ctrl)]

module_colors_treat <- module_color_map[unique_modules_treat]
module_colors_treat <- module_colors_treat[!is.na(module_colors_treat)]

annotation_colors <- list(
  module_ctrl = module_colors_ctrl,
  module_treat = module_colors_treat
)

pdf(file.path(output_dir, "07b_delta_heatmap_top_genes.pdf"), width = 10, height = 12)
tryCatch({
  pheatmap(delta_matrix_scaled,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = "Top 50 Genes by Connectivity Change\n(Z-scored deltas)",
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           fontsize_row = 6,
           angle_col = 45)
}, error = function(e) {
  cat("    Warning: Error creating annotated heatmap, creating simplified version\n")
  pheatmap(delta_matrix_scaled,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = "Top 50 Genes by Connectivity Change\n(Z-scored deltas)",
           cluster_cols = FALSE,
           show_rownames = TRUE,
           fontsize_row = 6,
           angle_col = 45)
})
dev.off()
cat("    Saved: 07b_delta_heatmap_top_genes.pdf\n")

# 12c. Module preservation heatmap
cat("  Creating module preservation heatmap...\n")
module_transition <- combined %>%
  filter(!is.na(module_ctrl) & !is.na(module_treat)) %>%
  group_by(module_ctrl, module_treat) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = module_treat, values_from = count, values_fill = 0) %>%
  column_to_rownames("module_ctrl") %>%
  as.matrix()

# Get all unique modules and sort them consistently
all_modules_sorted <- sort(unique(c(rownames(module_transition), colnames(module_transition))))

# Reorder matrix to have same order on both axes
module_transition_ordered <- matrix(0, 
                                    nrow = length(all_modules_sorted), 
                                    ncol = length(all_modules_sorted),
                                    dimnames = list(all_modules_sorted, all_modules_sorted))

# Fill in the values
for (i in rownames(module_transition)) {
  for (j in colnames(module_transition)) {
    module_transition_ordered[i, j] <- module_transition[i, j]
  }
}

# Calculate preservation percentage (row-wise)
module_preservation_pct <- sweep(module_transition_ordered, 1, rowSums(module_transition_ordered), "/") * 100

pdf(file.path(output_dir, "07c_module_preservation_heatmap.pdf"), width = 10, height = 10)
pheatmap(module_preservation_pct,
         color = colorRampPalette(c("white", "#4575B4", "#313695"))(100),
         main = "Module Preservation: Control (bottom)→ Treatment (right)\n(% of genes from ctrl module assigned to treat module)",
         display_numbers = matrix(sprintf("%.0f%%", module_preservation_pct), 
                                 nrow = nrow(module_preservation_pct)),
         number_color = "black",
         fontsize_number = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         angle_col = 45)
dev.off()
cat("    Saved: 07c_module_preservation_heatmap.pdf\n")

# 12d. Cross-condition correlation heatmap
cat("  Creating cross-condition correlation heatmap...\n")
# For each ctrl metric, correlate with each treat metric
cross_cor_matrix <- matrix(NA, nrow = length(pca_vars), ncol = length(pca_vars),
                           dimnames = list(paste0(pca_vars, "\n(", args$ctrl_label, ")"), 
                                         paste0(pca_vars, "\n(", args$treat_label, ")")))

for (i in 1:length(pca_vars)) {
  for (j in 1:length(pca_vars)) {
    ctrl_col <- paste0(pca_vars[i], "_ctrl")
    treat_col <- paste0(pca_vars[j], "_treat")
    
    valid_data <- combined %>%
      select(all_of(c(ctrl_col, treat_col))) %>%
      na.omit()
    
    if (nrow(valid_data) > 0) {
      cross_cor_matrix[i, j] <- cor(valid_data[[ctrl_col]], valid_data[[treat_col]])
    }
  }
}

pdf(file.path(output_dir, "07d_cross_condition_correlation_heatmap.pdf"), width = 10, height = 8)
pheatmap(cross_cor_matrix,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "Cross-Condition Metric Correlations\n(How do metrics in one condition relate to metrics in the other?)",
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 10,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         breaks = seq(-1, 1, length.out = 101),
         legend_breaks = c(-1, -0.5, 0, 0.5, 1),
         angle_col = 0,
         fontsize = 10)
dev.off()
cat("    Saved: 07d_cross_condition_correlation_heatmap.pdf\n")

# 12e. Module-specific summary heatmap
cat("  Creating module-specific metric summary heatmap...\n")
# Calculate mean values per module for each metric
module_summary_ctrl <- combined %>%
  filter(!is.na(module_ctrl)) %>%
  group_by(module_ctrl) %>%
  summarise(
    connectivity = mean(weighted_connectivity_ctrl, na.rm = TRUE),
    degree = mean(degree_ctrl, na.rm = TRUE),
    betweenness = mean(betweenness_centrality_ctrl, na.rm = TRUE),
    clustering = mean(clustering_coefficient_ctrl, na.rm = TRUE),
    eigenvector = mean(eigenvector_centrality_ctrl, na.rm = TRUE),
    expression = mean(mean_expression_ctrl, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  column_to_rownames("module_ctrl") %>%
  as.matrix()

module_summary_treat <- combined %>%
  filter(!is.na(module_treat)) %>%
  group_by(module_treat) %>%
  summarise(
    connectivity = mean(weighted_connectivity_treat, na.rm = TRUE),
    degree = mean(degree_treat, na.rm = TRUE),
    betweenness = mean(betweenness_centrality_treat, na.rm = TRUE),
    clustering = mean(clustering_coefficient_treat, na.rm = TRUE),
    eigenvector = mean(eigenvector_centrality_treat, na.rm = TRUE),
    expression = mean(mean_expression_treat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  column_to_rownames("module_treat") %>%
  as.matrix()

# Get all unique modules and ensure both matrices have the same rows
all_modules_for_summary <- sort(unique(c(rownames(module_summary_ctrl), rownames(module_summary_treat))))

# Reorder/expand matrices to have same rows
reorder_matrix <- function(mat, target_rows) {
  result <- matrix(NA, nrow = length(target_rows), ncol = ncol(mat),
                   dimnames = list(target_rows, colnames(mat)))
  common_rows <- intersect(rownames(mat), target_rows)
  result[common_rows, ] <- as.matrix(mat[common_rows, ])
  return(result)
}

module_summary_ctrl_ordered <- reorder_matrix(module_summary_ctrl, all_modules_for_summary)
module_summary_treat_ordered <- reorder_matrix(module_summary_treat, all_modules_for_summary)

# Scale each metric (z-score by column) - handle NAs
module_summary_ctrl_scaled <- scale(module_summary_ctrl_ordered)
module_summary_treat_scaled <- scale(module_summary_treat_ordered)

pdf(file.path(output_dir, "07e_module_summary_heatmap.pdf"), width = 12, height = 10)
par(mfrow = c(1, 2))

pheatmap(module_summary_ctrl_scaled,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = paste0(args$ctrl_label, " - Module Metric Averages\n(Z-scored)"),
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 8,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         angle_col = 45,
         na_col = "grey90")

pheatmap(module_summary_treat_scaled,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = paste0(args$treat_label, " - Module Metric Averages\n(Z-scored)"),
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 8,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         angle_col = 45,
         na_col = "grey90")

dev.off()
cat("    Saved: 07e_module_summary_heatmap.pdf\n")

# 12f. Expression-Network Integration bubble plot
cat("  Creating expression-network integration bubble plot...\n")
bubble_data <- diff_metrics %>%
  filter(!is.na(delta_expression) & !is.na(delta_connectivity) & !is.na(delta_variance)) %>%
  mutate(
    abs_delta_variance = abs(delta_variance),
    direction = case_when(
      delta_expression > 0 & delta_connectivity > 0 ~ "Both Increased",
      delta_expression > 0 & delta_connectivity < 0 ~ "Expr Up, Connect Down",
      delta_expression < 0 & delta_connectivity > 0 ~ "Expr Down, Connect Up",
      delta_expression < 0 & delta_connectivity < 0 ~ "Both Decreased",
      TRUE ~ "No Change"
    )
  )

# Control module coloring (reference modules)
p_bubble <- ggplot(bubble_data, aes(x = delta_expression, y = delta_connectivity, 
                                     size = abs_delta_variance, color = module_ctrl)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = module_color_map, name = "Control Module") +
  scale_size_continuous(name = "Abs Variance Change", range = c(1, 10)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(title = "Expression-Network Integration",
       subtitle = "Bubble size = magnitude of expression variance change",
       x = "Change in Mean Expression (delta)",
       y = "Change in Weighted Connectivity (delta)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "07f_expression_network_integration.pdf"), 
       p_bubble, width = 12, height = 8)
cat("    Saved: 07f_expression_network_integration.pdf\n")

# ----- 13. Summary statistics -----
cat("\n=== Summary Statistics ===\n")

cat("\nNetwork Metrics Comparison:\n")
for (var in pca_vars) {
  var_ctrl <- paste0(var, "_ctrl")
  var_treat <- paste0(var, "_treat")
  
  ctrl_vals <- combined[[var_ctrl]]
  treat_vals <- combined[[var_treat]]
  
  cat(sprintf("%-30s: Ctrl mean=%.4f, Treat mean=%.4f, delta=%.4f\n",
              var,
              mean(ctrl_vals, na.rm = TRUE),
              mean(treat_vals, na.rm = TRUE),
              mean(treat_vals - ctrl_vals, na.rm = TRUE)))
}

cat("\nExpression Statistics Comparison:\n")
expr_vars <- c("mean_expression", "variance_expression", "mad_expression")
for (var in expr_vars) {
  var_ctrl <- paste0(var, "_ctrl")
  var_treat <- paste0(var, "_treat")
  
  ctrl_vals <- combined[[var_ctrl]]
  treat_vals <- combined[[var_treat]]
  
  cat(sprintf("%-30s: Ctrl mean=%.4f, Treat mean=%.4f, delta=%.4f\n",
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
