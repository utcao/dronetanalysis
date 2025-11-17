Identifying global modifiers using co-expression network analysis

WGCNA pipeline script order:

1. Preprocess data
scr/preprocess/FinalProcessedDatasets.R

2. Subset dataset
scr/subset/01subset_dataset.R

3. Spearman correlation matrix generation
scr/spearman_corr/spearman_matrix.R

4. Adjacency threshold plot analysis
scr/soft_threshold/adjacency_matrix.R

5. Apply soft threshold after inspecting plots
scr/soft_threshold/apply_threshold.R

6. Calculate network-level features
scr/network_metrics/network_metrics.R

7. Calculate gene-level features
scr/gene_metrics/gene_level_metrics.R

8. Concatenate with expression data
scr/gene_metrics/add_expression_data.R

9. Compare gene metrics between diets
scr/analysis/compare_conditions.R

10. Perform GO enrichment analysis
scr/analysis/module_enrichment.R

11. Visualise the networks with igraph
scr/analysis/visualise_network.R