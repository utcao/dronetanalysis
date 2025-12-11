Identifying global modifiers using co-expression network analysis

WGCNA pipeline script order:

1. Preprocess data
src/scripts/preprocess/FinalProcessedDatasets.R

2. Subset dataset for testing (optional)
src/scripts/subset/01subset_dataset.R

3. Spearman correlation matrix generation
src/scripts/spearman_corr/spearman_matrix.R

4. Adjacency thresholding
src/scripts/soft_threshold/adjacency_matrix.R

6. Calculate network-level features and WGCNA modules
src/scripts/network_metrics/nm_and_modules.R

7. Calculate gene-level features
src/scripts/gene_metrics/gene_level_metrics.R

8. Concatenate with expression data
src/scripts/gene_metrics/add_expression_data.R

9. Compare gene metrics between diets
src/scripts/analysis/compare_conditions.R

10. Perform GO enrichment analysis
src/scripts/analysis/module_enrichment.R

11. Visualise the networks with igraph
src/scripts/analysis/visualise_network.R

12. NetRep module preservation
src/scripts/diff_matrix/module_preservation_netrep.R

Quintile analysis-specific scripts:

13 Subsetting of samples by gene expression
src/scripts/preprocess/subset_expr_by_gene_quintile.R

14.a) Adjacency difference matrix
src/scripts/diff_matrix/diff_adjacency.R
14.b) TOM difference matrix
src/scripts/diff_matrix/diff_TOM.R

15. Subset neighbours based on focal gene features
src/scripts/subset_features/subset_neighbour_features.R

16. Snakemake pipeline for quintile analysis
src/pipelines/Snakemake_template
