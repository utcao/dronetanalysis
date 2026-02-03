Identifying global modifiers using co-expression network analysis

WGCNA pipeline script order:

1. Preprocess data
src/scripts/00preprocess/FinalProcessedDatasets.R

2. Subset dataset for testing (optional)
src/scripts/01subset/02get_extreme_pop.R

3. Spearman correlation matrix generation
src/scripts/10spearman_corr/spearman_matrix.R

4. Adjacency thresholding
src/scripts/10wgcna/adjacency_matrix.R

6. Calculate network-level features and WGCNA modules
src/scripts/12network_metrics/network_modules_wgcna.R

7. Calculate gene-level features
src/scripts/13gene_metrics/gene_level_metrics.R

8. Concatenate with expression data
src/scripts/13gene_metrics/add_expression_stats.R

9. Compare gene metrics between diets
src/scripts/15analysis/compare_conditions.R

10. Perform GO enrichment analysis
src/scripts/15analysis/module_enrichment.R

11. Visualise the networks with igraph
src/scripts/15analysis/visualise_network.R

12. NetRep module preservation
src/scripts/15analysis/module_preservation_netrep.R

Quintile analysis-specific scripts:

13 Subsetting of samples by gene expression
src/scripts/00preprocess/subset_expr_by_gene_quintile.R

14.a) Adjacency difference matrix
src/scripts/11diff_matrix/diff_adjacency.R
14.b) TOM difference matrix
src/scripts/11diff_matrix/diff_tom.R

15. Subset neighbours based on focal gene features
src/scripts/14subset_features/subset_neighbour_features.R

Bootstrap co-expression pipeline (SGE two-stage):

1. src/scripts/01subset/01get_extreme_pop_bootstrap.py   (Stage 1 – generate indices)
2. src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py  (Stage 2 – per-gene array)
3. src/SGE_scripts/run_bootstrap_pipeline.sh              (submission script)

16. Snakemake pipeline for quintile analysis
src/pipelines/Snakemake_template
