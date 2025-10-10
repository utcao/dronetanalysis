# !/usr/bin/env Rscript
# ==============================================================================
# Weighted Gene Co-Expression Network Analysis (WGCNA) for drosophila transcriptomic data
#
# Script by Gabriel Thornes
#
# Last Updated: 09/10/2025
#
# This script::
#   1. Takes VST test datasets as input
#   2. Performs WGCNA analysis
#   3. Generates plots and saves results
# ==============================================================================

#################################
##### Packages and Setup ########
#################################

source("src/utils/utils_io.R")

# Load required packages
library(WGCNA)
library(yaml)

# Allow multi-threading 
enableWGCNAThreads()

# Load configuration
config <- yaml::read_yaml("config/config.yaml")

# Options for the analysis
options(stringsAsFactors = FALSE)

# Choose a set of soft-thresholding powers
powers <- c(1:20)

# Create output directories if they don't exist  
dir.create(config$output_dirs$wgcna_dir, recursive = TRUE, showWarnings = FALSE)

###############################
##### Data Loading ############
###############################

# Load test data for VST control
data_file <- file.path(config$project_dirs$subset_data_dir, "VSTdataCtrl_subset.txt")
data <- read.table(data_file, header=TRUE, sep="\t", row.names=1)

# Transpose data for WGCNA format (samples as rows)
datExpr <- as.data.frame(t(data))

# Store original gene names before any filtering
original_gene_names <- colnames(datExpr)

# Check for genes and samples with too many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  # Remove the offending genes and samples
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  # Update the gene names list to match filtered data
  original_gene_names <- original_gene_names[gsg$goodGenes]
}

###############################
##### Network Topology ########
###############################

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results - save to the WGCNA directory
pdf(file = file.path(config$output_dirs$wgcna_dir, "scale_independence.pdf"), width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.9, col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

# Choose an appropriate soft-thresholding power based on the plot (R^2 > 0.9)
softPower <- 5 # Change this based on your results

# Calculate adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average")

###############################
##### Module Construction #####
###############################

# Set the minimum module size
minModuleSize <- 30

# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)

# Calculate eigengenes for initial modules
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the result
pdf(file = file.path(config$output_dirs$wgcna_dir, "module_eigengene_clustering.pdf"), width = 12, height = 8)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Set the merge threshold (height at which to cut the tree)
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors <- merge$colors

# Eigengenes of the new merged modules
mergedMEs <- merge$newMEs

# Plot comparison of original vs merged modules
pdf(file = file.path(config$output_dirs$wgcna_dir, "module_comparison.pdf"), width = 12, height = 8)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                   c("Dynamic Tree Cut", "Merged dynamic"),
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05,
                   main = "Gene dendrogram and module colours before and after merging")
dev.off()

# Use the merged colors for subsequent analyses
moduleColors <- mergedColors

# Plot the dendrogram and colors underneath
pdf(file = file.path(config$output_dirs$wgcna_dir, "module_dendrogram_subset.pdf"), width = 12, height = 6)
plotDendroAndColors(geneTree, moduleColors, "Merged Modules",
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05,
                   main = "Gene dendrogram and final module colours")
dev.off()

# Save the results
output_dir <- file.path(config$output_dirs$wgcna_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save module colors and tree
save(dynamicMods, dynamicColors, moduleColors, mergedMEs, geneTree, 
     file = file.path(output_dir, "network_analysis_results.RData"))

# Create a data frame with gene information and module assignment
geneInfo <- data.frame(
  Gene_ID = original_gene_names,  # Use the preserved original gene names
  Original_Module = dynamicColors,
  Final_Module = moduleColors
)

# Write the gene module assignments to a file
write.table(geneInfo, 
           file = file.path(output_dir, "gene_module_assignment_test.txt"),
           sep = "\t", row.names = FALSE, quote = FALSE)

#################################
##### Hub Gene Identification ###
#################################

# Identify hub genes for each module
uniqueModules <- unique(moduleColors)

# Print module summary information
cat("\n=== MODULE SUMMARY ===\n")
cat("Number of modules identified (original):", length(unique(dynamicColors)), "\n")
cat("Number of modules after merging:", length(uniqueModules), "\n")
cat("Module colours:", paste(uniqueModules, collapse = ", "), "\n")

# Count genes per module
module_counts <- table(moduleColors)
cat("\nGenes per module:\n")
print(module_counts)

# Function to identify hub genes based on intramodular connectivity
hubGeneList <- sapply(uniqueModules, function(mod) {
  # Get genes in this module
  genes_in_module <- which(moduleColors == mod)
  if(length(genes_in_module) == 0) return(NA)
  
  # Get the module's expression data
  module_expr <- datExpr[, genes_in_module, drop = FALSE]
  
  # Calculate connectivity for each gene in the module
  gene_connectivity <- apply(module_expr, 2, function(x) {
    sum(abs(cor(x, module_expr, use = "complete.obs")))
  })
  
  # Return the gene with highest connectivity (hub gene)
  hub_gene_index <- which.max(gene_connectivity)
  hub_gene_name <- colnames(module_expr)[hub_gene_index]
  
  return(hub_gene_name)
})

# Set names for the hub gene list
names(hubGeneList) <- uniqueModules

hubGeneList <- na.omit(hubGeneList)
cat("\nModules with valid hub genes:", length(hubGeneList), "out of", length(uniqueModules), "\n")

hubGeneDF <- data.frame(
  Module = names(hubGeneList), 
  Hub_Gene = as.character(hubGeneList)
)

# Write the hub genes to a CSV file for ease of gene_id acquisition
write.table(hubGeneDF, 
           file = file.path(output_dir, "hub_genes_test.csv"),
           sep = ",", row.names = FALSE, quote = FALSE)

cat("WGCNA analysis complete. Results saved to:", output_dir, "\n")
