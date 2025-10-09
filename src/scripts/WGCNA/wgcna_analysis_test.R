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

# Load required packages
library(WGCNA)
library(yaml)

# Allow multi-threading 
enableWGCNAThreads()

# Load configuration
config <- yaml::read_yaml("config/config.yaml")

# Options for the analysis
options(stringsAsFactors = FALSE)

# Load test data for VST control
data_file <- file.path(config$project_dirs$test_data_dir, "VSTdataCtrl_test.txt")
data <- read.table(data_file, header=TRUE, sep="\t", row.names=1)

# Transpose data for WGCNA format (samples as rows)
datExpr <- as.data.frame(t(data))

# Check for genes and samples with too many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  # Remove the offending genes and samples
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Choose a set of soft-thresholding powers
powers <- c(1:20)

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results
pdf(file = file.path(config$output_dirs$figures_dir, "WGCNA", "scale_independence.pdf"), width = 9, height = 5)
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
softPower <- # Change this based on your results

# Calculate adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Set the minimum module size
minModuleSize <- 5

# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)

# Plot the dendrogram and colors underneath
pdf(file = file.path(config$output_dirs$figures_dir, "WGCNA", "module_dendrogram_test.pdf"), width = 12, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05,
                   main = "Gene dendrogram and module colors")
dev.off()

# Save the results
output_dir <- file.path(config$output_dirs$wgcna_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save module colors and tree
save(dynamicMods, dynamicColors, geneTree, 
     file = file.path(output_dir$wgcna_dir, "network_analysis_results.RData"))

# Create a data frame with gene information and module assignment
geneInfo <- data.frame(
  Gene_ID = colnames(datExpr),
  Module = dynamicColors
)

# Write the gene module assignments to a file
write.table(geneInfo, 
           file = file.path(output_dir, "gene_module_assignment_test.txt"),
           sep = "\t", row.names = FALSE, quote = FALSE)

cat("WGCNA analysis complete. Results saved to:", output_dirs$wgcna_dir, "\n")