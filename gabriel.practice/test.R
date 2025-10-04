#!/usr/bin/env Rscript
source("dronetanalysis/src/utils/utils_io.R")
# -----  Extract Key Settings from config.yaml -----
library(yaml)  # to read the YAML config
config <- yaml::read_yaml("config/config.yaml")
install.packages("data.table")
suppressPackageStartupMessages({
    library(purrr)
    library(data.table)
    library(glue)
    library(stringr)
    library(ggplot2)
    library(tibble)  # for rownames_to_column / column_to_rownames
})
library(data.table)

#read table (genes in first column)
dt <- fread("logCPM_Ctrl_Dros.csv")   # adjust filename
setnames(dt, 1, "gene")     # ensure first col named "gene"

#quick inspection
# per-sample ranges & how many negatives
summary_dt <- dt[, lapply(.SD, function(x) list(min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE))), .SDcols = -1]
print(summary_dt)
neg_counts <- dt[, lapply(.SD, function(x) sum(x < 0, na.rm=TRUE)), .SDcols = -1]
print(neg_counts)

#convert to matrix (genes x samples) for numeric ops
mat <- as.matrix(dt[, -1, with = FALSE])

#filter low-variance genes (important before correlation/network building)
gene_var <- rowVars <- apply(mat, 1, var, na.rm = TRUE)    # base R var per gene
keep <- which(gene_var > quantile(gene_var, 0.25))         # keep top 75% var genes (example)
mat_filt <- mat[keep, , drop = FALSE]

#optional: z-score per gene (centers each gene to mean 0, sd 1) — helps comparability
mat_z <- t(scale(t(mat_filt), center = TRUE, scale = TRUE))

#subset only a few genes
mat_z_sub<-mat_z[1:100,]

#compute pairwise Pearson correlations (genes x genes)
cor_mat <- cor(t(mat_z_sub), use = "pairwise.complete.obs", method = "pearson")

#convert to long data.table of unique pairs (upper triangle only)
library(data.table)
cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("geneA", "geneB", "correlation"))
# keep only one copy of each pair and remove self-self
cor_dt <- cor_dt[as.character(geneA) < as.character(geneB)]

#filter strong edges (e.g., abs(cor) > 0.7) and get top partners per gene
edges=cor_dt[,abs_corr :=abs(correlation)>0.7]
# top 3 partners for each geneA (example)
setorder(edges, geneA, -abs_corr)
top3 <- edges[, head(.SD, 3), by = geneA]