# voom/logCPM matrix and vst matrix acquisition for co-expression network analysis

# Script by James Tan
# Edited by Gabriel Thornes

# Last Updated: 7/10/25

#################################
##### Packages and Setup ########
#################################

rm(list = ls())

library(DESeq2)
library(dplyr)
library(edgeR)
library(limma)
library(tidyr) 
library(yaml)
config <- yaml::read_yaml("config/config.yaml")
rawcount_file <- config$raw_files$rawcount_file
rawmeta_file <- config$raw_files$rawmeta_file
voom_surrogate_files <- config$raw_files$voom_surrogate_files
vst_surrogate_files <- config$raw_files$vst_surrogate_files

#####################
##### Datasets ######
#####################

# Metadata on all samples, e.g., conditions and batches
head.info   <- read.table(rawmeta_file, h=T)

# Surrogate variables, estimated in the NvsHS_SurrogateVariables R scripts
tmm.voom.sv1 <- read.table(voom_surrogate_files[[1]])
tmm.voom.sv1 <- tmm.voom.sv1$x
tmm.voom.sv2 <- read.table(voom_surrogate_files[[2]])
tmm.voom.sv2 <- tmm.voom.sv2$x
tmm.voom.sv3 <- read.table(voom_surrogate_files[[3]])
tmm.voom.sv3 <- tmm.voom.sv3$x
tmm.voom.sv4 <- read.table(voom_surrogate_files[[4]])
tmm.voom.sv4 <- tmm.voom.sv4$x

vst.sv1 <- read.table(vst_surrogate_files[[1]])
vst.sv1 <- vst.sv1$x
vst.sv2 <- read.table(vst_surrogate_files[[2]])
vst.sv2 <- vst.sv2$x
vst.sv3 <- read.table(vst_surrogate_files[[3]])
vst.sv3 <- vst.sv3$x

# Based on PC inspection - Include 3 SVs for VST; 4 SVs for voom

# Expression - raw counts of filtered samples (see MakingGeneExpressionMatrix_head_HS&CTRL.R for filtering criteria)
raw.counts <- read.table(rawcount_file, h=T, check.names = F)

# Eliminate bimodal genes
Bimodal_genes <- read.csv('dataset/bimodal_genes/Bimodal_Genes.csv',row.names = 1) 
Non_bimodal_genes <- setdiff(rownames(raw.counts),rownames(Bimodal_genes))
raw.counts <- raw.counts[Non_bimodal_genes,]

# Only metadata on samples with counts
cov <- tibble::column_to_rownames(head.info, "id")
cov2 <- cov[c(colnames(raw.counts)),]

# Merge expression and metadata 
raw.counts.t <- data.frame(t(raw.counts)) 
raw.counts.t <- tibble::rownames_to_column(raw.counts.t, "id")
head.data <- merge(raw.counts.t, head.info, by = 'id')
head.data <- tibble::column_to_rownames(head.data, "id")

############################
##### Transformations ######
############################

# Preparing independent variable and batch variables
cov2$treatment <- as.factor(cov2$treatment)
cov2$RNAlibBatch <- as.factor(cov2$RNAlibBatch)
cov2$RNAseqBatch <- as.factor(cov2$RNAseqBatch)
cov2$egglayBatch<- as.factor(cov2$egglayBatch)
cov2$platingBatch <- as.factor(cov2$platingBatch)
cov2$well <- as.factor(cov2$well)

treatment <- cov2$treatment
rnalib <- cov2$RNAlibBatch
rnaseq <- cov2$RNAseqBatch
egg <-cov2$egglayBatch
plate <-cov2$platingBatch
well <-cov2$well

contrasts(rnalib) <- contr.sum(levels(rnalib)) # sums to zero contrast
contrasts(rnaseq) <- contr.sum(levels(rnaseq))
contrasts(egg) <- contr.sum(levels(egg))
contrasts(plate) <- contr.sum(levels(plate))
contrasts(well) <- contr.sum(levels(well))

# nullmodel - the biological factor of interest
nullmodel <- model.matrix(~treatment)

# Defining indices to subset the total dataset into Ctrl and HS downstream
conditions <- factor(t(head.data$treatment))
conditionsLevel <- levels(conditions)

############################
##### TMM Voom #############
############################

# Calculate TMM factors 
raw.counts.list <- DGEList(counts=raw.counts)
TMM.counts <- calcNormFactors(raw.counts.list, method = 'TMM') 

# Transform count data while keeping design
# First log2(cpm+0.5). Then it does voom keeping covariates separated
TMM.Voom.design <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4+treatment)
voom.object <- voom(TMM.counts, design = TMM.Voom.design, plot=T)

# For tests that do not accept batch effects as covariates:
# Correct counts for weights and batch effects but not for treatment
TMM.Voom.covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well+tmm.voom.sv1+tmm.voom.sv2+tmm.voom.sv3+tmm.voom.sv4)
voom.counts.bc <- limma::removeBatchEffect(voom.object, 
                                           covariates = TMM.Voom.covariates[,-1],  
                                           design = nullmodel) 
voom.counts.bc <- as.data.frame(voom.counts.bc)

# Split data by control or high sugar
voomdataN <- voom.counts.bc[,c(which(conditions==conditionsLevel[1]))]
voomdataHS <- voom.counts.bc[,c(which(conditions==conditionsLevel[2]))]

# Write VOOM data to tables
write.table(voomdataHS, 
           file=file.path(config$project_dirs$processed_data_dir, "VOOM", "voomdataHS.txt"),
           sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(voomdataN, 
           file=file.path(config$project_dirs$processed_data_dir, "VOOM", "voomdataN.txt"),
           sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#####################
##### VST ###########
#####################


# Apply variance-stabilizing transformation
data <- DESeqDataSetFromMatrix(countData = raw.counts,
                               colData = cov2,
                               design =  ~treatment) 
vst.data = vst(data, blind=F) 


# For tests that do not accept batch effects as covariates:
# Correct counts for batch effects but not for treatment
vst.counts <- assay(vst.data)

# Remove batch effects 
VST.covariates <- model.matrix(~rnalib+rnaseq+egg+plate+well+vst.sv1+vst.sv2+vst.sv3)
vst.counts.bc  <- limma::removeBatchEffect(vst.counts, covariates = VST.covariates[,-1],  design = nullmodel) 
vst.counts.bc.t <- t(vst.counts.bc)

# Split data by control or high sugar
VSTdataN <- vst.counts.bc[,c(which(conditions==conditionsLevel[1]))]
VSTdataHS <- vst.counts.bc[,c(which(conditions==conditionsLevel[2]))]

# Write VST data to tables
write.table(VSTdataHS, 
           file=file.path(config$project_dirs$processed_data_dir, "VST", "VSTdataHS.txt"),
           sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(VSTdataN, 
           file=file.path(config$project_dirs$processed_data_dir, "VST", "VSTdataN.txt"),
           sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)