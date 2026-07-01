library(openxlsx)
library(VennDiagram)
library(data.table)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

mad_shift_genes <- read.xlsx("results/variability/voomhs_all_genes_mad_summary.xlsx")
permu_shift_genes <- read.xlsx("results/variability/mad_permutation_hs.xlsx")

setDT(mad_shift_genes)
setDT(permu_shift_genes)

up_shift <- list(normal=mad_shift_genes[wilcoxon_p_adj < 0.05 & transcriptomic_mad_shift > 0, gene_id],
                permute = permu_shift_genes[ perm_p_adj < 0.05 & obs_delta_mean_mad > 0, gene_id])
down_shift <- list(normal=mad_shift_genes[wilcoxon_p_adj < 0.05 & transcriptomic_mad_shift < 0, gene_id],
                permute = permu_shift_genes[ perm_p_adj < 0.05 & obs_delta_mean_mad < 0, gene_id])

venn.diagram(down_shift, filename = "results/variability/venn_mad_downshift_genes_normal_vs_permute.png",
            fill=myCol[1:2], output=TRUE)
venn.diagram(up_shift, filename = "results/variability/venn_mad_upshift_genes_normal_vs_permute.png",
            fill=myCol[1:2], output=TRUE)