#!/usr/bin/env Rscript
# ==============================================================================
# Per-Gene Expression Statistics: Full CT + HS Matrix (condition-specific)
#
# Loads CT and HS VOOM expression matrices, computes per-gene statistics
# separately for each condition and for the combined full matrix:
#
#   CT-specific:   mean_ct, median_ct, sd_ct, mad_ct, cv2_ct
#   HS-specific:   mean_hs, median_hs, sd_hs, mad_hs, cv2_hs
#   Cross-condition differences (HS − CT): mad_hs_minus_ct, cv2_hs_minus_ct
#   Full CT+HS:    mean_full, median_full, sd_full, mad_full, cv2_full
#   Ranks (full):  rank_mean, rank_median, rank_mad, rank_cv2
#
# CT and HS annotation files are merged (union by ENSEMBL/gene_id) so that
# genes annotated in only one condition still receive a SYMBOL.
#
# Usage:
#   Rscript compute_full_mad_cv2_ranks.R \
#     --ct-file         data/processed/VOOM/voomdataCtrl.txt \
#     --hs-file         data/processed/VOOM/voomdataHS.txt \
#     --ct-mapping-file run_voomct/results_ct_voom/rewiring_hubs_ct_anno_0408_2026.tsv \
#     --hs-mapping-file run_voomhs/results_hs_voom/rewiring_hubs_hs_anno_0413_2026.tsv \
#     --output-file     results/variability/full_mad_cv2_ranks.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
  library(openxlsx)
})

# ----- CLI arguments -----
parser <- ArgumentParser(description = "Per-gene expression stats (CT, HS, and full matrix)")
parser$add_argument("--ct-file",          required = TRUE,
                    help = "CT VOOM expression matrix (tab-separated, genes x samples)")
parser$add_argument("--hs-file",          required = TRUE,
                    help = "HS VOOM expression matrix (tab-separated, genes x samples)")
parser$add_argument("--ct-mapping-file",  required = TRUE,
                    help = "CT annotation TSV with ENSEMBL and SYMBOL columns")
parser$add_argument("--hs-mapping-file",  required = TRUE,
                    help = "HS annotation TSV with ENSEMBL and SYMBOL columns")
parser$add_argument("--output-file",      required = TRUE,
                    help = "Output .xlsx file path (will NOT overwrite if exists)")
args <- parser$parse_args()

cat("=== Per-Gene Expression Statistics: CT, HS, and Full Matrix ===\n")
cat("CT file:          ", args$ct_file,          "\n")
cat("HS file:          ", args$hs_file,          "\n")
cat("CT mapping file:  ", args$ct_mapping_file,  "\n")
cat("HS mapping file:  ", args$hs_mapping_file,  "\n")
cat("Output file:      ", args$output_file,      "\n\n")

# ----- Guard: refuse to overwrite -----
if (file.exists(args$output_file))
  stop("Output file already exists: ", args$output_file,
       "\nRename or remove it before re-running.")

# ----- Validate inputs -----
for (f in c(args$ct_file, args$hs_file, args$ct_mapping_file, args$hs_mapping_file)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

out_dir <- dirname(args$output_file)
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

# ----- Load matrices -----
cat("Loading CT matrix...\n")
ct_dt <- fread(args$ct_file, data.table = TRUE)
gene_col <- colnames(ct_dt)[1]
cat("  Gene ID column: '", gene_col, "'\n", sep = "")
cat("  Dimensions:", nrow(ct_dt), "genes x", ncol(ct_dt) - 1, "samples\n")

cat("Loading HS matrix...\n")
hs_dt <- fread(args$hs_file, data.table = TRUE)
cat("  Dimensions:", nrow(hs_dt), "genes x", ncol(hs_dt) - 1, "samples\n\n")

# Ensure same gene order; merge on gene_id column
ct_genes <- ct_dt[[gene_col]]
hs_genes <- hs_dt[[gene_col]]
if (!identical(ct_genes, hs_genes)) {
  cat("Gene order differs — aligning by gene_id...\n")
  common_genes <- intersect(ct_genes, hs_genes)
  ct_dt <- ct_dt[match(common_genes, ct_genes)]
  hs_dt <- hs_dt[match(common_genes, hs_genes)]
  cat("  Common genes:", length(common_genes), "\n\n")
}

gene_ids   <- ct_dt[[gene_col]]
ct_samples <- setdiff(colnames(ct_dt), gene_col)
hs_samples <- setdiff(colnames(hs_dt), gene_col)

# Warn if any sample names overlap
overlap <- intersect(ct_samples, hs_samples)
if (length(overlap) > 0)
  warning(length(overlap), " sample name(s) appear in both CT and HS matrices.")

# ----- Helper: compute mean, median, sd, MAD, CV² from a matrix -----
cond_stats <- function(mat) {
  mn <- rowMeans(mat, na.rm = TRUE)
  sd_ <- apply(mat, 1, sd, na.rm = TRUE)
  list(
    mean   = mn,
    median = apply(mat, 1, median, na.rm = TRUE),
    sd     = sd_,
    mad    = apply(mat, 1, mad,    na.rm = TRUE),
    cv2    = (sd_ / mn)^2
  )
}

# ----- Compute per-condition and full-matrix statistics -----
cat("Computing CT-specific statistics...\n")
ct_mat <- as.matrix(ct_dt[, ..ct_samples])
ct_s   <- cond_stats(ct_mat)
cat("  CT:", nrow(ct_mat), "genes x", ncol(ct_mat), "samples\n")

cat("Computing HS-specific statistics...\n")
hs_mat <- as.matrix(hs_dt[, ..hs_samples])
hs_s   <- cond_stats(hs_mat)
cat("  HS:", nrow(hs_mat), "genes x", ncol(hs_mat), "samples\n")

cat("Computing full-matrix statistics...\n")
full_mat <- cbind(ct_mat, hs_mat)
rownames(full_mat) <- gene_ids
n_samples <- ncol(full_mat)
full_s <- cond_stats(full_mat)
cat("  Full:", nrow(full_mat), "genes x", n_samples, "samples\n\n")

# Ranks on full matrix: rank 1 = highest / most variable
rank_mean   <- rank(-full_s$mean,   ties.method = "average", na.last = "keep")
rank_median <- rank(-full_s$median, ties.method = "average", na.last = "keep")
rank_mad    <- rank(-full_s$mad,    ties.method = "average", na.last = "keep")
rank_cv2    <- rank(-full_s$cv2,    ties.method = "average", na.last = "keep")

results_dt <- data.table(
  gene_id         = gene_ids,
  # CT-specific
  mean_ct         = ct_s$mean,
  median_ct       = ct_s$median,
  sd_ct           = ct_s$sd,
  mad_ct          = ct_s$mad,
  cv2_ct          = ct_s$cv2,
  # HS-specific
  mean_hs         = hs_s$mean,
  median_hs       = hs_s$median,
  sd_hs           = hs_s$sd,
  mad_hs          = hs_s$mad,
  cv2_hs          = hs_s$cv2,
  # Cross-condition differences (HS − CT)
  mad_hs_minus_ct = hs_s$mad - ct_s$mad,
  cv2_hs_minus_ct = hs_s$cv2 - ct_s$cv2,
  # Full CT+HS matrix
  mean_full       = full_s$mean,
  median_full     = full_s$median,
  sd_full         = full_s$sd,
  mad_full        = full_s$mad,
  cv2_full        = full_s$cv2,
  # Ranks (full matrix)
  rank_mean       = rank_mean,
  rank_median     = rank_median,
  rank_mad        = rank_mad,
  rank_cv2        = rank_cv2
)

# ----- Join gene symbols (union of CT + HS annotations) -----
cat("Loading and merging annotation tables...\n")
load_mapping <- function(path, label) {
  dt <- fread(path, data.table = TRUE)
  if (!all(c("ENSEMBL", "SYMBOL") %in% colnames(dt)))
    stop(label, " annotation file must have 'ENSEMBL' and 'SYMBOL' columns.")
  dt[, .(gene_id = ENSEMBL, SYMBOL)]
}
ct_map <- load_mapping(args$ct_mapping_file, "CT")
hs_map <- load_mapping(args$hs_mapping_file, "HS")
# Union: prefer CT SYMBOL when both provide the same gene_id
mapping_dt <- unique(rbind(ct_map, hs_map), by = "gene_id")
cat("  CT annotation entries: ", nrow(ct_map), "\n")
cat("  HS annotation entries: ", nrow(hs_map), "\n")
cat("  Combined unique genes: ", nrow(mapping_dt), "\n\n")

results_dt <- merge(mapping_dt, results_dt, by = "gene_id", all.y = TRUE)

# Sort by full-matrix MAD rank
setorder(results_dt, rank_mad, na.last = TRUE)

# ----- Write xlsx -----
wb <- createWorkbook()
addWorksheet(wb, "Gene_Stats")
writeData(wb, "Gene_Stats", as.data.frame(results_dt),
          headerStyle = createStyle(textDecoration = "bold"))
setColWidths(wb, "Gene_Stats", cols = seq_len(ncol(results_dt)), widths = "auto")
saveWorkbook(wb, args$output_file, overwrite = FALSE)

# ----- Console summary -----
cat("=== Summary ===\n")
cat("Genes in full matrix:  ", nrow(results_dt), "\n")
cat("Total samples (CT+HS): ", n_samples,
    " (", length(ct_samples), " CT + ", length(hs_samples), " HS)\n", sep = "")
cat("Top 5 genes by mean expression (full matrix):\n")
top_mean <- copy(results_dt); setorder(top_mean, rank_mean, na.last = TRUE)
print(head(top_mean[, .(gene_id, SYMBOL, mean_full, rank_mean)], 5))
cat("Top 5 genes by MAD (full matrix):\n")
print(head(results_dt[, .(gene_id, SYMBOL, mad_full, rank_mad)], 5))
cat("Top 5 genes by CV² (full matrix):\n")
top_cv2 <- copy(results_dt); setorder(top_cv2, rank_cv2, na.last = TRUE)
print(head(top_cv2[, .(gene_id, SYMBOL, cv2_full, rank_cv2)], 5))
cat("Top 5 genes by MAD increase (HS > CT):\n")
top_delta <- copy(results_dt)
setorder(top_delta, -mad_hs_minus_ct, na.last = TRUE)
print(head(top_delta[, .(gene_id, SYMBOL, mad_ct, mad_hs, mad_hs_minus_ct)], 5))
cat("\nResults saved to: ", args$output_file, "\n")
