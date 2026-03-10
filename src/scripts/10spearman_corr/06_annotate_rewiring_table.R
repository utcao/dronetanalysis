#!/usr/bin/env Rscript
# ==============================================================================
# Annotate Rewiring Table with Gene Information
#
# Adds gene symbol, name, Entrez ID, and Ensembl ID columns to a rewiring
# (focus gene metrics) TSV file using the org.Dm.eg.db Bioconductor database.
#
# Usage:
#   Rscript 06_annotate_rewiring_table.R \
#       --input-tsv results/focus_gene_metrics.tsv \
#       --output-tsv results/focus_gene_metrics_annotated.tsv
# ==============================================================================

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(org.Dm.eg.db)
})

# Source annotation utility
# Source annotation utility
# This works for Rscript and interactive sessions
get_script_path <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- "--file="
    match <- grep(file_arg, cmd_args)
    if (length(match) > 0) {
        return(normalizePath(sub(file_arg, "", cmd_args[match])))
    } else {
        return(normalizePath(sys.frames()[[1]]$ofile))
    }
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
source(file.path(script_dir, "../../utils/utils_annotation.R"))

# ----- Command-line arguments -----
parser <- ArgumentParser(description = "Annotate rewiring table with gene names and IDs")
parser$add_argument("--input-tsv", required = TRUE,
                    help = "Path to input rewiring/focus gene metrics TSV")
parser$add_argument("--output-tsv", default = NULL,
                    help = "Path to output annotated TSV (default: input with _annotated suffix)")
parser$add_argument("--gene-id-col", default = "gene_id",
                    help = "Name of the column containing FBgn IDs (default: gene_id)")
args <- parser$parse_args()

# Default output path
if (is.null(args$output_tsv)) {
    args$output_tsv <- sub("\\.tsv$", "_annotated.tsv", args$input_tsv)
}

cat("=== Annotate Rewiring Table ===\n")
cat("Input:", args$input_tsv, "\n")
cat("Output:", args$output_tsv, "\n")
cat("Gene ID column:", args$gene_id_col, "\n\n")

# ----- Load data -----
dt <- fread(args$input_tsv)
cat("Loaded", nrow(dt), "rows,", ncol(dt), "columns\n")

if (!args$gene_id_col %in% colnames(dt)) {
    stop("Column '", args$gene_id_col, "' not found. Available: ",
         paste(colnames(dt), collapse = ", "))
}

fbgn_ids <- dt[[args$gene_id_col]]
cat("Unique gene IDs:", length(unique(fbgn_ids)), "\n")

# ----- Annotate -----
cat("Querying org.Dm.eg.db...\n")
annot <- annotate_fbgn(fbgn_ids)
cat("  Mapped:", sum(!is.na(annot$SYMBOL)), "/", nrow(annot), "genes to symbols\n")
cat("  Mapped:", sum(!is.na(annot$GENENAME)), "/", nrow(annot), "genes to names\n")
cat("  Mapped:", sum(!is.na(annot$ENTREZID)), "/", nrow(annot), "genes to Entrez IDs\n")
cat("  Mapped:", sum(!is.na(annot$ENSEMBL)), "/", nrow(annot), "genes to Ensembl IDs\n")

# ----- Merge -----
annot_dt <- as.data.table(annot)
setnames(annot_dt, "FLYBASE", args$gene_id_col)

# Left join to preserve all rows
dt_annotated <- merge(dt, annot_dt, by = args$gene_id_col, all.x = TRUE, sort = FALSE)

# Reorder columns: put annotation columns right after gene_id
id_col_pos <- which(colnames(dt_annotated) == args$gene_id_col)
annot_cols <- c("SYMBOL", "GENENAME", "ENTREZID", "ENSEMBL")
other_cols <- setdiff(colnames(dt_annotated), c(args$gene_id_col, annot_cols))

# Preserve original column order for non-annotation columns
orig_order <- colnames(dt)
other_cols_ordered <- orig_order[orig_order %in% other_cols]

setcolorder(dt_annotated, c(
    orig_order[1:id_col_pos],  # columns up to and including gene_id
    annot_cols,                 # annotation columns
    other_cols_ordered[!other_cols_ordered %in% orig_order[1:id_col_pos]]  # remaining original columns
))

# ----- Write output -----
dir.create(dirname(args$output_tsv), showWarnings = FALSE, recursive = TRUE)
fwrite(dt_annotated, args$output_tsv, sep = "\t")
cat("\nWrote", nrow(dt_annotated), "rows to", args$output_tsv, "\n")

# ----- Show sample -----
cat("\nFirst 5 rows (annotation columns):\n")
print(head(dt_annotated[, c(args$gene_id_col, annot_cols), with = FALSE], 5))
