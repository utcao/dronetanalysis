#!/usr/bin/env Rscript
# ==============================================================================
# Wrapper: Extract Top-N Genes by delta_mean_mad → Permutation Test
#
# Reads the MAD variability summary produced by
# compute_mad_transcriptomic_variability.R, ranks focus genes by delta_mean_mad,
# extracts the top N, and passes them as --focus-genes to
# compute_mad_permute_transcriptomic_noise.R.
#
# This replaces ad-hoc shell pipelines with a single, self-documented, logged
# command.  Use --dry-run to print the constructed command without running it.
# The selected gene list is always written to --gene-list-out for traceability.
#
# Usage (CT, top 200 by absolute delta, FDR pre-filter):
#   Rscript run_mad_permute_top_genes.R \
#     --mad-summary  results/variability/voomct_all_genes_mad_summary.xlsx \
#     --expr-file    data/processed/VOOM/voomdataCtrl.txt \
#     --mapping-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
#     --output-file  results/variability/mad_permutation_ct_top200.xlsx \
#     --top-n        200 \
#     --direction    abs \
#     --min-fdr      0.05 \
#     --n-perms      500 \
#     --condition-label "Control"
#
# Usage (HS, top 100 positive delta only, dry-run to inspect command):
#   Rscript run_mad_permute_top_genes.R \
#     --mad-summary  results/variability/voomhs_all_genes_mad_summary.xlsx \
#     --expr-file    data/processed/VOOM/voomdataHs.txt \
#     --mapping-file results/result_voomhs/rewiring_hubs_hs_anno_0408_2026.tsv \
#     --output-file  results/variability/mad_permutation_hs_top100.xlsx \
#     --top-n        100 \
#     --direction    positive \
#     --n-perms      500 \
#     --condition-label "HeatStress" \
#     --dry-run
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(argparse)
  library(openxlsx)
})

# ----- Resolve this script's own directory for the default permute-script path -----
get_script_dir <- function() {
  args      <- commandArgs(trailingOnly = FALSE)
  flag      <- "--file="
  path_flag <- args[grep(flag, args)]
  if (length(path_flag) == 0) return(getwd())
  dirname(normalizePath(sub(flag, "", path_flag)))
}

# ----- CLI arguments -----
parser <- ArgumentParser(
  description = paste(
    "Extract top-N focus genes by delta_mean_mad and run the permutation test.",
    "All --n-perms / --seed / --low-frac / --high-frac / --null-dist-dir args",
    "are forwarded to compute_mad_permute_transcriptomic_noise.R."
  )
)

# --- Gene selection ---
parser$add_argument("--mad-summary",
  required = TRUE,
  help = paste(
    "xlsx output of compute_mad_transcriptomic_variability.R",
    "(sheet 'MAD_variability'). Genes are ranked by delta_mean_mad."
  )
)
parser$add_argument("--top-n",
  type = "integer", default = 200,
  help = "Number of top genes to extract from the MAD summary (default: 200)"
)
parser$add_argument("--direction",
  default = "abs",
  help = paste(
    "Which delta_mean_mad direction to rank by.",
    "'abs'      = largest |delta| regardless of sign (default),",
    "'positive' = largest positive delta (HIGH more variable than LOW),",
    "'negative' = most negative delta (LOW more variable than HIGH)"
  )
)
parser$add_argument("--min-fdr",
  type = "double", default = NULL,
  help = paste(
    "Optional: pre-filter genes by wilcoxon_p_adj < min-fdr before ranking.",
    "Set to 0.05 to restrict to originally significant genes."
  )
)
parser$add_argument("--gene-list-out",
  default = NULL,
  help = paste(
    "Optional: path to write the selected gene list as a TSV (gene_id, SYMBOL,",
    "delta_mean_mad, wilcoxon_p_adj). Defaults to <output-file>_gene_list.tsv."
  )
)

# --- Pass-through: permutation script location ---
parser$add_argument("--permute-script",
  default = NULL,
  help = paste(
    "Path to compute_mad_permute_transcriptomic_noise.R.",
    "Defaults to the script in the same directory as this wrapper."
  )
)

# --- Pass-through: required permutation args ---
parser$add_argument("--expr-file",
  required = TRUE,
  help = "VOOM expression matrix (forwarded to permutation script)"
)
parser$add_argument("--mapping-file",
  required = TRUE,
  help = "TSV/CSV with gene_id and SYMBOL columns (forwarded)"
)
parser$add_argument("--output-file",
  required = TRUE,
  help = "Output .xlsx path for the permutation results (forwarded)"
)

# --- Pass-through: optional permutation args ---
parser$add_argument("--null-dist-dir",
  default = NULL,
  help = "Directory for per-gene null distribution TSVs (forwarded, optional)"
)
parser$add_argument("--low-frac",
  type = "double", default = 0.2,
  help = "LOW group fraction (forwarded, default: 0.2)"
)
parser$add_argument("--high-frac",
  type = "double", default = 0.2,
  help = "HIGH group fraction (forwarded, default: 0.2)"
)
parser$add_argument("--n-perms",
  type = "integer", default = 500,
  help = "Permutations per gene (forwarded, default: 500)"
)
parser$add_argument("--seed",
  type = "integer", default = 42,
  help = "Random seed (forwarded, default: 42)"
)
parser$add_argument("--condition-label",
  default = "Condition",
  help = "Condition label (forwarded)"
)

# --- Control ---
parser$add_argument("--dry-run",
  action = "store_true",
  help = "Print the constructed Rscript command without executing it"
)

args <- parser$parse_args()

# ----- Resolve permutation script path -----
if (is.null(args$permute_script)) {
  args$permute_script <- file.path(
    get_script_dir(),
    "compute_mad_permute_transcriptomic_noise.R"
  )
}

cat("=== run_mad_permute_top_genes ===\n")
cat("MAD summary:      ", args$mad_summary,     "\n")
cat("Top N:            ", args$top_n,            "\n")
cat("Direction:        ", args$direction,        "\n")
cat("Min FDR filter:   ",
    if (!is.null(args$min_fdr)) args$min_fdr else "(none)", "\n")
cat("Permute script:   ", args$permute_script,   "\n")
cat("Output file:      ", args$output_file,      "\n")
cat("Dry-run:          ", args$dry_run,          "\n\n")

# ----- Validate -----
if (!file.exists(args$mad_summary))
  stop("MAD summary not found: ", args$mad_summary)
if (!file.exists(args$permute_script))
  stop("Permutation script not found: ", args$permute_script)
if (!file.exists(args$expr_file))
  stop("Expression file not found: ", args$expr_file)
if (!file.exists(args$mapping_file))
  stop("Mapping file not found: ", args$mapping_file)
if (!args$direction %in% c("abs", "positive", "negative"))
  stop("--direction must be one of: abs, positive, negative")

# ----- Load MAD summary -----
cat("Loading MAD summary...\n")
mad_dt <- as.data.table(read.xlsx(args$mad_summary, sheet = "MAD_variability"))

required_cols <- c("gene_id", "delta_mean_mad", "wilcoxon_p_adj")
missing <- setdiff(required_cols, colnames(mad_dt))
if (length(missing) > 0)
  stop("MAD summary is missing columns: ", paste(missing, collapse = ", "))

cat("  Genes in summary: ", nrow(mad_dt), "\n")

# ----- Pre-filter by FDR -----
if (!is.null(args$min_fdr)) {
  before <- nrow(mad_dt)
  mad_dt <- mad_dt[!is.na(wilcoxon_p_adj) & wilcoxon_p_adj < args$min_fdr]
  cat("  After FDR <", args$min_fdr, "filter:", nrow(mad_dt),
      "(removed", before - nrow(mad_dt), ")\n")
  if (nrow(mad_dt) == 0)
    stop("No genes remain after FDR filter. Lower --min-fdr or remove it.")
}

# ----- Rank by delta_mean_mad -----
mad_dt <- mad_dt[!is.na(delta_mean_mad)]

if (args$direction == "abs") {
  mad_dt[, sort_key := abs(delta_mean_mad)]
  setorder(mad_dt, -sort_key)
  cat("  Ranking by |delta_mean_mad| (descending)\n")
} else if (args$direction == "positive") {
  mad_dt <- mad_dt[delta_mean_mad > 0]
  setorder(mad_dt, -delta_mean_mad)
  cat("  Ranking by delta_mean_mad > 0 (descending, HIGH more variable)\n")
} else {
  mad_dt <- mad_dt[delta_mean_mad < 0]
  setorder(mad_dt, delta_mean_mad)
  cat("  Ranking by delta_mean_mad < 0 (ascending, LOW more variable)\n")
}

if (nrow(mad_dt) == 0)
  stop("No genes with delta_mean_mad in the requested direction.")

# ----- Extract top N -----
n_avail  <- nrow(mad_dt)
n_select <- min(args$top_n, n_avail)
if (n_select < args$top_n)
  warning(sprintf(
    "--top-n %d requested but only %d genes available; using %d.",
    args$top_n, n_avail, n_select
  ))

selected <- mad_dt[seq_len(n_select)]
cat("  Selected top", n_select, "genes\n")
cat("  delta_mean_mad range: [",
    round(min(selected$delta_mean_mad), 4), ",",
    round(max(selected$delta_mean_mad), 4), "]\n\n")

# ----- Write gene list (always, for traceability) -----
gene_list_path <- args$gene_list_out
if (is.null(gene_list_path)) {
  gene_list_path <- sub("\\.xlsx$", "_gene_list.tsv", args$output_file)
}
dir.create(dirname(gene_list_path), recursive = TRUE, showWarnings = FALSE)

out_cols <- intersect(
  c("gene_id", "SYMBOL", "delta_mean_mad", "wilcoxon_p_adj"),
  colnames(selected)
)
fwrite(selected[, ..out_cols], gene_list_path, sep = "\t")
cat("Gene list written to:", gene_list_path, "\n\n")

# ----- Build focus-genes string -----
focus_genes_str <- paste(selected$gene_id, collapse = ",")

# ----- Construct Rscript call -----
script_args <- c(
  args$permute_script,
  "--expr-file",       args$expr_file,
  "--mapping-file",    args$mapping_file,
  "--output-file",     args$output_file,
  "--focus-genes",     focus_genes_str,
  "--low-frac",        args$low_frac,
  "--high-frac",       args$high_frac,
  "--n-perms",         args$n_perms,
  "--seed",            args$seed,
  "--condition-label", args$condition_label
)

if (!is.null(args$null_dist_dir)) {
  script_args <- c(script_args, "--null-dist-dir", args$null_dist_dir)
  dir.create(args$null_dist_dir, recursive = TRUE, showWarnings = FALSE)
}

dir.create(dirname(args$output_file), recursive = TRUE, showWarnings = FALSE)

# ----- Dry-run or execute -----
rscript_bin <- file.path(R.home("bin"), "Rscript")

cat("Command to run:\n")
cat(" ", rscript_bin, "\\\n")
for (i in seq(1, length(script_args), by = 2)) {
  if (i + 1 <= length(script_args)) {
    cat(sprintf("    %-20s %s \\\n", script_args[i], script_args[i + 1]))
  } else {
    cat("   ", script_args[i], "\n")
  }
}
cat("\n")

if (args$dry_run) {
  cat("Dry-run mode: command printed above, not executed.\n")
} else {
  cat("Executing...\n\n")
  exit_code <- system2(rscript_bin, args = script_args)
  if (exit_code != 0)
    stop("compute_mad_permute_transcriptomic_noise.R exited with code ", exit_code)
  cat("\nDone. Results saved to:", args$output_file, "\n")
}
