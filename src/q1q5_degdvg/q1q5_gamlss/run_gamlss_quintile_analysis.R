#!/usr/bin/env Rscript
# ==============================================================================
# run_gamlss_quintile_analysis.R
#
# Driver for DEG/DVG analysis comparing low-expression (Q1) vs high-expression
# (Q5) samples of a focal gene within a given condition.
#
# CLI uses argparse named flags (run with --help for full details):
#
# USAGE — interactive: edit PARAMETERS block below, then:
#   Rscript src/gamlss_q1q5/run_gamlss_quintile_analysis.R
#   (runs with defaults: FBgn0027560, 1 test gene, CT condition)
#
# USAGE — named flags (all optional; override the defaults below):
#   Rscript src/gamlss_q1q5/run_gamlss_quintile_analysis.R \
#     --focal-gene FBgn0027560 \
#     --condition-val 1 \
#     --full-run
#
#   --focal-gene GENE_ID   FlyBase gene ID for the Q1/Q5 split   (default: FBgn0027560)
#   --condition-val INT    1 = CT, 2 = HS                        (default: 1)
#   --full-run             analyse all ~8786 genes; omit for test mode
#   --n-test-genes N       limit to first N genes when not using --full-run (default: 1)
#
# Output files written to gamlss_dir:
#   {focal_gene}_{CT|HS}_Q1vsQ5_gamlss_dvgdeg.csv
#   {focal_gene}_{CT|HS}_Q1vsQ5_gamlss_dvgdeg_sum.csv
# ==============================================================================

# ==============================================================================
# PARAMETERS — defaults used when no CLI arguments are given
# ==============================================================================
focal_gene    <- "FBgn0027560"   # gene whose expression defines Q1/Q5 split
n_test_genes  <- 1L              # NULL for a full run; integer for test mode
condition_val <- 1L              # 1 = CT (control), 2 = HS (high sugar)

# Input files (paths relative to code/dronetanalysis/)
count_file <- "../../data/count/RawCounts_CPM1_head_hsctrl_onlyGEMMAsamples_Mar21.21.txt"
meta_file  <- "../../data/count/Info_RawCounts_CPM1_head_hsctrl_Jul9.20.txt"
sv_files   <- c(
  sv1 = "../../data/count/TMM_Voom_sv1_10.txt",
  sv2 = "../../data/count/TMM_Voom_sv2_10.txt",
  sv3 = "../../data/count/TMM_Voom_sv3_10.txt",
  sv4 = "../../data/count/TMM_Voom_sv4_10.txt"
)

# Sample grouping
condition_col <- "treatment"     # metadata column used to filter samples
trt_col_name  <- "expr_quintile" # name for the quintile column added to metadata
ref_group     <- "Q1"            # lower expression group (reference level)
frac_low      <- 0.20            # bottom fraction → Q1
frac_high     <- 0.80            # top fraction → Q5
label_low     <- "Q1"
label_high    <- "Q5"

# GAMLSS model formulas.
# fit_one wraps ALL formulas as ~ model.matrix(formula, meta), so GAMLSS always
# adds its own intercept on top. ~0 + causes collinearity and drops Q5 from the
# coefficient table. Use standard contrast coding (~expr_quintile, not ~0+).
formula_mu <- paste(
  "~expr_quintile +",
  "RNAlibBatch + RNAseqBatch + egglayBatch + platingBatch + well +",
  "sv1 + sv2 + sv3 + sv4"
)
formula_sigma  <- "~expr_quintile"
formula_nontrt <- paste(
  "~RNAlibBatch + RNAseqBatch + egglayBatch + platingBatch + well +",
  "sv1 + sv2 + sv3 + sv4"
)

# Sum-to-zero coding for categorical batch covariates only.
# SVs (sv1–sv4) are continuous — contr.sum applies only to factors.
# Set to NULL to use default treatment contrasts.
mu_sum_covars <- c("RNAlibBatch", "RNAseqBatch", "egglayBatch", "platingBatch", "well")

gamlss_dir    <- "../../results/gamlss"
n_workers     <- 1   # bplapply serializes full model objects back across IPC; use SerialParam to avoid OOM
fdr_threshold <- 0.05
# ==============================================================================

# --- CLI arguments (named flags; defaults come from the PARAMETERS block above) ---
library(argparse)
parser <- ArgumentParser(
  description = "GAMLSS DEG/DVG Q1-vs-Q5 analysis for a single focal gene"
)
parser$add_argument(
  "--focal-gene",
  default = focal_gene, metavar = "GENE_ID",
  help    = "FlyBase gene ID whose expression defines the Q1/Q5 split [default: %(default)s]"
)
parser$add_argument(
  "--n-test-genes",
  type = "integer", default = as.integer(n_test_genes), metavar = "N",
  help = "limit to first N genes for a quick test; ignored when --full-run is set [default: %(default)s]"
)
parser$add_argument(
  "--full-run",
  action = "store_true", default = FALSE,
  help   = "analyse all genes (sets n_test_genes = NULL, overrides --n-test-genes)"
)
parser$add_argument(
  "--condition-val",
  type = "integer", default = as.integer(condition_val), metavar = "INT",
  help = "condition code: 1 = CT, 2 = HS [default: %(default)s]"
)
parser$add_argument(
  "--out-dir",
  default = gamlss_dir, metavar = "DIR",
  help    = "output directory for result CSVs [default: %(default)s]"
)
.args         <- parser$parse_args()
focal_gene    <- .args$focal_gene
n_test_genes  <- if (.args$full_run) NULL else .args$n_test_genes
condition_val <- as.integer(.args$condition_val)
gamlss_dir    <- .args$out_dir

# --- derived values (always computed after CLI overrides) ---
cond_tag   <- c("1" = "CT", "2" = "HS")[[as.character(condition_val)]]
out_prefix <- glue::glue("{focal_gene}_{cond_tag}_{label_low}vs{label_high}")
# ==============================================================================

library(glue)
source("src/q1q5_degdvg/shared/preprocess_quintile.R")
source("src/q1q5_degdvg/q1q5_gamlss/tools_gamlss.R")
source("src/q1q5_degdvg/q1q5_gamlss/gamlss_dvgdeg_core.R")

if (!dir.exists(gamlss_dir)) dir.create(gamlss_dir, recursive = TRUE)

# --- Step 1: preprocess ---
message(glue("Preparing Q1/Q5 groups for focal gene: {focal_gene} ({cond_tag})"))
inputs <- prepare_quintile_groups(
  count_file    = count_file,
  meta_file     = meta_file,
  sv_files      = sv_files,
  focal_gene    = focal_gene,
  condition_col = condition_col,
  condition_val = condition_val,
  trt_col_name  = trt_col_name,
  frac_low      = frac_low,
  frac_high     = frac_high,
  label_low     = label_low,
  label_high    = label_high
)

if (!is.null(n_test_genes)) {
  inputs$count_mat <- inputs$count_mat[seq_len(min(n_test_genes, nrow(inputs$count_mat))), ]
  message(glue("TEST MODE: limiting to {nrow(inputs$count_mat)} genes"))
}

message(glue("Samples: n={ncol(inputs$count_mat)}, genes={nrow(inputs$count_mat)}"))
message("Group counts:")
print(table(inputs$meta_dt[[trt_col_name]]))

# --- Step 2: run GAMLSS DEG/DVG analysis ---
message("Fitting GAMLSS models...")
dvdg_tab <- run_gamlss_dvgdeg(
  count_mat      = inputs$count_mat,
  meta_dt        = inputs$meta_dt,
  trt_col        = trt_col_name,
  ref_group      = ref_group,
  formula_mu     = formula_mu,
  formula_sigma  = formula_sigma,
  formula_nontrt = formula_nontrt,
  mu_sum_covars  = mu_sum_covars,
  fdr_threshold  = fdr_threshold,
  n_workers      = n_workers,
  out_prefix     = out_prefix,
  gamlss_dir     = gamlss_dir
)

message(glue("Done. Results written to: {gamlss_dir}/{out_prefix}_gamlss_dvgdeg.csv"))
