# ==============================================================================
# preprocess_quintile.R
#
# Defines prepare_quintile_groups(): filters to a given condition, loads and
# merges surrogate variables, computes CPM-based expression quintile groups for
# a focal gene, and returns a count matrix + annotated metadata ready for
# run_gamlss_dvgdeg().
#
# Key file format notes:
#   count_file  — row 1 is sample IDs only (no "gene_id" header for col 1);
#                 data rows: gene_id <tab> count...
#   sv_files    — R write.table output; row labels are positional indices
#                 matching count_file column order; no sample ID column
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(tibble)
  library(edgeR)
  library(purrr)
})

# ------------------------------------------------------------------------------
# read_count_matrix
#
# Internal helper: reads the count file, handling the missing gene_id column
# header. Returns a matrix with rownames = gene IDs, colnames = sample IDs.
# Also returns the ordered sample ID vector (used for SV positional alignment).
# ------------------------------------------------------------------------------
.read_count_matrix <- function(count_file) {
  # The header line contains only sample IDs (no label for the gene_id column)
  header_sample_ids <- scan(count_file, what = "", nlines = 1,
                             sep = "\t", quiet = TRUE)
  count_dt <- fread(count_file, skip = 1, header = FALSE,
                    col.names    = c("gene_id", header_sample_ids),
                    colClasses   = c("character", rep("integer", length(header_sample_ids))),
                    showProgress = FALSE)
  count_mat <- as.matrix(
    tibble::column_to_rownames(as.data.frame(count_dt), "gene_id")
  )
  list(count_mat = count_mat, sample_ids = header_sample_ids)
}

# ------------------------------------------------------------------------------
# read_sv_files
#
# Internal helper: reads surrogate variable files (one SV per file) and assigns
# sample IDs by position (row i in file = column i in count_file).
# Returns a data.table with sample_id_col + one column per SV.
# ------------------------------------------------------------------------------
.read_sv_files <- function(sv_files, sample_ids, sample_id_col = "id") {
  n_samp <- length(sample_ids)

  sv_list <- imap(sv_files, function(sv_path, sv_name) {
    # R write.table format: header "x", rows "1" -0.002... (space-separated)
    sv_raw  <- read.table(sv_path, header = TRUE)   # rownames = positional indices
    sv_vals <- sv_raw[[1]]
    if (length(sv_vals) != n_samp) {
      stop(glue::glue(
        "SV file '{sv_path}' has {length(sv_vals)} rows ",
        "but count_file has {n_samp} samples. ",
        "SV row position must match count_file column order."
      ))
    }
    setNames(sv_vals, sample_ids)
  })

  sv_tab <- as.data.table(sv_list)
  sv_tab[[sample_id_col]] <- sample_ids
  setcolorder(sv_tab, c(sample_id_col, names(sv_files)))
  sv_tab
}

# ------------------------------------------------------------------------------
# prepare_quintile_groups
#
# Arguments:
#   count_file     path to the raw counts file (no gene_id header on row 1)
#   meta_file      path to the sample metadata file (tab-separated)
#   sv_files       named character vector: names become column names in meta,
#                  values are file paths; e.g. c(sv1 = "data/sv1.txt", ...)
#   focal_gene     gene ID whose expression defines the quintile split
#   condition_col  metadata column used for filtering (e.g. "treatment")
#   condition_val  value to retain (e.g. "CT")
#   trt_col_name   name for the new quintile column added to metadata
#   frac_low          lower quantile cutoff: samples <= this percentile → label_low
#   frac_high         upper quantile cutoff: samples >= this percentile → label_high
#                  samples strictly between frac_low and frac_high are dropped
#   label_low      factor label for the low-expression group (e.g. "Q1")
#   label_high     factor label for the high-expression group (e.g. "Q5")
#   sample_id_col  column in meta_file that matches colnames(count_mat)
#
# Returns:
#   list(count_mat, meta_dt) — both aligned to the same sample order,
#   with meta_dt containing trt_col_name as an ordered factor (label_low first)
# ------------------------------------------------------------------------------
prepare_quintile_groups <- function(
  count_file,
  meta_file,
  sv_files,
  focal_gene,
  condition_col,
  condition_val,
  trt_col_name,
  frac_low        = 0.20,
  frac_high       = 0.80,
  label_low    = "Q1",
  label_high   = "Q5",
  sample_id_col = "id"
) {
  # --- 1. read count matrix and get the sample order for SV alignment ---
  message("Reading count matrix...")
  count_info  <- .read_count_matrix(count_file)
  count_mat   <- count_info$count_mat
  sample_ids  <- count_info$sample_ids   # ordered as in count file columns

  # --- 2. load surrogate variables and assign sample IDs by position ---
  message("Loading surrogate variables...")
  sv_tab <- .read_sv_files(sv_files, sample_ids, sample_id_col)

  # --- 3. read metadata, merge SVs, filter to condition ---
  message("Loading metadata and merging SVs...")
  meta_dt <- fread(meta_file)
  meta_dt <- merge(meta_dt, sv_tab, by = sample_id_col, all.x = FALSE)
  meta_dt <- meta_dt[as.character(get(condition_col)) == as.character(condition_val)]

  if (nrow(meta_dt) == 0) {
    stop(glue::glue(
      "No samples remain after filtering '{condition_col}' == '{condition_val}'"
    ))
  }

  # --- 4. subset count matrix to kept samples ---
  kept_samples <- meta_dt[[sample_id_col]]
  missing      <- setdiff(kept_samples, colnames(count_mat))
  if (length(missing) > 0) {
    stop(glue::glue(
      "{length(missing)} sample(s) in metadata not found in count matrix: ",
      paste(head(missing, 5), collapse = ", ")
    ))
  }
  count_subset <- count_mat[, kept_samples, drop = FALSE]

  # check focal gene is present
  if (!focal_gene %in% rownames(count_subset)) {
    stop(glue::glue("focal_gene '{focal_gene}' not found in count matrix"))
  }

  # --- 5. compute CPM for focal gene and assign quintile groups ---
  # CPM normalises for library size so samples are comparable
  focal_cpm <- cpm(DGEList(count_subset))[focal_gene, ]

  low_cut  <- quantile(focal_cpm, frac_low)
  high_cut <- quantile(focal_cpm, frac_high)

  quintile_group <- dplyr::case_when(
    focal_cpm <= low_cut  ~ label_low,
    focal_cpm >= high_cut ~ label_high,
    TRUE                  ~ NA_character_
  )
  names(quintile_group) <- names(focal_cpm)

  n_low  <- sum(quintile_group == label_low,  na.rm = TRUE)
  n_high <- sum(quintile_group == label_high, na.rm = TRUE)
  n_drop <- sum(is.na(quintile_group))
  message(glue::glue(
    "Quintile assignment for {focal_gene} (condition: {condition_val}): ",
    "{label_low} n={n_low}, {label_high} n={n_high}, dropped n={n_drop}"
  ))

  if (n_low == 0 || n_high == 0) {
    stop("One or both quintile groups are empty — check frac_low/frac_high cutoffs")
  }

  # --- 6. add quintile column to metadata; drop mid-range samples ---
  meta_dt[[trt_col_name]] <- quintile_group[meta_dt[[sample_id_col]]]
  meta_dt <- meta_dt[!is.na(get(trt_col_name))]
  meta_dt[[trt_col_name]] <- factor(meta_dt[[trt_col_name]],
    levels = c(label_low, label_high))

  # --- 7. align count matrix to the final sample order in meta ---
  final_samples <- meta_dt[[sample_id_col]]
  count_final   <- count_subset[, final_samples, drop = FALSE]

  list(count_mat = count_final, meta_dt = meta_dt)
}
