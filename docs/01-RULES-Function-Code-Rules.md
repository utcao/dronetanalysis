# Function Code Rules

Rules for writing R functions in this project. Derived from code review of
`15analysis/plot_pca_l2l1_variability.R` and applicable project-wide.

---

## Variable Naming

**Use generic names that survive condition renaming. Use `ctrl`/`trt` to make
the reference vs treatment role explicit.**

Condition-specific prefixes (`ct_`, `hs_`) bleed into helper variables and
make functions brittle when conditions change. Use `ctrl` for the reference
condition and `trt` for the treatment/comparison condition. Pass the actual
suffixes as arguments so the function works for any pair.

| Avoid | Prefer |
|---|---|
| `ct_vals`, `hs_vals` | `ctrl_vals`, `trt_vals` |
| `ct_cols`, `hs_cols` | `ctrl_cols`, `trt_cols` |
| `cond1_suffix`, `cond2_suffix` (params) | `ctrl_suffix`, `trt_suffix` |
| `ct_r` (correlation result) | `cor_result`, `spearman_result` |
| `rho_lbl` | `cor_label` |
| `df_p` | `plot_df` (ggplot idiom) |
| `lbl_vec` | `gene_labels` (drop redundant `_vec`) |
| `sym_up`, `map_up` | `symbol_uc`, `map_uc` (suffix `_uc` = uppercase) |
| `ok` | `complete`, `both_observed` |

**Never shadow base R names.** `base`, `data`, `df`, `file`, `list`, `c`,
`table` are reserved. Use `feat_base`, `feature_names`, `result_dt`, etc.

---

## Scoping and Loop Structure

Separate **pre-loop setup** from **per-iteration work** with a blank line and
a comment. Variables computed once (gene annotation, label vectors) must not
sit unannounced next to per-iteration variables; the asymmetry confuses readers
who need to know what varies each iteration.

```r
# --- annotation: computed once, shared across all feature pairs ---
is_ann     <- symbol_uc %in% map_uc
gene_labels <- ifelse(is_ann, gene_map$label[match(symbol_uc, map_uc)], "")

for (i in seq_along(feature_names)) {
  # --- per-iteration: specific to this feature pair ---
  ctrl_vals <- data[[ctrl_cols[i]]]
  trt_vals  <- data[[trt_cols[i]]]
  ...
}
```

---

## Early Returns

Put each statement on its own line. A `;`-separated compound line is easy to
miss when scanning.

```r
# bad
if (!length(cols_c1)) { cat("no pairs\n"); return(invisible(NULL)) }

# good
if (!length(cols_c1)) {
  cat("  no matched column pairs found.\n")
  return(invisible(NULL))
}
```

---

## Logical Index Vectors

Prefer explicit assignment over the indirect `seq_along() %in% idx` pattern.

```r
# less clear
is_top <- seq_along(ctrl_vals) %in% top_idx

# clearer
is_top <- logical(length(ctrl_vals))
is_top[top_idx] <- TRUE
```

---

## Development Comments (Alternatives)

While a function is in active development, **keep alternative implementations
as commented-out code** so they can be swapped in without re-deriving them.
Mark them with `# ALT:` so they are easy to find and easy to remove at
stabilisation.

```r
# score by mean (default)
score <- (ctrl_vals + trt_vals) / 2
# ALT: distance from origin
# score <- sqrt(ctrl_vals^2 + trt_vals^2)
# ALT: max of the two conditions
# score <- pmax(ctrl_vals, trt_vals, na.rm = TRUE)
```

Remove `# ALT:` blocks before marking a function as stable.

---

## When to Add a Comment

Only when the **why** is non-obvious: a hidden constraint, a subtle invariant,
a workaround for a specific behaviour. Do not comment what the code already
says.

```r
# bad — obvious from the code
score <- (ctrl_vals + trt_vals) / 2  # compute the mean of the two conditions

# good — non-obvious choice
score <- (ctrl_vals + trt_vals) / 2  # arithmetic mean; Spearman rank is invariant to monotone transforms
```

---

## Intermediate Filter Variables

Avoid creating a named intermediate that is immediately subsetted away — it
suggests the intermediate is reused when it is not.

```r
# confusing: trt_match looks reusable but is discarded after one line
trt_match  <- sub("_ctrl$", "_trt", ctrl_cols)
keep       <- trt_match %in% all_cols
ctrl_cols  <- ctrl_cols[keep]
trt_cols   <- trt_match[keep]

# cleaner: keep both filtered vectors together, no orphan intermediate
trt_cols_cand <- sub("_ctrl$", "_trt", ctrl_cols)
pair_exists   <- trt_cols_cand %in% all_cols
ctrl_cols     <- ctrl_cols[pair_exists]
trt_cols      <- trt_cols_cand[pair_exists]
```

The rename from `keep` → `pair_exists` also makes the boolean's meaning
self-documenting.
