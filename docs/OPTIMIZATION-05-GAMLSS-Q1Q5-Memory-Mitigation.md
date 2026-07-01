# GAMLSS Q1/Q5 Memory Mitigation Guide

## Problem

HS jobs (larger sample size than CT — n=416 vs n=376 for the `FBgn0001233` example) in
`compute_degdvg_gamlss_q1q5.snakemake` intermittently die with:

```
Error in fwrite(dvdg_tab, ...) :
  Unable to allocate 8 MB * 1 thread buffers; '12: Cannot allocate memory'.
```

13 of 45 jobs failed this way in the 2026-06-23 run (`mem_gb=25` per job). Root cause and the
already-tried `rm()+gc()` fix are analyzed in
[FIX-04-GAMLSS-Q1Q5-CT-HS-Errors.md](FIX-04-GAMLSS-Q1Q5-CT-HS-Errors.md) — in short, three
full-size lists of ~8,786 GAMLSS model objects are held concurrently in memory during fitting,
and the existing cleanup only runs after peak memory has already been reached.

**This pipeline is known to run successfully on other, smaller datasets** — the fitting logic
itself is not broken, this dataset's larger sample counts just push some jobs over the
resource ceiling. Following that signal, the options below are ordered from **cheapest and
lowest-risk to most invasive**. Try them in order and stop as soon as jobs complete reliably —
don't reach for a bigger code change than the situation actually requires.

---

## Option 1 — Increase `mem_gb` (try this first)

**Cost:** zero code risk, one-line config change.
**When it's enough:** if the peak is only modestly above 25 GB, this alone resolves it.

Current setting (`src/q1q5_degdvg/q1q5_gamlss/compute_degdvg_gamlss_q1q5.snakemake:85-87`):

```python
resources:
    mem_gb   = 25,
    walltime = "240:0:0"
```

Bump this for all jobs (simplest):

```python
resources:
    mem_gb   = 45,
    walltime = "240:0:0"
```

Or, if only HS needs more (larger n), key it off the `cond` wildcard so CT jobs don't request
memory they don't need:

```python
resources:
    mem_gb   = lambda wildcards: 45 if wildcards.cond == "HS" else 25,
    walltime = "240:0:0"
```

No changes needed in `profiles/sge/config.yaml` — it already reads `{resources.mem_gb}` from
the rule.

---

## Option 2 — Tune `fwrite()` directly

**Cost:** one-line change in `gamlss_dvgdeg_core.R`.
**When to reach for this:** if Option 1 isn't available (cluster memory ceiling) or the crash
still occurs — the error message itself suggests this.

The failure is inside `data.table::fwrite`'s own buffer allocation (`gamlss_dvgdeg_core.R:307-309`
and `:321-323`). Passing an explicit single-thread, smaller buffer avoids the extra allocations
`fwrite` otherwise makes for multi-threaded writing:

```r
fwrite(dvdg_tab, file = ..., nThread = 1, buffMB = 4)
```

This doesn't reduce the peak memory used during fitting — it only lowers the amount `fwrite`
itself needs at the very end, which is useful if the process is already right at the edge of
its ceiling by that point.

---

## Option 3 — Free `m_sigma_int` / `m_mu_nontrt` immediately after use

**Cost:** small, localized change to `run_gamlss_dvgdeg()` in `gamlss_dvgdeg_core.R`.
**When to reach for this:** if Options 1-2 aren't sufficient or resources aren't available.

Currently (`gamlss_dvgdeg_core.R:133-159`) all three model-family lists are built up front and
kept alive until the final `rm()` at line 298 — even though `m_sigma_int` and `m_mu_nontrt`
are only read later for two scalars per gene each (`logLik()` and `df.fit`, used in the
`map2_dbl` LR-test calls at lines 220-226 and 246-252).

Instead of holding all three lists concurrently, extract the LR-test components for each
non-`m_full` family right after fitting it, and free the model list before fitting the next
one:

```r
m_full     <- bp_apply(seq_len(n_genes), fit_one, ...)
m_full_sum <- bp_apply(m_full, fit_sum_quiet)

m_sigma_int   <- bp_apply(seq_len(n_genes), fit_one, ..., sigma_formula = "~ 1")
sigma_int_ll  <- vapply(m_sigma_int, function(m) if (length(m)) logLik(m)[[1]] else NA_real_, numeric(1))
sigma_int_df  <- vapply(m_sigma_int, function(m) if (length(m)) m$df.fit else NA_real_, numeric(1))
rm(m_sigma_int); gc()

m_mu_nontrt  <- bp_apply(seq_len(n_genes), fit_one, ..., mu_formula = formula_nontrt)
mu_nontrt_ll <- vapply(m_mu_nontrt, function(m) if (length(m)) logLik(m)[[1]] else NA_real_, numeric(1))
mu_nontrt_df <- vapply(m_mu_nontrt, function(m) if (length(m)) m$df.fit else NA_real_, numeric(1))
rm(m_mu_nontrt); gc()
```

Then replace the `map2_dbl(m_full, m_sigma_int, ...)` / `map2_dbl(m_full, m_mu_nontrt, ...)`
LR-test calls inside the `foreach` block with vectorized arithmetic over the pre-extracted
`sigma_int_ll`/`sigma_int_df`/`mu_nontrt_ll`/`mu_nontrt_df` vectors and `m_full`'s own
`logLik`/`df.fit`. This bounds concurrent model-list memory to roughly 1.3× a single family
instead of 3×, since only `m_full` needs to stay alive for the full extraction phase (it feeds
`cpm_count_*`, `mu_*`, `BCV_*`, `sigma_*` in addition to the LR tests).

---

## Option 4 — Gene-batching (last resort)

**Cost:** largest change — restructures the fit/extract loop.
**When to reach for this:** only if cluster memory ceilings make Option 1 infeasible and
Options 2-3 still aren't enough.

Mirrors the batch-mode pattern already used elsewhere in this repo for the same problem shape
(see [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md), which processes genes in chunks
of 100-500 for the bootstrap pipeline instead of vectorizing over all genes at once).

Applied here: fit and extract genes in chunks (e.g. 500-1000 genes), write each chunk's result
rows via `fwrite(..., append = TRUE)`, then `rm()` + `gc()` that chunk's model objects before
moving to the next. Peak memory becomes O(batch_size) instead of O(n_genes) — at the cost of
some fitting throughput (TMM normalization and the treatment-level factor setup would need to
happen once up front, outside the per-chunk loop, since those depend on the full sample set,
not a gene subset).

Only build this if Options 1-3 are confirmed insufficient — it's a meaningfully larger change
to a script that already works correctly elsewhere.

---

## Option 5 — Diagnostic: quantify the actual per-object cost (do anytime, cheap)

Before deciding how far up this list to go, it's worth confirming what's actually driving
per-gene GAMLSS object size for this dataset. A high-cardinality categorical covariate in
`mu_sum_covars` (`RNAlibBatch`, `RNAseqBatch`, `egglayBatch`, `platingBatch`, `well`) would
inflate the sum-coded design matrix width, and GAMLSS objects retain their design matrices —
so this would compound across all ~8,786 × 3 fitted objects, and scale with HS's larger `n`.

Add a one-time diagnostic print near the top of `run_gamlss_dvgdeg()` (`gamlss_dvgdeg_core.R`,
right after the covariate coercion loop around line 93):

```r
message("Covariate cardinality: ", paste(
  vapply(mu_sum_covars, function(c) paste0(c, "=", nlevels(meta_dt[[c]])), character(1)),
  collapse = ", "
))
```

and, after the first `fit_one()` call succeeds, log `object.size(m_full[[1]])` to see the
per-object footprint directly. This quantifies whether Option 3 (moderate win) or Option 4
(bigger win, bigger effort) is actually worth pursuing, rather than guessing.

---

## Recommendation

Try **1 → 2 → 3 → 4** in that order, stopping as soon as HS jobs complete reliably. Run Option
5's diagnostic alongside whichever option is tried first — it costs nothing and tells you
whether the remaining options are likely to help before you invest time in them. Don't jump
straight to gene-batching (Option 4): the pipeline is already proven correct on other
datasets, so the fastest path to unblocking this dataset is more resources, not more code.

---

## Related Reading

- [FIX-04-GAMLSS-Q1Q5-CT-HS-Errors.md](FIX-04-GAMLSS-Q1Q5-CT-HS-Errors.md) — root cause analysis for both the CT and HS failures
- [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) — vectorized-vs-batch memory pattern (prior art for Option 4)
- [GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md](GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md) — pipeline overview and CLI usage
- [01-RULES-Function-Code-Rules.md](01-RULES-Function-Code-Rules.md) — R function code style rules to follow if Options 3 or 4 are implemented

---

**Last Updated:** 2026-07-01
**Status:** Documented — no option implemented yet; start with Option 1 on the next HPC run
