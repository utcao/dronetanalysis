# GAMLSS Q1/Q5 CT/HS Pipeline Errors

## Overview

`tmp/logs/log_snakemake_q1q5_degdevg_gamlss_0622.txt` records a 45-job SGE run of
`compute_degdvg_gamlss_q1q5.snakemake` (2026-06-23 → 2026-06-26, `mem_gb=25` per job). **19 of
45 jobs failed.** Two distinct failure modes appear — confirmed against
`tmp/logs/FBgn0001233_CT.log` and `tmp/logs/FBgn0001233_HS.log`, the two example per-job logs
for the same focal gene under each condition:

| Condition | Samples (n) | Failures | Failure mode |
|---|---|---|---|
| CT | 376 | 5 | `subscript out of bounds` inside `foreach %do%` |
| HS | 416 | 13 | `fwrite()` "Cannot allocate memory" |

`FBgn0001233` fails in **both** conditions, but for the two different reasons below — this
confirms they are separate bugs, not one root cause. HS has more samples than CT, which is
consistent with HS hitting memory limits more often (larger per-gene GAMLSS objects) while CT
mostly hits the fitting/indexing issue.

---

## Issue 1 — CT: `subscript out of bounds`

### Problem

```
ERROR [2026-06-23 12:47:04] Error in fitting model for gene 3027 FBgn0051014: NA's in the working vector or weights for parameter sigma

Error in { : task 1 failed - "subscript out of bounds"
Calls: run_gamlss_dvgdeg -> %do% -> <Anonymous>
```

(`tmp/logs/FBgn0001233_CT.log`)

Full run, 45-job log: CT failures were `FBgn0262739`, `FBgn0087035`, `FBgn0001233`,
`FBgn0004595`, `FBgn0016917`.

### Root Cause

`fit_one()` (`src/q1q5_degdvg/q1q5_gamlss/tools_gamlss.R:25-35`) already catches per-gene
GAMLSS fitting failures with `tryCatch` and returns `NULL` for that gene — the logged
`ERROR [...] NA's in the working vector...` line is *handled*, not the crash.

The crash itself is `foreach`'s standard error-wrapping behavior: both `%do%` and `%dopar%`
wrap any error raised inside the loop body as `"task N failed - <message>"`, regardless of
backend. So `"subscript out of bounds"` originates somewhere in the per-treatment extraction
block of `run_gamlss_dvgdeg()` (`gamlss_dvgdeg_core.R:180-293`) — specifically, code that
assumed every gene's fitted model has the same coefficient layout. When a gene's model matrix
has fewer columns (e.g. a batch covariate level absent after the Q1/Q5 subset drops a large
share of samples), indexing by a position derived from a different gene's model goes out of
bounds.

### Status: Already Fixed

`gamlss_dvgdeg_core.R:185-199` already replaces the global, gene-1-derived coefficient index
with a **per-gene** lookup:

```r
# single-model t-tests from the m_full summary table.
# coef_idx is looked up per gene rather than globally from gene-1's summary:
# when a gene's model matrix has fewer columns (dropped batch levels), its
# summary has fewer rows and a global index derived from gene 1 goes out of bounds.
p_mu_single <- vapply(m_full_sum, function(m) {
  if (!length(m)) return(NA_real_)
  idx <- grep(coef_trt, rownames(m))
  if (length(idx) < 1) return(NA_real_)
  m[idx[1], "Pr(>|t|)"]
}, numeric(1))
```

This matches the fix already applied ("I have applied the parameter index based on the
results of genes themself"). All other coefficient access in the extraction block
(`m$mu.coefficients[coef_trt]`, `m$sigma.coefficients[coef_trt]`, etc.) indexes by **name**,
not position, so a missing name safely returns `NA` rather than erroring.

### Residual Recommendation (optional follow-up — not applied)

The extraction block still runs for all genes in one `foreach` iteration with no per-gene
`tryCatch`. If a still-unknown edge case throws inside that block again, the entire job exits
via `Execution halted` and **every already-computed gene row for that job is lost** (results
are only written once, at the very end, via `fwrite`). If the subscript error recurs after
re-running the affected genes below, the next thing to try is wrapping the per-gene body of
the extraction (or the whole `foreach` iteration) in `tryCatch`, logging the offending gene
and continuing rather than aborting the whole 240h job.

### Verification

Re-run the previously-failed CT genes and confirm they complete without the earlier fix:

```bash
for gene in FBgn0262739 FBgn0087035 FBgn0001233 FBgn0004595 FBgn0016917; do
  /tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript \
    src/q1q5_degdvg/q1q5_gamlss/run_gamlss_quintile_analysis.R \
    --focal-gene "$gene" --full-run --condition-val 1 \
    --out-dir ../../results/gamlss_0609
done
```

---

## Issue 2 — HS: `fwrite()` "Cannot allocate memory"

### Problem

```
Fitting GAMLSS models...
Error in fwrite(dvdg_tab, file = file.path(gamlss_dir, glue("{out_prefix}_gamlss_dvgdeg.csv"))) :
  Unable to allocate 8 MB * 1 thread buffers; '12: Cannot allocate memory'. Please read ?fwrite for nThread, buffMB and verbose options.
Calls: run_gamlss_dvgdeg -> fwrite
Execution halted
Warning message:
system call failed: Cannot allocate memory
```

(`tmp/logs/FBgn0001233_HS.log`)

Full run, 45-job log: 13 HS jobs failed this way (vs. 5 CT failures from Issue 1), consistent
with HS's larger sample size (n=416 vs n=376 for CT) pushing per-gene GAMLSS objects closer to
the 25 GB job ceiling.

### Root Cause

All fitting for a job happens before any results are written:

```r
m_full      <- bp_apply(seq_len(n_genes), fit_one, ...)   # gamlss_dvgdeg_core.R:133-139
m_full_sum  <- bp_apply(m_full, fit_sum_quiet)             # :141
m_sigma_int <- bp_apply(seq_len(n_genes), fit_one, ...)    # :144-150
m_mu_nontrt <- bp_apply(seq_len(n_genes), fit_one, ...)    # :153-159
```

Three full-size lists of ~8,786 fitted GAMLSS objects (`m_full`, `m_sigma_int`,
`m_mu_nontrt`) — plus `m_full_sum` — are all held **simultaneously** in memory for the entire
fitting + extraction phase, even though `m_sigma_int` and `m_mu_nontrt` are only ever read for
two scalars per gene (`logLik()` and `df.fit`, inside `map2_dbl` at lines 220-226 and
246-252). Peak RSS occurs during this phase, not at `fwrite()` time — the `fwrite` failure is
just where the process finally has too little contiguous memory left to service even an 8 MB
buffer request.

### Why the Already-Tried Fix Didn't Help

`gamlss_dvgdeg_core.R:295-299` already does exactly what was tried ("I remove some
variables"):

```r
# Free the model lists — all values have been extracted into result_list.
# Calling gc() here reclaims the memory used by 3×n_genes GAMLSS objects
# before fwrite allocates its write buffers, preventing OOM on large sample sets.
rm(m_full, m_full_sum, m_sigma_int, m_mu_nontrt)
gc()
```

This runs **after** peak memory has already been reached — it frees memory before `fwrite`,
but by the time the process gets there under HS's larger sample size, the OS/cgroup ceiling
has likely already been hit (and/or heap fragmentation prevents even small allocations from
succeeding). Removing variables at the end doesn't lower the peak; it only affects what's live
at the very end of the job.

### Status: Open

See [OPTIMIZATION-05-GAMLSS-Q1Q5-Memory-Mitigation.md](OPTIMIZATION-05-GAMLSS-Q1Q5-Memory-Mitigation.md)
for a ranked series of mitigations, ordered from cheapest/lowest-risk (increase cluster
memory — this pipeline is already known to work correctly on other, smaller datasets) to most
invasive (restructure the fitting loop). Try them in order; stop as soon as jobs complete
reliably.

---

## Summary

| File | Line(s) | Issue | Status |
|---|---|---|---|
| `src/q1q5_degdvg/q1q5_gamlss/gamlss_dvgdeg_core.R` | 185-199 | CT: coefficient index derived from gene 1 caused out-of-bounds access for genes with fewer model-matrix columns | ✅ Fixed (per-gene `grep()` lookup) |
| `src/q1q5_degdvg/q1q5_gamlss/gamlss_dvgdeg_core.R` | 180-293 | CT: no per-gene error isolation in the extraction `foreach` body | ⚠️ Optional follow-up, not applied |
| `src/q1q5_degdvg/q1q5_gamlss/gamlss_dvgdeg_core.R` | 133-159, 295-299 | HS: three full-size model lists held concurrently; `rm()+gc()` runs after peak memory | ❌ Open — see mitigation doc |

---

## Related Reading

- [GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md](GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md) — pipeline overview, architecture, CLI usage
- [OPTIMIZATION-05-GAMLSS-Q1Q5-Memory-Mitigation.md](OPTIMIZATION-05-GAMLSS-Q1Q5-Memory-Mitigation.md) — ranked memory mitigation options for Issue 2
- [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) — prior art: vectorized-vs-batch memory pattern used elsewhere in this repo
- [01-RULES-Function-Code-Rules.md](01-RULES-Function-Code-Rules.md) — R function code style rules to follow if/when the residual/mitigation items above are implemented

---

**Last Updated:** 2026-07-01
**Status:** Issue 1 fixed and awaiting re-run verification; Issue 2 open, mitigations documented
