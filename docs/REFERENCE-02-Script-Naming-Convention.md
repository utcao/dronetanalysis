# Script Naming Convention Reference

## Overview

This document defines the naming convention applied to all R and Python scripts
in `src/scripts/15analysis/` as of 2026-04-27. The convention makes a script's
role immediately clear from its filename and groups related scripts naturally in
directory listings.

---

## Pattern

```
<verb>_<metric/method>_<description>.{R,py}
```

| Position | Role | Examples |
|---|---|---|
| **verb** | What the script *does* | `plot_`, `compute_`, `analyze_`, `run_` |
| **metric/method** | The primary data object or analytical method (anchor noun) | `mad_`, `itv_`, `pca_`, `hub_`, `conditions_`, `permutation_` |
| **description** | Further qualifier — what aspect or output | `_gene_comparison`, `_variability_ranks`, `_transcriptomic_variability` |

### Verb vocabulary

| Verb | When to use |
|---|---|
| `plot_` | Produces one or more figures / PDFs as primary output |
| `compute_` | Produces a numerical table (xlsx / TSV / CSV) as primary output, possibly with console stats |
| `analyze_` | Multi-output script: both figures and statistical tables |
| `run_` | Batch or pipeline execution wrapper — calls other scripts or tools |

### Metric/method as anchor

The second term is the consistent noun that anchors the filename to a specific
analytical context. Scripts belonging to the same analysis family share the same
anchor (`mad_`, `itv_`, `pca_`, …), so they sort together in directory listings.

---

## Applied Renames (2026-04-27)

The table below shows all files affected by this convention change. Files
starting with `module_` and `visualise_network.R` were left unchanged for
continuity.

| Old filename | New filename | Anchor | Change |
|---|---|---|---|
| `09_plot_permutation_null.R` | `plot_permutation_null_dist.R` | `permutation` | Removed stray numeric prefix; added `_dist` (histogram of null *distribution*) |
| `compare_conditions.R` | `analyze_conditions.R` | `conditions` | Added verb; multi-output → `analyze_` |
| `compute_full_mad_cv2_ranks.R` | `compute_mad_variability_ranks.R` | `mad` | Internal metric names (`full`, `cv2`) moved to comments |
| `pathway_enrichment_hubs.R` | `analyze_hub_enrichment.R` | `hub` | Verb-first; subject (`hub`) before method (`enrichment`) |
| `pca_l2l1_variability.R` | `plot_pca_l2l1_variability.R` | `pca` | Verb-first; `l2l1` is implementation detail |
| `plot_expr_violin.py` | `plot_expression_violin.py` | `expression` | Expanded abbreviation |
| `plot_gene_mad_variability.R` | `plot_mad_gene_comparison.R` | `mad` | `mad` is the anchor; `gene_comparison` = LOW vs HIGH group |
| `plot_sample_variability.R` | `plot_itv_sample_comparison.R` | `itv` | Named the metric (ITV); `sample_comparison` = LOW vs HIGH |
| `summarize_all_genes_mad_variability.R` | `compute_mad_transcriptomic_variability.R` | `mad` | `compute_` verb (produces table); `all_genes` removed as redundant |
| `run_variability_batch.py` | `run_variability_batch.py` | `variability` | No change — already consistent |

---

## Resulting Directory Layout

```
src/scripts/15analysis/
├── analyze_conditions.R                   # compare_conditions.R
├── analyze_hub_enrichment.R               # pathway_enrichment_hubs.R
├── compute_mad_transcriptomic_variability.R  # summarize_all_genes_mad_variability.R
├── compute_mad_variability_ranks.R        # compute_full_mad_cv2_ranks.R
├── module_enrichment.R                    (unchanged)
├── module_overlap_analysis.R              (unchanged)
├── module_preservation_netrep.R           (unchanged)
├── plot_expression_violin.py              # plot_expr_violin.py
├── plot_itv_sample_comparison.R           # plot_sample_variability.R
├── plot_mad_gene_comparison.R             # plot_gene_mad_variability.R
├── plot_pca_l2l1_variability.R                 # pca_l2l1_variability.R
├── plot_permutation_null_dist.R           # 09_plot_permutation_null.R
├── run_variability_batch.py               (unchanged)
├── visualise_network.R                    (unchanged)
└── figures/                               (new — ad-hoc post-pipeline vis scripts)
    └── plot_mad_rank_shift.R              (new)
```

### The `figures/` subfolder

Ad-hoc, post-pipeline visualization scripts that consume workflow outputs but do
not feed back into the pipeline. These are run manually to generate
publication-quality or exploratory figures. Future scripts of the same type go
here, following the same naming convention.

---

## Extending the Convention

When adding a new script to `15analysis/`:

1. **Choose verb** — does it produce figures (`plot_`), a table (`compute_`),
   both (`analyze_`), or wrap other tools (`run_`)?
2. **Choose anchor** — what is the primary metric or method? If MAD-related,
   use `mad_`; if ITV-related, use `itv_`; new methods get their own anchor.
3. **Add description** — brief qualifier: scope (gene, sample, module),
   output type (ranks, comparison, variability, enrichment), or target.
4. **Place ad-hoc vis scripts in `figures/`**, not the main folder.

---

## Related Reading

- [GUIDE-12-15Analysis-Overview.md](GUIDE-12-15Analysis-Overview.md) — Full inventory of all scripts and their purpose
- [GUIDE-05-Expression-Variability-Analysis.md](GUIDE-05-Expression-Variability-Analysis.md) — MAD and ITV scripts
- [GUIDE-11-PCA-Gene-Metrics.md](GUIDE-11-PCA-Gene-Metrics.md) — PCA and MAD variability ranks scripts
- [00-RULES-Documentation-Standards.md](00-RULES-Documentation-Standards.md) — Docs naming convention

---

**Last Updated:** 2026-04-27
**Status:** ✅ Active standard for `src/scripts/15analysis/`
