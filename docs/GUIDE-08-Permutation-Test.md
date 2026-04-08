# Permutation Test Guide

## Overview

`src/pipelines/Snakefile_permutation` orchestrates a permutation test that
justifies observed differential co-expression metrics by comparing them against
a **null distribution** generated from random sample assignments.

The core question: are the high-vs-low differences (rewiring, edge counts, effect
sizes) driven by the true **expression gradient**, or could they arise from any
random partition of samples?

---

## Prerequisites

- Stage 3 (`reconstruct_single`) must have completed for all target focus genes.
  The permutation test reads `networks/{gi}_{gene_id}.h5` and requires the
  `focus_gene/metrics/` attribute group written by Stage 3.
- `expression.h5`, `bootstrap_indices.h5`, `rewiring_hubs.tsv` from the
  Snakefile_bootstrap run.

---

## Method

### Algorithm per permutation replicate

1. Draw `k_low + k_high` samples uniformly at random from **all** samples
   (not restricted to the expression-extreme pool).
2. Split into `pseudo_low` (`k_low`) and `pseudo_high` (`k_high`).
3. Rank-normalize **all N genes** within each pseudo-group (required for
   correct Spearman). Cost: O(N × k).
4. Extract only the focus gene's row of correlations:
   ```python
   r_low  = (z_low[g]  @ z_low.T)  / (k_low - 1)   # (N,)
   r_high = (z_high[g] @ z_high.T) / (k_high - 1)
   ```
   This avoids computing the full N×N matrix — O(N × k) instead of O(N² × k).
5. Apply Fisher's Z-test + BH FDR on the N-1 focus-gene pairs.
6. Compute focus-gene metrics identical to Stage 3b.
7. (Optional) Compute L1×L1 sub-matrix for L2 expansion metrics.

### Metrics tested

| Metric | Description |
|--------|-------------|
| `n_diff_edges` | Number of significant differential edges |
| `focus_deg_low` | Focus gene degree in LOW network |
| `focus_deg_high` | Focus gene degree in HIGH network |
| `mean_abs_delta` | Mean \|Δr\| over differential edges |
| `max_abs_delta` | Max \|Δr\| over differential edges |
| `L1_n_nodes` | Number of direct differential partners |
| `L1_rewire` | L1 rewiring count (disappear + new + sign_change) |
| `L1_frac_rewire` | L1 rewiring fraction |
| `L2L1_deg` *(optional)* | Pure L2 nodes / L1 nodes |
| `L2L1_rewire` *(optional)* | L2 rewiring / L1 rewiring |

### Empirical p-value

Two-sided, mirroring `utils_permutation_net.R::permutate_p_two()`:

```
p = max( sum(|null| >= |obs|) / n_perm,  1/n_perm )
```

Floor at `1/n_perm` so p is never exactly 0.

### FDR note

The permutation applies BH correction on N-1 ≈ 8762 tests (focus-gene row
only). The real pipeline applies BH on N(N-1)/2 ≈ 38M tests. With `fdr_alpha=0.7`
(default) this asymmetry is negligible. For strict alphas the permutation is
conservative — the null n_diff_edges may be slightly larger than observed, making
it harder to reach significance.

---

## Scripts

| Script | Purpose |
|--------|---------|
| `src/scripts/10spearman_corr/07_permutation_test.py` | Per-gene: builds null distribution |
| `src/scripts/10spearman_corr/08_collect_permutation_pvals.py` | Aggregates → `permutation_pvals.tsv` |
| `src/scripts/15analysis/09_plot_permutation_null.R` | Histogram plots per gene |
| `src/pipelines/Snakefile_permutation` | Orchestration (independent of Snakefile_bootstrap) |

---

## Usage

### Step 1: Run the permutation test

```bash
# Dry run — check DAG without executing
snakemake -s code/dronetanalysis/src/pipelines/Snakefile_permutation \
    --config \
        expr_h5=results/expression.h5 \
        networks_dir=results/networks \
        indices_h5=results/bootstrap_indices.h5 \
        summary_tsv=results/rewiring_hubs.tsv \
        "perm_genes=[FBgn0001197,FBgn0262739]" \
        n_permutations=100 \
    -n

# Local run (4 parallel genes)
snakemake -s code/dronetanalysis/src/pipelines/Snakefile_permutation \
    --config \
        expr_h5=results/expression.h5 \
        networks_dir=results/networks \
        indices_h5=results/bootstrap_indices.h5 \
        summary_tsv=results/rewiring_hubs.tsv \
        "perm_genes=[FBgn0001197,FBgn0262739]" \
        n_permutations=100 \
    -j 4

# SGE cluster
snakemake -s code/dronetanalysis/src/pipelines/Snakefile_permutation \
    --configfile config/ct_voom_snakemake.yaml \
    --config \
        "perm_genes=[FBgn0001197,FBgn0262739]" \
        n_permutations=100 \
    --profile config/sge_profile
```

### Step 2: Collect p-values (automatic in Snakefile)

If running manually:

```bash
python src/scripts/10spearman_corr/08_collect_permutation_pvals.py \
    --null-dir results/permutation_null \
    --observed-tsv results/rewiring_hubs.tsv \
    --out-pvals-tsv results/permutation_pvals.tsv \
    --out-augmented-tsv results/rewiring_hubs_with_pvals.tsv \
    --alpha 0.05
```

### Step 3: Plot null distributions

```bash
Rscript src/scripts/15analysis/09_plot_permutation_null.R \
    --null-tsv-dir results/permutation_null \
    --pvals-tsv results/permutation_pvals.tsv \
    --output-dir results/visualization_data/perm_null \
    --condition-label "Control" \
    --metrics "n_diff_edges,L1_rewire,L1_frac_rewire,mean_abs_delta"
```

---

## Configuration Reference

All parameters are passed via `--config` or `--configfile`. Required keys must be
provided; optional keys have defaults shown.

### Required

| Key | Description |
|-----|-------------|
| `expr_h5` | Path to `expression.h5` from Snakefile_bootstrap Stage 0 |
| `networks_dir` | Path to `networks/` directory from Stage 3 |
| `indices_h5` | Path to `bootstrap_indices.h5` from Stage 1 |
| `summary_tsv` | Path to `rewiring_hubs.tsv` from Stage 3b |
| `perm_genes` | List of FBgn IDs to test, e.g. `[FBgn0001197, FBgn0262739]` |

### Optional

| Key | Default | Description |
|-----|---------|-------------|
| `perm_out_dir` | `<summary_tsv_dir>/permutation_null` | Output directory |
| `n_permutations` | `100` | Number of shuffles per gene |
| `perm_seed` | `9999` | Base seed; gene `i` uses `perm_seed + gene_index` |
| `fdr_alpha` | `0.7` | BH FDR threshold (should match pipeline run) |
| `corr_threshold` | `0.0001` | Minimum \|r\| to consider edge "present" |
| `compute_l2` | `False` | Compute L1×L1 sub-matrix for L2 metrics |
| `perm_alpha` | `0.05` | Significance threshold for `any_sig` column |
| `skip_perm_plot` | `False` | Skip R plotting step |
| `condition_label` | `"Condition"` | Label used in plot titles |
| `perm_plot_metrics` | `""` | Comma-separated metric names to plot (empty = all) |

---

## Outputs

### Directory layout

```
results/
├── permutation_null/
│   ├── {gi}_{gene_id}.h5                   ← null distributions (HDF5)
│   └── {gi}_{gene_id}_null_dist.tsv.gz     ← null distributions (TSV for R)
├── permutation_pvals.tsv                   ← empirical p-values, all genes
├── rewiring_hubs_with_pvals.tsv            ← rewiring_hubs.tsv + pval_* columns
└── visualization_data/perm_null/
    └── {gi}_{gene_id}_permutation_null.pdf ← one PDF per gene
```

### `permutation_pvals.tsv` columns

| Column | Description |
|--------|-------------|
| `gene_id` | Gene FlyBase ID |
| `gene_idx` | Gene index in expression matrix |
| `n_permutations` | Number of replicates run |
| `skipped` | True if gene had 0 significant edges in Stage 3 |
| `obs_{metric}` | Observed value for each metric |
| `pval_{metric}` | Empirical two-sided p-value |
| `null_mean_{metric}` | Mean of null distribution |
| `null_std_{metric}` | SD of null distribution |
| `null_all_zero_{metric}` | True if null is entirely zero |
| `n_metrics_sig` | Count of metrics with p < alpha |
| `any_sig` | Boolean: any metric significant |

### `{gi}_{gene_id}.h5` layout

```
meta/  (attrs)
    gene_index, gene_id, n_permutations, n_genes, n_samples
    k_low, k_high, fdr_alpha, corr_threshold, seed, compute_l2, skipped

observed/  (attrs — loaded from Stage 3 network file)
    n_diff_edges, focus_deg_low, focus_deg_high
    mean_abs_delta, max_abs_delta
    L1_n_nodes, L1_rewire, L1_frac_rewire
    [L2L1_deg, L2L1_rewire — if compute_l2]

null/  (datasets, shape: (n_permutations,), float32/int32)
    n_diff_edges, focus_deg_low, focus_deg_high, ...
```

---

## Resource Estimates

For N=8763 genes, k_low = k_high ≈ 280 samples:

| n_permutations | Time per gene | Memory per job |
|----------------|--------------|----------------|
| 100 | ~30 seconds | ~200 MB |
| 1000 | ~5–8 minutes | ~200 MB |

With 15–19 focus genes running in parallel on SGE:
- 100 permutations → ~1 minute wall-clock
- 1000 permutations → ~8 minutes wall-clock

---

## Interpreting Results

A gene passes the permutation test if its observed metric is **significantly
larger** than what random sample assignment would produce.

**Example interpretation:**

```
gene_id:     FBgn0001197
obs_L1_rewire:   45
null_mean_L1_rewire:  3.2
null_std_L1_rewire:   1.8
pval_L1_rewire:   0.001

→ The observed rewiring (45) is far above the null mean (3.2).
  Only 0.1% of random partitions produce rewiring ≥ 45.
  Conclusion: the high/low difference in this gene's neighbourhood is driven
  by the true expression gradient, not random chance.
```

**Null all-zero flag:**

When `null_all_zero_{m}=True`, all 100/1000 permutation replicates produced
zero for that metric — the null distribution is degenerate. This typically
means the metric is rarely non-zero under random assignment (strong signal),
so p = `1/n_perm` (the minimum possible p-value).

---

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|---------|
| Gene not found in networks/ | Stage 3 not run for that gene | Add gene to `gene_subset` in Snakefile_bootstrap config and re-run Stage 3 |
| `skipped=True` in output | Gene had 0 significant edges in Stage 3 | Check `rewiring_hubs.tsv` — gene may have very few differential edges at the chosen FDR threshold |
| All p-values = 1.0 | n_permutations too low | Increase `n_permutations` to 1000 |
| High null mean > observed | FDR asymmetry or low `fdr_alpha` | Expected behavior with strict alpha; documented asymmetry |
| Missing `focus_gene/metrics` in network HDF5 | Old Stage 3 output | Re-run Stage 3 (`reconstruct_single`) for affected genes |

---

## Related Reading

- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) - Full pipeline overview including Stage 7
- [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) - Snakefile_bootstrap config reference
- [GUIDE-04-Qualitative-Change-Metrics.md](GUIDE-04-Qualitative-Change-Metrics.md) - L1/L2 rewiring metric definitions
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) - Fisher's Z-test and FDR methodology

---

**Last Updated:** 2026-04-08
**Status:** ✅ Active
