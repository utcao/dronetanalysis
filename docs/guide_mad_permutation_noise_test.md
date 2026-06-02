# Guide: Permutation Validation of MAD Transcriptomic Noise

**Category:** Guide — Statistical validation of focus-gene-driven transcriptomic variability

---

## Problem

`compute_mad_transcriptomic_variability.R` tests whether a focus gene's expression gradient
(LOW vs HIGH expression samples) is associated with higher MAD (median absolute deviation) in
the rest of the transcriptome.  The test is a parametric Wilcoxon rank-sum followed by BH-FDR
correction.

In practice, this produces a very large number of significant genes even after FDR correction.
Two factors inflate the result:

1. The Wilcoxon test is applied to vectors of ~10 000 gene-level MAD values.  With that many
   data points the test is extremely powerful; even tiny, biologically irrelevant MAD differences
   yield p < 0.05.
2. The BH procedure corrects for multiple testing *across* focus genes, but the within-gene
   test is already over-powered.

A permutation null built from random sample partitions of equal size provides a more appropriate
reference: it quantifies how large a MAD difference arises under any arbitrary group split, so
significance can only be claimed when the expression-rank-based split is *more extreme* than
random.

---

## What "Breaking the Link" Means

In the original analysis the LOW/HIGH group assignment is determined entirely by the focus gene
G's expression rank:

```
LOW  = samples where G is expressed least (bottom k%)
HIGH = samples where G is expressed most  (top k%)
```

The permutation destroys this rank-to-sample link by **shuffling G's expression values across
samples** before re-ranking.  After the shuffle, the values are the same but the mapping
`value → sample` is broken.  When re-ranked, the derived LOW/HIGH split is effectively a
random partition of equal sizes (k_low, k_high) — not driven by G's biology.

```
expr_G_perm  <- sample(expr_G)           # break the value–sample identity
perm_low     <- order(expr_G_perm)[1:k_low]
perm_high    <- order(expr_G_perm, decreasing=TRUE)[1:k_high]
```

Repeating this N times gives the null distribution of `delta_mean_mad` under random group
assignment.  The two-sided permutation p-value is:

```
perm_p = #{|null_delta| >= |observed_delta|} / N
```

A focus gene is **validated** when it is significant under both the original parametric test
(BH-corrected) and the permutation test (BH-corrected).

---

## Script

`code/dronetanalysis/src/scripts/15analysis/compute_mad_permute_transcriptomic_noise.R`

### Key arguments

| Argument | Default | Description |
|---|---|---|
| `--expr-file` | required | VOOM expression matrix (genes × samples, tab-separated) |
| `--mapping-file` | required | TSV/CSV with `gene_id` and `SYMBOL` columns |
| `--output-file` | required | Output `.xlsx` path |
| `--null-dist-dir` | `NULL` | If set, writes per-gene `*_null_dist.tsv.gz` and `permutation_pvals.tsv` for `plot_permutation_null_dist.R` |
| `--focus-genes` | `NULL` | Comma-separated gene IDs to restrict permutation (saves compute) |
| `--low-frac` | `0.2` | Fraction of samples in LOW group |
| `--high-frac` | `0.2` | Fraction of samples in HIGH group |
| `--n-perms` | `500` | Permutations per gene |
| `--seed` | `42` | Random seed |
| `--condition-label` | `"Condition"` | Label added to output |

### Usage examples

**Full run on all genes (overnight):**
```bash
Rscript compute_mad_permute_transcriptomic_noise.R \
  --expr-file    data/processed/VOOM/voomdataCtrl.txt \
  --mapping-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --output-file  results/variability/mad_permutation_ct.xlsx \
  --n-perms      500 \
  --seed         42 \
  --condition-label "Control"
```

**Restricted to originally significant genes (fast, recommended first pass):**
```bash
# Extract sig gene IDs from the original xlsx first, then:
Rscript compute_mad_permute_transcriptomic_noise.R \
  --expr-file    data/processed/VOOM/voomdataCtrl.txt \
  --mapping-file results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --output-file  results/variability/mad_permutation_ct_sig.xlsx \
  --focus-genes  "FBgn001,FBgn002,FBgn003" \
  --n-perms      500
```

**With null distribution TSVs (for plotting):**
```bash
Rscript compute_mad_permute_transcriptomic_noise.R \
  --expr-file      data/processed/VOOM/voomdataCtrl.txt \
  --mapping-file   results/result_voomct/rewiring_hubs_ct_anno_0408_2026.tsv \
  --output-file    results/variability/mad_permutation_ct.xlsx \
  --null-dist-dir  results/variability/null_distributions_ct \
  --n-perms        500

# Then plot with the existing script:
Rscript plot_permutation_null_dist.R \
  --null-tsv-dir  results/variability/null_distributions_ct \
  --pvals-tsv     results/variability/null_distributions_ct/permutation_pvals.tsv \
  --output-dir    results/variability/null_distribution_plots \
  --metrics       "delta_mean_mad"
```

---

## Output Columns

| Column | Description |
|---|---|
| `gene_id`, `SYMBOL` | Gene identifiers |
| `obs_delta_mean_mad` | Observed mean(MAD_HIGH) − mean(MAD_LOW) |
| `obs_mean_mad_low` | Mean MAD across other genes in the LOW group |
| `obs_mean_mad_high` | Mean MAD across other genes in the HIGH group |
| `obs_wilcoxon_p` | Parametric Wilcoxon p-value (for cross-reference with original analysis) |
| `obs_p_adj` | BH-corrected parametric p-value |
| `perm_p` | Two-sided permutation p-value |
| `perm_p_adj` | BH-corrected permutation p-value |
| `n_perms` | Number of permutations used |
| `validated` | `TRUE` if both `obs_p_adj < 0.05` AND `perm_p_adj < 0.05` |
| `condition_label` | Condition label from `--condition-label` |

Results are sorted by `perm_p_adj` ascending.

---

## Interpreting the Results

**Expected outcome if inflation is real:**  
`n_sig` with `obs_p_adj < 0.05` will be much larger than `n_sig` with `perm_p_adj < 0.05`.
Many genes will have a low observed Wilcoxon p but a high permutation p, meaning any random
split of the same size would produce a comparable MAD difference.

**Diagnostic plot:**  
Plot `−log10(perm_p_adj)` on the y-axis vs `−log10(obs_p_adj)` on the x-axis.  Genes above
the diagonal (y > x) have stronger permutation evidence than parametric; genes below are
parametrically significant but not supported by the permutation null.  Only genes above a
`−log10` threshold on both axes should be considered confidently validated.

**Practical threshold:**  
`validated == TRUE` (both FDR < 0.05) is the recommended filter for downstream analysis.

---

## Computational Notes

- Each gene requires `n_perms` re-rankings and MAD computations across all other genes.
  With ~10 000 genes and 500 permutations, a full run takes several hours.
- If `matrixStats` is installed, the script uses `rowMads()` automatically (~10–50× faster
  than `apply(..., mad)`).  Install with `install.packages("matrixStats")`.
- Use `--focus-genes` to limit the run to originally significant genes for a quick first pass.
- The global seed is set once before the loop, so results are reproducible.  Using
  `--focus-genes` will produce different null distributions for the selected genes than a
  full run (different shuffle draws), but this is acceptable for an exploratory validation.
