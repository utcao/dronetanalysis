# Statistical Methods and Pipeline Efficiency Reference

## Overview

This document provides comprehensive technical documentation on the statistical methodology and pipeline design for differential co-expression analysis. It covers:

- **Part I**: Statistical testing methods for differential correlations
- **Part II**: Efficient pipeline design and implementation strategies

---

## Table of Contents

### Part I: Statistical Methods for Differential Co-Expression
1. [Problem Setup](#part-i-statistical-methods-for-differential-co-expression)
2. [Test Options Summary](#test-options-summary)
3. [Detailed Method Analysis](#detailed-analysis)
   - Bootstrap Difference Test (Recommended)
   - Fisher's Z-Transform Test
   - Permutation Test
   - Bayesian Credible Intervals
4. [Efficiency Comparison](#efficiency-comparison)
5. [Hypothesis Testing Framework](#hypothesis-testing-framework)
6. [Multiple Testing Correction](#multiple-testing-correction)
7. [Practical Recommendations](#practical-recommendations)
8. [Example Implementation](#example-implementation)
9. [Biological Interpretation](#biological-interpretation)

### Part II: Efficient Pipeline Design
1. [Questions Answered](#part-ii-efficient-pipeline-design-and-analysis)
2. [Pipeline Comparison](#proposed-efficient-pipeline)
3. [Efficiency Analysis](#efficiency-comparison-original-vs-proposed)
4. [Edge Selection Strategies](#selecting-significant-edges-strategy-options)
5. [Sparse Matrix Reconstruction](#reconstructing-correlation-matrices-from-significant-edges)
6. [Implementation Details](#recommended-pipeline-implementation)
7. [Implemented Pipeline](#implemented-pipeline-final)

---

# Part I: Statistical Methods for Differential Co-Expression

## Problem Setup

Given bootstrap correlation data for each edge:
- `r_low_base`: Base correlation in low expression subpopulation
- `r_low_boot`: [n_bootstrap] bootstrap correlations in low subpop
- `r_high_base`: Base correlation in high expression subpopulation
- `r_high_boot`: [n_bootstrap] bootstrap correlations in high subpop

**Goal**: Test if the correlation strength differs significantly between subpopulations.

---

## Test Options Summary

| Method | Null Hypothesis | Speed | Assumptions | Best For |
|--------|----------------|-------|-------------|----------|
| **1. Bootstrap Difference** | diff = 0 | ★★★★★ | Few | General use, robust |
| **2. Fisher's Z-transform** | Z_low = Z_high | ★★★★★ | Normality | Large samples, theoretical |
| **3. Permutation Test** | Exchangeable | ★★☆☆☆ | None | Gold standard, slow |
| **4. Bayesian Credible Intervals** | Overlapping CIs | ★★★☆☆ | Prior choice | Interpretation, uncertainty |

---

## Detailed Analysis

### 1. Bootstrap Difference Test (RECOMMENDED)

**Method**:
```
delta_boot[i] = r_high_boot[i] - r_low_boot[i]  # for i in 1..n_bootstrap
CI_95 = percentile(delta_boot, [2.5, 97.5])
p_value ≈ 2 * min(P(delta_boot > 0), P(delta_boot < 0))
```

**Pros**:
✅ **Fast**: O(n_bootstrap) per edge, ~50 operations
✅ **No assumptions**: Non-parametric, works for any distribution
✅ **Interpretable**: "95% CI of difference: [0.12, 0.28]"
✅ **Handles ties**: Works even if correlations near 0 or 1
✅ **Already computed**: You have the bootstrap samples!

**Cons**:
⚠️ Assumes bootstrap replicates are independent (they should be)
⚠️ With n_bootstrap=50, minimum p-value is ~0.02 (resolution limit)

**Implementation**:
```python
def bootstrap_difference_test(r_low_boot, r_high_boot, alpha=0.05):
    """
    Test if correlation differs between subpopulations using bootstrap.

    Returns
    -------
    delta_mean : float
        Mean difference (high - low)
    ci_low, ci_high : float
        95% confidence interval of difference
    p_value : float
        Two-tailed p-value
    significant : bool
        True if 0 not in CI (equivalent to p < alpha)
    """
    delta_boot = r_high_boot - r_low_boot
    delta_mean = np.mean(delta_boot)

    # Confidence interval (percentile method)
    ci_low, ci_high = np.percentile(delta_boot, [100*alpha/2, 100*(1-alpha/2)])

    # P-value: proportion of bootstrap samples crossing 0
    p_pos = np.mean(delta_boot > 0)
    p_neg = np.mean(delta_boot < 0)
    p_value = 2 * min(p_pos, p_neg)

    # Significant if CI excludes 0
    significant = (ci_low > 0) or (ci_high < 0)

    return delta_mean, ci_low, ci_high, p_value, significant
```

**When to use**: DEFAULT choice for most analyses

---

### 2. Fisher's Z-Transform Test

**Method**:
```
Z_low = 0.5 * log((1 + r_low) / (1 - r_low))   # Fisher's Z transform
Z_high = 0.5 * log((1 + r_high) / (1 - r_high))

SE_diff = sqrt(1/(n_low - 3) + 1/(n_high - 3))
z_score = (Z_high - Z_low) / SE_diff
p_value = 2 * (1 - norm.cdf(abs(z_score)))
```

**Pros**:
✅ **Fastest**: O(1) per edge, no bootstrap needed (can use base correlations)
✅ **Theoretical basis**: Well-established in statistics literature
✅ **Power**: Good statistical power for large n_samples

**Cons**:
⚠️ **Assumes normality**: Z-transform assumes correlations are normally distributed
⚠️ **Breaks near ±1**: Undefined for r = ±1 (log singularity)
⚠️ **Ignores bootstrap**: Doesn't use your valuable bootstrap data
⚠️ **Independent samples**: Assumes low/high groups are independent (they're from same genes)

**Implementation**:
```python
def fisher_z_test(r_low_base, r_high_base, n_low, n_high):
    """
    Fisher's Z-transform test for correlation difference.

    Parameters
    ----------
    r_low_base, r_high_base : float
        Base correlations in each subpopulation
    n_low, n_high : int
        Sample sizes in each subpopulation
    """
    # Fisher's Z transform
    z_low = 0.5 * np.log((1 + r_low_base) / (1 - r_low_base))
    z_high = 0.5 * np.log((1 + r_high_base) / (1 - r_high_base))

    # Standard error of difference
    se_diff = np.sqrt(1/(n_low - 3) + 1/(n_high - 3))

    # Z-score and p-value
    z_score = (z_high - z_low) / se_diff
    p_value = 2 * (1 - stats.norm.cdf(np.abs(z_score)))

    return z_score, p_value
```

**When to use**: Quick screening on base correlations only, large sample sizes (n > 30)

---

### 3. Permutation Test

**Method**:
```
# For each permutation:
1. Pool all bootstrap samples: [r_low_boot, r_high_boot]
2. Randomly split into two groups of size n_bootstrap
3. Compute mean difference for this split
4. Compare observed difference to permutation distribution
```

**Pros**:
✅ **Gold standard**: Exact test under null hypothesis
✅ **No assumptions**: Non-parametric, distribution-free
✅ **Handles dependencies**: Accounts for correlation structure

**Cons**:
⚠️ **SLOW**: O(n_perm * n_bootstrap) per edge, typically n_perm = 1000-10000
⚠️ **Overkill**: For ~1M edges, this is computationally prohibitive
⚠️ **Same result**: Usually agrees with bootstrap difference test

**Implementation**:
```python
def permutation_test(r_low_boot, r_high_boot, n_perm=10000):
    """
    Permutation test for correlation difference.

    WARNING: Very slow for large number of edges.
    """
    obs_diff = np.mean(r_high_boot) - np.mean(r_low_boot)
    pooled = np.concatenate([r_low_boot, r_high_boot])
    n = len(r_low_boot)

    perm_diffs = np.zeros(n_perm)
    for i in range(n_perm):
        np.random.shuffle(pooled)
        perm_diffs[i] = np.mean(pooled[n:]) - np.mean(pooled[:n])

    p_value = np.mean(np.abs(perm_diffs) >= np.abs(obs_diff))
    return p_value
```

**When to use**: Validation on a small subset of top edges, not for genome-wide testing

---

### 4. Bayesian Credible Intervals

**Method**:
```
# Fit posterior distributions for each subpopulation
posterior_low ~ Normal(mean(r_low_boot), std(r_low_boot))
posterior_high ~ Normal(mean(r_high_boot), std(r_high_boot))

# Sample from difference distribution
delta_samples = posterior_high - posterior_low
CI_95 = percentile(delta_samples, [2.5, 97.5])
prob_positive = P(delta_samples > 0)
```

**Pros**:
✅ **Interpretable**: "95% probability difference is positive"
✅ **Uncertainty quantification**: Full posterior distribution
✅ **Effect size**: Naturally incorporates magnitude of difference

**Cons**:
⚠️ **Complexity**: Requires MCMC or sampling framework
⚠️ **Prior choice**: Results depend on prior assumptions
⚠️ **Slower**: More computation than bootstrap difference
⚠️ **Subjective**: Different priors → different conclusions

**When to use**: High-stakes decisions, need for uncertainty quantification

---

## Efficiency Comparison

For **N_edges = 1,000,000** (typical gene-gene pairs):

| Method | Time per Edge | Total Time | Memory |
|--------|---------------|------------|--------|
| Bootstrap Difference | 10 μs | 10 seconds | O(n_bootstrap) |
| Fisher's Z | 1 μs | 1 second | O(1) |
| Permutation (1000x) | 10 ms | 10,000 seconds | O(n_perm) |
| Bayesian MCMC | 1 ms | 1,000 seconds | O(n_samples) |

**Recommendation**: **Bootstrap Difference** for initial analysis, **Fisher's Z** for screening

---

## Hypothesis Testing Framework

### Null Hypothesis Options

**Option A: No Difference in Correlation Strength**
```
H0: r_low = r_high  (correlations are equal)
H1: r_low ≠ r_high  (correlations differ)
```
**Use this if**: You care about any change in correlation magnitude

**Option B: No Difference in Absolute Correlation**
```
H0: |r_low| = |r_high|  (correlation strengths equal)
H1: |r_low| ≠ |r_high|  (correlation strengths differ)
```
**Use this if**: Sign doesn't matter (e.g., r=0.8 vs r=-0.8 both "strong")

**Option C: Gain/Loss of Correlation**
```
H0: (|r_low| < threshold AND |r_high| < threshold) OR (|r_low| ≥ threshold AND |r_high| ≥ threshold)
H1: Edge present in one group but not the other
```
**Use this if**: You want edges that "appear" or "disappear" (most biologically interesting)

---

## Multiple Testing Correction

With ~1M edges, you MUST correct for multiple testing:

### 1. False Discovery Rate (FDR) - RECOMMENDED
```python
from statsmodels.stats.multitest import multipletests
reject, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
```
**Controls**: Proportion of false discoveries among all discoveries
**Threshold**: Typically FDR < 0.05 or 0.10

### 2. Bonferroni Correction - TOO CONSERVATIVE
```python
alpha_corrected = 0.05 / n_tests
```
**Controls**: Family-wise error rate
**Problem**: Too stringent for exploratory biology (almost no discoveries)

### 3. Permutation-based FDR
```python
# Permute sample labels, recalculate all p-values, estimate null distribution
# Time-consuming but most accurate
```

---

## Practical Recommendations

### **Two-Stage Strategy (RECOMMENDED)**

**Stage 1: Fast Screening (Fisher's Z)**
- Use base correlations only
- Apply Fisher's Z test to all edges
- Threshold: p < 0.001 (lenient, for screening)
- **Result**: ~1,000-10,000 candidate edges

**Stage 2: Robust Testing (Bootstrap Difference)**
- Use full bootstrap data on candidates only
- Apply bootstrap difference test
- Apply FDR correction (BH method)
- Threshold: FDR < 0.05
- **Result**: ~100-1,000 significant edges

**Benefits**:
✅ Fast initial screening (1 second)
✅ Rigorous testing on candidates (10 seconds)
✅ Leverages bootstrap data where it matters
✅ Computational feasible for 1M edges

---

## Example Implementation

```python
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def analyze_differential_coexpression(
    low_base, low_boot,    # (n_edges,), (n_edges, n_bootstrap)
    high_base, high_boot,  # (n_edges,), (n_edges, n_bootstrap)
    n_low, n_high,         # sample sizes
    fdr_threshold=0.05
):
    """
    Two-stage analysis: Fisher's Z screening + bootstrap testing.

    Returns
    -------
    results : dict
        - candidate_idx: indices passing screening
        - pvals: bootstrap p-values for candidates
        - delta_mean: mean differences
        - ci_low, ci_high: confidence intervals
        - significant: boolean array (FDR < threshold)
    """
    n_edges = len(low_base)

    # ===== Stage 1: Fisher's Z screening =====
    z_low = 0.5 * np.log((1 + low_base) / (1 - low_base))
    z_high = 0.5 * np.log((1 + high_base) / (1 - high_base))
    se_diff = np.sqrt(1/(n_low - 3) + 1/(n_high - 3))
    z_score = (z_high - z_low) / se_diff
    pvals_screen = 2 * (1 - stats.norm.cdf(np.abs(z_score)))

    # Select candidates (lenient threshold for screening)
    candidate_idx = np.where(pvals_screen < 0.001)[0]
    print(f"Stage 1: {len(candidate_idx)} / {n_edges} candidates")

    # ===== Stage 2: Bootstrap testing on candidates =====
    n_candidates = len(candidate_idx)
    pvals_boot = np.zeros(n_candidates)
    delta_mean = np.zeros(n_candidates)
    ci_low = np.zeros(n_candidates)
    ci_high = np.zeros(n_candidates)

    for i, edge_idx in enumerate(candidate_idx):
        delta_boot = high_boot[edge_idx] - low_boot[edge_idx]
        delta_mean[i] = np.mean(delta_boot)
        ci_low[i], ci_high[i] = np.percentile(delta_boot, [2.5, 97.5])

        # Two-tailed p-value
        p_pos = np.mean(delta_boot > 0)
        p_neg = np.mean(delta_boot < 0)
        pvals_boot[i] = 2 * min(p_pos, p_neg)

    # FDR correction
    reject, pvals_corrected, _, _ = multipletests(
        pvals_boot, alpha=fdr_threshold, method='fdr_bh'
    )

    print(f"Stage 2: {np.sum(reject)} / {n_candidates} significant (FDR < {fdr_threshold})")

    return {
        'candidate_idx': candidate_idx,
        'pvals': pvals_boot,
        'pvals_corrected': pvals_corrected,
        'delta_mean': delta_mean,
        'ci_low': ci_low,
        'ci_high': ci_high,
        'significant': reject,
        'fdr_threshold': fdr_threshold
    }
```

---

## Biological Interpretation

### Effect Size Thresholds

**Weak difference**: |Δr| < 0.1
**Moderate difference**: 0.1 ≤ |Δr| < 0.3
**Strong difference**: |Δr| ≥ 0.3

**Recommendation**: Filter by both statistical significance AND effect size:
```python
# Significant AND moderate effect
significant_edges = (results['significant']) & (np.abs(results['delta_mean']) > 0.1)
```

### Direction of Change

**Gain of correlation**: r_low ≈ 0, r_high >> 0 (or << 0)
**Loss of correlation**: r_low >> 0, r_high ≈ 0
**Sign flip**: r_low > 0, r_high < 0 (or vice versa) - RARE but very interesting!
**Magnitude change**: Both significant, but different strengths

---

## Final Recommendation

**Default Pipeline**:
1. **Screen** with Fisher's Z (all edges, 1 second)
2. **Test** with Bootstrap Difference (top candidates, 10 seconds)
3. **Correct** with FDR BH (q < 0.05)
4. **Filter** by effect size (|Δr| > 0.1)

**Validation** (optional):
- Run permutation test on top 100 edges
- Compare with Fisher's Z + Bootstrap results
- Should agree ~95%+

**Output**: Ranked list of edges by:
- FDR-corrected p-value (significance)
- Mean difference |Δr| (effect size)
- Confidence interval width (precision)

---

# Part II: Efficient Pipeline Design and Analysis

## Questions Answered

### Q1: Fisher's Z vs t-distribution for testing single correlations

**Both test the same null hypothesis (r = 0), but with different approaches:**

#### Fisher's Z-Transform Test
```python
# Test H0: r = 0
z = 0.5 * np.log((1 + r) / (1 - r))
SE = 1 / np.sqrt(n - 3)
z_score = z / SE
p_value = 2 * (1 - norm.cdf(abs(z_score)))
```

**Assumptions:**
- Large sample (n > 30 recommended)
- Correlation approximately normally distributed after transform
- Works better for r far from 0

#### t-Distribution Test
```python
# Test H0: r = 0
t = r * np.sqrt(n - 2) / np.sqrt(1 - r**2)
p_value = 2 * (1 - t_dist.cdf(abs(t), df=n-2))
```

**Assumptions:**
- Data are bivariate normal
- Exact test for small samples
- Standard approach in most statistics texts

**Which to use?**
- **t-test**: More common, exact for normally distributed data
- **Fisher's Z**: Better for comparing TWO correlations (your Hypothesis 3)

**Recommendation**: Use **t-test** for single correlation significance (Hypotheses 1 & 2), **Fisher's Z** for comparing correlations (Hypothesis 3)

---

### Q2: Is it common to use bootstrap after Fisher's Z?

**Yes, absolutely!** This is a standard two-stage approach:

**Stage 1 (Fisher's Z)**: Fast parametric screening
- Assumes normality after Z-transform
- Very fast: O(1) per edge
- Good for reducing search space

**Stage 2 (Bootstrap)**: Robust non-parametric confirmation
- No distributional assumptions
- Provides confidence intervals
- Gold standard for final inference

**Why this is good practice:**
1. Fisher's Z is **optimistic** (assumes perfect normality)
2. Bootstrap is **conservative** (makes no assumptions)
3. Edges passing both tests are very robust
4. Computationally efficient: test millions → confirm thousands

**Similar workflows in genomics:**
- GWAS: Chi-square screening → permutation testing
- RNA-seq: t-test screening → exact tests on DE genes
- ChIP-seq: Poisson screening → FDR on peaks

---

## Proposed Efficient Pipeline

### Original Pipeline (Inefficient)

```
Stage 1: Generate indices (1 job, ~30s)
    → bootstrap_indices.h5

Stage 2: Compute base + bootstrap for ALL edges (20,000 jobs)
    For each gene:
        low_base:  Compute correlations → (200M edges)
        high_base: Compute correlations → (200M edges)
        low_boot:  50 bootstrap replicates → (50 × 200M)
        high_boot: 50 bootstrap replicates → (50 × 200M)

    Total computations: 20,000 × (2 + 100) correlations = 2,040,000 full correlation matrices
    Time per gene: ~10 min
    Total time: 20,000 × 10 min = 138 days (with parallelization: hours)
    Storage: 20,000 × 40 GB = 800 TB
```

**Problems:**
❌ Computes bootstrap for ALL 200M edges (most non-significant)
❌ Massive storage (800 TB)
❌ Stage 3 must process 200M × 20,000 genes

---

### Proposed Efficient Pipeline

```
Stage 1: Generate indices (1 job, ~30s)
    → bootstrap_indices.h5

Stage 2a: Compute BASE correlations + significance (20,000 jobs OR 1 merged job)
    For each gene:
        low_base:  Compute correlations → (200M edges)
        high_base: Compute correlations → (200M edges)

        # Test single correlations
        low_pval:  t-test for r_low ≠ 0 → (200M p-values)
        high_pval: t-test for r_high ≠ 0 → (200M p-values)

        # Test differential
        diff_z:    Fisher's Z for r_low ≠ r_high → (200M Z-scores)
        diff_pval: p-value from Z-score → (200M p-values)

    Apply FDR correction:
        low_qval:  FDR on low_pval → (200M q-values)
        high_qval: FDR on high_pval → (200M q-values)
        diff_qval: FDR on diff_pval → (200M q-values)

    Output: base_correlations.h5
        low/corr_triu      (n_tests,)
        low/pval           (n_tests,)
        low/qval           (n_tests,)
        high/corr_triu     (n_tests,)
        high/pval          (n_tests,)
        high/qval          (n_tests,)
        diff/fisher_z      (n_tests,)
        diff/pval          (n_tests,)
        diff/qval          (n_tests,)

    Time per gene: ~2 min (no bootstrap!)
    Total time: 20,000 × 2 min = 27 days (parallelized: ~1 hour)
    Storage: 20,000 × 7.2 GB = 144 GB (much better!)

Stage 2b: Bootstrap ONLY significant differential edges
    # Select edges where diff_qval < 0.05
    # Typical: ~1,000-10,000 edges out of 200M

    For each significant edge (can parallelize):
        For bootstrap rep in 1..50:
            low_boot_corr[edge]:  Compute correlation on low_boot samples
            high_boot_corr[edge]: Compute correlation on high_boot samples

    Output: bootstrap_differential.h5
        edges/indices         (n_sig, 2) - gene pair indices
        edges/low_boot        (n_sig, 50) - bootstrap correlations
        edges/high_boot       (n_sig, 50) - bootstrap correlations
        edges/delta_boot      (n_sig, 50) - high - low
        edges/ci_low          (n_sig,) - 95% CI lower
        edges/ci_high         (n_sig,) - 95% CI upper

    Time: 10,000 edges × 50 bootstrap × 2 subpops × 0.001s = 1000s (~15 min)
    Storage: 10,000 × 50 × 2 × 4 bytes = 4 MB (tiny!)

Stage 3: Final analysis and matrix reconstruction
    # Load significant edges from Stage 2a
    # Load bootstrap CI from Stage 2b
    # Reconstruct sparse correlation matrices
```

**Benefits:**
✅ Computes bootstrap only for significant edges (~0.005% of edges)
✅ Storage: 144 GB (vs 800 TB = 5,500× reduction)
✅ Stage 2b runs in 15 minutes (vs days of bootstrap computation)
✅ Clear separation: significance → bootstrap confirmation

---

## Efficiency Comparison: Original vs Proposed

| Metric | Original Pipeline | Proposed Pipeline | Speedup |
|--------|------------------|-------------------|---------|
| **Bootstrap computations** | 200M × 50 × 2 = 20B | 10K × 50 × 2 = 1M | 20,000× |
| **Storage** | 800 TB | 144 GB + 4 MB | 5,500× |
| **Stage 2 time** | Days (parallelized) | ~1 hour + 15 min | 50-100× |
| **I/O in Stage 3** | Read 800 TB | Read 144 GB | 5,500× |

---

## Selecting Significant Edges: Strategy Options

### Your Question: "Should I require both significant OR either one?"

This depends on your biological question:

#### Option A: Union (Either Low OR High Significant)
```python
sig_edges = (low_qval < 0.05) | (high_qval < 0.05)
```

**Biological interpretation:**
"Edges that are functionally relevant in AT LEAST ONE subpopulation"

**Use when:**
- Interested in edges that "emerge" or "disappear"
- Example: Edge present in high expression but not low

**Risk:** More edges, some may be false positives

---

#### Option B: Intersection (Both Low AND High Significant)
```python
sig_edges = (low_qval < 0.05) & (high_qval < 0.05)
```

**Biological interpretation:**
"Edges that are robustly present in BOTH subpopulations"

**Use when:**
- Interested in context-independent core interactions
- Testing if STRENGTH differs (not presence/absence)

**Risk:** Miss interesting differential edges that are strong in one condition

---

#### Option C: Differential Priority (Significant Difference, Regardless of Individual Significance)
```python
sig_edges = (diff_qval < 0.05)
# Don't filter by low_qval or high_qval
```

**Biological interpretation:**
"Edges where correlation CHANGES between conditions"

**Use when:**
- Interested in context-dependent rewiring
- Don't care if edge is weak in both, just that it differs

**Risk:** May include edges that are weak in both (low power)

---

#### Option D: Hybrid (Recommended for Most Analyses)
```python
# Require: significant in at least one + significant difference
sig_individual = (low_qval < 0.05) | (high_qval < 0.05)
sig_differential = (diff_qval < 0.05)
sig_edges = sig_individual & sig_differential
```

**Biological interpretation:**
"Edges that are functionally relevant AND change between conditions"

**Use when:**
- Want to avoid testing noise (both correlations near 0)
- Focus on biologically meaningful rewiring

**This is the most common approach in differential co-expression studies.**

---

## Reconstructing Correlation Matrices from Significant Edges

### Sparse Matrix Reconstruction

```python
import scipy.sparse as sp

def reconstruct_sparse_correlation_matrix(
    edge_indices: np.ndarray,  # (n_sig, 2) - gene pairs
    edge_values: np.ndarray,   # (n_sig,) - correlations
    n_genes: int
) -> sp.csr_matrix:
    """
    Reconstruct sparse correlation matrix from significant edges.

    Returns
    -------
    corr_matrix : (n_genes, n_genes) sparse matrix
        Symmetric correlation matrix with:
        - Diagonal = 1.0
        - Significant edges = correlation value
        - Non-significant edges = 0 (implicit)
    """
    # Initialize with diagonal
    row = []
    col = []
    data = []

    # Add diagonal (self-correlations = 1)
    row.extend(range(n_genes))
    col.extend(range(n_genes))
    data.extend([1.0] * n_genes)

    # Add significant edges (symmetric)
    for (i, j), val in zip(edge_indices, edge_values):
        row.extend([i, j])
        col.extend([j, i])
        data.extend([val, val])

    # Create sparse matrix
    corr_matrix = sp.csr_matrix(
        (data, (row, col)),
        shape=(n_genes, n_genes)
    )

    return corr_matrix

# Usage
low_corr_sparse = reconstruct_sparse_correlation_matrix(
    sig_edge_indices, low_corr_values, n_genes
)
high_corr_sparse = reconstruct_sparse_correlation_matrix(
    sig_edge_indices, high_corr_values, n_genes
)
```

**Storage:**
- Dense matrix: 20,000 × 20,000 × 4 bytes = 1.6 GB
- Sparse matrix (10K edges): 10,000 × 3 × 4 bytes = 120 KB (13,000× smaller)

---

## Recommended Pipeline Implementation

### Stage 2a: Base Correlations + Significance

```python
#!/usr/bin/env python3
"""
Stage 2a: Compute base correlations and test significance.

For each gene's low/high subpopulations:
1. Compute full correlation matrix
2. Test H0: r = 0 (t-test)
3. Test H0: r_low = r_high (Fisher's Z)
4. Apply FDR correction
"""

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def test_correlation_significance(r: np.ndarray, n: int) -> np.ndarray:
    """
    Test H0: correlation = 0 using t-distribution.

    Parameters
    ----------
    r : array of correlations
    n : sample size

    Returns
    -------
    pval : array of p-values
    """
    # Handle edge cases (r = ±1)
    r_clip = np.clip(r, -0.9999, 0.9999)

    # t-statistic: t = r * sqrt(n - 2) / sqrt(1 - r^2)
    t_stat = r_clip * np.sqrt(n - 2) / np.sqrt(1 - r_clip**2)

    # Two-tailed p-value
    pval = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=n-2))

    return pval


def test_correlation_difference(
    r_low: np.ndarray,
    r_high: np.ndarray,
    n_low: int,
    n_high: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Test H0: r_low = r_high using Fisher's Z-transform.

    Returns
    -------
    z_score : array of Z-scores
    pval : array of p-values
    """
    # Fisher's Z transform
    r_low_clip = np.clip(r_low, -0.9999, 0.9999)
    r_high_clip = np.clip(r_high, -0.9999, 0.9999)

    z_low = 0.5 * np.log((1 + r_low_clip) / (1 - r_low_clip))
    z_high = 0.5 * np.log((1 + r_high_clip) / (1 - r_high_clip))

    # Standard error of difference
    se_diff = np.sqrt(1.0 / (n_low - 3) + 1.0 / (n_high - 3))

    # Z-score
    z_score = (z_high - z_low) / se_diff

    # Two-tailed p-value
    pval = 2 * (1 - stats.norm.cdf(np.abs(z_score)))

    return z_score, pval


def compute_base_with_significance(
    expr: np.ndarray,  # (n_genes, n_samples)
    low_indices: np.ndarray,
    high_indices: np.ndarray,
    fdr_alpha: float = 0.05
) -> dict:
    """
    Compute base correlations and test significance.

    Returns
    -------
    results : dict with keys:
        low/corr_triu, low/pval, low/qval
        high/corr_triu, high/pval, high/qval
        diff/z_score, diff/pval, diff/qval
    """
    n_genes = expr.shape[0]

    # Compute correlations
    expr_low = expr[:, low_indices]
    expr_high = expr[:, high_indices]

    corr_low = np.corrcoef(expr_low)
    corr_high = np.corrcoef(expr_high)

    # Extract upper triangle
    triu_idx = np.triu_indices(n_genes, k=1)
    corr_low_triu = corr_low[triu_idx]
    corr_high_triu = corr_high[triu_idx]

    # Test individual correlations
    low_pval = test_correlation_significance(corr_low_triu, len(low_indices))
    high_pval = test_correlation_significance(corr_high_triu, len(high_indices))

    # Test differential
    diff_z, diff_pval = test_correlation_difference(
        corr_low_triu, corr_high_triu,
        len(low_indices), len(high_indices)
    )

    # FDR correction
    _, low_qval, _, _ = multipletests(low_pval, alpha=fdr_alpha, method='fdr_bh')
    _, high_qval, _, _ = multipletests(high_pval, alpha=fdr_alpha, method='fdr_bh')
    _, diff_qval, _, _ = multipletests(diff_pval, alpha=fdr_alpha, method='fdr_bh')

    return {
        'low/corr_triu': corr_low_triu.astype(np.float32),
        'low/pval': low_pval.astype(np.float32),
        'low/qval': low_qval.astype(np.float32),
        'high/corr_triu': corr_high_triu.astype(np.float32),
        'high/pval': high_pval.astype(np.float32),
        'high/qval': high_qval.astype(np.float32),
        'diff/z_score': diff_z.astype(np.float32),
        'diff/pval': diff_pval.astype(np.float32),
        'diff/qval': diff_qval.astype(np.float32),
    }
```

---

### Stage 2b: Bootstrap for Significant Edges Only

```python
def bootstrap_significant_edges(
    expr: np.ndarray,
    edge_indices: np.ndarray,  # (n_sig, 2) - gene pairs
    low_boot_indices: np.ndarray,  # (n_bootstrap, k_low)
    high_boot_indices: np.ndarray,  # (n_bootstrap, k_high)
) -> dict:
    """
    Compute bootstrap correlations for significant edges only.

    Returns
    -------
    results : dict with keys:
        edges/low_boot  (n_sig, n_bootstrap)
        edges/high_boot (n_sig, n_bootstrap)
        edges/delta_boot (n_sig, n_bootstrap)
        edges/ci_low, edges/ci_high (n_sig,)
    """
    n_sig = len(edge_indices)
    n_bootstrap = len(low_boot_indices)

    low_boot = np.zeros((n_sig, n_bootstrap), dtype=np.float32)
    high_boot = np.zeros((n_sig, n_bootstrap), dtype=np.float32)

    for b in range(n_bootstrap):
        # Subset expression
        expr_low_b = expr[:, low_boot_indices[b]]
        expr_high_b = expr[:, high_boot_indices[b]]

        # Compute correlations for significant edges only
        for i, (gene_i, gene_j) in enumerate(edge_indices):
            # Low bootstrap
            x_low = expr_low_b[gene_i, :]
            y_low = expr_low_b[gene_j, :]
            low_boot[i, b] = np.corrcoef(x_low, y_low)[0, 1]

            # High bootstrap
            x_high = expr_high_b[gene_i, :]
            y_high = expr_high_b[gene_j, :]
            high_boot[i, b] = np.corrcoef(x_high, y_high)[0, 1]

    # Compute delta and CI
    delta_boot = high_boot - low_boot
    ci_low = np.percentile(delta_boot, 2.5, axis=1)
    ci_high = np.percentile(delta_boot, 97.5, axis=1)

    return {
        'edges/low_boot': low_boot,
        'edges/high_boot': high_boot,
        'edges/delta_boot': delta_boot,
        'edges/ci_low': ci_low.astype(np.float32),
        'edges/ci_high': ci_high.astype(np.float32),
    }
```

---

## Summary & Recommendations

### Statistical Tests Summary

| Test | Purpose | Method | Use For |
|------|---------|--------|---------|
| **t-test** | r ≠ 0 | `t = r√(n-2)/√(1-r²)` | Single correlation significance |
| **Fisher's Z** | r₁ ≠ r₂ | `Z = (Z₁ - Z₂)/SE` | Comparing two correlations |
| **Bootstrap** | Δr ≠ 0 | Percentile CI | Robust confirmation |

### Pipeline Recommendation

**Use the proposed efficient pipeline:**

1. **Stage 2a**: Base + significance (fast, all edges)
   - Compute low/high correlations
   - Test r ≠ 0 (t-test)
   - Test r_low ≠ r_high (Fisher's Z)
   - FDR correction

2. **Stage 2b**: Bootstrap confirmation (slow, but only ~10K edges)
   - Bootstrap only edges where `diff_qval < 0.05`
   - Get confidence intervals

3. **Stage 3**: Select edges using hybrid filter:
   ```python
   sig_individual = (low_qval < 0.05) | (high_qval < 0.05)
   sig_differential = (diff_qval < 0.05)
   final_edges = sig_individual & sig_differential & (ci_excludes_zero)
   ```

### Efficiency Gains

- **Storage**: 800 TB → 144 GB (5,500× reduction)
- **Computation**: 20 billion bootstraps → 1 million (20,000× reduction)
- **Time**: Days → ~1.5 hours total

This is a **much better pipeline design**!

---

## IMPLEMENTED PIPELINE (Final)

The following scripts have been implemented based on the analysis above:

### Scripts

| Stage | Script | Output |
|-------|--------|--------|
| 0 | `00preprocess/00convert_expr_to_hdf5.py` | `expression.h5` |
| 1 | `01subset/01get_extreme_pop_bootstrap.py` | `bootstrap_indices.h5` |
| 2a | `10spearman_corr/02a_calc_base_correlations.py` | `base_correlations.h5` |
| 2b | `10spearman_corr/02b_bootstrap_significant_edges.py` | `bootstrap_significant.h5` |
| 3 | `10spearman_corr/03_reconstruct_diff_network.py` | `differential_network.h5` |
| 4 | `10spearman_corr/04_collect_focus_gene_topology.py` | `focus_gene_topology.h5` |
| 5 | `10spearman_corr/05_prepare_visualization_data.py` | `visualization_data/` |

### Quick Test (Toy Data)

```bash
# Run complete pipeline with toy data (5 genes, 50 samples)
bash src/SGE_scripts/run_bootstrap_pipeline.sh --toy --local --out-dir results_test
```

### Full Usage

```bash
# Stage 0: Convert TSV to HDF5 (optional, for faster loading)
python src/scripts/00preprocess/00convert_expr_to_hdf5.py \
    --expr-tsv data/expression.tsv \
    --out-h5 data/expression.h5

# Stage 1: Generate bootstrap indices
python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
    --in-h5 data/expression.h5 \
    --out-h5 results/bootstrap_indices.h5 \
    --low-frac 0.2 --high-frac 0.2 \
    --n-bootstrap 50

# Stage 2a: Compute base correlations + significance tests
python src/scripts/10spearman_corr/02a_calc_base_correlations.py \
    --expr-h5 data/expression.h5 \
    --indices-h5 results/bootstrap_indices.h5 \
    --out-h5 results/base_correlations.h5 \
    --fdr-alpha 0.05

# Stage 2b: Bootstrap only significant edges
python src/scripts/10spearman_corr/02b_bootstrap_significant_edges.py \
    --expr-h5 data/expression.h5 \
    --indices-h5 results/bootstrap_indices.h5 \
    --base-h5 results/base_correlations.h5 \
    --out-h5 results/bootstrap_significant.h5 \
    --edge-selection sig_edges

# Stage 3: Reconstruct differential network
python src/scripts/10spearman_corr/03_reconstruct_diff_network.py \
    --base-h5 results/base_correlations.h5 \
    --boot-h5 results/bootstrap_significant.h5 \
    --out-h5 results/differential_network.h5 \
    --edge-selection sig_edges

# Stage 5: Prepare visualization data
python src/scripts/10spearman_corr/05_prepare_visualization_data.py \
    --diff-h5 results/differential_network.h5 \
    --out-dir results/visualization_data \
    --top-n 10

# Stage 6: Run R visualization
cd results/visualization_data && Rscript visualize_networks.R
```

### Edge Selection Modes

The `--edge-selection` parameter controls which edges are analyzed:

| Mode | Formula | Description |
|------|---------|-------------|
| `sig_edges` | `(low_qval < 0.05 \| high_qval < 0.05) & (diff_qval < 0.05)` | Edges significant in at least one condition AND significantly different |
| `sig_differential` | `diff_qval < 0.05` | Edges with significant differential correlation only |

### Output Files

**`base_correlations.h5`** (Stage 2a):
```
meta/           - n_genes, k_low, k_high, fdr_alpha
low/corr_triu   - (n_tests,) low subpop correlations
low/pval_triu   - (n_tests,) p-values from t-test
low/qval_triu   - (n_tests,) FDR-corrected q-values
high/corr_triu  - (n_tests,) high subpop correlations
high/pval_triu  - (n_tests,) p-values
high/qval_triu  - (n_tests,) q-values
diff/fisher_z   - (n_tests,) Fisher's Z statistic
diff/pval_triu  - (n_tests,) p-values
diff/qval_triu  - (n_tests,) q-values
diff/delta_triu - (n_tests,) r_high - r_low
significant/    - boolean masks and indices
```

**`bootstrap_significant.h5`** (Stage 2b):
```
meta/           - n_sig_edges, n_bootstrap, edge_selection_mode
edges/indices   - (n_sig,) flat triu indices
edges/gene_i    - (n_sig,) gene i for each edge
edges/gene_j    - (n_sig,) gene j for each edge
base/delta      - (n_sig,) original delta (high - low)
base/r_low      - (n_sig,) original low correlation
base/r_high     - (n_sig,) original high correlation
boot/delta      - (n_sig, n_bootstrap) bootstrap deltas
boot/delta_mean - (n_sig,) mean bootstrap delta
boot/delta_std  - (n_sig,) std of bootstrap delta
boot/ci_low     - (n_sig,) 95% CI lower bound
boot/ci_high    - (n_sig,) 95% CI upper bound
boot/bias       - (n_sig,) delta_mean - delta_base (bootstrap bias)
pval/bootstrap_pval - (n_sig,) bootstrap p-value
```

**`differential_network.h5`** (Stage 3):
```
meta/           - n_genes, n_significant, edge_selection, min_effect
edges/          - filtered edge data (gene_i, gene_j, delta, CI, p-values)
matrices/       - sparse matrix representation for reconstruction
```

### Key Features

1. **Three significance tests**:
   - Individual low correlation (t-test)
   - Individual high correlation (t-test)
   - Differential correlation (Fisher's Z)

2. **Bootstrap only for significant edges** (~100× speedup)

3. **Bootstrap bias estimation**: `bias = delta_boot_mean - delta_base`

4. **Sparse matrix reconstruction** for downstream network analysis

---

## Related Reading

- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) - Full pipeline workflow guide
- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) - Network topology metrics
- [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) - Memory optimization strategies
- [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) - Storage optimization guide
- [FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md) - Summary of critical fixes

---

**Last Updated:** 2026-02-13
**Status:** ✅ Complete technical reference
