# Statistical Tests for Differential Co-Expression Analysis

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
