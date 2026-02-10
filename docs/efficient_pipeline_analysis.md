# Efficient Pipeline Design: Base Significance + Bootstrap for Differential Edges

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

## Proposed Efficient Pipeline Analysis

Your proposed pipeline is **excellent**! Let me refine it and compare efficiency.

### Original Pipeline (Current 02calc_corr_edge_bootstrap_corr.py)

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

## Comparison with 01calc_corr_edge.py

Let me check what `01calc_corr_edge.py` does (assuming it computes base correlations):

**Typical 01calc_corr_edge.py workflow:**
```python
# Compute full correlation matrix on entire dataset
expr = load_expression()  # (n_genes, n_samples)
corr_matrix = np.corrcoef(expr)  # (n_genes, n_genes)
# Extract upper triangle
corr_triu = corr_matrix[np.triu_indices(n_genes, k=1)]
# Test significance (t-test)
t_stat = r * np.sqrt(n - 2) / np.sqrt(1 - r**2)
pval = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-2))
```

**Your proposed Stage 2a does the same, but for TWO subpopulations:**
```python
# Low subpopulation
expr_low = expr[:, low_sample_indices]
corr_low = np.corrcoef(expr_low)
# ... test significance ...

# High subpopulation  
expr_high = expr[:, high_sample_indices]
corr_high = np.corrcoef(expr_high)
# ... test significance ...

# Differential test
fisher_z_low = 0.5 * np.log((1 + corr_low) / (1 - corr_low))
fisher_z_high = 0.5 * np.log((1 + corr_high) / (1 - corr_high))
se_diff = np.sqrt(1/(n_low - 3) + 1/(n_high - 3))
z_score = (fisher_z_high - fisher_z_low) / se_diff
```

**Efficiency comparison:**
- `01calc_corr_edge.py`: 1 correlation matrix + 1 significance test
- **Stage 2a**: 2 correlation matrices + 2 significance tests + 1 differential test
- **Ratio**: ~3× slower than 01calc_corr_edge.py (still very fast!)

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
