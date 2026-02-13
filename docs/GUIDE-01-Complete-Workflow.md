# Complete Differential Co-expression Pipeline Workflow

## Overview

This pipeline identifies **differential co-expression** between low and high expression subpopulations, using bootstrap confidence intervals for robust statistical inference.

---

## Pipeline Stages

| Stage | Script | Output | Description |
|-------|--------|--------|-------------|
| 0 | `00preprocess/00convert_expr_to_hdf5.py` | `expression.h5` | Convert TSV to HDF5 |
| 1 | `01subset/01get_extreme_pop_bootstrap.py` | `bootstrap_indices.h5` | Generate bootstrap sample indices |
| 2a | `10spearman_corr/02a_calc_base_correlations.py` | `base_correlations.h5` | Compute correlations + significance tests |
| 2b | `10spearman_corr/02b_bootstrap_significant_edges.py` | `bootstrap_significant.h5` | Bootstrap only significant edges |
| 3 | `10spearman_corr/03_reconstruct_diff_network.py` | `differential_network.h5` | Reconstruct network with topology |
| 4 | `10spearman_corr/04_collect_focus_gene_topology.py` | `focus_gene_topology.h5` | (Optional) Aggregate per-gene networks |
| 5 | `10spearman_corr/05_prepare_visualization_data.py` | `visualization_data/` | Export for R visualization |

---

## Quick Start

### Option 1: Local Test with Toy Data
```bash
# Run complete pipeline locally (no SGE needed)
bash src/SGE_scripts/run_bootstrap_pipeline.sh --toy --local --out-dir results_toy

# Outputs in results_toy/:
#   bootstrap_indices.h5
#   base_correlations.h5
#   bootstrap_significant.h5
#   differential_network.h5
#   visualization_data/
```

### Option 2: Real Data with SGE
```bash
# Convert TSV to HDF5 and run full pipeline
bash src/SGE_scripts/run_bootstrap_pipeline.sh \
    --expr-tsv dataset/voomdataCtrl.txt \
    --out-dir results \
    --n-bootstrap 50 \
    --seed 42

# Or with existing HDF5
bash src/SGE_scripts/run_bootstrap_pipeline.sh \
    --expr-h5 dataset/expression.h5 \
    --out-dir results
```

---

## Detailed Stage Commands

### Stage 0: Convert Expression TSV to HDF5 (Optional)
```bash
python src/scripts/00preprocess/00convert_expr_to_hdf5.py \
    --expr-tsv data/expression.tsv \
    --out-h5 data/expression.h5 \
    --compression gzip
```

### Stage 1: Generate Bootstrap Indices
```bash
# From HDF5
python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
    --in-h5 data/expression.h5 \
    --out-h5 results/bootstrap_indices.h5 \
    --low-frac 0.2 --high-frac 0.2 \
    --n-bootstrap 50 --bootstrap-frac 0.8 \
    --seed 42

# Or with toy data
python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
    --toy \
    --out-h5 results/bootstrap_indices.h5
```

### Stage 2a: Compute Base Correlations + Significance Tests
```bash
python src/scripts/10spearman_corr/02a_calc_base_correlations.py \
    --expr-h5 data/expression.h5 \
    --indices-h5 results/bootstrap_indices.h5 \
    --out-h5 results/base_correlations.h5 \
    --fdr-alpha 0.05

# Or with toy data
python src/scripts/10spearman_corr/02a_calc_base_correlations.py \
    --toy \
    --indices-h5 results/bootstrap_indices.h5 \
    --out-h5 results/base_correlations.h5
```

**Significance tests performed:**
- Individual correlation significance (t-test for r ≠ 0)
- Differential correlation significance (Fisher's Z for r_low ≠ r_high)
- FDR correction (Benjamini-Hochberg)

### Stage 2b: Bootstrap Significant Edges Only
```bash
python src/scripts/10spearman_corr/02b_bootstrap_significant_edges.py \
    --expr-h5 data/expression.h5 \
    --indices-h5 results/bootstrap_indices.h5 \
    --base-h5 results/base_correlations.h5 \
    --out-h5 results/bootstrap_significant.h5 \
    --edge-selection sig_edges

# Or with toy data
python src/scripts/10spearman_corr/02b_bootstrap_significant_edges.py \
    --toy \
    --indices-h5 results/bootstrap_indices.h5 \
    --base-h5 results/base_correlations.h5 \
    --out-h5 results/bootstrap_significant.h5
```

**Key insight:** Only bootstraps significant edges (~1% of all edges), providing ~100x speedup.

### Stage 3: Reconstruct Differential Network
```bash
python src/scripts/10spearman_corr/03_reconstruct_diff_network.py \
    --base-h5 results/base_correlations.h5 \
    --boot-h5 results/bootstrap_significant.h5 \
    --out-h5 results/differential_network.h5 \
    --edge-selection sig_edges

# With per-gene metrics and TSV output
python src/scripts/10spearman_corr/03_reconstruct_diff_network.py \
    --base-h5 results/base_correlations.h5 \
    --boot-h5 results/bootstrap_significant.h5 \
    --out-h5 results/differential_network.h5 \
    --out-tsv results/rewiring_hubs.tsv \
    --calc-per-gene-metrics
```

**Features:**
- Global topology for LOW, HIGH, and DIFF networks
- Qualitative edge classification (disappear, new, sign_change, strengthen, weaken)
- Focus gene neighborhood analysis (1st and 2nd layer partners)
- Optional per-gene metrics (--calc-per-gene-metrics)

### Stage 4: Collect Focus Gene Topology (Optional)
```bash
# Only needed for per-gene network analysis
python src/scripts/10spearman_corr/04_collect_focus_gene_topology.py \
    --network-dir results/networks \
    --focus-genes top:50 \
    --n-genes 20000 \
    --out-h5 results/focus_gene_topology.h5 \
    --n-jobs 8
```

### Stage 5: Prepare Visualization Data
```bash
python src/scripts/10spearman_corr/05_prepare_visualization_data.py \
    --diff-h5 results/differential_network.h5 \
    --out-dir results/visualization_data \
    --top-n 10
```

**Outputs:**
```
visualization_data/
├── edges_low.tsv          - Low network edges
├── edges_high.tsv         - High network edges
├── edges_diff.tsv         - Differential edges
├── nodes.tsv              - Node attributes (3 degree types)
├── ego_networks/          - GraphML files for Cytoscape
│   ├── ego_gene_0.graphml
│   └── ...
└── visualize_networks.R   - R script template
```

### Stage 6: R Visualization
```bash
cd results/visualization_data
Rscript visualize_networks.R
```

---

## Three Network Types

For each analysis, you get **3 networks**:

1. **Low network** (`corr_low`, `qval_low`)
   - Co-expression in LOW expression subpopulation
   - Degree = hub status in baseline state

2. **High network** (`corr_high`, `qval_high`)
   - Co-expression in HIGH expression subpopulation
   - Degree = hub status in altered state

3. **Differential network** (`delta = corr_high - corr_low`, `qval_diff`)
   - CHANGES in co-expression
   - Degree = rewiring magnitude

**Biological interpretation:**
```
Gene A:
  - degree_low = 5   (minor player in baseline)
  - degree_high = 50  (becomes hub in high condition)
  - degree_diff = 45  (massive rewiring)

Conclusion: Gene A is a condition-specific regulatory hub
```

---

## Edge Selection Modes

The `--edge-selection` parameter controls which edges are analyzed:

| Mode | Formula | Description |
|------|---------|-------------|
| `sig_edges` | `(low_qval < α \| high_qval < α) & (diff_qval < α)` | Edges significant in at least one condition AND significantly different |
| `sig_differential` | `diff_qval < α` | Edges with significant differential correlation only |

**Recommendation:** Use `sig_edges` (default) to focus on biologically meaningful rewiring.

---

## Qualitative Edge Classification

Edges are classified based on how they change from low to high:

| Category | Definition | Biological Meaning |
|----------|------------|-------------------|
| `disappear` | Present in low, absent in high | Edge lost in high expression |
| `new` | Absent in low, present in high | Edge gained in high expression |
| `sign_change` | Sign flips (+ → - or - → +) | Regulatory inversion |
| `strengthen` | Same sign, \|r_high\| > \|r_low\| | Edge strengthens |
| `weaken` | Same sign, \|r_high\| < \|r_low\| | Edge weakens |
| `unchanged` | No significant change | Stable edge |

---

## HDF5 Output Schemas

### base_correlations.h5 (Stage 2a)
```
meta/
    n_genes, n_samples, k_low, k_high, fdr_alpha
low/
    corr_triu      (n_tests,) - Spearman correlations
    pval_triu      (n_tests,) - p-values from t-test
    qval_triu      (n_tests,) - FDR-corrected q-values
high/
    corr_triu, pval_triu, qval_triu
diff/
    fisher_z       (n_tests,) - Fisher's Z statistic
    pval_triu, qval_triu
    delta_triu     (n_tests,) - r_high - r_low
significant/
    sig_low, sig_high, sig_individual, sig_differential, sig_edges
    indices        (n_sig,) - flat indices of sig_edges
```

### bootstrap_significant.h5 (Stage 2b)
```
meta/
    n_sig_edges, n_bootstrap, edge_selection_mode, ci_alpha
edges/
    indices        (n_sig,) - flat triu indices
    gene_i, gene_j (n_sig,) - gene pair indices
base/
    delta, r_low, r_high
boot/
    delta          (n_sig, n_bootstrap) - bootstrap deltas
    delta_mean, delta_std, ci_low, ci_high, bias
pval/
    bootstrap_pval (n_sig,) - bootstrap p-value for Δr ≠ 0
```

### differential_network.h5 (Stage 3)
```
meta/
    n_genes, n_significant, edge_selection, corr_threshold
edges/
    gene_i, gene_j, delta_base, r_low, r_high
    qual_score, qual_label
    ci_low, ci_high, pval_boot
topology/
    global_low/    - degrees, n_edges, density, ...
    global_high/   - degrees, n_edges, density, ...
    global_diff/   - degrees, n_edges, density, ...
    per_gene/      - (optional) per-gene metrics
focus_gene/
    direct_partners, indirect_partners
    direct_stats, two_layer_stats
matrices/
    delta_data, delta_indices, delta_indptr  - CSR sparse matrix
```

---

## Efficiency Summary

| Approach | Time | Storage | Use When |
|----------|------|---------|----------|
| **Global (single reference)** | ~30 min | ~8 GB | Standard differential co-expression |
| **Per-gene networks** | ~3 hours | ~140 TB | Context-specific analysis |
| **Focus gene collection** | ~2 hours | ~4 GB | Top rewiring genes deep dive |

**Recommendation:** Start with global analysis, use per-gene only for targeted follow-up.

---

## Testing the Pipeline

```bash
# Quick smoke test (5 genes, 50 samples, runs locally)
bash src/SGE_scripts/run_bootstrap_pipeline.sh --toy --local --out-dir results_test

# Verify outputs
ls -la results_test/
# Should see:
#   bootstrap_indices.h5
#   base_correlations.h5
#   bootstrap_significant.h5
#   differential_network.h5
#   visualization_data/

# Check HDF5 contents
python -c "import h5py; h5 = h5py.File('results_test/differential_network.h5', 'r'); print(list(h5.keys()))"
```

---

## Related Reading

- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) - Understanding network topology metrics
- [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) - Memory optimization strategies (vectorized vs batch mode)
- [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) - Storage optimization (common, minimal, full modes)
- [FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md) - Summary of critical pipeline fixes
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) - Statistical methodology and pipeline design
