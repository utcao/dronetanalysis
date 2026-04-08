# Complete Differential Co-expression Pipeline Workflow

## Overview

This pipeline identifies **differential co-expression** between low and high expression subpopulations, using bootstrap confidence intervals for robust statistical inference.

---

## Pipeline Stages

| Stage | Script | Output | Description |
|-------|--------|--------|-------------|
| 0 | `00preprocess/00convert_expr_to_hdf5.py` | `expression.h5` | Convert TSV to HDF5 |
| 1 | `01subset/01get_extreme_pop_bootstrap.py` | `bootstrap_indices.h5` | Generate bootstrap sample indices |
| 2a | `10spearman_corr/02a_calc_base_correlations.py` | `base_correlations/{gi}_{gene}.h5` | Compute correlations + significance tests (per gene) |
| 2b | `10spearman_corr/02b_bootstrap_significant_edges.py` | `bootstrap_significant/{gi}_{gene}.h5` | Bootstrap only significant edges (per gene) |
| 3 | `10spearman_corr/03_reconstruct_diff_network.py` | `networks/{gi}_{gene}.h5` | Reconstruct full topology + focus gene metrics (per gene) |
| 3b | `10spearman_corr/03b_collect_networks.py` | `differential_network_summary.h5`, `rewiring_hubs.tsv` | Aggregate per-gene networks into summary |
| 4 | `10spearman_corr/04_collect_focus_gene_topology.py` | `focus_gene_topology.h5` | Collect per-gene topology for focus genes (optional) |
| 5 | `10spearman_corr/05_prepare_visualization_data.py` | `visualization_data/` | (Optional) Export for R/Cytoscape visualization |
| 6 | `10spearman_corr/06_annotate_rewiring_table.R` | `*_annotated.tsv` | Annotate gene IDs with symbols/names |
| 7 | `10spearman_corr/07_permutation_test.py` | `permutation_null/{gi}_{gene}.h5` | Permutation test: null distribution per focus gene (independent Snakefile) |
| 7b | `10spearman_corr/08_collect_permutation_pvals.py` | `permutation_pvals.tsv`, `rewiring_hubs_with_pvals.tsv` | Aggregate null distributions → empirical p-values |

---

## Quick Start

### Option 0: Snakemake (recommended)

```bash
# Dry run (shows what will execute)
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml --config out_dir=results_test -n

# Full run (15-gene subset from config)
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml --config out_dir=results_test -j 15
```

See [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) for full Snakemake
documentation including config reference, stage descriptions, and troubleshooting.

### Option 1: Local Test with Toy Data (bash script)
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

### Stage 3: Reconstruct Single-Gene Differential Network (per gene)

Run once per reference gene (Snakemake schedules these as an array job via the
`reconstruct_single` rule):

```bash
python src/scripts/10spearman_corr/03_reconstruct_diff_network.py \
    --base-h5 results/base_correlations/0001_FBgn0001234.h5 \
    --boot-h5 results/bootstrap_significant/0001_FBgn0001234.h5 \
    --out-h5 results/networks/0001_FBgn0001234.h5 \
    --edge-selection sig_differential
```

**Features:**
- Global topology for LOW, HIGH, and DIFF networks (with degree arrays)
- Qualitative edge classification (disappear, new, sign_change, strengthen, weaken)
- Focus gene neighborhood analysis (direct L1 partners + outer-layer L1↔L1 / L1→L2)
- Two-layer edge definitions:
  - `direct_edges` = focus→L1 edges only
  - `two_layer_edges` = L1↔L1 + L1→L2 edges (outer layer, excludes ego)
  - `full_two_layer_edges` = all three (complete 2-hop ego network)
- L1/L2/ratio metrics stored in `focus_gene/metrics/` for downstream collection

### Stage 3b: Collect Networks (single aggregation job)

Reads all per-gene network files from Stage 3 and aggregates them:

```bash
# Preferred: reads Stage 3 output (no re-computation)
python src/scripts/10spearman_corr/03b_collect_networks.py \
    --networks-dir results/networks \
    --out-h5 results/differential_network_summary.h5 \
    --out-focus-tsv results/rewiring_hubs.tsv \
    --annotate

# Fallback: reads raw base/boot dirs (no Stage 3 output required)
python src/scripts/10spearman_corr/03b_collect_networks.py \
    --base-dir results/base_correlations \
    --boot-dir results/bootstrap_significant \
    --out-h5 results/differential_network_summary.h5 \
    --out-focus-tsv results/rewiring_hubs.tsv
```

Per-gene metrics written to `rewiring_hubs.tsv`, see column reference in
[GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md).
Optional `--annotate` adds gene symbols/names (requires R).

### Stage 4: Collect Focus Gene Topology

```bash
# Reads per-gene network files from Stage 3 (networks/ directory)
python src/scripts/10spearman_corr/04_collect_focus_gene_topology.py \
    --network-dir results/networks \
    --focus-genes top:50 \
    --n-genes 8763 \
    --n-jobs 4 \
    --out-h5 results/focus_gene_topology.h5
```

Enable via Snakemake with `skip_focus_topology: false` in config (see
[GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md)).

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

### Stage 6: Annotate Rewiring Table
```bash
# Standalone: annotate any TSV containing FBgn gene IDs
Rscript src/scripts/10spearman_corr/06_annotate_rewiring_table.R \
    --input-tsv results/focus_gene_metrics.tsv \
    --output-tsv results/focus_gene_metrics_annotated.tsv

# Or with a custom gene ID column name
Rscript src/scripts/10spearman_corr/06_annotate_rewiring_table.R \
    --input-tsv results/focus_gene_metrics.tsv \
    --gene-id-col gene_id
```

**Adds these columns** (using `org.Dm.eg.db` local database):
- `SYMBOL` - Gene symbol (e.g., Hsp70Aa)
- `GENENAME` - Full gene name/description
- `ENTREZID` - NCBI Entrez Gene ID
- `ENSEMBL` - Ensembl gene ID

**Notes:**
- Uses `src/utils/utils_annotation.R` which provides a reusable `annotate_fbgn()` function
- Can also be triggered automatically from Stage 3 with the `--annotate` flag
- No network access required (local SQLite database)

### Stage 7: R Visualization
```bash
cd results/visualization_data
Rscript visualize_networks.R
```

---

## Permutation Test (Stage 7 — Independent Snakefile)

Justifies observed differential co-expression metrics by comparing them against
a null distribution from random sample assignment. Runs via `Snakefile_permutation`
independently of the main pipeline.

```bash
# Run permutation test for selected focus genes
snakemake -s code/dronetanalysis/src/pipelines/Snakefile_permutation \
    --config \
        expr_h5=results/expression.h5 \
        networks_dir=results/networks \
        indices_h5=results/bootstrap_indices.h5 \
        summary_tsv=results/rewiring_hubs.tsv \
        "perm_genes=[FBgn0001197,FBgn0262739]" \
        n_permutations=100 \
    -j 4
```

**Outputs:**
- `permutation_null/{gi}_{gene_id}.h5` — null distributions per gene
- `permutation_pvals.tsv` — empirical p-values for all metrics
- `rewiring_hubs_with_pvals.tsv` — `rewiring_hubs.tsv` augmented with `pval_*` columns
- `visualization_data/perm_null/*.pdf` — null distribution histograms with observed values

See [GUIDE-08-Permutation-Test.md](GUIDE-08-Permutation-Test.md) for full documentation.

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

## Running a Gene Subset (Testing / Re-running Specific Genes)

To re-run stages 2a and 2b for only a specific set of genes (e.g. to test code changes without processing all ~8,700 genes), use the `gene_subset` config parameter.

### Quick start

```bash
# Dry run — shows what would be executed for the 15 ct_voom focus genes:
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml -n

# Local execution (runs 02a then 02b for each gene in the subset):
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml -j 15

# Run only stage 2a (skip 2b and stage 3):
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --until base_correlations -j 15
```

Snakemake skips output files that already exist, so re-running is safe and only processes missing or outdated files.

### How it works

The `gene_subset` key in the YAML config limits per-gene jobs to named genes:

```yaml
gene_subset:
  - FBgn0001233
  - FBgn0002563
```

- **Stages 0 and 1** always run for the full dataset (bootstrap indices cover all genes).
- **Stages 2a, 2b, 3 (reconstruct_single), 3b (collect_networks), and 4** respect the
  subset — Snakemake only schedules jobs for genes in the subset, via the `_filter_genes()`
  aggregate input function.
- `collect_networks` (Stage 3b) reads only files that exist in `networks/`, so it
  automatically handles partial subsets without any special flags.
- To process all genes, remove or set `gene_subset: []`.

### Config files

| Config | Purpose |
|--------|---------|
| `config/ct_voom_snakemake.yaml` | ct_voom — 15 focus genes (subset testing) |
| `config/hs_voom_snakemake.yaml` | hs_voom — all 8763 genes (full production run) |

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

- [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) - Snakemake workflow: config reference, stage descriptions, troubleshooting
- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) - Understanding network topology metrics
- [GUIDE-08-Permutation-Test.md](GUIDE-08-Permutation-Test.md) - Permutation test: null distribution method, config, and output interpretation
- [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) - Memory optimization strategies (vectorized vs batch mode)
- [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) - Storage optimization (common, minimal, full modes)
- [FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md) - Summary of critical pipeline fixes
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) - Statistical methodology and pipeline design
