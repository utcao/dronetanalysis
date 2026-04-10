# Documentation Index

Welcome to the dronetanalysis documentation! This index provides quick navigation to all available guides and references.

---

## Quick Start

**New to the pipeline?** Start here:
1. [Complete Workflow Guide](GUIDE-01-Complete-Workflow.md) - Step-by-step pipeline execution guide
2. [Memory Optimization](OPTIMIZATION-01-Memory.md) - Choose vectorized or batch mode for your dataset
3. [Storage Optimization](OPTIMIZATION-02-Storage.md) - Reduce disk usage by 87-97%

---

## Documentation Categories

### 📘 User Guides

| Document | Description |
|----------|-------------|
| [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) | Complete pipeline workflow from expression data to differential network analysis |
| [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) | Understanding network topology metrics (degree, betweenness, clustering, etc.) |
| [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) | Snakemake workflow: config reference, all parameters, gene subset filtering, Stage 4 |
| [GUIDE-04-Qualitative-Change-Metrics.md](GUIDE-04-Qualitative-Change-Metrics.md) | Qualitative change categories (DISAPPEAR/NEW/SIGN_CHANGE/STRENGTHEN/WEAKEN), partition identities, L1/L2 focus-gene metrics, sanity checks |
| [GUIDE-05-Expression-Variability-Analysis.md](GUIDE-05-Expression-Variability-Analysis.md) | Gene-level MAD variability and sample-level ITV across LOW/HIGH expression groups; single-gene violin plot + genome-wide batch xlsx summary with BH-FDR |
| [GUIDE-06-Expression-Violin-Plot.md](GUIDE-06-Expression-Violin-Plot.md) | Violin/boxplot with individual selection (top/bottom X%), dot highlighting, and Wilcoxon rank-sum + BH-FDR annotation; exports PDF/SVG/PNG |
| [GUIDE-07-Pathway-Enrichment.md](GUIDE-07-Pathway-Enrichment.md) | GO/KEGG/GSEA enrichment on top or bottom N rewiring hub genes ranked by L2L1 ratio metrics; Excel + PDF output |
| [GUIDE-08-Permutation-Test.md](GUIDE-08-Permutation-Test.md) | Permutation test validating that observed differential co-expression metrics are driven by the expression gradient; null distributions, empirical p-values, histogram plots |

### ⚡ Optimization Guides

| Document | Description |
|----------|-------------|
| [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Vectorized vs batch processing modes to manage memory (12-30 GB) |
| [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Sparse storage modes reducing disk usage from 7 TB to 0.9 TB (common) or 0.3 TB (minimal) |
| [OPTIMIZATION-03-Network-Reconstruction.md](OPTIMIZATION-03-Network-Reconstruction.md) | Stage 3 performance optimization providing 200-800x speedup for network reconstruction |

### 🔧 Critical Fixes

| Document | Description |
|----------|-------------|
| [FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md) | Executive summary of all three critical pipeline fixes |
| [FIX-02-HDF5-Attributes.md](FIX-02-HDF5-Attributes.md) | Fix for HDF5 64 KB attribute limit by storing gene names as datasets |
| [FIX-03-HPC-SGE-Pipeline.md](FIX-03-HPC-SGE-Pipeline.md) | Fix for HDF5 concurrent read failure on NFS and inflated mem_mb values on SGE HPC |

### 📊 Technical References

| Document | Description |
|----------|-------------|
| [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) | Statistical testing methods (bootstrap, Fisher's Z) and efficient pipeline design |

### 📋 Documentation Standards

| Document | Description |
|----------|-------------|
| [00-RULES-Documentation-Standards.md](00-RULES-Documentation-Standards.md) | Standards for creating and maintaining documentation (naming, templates, cross-references) |

### 🗂️ Archived Documents

Historical development logs and planning documents:

| Document | Description |
|----------|-------------|
| [archive/DEV-LOG-Bootstrap-Implementation.md](archive/DEV-LOG-Bootstrap-Implementation.md) | Development log for bootstrap correlation reconstruction |
| [archive/DEV-LOG-Rewiring-Implementation.md](archive/DEV-LOG-Rewiring-Implementation.md) | Development log for network rewiring implementation |
| [archive/DEV-PLAN-Focus-Gene-Collection.md](archive/DEV-PLAN-Focus-Gene-Collection.md) | Planning document for focus gene TSV collection feature |

---

## Recommended Reading Order

### For New Users

1. **[GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md)** - Understand the full pipeline
2. **[OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md)** - Choose your processing mode
3. **[OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md)** - Configure storage settings

### For Troubleshooting

1. **[FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md)** - Check if your issue is a known fix
2. **Specific fix guides** - FIX-02-HDF5-Attributes.md, OPTIMIZATION-01-Memory.md, OPTIMIZATION-02-Storage.md
3. **[REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md)** - Understand statistical methodology

### For Advanced Users

1. **[GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md)** - Interpret network topology results
2. **[GUIDE-04-Qualitative-Change-Metrics.md](GUIDE-04-Qualitative-Change-Metrics.md)** - Understand qualitative change categories and sanity-check identities
3. **[REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md)** - Deep dive into statistical methods and pipeline design

---

## Troubleshooting Index

### Common Errors and Solutions

| Error | Document | Section |
|-------|----------|---------|
| `OSError: Unable to synchronously create attribute (object header message is too large)` | [FIX-02-HDF5-Attributes.md](FIX-02-HDF5-Attributes.md) | Root Cause |
| `OSError: Unable to synchronously open file (bad object header version number)` on HPC | [FIX-03-HPC-SGE-Pipeline.md](FIX-03-HPC-SGE-Pipeline.md) | Fix 1 |
| `numpy._core._exceptions._ArrayMemoryError: Unable to allocate X GiB` | [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Solution |
| `ValueError: All chunk dimensions must be positive` | [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Edge Cases |
| Out of disk space (files ~1.3 GB each) | [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Problem Solved |
| Stage 1 memory allocation insufficient | [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Memory Requirements by Stage |

### Performance Issues

| Issue | Document | Recommendation |
|-------|----------|----------------|
| Slow Stage 1 execution | [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Use vectorized mode (default) |
| Slow Stage 3 execution (hub genes) | [OPTIMIZATION-03-Network-Reconstruction.md](OPTIMIZATION-03-Network-Reconstruction.md) | Optimizations auto-enabled (200-800x faster) |
| High memory usage (>30 GB) | [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Use batch mode with --batch-size 100 |
| Excessive disk usage (>1 TB) | [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Use --storage-mode common (default) |
| Need maximum space savings | [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Use --storage-mode minimal |

---

## Pipeline Stages Overview

| Stage | Script | Input | Output | Key Parameters |
|-------|--------|-------|--------|----------------|
| **0** | 00convert_expr_to_hdf5.py | expression.tsv | expression.h5 | - |
| **1** | 01get_extreme_pop_bootstrap.py | expression.h5 | bootstrap_indices.h5 | --batch-size (optional) |
| **2a** | 02a_calc_base_correlations.py | expression.h5, bootstrap_indices.h5 | base_correlations/{gi}_{gene}.h5 | --storage-mode |
| **2b** | 02b_bootstrap_significant_edges.py | expression.h5, base_correlations/ | bootstrap_significant/{gi}_{gene}.h5 | --storage-mode |
| **3** | 03_reconstruct_diff_network.py | base_correlations/{gi}.h5, bootstrap_significant/{gi}.h5 | networks/{gi}_{gene}.h5 | --edge-selection |
| **3b** | 03b_collect_networks.py | networks/ | differential_network_summary.h5, rewiring_hubs.tsv | --annotate |
| **4** | 04_collect_focus_gene_topology.py | networks/ | focus_gene_topology.h5 | --focus-genes, --n-jobs |
| **5** | 05_prepare_visualization_data.py | differential_network_summary.h5 | visualization_data/ | - |
| **6** | 06_annotate_rewiring_table.R | rewiring_hubs.tsv | rewiring_hubs_annotated.tsv | - |
| **7** | 07_permutation_test.py | expression.h5, networks/ | permutation_null/{gi}_{gene}.h5 | --n-permutations, --perm-genes |
| **7b** | 08_collect_permutation_pvals.py | permutation_null/ | permutation_pvals.tsv | --alpha |

Stages 2a, 2b, 3, 3b, and 4 respect `gene_subset` config — only named genes are processed.
Stage 4 is optional; enable with `skip_focus_topology: false` in config.
Stages 7 and 7b run via `Snakefile_permutation` (independent of `Snakefile_bootstrap`).

For detailed stage information, see [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md).

---

## Key Features Implemented

### ✅ Memory Optimization (2026-02-13)
- **Dual-mode processing**: Vectorized (fast, 30 GB) or batch (memory-efficient, 12 GB)
- **Automatic memory allocation**: Pipeline adjusts SGE job requirements based on mode
- **50% memory reduction**: Changed int64 to int32 for bootstrap indices
- **Status**: Implemented and tested

### ✅ Storage Optimization (2026-02-13)
- **Three storage modes**: common (default, 87% reduction), minimal (97% reduction), full (legacy)
- **Sparse storage**: Only significant edges stored, not all 38M edges
- **Removed unused arrays**: Eliminated p-values, q-values, and boot/delta arrays never read downstream
- **Backwards compatible**: Stage 3 auto-detects and reads both old and new formats
- **Status**: Implemented and tested

### ✅ HDF5 Attribute Fix (2026-02-13)
- **Gene names as datasets**: Avoids 64 KB attribute size limit
- **No downstream changes required**: All scripts already read from datasets
- **Supports unlimited genes**: Can handle 20,000+ genes without errors
- **Status**: Implemented and tested

---

## Getting Help

### Pipeline Execution

For complete command-line examples and SGE cluster usage:
```bash
# See the complete workflow guide
less docs/GUIDE-01-Complete-Workflow.md

# Or run the pipeline script with --help
bash src/SGE_scripts/run_bootstrap_pipeline.sh --help
```

### Test Dataset

Run the pipeline on toy data to verify installation:
```bash
bash src/SGE_scripts/run_bootstrap_pipeline.sh \
    --toy --local \
    --out-dir test_output \
    --n-bootstrap 10
```

### Parameter Selection

- **Memory mode**: See [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) → "Recommended Batch Sizes"
- **Storage mode**: See [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) → "Storage Mode Details"
- **Network metrics**: See [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) → "Metric Definitions"

---

## Recent Updates

### 2026-04-08
- ✅ **New Stage 7**: Added permutation test pipeline (`Snakefile_permutation`) to validate that observed differential co-expression metrics are driven by the expression gradient, not random chance
  - `src/scripts/10spearman_corr/07_permutation_test.py` — per-gene null distribution via focus-row-only computation (O(N×k) not O(N²×k))
  - `src/scripts/10spearman_corr/08_collect_permutation_pvals.py` — aggregate to `permutation_pvals.tsv` and `rewiring_hubs_with_pvals.tsv`
  - `src/scripts/15analysis/09_plot_permutation_null.R` — null distribution histograms with observed values
  - `src/pipelines/Snakefile_permutation` — independent Snakefile (no changes to existing pipeline)
- ✅ **Documentation**: Added [GUIDE-08-Permutation-Test.md](GUIDE-08-Permutation-Test.md)
- ✅ **New Script**: Added `src/scripts/15analysis/pathway_enrichment_hubs.R` — GO/KEGG/GSEA enrichment on top/bottom N rewiring hub genes; ranked by any numeric column (default `L2L1_rewire`); cascading p-value cutoffs; optional GO simplification and GSEA; Excel + PDF output
- ✅ **Documentation**: Added [GUIDE-07-Pathway-Enrichment.md](GUIDE-07-Pathway-Enrichment.md)

### 2026-03-22
- ✅ **New Script**: Added `src/scripts/plot_expr_violin.py` — violin/boxplot for expression matrices with individual selection (top/bottom X%), dot highlighting, Wilcoxon rank-sum + BH-FDR, PDF/SVG/PNG export
- ✅ **`--gene_name` parameter**: human-readable gene symbol auto-appended to output file names (`{prefix}_{gene_id}_{gene_name}.pdf`)
- ✅ **Output folder** renamed from `boxplot_violinplot/` to `results/expr_violin/`
- ✅ **Documentation**: Added [GUIDE-06-Expression-Violin-Plot.md](GUIDE-06-Expression-Violin-Plot.md)

### 2026-04-10
- ✅ **New Script**: Added `src/scripts/15analysis/summarize_all_genes_mad_variability.R` — genome-wide batch MAD variability analysis; loads expression matrix once, iterates over all genes, applies BH-FDR correction, joins gene symbols from a mapping file, writes ranked xlsx table via `openxlsx` (bold headers, auto column widths) with console summary (% significant, % increased/decreased variability)
- ✅ **Documentation**: Updated [GUIDE-05-Expression-Variability-Analysis.md](GUIDE-05-Expression-Variability-Analysis.md) — new Script 2 section with argument table, output column reference, spot-check instructions, and `writexl` prerequisite note

### 2026-03-21 (variability analysis)
- ✅ **New Analysis Scripts**: Added two standalone R scripts for expression variability analysis:
  - `plot_gene_mad_variability.R` — per-gene MAD across LOW/HIGH sample groups (one dot per gene); Wilcoxon test; `--partner-type` stub for future direct/indirect partner restriction
  - `plot_sample_variability.R` — per-sample ITV (median/mean/sum of |log-expr − population mean|) across LOW/HIGH groups (one dot per sample); `--summary-metric` embedded in output filename
- ✅ **Documentation**: Added [GUIDE-05-Expression-Variability-Analysis.md](GUIDE-05-Expression-Variability-Analysis.md)

### 2026-03-21
- ✅ **Qualitative Metric Fixes** (Stage 3 + 3b):
  - Fixed `focus_deg_low` / `focus_deg_high` always being 0 in summary output (were not saved to `focus_gene/metrics/` HDF5 group due to missing prefix in `save_results`)
  - Fixed inconsistent STRENGTHEN/WEAKEN definition in `compute_edge_stats`: now uses pure `qual_score` codes (consistent with global `qual_summary` and mode-B counts); previously double-counted SIGN_CHANGE edges
  - Added `L1_n_edges_low`, `L1_n_edges_high`, `L2_n_edges_low`, `L2_n_edges_high` metrics throughout (Stage 3 HDF5, Stage 3b summary HDF5, TSV output)
- ✅ **Documentation**: Added [GUIDE-04-Qualitative-Change-Metrics.md](GUIDE-04-Qualitative-Change-Metrics.md) with full documentation of the 6 qualitative categories, partition identities, and L1/L2 sanity checks; updated [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) with differential network section and cross-reference

### 2026-03-20
- ✅ **HPC SGE Fixes**: Fixed `OSError: bad object header version number` on NFS by adding
  `HDF5_USE_FILE_LOCKING=FALSE` to SGE submit command; corrected 10× inflated `mem_mb`
  values in Snakefile; increased `latency-wait` to 120s
- ✅ **Documentation**: Added [FIX-03-HPC-SGE-Pipeline.md](FIX-03-HPC-SGE-Pipeline.md)

### 2026-03-11
- ✅ **Pipeline Refactoring**: Stage 3 split into single-gene `reconstruct_single` (per-gene
  array job) + `03b_collect_networks.py` (aggregation). Stage 4 now fully wired in Snakemake.
- ✅ **Gene Subset Propagation**: `gene_subset` config now propagates through Stages 2a, 2b,
  3, 3b, and 4 via `_filter_genes()` aggregate inputs. Stages 0 and 1 unaffected.
- ✅ **Stage 4 Enabled**: `collect_focus_gene_topology` rule active; enable with
  `skip_focus_topology: false`. Fixed `mean_delta` attribute bug in Stage 4 reader.
- ✅ **Snakemake Pipeline**: Added `GUIDE-03-Snakemake-Pipeline.md` with full config
  reference, stage descriptions, checkpoint behaviour explanation, and Gene Subset Filtering section
- ✅ **Network Metrics Renamed**: `str_*` → `conn_*` throughout Stage 3 and TSV output
  (connectivity = sum/mean \|r\|); fixed stale key inconsistency bug in `compute_edge_stats`
- ✅ **New Focus Gene Metrics**: `focus_deg_low/high`, `L1/L2_conn_mean_low/high`,
  `L2_n_disappear/new/sign_chg/strengthen/weaken`, `L2_mean_abs_dr`,
  `n_direct_edges`, `n_l1_to_l1_edges`, `n_l1_to_l2_edges`
- ✅ **Two-layer Edge Definition**: `two_layer_edges` now = L1↔L1 + L1→L2 (outer layers
  only); `full_two_layer_edges` = complete 2-hop ego network including focus→L1

### 2026-02-14
- ✅ **Network Reconstruction Optimization**: Implemented 200-800x speedup for Stage 3
  - Edge filtering replaces O(L1²) nested loops (10-100x faster for hub genes)
  - Multiprocessing parallelization (16x speedup on 16-core machines)
  - Vectorized operations and early filtering (20-30% additional gains)
- ✅ **Documentation**: Added OPTIMIZATION-03-Network-Reconstruction.md

### 2026-02-13
- ✅ **Storage Optimization**: Implemented sparse storage with 87-97% disk savings
- ✅ **Memory Optimization**: Added batch processing mode for large datasets
- ✅ **HDF5 Fix**: Fixed attribute size limit for datasets with 20,000+ genes
- ✅ **Documentation**: Created comprehensive fix summaries and optimization guides

---

## Contributing

When adding new documentation:
1. Follow the naming convention: `[CATEGORY]-[NUMBER]-[Name].md` (see [00-RULES-Documentation-Standards.md](00-RULES-Documentation-Standards.md))
2. Add entry to this README.md index
3. Include "Related Reading" section with cross-references
4. Update FIX-01-Critical-Issues-Summary.md if documenting a fix

---

**Last Updated:** 2026-04-10
**Pipeline Version:** 2.3 (permutation test: Stage 7 + Snakefile_permutation)
**Status:** ✅ All documentation current
