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

### ⚡ Optimization Guides

| Document | Description |
|----------|-------------|
| [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Vectorized vs batch processing modes to manage memory (12-30 GB) |
| [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Sparse storage modes reducing disk usage from 7 TB to 0.9 TB (common) or 0.3 TB (minimal) |

### 🔧 Critical Fixes

| Document | Description |
|----------|-------------|
| [FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md) | Executive summary of all three critical pipeline fixes |
| [FIX-02-HDF5-Attributes.md](FIX-02-HDF5-Attributes.md) | Fix for HDF5 64 KB attribute limit by storing gene names as datasets |

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
2. **[REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md)** - Deep dive into statistical methods and pipeline design

---

## Troubleshooting Index

### Common Errors and Solutions

| Error | Document | Section |
|-------|----------|---------|
| `OSError: Unable to synchronously create attribute (object header message is too large)` | [FIX-02-HDF5-Attributes.md](FIX-02-HDF5-Attributes.md) | Root Cause |
| `numpy._core._exceptions._ArrayMemoryError: Unable to allocate X GiB` | [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Solution |
| `ValueError: All chunk dimensions must be positive` | [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Edge Cases |
| Out of disk space (files ~1.3 GB each) | [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Problem Solved |
| Stage 1 memory allocation insufficient | [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Memory Requirements by Stage |

### Performance Issues

| Issue | Document | Recommendation |
|-------|----------|----------------|
| Slow Stage 1 execution | [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Use vectorized mode (default) |
| High memory usage (>30 GB) | [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) | Use batch mode with --batch-size 100 |
| Excessive disk usage (>1 TB) | [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Use --storage-mode common (default) |
| Need maximum space savings | [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) | Use --storage-mode minimal |

---

## Pipeline Stages Overview

| Stage | Script | Input | Output | Key Parameters |
|-------|--------|-------|--------|----------------|
| **0** | 00convert_expr_to_hdf5.py | expression.tsv | expression.h5 | - |
| **1** | 01get_extreme_pop_bootstrap.py | expression.h5 | bootstrap_indices.h5 | --batch-size (optional) |
| **2a** | 02a_calc_base_correlations.py | expression.h5, bootstrap_indices.h5 | base_correlations/ | --storage-mode |
| **2b** | 02b_bootstrap_significant_edges.py | expression.h5, base_correlations/ | bootstrap_significant/ | --storage-mode |
| **3** | 03_reconstruct_diff_network.py | base_correlations/, bootstrap_significant/ | differential_network_summary.h5 | - |
| **4** | 04_collect_focus_gene_topology.py | differential_network_summary.h5 | focus_gene_topology/ | --focus-genes |
| **5** | 05_prepare_visualization_data.py | differential_network_summary.h5 | visualization_data/ | - |

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

**Last Updated:** 2026-02-13
**Pipeline Version:** 2.0 (with optimizations)
**Status:** ✅ All documentation current
