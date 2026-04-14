# Multi-Input Workflow Guide

## Overview

The pipeline (`Snakefile_bootstrap`) processes **one expression input at a time**. Each run
is tied to a single configfile, a single `out_dir`, and a single Snakemake state directory
(`.snakemake/`). When two or more expression matrices must be processed independently (e.g.,
cell-type voom-normalised data **and** host-species voom-normalised data), a separate run
directory is created per input instead of running everything from the project root.

This guide documents the convention used in this project and explains how to replicate it for
new inputs.

---

## Limitation

Running multiple inputs from the same working directory causes conflicts:

| Resource | Conflict |
|----------|----------|
| `.snakemake/` state directory | Snakemake tracks job status per working directory; two concurrent runs overwrite each other's state |
| `out_dir` in config | Both runs would write to the same output path unless manually changed every time |
| Log files | SGE and Snakemake log files collide |

There is currently **no native multi-input support** in `Snakefile_bootstrap`. Each input
requires its own invocation.

---

## Solution: Per-Input Run Directories with Symlinks

Create a dedicated subdirectory for each input. Inside each directory, symlink `code/` and
`data/` back to the shared project-level copies to avoid duplication.

### Directory Layout

```
modifiers_dronet/          ← project root
├── code/                  ← shared pipeline code (one copy)
├── data/                  ← shared reference/input data (one copy)
├── run_voomct/            ← run directory for cell-type voom input
│   ├── code -> ../code    (symlink)
│   ├── data -> ../data    (symlink)
│   ├── results_ct_voom/   ← pipeline output for this input
│   └── tmp/
└── run_voomhs/            ← run directory for host-species voom input
    ├── code -> ../code    (symlink)
    ├── data -> ../data    (symlink)
    ├── results_hs_voom/   ← pipeline output for this input
    └── tmp/
```

Each run directory has its own independent:
- Snakemake state (`.snakemake/`)
- Output directory (`results_*/`)
- Temporary files (`tmp/`)
- Log files (stored at the project root as `log_dronet_<tag>_*.txt`)

The symlinked `code/` and `data/` are read-only during normal pipeline execution, so sharing
them is safe even when both runs execute concurrently.

---

## Step-by-Step: Adding a New Input

### Step 1 — Create the run directory

```bash
cd /fml/ag-pallares_projects/modifiers_dronet
mkdir run_<tag>          # e.g. run_voomct, run_voomhs, run_newexpr
```

### Step 2 — Create symlinks

```bash
cd run_<tag>
ln -s ../code code
ln -s ../data data
```

### Step 3 — Create a config file

Copy an existing config and adjust the input path and output directory:

```bash
cp code/dronetanalysis/config/ct_voom_snakemake.yaml \
   code/dronetanalysis/config/<tag>_snakemake.yaml
```

Key fields to update in the new config:

```yaml
expr_h5: data/<new_input>.h5          # or expr_tsv: data/<new_input>.tsv
out_dir: results_<tag>/               # unique output dir within this run directory
gene_subset: [...]                    # adjust as needed
```

### Step 4 — Run the pipeline from inside the run directory

```bash
cd run_<tag>
snakemake -s code/dronetanalysis/src/pipelines/Snakefile_bootstrap \
    --configfile code/dronetanalysis/config/<tag>_snakemake.yaml \
    --profile code/dronetanalysis/config/sge_profile \
    -j 100
```

Snakemake writes its state to `run_<tag>/.snakemake/`, fully isolated from other runs.

---

## Running Multiple Inputs Concurrently

Both run directories can be active at the same time because:

- Snakemake state is isolated per directory.
- Output paths are separate (`results_ct_voom/` vs `results_hs_voom/`).
- `code/` and `data/` symlinks are read-only during execution.

Open two terminals (or submit both via SGE) and launch each from its own run directory.

---

## Existing Run Directories in This Project

| Directory | Input | Config | Output |
|-----------|-------|--------|--------|
| `run_voomct/` | Cell-type voom-normalised expression | `ct_voom_snakemake.yaml` | `results_ct_voom/` |
| `run_voomhs/` | Host-species voom-normalised expression | `hs_voom_snakemake.yaml` | `results_hs_voom/` |

## Available Data Sources

In addition to the expression matrices, the project contains a SNP genotype file that can be integrated with MAD/ITV results or used as a standalone input for genotype-based analyses:

| File | Location | Samples | SNPs | Notes |
|------|----------|---------|------|-------|
| `Dmel_head_hs_ct_Miss80_MAF5_LD8_HWE_1975ind.vcf` | `data/snp/` | 1,975 | 413,348 | Same `{lineID}_{wellID}` sample names as expression matrices |

The SNP file uses the same `{lineID}_{wellID}` naming convention (e.g. `106_A10`) as both `voomdataCtrl.txt` and `voomdataHS.txt`. This means the sample axis can be joined directly between expression MAD results and SNP-derived metrics without any name remapping.

> **Important for multi-input analyses:** 31 inbred lines contribute multiple replicate samples each. When combining expression and SNP data, always group by `lineID` (the numeric prefix) — individual replicates within the same line are near-genetic clones and must not be treated as independent observations. See [`docs/dataset_snp_structure.md`](../../../../docs/dataset_snp_structure.md) for the full sample breakdown per line.

For SNP-specific operations (subsetting, numeric conversion, sequence reconstruction), see [`docs/guide_snp_operations.md`](../../../../docs/guide_snp_operations.md).

---

## Troubleshooting

### Snakemake says jobs are already running / locked

Each run directory maintains its own lock. If you see a lock error, ensure you are running
from inside the correct `run_<tag>/` directory, not from the project root.

```bash
# Unlock if the previous run was interrupted:
cd run_<tag>
snakemake -s code/dronetanalysis/src/pipelines/Snakefile_bootstrap \
    --configfile code/dronetanalysis/config/<tag>_snakemake.yaml \
    --unlock
```

### Output files appear in the wrong directory

Verify `out_dir` in the configfile points to a path **relative to the run directory**, not
an absolute path. For example `out_dir: results_ct_voom/` resolves to
`run_voomct/results_ct_voom/` when Snakemake is invoked from `run_voomct/`.

### Symlink points to wrong location

```bash
ls -la run_<tag>/code run_<tag>/data
# Should show: code -> ../code   data -> ../data
```

If a symlink is broken, recreate it:

```bash
cd run_<tag>
rm code && ln -s ../code code
rm data && ln -s ../data data
```

---

## Related Reading

- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) - Full pipeline workflow from expression data to differential network analysis
- [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) - Config reference, all parameters, SGE profile setup
- [FIX-03-HPC-SGE-Pipeline.md](FIX-03-HPC-SGE-Pipeline.md) - HDF5 locking fixes relevant when running on NFS with concurrent jobs
- [docs/dataset_snp_structure.md](../../../../docs/dataset_snp_structure.md) — SNP VCF structure, sample naming, and joining with expression data
- [docs/guide_snp_operations.md](../../../../docs/guide_snp_operations.md) — Practical SNP file operations for subsetting and conversion

---

**Last Updated:** 2026-04-14
**Status:** Active workaround (no native multi-input support in Snakefile_bootstrap)
