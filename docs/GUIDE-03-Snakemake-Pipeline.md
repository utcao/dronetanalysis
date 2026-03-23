# Snakemake Pipeline Guide

## Overview

This guide explains how to run the differential co-expression pipeline using Snakemake
(`src/pipelines/Snakefile_bootstrap`). Snakemake manages job dependencies, parallelism,
and incremental re-runs automatically, making it the preferred method over the manual bash
script for most use cases.

---

## Prerequisites

- Snakemake ≥ 8.0 installed (`conda install -c bioconda snakemake`)
- Python dependencies installed (h5py, numpy, scipy)
- Expression data as TSV or HDF5
- For SGE cluster submission: `pip install snakemake-executor-plugin-cluster-generic`

---

## Configuration

All pipeline parameters are set in a YAML config file. Two example configs are provided:

| Config | Purpose |
|--------|---------|
| `config/ct_voom_snakemake.yaml` | Control (ct_voom), 15-gene subset — Stage 6 enabled |
| `config/hs_voom_snakemake.yaml` | Heat-shock (hs_voom), 15-gene subset — Stage 6 enabled |

### Complete Config Reference

```yaml
# --- Input (provide one of the three) ---
expr_tsv: "dataset/processed/VOOM/voomdataCtrl.txt"  # TSV → auto-converted to HDF5
expr_h5:  "results/expression.h5"                    # existing HDF5, skips conversion
# toy: true                                           # built-in synthetic data

# --- Output ---
out_dir: "results_ct_voom"   # all outputs written here

# --- Bootstrap parameters ---
low_frac:       0.2     # fraction of samples in LOW expression group
high_frac:      0.2     # fraction of samples in HIGH expression group
n_bootstrap:    1000    # number of bootstrap resamples
bootstrap_frac: 0.8     # fraction of group used per resample
seed:           42      # random seed for reproducibility
batch_size:     null    # null = vectorized (fast, ~30 GB RAM);
                        # integer = batch mode (slower, ~12 GB RAM)

# --- Significance ---
fdr_alpha:      0.7               # FDR threshold for edge significance
edge_selection: "sig_differential" # "sig_differential" or "sig_edges"
no_ci_filter:   false             # false = require CI to exclude zero (default)

# --- Storage ---
storage_mode: "common"   # "common" (default), "minimal", or "full"
                         # see OPTIMIZATION-02-Storage.md

# --- Thresholds (Stage 3) ---
min_effect:      0.0    # minimum |Δr| to retain an edge
corr_threshold:  0.0001 # minimum |r| to count as "present"
annotate:        false  # true = annotate output TSV with gene symbols (requires R)

# --- Optional stages ---
skip_visualization:  true  # false = run Stage 5 (visualization export)

# --- Stage 4: Focus gene topology ---
# Set skip_focus_topology: false to enable.
# Requires Stage 3 (reconstruct_single) to complete for all relevant genes first.
skip_focus_topology: true
focus_genes: "top:50"    # "top:N", "0,1,2,...", or "range:start:end"
focus_n_jobs: 4          # parallel workers for Stage 4

# --- Stage 6: Variability analysis (gene_subset only) ---
# Runs per-gene MAD and sample-level ITV plots for every gene in gene_subset.
# Requires expr_tsv to be set (R scripts read the TSV directly, not the HDF5).
skip_variability: false          # false = enable; true = skip (default)
condition_label: "Control"       # label shown in plot titles
variability_save_csv: false      # true = also write per-gene CSV files

# --- Gene subset (optional) ---
# Restrict per-gene jobs (Stages 2a, 2b, 3, 3b, 4) to a named subset of genes.
# Also drives Stage 6 (variability) when skip_variability: false.
# Stages 0 and 1 always process the full expression matrix.
# Remove gene_subset (or set to []) to process all genes.
gene_subset:
  - FBgn0002563
  - FBgn0027844
```

---

## Running the Pipeline

### Dry run (always do this first)

```bash
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    -n
```

> **Note on checkpoint behaviour:** The dry run will show only a small number of jobs
> (`preprocess`, `collect_networks`, `all`) even when dozens of per-gene jobs will run.
> This is expected — Snakemake cannot enumerate per-gene Stage 2a/2b/3 jobs until the
> `preprocess` checkpoint completes and `expression.h5` is created. All stages will
> execute correctly in a real run.

### Full run from scratch (new output folder)

```bash
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --config out_dir=results_ct_voom_v2 \
    -j 15
```

### Re-run only Stage 3b collection (existing intermediate files)

```bash
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --forcerun collect_networks \
    -j 1
```

### Run only through Stage 2b (gene subset testing)

```bash
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --until bootstrap_significant \
    -j 15
```

### Run through Stage 3 including Stage 4 (focus gene topology)

```bash
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --config skip_focus_topology=false \
    -j 15
```

### Run Stage 6 variability analysis (gene_subset configs)

```bash
# Enable variability plots for all genes in gene_subset
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --config skip_variability=false \
    -j 15
```

> **Requirements:** `expr_tsv` must be set in the config (the R scripts read the TSV
> directly). `gene_subset` must be non-empty, otherwise Stage 6 is silently skipped.

### Override config values on the command line

```bash
# Run with a new output folder and a different storage mode
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --config out_dir=results_new storage_mode=minimal annotate=true \
    -j 15
```

### SGE cluster submission (Snakemake 8+/9+)

Snakemake 8+ removed the old `--cluster` flag. Use the `cluster-generic` executor plugin instead.

**Option A — inline (one-off run):**

```bash
mkdir -p logs
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/hs_voom_snakemake.yaml \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -cwd -j y -q long.q -l h_vmem={resources.mem_mb}M -l h_rt={resources.runtime}:0:0 -o logs/ -e logs/ -N {rule}" \
    --default-resources mem_mb=8000 runtime=60 \
    --jobs 200 \
    --latency-wait 60
```

**Option B — profile (recommended, reusable):**

A ready-made profile is provided at `config/sge_profile/config.yaml`. Just run:

```bash
mkdir -p logs
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/hs_voom_snakemake.yaml \
    --profile config/sge_profile
```

**Resource mapping:** Each rule's `resources:` block sets `mem_mb` (MB) and `runtime`
(minutes). These map to SGE's `-l h_vmem` and `-l h_rt` in the submit command.

| Rule | mem_mb | runtime |
|------|--------|---------|
| `base_correlations` | 32000 | 240 min |
| `bootstrap_significant` | 16000 | 120 min |
| `reconstruct_single` | 4000 | 30 min |
| `collect_networks` | 8000 | 120 min |
| `collect_focus_gene_topology` | 8000 | 120 min |
| `gene_mad_variability` | 8000 | 30 min |
| `gene_sample_variability` | 8000 | 30 min |

> **Queue name:** Replace `long.q` in `config/sge_profile/config.yaml` with your cluster's
> actual queue (check with `qstat -g c`).

> **Latency wait:** `--latency-wait 60` (60 seconds) is important on NFS/shared filesystems
> so Snakemake waits for output files to appear after a job finishes.

---

## Pipeline DAG and Stages

```
[expr_tsv / expr_h5 / --toy]
        │
        ▼
   Stage 0: preprocess (CHECKPOINT)
   → expression.h5
        │
        ▼
   Stage 1: bootstrap_indices
   → bootstrap_indices.h5
        │
        ├────────────────────────────────────────────┐
        ▼                                            ▼
   Stage 2a: base_correlations               (per gene, ×N)
   → base_correlations/{gi}_{name}.h5
        │
        ▼
   Stage 2b: bootstrap_significant           (per gene, ×N)
   → bootstrap_significant/{gi}_{name}.h5
        │
        ▼
   Stage 3: reconstruct_single              (per gene, ×N)
   → networks/{gi}_{name}.h5
        │
        ├──────────────────────────────────┐
        ▼                                  ▼
   Stage 3b: collect_networks         Stage 4: collect_focus_gene_topology
   → differential_network_summary.h5  → focus_gene_topology.h5
   → rewiring_hubs.tsv                  (if skip_focus_topology: false)
        │
        └── [Stage 5: visualization, if skip_visualization: false]
            → visualization_data/

[Stage 6: gene_mad_variability + gene_sample_variability, if gene_subset set
 and skip_variability: false — reads expr_tsv directly, one job per subset gene]
→ variability/{gene_id}_mad_variability.pdf
→ variability/{gene_id}_sample_variability_median.pdf
```

### Why Stage 2 shows as a checkpoint

`preprocess` is a Snakemake **checkpoint** because the number and names of per-gene jobs
in Stages 2a, 2b, and 3 depend on the gene names inside `expression.h5`. Snakemake cannot
know how many `{gi}_{gene_id}` targets to create until the file exists. Once `preprocess`
completes, the DAG is re-evaluated and all per-gene jobs are scheduled automatically.

---

## Gene Subset Filtering

Setting `gene_subset` in the config restricts per-gene jobs to only the named genes:

```yaml
gene_subset:
  - FBgn0002563
  - FBgn0027844
```

**Which stages are affected:**

| Stage | Subset filtering |
|-------|-----------------|
| Stage 0 (preprocess) | Always runs on full data |
| Stage 1 (bootstrap_indices) | Always runs on full data |
| Stage 2a (base_correlations) | ✅ Only subset genes scheduled |
| Stage 2b (bootstrap_significant) | ✅ Only subset genes scheduled |
| Stage 3 (reconstruct_single) | ✅ Only subset genes scheduled |
| Stage 3b (collect_networks) | ✅ Reads only existing per-gene files |
| Stage 4 (collect_focus_topology) | ✅ Reads only existing per-gene files |
| Stage 5 (visualization) | Reads aggregated summary, unaffected |
| Stage 6 (variability) | ✅ Only runs when gene_subset is non-empty |

This uses the `_filter_genes()` function and Snakemake's **checkpoint + aggregate input**
pattern: each downstream stage's input is resolved dynamically from the gene names in
`expression.h5`, filtered by `GENE_SUBSET`. Only files whose gene IDs are in the subset
are ever requested, so Snakemake never schedules jobs for excluded genes.

Remove `gene_subset` (or set to `[]`) to process all genes.

---

## Stage Descriptions

### Stage 0 — Preprocess

Converts expression TSV to HDF5, or generates toy data.
- **Script:** `src/scripts/00preprocess/00convert_expr_to_hdf5.py`
- **Output:** `{out_dir}/expression.h5`
- **Triggered by:** `expr_tsv`, `expr_h5`, or `toy: true` in config

### Stage 1 — Bootstrap Indices

Partitions samples into LOW/HIGH groups and generates bootstrap resamples.
- **Script:** `src/scripts/01subset/01get_extreme_pop_bootstrap.py`
- **Output:** `{out_dir}/bootstrap_indices.h5`
- **Key params:** `low_frac`, `high_frac`, `n_bootstrap`, `bootstrap_frac`, `seed`, `batch_size`

### Stage 2a — Base Correlations (per gene, array job)

Computes Spearman correlations for each reference gene's LOW/HIGH subsets with FDR correction.
- **Script:** `src/scripts/10spearman_corr/02a_calc_base_correlations.py`
- **Output:** `{out_dir}/base_correlations/{gi}_{gene_id}.h5` (one per gene)
- **Key params:** `fdr_alpha`, `storage_mode`

### Stage 2b — Bootstrap Significant Edges (per gene, array job)

Bootstraps only the significant edges identified in Stage 2a (~100× speedup).
- **Script:** `src/scripts/10spearman_corr/02b_bootstrap_significant_edges.py`
- **Output:** `{out_dir}/bootstrap_significant/{gi}_{gene_id}.h5` (one per gene)
- **Key params:** `edge_selection`, `storage_mode`

### Stage 3 — Reconstruct Single Network (per gene, array job)

Reconstructs the full differential network for each reference gene: qualitative edge
classification, global topology (LOW/HIGH/DIFF networks with degree distributions), and
focus gene neighborhood analysis (L1/L2 partners with connectivity metrics).
- **Script:** `src/scripts/10spearman_corr/03_reconstruct_diff_network.py`
- **Output:** `{out_dir}/networks/{gi}_{gene_id}.h5` (one per gene, ~4 GB RAM per job)
- **Key params:** `edge_selection`, `min_effect`, `corr_threshold`, `no_ci_filter`
- **Note:** Scheduled as one Snakemake job per gene (like Stages 2a/2b), enabling native
  Snakemake parallelism across a cluster.

### Stage 3b — Collect Networks (single job)

Aggregates per-gene network files from Stage 3 into a summary HDF5 and the
`rewiring_hubs.tsv` ranking table. Reads from `networks/` (Stage 3 output).
- **Script:** `src/scripts/10spearman_corr/03b_collect_networks.py`
- **Output:** `{out_dir}/differential_network_summary.h5`, `{out_dir}/rewiring_hubs.tsv`
- **Key params:** `annotate`
- **Fallback mode:** Can also read raw `--base-dir`/`--boot-dir` directly (standalone use
  without running Stage 3 first).

### Stage 4 — Focus Gene Topology (optional)

Collects per-gene topology matrices (degree_low, degree_high, degree_diff) and L1/L2
metrics across all reference gene contexts for a selected set of focus genes.
- **Script:** `src/scripts/10spearman_corr/04_collect_focus_gene_topology.py`
- **Output:** `{out_dir}/focus_gene_topology.h5`
- **Enable with:** `skip_focus_topology: false` in config
- **Key params:** `focus_genes` ("top:N", "0,1,2,...", or "range:start:end"), `focus_n_jobs`

### Stage 5 — Visualization Export (optional)

Exports the differential network to TSV and GraphML files for R/Cytoscape visualization.
- **Script:** `src/scripts/10spearman_corr/05_prepare_visualization_data.py`
- **Output:** `{out_dir}/visualization_data/`
- **Enable with:** `skip_visualization: false` in config

### Stage 6 — Variability Analysis (optional; gene_subset only)

For each gene in `gene_subset`, produces two variability plots comparing the LOW vs HIGH
sample groups split by that gene's expression rank:

| Rule | Script | Output |
|------|--------|--------|
| `gene_mad_variability` | `plot_gene_mad_variability.R` | `{gene_id}_mad_variability.pdf` — MAD per gene |
| `gene_sample_variability` | `plot_sample_variability.R` | `{gene_id}_sample_variability_median.pdf` — ITV per sample |

- **Scripts:** `src/scripts/15analysis/plot_gene_mad_variability.R`, `src/scripts/15analysis/plot_sample_variability.R`
- **Output:** `{out_dir}/variability/` (one PDF per gene per script)
- **Enable with:** `skip_variability: false` in config (requires `gene_subset` and `expr_tsv`)
- **Key params:** `condition_label`, `variability_save_csv`, `low_frac`, `high_frac`
- **Note:** Scripts read `expr_tsv` directly; no pipeline HDF5 outputs are required.
  Each rule is one Snakemake job per gene, parallelised natively across available cores.

---

## Output Files

After a successful run:

```
{out_dir}/
├── expression.h5                      # Stage 0
├── bootstrap_indices.h5               # Stage 1
├── base_correlations/                 # Stage 2a (per gene) — temp(), auto-deleted after Stage 3
│   ├── 0000_FBgn0000001.h5
│   └── ...
├── bootstrap_significant/             # Stage 2b (per gene) — temp(), auto-deleted after Stage 3
│   ├── 0000_FBgn0000001.h5
│   └── ...
├── networks/                          # Stage 3 (per gene)
│   ├── 0000_FBgn0000001.h5
│   └── ...
├── differential_network_summary.h5    # Stage 3b — aggregated metrics
├── rewiring_hubs.tsv                  # Stage 3b — per-gene TSV (sortable)
├── focus_gene_topology.h5             # Stage 4 (if enabled)
├── visualization_data/                # Stage 5 (if enabled)
│   ├── edges_low.tsv
│   ├── edges_high.tsv
│   ├── edges_diff.tsv
│   └── nodes.tsv
└── variability/                       # Stage 6 (if enabled + gene_subset set)
    ├── FBgn0001197_mad_variability.pdf
    ├── FBgn0001197_sample_variability_median.pdf
    └── ...
```

### rewiring_hubs.tsv columns

Sorted by `L2L1_deg` (descending). Key columns:

| Column | Type | Description |
|--------|------|-------------|
| `gene_idx` | int | Gene index in expression matrix |
| `gene_id` | str | Gene ID (e.g. FBgn identifier) |
| `n_sig_total` | int | Total significant edges |
| `focus_deg_low` | int | Degree of focus gene in LOW network |
| `focus_deg_high` | int | Degree of focus gene in HIGH network |
| `L1_deg_diff` | int | Degree in DIFFERENTIAL network (direct partners) |
| `L1_n_disappear` | int | L1 edges lost (low→high) |
| `L1_n_new` | int | L1 edges gained |
| `L1_n_sign_chg` | int | L1 edges with sign change |
| `L1_n_strengthen` | int | L1 edges strengthened |
| `L1_n_weaken` | int | L1 edges weakened |
| `L1_rewire` | int | L1 rewiring score (disappear+new+sign_chg) |
| `L1_conn_low` | float | Sum \|r\| of L1 edges in LOW condition |
| `L1_conn_high` | float | Sum \|r\| of L1 edges in HIGH condition |
| `L1_conn_diff` | float | Sum \|Δr\| of L1 edges |
| `L1_conn_mean_low` | float | Mean \|r\| per L1 edge in LOW |
| `L1_conn_mean_high` | float | Mean \|r\| per L1 edge in HIGH |
| `L1_mean_abs_dr` | float | Mean \|Δr\| per L1 edge |
| `L2_deg_diff` | int | Indirect partners count (L1∪L2 nodes) |
| `L2_n_disappear` | int | Outer-layer edges lost |
| `L2_n_new` | int | Outer-layer edges gained |
| `L2_n_sign_chg` | int | Outer-layer sign changes |
| `L2_n_strengthen` | int | Outer-layer edges strengthened |
| `L2_n_weaken` | int | Outer-layer edges weakened |
| `L2_rewire` | int | L2 rewiring score |
| `L2_conn_low` | float | Sum \|r\| of outer-layer edges in LOW |
| `L2_conn_high` | float | Sum \|r\| of outer-layer edges in HIGH |
| `L2_conn_diff` | float | Sum \|Δr\| of outer-layer edges |
| `L2_conn_mean_low` | float | Mean \|r\| per outer-layer edge in LOW |
| `L2_conn_mean_high` | float | Mean \|r\| per outer-layer edge in HIGH |
| `L2_mean_abs_dr` | float | Mean \|Δr\| per outer-layer edge |
| `n_direct_edges` | int | Count of focus→L1 edges |
| `n_l1_to_l1_edges` | int | Count of L1↔L1 edges |
| `n_l1_to_l2_edges` | int | Count of L1→L2 edges |
| `L2L1_deg` | float | Ratio: indirect/direct partner count |
| `L2L1_rewire` | float | Ratio: L2/L1 rewiring score |
| `L2L1_conn` | float | Ratio: L2/L1 connectivity change |
| `HL_conn_L1` | float | Ratio: L1 high/low connectivity |
| `HL_conn_L2` | float | Ratio: L2 high/low connectivity |

> **Connectivity (`conn`) vs strength (`str`):** All metrics use `conn` (connectivity),
> defined as sum or mean of absolute correlation values (\|r\|). This replaces the earlier
> `str_` prefix.

---

## Troubleshooting

### Dry run shows only a few jobs

Expected. See [checkpoint behaviour note](#dry-run-always-do-this-first) above.

### "Nothing to be done" on re-run

All output files already exist. To force a specific stage:
```bash
snakemake ... --forcerun collect_networks -j 1
```
To force everything from scratch in a new directory:
```bash
snakemake ... --config out_dir=results_new -j 15
```

### Stage 3b fails with KeyError / missing metrics

Likely a stale `networks/` file from a previous version. Re-run Stage 3:
```bash
snakemake ... --forcerun reconstruct_single -j 15
```

### storage_mode: "common" not recognised

Check that `02a_calc_base_correlations.py` and `02b_bootstrap_significant_edges.py`
support `"common"` as a valid `--storage-mode` value.
See [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) for storage mode details.

---

## Related Reading

- [GUIDE-01-Complete-Workflow.md](GUIDE-01-Complete-Workflow.md) - Full pipeline stage commands and HDF5 schemas
- [GUIDE-02-Network-Metrics.md](GUIDE-02-Network-Metrics.md) - Understanding focus gene connectivity metrics
- [OPTIMIZATION-02-Storage.md](OPTIMIZATION-02-Storage.md) - Storage mode selection (`common`, `minimal`, `full`)
- [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) - Memory mode selection (`batch_size` parameter)
- [REFERENCE-01-Statistical-Methods.md](REFERENCE-01-Statistical-Methods.md) - Statistical methodology

---

**Last Updated:** 2026-03-23
**Status:** ✅ Active — Stage 4 enabled (conditional on `skip_focus_topology: false`); Stage 6 enabled (conditional on `skip_variability: false` + `gene_subset` non-empty)
