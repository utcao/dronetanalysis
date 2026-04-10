# GUIDE-NM-02: Multi-Dataset Support (Future To-Do)

> **Status:** Not yet implemented. This document describes the planned refactor to allow
> the pipeline to process multiple expression datasets in a single Snakemake invocation.

---

## Overview

The current pipeline processes **one expression dataset at a time**. All output paths are
module-level string constants derived from a single `out_dir` config key, making it
impossible to fan out over multiple datasets in one run.

The planned refactor introduces a `{dataset}` wildcard throughout all rules, a `datasets:`
dict in the config, and a per-dataset config helper — a standard Snakemake multi-sample
pattern. Shared parameters (thresholds, compression) remain at the top level and can be
overridden per-dataset.

**Is this reasonable for a pipeline?** Yes — this is the idiomatic Snakemake pattern for
running the same workflow over multiple samples. The main complexity comes from the
checkpoint + dynamic per-gene parallelization, which requires careful wildcard handling but
is well-supported by Snakemake.

---

## Future To-Do

### FTD-NM-01: Introduce `{dataset}` wildcard throughout all rules

#### 1. Config — new structure (`config/example.yaml`)

```yaml
# Base output directory — each dataset lands at {out_dir}/{dataset}/
out_dir: "results/network_metrics"

# Dataset definitions (dict keys become {dataset} wildcard values)
datasets:
  ctrl:
    expr_tsv: "data/VOOM/ctrl.txt"
  treated:
    expr_tsv: "data/VOOM/treated.txt"
    skip_global_metrics: false        # per-dataset override of global default
    compute_betweenness: false

# Shared parameters (apply to all datasets unless overridden in their block)
fdr_alpha:             0.05
corr_threshold:        0.0001
compression_level:     6
# chunk_size: null
skip_global_metrics:   true
compute_betweenness:   false
betweenness_n_samples: 500
```

Settings that are **per-dataset** (must live inside the dataset block):
`expr_tsv`, `expr_h5`, `toy`, `gene_subset`

Settings that are **global or per-dataset** (top-level default, per-dataset override):
`fdr_alpha`, `corr_threshold`, `chunk_size`, `skip_global_metrics`,
`compute_betweenness`, `betweenness_n_samples`

Settings that are **global only** (no per-dataset override makes sense):
`out_dir`, `compression_level`

#### 2. Snakefile — header changes

Replace all module-level path constants with `{dataset}`-parameterized templates.
The critical syntax: double-braces inside an f-string produce a literal `{...}` wildcard:

```python
DATASETS = config.get("datasets", {})
OUT      = config.get("out_dir", "results/network_metrics")

# Backwards-compat shim: wrap legacy flat config as single "default" dataset
if not config.get("datasets"):
    ds_keys = {"expr_tsv", "expr_h5", "toy", "gene_subset",
               "skip_global_metrics", "compute_betweenness", "betweenness_n_samples"}
    config["datasets"] = {"default": {k: config[k] for k in ds_keys if k in config}}
DATASETS = config["datasets"]

def _ds_cfg(ds, key, default=None):
    """Resolve a setting: per-dataset override → global config → hardcoded default."""
    return DATASETS[ds].get(key, config.get(key, default))

# Path templates — {{dataset}} → literal "{dataset}" wildcard after f-string eval
EXPR_H5_PAT     = f"{OUT}/{{dataset}}/expression.h5"
CORR_H5_PAT     = f"{OUT}/{{dataset}}/corr_significant.h5"
NETWORK_H5_PAT  = f"{OUT}/{{dataset}}/network_edges.h5"
PER_GENE_PAT    = f"{OUT}/{{dataset}}/per_gene/{{gi}}_{{gene_id}}.h5"
SUMMARY_H5_PAT  = f"{OUT}/{{dataset}}/network_metrics_summary.h5"
SUMMARY_TSV_PAT = f"{OUT}/{{dataset}}/network_metrics.tsv"
GLOBAL_H5_PAT   = f"{OUT}/{{dataset}}/global_topology.h5"
```

Add `dataset` to `wildcard_constraints` to prevent path-separator ambiguity:

```python
wildcard_constraints:
    gi      = r"\d{4}",
    dataset = r"[^/]+",
```

#### 3. Gene discovery functions

```python
def _filter_genes(names, gene_subset):
    pairs = list(enumerate(names))
    if gene_subset:
        pairs = [(i, n) for i, n in pairs if n in gene_subset]
    return pairs

def all_per_gene_files(wildcards):
    # Pass dataset= explicitly to locate the right checkpoint output
    h5_path = checkpoints.preprocess.get(dataset=wildcards.dataset).output[0]
    names = _read_gene_names(h5_path)
    gene_subset = set(_ds_cfg(wildcards.dataset, "gene_subset", []))
    per_gene_dir = f"{OUT}/{wildcards.dataset}/per_gene"
    return [f"{per_gene_dir}/{i:04d}_{name}.h5"
            for i, name in _filter_genes(names, gene_subset)]
```

#### 4. `rule all`

```python
def _all_targets(wildcards):
    targets = (expand(SUMMARY_H5_PAT,  dataset=list(DATASETS)) +
               expand(SUMMARY_TSV_PAT, dataset=list(DATASETS)))
    for ds in DATASETS:
        if not _ds_cfg(ds, "skip_global_metrics", True):
            targets.append(GLOBAL_H5_PAT.format(dataset=ds))
    return targets

rule all:
    input: _all_targets
```

#### 5. `checkpoint preprocess`

```python
checkpoint preprocess:
    output: EXPR_H5_PAT
    params:
        toy      = lambda wc: _ds_cfg(wc.dataset, "toy", False),
        expr_tsv = lambda wc: _ds_cfg(wc.dataset, "expr_tsv", ""),
        expr_h5  = lambda wc: _ds_cfg(wc.dataset, "expr_h5", ""),
    resources:
        mem_mb=16000, runtime=30, threads=4,
    run:
        if params.toy:
            shell(f"python {PREPROCESS_SCR} --toy --out-h5 {{output}}")
        elif params.expr_h5:
            # Symlink external file so output matches EXPR_H5_PAT
            shell(f"ln -sfn {params.expr_h5} {{output}}")
        elif params.expr_tsv:
            shell(f"python {PREPROCESS_SCR} --expr-tsv {params.expr_tsv} "
                  f"--out-h5 {{output}} --compression gzip "
                  f"--compression-level {COMPRESSION}")
        else:
            raise ValueError(f"Dataset '{wildcards.dataset}': provide expr_tsv, expr_h5, or toy=true.")
```

> **Note on `expr_h5`:** The original pipeline allowed the checkpoint output to *be* the
> external file (by setting `EXPR_H5_FINAL = EXPR_H5`). With the `{dataset}` wildcard, the
> output path must match `EXPR_H5_PAT`. A symlink avoids data duplication. If the external
> file is later moved or deleted, the symlink will break.

#### 6. `compute_correlations`

Per-dataset params via lambdas; `COMPRESSION` remains global:

```python
rule compute_correlations:
    input:  expr = EXPR_H5_PAT
    output: CORR_H5_PAT
    params:
        fdr_alpha      = lambda wc: _ds_cfg(wc.dataset, "fdr_alpha", FDR_ALPHA),
        corr_threshold = lambda wc: _ds_cfg(wc.dataset, "corr_threshold", CORR_THRESHOLD),
        chunk_flag     = lambda wc: (
            f"--chunk-size {_ds_cfg(wc.dataset, 'chunk_size', CHUNK_SIZE)}"
            if _ds_cfg(wc.dataset, "chunk_size", CHUNK_SIZE) else ""
        ),
    resources:
        mem_mb=8000, runtime=60, threads=8,
    shell: ...
```

#### 7. `build_network` and `per_gene_metrics`

Update `input:` / `output:` to `*_PAT`. No logic changes needed.
`per_gene_metrics` gains `{dataset}` automatically via `NETWORK_H5_PAT` and `PER_GENE_PAT`.

#### 8. `collect_metrics`

```python
rule collect_metrics:
    input:
        per_gene_files = all_per_gene_files,
        network        = NETWORK_H5_PAT,
    output:
        h5  = SUMMARY_H5_PAT,
        tsv = SUMMARY_TSV_PAT,
    params:
        per_gene_dir = lambda wc: f"{OUT}/{wc.dataset}/per_gene",
    ...
```

#### 9. `global_metrics`

**Remove** the `if not SKIP_GLOBAL_METRICS:` guard. Define unconditionally — whether it
runs is controlled entirely by `rule all` requesting its output per-dataset:

```python
rule global_metrics:
    input:  network = NETWORK_H5_PAT
    output: GLOBAL_H5_PAT
    params:
        betw_flag     = lambda wc: "--compute-betweenness" if _ds_cfg(wc.dataset, "compute_betweenness", False) else "",
        betw_nsamples = lambda wc: _ds_cfg(wc.dataset, "betweenness_n_samples", BETWEENNESS_NSAMPLES),
    ...
```

---

## Files to Modify

| File | Change |
|------|--------|
| `Snakefile` | Full header rewrite + all rules updated with `{dataset}` wildcard |
| `config/example.yaml` | Restructure to `datasets:` dict |

Scripts (`scripts/01_*.py` … `04_*.py`) do **not** need changes — they operate on individual
HDF5 files and are unaware of the dataset concept.

---

## Verification

```bash
# Dry-run — should show separate branches per dataset in the DAG
snakemake -s network_metrics/Snakefile \
    --configfile network_metrics/config/example.yaml -n

# Visualize DAG
snakemake -s network_metrics/Snakefile \
    --configfile network_metrics/config/example.yaml \
    --dag | dot -Tpng > dag_multi.png

# Toy smoke test with two datasets
snakemake -s network_metrics/Snakefile \
    --config 'datasets={"ctrl": {"toy": true}, "treated": {"toy": true}}' \
    out_dir=tmp/test_multi -j 8
# Expected: tmp/test_multi/ctrl/network_metrics.tsv
#           tmp/test_multi/treated/network_metrics.tsv
```

---

## Pitfalls

| # | Issue | Mitigation |
|---|-------|------------|
| 1 | `{{dataset}}` vs `{dataset}` in f-strings | Double braces escape to literal `{dataset}`. Single braces would try to evaluate a Python variable and raise `NameError`. |
| 2 | `checkpoints.preprocess.get(**wildcards)` → `.get(dataset=wildcards.dataset)` | Explicit keyword prevents subtle bugs if additional wildcards appear on the call site. |
| 3 | `COMPRESSION` is parse-time global in shell strings | Fine as long as compression is not per-dataset. If per-dataset compression is ever needed, move to a `params:` lambda. |
| 4 | `expr_h5` symlink breaks if source moves | Document in config comments. |
| 5 | Backwards compat: existing flat configs produce output at `{out_dir}/default/...` | Acceptable change. The shim avoids a hard failure; update configs when convenient. |

---

*Document version: 2026-04-10*
*Pipeline: `network_metrics/`*
*Status: Future To-Do — not yet implemented*
