# Snakemake Workflow Efficiency and Modular Organization Guide

## Problem

Naïve Snakemake workflows grow brittle as analyses scale:

- Paths and tool locations hardcoded inside every rule — one rename touches every file
- One monolithic Snakefile becomes hundreds of lines and hard to navigate
- Cluster resources duplicated or inconsistent across rules
- The `--cluster` flag was removed in Snakemake 9, silently breaking SGE profiles
- Reserved resource names (e.g. `runtime`) reject valid `HH:MM:SS` strings
- Long analyses submitted interactively die when the SSH connection drops
- Unclear how to add focal genes mid-project without accidentally re-running existing results

This guide covers the patterns that address each of these issues, using the Q1/Q5 DEG/DVG
workflows (`src/q1q5_degdvg/`) as worked examples throughout.

---

## Solution

### Pattern 1: Centralized path variables

**Naïve** — hardcode paths inside every rule:

```python
rule fit_gamlss:
    output: "../../results/gamlss_0608_2026/{gene}_{cond}_Q1vsQ5_gamlss_dvgdeg.csv"
    log:    "../../logs/gamlss_0608_2026/{gene}_{cond}.log"
    shell:
        "/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript ..."
```

**Efficient** — declare all paths as top-level variables; reference them in rules:

```python
RSCRIPT    = "/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript"
GAMLSS_DIR = "../../results/gamlss_0608_2026"
LOG_DIR    = "../../logs/gamlss_0608_2026"

rule fit_gamlss:
    output: GAMLSS_DIR + "/{gene}_{cond}_Q1vsQ5_gamlss_dvgdeg.csv"
    log:    LOG_DIR    + "/{gene}_{cond}.log"
    shell:  "{RSCRIPT} src/gamlss_q1q5/run_gamlss_quintile_analysis.R ..."
```

Changing an output directory now requires editing one line, not every rule.

---

### Pattern 2: Read gene lists from external files

**Naïve** — hardcode the input list in the workflow:

```python
FOCAL_GENES = ["FBgn0027560", "FBgn0039752"]
```

**Efficient** — read from an external file that can be edited without touching the
Snakefile. Using a CSV allows you to carry annotation columns (symbol, description, etc.)
alongside the IDs the workflow actually uses:

```python
import csv
with open("src/gamlss_q1q5/focal_genes.csv") as fh:
    FOCAL_GENES = [row["gene_id"] for row in csv.DictReader(fh)]
```

`focal_genes.csv` (`src/gamlss_q1q5/focal_genes.csv`):
```csv
gene_id,symbol,description
FBgn0027560,ken,ken and barbie
FBgn0039752,CG12345,CG12345
FBgn0001234,sna,snail
```

The workflow reads only `gene_id`; annotation columns are ignored by the pipeline but
available for manual reference or downstream scripts. Adding or removing a gene is a
one-row CSV edit — no Snakefile change, no dry-run surprises.

---

### Pattern 3: Use `executor: cluster-generic` instead of `--cluster`

Snakemake 9 **removed** the `--cluster` command-line flag entirely. The replacement is
the `cluster-generic` executor plugin:

```bash
# install once, in the env used to run snakemake
pip install snakemake-executor-plugin-cluster-generic
```

**Old profile `config.yaml` (broken in Snakemake 9)**:

```yaml
cluster: "qsub -pe parallel {threads} -l h_vmem={resources.mem_mb}M -l h_rt=48:0:0"
jobs: 50
```

Snakemake 9 ignores the `cluster:` key and treats the qsub string as a file target,
producing `MissingRuleException: No rule to produce qsub -A ...`.

**New profile `config.yaml`**:

```yaml
executor: cluster-generic

# >- joins lines with spaces and strips the trailing newline → one-line qsub command
cluster-generic-submit-cmd: >-
  qsub
  -A child_gamlss_q1q5
  -pe parallel {threads}
  -l h_vmem={resources.mem_gb}G
  -l h_rt={resources.walltime}
  -N {rule}_{wildcards.gene}_{wildcards.cond}
  -j y
  -cwd
  -o /dev/null

jobs: 50
latency-wait: 120    # seconds to wait for NFS output after job completes
keep-going: true
rerun-incomplete: true
printshellcmds: true
```

Snakemake expands `{threads}`, `{resources.*}`, `{rule}`, `{wildcards.*}`, and `{jobid}`
in the submit command. Wildcards let you embed meaningful job names in the SGE queue
(e.g. `gamlss_quintile_FBgn0027560_CT`).

---

### Pattern 4: Avoid the reserved `runtime` resource

Snakemake ≥7 reserves `resources.runtime` as **integer minutes**. Passing a `HH:MM:SS`
string causes:

```
WorkflowError: Resource 'runtime' with value '48:0:0' could not be parsed as minutes
```

**Fix**: use any other name. The convention here is `walltime`:

```python
rule fit_gamlss:
    resources:
        mem_gb   = 25,
        walltime = "48:0:0"   # passed as {resources.walltime} in the cluster submit cmd
```

```yaml
cluster-generic-submit-cmd: >-
  qsub
  -l h_rt={resources.walltime}
  ...
```

`walltime` is an arbitrary user resource — Snakemake passes it through unvalidated.

---

### Pattern 5: Submit Snakemake as a master qsub job

Running Snakemake interactively means the workflow dies if your connection drops.
A lightweight master job keeps Snakemake alive for the full analysis duration:

```bash
#!/bin/sh
#$ -A gamlss_q1q5
#$ -pe parallel 1
#$ -l h_vmem=4G
#$ -l h_rt=560:0:0     # must exceed longest child job (48 h) + queue wait time
#$ -N gamlss_master
#$ -j y
#$ -cwd

set -e
set -x

SNAKEMAKE=/tmp/global2/caoyt/miniforge3/bin/snakemake

mkdir -p ../../results/gamlss_0608_2026
mkdir -p ../../logs/gamlss_0608_2026

$SNAKEMAKE \
  -s src/gamlss_q1q5/compute_degdvg_gamlss_q1q5.snakemake \
  --profile src/gamlss_q1q5/profiles/sge
```

The master job runs only Snakemake (4 GB, 1 slot). All heavy work is in child jobs.
Email notifications can be added with `-m beas -M user@example.com`.

---

### Pattern 6: Two-environment split via full interpreter paths

A common pattern: Snakemake needs one conda environment (with the executor plugin),
but the R analysis scripts need another (`ganlss` env with GAMLSS packages).

**Solution**: never `conda activate` inside child jobs. Call the interpreter by its
full path from the target environment:

```python
# In the Snakefile — RSCRIPT points to the ganlss env's Rscript
RSCRIPT = "/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript"
```

```bash
# master job: uses base env (snakemake + executor plugin)
# child jobs: use {RSCRIPT} directly — no conda activate needed
```

This avoids conda activation inside qsub child jobs, which can fail or inherit the
wrong environment depending on SGE's shell initialization order.

---

## Multi-File Rule Organization

For complex pipelines, a single Snakefile becomes unwieldy. Snakemake's `include:`
directive merges separate `.smk` rule files at parse time — identical to writing
everything in one file, but organized by concern.

### When to split

| Split into separate files | Keep in one file |
|---|---|
| Rules with substantially different concerns (preprocess / fit / aggregate) | Pipeline is under ~100 lines |
| Rules conditionally included (`if config.get("run_permutation"): include: ...`) | Rules share many interdependent `expand()` calls |
| Rules shared across multiple pipelines | You want the whole workflow readable in one pass |

### Directory layout

```
src/
  pipelines/
    Snakefile_main             ← master Snakefile: config, rule all, include: calls
    rules/
      preprocess.smk           ← condition filtering and quintile assignment
      fit_models.smk           ← GAMLSS model fitting rules
      aggregate.smk            ← result collection and summary
  gamlss_q1q5/
    profiles/
      sge/config.yaml          ← cluster executor profile
    focal_genes.txt            ← one gene ID per line
```

### Master Snakefile (`Snakefile_main`)

Declares all globals first, then includes rule files:

```python
# ==== globals ====
RSCRIPT    = "/tmp/global2/caoyt/miniforge3/envs/ganlss/bin/Rscript"
RESULTS    = "../../results/gamlss"
LOG_DIR    = "../../logs/gamlss"
CONDITIONS = {"CT": 1, "HS": 2}

with open("src/gamlss_q1q5/focal_genes.txt") as fh:
    FOCAL_GENES = [l.strip() for l in fh if l.strip() and not l.startswith("#")]

# ==== include rule files (order matters if rules reference each other) ====
include: "rules/preprocess.smk"
include: "rules/fit_models.smk"
include: "rules/aggregate.smk"

# ==== target rule — always defined in the master ====
rule all:
    input:
        expand(RESULTS + "/{gene}_{cond}_Q1vsQ5_gamlss_dvgdeg.csv",
               gene=FOCAL_GENES, cond=CONDITIONS.keys())
```

### Rule file (`rules/fit_models.smk`)

Rule files inherit all globals declared before the `include:` statement — no imports
or re-declarations needed:

```python
# RSCRIPT, RESULTS, LOG_DIR, CONDITIONS, FOCAL_GENES are all available here.

rule fit_gamlss:
    output:
        main    = RESULTS + "/{gene}_{cond}_Q1vsQ5_gamlss_dvgdeg.csv",
        summary = RESULTS + "/{gene}_{cond}_Q1vsQ5_gamlss_dvgdeg_sum.csv"
    params:
        condition_val = lambda wildcards: CONDITIONS[wildcards.cond]
    threads: 4
    resources:
        mem_gb   = 25,
        walltime = "48:0:0"
    log:
        LOG_DIR + "/{gene}_{cond}.log"
    shell:
        """
        mkdir -p {LOG_DIR}
        {RSCRIPT} src/gamlss_q1q5/run_gamlss_quintile_analysis.R \
            --focal-gene {wildcards.gene} \
            --full-run \
            --condition-val {params.condition_val} \
            > {log} 2>&1
        """
```

### How `include:` works

`include:` is **textual inclusion at parse time**, not a Python import:

- ✅ All globals declared before the `include:` line are available in the included file
- ✅ Rules in included files are visible to the master's `rule all` input block
- ✅ Included files can themselves include other files
- ❌ Included files cannot be run standalone (`snakemake -s rules/fit_models.smk`
  would fail — the globals are missing)
- ❌ Circular includes are not supported

### Conditional includes

Use Python `if` to optionally enable a rule set:

```python
if config.get("run_permutation", False):
    include: "rules/permutation.smk"
```

Controlled from the command line with `--config run_permutation=true` or from a
`config.yaml` file passed via `--configfile`.

### Sharing rules across pipelines

When two pipelines (e.g. `Snakefile_CT` and `Snakefile_HS`) need the same fitting
rules, place shared logic in a single `.smk` file and include it from both:

```python
# Snakefile_CT
CONDITION = "CT"
COND_VAL  = 1
include: "rules/fit_models.smk"   # uses CONDITION and COND_VAL

# Snakefile_HS
CONDITION = "HS"
COND_VAL  = 2
include: "rules/fit_models.smk"   # same file, different globals
```

The included file uses whatever `CONDITION` was set by the master before the
`include:` call.

---

### Pattern 7: Incremental runs — add genes without re-running existing results

Snakemake decides whether to run a job by checking whether its **output files already
exist on disk**. It does not re-run a job simply because the input gene list changed.
This means adding genes to `focal_genes.csv` and re-submitting will only run the new
genes — existing result CSVs are untouched.

**How it works:**

The `with open(...)` block that reads `focal_genes.csv` runs as plain Python at
**parse time**, outside any rule. Snakemake does not track it as a rule dependency.
Editing the CSV expands `FOCAL_GENES`, which adds new expected output paths. For each
path Snakemake checks: does the file exist? If yes → skip. If no → run.

```python
# This block is NOT a tracked Snakemake input — editing the CSV never
# invalidates existing outputs.
with open("src/q1q5_degdvg/shared/focal_genes.csv") as fh:
    FOCAL_GENES = [row["gene_id"] for row in csv.DictReader(fh)]
```

**Rules for safe incremental use:**

1. **Keep the output directory name unchanged.**
   Output paths are built from `GAMLSS_DIR` / `DESEQ2_DIR` / `LIMMA_DIR`. If you
   rename those variables (e.g. `deseq2_0609` → `deseq2_0610`), Snakemake won't
   find the old files and will try to re-run every gene. Change the directory name
   only when you intentionally want a fresh run with new parameters.

2. **Always dry-run first after editing the CSV.**
   ```bash
   /tmp/global2/caoyt/miniforge3/bin/snakemake \
     -s src/q1q5_degdvg/q1q5_deseq2/compute_deg_deseq2_q1q5.snakemake -n -p
   ```
   The dry-run lists every job that *would* be submitted. Verify it shows only the
   new genes × conditions before submitting.

3. **Never use `-F` / `--forceall` when protecting existing results.**
   `--forceall` bypasses all existence checks and re-runs every job unconditionally.
   It is the one flag that can overwrite existing CSVs.

4. **`rerun-incomplete: true` is safe but worth understanding.**
   Set in all three SGE profiles, this flag tells Snakemake to re-run any job whose
   output file was only partially written (e.g. the cluster job was killed mid-run).
   It applies to new genes as well — if a new gene's job is killed, the incomplete
   output is detected and re-run automatically on the next submission.

5. **Removing a gene from the CSV does not delete its results.**
   Snakemake never removes output files. The old CSV can be trimmed freely; the
   corresponding result files remain on disk until you delete them manually.

**Worked example — adding two genes to an existing run:**

```bash
# 1. Edit the CSV (append two rows)
echo "FBgn0040765,InR,Insulin receptor"   >> src/q1q5_degdvg/shared/focal_genes.csv
echo "FBgn0000137,aos,argos"              >> src/q1q5_degdvg/shared/focal_genes.csv

# 2. Dry-run — should show exactly 4 new jobs (2 genes × 2 conditions)
/tmp/global2/caoyt/miniforge3/bin/snakemake \
  -s src/q1q5_degdvg/q1q5_deseq2/compute_deg_deseq2_q1q5.snakemake -n -p

# 3. Submit (master qsub; only the 4 new jobs are dispatched)
qsub src/q1q5_degdvg/q1q5_deseq2/submit_deseq2_q1q5.sh
```

---

## Quick Reference: Naïve vs Efficient

| Issue | Naïve | Efficient |
|---|---|---|
| Change output directory | Edit every rule | Change one variable at the top |
| Add focal gene | Edit Snakefile | Add one row to `focal_genes.csv` |
| Re-run only new genes | Delete old results or `--forceall` | Just submit — Snakemake skips existing outputs |
| Accidentally re-run everything | Use `--forceall` or rename output dir | Never rename dir mid-analysis; never use `-F` |
| Snakemake 9 + SGE | `--cluster` → `MissingRuleException` | `executor: cluster-generic` + plugin |
| `HH:MM:SS` walltime | `WorkflowError: runtime could not be parsed` | Use `walltime` resource |
| Connection drop mid-run | Snakemake process killed, jobs orphaned | Master qsub job continues unattended |
| Workflow > 200 lines | Hard to navigate | `include:` splits into focused `.smk` files |
| Two conda envs | `conda activate` in child job (fragile) | Full interpreter path in `RSCRIPT` variable |
| Partial output from crashed job | Stale file blocks re-run | `rerun-incomplete: true` in profile detects and re-runs |

---

## Related Reading

- [GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md](GUIDE-19-GAMLSS-DEG-DVG-Q1Q5.md) — Q1/Q5 analysis workflow (worked example of all patterns above)
- [GUIDE-20-DESeq2-limma-voom-DEG-Q1Q5.md](GUIDE-20-DESeq2-limma-voom-DEG-Q1Q5.md) — DESeq2 and limma-voom Q1/Q5 workflow
- [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) — Snakemake config reference and SGE submission for the bootstrap pipeline
- [FIX-03-HPC-SGE-Pipeline.md](FIX-03-HPC-SGE-Pipeline.md) — SGE-specific fixes: NFS locking, `latency-wait`, memory allocation

---

**Last Updated:** 2026-06-09
**Status:** ✅ Active — patterns validated with Snakemake 9.13.2 + snakemake-executor-plugin-cluster-generic 0.3.x
