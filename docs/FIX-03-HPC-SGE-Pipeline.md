# HPC / SGE Pipeline Fixes

## Overview

Three issues observed when running `Snakefile_bootstrap` on an SGE HPC cluster:

1. **Fatal**: `OSError: Unable to synchronously open file (bad object header version number)` — HDF5 concurrent read failure on NFS
2. **Resource**: `mem_mb` values in the Snakefile were 10× larger than required
3. **Resource**: `latency-wait` was insufficient for NFS filesystem propagation

---

## Fix 1 — HDF5 Concurrent Read Failure on NFS [CRITICAL]

### Problem

Stage 2a (`base_correlations`) fails with:

```
File "/tmp/global2/caoyt/miniforge3/lib/python3.12/site-packages/h5py/_hl/files.py"
  ...
OSError: Unable to synchronously open file (bad object header version number)
```

The job reads `bootstrap_indices.h5` (written once by Stage 1) and `expression.h5`
(written once by Stage 0). Both fail when many parallel jobs access them simultaneously.

### Root Cause

**This is NOT a file-incompleteness problem.** Snakemake's DAG guarantees Stage 1 fully
exits (exit code 0, file closed and flushed) before any Stage 2a job is submitted.
The `latency-wait` setting then adds an NFS propagation buffer on top.

The error arises from **NFS + HDF5 file locking interaction**:

1. Stage 2a submits up to ~8,763 jobs to SGE nearly simultaneously
2. Many jobs start at similar times and all call `h5py.File(bootstrap_indices_h5, "r")`
3. HDF5 internally calls `fcntl(F_SETLK)` to acquire a POSIX read lock
4. On NFS, POSIX locks are routed through the **Network Lock Manager (NLM)** daemon
5. Hundreds of simultaneous `fcntl()` calls to the same NFS file overwhelm NLM, which
   returns errors, timeouts, or stale lock state to some clients
6. A client that received a bad NLM response reads a stale NFS page for the HDF5
   object header region → the version byte is 0x00 or garbage instead of 0x01/0x02
7. HDF5 raises `"bad object header version number"`

The file itself is not corrupt — only the client's NFS page cache view is inconsistent
due to the failed lock operation.

### Solution

Set `HDF5_USE_FILE_LOCKING=FALSE` in every SGE job environment. This bypasses the
broken `fcntl()` / NLM path entirely. Each job opens the file directly via normal NFS
I/O, which is consistent because no writes are occurring concurrently.

**Why this is safe here:**

| File | Writer status when Stage 2a runs | Stage 2a access |
|------|----------------------------------|-----------------|
| `bootstrap_indices.h5` | Stage 1 fully done (Snakemake DAG) | Read-only, 8763 concurrent readers |
| `expression.h5` | Stage 0 fully done (checkpoint) | Read-only, 8763 concurrent readers |
| `base_correlations/{gi}.h5` | Written by this job only | Single writer, no concurrent readers |

Since there is no concurrent writer for the shared files, locking provides no safety
benefit — it only introduces the NLM contention that causes the error.

**When `HDF5_USE_FILE_LOCKING=FALSE` would NOT be safe:**
- A job writes to a file while another job reads the same file simultaneously
- SWMR (Single Writer Multiple Reader) mode is in use (not applicable here)

### Change Applied

**File:** `config/sge_profile/config.yaml`

```yaml
# Before
cluster-generic-submit-cmd: "qsub ... -N {rule}"

# After
cluster-generic-submit-cmd: "qsub ... -v HDF5_USE_FILE_LOCKING=FALSE -N {rule}"
```

The `-v` flag passes the environment variable to every submitted SGE job.

---

## Fix 2 — Inflated `mem_mb` Values in Snakefile [RESOURCE]

### Problem

All `mem_mb` resource values in `src/pipelines/Snakefile_bootstrap` were 10× larger
than documented requirements. This caused either:
- Jobs queued indefinitely waiting for high-memory nodes
- Unnecessary slot consumption on shared clusters

### Root Cause

Extra zero appended to each value (e.g., `32000` → `320000`).

### Solution

**File:** `src/pipelines/Snakefile_bootstrap`

| Rule | Before (MB) | After (MB) | Reference |
|------|-------------|------------|-----------|
| `preprocess` | 160,000 | 16,000 | Stage 0 peak ~0.5 GB |
| `bootstrap_indices` (batch mode) | 160,000 | 16,000 | Batch peak ~3–15 GB |
| `bootstrap_indices` (vectorized) | 30,000 | 30,000 | ✅ already correct |
| `base_correlations` | 320,000 | 32,000 | Stage 2a peak ~1.5 GB |
| `bootstrap_significant` | 160,000 | 16,000 | Stage 2b peak ~14 GB |
| `reconstruct_single` | 80,000 | 8,000 | Stage 3 peak ~4 GB |
| `collect_networks` | 80,000 | 8,000 | Stage 3b peak ~2 GB |
| `collect_focus_gene_topology` | 80,000 | 8,000 | Stage 4 peak ~2 GB |

> For detailed memory requirements by stage, see
> [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md).

---

## Fix 3 — Increase `latency-wait` [RESOURCE]

### Problem

`latency-wait: 60` (60 seconds) may not be sufficient for NFS filesystems under heavy
load when ~8,763 jobs are submitted in a short window. If Snakemake checks for output
files before NFS propagates them from the writer node, it prematurely declares a job
failed.

### Solution

**File:** `config/sge_profile/config.yaml`

```yaml
# Before
latency-wait: 60

# After
latency-wait: 120
```

This doubles the NFS propagation window. Combined with Fix 1 (no lock contention),
Stage 2a jobs will reliably find and read `bootstrap_indices.h5` after Stage 1 completes.

---

## Verification

After applying fixes, test Stage 1 → Stage 2a transition on a small gene subset:

```bash
# 1. Run through Stage 1 only
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --until bootstrap_indices \
    --profile config/sge_profile

# 2. Confirm the HDF5 file is readable
python -c "
import h5py
with h5py.File('results_ct_voom/bootstrap_indices.h5', 'r') as f:
    print('Keys:', list(f.keys()))
    print('n_bootstrap:', f['meta'].attrs['n_bootstrap'])
print('OK')
"

# 3. Run Stage 2a on the test subset (15 genes in ct_voom config)
snakemake -s src/pipelines/Snakefile_bootstrap \
    --configfile config/ct_voom_snakemake.yaml \
    --until base_correlations \
    --profile config/sge_profile

# 4. Confirm no "bad object header" errors in job logs
grep -r "bad object header" tmp/logs/
```

---

## Summary of Changes

| File | Change |
|------|--------|
| `config/sge_profile/config.yaml` | Added `-v HDF5_USE_FILE_LOCKING=FALSE` to `qsub` command |
| `config/sge_profile/config.yaml` | `latency-wait`: 60 → 120 |
| `src/pipelines/Snakefile_bootstrap` | `mem_mb` corrected for 7 rules (removed extra zero) |

---

## Related Reading

- [GUIDE-03-Snakemake-Pipeline.md](GUIDE-03-Snakemake-Pipeline.md) - SGE cluster submission setup, resource mapping table
- [OPTIMIZATION-01-Memory.md](OPTIMIZATION-01-Memory.md) - Detailed memory requirements per stage
- [FIX-01-Critical-Issues-Summary.md](FIX-01-Critical-Issues-Summary.md) - Other critical pipeline fixes

---

**Last Updated:** 2026-03-20
**Status:** ✅ Applied — fixes committed to `config/sge_profile/config.yaml` and `src/pipelines/Snakefile_bootstrap`
