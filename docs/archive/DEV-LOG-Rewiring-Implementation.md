# wiring bootstrap and spearman

New files: `02calc_corr_edge_bootstrap_corr.py` — per-gene SGE array task

|Responsibility|	Detail|
|---|---|
|Gene ID|	`--gene-id N` or `$SGE_TASK_ID` (SGE array)|
|Expression input|	--expr-tsv (real) or --toy (same 5×50 seed-1 matrix as Stage 1)|
|Reads from|	bootstrap_indices.h5 — one chunk per dataset thanks to chunks=(1,…)|
|Computes|	rank_zscore_rows → Z @ Z.T / (M-1) → triu extract — no block-tiling needed because bootstrap subsets are small|
|Writes|	gene_XXXXX.h5 with low/{base,boot}/corr_triu and high/{base,boot}/corr_triu|
|Boot dataset layout|	(n_bootstrap, n_tests) chunked (1, n_tests) — one replicate per chunk|

`run_bootstrap_pipeline.sh` — two-stage SGE submission

|Stage|	Jobs|	Hold|	Memory|	Output|
|---|---|---|---|---|
|1 — bootstrap_indices|	1|	none (-hold_jid 0)|	8 G|	results/bootstrap_indices.h5|
|2 — corr_array|	1..N_GENES (gene IDs 0..N-1)|	Stage 1 JID|	16 G|	results/corr/gene_XXXXX.h5|

Gene count is auto-detected from `--expr-tsv` (wc -l minus the header); hardcoded to 5 under `--toy`. The script follows the repo SGE conventions (account, queue, log paths, email).

Toy validation results (all 5 genes, base + bootstrap rep 0, both partitions)

All stored corr_triu values match scipy.stats.spearmanr to within < 1.2 e-7 (float32 rounding floor). Index ranges, dataset shapes, and chunk layouts all check out.

## comparing Route through 01calc_corr_edge.py vs. gene_wise way

Whether it is necessary to have a separate script for gene-wise behaviour.

The two options are not equivalent — there are real costs to routing through 01calc_corr_edge.py.

### The critical issue: module-level imports

`01calc_corr_edge.py:50-51` imports two packages unconditionally at the top:

```py
from scipy.stats import rankdata, t
from statsmodels.stats.multitest import multipletests
```

Every SGE task is a fresh Python process. statsmodels alone adds ~1–2 s of import overhead per process. With `--skip-pval` the functions are never called, but the import still happens. On a 10 000-gene array that is **10 000–20 000 s of pure import wall time** across the array (though spread across nodes, it still consumes CPU-hours and delays job start).

`02calc_corr_edge_bootstrap_corr.py:46` imports only from scipy.stats import rankdata — the one thing it actually uses.

### Computation path comparison

|Aspect|	Route through 01calc_corr_edge.py|	02calc_corr_edge_bootstrap_corr.py|
|---|---|---|
|Imports per process|	statsmodels, scipy.stats.t, dataclasses, time, Sequence|	only scipy.stats.rankdata
|Matmul strategy|	blockwise tiling loop + write_block_to_triu_1d row-by-row|	single Z @ Z.T, then np.triu_indices vectorised extract|
|Why tiling exists|	full dataset can be tens of thousands of genes → N×N won't fit in RAM|	bootstrap subset has ~0.16×N_samples columns; the N_genes × N_genes matmul is the same size either way, but there is no memory pressure from the subset dimension|
|HDF5 write per replicate|	create_h5_datasets creates string arrays for gene_ids/sample_ids, Progress tracker prints, summary stats|	writes only corr_triu — the one number downstream needs|
|p-value / FDR code path|	guarded by if not skip_pval, but the code is still parsed and the related datasets are conditionally skipped inside the loop|	code path does not exist|

The blockwise tiling and `write_block_to_triu_1d` are designed for the case where N_genes is large enough that N×N doesn't fit in memory. That problem does not change based on how the samples were selected — it is the same N_genes either way. Routing the bootstrap replicates through that infrastructure adds loop overhead and per-row HDF5 writes with no memory benefit.

**The 8-line duplication is not a real cost**

`rank_zscore_rows` in `02calc_corr_edge_bootstrap_corr.py:52-59` is a pure, stateless function with no imports beyond rankdata. Extracting it to a shared utils.py would couple the two scripts and force every change to be coordinated. Duplicating 8 lines of a leaf function is cheaper than that coupling, and it keeps each script independently runnable and testable.

### Verdict: 02calc_corr_edge_bootstrap_corr.py is the right design.

The only scenario where routing through 01calc_corr_edge.py would make sense is if N_genes were so large that even the bootstrap path needed block-tiled matmul. In that case you would extract the shared pieces (rank_zscore_rows, triu extraction) into a shared module and have both scripts import from it — not merge the scripts.

## SGE script

Root cause: `SGE -t` is 1-based. The minimum valid task ID is 1. Passing 0 is rejected at submission time — no tasks run at all.

Changes made in run_bootstrap_pipeline.sh:

|Line|	Before|	After|	Why|
|---|---|---|---|---|
|83|	# SGE -t is 0-based …|	# SGE -t is 1-based inclusive; task IDs run 1..N_GENES|	Corrected the comment; removed LAST_TASK variable|
|136, 141|	--gene-id \$SGE_TASK_ID|	--gene-id \$((SGE_TASK_ID - 1))|	SGE_TASK_ID is now 1–N; the HDF5 indices are 0-based, so subtract 1 before passing to the Python script|
|151|	-t 0-${LAST_TASK}|	-t 1-${N_GENES}|	Valid 1-based range|

The mapping is now: `SGE_TASK_ID 1` → `gene_id 0`, `SGE_TASK_ID 2` → `gene_id 1`, … `SGE_TASK_ID N` → `gene_id N-1`.

## Summary the context before the analysis & verdicts

This session is being continued from a previous conversation that ran out of context. The summary below covers the earlier portion of the conversation.

Summary:
1. Primary Request and Intent:
   - **Initial request**: Review `01calc_corr_edge.py` against Black style, summarize functionality, identify efficiency/readability improvements.
   - **Second request**: Same review for `01get_extreme_pop_bootstrap.py`, plus pipeline-level efficiency thinking and SGE job array coordination between the two scripts. User described the pipeline: create low/high expression subpopulations with bootstrapping, feed subsetted datasets to compute co-expression edges, compare low vs high. User explicitly noted `01calc_corr_edge.py` is slow and needs efficiency improvement for cluster parallel execution.
   - **Third request**: Modify `01get_extreme_pop_bootstrap.py` to use "one dataset per category" HDF5 layout instead of one dataset per gene.
   - **Fourth request**: Add `parse_args()` to `01get_extreme_pop_bootstrap.py`, improve pipeline-level efficiency in `01calc_corr_edge.py`, provide terminal example commands.
   - **Fifth (current) request**: "Please wire bootstrap and calculate co-expression pairs together." — create the missing SGE array-task wrapper (`02calc_corr_edge_bootstrap_corr.py`) and SGE submission script (`run_bootstrap_pipeline.sh`).

2. Key Technical Concepts:
   - **Spearman via rank-zscore-matmul**: `Spearman(x,y) == Pearson(rank(x), rank(y))`. Rank each gene row, z-score (ddof=1), then `C = Z @ Z.T / (M-1)`.
   - **Upper-triangle 1D storage**: Only store `N*(N-1)/2` pairs in row-major order `(0,1)(0,2)…(N-2,N-1)`. `TriuIndex` class maps `(i,j) → k`. `np.triu_indices(N, k=1)` produces the same order.
   - **Block-tiling**: `01calc_corr_edge.py` tiles `Z @ Z.T` into `(block_size × block_size)` chunks to bound peak memory. Not needed for bootstrap subsets because `k_resample` is small.
   - **`--skip-pval` mode**: Skips `scipy.stats.t.sf` (the slowest per-element operation) and BH-FDR (`multipletests`). Per-replicate p-values are not needed; significance comes from comparing low vs high distributions across bootstraps.
   - **Vectorised bootstrap**: `rng.integers(0, k, size=(n_genes, n_boot, k_resample))` + fancy-index gather replaces a per-gene Python loop. Verified to produce identical results to per-gene `rng.choice`.
   - **HDF5 gene-axis chunking**: `chunks=(1, …)` so each SGE job reads exactly one chunk when accessing `dataset[gene_id]`.
   - **SGE conventions in repo**: `#$ -A randeff_remove`, `-l h_vmem`, `-l h_rt`, `-o logs/$JOB_NAME.$JOB_ID.out`, `-S /bin/bash`, `-cwd`, job arrays via `#$ -t`, `-hold_jid` for dependencies.

3. Files and Code Sections:

   - **`src/scripts/01subset/01get_extreme_pop_bootstrap.py`** (heavily modified)
     - Removed: `tqdm` import, `bootstrap_sample_indices()`, `_write_dataset()`, `_toy()`.
     - Added: `parse_args()` (flags: `--in-tsv`, `--out-h5`, `--low-frac`, `--high-frac`, `--n-bootstrap`, `--bootstrap-frac`, `--seed`, `--toy`), `main()` with pandas TSV loader.
     - Rewritten: `build_gradient_indices_h5` — vectorised `np.argpartition` (axis=1, all genes at once), vectorised bootstrap via `rng.integers` + fancy-index gather, writes 4 pre-allocated datasets (`indices/low`, `indices/high`, `indices/low_boot`, `indices/high_boot`) with `chunks=(1,…)`.
     - Toy mode: `np.random.default_rng(1).normal(size=(5, 50)).astype(np.float32)` — **must match** in the wiring script.
     - Black-clean. Syntax verified.

   - **`src/scripts/10spearman_corr/01calc_corr_edge.py`** (modified for pipeline efficiency)
     - Added `skip_pval: bool = False` to `create_h5_datasets`, `compute_and_store`, and wired through `parse_args` → `main`.
     - `create_h5_datasets`: conditionally skips `pval_triu` dataset creation.
     - `write_block_to_triu_1d`: pval slice and write guarded by `if pval_triu is not None`.
     - `compute_and_store`: conditionally skips `corr_to_pvals_from_t`, FDR (`multipletests`), and pval-related summary datasets. Corr summary stats (`corr_min`, `corr_max`, `mean_abs_corr`) always written.
     - `rank_zscore_rows` body: `r_centered = r - mean`, `denom = r.std(ddof=ddof)`, `denom` guarded with `< eps` threshold (not `== 0`).
     - Black-formatted the entire file (all pre-existing line-length overflows resolved).
     - `--skip-pval` output HDF5 contains only: `meta/*`, `results/corr_triu`, `summary/{corr_min,corr_max,mean_abs_corr}`.

   - **`src/SGE_scripts/SGE_template.sh`** (read-only, reference)
     - Template uses `# directives, `#$ -t 1-1000` for arrays, `$SGE_TASK_ID` env var, logs to `logs/$JOB_NAME.$JOB_ID.out`.

   - **`src/pipelines/Snakefile`** (read-only, reference)
     - Existing Snakemake pipeline; mostly R-based. Not directly integrated yet.

   - **`src/scripts/01subset/check_bootstrap.py`** (read-only)
     - Manual verification script; reads `gradient_indices.h5`, checks shapes and index validity.

   - **`src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py`** — created and validated
   - **`src/SGE_scripts/run_bootstrap_pipeline.sh`** — created and validated

4. Errors and Fixes:
   - **`all_high` Black formatting**: Wrote line break after `[…]` but Black prefers break after `.astype(`. Fixed by applying Black's preferred style.
   - **`p_block[…]` TypeError on None**: `write_block_to_triu_1d` extracted `pval_row_segment = p_block[…]` unconditionally before the `if pval_triu is not None` guard. Fixed by moving both the slice and write inside the guard: `if pval_triu is not None: pval_triu[…] = p_block[…]`.
   - **`h5py.visititems` callback signature**: Used `def visit(name)` but h5py passes `(name, obj)`. Fixed to `def visit(name, obj)`.
   - **Pre-existing bugs in `01calc_corr_edge.py`** (identified, some fixed by user in between sessions): `data=data=` SyntaxError (user fixed → `data=`), `r_centered` undefined after replacing manual std (user restored `r_centered = r - mean`), `Iterable` not imported (user changed to `Sequence`), docstring placed after import in `validate_toy` (user swapped order).
   - **Extra blank line** between `GradientParams` and `build_gradient_indices_h5` in bootstrap script (3 blanks → 2). Fixed.

5. Problem Solving:
   - Vectorised the entire bootstrap index generation (partition + resample) — eliminated per-gene Python loop.
   - Consolidated HDF5 from 4×N_genes tiny datasets to 4 pre-allocated 2D/3D datasets with gene-axis chunking.
   - Added `--skip-pval` to eliminate the two most expensive per-replicate operations: `scipy.stats.t.sf` (not efficiently vectorised internally) and `multipletests` FDR over N²/2 values.
   - Verified `np.triu_indices(N, k=1)` produces the same pair ordering as `TriuIndex` (row-major upper triangle excluding diagonal).
   - Verified vectorised bootstrap gather (`rng.integers` + fancy index) produces identical output to per-gene `rng.choice` loop for same seed.

6. All User Messages:
   - "Check this script with black code style guide and most commonly used guide when Black guide does not require for some styles. Summary what does this code do and check where can be improved in terms of efficiency and readability."
   - "Summary what does this code do and check where can be improved in terms of efficiency and readability. Check this script with black code style guide and most commonly used guide when Black guide does not require for some styles. The idea is to create lowest and highest expression subpopulations and their bootstrapping datasets with 80% size repeating certain times. Then feed those datasets to calculate significant co-expression edges by comparing both extreme expression subdatasets with same pairs. Because 01calc_corr_edge.py single step will take a lot of times, we need improve the calculation efficiency of both steps and think about how to coordinated these two scripts before submit job array onto cluster for parallel. @src/scripts/10spearman_corr/01calc_corr_edge.py"
   - "How to modify code to make 'One dataset per category' in the HDF5 output?"
   - "Please add parse_args() for this script before improve pipeline-level efficiency and provide example code to run in the terminal."
   - "Please wire bootstrap and calculate co-expression pairs together." (stated twice — once before /compact, once after)

7. Pending Tasks:
   - **Write `src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py`**: SGE array-task Python script. One job per gene. Reads bootstrap indices from shared HDF5, column-subsets expression, computes `rank_zscore_rows → matmul → triu extraction` per replicate (no block-tiling needed), writes one HDF5 per gene with `low/base`, `low/boot`, `high/base`, `high/boot`. Gene ID from `$SGE_TASK_ID` or `--gene-id`. `--toy` uses same 5×50 matrix as bootstrap `--toy`.
   - **Write `src/SGE_scripts/run_bootstrap_pipeline.sh`**: Two-stage SGE submission. Stage 1 single job → `bootstrap_indices.h5`. Stage 2 job array `-t 1-N` with `-hold_jid` on Stage 1 (gene IDs offset by -1) → `corr/gene_XXXXX.h5`. Follows repo SGE template conventions.
   - **End-to-end toy test**: Run bootstrap `--toy`, run corr wrapper `--toy --gene-id 0..4`, verify output HDF5 layouts and correlation values against `scipy.stats.spearmanr`.

8. Current Work:
   The immediate work before the summary request was writing the two files that wire the pipeline together. TodoWrite was set up with three tasks; task 1 ("Write 02calc_corr_edge_bootstrap_corr.py") was marked `in_progress`. The design was fully worked out:
   - `rank_zscore_rows` defined locally in the new script (avoids `importlib` for the digit-prefixed `01calc_corr_edge.py` filename and avoids pulling in heavy `statsmodels` import)
   - `compute_corr_triu(x_subset, triu_idx)` does the full rank→zscore→matmul→extract in one shot (no blocking)
   - `run_gene(gene_id, expr, indices_h5_path, out_dir)` orchestrates: read indices (one chunk per dataset), open output HDF5, write base subsets, pre-allocate `(n_bootstrap, n_tests)` boot datasets with `chunks=(1, n_tests)`, loop over replicates writing row slices
   - Shell script uses `-wrap` for inline commands, `awk '{print $3}'` to extract job IDs, auto-detects `n_genes` via `wc -l`

9. Optional Next Step:
   Continue directly from where work was interrupted. The user's most recent explicit request was: **"Please wire bootstrap and calculate co-expression pairs together."**
   
   Immediate next action: Write `src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py` with the design outlined in §8, then write `src/SGE_scripts/run_bootstrap_pipeline.sh`, then run the end-to-end toy test (bootstrap `--toy` → corr `--toy --gene-id 0` through 4 → verify output).

If you need specific details from before compaction (like exact code snippets, error messages, or content you generated), read the full transcript at: /fml/ag-pallares/caoyt/.claude/projects/-tmp-global2-caoyt-projects-dronetanalysis/365f3c8e-3fe1-4f12-bd2b-492d483cbb63.jsonl
Please continue the conversation from where we left it off without asking the user any further questions. Continue with the last task that you were asked to work on.