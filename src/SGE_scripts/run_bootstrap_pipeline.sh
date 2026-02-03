#!/bin/bash
# =============================================================================
# run_bootstrap_pipeline.sh
# Two-stage SGE submission for the bootstrap co-expression pipeline.
#
# Stage 1 (single job)   – generate shared bootstrap indices  → bootstrap_indices.h5
# Stage 2 (job array)    – per-gene Spearman correlations     → corr/gene_XXXXX.h5
#
# Usage
# -----
#   # Real data
#   bash src/SGE_scripts/run_bootstrap_pipeline.sh \
#       --expr-tsv   dataset/voomdataCtrl.txt \
#       --out-dir    results \
#       --n-bootstrap 50 --seed 42
#
#   # Toy smoke test (no input file needed; 5 genes, 50 samples)
#   bash src/SGE_scripts/run_bootstrap_pipeline.sh --toy --out-dir results
#
# Notes
# -----
# * Requires qsub (SGE / UGE).
# * For --toy, gene count is hard-coded to 5.  For real data, n_genes is
#   auto-detected from the header row of --expr-tsv (number of columns minus 1
#   would give samples; row count minus 1 gives genes).
# * Stage 2 holds on the Stage 1 job ID so it cannot start until indices exist.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# defaults
# ---------------------------------------------------------------------------
EXPR_TSV=""
OUT_DIR="results"
LOW_FRAC=0.2
HIGH_FRAC=0.2
N_BOOTSTRAP=50
BOOTSTRAP_FRAC=0.8
SEED=42
TOY=false

# ---------------------------------------------------------------------------
# parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --expr-tsv)       EXPR_TSV="$2";        shift 2 ;;
        --out-dir)        OUT_DIR="$2";          shift 2 ;;
        --low-frac)       LOW_FRAC="$2";         shift 2 ;;
        --high-frac)      HIGH_FRAC="$2";        shift 2 ;;
        --n-bootstrap)    N_BOOTSTRAP="$2";      shift 2 ;;
        --bootstrap-frac) BOOTSTRAP_FRAC="$2";   shift 2 ;;
        --seed)           SEED="$2";             shift 2 ;;
        --toy)            TOY=true;              shift   ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------
if [[ "$TOY" == false && -z "$EXPR_TSV" ]]; then
    echo "ERROR: provide --expr-tsv or use --toy." >&2
    exit 1
fi

INDICES_H5="${OUT_DIR}/bootstrap_indices.h5"
CORR_DIR="${OUT_DIR}/corr"

# ---------------------------------------------------------------------------
# detect gene count
# ---------------------------------------------------------------------------
if [[ "$TOY" == true ]]; then
    N_GENES=5
else
    # expression TSV: first column is gene ID; every subsequent line is one gene
    # subtract 1 for the header row
    N_GENES=$(( $(wc -l < "$EXPR_TSV") - 1 ))
fi
# SGE -t is 1-based inclusive; task IDs run 1..N_GENES

echo "=============================================="
echo "  Bootstrap co-expression pipeline"
echo "  toy=$TOY  n_genes=$N_GENES  n_bootstrap=$N_BOOTSTRAP  seed=$SEED"
echo "  indices  → $INDICES_H5"
echo "  corr     → $CORR_DIR/"
echo "=============================================="

# ---------------------------------------------------------------------------
# Stage 1 – bootstrap indices (single job, no dependency)
# ---------------------------------------------------------------------------
if [[ "$TOY" == true ]]; then
    STAGE1_CMD="python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
        --toy \
        --out-h5 $INDICES_H5 \
        --low-frac $LOW_FRAC --high-frac $HIGH_FRAC \
        --n-bootstrap $N_BOOTSTRAP --bootstrap-frac $BOOTSTRAP_FRAC \
        --seed $SEED"
else
    STAGE1_CMD="python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
        --in-tsv $EXPR_TSV \
        --out-h5 $INDICES_H5 \
        --low-frac $LOW_FRAC --high-frac $HIGH_FRAC \
        --n-bootstrap $N_BOOTSTRAP --bootstrap-frac $BOOTSTRAP_FRAC \
        --seed $SEED"
fi

echo "[Stage 1] Submitting bootstrap index job..."
STAGE1_JID=$(qsub \
    -N bootstrap_indices \
    -pe parallel 4 \
    -l h_vmem=8G \
    -l h_rt=1:0:0 \
    -q long.q \
    -j y \
    -m beas \
    -M yutao.cao@tuebingen.mpg.de \
    -o logs/bootstrap_indices.\$JOB_ID.out \
    -e logs/bootstrap_indices.\$JOB_ID.err \
    -S /bin/bash \
    -cwd \
    -b y "$STAGE1_CMD" | awk '{print $3}')

echo "  Stage 1 JID = $STAGE1_JID"

# ---------------------------------------------------------------------------
# Stage 2 – per-gene correlation array (holds on Stage 1)
# ---------------------------------------------------------------------------
if [[ "$TOY" == true ]]; then
    STAGE2_CMD="python src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py \
        --toy \
        --gene-id \$((SGE_TASK_ID - 1)) \
        --indices-h5 $INDICES_H5 \
        --out-dir $CORR_DIR"
else
    STAGE2_CMD="python src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py \
        --gene-id \$((SGE_TASK_ID - 1)) \
        --indices-h5 $INDICES_H5 \
        --expr-tsv $EXPR_TSV \
        --out-dir $CORR_DIR"
fi

# -t 1-${N_GENES}: Job array, which is used to Creates N_GENES separate tasks, numbered from 1 to N_GENES.
# -l h_vmem=16G: Hard memory limit: Requests 16 GB of virtual memory per task.
# -l h_rt=10:0:0: Hard runtime limit: Maximum 10 hours per task.
# -q long.q: Queue name
# -j y: Merge stdout and stderr: Combines both output and error streams into the .out file.
# -S /bin/bash: Shell to use: Run the job using bash.

# | awk '{print $3}': awk '{print $3}' extracts the job ID from the output info after submission

echo "[Stage 2] Submitting per-gene correlation array (1–${N_GENES}, gene IDs 0–$((N_GENES - 1)))..."
STAGE2_JID=$(qsub \
    -N corr_array \
    -pe parallel 4 \
    -t 1-${N_GENES} \
    -l h_vmem=16G \
    -l h_rt=10:0:0 \
    -q long.q \
    -j y \
    -m beas \
    -M yutao.cao@tuebingen.mpg.de \
    -o logs/corr_array.\$JOB_ID.\$SGE_TASK_ID.out \
    -e logs/corr_array.\$JOB_ID.\$SGE_TASK_ID.err \
    -S /bin/bash \
    -cwd \
    -hold_jid "$STAGE1_JID" \
    -b y "$STAGE2_CMD" | awk '{print $3}')

echo "  Stage 2 JID = $STAGE2_JID"
echo ""
echo "Monitor with:  qstat -u $(whoami)"
echo "Done when all tasks finish.  Outputs in $CORR_DIR/"
