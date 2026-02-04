#!/bin/bash
# =============================================================================
# run_bootstrap_pipeline.sh
# Three-stage SGE submission for the bootstrap co-expression pipeline.
#
# Stage 0 (optional)     - convert TSV to HDF5               -> expression.h5
# Stage 1 (single job)   - generate shared bootstrap indices -> bootstrap_indices.h5
# Stage 2 (job array)    - per-gene Spearman correlations    -> corr/gene_XXXXX.h5
#
# Usage
# -----
#   # Real data (automatic HDF5 conversion from TSV)
#   bash src/SGE_scripts/run_bootstrap_pipeline.sh \
#       --expr-tsv   dataset/voomdataCtrl.txt \
#       --out-dir    results \
#       --n-bootstrap 50 --seed 42
#
#   # Use existing HDF5 (skip conversion)
#   bash src/SGE_scripts/run_bootstrap_pipeline.sh \
#       --expr-h5    dataset/expression.h5 \
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
#   auto-detected from the expression input file (row count minus 1 header).
# * If --expr-tsv is provided, Stage 0 converts it to HDF5 first (~100x faster loading).
# * If --expr-h5 is provided, Stage 0 is skipped entirely.
# * Stage 1 holds on Stage 0 (if run); Stage 2 holds on Stage 1.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# defaults
# ---------------------------------------------------------------------------
EXPR_TSV=""
EXPR_H5=""
OUT_DIR="results"
LOW_FRAC=0.2
HIGH_FRAC=0.2
N_BOOTSTRAP=50
BOOTSTRAP_FRAC=0.8
SEED=42
TOY=false
SKIP_PREPROCESS=false  # Set to true if --expr-h5 is provided

# ---------------------------------------------------------------------------
# parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --expr-tsv)       EXPR_TSV="$2";        shift 2 ;;
        --expr-h5)        EXPR_H5="$2";         shift 2 ;;  # NEW: Accept pre-converted HDF5
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
# validate inputs and determine workflow
# ---------------------------------------------------------------------------
if [[ "$TOY" == false && -z "$EXPR_TSV" && -z "$EXPR_H5" ]]; then
    echo "ERROR: provide --expr-tsv, --expr-h5, or use --toy." >&2
    exit 1
fi

# Check for conflicting inputs
if [[ -n "$EXPR_TSV" && -n "$EXPR_H5" ]]; then
    echo "ERROR: provide either --expr-tsv or --expr-h5, not both." >&2
    exit 1
fi

# Determine if we need Stage 0 preprocessing
if [[ -n "$EXPR_H5" ]]; then
    SKIP_PREPROCESS=true
    EXPR_H5_FINAL="$EXPR_H5"
elif [[ -n "$EXPR_TSV" ]]; then
    SKIP_PREPROCESS=false
    EXPR_H5_FINAL="${OUT_DIR}/expression.h5"
fi

INDICES_H5="${OUT_DIR}/bootstrap_indices.h5"
CORR_DIR="${OUT_DIR}/corr"

# ---------------------------------------------------------------------------
# detect gene count
# ---------------------------------------------------------------------------
if [[ "$TOY" == true ]]; then
    N_GENES=5
elif [[ -n "$EXPR_H5" ]]; then
    # For HDF5 input, we need to query the file for gene count
    # Using h5ls or python to get the shape
    N_GENES=$(python3 -c "import h5py; print(h5py.File('$EXPR_H5', 'r')['expr'].shape[0])")
else
    # TSV: first column is gene ID; every subsequent line is one gene
    # subtract 1 for the header row
    N_GENES=$(( $(wc -l < "$EXPR_TSV") - 1 ))
fi
# SGE -t is 1-based inclusive; task IDs run 1..N_GENES

echo "=============================================="
echo "  Bootstrap co-expression pipeline"
echo "  toy=$TOY  n_genes=$N_GENES  n_bootstrap=$N_BOOTSTRAP  seed=$SEED"
if [[ "$SKIP_PREPROCESS" == true ]]; then
    echo "  Using existing HDF5: $EXPR_H5_FINAL"
else
    echo "  Converting TSV -> HDF5: $EXPR_TSV -> $EXPR_H5_FINAL"
fi
echo "  indices  -> $INDICES_H5"
echo "  corr     -> $CORR_DIR/"
echo "=============================================="

# ---------------------------------------------------------------------------
# Stage 0 - TSV to HDF5 conversion (optional, only if --expr-tsv provided)
# ---------------------------------------------------------------------------
STAGE0_JID=""
HOLD_STAGE0=""

if [[ "$SKIP_PREPROCESS" == false && "$TOY" == false ]]; then
    echo "[Stage 0] Submitting TSV->HDF5 conversion job..."
    
    STAGE0_CMD="python src/scripts/00preprocess/00convert_expr_to_hdf5.py \
        --expr-tsv $EXPR_TSV \
        --out-h5 $EXPR_H5_FINAL \
        --compression gzip \
        --compression-level 4"
    
    STAGE0_JID=$(qsub \
        -N expr_to_hdf5 \
        -pe parallel 2 \
        -l h_vmem=16G \
        -l h_rt=2:0:0 \
        -q long.q \
        -j y \
        -m beas \
        -M yutao.cao@tuebingen.mpg.de \
        -o logs/expr_to_hdf5.\$JOB_ID.out \
        -e logs/expr_to_hdf5.\$JOB_ID.err \
        -S /bin/bash \
        -cwd \
        -b y "$STAGE0_CMD" | awk '{print $3}')
    
    echo "  Stage 0 JID = $STAGE0_JID"
    HOLD_STAGE0="-hold_jid $STAGE0_JID"
else
    echo "[Stage 0] Skipped (using existing HDF5 or toy data)"
fi

# ---------------------------------------------------------------------------
# Stage 1 - bootstrap indices (single job, holds on Stage 0 if applicable)
# ---------------------------------------------------------------------------
if [[ "$TOY" == true ]]; then
    STAGE1_CMD="python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
        --toy \
        --out-h5 $INDICES_H5 \
        --low-frac $LOW_FRAC --high-frac $HIGH_FRAC \
        --n-bootstrap $N_BOOTSTRAP --bootstrap-frac $BOOTSTRAP_FRAC \
        --seed $SEED"
else
    # Use HDF5 input (either from Stage 0 or user-provided)
    STAGE1_CMD="python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
        --in-h5 $EXPR_H5_FINAL \
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
    $HOLD_STAGE0 \
    -b y "$STAGE1_CMD" | awk '{print $3}')

echo "  Stage 1 JID = $STAGE1_JID"

# ---------------------------------------------------------------------------
# Stage 2 - per-gene correlation array (holds on Stage 1)
# ---------------------------------------------------------------------------
if [[ "$TOY" == true ]]; then
    STAGE2_CMD="python src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py \
        --toy \
        --gene-id \$((SGE_TASK_ID - 1)) \
        --indices-h5 $INDICES_H5 \
        --out-dir $CORR_DIR"
else
    # Use HDF5 expression file
    STAGE2_CMD="python src/scripts/10spearman_corr/02calc_corr_edge_bootstrap_corr.py \
        --gene-id \$((SGE_TASK_ID - 1)) \
        --indices-h5 $INDICES_H5 \
        --expr-h5 $EXPR_H5_FINAL \
        --out-dir $CORR_DIR"
fi

# -t 1-${N_GENES}: Job array, creates N_GENES separate tasks, numbered from 1 to N_GENES.
# -l h_vmem=16G: Hard memory limit: Requests 16 GB of virtual memory per task.
# -l h_rt=10:0:0: Hard runtime limit: Maximum 10 hours per task.
# -q long.q: Queue name
# -j y: Merge stdout and stderr: Combines both output and error streams into the .out file.
# -S /bin/bash: Shell to use: Run the job using bash.
# -hold_jid: Wait for the specified job to complete before starting

# | awk '{print $3}': awk '{print $3}' extracts the job ID from the output info after submission

echo "[Stage 2] Submitting per-gene correlation array (1-${N_GENES}, gene IDs 0-$((N_GENES - 1)))..."
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
echo "Pipeline submitted successfully!"
echo ""
echo "Job summary:"
if [[ -n "$STAGE0_JID" ]]; then
    echo "  Stage 0 (TSV->HDF5):  Job $STAGE0_JID"
fi
echo "  Stage 1 (indices):    Job $STAGE1_JID"
echo "  Stage 2 (corr array): Job $STAGE2_JID (${N_GENES} tasks)"
echo ""
echo "Monitor with:  qstat -u $(whoami)"
echo "View logs in:  logs/"
echo "Done when all tasks finish. Outputs in $CORR_DIR/"
