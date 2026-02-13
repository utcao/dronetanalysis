#!/bin/bash
# =============================================================================
# run_bootstrap_pipeline.sh
# Complete differential co-expression pipeline with SGE job submission.
#
# Pipeline Stages:
# Stage 0 (optional)  - Convert TSV to HDF5              -> expression.h5
# Stage 1 (single)    - Generate bootstrap indices       -> bootstrap_indices.h5
# Stage 2a (single)   - Compute base correlations + FDR  -> base_correlations.h5
# Stage 2b (single)   - Bootstrap significant edges only -> bootstrap_significant.h5
# Stage 3 (single)    - Reconstruct differential network -> differential_network.h5
# Stage 4 (optional)  - Collect focus gene topology      -> focus_gene_topology.h5
# Stage 5 (single)    - Prepare visualization data       -> visualization_data/
#
# Usage
# -----
#   # Real data (automatic HDF5 conversion from TSV)
#   bash src/SGE_scripts/run_bootstrap_pipeline.sh \
#       --expr-tsv   dataset/test/voomdataCtrl_test.txt \
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
#   # Local execution (no SGE, for testing)
#   bash src/SGE_scripts/run_bootstrap_pipeline.sh --toy --local --out-dir results
#
# Notes
# -----
# * Requires qsub (SGE / UGE) unless --local is specified.
# * For --toy, gene count is hard-coded to 5.
# * Stage 4 is for per-gene networks (many files) - skipped by default.
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
N_BOOTSTRAP=500
BOOTSTRAP_FRAC=0.8
FDR_ALPHA=0.9
SEED=42
TOY=false
LOCAL=false
SKIP_PREPROCESS=false
SKIP_STAGE4=true  # Stage 4 is for per-gene networks, skip by default
CALC_PER_GENE_METRICS=false
NO_CI_FILTER=false
BATCH_SIZE=""  # Batch size for Stage 1 (empty = vectorized mode)

# ---------------------------------------------------------------------------
# parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --expr-tsv)       EXPR_TSV="$2";        shift 2 ;;
        --expr-h5)        EXPR_H5="$2";         shift 2 ;;
        --out-dir)        OUT_DIR="$2";         shift 2 ;;
        --low-frac)       LOW_FRAC="$2";        shift 2 ;;
        --high-frac)      HIGH_FRAC="$2";       shift 2 ;;
        --n-bootstrap)    N_BOOTSTRAP="$2";     shift 2 ;;
        --bootstrap-frac) BOOTSTRAP_FRAC="$2";  shift 2 ;;
        --fdr-alpha)      FDR_ALPHA="$2";       shift 2 ;;
        --seed)           SEED="$2";            shift 2 ;;
        --batch-size)     BATCH_SIZE="$2";      shift 2 ;;
        --toy)            TOY=true;             shift   ;;
        --local)          LOCAL=true;           shift   ;;
        --with-stage4)    SKIP_STAGE4=false;    shift   ;;
        --calc-per-gene-metrics) CALC_PER_GENE_METRICS=true; shift ;;
        --no-ci-filter) NO_CI_FILTER=true; shift ;;
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

if [[ -n "$EXPR_TSV" && -n "$EXPR_H5" ]]; then
    echo "ERROR: provide either --expr-tsv or --expr-h5, not both." >&2
    exit 1
fi

# Determine if we need Stage 0 preprocessing
if [[ "$TOY" == true ]]; then
    # TOY mode: Stage 0 generates structured toy expression.h5
    SKIP_PREPROCESS=false
    EXPR_H5_FINAL="${OUT_DIR}/expression.h5"
elif [[ -n "$EXPR_H5" ]]; then
    SKIP_PREPROCESS=true
    EXPR_H5_FINAL="$EXPR_H5"
elif [[ -n "$EXPR_TSV" ]]; then
    SKIP_PREPROCESS=false
    EXPR_H5_FINAL="${OUT_DIR}/expression.h5"
fi

# Output file paths
INDICES_H5="${OUT_DIR}/bootstrap_indices.h5"
BASE_CORR_DIR="${OUT_DIR}/base_correlations"
BOOT_SIG_DIR="${OUT_DIR}/bootstrap_significant"
DIFF_NET_H5="${OUT_DIR}/differential_network_summary.h5"
DIFF_NET_TSV="${OUT_DIR}/rewiring_hubs.tsv"
FOCUS_TOPO_H5="${OUT_DIR}/focus_gene_topology.h5"
VIS_DIR="${OUT_DIR}/visualization_data"

# Create output directories
mkdir -p "${OUT_DIR}"
mkdir -p "${BASE_CORR_DIR}"
mkdir -p "${BOOT_SIG_DIR}"
mkdir -p logs

# ---------------------------------------------------------------------------
# detect gene count (toy mode detected after Stage 0 runs)
# ---------------------------------------------------------------------------
if [[ "$TOY" == true ]]; then
    N_GENES=10  # toy data: 10 genes, 60 samples
elif [[ -n "$EXPR_H5" ]]; then
    N_GENES=$(python3 -c "import h5py; print(h5py.File('$EXPR_H5', 'r')['expr'].shape[0])")
else
    N_GENES=$(( $(wc -l < "$EXPR_TSV") - 1 ))
fi

echo "=============================================="
echo "  Differential Co-expression Pipeline"
echo "=============================================="
echo "  Mode:        $([ "$TOY" == true ] && echo 'TOY' || echo 'REAL DATA')"
echo "  Execution:   $([ "$LOCAL" == true ] && echo 'LOCAL' || echo 'SGE')"
echo "  n_genes:     $N_GENES"
echo "  n_bootstrap: $N_BOOTSTRAP"
echo "  fdr_alpha:   $FDR_ALPHA"
echo "  seed:        $SEED"
echo ""
if [[ "$TOY" == true ]]; then
    echo "  Input:       TOY DATA (10 genes, 60 samples)"
elif [[ "$SKIP_PREPROCESS" == true ]]; then
    echo "  Input HDF5:  $EXPR_H5_FINAL"
else
    echo "  Input TSV:   $EXPR_TSV"
    echo "  Output HDF5: $EXPR_H5_FINAL"
fi
echo ""
echo "  Outputs:"
echo "    Stage 1:  $INDICES_H5"
echo "    Stage 2a: $BASE_CORR_DIR/ (per-gene)"
echo "    Stage 2b: $BOOT_SIG_DIR/ (per-gene)"
echo "    Stage 3:  $DIFF_NET_H5"
echo "    Stage 5:  $VIS_DIR/"
echo "=============================================="

# ---------------------------------------------------------------------------
# Helper function: run command locally or via SGE
# ---------------------------------------------------------------------------
run_job() {
    local job_name="$1"
    local hold_jid="$2"
    local cmd="$3"
    local vmem="${4:-8G}"
    local runtime="${5:-2:0:0}"

    if [[ "$LOCAL" == true ]]; then
        echo "[LOCAL] Running $job_name..." >&2
        if ! eval "$cmd"; then
            echo "ERROR: $job_name failed!" >&2
            exit 1
        fi
        echo "" >&2
        # Return empty string — local mode has no JID
        echo ""
    else
        local qsub_opts="-N $job_name \
            -pe parallel 4 \
            -l h_vmem=$vmem \
            -l h_rt=$runtime \
            -q long.q \
            -j y \
            -o logs/${job_name}.\$JOB_ID.out \
            -e logs/${job_name}.\$JOB_ID.err \
            -S /bin/bash \
            -cwd"

        if [[ -n "$hold_jid" ]]; then
            qsub_opts="$qsub_opts -hold_jid $hold_jid"
        fi

        local qsub_output
        qsub_output=$(qsub $qsub_opts -b y "$cmd" 2>&1)
        local jid=$(echo "$qsub_output" | awk '{print $3}')

        if [[ -z "$jid" || ! "$jid" =~ ^[0-9]+$ ]]; then
            echo "ERROR: $job_name failed to submit. qsub output:" >&2
            echo "  $qsub_output" >&2
            exit 1
        fi

        echo "  $job_name: Job $jid" >&2
        echo "$jid"
    fi
}

# Helper: run array job (SGE -t) or loop locally
run_array_job() {
    local job_name="$1"
    local hold_jid="$2"
    local cmd_template="$3"  # Must contain GENE_INDEX placeholder
    local n_tasks="$4"
    local vmem="${5:-8G}"
    local runtime="${6:-2:0:0}"

    if [[ "$LOCAL" == true ]]; then
        echo "[LOCAL] Running $job_name (${n_tasks} tasks)..." >&2
        for g in $(seq 0 $((n_tasks - 1))); do
            local cmd="${cmd_template//GENE_INDEX/$g}"
            if ! eval "$cmd"; then
                echo "ERROR: $job_name task $g failed!" >&2
                exit 1
            fi
        done
        echo "" >&2
        # Return empty string — local mode has no JID
        echo ""
    else
        local qsub_opts="-N $job_name \
            -t 1-${n_tasks} \
            -pe parallel 4 \
            -l h_vmem=$vmem \
            -l h_rt=$runtime \
            -q long.q \
            -j y \
            -o logs/${job_name}.\$JOB_ID.\$TASK_ID.out \
            -e logs/${job_name}.\$JOB_ID.\$TASK_ID.err \
            -S /bin/bash \
            -cwd"

        if [[ -n "$hold_jid" ]]; then
            qsub_opts="$qsub_opts -hold_jid $hold_jid"
        fi

        # Replace GENE_INDEX with $(($SGE_TASK_ID - 1)) for 0-based gene index
        local sge_cmd="${cmd_template//GENE_INDEX/\$((\$SGE_TASK_ID - 1))}"
        local qsub_output
        qsub_output=$(qsub $qsub_opts -b y "$sge_cmd" 2>&1)
        # Array job output: "Your job-array 12345.1-10:1 ..." — extract numeric JID before the dot
        local jid=$(echo "$qsub_output" | awk '{print $3}' | cut -d. -f1)

        if [[ -z "$jid" || ! "$jid" =~ ^[0-9]+$ ]]; then
            echo "ERROR: $job_name failed to submit. qsub output:" >&2
            echo "  $qsub_output" >&2
            exit 1
        fi

        echo "  $job_name: Array Job $jid (${n_tasks} tasks)" >&2
        echo "$jid"
    fi
}

# ---------------------------------------------------------------------------
# Stage 0 - TSV to HDF5 conversion (optional)
# ---------------------------------------------------------------------------
STAGE0_JID=""

if [[ "$TOY" == true ]]; then
    echo ""
    echo "[Stage 0] Generating toy data with known differential edges..."

    STAGE0_CMD="python src/scripts/00preprocess/00convert_expr_to_hdf5.py \
        --toy \
        --out-h5 $EXPR_H5_FINAL"

    STAGE0_JID=$(run_job "stage0_toy" "" "$STAGE0_CMD" "4G" "0:10:0")
elif [[ "$SKIP_PREPROCESS" == false ]]; then
    echo ""
    echo "[Stage 0] TSV -> HDF5 conversion..."

    STAGE0_CMD="python src/scripts/00preprocess/00convert_expr_to_hdf5.py \
        --expr-tsv $EXPR_TSV \
        --out-h5 $EXPR_H5_FINAL \
        --compression gzip \
        --compression-level 4"

    STAGE0_JID=$(run_job "stage0_convert" "" "$STAGE0_CMD" "16G" "2:0:0")
else
    echo ""
    echo "[Stage 0] Skipped (using existing HDF5)"
fi

# ---------------------------------------------------------------------------
# Stage 1 - Bootstrap indices
# ---------------------------------------------------------------------------
echo ""
echo "[Stage 1] Generate bootstrap indices..."

STAGE1_CMD="python src/scripts/01subset/01get_extreme_pop_bootstrap.py \
    --in-h5 $EXPR_H5_FINAL \
    --out-h5 $INDICES_H5 \
    --low-frac $LOW_FRAC --high-frac $HIGH_FRAC \
    --n-bootstrap $N_BOOTSTRAP --bootstrap-frac $BOOTSTRAP_FRAC \
    --seed $SEED"

# Add batch-size parameter if specified
if [[ -n "$BATCH_SIZE" ]]; then
    STAGE1_CMD="$STAGE1_CMD --batch-size $BATCH_SIZE"
    STAGE1_MEM="12G"  # Batch mode needs less memory
    echo "  (Batch mode enabled: $BATCH_SIZE genes per batch, memory: $STAGE1_MEM)"
else
    STAGE1_MEM="30G"  # Vectorized mode needs more memory
    echo "  (Vectorized mode: faster but needs more memory: $STAGE1_MEM)"
fi

STAGE1_JID=$(run_job "stage1_indices" "$STAGE0_JID" "$STAGE1_CMD" "$STAGE1_MEM" "1:0:0")

# ---------------------------------------------------------------------------
# Stage 2a - Base correlations + significance tests (per-gene array job)
# ---------------------------------------------------------------------------
echo ""
echo "[Stage 2a] Compute base correlations (${N_GENES} genes, per-gene array)..."

STAGE2A_CMD="python src/scripts/10spearman_corr/02a_calc_base_correlations.py \
    --expr-h5 $EXPR_H5_FINAL \
    --indices-h5 $INDICES_H5 \
    --out-dir $BASE_CORR_DIR \
    --gene-index GENE_INDEX \
    --fdr-alpha $FDR_ALPHA"

STAGE2A_JID=$(run_array_job "stage2a_base" "$STAGE1_JID" "$STAGE2A_CMD" "$N_GENES" "32G" "4:0:0")

# ---------------------------------------------------------------------------
# Stage 2b - Bootstrap significant edges only (per-gene array job)
# ---------------------------------------------------------------------------
echo ""
echo "[Stage 2b] Bootstrap significant edges (${N_GENES} genes, per-gene array)..."

STAGE2B_CMD="python src/scripts/10spearman_corr/02b_bootstrap_significant_edges.py \
    --expr-h5 $EXPR_H5_FINAL \
    --indices-h5 $INDICES_H5 \
    --base-dir $BASE_CORR_DIR \
    --out-dir $BOOT_SIG_DIR \
    --gene-index GENE_INDEX \
    --edge-selection sig_differential"

STAGE2B_JID=$(run_array_job "stage2b_boot" "$STAGE2A_JID" "$STAGE2B_CMD" "$N_GENES" "20G" "2:0:0")

# ---------------------------------------------------------------------------
# Stage 3 - Collect per-gene networks into summary
# ---------------------------------------------------------------------------
echo ""
echo "[Stage 3] Collect per-gene networks into summary..."

STAGE3_CMD="python src/scripts/10spearman_corr/03_reconstruct_diff_network.py \
    --base-dir $BASE_CORR_DIR \
    --boot-dir $BOOT_SIG_DIR \
    --out-h5 $DIFF_NET_H5 \
    --out-focus-tsv $DIFF_NET_TSV \
    --edge-selection sig_differential"

if [[ "$NO_CI_FILTER" == true ]]; then
    STAGE3_CMD="$STAGE3_CMD --no-ci-filter"
fi

STAGE3_JID=$(run_job "stage3_collect" "$STAGE2B_JID" "$STAGE3_CMD" "16G" "1:0:0")

# ---------------------------------------------------------------------------
# Stage 4 - Collect focus gene topology (optional, for per-gene networks)
# ---------------------------------------------------------------------------
STAGE4_JID=""
if [[ "$SKIP_STAGE4" == false ]]; then
    echo ""
    echo "[Stage 4] Collect focus gene topology..."

    STAGE4_CMD="python src/scripts/10spearman_corr/04_collect_focus_gene_topology.py \
        --network-dir ${OUT_DIR}/networks \
        --focus-genes top:50 \
        --n-genes $N_GENES \
        --out-h5 $FOCUS_TOPO_H5 \
        --n-jobs 4"

    STAGE4_JID=$(run_job "stage4_topology" "$STAGE3_JID" "$STAGE4_CMD" "8G" "2:0:0")
else
    echo ""
    echo "[Stage 4] Skipped (per-gene networks not requested, use --with-stage4 to enable)"
    STAGE4_JID="$STAGE3_JID"  # Use Stage 3 as dependency for Stage 5
fi

# ---------------------------------------------------------------------------
# Stage 5 - Prepare visualization data (TODO: update for per-gene mode)
# ---------------------------------------------------------------------------
echo ""
echo "[Stage 5] Skipped (requires per-gene visualization update)"
STAGE5_JID="$STAGE4_JID"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Pipeline submitted successfully!"
echo "=============================================="

if [[ "$LOCAL" == true ]]; then
    echo ""
    echo "Execution complete. Outputs in $OUT_DIR/"
else
    echo ""
    echo "Job summary:"
    if [[ -n "$STAGE0_JID" ]]; then
        echo "  Stage 0 (TSV->HDF5):     Job $STAGE0_JID"
    fi
    echo "  Stage 1 (indices):       Job $STAGE1_JID"
    echo "  Stage 2a (base corr):    Array Job $STAGE2A_JID (${N_GENES} tasks)"
    echo "  Stage 2b (bootstrap):    Array Job $STAGE2B_JID (${N_GENES} tasks)"
    echo "  Stage 3 (collect):       Job $STAGE3_JID"
    if [[ "$SKIP_STAGE4" == false ]]; then
        echo "  Stage 4 (topology):      Job $STAGE4_JID"
    fi
    echo "  Stage 5 (visualization): Job $STAGE5_JID"
    echo ""
    echo "Monitor with:  qstat -u $(whoami)"
    echo "View logs in:  logs/"
    echo "Outputs in:    $OUT_DIR/"
fi

echo ""
echo "Next steps after completion:"
echo "  1. Examine network summary: $DIFF_NET_H5"
echo "  2. Examine rewiring hubs:   $DIFF_NET_TSV"
echo "  3. Per-gene files:          $BASE_CORR_DIR/ and $BOOT_SIG_DIR/"
echo "  4. Run R visualization:     cd $VIS_DIR && Rscript visualize_networks.R"
