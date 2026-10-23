#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
# submit_pearson_repro.sh — Slim whole-cell Pearson pipeline (GFP-R3HDM4)
# ═══════════════════════════════════════════════════════════════════════════
#
# Reproduces the Fig 7D swarm plot (pearson_L29_swarm_DMSO_vs_PT179_23h.pdf).
#
#   Step 1: Cellpose-SAM cell segmentation (segment_cells.py)          [GPU array]
#   Step 2: Per-cell whole-cell Pearson (wholecell_pearson.py)         [CPU array]
#   Step 3: Swarm plot + Mann-Whitney stats (plot_pearson_swarm.py)    [CPU single]
#
# No ER masking / Manders: whole-cell Pearson does not depend on them.
#
# Usage:
#   bash submit_pearson_repro.sh [--start-step N] [--stop-step N] [--seg-out DIR] [--dry-run]
#
# Default starts at step 2, reusing the existing masks in
# ${TLSD_CONFOCAL_DATA}/output/cell_masks. To re-segment without overwriting
# the original masks, pass e.g.
#   --start-step 1 --seg-out /path/to/new_output
#
# Environment:
#   TLSD_CONFOCAL_DATA    experiment folder (Denoised/, output/)        [required]
#   TLSD_OUT              output root; results go to confocal_colocalization/ [required]
#   TLSD_PYTHON           interpreter of the nd2_coloc env (default: python3)
#   TLSD_GPU_PARTITION    GPU partition for step 1                      [step 1]
#   TLSD_SLURM_PARTITION  CPU partition for steps 2-3 (default: cluster default)
#   TLSD_SLURM_ACCOUNT    Slurm account (default: none)
#   TLSD_JOB_PREFIX       job-name prefix (default: ${USER}_tlsd_confocal)
# ═══════════════════════════════════════════════════════════════════════════

set -euo pipefail

REPRO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENT_DIR="${TLSD_CONFOCAL_DATA:?set TLSD_CONFOCAL_DATA to the confocal experiment folder}"
OUT_DIR="${TLSD_OUT:?set TLSD_OUT to the output root}/confocal_colocalization"

START_STEP=2
STOP_STEP=3
DRY_RUN=false
SEG_OUT="${EXPERIMENT_DIR}/output"   # holds cell_masks/ and channels/

while [[ $# -gt 0 ]]; do
    case "$1" in
        --start-step) START_STEP="$2"; shift 2 ;;
        --stop-step)  STOP_STEP="$2"; shift 2 ;;
        --seg-out)    SEG_OUT="$2"; shift 2 ;;
        --dry-run)    DRY_RUN=true; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════

INPUT_DIR="${EXPERIMENT_DIR}/Denoised"      # NIS Batch Denoise ND2s
PATTERN="*GFP-R3HDM4*"                      # 260402 acquisition (10 images)
SCRIPTS_DIR="${REPRO_DIR}"
PER_CELL_DIR="${OUT_DIR}/per_cell"
FIG_DIR="${OUT_DIR}/figures"
LOG_DIR="${OUT_DIR}/logs"

# SLURM (site-specific; set through the environment)
ACCOUNT="${TLSD_SLURM_ACCOUNT:-}"
GPU_PARTITION="${TLSD_GPU_PARTITION:-}"
CPU_PARTITION="${TLSD_SLURM_PARTITION:-}"
PY="${TLSD_PYTHON:-python3}"                # nd2_coloc env (envs/nd2_coloc.yml)
JOB_PREFIX="${TLSD_JOB_PREFIX:-${USER}_tlsd_confocal}"

# Channels (ND2 order: 488, 545, 640 nm)
CHANNELS="0,1,2"
CHANNEL_NAMES="GFP-R3HDM4,Sec61B,L29"      # GFP, JF552, JFX650

# Segmentation (values used for the original figure; see cell_masks/*_seg_meta.json)
CELLPOSE_DIAM=260
MIN_CELL_VOL=500
FLOW_THRESHOLD=0.4
CELLPROB_THRESHOLD=0.0
STITCH_THRESHOLD=0.5

# Quantification
MIN_CELL_VOXELS=500
OUTLIER_CLIP="0.1,99.9"

# Plot
LINE="L29"
TIMEPOINT="23h"
CONTROL="DMSO"
TREATMENT="PT179"

# ═══════════════════════════════════════════════════════════════════════════

if [ "$START_STEP" -le 1 ] && [ "$STOP_STEP" -ge 1 ] && [ -z "$GPU_PARTITION" ]; then
    echo "ERROR: step 1 needs a GPU partition: export TLSD_GPU_PARTITION=<gpu partition>"
    exit 1
fi
ACCOUNT_LINE=""
[ -n "$ACCOUNT" ] && ACCOUNT_LINE="#SBATCH --account=${ACCOUNT}"
CPU_PARTITION_LINE=""
[ -n "$CPU_PARTITION" ] && CPU_PARTITION_LINE="#SBATCH --partition=${CPU_PARTITION}"
# Job preamble: the env's own lib/ first (CUDA/torch shared libraries)
ENV_SETUP="export LD_LIBRARY_PATH=\"\$(\"${PY}\" -c 'import sys; print(sys.prefix)')/lib:\${LD_LIBRARY_PATH:-}\""

mkdir -p "$LOG_DIR" "$PER_CELL_DIR" "$FIG_DIR"

N_FILES=$(ls -1 "${INPUT_DIR}"/${PATTERN}.nd2 2>/dev/null | wc -l)
if [ "$N_FILES" -eq 0 ]; then
    echo "ERROR: No ${PATTERN}.nd2 files in ${INPUT_DIR}"
    exit 1
fi
LAST_IDX=$((N_FILES - 1))

echo "============================================================"
echo "  Whole-cell Pearson pipeline"
echo "  Input:     ${INPUT_DIR}/${PATTERN}.nd2 (${N_FILES} files)"
echo "  Seg out:   ${SEG_OUT}"
echo "  Channels:  ${CHANNELS} (${CHANNEL_NAMES})"
echo "  Steps:     ${START_STEP}-${STOP_STEP}"
echo "============================================================"

submit_or_echo() {
    if $DRY_RUN; then
        echo "[DRY RUN] Would submit:" >&2
        echo "$1" >&2
        echo "---" >&2
        echo "dry_run_placeholder"
    else
        echo "$1" | sbatch --parsable
    fi
}

# ── Step 1: Cellpose-SAM segmentation ──
JOB1=""
if [ "$START_STEP" -le 1 ] && [ "$STOP_STEP" -ge 1 ]; then
    JOB1=$(submit_or_echo "#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_1_segment
#SBATCH --output=${LOG_DIR}/segment_%A_%a.out
#SBATCH --error=${LOG_DIR}/segment_%A_%a.err
#SBATCH --array=0-${LAST_IDX}%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=${GPU_PARTITION}
${ACCOUNT_LINE}

${ENV_SETUP}

\"${PY}\" ${SCRIPTS_DIR}/segment_cells.py \\
    --input \"${INPUT_DIR}\" \\
    --output \"${SEG_OUT}\" \\
    --pattern \"${PATTERN}.nd2\" \\
    --index \${SLURM_ARRAY_TASK_ID} \\
    --channels ${CHANNELS} \\
    --cellpose-diam ${CELLPOSE_DIAM} \\
    --min-cell-vol ${MIN_CELL_VOL} \\
    --flow-threshold ${FLOW_THRESHOLD} \\
    --cellprob-threshold ${CELLPROB_THRESHOLD} \\
    --stitch-threshold ${STITCH_THRESHOLD}
")
    echo "  Step 1 (segment):  job ${JOB1}"
fi

# ── Step 2: Whole-cell Pearson ──
JOB2=""
if [ "$START_STEP" -le 2 ] && [ "$STOP_STEP" -ge 2 ]; then
    DEP=""
    [ -n "$JOB1" ] && DEP="#SBATCH --dependency=afterok:${JOB1}"
    JOB2=$(submit_or_echo "#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_2_pearson
#SBATCH --output=${LOG_DIR}/pearson_%A_%a.out
#SBATCH --error=${LOG_DIR}/pearson_%A_%a.err
#SBATCH --array=0-${LAST_IDX}%10
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
${CPU_PARTITION_LINE}
${ACCOUNT_LINE}
${DEP}

${ENV_SETUP}

\"${PY}\" ${REPRO_DIR}/wholecell_pearson.py \\
    --input \"${SEG_OUT}\" \\
    --output \"${PER_CELL_DIR}\" \\
    --pattern \"${PATTERN}\" \\
    --index \${SLURM_ARRAY_TASK_ID} \\
    --channels ${CHANNELS} \\
    --channel-names ${CHANNEL_NAMES} \\
    --min-cell-voxels ${MIN_CELL_VOXELS} \\
    --outlier-clip ${OUTLIER_CLIP}
")
    echo "  Step 2 (pearson):  job ${JOB2}"
fi

# ── Step 3: Plot ──
JOB3=""
if [ "$START_STEP" -le 3 ] && [ "$STOP_STEP" -ge 3 ]; then
    DEP=""
    [ -n "$JOB2" ] && DEP="#SBATCH --dependency=afterok:${JOB2}"
    JOB3=$(submit_or_echo "#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_3_plot
#SBATCH --output=${LOG_DIR}/plot_%j.out
#SBATCH --error=${LOG_DIR}/plot_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
${CPU_PARTITION_LINE}
${ACCOUNT_LINE}
${DEP}

${ENV_SETUP}

\"${PY}\" ${REPRO_DIR}/plot_pearson_swarm.py \\
    --input \"${PER_CELL_DIR}\" \\
    --output \"${FIG_DIR}\" \\
    --line ${LINE} --timepoint ${TIMEPOINT} \\
    --control ${CONTROL} --treatment ${TREATMENT}
")
    echo "  Step 3 (plot):     job ${JOB3}"
fi

echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Logs:    ${LOG_DIR}/"
