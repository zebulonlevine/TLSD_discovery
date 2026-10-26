#!/bin/bash
#SBATCH --job-name=rezone
#SBATCH --output=logs/REZONE_%A_%a.log
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
PY="${TLSD_PYTHON:-python3}"
ER=${TLSD_REPO:?set TLSD_REPO to the repository root}/single_molecule/er
ANALYSIS=${TLSD_SPT_DATA:?set TLSD_SPT_DATA}/all_data/2025_09_15_zeb_ER/tracking_output_q=15p0
ERDIR=${TLSD_SPT_DATA:?set TLSD_SPT_DATA}/all_data/2025_09_15_zeb_ER/preprocessed_corrected
$PY $ER/rerun_er_zoning.py --analysis-dir "$ANALYSIS" --er-dir "$ERDIR" \
    --pct-floor 60 --require-broad false --min-frames 2 --fov-index $SLURM_ARRAY_TASK_ID
