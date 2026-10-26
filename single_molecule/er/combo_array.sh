#!/bin/bash
#SBATCH --job-name=cwER
#SBATCH --output=logs/CW_%A_%a.log
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
PY="${TLSD_PYTHON:-python3}"
ER=${TLSD_REPO:?set TLSD_REPO to the repository root}/single_molecule/er
CODE=${TLSD_REPO:?set TLSD_REPO to the repository root}/single_molecule/pipeline
ANALYSIS=${TLSD_SPT_DATA:?set TLSD_SPT_DATA}/all_data/2025_09_15_zeb_ER/tracking_output_q=15p0
export PYTHONPATH=$CODE:$PYTHONPATH
PKL=$ANALYSIS/settings.pkl
COMBO=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $ANALYSIS/logs/zone_combos.txt)
echo "=== combinationwise: $COMBO  (PYTHONPATH=$PYTHONPATH) ==="
$PY $ER/run_combinationwise.py -s $PKL -c "$COMBO" --force
