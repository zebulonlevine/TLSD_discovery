#!/bin/bash
## call_megapool_operations.slurm
#SBATCH --job-name=saspt_megapool
#SBATCH --output=logs/MEGAPOOL_output_%j.log
##SBATCH --error=logs/error_%j.log
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# Required exports (via sbatch --export):
#   MEGAPOOL       — megapool name (must appear in megapools.yaml)
#   PKL_PATH       — path to settings.pkl
#   CODE_DIRECTORY — code dir containing run_megapool.py
# Optional:
#   FORCE          — set to "1" to force re-pooling/re-saSPT

echo "Settings file: $PKL_PATH"
echo "Megapool:      $MEGAPOOL"
echo "Code from:     $CODE_DIRECTORY"
echo "Force:         ${FORCE:-0}"

EXTRA_ARGS=""
if [[ "${FORCE:-0}" == "1" ]]; then
    EXTRA_ARGS="--force"
fi

"${TLSD_PYTHON:-python3}" $CODE_DIRECTORY/run_megapool.py --megapool "$MEGAPOOL" --settingsPKL "$PKL_PATH" $EXTRA_ARGS
