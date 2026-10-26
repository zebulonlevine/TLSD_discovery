#!/bin/bash
## call_combinationwise_operations.slurm
#SBATCH --job-name=saspt_combo
#SBATCH --output=logs/COMBINATIONWISE_output_%A_%a.log
##SBATCH --error=logs/error_%A_%a.log
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4             # number of cores/threads requested.

# Access environment variables exported from Python
echo "Settings file: $PKL_PATH"
echo "Combination:   $COMBINATION"
echo "Code from:     $CODE_DIRECTORY"
echo "Force:         ${FORCE:-0}"

EXTRA_ARGS=""
if [[ "${FORCE:-0}" == "1" ]]; then
    EXTRA_ARGS="--force"
fi
# SKIP_ROLLING=1 drops the rolling-window saSPT fit (a second full fit over the
# pooled CSV). No paper-figure output reads rollingMLE; only the dwell and
# transition analyses do.
if [[ "${SKIP_ROLLING:-0}" == "1" ]]; then
    EXTRA_ARGS="$EXTRA_ARGS --skip-rolling-windows"
fi

# Example: pass to your processing script
"${TLSD_PYTHON:-python3}" $CODE_DIRECTORY/run_combinationwise.py --combination "$COMBINATION" --settingsPKL "$PKL_PATH" $EXTRA_ARGS
