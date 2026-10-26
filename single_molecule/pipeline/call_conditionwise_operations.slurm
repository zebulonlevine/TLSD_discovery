#!/bin/bash
## call_filewise_operations.slurm
#SBATCH --job-name=saspt_condition
#SBATCH --output=logs/CONDITIONWISE_output_%A_%a.log
##SBATCH --error=logs/error_%A_%a.log
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4             # number of cores/threads requested.

# Access environment variables exported from Python

echo "Settings file: $PKL_PATH"
echo "Condition:     $CONDITION"
echo "Code from:     $CODE_DIRECTORY"

# Example: pass to your processing script
"${TLSD_PYTHON:-python3}" $CODE_DIRECTORY/run_conditionwise.py --condition "$CONDITION" --settingsPKL "$PKL_PATH"
