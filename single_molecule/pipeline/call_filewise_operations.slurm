#!/bin/bash
## call_filewise_operations.slurm
#SBATCH --job-name=process_file
#SBATCH --output=logs/FILEWISE_output_%A_%a.log
##SBATCH --error=logs/error_%A_%a.log
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1             # number of cores/threads requested.

# Access environment variables exported from Python

IFS=":" read -r -a ITEMS_ARR <<< "$ITEM_STRING"
item=${ITEMS_ARR[$SLURM_ARRAY_TASK_ID]}

echo "ND2 file:      $item"
echo "Settings file: $PKL_PATH"
echo "Condition:     $CONDITION"
echo "Code from:     $CODE_DIRECTORY"

# Example: pass to your processing script
"${TLSD_PYTHON:-python3}" $CODE_DIRECTORY/run_filewise.py --condition "$CONDITION" --inputfile "$item" --settingsPKL "$PKL_PATH"
