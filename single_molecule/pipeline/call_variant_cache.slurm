#!/bin/bash
#SBATCH --job-name=tlsd_varcache
#SBATCH --output=logs/VARCACHE_%j.log
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
# --mem supplied per job on the command line.
set -euo pipefail
CODE=${TLSD_REPO:?set TLSD_REPO to the repository root}/single_molecule/pipeline
export PYTHONPATH="$CODE:${PYTHONPATH:-}"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
"${TLSD_PYTHON:-python3}" -u \
  "$CODE/build_variant_cache_one.py" --dataset "$DATASET" --condition "$COND" --variant "$VARIANT"
