#!/bin/bash
#SBATCH --job-name=tlsd_focaldepth
#SBATCH --output=logs/FOCALDEPTH_%j.log
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
# Empirical focal-depth fit + dz sweep. Measured 1.48 GB (job 10938892,
# 100 bootstrap draws over the 60S_bothgeno cache). 8G at ~5x headroom
# for the larger claim-2 caches.
set -euo pipefail
CODE_DIRECTORY=${TLSD_REPO:?set TLSD_REPO to the repository root}/single_molecule/pipeline
export PYTHONPATH="$CODE_DIRECTORY:${PYTHONPATH:-}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
"${TLSD_PYTHON:-python3}" -u \
    "$CODE_DIRECTORY/calibrate_focal_depth.py" "$@"
