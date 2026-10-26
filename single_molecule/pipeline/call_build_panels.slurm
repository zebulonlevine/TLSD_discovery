#!/bin/bash
#SBATCH --job-name=tlsd_buildpanels
#SBATCH --output=logs/BUILDPANELS_%j.log
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
# Rebuild paper SPT panels at saSPT's own fit settings.
# Memory: dominated by the float32 likelihood cache held in RAM
# (60S_bothgeno: 192k subtracks x 150 x 5 x 4 B ~ 0.6 GB, plus the decompressed
# npz copy during load). 24G is a first untested request; re-size from `seff`.
set -euo pipefail
FIGDIR=${TLSD_REPO:?set TLSD_REPO to the repository root}/single_molecule/figures
CODE_DIRECTORY=${TLSD_REPO:?set TLSD_REPO to the repository root}/single_molecule/pipeline
export PYTHONPATH="$FIGDIR:$CODE_DIRECTORY:${PYTHONPATH:-}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
SCRIPT="$1"; shift
"${TLSD_PYTHON:-python3}" -u "$FIGDIR/$SCRIPT" "$@"
