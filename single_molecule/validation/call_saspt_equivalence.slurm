#!/bin/bash
#SBATCH --job-name=tlsd_saspt_equiv
#SBATCH --output=logs/SASPT_EQUIV_%j.log
#SBATCH --mem=12G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
# Numerical equivalence suite: real saspt 0.4.0 vs the Fork-B VBEM port.
# Memory: measured 1.87 GB at n_traj=2000 (job 10938888). Dominated by the
# (n_D, n_loc, n_sub) float64 arrays saspt's own infer_posterior
# materializes, so it scales ~linearly in n_traj: 20000 -> ~8 GB peak.
# 12G at ~1.5x headroom.
set -euo pipefail
CODE_DIRECTORY=${TLSD_REPO:?set TLSD_REPO to the repository root}/single_molecule/pipeline
export PYTHONPATH="$CODE_DIRECTORY:${PYTHONPATH:-}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
"${TLSD_PYTHON:-python3}" -u \
    "$CODE_DIRECTORY/test_saspt_equivalence.py" "$@"
