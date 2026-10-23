#!/bin/bash
# Render every screen-derived panel. Needs the pipeline, scoring and GO outputs
# under the data root. Environment: envs/crispri-screen.yml.
#   srun --ntasks=1 -c2 --mem=16G bash screen_figures/run_all_figures.sh
# (--ntasks=1 matters: parallel writers corrupt PDFs.)
set -euo pipefail
cd "$(dirname "$0")/panels"
PY=${PY:-python}
export MPLBACKEND=Agg
for s in fig1_figS1_figS5A figS2B_sublibrary_selection fig2C_synth_vs_deg \
         fig2D_figS3B_category_violins figS2C_score_construction figS2D_figS7A_lfc_correlation \
         figS2E_go_scores figS2F_growth_vs_scores figS3A_subunit_scatter \
         fig6B_figS7BC_suppressor_scatter figS7EF_growth_epistasis; do
    echo "== $s"; "$PY" "$s.py"
done
