#!/bin/bash
# Per-sgRNA synthesis / degradation line scores for the 8 sublibrary screens
# (60S = L29, 40S = S3; sgNT = "wt" and sgUBA5 arms).
#
# Scripts the screen_viz.py line-score-pipeline runs that produced the published
# Sublibrary_results/scores/ (run interactively 2026-01-15 for sgNT and
# 2026-01-26 for sgUBA5). x = the ratioON comparison, y = new (synthesis, line
# y = x) or old (degradation, line y = -x); score = parallel*|parallel| /
# (|parallel| + perp), with defaults --min-count 20 --top-n 3 --fdr 0.05.
#
# Reads  ${TLSD_DATA}/Sublibrary_results/results/<comparison>/  (MAGeCK output)
# Writes ${TLSD_OUT}/Sublibrary_results/scores/<screen>_<score>[_sgUBA5]/
# The per-sgRNA file is renamed sgrna_scores.csv -> <dir>_sgrna_scores.csv. (By
# hand, two published directories got irregular names; the reader in
# final_pvalues.py accepts any *sgrna_scores*.csv, so these regular names work.)
# The gene-level files this also writes (gene_scores_*.csv) are superseded by
# export_scores.py and are not used downstream.
#
# Usage: bash run_line_scores.sh      (roots from TLSD_DATA/TLSD_OUT or config/paths.yaml)
# Heavy enough to run under Slurm:  srun --ntasks=1 -c1 --mem=4G bash run_line_scores.sh
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
PY=${PY:-python}
root() { "$PY" -c "import sys; sys.path.insert(0, '$REPO'); from common import paths; print(paths.$1())"; }
DATA=$(root data_root); OUT=$(root out_root)
VIZ="$REPO/screen_figures/lib/screen_viz.py"
RES="$DATA/Sublibrary_results/results"

run() {   # screen score cond
    local screen=$1 score=$2 cond=$3 ycomp line sub lab tag
    case $score in
        Synthesis)   ycomp=new; line="1,0";  lab="New" ;;
        Degradation) ycomp=old; line="-1,0"; lab="Old" ;;
    esac
    sub=$([ "$screen" = L29 ] && echo 60S || echo 40S)
    tag=$([ "$cond" = wt ] && echo "" || echo " sgUBA5")
    local x="${screen}_ratioON_top_v_bottom_${cond}" y="${screen}_${ycomp}_top_v_bottom_${cond}"
    local dir="${screen}_${score}$([ "$cond" = wt ] || echo _sgUBA5)"
    local o="$OUT/Sublibrary_results/scores/$dir"
    echo "== $dir"
    "$PY" "$VIZ" line-score-pipeline \
        --sgrna1 "$RES/$x/$x.sgrna_summary.txt" --sgrna2 "$RES/$y/$y.sgrna_summary.txt" \
        --results1 "$RES/$x/${x}_all_genes.csv" --results2 "$RES/$y/${y}_all_genes.csv" \
        --line="$line" --names "$sub Turnover (New/Old)$tag,$lab $sub Level$tag" \
        --fdr 0.05 -o "$o"
    mv "$o/sgrna_scores.csv" "$o/${dir}_sgrna_scores.csv"
}

for cond in wt sgUBA5; do
    for screen in L29 S3; do
        for score in Synthesis Degradation; do
            run "$screen" "$score" "$cond"
        done
    done
done
