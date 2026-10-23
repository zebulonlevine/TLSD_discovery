#!/usr/bin/env bash
# Orchestrator for the Morgens 2017 ncomms15178 ricin screen reanalysis.
#
# Pipeline:
#   00_prep_counts.py        xlsx -> MAGeCK count.txt
#   01_run_mageck.sh         mageck test (median normalization + sgRNA LFCs)
#   02_gene_scores.py        MWU on top-N sgRNAs vs NTC pool -> gene-level table
#   03_waterfall.py          waterfall plot with --highlight controls
#
# Edit CONDA_ENV below to match the env where mageck is installed.
# Edit the highlight blocks at the bottom to control the plot.
#
# Usage:  bash run_all.sh       (with the mageck env, envs/mageck.yml, active)
# Reads  ${TLSD_DATA}/ricin_data_reanalysis/ncomms15178-s4.xlsx  (Morgens 2017 Supplementary Data 4)
# Writes ${TLSD_OUT}/ricin_data_reanalysis/results/

set -euo pipefail
cd "$(dirname "$0")"          # scripts and gene lists are referenced relative to this folder
REPO=$(cd .. && pwd)
PY=${PY:-python}
root() { "$PY" -c "import sys; sys.path.insert(0, '$REPO'); from common import paths; print(paths.$1())"; }
DATA=$(root data_root); OUT=$(root out_root)

# ----------------------------------------------------------------------------
# Configure
# ----------------------------------------------------------------------------
CONDA_ENV="mageck"                                # <-- EDIT this if your env name differs
XLSX="${DATA}/ricin_data_reanalysis/ncomms15178-s4.xlsx"   # Morgens 2017 supplementary
OUTDIR="${OUT}/ricin_data_reanalysis/results"
PREFIX="${OUTDIR}/K562_vs_Ricin"
TOPN=3                                            # top-N sgRNAs per gene for MWU
NTC_POOL=2000                                     # subsampled NTC pool size for MWU
FDR=0.05

# Replicate IDs as they appear in the count file (set by 00_prep_counts.py)
CTRL_IDS="K562_rep1,K562_rep2"
TRT_IDS="Ricin_rep1,Ricin_rep2"

mkdir -p "$OUTDIR"

# ----------------------------------------------------------------------------
# 00. Prep MAGeCK count table from the Morgens xlsx
# ----------------------------------------------------------------------------
if [[ ! -f "${OUTDIR}/counts.txt" ]]; then
    echo "[$(date +%H:%M:%S)] Step 00: building count table from $XLSX"
    "$PY" 00_prep_counts.py \
        --xlsx "$XLSX" \
        --ctrl-sheets "K562 replicate 1" "K562 replicate 2" \
        --trt-sheets  "Ricin replicate 1" "Ricin replicate 2" \
        --ctrl-names  K562_rep1 K562_rep2 \
        --trt-names   Ricin_rep1 Ricin_rep2 \
        --output      "${OUTDIR}/counts.txt"
else
    echo "[$(date +%H:%M:%S)] Step 00: ${OUTDIR}/counts.txt exists, skipping"
fi

# ----------------------------------------------------------------------------
# 01. Run MAGeCK test (gives median-normalized sgRNA LFCs)
# ----------------------------------------------------------------------------
if [[ ! -f "${PREFIX}.sgrna_summary.txt" ]]; then
    echo "[$(date +%H:%M:%S)] Step 01: running MAGeCK in conda env '$CONDA_ENV'"
    bash 01_run_mageck.sh \
        "${OUTDIR}/counts.txt" \
        "$PREFIX" \
        "$CTRL_IDS" \
        "$TRT_IDS" \
        "$CONDA_ENV"
else
    echo "[$(date +%H:%M:%S)] Step 01: ${PREFIX}.sgrna_summary.txt exists, skipping"
fi

# ----------------------------------------------------------------------------
# 02. Gene-level MWU pipeline (epsilon, p, product, FDR)
# ----------------------------------------------------------------------------
echo "[$(date +%H:%M:%S)] Step 02: gene-level MWU scoring"
"$PY" 02_gene_scores.py \
    --sgrna-summary "${PREFIX}.sgrna_summary.txt" \
    --output-prefix "${PREFIX}" \
    --top-n "$TOPN" \
    --ntc-pool-size "$NTC_POOL" \
    --fdr "$FDR"

# ----------------------------------------------------------------------------
# 03. Waterfall plot.  Add/remove --highlight blocks freely.
# ----------------------------------------------------------------------------
echo "[$(date +%H:%M:%S)] Step 03: waterfall plot"
"$PY" 03_waterfall.py \
    --gene-scores "${PREFIX}_all_genes.csv" \
    --rank-by epsilon \
    --output "${PREFIX}_waterfall.pdf" \
    --background-known \
    --highlight 'label=R3HDM4;color=#ffcc00;size=450;marker=*;genes=R3HDM4;annotate=true' \
    --highlight 'label=60S ribosomal proteins;color=#9467bd;size=45;genes=@rpl_genes.txt' \
    --highlight 'label=Ricin trafficking;color=#e41a1c;size=35;genes=@ricin_trafficking.txt' \
    --highlight 'label=mTORC1 negative regulators;color=#377eb8;size=35;genes=NPRL2,NPRL3,DEPDC5,TSC1,TSC2,SZT2,KPTN,ITFG2,C12orf66'

echo "[$(date +%H:%M:%S)] Done. Outputs in $OUTDIR/"
