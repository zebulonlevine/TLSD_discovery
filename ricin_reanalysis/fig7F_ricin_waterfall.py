#!/usr/bin/env python3
"""Fig 7F: waterfall of gene effects in the ricin screen (Morgens et al. 2017,
K562 ricin vs untreated, re-analysed with run_all.sh), with TLSD1/R3HDM4 highlighted.

From run_waterfall.ipynb cells 2 (setup) and 9 (plot). Cell 9 as saved passed
`size=20`, an argument plot_waterfall() has never had (TypeError). The published
PDF was rendered before that line was added, so it is dropped here. Two unused
gene lists in the same cell are dropped too.

Reads   ${TLSD_DATA}/ricin_data_reanalysis/results/K562_vs_Ricin_all_genes.csv  (run_all.sh)
Writes  ${TLSD_OUT}/ricin_data_reanalysis/waterfall_ricin_R3HDM4.pdf
Run with envs/mageck.yml or envs/crispri-screen.yml:  python fig7F_ricin_waterfall.py
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from common import paths  # noqa: E402
from waterfall import plot_waterfall, Highlight  # noqa: E402

GENE_SCORES = paths.inp('ricin_data_reanalysis/results/K562_vs_Ricin_all_genes.csv')

fig, ax = plot_waterfall(
    gene_scores=GENE_SCORES,
    rank_by="epsilon",
    background_known=False,
    figsize=(6.5, 3.25),
    highlights=[
        Highlight(label="R3HDM4", color="#2ca02c", size=450, marker=".",
                  genes=["R3HDM4"], annotate=True, zorder=100),
    ],
)
plt.savefig(paths.outp('ricin_data_reanalysis/waterfall_ricin_R3HDM4.pdf', mkdir=True))
