#!/usr/bin/env python3
"""Fig 6B (degradation) and Fig S7B/C (synthesis): suppressor-screen gene scores,
sgNT (x) vs sgUBA5 (y), 60S and 40S.

Extracted from Fig5.ipynb cell(s) 1, 3, 5, 7, 10 (cell 1 = path setup). Cell code is verbatim except for paths:
inputs come from the re-scored tree via common.paths, categories from
screen_figures/categories/, and outputs go to fixed_NotI/plots/*_fixedNotI.pdf
under the output root. These are the same substitutions the original
regenerate_pRanked_figures.py harness applied when it exec'd these cells.

Run with envs/crispri-screen.yml:  python fig6B_figS7BC_suppressor_scatter.py
"""
import os
import sys

import matplotlib
if 'ipykernel' not in sys.modules:            # headless when run from the CLI
    matplotlib.use('Agg')
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'screen_figures', 'lib'))
from common import paths  # noqa: E402
import screen_viz  # noqa: E402

os.makedirs(paths.outp('fixed_NotI', 'plots'), exist_ok=True)
os.makedirs(paths.outp('fixed_NotI', 'compare_csv'), exist_ok=True)


# ---- Fig5.ipynb cell 1 --------------------------------------------------
pathDeg60S_sgUBA5 = paths.inp('fixed_NotI/scores/L29_Degradation_sgUBA5')
pathDeg40S_sgUBA5 = paths.inp('fixed_NotI/scores/S3_Degradation_sgUBA5')
pathSynth60S_sgUBA5 = paths.inp('fixed_NotI/scores/L29_Synthesis_sgUBA5')
pathSynth40S_sgUBA5 = paths.inp('fixed_NotI/scores/S3_Synthesis_sgUBA5')
pathDeg60S = paths.inp('fixed_NotI/scores/L29_Degradation')
pathDeg40S = paths.inp('fixed_NotI/scores/S3_Degradation')
pathSynth60S = paths.inp('fixed_NotI/scores/L29_Synthesis')
pathSynth40S = paths.inp('fixed_NotI/scores/S3_Synthesis')

growth40S = paths.inp('fixed_NotI/results/S3_unsorted_v_early_wt')
growth40S_sgUBA5 = paths.inp('fixed_NotI/results/S3_unsorted_v_early_sgUBA5')
growth60S = paths.inp('fixed_NotI/results/L29_unsorted_v_early_wt')
growth60S_sgUBA5 = paths.inp('fixed_NotI/results/L29_unsorted_v_early_sgUBA5')

# (pathCats is set by each plotting cell below)




# ---- Fig5.ipynb cell 3 --------------------------------------------------
comp1 = screen_viz.load_comparison(os.path.join(pathDeg60S,'gene_scores_all_genes.csv'), os.path.join(pathDeg60S,'gene_scores_hits_fdr0.05.csv'), '60S Degradation Score')
comp2 = screen_viz.load_comparison(os.path.join(pathDeg60S_sgUBA5,'gene_scores_all_genes.csv'), os.path.join(pathDeg60S_sgUBA5,'gene_scores_hits_fdr0.05.csv'), '60S Degradation Score (sgUBA5)')


pathCats = paths.repo('screen_figures', 'categories', 'CatsSublib_R3HDM4_40S.csv')
categories = screen_viz.load_categories(pathCats)

screen_viz.plot_comparison_scatter(
    comp1, comp2, paths.outp('fixed_NotI/plots/compare_60SDegscores_sgUBA5_R3HDM4cats_fixedNotI.pdf'),
    categories=categories, title = "60S Degradation vs. 60S Degradation sgUBA5",
    show_diagonal = True,show_correlation = False,figsize = (7,4),
    legend_outside = True,colors = {
        '40S Biogenesis/Proteins' :  '#74c476',
        'Autophagy/Lysosome' : '#c51b7d',
        'R3HDM4' : '#231F20',
        'UFMylation' : '#E69F00',
        'Translation initiation': '#5CB1E3'
    }, size_by_hits = False, point_size = 70,
    category_alpha={'R3HDM4': 1.0}
)
deg60Scomp = comp1.results.merge(comp2.results, suffixes = ['_sgNT', '_sgUBA5'], on=['gene_id','gene'])
deg60Scomp['diff_epsilon_sgUBA5-sgNT'] = deg60Scomp['epsilon_sgUBA5'] - deg60Scomp['epsilon_sgNT']
deg60Scomp['category'] = ''
for i in categories.category.unique():
    deg60Scomp.loc[deg60Scomp.gene.isin(categories.gene[categories.category == i]),'category'] = i
deg60Scomp.sort_values('epsilon_sgNT').to_csv(paths.outp('fixed_NotI/compare_csv/degradation60S_compare_fixedNotI.csv'))
deg60Scomp.sort_values('diff_epsilon_sgUBA5-sgNT')


# ---- Fig5.ipynb cell 5 --------------------------------------------------
comp1 = screen_viz.load_comparison(os.path.join(pathDeg40S,'gene_scores_all_genes.csv'), os.path.join(pathDeg40S,'gene_scores_hits_fdr0.05.csv'), '40S Degradation Score')
comp2 = screen_viz.load_comparison(os.path.join(pathDeg40S_sgUBA5,'gene_scores_all_genes.csv'), os.path.join(pathDeg40S_sgUBA5,'gene_scores_hits_fdr0.05.csv'), '40S Degradation Score (sgUBA5)')


pathCats = paths.repo('screen_figures', 'categories', 'CatsSublib_R3HDM4_40S.csv')
categories = screen_viz.load_categories(pathCats)

screen_viz.plot_comparison_scatter(
    comp1, comp2, paths.outp('fixed_NotI/plots/compare_40SDegscores_sgUBA5_R3HDM4cats_fixedNotI.pdf'),
    categories=categories, title = "40S Degradation vs. 40S Degradation sgUBA5",
    show_diagonal = True,show_correlation = False,figsize = (7,4),
    legend_outside = True,colors = {
        '40S Biogenesis/Proteins' :  '#74c476',
        'Autophagy/Lysosome' : '#c51b7d',
        'R3HDM4' : '#231F20',
        'UFMylation' : '#E69F00',
        'Translation initiation': '#5CB1E3'
    }, size_by_hits = False, point_size = 70,
    category_alpha={'R3HDM4': 1.0}
)
deg40Scomp = comp1.results.merge(comp2.results, suffixes = ['_sgNT', '_sgUBA5'], on=['gene_id','gene'])
deg40Scomp['diff_epsilon_sgUBA5-sgNT'] = deg40Scomp['epsilon_sgUBA5'] - deg40Scomp['epsilon_sgNT']
deg40Scomp['category'] = ''
for i in categories.category.unique():
    deg40Scomp.loc[deg40Scomp.gene.isin(categories.gene[categories.category == i]),'category'] = i
deg40Scomp.sort_values('epsilon_sgNT').to_csv(paths.outp('fixed_NotI/compare_csv/degradation40S_compare_fixedNotI.csv'))

deg40Scomp.sort_values('epsilon_sgNT')


# ---- Fig5.ipynb cell 7 --------------------------------------------------
comp1 = screen_viz.load_comparison(os.path.join(pathSynth60S,'gene_scores_all_genes.csv'), os.path.join(pathSynth60S,'gene_scores_hits_fdr0.05.csv'), '60S Synthesis Score')
comp2 = screen_viz.load_comparison(os.path.join(pathSynth60S_sgUBA5,'gene_scores_all_genes.csv'), os.path.join(pathSynth60S_sgUBA5,'gene_scores_hits_fdr0.05.csv'), '60S Synthesis Score (sgUBA5)')


pathCats = paths.repo('screen_figures', 'categories', 'CatsSublib_R3HDM4_40S_2.csv')
categories = screen_viz.load_categories(pathCats)

screen_viz.plot_comparison_scatter(
    comp1, comp2, paths.outp('fixed_NotI/plots/compare_60SSynthscores_sgUBA5_R3HDM4cats_fixedNotI.pdf'),
    categories=categories, title = "60S Synthesis vs. 60S Synthesis sgUBA5",
    show_diagonal = True,show_correlation = False,figsize = (7,4),
    legend_outside = True,colors = {
        '40S Biogenesis/Proteins' :  '#74c476',
        '60S Biogenesis/Proteins' :  '#238b45',
        'Both Biogenesis' : '#00441b',
        'Autophagy/Lysosome' : '#c51b7d',
        'R3HDM4' : '#231F20',
        'UFMylation' : '#E69F00',
        'Translation initiation': '#5CB1E3'
    }, size_by_hits = False, point_size = 70,
    category_alpha={'R3HDM4': 1.0}
)
Synth60Scomp = comp1.results.merge(comp2.results, suffixes = ['_sgNT', '_sgUBA5'], on=['gene_id','gene'])
Synth60Scomp['diff_epsilon_sgUBA5-sgNT'] = Synth60Scomp['epsilon_sgUBA5'] - Synth60Scomp['epsilon_sgNT']
Synth60Scomp['category'] = ''
for i in categories.category.unique():
    Synth60Scomp.loc[Synth60Scomp.gene.isin(categories.gene[categories.category == i]),'category'] = i
Synth60Scomp.sort_values('epsilon_sgUBA5').to_csv(paths.outp('fixed_NotI/compare_csv/synthesis60S_compare_fixedNotI.csv'))
Synth60Scomp.sort_values('epsilon_sgUBA5')


# ---- Fig5.ipynb cell 10 --------------------------------------------------
comp1 = screen_viz.load_comparison(os.path.join(pathSynth40S,'gene_scores_all_genes.csv'), os.path.join(pathSynth40S,'gene_scores_hits_fdr0.05.csv'), '40S Synthesis Score')
comp2 = screen_viz.load_comparison(os.path.join(pathSynth40S_sgUBA5,'gene_scores_all_genes.csv'), os.path.join(pathSynth40S_sgUBA5,'gene_scores_hits_fdr0.05.csv'), '40S Synthesis Score (sgUBA5)')


pathCats = paths.repo('screen_figures', 'categories', 'CatsSublib_R3HDM4_40S_2.csv')
categories = screen_viz.load_categories(pathCats)

screen_viz.plot_comparison_scatter(
    comp1, comp2, paths.outp('fixed_NotI/plots/compare_40SSynthscores_sgUBA5_R3HDM4cats_fixedNotI.pdf'),
    categories=categories, title = "40S Synthesis vs. 40S Synthesis sgUBA5",
    show_diagonal = True,show_correlation = False,figsize = (7,4),
    legend_outside = True,colors = {
        '40S Biogenesis/Proteins' :  '#74c476',
        '60S Biogenesis/Proteins' :  '#238b45',
        'Both Biogenesis' : '#00441b',
        'Autophagy/Lysosome' : '#c51b7d',
        'R3HDM4' : '#231F20',
        'UFMylation' : '#E69F00',
        'Translation initiation': '#5CB1E3'
    }, size_by_hits = False, point_size = 70,
    category_alpha={'R3HDM4': 1.0}
)
Synth40Scomp = comp1.results.merge(comp2.results, suffixes = ['_sgNT', '_sgUBA5'], on=['gene_id','gene'])
Synth40Scomp['diff_epsilon_sgUBA5-sgNT'] = Synth40Scomp['epsilon_sgUBA5'] - Synth40Scomp['epsilon_sgNT']
Synth40Scomp['category'] = ''
for i in categories.category.unique():
    Synth40Scomp.loc[Synth40Scomp.gene.isin(categories.gene[categories.category == i]),'category'] = i
Synth40Scomp.sort_values('diff_epsilon_sgUBA5-sgNT').to_csv(paths.outp('fixed_NotI/compare_csv/synthesis40S_compare_fixedNotI.csv'))
Synth40Scomp.sort_values('diff_epsilon_sgUBA5-sgNT')
