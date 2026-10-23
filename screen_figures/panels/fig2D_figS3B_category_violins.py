#!/usr/bin/env python3
"""Fig 2D (degradation) and Fig S3B (synthesis): category violins of gene
scores, 60S and 40S; outlined points pass 5% FDR.

Extracted from Fig2.ipynb cell(s) 41, 42, 43, 44. Cell code is verbatim except for paths:
inputs come from the re-scored tree via common.paths, categories from
screen_figures/categories/, and outputs go to fixed_NotI/plots/*_fixedNotI.pdf
under the output root. These are the same substitutions the original
regenerate_pRanked_figures.py harness applied when it exec'd these cells.

Run with envs/crispri-screen.yml:  python fig2D_figS3B_category_violins.py
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


# ---- Fig2.ipynb cell 41 --------------------------------------------------
catsMany = screen_viz.load_categories(paths.repo('screen_figures', 'categories', 'manyCats.csv'))
screen_viz.plot_category_enrichment(
    results_file = paths.inp('fixed_NotI/scores/L29_Degradation/gene_scores_all_genes.csv'), 
    output_file = paths.outp('fixed_NotI/plots/60Sdeg_violinPlots_fixedNotI.pdf'),
    categories = catsMany,
    hits_file =  paths.inp('fixed_NotI/scores/L29_Degradation/gene_scores_hits_fdr0.05.csv'),
    title = '60S Degradation',
    figsize = (13,3.5),
    orientation = 'vertical',
    show_legend = False,
    colors = {x : 'darkgrey' for x in catsMany.category.unique()},
    category_order = ['Both Biogenesis','40S Biogenesis/Proteins','60S Biogenesis/Proteins',
                      'RQC','Autophagy','Lysosome','Proteasome','E3','UFMylation'],
    show_ntc = True
)


# ---- Fig2.ipynb cell 42 --------------------------------------------------
catsMany = screen_viz.load_categories(paths.repo('screen_figures', 'categories', 'manyCats.csv'))
screen_viz.plot_category_enrichment(
    results_file = paths.inp('fixed_NotI/scores/S3_Degradation/gene_scores_all_genes.csv'), 
    output_file = paths.outp('fixed_NotI/plots/40Sdeg_violinPlots_fixedNotI.pdf'),
    categories = catsMany,
    hits_file =  paths.inp('fixed_NotI/scores/S3_Degradation/gene_scores_hits_fdr0.05.csv'),
    title = '40S Degradation',
    figsize = (13,3.5),
    orientation = 'vertical',
    show_legend = False,
    colors = {x : 'darkgrey' for x in catsMany.category.unique()},
    category_order = ['Both Biogenesis','40S Biogenesis/Proteins','60S Biogenesis/Proteins',
                      'RQC','Autophagy','Lysosome','Proteasome','E3','UFMylation'],
    show_ntc = True
)


# ---- Fig2.ipynb cell 43 --------------------------------------------------
catsMany = screen_viz.load_categories(paths.repo('screen_figures', 'categories', 'manyCats.csv'))
screen_viz.plot_category_enrichment(
    results_file = paths.inp('fixed_NotI/scores/L29_Synthesis/gene_scores_all_genes.csv'), 
    output_file = paths.outp('fixed_NotI/plots/60Ssynth_violinPlots_fixedNotI.pdf'),
    categories = catsMany,
    hits_file =  paths.inp('fixed_NotI/scores/L29_Synthesis/gene_scores_hits_fdr0.05.csv'),
    title = '60S Synth',
    figsize = (13,3.5),
    orientation = 'vertical',
    show_legend = False,
    colors = {x : 'darkgrey' for x in catsMany.category.unique()},
    category_order = ['Both Biogenesis','40S Biogenesis/Proteins','60S Biogenesis/Proteins',
                      'RQC','Autophagy','Lysosome','Proteasome','E3','UFMylation'],
    show_ntc = True
)


# ---- Fig2.ipynb cell 44 --------------------------------------------------
catsMany = screen_viz.load_categories(paths.repo('screen_figures', 'categories', 'manyCats.csv'))
screen_viz.plot_category_enrichment(
    results_file = paths.inp('fixed_NotI/scores/S3_Synthesis/gene_scores_all_genes.csv'), 
    output_file = paths.outp('fixed_NotI/plots/40Ssynth_violinPlots_fixedNotI.pdf'),
    categories = catsMany,
    hits_file =  paths.inp('fixed_NotI/scores/S3_Synthesis/gene_scores_hits_fdr0.05.csv'),
    title = '40S Synth',
    figsize = (13,3.5),
    orientation = 'vertical',
    show_legend = False,
    colors = {x : 'darkgrey' for x in catsMany.category.unique()},
    category_order = ['Both Biogenesis','40S Biogenesis/Proteins','60S Biogenesis/Proteins',
                      'RQC','Autophagy','Lysosome','Proteasome','E3','UFMylation'],
    show_ntc = True
)
