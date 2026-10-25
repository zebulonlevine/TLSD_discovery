#!/usr/bin/env python3
"""Fig S3D (genome-wide vs sublibrary) and Fig S8A (sgNT vs sgUBA5): Pearson
correlation heatmaps of gene-level log2 fold changes.

Extracted from Fig2.ipynb cell(s) 6, 7. Cell code is verbatim except for paths:
inputs come from the re-scored tree via common.paths, categories from
screen_figures/categories/, and outputs go to fixed_NotI/plots/*_fixedNotI.pdf
under the output root. These are the same substitutions the original
regenerate_pRanked_figures.py harness applied when it exec'd these cells.

Run with envs/crispri-screen.yml:  python figS3D_figS8A_lfc_correlation.py
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


# ---- Fig2.ipynb cell 6 --------------------------------------------------
#compare genomewide to 
fileBaseSublibrary =  paths.inp('fixed_NotI/results')
fileBaseGenomewide = paths.inp('GenomeWide_results/results')
results_files = [
    os.path.join(fileBaseGenomewide,'L29_TopvBottom','L29_TopvBottom_all_genes.csv'),
    os.path.join(fileBaseGenomewide,'S3_TopvBottom','S3_TopvBottom_all_genes.csv'),
    os.path.join(fileBaseGenomewide,'L29_UnsortedvEarly','L29_UnsortedvEarly_all_genes.csv'),
    os.path.join(fileBaseGenomewide,'S3_UnsortedvEarly','S3_UnsortedvEarly_all_genes.csv'),
    
    os.path.join(fileBaseSublibrary,'L29_ratioON_top_v_bottom_wt','L29_ratioON_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_ratioON_top_v_bottom_wt','S3_ratioON_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'L29_unsorted_v_early_wt','L29_unsorted_v_early_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_unsorted_v_early_wt','S3_unsorted_v_early_wt_all_genes.csv'),
    
    os.path.join(fileBaseSublibrary,'L29_old_top_v_bottom_wt','L29_old_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_old_top_v_bottom_wt','S3_old_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'L29_new_top_v_bottom_wt','L29_new_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_new_top_v_bottom_wt','S3_new_top_v_bottom_wt_all_genes.csv'),
]
output_file = paths.outp('fixed_NotI/plots/heatmap_correlation_wt_fixedNotI.pdf')

screen_viz.plot_foldchange_heatmap(
    results_files,               # list[str]: _all_genes.csv files
    output_file,                 # str: output path
    comparison_names = None,     # list[str]: names for each file
    method           = 'pearson',# str: 'pearson' or 'spearman'
    cluster          = True,     # bool
    sample_order     = None,     # list[str]: manual order (disables clustering)
    title            = None,     # str
    figsize          = (10, 8),  # tuple[float, float]
    epsilon_col      = 'epsilon',# str
    merge_on         = 'gene'      # str: 'gene', 'gene_id', or auto-detect if None
)


# ---- Fig2.ipynb cell 7 --------------------------------------------------
#compare sgUBA5, wt
fileBaseSublibrary =  paths.inp('fixed_NotI/results')
fileBaseGenomewide = paths.inp('GenomeWide_results/results')
results_files = [
    
    
    os.path.join(fileBaseSublibrary,'L29_ratioON_top_v_bottom_wt','L29_ratioON_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_ratioON_top_v_bottom_wt','S3_ratioON_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'L29_unsorted_v_early_wt','L29_unsorted_v_early_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_unsorted_v_early_wt','S3_unsorted_v_early_wt_all_genes.csv'),
    
    os.path.join(fileBaseSublibrary,'L29_old_top_v_bottom_wt','L29_old_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_old_top_v_bottom_wt','S3_old_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'L29_new_top_v_bottom_wt','L29_new_top_v_bottom_wt_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_new_top_v_bottom_wt','S3_new_top_v_bottom_wt_all_genes.csv'),
    
    os.path.join(fileBaseSublibrary,'L29_ratioON_top_v_bottom_sgUBA5','L29_ratioON_top_v_bottom_sgUBA5_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_ratioON_top_v_bottom_sgUBA5','S3_ratioON_top_v_bottom_sgUBA5_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'L29_unsorted_v_early_sgUBA5','L29_unsorted_v_early_sgUBA5_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_unsorted_v_early_sgUBA5','S3_unsorted_v_early_sgUBA5_all_genes.csv'),
    
    os.path.join(fileBaseSublibrary,'L29_old_top_v_bottom_sgUBA5','L29_old_top_v_bottom_sgUBA5_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_old_top_v_bottom_sgUBA5','S3_old_top_v_bottom_sgUBA5_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'L29_new_top_v_bottom_sgUBA5','L29_new_top_v_bottom_sgUBA5_all_genes.csv'),
    os.path.join(fileBaseSublibrary,'S3_new_top_v_bottom_sgUBA5','S3_new_top_v_bottom_sgUBA5_all_genes.csv'),
]
output_file = paths.outp('fixed_NotI/plots/heatmap_correlation_sgUBA5_fixedNotI.pdf')

screen_viz.plot_foldchange_heatmap(
    results_files,               # list[str]: _all_genes.csv files
    output_file,                 # str: output path
    comparison_names = None,     # list[str]: names for each file
    method           = 'pearson',# str: 'pearson' or 'spearman'
    cluster          = True,     # bool
    sample_order     = None,     # list[str]: manual order (disables clustering)
    title            = None,     # str
    figsize          = (10, 8),  # tuple[float, float]
    epsilon_col      = 'epsilon',# str
    merge_on         = 'gene'      # str: 'gene', 'gene_id', or auto-detect if None
)
