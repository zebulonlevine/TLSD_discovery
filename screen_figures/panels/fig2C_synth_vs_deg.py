#!/usr/bin/env python3
"""Fig 2C: gene-level synthesis (y) vs degradation (x) scores, 60S and 40S,
biogenesis split by subunit.

Extracted from Fig2.ipynb cell(s) 47, 48. Cell code is verbatim except for paths:
inputs come from the re-scored tree via common.paths, categories from
screen_figures/categories/, and outputs go to fixed_NotI/plots/*_fixedNotI.pdf
under the output root. These are the same substitutions the original
regenerate_pRanked_figures.py harness applied when it exec'd these cells.

Run with envs/crispri-screen.yml:  python fig2C_synth_vs_deg.py
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


# ---- Fig2.ipynb cell 47 --------------------------------------------------
# =========================================================================
# NEW VERSION of compare_60Sscores.pdf -- ribosome biogenesis split into
# multiple shades of green (Both subunits / 40S / 60S Biogenesis/Proteins),
# with RNA Pol I placed in "Both subunits". Points are rasterized (dpi=600 raster
# layer). Writes NEW files; the original compare_60Sscores.pdf is untouched.
# =========================================================================
pathDeg60S   = paths.inp('fixed_NotI/scores/L29_Degradation')
pathSynth60S = paths.inp('fixed_NotI/scores/L29_Synthesis')

# Build a category table where "Biogenesis" is split by subunit into three
# shades-of-green groups, using GenomeWideCats_synth_bySubunit.csv. RNA Pol I
# subunits are forced into "Both subunits". Lysosome/Autophagy + UFMylation are
# carried over from GenomeWideCats_threeOptions.csv unchanged.
pathCats      = paths.repo('screen_figures', 'categories', 'GenomeWideCats_threeOptions.csv')
pathBySubunit = paths.repo('screen_figures', 'categories', 'GenomeWideCats_synth_bySubunit.csv')

rna_pol_I = ['POLR1A','POLR1B','POLR1C','POLR1D','POLR1E','POLR1F','POLR1G','POLR1H']

bio_bysub = screen_viz.load_categories(pathBySubunit)[['gene','category']].copy()
bio_bysub['category'] = bio_bysub['category'].replace({'Both Biogenesis':'Both subunits'})
bio_bysub.loc[bio_bysub['gene'].isin(rna_pol_I), 'category'] = 'Both subunits'

_base    = screen_viz.load_categories(pathCats)[['gene','category']]
_non_bio = _base.loc[_base['category'] != 'Biogenesis']
categories_bysub = pd.concat([bio_bysub, _non_bio], ignore_index=True)

# Manual curation: assign these genes to Lysosome/Autophagy
# (not present in the source category files).
_lyso_manual = ['RNASEK', 'ATG14']
categories_bysub = categories_bysub[~categories_bysub['gene'].isin(_lyso_manual)]
categories_bysub = pd.concat(
    [categories_bysub, pd.DataFrame({'gene': _lyso_manual,
                                     'category': ['Lysosome/Autophagy'] * len(_lyso_manual)})],
    ignore_index=True)

# Three shades of green for ribosome biogenesis; other categories keep their colors.
green_colors = {
    'Both subunits'           : '#00441b',  # dark green (includes RNA Pol I)
    '60S Biogenesis/Proteins' : '#238b45',  # medium green
    '40S Biogenesis/Proteins' : '#74c476',  # light green
    'Lysosome/Autophagy'      : '#c51b7d',
    'UFMylation'              : '#E69F00',
}
bio_order = ['Both subunits','60S Biogenesis/Proteins','40S Biogenesis/Proteins',
             'Lysosome/Autophagy','UFMylation']

comp1 = screen_viz.load_comparison(os.path.join(pathDeg60S,'gene_scores_all_genes.csv'), os.path.join(pathDeg60S,'gene_scores_hits_fdr0.05.csv'), '60S Degradation Score')
comp2 = screen_viz.load_comparison(os.path.join(pathSynth60S,'gene_scores_all_genes.csv'), os.path.join(pathSynth60S,'gene_scores_hits_fdr0.05.csv'), '60S Synthesis Score')

# Label the 5 most extreme genes in EACH direction (top/bottom of x, top/bottom of y).
def _top_each_direction(cx, cy, n=5):
    dx = screen_viz.deduplicate_by_max_abs_epsilon(cx.results, 'gene_id', 'epsilon')
    dy = screen_viz.deduplicate_by_max_abs_epsilon(cy.results, 'gene_id', 'epsilon')
    m = dx[['gene_id','gene','epsilon']].merge(dy[['gene_id','epsilon']], on='gene_id', suffixes=('_x','_y'))
    genes = set()
    for col in ['epsilon_x','epsilon_y']:
        genes |= set(m.nlargest(n, col)['gene'])
        genes |= set(m.nsmallest(n, col)['gene'])
    return sorted(genes)

label_list = _top_each_direction(comp1, comp2, 5)

screen_viz.plot_comparison_scatter(
    comp1, comp2, paths.outp('fixed_NotI/plots/compare_60Sscores_bioBySubunit_fixedNotI.pdf'),
    categories=categories_bysub, title = None,
    label_genes=label_list, label_top_n_x=0, label_top_n_y=0, adjust_labels=True,
    show_diagonal = False, show_correlation = False, figsize = (4.2,2.8),
    legend_outside = True, size_by_hits = False, point_size = 30, fontsize = 8, show_axis_labels = False, dpi = 600,
    colors = green_colors, category_order = bio_order,
)

# --- No-legend version (same plot, legend suppressed) ---
screen_viz.plot_comparison_scatter(
    comp1, comp2, paths.outp('fixed_NotI/plots/compare_60Sscores_bioBySubunit_noLegend_fixedNotI.pdf'),
    categories=categories_bysub, title = None,
    label_genes=label_list, label_top_n_x=0, label_top_n_y=0, adjust_labels=True,
    show_diagonal = False, show_correlation = False, figsize = (4.2,2.8),
    legend_outside = True, size_by_hits = False, point_size = 30, fontsize = 8, show_axis_labels = False, dpi = 600,
    colors = green_colors, category_order = bio_order, show_legend = False,
)


# ---- Fig2.ipynb cell 48 --------------------------------------------------
# =========================================================================
# NEW VERSION of compare_40Sscores.pdf -- ribosome biogenesis split into
# multiple shades of green (Both subunits / 40S / 60S Biogenesis/Proteins),
# with RNA Pol I placed in "Both subunits". Points are rasterized (dpi=600 raster
# layer). Writes NEW files; the original compare_40Sscores.pdf is untouched.
# =========================================================================
pathDeg40S   = paths.inp('fixed_NotI/scores/S3_Degradation')
pathSynth40S = paths.inp('fixed_NotI/scores/S3_Synthesis')

# Build a category table where "Biogenesis" is split by subunit into three
# shades-of-green groups, using GenomeWideCats_synth_bySubunit.csv. RNA Pol I
# subunits are forced into "Both subunits". Lysosome/Autophagy + UFMylation are
# carried over from GenomeWideCats_threeOptions.csv unchanged.
pathCats      = paths.repo('screen_figures', 'categories', 'GenomeWideCats_threeOptions.csv')
pathBySubunit = paths.repo('screen_figures', 'categories', 'GenomeWideCats_synth_bySubunit.csv')

rna_pol_I = ['POLR1A','POLR1B','POLR1C','POLR1D','POLR1E','POLR1F','POLR1G','POLR1H']

bio_bysub = screen_viz.load_categories(pathBySubunit)[['gene','category']].copy()
bio_bysub['category'] = bio_bysub['category'].replace({'Both Biogenesis':'Both subunits'})
bio_bysub.loc[bio_bysub['gene'].isin(rna_pol_I), 'category'] = 'Both subunits'

_base    = screen_viz.load_categories(pathCats)[['gene','category']]
_non_bio = _base.loc[_base['category'] != 'Biogenesis']
categories_bysub = pd.concat([bio_bysub, _non_bio], ignore_index=True)

# Manual curation: assign these genes to Lysosome/Autophagy
# (not present in the source category files).
_lyso_manual = ['RNASEK', 'ATG14']
categories_bysub = categories_bysub[~categories_bysub['gene'].isin(_lyso_manual)]
categories_bysub = pd.concat(
    [categories_bysub, pd.DataFrame({'gene': _lyso_manual,
                                     'category': ['Lysosome/Autophagy'] * len(_lyso_manual)})],
    ignore_index=True)

# Three shades of green for ribosome biogenesis; other categories keep their colors.
green_colors = {
    'Both subunits'           : '#00441b',  # dark green (includes RNA Pol I)
    '60S Biogenesis/Proteins' : '#238b45',  # medium green
    '40S Biogenesis/Proteins' : '#74c476',  # light green
    'Lysosome/Autophagy'      : '#c51b7d',
    'UFMylation'              : '#E69F00',
}
bio_order = ['Both subunits','60S Biogenesis/Proteins','40S Biogenesis/Proteins',
             'Lysosome/Autophagy','UFMylation']

comp1 = screen_viz.load_comparison(os.path.join(pathDeg40S,'gene_scores_all_genes.csv'), os.path.join(pathDeg40S,'gene_scores_hits_fdr0.05.csv'), '40S Degradation Score')
comp2 = screen_viz.load_comparison(os.path.join(pathSynth40S,'gene_scores_all_genes.csv'), os.path.join(pathSynth40S,'gene_scores_hits_fdr0.05.csv'), '40S Synthesis Score')

# Label the 5 most extreme genes in EACH direction (top/bottom of x, top/bottom of y).
def _top_each_direction(cx, cy, n=5):
    dx = screen_viz.deduplicate_by_max_abs_epsilon(cx.results, 'gene_id', 'epsilon')
    dy = screen_viz.deduplicate_by_max_abs_epsilon(cy.results, 'gene_id', 'epsilon')
    m = dx[['gene_id','gene','epsilon']].merge(dy[['gene_id','epsilon']], on='gene_id', suffixes=('_x','_y'))
    genes = set()
    for col in ['epsilon_x','epsilon_y']:
        genes |= set(m.nlargest(n, col)['gene'])
        genes |= set(m.nsmallest(n, col)['gene'])
    return sorted(genes)

label_list = _top_each_direction(comp1, comp2, 5)

screen_viz.plot_comparison_scatter(
    comp1, comp2, paths.outp('fixed_NotI/plots/compare_40Sscores_bioBySubunit_fixedNotI.pdf'),
    categories=categories_bysub, title = None,
    label_genes=label_list, label_top_n_x=0, label_top_n_y=0, adjust_labels=True,
    show_diagonal = False, show_correlation = False, figsize = (4.2,2.8),
    legend_outside = True, size_by_hits = False, point_size = 30, fontsize = 8, show_axis_labels = False, dpi = 600,
    colors = green_colors, category_order = bio_order,
)

# --- No-legend version (same plot, legend suppressed) ---
screen_viz.plot_comparison_scatter(
    comp1, comp2, paths.outp('fixed_NotI/plots/compare_40Sscores_bioBySubunit_noLegend_fixedNotI.pdf'),
    categories=categories_bysub, title = None,
    label_genes=label_list, label_top_n_x=0, label_top_n_y=0, adjust_labels=True,
    show_diagonal = False, show_correlation = False, figsize = (4.2,2.8),
    legend_outside = True, size_by_hits = False, point_size = 30, fontsize = 8, show_axis_labels = False, dpi = 600,
    colors = green_colors, category_order = bio_order, show_legend = False,
)
