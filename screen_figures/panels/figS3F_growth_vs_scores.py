#!/usr/bin/env python3
"""Fig S3F: Pearson correlation of growth, synthesis and degradation gene scores
(sgNT; cell 50 also writes the sgNT + sgUBA5 variant).

Extracted from Fig2.ipynb cell(s) 49, 50. Cell code is verbatim except for paths:
inputs come from the re-scored tree via common.paths, categories from
screen_figures/categories/, and outputs go to fixed_NotI/plots/*_fixedNotI.pdf
under the output root. These are the same substitutions the original
regenerate_pRanked_figures.py harness applied when it exec'd these cells.

Run with envs/crispri-screen.yml:  python figS3F_growth_vs_scores.py
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


# ---- Fig2.ipynb cell 49 --------------------------------------------------
# =========================================================================
# Correlation heatmap (WT only): growth effects vs. synthesis/degradation
# scores, across both sublibraries (60S = L29, 40S = S3).
# =========================================================================
# Sublibraries: L29 = 60S, S3 = 40S. "Growth" = unsorted_v_early (the growth/fitness
# comparison used as "Growth" elsewhere in this notebook). Synthesis/Degradation are the
# computed scores under Sublibrary_results/scores/. Correlation = pearson on gene epsilon.
sublibResults = paths.inp('fixed_NotI/results')
sublibScores  = paths.inp('fixed_NotI/scores')
sublibPlots   = paths.outp('fixed_NotI/plots')

def _growth_file(sub, cond):
    d = f'{sub}_unsorted_v_early_{cond}'
    return os.path.join(sublibResults, d, d + '_all_genes.csv')

def _score_file(sub, score, cond):
    sfx = '' if cond == 'wt' else '_sgUBA5'
    return os.path.join(sublibScores, f'{sub}_{score}{sfx}', 'gene_scores_all_genes.csv')

_subs = [('60S', 'L29'), ('40S', 'S3')]   # (display label, MAGeCK dir prefix)

def _inputs(cond, suffix=''):
    """Return (files, names) for growth + synthesis + degradation of both sublibraries."""
    files, names = [], []
    for measure, getter in [('Growth',      lambda sub: _growth_file(sub, cond)),
                            ('Synthesis',    lambda sub: _score_file(sub, 'Synthesis', cond)),
                            ('Degradation',  lambda sub: _score_file(sub, 'Degradation', cond))]:
        for lbl, sub in _subs:
            files.append(getter(sub))
            names.append(f'{lbl} {measure}{suffix}')
    return files, names

files_wt, names_wt = _inputs('wt')

screen_viz.plot_foldchange_heatmap(
    files_wt,
    paths.outp('fixed_NotI/plots/heatmap_growth_vs_synthdeg_wt_fixedNotI.pdf'),
    comparison_names = names_wt,
    method      = 'pearson',
    cluster     = True,
    epsilon_col = 'epsilon',
    merge_on    = 'gene',
    figsize     = (8, 6),
    title       = 'WT: growth vs. synthesis / degradation scores',
)


# ---- Fig2.ipynb cell 50 --------------------------------------------------
# =========================================================================
# Same correlation heatmap, now including sgUBA5 alongside WT:
# growth vs. synthesis/degradation scores for both sublibraries x {WT, sgUBA5}.
# =========================================================================
# Sublibraries: L29 = 60S, S3 = 40S. "Growth" = unsorted_v_early (the growth/fitness
# comparison used as "Growth" elsewhere in this notebook). Synthesis/Degradation are the
# computed scores under Sublibrary_results/scores/. Correlation = pearson on gene epsilon.
sublibResults = paths.inp('fixed_NotI/results')
sublibScores  = paths.inp('fixed_NotI/scores')
sublibPlots   = paths.outp('fixed_NotI/plots')

def _growth_file(sub, cond):
    d = f'{sub}_unsorted_v_early_{cond}'
    return os.path.join(sublibResults, d, d + '_all_genes.csv')

def _score_file(sub, score, cond):
    sfx = '' if cond == 'wt' else '_sgUBA5'
    return os.path.join(sublibScores, f'{sub}_{score}{sfx}', 'gene_scores_all_genes.csv')

_subs = [('60S', 'L29'), ('40S', 'S3')]   # (display label, MAGeCK dir prefix)

def _inputs(cond, suffix=''):
    """Return (files, names) for growth + synthesis + degradation of both sublibraries."""
    files, names = [], []
    for measure, getter in [('Growth',      lambda sub: _growth_file(sub, cond)),
                            ('Synthesis',    lambda sub: _score_file(sub, 'Synthesis', cond)),
                            ('Degradation',  lambda sub: _score_file(sub, 'Degradation', cond))]:
        for lbl, sub in _subs:
            files.append(getter(sub))
            names.append(f'{lbl} {measure}{suffix}')
    return files, names

files_wt,   names_wt   = _inputs('wt',     ' WT')
files_uba5, names_uba5 = _inputs('sgUBA5', ' sgUBA5')

screen_viz.plot_foldchange_heatmap(
    files_wt + files_uba5,
    paths.outp('fixed_NotI/plots/heatmap_growth_vs_synthdeg_wt_sgUBA5_fixedNotI.pdf'),
    comparison_names = names_wt + names_uba5,
    method      = 'pearson',
    cluster     = True,
    epsilon_col = 'epsilon',
    merge_on    = 'gene',
    figsize     = (10, 8),
    title       = 'WT + sgUBA5: growth vs. synthesis / degradation scores',
)
