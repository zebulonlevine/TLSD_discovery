#!/usr/bin/env python3
"""Fig S2C: sgRNA-level score construction scatters (flat grey and coloured by
line score), 60S and 40S, synthesis and degradation.

Extracted from Fig2.ipynb cell(s) 51, 52. Cell code is verbatim except for paths:
inputs come from the re-scored tree via common.paths, categories from
screen_figures/categories/, and outputs go to fixed_NotI/plots/*_fixedNotI.pdf
under the output root. These are the same substitutions the original
regenerate_pRanked_figures.py harness applied when it exec'd these cells.

Run with envs/crispri-screen.yml:  python figS2C_score_construction.py
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


# ---- Fig2.ipynb cell 51 --------------------------------------------------
# =========================================================================
# Square, rasterized-PDF versions of the sgRNA line-score deviation "flat"
# plots (plain grey), matching the core-plot aspect ratio. New filenames;
# the original *_score_dev_flat.pdf files are left untouched.
# =========================================================================
_R = paths.inp('fixed_NotI/results')
def _sg(name):
    return os.path.join(_R, name, name + '.sgrna_summary.txt')

# (out_prefix, sgrna1_dir, sgrna2_dir, slope, comp1_name, comp2_name)
_flat_specs = [
    ('60S_synth_score_dev_flat_square',
     'L29_new_top_v_bottom_wt', 'L29_ratioON_top_v_bottom_wt',  1,
     'New 60S log2 fold change', '60S New/Old Ratio log2 fold change'),
    ('60S_degrad_score_dev_flat_square',
     'L29_old_top_v_bottom_wt', 'L29_ratioON_top_v_bottom_wt', -1,
     'Old 60S log2 fold change', '60S New/Old Ratio log2 fold change'),
    ('40S_synth_score_dev_flat_square',
     'S3_new_top_v_bottom_wt', 'S3_ratioON_top_v_bottom_wt',  1,
     'New 40S log2 fold change', '40S New/Old Ratio log2 fold change'),
    ('40S_degrad_score_dev_flat_square',
     'S3_old_top_v_bottom_wt', 'S3_ratioON_top_v_bottom_wt', -1,
     'Old 40S log2 fold change', '40S New/Old Ratio log2 fold change'),
]

for out, sg1, sg2, slope, c1, c2 in _flat_specs:
    screen_viz.plot_sgrna_score_scatter(
        _sg(sg1), _sg(sg2), slope, 0, paths.outp('fixed_NotI/plots', out + '_fixedNotI.pdf'),
        comp1_name=c1, comp2_name=c2,
        show_ref_line=True, figsize=(2.8, 2.8), show_legend=False,
        plain=True, show_axis_lines=False, title='', dpi=600,
    )


# ---- Fig2.ipynb cell 52 --------------------------------------------------
# =========================================================================
# Non-flat (colored by weighted line score) versions of the sgRNA score
# development plots, at the SAME square sizing as the *_flat_square plots.
# Uses the perp-weighted score (score_mode='perp_weighted'), not fully-parallel.
# Colorbar has its own axes (plot region identical to the flat squares)
# and a fixed shared scale (vmin/vmax). New filenames; nothing overwritten.
# =========================================================================
_R = paths.inp('fixed_NotI/results')
def _sg(name):
    return os.path.join(_R, name, name + '.sgrna_summary.txt')

# (out_prefix, sgrna1_dir, sgrna2_dir, slope, comp1_name, comp2_name)
_dev_specs = [
    ('60S_synth_score_dev_square',
     'L29_new_top_v_bottom_wt', 'L29_ratioON_top_v_bottom_wt',  1,
     'New 60S log2 fold change', '60S New/Old Ratio log2 fold change'),
    ('60S_degrad_score_dev_square',
     'L29_old_top_v_bottom_wt', 'L29_ratioON_top_v_bottom_wt', -1,
     'Old 60S log2 fold change', '60S New/Old Ratio log2 fold change'),
    ('40S_synth_score_dev_square',
     'S3_new_top_v_bottom_wt', 'S3_ratioON_top_v_bottom_wt',  1,
     'New 40S log2 fold change', '40S New/Old Ratio log2 fold change'),
    ('40S_degrad_score_dev_square',
     'S3_old_top_v_bottom_wt', 'S3_ratioON_top_v_bottom_wt', -1,
     'Old 40S log2 fold change', '40S New/Old Ratio log2 fold change'),
]

for out, sg1, sg2, slope, c1, c2 in _dev_specs:
    screen_viz.plot_sgrna_score_scatter(
        _sg(sg1), _sg(sg2), slope, 0, paths.outp('fixed_NotI/plots', out + '_fixedNotI.pdf'),
        comp1_name=c1, comp2_name=c2,
        show_ref_line=True, figsize=(2.8, 2.8), show_legend=False,
        plain=False, score_mode='perp_weighted', show_axis_lines=False,
        title='', dpi=600,
        cbar_separate=True,     # colorbar in its own region; plot area unchanged
        vmin=-3.1, vmax=3.1,    # fixed line-score scale shared across all four
    )
