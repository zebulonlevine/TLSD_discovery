#!/usr/bin/env python3
"""
Per-sgRNA noise model for the synthesis / degradation line scores, using
pooled within-group replicate residuals instead of MAGeCK's control-only
mean-variance model. `build()` is used by final_pvalues.sgrna_table().

Why this is possible at all
---------------------------
Both scores are *linear* in two MAGeCK log-fold-changes (recovered from the
stored distances: slope +1 for synthesis, -1 for degradation, intercept 0):

    synthesis  parallel = (LFC_ratio + LFC_new) / sqrt(2)
    degradation parallel = (LFC_ratio - LFC_old) / sqrt(2)

so the variance propagates in closed form -- no delta method needed. The
`perp_weighted` score is NOT linear, so the test is built on `parallel_dist`
(the signed projection along the reference line, i.e. the actual effect axis)
and `perp_dist` is kept only as a QC/consistency measure.

Variance model
--------------
Every group here has 2 replicates (4 groups per score: x_Top, x_Bottom,
y_Top, y_Bottom). Subtract each guide's per-group mean and pool the residuals
across all 4 groups:

    s2_i = sum_g (l_i,g,1 - l_i,g,2)^2 / 2  / 4        (4 df, signal-free)

With 2 replicates per group, Var(mean_g) = s2/2, so
Var(LFC) = s2/2 + s2/2 = s2 and

    Var(parallel) = [Var(LFC_x) + Var(LFC_y)] / 2 = s2

i.e. z = parallel / sqrt(s2_adj). 4 df is still noisy per guide, so s2 is
shrunk toward a fitted log-log mean-variance trend (MAGeCK's idea, but fit to
signal-free within-group residuals rather than to control samples, and never
to control+treatment pooled -- which is what MAGeCK falls back to when a
comparison has one replicate, and which absorbs real signal into the null).

Calibration is checked directly: NTC sgRNA p-values must be ~Uniform(0,1).
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from common import paths  # noqa: E402

COUNTS = paths.inp('Sublibrary_results/results/counts_for_mageck.txt')

NTC_LABELS = frozenset({'negative_control', 'negative control', 'non-targeting',
                        'nontargeting', 'ntc', 'neg_control', 'negative', 'safe_harbor'})


def is_ntc(names):
    head = names.astype(str).str.strip().str.split('_').str[0].str.lower()
    return head.isin(NTC_LABELS)


# score -> (y group stem, slope). x is always ratioON.
SCORES = {'Synthesis': ('new', +1.0), 'Degradation': ('old', -1.0)}
SCREENS = {'L29': '60S', 'S3': '40S'}
REPS = ('rep1', 'rep2')
PSEUDO = 32.0          # log2 pseudocount; damps low-count guides


def median_size_factors(mat):
    """DESeq/MAGeCK-style median-of-ratios size factors."""
    logm = np.log(mat.replace(0, np.nan))
    ref = logm.mean(axis=1)                       # geometric mean per guide
    sf = np.exp((logm.sub(ref, axis=0)).median(axis=0, skipna=True))
    return sf / np.exp(np.log(sf).mean())


COMPARISONS = paths.repo('crispr_screens/pipeline/comparisons_sublibrary.csv')


def _arms(comp_name):
    """(treatment samples, control samples) for a comparison, from the table.

    Do NOT infer orientation from the comparison name: '*_ratioON_top_v_bottom_wt'
    is actually treatment=Bottom, control=Top -- reversed relative to the
    '*_new/old_top_v_bottom_wt' comparisons. Assuming Top=treatment flips the
    sign of the x axis and destroys the score (corr with the published score
    drops to -0.05).
    """
    tab = pd.read_csv(COMPARISONS).set_index('comparison_name')
    row = tab.loc[comp_name]
    split = lambda v: [x.strip() for x in str(v).split(',') if x.strip()]
    return split(row['treatment']), split(row['control'])


def build(screen, score, cond='wt', min_count=20):
    """cond: 'wt' or 'sgUBA5'. Both have 2 replicates in all 4 groups and the
    same reversed ratioON orientation."""
    ystem, slope = SCORES[score]
    x_treat, x_ctrl = _arms(f'{screen}_ratioON_top_v_bottom_{cond}')
    y_treat, y_ctrl = _arms(f'{screen}_{ystem}_top_v_bottom_{cond}')
    groups = {'x_top': x_treat, 'x_bot': x_ctrl, 'y_top': y_treat, 'y_bot': y_ctrl}
    cols = [c for v in groups.values() for c in v]

    raw = pd.read_table(COUNTS, usecols=['sgRNA', 'Gene'] + cols)
    sf = median_size_factors(raw[cols])
    norm = raw[cols].div(sf, axis=1)

    # depth filter, mirroring the existing min_count=20 on treat_mean
    keep = norm.mean(axis=1) > min_count
    raw, norm = raw[keep].reset_index(drop=True), norm[keep].reset_index(drop=True)

    l = np.log2(norm + PSEUDO)

    # per-group means and signal-free pooled residual variance
    gmean, ss = {}, np.zeros(len(l))
    for g, (c1, c2) in groups.items():
        gmean[g] = (l[c1] + l[c2]) / 2.0
        ss += (l[c1] - l[c2]).values ** 2 / 2.0        # 1 df per group
    n_df = len(groups)                                  # 4 groups x 1 df
    s2_raw = ss / n_df

    lfc_x = (gmean['x_top'] - gmean['x_bot']).values
    lfc_y = (gmean['y_top'] - gmean['y_bot']).values
    norm_f = np.sqrt(1 + slope ** 2)
    parallel = (lfc_x + slope * lfc_y) / norm_f
    perp = np.abs(lfc_y - slope * lfc_x) / norm_f

    # shrink s2 toward a fitted log-log trend against expression
    expr = np.log2(norm.mean(axis=1).values + PSEUDO)
    ok = np.isfinite(s2_raw) & (s2_raw > 0) & np.isfinite(expr)
    k, b = np.polyfit(expr[ok], np.log2(s2_raw[ok]), 1)
    s2_fit = 2 ** (k * expr + b)
    # conservative: never let a guide claim less noise than the trend predicts
    s2_adj = np.maximum(s2_raw, s2_fit)

    var_parallel = s2_adj                               # see module docstring
    z = parallel / np.sqrt(var_parallel)
    p = 2 * stats.norm.sf(np.abs(z))

    out = pd.DataFrame({
        'sgrna': raw['sgRNA'], 'gene': raw['Gene'],
        'lfc_x_ratio': lfc_x, f'lfc_y_{ystem}': lfc_y,
        'parallel': parallel, 'perp': perp,
        'expr_log2': expr, 's2_raw': s2_raw, 's2_fit': s2_fit, 's2_adj': s2_adj,
        'z': z, 'pvalue': p,
    })
    out['is_ntc'] = is_ntc(out['gene']).values
    out.attrs['fit'] = (k, b)
    out.attrs['slope'] = slope
    out.attrs['cond'] = cond
    return out
