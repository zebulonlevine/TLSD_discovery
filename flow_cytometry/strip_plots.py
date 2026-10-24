"""Shared helpers for the RiboHalo knockdown strip plots (Fig 3H, 4B, 4D, 7A; Fig S4I, S8B–D).

Each well is one point: the FlowJo median fluorescence of HaloTag-positive, GFP+/BFP+ singlets
(old ribosomes = first label, JF549, channel YL1-A; new ribosomes = chase label, JF646, RL1-A),
divided by the mean of the control (sgNT) wells of the same cell line and condition. The red line is
each group's mean.

Statistics: Welch's t-test between pairs of groups; Benjamini-Hochberg across the comparisons drawn in
one panel.
"""
import itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

GROUP_COL = 'SampleID'
JITTER_SEED = 0     # seaborn jitters strip points with the global NumPy RNG; seeded so reruns match

# group mean as a red line; box, whiskers, caps, median and fliers hidden
MEAN_LINE = dict(showmeans=True, meanline=True, meanprops={'color': 'red', 'ls': '-', 'lw': 2},
                 medianprops={'visible': False}, whiskerprops={'visible': False}, zorder=10,
                 showfliers=False, showbox=False, showcaps=False)


def normalize_to_control(df, value_col, out_col, key, is_control):
    """df[out_col] = value / mean(value of the control wells with the same key).

    key: Series giving each well's normalization group (e.g. cell line, or line + condition);
    is_control: boolean Series marking the control (sgNT) wells.
    """
    df[out_col] = np.nan
    for k in key[is_control].unique():
        rows = key == k
        df.loc[rows, out_col] = df.loc[rows, value_col] / df.loc[rows & is_control, value_col].mean()
    return df


def strip_with_mean(data, y, path, title=None, xlabel=None, ylabel=None):
    """One point per well, red line at the group mean; groups in the order of `data`
    (or of its categorical SampleID)."""
    sns.boxplot(x=GROUP_COL, y=y, data=data, **MEAN_LINE)
    np.random.seed(JITTER_SEED)
    sns.stripplot(data=data, x=GROUP_COL, y=y, color='black')
    if title is not None:
        plt.title(title)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    plt.xticks(rotation=90)
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close('all')


def group_summary(df, value_col, panel=''):
    """n, mean, median, SD and SEM per group."""
    groups = sorted(df[GROUP_COL].unique())
    desc = (df.groupby(GROUP_COL, observed=True)[value_col]
              .agg(n='count', mean='mean', median='median', std='std', sem='sem')
              .reindex(groups).reset_index())
    desc.insert(0, 'panel', panel)
    return desc


def pairwise_welch_bh(df, value_col, pairs=None, panel=''):
    """Welch's t-test for each pair of groups (all pairs by default), BH across the pairs."""
    if pairs is None:
        pairs = itertools.combinations(sorted(df[GROUP_COL].unique()), 2)
    rows = []
    for a, b in pairs:
        x = df.loc[df[GROUP_COL] == a, value_col].dropna().values
        y = df.loc[df[GROUP_COL] == b, value_col].dropna().values
        t, p = stats.ttest_ind(x, y, equal_var=False)
        v1, v2 = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
        rows.append(dict(panel=panel, metric=value_col, group1=a, group2=b,
                         n1=len(x), n2=len(y), mean1=x.mean(), mean2=y.mean(), t=t,
                         df_welch=(v1 + v2) ** 2 / (v1 ** 2 / (len(x) - 1) + v2 ** 2 / (len(y) - 1)),
                         p=p))
    out = pd.DataFrame(rows)
    out['p_BH'] = multipletests(out.p, method='fdr_bh')[1]
    return out
