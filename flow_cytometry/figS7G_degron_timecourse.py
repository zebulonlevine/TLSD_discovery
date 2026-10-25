"""Fig S7G: old, new and total ribosomes after UBA5 degradation (SD40-UBA5 degron, PT-179).

RPE1 HaloTag-RPS24 (40S, line 'S24') and RPL29-HaloTag (60S, line 'L29') cells carrying
SNAP-Sec61b-2A-SD40-UBA5 over sgUBA5 were pulse-labelled with JF549 (old ribosomes) and chased with
JF646 (new ribosomes), with 20 uM PT-179 or DMSO from the time of labelling; 4 wells per condition at
days 1, 4, 7 and 10. Single-dye wells (JF549 only, JF646 only) give each dye's total signal.
Input: FlowJo export of per-well medians (YL1-A = JF549, RL1-A = JF646) after singlet and
GFP/BFP gating.

  1. Background: the mean JF646-channel signal of JF549-only wells (and vice versa) is subtracted
     from the single-dye totals.
  2. One S24 Base JF646-only DMSO group (day 10) is excluded (failed chase labelling).
  3. Normalization per cell line and genotype, pooled over days: old = JF549 / mean DMSO JF549-only
     signal, new = JF646 / mean DMSO JF646-only signal (the denominators are background-subtracted,
     the per-well numerators are not); total = old + new.
  4. Mean per condition (old, new: pulse-chase wells, n = 4; total: single-dye wells, n = 8);
     error = replicate SEM and the SEM of the normalizing mean, combined in quadrature.
  5. DMSO vs PT-179 at each day: Welch's t-test, BH across all 96 comparisons of the experiment
     (2 lines x 4 genotypes x old/new/total x 4 days). The total test pools single-dye and
     pulse-chase wells (n = 12). Fig S7G shows the 'Base' genotype (degron only); the stars on the
     panel are the BH-adjusted values.

Outputs (output root, flow_cytometry/):
  figS7G_S24_Base.pdf, figS7G_L29_Base.pdf     40S (top row) and 60S (bottom row) of the panel
  figS7G_summary.csv                           plotted means and errors, all genotypes
  figS7G_welch_tests.csv                       all 96 tests, raw and BH-adjusted p

    python flow_cytometry/figS7G_degron_timecourse.py
"""
import os
import sys
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib as mpl
if 'ipykernel' not in sys.modules:
    mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402

DATA_FILE = 'flow_timecourse_2603.xlsx'
OUT_DIR = 'flow_cytometry'
PANEL_GENETICS = ['Base']              # genotype(s) drawn for the paper panel

GROUP_COLS = ['Cell Line', 'Genetics']             # normalization groups (pooled over days)
CONDITION_COLS = ['Day', 'Cell Line', 'Genetics', 'PT179']
METRICS = ('old', 'new', 'total')
VALUE_COL = {'old': 'old_frac', 'new': 'new_frac', 'total': 'total_frac'}
YLABEL = {'old': 'Old / Total (549)', 'new': 'New / Total (646)', 'total': 'Normalized total'}
LOG_SCALE = ('old', 'new')

GENETICS_COLORS = {'Base': '#333333', 'CJZ4': '#1f77b4', 'sgRIOK3': '#e69f00', 'sgR3HDM4': '#2ca02c'}
PT179_STYLE = {0: {'ls': '-', 'marker': 'o'}, 20: {'ls': '--', 'marker': 's'}}
STYLE = {
    'figure.dpi': 130, 'figure.facecolor': 'white',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.labelsize': 11, 'axes.titlesize': 12,
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 9, 'legend.frameon': False,
}


def out(name):
    return paths.outp(OUT_DIR, name, mkdir=True)


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------

def load():
    df = pd.read_excel(paths.inp('flow_cytometry', 'S7G', DATA_FILE))
    df['PT179'] = df['PT179'].astype(int)
    df['SampleType'] = df['SampleType'].astype(str)
    return df


def subtract_background(df):
    """Spillover/background of each dye's channel, measured in wells stained with the other dye only."""
    backgrounds = {'RL1A_Median': df[df['SampleType'] == '549'].RL1A_Median.mean(),
                   'YL1A_Median': df[df['SampleType'] == '646'].YL1A_Median.mean()}
    for col, bg in backgrounds.items():
        df[col + '_corr'] = df[col] - bg
    return df


def drop_outlier(df):
    """S24 Base, JF646-only, DMSO, day 10: no usable chase labelling."""
    return df[~((df['Cell Line'] == 'S24') & (df['Genetics'] == 'Base') & (df['SampleType'] == '646')
                & (df['PT179'] == 0) & (df['Day'] == 10))]


def normalization_table(df):
    """Mean and SEM of the DMSO single-dye totals per cell line and genotype, all days pooled."""
    s549 = df[(df['SampleType'] == '549') & (df['PT179'] == 0)]
    norm_old = (s549.groupby(GROUP_COLS)['YL1A_Median_corr']
                .agg(total_old_mean='mean', total_old_sem='sem', total_old_n='count').reset_index())
    s646 = df[(df['SampleType'] == '646') & (df['PT179'] == 0)]
    norm_new = (s646.groupby(GROUP_COLS)['RL1A_Median_corr']
                .agg(total_new_mean='mean', total_new_sem='sem', total_new_n='count').reset_index())
    return norm_old.merge(norm_new, on=GROUP_COLS, how='outer')


def normalize(df, norm):
    """Old, new and total fractions per well, with the relative error of the normalizing mean."""
    d = df.copy().merge(norm, on=GROUP_COLS, how='left')
    d['old_frac'] = d['YL1A_Median'] / d['total_old_mean']
    d['new_frac'] = d['RL1A_Median'] / d['total_new_mean']
    d['old_frac_denom_err'] = d['old_frac'] * (d['total_old_sem'] / d['total_old_mean'])
    d['new_frac_denom_err'] = d['new_frac'] * (d['total_new_sem'] / d['total_new_mean'])
    d['total_frac'] = d['old_frac'] + d['new_frac']
    d['total_frac_err'] = np.sqrt(d['old_frac_denom_err'] ** 2 + d['new_frac_denom_err'] ** 2)
    return d


def summarize_with_combined_err(data, group_cols, value_col, denom_err_col):
    """Mean +/- combined error (replicate SEM and denominator uncertainty in quadrature) per group."""
    records = []
    for keys, grp in data.groupby(group_cols):
        row = dict(zip(group_cols, keys)) if isinstance(keys, tuple) else {group_cols[0]: keys}
        vals = grp[value_col].values
        denom_errs = grp[denom_err_col].values
        row['n'] = len(vals)
        row['mean'] = np.mean(vals)
        row['replicate_sem'] = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
        row['denom_err'] = np.mean(denom_errs)
        row['combined_err'] = np.sqrt(row['replicate_sem'] ** 2 + row['denom_err'] ** 2)
        records.append(row)
    return pd.DataFrame(records)


def summaries(d):
    """Old and new from the pulse-chase wells; total from the single-dye wells."""
    test = d[d['SampleType'] == 'test']
    single = d[d['SampleType'].isin(['549', '646'])]
    return {
        'old': summarize_with_combined_err(test, CONDITION_COLS, 'old_frac', 'old_frac_denom_err').assign(metric='old'),
        'new': summarize_with_combined_err(test, CONDITION_COLS, 'new_frac', 'new_frac_denom_err').assign(metric='new'),
        'total': summarize_with_combined_err(single, CONDITION_COLS, 'total_frac', 'total_frac_err').assign(metric='total'),
    }


def test_wells(d, cl, gen, metric):
    """Wells tested for one panel: pulse-chase wells for old/new; all wells for total."""
    sub = d[(d['Cell Line'] == cl) & (d['Genetics'] == gen)]
    return sub[sub['SampleType'] == 'test'] if metric in ('old', 'new') else sub


def welch_tests(d):
    """DMSO vs PT-179 per day for every line, genotype and measure; BH across all of them."""
    records = []
    for cl in sorted(d['Cell Line'].unique()):
        for gen in sorted(d['Genetics'].unique()):
            for metric in METRICS:
                sub_r = test_wells(d, cl, gen, metric)
                vcol = VALUE_COL[metric]
                for day in sorted(sub_r['Day'].unique()):
                    day_data = sub_r[sub_r['Day'] == day]
                    for g1, g2 in combinations(sorted(day_data['PT179'].unique()), 2):
                        v1 = day_data[day_data['PT179'] == g1][vcol].dropna().values
                        v2 = day_data[day_data['PT179'] == g2][vcol].dropna().values
                        if len(v1) < 2 or len(v2) < 2:
                            continue
                        t_stat, p_val = stats.ttest_ind(v1, v2, equal_var=False)
                        records.append({'Cell Line': cl, 'Genetics': gen, 'Metric': metric, 'Day': day,
                                        'PT179 1': g1, 'n1': len(v1), 'Mean 1': np.mean(v1),
                                        'PT179 2': g2, 'n2': len(v2), 'Mean 2': np.mean(v2),
                                        't_stat': t_stat, 'p_raw': p_val})
    t = pd.DataFrame(records)
    t['p_bh'] = multipletests(t['p_raw'].values, alpha=0.05, method='fdr_bh')[1]
    t['n_comparisons'] = len(t)
    t['sig'] = np.where(t['p_bh'] < 0.001, '***', np.where(t['p_bh'] < 0.01, '**',
                        np.where(t['p_bh'] < 0.05, '*', 'ns')))
    return t


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------

def plot_line(summ, cl, gen, path):
    """Old, new and total over time, DMSO (solid) vs PT-179 (dashed), for one line and genotype."""
    with mpl.rc_context(STYLE):
        fig, axes = plt.subplots(1, len(METRICS), figsize=(5.5 * len(METRICS), 4.5))
        for ax, metric in zip(axes, METRICS):
            sub_s = summ[metric]
            sub_s = sub_s[(sub_s['Cell Line'] == cl) & (sub_s['Genetics'] == gen)]
            for pt, grp in sub_s.groupby('PT179'):
                grp = grp.sort_values('Day')
                style = PT179_STYLE.get(pt, {'ls': ':', 'marker': 'D'})
                pt_label = '' if pt == 0 else ' +PT179'
                ax.errorbar(grp['Day'], grp['mean'], yerr=grp['combined_err'],
                            color=GENETICS_COLORS.get(gen, 'gray'), capsize=3, lw=1.5, zorder=2,
                            label=f'{pt_label.strip() or "untreated"}', **style)
            if metric in LOG_SCALE:
                ax.set_yscale('log')
            else:
                ax.axhline(1.0, color='gray', ls=':', lw=0.8, alpha=0.5)
            ax.set_xlabel('Day')
            ax.set_ylabel(YLABEL[metric])
            ax.set_title(metric.capitalize())
            ax.legend(fontsize=7)
        fig.suptitle(f'{cl} — {gen}', fontsize=13, y=1.02)
        plt.tight_layout()
        fig.savefig(path, dpi=300, bbox_inches='tight')
        plt.close(fig)


def main():
    df = drop_outlier(subtract_background(load()))
    d = normalize(df, normalization_table(df))
    summ = summaries(d)
    pd.concat(summ.values(), ignore_index=True).to_csv(out('figS7G_summary.csv'), index=False)
    tests = welch_tests(d)
    tests.to_csv(out('figS7G_welch_tests.csv'), index=False)
    for cl in sorted(d['Cell Line'].unique()):
        for gen in PANEL_GENETICS:
            plot_line(summ, cl, gen, out(f'figS7G_{cl}_{gen}.pdf'))
    shown = tests[tests['Genetics'].isin(PANEL_GENETICS)]
    print(f'{len(tests)} DMSO vs PT-179 tests; panel ({", ".join(PANEL_GENETICS)}):')
    print(shown.pivot_table(index=['Cell Line', 'Metric'], columns='Day', values='sig', aggfunc='first').to_string())
    print('wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
