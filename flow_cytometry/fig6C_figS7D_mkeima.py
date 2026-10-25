"""RPL29-mKeima ribophagy assay (Fig 6C, Fig S7D).

RPE1 CRISPRi cells expressing RPL29-mKeima in three backgrounds:
  - polyclonal (L29_mKeima);
  - the same cells also carrying the degradable SD40-UBA5 construct (CJZ18);
  - those degron cells treated with PT-179 for 7 days to degrade UBA5 (CJZ18_PT).
Each background was grown in full medium (df) or starved in EBSS for 20 h, with or without overnight
bafilomycin A1 (BafA1), 6 wells per arm. The input is the FlowJo export: per-well median neutral
(VL3) and acidic (YL2) mKeima fluorescence of BFP+ single cells.

  1. Each well's median is divided by the mean of the untreated wells of the same background and
     medium.
  2. BafA1 effect per background x medium = mean(BafA1) - mean(untreated) of the normalized values;
     error bar = combined SD, sqrt(var_BafA1 + var_untreated).
  3. BafA1 vs untreated: Welch's t-test for each background x medium, Benjamini-Hochberg-adjusted
     across the 6.
  4. Between background x medium groups: Welch's t-test on the BafA1 effects from their summary
     statistics, Benjamini-Hochberg-adjusted across the 15 pairs.

Outputs (output root, flow_cytometry/):
  fig6C_acidic_bafA1_effect.pdf      acidic signal (Fig 6C)
  figS7D_neutral_bafA1_effect.pdf    neutral signal (Fig S7D)
  fig6C_figS7D_bafA1_effect.csv      effects, combined SDs, BafA1 vs untreated tests
  fig6C_figS7D_between_groups.csv    between-group tests

    python flow_cytometry/fig6C_figS7D_mkeima.py
"""
import os
import sys
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import paths  # noqa: E402

DATA = ('flow_cytometry', '6C_S7D', 'DataWorkup_L29_mKeima_CJZ18.xlsx')
OUT_DIR = 'flow_cytometry'
SIGNALS = {'acidic': 'Acidic', 'neutral': 'Neutral'}       # column prefix -> plot label
UNTREATED, TREATED = 0, '16 nM'                             # BafA column values in the export
PANELS = {'acidic': 'fig6C_acidic_bafA1_effect.pdf', 'neutral': 'figS7D_neutral_bafA1_effect.pdf'}


def out(name):
    return paths.outp(OUT_DIR, name, mkdir=True)


def load():
    """Per-well medians, normalized to the untreated mean of the same background and medium."""
    df = pd.read_excel(paths.inp(*DATA))
    df['group'] = df['condition'] + '_' + df['Media']
    untreated = df['BafA'] == UNTREATED
    for sig in SIGNALS:
        ref = df[untreated].groupby('group')[f'{sig}_median'].mean()
        df[f'{sig}_norm'] = df[f'{sig}_median'] / df['group'].map(ref)
    return df


def bafA1_effects(df):
    """BafA1 effect, combined SD and BafA1-vs-untreated Welch test per background x medium."""
    treated, untreated = df[df['BafA'] == TREATED], df[df['BafA'] == UNTREATED]
    cols = [f'{sig}_norm' for sig in SIGNALS]
    t_agg = treated.groupby('group')[cols].agg(['mean', 'var', 'count'])
    u_agg = untreated.groupby('group')[cols].agg(['mean', 'var', 'count'])
    eff = pd.DataFrame({'group': t_agg.index,
                        'n_bafA1': t_agg[(cols[0], 'count')].values,
                        'n_untreated': u_agg.loc[t_agg.index, (cols[0], 'count')].values})
    for sig in SIGNALS:
        c = f'{sig}_norm'
        eff[f'{sig}_effect'] = (t_agg[(c, 'mean')] - u_agg.loc[t_agg.index, (c, 'mean')]).values
        eff[f'{sig}_effect_var'] = (t_agg[(c, 'var')] + u_agg.loc[t_agg.index, (c, 'var')]).values
        eff[f'{sig}_effect_sd'] = np.sqrt(eff[f'{sig}_effect_var'])
        tests = [stats.ttest_ind(treated.loc[treated.group == g, c], untreated.loc[untreated.group == g, c],
                                 equal_var=False) for g in eff.group]
        eff[f'{sig}_t'] = [r.statistic for r in tests]
        eff[f'{sig}_p'] = [r.pvalue for r in tests]
        eff[f'{sig}_p_bh'] = multipletests(eff[f'{sig}_p'], method='fdr_bh')[1]
    return eff


def between_groups(eff):
    """Welch's t-test between the BafA1 effects of every pair of groups (each effect's SD is the
    combined SD, n = wells per arm), Benjamini-Hochberg-adjusted per signal."""
    rows = []
    for (_, a), (_, b) in combinations(eff.iterrows(), 2):
        row = {'group1': a.group, 'group2': b.group}
        for sig in SIGNALS:
            m1, v1, m2, v2 = a[f'{sig}_effect'], a[f'{sig}_effect_var'], b[f'{sig}_effect'], b[f'{sig}_effect_var']
            r = stats.ttest_ind_from_stats(m1, np.sqrt(v1), a.n_bafA1, m2, np.sqrt(v2), b.n_bafA1, equal_var=False)
            row.update({f'{sig}_effect_diff': m1 - m2, f'{sig}_effect_diff_var': v1 + v2,
                        f'{sig}_t': r.statistic, f'{sig}_p': r.pvalue})
        rows.append(row)
    pairs = pd.DataFrame(rows)
    for sig in SIGNALS:
        pairs[f'{sig}_p_bh'] = multipletests(pairs[f'{sig}_p'], method='fdr_bh')[1]
    return pairs


def effect_bars(eff, sig, path):
    """BafA1 - untreated difference of normalized medians per group; error bars = combined SD."""
    label = SIGNALS[sig]
    plt.figure(figsize=(10, 6))
    plt.bar(eff['group'], eff[f'{sig}_effect'], yerr=eff[f'{sig}_effect_sd'], capsize=5)
    plt.xticks(rotation=90)
    plt.xlabel('Condition')
    plt.ylabel(f'Difference in Mean {label} mKeima Signal (BafA - no BafA)')
    plt.title(f'Difference in Mean {label} mKeima Signal Between BafA and no BafA Conditions')
    plt.savefig(path, bbox_inches='tight')
    plt.close('all')


def main():
    eff = bafA1_effects(load())
    eff.to_csv(out('fig6C_figS7D_bafA1_effect.csv'), index=False)
    between_groups(eff).to_csv(out('fig6C_figS7D_between_groups.csv'), index=False)
    for sig, name in PANELS.items():
        effect_bars(eff, sig, out(name))
    print(eff[['group', 'acidic_effect', 'acidic_p_bh', 'neutral_effect', 'neutral_p_bh']].to_string(index=False))
    print('wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
