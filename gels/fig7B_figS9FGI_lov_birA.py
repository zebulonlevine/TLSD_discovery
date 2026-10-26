"""Fig 7B, S9F, S9G, S9I: ER-ribosome proximity labelling with LOV-BirA, read out on gels.

RPL29-HaloTag-AviTag cells with ER-targeted LOV-BirA (six genotypes: sgNT, sgUBA5, sgTLSD1,
sgUBA5 sgTLSD1, and TLSD1 overexpression over sgNT or sgUBA5) were labelled with JF549 for 1 h and
biotinylated for 15 min in blue light (L) or kept dark (D, biotin background), then lysed at once
(day 0) or after 5 days (light only). Lanes were spread over 8 gels (one replicate of 3 conditions
per genotype per gel). Input: ImageStudio band intensities with local background (sheet `gelRedo`):
in-gel JF549 fluorescence, streptavidin (biotin) and anti-HaloTag blots, and tubulin from the
re-probed membrane. Each signal has an intact upper band, a ~40 kDa intermediate band (post-lysis
cleavage) and, in day-5 lanes, a ~30 kDa lysosomal HaloTag band (lower); the sheet also holds each
band divided by the lane's tubulin (`*_per_tubulin`).

  1. Intact signal = upper + intermediate band, per tubulin.
  2. Per-gel normalization: each tubulin-normalized signal is divided by its mean over all lanes of
     the same gel.
  3. Panel quantities:
       7B   biotin / JF549 of the same lane, / mean of sgNT day-0 light lanes
       S9F  JF549, / mean of all sgNT day-0 lanes (light and dark)
       S9G  biotin, / mean of sgNT day-0 light lanes
       S9I  JF549 lower band / JF549 intact (raw intensities of the same lane), day 5 only,
            / mean of sgNT day 5
  4. Mean and 95% CI (t distribution, n - 1 df) per group.
  5. Welch's t-tests, Benjamini-Hochberg adjusted within each family:
       7B         day 0 vs day 5 (light) within each genotype                            6 tests
       S9F, S9G   within genotype: day-0 light vs day 5, day 5 vs day-0 dark             12 tests
                  across genotypes, day-0 light lanes                                     15 tests
                  across genotypes, day-5 lanes                                           15 tests
       S9I        all pairs of genotypes, day 5                                           15 tests

The sheet includes a repeat run of one sgNT and one sgUBA5 replicate (gel 8, column `extra` = Y),
so those two genotypes have n = 7 and the others n = 6. The original analysis and the placed panels
include these lanes. `--exclude-repeat-lanes` drops them before step 2 and writes to
gels/exclude_repeat_lanes/ instead.

Points are jittered with a fixed seed. Significance brackets are not drawn; the stars are in the
*_welch.csv tables.

Outputs (output root, gels/):
  fig7B_er_fraction.pdf, figS9F_jf549.pdf, figS9G_biotin.pdf, figS9I_ribophagy.pdf
  fig7B_figS9FGI_values.csv    per-lane values of the four panel quantities
  fig7B_figS9FGI_summary.csv   mean, SEM, 95% CI and n per panel and group
  fig7B_figS9FGI_welch.csv     all tests: raw and BH-adjusted p, stars

    python gels/fig7B_figS9FGI_lov_birA.py [--exclude-repeat-lanes]
"""
import argparse
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

DATA_FILE = ('gels', '7B_S9FGI', 'BirAgels.xlsx')
OUT_DIR = 'gels'

# genotype names in the sheet (R3HDM4 = TLSD1; "UBC R3HDM4" = TLSD1 overexpression) -> panel labels
GENOTYPES = {
    'sgNT': 'sgNT',
    'sgUBA5': 'sgUBA5',
    'sgR3HDM4': 'sgTLSD1',
    'sgUBA5/sgR3HDM4': 'sgUBA5 sgTLSD1',
    'UBC R3HDM4 sgNT': 'sgNT oeTLSD1',
    'UBC R3HDM4 sgUBA5': 'sgUBA5 oeTLSD1',
}
TIMING_LABEL = {'L': 'L 0 d', '5d': 'L 5 d', 'D': 'D 0 d'}

# panel -> (value column, timings shown, y label, y limits, point colour)
PANELS = {
    '7B': ('er_fraction', ['L', '5d'], 'Relative Fraction ER 60S\n(Biotin/JF549)', (0.4, 1.8), 'grey'),
    'S9F': ('jf549', ['L', '5d', 'D'], 'JF549 Levels (normalized to sgNT)', (0, 1.8), '#EE3B3B'),
    'S9G': ('biotin', ['L', '5d', 'D'], 'Biotinylated RPL29-HT-Avi Levels', (0, 1.8), 'grey'),
    'S9I': ('ribophagy', ['5d'], 'Relative Ribophagy Flux\n(fraction lower/upper)', None, 'grey'),
}
FILES = {'7B': 'fig7B_er_fraction.pdf', 'S9F': 'figS9F_jf549.pdf',
         'S9G': 'figS9G_biotin.pdf', 'S9I': 'figS9I_ribophagy.pdf'}
STYLE = {'figure.dpi': 100, 'figure.facecolor': 'white', 'pdf.fonttype': 42,
         'axes.spines.top': False, 'axes.spines.right': False}


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------

def load(exclude_repeat_lanes=False):
    df = pd.read_excel(paths.inp(*DATA_FILE), sheet_name='gelRedo')
    if exclude_repeat_lanes:
        df = df[df['extra'] != 'Y'].reset_index(drop=True)
    return df


def quantify(df):
    """Steps 1-3: per-lane panel quantities."""
    d = df.copy()
    for sig in ('Biotin', 'JF549'):
        d[f'{sig}_intact_perTubulin'] = d[f'{sig}Upper_per_tubulin'] + d[f'{sig}Middle_per_tubulin']
        col = f'{sig}_intact_perTubulin'
        d[col] = d[col] / d.groupby('gel')[col].transform('mean')
    d['JF549_intact'] = d['JF549Upper'] + d['JF549Middle']
    d['group'] = d['genetics'] + '_' + d['timing']

    def rel(col, ref_groups):
        return d[col] / d[col][d['group'].isin(ref_groups)].mean()

    d['er_fraction'] = d['Biotin_intact_perTubulin'] / d['JF549_intact_perTubulin']
    d['er_fraction'] = rel('er_fraction', ['sgNT_L'])
    d['jf549'] = rel('JF549_intact_perTubulin', ['sgNT_L', 'sgNT_D'])
    d['biotin'] = rel('Biotin_intact_perTubulin', ['sgNT_L'])
    d['ribophagy'] = d['JF549Lower'] / d['JF549_intact']
    d['ribophagy'] = rel('ribophagy', ['sgNT_5d'])
    d['genotype'] = d['genetics'].map(GENOTYPES)
    return d


def summarize(d):
    rows = []
    for panel, (col, timings, *_rest) in PANELS.items():
        for gen in GENOTYPES:
            for t in timings:
                v = d.loc[(d['genetics'] == gen) & (d['timing'] == t), col]
                n, m, se = v.count(), v.mean(), v.sem()
                half = se * stats.t.ppf(0.975, df=n - 1) if n > 1 else 0.0
                rows.append({'panel': panel, 'genotype': GENOTYPES[gen], 'timing': TIMING_LABEL[t], 'n': n,
                             'mean': m, 'sem': se, 'ci95_lo': m - half, 'ci95_hi': m + half})
    return pd.DataFrame(rows)


def _stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def families():
    """panel -> {family: [((genetics, timing), (genetics, timing)), ...]} (step 5)."""
    gens = list(GENOTYPES)
    within = [p for g in gens for p in (((g, 'L'), (g, '5d')), ((g, '5d'), (g, 'D')))]
    across = {t: [((a, t), (b, t)) for a, b in combinations(gens, 2)] for t in ('L', '5d')}
    fam_FG = {'within genotype': within, 'day 0 across genotypes': across['L'],
              'day 5 across genotypes': across['5d']}
    return {'7B': {'day 0 vs day 5': [((g, 'L'), (g, '5d')) for g in gens]},
            'S9F': fam_FG, 'S9G': fam_FG,
            'S9I': {'day 5 across genotypes': across['5d']}}


def welch_tests(d):
    rows = []
    for panel, fams in families().items():
        col = PANELS[panel][0]
        for fam, pairs in fams.items():
            block = []
            for (g1, t1), (g2, t2) in pairs:
                a = d.loc[(d['genetics'] == g1) & (d['timing'] == t1), col].dropna()
                b = d.loc[(d['genetics'] == g2) & (d['timing'] == t2), col].dropna()
                r = stats.ttest_ind(a, b, equal_var=False)
                block.append({'panel': panel, 'family': fam,
                              'group1': f'{GENOTYPES[g1]} {TIMING_LABEL[t1]}',
                              'group2': f'{GENOTYPES[g2]} {TIMING_LABEL[t2]}',
                              'n1': len(a), 'n2': len(b), 'mean1': a.mean(), 'mean2': b.mean(),
                              't': r.statistic, 'df': r.df, 'p': r.pvalue})
            block = pd.DataFrame(block)
            block['p_bh'] = multipletests(block['p'], method='fdr_bh')[1]
            block['sig'] = block['p_bh'].map(_stars)
            rows.append(block)
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def draw_strip(data, x, y, order, ax, point_color='grey', jitter=0.25, seed=0):
    """Points (jittered), black mean and 95% CI (t) per category of `x`, in `order`."""
    pos = data[x].map({g: i for i, g in enumerate(order)}).astype(float)
    rng = np.random.default_rng(seed)
    ax.scatter(pos + rng.uniform(-jitter, jitter, size=len(data)), data[y], c=point_color,
               edgecolors='none', marker='o', alpha=0.6, s=25, zorder=1, linewidth=0.5)
    summary = data.groupby(x)[y].agg(['mean', 'sem', 'count']).reindex(order)
    for i, grp in enumerate(order):
        m, sem, n = summary.loc[grp, 'mean'], summary.loc[grp, 'sem'], summary.loc[grp, 'count']
        ax.hlines(m, i - 0.25, i + 0.25, color='black', linewidth=2.5, zorder=3)
        tcrit = stats.t.ppf(0.975, df=n - 1) if n > 1 else 0.0
        ax.errorbar(i, m, yerr=sem * tcrit, fmt='none', ecolor='black', elinewidth=2.5, capsize=4, zorder=3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_xlim(-0.5 - 0.15, len(order) - 0.5 + 0.15)


def panel_data(d, panel):
    """Rows and category order for one panel (x = 'group', i.e. genetics_timing)."""
    col, timings, *_rest = PANELS[panel]
    order = [f'{g}_{t}' for g in GENOTYPES for t in timings]
    return d[d['timing'].isin(timings)], order


FIGSIZE = {'7B': (6.4, 3.2), 'S9F': (7.2, 3.4), 'S9G': (7.2, 3.4), 'S9I': (3.6, 3.4)}


def new_figure(panel):
    return plt.subplots(figsize=FIGSIZE[panel])


def decorate(ax, panel):
    col, timings, ylabel, ylim, _color = PANELS[panel]
    labels = [f'{GENOTYPES[g]}\n{TIMING_LABEL[t]}' if len(timings) > 1 else GENOTYPES[g]
              for g in GENOTYPES for t in timings]
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)


def plot_panel(d, panel, path):
    data, order = panel_data(d, panel)
    with mpl.rc_context(STYLE):
        fig, ax = new_figure(panel)
        draw_strip(data, 'group', PANELS[panel][0], order, ax, point_color=PANELS[panel][4])
        decorate(ax, panel)
        fig.savefig(path, bbox_inches='tight')
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--exclude-repeat-lanes', action='store_true',
                    help='drop the repeat-run sgNT/sgUBA5 lanes (extra = Y); writes to gels/exclude_repeat_lanes/')
    args = ap.parse_args()
    sub = (OUT_DIR, 'exclude_repeat_lanes') if args.exclude_repeat_lanes else (OUT_DIR,)

    def out(name):
        return paths.outp(*sub, name, mkdir=True)

    d = quantify(load(args.exclude_repeat_lanes))
    cols = ['number', 'gel', 'rep', 'extra', 'genotype', 'timing'] + [p[0] for p in PANELS.values()]
    d[cols].to_csv(out('fig7B_figS9FGI_values.csv'), index=False)
    summarize(d).to_csv(out('fig7B_figS9FGI_summary.csv'), index=False)
    tests = welch_tests(d)
    tests.to_csv(out('fig7B_figS9FGI_welch.csv'), index=False)
    for panel in PANELS:
        plot_panel(d, panel, out(FILES[panel]))
    sig = tests[tests['sig'] != 'ns']
    print(sig[['panel', 'family', 'group1', 'group2', 'p_bh', 'sig']].to_string(index=False))
    print('wrote', paths.outp(*sub))


if __name__ == '__main__':
    main()
