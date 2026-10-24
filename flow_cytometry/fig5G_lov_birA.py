"""Optogenetic proximity labeling of ER ribosomes, flow cytometry (Fig 5G).

RPE1 CRISPRi cells (sgNT or sgUBA5) expressing RPL29-HaloTag-AviTag and ER-tethered LOV-BirA-Sec61B
were labeled with HaloTag dye (JF549 and/or JF646) and biotin under 20 min blue light (or kept dark),
then collected immediately or after a 4-day chase. Streptavidin-AF568 stains the biotinylated (ER)
ribosomes. Plate LOV-BirA_test1; the over-fixed day-3 samples are not analyzed.

  1. Gating (FCS files; --fcs-dir): compensation with the acquisition spillover matrix; cells
     (FSC-H 150,000-600,000, SSC-H 40,000-600,000); singlets (FSC-W 5th-90th percentile); GFP+
     (guide vector, Comp-BL1-A >= 2,000) and BFP+ (Comp-VL1-A >= 400). Per-well medians of
     Comp-YL1-A (JF549), Comp-YL2-A (streptavidin-AF568) and Comp-RL1-A (JF646) for GFP+BFP+
     singlets. Wells with <= 1,000 such cells are dropped (none were).
  2. Biotin, per genotype and HaloTag-dye condition: B = (M - D) / (L - D), D and L the mean
     medians of the matched no-light and immediate-collection wells. The JF549-pulse/JF646-chase
     wells (day 4 only) use the mean over all dye conditions.
  3. HaloTag, per genotype and illumination: background = mean of wells not stained with the dye;
     the background-subtracted signal is divided by the mean of wells stained with that dye alone
     (day 0 and day 4 pooled). The JF549 panel shows illuminated wells only.
  4. sgNT vs sgUBA5 at the 4d chase: t-test with the error propagated (delta method) from every
     plotted well, Welch-Satterthwaite df (the stars in the panel); the per-well Welch test on the
     4d-chase values is reported alongside.

Input (data root, flow_cytometry/5G/): channel_statistics_dark.csv, the step-1 per-well table
(columns filename, type, staining, Genetics, *_count, *_pct, Comp-YL1-A_median, Comp-YL2-A_median,
Comp-RL1-A_median). With --fcs-dir and --metadata the table is rebuilt from the FCS files and
written to fig5G_channel_statistics.csv first.

Outputs (output root, flow_cytometry/):
  fig5G_flow.pdf / .png                both panels; fig5G_JF549.pdf, fig5G_biotin.pdf single panels
  fig5G_per_well.csv                   per-well raw medians and normalized values
  fig5G_summary_n.csv                  n, mean, SD, SEM per plotted group
  fig5G_welch_tests.csv                tests
  fig5G_normalization_checks.csv       illumination effect on HaloTag; day-4 vs day-0 fully labeled

    python flow_cytometry/fig5G_lov_birA.py
    python flow_cytometry/fig5G_lov_birA.py --fcs-dir <Zeb251125/dark> --metadata <dark_metadata.csv>
"""
import argparse
import contextlib
import io
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
if 'ipykernel' not in sys.modules:
    mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402

DATA = ('flow_cytometry', '5G', 'channel_statistics_dark.csv')
OUT_DIR = 'flow_cytometry'

GENOTYPES = {203: 'sgNT', 198: 'sgUBA5'}          # plasmid numbers in metadata 'Genetics'
COLORS = {'sgNT': '#6d6e71', 'sgUBA5': '#d96227'}  # Figure 5 colors
MIN_CELLS = 1000                                   # minimum GFP+BFP+ singlets per well

# HaloTag dye condition encoded in the metadata 'staining' column (all wells also get streptavidin-AF568)
PC = 'JF549 pulse / JF646 chase'
DYE = {'SA-AF568': 'none', 'SA-AF568_JF549': 'JF549', 'SA-AF568_JF646': 'JF646', 'SA-AF568_4dJF549_JF646New': PC}
# collection day and light condition encoded in 'type'
TIMING = {'Light20min': ('d0', 'light'), 'Dark20min': ('d0', 'dark'), '4dChase': ('d4', 'light')}
CHANNELS = {'Comp-YL1-A_median': 'JF549_raw', 'Comp-YL2-A_median': 'biotin_raw',
            'Comp-RL1-A_median': 'JF646_raw', 'GFP_BFP+_count': 'n_cells'}

PANELS = {
    'JF549': dict(value='JF549_norm', group='JF549_group', order=['fully labeled', 'unlabeled', '4d chase'],
                  title='Total old 60S ribosomes', ylabel='Normalized Halo-JF549', primary='propagated'),
    'biotin': dict(value='biotin_norm', group='biotin_group', order=['immediate collection', 'no light', '4d chase'],
                   title='ER-biotinylated old 60S ribosomes', ylabel='Normalized streptavidin staining', primary='propagated'),
}
REF_DYES = ['none', 'JF549', 'JF646']
YLIM = {'JF549': (-0.15, 1.25), 'biotin': (-0.45, 1.6)}
# Illustrator-editable output: TrueType (live) text, no clipping masks, no transparency, no background boxes
STYLE = {'font.family': 'Arial', 'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
         'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
         'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False,
         'axes.spines.right': False}


def out(name):
    return paths.outp(OUT_DIR, name, mkdir=True)


# ---------------------------------------------------------------------------
# step 1: FCS -> per-well medians (needs the raw FCS files)
# ---------------------------------------------------------------------------

def channel_statistics_from_fcs(fcs_dir, metadata_csv):
    """Compensate, gate and return per-well medians of GFP+BFP+ singlets for the analyzed wells."""
    from fcs_gating import FCAnalyzer
    fcs_dir = Path(fcs_dir)
    meta = pd.read_csv(metadata_csv)
    meta = meta.loc[meta.Genetics.isin(GENOTYPES) & meta.type.isin(TIMING) & meta.staining.isin(DYE)]
    meta = meta.loc[[(fcs_dir / f).exists() for f in meta.filename]]
    print(len(meta), 'wells selected')

    analyzer = FCAnalyzer(output_dir=paths.outp(OUT_DIR, mkdir=True))
    with contextlib.redirect_stdout(io.StringIO()):   # silence per-file load messages
        analyzer.load_files([str(fcs_dir / f) for f in meta.filename], metadata_csv=metadata_csv)
    # all wells carry the same acquisition spillover matrix, so the first file's is used
    spill = {fa.metadata['_fcs_meta'].get('$SPILLOVER') for fa in analyzer.files.values()}
    assert len(spill) == 1, 'spillover matrix differs between files'
    analyzer.apply_compensation()

    analyzer.gate_cells(fsc='FSC-H', ssc='SSC-H', fsc_min=150000, fsc_max=600000,
                        ssc_min=40000, ssc_max=600000, method='rectangle', name='Cells')
    analyzer.gate_singlets(area='FSC-A', height='FSC-H', width='FSC-W', method='width',
                           width_percentile_low=5, width_percentile_high=90, name='Singlets')
    analyzer.gate_fluorophore(channel='Comp-BL1-A', threshold=2000, positive_name='GFP+', negative_name='GFP-')
    analyzer.gate_fluorophore(channel='Comp-VL1-A', threshold=400, positive_name='BFP+', negative_name='BFP-')
    analyzer.combine_gates(['GFP+', 'BFP+'], operation='AND', name='GFP_BFP+')
    return analyzer.get_channel_statistics(
        channels=['Comp-YL1-A', 'Comp-YL2-A', 'Comp-RL1-A'], gate=['GFP_BFP+'], stats=['median'],
        include_metadata=['type', 'staining', 'Genetics'], include_counts=True)


# ---------------------------------------------------------------------------
# steps 2-3: per-well normalization
# ---------------------------------------------------------------------------

def load_wells(table):
    """Analyzed wells of the per-well table, with their design annotations."""
    table = table.loc[table.Genetics.isin(GENOTYPES) & table.type.isin(TIMING) & table.staining.isin(DYE)]
    for g in ['Cells', 'Singlets']:
        print(f"{g}: {table[g + '_pct'].min():.1f}-{table[g + '_pct'].max():.1f}% of parent")
    wells = table[['filename', 'type', 'staining', 'Genetics', 'total_events'] + list(CHANNELS)]
    wells = wells.rename(columns=CHANNELS)
    wells = wells[['filename', 'type', 'staining', 'Genetics', 'total_events', 'n_cells',
                   'JF549_raw', 'biotin_raw', 'JF646_raw']].reset_index(drop=True)
    wells['genotype'] = wells.Genetics.map(GENOTYPES)
    wells['dye'] = wells.staining.map(DYE)
    wells['day'] = wells.type.map(lambda t: TIMING[t][0])
    wells['light'] = wells.type.map(lambda t: TIMING[t][1])
    wells['has_JF549'] = wells.dye.isin(['JF549', PC])
    wells['has_JF646'] = wells.dye.isin(['JF646', PC])

    low = wells.n_cells <= MIN_CELLS
    print(f'{len(wells)} wells; GFP+BFP+ singlets per well: {wells.n_cells.min():,}-{wells.n_cells.max():,} '
          f'(median {wells.n_cells.median():,.0f}); {low.sum()} wells with <= {MIN_CELLS} cells excluded')
    return wells.loc[~low].reset_index(drop=True)


def normalize_biotin(wells):
    """B = (M - D) / (L - D) with the matched no-light (D) and immediate-collection (L) means."""
    d0 = wells.loc[wells.day == 'd0']
    refs = {}   # (genotype, dye) -> (no-light mean D, immediate-collection mean L)
    for (g, dye), sub in d0.groupby(['genotype', 'dye']):
        refs[(g, dye)] = (sub.loc[sub.light == 'dark', 'biotin_raw'].mean(), sub.loc[sub.light == 'light', 'biotin_raw'].mean())
    for g, sub in d0.groupby('genotype'):
        refs[(g, PC)] = (sub.loc[sub.light == 'dark', 'biotin_raw'].mean(), sub.loc[sub.light == 'light', 'biotin_raw'].mean())

    D, L = np.array([refs[k] for k in zip(wells.genotype, wells.dye)]).T
    wells['biotin_bg'] = D
    wells['biotin_norm'] = (wells.biotin_raw - D) / (L - D)
    wells['biotin_group'] = wells.type.map({'Light20min': 'immediate collection', 'Dark20min': 'no light',
                                            '4dChase': '4d chase'})
    return wells


def normalize_halotag(wells):
    """Background-subtract and normalize JF549 and JF646, matched for genotype and illumination."""
    for dye, flag in [('JF549', 'has_JF549'), ('JF646', 'has_JF646')]:
        key = pd.MultiIndex.from_frame(wells[['genotype', 'light']])
        bg = wells.loc[~wells[flag]].groupby(['genotype', 'light'])[f'{dye}_raw'].mean()
        wells[f'{dye}_bgsub'] = wells[f'{dye}_raw'] - bg.reindex(key).values
        ref = wells.loc[wells.dye == dye].groupby(['genotype', 'light'])[f'{dye}_bgsub'].mean()
        wells[f'{dye}_norm'] = wells[f'{dye}_bgsub'] / ref.reindex(key).values

    wells['JF549_group'] = np.where(
        wells.light == 'light',
        np.select([wells.dye == PC, wells.has_JF549], ['4d chase', 'fully labeled'], 'unlabeled'),
        None)
    return wells


def normalization_checks(wells):
    """The checks behind step 3: blue light lowers day-0 HaloTag signal (so references are matched for
    illumination), and fully labeled JF549 signal is stable from day 0 to day 4 (so both are pooled)."""
    rows = []
    d0 = wells.loc[wells.day == 'd0']
    for dye in ['JF549', 'JF646']:
        bg = d0.loc[d0.dye == 'none'].groupby('genotype')[f'{dye}_raw'].mean()
        for g, sub in d0.loc[d0.dye == dye].groupby('genotype'):
            v = sub[f'{dye}_raw'] - bg[g]
            lit, dark = v[sub.light == 'light'], v[sub.light == 'dark']
            rows.append(dict(check=f'{dye}, day 0: illuminated / dark', genotype=g, ratio=lit.mean() / dark.mean(),
                             p_welch=stats.ttest_ind(lit, dark, equal_var=False).pvalue, n1=len(lit), n2=len(dark)))
    fl = wells.loc[wells.JF549_group == 'fully labeled']
    for g, sub in fl.groupby('genotype'):
        a, b = sub.loc[sub.day == 'd0', 'JF549_bgsub'], sub.loc[sub.day == 'd4', 'JF549_bgsub']
        rows.append(dict(check='JF549 fully labeled: day 4 / day 0', genotype=g, ratio=b.mean() / a.mean(),
                         p_welch=stats.ttest_ind(a, b, equal_var=False).pvalue, n1=len(a), n2=len(b)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# step 4: statistics
# ---------------------------------------------------------------------------

def welch(a, b):
    r = stats.ttest_ind(a, b, equal_var=False)
    return dict(estimate=a.mean() - b.mean(), se=np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)),
                t=r.statistic, df=r.df, p=r.pvalue)


def propagated(F, groups):
    # delta method: F at the group means; each group's squared SEM is weighted by (dF/d mean)^2
    means = {k: np.mean(v) for k, v in groups.items()}
    parts = []
    for k, v in groups.items():
        h = 1e-6 * max(abs(means[k]), 1.0)
        grad = (F({**means, k: means[k] + h}) - F({**means, k: means[k] - h})) / (2 * h)
        parts.append((grad ** 2 * np.var(v, ddof=1) / len(v), len(v) - 1))
    return F(means), parts


def propagated_welch(fit1, fit2):
    # t-test on F1 - F2 with Welch-Satterthwaite df over all variance components of both genotypes
    (f1, parts1), (f2, parts2) = fit1, fit2
    u, dof = np.array(parts1 + parts2).T
    se = np.sqrt(u.sum())
    t = (f1 - f2) / se
    df = u.sum() ** 2 / np.sum(u ** 2 / dof)
    return dict(estimate=f1 - f2, se=se, t=t, df=df, p=2 * stats.t.sf(abs(t), df))


def jf549_groups(wells, g):
    # every plotted JF549 well: 4d chase, fully labeled (reference), unlabeled (background)
    sub = wells.loc[wells.genotype == g]
    return {grp: sub.loc[sub.JF549_group == grp, 'JF549_raw'].to_numpy() for grp in PANELS['JF549']['order']}


def jf549_F(groups):
    return lambda m: (m['4d chase'] - m['unlabeled']) / (m['fully labeled'] - m['unlabeled'])


def biotin_groups(wells, g):
    # every plotted biotin well, split by the dye condition that sets its reference
    sub = wells.loc[wells.genotype == g]
    return {(grp, d): v.to_numpy() for (grp, d), v in sub.groupby(['biotin_group', 'dye']).biotin_raw}


def biotin_F(groups):
    n = {d: len(v) for (grp, d), v in groups.items() if grp == '4d chase'}
    def F(m):
        D = {d: m[('no light', d)] for d in REF_DYES}
        L = {d: m[('immediate collection', d)] for d in REF_DYES}
        D[PC], L[PC] = np.mean([D[d] for d in REF_DYES]), np.mean([L[d] for d in REF_DYES])
        return sum(n[d] * (m[('4d chase', d)] - D[d]) / (L[d] - D[d]) for d in n) / sum(n.values())
    return F


GROUPS = {'JF549': (jf549_groups, jf549_F), 'biotin': (biotin_groups, biotin_F)}


def statistics(wells):
    """Per-group n / mean / SD / SEM, and the sgNT vs sgUBA5 tests at the 4d chase."""
    summary, tests = [], []
    for panel, p in PANELS.items():
        for grp in p['order']:
            sub = wells.loc[wells[p['group']] == grp]
            for g in GENOTYPES.values():
                v = sub.loc[sub.genotype == g, p['value']]
                comp = sub.loc[sub.genotype == g].groupby(['day', 'light', 'dye']).size()
                summary.append(dict(panel=panel, condition=grp, genotype=g, n_wells=len(v), mean=v.mean(), sd=v.std(),
                                    sem=v.sem(), wells=', '.join(f'{d}/{l}/{dy}: {k}' for (d, l, dy), k in comp.items())))

        chase = wells.loc[wells[p['group']] == '4d chase']
        a, b = (chase.loc[chase.genotype == g, p['value']] for g in ['sgNT', 'sgUBA5'])
        tests.append(dict(panel=panel, comparison='4d chase: sgNT vs sgUBA5', method='Welch, per-well 4d-chase values',
                          primary=p['primary'] == 'per-well', n1=len(a), n2=len(b), **welch(a, b)))

        get_groups, make_F = GROUPS[panel]
        fits, n_used = [], []
        for g, plotted in [('sgNT', a), ('sgUBA5', b)]:
            groups = get_groups(wells, g)
            fits.append(propagated(make_F(groups), groups))
            n_used.append(sum(len(v) for v in groups.values()))
            assert np.isclose(fits[-1][0], plotted.mean())   # F from group means = mean of the plotted 4d-chase points
        tests.append(dict(panel=panel, comparison='4d chase: sgNT vs sgUBA5',
                          method='Welch-Satterthwaite, error propagated from all plotted wells',
                          primary=p['primary'] == 'propagated', n1=n_used[0], n2=n_used[1], **propagated_welch(*fits)))

    # supporting: does 4d-chase biotin differ from the no-light baseline within each genotype?
    for g in GENOTYPES.values():
        sub = wells.loc[wells.genotype == g]
        a, b = sub.loc[sub.biotin_group == '4d chase', 'biotin_norm'], sub.loc[sub.biotin_group == 'no light', 'biotin_norm']
        tests.append(dict(panel='biotin', comparison=f'{g}: 4d chase vs no light', method='Welch, per-well values',
                          primary=False, n1=len(a), n2=len(b), **welch(a, b)))
    return pd.DataFrame(summary), pd.DataFrame(tests)


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------

def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'


def dot_panel(ax, wells, tests, panel, ylim, seed=0):
    """Per-well values (jittered), group means (black lines), n per column, and the primary 4d-chase test."""
    p, rng = PANELS[panel], np.random.default_rng(seed)
    ax.patch.set_visible(False)
    ax.figure.patch.set_visible(False)
    offsets = dict(zip(GENOTYPES.values(), [-0.19, 0.19]))
    for i, grp in enumerate(p['order']):
        for g, dx in offsets.items():
            v = wells.loc[(wells[p['group']] == grp) & (wells.genotype == g), p['value']].to_numpy()
            x = i + dx + rng.uniform(-0.06, 0.06, len(v))
            ax.scatter(x, v, s=16, color=COLORS[g], edgecolor='white', linewidth=0.4, zorder=2, clip_on=False)
            ax.plot([i + dx - 0.13, i + dx + 0.13], [v.mean()] * 2, color='black', lw=1, zorder=3, clip_on=False)
            ax.text(i + dx, ylim[0] + 0.015 * (ylim[1] - ylim[0]), f'{len(v)}', ha='center', va='bottom',
                    fontsize=5.5, color='#555555')
    ax.text(-0.58, ylim[0] + 0.015 * (ylim[1] - ylim[0]), 'n =', ha='left', va='bottom', fontsize=5.5, color='#555555')
    # 4d chase: sgNT vs sgUBA5, primary test for this panel
    t = tests.loc[(tests.panel == panel) & tests.primary].iloc[0]
    top = wells.loc[wells[p['group']] == '4d chase', p['value']].max()
    y = top + 0.06 * (ylim[1] - ylim[0])
    i = p['order'].index('4d chase')
    ax.plot([i + offsets['sgNT'], i + offsets['sgUBA5']], [y, y], color='black', lw=0.6, clip_on=False)
    ax.text(i, y, stars(t.p), ha='center', va='bottom', fontsize=7 if t.p < 0.05 else 6)
    ax.axhline(0, color='#bbbbbb', lw=0.5, zorder=1, clip_on=False)
    ax.set_xticks(range(len(p['order'])), [s.replace(' ', '\n', 1) if len(s) > 10 else s for s in p['order']])
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.6, len(p['order']) - 0.4)
    ax.set_ylim(*ylim)
    ax.set_ylabel(p['ylabel'])
    ax.set_title(p['title'], fontsize=7)
    handles = [Line2D([], [], ls='', marker='o', ms=4, mfc=c, mec='white', mew=0.4, label=g)
               for g, c in COLORS.items()]
    ax.legend(handles=handles, frameon=False, loc='upper center', handletextpad=0.1, borderaxespad=0.2)


def figure(wells, tests):
    mpl.rcParams.update(STYLE)
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.3), constrained_layout=True)
    for ax, panel in zip(axes, PANELS):
        dot_panel(ax, wells, tests, panel, YLIM[panel])
    fig.savefig(out('fig5G_flow.pdf'))
    fig.savefig(out('fig5G_flow.png'), dpi=300)
    for panel in PANELS:
        f, ax = plt.subplots(figsize=(2.3, 2.3), constrained_layout=True)
        dot_panel(ax, wells, tests, panel, YLIM[panel])
        f.savefig(out(f'fig5G_{panel}.pdf'))
        plt.close(f)
    plt.close('all')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fcs-dir', help='raw FCS files (rebuild the per-well table from them)')
    ap.add_argument('--metadata', help='well metadata CSV (filename, type, staining, Genetics); with --fcs-dir')
    a = ap.parse_args()
    if a.fcs_dir:
        table = channel_statistics_from_fcs(a.fcs_dir, a.metadata)
        table.to_csv(out('fig5G_channel_statistics.csv'), index=False)
    else:
        table = pd.read_csv(paths.inp(*DATA))

    wells = normalize_halotag(normalize_biotin(load_wells(table)))
    checks = normalization_checks(wells)
    summary, tests = statistics(wells)
    checks.to_csv(out('fig5G_normalization_checks.csv'), index=False)
    summary.to_csv(out('fig5G_summary_n.csv'), index=False)
    tests.to_csv(out('fig5G_welch_tests.csv'), index=False)
    wells.to_csv(out('fig5G_per_well.csv'), index=False)
    figure(wells, tests)
    print(tests[['panel', 'comparison', 'method', 'estimate', 'p']].round(4).to_string(index=False))
    print('wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
