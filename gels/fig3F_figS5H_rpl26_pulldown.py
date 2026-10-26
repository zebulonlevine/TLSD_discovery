"""Fig 3F, S5H: new/old ratio of RPL26-3xFLAG-tagged vs untagged ribosomes (wildtype vs 5K5R).

Quiescent RPL29-HaloTag cells carrying no tagged RPL26 (2 replicates), RPL26-3xFLAG (wildtype, 3) or
non-UFMylatable RPL26-3xFLAG(5K5R) (3) were pulse-labelled with JF549 (old ribosomes) and chased with
JF646 (new ribosomes) for 5 days. Ribosomes were pelleted through a sucrose cushion (RP) and
anti-FLAG pulled down into flowthrough (FT, untagged) and elute (tagged). Input: ImageStudio band
intensities (local background) of in-gel JF549 / JF646 fluorescence and of the anti-RPL26 blot
(tagged and untagged bands) for lysate, RP, FT and elute, one row per sample (sheet `data`).

  1. new/old = (JF646 / RPL26) / (JF549 / RPL26), each fluorescence column first scaled to its
     maximum. The RPL26 and maximum scalings cancel in the ratio; they are kept so the values equal
     the original notebook's exactly.
  2. Normalize to the ribosome pellet: divide by the median RP new/old of the same construct
     (the mean of the two RP values for the no-tag control), then log2.
  3. Elute (tagged) vs FT (untagged) within each construct: Welch's t-test on the log2 values,
     unadjusted (n = 3 per group).

Fig 3F plots FT and elute for wildtype and 5K5R (points and mean). Fig S5H adds the ribosome
pellet and the no-tag control (points and range; the original seaborn bootstrap 95% CI equals the
range for n <= 3). Replicates are spread by fixed offsets, not random jitter.

Outputs (output root, gels/):
  fig3F_rpl26_pulldown.pdf, figS5H_rpl26_pulldown.pdf
  fig3F_figS5H_values.csv    per-sample log2 new/old (RP-normalized), all fractions
  fig3F_figS5H_summary.csv   mean, min, max and n per construct and fraction
  fig3F_figS5H_welch.csv     elute vs FT per construct

    python gels/fig3F_figS5H_rpl26_pulldown.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
if 'ipykernel' not in sys.modules:
    mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402

DATA_FILE = ('gels', '3F_S5H', 'QuiescentPulldownFluorGels.xlsx')
OUT_DIR = 'gels'
FRACTIONS = ['Lysate', 'RP', 'FT', 'elute']
CONSTRUCTS = {0: 'no tagged\nRPL26', 101: 'wildtype', 193: '5K5R'}   # plasmid IDs in the sheet
COLORS = {'RP': '#000000', 'FT': '#BCBEC0', 'elute': '#00AEEF'}
LABELS = {'RP': 'all ribosomes (sucrose pellet)', 'FT': 'untagged (flowthrough)', 'elute': 'tagged (elute)'}
STYLE = {
    'figure.dpi': 130, 'figure.facecolor': 'white', 'font.size': 9,
    'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'pdf.fonttype': 42,
}


def out(name):
    return paths.outp(OUT_DIR, name, mkdir=True)


def load():
    data = pd.read_excel(paths.inp(*DATA_FILE), sheet_name='data')
    data['Construct'] = data['Construct'].fillna(0)                 # blank = no tagged RPL26
    long = []
    for frac in FRACTIONS:
        d = data[['Construct', 'Rep']].copy()
        d['fraction'] = frac
        d['549'] = data[f'{frac}_549']
        d['646'] = data[f'{frac}_646']
        d['L26_all'] = data[f'{frac}_L26_tagged'] + data[f'{frac}_L26_untagged']
        long.append(d)
    return pd.concat(long, ignore_index=True)


def new_per_old(d):
    d = d.copy()
    n549 = d['549'] / d['L26_all']
    n646 = d['646'] / d['L26_all']
    d['newperold'] = (n646 / n646.max()) / (n549 / n549.max())
    rp_median = d[d['fraction'] == 'RP'].groupby('Construct')['newperold'].median()
    d['newperold_RPnorm'] = d['newperold'] / d['Construct'].map(rp_median)
    d['log2_newperold_RPnorm'] = d['newperold_RPnorm'].apply(lambda x: np.log2(x) if x > 0 else np.nan)
    return d.dropna(subset=['log2_newperold_RPnorm'])


def welch(d):
    rows = []
    for c in (101, 193):
        ft = d.loc[(d['Construct'] == c) & (d['fraction'] == 'FT'), 'log2_newperold_RPnorm']
        el = d.loc[(d['Construct'] == c) & (d['fraction'] == 'elute'), 'log2_newperold_RPnorm']
        r = stats.ttest_ind(el, ft, equal_var=False)
        p = r.pvalue
        rows.append({'construct': CONSTRUCTS[c], 'n_elute': len(el), 'n_FT': len(ft),
                     'mean_elute': el.mean(), 'mean_FT': ft.mean(), 't': r.statistic, 'df': r.df,
                     'p_welch': p, 'sig': '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'})
    return pd.DataFrame(rows)


def summary(d):
    g = d.groupby(['Construct', 'fraction'])['log2_newperold_RPnorm']
    s = g.agg(['mean', 'min', 'max', 'count']).reset_index()
    s['construct'] = s['Construct'].map(CONSTRUCTS).str.replace('\n', ' ')
    return s


def _spread(n, width):
    """Fixed horizontal offsets so replicate points with equal values stay visible."""
    return np.linspace(-width, width, n) if n > 1 else np.zeros(1)


def _xpos(constructs, fractions, width):
    offs = np.linspace(-width, width, len(fractions)) if len(fractions) > 1 else [0.0]
    return {(c, f): i + o for i, c in enumerate(constructs) for f, o in zip(fractions, offs)}


def plot_3F(d, path):
    constructs, fractions = [101, 193], ['FT', 'elute']
    pos = _xpos(constructs, fractions, 0.17)
    with mpl.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(2.4, 1.9))
        for (c, f), x in pos.items():
            v = d.loc[(d['Construct'] == c) & (d['fraction'] == f), 'log2_newperold_RPnorm']
            ax.scatter(x + _spread(len(v), 0.06), v, s=22, color=COLORS[f], lw=0, zorder=2,
                       label={'FT': 'untagged', 'elute': '3xFLAG'}[f] if c == constructs[0] else None)
            ax.hlines(v.mean(), x - 0.11, x + 0.11, color=COLORS[f], lw=1.5, zorder=3)
        ax.set_xticks(range(len(constructs)))
        ax.set_xticklabels([CONSTRUCTS[c] for c in constructs])
        ax.set_xlim(-0.5, len(constructs) - 0.5)
        ax.set_ylim(-0.42, 0.62)
        ax.set_ylabel('Log$_2$ FC new/old ratio\n(ribosome subset/total ribosomes)')
        ax.set_xlabel('tagged RPL26')
        ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=7, handletextpad=0.2)
        fig.savefig(path, dpi=300, bbox_inches='tight')
        plt.close(fig)


def plot_S5H(d, path):
    constructs, fractions = [0, 101, 193], ['RP', 'FT', 'elute']
    pos = _xpos(constructs, fractions, 0.22)
    with mpl.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(2.8, 2.0))
        for (c, f), x in pos.items():
            v = d.loc[(d['Construct'] == c) & (d['fraction'] == f), 'log2_newperold_RPnorm']
            if v.empty:
                continue
            ax.vlines(x, v.min(), v.max(), color='#666666', lw=1.5, zorder=1)
            ax.scatter(x + _spread(len(v), 0.05), v, s=26, color=COLORS[f], lw=0, zorder=2,
                       label=LABELS[f] if c == 101 else None)
        ax.set_xticks(range(len(constructs)))
        ax.set_xticklabels([CONSTRUCTS[c] if c == 0 else f'{CONSTRUCTS[c]}\ntagged\nRPL26' for c in constructs])
        ax.set_xlim(-0.5, len(constructs) - 0.5)
        ax.set_ylim(-0.42, 0.62)
        ax.set_ylabel('Log$_2$ FC new/old ratio\n(ribosome subset/total ribosomes)')
        ax.legend(loc='lower left', bbox_to_anchor=(0.1, 1.0), fontsize=7, handletextpad=0.2)
        fig.savefig(path, dpi=300, bbox_inches='tight')
        plt.close(fig)


def main():
    d = new_per_old(load())
    d.assign(construct=d['Construct'].map(CONSTRUCTS).str.replace('\n', ' ')).to_csv(
        out('fig3F_figS5H_values.csv'), index=False)
    summary(d).to_csv(out('fig3F_figS5H_summary.csv'), index=False)
    tests = welch(d)
    tests.to_csv(out('fig3F_figS5H_welch.csv'), index=False)
    plot_3F(d, out('fig3F_rpl26_pulldown.pdf'))
    plot_S5H(d, out('figS5H_rpl26_pulldown.pdf'))
    print(tests.to_string(index=False))
    print('wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
