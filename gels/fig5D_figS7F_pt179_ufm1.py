"""Fig 5D (right), S7F: loss of RPL26 UFMylation after PT-179 degrades SD40-UBA5.

60S (RPL29-HaloTag) and 40S (HaloTag-RPS24) reporter lines expressing SD40-UBA5 over sgUBA5 were
treated with 20 uM PT-179 for 0-24 h, with a 10 min anisomycin pulse before lysis to induce RPL26
UFMylation, and immunoblotted for UFM1 and GAPDH (n = 3 per time point). Input: per-lane UFM1
signal from ImageStudio, already normalized in the spreadsheet (relative to cells with endogenous
UBA5; the normalization step itself is not part of this repository); columns `norm_FullLane_1`
(60S reporter line) and `norm_FullLane_4` (40S reporter line). The sheet codes the 24 h sample as
time 3.5, its position on the broken axis.

  1. Mean of the 3 replicates per line and time point.
  2. Plot replicates and the mean line on an axis broken between 2 h and 24 h.

Fig 5D plots the 40S reporter line; Fig S7F plots both lines. No statistics.

Outputs (output root, gels/):
  fig5D_pt179_ufm1.pdf, figS7F_pt179_ufm1.pdf
  fig5D_figS7F_values.csv   per-replicate values (time in hours and axis position)
  fig5D_figS7F_means.csv    mean per line and time point

    python gels/fig5D_figS7F_pt179_ufm1.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
if 'ipykernel' not in sys.modules:
    mpl.use('Agg')
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402

DATA_FILE = ('gels', '5D_S7F', 'blotData_fullLane_UFM1.xlsx')
OUT_DIR = 'gels'
LINES = {'norm_FullLane_1': '60S reporter', 'norm_FullLane_4': '40S reporter'}
COLORS = {'60S reporter': '#000000', '40S reporter': '#A7A9AC'}
AXIS_POS_24H = 3.5           # the sheet's time value for the 24 h sample
BREAK_AT = 2.95              # x position of the axis-break marks
STYLE = {
    'figure.dpi': 130, 'figure.facecolor': 'white', 'font.size': 9,
    'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'pdf.fonttype': 42,
}


def out(name):
    return paths.outp(OUT_DIR, name, mkdir=True)


def load():
    df = pd.read_excel(paths.inp(*DATA_FILE))
    long = df.melt(id_vars='time', value_vars=list(LINES), var_name='column', value_name='ufm1')
    long['line'] = long['column'].map(LINES)
    long = long.rename(columns={'time': 'axis_pos'})
    long['time_h'] = long['axis_pos'].where(long['axis_pos'] != AXIS_POS_24H, 24.0)
    return long[['line', 'time_h', 'axis_pos', 'ufm1']]


def means(long):
    return long.groupby(['line', 'time_h', 'axis_pos'], sort=False)['ufm1'].agg(['mean', 'count']).reset_index()


def _axis_break(ax):
    kw = dict(color='black', lw=1.0, clip_on=False, transform=ax.get_xaxis_transform())
    for dx in (-0.06, 0.06):
        ax.plot([BREAK_AT + dx - 0.05, BREAK_AT + dx + 0.05], [-0.025, 0.025], **kw)


def plot(long, m, lines, path, ylabel, figsize):
    with mpl.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=figsize)
        for line in lines:
            v = long[long['line'] == line]
            mm = m[m['line'] == line]
            ax.scatter(v['axis_pos'], v['ufm1'], s=14, color=COLORS[line], lw=0, zorder=2,
                       label=f'{line} UFM1' if len(lines) > 1 else 'replicate')
            ax.plot(mm['axis_pos'], mm['mean'], color=COLORS[line], lw=1.5, zorder=1,
                    label='(mean)' if len(lines) > 1 else 'mean')
        ax.set_xticks([0, 0.5, 1, 1.5, 2, AXIS_POS_24H])
        ax.set_xticklabels(['0', '0.5', '1', '1.5', '2', '24'])
        ax.set_xlim(-0.2, AXIS_POS_24H + 0.2)
        ax.set_ylim(0, 0.4)
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
        _axis_break(ax)
        ax.set_xlabel('Time Since PT-179 addition (h)')
        ax.set_ylabel(ylabel)
        ax.legend(loc='upper right', fontsize=7)
        fig.savefig(path, dpi=300, bbox_inches='tight')
        plt.close(fig)


def main():
    long = load()
    m = means(long)
    long.to_csv(out('fig5D_figS7F_values.csv'), index=False)
    m.to_csv(out('fig5D_figS7F_means.csv'), index=False)
    plot(long, m, ['40S reporter'], out('fig5D_pt179_ufm1.pdf'),
         'Fraction endogenous\nUFMylation', (2.0, 1.9))
    plot(long, m, ['60S reporter', '40S reporter'], out('figS7F_pt179_ufm1.pdf'),
         'Full-lane UFMylation\n(fraction of cells without\nsgUBA5 + SD40-UBA5)', (2.3, 1.9))
    print(m.pivot(index='time_h', columns='line', values='mean').round(4).to_string())
    print('wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
