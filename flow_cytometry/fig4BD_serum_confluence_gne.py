"""Fig 4B, 4D: old 60S after UBA5 knockdown by serum and confluence, and with GNE-7883 (RiboHalo, 5 days).

RPL29-HaloTag RPE1 cells, sgNT or sgUBA5, 3 wells each, chased for 5 days at confluent or
non-confluent density in 1% or 10% FBS, with or without GNE-7883 (added at the start of the chase).
Per-well FlowJo medians of old (first-label) 60S are normalized to the mean of the sgNT wells of the
same condition (strip_plots.py). One plot holds both panels: Fig 4B is the serum x confluence
groups, Fig 4D the non-confluent 1% FBS groups with and without GNE-7883.

Input (data root, flow_cytometry/4B_4D/): Confluence_GNE_5d.xlsx (FlowJo table export)
Outputs (output root, flow_cytometry/):
  fig4BD_old_L29.pdf
  fig4BD_stats_groups.csv     n, mean, median, SD, SEM per group
  fig4BD_stats_pairwise.csv   sgUBA5 vs sgNT per condition, Welch's t-test, BH across the 4 conditions

    python flow_cytometry/fig4BD_serum_confluence_gne.py
"""
import os
import sys

import pandas as pd
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402
import strip_plots as spl  # noqa: E402

XLSX = 'Confluence_GNE_5d.xlsx'
CONTROL, KNOCKDOWN = 'NT2', 'NT_UBA5'           # sgNT, sgUBA5 suffixes of the sample IDs
CONDITIONS = ['L29_confluent_1FBS', 'L29_confluent_10FBS', 'L29_nonconfluent_1FBS',   # Fig 4B
              'L29_nonconfluent_1FBS_GNE']                                           # Fig 4D (+ its control)
ORDER = [f'{c}_{g}' for c in CONDITIONS for g in (CONTROL, KNOCKDOWN)]


def out(name):
    return paths.outp('flow_cytometry', name, mkdir=True)


def load():
    """Plotted wells, old median normalized to the sgNT mean of the same condition."""
    df = pd.read_excel(paths.inp('flow_cytometry', '4B_4D', XLSX), sheet_name='Sheet1')
    condition = df.SampleID.str.replace(f'_({CONTROL}|{KNOCKDOWN})$', '', regex=True)
    spl.normalize_to_control(df, 'old_median', 'old_norm', key=condition,
                             is_control=df.SampleID.str.endswith('_' + CONTROL))
    sub = df.loc[df.SampleID.isin(ORDER)].copy()
    sub['SampleID'] = pd.Categorical(sub['SampleID'], categories=ORDER, ordered=True)
    return sub


def main():
    sub = load()
    spl.strip_with_mean(sub, 'old_norm', out('fig4BD_old_L29.pdf'))

    spl.group_summary(sub, 'old_norm', panel='fig4BD old (L29)').to_csv(
        out('fig4BD_stats_groups.csv'), index=False)
    pairs = [(f'{c}_{CONTROL}', f'{c}_{KNOCKDOWN}') for c in CONDITIONS]
    spl.pairwise_welch_bh(sub, 'old_norm', pairs=pairs, panel='fig4BD old (L29)').to_csv(
        out('fig4BD_stats_pairwise.csv'), index=False)
    print('wrote', paths.outp('flow_cytometry'))


if __name__ == '__main__':
    main()
