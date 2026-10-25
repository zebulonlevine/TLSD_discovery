"""Fig 7A, S9B–D: old and new ribosomes after UBA5, TLSD1 (R3HDM4) or double knockdown (RiboHalo).

RPL29-HaloTag (60S) and RPS3-HaloTag (40S) RPE1 lines, sgNT, sgUBA5, sgTLSD1 or sgUBA5 + sgTLSD1,
4 wells each. Per-well FlowJo medians are normalized to the mean of the sgNT wells of the same line
(strip_plots.py). The 3xHA-TLSD1 rescue lines (samples *_240, *_241) were collected in the same plate
but are not part of these panels.
Fig 7A: old 60S; Fig S9B: old 40S; Fig S9C: new 60S; Fig S9D: new 40S.

Input (data root, flow_cytometry/7A_S9B-D/): R3HDM4_L29_knockdown_cover_3xHA.xlsx (FlowJo table export)
Outputs (output root, flow_cytometry/):
  fig7A_old_L29.pdf, figS9B_old_S3.pdf, figS9C_new_L29.pdf, figS9D_new_S3.pdf
  fig7A_figS9BD_stats_groups.csv     n, mean, median, SD, SEM per group
  fig7A_figS9BD_stats_pairwise.csv   Welch's t-test for every pair of groups, BH within each panel

    python flow_cytometry/fig7A_figS9BD_uba5_tlsd1.py
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

XLSX = 'R3HDM4_L29_knockdown_cover_3xHA.xlsx'
PANELS = {('old', 'L29'): 'fig7A', ('old', 'S3'): 'figS9B', ('new', 'L29'): 'figS9C', ('new', 'S3'): 'figS9D'}
CONTROL = 'NT2'                         # sgNT wells
RESCUE = ('240', '241')                 # 3xHA-TLSD1 rescue lines, not plotted


def out(name):
    return paths.outp('flow_cytometry', name, mkdir=True)


def panel_data(data, m, line):
    """Wells of one line (rescue lines dropped), `m` median normalized to the sgNT mean."""
    sub = data.loc[data.SampleID.str.contains(line)
                   & ~data.SampleID.str.contains('|'.join(RESCUE)), :].copy()
    spl.normalize_to_control(sub, m + '_median', m + '_norm', key=pd.Series(line, index=sub.index),
                             is_control=sub.SampleID.str.contains(CONTROL))
    return sub


def main():
    data = pd.read_excel(paths.inp('flow_cytometry', '7A_S9B-D', XLSX), sheet_name='Sheet1')
    summaries, tests = [], []
    for (m, line), panel in PANELS.items():
        sub = panel_data(data, m, line)
        spl.strip_with_mean(sub.sort_values('SampleID'), m + '_norm', out(f'{panel}_{m}_{line}.pdf'),
                            title=f'{m} ({line})', xlabel='Sample', ylabel=f'{m} Median')
        summaries.append(spl.group_summary(sub, m + '_norm', panel=f'{panel} {m} ({line})'))
        tests.append(spl.pairwise_welch_bh(sub, m + '_norm', panel=f'{panel} {m} ({line})'))
    pd.concat(summaries, ignore_index=True).to_csv(out('fig7A_figS9BD_stats_groups.csv'), index=False)
    pd.concat(tests, ignore_index=True).to_csv(out('fig7A_figS9BD_stats_pairwise.csv'), index=False)
    print('wrote', paths.outp('flow_cytometry'))


if __name__ == '__main__':
    main()
