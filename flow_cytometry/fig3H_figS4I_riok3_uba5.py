"""Fig 3H, S4I: old and new ribosomes after RIOK3, UBA5 or double knockdown (RiboHalo, 5-day chase).

RPS3-HaloTag (40S) and RPL29-HaloTag (60S) RPE1 lines, sgNT, sgUBA5, sgRIOK3 or sgUBA5 + sgRIOK3,
4 wells each; two RIOK3 sgRNAs, each run as its own experiment (sheets 'guide1', 'guide2').
Per-well FlowJo medians are normalized to the mean of the sgNT wells of the same line (strip_plots.py).
Fig 3H: old ribosomes (40S left, 60S right); Fig S4I: new ribosomes.

Input (data root, flow_cytometry/3H_S4I/): RIOK3.xlsx (FlowJo table export, one sheet per guide)
Outputs (output root, flow_cytometry/):
  fig3H_old_{S3,L29}_{guide1,guide2}.pdf, figS4I_new_{S3,L29}_{guide1,guide2}.pdf
  fig3H_figS4I_stats_groups.csv     n, mean, median, SD, SEM per group
  fig3H_figS4I_stats_pairwise.csv   Welch's t-test for every pair of groups, BH within each panel

    python flow_cytometry/fig3H_figS4I_riok3_uba5.py
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

XLSX = 'RIOK3.xlsx'
GUIDES = ['guide1', 'guide2']           # sheets; RIOK3 sgRNA 2 and sgRNA 3
LINES = ['L29', 'S3']                   # RPL29-HaloTag (60S), RPS3-HaloTag (40S)
PANELS = {'old': 'fig3H', 'new': 'figS4I'}
COLUMNS = ['fileName', 'SampleID', 'Count', 'FreqCells', 'FreqSinglets', 'FreqPos', 'oldpernew_median',
           'oldpernew_CV', 'new_median', 'new_CV', 'old_median', 'old_CV']   # new = RL1-A, old = YL1-A


def out(name):
    return paths.outp('flow_cytometry', name, mkdir=True)


def load(guide):
    df = pd.read_excel(paths.inp('flow_cytometry', '3H_S4I', XLSX), sheet_name=guide)
    df.columns = COLUMNS
    line = df.SampleID.str.split('_').str[0]
    for m in PANELS:
        spl.normalize_to_control(df, m + '_median', m + '_norm', key=line,
                                 is_control=df.SampleID.str.contains('_NT$'))
    return df


def main():
    summaries, tests = [], []
    for guide in GUIDES:
        df = load(guide)
        for m, panel in PANELS.items():
            for line in LINES:
                d = df.loc[df.SampleID.str.contains(line), :]
                spl.strip_with_mean(d.sort_values('SampleID'), m + '_norm',
                                    out(f'{panel}_{m}_{line}_{guide}.pdf'),
                                    title=f'{m} ({line})', xlabel='Sample', ylabel=f'{m} Median')
                name = f'{panel} {m} ({line}) {guide}'
                summaries.append(spl.group_summary(d, m + '_norm', panel=name))
                tests.append(spl.pairwise_welch_bh(d, m + '_norm', panel=name))
    pd.concat(summaries, ignore_index=True).to_csv(out('fig3H_figS4I_stats_groups.csv'), index=False)
    pd.concat(tests, ignore_index=True).to_csv(out('fig3H_figS4I_stats_pairwise.csv'), index=False)
    print('wrote', paths.outp('flow_cytometry'))


if __name__ == '__main__':
    main()
