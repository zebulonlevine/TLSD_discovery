"""Fig 1C: RiboHalo turnover curves in quiescent RPE1 cells (exponential fit, t1/2).

Quiescent (confluent, 1% serum, 20 uM STLC) RPE1 cells with HaloTag knocked into RPL29 (60S, clone
L29A2) or RPS3 (40S, clone S3F10) were pulse-labelled with JF549-Halo or JF646-Halo, chased, and at
each time point (days 1-15) counter-stained with the other dye. Singly-labelled wells (same dye at
start and end) give the full signal for each dye. Gating (single viable HaloTag+ cells) and per-well
median fluorescence were done in FlowJo; this script starts from the FlowJo exports
(JF646 = RL1-A / Comp-APC-A, JF549 = YL1-A / Comp-PE-A).

  1. Each pulse-chase well's median is divided by the mean of the same-day singly-labelled
     controls for that dye (6 wells per dye); the control SEM is propagated (uncertainties).
  2. Fraction old per well = mean of the two estimates: pulse dye / control, and
     1 - chase dye / control.
  3. Per day: mean of the 6 pulse-chase wells (3 per dye order); error = root-sum-square of the
     replicate SEM and the mean propagated control error.
  4. A * exp(-k t) fitted to the daily means, weighted by that error (absolute sigma);
     t1/2 = ln 2 / k, shown with its 95% confidence interval (1.96 SE), and a 95% band on the curve.

The per-cell new/old histograms in Fig 1C were drawn in FlowJo (no code).

Inputs (data root, flow_cytometry/1C/): day{1,3,5,7,9,11,13,15}.csv (FlowJo exports; day 5 was
acquired in tubes and is annotated from the file names), PlateMap_Feb.xlsx (plate layout).

Outputs (output root, flow_cytometry/):
  fig1C_RPL29_turnover_fit.pdf, fig1C_RPS3_turnover_fit.pdf   large (left) and small (right) subunit
  fig1C_turnover_fits.csv      A, k, t1/2 (SE and 95% CI) per line
  fig1C_turnover_points.csv    daily means and errors plotted
  fig1C_turnover_wells.csv     per-well normalized values

    python flow_cytometry/fig1C_ribohalo_turnover.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import paths  # noqa: E402

DATA = ('flow_cytometry', '1C')
DAYS = [1, 3, 5, 7, 9, 11, 13, 15]       # days 17 and 19 were collected but most cells were dead
PLATE_MAP = 'PlateMap_Feb.xlsx'
LINES = {'L29A2': 'RPL29', 'S3F10': 'RPS3'}   # clone -> tagged protein (60S, 40S)
PULSE_CHASE = ['549first', '646first']
FLOWJO_COLUMNS = {                         # FlowJo statistic -> column
    'Count': 'viable_count',
    'Freq. of Parent (%)': 'viable_percent',
    'Median (RL1-A)': 'median_RL1A', 'Median (Comp-APC-A :: null)': 'median_RL1A',    # JF646
    'Median (YL1-A)': 'median_YL1A', 'Median (Comp-PE-A :: null)': 'median_YL1A',     # JF549
}
TUBE_TREATMENTS = {'549_646': '549first', '646_549': '646first', '549ctrl': '549ctrl', '646ctrl': '646ctrl'}


def out(name):
    return paths.outp('flow_cytometry', name, mkdir=True)


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

def load_day(day, plate_map):
    """One day's FlowJo export, annotated with cell line and treatment."""
    t = pd.read_csv(paths.inp(*DATA, f'day{day}.csv'))
    t.columns = ['FCS'] + [FLOWJO_COLUMNS[c.split(' | ')[-1]] for c in t.columns[1:]]
    if day == 5:                           # acquired in tubes: annotation is in the file name
        name = t.FCS.str.slice(None, -4)
        t['Cell Line'] = name.str.slice(0, 5).str.strip('_')
        t['Treatment'] = ''
        for token, treatment in TUBE_TREATMENTS.items():
            t.loc[name.str.contains(token), 'Treatment'] = treatment
    else:                                  # 'A1.fcs' or 'Specimen_001_C1_C01_001.fcs' -> well
        t['Well'] = t.FCS.str.extract(r'(?:^|_)([A-H]\d{1,2})(?:_|\.fcs$)', expand=False)
        t = pd.merge(t, plate_map.loc[plate_map.Day == day, ['Well', 'Cell Line', 'Treatment']], on='Well', how='left')
    t['Day'] = day
    return t


def load_wells():
    plate_map = pd.read_excel(paths.inp(*DATA, PLATE_MAP))
    plate_map['Well'] = plate_map.Row.str.cat(plate_map.Column.astype(str))
    return pd.concat([load_day(d, plate_map) for d in DAYS], ignore_index=True)


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------

def control(grouped, day, treatment, col):
    """Mean +/- SEM of the same-day singly-labelled control wells for one dye."""
    row = grouped.loc[(grouped.Treatment == treatment) & (grouped.Day == day)]
    return unc.ufloat(row[(col, 'mean')].values[0], row[(col, 'std')].values[0] / np.sqrt(row[(col, 'count')].values[0]))


def fraction_old(wells):
    """Per-well fraction of old (pulse-labelled) protein, with the propagated control error."""
    wells = wells.copy()
    grouped = wells[['Day', 'Treatment', 'median_YL1A', 'median_RL1A']].groupby(
        ['Day', 'Treatment']).agg(['mean', 'std', 'count']).reset_index()
    for col in ['FractionMax_RL1A', 'FractionMax_RL1A_err', 'FractionMax_YL1A', 'FractionMax_YL1A_err',
                'FractionRemain', 'FractionRemain_err']:
        wells[col] = np.nan
    for day in wells.Day.unique():
        ctrl_rl = control(grouped, day, '646ctrl', 'median_RL1A')
        ctrl_yl = control(grouped, day, '549ctrl', 'median_YL1A')
        for idx in wells.index[wells.Day == day]:
            norm_rl = wells.loc[idx, 'median_RL1A'] / ctrl_rl
            norm_yl = wells.loc[idx, 'median_YL1A'] / ctrl_yl
            wells.loc[idx, ['FractionMax_RL1A', 'FractionMax_RL1A_err']] = norm_rl.nominal_value, norm_rl.std_dev
            wells.loc[idx, ['FractionMax_YL1A', 'FractionMax_YL1A_err']] = norm_yl.nominal_value, norm_yl.std_dev
            treatment = wells.loc[idx, 'Treatment']
            if treatment == '549first':          # JF549 = old, JF646 = new
                old_yl, old_rl = norm_yl, 1 - norm_rl
            elif treatment == '646first':        # JF646 = old, JF549 = new
                old_yl, old_rl = 1 - norm_yl, norm_rl
            else:
                continue
            combined = (old_rl + old_yl) / 2
            wells.loc[idx, 'FractionRemain'] = combined.nominal_value
            wells.loc[idx, 'FractionRemain_err'] = combined.std_dev
    return wells


def daily_points(wells):
    """Mean of the pulse-chase wells per day; error = replicate SEM (+) mean propagated error."""
    t = wells.loc[wells.Treatment.isin(PULSE_CHASE), ['Day', 'FractionRemain', 'FractionRemain_err']].groupby('Day').agg(
        mean=('FractionRemain', 'mean'), n=('FractionRemain', 'count'), std=('FractionRemain', 'std'),
        mean_prop_err=('FractionRemain_err', 'mean'))
    t['rep_sem'] = t['std'] / np.sqrt(t['n'])
    t['combined_sem'] = np.sqrt(t.rep_sem ** 2 + t.mean_prop_err ** 2)
    return t


def exp_decay(x, A, k):
    """Exponential decay from A with rate k."""
    return A * np.exp(-k * x)


def fit_and_plot(points, clone, path):
    """Weighted exponential fit; curve with 95% band, daily means with error bars."""
    popt, pcov = curve_fit(exp_decay, points.index, points['mean'].values,
                           sigma=points.combined_sem.values, absolute_sigma=True)
    A, k = unc.correlated_values(popt, pcov)
    px = np.linspace(0, 16, 100)
    py = A * unp.exp(-k * px)
    nom, std = unp.nominal_values(py), unp.std_devs(py)
    half_life = unp.log(2) / k

    plt.figure()
    plt.plot(px, nom, c='darkblue')
    plt.fill_between(px, nom - 1.96 * std, nom + 1.96 * std, color='lightblue', alpha=0.5)
    plt.errorbar(x=points.index, y=points['mean'], color='k', xerr=None, yerr=points['combined_sem'],
                 fmt='o', capsize=3)
    plt.title(f'{clone} exponential decay')
    plt.xlabel('Days')
    plt.ylabel('Fraction ribosomes old')
    plt.ylim(0, nom.max() + 0.05)
    plt.savefig(path)
    plt.close('all')
    return {'A': A.nominal_value, 'A_se': A.std_dev, 'k_per_day': k.nominal_value, 'k_se': k.std_dev,
            't_half_days': half_life.nominal_value, 't_half_se': half_life.std_dev,
            't_half_ci95_low': half_life.nominal_value - 1.96 * half_life.std_dev,
            't_half_ci95_high': half_life.nominal_value + 1.96 * half_life.std_dev}


def main():
    wells = load_wells()
    fits, points_all, wells_all = [], [], []
    for clone, protein in LINES.items():
        w = fraction_old(wells.loc[(wells['Cell Line'] == clone) & (wells.Day < 16)].reset_index(drop=True))
        points = daily_points(w)
        fit = fit_and_plot(points, clone, out(f'fig1C_{protein}_turnover_fit.pdf'))
        fits.append({'clone': clone, 'protein': protein, 'n_days': len(points), **fit})
        points_all.append(points.reset_index().assign(clone=clone, protein=protein))
        wells_all.append(w.assign(protein=protein))
        print(f'{protein} ({clone}): t1/2 = {fit["t_half_days"]:.3f} +/- {1.96 * fit["t_half_se"]:.3f} days (95% CI)')
    pd.DataFrame(fits).to_csv(out('fig1C_turnover_fits.csv'), index=False)
    pd.concat(points_all, ignore_index=True).to_csv(out('fig1C_turnover_points.csv'), index=False)
    pd.concat(wells_all, ignore_index=True).to_csv(out('fig1C_turnover_wells.csv'), index=False)
    print('wrote', paths.outp('flow_cytometry'))


if __name__ == '__main__':
    main()
