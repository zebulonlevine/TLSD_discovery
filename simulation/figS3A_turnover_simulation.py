"""Fig S3A: increased synthesis and increased degradation raise new/old alike.

Two perturbations from steady state (t1/2 = 7 d, synthesis 100 molecules/day):

  synthesis    k_syn x 3
  degradation  k_deg x f, with f found (brentq) so that new/old equals the
               synthesis scenario's new/old 14 days after the change (f ~ 1.46)

Molecules are "new" if no older than 7 days (the initial half-life) and "old"
otherwise. Each run settles for 120 days at the initial rates first; time is
reported relative to the rate change.

Writes, under the output root (common/paths.py), in simulation/:
  figS3A_ratio_and_abundance.pdf  the two panels of Fig S3A (new/old ratio; log2 new and old)
  figS3A_curves.csv               source data for both scenarios
  increased_ratio_summary.pdf     diagnostic figure with every quantity for both scenarios

    python figS3A_turnover_simulation.py        # ~1 min
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
if 'ipykernel' not in sys.modules:            # headless when run from the CLI
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402
from biomolecule_turnover import (  # noqa: E402
    SimulationParams, simulate_turnover, calculate_steady_state,
)

# baseline
T_HALF = 7.0                       # initial half-life = the new/old age cutoff (days)
K_DEG_BASE = np.log(2) / T_HALF
K_SYN_BASE = 100.0
N_SS = K_SYN_BASE / K_DEG_BASE

# timing (days)
T_SETTLE = 120.0                   # hidden settling period at the initial rates
T_PRE = 5.0                        # visible window before the change
T_MATCH = 14.0                     # where the two scenarios' new/old ratios are matched
T_POST = 50.0                      # visible window after the change
DT = 0.05
T_CHANGE = T_SETTLE + T_PRE
T_MAX = T_SETTLE + T_PRE + T_POST

SYN_FACTOR = 3.0

# Fig S3A colours (degradation as the lysosome/autophagy magenta used throughout)
DEG_COLOR = '#C51B7D'
SYN_COLOR = '#4E9442'


def run_scenario(k_syn_final=None, k_deg_final=None):
    """Simulate a rate change at T_CHANGE; time is returned relative to the change."""
    params = SimulationParams(
        k_syn_initial=K_SYN_BASE, k_deg_initial=K_DEG_BASE, N0=N_SS,
        t_max=T_MAX, dt=DT, t_change=T_CHANGE,
        k_syn_final=k_syn_final, k_deg_final=k_deg_final)
    results = simulate_turnover(params)
    results.time = results.time - T_CHANGE
    return results, params


def value_at(results, series, t_rel):
    return series[np.argmin(np.abs(results.time - t_rel))]


def run_both():
    """The two scenarios of Fig S3A, plus the matched degradation factor."""
    syn, syn_p = run_scenario(k_syn_final=K_SYN_BASE * SYN_FACTOR)
    target = value_at(syn, syn.ratio_new_to_old, T_MATCH)

    def mismatch(k_deg_new):
        r, _ = run_scenario(k_deg_final=k_deg_new)
        return value_at(r, r.ratio_new_to_old, T_MATCH) - target

    k_deg_matched = brentq(mismatch, K_DEG_BASE * 1.01, K_DEG_BASE * 5, xtol=1e-8)
    deg, deg_p = run_scenario(k_deg_final=k_deg_matched)
    return dict(syn=syn, syn_p=syn_p, deg=deg, deg_p=deg_p, target=target,
                deg_factor=k_deg_matched / K_DEG_BASE)


def write_curves(sim, path):
    rows = []
    for name, label in (('deg', f'degradation x{sim["deg_factor"]:.3f}'),
                        ('syn', f'synthesis x{SYN_FACTOR:g}')):
        r = sim[name]
        m = (r.time >= -T_PRE) & (r.time <= T_POST)
        rows.append(pd.DataFrame({'scenario': label, 'day': r.time[m], 'new_old_ratio': r.ratio_new_to_old[m],
                                  'total': r.total[m], 'new': r.new[m], 'old': r.old[m]}))
    pd.concat(rows).to_csv(path, index=False, float_format='%.6g')


def plot_figS3A(sim, path):
    """The two panels as laid out in Fig S3A."""
    fig, (ax_r, ax_a) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    deg_lab = f'↑ degradation ({sim["deg_factor"]:.2f}x)'
    syn_lab = f'↑ synthesis ({SYN_FACTOR:g}x)'
    for name, color, label in (('deg', DEG_COLOR, deg_lab), ('syn', SYN_COLOR, syn_lab)):
        r = sim[name]
        m = (r.time >= -T_PRE) & (r.time <= T_POST)
        ax_r.plot(r.time[m], r.ratio_new_to_old[m], color=color, lw=2, label=label)
        ax_a.plot(r.time[m], np.log2(r.new[m]), color=color, lw=2, ls='--')
        ax_a.plot(r.time[m], np.log2(r.old[m]), color=color, lw=2)
        for series in (r.new, r.old):
            ax_a.plot(T_MATCH, np.log2(value_at(r, series, T_MATCH)), 'o', color=color, ms=6, zorder=5)
    ax_r.plot(T_MATCH, sim['target'], 'ko', ms=6, zorder=5)
    ax_r.set_ylabel('New/Old Ratio')
    ax_r.legend(frameon=False, fontsize=8, loc='upper right')
    ax_a.set_ylabel('log$_2$ (Abundance)')
    ax_a.plot([], [], 'k--', lw=2, label=f'new molecules\n(< {T_HALF:g} days)')
    ax_a.plot([], [], 'k-', lw=2, label=f'old molecules\n(> {T_HALF:g} days)')
    ax_a.legend(frameon=False, fontsize=8, loc='center right')
    for ax in (ax_r, ax_a):
        ax.axvline(0, color='0.8', ls='--', lw=1.5)
        ax.set_xlim(-T_PRE, T_POST)
        ax.set_xlabel('Days')
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Simulated Data')
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def plot_summary(sim, path):
    """Diagnostic figure: total, old and new for each scenario, then the overlays."""
    results_inc_deg, params_inc_deg = sim['deg'], sim['deg_p']
    results_inc_syn, params_inc_syn = sim['syn'], sim['syn_p']
    k_deg_increase_factor, ratio_target = sim['deg_factor'], sim['target']
    x_lo, x_hi, t_match = -T_PRE, T_POST, T_MATCH
    colors = {'deg': '#F15A29', 'syn': '#1C75BC', 'total': '#1f77b4', 'old': '#d62728', 'new': '#2ca02c'}

    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.35, height_ratios=[1, 1, 1.1])

    # rows 0-1: per-scenario total, old, new
    for results, params, title, col in (
            (results_inc_deg, params_inc_deg, f'↑ Degradation (k_deg ×{k_deg_increase_factor:.2f})', 0),
            (results_inc_syn, params_inc_syn, f'↑ Synthesis (k_syn ×{SYN_FACTOR:.1f})', 2)):
        ax_total = fig.add_subplot(gs[0, col:col + 2])
        ax_old = fig.add_subplot(gs[1, col])
        ax_new = fig.add_subplot(gs[1, col + 1])
        mask = (results.time >= x_lo) & (results.time <= x_hi)
        t = results.time[mask]

        ax_total.plot(t, results.total[mask], color=colors['total'], lw=2)
        ax_total.axhline(N_SS, color='black', ls=':', alpha=0.4, label=f'Initial SS = {N_SS:.0f}')
        k_s = params.k_syn_final if params.k_syn_final is not None else params.k_syn_initial
        k_d = params.k_deg_final if params.k_deg_final is not None else params.k_deg_initial
        ss_new = calculate_steady_state(k_s, k_d)
        ax_total.axhline(ss_new, color='orange', ls=':', alpha=0.7, label=f'New SS = {ss_new:.0f}')
        ax_total.set_ylabel('Total Abundance', fontsize=11)
        ax_total.set_title(f'{title}\nTotal Abundance', fontsize=12, fontweight='bold')
        ax_total.legend(loc='best', fontsize=9)

        ax_old.plot(t, results.old[mask], color=colors['old'], lw=2)
        ax_old.set_ylabel('Old Molecules', fontsize=11)
        ax_old.set_title(f'Old (age > {T_HALF:.0f} d)', fontsize=11, fontweight='bold')
        ax_old.set_xlabel('Days relative to perturbation', fontsize=10)

        ax_new.plot(t, results.new[mask], color=colors['new'], lw=2)
        ax_new.set_ylabel('New Molecules', fontsize=11)
        ax_new.set_title(f'New (age ≤ {T_HALF:.0f} d)', fontsize=11, fontweight='bold')
        ax_new.set_xlabel('Days relative to perturbation', fontsize=10)

        for ax in (ax_total, ax_old, ax_new):
            ax.axvline(0, color='gray', ls='--', alpha=0.5, lw=1.5)
            ax.axvline(t_match, color='red', ls='-', alpha=0.4, lw=2)
            ax.set_xlim(x_lo, x_hi)
            ax.grid(True, alpha=0.25)

    # row 2: overlays
    ax_ratio = fig.add_subplot(gs[2, 0])
    ax_abund = fig.add_subplot(gs[2, 1])
    ax_old_ov = fig.add_subplot(gs[2, 2])
    ax_new_ov = fig.add_subplot(gs[2, 3])
    mask_deg = (results_inc_deg.time >= x_lo) & (results_inc_deg.time <= x_hi)
    mask_syn = (results_inc_syn.time >= x_lo) & (results_inc_syn.time <= x_hi)
    t_deg, t_syn = results_inc_deg.time[mask_deg], results_inc_syn.time[mask_syn]
    deg_label = f'↑ Deg (k_deg ×{k_deg_increase_factor:.2f})'
    syn_label = f'↑ Syn (k_syn ×{SYN_FACTOR:.1f})'
    idx_deg_m = np.argmin(np.abs(results_inc_deg.time - t_match))
    idx_syn_m = np.argmin(np.abs(results_inc_syn.time - t_match))

    def fmt(ax):
        ax.axvline(0, color='gray', ls='--', alpha=0.5, lw=1.5)
        ax.axvline(t_match, color='red', ls='-', alpha=0.4, lw=2)
        ax.set_xlim(x_lo, x_hi)
        ax.grid(True, alpha=0.25)
        ax.set_xlabel('Days', fontsize=10)

    ax_ratio.plot(t_deg, results_inc_deg.ratio_new_to_old[mask_deg], color=colors['deg'], lw=2.5, ls='-', label=deg_label)
    ax_ratio.plot(t_syn, results_inc_syn.ratio_new_to_old[mask_syn], color=colors['syn'], lw=2.5, ls='--', label=syn_label)
    ax_ratio.plot(t_match, ratio_target, 'ko', ms=9, zorder=5)
    ax_ratio.axhline(1.0, color='black', ls=':', alpha=0.4)
    ax_ratio.set_ylabel('New / Old Ratio', fontsize=11)
    ax_ratio.set_title('New / Old Ratio', fontsize=12, fontweight='bold')
    ax_ratio.legend(loc='upper right', fontsize=8, framealpha=0.9)
    fmt(ax_ratio)

    ax_abund.plot(t_deg, np.log2(results_inc_deg.total[mask_deg]), color=colors['deg'], lw=2.5, ls='-', label=deg_label)
    ax_abund.plot(t_syn, np.log2(results_inc_syn.total[mask_syn]), color=colors['syn'], lw=2.5, ls='--', label=syn_label)
    ax_abund.axhline(np.log2(N_SS), color='black', ls=':', alpha=0.4, label='Initial SS')
    ax_abund.plot(t_match, np.log2(results_inc_deg.total[idx_deg_m]), 'o', color=colors['deg'], ms=9, zorder=5)
    ax_abund.plot(t_match, np.log2(results_inc_syn.total[idx_syn_m]), 'o', color=colors['syn'], ms=9, zorder=5)
    ax_abund.set_ylabel('log₂(Total Abundance)', fontsize=11)
    ax_abund.set_title('Total Abundance (log₂)', fontsize=12, fontweight='bold')
    ax_abund.legend(loc='center right', fontsize=8, framealpha=0.9)
    fmt(ax_abund)

    for ax, attr, name in ((ax_old_ov, 'old', f'Old (age > {T_HALF:.0f} d) — log₂'),
                           (ax_new_ov, 'new', f'New (age ≤ {T_HALF:.0f} d) — log₂')):
        ax.plot(t_deg, np.log2(getattr(results_inc_deg, attr)[mask_deg]), color=colors['deg'], lw=2.5, ls='-', label=deg_label)
        ax.plot(t_syn, np.log2(getattr(results_inc_syn, attr)[mask_syn]), color=colors['syn'], lw=2.5, ls='--', label=syn_label)
        ax.plot(t_match, np.log2(getattr(results_inc_deg, attr)[idx_deg_m]), 'o', color=colors['deg'], ms=9, zorder=5)
        ax.plot(t_match, np.log2(getattr(results_inc_syn, attr)[idx_syn_m]), 'o', color=colors['syn'], ms=9, zorder=5)
        ax.set_ylabel(f'log₂({attr.capitalize()} Molecules)', fontsize=11)
        ax.set_title(name, fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=8, framealpha=0.9)
        fmt(ax)

    fig.suptitle('Increased New/Old Ratio: Degradation vs Synthesis\n'
                 f'Age cutoff = {T_HALF:.0f} d (initial t½)  ·  Ratios cross at +{t_match:.0f} d  ·  '
                 f'{T_SETTLE:.0f} d settling (hidden)', fontsize=14, fontweight='bold', y=0.98)
    fig.savefig(path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    sim = run_both()
    deg_ss = calculate_steady_state(K_SYN_BASE, K_DEG_BASE * sim['deg_factor'])
    syn_ss = calculate_steady_state(K_SYN_BASE * SYN_FACTOR, K_DEG_BASE)
    print(f'half-life (age cutoff) = {T_HALF:.0f} d; steady state = {N_SS:.1f}')
    print(f'synthesis x{SYN_FACTOR:.2f}: new/old at +{T_MATCH:.0f} d = {sim["target"]:.3f}, new steady state = {syn_ss:.0f}')
    print(f'degradation x{sim["deg_factor"]:.3f}: new/old at +{T_MATCH:.0f} d = '
          f'{value_at(sim["deg"], sim["deg"].ratio_new_to_old, T_MATCH):.3f}, new steady state = {deg_ss:.0f}')

    write_curves(sim, paths.outp('simulation', 'figS3A_curves.csv', mkdir=True))
    plot_figS3A(sim, paths.outp('simulation', 'figS3A_ratio_and_abundance.pdf'))
    plot_summary(sim, paths.outp('simulation', 'increased_ratio_summary.pdf'))
    print('wrote', paths.outp('simulation'))


if __name__ == '__main__':
    main()
