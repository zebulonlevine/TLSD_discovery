#!/usr/bin/env python
"""
replot_refit_panels.py
========================

Redraw the quantification panels from the VBEM-refit statistics.

Which panels this touches, and which it deliberately does not:

  NEEDS redrawing -- error bars and significance stars change, because
  refitting widens the intervals (1.2-2.7x measured):
      Fig S7C    claim-1 60S per-bin difference
      Fig S7A    claim-1 60S/40S W1 against the null
      Fig S7D    claim-2 drug per-bin effects
      Fig 5F-r   zone per-bin effect

  DOES NOT need redrawing -- these show point estimates only, and the refit
  reproduces every point estimate exactly (verified to <1e-5), so the
  committed PDFs are already final as far as the STATISTICS go:
      Fig 5C, S7B, S7E, 5F-left        density overlays
      Fig S7E table                    drug deltas are point estimates
  (Presentation is separate: 5C, S7B, S7E and S7J were rebuilt in Sep 2026 to
  fill every curve rather than one. No number moved.)

Outputs carry a `_refit` suffix so the no-refit panels stay on disk for
comparison rather than being overwritten.
"""
import argparse, glob, json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
import paperfig_style as ps
from submit_panel_refit import CONTRASTS, OUT_ROOT

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
OUT = tlsd_paths.outp('PaperFigures_updated')


def as_bin_frame(key):
    """Refit stats CSV -> the column contract plot_bin_effect expects."""
    p = os.path.join(OUT_ROOT, key, f'{key}_refit_stats.csv')
    if not os.path.exists(p):
        return None
    d = pd.read_csv(p)
    d = d[d.quantity.str.startswith('bin_')].reset_index(drop=True)
    return pd.DataFrame(dict(
        label=[q[4:] for q in d.quantity],
        observed_diff=d.observed, ci_lo=d.ci_lo, ci_hi=d.ci_hi,
        stars_bh=[ps.stars(x) for x in d.p_bh]))


def draws(key):
    p = os.path.join(OUT_ROOT, key, f'{key}_refit_draws.npz')
    return np.load(p) if os.path.exists(p) else None


def w1_panel(keys_labels, out_pdf):
    series = []
    for key, label in keys_labels:
        d = draws(key)
        if d is None:
            print(f'  skip {out_pdf}: {key} not aggregated'); return
        st = pd.read_csv(os.path.join(OUT_ROOT, key, f'{key}_refit_stats.csv'))
        w1 = float(st.loc[st.quantity == 'W1', 'observed'].iloc[0])
        series.append(dict(label=label, w1=w1,
                           boot=d['boot_w1'], null=d['null_w1']))
    pc.plot_w1_vs_null(series, out_pdf)
    print(f'  wrote {pc.tag(out_pdf)}')


def bin_panel(specs, out_pdf, ylabel, title=None, bracket_csv=None):
    series = []
    for key, label, color in specs:
        df = as_bin_frame(key)
        if df is None:
            print(f'  skip {out_pdf}: {key} not aggregated'); return
        series.append(dict(df=df, label=label, color=color, marker='^'))
    brackets = None
    if bracket_csv and os.path.exists(bracket_csv):
        b = pd.read_csv(bracket_csv)
        brackets = {i: s for i, s in enumerate(b.stars.fillna(''))
                    if isinstance(s, str) and s}
    pc.plot_bin_effect(series, out_pdf, ylabel=ylabel, title=title,
                       brackets=brackets)
    print(f'  wrote {pc.tag(out_pdf)}')


def figS7D_combined():
    """S7D, all four drug contrasts on ONE axis.

    Colour separates genotype (grey family sgNT, orange family sgUBA5) so the
    genotype comparison is direct rather than split across two panels; shade
    and marker separate drug within each genotype. plot_bin_effect spaces k
    series across +-0.16 bin widths, so four fit without overlapping.
    """
    specs = [('claim2_NT_HHT',   'sgNT + HHT',    ps.CONTROL,        '^'),
             ('claim2_NT_eme',   'sgNT + eme',    ps.CONTROL_LIGHT,  'o'),
             ('claim2_UBA5_HHT', 'sgUBA5 + HHT',  ps.PERTURBED,      '^'),
             ('claim2_UBA5_eme', 'sgUBA5 + eme',  ps.PERTURBED_LIGHT,'o')]
    series = []
    for key, label, color, marker in specs:
        df = as_bin_frame(key)
        if df is None:
            print(f'  skip figS7D_combined: {key} not aggregated'); return
        series.append(dict(df=df, label=label, color=color, marker=marker))
    out = f'{OUT}/figureS7/D/figS7D_drug_shift_both_genotypes_refit.pdf'
    pc.plot_bin_effect(series, out, ylabel=ps.YLABEL_DELTA,
                       size='effect_wide_tall', legend_ncol=4,
                       title='drug effect on 60S diffusion, both genotypes')
    print(f'  wrote {pc.tag(out)}')


def figS7E():
    """S7E density overlays, coloured by GENOTYPE family rather than by drug.

    The manuscript's convention is grey for the control genotype and orange for
    sgUBA5, so colouring these two panels by drug (black / orange / blue in
    both) fought the rest of the figure. Each panel now uses one hue family
    with three shades, darkest untreated, so the two panels are visually
    parallel and the genotype is readable from colour alone.

    Redrawn from the pmfs stored by aggregate_panel_refit, so no refitting:
    every condition appears as the reference arm or the perturbed arm of one of
    the claim-2 contrasts.
    """
    from submit_panel_refit import NSUB
    from region_schemes import SCHEME_60S
    # condition -> (contrast, which arm of it)
    SRC = {'S3_NT_none':   ('claim2_NT_HHT', 'pmf_a'),
           'S3_NT_HT10':   ('claim2_NT_HHT', 'pmf_b'),
           'S3_NT_eme360': ('claim2_NT_eme', 'pmf_b'),
           'S3_UBA5_none':   ('claim2_UBA5_HHT', 'pmf_a'),
           'S3_UBA5_HT10':   ('claim2_UBA5_HHT', 'pmf_b'),
           'S3_UBA5_eme360': ('claim2_UBA5_eme', 'pmf_b')}
    FAMILY = {'NT':   [ps.BLACK, ps.CONTROL, ps.CONTROL_LIGHT],
              'UBA5': [ps.PERTURBED_DARK, ps.PERTURBED, ps.PERTURBED_LIGHT]}
    LABEL = {'none': 'untreated', 'HT10': 'harringtonine', 'eme360': 'emetine'}
    GENO = {'NT': 'sgNT', 'UBA5': 'sgUBA5'}

    dc = None
    for geno in ('NT', 'UBA5'):
        curves = []
        for i, drug in enumerate(('none', 'HT10', 'eme360')):
            cond = f'S3_{geno}_{drug}'
            key, arm = SRC[cond]
            f = os.path.join(OUT_ROOT, key, f'{key}_refit_draws.npz')
            if not os.path.exists(f):
                print(f'  skip figS7E: {key} not aggregated'); return
            z = np.load(f)
            curves.append(dict(pmf=z[arm], color=FAMILY[geno][i],
                               label=f'{LABEL[drug]} (n={NSUB[cond]:,})',
                               fill=True))
        if dc is None:
            # read the grid from the deployed settings rather than
            # reconstructing it: logspace(-3, 2, 150) happens to match today,
            # but a grid change would then silently mis-plot the x-axis
            import pickle
            from submit_panel_refit import TIERS
            _s = pickle.load(open(pc.settings_path(TIERS['claim2']), 'rb'))
            dc = np.asarray((_s.get('saspt') or _s)['diff_coefs'])
            assert len(dc) == len(curves[0]['pmf']), (
                f'grid is {len(dc)} states but the stored pmf has '
                f'{len(curves[0]["pmf"])}')
        out = f'{OUT}/figureS7/E/figS7E_overlay_{GENO[geno]}_refit.pdf'
        pc.plot_overlay(curves, dc, SCHEME_60S, out,
                        xlabel=r'Estimated Diffusion Coefficient 60S (µm$^2$/sec)')
        print(f'  wrote {pc.tag(out)}')


def figS7D_genotype():
    """S7D third sub-panel: sgUBA5 - sgNT at each drug level, one axis.

    Colour separates drug, matching DRUG_COLOR in build_claim2_panels
    (untreated black, harringtonine orange, emetine blue), so this reads as
    "does the genotype effect survive each perturbation".
    """
    specs = [('claim2_geno_none', 'untreated',     ps.BLACK,     'o'),
             ('claim2_geno_HHT',  'harringtonine', ps.PERTURBED, '^'),
             ('claim2_geno_eme',  'emetine',       ps.BLUE,      's')]
    series = []
    for key, label, color, marker in specs:
        df = as_bin_frame(key)
        if df is None:
            print(f'  skip figS7D_genotype: {key} not aggregated'); return
        series.append(dict(df=df, label=label, color=color, marker=marker))
    out = f'{OUT}/figureS7/D/figS7D_genotype_effect_by_drug_refit.pdf'
    pc.plot_bin_effect(series, out, ylabel=ps.YLABEL_DELTA,
                       size='effect_wide_tall', legend_ncol=3,
                       title='sgUBA5 - sgNT, by drug')
    print(f'  wrote {pc.tag(out)}')


def fig5F():
    """5F-right: zone effect per treatment, with the interaction as brackets.

    Reads the 5F stats frame rather than a per-contrast one -- the zone
    resampling has its own layout (aggregate_zone_refit.py), because its
    permutation is trajectory-level and refit-free while its bootstrap refits.
    """
    csv = pc.tag(f'{OUT}/figure5/F/fig5F_refit_stats.csv')
    if not os.path.exists(csv):
        print('  skip fig5F: not aggregated'); return
    d = pd.read_csv(csv)
    series = []
    for treat, color in (('DMSO', ps.CONTROL), ('PT-179 (24h)', ps.PERTURBED)):
        z = d[(d.panel == '5F_zone_effect') & (d.treatment == treat)]
        series.append(dict(
            df=pd.DataFrame(dict(label=z.bin.to_numpy(),
                                 observed_diff=z.observed.to_numpy(),
                                 ci_lo=z.ci_lo.to_numpy(),
                                 ci_hi=z.ci_hi.to_numpy(),
                                 stars_bh=[ps.stars(x) for x in z.p_bh])),
            # the contrast direction is in the title; repeating it per series
            # overflows the legend at this panel width
            label=treat, color=color, marker='^'))
    inter = d[d.panel == '5F_interaction'].reset_index(drop=True)
    brackets = {i: ps.stars(r.p_bh) for i, r in inter.iterrows()
                if ps.stars(r.p_bh)}
    out = f'{OUT}/figure5/F/fig5F_right_zone_bin_effect_refit.pdf'
    pc.plot_bin_effect(series, out, ylabel=ps.YLABEL_DELTA,
                       title='60S on ER - off ER, by treatment',
                       brackets=brackets)
    print(f'  wrote {pc.tag(out)}')


PANELS = {
    'fig5F': fig5F,
    'figS7D_genotype': figS7D_genotype,
    'figS7E': figS7E,
    'figS7A': lambda: w1_panel(
        [('claim1_60S', 'RPL29-Halo (60S)'), ('claim1_40S', 'Halo-RPS24 (40S)')],
        f'{OUT}/figureS7/A/figS7A_W1_60S_40S_refit.pdf'),
    'figS7C': lambda: bin_panel(
        [('claim1_60S', 'sgUBA5 - sgNT', ps.PERTURBED)],
        f'{OUT}/figureS7/C/figS7C_60S_bin_difference_refit.pdf',
        ps.YLABEL_DELTA),
    'figS7D': lambda: [
        bin_panel([(k, lab, col)],
                  f'{OUT}/figureS7/D/figS7D_{name}_refit.pdf', ps.YLABEL_DELTA,
                  title=title)
        for k, lab, col, name, title in [
            ('claim2_NT_HHT', 'HHT - untreated', ps.PERTURBED,
             'drug_shift_sgNT_HHT', 'sgNT: harringtonine'),
            ('claim2_NT_eme', 'emetine - untreated', ps.BLUE,
             'drug_shift_sgNT_eme', 'sgNT: emetine'),
            ('claim2_UBA5_HHT', 'HHT - untreated', ps.PERTURBED,
             'drug_shift_sgUBA5_HHT', 'sgUBA5: harringtonine'),
            ('claim2_UBA5_eme', 'emetine - untreated', ps.BLUE,
             'drug_shift_sgUBA5_eme', 'sgUBA5: emetine')]]
             + [figS7D_combined()],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', action='append', choices=sorted(PANELS))
    ap.add_argument('--all', action='store_true')
    a = ap.parse_args()
    ps.apply()
    for name in (sorted(PANELS) if a.all else (a.panel or [])):
        print(f'{name}:')
        PANELS[name]()
    return 0


if __name__ == '__main__':
    sys.exit(main())
