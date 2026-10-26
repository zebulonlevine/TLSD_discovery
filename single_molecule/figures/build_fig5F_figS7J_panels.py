"""
build_fig5F_figS7J_panels.py -- Figure 5F and Supplementary Figure S7J.

Dataset `2025_09_15_zeb_ER`, tier `q=15p0`, 15.091 ms / 0.11 um.
RPL29-Halo (60S) co-imaged with SNAP-Sec61B. `noDegron` = DMSO (10 FOV),
`Degron` = 24 h PT-179 (14 FOV). "on ER" is the `on_ER_broad` zone
(on_structure + near_ER); "off ER" is the positively-defined `off_ER` zone
(>=95% spatial purity outside the broad mask). Trajectories crossing the mask
boundary are `intermediate`/`ambiguous` and enter neither zone.

Panels
------
  5F left   density overlay, 4 curves: {DMSO, PT-179} x {on ER, off ER}
  5F right  per-bin (on ER - off ER) difference within each treatment, with a
            FOV-bootstrap 95% CI
  S7J       density overlay ignoring zone, DMSO vs PT-179 (the same
            sub-trajectories pooled across zones)

Why the null here is not the FOV-relabelling null used elsewhere
---------------------------------------------------------------
Every FOV contains both zones, so relabelling whole FOVs cannot generate a
zone contrast. The null is instead a within-FOV, TRAJECTORY-level zone-label
permutation (each FOV's zone counts preserved, all of a trajectory's
sub-trajectories moving together) -- `forkb_er_zone.zone_refit_permutation`.
Permuting at the sub-trajectory level is anticonservative and produced a real
calibration failure that `run_v4_forkb_zone_negctrl.py` caught.

Run:  sbatch call_build_panels.slurm build_fig5F_figS7J_panels.py
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
import paperfig_style as ps
import er_zone_naive_cache as ezc
import forkb_er_zone as fbz
from occupancy_compare import _subband_summary, _make_likelihood_from_cache, _vbem_infer
from softpost_region_panels import did_from_occ
from saspt_likelihood_cache import _load_saspt_settings, _grid_checksum
from multiplicity import bh, holm, perm_p_addone
from region_schemes import SCHEME_60S

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
ROOT = tlsd_paths.spt('all_data', '2025_09_15_zeb_ER', 'tracking_output_q=15p0')
# Resolved through panel_core so PANEL_VARIANT selects the matching settings
# block and cache directory, rather than being hard-coded here.
SETTINGS_PKL = pc.settings_path(ROOT)
CACHE_DIR = os.path.join(ROOT, pc.CACHE_SUBDIR)
# The variant caches are named by condition; the default ones carry the
# _standalone_ prefix from the batch scripts that built them.
CACHE_NAME = (lambda d: d) if pc.VARIANT else (lambda d: f'_standalone_{d}')
OUT = tlsd_paths.outp('PaperFigures_updated')

ZONE_ON, ZONE_OFF = 'on_ER_broad', 'off_ER'
# Display names: the figure speaks in treatments, the data in degron state.
TREATMENT = {'noDegron': 'DMSO', 'Degron': 'PT-179 (24h)'}
COLOR = {'noDegron': ps.CONTROL, 'Degron': ps.PERTURBED}
COLOR_OFF = {'noDegron': ps.CONTROL_LIGHT, 'Degron': ps.PERTURBED_LIGHT}


def zone_masks(cache, traj_zone):
    traj_id = np.asarray(cache['orig_trajectory']).astype(int)
    zone = np.array([traj_zone.get(t, 'unassigned') for t in traj_id])
    return zone


def per_fov_occ_by_zone(ctx, cache, zone):
    """Per-FOV soft occupancy for each zone, on a shared FOV index.

    One VBEM fit over the POOLED sub-trajectories of both zones (so the two
    zones share a single posterior tilt and are directly comparable), then the
    per-FOV mean of each sub-trajectory's marginal, taken over that zone's
    sub-trajectories only.
    """
    keep_pool = np.isin(zone, [ZONE_ON, ZONE_OFF])
    L_sub = cache['L'][keep_pool]
    jumps = np.asarray(cache['jumps_per_track'])[keep_pool].astype(np.float32)
    fov = cache['fov_of_subtrack'][keep_pool]
    zsub = zone[keep_pool]

    prior = np.full(L_sub.shape[1:], pc.CONC_PARAM, dtype=np.float32)
    w = np.ones_like(jumps) if pc.EQUAL_WEIGHT else None
    _, _, tau = _vbem_infer(L_sub, jumps, prior, ctx['likelihood'],
                            pc._VBEM_ITERS, tol=pc._VBEM_TOL, weights=w,
                            apply_defoc=False)

    n_sub, n_D, n_loc = L_sub.shape
    norm = L_sub.reshape(n_sub, -1) @ tau.reshape(-1)
    norm[norm == 0] = 1.0
    marg = np.zeros((n_sub, n_D), dtype=np.float64)
    for j in range(n_loc):
        marg += L_sub[:, :, j] * tau[None, :, j]
    marg /= norm[:, None]

    # Index BOTH zones on the same FOV list. The zones are paired within a FOV
    # (every FOV contains both), and the downstream contrast is a per-FOV
    # difference, so the two matrices must be row-aligned. Building each zone's
    # FOV list independently would silently misalign them whenever a FOV has
    # sub-trajectories in one zone but not the other.
    fovs = sorted(set(fov))
    out = {}
    for z in (ZONE_ON, ZONE_OFF):
        sel = zsub == z
        occ = np.full((len(fovs), n_D), np.nan)
        for i, f in enumerate(fovs):
            rows = sel & (fov == f)
            if rows.any():
                occ[i] = marg[rows].mean(axis=0)
        out[z] = (occ, int(sel.sum()))
    return fovs, out


def build(degron):
    # build_zone_labeled_cache returns the LOADED CACHE (not a path) plus the
    # trajectory -> zone map; it loads without a checksum assert, so the grid
    # is checked against the deployed settings.pkl explicitly here.
    cache, traj_zone = ezc.build_zone_labeled_cache(
        ROOT, SETTINGS_PKL, degron, [ZONE_ON, ZONE_OFF], CACHE_DIR,
        cache_name=CACHE_NAME(degron), verbose=False)
    settings = _load_saspt_settings(SETTINGS_PKL)
    live = _grid_checksum(settings)
    if live != cache['grid_checksum']:
        raise RuntimeError(
            f'grid checksum mismatch for {degron}: cache={cache["grid_checksum"][:12]} '
            f'vs {os.path.basename(SETTINGS_PKL)}={live[:12]} -- rebuild the cache')
    ctx = dict(cache=cache, likelihood=_make_likelihood_from_cache(cache),
               settings=settings, diff_coefs=np.asarray(cache['diff_coefs']))
    zone = zone_masks(ctx['cache'], traj_zone)
    fov_list, by_zone = per_fov_occ_by_zone(ctx, ctx['cache'], zone)
    occ_on, n_on = by_zone[ZONE_ON]
    occ_off, n_off = by_zone[ZONE_OFF]
    # zone-agnostic curve for S7J: every sub-trajectory of this condition,
    # including the intermediate/ambiguous ones the zone contrast excludes
    occ_all, _, n_all_per_fov = pc.per_fov_occupancy(
        ctx, list(ctx['cache']['fov_metadata']['fov_basename']))
    n_all = int(n_all_per_fov.sum())
    return dict(ctx=ctx, occ_on=occ_on, occ_off=occ_off, occ_all=occ_all,
                n_on=n_on, n_off=n_off, n_all=n_all,
                n_fov=len(fov_list), degron=degron)


def zone_bin_effect(ctx, occ_on, occ_off, degron, seed=0):
    """Per-bin (on ER - off ER) with a FOV bootstrap. The permutation column
    is left to `forkb_er_zone.zone_refit_permutation`, which carries the
    correct within-FOV trajectory-level null; here the FOV-paired bootstrap
    supplies the plotted interval."""
    rows = _subband_summary(ctx['diff_coefs'], occ_off, occ_on, SCHEME_60S,
                            None, None, seed=seed, n_boot=pc.N_BOOT,
                            n_perm=pc.N_PERM)
    df = pd.DataFrame(rows)
    df['p_perm'] = [perm_p_addone(p, pc.N_PERM) for p in df['p_value_two_sided']]
    df['p_bh'] = bh(df['p_perm'].to_numpy())
    df['p_holm'] = holm(df['p_perm'].to_numpy())
    df['stars_bh'] = [ps.stars(p) for p in df['p_bh']]
    df['treatment'] = TREATMENT[degron]
    df['contrast'] = f'{ZONE_ON} - {ZONE_OFF}'
    df['n_fov'] = occ_on.shape[0]
    df['equal_weight'] = pc.EQUAL_WEIGHT
    df['null_note'] = ('FOV bootstrap shown; inferential null is the within-FOV '
                       'trajectory-level zone permutation (forkb_er_zone)')
    return df


def zone_treatment_interaction(ctx, fits, seed=0):
    """Per-bin (on ER - off ER) difference between treatments, from per-FOV
    paired zone differences.

    The zone contrast is taken within each FOV, so the FOV is the unit and the
    zone main effect is differenced away before any resampling. The null then
    relabels treatment across FOVs, which is exact under "treatment does not
    change the zone difference".
    """
    rng = np.random.default_rng(seed)
    diff = ctx['diff_coefs']
    d = {k: fits[k]['occ_on'] - fits[k]['occ_off'] for k in ('noDegron', 'Degron')}
    rows = []
    for label, (lo, hi) in SCHEME_60S.items():
        m = (diff >= lo) & (diff < hi)
        a = np.nansum(d['noDegron'][:, m], axis=1)
        b = np.nansum(d['Degron'][:, m], axis=1)
        obs = float(b.mean() - a.mean())

        pooled, na = np.concatenate([a, b]), len(a)
        null = np.empty(pc.N_PERM)
        for i in range(pc.N_PERM):
            pm = rng.permutation(len(pooled))
            null[i] = pooled[pm[na:]].mean() - pooled[pm[:na]].mean()
        boot = np.empty(pc.N_BOOT)
        for i in range(pc.N_BOOT):
            boot[i] = (b[rng.integers(0, len(b), len(b))].mean()
                       - a[rng.integers(0, len(a), len(a))].mean())
        rows.append(dict(
            label=label, d_lo=lo, d_hi=hi, observed_diff=obs,
            zone_diff_DMSO=float(a.mean()), zone_diff_PT179=float(b.mean()),
            ci_lo=float(np.percentile(boot, 2.5)),
            ci_hi=float(np.percentile(boot, 97.5)),
            p_perm=float(perm_p_addone((np.abs(null) >= abs(obs)).mean(), pc.N_PERM)),
            null_sd=float(null.std()), boot_sd=float(boot.std()),
            n_fov_DMSO=len(a), n_fov_PT179=len(b)))
    df = pd.DataFrame(rows)
    df['p_bh'] = bh(df['p_perm'].to_numpy())
    df['p_holm'] = holm(df['p_perm'].to_numpy())
    df['stars_bh'] = [ps.stars(x) for x in df['p_bh']]
    df['contrast'] = '(on ER - off ER), PT-179 minus DMSO'
    df['null'] = 'treatment relabelling on per-FOV paired zone differences'
    return df


def main():
    ps.apply()
    fits = {d: build(d) for d in ('noDegron', 'Degron')}
    ctx0 = fits['noDegron']['ctx']
    dc = ctx0['diff_coefs']

    # --- Fig 5F left: four curves, colour = treatment, dash = zone ---
    curves = []
    for d in ('noDegron', 'Degron'):
        f = fits[d]
        curves.append(dict(pmf=np.nanmean(f['occ_on'], axis=0),
                           label=f"{TREATMENT[d]} on ER (n={f['n_on']:,})",
                           color=COLOR[d], ls='-'))
        curves.append(dict(pmf=np.nanmean(f['occ_off'], axis=0),
                           label=f"{TREATMENT[d]} off ER (n={f['n_off']:,})",
                           color=COLOR_OFF[d], ls='--'))
    pc.plot_overlay(curves, dc, SCHEME_60S,
                    f'{OUT}/figure5/F/fig5F_left_zone_posterior_overlay.pdf',
                    xlabel=r'Estimated Diffusion Coefficient 60S (µm$^2$/sec)')

    # --- Fig 5F right: per-bin (on ER - off ER), one series per treatment ---
    dfs, series = [], []
    for d in ('noDegron', 'Degron'):
        f = fits[d]
        df = zone_bin_effect(f['ctx'], f['occ_on'], f['occ_off'], d)
        dfs.append(df)
        series.append(dict(df=df, label=TREATMENT[d], color=COLOR[d], marker='o'))
    # Stars on this panel test the TREATMENT x ZONE interaction -- per the draft
    # legend, "p-values from FOV-resampling between conditions" -- not each
    # series against zero. That is the quantity the text actually claims: bin III
    # accumulates on-ER only under PT-179, with no difference untreated.
    #
    # Because every FOV contains both zones, the zone contrast can be taken
    # WITHIN a FOV first: d_f = occ_on(f) - occ_off(f). That per-FOV difference
    # already removes the zone main effect, so relabelling treatment on d_f is a
    # clean two-sample permutation with nothing left to contaminate the null --
    # the problem that forces a Freedman-Lane residual null when the four cells
    # of a 2x2 are separate FOVs does not arise here.
    did = zone_treatment_interaction(ctx0, fits)
    did.to_csv(pc.tag(f'{OUT}/figure5/F/fig5F_right_zone_treatment_DiD_stats.csv'),
               index=False)

    pc.plot_bin_effect(series, f'{OUT}/figure5/F/fig5F_right_zone_bin_effect.pdf',
                       ylabel=ps.YLABEL_DELTA + '\n(on ER - off ER)',
                       brackets={i: st for i, st in enumerate(did['stars_bh']) if st})
    print('\nzone x treatment interaction (the stars on 5F right):')
    print(did[['label', 'zone_diff_DMSO', 'zone_diff_PT179', 'observed_diff',
               'ci_lo', 'ci_hi', 'p_perm', 'p_bh', 'stars_bh']].round(4).to_string(index=False))
    pd.concat(dfs).to_csv(pc.tag(f'{OUT}/figure5/F/fig5F_right_zone_bin_effect_stats.csv'),
                          index=False)

    # --- Fig S7J: zone-agnostic overlay ---
    pc.plot_overlay(
        [dict(pmf=np.nanmean(fits['noDegron']['occ_all'], axis=0),
              label=f"DMSO (n={fits['noDegron']['n_all']:,})",
              color=ps.CONTROL, fill=True),
         dict(pmf=np.nanmean(fits['Degron']['occ_all'], axis=0),
              label=f"PT-179 24h (n={fits['Degron']['n_all']:,})",
              color=ps.PERTURBED, fill=True)],
        dc, SCHEME_60S, f'{OUT}/figureS7/J/figS7J_total_posterior_overlay.pdf',
        xlabel=r'Estimated Diffusion Coefficient 60S (µm$^2$/sec)')

    for df in dfs:
        print(f"\n{df['treatment'].iloc[0]}  (on ER - off ER)")
        print(df[['label', 'observed_diff', 'ci_lo', 'ci_hi',
                  'p_perm', 'p_bh', 'stars_bh']].to_string(index=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
