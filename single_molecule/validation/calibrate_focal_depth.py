"""
calibrate_focal_depth.py
========================

Pick a defensible `focal_depth` (dz) for the saSPT defocalization correction,
and quantify how much any published conclusion depends on that choice.

Why this is needed
------------------
saSPT applies the defocalization correction unconditionally
(`saspt/sa.py:374`): each diffusion state's occupancy is divided by the
fraction of particles of that D still inside the focal slab one frame later.
The size of that correction is set entirely by `focal_depth`. All four datasets
here inherited `focal_depth = 0.1` um from an `init_settings.py` default, which
is a TIRF-like number; the imaging is HiLo on a Nikon Ti2 / Ilas2 with a Prime
95B at 0.11 um/px. At dz = 0.1 um the correction inflates the D = 10 um^2/s
occupancy ~14-fold and D = 100 ~33-fold; at dz = 0.7 um those are 2.2x and
4.7x. So the number matters, and it was never measured.

Three independent handles are produced here:

1. `optical_dof()`   -- Berek's depth-of-field formula. A physical UPPER bound
                        on the diffraction-limited focal slab.
2. `empirical_dz()`  -- fit dz from the data itself: fast particles leave the
                        slab, so observed track length falls with D. Taking the
                        ratio against the slowest D bin cancels photobleaching
                        (which is D-independent). This is a LOWER bound, since
                        linking failure and motion blur also shorten fast
                        tracks and are absorbed into the same decay.
3. `sweep_dz()`      -- recompute every dz-sensitive statistic across a grid of
                        dz and report which conclusions are invariant.

The sweep is cheap because of a structural fact verified by
`test_saspt_equivalence.py` T2(iii): saSPT applies the correction AFTER the
VBEM loop, so the converged `exp_log_tau` -- and therefore the whole fit -- is
independent of dz. One VBEM fit per group serves the entire sweep; only the
final division by `f_remain` changes. Per-FOV soft occupancy (the bin panels)
is built from `exp_log_tau` and so does not move with dz at all; only the
population posterior (overlay curves) and W1 do.

Run via SLURM -- `empirical_dz`/`sweep_dz` load caches and run VBEM:
    sbatch call_calibrate_focal_depth.slurm
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from saspt.defoc import f_remain_rbm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pipeline')))  # the modules under test live in ../pipeline
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
import occupancy_compare as oc
from saspt_likelihood_cache import load_likelihood_cache, _load_saspt_settings

DZ_GRID = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, np.inf]

# Nikon Ti2 / Ilas2 / Photometrics Prime 95B, 0.11 um/px => 100x with 11 um
# camera pixels. NA is NOT recorded in the ND2 metadata (only "Ilas2 TIRF
# Preset: TIRF 561" and "Modality: Widefield Fluorescence"), so the objective
# NA is the one number here that must be confirmed by hand. 1.49 is the
# standard 100x Apo TIRF pairing for this setup.
OPTICS = dict(NA=1.49, wavelength_um=0.58, n_immersion=1.515,
              magnification=100.0, camera_pixel_um=11.0)


def optical_dof(NA, wavelength_um, n_immersion, magnification, camera_pixel_um):
    """Berek's depth of field: the diffraction term plus the detector term.

        DOF = lambda * n / NA^2  +  n * e / (M * NA)

    Returns um. This is the depth over which the image stays diffraction-
    limited; the axial range over which a single molecule is still detectable
    and localizable is somewhat larger, so treat this as a soft lower bound on
    a reasonable `focal_depth` rather than a hard ceiling.
    """
    diffraction = wavelength_um * n_immersion / NA ** 2
    detector = n_immersion * camera_pixel_um / (magnification * NA)
    return dict(dof_um=diffraction + detector, diffraction_term_um=diffraction,
                detector_term_um=detector, **dict(
                    NA=NA, wavelength_um=wavelength_um,
                    n_immersion=n_immersion, magnification=magnification,
                    camera_pixel_um=camera_pixel_um))


# ---------------------------------------------------------------------------
# empirical dz from the track-length / D relationship
# ---------------------------------------------------------------------------

def _map_diff_coef(cache):
    """Per-subtrack MAP diffusion coefficient from the label-blind likelihood.

    Uses the cached naive (prior-free) likelihood marginalized over the
    localization-error axis, so the D assignment does not depend on any fit or
    on which group the subtrack belongs to.
    """
    L = cache['L']                       # (n_sub, n_D, n_loc)
    marg = L.sum(axis=2)
    return np.asarray(cache['diff_coefs'])[marg.argmax(axis=1)]


def empirical_dz(cache_path, settings_pkl, n_d_bins=6, min_per_bin=2000,
                 n_boot=200, seed=0, max_jumps=10):
    """Fit dz from the D-dependence of observed track length.

    Model
    -----
    A subtrack of diffusion coefficient D survives to n jumps only if it stays
    both photostable and inside the focal slab:

        S_obs(n | D) = S_photo(n) * f_remain(D, n, dt, dz)

    Photobleaching and blinking are D-independent, so dividing by the slowest D
    bin (where f_remain ~ 1) cancels S_photo:

        S_obs(n | D) / S_obs(n | D_ref)  ~  f_remain(D, n) / f_remain(D_ref, n)

    dz is fit by least squares on the log of that ratio over all (D bin, n).

    Caveat, stated in the output
    ----------------------------
    Fast particles are ALSO lost to linking failure (a jump beyond the tracker's
    search radius) and to motion blur. Both shorten fast tracks and are absorbed
    into the same decay, so the fitted dz is biased SHALLOW -- it is a lower
    bound on the true focal depth, i.e. an upper bound on the correction.
    """
    s = _load_saspt_settings(settings_pkl)
    dt = float(s['frame_interval'])
    cache = load_likelihood_cache(cache_path, settings_pkl=settings_pkl)

    D_map = _map_diff_coef(cache)
    jumps = np.asarray(cache['jumps_per_track']).astype(int)
    ok = jumps >= 1
    D_map, jumps = D_map[ok], jumps[ok]

    # Log-spaced D bins over the populated range, merged until each holds
    # enough subtracks for a stable survival curve.
    logD = np.log10(D_map)
    qs = np.linspace(0, 100, n_d_bins + 1)
    edges = np.unique(np.percentile(logD, qs))
    idx = np.clip(np.digitize(logD, edges[1:-1]), 0, len(edges) - 2)

    bins = []
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() < min_per_bin:
            continue
        j = jumps[m]
        surv = np.array([(j >= n).mean() for n in range(1, max_jumps + 1)])
        bins.append(dict(D=float(np.median(D_map[m])), n=int(m.sum()),
                         surv=surv, mask=m))
    if len(bins) < 2:
        raise RuntimeError('not enough populated D bins to fit dz')

    ref = bins[0]                       # slowest bin is the photobleaching ref
    ns = np.arange(1, max_jumps + 1)

    def _resid(dz):
        r = []
        for b in bins[1:]:
            obs = np.log(np.maximum(b['surv'], 1e-12)
                         / np.maximum(ref['surv'], 1e-12))
            fr_b = f_remain_rbm(b['D'], max_jumps, dt, dz)
            fr_r = f_remain_rbm(ref['D'], max_jumps, dt, dz)
            pred = np.log(np.maximum(fr_b, 1e-12) / np.maximum(fr_r, 1e-12))
            good = (b['surv'] > 0.01) & (ref['surv'] > 0.01)
            r.append(obs[good] - pred[good])
        return float(np.sum(np.concatenate(r) ** 2))

    fit = minimize_scalar(_resid, bounds=(0.05, 5.0), method='bounded')
    dz_hat = float(fit.x)

    # Subtrack bootstrap for a CI on dz.
    rng = np.random.default_rng(seed)
    n_sub = D_map.size
    draws = []
    for _ in range(n_boot):
        take = rng.integers(0, n_sub, n_sub)
        Db, jb = D_map[take], jumps[take]
        lb = np.log10(Db)
        ib = np.clip(np.digitize(lb, edges[1:-1]), 0, len(edges) - 2)
        bb = []
        for b in range(len(edges) - 1):
            m = ib == b
            if m.sum() < min_per_bin:
                continue
            j = jb[m]
            bb.append(dict(D=float(np.median(Db[m])),
                           surv=np.array([(j >= n).mean() for n in ns])))
        if len(bb) < 2:
            continue
        rf = bb[0]

        def _rb(dz, bb=bb, rf=rf):
            r = []
            for b in bb[1:]:
                obs = np.log(np.maximum(b['surv'], 1e-12)
                             / np.maximum(rf['surv'], 1e-12))
                pred = np.log(np.maximum(f_remain_rbm(b['D'], max_jumps, dt, dz), 1e-12)
                              / np.maximum(f_remain_rbm(rf['D'], max_jumps, dt, dz), 1e-12))
                good = (b['surv'] > 0.01) & (rf['surv'] > 0.01)
                r.append(obs[good] - pred[good])
            return float(np.sum(np.concatenate(r) ** 2))

        draws.append(float(minimize_scalar(_rb, bounds=(0.05, 5.0),
                                           method='bounded').x))

    draws = np.asarray(draws)
    return dict(
        cache=os.path.basename(cache_path),
        dz_hat_um=dz_hat,
        dz_ci95=[float(np.percentile(draws, 2.5)),
                 float(np.percentile(draws, 97.5))] if draws.size else None,
        n_boot_ok=int(draws.size),
        frame_interval_s=dt,
        n_subtracks=int(D_map.size),
        d_bins=[dict(D_median=b['D'], n=b['n'], survival=b['surv'].tolist())
                for b in bins],
        interpretation=('LOWER bound on focal depth: linking failure and motion '
                        'blur also shorten fast tracks and are absorbed into '
                        'the same decay, biasing the fit shallow.'))


# ---------------------------------------------------------------------------
# dz sweep -- one VBEM fit per group serves the whole grid
# ---------------------------------------------------------------------------

def _n_raw_for_group(cache, likelihood, fovs, conc_param=1.0, max_iter=200,
                     equal_weight=False):
    """Uncorrected Dirichlet sufficient statistic for one group of FOVs.

    Run with apply_defoc=False so `n` is the raw statistic; the correction is
    a pure post-hoc division (verified bitwise in test_saspt_equivalence T2),
    so every dz in the sweep reuses this single fit.
    """
    mask = np.isin(cache['fov_of_subtrack'], fovs)
    L_sub = cache['L'][mask]
    jumps = np.asarray(cache['jumps_per_track'])[mask].astype(np.float32)
    prior = np.full(L_sub.shape[1:], conc_param, dtype=np.float32)
    w = np.ones_like(jumps) if equal_weight else None
    n_raw, _, _ = oc._vbem_infer(L_sub, jumps, prior, likelihood, max_iter,
                                 weights=w, apply_defoc=False)
    return np.asarray(n_raw, dtype=np.float64), int(mask.sum())


def _pmf_at_dz(n_raw, diff_coefs, dt, dz):
    """Apply saSPT's own correction at an arbitrary dz, then marginalize.

    Reproduces `defoc_corr_rbme` exactly: divide along the diffusion-coefficient
    axis by the one-frame focal-slab survival, then normalize.
    """
    if np.isinf(dz):
        n = n_raw.copy()
    else:
        fr = np.array([f_remain_rbm(d, 1, dt, dz)[0] for d in diff_coefs])
        n = (n_raw.T / fr).T
    p = n.sum(axis=1)
    return p / p.sum()


def sweep_dz(cache_path, settings_pkl, group_a, group_b, name_a='A', name_b='B',
             dz_grid=DZ_GRID, equal_weight=False, max_iter=200):
    """W1 between two groups' population posteriors across a grid of dz.

    Only two VBEM fits total, regardless of how many dz values are swept.
    """
    s = _load_saspt_settings(settings_pkl)
    dt = float(s['frame_interval'])
    cache = load_likelihood_cache(cache_path, settings_pkl=settings_pkl)
    likelihood = oc._make_likelihood_from_cache(cache)
    D = np.asarray(cache['diff_coefs'])
    log_grid = np.log10(D)

    fa = oc._resolve_group(cache, group_a)
    fb = oc._resolve_group(cache, group_b)
    n_a, ns_a = _n_raw_for_group(cache, likelihood, fa, max_iter=max_iter,
                                 equal_weight=equal_weight)
    n_b, ns_b = _n_raw_for_group(cache, likelihood, fb, max_iter=max_iter,
                                 equal_weight=equal_weight)

    rows = []
    for dz in dz_grid:
        pa = _pmf_at_dz(n_a, D, dt, dz)
        pb = _pmf_at_dz(n_b, D, dt, dz)
        rows.append(dict(
            dz_um=(None if np.isinf(dz) else float(dz)),
            w1=float(oc._w1_logD(log_grid, pa, pb)),
            median_logD_a=float(np.interp(0.5, np.cumsum(pa), log_grid)),
            median_logD_b=float(np.interp(0.5, np.cumsum(pb), log_grid))))
    return dict(contrast=f'{name_a}_vs_{name_b}',
                n_fov=[len(fa), len(fb)], n_subtracks=[ns_a, ns_b],
                equal_weight=bool(equal_weight), sweep=rows)


# ---------------------------------------------------------------------------

DATASETS = {
    'fig5C_60S': dict(
        tier=tlsd_paths.spt('all_data', '2025_06_06_zeb', 'tracking_output_q=15p0_splitsize10'),
        cache='occupancy_compare_cache/60S_bothgeno_likelihood_cache.npz',
        group_a={'source_member': 'sgNT_60S'}, group_b={'source_member': 'sgUBA5_60S'},
        name_a='sgNT_60S', name_b='sgUBA5_60S'),
    'fig5F_ER': dict(
        tier=tlsd_paths.spt('all_data', '2025_09_15_zeb_ER', 'tracking_output_q=15p0'),
        cache='occupancy_compare_cache/_standalone_noDegron_likelihood_cache.npz',
        group_a=None, group_b=None, name_a='noDegron', name_b=None),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--datasets', nargs='*', default=list(DATASETS))
    ap.add_argument('--n-boot', type=int, default=200)
    ap.add_argument('--skip-sweep', action='store_true')
    ap.add_argument('--out', default=tlsd_paths.outp('PaperFigures_updated', 'verification',
                                                     'focal_depth_calibration.json'))
    args = ap.parse_args()

    out = dict(optical_bound=optical_dof(**OPTICS),
               optics_note=('NA is NOT in the ND2 metadata (only "Ilas2 TIRF '
                            'Preset: TIRF 561", Nikon Ti2, Prime 95B, '
                            '0.11 um/px => 100x). NA=1.49 assumed; confirm.'),
               empirical={}, sweeps={})

    for key in args.datasets:
        cfg = DATASETS[key]
        spkl = os.path.join(cfg['tier'], 'settings.pkl')
        cpath = os.path.join(cfg['tier'], cfg['cache'])
        try:
            print(f'[{key}] empirical dz ...', flush=True)
            out['empirical'][key] = empirical_dz(cpath, spkl, n_boot=args.n_boot)
        except Exception as exc:
            out['empirical'][key] = dict(error=f'{type(exc).__name__}: {exc}')
            print(f'  FAILED: {exc}', file=sys.stderr, flush=True)

        if not args.skip_sweep and cfg.get('group_b'):
            try:
                print(f'[{key}] dz sweep ...', flush=True)
                out['sweeps'][key] = sweep_dz(
                    cpath, spkl, cfg['group_a'], cfg['group_b'],
                    cfg['name_a'], cfg['name_b'])
            except Exception as exc:
                out['sweeps'][key] = dict(error=f'{type(exc).__name__}: {exc}')
                print(f'  FAILED: {exc}', file=sys.stderr, flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as fh:
        json.dump(out, fh, indent=2, default=oc._json_default)

    b = out['optical_bound']
    print(f"\noptical DOF (Berek, NA={b['NA']}) = {b['dof_um']:.3f} um "
          f"(diffraction {b['diffraction_term_um']:.3f} + "
          f"detector {b['detector_term_um']:.3f})")
    for k, v in out['empirical'].items():
        if 'error' in v:
            print(f"  {k}: ERROR {v['error']}")
        else:
            ci = v['dz_ci95']
            print(f"  {k}: empirical dz = {v['dz_hat_um']:.3f} um "
                  f"(95% CI {ci[0]:.3f}-{ci[1]:.3f}) [LOWER bound]")
    for k, v in out['sweeps'].items():
        if 'error' in v:
            print(f"  sweep {k}: ERROR {v['error']}")
            continue
        print(f"  sweep {k} ({v['contrast']}):")
        for r in v['sweep']:
            dz = 'inf' if r['dz_um'] is None else f"{r['dz_um']:.1f}"
            print(f"      dz={dz:>4}  W1={r['w1']:.4f}")
    print(f'\nwrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
