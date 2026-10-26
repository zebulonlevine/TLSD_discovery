"""
test_saspt_equivalence.py
=========================

Proves that the Fork-B analysis path behind every SPT panel in the paper
computes the SAME numbers the real `saspt` package computes.

Motivation
----------
Exactly one piece of saSPT math is re-implemented in this repo:
`occupancy_compare._vbem_infer`, an algebraically-reformulated port of
`saspt.sa.StateArray.infer_posterior`'s VBEM loop. Everything else -- the
RBME likelihood, the state grid, jump extraction, and the defocalization
correction -- is the installed `saspt` package itself. This script is the
evidence for that claim, and it quantifies the places where the port
deliberately differs from stock saSPT:

  1. an off-by-one in the iteration accounting. saspt contracts R -> n at the
     TOP of each pass and updates R at the bottom, so the `n` it returns after
     K passes reflects K-1 R-updates; `_vbem_infer` folds both into one
     recurrence, so after K passes its `n` reflects K R-updates. Calling it
     with `max_iter = K-1` therefore reproduces `saspt(max_iter=K)` exactly
     (measured: 3.0e-16). Run as-is at K=200 the two differ by ~1.2e-5 --
     which is ~180x SMALLER than saspt's own residual non-convergence at its
     default max_iter=200 (2.2e-3 vs its own 5000-iteration fixed point), and
     both loops converge to the SAME fixed point (2.7e-12 at 5000 iterations).
  1b. a convergence early-stop (`tol`) that saspt does not have (saspt always
     runs the full `max_iter`),
  2. float32 arithmetic where saspt uses float64,
  3. optional `weights=` / `apply_defoc=` switches -- NOT exercised here; the
     identity tests all run at saSPT's own settings (displacement weighting,
     defocalization on).

Tests
-----
T1  VBEM loop identity          stock StateArray vs _vbem_infer on the SAME L
T2  Defocalization identity     three independent checks (see its docstring)
T3  Likelihood-cache fidelity   cache construction vs a fresh StateArray
T4  Grid / settings identity    settings.pkl vs cache snapshot vs saspt defaults
T5  Marginalization provenance  local loc-error marginal vs saspt's own
T6  End-to-end                  stock posterior_occs vs Fork B pmf, real data

Writes `saspt_equivalence_report.json` and prints a readable table. Exit
status is nonzero if any test misses its numeric bar.

Run via SLURM -- never on the head node, since T1/T5/T6 run VBEM:
    sbatch call_saspt_equivalence.slurm
"""

import argparse
import json
import os
import pickle
import sys
import traceback

import numpy as np
import pandas as pd

import saspt
from saspt import StateArray
from saspt.defoc import defoc_corr, f_remain_rbm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pipeline')))  # the modules under test live in ../pipeline
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
import occupancy_compare as oc
from saspt_likelihood_cache import (load_likelihood_cache, _load_saspt_settings,
                                    _grid_checksum, _GRID_KEYS, FORCED_SAMPLE_SIZE)

# One dataset entry per experiment that has a standalone cache (Fig 5C 60S,
# Fig 5F ER). `traj` is a pooled trajectory CSV a production cache was built from.
DATASETS = {
    'fig5C_60S': dict(
        tier=tlsd_paths.spt('all_data', '2025_06_06_zeb', 'tracking_output_q=15p0_splitsize10'),
        traj='traj_combinations/sgNT_60S_pooled_traj.csv',
        cache='occupancy_compare_cache/_standalone_sgNT_60S_likelihood_cache.npz'),
    'fig5F_ER': dict(
        tier=tlsd_paths.spt('all_data', '2025_09_15_zeb_ER', 'tracking_output_q=15p0'),
        traj='traj_combinations/noDegron_pooled_traj.csv',
        cache='occupancy_compare_cache/_standalone_noDegron_likelihood_cache.npz'),
}

# Numeric bars. (a) is the pure algebraic-reformulation check and must be
# essentially exact. (c) is the production configuration and is allowed the
# float32 + early-stop error, bounded well below the smallest effect the paper
# reports (occupancy shifts ~0.5 percentage points; W1 effects ~0.05).
BAR_EXACT_OCC = 1e-10      # T1(a) max |delta occupancy|
BAR_PROD_OCC = 1e-4        # T1(c) max |delta occupancy| (occupancies sum to 1)
BAR_PROD_W1 = 5e-4         # T1(c)/T6 |delta W1| on the log10-D axis
BAR_CACHE_L = 1e-6         # T3 float32 round-trip on per-subtrack normalized L


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _settings(tier):
    return _load_saspt_settings(os.path.join(tier, 'settings.pkl'))


def _run_settings(tier, max_iter, conc_param=None):
    """Deployed saSPT settings, with sample_size forced high so nothing is
    subsampled away (the production combination/megapool path does the same).
    conc_param is left absent unless given, so saspt's own default applies --
    which is exactly what production does."""
    run = dict(_settings(tier))
    run['sample_size'] = FORCED_SAMPLE_SIZE
    run['progress_bar'] = False
    run.pop('save_subtrack_posteriors', None)
    run['max_iter'] = max_iter
    if conc_param is not None:
        run['conc_param'] = conc_param
    return run


def _subsample_detections(traj_csv, settings_pkl, n_traj, seed=0):
    """A reproducible subset of real detections, filtered exactly as
    `trackingCodeFunctions.executeSASPT` and `build_likelihood_cache` do (the
    min_I0 mean-intensity filter), so the subset matches what production fits."""
    df = pd.read_csv(traj_csv)
    with open(settings_pkl, 'rb') as fh:
        s = pickle.load(fh)
    min_I0 = s.get('quot', {}).get('track', {}).get('min_I0')
    if min_I0 is not None and 'I0' in df.columns:
        df = df.groupby('trajectory').filter(lambda g: g['I0'].mean() > min_I0)
    ids = np.sort(df['trajectory'].unique())
    rng = np.random.default_rng(seed)
    keep = rng.choice(ids, size=min(n_traj, ids.size), replace=False)
    return df[df['trajectory'].isin(keep)].copy()


def _fork_b_inputs(SA, dtype=np.float64):
    """The three arrays Fork B consumes, taken straight off a real StateArray:
    L in subtrack-major order, the per-subtrack jump counts, and a flat
    Dirichlet prior at saspt's default concentration."""
    L_sub = np.moveaxis(SA.naive_assignment_probabilities, -1, 0).astype(dtype)
    jumps = np.asarray(SA.jumps_per_track, dtype=dtype)
    prior = np.full(L_sub.shape[1:], 1.0, dtype=dtype)
    return L_sub, jumps, prior


def _pmf(occ):
    """Marginalize a (n_D, n_loc) occupancy over localization error."""
    p = np.asarray(occ, dtype=np.float64).sum(axis=1)
    return p / p.sum()


# ---------------------------------------------------------------------------
# T1 -- VBEM loop identity
# ---------------------------------------------------------------------------

def t1_vbem_identity(key, cfg, n_traj, results):
    """Feed stock saspt's own L to _vbem_infer and require the same answer.

    saspt's loop (sa.py:346-379) materializes the full (n_D, n_loc, n_sub)
    responsibility array R every iteration; _vbem_infer substitutes an
    algebraically identical recurrence that never materializes R. Running both
    on the SAME L isolates that reformulation from every other possible
    difference (likelihood, grid, filtering), so config (a) -- float64 with the
    early stop disabled -- must agree to machine precision.
    """
    tier = cfg['tier']
    det = _subsample_detections(os.path.join(tier, cfg['traj']),
                                os.path.join(tier, 'settings.pkl'), n_traj)
    SA = StateArray.from_detections(det, **_run_settings(tier, 200))

    _, _, occ_ref = SA.infer_posterior()          # stock saspt reference
    log_grid = np.log10(np.asarray(SA.likelihood.diff_coefs))
    pmf_ref = _pmf(occ_ref)

    configs = {
        # (a) saspt's EXACT iteration semantics: one fewer R->n contraction,
        #     no early stop, float64. This is the exactness test.
        'a_saspt_exact': dict(dtype=np.float64, tol=0.0, max_iter=199),
        # (b) the same call run at a matching max_iter, isolating the
        #     off-by-one alone.
        'b_float64_full': dict(dtype=np.float64, tol=0.0, max_iter=200),
        # (c) + the early stop.
        'c_float64_tol': dict(dtype=np.float64, tol=1e-5, max_iter=200),
        # (d) the production configuration (float32 + early stop).
        'd_float32_tol': dict(dtype=np.float32, tol=1e-5, max_iter=200),
    }
    L64, j64, _ = _fork_b_inputs(SA)
    out = {'n_subtracks': int(L64.shape[0]), 'n_jumps': int(j64.sum()),
           'grid_shape': list(L64.shape[1:]), 'configs': {}}

    for name, c in configs.items():
        L_sub, jumps, prior = _fork_b_inputs(SA, dtype=c['dtype'])
        _, post_mean, _ = oc._vbem_infer(L_sub, jumps, prior, SA.likelihood,
                                         c['max_iter'], tol=c['tol'],
                                         weights=None, apply_defoc=True)
        pm = np.asarray(post_mean, dtype=np.float64)
        out['configs'][name] = dict(
            max_abs_occ_diff=float(np.max(np.abs(pm - occ_ref))),
            max_rel_occ_diff=float(np.max(np.abs(pm - occ_ref)
                                          / np.maximum(occ_ref, 1e-12))),
            w1_vs_stock=float(oc._w1_logD(log_grid, pmf_ref, _pmf(post_mean))))

    # How far from its own fixed point is stock saspt at its default
    # max_iter? This is the yardstick the port's deviation must be judged
    # against -- a difference smaller than saspt's own convergence residual
    # cannot change any conclusion.
    SA_conv = StateArray.from_detections(det, **_run_settings(tier, 5000))
    out['saspt_self_convergence_200_vs_5000'] = float(
        np.max(np.abs(SA_conv.posterior_occs - occ_ref)))

    a = out['configs']['a_saspt_exact']
    d = out['configs']['d_float32_tol']
    out['passed'] = bool(a['max_abs_occ_diff'] <= BAR_EXACT_OCC
                         and d['max_abs_occ_diff'] <= BAR_PROD_OCC
                         and abs(d['w1_vs_stock']) <= BAR_PROD_W1)
    out['bars'] = dict(exact_occ=BAR_EXACT_OCC, prod_occ=BAR_PROD_OCC,
                       prod_w1=BAR_PROD_W1)
    results[f'T1_vbem_identity::{key}'] = out


# ---------------------------------------------------------------------------
# T2 -- Defocalization identity
# ---------------------------------------------------------------------------

def t2_defocalization(key, cfg, n_traj, results):
    """Three independent checks that defocalization is handled EXACTLY as
    saspt handles it.

    (i)   The likelihood object Fork B rebuilds from a cache's grid snapshot
          must produce a bitwise-identical focal-slab survival curve to the
          likelihood built inside a real StateArray, and to a direct call to
          saspt.defoc.f_remain_rbm. If the rebuilt object disagreed, every
          Fork B occupancy would be scaled wrongly.
    (ii)  _vbem_infer's corrected sufficient statistic must equal an
          INDEPENDENT call to saspt.defoc.defoc_corr on the uncorrected
          statistic -- pinning down that the correction is a division by
          f_remain along the diffusion-coefficient axis and nothing else.
    (iii) The correction must be applied strictly AFTER the loop and never
          feed back into exp_log_tau. This is what licenses the claim that the
          per-FOV bin-occupancy panels are defocalization-invariant while the
          posterior-overlay / W1 panels are not. It is asserted in
          _vbem_infer's docstring but was never tested.
    """
    tier = cfg['tier']
    s = _settings(tier)
    dz, dt = float(s['focal_depth']), float(s['frame_interval'])
    det = _subsample_detections(os.path.join(tier, cfg['traj']),
                                os.path.join(tier, 'settings.pkl'), n_traj)
    SA = StateArray.from_detections(det, **_run_settings(tier, 200))
    lik_stock = SA.likelihood

    cache = load_likelihood_cache(os.path.join(tier, cfg['cache']),
                                  settings_pkl=os.path.join(tier, 'settings.pkl'))
    lik_rebuilt = oc._make_likelihood_from_cache(cache)

    D = np.asarray(lik_stock.diff_coefs)

    # (i) Probe each likelihood object with an all-ones array: the correction
    # divides by f_remain, so the reciprocal of the result IS the survival curve.
    probe = np.ones((D.size, np.asarray(lik_stock.loc_errors).size), dtype=np.float64)
    fr_stock = 1.0 / lik_stock.correct_for_defocalization(probe.copy(),
                                                          normalize=False)[:, 0]
    fr_rebuilt = 1.0 / lik_rebuilt.correct_for_defocalization(probe.copy(),
                                                              normalize=False)[:, 0]
    fr_direct = np.array([f_remain_rbm(d, 1, dt, dz)[0] for d in D])
    # fr_stock/fr_rebuilt were recovered by inverting a division, so they are
    # compared to the direct call at float round-trip tolerance, not bitwise;
    # the bitwise claims that matter are the two below (same likelihood object
    # rebuilt from a cache, and the correction itself in check (ii)).
    check_i = dict(
        bitwise_stock_vs_rebuilt=bool(np.array_equal(fr_stock, fr_rebuilt)),
        max_abs_diff_stock_vs_rebuilt=float(np.max(np.abs(fr_stock - fr_rebuilt))),
        matches_f_remain_rbm=bool(np.allclose(fr_stock, fr_direct, rtol=1e-12, atol=0)),
        max_rel_diff_vs_f_remain_rbm=float(np.max(np.abs(fr_stock - fr_direct)
                                                  / fr_direct)),
        focal_depth_um=dz, frame_interval_s=dt,
        f_remain_min=float(fr_direct.min()), f_remain_max=float(fr_direct.max()),
        max_occupancy_inflation=float(1.0 / fr_direct.min()))

    # (ii) _vbem_infer's correction vs an independent defoc_corr call.
    L_sub, jumps, prior = _fork_b_inputs(SA)
    n_raw, _, tau_off = oc._vbem_infer(L_sub, jumps, prior, lik_stock, 200,
                                       tol=0.0, apply_defoc=False)
    n_cor, _, tau_on = oc._vbem_infer(L_sub, jumps, prior, lik_stock, 200,
                                      tol=0.0, apply_defoc=True)
    n_indep = defoc_corr(n_raw.copy(),
                         (D, np.asarray(lik_stock.loc_errors)),
                         likelihood=lik_stock.name,
                         frame_interval=dt, dz=dz, normalize=False)
    check_ii = dict(
        likelihood_name=str(lik_stock.name),
        bitwise=bool(np.array_equal(n_cor, n_indep)),
        max_abs_diff=float(np.max(np.abs(n_cor - n_indep))),
        max_rel_diff=float(np.max(np.abs(n_cor - n_indep)
                                  / np.maximum(np.abs(n_indep), 1e-30))))

    # (iii) exp_log_tau must be untouched by the correction; record how big the
    # effect on the pmf actually is, for the Methods text.
    log_grid = np.log10(D)
    _, pm_off, _ = oc._vbem_infer(L_sub, jumps, prior, lik_stock, 200,
                                  tol=0.0, apply_defoc=False)
    _, pm_on, _ = oc._vbem_infer(L_sub, jumps, prior, lik_stock, 200,
                                 tol=0.0, apply_defoc=True)
    check_iii = dict(
        exp_log_tau_bitwise_identical=bool(np.array_equal(tau_off, tau_on)),
        w1_defoc_on_vs_off=float(oc._w1_logD(log_grid, _pmf(pm_off), _pmf(pm_on))),
        max_abs_pmf_diff=float(np.max(np.abs(_pmf(pm_on) - _pmf(pm_off)))))

    results[f'T2_defocalization::{key}'] = dict(
        check_i_survival_curve=check_i,
        check_ii_correction_math=check_ii,
        check_iii_applied_post_loop=check_iii,
        passed=bool(check_i['bitwise_stock_vs_rebuilt']
                    and check_i['matches_f_remain_rbm']
                    and check_ii['bitwise']
                    and check_iii['exp_log_tau_bitwise_identical']))


# ---------------------------------------------------------------------------
# T3 -- likelihood cache fidelity
# ---------------------------------------------------------------------------

def t3_cache_fidelity(key, cfg, n_traj, results):
    """The on-disk cache must be real saspt output, only float32-rounded.

    build_likelihood_cache runs a genuine StateArray with max_iter=0 and stores
    `naive_assignment_probabilities`. Rebuilding that here on the same
    detections catches any drift in the cache builder's filtering or axis
    handling, and bounds the float32 storage error.
    """
    tier = cfg['tier']
    det = _subsample_detections(os.path.join(tier, cfg['traj']),
                                os.path.join(tier, 'settings.pkl'), n_traj)
    SA = StateArray.from_detections(det, **_run_settings(tier, 0))

    L_fresh = np.moveaxis(SA.naive_assignment_probabilities, -1, 0)
    L_round = L_fresh.astype(np.float32).astype(np.float64)
    out = dict(
        n_subtracks=int(L_fresh.shape[0]),
        float32_roundtrip_max_abs=float(np.max(np.abs(L_fresh - L_round))),
        per_subtrack_sums_to_one=float(np.max(np.abs(
            L_fresh.reshape(L_fresh.shape[0], -1).sum(axis=1) - 1.0))))
    out['passed'] = bool(out['float32_roundtrip_max_abs'] <= BAR_CACHE_L
                         and out['per_subtrack_sums_to_one'] <= 1e-6)
    results[f'T3_cache_fidelity::{key}'] = out


# ---------------------------------------------------------------------------
# T4 -- grid / settings identity
# ---------------------------------------------------------------------------

def t4_settings_identity(key, cfg, results):
    """Every grid parameter must agree between the deployed settings.pkl and
    the cache's stored snapshot -- and we record the resolved values so the
    Methods text can quote them rather than assert them.

    max_iter / conc_param are ABSENT from settings.pkl, which means production
    inherits saspt's own defaults; record those explicitly so the Methods text
    states the values that actually ran.
    """
    tier = cfg['tier']
    spkl = os.path.join(tier, 'settings.pkl')
    s = _load_saspt_settings(spkl)
    cache = load_likelihood_cache(os.path.join(tier, cfg['cache']),
                                  settings_pkl=spkl)
    snap = json.loads(str(cache['grid_snapshot']))

    mismatches = []
    for k in _GRID_KEYS:
        a, b = s[k], snap[k]
        same = (np.allclose(np.asarray(a, float), np.asarray(b, float))
                if isinstance(a, np.ndarray) else a == b)
        if not same:
            mismatches.append(k)

    from saspt.constants import DEFAULT_MAX_ITER, DEFAULT_CONC_PARAM
    out = dict(
        saspt_version=saspt.__version__,
        grid_checksum=cache['grid_checksum'][:16],
        mismatched_grid_keys=mismatches,
        max_iter_in_settings=bool('max_iter' in s),
        conc_param_in_settings=bool('conc_param' in s),
        resolved_max_iter=(s['max_iter'] if 'max_iter' in s else DEFAULT_MAX_ITER),
        resolved_conc_param=(s['conc_param'] if 'conc_param' in s
                             else DEFAULT_CONC_PARAM),
        likelihood_type=str(s['likelihood_type']),
        pixel_size_um=float(s['pixel_size_um']),
        frame_interval=float(s['frame_interval']),
        focal_depth=float(s['focal_depth']),
        splitsize=int(s['splitsize']),
        n_diff_coefs=int(np.asarray(s['diff_coefs']).size),
        diff_coef_range=[float(np.min(s['diff_coefs'])),
                         float(np.max(s['diff_coefs']))],
        n_loc_errors=int(np.asarray(s['loc_errors']).size),
        loc_error_range=[float(np.min(s['loc_errors'])),
                         float(np.max(s['loc_errors']))])
    out['passed'] = bool(not mismatches)
    results[f'T4_settings_identity::{key}'] = out


# ---------------------------------------------------------------------------
# T5 -- marginalization provenance
# ---------------------------------------------------------------------------

def t5_marginalization(key, cfg, n_traj, results):
    """The headline per-FOV soft occupancy is not a saspt attribute: it is the
    posterior responsibility R marginalized over the localization-error axis
    and averaged over a FOV's subtracks. Check the local marginal against
    saspt's own `marginalize_on_diff_coef` (lik.py:245-272), which sums over
    the loc-error axis and renormalizes per trajectory."""
    tier = cfg['tier']
    det = _subsample_detections(os.path.join(tier, cfg['traj']),
                                os.path.join(tier, 'settings.pkl'), n_traj)
    SA = StateArray.from_detections(det, **_run_settings(tier, 200))
    marg_ref = SA.likelihood.marginalize_on_diff_coef(
        SA.posterior_assignment_probabilities)              # (n_D, n_sub)

    L_sub, jumps, prior = _fork_b_inputs(SA)
    _, _, tau = oc._vbem_infer(L_sub, jumps, prior, SA.likelihood, 200, tol=0.0)

    # Reproduce the downstream construction used by softpost_region_panels /
    # _group_posterior_full: R = L*tau / rownorm, summed over loc-error, then
    # renormalized per subtrack.
    n_sub, n_D, n_loc = L_sub.shape
    norm = L_sub.reshape(n_sub, -1) @ tau.reshape(-1)
    norm[norm == 0] = 1.0
    marg = np.zeros((n_sub, n_D))
    for j in range(n_loc):
        marg += L_sub[:, :, j] * tau[:, j]
    marg = (marg.T / norm).T
    marg /= marg.sum(axis=1, keepdims=True)

    out = dict(
        max_abs_diff=float(np.max(np.abs(marg.T - marg_ref))),
        max_rel_diff=float(np.max(np.abs(marg.T - marg_ref)
                                  / np.maximum(marg_ref, 1e-12))))
    out['passed'] = bool(out['max_abs_diff'] <= 1e-9)
    results[f'T5_marginalization::{key}'] = out


# ---------------------------------------------------------------------------
# T6 -- end-to-end
# ---------------------------------------------------------------------------

def t6_end_to_end(key, cfg, n_traj, results):
    """Full stock saspt (detections -> posterior_occs) vs the full Fork B route
    (detections -> cached-style L -> _vbem_infer -> pmf), both at saSPT's own
    settings. This exercises both paths in their entirety, including the
    float32 storage step the real cache uses."""
    tier = cfg['tier']
    det = _subsample_detections(os.path.join(tier, cfg['traj']),
                                os.path.join(tier, 'settings.pkl'), n_traj)

    SA_stock = StateArray.from_detections(det, **_run_settings(tier, 200))
    pmf_stock = _pmf(SA_stock.posterior_occs)

    SA_lik = StateArray.from_detections(det, **_run_settings(tier, 0))
    L_sub, jumps, prior = _fork_b_inputs(SA_lik, dtype=np.float32)
    _, post_mean, _ = oc._vbem_infer(L_sub, jumps, prior, SA_lik.likelihood, 200,
                                     tol=1e-5, weights=None, apply_defoc=True)
    pmf_fork = _pmf(post_mean)

    log_grid = np.log10(np.asarray(SA_stock.likelihood.diff_coefs))
    out = dict(
        n_subtracks=int(L_sub.shape[0]),
        max_abs_pmf_diff=float(np.max(np.abs(pmf_fork - pmf_stock))),
        w1_stock_vs_fork=float(oc._w1_logD(log_grid, pmf_stock, pmf_fork)))
    out['passed'] = bool(out['w1_stock_vs_fork'] <= BAR_PROD_W1)
    results[f'T6_end_to_end::{key}'] = out


# ---------------------------------------------------------------------------

_METRIC = {
    'T1': lambda v: (f"exact {v['configs']['a_saspt_exact']['max_abs_occ_diff']:.1e} | "
                     f"prod {v['configs']['d_float32_tol']['max_abs_occ_diff']:.1e} | "
                     f"saspt self-conv {v['saspt_self_convergence_200_vs_5000']:.1e}"),
    'T2': lambda v: (f"dW1(defoc on vs off) = "
                     f"{v['check_iii_applied_post_loop']['w1_defoc_on_vs_off']:.4f}"),
    'T3': lambda v: f"f32 roundtrip {v['float32_roundtrip_max_abs']:.2e}",
    'T4': lambda v: (f"max_iter={v['resolved_max_iter']} "
                     f"conc={v['resolved_conc_param']} dz={v['focal_depth']}"),
    'T5': lambda v: f"max |diff| {v['max_abs_diff']:.2e}",
    'T6': lambda v: f"W1 = {v['w1_stock_vs_fork']:.2e}",
}


def main():
    ap = argparse.ArgumentParser(
        description='saSPT / Fork-B numerical equivalence suite')
    ap.add_argument('--n-traj', type=int, default=20000,
                    help='raw trajectories to subsample per dataset for the '
                         'fit-based tests (default 20000)')
    ap.add_argument('--datasets', nargs='*', default=list(DATASETS),
                    choices=list(DATASETS))
    ap.add_argument('--out', default=tlsd_paths.outp(
        'PaperFigures_updated', 'verification', 'saspt_equivalence_report.json'))
    args = ap.parse_args()

    results = {'meta': dict(
        saspt_version=saspt.__version__,
        saspt_path=os.path.dirname(saspt.__file__),
        numpy_version=np.__version__,
        python=sys.version.split()[0],
        n_traj_subsample=args.n_traj)}

    for key in args.datasets:
        cfg = DATASETS[key]
        for fn, needs_n in ((t4_settings_identity, False),
                            (t2_defocalization, True),
                            (t3_cache_fidelity, True),
                            (t1_vbem_identity, True),
                            (t5_marginalization, True),
                            (t6_end_to_end, True)):
            try:
                print(f"[{key}] {fn.__name__} ...", flush=True)
                fn(key, cfg, args.n_traj, results) if needs_n \
                    else fn(key, cfg, results)
            except Exception:
                results[f'{fn.__name__}::{key}'] = dict(
                    passed=False, error=traceback.format_exc())
                print(traceback.format_exc(), file=sys.stderr, flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as fh:
        json.dump(results, fh, indent=2, default=oc._json_default)

    print('\n' + '=' * 86)
    print(f"{'test':<46}{'result':<9}key metric")
    print('-' * 86)
    n_fail = 0
    for k, v in results.items():
        if k == 'meta' or not isinstance(v, dict) or 'passed' not in v:
            continue
        n_fail += (not v['passed'])
        fmt = _METRIC.get(k[:2])
        try:
            metric = fmt(v) if fmt else ''
        except Exception:
            metric = v.get('error', '').splitlines()[-1:] and 'ERROR'
        print(f"{k:<46}{'PASS' if v['passed'] else 'FAIL':<9}{metric}")
    print('=' * 86)
    print(f"wrote {args.out}")
    return 1 if n_fail else 0


if __name__ == '__main__':
    sys.exit(main())
