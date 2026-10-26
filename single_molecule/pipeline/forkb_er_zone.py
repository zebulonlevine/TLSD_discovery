"""
forkb_er_zone.py
==================

Fork-B-style (VBEM-refit) zone comparator for the ER on-ER-vs-off-ER
contrast (Fig 5F). Replaces er_zone_naive_cache.py's naive/
no-refit approach per Zeb's correction (Fork B is the default statistic,
not naive occupancy) -- and independently fixes a real calibration bug
V4's negative control caught in that module: zone labels must be permuted
at the TRAJECTORY level (all of a trajectory's subtracks move together),
never at the subtrack level, which is anticonservative (CLAUDE.md already
documents this exact failure mode for FOV-level permutation; the same
principle applies one level down here).

Generalizes occupancy_compare._group_posterior_full (which selects
subtracks by FOV membership) to select by an arbitrary subtrack-level
boolean mask instead -- zone membership is a per-trajectory property, not
an FOV-level one, so the FOV-based selector doesn't apply here.
"""

import json
import os
import time

import numpy as np
import pandas as pd

from occupancy_compare import (
    _vbem_infer, _make_likelihood_from_cache, _w1_logD, _subband_summary,
    plot_region_effect_vs_null, load_likelihood_cache,
)


def _group_posterior_full_mask(cache, likelihood, subtrack_mask, conc_param, max_iter,
                               equal_weight=False, apply_defoc=True):
    """VBEM-refit posterior for an arbitrary subtrack-level boolean mask
    (not necessarily FOV-aligned). Returns (pmf, occ_by_fov, fovs, n_sub_by_fov):
    pmf = population-level soft posterior; occ_by_fov = per-FOV posterior-
    MLE_D histogram over the FOVs actually present in this mask (real
    per-FOV replicates, for bootstrap CIs)."""
    L_sub = cache['L'][subtrack_mask]
    jumps_sub = cache['jumps_per_track'][subtrack_mask].astype(np.float32)
    prior = np.full(L_sub.shape[1:], conc_param, dtype=np.float32)
    w = np.ones_like(jumps_sub) if equal_weight else None
    _, post_mean, exp_log_tau_2d = _vbem_infer(L_sub, jumps_sub, prior, likelihood,
                                               max_iter, weights=w,
                                               apply_defoc=apply_defoc)
    pmf = post_mean.sum(axis=1)

    n_sub_local, n_D, n_loc = L_sub.shape
    L_flat = L_sub.reshape(n_sub_local, n_D * n_loc)
    norm = L_flat @ exp_log_tau_2d.reshape(-1)
    norm[norm == 0] = 1.0
    marginal = np.zeros((n_sub_local, n_D), dtype=np.float32)
    for j in range(n_loc):
        marginal += L_sub[:, :, j] * exp_log_tau_2d[None, :, j]
    marginal /= norm[:, None]
    mle_d_idx = np.argmax(marginal, axis=1)
    one_hot = np.zeros((n_sub_local, n_D), dtype=np.float64)
    one_hot[np.arange(n_sub_local), mle_d_idx] = 1.0

    fov_local = cache['fov_of_subtrack'][subtrack_mask]
    fovs = np.unique(fov_local)
    occ = np.full((len(fovs), n_D), np.nan)
    n_sub = np.zeros(len(fovs), dtype=int)
    for i, fov in enumerate(fovs):
        rows = fov_local == fov
        n = int(rows.sum())
        n_sub[i] = n
        if n > 0:
            occ[i] = one_hot[rows].mean(axis=0)
    return pmf, occ, fovs, n_sub


def zone_refit_permutation(cache_path, traj_zone, zone_a, zone_b, name_a, name_b,
                            out_dir, n_perm=100, conc_param=1.0, max_iter=200,
                            seed=0, verbose=True, force=False, cond_label='',
                            equal_weight=False, apply_defoc=True):
    """Fork-B (VBEM-refit) zone comparator with a within-FOV,
    TRAJECTORY-level zone-relabeling permutation null (fixes the subtrack-
    level anticonservative bug in er_zone_naive_cache.zone_within_fov_w1).

    traj_zone: dict {raw trajectory id (int) -> zone label str}, as built
    by er_zone_naive_cache.build_zone_labeled_cache (reused unchanged --
    only the statistical engine below is new).

    Idempotent: if the output json+npz already exist and not `force`, the
    (expensive) VBEM refit + permutation is skipped and the cached result
    dict + npz path are returned -- so a re-run only redoes the cheap
    downstream plotting.
    """
    os.makedirs(out_dir, exist_ok=True)
    tag = f'{cond_label}_{name_a}_vs_{name_b}' if cond_label else f'{name_a}_vs_{name_b}'
    if equal_weight or not apply_defoc:
        tag += '_eqfit'
    json_path = os.path.join(out_dir, f'zone_refit_permutation_{tag}.json')
    npz_path_cached = os.path.join(out_dir, f'zone_refit_permutation_{tag}.npz')
    if not force and os.path.exists(json_path) and os.path.exists(npz_path_cached):
        if verbose:
            print(f'[{tag}] cached -> reusing {npz_path_cached} (skipping VBEM refit)')
        with open(json_path) as f:
            return json.load(f), npz_path_cached

    cache = load_likelihood_cache(cache_path)
    likelihood = _make_likelihood_from_cache(cache)
    diff_coefs = cache['diff_coefs']
    log_grid = np.log10(diff_coefs)

    traj_id = np.asarray(cache['orig_trajectory']).astype(int)
    zone = np.array([traj_zone.get(t, 'unassigned') for t in traj_id])
    fov_of_subtrack = cache['fov_of_subtrack']

    keep = np.isin(zone, [zone_a, zone_b])
    traj_id, zone, fov_of_subtrack = traj_id[keep], zone[keep], fov_of_subtrack[keep]
    L_all = cache['L'][keep]
    jumps_all = cache['jumps_per_track'][keep]
    sub_cache = {'L': L_all, 'jumps_per_track': jumps_all,
                 'fov_of_subtrack': fov_of_subtrack}

    # ---- observed ----
    t0 = time.time()
    mask_a_obs = zone == zone_a
    mask_b_obs = zone == zone_b
    pmf_a, occ_a_mled, fovs_a, n_sub_a = _group_posterior_full_mask(
        sub_cache, likelihood, mask_a_obs, conc_param, max_iter,
        equal_weight=equal_weight, apply_defoc=apply_defoc)
    pmf_b, occ_b_mled, fovs_b, n_sub_b = _group_posterior_full_mask(
        sub_cache, likelihood, mask_b_obs, conc_param, max_iter,
        equal_weight=equal_weight, apply_defoc=apply_defoc)
    w1_obs = _w1_logD(log_grid, pmf_a, pmf_b)
    if verbose:
        print(f'[{name_a} vs {name_b}] observed fits: {time.time()-t0:.1f}s  '
              f'W1={w1_obs:.4f}  n_fov_a={len(fovs_a)} n_fov_b={len(fovs_b)}')

    # trajectory-level table (one row per trajectory, not per subtrack) --
    # the permutation unit. A trajectory's zone label is shuffled as one
    # atomic unit, never split across its own subtracks.
    traj_table = (pd.DataFrame({'traj_id': traj_id, 'fov': fov_of_subtrack, 'zone': zone})
                    .drop_duplicates('traj_id')
                    .reset_index(drop=True))
    fov_groups = {f: g.index.values for f, g in traj_table.groupby('fov')}
    z0 = traj_table['zone'].values

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    null_curves = np.empty((n_perm, len(diff_coefs)))
    t0 = time.time()
    for k in range(n_perm):
        zp = z0.copy()
        for f, idx in fov_groups.items():
            zp[idx] = rng.permutation(zp[idx])
        traj_to_pseudo = dict(zip(traj_table['traj_id'].values, zp))
        pseudo_zone = pd.Series(traj_id).map(traj_to_pseudo).values
        mask_a = pseudo_zone == zone_a
        mask_b = pseudo_zone == zone_b
        pa, _, _, _ = _group_posterior_full_mask(sub_cache, likelihood, mask_a, conc_param, max_iter,
                                       equal_weight=equal_weight, apply_defoc=apply_defoc)
        pb, _, _, _ = _group_posterior_full_mask(sub_cache, likelihood, mask_b, conc_param, max_iter,
                                       equal_weight=equal_weight, apply_defoc=apply_defoc)
        null[k] = _w1_logD(log_grid, pa, pb)
        null_curves[k] = pb - pa
        if verbose and (k + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f'  permutation {k+1}/{n_perm}  ({elapsed/(k+1):.1f}s/perm avg)')

    null_p95 = float(np.percentile(null, 95))
    p_value = float((null >= w1_obs).mean())
    verdict = 'separation IS above the relabeling null' if p_value < 0.05 else 'within the relabeling null'

    tag = f'{cond_label}_{name_a}_vs_{name_b}' if cond_label else f'{name_a}_vs_{name_b}'
    result = {
        'name_a': name_a, 'name_b': name_b, 'w1_observed': w1_obs,
        'null_p95': null_p95, 'null_mean': float(null.mean()), 'p_value': p_value,
        'verdict': verdict, 'n_perm': n_perm, 'seed': seed,
        'n_fov_a': len(fovs_a), 'n_fov_b': len(fovs_b),
        'n_sub_a': int(mask_a_obs.sum()), 'n_sub_b': int(mask_b_obs.sum()),
    }
    with open(os.path.join(out_dir, f'zone_refit_permutation_{tag}.json'), 'w') as f:
        json.dump(result, f, indent=2)
    np.savez_compressed(
        os.path.join(out_dir, f'zone_refit_permutation_{tag}.npz'),
        diff_coefs=diff_coefs, pmf_a=pmf_a, pmf_b=pmf_b, null=null,
        null_curves=null_curves, occ_a_mled=occ_a_mled, occ_b_mled=occ_b_mled,
        fovs_a=fovs_a, fovs_b=fovs_b, **{k: v for k, v in result.items()
                                          if k not in ('name_a', 'name_b', 'verdict')})

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(null, bins=30, color='gray', alpha=0.6, label=f'trajectory-level within-FOV null (n_perm={n_perm})')
    ax.axvline(null_p95, color='k', ls='--', lw=1.2, label=f'null p95={null_p95:.4f}')
    ax.axvline(w1_obs, color='#cc3333', lw=2, label=f'observed W1={w1_obs:.4f}')
    ax.set_xlabel(r'W1(log$_{10}$ D)')
    ax.set_title(f'Fork B zone comparator: {name_a} vs {name_b} (p={p_value:.4f})')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'zone_refit_null_distribution_{tag}.pdf'))
    plt.close(fig)

    return result, os.path.join(out_dir, f'zone_refit_permutation_{tag}.npz')
