#!/usr/bin/env python
"""
run_zone_refit_batch.py
=========================

Resampling for Fig 5F's zone contrasts, one condition per task.

Two schemes, and they cost wildly different amounts because of how the panel's
estimator is built:

  permutation  Zone labels are shuffled among TRAJECTORIES within each field
               of view (all of a trajectory's sub-trajectories move together).
               per_fov_occ_by_zone fits ONE VBEM over the pooled
               sub-trajectories of both zones and then partitions the
               marginals, so relabelling zones leaves the pooled set -- and
               therefore exp_log_tau and every marginal -- bit-identical. The
               fit is invariant under this null, so it is computed once and
               every permutation is just a re-partition. No refit is needed
               and the draw count can be large enough to escape the p floor.
               (forkb_er_zone.zone_refit_permutation DOES refit per draw, but
               only because it fits each zone separately -- a different
               estimator from the one 5F plots.)

  bootstrap    Fields of view are resampled with replacement, which DOES
               change the pooled set, so each draw needs a real refit. This is
               the expensive half, and the one that propagates fit
               uncertainty into the interval.

Writes per-FOV zone differences per draw, so the zone x treatment interaction
can be formed downstream from the two conditions' draws.
"""
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
from occupancy_compare import _make_likelihood_from_cache
import er_zone_naive_cache as ezc
import build_claim3_panels as c3

ZONE_ON, ZONE_OFF = c3.ZONE_ON, c3.ZONE_OFF


def pooled_fit(ctx, L_sub, jumps):
    prior = np.full(L_sub.shape[1:], pc.CONC_PARAM, dtype=np.float32)
    w = np.ones_like(jumps) if pc.EQUAL_WEIGHT else None
    from occupancy_compare import _vbem_infer
    _, _, tau = _vbem_infer(L_sub, jumps, prior, ctx['likelihood'],
                            pc._VBEM_ITERS, tol=pc._VBEM_TOL, weights=w,
                            apply_defoc=False)
    return tau


def zone_diff_from_partition(marg, fov, zsub, fov_slots):
    """Per-FOV (on - off) occupancy difference, one row per FOV slot."""
    out = np.full((len(fov_slots), marg.shape[1]), np.nan)
    for i, f in enumerate(fov_slots):
        in_f = fov == f
        on, off = in_f & (zsub == ZONE_ON), in_f & (zsub == ZONE_OFF)
        if on.any() and off.any():
            out[i] = marg[on].mean(axis=0) - marg[off].mean(axis=0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--condition', required=True, choices=['noDegron', 'Degron'])
    ap.add_argument('--scheme', required=True, choices=['permutation', 'bootstrap'])
    ap.add_argument('--n-draws', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out-dir', required=True)
    # Batching: the bootstrap refits per draw and is far too slow to run 500
    # serially (measured 173-204 s/draw, i.e. 24-28 h against any sane
    # walltime), so both schemes are sliced across array tasks. Draw k is a
    # pure function of (seed, k), so slicing changes wall-clock only.
    ap.add_argument('--batch-index', type=int, default=0)
    ap.add_argument('--n-batches', type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    bs = -(-a.n_draws // a.n_batches)
    k_lo, k_hi = a.batch_index * bs, min((a.batch_index + 1) * bs, a.n_draws)
    if k_lo >= a.n_draws:
        print(f'batch {a.batch_index}: nothing to do'); return
    out_path = os.path.join(
        a.out_dir,
        f'zone_{a.condition}_{a.scheme}_{a.batch_index:03d}.npz')
    if os.path.exists(out_path):
        print(f'already complete: {out_path}')
        return

    cache, traj_zone = ezc.build_zone_labeled_cache(
        c3.ROOT, c3.SETTINGS_PKL, a.condition, [ZONE_ON, ZONE_OFF],
        c3.CACHE_DIR, cache_name=c3.CACHE_NAME(a.condition), verbose=False)
    zone = c3.zone_masks(cache, traj_zone)
    keep = np.isin(zone, [ZONE_ON, ZONE_OFF])
    L = cache['L'][keep]
    jumps = np.asarray(cache['jumps_per_track'])[keep].astype(np.float32)
    fov = cache['fov_of_subtrack'][keep]
    traj = np.asarray(cache['orig_trajectory']).astype(int)[keep]
    zsub = zone[keep]
    ctx = dict(likelihood=_make_likelihood_from_cache(cache),
               diff_coefs=cache['diff_coefs'])
    fovs = sorted(set(fov))
    print(f'[{a.condition}] {L.shape[0]:,} sub-trajectories in ON/OFF, '
          f'{len(fovs)} FOVs, scheme={a.scheme}', flush=True)

    t0 = time.time()
    if a.scheme == 'permutation':
        # one fit, reused by every draw -- see module docstring
        tau = pooled_fit(ctx, L, jumps)
        marg = pc._marginals(L, tau)
        obs = zone_diff_from_partition(marg, fov, zsub, fovs)
        # trajectory is the permutation unit: one label per trajectory,
        # shuffled only among trajectories sharing a field of view
        # One label per TRAJECTORY, permuted only among trajectories sharing a
        # field of view, then expanded to sub-trajectories through the inverse
        # index. Assigning through a per-trajectory boolean mask instead makes
        # each draw O(n_traj * n_subtracks) -- ~1e9 ops here.
        _, first, inv = np.unique(traj, return_index=True, return_inverse=True)
        tr_zone, tr_fov = zsub[first], fov[first]
        fov_groups = [np.flatnonzero(tr_fov == f) for f in fovs]
        seeds = np.random.SeedSequence(a.seed).spawn(a.n_draws)   # hoisted
        draws = np.empty((k_hi - k_lo, len(fovs), marg.shape[1]))
        for n, k in enumerate(range(k_lo, k_hi)):
            rng = np.random.default_rng(seeds[k])
            tz = tr_zone.copy()
            for g in fov_groups:
                tz[g] = rng.permutation(tz[g])
            draws[n] = zone_diff_from_partition(marg, fov, tz[inv], fovs)
            if n == 0:
                print(f'  first draw {time.time()-t0:.1f}s', flush=True)
    else:
        obs_tau = pooled_fit(ctx, L, jumps)
        obs = zone_diff_from_partition(pc._marginals(L, obs_tau), fov, zsub, fovs)
        rows_by_fov = {f: np.flatnonzero(fov == f) for f in fovs}
        seeds = np.random.SeedSequence(a.seed).spawn(a.n_draws)   # hoisted
        draws = np.empty((k_hi - k_lo, len(fovs), len(ctx['diff_coefs'])))
        for n, k in enumerate(range(k_lo, k_hi)):
            rng = np.random.default_rng(seeds[k])
            pick = [fovs[i] for i in rng.integers(0, len(fovs), len(fovs))]
            blocks = [rows_by_fov[f] for f in pick]
            idx = np.concatenate(blocks)
            tau = pooled_fit(ctx, L[idx], jumps[idx])
            marg = pc._marginals(L[idx], tau)
            fk, zk = fov[idx].copy(), zsub[idx]
            # relabel each draw slot uniquely so a FOV drawn twice counts twice
            off = 0
            for i, blk in enumerate(blocks):
                fk[off:off + len(blk)] = f'__slot{i}'
                off += len(blk)
            draws[n] = zone_diff_from_partition(
                marg, fk, zk, [f'__slot{i}' for i in range(len(blocks))])
            if n == 0:
                print(f'  first draw {time.time()-t0:.1f}s', flush=True)

    el = time.time() - t0
    np.savez_compressed(out_path, obs=obs, draws=draws,
                        k_lo=k_lo, k_hi=k_hi,
                        diff_coefs=ctx['diff_coefs'], fovs=np.array(fovs, dtype=object),
                        n_draws=a.n_draws, seed=a.seed, elapsed_sec=el)
    nd = k_hi - k_lo
    print(f'{nd} draws [{k_lo},{k_hi}) in {el:.0f}s '
          f'({el/nd:.2f}s/draw) -> {out_path}')


if __name__ == '__main__':
    main()
