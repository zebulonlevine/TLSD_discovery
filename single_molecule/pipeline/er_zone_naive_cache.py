"""
er_zone_naive_cache.py
========================

Zone-mapping utility for the ER dataset's on-ER-vs-off-ER contrast
(Fig 5F). Reuses saspt_likelihood_cache.build_likelihood_cache
on the all-zone (unfiltered) combinationwise pool, then derives each
trajectory's zone label by joining against the zone-specific combos
(*_on_structure, *_on_ER_broad, *_off_ER) on (fov_basename, x, y, frame) --
the same cross-reference join strategy already validated elsewhere in this
project (CLAUDE.md: "verified by (x,y,frame) identity matched 100%").
Verified here too: 100% match, zero on_structure/off_ER overlap.

WHY a new path instead of merge_likelihood_caches / CombinedFit:
  - merge_likelihood_caches hard-rejects the same FOV appearing in more
    than one member -- correct for its intended use (distinct conditions),
    but ER zone combos deliberately violate that (one physical FOV
    contributes different trajectories to different zone pools).
  - CombinedFit needs a subtrack_posteriors/*.npz, which is empty for this
    dataset (never populated -- optional output, opt-in emission never
    wired up here).

So: ONE naive cache from the all-zone pool (no FOV collision, since it's a
single unfiltered pool), zone assigned to each of ITS subtracks via the
join above. build_zone_labeled_cache is reused by forkb_er_zone.py, which
holds the actual Fork-B (VBEM-refit) statistical engine for this contrast.

This module used to ALSO contain a naive (no-refit) statistical test,
zone_within_fov_w1 -- removed after V4's negative control caught it
anticonservative (permuted zone labels at the SUBTRACK level; a trajectory's
own subtracks could get split across pseudo-zones, which is exactly the
"subtrack-level shuffling is anticonservative" failure mode CLAUDE.md
already documents for FOV-level permutation, one level down). The fix
(permute at the TRAJECTORY level) lives in forkb_er_zone.zone_refit_permutation,
which is also VBEM-refit-based per the project's "Fork B is the default,
not naive occupancy" standard -- so it wasn't reintroduced here fixed-but-naive.
"""

import os

import pandas as pd

import saspt_likelihood_cache as slc


def _traj_zone_map(root, combo_all, zone_combo, key_round=3):
    """Set of (fov_basename, x, y, frame) tuples for one zone combo, keyed
    to the SAME rounding as the all-pool join."""
    traj = pd.read_csv(os.path.join(root, 'traj_combinations', f'{zone_combo}_pooled_traj.csv'))
    fov = pd.read_csv(os.path.join(root, 'traj_combinations', f'{zone_combo}_fov_map.csv'))
    m = traj.merge(fov[['trajectory', 'fov_basename']], on='trajectory', how='left')
    keys = set(map(tuple, m[['fov_basename', 'x', 'y', 'frame']].round(key_round).values))
    return keys


def build_zone_labeled_cache(root, settings_pkl, combo_all, zones, out_dir,
                              cache_name=None, force=False, verbose=True):
    """Builds (or loads) the all-pool naive cache, and returns
    (cache, trajectory_zone) where trajectory_zone is a dict mapping the
    cache's own local trajectory id (matching parse_pooled_id(orig_trajectory))
    -> zone label string (one of `zones`, or 'unassigned' if it matched no
    zone combo -- e.g. ambiguous/intermediate trajectories per CLAUDE.md's
    5-label ER scheme, which this contrast doesn't use).
    """
    cache_name = cache_name or f'_standalone_{combo_all}'
    traj_csv = os.path.join(root, 'traj_combinations', f'{combo_all}_pooled_traj.csv')
    fov_csv = os.path.join(root, 'traj_combinations', f'{combo_all}_fov_map.csv')
    cache_path = slc.build_likelihood_cache(
        traj_csv, settings_pkl, out_dir, fov_map_csv=fov_csv,
        cache_name=cache_name, force=force, verbose=verbose)
    cache = slc.load_likelihood_cache(cache_path)

    all_traj = pd.read_csv(traj_csv)
    all_fov = pd.read_csv(fov_csv)
    all_m = all_traj.merge(all_fov[['trajectory', 'fov_basename']], on='trajectory', how='left')
    all_m['key'] = list(map(tuple, all_m[['fov_basename', 'x', 'y', 'frame']].round(3).values))

    zone_of_key = {}
    for zone in zones:
        keys = _traj_zone_map(root, combo_all, f'{combo_all}_{zone}')
        for k in keys:
            zone_of_key[k] = zone
        if verbose:
            print(f'  zone={zone}: {len(keys):,} localizations')

    all_m['zone'] = all_m['key'].map(zone_of_key).fillna('unassigned')
    traj_zone = all_m.groupby('trajectory')['zone'].agg(
        lambda s: s.value_counts().idxmax()).to_dict()

    if verbose:
        vc = pd.Series(list(traj_zone.values())).value_counts()
        print(f'  trajectory-level zone assignment: {vc.to_dict()}')

    return cache, traj_zone
