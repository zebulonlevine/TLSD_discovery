"""
build_variant_cache_one.py -- build ONE per-condition likelihood cache under a
named inference variant, into its own cache directory so existing caches are
untouched.

Variant
-------
36state      saSPT's default 36-state localization-error grid (0-0.070 um),
             with the localization-error axis marginalized as saSPT
             intends. This is the configuration every panel uses
             (PANEL_VARIANT=36state).

Why per-condition rather than the merged multi-condition caches used before:
at 36 localization-error states the merged claim-2 cache would be ~25 GB
stored and ~50 GB during saSPT's float64 build. Nothing is lost -- every
downstream statistic fits each condition over its own sub-trajectories only
(panel_core.per_fov_occupancy), so a merged cache buys only a shared-grid
guarantee that the stored checksum already provides.

Usage (one condition per SLURM job, memory sized per condition):
    python build_variant_cache_one.py --dataset claim1 --condition sgNT_60S --variant 36state
"""
import argparse
import os
import pickle
import sys

import numpy as np

CODE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODE)
from saspt_likelihood_cache import build_likelihood_cache

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
A = tlsd_paths.spt('all_data')
FOCAL_DEPTH = 0.7

VARIANTS = {
    '36state':    lambda ds: np.arange(0.0, 0.072, 0.002),
}

DATASETS = {
    'claim1': dict(tier=f'{A}/2025_06_06_zeb/tracking_output_q=15p0_splitsize10',
                   tc='traj_combinations',
                   conditions=['sgNT_60S', 'sgUBA5_60S', 'sgNT_40S', 'sgUBA5_40S']),
    'claim2': dict(tier=f'{A}/2026_07_08/tracking_output_q=15p0',
                   tc='traj_combinations',
                   conditions=['S3_NT_none', 'S3_NT_HT10', 'S3_NT_eme360',
                               'S3_UBA5_none', 'S3_UBA5_HT10', 'S3_UBA5_eme360']),
    'claim3': dict(tier=f'{A}/2025_09_15_zeb_ER/tracking_output_q=15p0',
                   tc='traj_combinations',
                   conditions=['noDegron', 'Degron']),
}


def variant_settings(tier, out_dir, variant, dataset):
    """A settings.pkl carrying the variant's grid and focal_depth, written into
    the variant cache dir so load_likelihood_cache's checksum assertion has the
    matching reference and cannot pair a cache with the wrong grid."""
    os.makedirs(out_dir, exist_ok=True)
    p = os.path.join(out_dir, f'settings_{variant}.pkl')
    grid = VARIANTS[variant](dataset)
    if grid is None or (hasattr(grid, '__len__') and None in np.atleast_1d(grid)):
        raise SystemExit(f'no localization error defined for {dataset}/{variant}')
    if not os.path.exists(p):
        s = pickle.load(open(os.path.join(tier, 'settings.pkl'), 'rb'))
        s['saspt']['loc_errors'] = grid
        s['saspt']['focal_depth'] = FOCAL_DEPTH
        s['saspt']['progress_bar'] = False
        pickle.dump(s, open(p, 'wb'))
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dataset', required=True, choices=list(DATASETS))
    ap.add_argument('--condition', required=True)
    ap.add_argument('--variant', required=True, choices=list(VARIANTS))
    args = ap.parse_args()
    cfg = DATASETS[args.dataset]
    tier, cond = cfg['tier'], args.condition
    out_dir = os.path.join(tier, f'occupancy_compare_cache_{args.variant}')
    spkl = variant_settings(tier, out_dir, args.variant, args.dataset)
    tc = os.path.join(tier, cfg['tc'])
    grid = VARIANTS[args.variant](args.dataset)
    print(f'{args.dataset}/{cond}: variant {args.variant}, '
          f'{len(grid)} loc-error state(s) '
          f'{np.min(grid):.4f}-{np.max(grid):.4f} um, focal_depth {FOCAL_DEPTH}')
    p = build_likelihood_cache(
        os.path.join(tc, f'{cond}_pooled_traj.csv'), spkl, out_dir,
        fov_map_csv=os.path.join(tc, f'{cond}_fov_map.csv'),
        cache_name=cond, force=True, verbose=True)
    import resource
    print(f'  wrote {p}')
    print(f'  peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.2f} GB')
    return 0


if __name__ == '__main__':
    sys.exit(main())
