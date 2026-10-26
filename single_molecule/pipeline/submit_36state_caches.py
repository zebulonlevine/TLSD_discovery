"""
submit_36state_caches.py -- submit one cache-build job per condition, each with
its memory request computed from that condition's sub-trajectory count.

A single array-wide request would have to cover the largest condition, leaving
the small ones at ~20% utilisation. Sizing per condition keeps every job near
the 70% target.

Calibration: a measured build of 83,158 sub-trajectories at 36 loc-error states
peaked at 11.1 GB, i.e. ~6.2x the stored float32 L array (saSPT builds in
float64 and takes an exp() copy). Request = 1.4x that estimate, floored at 8 GB.
"""
import argparse
import os
import subprocess
import sys


CODE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODE)
from build_variant_cache_one import DATASETS

# Sub-trajectories per condition, read from the 36-state cache headers (the
# same numbers as NSUB in figures/submit_panel_refit.py). Used only for sizing.
NSUB = {
    'sgNT_60S': 83158, 'sgUBA5_60S': 109543,
    'sgNT_40S': 166062, 'sgUBA5_40S': 91480,
    'S3_NT_none': 204998, 'S3_NT_HT10': 168080, 'S3_NT_eme360': 134970,
    'S3_UBA5_none': 271998, 'S3_UBA5_HT10': 132139, 'S3_UBA5_eme360': 264176,
    'noDegron': 115656, 'Degron': 162860,
}
PEAK_PER_L = 6.2       # measured ratio of peak RSS to stored float32 L size
HEADROOM = 1.4         # -> ~71% utilisation


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--datasets', nargs='*', default=list(DATASETS))
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    print(f"{'dataset':8} {'condition':16} {'n_sub':>9} {'L GB':>6} "
          f"{'est peak':>9} {'--mem':>7}")
    total = 0
    for ds in args.datasets:
        for cond in DATASETS[ds]['conditions']:
            n = NSUB[cond]
            l_gb = n * 150 * 36 * 4 / 1e9
            mem = max(8, int(l_gb * PEAK_PER_L * HEADROOM + 1))
            print(f'{ds:8} {cond:16} {n:>9,} {l_gb:>6.2f} '
                  f'{l_gb*PEAK_PER_L:>8.1f}G {mem:>6}G')
            total += mem
            cmd = ['sbatch', f'--mem={mem}G', '--parsable',
                   '--export', f'ALL,DATASET={ds},COND={cond},VARIANT=36state',
                   os.path.join(CODE, 'call_variant_cache.slurm')]
            if args.dry_run:
                continue
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                print(f'   sbatch FAILED: {r.stderr.strip()}')
            else:
                print(f'   job {r.stdout.strip()}')
    print(f'\n{"DRY RUN -- nothing submitted" if args.dry_run else "submitted"}; '
          f'peak aggregate memory if all ran at once: {total} GB '
          f'(they will not -- the partition schedules them)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
