#!/usr/bin/env python
"""
submit_panel_refit_all.py
===========================

Submit every remaining contrast's refit resampling as ONE throttled array job.

Why unified rather than one array per contrast: the cost here is memory
bandwidth, not CPU. A draw takes 121 s alone on a node and 520 s at ~12 array
tasks per node -- a 4.3x penalty -- so saturating the per-user CPU cap costs
about 4x the core-hours AND finishes later. Bounding total concurrency is
therefore the single most important knob, and separate array jobs cannot share
a `%throttle`. One array over (contrast x scheme x batch) can.

    python submit_panel_refit_all.py --dry-run
    python submit_panel_refit_all.py --throttle 32

Skips contrasts that already have aggregated stats, so it is safe to re-run.
"""
import argparse, json, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from submit_panel_refit import (CONTRASTS, TIERS, NSUB, OUT_ROOT, CODE, FIG,
                                SEEDS, SEC_PER_SUBTRACK, THREADS, sizing, sbatch)

# Re-measured after the einsum change, alone on a node (see
# diagnostics/check_einsum_regression.py): 121 s/draw at 192,701 pooled
# sub-trajectories.
SEC_PER_SUBTRACK_UNCONTENDED = 121.1 / 192701
PY = sys.executable


# Figure order, not alphabetical: tasks run in manifest order, so this is
# what decides which panel's statistics are finishable first. The Fig S7A/S7C
# 60S contrast (figS7AC_60S) is already complete; 5F needs the trajectory-level
# zone machinery and is queued separately.
PRIORITY = [
    'figS7A_40S',                              # Fig S7A (60S done)
    'figS7D_NT_HHT', 'figS7D_NT_eme',          # Fig S7D/E
    'figS7D_UBA5_HHT', 'figS7D_UBA5_eme',
    'figS7D_geno_none', 'figS7D_geno_HHT', 'figS7D_geno_eme',
]


def seed_for(key):
    """Stable per-contrast seed.

    A shared seed makes two contrasts' resamples correlated whenever their
    group sizes match: SeedSequence(seed).spawn(n)[k] hands both the same RNG
    state, so rng.integers(0, 20, 20) returns the identical index pattern.
    Each contrast is individually valid, but any statistic formed ACROSS
    contrasts -- a difference-in-differences built from two arms, say --
    inherits a correlation that is an artefact of the RNG, not of the data.

    The seeds are fixed in submit_panel_refit.SEEDS: each is the sha256 of the
    contrast's original key, fixed when the deposited draws were made, so
    renaming a contrast does not change its draws.
    """
    return SEEDS[key]


def already_done(key):
    return os.path.exists(os.path.join(OUT_ROOT, key, f'{key}_refit_stats.csv'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-draws', type=int, default=500)
    ap.add_argument('--batch', type=int, default=10)
    ap.add_argument('--throttle', type=int, default=32,
                    help='max concurrently RUNNING array tasks (default 32; '
                         'the measured knee is well below the CPU cap)')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--include-done', action='store_true')
    ap.add_argument('--merge-job', help='reuse an in-flight merge job instead '
                                        'of resubmitting the merges')
    # Partition LIST, most-idle first. This work is memory-bandwidth bound, so
    # the number of distinct NODES the tasks spread over matters more than the
    # per-user CPU quota: 24 currently has 13 fully idle nodes (928 idle CPUs)
    # and a 240 CPU / 2T cap, against 160 / 1T on 20. `weissman` is only 3
    # nodes / 200 CPUs and is often fully allocated, but costs nothing to list.
    ap.add_argument('--partition', default='24,20,weissman')
    args = ap.parse_args()

    ordered = PRIORITY + [k for k in sorted(CONTRASTS) if k not in PRIORITY]
    keys = [k for k in ordered
            if args.include_done or not already_done(k)]
    if not keys:
        print('nothing to do -- every contrast already has aggregated stats')
        return

    tasks, merges, core_h, max_mem = [], [], 0.0, 0
    for key in keys:
        tier_key, ca, cb, label = CONTRASTS[key]
        tier, out_dir = TIERS[tier_key], os.path.join(OUT_ROOT, key)
        pooled = NSUB[ca] + NSUB[cb]
        mem, _, n_batches = sizing(ca, cb, args.n_draws, args.batch)
        max_mem = max(max_mem, mem)
        core_h += pooled * SEC_PER_SUBTRACK_UNCONTENDED * args.n_draws * 2 * THREADS / 3600
        merges.append((tier, ca, cb))
        merged_rel = os.path.join('occupancy_compare_cache_36state',
                                  f'_joint_{ca}_vs_{cb}_likelihood_cache.npz')
        os.makedirs(out_dir, exist_ok=True)
        for scheme in ('permutation', 'bootstrap'):
            spec_file = os.path.join(out_dir, f'_spec_{scheme}.json')
            if not args.dry_run:
                with open(spec_file, 'w') as f:
                    json.dump(dict(tier=tier, cache_rel=merged_rel, cond_a=ca,
                                   cond_b=cb, label=label, scheme=scheme,
                                   n_draws=args.n_draws,
                                   seed=seed_for(key),
                                   out_dir=out_dir, key=key), f, indent=2)
            for b in range(n_batches):
                tasks.append(dict(key=key, scheme=scheme, batch_index=b,
                                  n_batches=n_batches, spec_file=spec_file))
    # tasks are emitted contrast-by-contrast, so the throttle drains them in
    # figure order and each panel's numbers land as early as they can rather
    # than everything finishing at once at the end

    per_task_min = args.batch * (max(NSUB[CONTRASTS[k][1]] + NSUB[CONTRASTS[k][2]]
                                     for k in keys)
                                 * SEC_PER_SUBTRACK_UNCONTENDED) / 60
    print(f'{len(keys)} contrast(s), {len(tasks)} array tasks, '
          f'throttle {args.throttle}')
    print(f'  ~{core_h:.0f} core-hours if uncontended; '
          f'~{len(tasks)*per_task_min/args.throttle/60:.1f} h wall at that throttle')
    print(f'  --mem={max_mem}G --cpus-per-task={THREADS}; '
          f'worst task ~{per_task_min:.0f} min')
    print('  contrasts: ' + ', '.join(keys))
    MAX_ARRAY = 1001   # scontrol show config | grep MaxArraySize
    if len(tasks) > MAX_ARRAY:
        need = -(-args.n_draws * 2 * len(keys) // MAX_ARRAY)
        print(f'\nERROR: {len(tasks)} tasks exceeds MaxArraySize={MAX_ARRAY}. '
              f'Raise --batch to at least {need} (currently {args.batch}).')
        sys.exit(1)
    if args.dry_run:
        return

    os.makedirs(OUT_ROOT, exist_ok=True)
    manifest = os.path.join(OUT_ROOT, '_manifest.json')
    with open(manifest, 'w') as f:
        json.dump({'tasks': tasks}, f, indent=2)

    env = (f'export PYTHONPATH={FIG}:{CODE}:$PYTHONPATH PANEL_VARIANT=36state '
           f'OMP_NUM_THREADS={THREADS} OPENBLAS_NUM_THREADS={THREADS} '
           f'MKL_NUM_THREADS={THREADS}; ')
    merge_lines = ' '.join(
        f'{PY} -u {FIG}/merge_contrast_cache.py --tier {t} '
        f'--cond-a {a} --cond-b {b} &&' for t, a, b in merges)
    merge_id = args.merge_job or sbatch(
        ['--job-name', 'tlsd_refit_merges', '--mem',
                       f'{max_mem}G', '--cpus-per-task', '2', '--partition',
                       '20', '--account', 'weissman', '--time', '6:00:00',
                       '--output', f'{OUT_ROOT}/merges_%j.log',
                       '--wrap', env + merge_lines + ' echo ALL_MERGES_DONE'])

    array_id = sbatch([
        '--job-name', 'tlsd_refit_all', '--mem', f'{max_mem}G',
        '--cpus-per-task', str(THREADS), '--partition', args.partition,
        '--account', 'weissman', '--time', '8:00:00',
        '--array', f'0-{len(tasks)-1}%{args.throttle}',
        '--dependency', f'afterok:{merge_id}',
        '--output', f'{OUT_ROOT}/refit_%A_%a.log',
        '--wrap', env + f'{PY} -u {FIG}/run_panel_refit_batch.py '
                        f'--manifest {manifest} --task-index $SLURM_ARRAY_TASK_ID'])

    agg_ids = []
    for key in keys:
        # The aggregator reloads the MERGED cache to recompute the observed
        # statistic, so it needs memory on the same scale as a refit task, not
        # a flat 12G. Two Fig S7D aggregates were OOM-killed at 12G: pooled L
        # is 8-11 GiB there against 3.9 GiB for Fig S7A. Peak is about
        # 1.5 x L (the array plus one group's slice copy) plus the marginals.
        _, _ca, _cb, _ = CONTRASTS[key]
        _l = (NSUB[_ca] + NSUB[_cb]) * 150 * 36 * 4 / 2 ** 30
        _mem = max(12, int(-(-(_l * 1.5 + 0.7) * 1.4 // 4) * 4))
        agg_ids.append(sbatch([
            '--job-name', f'tlsd_agg_{key}', '--mem', f'{_mem}G',
            '--cpus-per-task', '1', '--partition', args.partition,
            '--account', 'weissman', '--time', '2:00:00',
            '--dependency', f'afterany:{array_id}',
            '--output', f'{OUT_ROOT}/{key}/aggregate_%j.log',
            '--wrap', env + f'{PY} -u {FIG}/aggregate_panel_refit.py --key {key}']))
    print(f'\n  merges={merge_id}  array={array_id} ({len(tasks)} tasks)')
    print(f'  aggregates: {" ".join(agg_ids)}')


if __name__ == '__main__':
    main()
