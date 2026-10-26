#!/usr/bin/env python
"""
submit_panel_refit.py
=======================

Submit the VBEM-refit resampling for one or more paper-panel contrasts:
a merge job, then a permutation array (p-values) and a bootstrap array (CIs),
then an aggregation job.

    python submit_panel_refit.py --contrast claim1_60S
    python submit_panel_refit.py --all --dry-run

Sizing (measured, diagnostics/time_vbem_refit.py): 0.93 ms per sub-trajectory
per VBEM fit at 2 threads. Two threads beats eight on CORE-seconds (162 vs
247 for claim 1), so under the per-user CPU cap narrow tasks give ~50% more
throughput; memory, not CPU, is what actually binds concurrency here.

Dependencies use `afterany` plus an explicit completeness check in the
aggregator, not `afterok`: a single failed array task under `afterok` leaves
the aggregate permanently blocked as DependencyNeverSatisfied, which is both
silent and annoying to unpick. `afterany` + a check fails loudly instead.
"""

import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
A = tlsd_paths.spt('all_data')
CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pipeline'))
FIG = os.path.dirname(os.path.abspath(__file__))           # figure code
OUT_ROOT = tlsd_paths.outp('PaperFigures_updated', 'refit_resampling')
PY = sys.executable

TIERS = {
    'claim1': f'{A}/2025_06_06_zeb/tracking_output_q=15p0_splitsize10',
    'claim2': f'{A}/2026_07_08/tracking_output_q=15p0',
    'claim3': f'{A}/2025_09_15_zeb_ER/tracking_output_q=15p0',
}

# n_sub per condition, read from the cache headers -- used only for sizing.
NSUB = {
    'sgNT_60S': 83158, 'sgUBA5_60S': 109543,
    'sgNT_40S': 166062, 'sgUBA5_40S': 91480,
    'S3_NT_none': 204998, 'S3_NT_HT10': 168080, 'S3_NT_eme360': 134970,
    'S3_UBA5_none': 271998, 'S3_UBA5_HT10': 132139, 'S3_UBA5_eme360': 264176,
    'noDegron': 115656, 'Degron': 162860,
}

CONTRASTS = {
    'claim1_60S':      ('claim1', 'sgNT_60S', 'sgUBA5_60S', 'Fig S7A/S7C 60S'),
    'claim1_40S':      ('claim1', 'sgNT_40S', 'sgUBA5_40S', 'Fig S7A 40S'),
    'claim2_NT_HHT':   ('claim2', 'S3_NT_none', 'S3_NT_HT10', 'Fig S7D sgNT HHT'),
    'claim2_NT_eme':   ('claim2', 'S3_NT_none', 'S3_NT_eme360', 'Fig S7D sgNT emetine'),
    'claim2_UBA5_HHT': ('claim2', 'S3_UBA5_none', 'S3_UBA5_HT10', 'Fig S7D sgUBA5 HHT'),
    'claim2_UBA5_eme': ('claim2', 'S3_UBA5_none', 'S3_UBA5_eme360', 'Fig S7D sgUBA5 emetine'),
    'claim2_geno_none':('claim2', 'S3_NT_none', 'S3_UBA5_none', 'Fig S7D genotype, untreated'),
    'claim2_geno_HHT': ('claim2', 'S3_NT_HT10', 'S3_UBA5_HT10', 'Fig S7D genotype, HHT'),
    'claim2_geno_eme': ('claim2', 'S3_NT_eme360', 'S3_UBA5_eme360', 'Fig S7D genotype, emetine'),
}

SEC_PER_SUBTRACK = 0.93e-3
THREADS = 2


def sizing(cond_a, cond_b, n_draws, batch):
    """(mem_gb, minutes_per_task, n_batches) for one scheme's array."""
    pooled = NSUB[cond_a] + NSUB[cond_b]
    # peak RSS ~ pooled L held in memory + one group's slice copied by _fit_idx
    l_gib = pooled * 150 * 36 * 4 / 2 ** 30
    mem = int(max(12, -(-(l_gib * 1.7 + 3) // 4) * 4))
    per_draw = pooled * SEC_PER_SUBTRACK
    n_batches = -(-n_draws // batch)
    minutes = per_draw * batch / 60
    return mem, minutes, n_batches


def sbatch(args_list):
    r = subprocess.run(['sbatch', '--parsable'] + args_list,
                       capture_output=True, text=True, cwd=CODE)
    if r.returncode != 0:
        raise RuntimeError(f'sbatch failed: {r.stderr}')
    return r.stdout.strip().split(';')[0]


def submit_contrast(key, n_draws, batch, dry_run=False):
    tier_key, cond_a, cond_b, label = CONTRASTS[key]
    tier = TIERS[tier_key]
    out_dir = os.path.join(OUT_ROOT, key)
    merged_rel = os.path.join('occupancy_compare_cache_36state',
                              f'_joint_{cond_a}_vs_{cond_b}_likelihood_cache.npz')
    mem, minutes, n_batches = sizing(cond_a, cond_b, n_draws, batch)
    walltime = max(2, int(minutes / 60 * 3) + 1)   # 3x headroom

    print(f'\n=== {key}  ({label}) ===')
    print(f'  {cond_a} vs {cond_b}   pooled n_sub='
          f'{NSUB[cond_a]+NSUB[cond_b]:,}')
    print(f'  {n_draws} draws x 2 schemes, batch={batch} -> {n_batches} tasks each')
    print(f'  --mem={mem}G --cpus-per-task={THREADS} --time={walltime}:00:00 '
          f'(~{minutes:.0f} min/task)')
    core_h = (NSUB[cond_a]+NSUB[cond_b]) * SEC_PER_SUBTRACK * n_draws * 2 * THREADS / 3600
    print(f'  ~{core_h:.0f} core-hours')
    if dry_run:
        return None

    os.makedirs(out_dir, exist_ok=True)
    merge_id = sbatch([
        '--job-name', f'tlsd_merge_{key}', '--mem', f'{mem}G',
        '--cpus-per-task', '2', '--partition', '20', '--account', 'weissman',
        '--time', '2:00:00', '--output', f'{out_dir}/merge_%j.log',
        '--wrap',
        f'export PYTHONPATH={FIG}:{CODE}:$PYTHONPATH PANEL_VARIANT=36state; '
        f'{PY} -u '
        f'{FIG}/merge_contrast_cache.py --tier {tier} --cond-a {cond_a} '
        f'--cond-b {cond_b}'])

    ids = {'merge': merge_id}
    for scheme in ('permutation', 'bootstrap'):
        spec = dict(tier=tier, cache_rel=merged_rel, cond_a=cond_a,
                    cond_b=cond_b, label=label, scheme=scheme,
                    n_draws=n_draws, seed=0, out_dir=out_dir, key=key)
        spec_file = os.path.join(out_dir, f'_spec_{scheme}.json')
        with open(spec_file, 'w') as f:
            json.dump(spec, f, indent=2)
        ids[scheme] = sbatch([
            '--job-name', f'tlsd_{scheme[:4]}_{key}', '--mem', f'{mem}G',
            '--cpus-per-task', str(THREADS), '--partition', '20',
            '--account', 'weissman', '--time', f'{walltime}:00:00',
            '--array', f'0-{n_batches-1}',
            '--dependency', f'afterok:{merge_id}',
            '--output', f'{out_dir}/{scheme}_%A_%a.log',
            '--wrap',
            f'export PYTHONPATH={FIG}:{CODE}:$PYTHONPATH PANEL_VARIANT=36state '
            f'OMP_NUM_THREADS={THREADS} OPENBLAS_NUM_THREADS={THREADS} '
            f'MKL_NUM_THREADS={THREADS}; '
            f'{PY} -u '
            f'{FIG}/run_panel_refit_batch.py --spec-file {spec_file} '
            f'--batch-index $SLURM_ARRAY_TASK_ID --n-batches {n_batches}'])

    ids['aggregate'] = sbatch([
        '--job-name', f'tlsd_agg_{key}', '--mem', '8G',
        '--cpus-per-task', '1', '--partition', '20', '--account', 'weissman',
        '--time', '1:00:00',
        '--dependency', f'afterany:{ids["permutation"]}:{ids["bootstrap"]}',
        '--output', f'{out_dir}/aggregate_%j.log',
        '--wrap',
        f'export PYTHONPATH={FIG}:{CODE}:$PYTHONPATH PANEL_VARIANT=36state; '
        f'{PY} -u '
        f'{FIG}/aggregate_panel_refit.py --key {key}'])

    print('  jobs: ' + '  '.join(f'{k}={v}' for k, v in ids.items()))
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--contrast', action='append', choices=sorted(CONTRASTS))
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--n-draws', type=int, default=500)
    ap.add_argument('--batch', type=int, default=10)
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    keys = sorted(CONTRASTS) if args.all else (args.contrast or [])
    if not keys:
        ap.error('pass --contrast KEY (repeatable) or --all')
    total = 0.0
    for k in keys:
        _, a, b, _ = CONTRASTS[k]
        total += (NSUB[a]+NSUB[b]) * SEC_PER_SUBTRACK * args.n_draws * 2 * THREADS / 3600
        submit_contrast(k, args.n_draws, args.batch, args.dry_run)
    print(f'\nTOTAL across {len(keys)} contrast(s): ~{total:.0f} core-hours')


if __name__ == '__main__':
    main()
