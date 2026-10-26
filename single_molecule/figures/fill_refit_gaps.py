#!/usr/bin/env python
"""
fill_refit_gaps.py
====================

Find and resubmit refit batches that never landed -- tasks killed by the
walltime, which write nothing because the npz is only saved after all of a
batch's draws finish.

IMPORTANT: resubmits at the SAME --batch size as the original run. Batch index
k covers draws [k*batch, (k+1)*batch), so a different batch size remaps every
index onto a different draw range while the filenames stay batch_NNNN.npz. The
resume check tests filename existence, so a changed batch size would silently
skip the wrong ranges. (The aggregator's completeness check over draw indices
would catch the resulting gap, but only after the fact.)

Contention is what caused the timeouts -- 121 s/draw alone on a node against
529-1384 s under a saturated array -- so gap-filling once the array has
drained runs at close to full speed.
"""
import argparse, glob, json, os, subprocess, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from submit_panel_refit import CONTRASTS, TIERS, NSUB, OUT_ROOT, CODE, FIG, sbatch

PY = sys.executable
SEC_PER_SUBTRACK = 121.1 / 192701


def infer_batch(key, scheme, default):
    """Read the batch size back off an existing npz.

    Contrasts were not all run at the same batch size -- figS7AC_60S was the
    pilot at 10, the rest at 20 -- and assuming one global value flags every
    pilot batch as inconsistent. Batch index k covers [k*batch, (k+1)*batch),
    so `batch` is recoverable as the span of any stored k_indices.
    """
    for p in sorted(glob.glob(os.path.join(
            OUT_ROOT, key, f'_{scheme}_batches', 'batch_*.npz'))):
        try:
            k = np.load(p)['k_indices']
            if len(k) > 1:
                return int(k.max() - k.min() + 1)
        except Exception:
            continue
    return default


def missing(key, scheme, n_draws, batch):
    """Batch indices with no npz, plus any whose stored k_indices disagree."""
    d = os.path.join(OUT_ROOT, key, f'_{scheme}_batches')
    n_batches = -(-n_draws // batch)
    gaps = []
    for b in range(n_batches):
        p = os.path.join(d, f'batch_{b:04d}.npz')
        if not os.path.exists(p):
            # a .partial.npz means a task was killed mid-batch; the refill
            # resumes from it rather than recomputing the draws it did finish
            part = p + '.partial.npz'
            if os.path.exists(part):
                try:
                    n = len(np.load(part)['k_done'])
                    print(f'  ~  {key}/{scheme} batch {b}: partial, '
                          f'{n} draw(s) salvageable')
                except Exception:
                    pass
            gaps.append(b); continue
        try:
            k = np.load(p)['k_indices']
            if int(k.min()) != b * batch:
                print(f'  !! {key}/{scheme} batch {b}: k starts at {k.min()}, '
                      f'expected {b*batch} -- batch size changed between runs')
                gaps.append(b)
        except Exception as e:
            print(f'  !! {key}/{scheme} batch {b} unreadable ({e}); refilling')
            gaps.append(b)
    return gaps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-draws', type=int, default=500)
    ap.add_argument('--batch', type=int, default=20,
                    help='MUST match the original run')
    ap.add_argument('--skip', action='append', default=[],
                    help="contrast keys to ignore (e.g. dropped contrasts)")
    ap.add_argument('--throttle', type=int, default=24,
                    help='low, since timeouts were caused by contention')
    ap.add_argument('--time', default='8:00:00')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--force', action='store_true',
                    help='run even with refit jobs still in flight')
    a = ap.parse_args()

    inflight = subprocess.run(
        ['squeue', '-u', os.environ.get('USER', ''), '-h', '-t',
         'RUNNING,PENDING', '-o', '%j'], capture_output=True, text=True).stdout
    n_live = sum(1 for ln in inflight.split('\n') if 'refit' in ln)
    if n_live and not a.force:
        print(f'REFUSING: {n_live} refit job(s) still RUNNING/PENDING. Batches '
              f'in flight have written nothing yet and would be counted as '
              f'gaps, duplicating work. Wait for the array to drain, or pass '
              f'--force if you know better.')
        sys.exit(1)

    tasks = []
    for key in sorted(CONTRASTS):
        if key in a.skip:
            continue
        for scheme in ('permutation', 'bootstrap'):
            spec = os.path.join(OUT_ROOT, key, f'_spec_{scheme}.json')
            if not os.path.exists(spec):
                continue
            batch = infer_batch(key, scheme, a.batch)
            if batch != a.batch:
                print(f'  {key}/{scheme}: batch size {batch} (not {a.batch}) '
                      f'-- inferred from stored k_indices')
            g = missing(key, scheme, a.n_draws, batch)
            if g:
                print(f'  {key:<18} {scheme:<12} missing {len(g)}: {g}')
                tasks += [dict(spec_file=spec, batch_index=b,
                               n_batches=-(-a.n_draws // batch)) for b in g]
    if not tasks:
        print('no gaps -- every batch present and consistent')
        return
    keys_used = {json.load(open(t['spec_file']))['key'] for t in tasks}
    l_max = max((NSUB[CONTRASTS[k][1]] + NSUB[CONTRASTS[k][2]])
                * 150 * 36 * 4 / 2 ** 30 for k in keys_used)
    mem = max(12, int(-(-(l_max * 1.7 + 3) // 4) * 4))
    print(f'\n{len(tasks)} gap task(s), --mem={mem}G, throttle {a.throttle}')
    if a.dry_run:
        return
    man = os.path.join(OUT_ROOT, '_manifest_gapfill.json')
    with open(man, 'w') as f:
        json.dump({'tasks': tasks}, f, indent=2)
    env = (f'export PYTHONPATH={FIG}:{CODE}:$PYTHONPATH PANEL_VARIANT=36state '
           f'OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2; ')
    jid = sbatch(['--job-name', 'tlsd_refit_gapfill', '--mem', f'{mem}G',
                  '--cpus-per-task', '2', '--partition', '24,20,weissman',
                  '--account', 'weissman', '--time', a.time,
                  '--array', f'0-{len(tasks)-1}%{a.throttle}',
                  '--output', f'{OUT_ROOT}/gapfill_%A_%a.log',
                  '--wrap', env + f'{PY} -u {FIG}/run_panel_refit_batch.py '
                                  f'--manifest {man} '
                                  f'--task-index $SLURM_ARRAY_TASK_ID'])
    print(f'  gapfill array = {jid}')


if __name__ == '__main__':
    main()
