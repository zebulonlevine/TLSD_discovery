#!/usr/bin/env python
"""
run_panel_refit_batch.py
==========================

One SLURM array task's worth of VBEM-REFIT resampling draws for a paper panel
contrast -- permutation (for p-values) or bootstrap (for CIs), selected by the
spec's `scheme`.

Why this exists rather than occupancy_compare.refit_permutation_null_batch:
that path's population statistic is `post_mean.sum(axis=1)`, the Dirichlet
posterior mean, which is DISPLACEMENT-weighted. Every panel in this figure set
reports the PER-FOV MEAN of sub-trajectory posterior marginals, weighting each
field of view equally, because the field of view is the unit of analysis. A
null built on the first estimator does not test the second. This worker calls
panel_core.occupancy_from_draw, so the resampled statistic is computed by the
same code as the observed one and cannot drift from it.

Draw k depends only on (seed, k) -- never on which task computed it or in what
order -- so batching across array tasks changes wall-clock only, never the
statistics.

Outputs one npz per batch: w1 per draw, and `curves` (pmf_b - pmf_a per draw)
from which every per-bin interval and p-value is an integral.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
from occupancy_compare import _w1_logD, _resolve_group

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots


def resolve_draw(scheme, fovs_a, fovs_b, seed, n_draws, k):
    """The FOV draw for index k. Keyed only by (seed, k) via SeedSequence."""
    rng = np.random.default_rng(np.random.SeedSequence(seed).spawn(n_draws)[k])
    a, b = np.asarray(fovs_a), np.asarray(fovs_b)
    if scheme == 'permutation':
        pooled = np.concatenate([a, b])
        perm = rng.permutation(len(pooled))
        return pooled[perm[:len(a)]], pooled[perm[len(a):]]
    if scheme == 'bootstrap':
        return (a[rng.integers(0, len(a), len(a))],
                b[rng.integers(0, len(b), len(b))])
    raise ValueError(f'unknown scheme {scheme!r}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--spec-file')
    ap.add_argument('--batch-index', type=int)
    ap.add_argument('--n-batches', type=int)
    # Unified mode: one array job covers every contrast x scheme x batch, so a
    # single %throttle bounds TOTAL concurrency. Separate array jobs cannot
    # share a throttle, and letting them all run saturated memory bandwidth --
    # measured 4.3x slower per draw at ~12 tasks/node than alone on a node.
    ap.add_argument('--manifest')
    ap.add_argument('--task-index', type=int)
    args = ap.parse_args()

    if args.manifest is not None:
        with open(args.manifest) as f:
            entry = json.load(f)['tasks'][args.task_index]
        args.spec_file = entry['spec_file']
        args.batch_index = entry['batch_index']
        args.n_batches = entry['n_batches']
    if args.spec_file is None:
        ap.error('pass --spec-file/--batch-index/--n-batches, '
                 'or --manifest/--task-index')

    with open(args.spec_file) as f:
        spec = json.load(f)
    scheme, n_draws = spec['scheme'], spec['n_draws']
    # The deposited specs give the tier relative to the SPT data root and the
    # output directory relative to the output root.
    if not os.path.isabs(spec['tier']):
        spec['tier'] = tlsd_paths.spt(spec['tier'])
    if not os.path.isabs(spec['out_dir']):
        spec['out_dir'] = tlsd_paths.outp(spec['out_dir'])

    batch_size = -(-n_draws // args.n_batches)
    start = args.batch_index * batch_size
    end = min(start + batch_size, n_draws)
    if start >= n_draws:
        print(f'batch {args.batch_index}: nothing to do')
        return
    k_indices = list(range(start, end))

    out_dir = os.path.join(spec['out_dir'], f'_{scheme}_batches')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f'batch_{args.batch_index:04d}.npz')
    # Resume: a completed batch is a pure function of (seed, k), so an existing
    # npz is exactly what this task would recompute. Makes re-throttling the
    # array free -- the concurrency knob matters a lot here (memory bandwidth,
    # not CPU) and may need a second attempt to get right.
    if os.path.exists(out_path):
        print(f'batch {args.batch_index} already complete: {out_path}')
        return

    ctx = pc.load(spec['tier'], spec['cache_rel'])
    fovs_a = _resolve_group(ctx['cache'], {'source_member': spec['cond_a']})
    fovs_b = _resolve_group(ctx['cache'], {'source_member': spec['cond_b']})
    lg = np.log10(ctx['diff_coefs'])
    print(f"[{spec['label']}] {scheme} batch {args.batch_index}/"
          f"{args.n_batches}: k in [{start},{end}) of {n_draws}; "
          f"n_fov {len(fovs_a)} vs {len(fovs_b)}", flush=True)

    import time
    # Checkpoint after every draw. A draw here costs 120 s alone on a node but
    # 500-1400 s under a saturated array, so an 8 h walltime is not a safe
    # margin for a 20-draw batch -- and without a checkpoint a task killed at
    # the limit writes NOTHING and its whole batch has to be recomputed.
    # Checkpointing bounds the loss to the draw in progress.
    ckpt = out_path + '.partial.npz'
    done_k, done_w1, done_curves = [], [], []
    if os.path.exists(ckpt):
        try:
            z = np.load(ckpt)
            done_k = z['k_done'].tolist()
            done_w1, done_curves = list(z['w1']), list(z['curves'])
            print(f'  resuming from checkpoint: {len(done_k)} draw(s) already '
                  f'done ({min(done_k)}..{max(done_k)})', flush=True)
        except Exception as e:
            print(f'  checkpoint unreadable ({e}); starting the batch over')
            done_k, done_w1, done_curves = [], [], []

    t0 = time.time()
    todo = [k for k in k_indices if k not in set(done_k)]
    for i, k in enumerate(todo):
        da, db = resolve_draw(scheme, fovs_a, fovs_b, spec['seed'], n_draws, k)
        pa = np.nanmean(pc.occupancy_from_draw(ctx, da), axis=0)
        pb = np.nanmean(pc.occupancy_from_draw(ctx, db), axis=0)
        done_k.append(k)
        done_w1.append(_w1_logD(lg, pa, pb))
        done_curves.append(pb - pa)
        np.savez_compressed(ckpt, k_done=np.array(done_k),
                            w1=np.array(done_w1),
                            curves=np.array(done_curves))
        if i == 0:
            print(f'  first draw took {time.time()-t0:.1f}s', flush=True)
    elapsed = time.time() - t0
    order = np.argsort(done_k)
    k_indices = list(np.array(done_k)[order])
    w1 = np.array(done_w1)[order]
    curves = np.array(done_curves)[order]

    os.makedirs(out_dir, exist_ok=True)
    np.savez_compressed(out_path, k_indices=np.array(k_indices), w1=w1,
                        curves=curves, elapsed_sec=elapsed)
    if os.path.exists(ckpt):
        os.remove(ckpt)      # batch is complete; the checkpoint is redundant
    print(f'{len(k_indices)} draws in {elapsed:.1f}s '
          f'({elapsed/len(k_indices):.1f}s/draw) -> {out_path}')


if __name__ == '__main__':
    main()
