#!/usr/bin/env python
"""Merge two per-condition 36-state likelihood caches into one pooled cache.

The variant cache layout stores one cache per condition, but a permutation
null has to relabel fields of view ACROSS the two conditions, which needs both
in a single L array. Merging is exact rather than an approximation: the naive
per-sub-trajectory likelihood depends only on that sub-trajectory's own jumps
and the shared grid, so concatenation gives byte-identical L to one joint
saSPT call (saspt_likelihood_cache.merge_likelihood_caches).

Idempotent -- skips if the merged cache already exists.
"""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saspt_likelihood_cache import merge_likelihood_caches
import panel_core as pc

ap = argparse.ArgumentParser()
ap.add_argument('--tier', required=True)
ap.add_argument('--cond-a', required=True)
ap.add_argument('--cond-b', required=True)
a = ap.parse_args()

cdir = os.path.join(a.tier, pc.CACHE_SUBDIR)
name = f'_joint_{a.cond_a}_vs_{a.cond_b}'
out = os.path.join(cdir, f'{name}_likelihood_cache.npz')
if os.path.exists(out):
    print(f'already present, skipping: {out}')
    sys.exit(0)
p = merge_likelihood_caches(
    {a.cond_a: os.path.join(cdir, f'{a.cond_a}_likelihood_cache.npz'),
     a.cond_b: os.path.join(cdir, f'{a.cond_b}_likelihood_cache.npz')},
    cdir, name, force=False, verbose=True)
print(f'wrote {p}  ({os.path.getsize(p)/2**30:.2f} GiB)')
