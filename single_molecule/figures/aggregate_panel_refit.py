#!/usr/bin/env python
"""
aggregate_panel_refit.py
==========================

Combine one contrast's permutation and bootstrap batches into the reported
statistics: W1 with a bootstrap CI and a permutation p-value, and the same
per diffusion bin.

Per-bin quantities are integrals of the stored per-draw pmf-difference curves,
so the bins can be changed without refitting anything.

Checks completeness explicitly rather than relying on `afterok`: every draw
index in range(n_draws) must be present exactly once, for both schemes, or
this fails loudly. A silently short null would just quietly shrink p.
"""
import argparse, glob, json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
from occupancy_compare import _w1_logD
from region_schemes import SCHEME_60S, SCHEME_40S
from multiplicity import bh, holm, perm_p_addone
from submit_panel_refit import CONTRASTS, TIERS, OUT_ROOT

SCHEME_FOR = {'claim1_40S': SCHEME_40S}


def load_scheme_batches(out_dir, scheme, n_draws):
    files = sorted(glob.glob(os.path.join(out_dir, f'_{scheme}_batches', '*.npz')))
    if not files:
        raise SystemExit(f'FAIL: no {scheme} batches in {out_dir}')
    ks, w1s, curves = [], [], []
    for f in files:
        d = np.load(f)
        ks.append(d['k_indices']); w1s.append(d['w1']); curves.append(d['curves'])
    k = np.concatenate(ks)
    order = np.argsort(k)
    k = k[order]
    missing = sorted(set(range(n_draws)) - set(k.tolist()))
    dupes = sorted({int(x) for x in k[:-1][k[1:] == k[:-1]]})
    if missing or dupes:
        raise SystemExit(
            f'FAIL: {scheme} incomplete -- {len(missing)} missing '
            f'{missing[:8]}{"..." if len(missing) > 8 else ""}, '
            f'{len(dupes)} duplicated {dupes[:8]}')
    return (np.concatenate(w1s)[order], np.concatenate(curves)[order])


def bin_masses(pmf, diff_coefs, scheme):
    return np.array([pmf[(diff_coefs >= lo) & (diff_coefs < hi)].sum()
                     for lo, hi in scheme.values()])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--key', required=True, choices=sorted(CONTRASTS))
    args = ap.parse_args()
    key = args.key
    tier_key, cond_a, cond_b, label = CONTRASTS[key]
    tier, out_dir = TIERS[tier_key], os.path.join(OUT_ROOT, key)
    spec = json.load(open(os.path.join(out_dir, '_spec_permutation.json')))
    n_draws = spec['n_draws']
    scheme = SCHEME_FOR.get(key, SCHEME_60S)

    null_w1, null_curves = load_scheme_batches(out_dir, 'permutation', n_draws)
    boot_w1, boot_curves = load_scheme_batches(out_dir, 'bootstrap', n_draws)
    print(f'[{key}] {n_draws} permutation + {n_draws} bootstrap draws, complete')

    # observed: the same code path the panels use
    ctx = pc.load(tier, spec['cache_rel'])
    dc = ctx['diff_coefs']
    occ_a, _, _ = pc.per_fov_occupancy(ctx, {'source_member': cond_a})
    occ_b, _, _ = pc.per_fov_occupancy(ctx, {'source_member': cond_b})
    pmf_a, pmf_b = np.nanmean(occ_a, axis=0), np.nanmean(occ_b, axis=0)
    lg = np.log10(dc)
    w1_obs = _w1_logD(lg, pmf_a, pmf_b)

    p_w1 = perm_p_addone(float((null_w1 >= w1_obs).mean()), n_draws)
    rows = [dict(quantity='W1', observed=w1_obs,
                 ci_lo=float(np.percentile(boot_w1, 2.5)),
                 ci_hi=float(np.percentile(boot_w1, 97.5)),
                 null_p95=float(np.percentile(null_w1, 95)),
                 p_perm=p_w1, p_bh=np.nan, stars='')]

    edges = list(scheme.values())
    idx = [np.flatnonzero((dc >= lo) & (dc < hi)) for lo, hi in edges]
    obs_delta = bin_masses(pmf_b, dc, scheme) - bin_masses(pmf_a, dc, scheme)
    null_delta = np.stack([null_curves[:, i].sum(axis=1) for i in idx], axis=1)
    boot_delta = np.stack([boot_curves[:, i].sum(axis=1) for i in idx], axis=1)
    p_bins = np.array([perm_p_addone(
        float((np.abs(null_delta[:, j]) >= abs(obs_delta[j])).mean()), n_draws)
        for j in range(len(edges))])
    p_bins_bh = bh(p_bins)
    import paperfig_style as ps
    for j, name in enumerate(scheme):
        rows.append(dict(
            quantity=f'bin_{name}', observed=obs_delta[j],
            ci_lo=float(np.percentile(boot_delta[:, j], 2.5)),
            ci_hi=float(np.percentile(boot_delta[:, j], 97.5)),
            null_p95=float(np.percentile(np.abs(null_delta[:, j]), 95)),
            p_perm=p_bins[j], p_bh=p_bins_bh[j], stars=ps.stars(p_bins_bh[j])))

    df = pd.DataFrame(rows)
    df['contrast'], df['label'] = key, label
    df['cond_a'], df['cond_b'] = cond_a, cond_b
    df['n_draws'], df['refit'] = n_draws, True
    csv = os.path.join(out_dir, f'{key}_refit_stats.csv')
    df.to_csv(csv, index=False)
    np.savez_compressed(os.path.join(out_dir, f'{key}_refit_draws.npz'),
                        null_w1=null_w1, boot_w1=boot_w1,
                        null_delta=null_delta, boot_delta=boot_delta,
                        obs_delta=obs_delta, pmf_a=pmf_a, pmf_b=pmf_b)
    print(df[['quantity', 'observed', 'ci_lo', 'ci_hi', 'null_p95',
              'p_perm', 'p_bh', 'stars']].to_string(index=False))
    print(f'wrote {csv}')


if __name__ == '__main__':
    main()
