#!/usr/bin/env python
"""
aggregate_zone_refit.py
=========================

Fig 5F from the zone resampling draws.

Three quantities, each with the null that actually matches it:

  zone effect, per treatment   (on ER - off ER), averaged over fields of view.
        p from the 2000-draw trajectory-level zone permutation. That null
        needs no refit -- per_fov_occ_by_zone fits one VBEM over the pooled
        sub-trajectories of both zones and then partitions, so relabelling
        zones leaves the fit bit-identical. Being refit-free is why it runs at
        2000 draws and is NOT capped at the 1/501 floor the FOV-shuffle
        contrasts hit.
        CI from the 500-draw refit bootstrap, which does refit, because
        resampling fields of view changes the pooled set.

  zone x treatment interaction   the difference of the two zone effects.
        This does NOT have the 2x2 Freedman-Lane problem: the zone contrast is
        taken WITHIN each field of view, so the zone main effect is differenced
        away before any resampling and the response is one number per FOV.
        Permuting treatment across FOVs is then an exact two-sample test.
        CI by all-pairs from the two treatments' independent bootstraps; the
        pairs are not independent, so the CI's resolution is the per-arm
        resolution (recorded as p_floor). The permutation p uses the
        treatment-permutation count (add-one, floor 1/(n_perm+1)).
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
import paperfig_style as ps
from multiplicity import bh, perm_p_addone
from region_schemes import SCHEME_60S

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
OUT = tlsd_paths.outp('PaperFigures_updated')
D = os.path.join(OUT, 'refit_resampling', 'fig5F')
TREAT = {'noDegron': 'DMSO', 'Degron': 'PT-179 (24h)'}


def load(cond, scheme):
    fs = sorted(glob.glob(os.path.join(D, f'zone_{cond}_{scheme}_*.npz')))
    if not fs:
        raise SystemExit(f'FAIL: no {cond} {scheme} draws in {D}')
    parts, obs, dc = [], None, None
    for f in fs:
        z = np.load(f, allow_pickle=True)
        parts.append(z['draws']); obs, dc = z['obs'], z['diff_coefs']
    return obs, np.concatenate(parts, axis=0), dc


def per_bin(arr, dc):
    """(..., n_fov, n_D) -> (..., n_bins), integrating each bin then averaging
    over fields of view."""
    out = []
    for lo, hi in SCHEME_60S.values():
        m = (dc >= lo) & (dc < hi)
        out.append(np.nanmean(np.nansum(arr[..., m], axis=-1), axis=-1))
    return np.stack(out, axis=-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-perm-treat', type=int, default=10000)
    a = ap.parse_args()

    obs_b, boot_b, dc = {}, {}, None
    rows = []
    per_fov_obs = {}
    for cond in ('noDegron', 'Degron'):
        o, perm, dc = load(cond, 'permutation')
        _, boot, _ = load(cond, 'bootstrap')
        per_fov_obs[cond] = o
        obs = per_bin(o, dc)
        null = per_bin(perm, dc)
        bt = per_bin(boot, dc)
        obs_b[cond], boot_b[cond] = obs, bt
        p = np.array([perm_p_addone(
            float((np.abs(null[:, j]) >= abs(obs[j])).mean()), null.shape[0])
            for j in range(len(obs))])
        pb = bh(p)
        for j, name in enumerate(SCHEME_60S):
            rows.append(dict(
                panel='5F_zone_effect', treatment=TREAT[cond], bin=name,
                observed=obs[j], ci_lo=float(np.percentile(bt[:, j], 2.5)),
                ci_hi=float(np.percentile(bt[:, j], 97.5)),
                p_perm=p[j], p_bh=pb[j], stars=ps.stars(pb[j]),
                n_perm=null.shape[0], n_boot=bt.shape[0]))
        print(f'[{TREAT[cond]}] {null.shape[0]} zone permutations, '
              f'{bt.shape[0]} refit bootstrap draws')

    # interaction: exact treatment permutation on per-FOV zone differences
    A = per_bin(per_fov_obs['noDegron'][None], dc)[0]
    B = per_bin(per_fov_obs['Degron'][None], dc)[0]
    a_fov = np.stack([np.nansum(per_fov_obs['noDegron'][:, (dc >= lo) & (dc < hi)], axis=1)
                      for lo, hi in SCHEME_60S.values()], axis=1)
    b_fov = np.stack([np.nansum(per_fov_obs['Degron'][:, (dc >= lo) & (dc < hi)], axis=1)
                      for lo, hi in SCHEME_60S.values()], axis=1)
    rng = np.random.default_rng(0)
    pooled = np.concatenate([a_fov, b_fov], axis=0)
    na = a_fov.shape[0]
    did_obs = np.nanmean(b_fov, axis=0) - np.nanmean(a_fov, axis=0)
    null_t = np.empty((a.n_perm_treat, len(did_obs)))
    for i in range(a.n_perm_treat):
        pm = rng.permutation(pooled.shape[0])
        null_t[i] = (np.nanmean(pooled[pm[na:]], axis=0)
                     - np.nanmean(pooled[pm[:na]], axis=0))
    # all-pairs CI from the two refit bootstraps
    pairs = (boot_b['Degron'][None, :, :] - boot_b['noDegron'][:, None, :]
             ).reshape(-1, len(did_obs))
    floor = 1.0 / (boot_b['Degron'].shape[0] + 1)
    p_i = np.array([perm_p_addone(
        float((np.abs(null_t[:, j]) >= abs(did_obs[j])).mean()), a.n_perm_treat)
        for j in range(len(did_obs))])
    p_ibh = bh(p_i)
    for j, name in enumerate(SCHEME_60S):
        rows.append(dict(
            panel='5F_interaction', treatment='PT-179 - DMSO', bin=name,
            observed=did_obs[j],
            ci_lo=float(np.percentile(pairs[:, j], 2.5)),
            ci_hi=float(np.percentile(pairs[:, j], 97.5)),
            p_perm=p_i[j], p_bh=p_ibh[j], stars=ps.stars(p_ibh[j]),
            n_perm=a.n_perm_treat, n_boot=pairs.shape[0], p_floor=floor))

    df = pd.DataFrame(rows)
    out = os.path.join(OUT, 'figure5/F/fig5F_refit_stats.csv')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(pc.tag(out), index=False)
    for panel in ('5F_zone_effect', '5F_interaction'):
        print(f'\n=== {panel} ===')
        print(df[df.panel == panel][['treatment', 'bin', 'observed', 'ci_lo',
                                     'ci_hi', 'p_perm', 'p_bh', 'stars']]
              .to_string(index=False))
    print(f'\nwrote {pc.tag(out)}')


if __name__ == '__main__':
    main()
