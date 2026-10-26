#!/usr/bin/env python
"""
write_refit_sidecars.py
=========================

Consolidate the per-contrast refit statistics into one sidecar CSV per panel
directory, so each panel ships with the numbers that are actually plotted on
it. Without this the refit stats live only under refit_resampling/ keyed by
contrast, while the panel directory still carries the superseded no-refit CSV.

Adds two columns the raw per-contrast files do not carry:

  p_at_floor   True when p equals the smallest value the draw count can
               express, 1/(B+1). At B=500 that is 0.002, so *** (p<0.001) is
               UNREACHABLE for these contrasts. The superseded no-refit files
               show *** off 2000 draws, so a reader comparing the two would
               otherwise read the refit as weaker when it is only
               coarser-grained.
  supersedes   the no-refit file this replaces, so the provenance is explicit
               rather than implied by a filename suffix. Empty for S7A and S7C:
               their no-refit builds are not part of this repository.
"""
import glob, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
from submit_panel_refit import CONTRASTS, OUT_ROOT

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
OUT = tlsd_paths.outp('PaperFigures_updated')

# panel dir -> (sidecar basename, contrasts, the no-refit file it supersedes)
SIDECARS = {
    'figureS7/A': ('figS7A_W1_refit_stats', ['figS7AC_60S', 'figS7A_40S'], ''),
    'figureS7/C': ('figS7C_60S_bin_difference_refit_stats', ['figS7AC_60S'], ''),
    'figureS7/D': ('figS7D_refit_stats',
                   ['figS7D_NT_HHT', 'figS7D_NT_eme', 'figS7D_UBA5_HHT',
                    'figS7D_UBA5_eme', 'figS7D_geno_none', 'figS7D_geno_HHT',
                    'figS7D_geno_eme'],
                   'figS7D_stats_36state.csv'),
}


def main():
    for panel_dir, (base, keys, supersedes) in SIDECARS.items():
        # Partial rather than all-or-nothing: a multi-contrast panel should
        # still ship a sidecar for the contrasts that are done, with the
        # missing ones named in the file, instead of having none at all.
        frames, absent = [], []
        for k in keys:
            p = os.path.join(OUT_ROOT, k, f'{k}_refit_stats.csv')
            if not os.path.exists(p):
                absent.append(k)
                continue
            d = pd.read_csv(p)
            _, ca, cb, label = CONTRASTS[k]
            # aggregate_panel_refit already writes `contrast`/`label`; reorder
            # rather than insert so this is idempotent
            d['contrast'], d['contrast_label'] = k, label
            front = ['contrast', 'contrast_label', 'quantity']
            d = d[front + [c for c in d.columns if c not in front]]
            n = int(d['n_draws'].iloc[0]) if 'n_draws' in d else 500
            floor = 1.0 / (n + 1)
            d['p_at_floor'] = np.isclose(d['p_perm'], floor, rtol=1e-6)
            frames.append(d)
        if not frames:
            print(f'  {panel_dir}: nothing aggregated yet -- skipped')
            continue
        df = pd.concat(frames, ignore_index=True)
        df['supersedes'] = supersedes
        df['contrasts_missing'] = ';'.join(absent) if absent else ''
        df['refit'] = True
        out = pc.tag(os.path.join(OUT, panel_dir, f'{base}.csv'))
        os.makedirs(os.path.dirname(out), exist_ok=True)
        df.to_csv(out, index=False)
        nf = int(df['p_at_floor'].sum())
        msg = (f'  wrote {os.path.relpath(out, OUT)}  '
               f'({len(keys) - len(absent)}/{len(keys)} contrast(s), '
               f'{len(df)} rows, {nf} p-value(s) at the '
               f'{1/(int(df["n_draws"].iloc[0])+1):.4f} floor)')
        if absent:
            msg += f'  MISSING: {",".join(absent)}'
        print(msg)


if __name__ == '__main__':
    main()
