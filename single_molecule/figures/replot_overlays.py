#!/usr/bin/env python
"""
replot_overlays.py
==================

Redraw any posterior-density overlay from its `_curves.csv` + `_plotspec.json`
sidecars. No VBEM, no likelihood cache, no cluster job -- it reads the pmf that
was plotted and hands it back to `panel_core.plot_overlay`.

Use it for anything presentational: fills, alpha, colours, axis limits, the
y-cap. Use a builder only when a NUMBER has to change, because only the builder
refits.

    # restyle every overlay, faithfully (should be a no-op)
    python replot_overlays.py --all --check

    # what the Sep 2026 fill change would have cost with this in place
    python replot_overlays.py --all --fill all

    # one panel, fainter tint
    python replot_overlays.py figure5/C/fig5C_60S_posterior_overlay_36state \
        --fill all --fill-alpha 0.10

`--check` re-emits to a temporary directory and rasterizes both copies at
300 dpi instead of overwriting, so "is the sidecar actually sufficient to
reproduce the committed panel" is a question you can answer rather than assume.

The sidecars are written unconditionally by plot_overlay, so a panel has them
from its first build onwards. A panel built before Sep 2026 has none until its
builder is run once more; `--all` reports those rather than skipping silently.
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
import paperfig_style as ps

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
HERE = tlsd_paths.outp('PaperFigures_updated')


def out_pdf_for(spec_path, suffix):
    """The untagged PDF path plot_overlay expects, derived from where the spec
    sits rather than from anything stored in it -- so the sidecars survive the
    tree being moved or exported."""
    stem = spec_path[:-len('_plotspec.json')]
    if suffix and stem.endswith(suffix):
        stem = stem[:-len(suffix)]
    return stem + '.pdf'


def load(spec_path):
    with open(spec_path) as fh:
        spec = json.load(fh)
    spec['out_pdf'] = out_pdf_for(spec_path, spec.get('suffix', ''))
    df = pd.read_csv(os.path.join(os.path.dirname(spec_path), spec['curves_csv']))
    curves = [dict(pmf=df[f"pmf_{c['key']}"].to_numpy(), label=c['label'],
                   color=c['color'], ls=c['ls'], lw=c['lw'], fill=c['fill'],
                   **({'alpha': c['alpha']} if c.get('alpha') is not None else {}))
              for c in spec['curves']]
    return spec, df['D'].to_numpy(), curves


def emit(spec, D, curves, out_pdf, fill=None, fill_alpha=None, ymax='saved'):
    if fill == 'all':
        for c in curves:
            c['fill'] = True
    elif fill == 'none':
        for c in curves:
            c['fill'] = False
    return pc.plot_overlay(
        curves, D, {k: tuple(v) for k, v in spec['scheme'].items()}, out_pdf,
        size=spec['size'], xlabel=spec['xlabel'], ylabel=spec['ylabel'],
        xlim=tuple(spec['xlim']), legend_loc=spec['legend_loc'],
        ymax=spec['ymax'] if ymax == 'saved' else ymax,
        fill_alpha=fill_alpha if fill_alpha is not None else spec['fill_alpha'])


def raster(pdf, tmp):
    subprocess.run(['pdftoppm', '-r', '300', '-png', '-f', '1', '-l', '1',
                    pdf, os.path.join(tmp, 'p')], check=True)
    from PIL import Image
    f = [p for p in os.listdir(tmp) if p.endswith('.png')][0]
    return np.asarray(Image.open(os.path.join(tmp, f)).convert('RGB')).astype(int)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('panels', nargs='*',
                    help='panel paths relative to PaperFigures_updated, with or '
                         'without the _plotspec.json / .pdf suffix')
    ap.add_argument('--all', action='store_true', help='every panel with sidecars')
    ap.add_argument('--fill', choices=('saved', 'all', 'none'), default='saved')
    ap.add_argument('--fill-alpha', type=float)
    ap.add_argument('--ymax', type=float)
    ap.add_argument('--check', action='store_true',
                    help='re-emit to a temp dir and raster-compare instead of '
                         'overwriting')
    a = ap.parse_args()
    ps.apply()

    if a.all:
        specs = sorted(glob.glob(os.path.join(HERE, 'figure*', '**',
                                              '*_plotspec.json'), recursive=True))
    else:
        specs = []
        for p in a.panels:
            base = os.path.join(HERE, p)
            for cand in (base, base.replace('.pdf', ''),
                         base.replace('_plotspec.json', '')):
                if os.path.exists(cand + '_plotspec.json'):
                    specs.append(cand + '_plotspec.json'); break
            else:
                print(f'  no sidecar for {p} -- run its builder once to write one')
    if not specs:
        print('nothing to do'); return 1

    rc = 0
    for sp in specs:
        rel = os.path.relpath(sp, HERE)
        spec, D, curves = load(sp)
        if spec['suffix'] != pc.SUFFIX:
            print(f'  SKIP {rel}: written under suffix "{spec["suffix"]}" but '
                  f'this process has "{pc.SUFFIX}" -- set PANEL_VARIANT to match')
            rc = 1
            continue
        if not a.check:
            out = emit(spec, D, curves, spec['out_pdf'], a.fill, a.fill_alpha,
                       a.ymax if a.ymax is not None else 'saved')
            print(f'  wrote {os.path.relpath(out, HERE)}')
            continue
        committed = pc.tag(spec['out_pdf'])
        with tempfile.TemporaryDirectory() as td:
            probe = os.path.join(td, os.path.basename(spec['out_pdf']))
            emit(spec, D, curves, probe, a.fill, a.fill_alpha,
                 a.ymax if a.ymax is not None else 'saved')
            with tempfile.TemporaryDirectory() as ta, \
                    tempfile.TemporaryDirectory() as tb:
                x, y = raster(committed, ta), raster(pc.tag(probe), tb)
        if x.shape != y.shape:
            print(f'  SHAPE  {rel}: {x.shape} vs {y.shape}'); rc = 1; continue
        n = int((x != y).any(axis=2).sum())
        print(f'  {"reproduces" if n == 0 else f"DIFFERS {n}px"}  {rel}')
        rc = rc or (0 if n == 0 else 1)
    return rc


if __name__ == '__main__':
    sys.exit(main())
