"""
build_fig5B_snapshot.py -- Figure 5B images: two raw SPT frames (RPL29-Halo /
60S), zoomed in, with every molecule present at that frame drawn as its full
track so far.

The placed 5B shows two unlabelled sgNT crops (fig5B_snapshot_spec.json). It was
drafted as sgNT vs sgUBA5, but the right-hand crop came from an NT-guide FOV, so
the condition labels were dropped (2026-10-26); the panel only illustrates raw
data. --scout still ranks candidates per condition.

Dataset `2025_06_06_zeb`, tier `q=15p0_splitsize10`, 0p015 arm (15.091 ms,
0.11 um/px) -- the same files that feed the 5C posterior overlay.

Illustrative only (one FOV per crop). Tracks are drawn in ONE color: the
panel shows that long-lived and short-lived trajectories coexist, it makes no
diffusion-class claim. Long vs short is used only to CHOOSE the crop.

Caption caveat: crops are chosen for tracks whose path spans >= 0.44 um, so
they over-represent MOBILE molecules -- most long-lived 60S are near-immobile
at this scale and draw as a dot. The image shows the kinds of track present,
not their proportions.

Two modes
---------
  --scout   score every (FOV, frame, crop) for having both long and short tracks
            present, write the top candidates per condition to a contact sheet
            (detections at T ringed, as a frame-registration check) and to
            fig5B_snapshot_candidates.json. The top-ranked pair is written as
            the default fig5B_snapshot_spec.json unless one already exists.
  --render  draw the final 1x2 panel from fig5B_snapshot_spec.json in the
            output folder, or, if there is none, from the copy shipped next to
            this script (the FOVs, frames and crops placed in Fig 5B). Edit the
            spec (stem / frame_T / crop) to swap in another candidate; an
            optional per-panel "title" and a top-level "show_time" add the
            condition label and the frame time.

The image is ONE instantaneous raw frame (no projection, no background
subtraction), so single molecules sit only ~50-100 ADU above a ~160 ADU
background. Both panels share one grey scale (percentiles of the two crops
pooled) so brightness is comparable between them.

Run:  sbatch call_build_panels.slurm build_fig5B_snapshot.py --scout
      sbatch call_build_panels.slurm build_fig5B_snapshot.py --render
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import pims
import matplotlib.patheffects as pe
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'er')))
import paperfig_style as ps
import viz_tools as vt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
EXP = tlsd_paths.spt('all_data', '2025_06_06_zeb')
TIER = os.path.join(EXP, 'tracking_output_q=15p0_splitsize10')
OUT = tlsd_paths.outp('PaperFigures_updated', 'figure5', 'B')
CAND_JSON = os.path.join(OUT, 'fig5B_snapshot_candidates.json')
SPEC_JSON = os.path.join(OUT, 'fig5B_snapshot_spec.json')
SPEC_SHIPPED = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'fig5B_snapshot_spec.json')

PIXEL_SIZE = 0.11          # um / px (ND2 metadata, METHODS table)
FRAME_INTERVAL = 0.015091  # s / frame

# Same member conditions as the sgNT_60S / sgUBA5_60S combos behind Fig 5C.
CONDITIONS = {
    'sgNT': 'halo=RPL29_cond=NT-guide_exp=0p015',
    'sgUBA5': 'halo=RPL29_cond=sgUBA5-LT_exp=0p015',
}

# Crop selection (scout only)
CROP_PX = 80               # ~8.8 um square: small enough that tails read
GRID_PX = 20               # crops are scored on a CROP_PX/4 grid
T_STRIDE = 10              # candidate frames every 10 frames
LONG_ELAPSED = 30          # "long": >= 30 frames (0.45 s) of track drawn by T
SHORT_TOTAL = 8            # "short": the whole trajectory is <= 8 frames
MIN_EXTENT_PX = 4          # a track counts only if its drawn path spans
                           # >= 4 px (0.44 um) -- most long-lived 60S are
                           # near-immobile and would draw as a bare dot
MAX_TRACKS = 10            # skip crops too crowded to read
N_CANDIDATES = 6
MIN_SEPARATION_T = 100     # candidates from one FOV must differ by >= 100 frames

# One colour for every track, not a condition colour: the magenta of the placed
# Fig 5B (#EC008C, sampled from the draft). The panel was first rendered in
# yellow #ffe600 and recoloured by hand; the legend ("past trajectory",
# "current position") was also added by hand.
TRACK_COLOR = '#EC008C'
TAIL_LW = 0.8
TIP_SIZE = 2.5
DISPLAY_PCT = (1.0, 99.8)  # shared grey scale over the two crops' pixels


def stems(condition):
    tdir = os.path.join(TIER, 'trajectories')
    # Trailing-underscore anchor so a condition never matches a sibling.
    return sorted(f[:-len('_traj.csv')] for f in os.listdir(tdir)
                  if f.startswith(condition + '_') and f.endswith('_traj.csv'))


def load_traj(stem):
    df = pd.read_csv(os.path.join(TIER, 'trajectories', f'{stem}_traj.csv'))
    df = df.dropna(subset=['x', 'y'])
    g = df.groupby('trajectory')['frame']
    df['t_start'] = g.transform('min')
    df['n_frames'] = g.transform('size')
    return df


def raw_frame(stem, frame_T):
    with pims.open(os.path.join(EXP, 'raw_data', f'{stem}.nd2')) as r:
        return np.asarray(r[frame_T], dtype=float)


def registration_offset(stem, df, frame_T):
    """Frame offset (-1, 0, +1) whose image is brightest at the detections of
    frame_T. 0 means trajectory frame indices address the ND2 directly."""
    det = df[df['frame'] == frame_T]
    yy = det['y_detect'].astype(int).values
    xx = det['x_detect'].astype(int).values
    score = {}
    with pims.open(os.path.join(EXP, 'raw_data', f'{stem}.nd2')) as r:
        for off in (-1, 0, 1):
            t = frame_T + off
            if 0 <= t < len(r):
                img = np.asarray(r[t], dtype=float)
                score[off] = float(np.mean(img[yy, xx]) - np.median(img))
    return max(score, key=score.get), score


def classify_at(df, frame_T):
    """Rows of frame_T, one per present trajectory, with `long` / `short`
    flags. Both require the path drawn up to T to span >= MIN_EXTENT_PX
    (max distance of any earlier point from the T position)."""
    at = df[df['frame'] == frame_T].set_index('trajectory')
    hist = df[(df['frame'] <= frame_T) & df['trajectory'].isin(at.index)]
    d = np.hypot(hist['x'].values - at.loc[hist['trajectory'], 'x'].values,
                 hist['y'].values - at.loc[hist['trajectory'], 'y'].values)
    extent = pd.Series(d, index=hist['trajectory'].values).groupby(level=0).max()
    at['extent'] = extent.reindex(at.index).values
    visible = at['extent'] >= MIN_EXTENT_PX
    at['long'] = visible & ((frame_T - at['t_start']) >= LONG_ELAPSED)
    at['short'] = visible & (at['n_frames'] <= SHORT_TOTAL)
    return at.reset_index()


def score_fov(stem, df, H, W):
    """Rows of (stem, T, crop, n_long, n_short, n_total) for every crop that
    holds at least one long and one short track at T."""
    ny, nx = H // GRID_PX, W // GRID_PX
    k = CROP_PX // GRID_PX
    rows = []
    for T in range(0, int(df['frame'].max()) + 1, T_STRIDE):
        if not (df['frame'] == T).any():
            continue
        at = classify_at(df, T)
        gy = np.clip((at['y'].values // GRID_PX).astype(int), 0, ny - 1)
        gx = np.clip((at['x'].values // GRID_PX).astype(int), 0, nx - 1)
        is_long = at['long'].values
        is_short = at['short'].values
        grids = []
        for m in (is_long, is_short, np.ones(len(at), bool)):
            g = np.zeros((ny, nx))
            np.add.at(g, (gy[m], gx[m]), 1)
            # k x k window sums via a 2-D cumulative sum
            c = np.pad(g.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
            grids.append(c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k])
        n_long, n_short, n_tot = grids
        ok = (n_long >= 1) & (n_short >= 1) & (n_tot <= MAX_TRACKS)
        for iy, ix in zip(*np.nonzero(ok)):
            y0, x0 = int(iy) * GRID_PX, int(ix) * GRID_PX
            rows.append(dict(stem=stem, frame_T=T,
                             crop=[y0, y0 + CROP_PX, x0, x0 + CROP_PX],
                             n_long=int(n_long[iy, ix]),
                             n_short=int(n_short[iy, ix]),
                             n_total=int(n_tot[iy, ix])))
    return rows


def rank(rows):
    """Best = most balanced long/short mix, then most tracks; keep candidates
    spread out (one per FOV-and-time neighbourhood, non-overlapping crops)."""
    rows = sorted(rows, key=lambda r: (min(r['n_long'], r['n_short']),
                                       r['n_long'] + r['n_short'],
                                       -r['n_total']), reverse=True)
    picked = []
    for r in rows:
        clash = any(p['stem'] == r['stem']
                    and abs(p['frame_T'] - r['frame_T']) < MIN_SEPARATION_T
                    for p in picked)
        if not clash:
            picked.append(r)
        if len(picked) == N_CANDIDATES:
            break
    return picked


def draw_snapshot(ax, df, img, frame_T, crop, vmin, vmax, rings=False):
    """Grey raw crop + full tracks (start..T) of molecules present at T in
    the crop. Returns (n_long, n_short, n_drawn)."""
    y0, y1, x0, x1 = crop
    ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, origin='upper',
              interpolation='nearest')
    keep = vt._present_in_crop(df, frame_T, 0, crop, 1.0)
    vt._draw_tails(ax, df, frame_T, keep, 'fixed', None, 1.0, None, None,
                   TIP_SIZE, TAIL_LW, fixed_color=TRACK_COLOR,
                   show_tip_markers=False)
    at = classify_at(df, frame_T)
    at = at[at['trajectory'].isin(keep)]
    # Small tip dot without vt's white halo, which swamps a short tail.
    ax.scatter(at['x'], at['y'], s=TIP_SIZE, color=TRACK_COLOR, lw=0,
               zorder=6)
    if rings:
        det = df[df['frame'] == frame_T]
        ax.scatter(det['x_detect'], det['y_detect'], s=30, facecolors='none',
                   edgecolors='#23AAE1', linewidths=0.6, zorder=7)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return int(at['long'].sum()), int(at['short'].sum()), len(keep)


def scale_bar(ax, crop, bar_um=2):
    y0, y1, x0, x1 = crop
    bar_px = bar_um / PIXEL_SIZE
    bx = x1 - bar_px - 0.06 * (x1 - x0)
    by = y1 - 0.07 * (y1 - y0)
    ax.plot([bx, bx + bar_px], [by, by], '-', color='white', lw=1.5,
            solid_capstyle='butt',
            path_effects=[pe.withStroke(linewidth=2.5, foreground='black')])
    ax.text(bx + bar_px / 2, by - 0.025 * (y1 - y0), f'{bar_um} µm',
            color='white', ha='center', va='bottom', fontsize=6.5,
            path_effects=[pe.withStroke(linewidth=1.5, foreground='black')])


def crop_pixels(img, crop):
    y0, y1, x0, x1 = crop
    return img[y0:y1, x0:x1].ravel()


def scout():
    os.makedirs(OUT, exist_ok=True)
    cands = {}
    for label, cond in CONDITIONS.items():
        rows = []
        for stem in stems(cond):
            df = load_traj(stem)
            with pims.open(os.path.join(EXP, 'raw_data', f'{stem}.nd2')) as r:
                H, W = r.frame_shape
            rows += score_fov(stem, df, H, W)
        cands[label] = rank(rows)
        print(f'{label}: {len(rows)} qualifying crops, top {len(cands[label])}')

    ps.apply()
    pdf_path = os.path.join(OUT, 'fig5B_snapshot_candidates.pdf')
    with PdfPages(pdf_path) as pdf:
        for label, picked in cands.items():
            fig, axes = ps.plt.subplots(2, 3, figsize=(7.5, 5.4))
            for i, (ax, c) in enumerate(zip(axes.ravel(), picked)):
                df = load_traj(c['stem'])
                off, sc = registration_offset(c['stem'], df, c['frame_T'])
                c['registration_offset'] = off
                c['registration_score'] = sc
                img = raw_frame(c['stem'], c['frame_T'])
                lo, hi = np.percentile(crop_pixels(img, c['crop']), DISPLAY_PCT)
                draw_snapshot(ax, df, img, c['frame_T'], c['crop'], lo, hi,
                              rings=True)
                scale_bar(ax, c['crop'])
                fov = c['stem'].rsplit('_', 1)[-1]
                ax.set_title(f"#{i + 1}  FOV {fov}  T={c['frame_T']}  "
                             f"crop={c['crop']}\nlong {c['n_long']} / short "
                             f"{c['n_short']} / all {c['n_total']}   "
                             f"reg. offset {off:+d}", fontsize=6)
                print(label, i + 1, c['stem'], c['frame_T'], c['crop'],
                      'reg', {k: round(v, 1) for k, v in sc.items()})
            for ax in axes.ravel()[len(picked):]:
                ax.axis('off')
            fig.suptitle(f'{label} candidates (blue rings = detections at T; '
                         f'they should sit on spots)', fontsize=8)
            fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.88,
                                wspace=0.08, hspace=0.35)
            pdf.savefig(fig)
            ps.plt.close(fig)
    with open(CAND_JSON, 'w') as f:
        json.dump(cands, f, indent=2)
    print('wrote', pdf_path, CAND_JSON)

    if not os.path.exists(SPEC_JSON):
        spec = {'display_pct': list(DISPLAY_PCT),
                'panels': {lab: {k: picked[0][k]
                                 for k in ('stem', 'frame_T', 'crop')}
                           for lab, picked in cands.items()}}
        with open(SPEC_JSON, 'w') as f:
            json.dump(spec, f, indent=2)
        print('wrote default spec (rank 1 per condition):', SPEC_JSON)


def render():
    with open(SPEC_JSON if os.path.exists(SPEC_JSON) else SPEC_SHIPPED) as f:
        spec = json.load(f)
    panels = spec['panels']

    data = {}
    for label, p in panels.items():
        df = load_traj(p['stem'])
        off, sc = registration_offset(p['stem'], df, p['frame_T'])
        assert off == 0, (f"{label}: detections at T line up with ND2 frame "
                          f"T{off:+d}, not T -- fix indexing first ({sc})")
        data[label] = (df, raw_frame(p['stem'], p['frame_T']))

    # One grey scale for both panels.
    pooled = np.concatenate([crop_pixels(img, panels[lab]['crop'])
                             for lab, (_, img) in data.items()])
    vmin, vmax = np.percentile(pooled, spec.get('display_pct', DISPLAY_PCT))

    ps.apply()
    fig, axes = ps.plt.subplots(1, 2, figsize=(3.25, 1.78))
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.86,
                        wspace=0.04)
    for ax, (label, (df, img)) in zip(axes, data.items()):
        p = panels[label]
        n_long, n_short, n_drawn = draw_snapshot(ax, df, img, p['frame_T'],
                                                 p['crop'], vmin, vmax)
        assert n_long >= 1 and n_short >= 1, (label, n_long, n_short)
        scale_bar(ax, p['crop'])
        t_s = p['frame_T'] * FRAME_INTERVAL
        if p.get('title'):                 # the placed 5B is unlabelled
            ax.set_title(p['title'], fontsize=7, pad=2,
                         color=ps.COLORS.get(p['title'], 'black'))
        if spec.get('show_time', False):
            ax.text(0.03, 0.97, f't = {t_s:.2f} s', transform=ax.transAxes,
                    ha='left', va='top', fontsize=6, color='white',
                    path_effects=[pe.withStroke(linewidth=1.5,
                                                foreground='black')])
        print(f"{label}: {p['stem']} T={p['frame_T']} crop={p['crop']} "
              f"drawn={n_drawn} long={n_long} short={n_short}")
    print(f'grey scale {vmin:.1f}-{vmax:.1f} ADU '
          f'(pooled crop percentiles {spec.get("display_pct", DISPLAY_PCT)})')

    out = os.path.join(OUT, 'fig5B_snapshot.pdf')
    ps.save(fig, out, metadata={'CreationDate': None})
    print('wrote', out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--scout', action='store_true')
    g.add_argument('--render', action='store_true')
    a = ap.parse_args()
    scout() if a.scout else render()
