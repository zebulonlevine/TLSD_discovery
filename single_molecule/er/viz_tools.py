"""
viz_tools.py — ER+SPT composite visualization + zoom-selection toolkit.

Consolidates four cell-scripts (composite_snapshot, spt_projection, comet_movie,
scout_loop) into one importable module, and adds:
  - composite_movie : two-channel movie wrapper (sweeps frame_T, builds a
                      composite per frame, mp4 -> gif -> PNG-seq fallback).
                      Replaces the old comet_movie -> comet_snapshot call.
  - build_zoom_table: per-trajectory (Level 1) and per-spatial-bin (Level 2)
                      tables for quantitatively choosing zoom crops, reading
                      the pre-computed classification + combo MLE files.

USAGE (in a notebook — survives editing without stale-def errors):

    import sys; sys.path.insert(0, '/path/holding/viz_tools')
    import importlib, viz_tools; importlib.reload(viz_tools)
    import viz_tools as vt          # call vt.composite_snapshot(...), etc.

Use the ``vt.`` prefix, NOT ``from viz_tools import *``. With ``import *`` the
names in your namespace keep pointing at the OLD functions after reload until
you re-run the import line — the exact silent staleness this module exists to
kill. With ``vt.`` prefix, ``importlib.reload(viz_tools)`` alone is enough.

Paths are passed explicitly via ExpPaths so the module is dataset-agnostic
(no module-level ANALYSIS_DIR / ER_SOURCE / MASK globals).

CROP CONVENTION: every crop in this module is (y0, y1, x0, x1) — rows first.
build_zoom_table emits crops in exactly that order, so its output pastes
straight into composite_snapshot / composite_movie with nothing to transpose.
"""

import os
import re
import glob
from dataclasses import dataclass
from collections import namedtuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
import matplotlib.patheffects as pe
import tifffile


# =============================================================================
# Constants (single source of truth)
# =============================================================================
PIXEL_SIZE = 0.11          # um / px
FRAME_INTERVAL = 0.015091  # s / frame
COORD_SCALE = 1.0          # trajectory coords are already in PIXELS

D_BREAKS = [0, 0.003, 0.04, 0.16, 0.7, 2.3, 30]
D_LABELS = ['immobile', 'very_slow', 'slow', 'medium', 'free 60S', 'free protein']
D_COLORS = {'immobile': '#3b4cc0', 'very_slow': '#6c8ed4', 'slow': '#9dc3d8',
            'medium': '#d6c89a', 'free 60S': '#e8954e', 'free protein': '#d6453b'}
# Zone palette (2026-08-02 revision). Previous scheme was
#   on_*=#1a9c1a (mid green), off_ER=#ffe600 (yellow), intermediate/ambiguous=#cccccc
# Changed because bright-mid green tracks were indistinguishable from the green
# ER-mask fill they are drawn on. Now: a LIGHT SPRING green for on-ER (reads
# against the dim green mask), grey for off_ER, and a warm tan neutral for the
# undecided classes. Tan rather than a slate blue-grey (#7a8ba3, tried and
# rejected): slate is separable from the off_ER grey in zoom panels but the
# two collapse together at full-FOV track linewidths, and the full FOV is the
# panel most likely to be printed small.
# WARNING: grey meant intermediate/ambiguous under the OLD scheme and means
# off_ER under this one. Re-render any zone figure still on the old palette
# rather than circulating both.
ZONE_COLORS = {'on_structure': '#3dfc8f', 'near_ER': '#c79a1f',
               'on_ER_broad': '#3dfc8f', 'off_ER': '#cfcfcf',
               'intermediate': '#b9a074', 'ambiguous': '#b9a074'}


def diff_class(d):
    """Map an MLE_D value to its diffusion-state label via D_BREAKS.

    Returns np.nan for NaN/None input so it survives missing MLE.
    """
    if d is None or (isinstance(d, float) and np.isnan(d)):
        return np.nan
    idx = int(np.digitize(d, D_BREAKS[1:-1]))
    idx = min(idx, len(D_LABELS) - 1)
    return D_LABELS[idx]


# =============================================================================
# Paths
# =============================================================================
@dataclass
class ExpPaths:
    """Experiment path bundle.

    Pass exp_dir (raw/preprocessed root) and analysis_dir (tracking_output_*);
    the rest are derived. Override any derived dir explicitly if your layout
    differs.
    """
    exp_dir: str
    analysis_dir: str

    @property
    def er_seg_dir(self):
        return os.path.join(self.analysis_dir, 'er_segmentation')

    @property
    def traj_dir(self):
        return os.path.join(self.analysis_dir, 'trajectories')

    @property
    def preprocessed_dir(self):
        return os.path.join(self.exp_dir, 'preprocessed')

    @property
    def preprocessed_corrected_dir(self):
        return os.path.join(self.exp_dir, 'preprocessed_corrected')

    @property
    def traj_combinations_dir(self):
        return os.path.join(self.analysis_dir, 'traj_combinations')

    @property
    def mle_combinations_dir(self):
        return os.path.join(self.analysis_dir, 'MLE_combinations')


FOVData = namedtuple('FOVData',
                     ['traj', 'er', 'mask', 'zone_dict', 'masks', 'spt_tif_path'])


# =============================================================================
# Image ops
# =============================================================================
def _norm(img, lo_pct, hi_pct, gain=1.0):
    """Percentile-normalize a grayscale image to [0,1] with a brightness gain."""
    lo, hi = np.percentile(img, [lo_pct, hi_pct])
    if hi <= lo:
        return np.zeros_like(img, dtype=float)
    out = np.clip((img.astype(float) - lo) / (hi - lo), 0, 1) * gain
    return np.clip(out, 0, 1)


def build_composite(er_image, spt_frame, mask=None,
                    er_gain=1.0, spt_gain=1.4, mask_alpha=0.30,
                    er_pct=(1, 99), spt_pct=(50, 99.9),
                    mask_color=(1.0, 1.0, 0.0)):
    """Additive RGB composite: ER->green, ribosome->magenta, mask->mask_color."""
    g = _norm(er_image, *er_pct, gain=er_gain)
    m = _norm(spt_frame, *spt_pct, gain=spt_gain)

    rgb = np.zeros((*er_image.shape, 3), dtype=float)
    rgb[..., 1] += g            # green from ER
    rgb[..., 0] += m            # red from ribosome
    rgb[..., 2] += m            # blue from ribosome (R+B = magenta)

    if mask is not None:
        mk = mask.astype(bool)
        rgb[mk, 0] += mask_alpha * mask_color[0]
        rgb[mk, 1] += mask_alpha * mask_color[1]
        rgb[mk, 2] += mask_alpha * mask_color[2]

    return np.clip(rgb, 0, 1)


# =============================================================================
# Zone / color helpers
# =============================================================================
def _zone_display(label, pool_on_near):
    """Map a raw zone label to its display class (collapses on+near if asked)."""
    if pool_on_near and label in ('on_structure', 'near_ER'):
        return 'on_ER_broad'
    return label


def _zone_color(label, pool_on_near, pooled_color=None):
    """Color for a zone label, honoring pooling + optional pooled color."""
    disp = _zone_display(label, pool_on_near)
    if disp == 'on_ER_broad' and pooled_color is not None:
        return pooled_color
    return ZONE_COLORS.get(disp, '#cccccc')


# =============================================================================
# SPT frame loader (single frame or short projection)
# =============================================================================
def _project_pages(tf, frame_T, projection=None, mode='max'):
    """Read frame_T (or a projection over [T, T+N)) from an OPEN TiffFile."""
    n_pages = len(tf.pages)
    if projection is None or projection <= 1:
        return tf.pages[frame_T].asarray().astype(float)
    hi = min(frame_T + projection, n_pages)
    stack = np.stack([tf.pages[i].asarray().astype(float)
                      for i in range(frame_T, hi)])
    return stack.mean(axis=0) if mode == 'mean' else stack.max(axis=0)


def load_spt_frame(spt_tif_path, frame_T, projection=None, mode='max'):
    """Load the magenta-channel image for a snapshot (opens + closes the tif).

    projection : None -> single instantaneous frame (faint, honest).
                 int N -> project [frame_T, frame_T+N) (lifts faint spots but
                          smears moving particles into streaks; NOT
                          instantaneous — caption accordingly).
    """
    with tifffile.TiffFile(spt_tif_path) as tf:
        return _project_pages(tf, frame_T, projection, mode)


# =============================================================================
# FOV loader
# =============================================================================
def load_fov(stem, paths, er_source='corrected', mask_name='on_ER_broad'):
    """Load + bundle the per-FOV inputs.

    Returns FOVData(traj, er, mask, zone_dict, masks, spt_tif_path):
      traj   : trajectory DataFrame (NaN x/y dropped)
      er     : ER image (corrected illum-corrected .npy, or processed deconv)
      mask   : the single boolean mask selected by mask_name (or None)
      zone_dict : {trajectory_id: er_zone_label} or None
      masks  : {'on_structure','near_ER','off_ER','on_ER_broad'} bool dict
      spt_tif_path : path to the SPT tif (caller opens it; the movie opens once)
    """
    traj = pd.read_csv(os.path.join(paths.traj_dir, f'{stem}_traj.csv')
                       ).dropna(subset=['x', 'y'])

    if er_source == 'processed':
        er = np.load(os.path.join(paths.er_seg_dir, f'{stem}_er_processed.npy'))
    else:
        er = np.load(os.path.join(paths.preprocessed_corrected_dir,
                                  f'{stem}_er.npy'))

    masks = {}
    npz_path = os.path.join(paths.er_seg_dir, f'{stem}_er_masks.npz')
    if os.path.exists(npz_path):
        with np.load(npz_path) as mm:
            on = mm['on_structure'].astype(bool)
            near = mm['near_ER'].astype(bool)
            off = mm['off_ER'].astype(bool)
        masks = {'on_structure': on, 'near_ER': near, 'off_ER': off,
                 'on_ER_broad': on | near}

    mask = masks.get(mask_name) if (mask_name and masks) else None

    zone_dict = None
    cls_path = os.path.join(paths.er_seg_dir, f'{stem}_traj_classification.csv')
    if os.path.exists(cls_path):
        cls = pd.read_csv(cls_path)
        zone_dict = dict(zip(cls['trajectory'], cls['er_zone_label']))

    spt_tif_path = os.path.join(paths.preprocessed_dir, f'{stem}_spt.tif')
    return FOVData(traj, er, mask, zone_dict, masks, spt_tif_path)


# =============================================================================
# Drawing primitives
# =============================================================================
def _draw_tails(ax, traj_df, frame_T, keep_ids, color_by,
                tail_length, coord_scale,
                traj_mle_dict, traj_zone_dict, marker_size, tail_lw,
                pool_on_near=False, pooled_color=None, show_tip_markers=True,
                fixed_color=None):
    """Draw comet tails (+ optional tip markers) for the given trajectory ids.

    color_by='fixed' draws every tail in `fixed_color` — use when the panel
    shows tracks WITHOUT a classification claim (e.g. the raw co-stain panel,
    where uniform magenta says "these are ribosome tracks", not "this one is
    on-ER"). Zone/MLE coloring belongs in a separate, explicitly labelled panel.
    """
    df = traj_df
    if color_by == 'trajectory':
        cmap = plt.get_cmap('tab20')
        cat_color = {tid: cmap(i % 20) for i, tid in enumerate(sorted(keep_ids))}

    for tid in keep_ids:
        g = df[df['trajectory'] == tid].sort_values('frame')
        lo = frame_T - tail_length if tail_length is not None else -np.inf
        seg = g[(g['frame'] <= frame_T) & (g['frame'] >= lo)]
        if len(seg) == 0:
            continue
        xs = (seg['x'] / coord_scale).values
        ys = (seg['y'] / coord_scale).values

        if color_by == 'time':
            pts = np.array([xs, ys]).T.reshape(-1, 1, 2)
            segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
            tnorm = np.linspace(0.25, 1.0, len(segs)) if len(segs) else np.array([])
            lc = LineCollection(segs, cmap='cool', array=tnorm,
                                linewidth=tail_lw, alpha=0.95)
            ax.add_collection(lc)
            tip_c = plt.get_cmap('cool')(1.0)
        else:
            if color_by == 'trajectory':
                tip_c = cat_color[tid]
            elif color_by == 'mle':
                d = (traj_mle_dict or {}).get(tid, np.nan)
                if np.isnan(d):
                    tip_c = '#dddddd'
                else:
                    cat = pd.cut([d], bins=D_BREAKS, labels=D_LABELS,
                                 include_lowest=True, right=False)[0]
                    tip_c = D_COLORS.get(str(cat), '#dddddd')
            elif color_by == 'zone':
                z = (traj_zone_dict or {}).get(tid, 'off_ER')
                tip_c = _zone_color(z, pool_on_near, pooled_color)
            elif color_by == 'fixed':
                tip_c = fixed_color or '#ffffff'
            else:
                tip_c = '#ffffff'
            ax.plot(xs, ys, '-', color=tip_c, linewidth=tail_lw, alpha=0.9,
                    solid_capstyle='round',
                    path_effects=[pe.withStroke(linewidth=tail_lw + 1.2,
                                                foreground='black', alpha=0.4)])
        if show_tip_markers:
            tip = seg.iloc[-1]
            ax.scatter([(tip['x'] / coord_scale)], [(tip['y'] / coord_scale)],
                       s=marker_size, color=tip_c, edgecolors='white',
                       linewidths=0.6, zorder=6)


def _draw_rings(ax, traj_df, frame_T, tail_length, coord_scale, traj_zone_dict,
                loc_present_only=True, loc_color_by_zone=False,
                loc_marker_size=40, pool_on_near=False, pooled_color=None,
                fixed_color=None):
    """Hollow rings at tracked localizations (frame_T, or the tail window)."""
    if loc_present_only:
        locs = traj_df[traj_df['frame'] == frame_T]
    else:
        lo = frame_T - tail_length if tail_length else frame_T
        locs = traj_df[(traj_df['frame'] <= frame_T) & (traj_df['frame'] >= lo)]
    lx = (locs['x'] / coord_scale).values
    ly = (locs['y'] / coord_scale).values
    if loc_color_by_zone and traj_zone_dict is not None:
        ring_colors = [_zone_color((traj_zone_dict or {}).get(t, 'off_ER'),
                                   pool_on_near, pooled_color)
                       for t in locs['trajectory'].values]
        ax.scatter(lx, ly, s=loc_marker_size, facecolors='none',
                   edgecolors=ring_colors, linewidths=1.1, alpha=0.95, zorder=7)
    else:
        ax.scatter(lx, ly, s=loc_marker_size, facecolors='none',
                   edgecolors=(fixed_color or '#00ffff'),
                   linewidths=0.9, alpha=0.9, zorder=7)


def _present_in_crop(traj_df, frame_T, present_window, crop, coord_scale):
    """IDs of trajectories present at T whose T-position is inside the crop."""
    df = traj_df
    pm = (df['frame'] >= frame_T - present_window) & \
         (df['frame'] <= frame_T + present_window)
    ids = df.loc[pm, 'trajectory'].unique()
    if crop is None:
        return list(ids)
    y0, y1, x0, x1 = crop
    keep = []
    for tid in ids:
        g = df[df['trajectory'] == tid]
        at_T = g.iloc[(g['frame'] - frame_T).abs().argmin()]
        xx, yy = at_T['x'] / coord_scale, at_T['y'] / coord_scale
        if (y0 <= yy < y1) and (x0 <= xx < x1):
            keep.append(tid)
    return keep


def _scale_bar(ax, x0, x1, y0, y1, pixel_size=PIXEL_SIZE, bar_um=5):
    bar_px = bar_um / pixel_size
    if (x1 - x0) <= bar_px * 1.5:
        return
    bx = x1 - bar_px - 0.05 * (x1 - x0)
    by = y1 - 0.07 * (y1 - y0)
    ax.plot([bx, bx + bar_px], [by, by], '-', color='white', linewidth=3,
            path_effects=[pe.withStroke(linewidth=5, foreground='black')])
    ax.text(bx + bar_px / 2, by - 0.02 * (y1 - y0), f'{bar_um} \u00b5m',
            color='white', ha='center', va='bottom', fontsize=9,
            path_effects=[pe.withStroke(linewidth=2, foreground='black')])


def draw_zoom_box(ax, crop, color='white', lw=1.8, ls='--', label=None):
    """Draw a rectangle marking `crop` (y0, y1, x0, x1) on a full-image axis.

    Use on the FULL-FOV view to show where a zoom panel/movie is taken from.
    Pass the same crop tuple you hand the zoom call.
    """
    from matplotlib.patches import Rectangle
    y0, y1, x0, x1 = crop
    rect = Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                     edgecolor=color, linewidth=lw, linestyle=ls, zorder=8,
                     path_effects=[pe.withStroke(linewidth=lw + 1.6,
                                                 foreground='black', alpha=0.5)])
    ax.add_patch(rect)
    if label:
        ax.text(x0, y0 - 4, label, color=color, fontsize=9, va='bottom',
                path_effects=[pe.withStroke(linewidth=2, foreground='black')])
    return rect


# =============================================================================
# Public: composite snapshot
# =============================================================================
def composite_snapshot(traj_df, er_image, spt_frame, frame_T,
                       crop=None, color_by='zone', tail_length=60,
                       mask=None, mask2=None,
                       layout='composite',     # 'composite'|'side_by_side'|'both'
                       mask_display='tint',    # 'off'|'tint'|'panel'
                       mask_color=(1.0, 1.0, 0.0),
                       mask2_color=(0.0, 1.0, 1.0),
                       er_gain=1.0, spt_gain=1.4, mask_alpha=0.30,
                       er_pct=(1, 99), spt_pct=(50, 99.9),
                       traj_mle_dict=None, traj_zone_dict=None,
                       coord_scale=COORD_SCALE, present_window=1,
                       marker_size=26, tail_lw=1.3,
                       show_localizations=False, loc_marker_size=40,
                       loc_present_only=True, loc_color_by_zone=False,
                       show_tip_markers=True,
                       pool_on_near=False, pooled_color=None,
                       zoom_box=None,
                       axes=None, title=None, pixel_size=PIXEL_SIZE):
    """Composite snapshot at frame_T. See module docstring for crop convention.

    Returns the list of Axes.
    """
    if crop is not None:
        y0, y1, x0, x1 = crop
    else:
        y0, y1, x0, x1 = 0, er_image.shape[0], 0, er_image.shape[1]

    keep_ids = _present_in_crop(traj_df, frame_T, present_window, crop, coord_scale)

    tint_mask = mask if (mask is not None and mask_display == 'tint') else None
    tint_mask2 = mask2 if (mask2 is not None and mask_display == 'tint') else None

    # Green ER panel (with optional mask tints)
    g_only = np.zeros((*er_image.shape, 3))
    g_only[..., 1] = _norm(er_image, *er_pct, gain=er_gain)
    for tm, tc in [(tint_mask, mask_color), (tint_mask2, mask2_color)]:
        if tm is not None:
            mk = tm.astype(bool)
            g_only[mk, 0] += mask_alpha * tc[0]
            g_only[mk, 1] += mask_alpha * tc[1]
            g_only[mk, 2] += mask_alpha * tc[2]
    g_only = np.clip(g_only, 0, 1)

    # Magenta ribosome panel
    m_only = np.zeros((*er_image.shape, 3))
    mc = _norm(spt_frame, *spt_pct, gain=spt_gain)
    m_only[..., 0] = mc
    m_only[..., 2] = mc

    composite = build_composite(er_image, spt_frame, tint_mask,
                                er_gain, spt_gain, mask_alpha, er_pct, spt_pct,
                                mask_color=mask_color)
    if tint_mask2 is not None:
        mk = tint_mask2.astype(bool)
        composite[mk, 0] = np.clip(composite[mk, 0] + mask_alpha * mask2_color[0], 0, 1)
        composite[mk, 1] = np.clip(composite[mk, 1] + mask_alpha * mask2_color[1], 0, 1)
        composite[mk, 2] = np.clip(composite[mk, 2] + mask_alpha * mask2_color[2], 0, 1)

    mask_panel = None
    if (mask is not None or mask2 is not None) and mask_display == 'panel':
        mask_panel = np.zeros((*er_image.shape, 3))
        for mm, cc in [(mask, mask_color), (mask2, mask2_color)]:
            if mm is not None:
                mk = mm.astype(bool)
                mask_panel[mk, 0] = cc[0]
                mask_panel[mk, 1] = cc[1]
                mask_panel[mk, 2] = cc[2]

    panels = {'composite': composite, 'ribosome': m_only, 'er': g_only,
              'mask': mask_panel}
    if layout == 'composite':
        order = ['composite']
    elif layout == 'side_by_side':
        order = ['ribosome', 'er']
    elif layout == 'both':
        order = ['composite', 'ribosome', 'er']
    else:
        raise ValueError(f'unknown layout {layout!r}')
    if mask_panel is not None:
        order = order + ['mask']

    if axes is None:
        fig, axs = plt.subplots(1, len(order),
                                figsize=(8.5 * len(order), 8.5), squeeze=False)
        axs = list(axs[0])
    else:
        axs = axes if isinstance(axes, (list, tuple)) else [axes]

    panel_titles = {'composite': 'composite', 'ribosome': 'ribosome (frame T)',
                    'er': 'ER', 'mask': 'ER mask'}
    for ax, name in zip(axs, order):
        ax.imshow(panels[name], origin='upper')
        if name in ('composite', 'ribosome'):
            _draw_tails(ax, traj_df, frame_T, keep_ids, color_by, tail_length,
                        coord_scale, traj_mle_dict, traj_zone_dict,
                        marker_size, tail_lw,
                        pool_on_near=pool_on_near, pooled_color=pooled_color,
                        show_tip_markers=show_tip_markers)
            if show_localizations:
                _draw_rings(ax, traj_df, frame_T, tail_length, coord_scale,
                            traj_zone_dict, loc_present_only, loc_color_by_zone,
                            loc_marker_size, pool_on_near, pooled_color)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(panel_titles[name], fontsize=10)
        _scale_bar(ax, x0, x1, y0, y1, pixel_size=pixel_size)
        if zoom_box is not None:
            draw_zoom_box(ax, zoom_box)

    if title:
        fig = axs[0].get_figure()
        fig.suptitle(title, y=1.0, fontsize=12)

    tail_ax = axs[0]
    if color_by == 'mle':
        from matplotlib.patches import Patch
        h = [Patch(color=D_COLORS[l], label=l) for l in D_LABELS]
        tail_ax.legend(handles=h, fontsize=7, loc='upper left', framealpha=0.8)
    elif color_by == 'zone':
        from matplotlib.patches import Patch
        seen = sorted(set(_zone_display((traj_zone_dict or {}).get(t, 'off_ER'),
                                        pool_on_near) for t in keep_ids))
        h = [Patch(color=_zone_color(z, pool_on_near, pooled_color), label=z)
             for z in seen]
        tail_ax.legend(handles=h, fontsize=7, loc='upper left', framealpha=0.8)

    return axs


# =============================================================================
# Public: ER-only panel (no ribosome channel) with optional multi-color masks
# =============================================================================
def er_panel(ax, er_image, masks=None, mask_colors=None,
             er_pct=(1, 99), er_gain=1.0, mask_alpha=0.30,
             channel='green', crop=None, zoom_box=None,
             pixel_size=PIXEL_SIZE, title=None,
             traj_df=None, frame_T=None, color_by='zone', tail_length=60,
             traj_mle_dict=None, traj_zone_dict=None,
             coord_scale=COORD_SCALE, present_window=1,
             marker_size=26, tail_lw=1.3, show_tip_markers=True,
             show_localizations=False, loc_marker_size=40,
             loc_present_only=True, loc_color_by_zone=False,
             pool_on_near=False, pooled_color=None, fixed_color=None):
    """Render ONE ER image (no ribosome channel) on `ax`, optionally tinting
    one or more boolean masks each in its own color, and OPTIONALLY overlaying
    trajectory tails/tips/rings (no raw SPT data).

    masks / mask_colors : parallel lists, e.g.
        masks=[m['on_structure'], m['near_ER']],
        mask_colors=[(1, 1, 0), (0, 1, 1)]      # yellow, cyan
      (mask_colors defaults to yellow/cyan/magenta if omitted.)
    channel : 'green' (matches the composite) or 'gray'.

    Track overlay: pass traj_df AND frame_T to draw comet tails (colored by
    `color_by`) on top of the masked ER, without the raw SPT channel — i.e.
    tracks on the ER zone map. show_tip_markers / show_localizations control
    the tip dots and localization rings as in composite_snapshot.

    For "ER with and without the overlay", call this twice into two axes —
    once with masks=None, once with masks set.
    """
    base = _norm(er_image, *er_pct, gain=er_gain)
    rgb = np.zeros((*er_image.shape, 3))
    if channel == 'gray':
        rgb[..., 0] = rgb[..., 1] = rgb[..., 2] = base
    else:
        rgb[..., 1] = base
    if masks:
        if mask_colors is None:
            mask_colors = [(1, 1, 0), (0, 1, 1), (1, 0, 1)][:len(masks)]
        for m, c in zip(masks, mask_colors):
            mk = m.astype(bool)
            for i in range(3):
                rgb[mk, i] = np.clip(rgb[mk, i] + mask_alpha * c[i], 0, 1)
    ax.imshow(np.clip(rgb, 0, 1), origin='upper')

    # Optional track overlay (tails/tips/rings) — NO raw SPT channel.
    if traj_df is not None and frame_T is not None:
        keep_ids = _present_in_crop(traj_df, frame_T, present_window, crop,
                                    coord_scale)
        _draw_tails(ax, traj_df, frame_T, keep_ids, color_by, tail_length,
                    coord_scale, traj_mle_dict, traj_zone_dict, marker_size,
                    tail_lw, pool_on_near=pool_on_near, pooled_color=pooled_color,
                    show_tip_markers=show_tip_markers, fixed_color=fixed_color)
        if show_localizations:
            _draw_rings(ax, traj_df, frame_T, tail_length, coord_scale,
                        traj_zone_dict, loc_present_only, loc_color_by_zone,
                        loc_marker_size, pool_on_near, pooled_color,
                        fixed_color=fixed_color)

    if crop is not None:
        y0, y1, x0, x1 = crop
    else:
        y0, y1, x0, x1 = 0, er_image.shape[0], 0, er_image.shape[1]
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_xticks([]); ax.set_yticks([])
    _scale_bar(ax, x0, x1, y0, y1, pixel_size=pixel_size)
    if zoom_box is not None:
        draw_zoom_box(ax, zoom_box)
    if title:
        ax.set_title(title, fontsize=10)
    return ax


# =============================================================================
# Public: mask-only panel (colored zones, NO grayscale/green ER base image)
# =============================================================================
def mask_only_panel(ax, shape, masks, mask_colors=None, background=(0, 0, 0),
                     crop=None, zoom_box=None, pixel_size=PIXEL_SIZE, title=None,
                     traj_df=None, frame_T=None, color_by='zone', tail_length=60,
                     traj_mle_dict=None, traj_zone_dict=None,
                     coord_scale=COORD_SCALE, present_window=1,
                     marker_size=26, tail_lw=1.3, show_tip_markers=True,
                     show_localizations=False, loc_marker_size=40,
                     loc_present_only=True, loc_color_by_zone=False,
                     pool_on_near=False, pooled_color=None, fixed_color=None,
                     spt_frame=None, spt_pct=(95, 99.95), spt_gain=1.5,
                     spt_channel='magenta'):
    """Render ONLY the colored zone mask(s) on `ax` -- no ER image underneath
    at all (unlike er_panel's `masks=`, which tints a grayscale/green ER
    base). `shape` is the (H, W) image shape (e.g. `er_image.shape`), used
    only to size the canvas -- the ER image itself is never read.

    masks / mask_colors : parallel lists, e.g.
        masks=[m['on_structure'], m['near_ER']],
        mask_colors=[(0.10, 0.61, 0.10), (0, 1, 1)]   # green, cyan
      (mask_colors defaults to yellow/cyan/magenta if omitted, matching
      er_panel's default.) Pixels outside every mask are `background`
      (default solid black).

    spt_frame : optional raw ribosome image (single frame or projection) added
      ADDITIVELY on top of the mask canvas, in `spt_channel` ('magenta' or
      'gray'). Use to show the mask together with the raw particles it is
      classifying, rather than the mask as a bare cartoon. The mask stays the
      green channel and magenta stays the ribosome channel, so this panel
      keeps the same color grammar as the green/magenta co-stain merge.

    Same optional track overlay (tails/tips/rings) as er_panel, drawn on
    top of the mask-only canvas. Pass color_by='fixed' + fixed_color for
    uniformly-colored tracks (no per-track classification claim).
    """
    rgb = np.zeros((*shape, 3))
    rgb[..., 0], rgb[..., 1], rgb[..., 2] = background
    if masks:
        if mask_colors is None:
            mask_colors = [(1, 1, 0), (0, 1, 1), (1, 0, 1)][:len(masks)]
        for m, c in zip(masks, mask_colors):
            mk = m.astype(bool)
            for i in range(3):
                rgb[mk, i] = np.clip(c[i], 0, 1)
    if spt_frame is not None:
        s = _norm(spt_frame, *spt_pct, gain=spt_gain)
        if spt_channel == 'magenta':
            rgb[..., 0] += s
            rgb[..., 2] += s
        else:
            rgb[..., 0] += s
            rgb[..., 1] += s
            rgb[..., 2] += s
    ax.imshow(np.clip(rgb, 0, 1), origin='upper')

    if traj_df is not None and frame_T is not None:
        keep_ids = _present_in_crop(traj_df, frame_T, present_window, crop,
                                    coord_scale)
        _draw_tails(ax, traj_df, frame_T, keep_ids, color_by, tail_length,
                    coord_scale, traj_mle_dict, traj_zone_dict, marker_size,
                    tail_lw, pool_on_near=pool_on_near, pooled_color=pooled_color,
                    show_tip_markers=show_tip_markers, fixed_color=fixed_color)
        if show_localizations:
            _draw_rings(ax, traj_df, frame_T, tail_length, coord_scale,
                        traj_zone_dict, loc_present_only, loc_color_by_zone,
                        loc_marker_size, pool_on_near, pooled_color,
                        fixed_color=fixed_color)

    if crop is not None:
        y0, y1, x0, x1 = crop
    else:
        y0, y1, x0, x1 = 0, shape[0], 0, shape[1]
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_xticks([]); ax.set_yticks([])
    _scale_bar(ax, x0, x1, y0, y1, pixel_size=pixel_size)
    if zoom_box is not None:
        draw_zoom_box(ax, zoom_box)
    if title:
        ax.set_title(title, fontsize=10)
    return ax


# =============================================================================
# Public: SPT-only panel (ribosome channel + tails, grayscale by default)
# =============================================================================
def spt_panel(ax, traj_df, spt_frame, frame_T, crop=None,
              color_by='zone', tail_length=60,
              traj_mle_dict=None, traj_zone_dict=None,
              spt_pct=(99.5, 99.99), spt_gain=1.0, channel='gray',
              coord_scale=COORD_SCALE, present_window=1,
              marker_size=26, tail_lw=1.3,
              show_localizations=False, loc_marker_size=40,
              loc_present_only=True, loc_color_by_zone=False,
              show_tip_markers=True, pool_on_near=False, pooled_color=None,
              fixed_color=None,
              zoom_box=None, pixel_size=PIXEL_SIZE, title=None):
    """Render the ribosome SPT channel alone (+ comet tails) on `ax`.

    channel : 'gray' (default — for split-out panels) or 'magenta' (matches
    the composite). Tails/rings still color by `color_by` regardless.
    """
    base = _norm(spt_frame, *spt_pct, gain=spt_gain)
    rgb = np.zeros((*spt_frame.shape, 3))
    if channel == 'magenta':
        rgb[..., 0] = base
        rgb[..., 2] = base
    else:
        rgb[..., 0] = rgb[..., 1] = rgb[..., 2] = base
    ax.imshow(rgb, origin='upper')

    keep_ids = _present_in_crop(traj_df, frame_T, present_window, crop,
                                coord_scale)
    _draw_tails(ax, traj_df, frame_T, keep_ids, color_by, tail_length,
                coord_scale, traj_mle_dict, traj_zone_dict, marker_size, tail_lw,
                pool_on_near=pool_on_near, pooled_color=pooled_color,
                show_tip_markers=show_tip_markers, fixed_color=fixed_color)
    if show_localizations:
        _draw_rings(ax, traj_df, frame_T, tail_length, coord_scale,
                    traj_zone_dict, loc_present_only, loc_color_by_zone,
                    loc_marker_size, pool_on_near, pooled_color,
                    fixed_color=fixed_color)

    if crop is not None:
        y0, y1, x0, x1 = crop
    else:
        y0, y1, x0, x1 = 0, spt_frame.shape[0], 0, spt_frame.shape[1]
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_xticks([]); ax.set_yticks([])
    _scale_bar(ax, x0, x1, y0, y1, pixel_size=pixel_size)
    if zoom_box is not None:
        draw_zoom_box(ax, zoom_box)
    if title:
        ax.set_title(title, fontsize=10)
    return ax


# =============================================================================
# Public: composite MOVIE (replaces comet_movie's comet_snapshot call)
# =============================================================================
def composite_movie(traj_df, er_image, spt_tif_path, frames, out_path,
                    crop=None, color_by='zone', tail_length=60,
                    mask=None, mask2=None, mask_display='tint',
                    mask_color=(1.0, 1.0, 0.0), mask2_color=(0.0, 1.0, 1.0),
                    traj_mle_dict=None, traj_zone_dict=None,
                    er_pct=(1, 99), er_gain=1.0,
                    spt_pct=(99.5, 99.99), spt_gain=1.0,
                    projection=None, projection_mode='max',
                    show_localizations=True, loc_marker_size=60,
                    loc_present_only=True, loc_color_by_zone=True,
                    present_window=1, show_tip_markers=True,
                    pool_on_near=False, pooled_color=None,
                    coord_scale=COORD_SCALE, pixel_size=PIXEL_SIZE,
                    zoom_box=None,
                    fps=15, dpi=120, figsize=(9, 9),
                    title_prefix='', progress_every=25):
    """Two-channel composite movie over a frame sequence.

    DEFAULTS encode the agreed video strategy: INSTANTANEOUS magenta channel
    (projection=None) with raw haze suppressed (spt_pct=(99.5,99.99)) plus
    localization rings (show_localizations=True). Rationale: a movie already
    advances in time on playback, so a within-frame projection would smear
    motion a second time and the streak would LEAD the true position (it
    integrates future frames T..T+N). Rings at exactly T are unambiguous; the
    tails carry trajectory history.

    `projection` is available but discouraged for movies; if detections strobe
    frame-to-frame, raise `present_window` first (anti-flicker) before
    reaching for a short projection.

    color_by 'zone'/'mle'/'time' are stable frame-to-frame; 'trajectory'
    colors may shift as the present-set changes.

    out_path extension chooses format: .mp4 (ffmpeg) -> .gif (Pillow) ->
    directory of PNGs. Returns the path actually written.

    One figure/axes is reused and cleared each frame, and the SPT tif is opened
    ONCE for the whole sweep.
    """
    frames = list(frames)
    if not frames:
        raise ValueError("frames is empty")

    fig, ax = plt.subplots(figsize=figsize)
    tf = tifffile.TiffFile(spt_tif_path)

    def draw(frame_T):
        ax.clear()
        spt = _project_pages(tf, frame_T, projection, projection_mode)
        composite_snapshot(
            traj_df, er_image, spt, frame_T,
            crop=crop, color_by=color_by, tail_length=tail_length,
            mask=mask, mask2=mask2, mask_display=mask_display,
            mask_color=mask_color, mask2_color=mask2_color,
            er_pct=er_pct, er_gain=er_gain, spt_pct=spt_pct, spt_gain=spt_gain,
            traj_mle_dict=traj_mle_dict, traj_zone_dict=traj_zone_dict,
            coord_scale=coord_scale, present_window=present_window,
            show_localizations=show_localizations,
            loc_marker_size=loc_marker_size, loc_present_only=loc_present_only,
            loc_color_by_zone=loc_color_by_zone,
            show_tip_markers=show_tip_markers,
            pool_on_near=pool_on_near, pooled_color=pooled_color,
            layout='composite', axes=[ax], pixel_size=pixel_size,
            zoom_box=zoom_box,
            title=f'{title_prefix}frame {frame_T}')

    def _run(writer, path):
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        with writer.saving(fig, path, dpi=dpi):
            for i, fT in enumerate(frames):
                draw(fT)
                writer.grab_frame()
                if progress_every and i % progress_every == 0:
                    print(f'  frame {i+1}/{len(frames)} (T={fT})')

    try:
        ext = os.path.splitext(out_path)[1].lower()

        if ext == '.mp4':
            if animation.writers.is_available('ffmpeg'):
                _run(animation.FFMpegWriter(fps=fps, bitrate=4000), out_path)
                print(f'Wrote mp4: {out_path}')
                return out_path
            print('  ffmpeg not available; falling back to GIF.')
            out_path = os.path.splitext(out_path)[0] + '.gif'
            ext = '.gif'

        if ext == '.gif':
            if animation.writers.is_available('pillow'):
                _run(animation.PillowWriter(fps=fps), out_path)
                print(f'Wrote gif: {out_path}')
                return out_path
            print('  Pillow writer not available; falling back to PNG sequence.')
            out_path = os.path.splitext(out_path)[0] + '_frames'

        # PNG sequence fallback
        seq_dir = out_path
        os.makedirs(seq_dir, exist_ok=True)
        for i, fT in enumerate(frames):
            draw(fT)
            fig.savefig(os.path.join(seq_dir, f'frame_{i:05d}_T{fT:05d}.png'),
                        dpi=dpi, bbox_inches='tight')
            if progress_every and i % progress_every == 0:
                print(f'  frame {i+1}/{len(frames)} (T={fT})')
        print(f'Wrote {len(frames)} PNGs to: {seq_dir}')
        print('Assemble with ffmpeg, e.g.:')
        print(f"  ffmpeg -framerate {fps} -pattern_type glob "
              f"-i '{seq_dir}/frame_*.png' -c:v libx264 -pix_fmt yuv420p "
              f'{seq_dir}.mp4')
        return seq_dir
    finally:
        tf.close()
        plt.close(fig)


# =============================================================================
# Combo MLE join (pooled_id -> per-FOV original trajectory id -> MLE_D)
# =============================================================================
def _load_combo_mle(paths, combo):
    """Return {(stem, original_traj_id): (median_MLE_D, n_subtracks)} for one
    ALL-zone combo, by joining its pooled MLE file to its fov_map sidecar.

    Uses contiguous per-FOV pooled-id blocks: for each fov_basename the pooled
    ids form a contiguous run, so original_id = pooled_id - min(pooled_id in
    that FOV). For the all-zone combos that original_id equals the canonical
    per-FOV trajectory id, which is also the classification CSV key.
    """
    mle_path = os.path.join(paths.mle_combinations_dir, f'{combo}_pooled_MLE.csv')
    fov_path = os.path.join(paths.traj_combinations_dir, f'{combo}_fov_map.csv')
    if not (os.path.exists(mle_path) and os.path.exists(fov_path)):
        print(f"  [MLE] combo '{combo}': missing "
              f"{'MLE' if not os.path.exists(mle_path) else 'fov_map'} file; "
              f"skipping.")
        return {}

    mle = pd.read_csv(mle_path)
    mle['pooled_id'] = mle['orig_trajectory'].str.extract(r'::(\d+)$').astype(int)

    fov_map = pd.read_csv(fov_path)
    # offset per FOV = min pooled id in that FOV's contiguous block
    offset = fov_map.groupby('fov_basename')['trajectory'].min().to_dict()
    pooled_to_fov = dict(zip(fov_map['trajectory'], fov_map['fov_basename']))

    mle['fov_basename'] = mle['pooled_id'].map(pooled_to_fov)
    mle = mle.dropna(subset=['fov_basename'])
    mle['offset'] = mle['fov_basename'].map(offset)
    mle['orig_id'] = (mle['pooled_id'] - mle['offset']).astype(int)
    mle['stem'] = mle['fov_basename'].str.replace(r'_traj\.csv$', '', regex=True)

    out = {}
    grp = mle.groupby(['stem', 'orig_id'])['MLE_D']
    med = grp.median()
    cnt = grp.size()
    for key in med.index:
        out[key] = (float(med[key]), int(cnt[key]))
    return out


def _mle_map_for_combos(paths, mle_combos):
    """Merge per-combo MLE maps into one {(stem, orig_id): (MLE_D, n_subtracks)}."""
    if not mle_combos:
        return {}
    merged = {}
    for combo in mle_combos:
        merged.update(_load_combo_mle(paths, combo))
    return merged


# =============================================================================
# Zoom-selection tables (Level 1 per-trajectory, Level 2 per-spatial-bin)
# =============================================================================
def zone_diff_predicate(zones=None, diffs=None):
    """Build a population predicate f(zone_label, diff_class) -> bool.

    zones / diffs : optional allow-lists. A trajectory matches if its zone is
    in `zones` (when given) AND its diff_class is in `diffs` (when given).
    """
    zset = set(zones) if zones else None
    dset = set(diffs) if diffs else None

    def f(zone_label, dclass):
        if zset is not None and zone_label not in zset:
            return False
        if dset is not None:
            if dclass is None or (isinstance(dclass, float) and np.isnan(dclass)):
                return False
            if dclass not in dset:
                return False
        return True
    return f


def _fov_shape(paths, stem):
    """Cheap (H, W) for a FOV via memmap; falls back to None if unavailable."""
    for p in (os.path.join(paths.preprocessed_corrected_dir, f'{stem}_er.npy'),
              os.path.join(paths.er_seg_dir, f'{stem}_er_processed.npy')):
        if os.path.exists(p):
            try:
                return np.load(p, mmap_mode='r').shape[:2]
            except Exception:
                pass
    return None


def _discover_stems(paths, conditions=None):
    stems = [os.path.basename(f)[:-len('_traj.csv')]
             for f in sorted(glob.glob(os.path.join(paths.traj_dir,
                                                     '*_traj.csv')))]
    if conditions:
        stems = [s for s in stems if any(s.startswith(c) for c in conditions)]
    return stems


def _per_traj_for_stem(stem, paths, mle_map, stride):
    """Level-1 rows for one FOV; also returns the localization df for binning."""
    traj = pd.read_csv(os.path.join(paths.traj_dir, f'{stem}_traj.csv')
                       ).dropna(subset=['x', 'y'])
    if len(traj) == 0:
        return pd.DataFrame(), traj

    agg = traj.groupby('trajectory').agg(
        n_frames=('frame', 'count'),
        frame_start=('frame', 'min'), frame_end=('frame', 'max'),
        x_min=('x', 'min'), x_max=('x', 'max'),
        y_min=('y', 'min'), y_max=('y', 'max'),
        x_cen=('x', 'mean'), y_cen=('y', 'mean'),
    ).reset_index()
    agg['span'] = np.maximum(agg['x_max'] - agg['x_min'],
                             agg['y_max'] - agg['y_min'])
    agg['stem'] = stem

    # classification merge
    cls_path = os.path.join(paths.er_seg_dir, f'{stem}_traj_classification.csv')
    if os.path.exists(cls_path):
        cls = pd.read_csv(cls_path)
        keep = [c for c in ['trajectory', 'er_zone_label', 'frac_on',
                            'frac_near', 'enrichment_p', 'enrichment_p_near',
                            'n_on', 'n_near'] if c in cls.columns]
        agg = agg.merge(cls[keep], on='trajectory', how='left')
    if 'er_zone_label' not in agg.columns:
        agg['er_zone_label'] = np.nan

    # MLE merge
    agg['MLE_D'] = [mle_map.get((stem, int(t)), (np.nan, 0))[0]
                    for t in agg['trajectory']]
    agg['n_subtracks'] = [mle_map.get((stem, int(t)), (np.nan, 0))[1]
                          for t in agg['trajectory']]
    agg['diff_class'] = agg['MLE_D'].map(diff_class)

    # home bin (centroid)
    agg['hb_row'] = (agg['y_cen'] // stride).astype(int)
    agg['hb_col'] = (agg['x_cen'] // stride).astype(int)
    return agg, traj


def build_zoom_table(paths, stems=None, conditions=None, mle_combos=None,
                     bin_size=120, stride=None, min_dwell=4,
                     populations=None, exclude_ambiguous=True):
    """Build per-trajectory (Level 1) and per-spatial-bin (Level 2) tables for
    quantitative zoom selection.

    Parameters
    ----------
    paths : ExpPaths
    stems : list of FOV stems, or None to discover all in traj_dir.
    conditions : optional list; keeps stems starting with any of them.
    mle_combos : list of ALL-zone combo names (e.g. ['Degron','noDegron'])
        providing MLE_D via the pooled-MLE + fov_map join. None -> zone-only
        (MLE_D / diff_class become NaN). Pass the all-zone combos, NOT the
        per-zone ones, so original_id reconstruction matches the canonical
        per-FOV trajectory ids. NOTE: these D's carry the joint-fit ~+0.16
        log10 offset — valid for RELATIVE zoom ranking, not as per-FOV MLE_D.
    bin_size : square bin (px) = the crop size you intend to use.
    stride : bin step (px); None -> stride=bin_size (non-overlapping). Set
        stride < bin_size for overlapping bins that catch edge-straddlers.
    min_dwell : a trajectory counts toward a bin only if >= this many of its
        localizations fall inside the bin (curbs corner-clip inflation from
        long fast tracks; also a legibility filter).
    populations : dict {name: predicate(zone_label, diff_class)->bool}, e.g.
        {'A': zone_diff_predicate(['on_structure'], ['immobile','confined']),
         'B': zone_diff_predicate(['off_ER'], ['fast','super-fast'])}.
        contrast_score = min over the populations' per-bin counts (high only
        when ALL requested populations are present together). None -> no
        contrast_score / cooccurrence columns.
    exclude_ambiguous : drop 'ambiguous' (and NaN-label) trajectories from the
        counts/score so you zoom on confident calls. n_ambiguous is reported
        per bin regardless, so power-dead regions are visible.

    Returns
    -------
    (level1_df, level2_df)
    """
    if stride is None:
        stride = bin_size
    if stems is None:
        stems = _discover_stems(paths, conditions)
    elif conditions:
        stems = [s for s in stems if any(s.startswith(c) for c in conditions)]

    mle_map = _mle_map_for_combos(paths, mle_combos)
    pop_names = list(populations.keys()) if populations else []

    level1_parts = []
    level2_rows = []

    for stem in stems:
        per_traj, locs = _per_traj_for_stem(stem, paths, mle_map, stride)
        if len(per_traj) == 0:
            continue
        level1_parts.append(per_traj)

        # per-traj attribute lookups for binning
        zone_of = dict(zip(per_traj['trajectory'], per_traj['er_zone_label']))
        diff_of = dict(zip(per_traj['trajectory'], per_traj['diff_class']))
        nsub_of = dict(zip(per_traj['trajectory'], per_traj['n_subtracks']))
        pop_of = {name: {t for t in per_traj['trajectory']
                         if populations[name](zone_of.get(t), diff_of.get(t))}
                  for name in pop_names}

        shape = _fov_shape(paths, stem)
        if shape is not None:
            H, W = shape
        else:
            H = int(np.ceil(locs['y'].max())) + 1
            W = int(np.ceil(locs['x'].max())) + 1

        xv = locs['x'].values
        yv = locs['y'].values
        tv = locs['trajectory'].values
        fv = locs['frame'].values

        for y0 in range(0, H, stride):
            y1 = min(y0 + bin_size, H)
            for x0 in range(0, W, stride):
                x1 = min(x0 + bin_size, W)
                inb = (yv >= y0) & (yv < y1) & (xv >= x0) & (xv < x1)
                if not inb.any():
                    continue
                bt = tv[inb]
                # dwell count per trajectory in this bin
                uniq, counts = np.unique(bt, return_counts=True)
                touching = uniq[counts >= min_dwell]
                if len(touching) == 0:
                    continue

                labels = [zone_of.get(t) for t in touching]
                n_ambiguous = sum(1 for l in labels
                                  if l == 'ambiguous' or l is None
                                  or (isinstance(l, float) and np.isnan(l)))

                if exclude_ambiguous:
                    consider = [t for t, l in zip(touching, labels)
                                if not (l == 'ambiguous' or l is None
                                        or (isinstance(l, float) and np.isnan(l)))]
                else:
                    consider = list(touching)

                row = {'stem': stem,
                       'bin_y0': y0, 'bin_y1': y1, 'bin_x0': x0, 'bin_x1': x1,
                       'crop': (y0, y1, x0, x1),
                       'n_traj_touching': len(touching),
                       'n_ambiguous': n_ambiguous,
                       'total_locs': int(inb.sum())}

                czone = {z: 0 for z in ('on_structure', 'near_ER', 'off_ER')}
                czone_sub = {z: 0 for z in ('on_structure', 'near_ER', 'off_ER')}
                cdiff = {d: 0 for d in D_LABELS}
                mle_vals = []
                for t in consider:
                    z = zone_of.get(t)
                    if z in czone:
                        czone[z] += 1
                        czone_sub[z] += int(nsub_of.get(t, 0) or 0)
                    d = diff_of.get(t)
                    if isinstance(d, str) and d in cdiff:
                        cdiff[d] += 1
                    mv = per_traj.loc[per_traj['trajectory'] == t, 'MLE_D']
                    if len(mv) and not np.isnan(mv.values[0]):
                        mle_vals.append(float(mv.values[0]))

                for z in czone:
                    row[f'n_{z}'] = czone[z]
                    row[f'nsub_{z}'] = czone_sub[z]
                for d in D_LABELS:
                    row[f'n_{d.replace("-", "_").replace(" ", "_")}'] = cdiff[d]

                row['dominant_zone'] = (max(czone, key=czone.get)
                                        if sum(czone.values()) else np.nan)
                row['dominant_diff_class'] = (max(cdiff, key=cdiff.get)
                                              if sum(cdiff.values()) else np.nan)
                row['med_MLE_D'] = float(np.median(mle_vals)) if mle_vals else np.nan
                if len(mle_vals) >= 2:
                    q75, q25 = np.percentile(mle_vals, [75, 25])
                    row['iqr_MLE_D'] = float(q75 - q25)
                else:
                    row['iqr_MLE_D'] = np.nan

                # contrast + cooccurrence window
                if populations:
                    consider_set = set(consider)
                    pop_counts = {}
                    pop_frames = {}
                    for name in pop_names:
                        members = pop_of[name] & consider_set
                        pop_counts[name] = len(members)
                        # in-bin frames for this population's members
                        msk = inb & np.isin(tv, list(members)) if members else \
                            np.zeros_like(inb)
                        pf = fv[msk]
                        pop_frames[name] = (int(pf.min()), int(pf.max())) \
                            if len(pf) else None
                        row[f'count_{name}'] = pop_counts[name]
                    row['contrast_score'] = min(pop_counts.values()) \
                        if pop_counts else 0
                    # co-occurrence window: overlap of every population's range
                    if all(pop_frames[n] is not None for n in pop_names) and pop_names:
                        lo = max(pop_frames[n][0] for n in pop_names)
                        hi = min(pop_frames[n][1] for n in pop_names)
                        if lo <= hi:
                            row['cooccur_frame_start'] = lo
                            row['cooccur_frame_end'] = hi
                        else:
                            row['cooccur_frame_start'] = np.nan
                            row['cooccur_frame_end'] = np.nan
                    else:
                        row['cooccur_frame_start'] = np.nan
                        row['cooccur_frame_end'] = np.nan

                level2_rows.append(row)

    level1_df = (pd.concat(level1_parts, ignore_index=True)
                 if level1_parts else pd.DataFrame())
    level2_df = pd.DataFrame(level2_rows)
    if populations and len(level2_df):
        level2_df = level2_df.sort_values('contrast_score', ascending=False
                                          ).reset_index(drop=True)
    return level1_df, level2_df


def top_crops(level2_df, k=10, by='contrast_score'):
    """Return the top-k bins as paste-ready crops.

    Output columns: stem, crop (y0,y1,x0,x1), crop_str, the ranking value, and
    the co-occurrence frame window if present. Also prints them.
    """
    if len(level2_df) == 0:
        print('(no bins)')
        return level2_df
    cols = ['stem', 'crop', by]
    for extra in ('cooccur_frame_start', 'cooccur_frame_end',
                  'dominant_zone', 'dominant_diff_class', 'n_traj_touching'):
        if extra in level2_df.columns:
            cols.append(extra)
    out = level2_df.sort_values(by, ascending=False).head(k)[cols].copy()
    out['crop_str'] = out['crop'].map(lambda c: f'({c[0]}, {c[1]}, {c[2]}, {c[3]})')
    print(f'Top {len(out)} bins by {by} (crop = y0, y1, x0, x1):')
    for _, r in out.iterrows():
        extra = ''
        if 'cooccur_frame_start' in out.columns and \
                not pd.isna(r.get('cooccur_frame_start')):
            extra = (f"  frames {int(r['cooccur_frame_start'])}-"
                     f"{int(r['cooccur_frame_end'])}")
        print(f"  {r['stem']}  crop={r['crop_str']}  {by}={r[by]}{extra}")
    return out
