"""
build_er_costain_figure.py
==========================

Main-text ER / ribosome co-stain + masking figure (Fig 5E: FOV 001, the full-FOV
merge and the zoom-1 merge, mask + tracks and classification panels).

This script builds the explanatory panel that walks a reader from the two raw channels to the ER mask used by the zone
classification:

  col 1  ER                    green channel alone
  col 2  ribosome (SPT)        magenta channel alone, max-projected
  col 3  merge                 green ER + magenta ribosome, NO mask
  col 4  ER mask + tracks      on_ER_broad mask (green) + the raw magenta
                               particles + comet tails, ALL TAILS ONE COLOUR
                               (magenta) -- this panel makes no per-track
                               classification claim, it just shows that tracks
                               and mask are being compared in the same frame
  col 5  zone classification   same mask + particles, tails now coloured by
                               er_zone_label (blue = on_ER_broad,
                               red = off_ER, purple = intermediate).
                               `ambiguous` tracks are OMITTED here -- see
                               EXCLUDE_ZONES.

Splitting 4 and 5 is deliberate: the classification colours are a *result*,
so they get their own panel rather than being smuggled into the panel that
introduces the mask.

Rows: row 0 = full FOV (dashed boxes mark the zooms); rows 1.. = zoom crops.

DISPLAY SCALING. Both channels are normalized ONCE into display arrays that
every panel then shares (passed through with pct=(0,100), gain=1). This
guarantees the green in col 1 is the same green in cols 3-5 -- a per-panel
re-normalization would silently change the apparent ER brightness between the
single-channel and merge panels. ER additionally gets a gamma (>1) because
its intensity is broadly distributed (p25~387, p75~525 within a 297-1108
range), so a plain percentile stretch renders as a flat green wash instead of
visible ER structure -- this is why the older notebook composite looked
saturated.

Crops: taken from CROPS below when present (FOV 001 reuses the two crops the
existing figure4_claim3 panels already use, so the new figure lines up with
what is already in the deck). Any FOV without an entry falls back to
viz_tools.build_zoom_table contrast_score selection.
"""

import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'er')))
import viz_tools as vt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
OUT = tlsd_paths.outp('PaperFigures', 'figure4_claim3', 'costain')
os.makedirs(OUT, exist_ok=True)

paths = vt.ExpPaths(
    exp_dir=tlsd_paths.spt('all_data', '2025_09_15_zeb_ER'),
    analysis_dir=os.path.join(tlsd_paths.spt('all_data', '2025_09_15_zeb_ER'), 'tracking_output_q=15p0'),
)

FOVS = [
    ('line=L29_degron=000_HHT=000_001', 'noDegron, HHT=000'),
    ('line=L29_degron=000_HHT=005_t=000_023', 'noDegron, HHT=005 (treated)'),
    ('line=L29_degron=001_HHT=000_017', 'Degron, HHT=000'),
    ('line=L29_degron=001_HHT=005_t=000_037', 'Degron, HHT=005 (treated)'),
]

# Fixed crops for the placed FOV (originally picked for the FOV QC grids).
# Omit a stem to auto-pick via build_zoom_table.
CROPS = {
    'line=L29_degron=000_HHT=000_001': [(0, 160, 160, 320), (320, 480, 0, 160)],
}

MLE_COMBOS = ['noDegron', 'Degron']
BIN_SIZE = 160
N_ZOOMS = 2

FRAME_T = 300
TAIL_LENGTH = 80
SNAP_PROJ_N, SNAP_PROJ_MODE = 8, 'max'

# --- display scaling (shared by every panel; see module docstring) ----------
ER_PCT, ER_GAMMA = (1, 99.9), 1.6
SPT_PCT, SPT_GAIN = (95, 99.95), 1.5

# --- colours ---------------------------------------------------------------
# Mask fill is a DIM green: the bright ZONE_COLORS green is reserved for
# on_ER_broad tracks in col 5, and bright-green tracks on a bright-green mask
# is the low-contrast problem flagged in docs/HANDOFF_er_viz.md.
MASK_COLOR = (0.06, 0.34, 0.06)
TRACE_MAGENTA = '#ff45d0'
# Classification colours of the published Fig 5E (on ER blue, off ER red,
# intermediate = "undetermined" purple). The placed panel was recoloured from
# the viz_tools palette by hand; this sets the same colours at render time.
# Only this script's process sees the change.
vt.ZONE_COLORS.update({'on_structure': '#23AAE1', 'on_ER_broad': '#23AAE1',
                       'off_ER': '#ED1C24', 'intermediate': '#662D91'})
ON_BROAD_COLOR = vt.ZONE_COLORS['on_ER_broad']      # blue
OFF_ER_COLOR = vt.ZONE_COLORS['off_ER']             # red
ZOOM_BOX_COLORS = ['white', '#00d0ff']

# Zones dropped from the col-5 classification panel ONLY. `ambiguous` is the
# untestable-short class: in FOV 001 it is 54% of trajectories but just 9.7%
# of localizations (14,812 / 152,862), so omitting it removes half the tracks
# while costing ~a tenth of observed particle-time, and stops an unclassified
# majority from swamping the panel whose job is to show the on/off-ER split.
# Col 4 still draws EVERY track, so the pair reads as
# "all tracks" -> "the classifiable subset". Caption it that way.
EXCLUDE_ZONES = ('ambiguous',)

SHOW_TITLES = True      # False for the final main-text version (label in AI)


def pick_zoom_crops(stem, n=N_ZOOMS):
    """contrast_score crop selection (only for FOVs missing a CROPS entry)."""
    pops = {
        'A': vt.zone_diff_predicate(zones=['on_structure'],
                                    diffs=['immobile', 'very_slow', 'slow']),
        'B': vt.zone_diff_predicate(zones=['off_ER'],
                                    diffs=['free 60S', 'free protein']),
    }
    l1, l2 = vt.build_zoom_table(
        paths, stems=[stem], mle_combos=MLE_COMBOS,
        bin_size=BIN_SIZE, stride=BIN_SIZE, min_dwell=2,
        populations=pops, exclude_ambiguous=True)
    top = vt.top_crops(l2, k=n, by='contrast_score')
    return list(top['crop']) if len(top) else []


def _title(text):
    return text if SHOW_TITLES else None


def drop_zones(traj, zone_dict, zones=EXCLUDE_ZONES):
    """Trajectory rows whose er_zone_label is not in `zones`.

    Applied to the classification panel only. load_fov guarantees every
    trajectory id is present in zone_dict (verified on this dataset), so an
    unmapped id would be a real upstream problem rather than something to
    silently pass through as the _zone_color 'off_ER' default -- hence
    dropping NaN here too.
    """
    if not zones:
        return traj
    lab = traj['trajectory'].map(zone_dict)
    return traj[lab.notna() & ~lab.isin(zones)]


def panel_row(axes_row, traj, er_disp, spt_disp, on_broad, zone_dict,
              crop, show_rings, present_window, zoom_boxes=None,
              row_label=None):
    """Render one 5-column row for `crop` (None = full FOV)."""
    prefix = f'{row_label} | ' if row_label else ''

    # 1 — ER alone (green)
    vt.er_panel(axes_row[0], er_disp, channel='green', crop=crop,
                er_pct=(0, 100), er_gain=1.0,
                title=_title(f'{prefix}ER'))

    # 2 — ribosome alone (magenta), no tracks
    vt.spt_panel(axes_row[1], traj, spt_disp, FRAME_T, crop=crop,
                 channel='magenta', spt_pct=(0, 100), spt_gain=1.0,
                 tail_length=0, show_tip_markers=False,
                 show_localizations=False, present_window=present_window,
                 title=_title('ribosome (60S, max-proj)'))

    # 3 — merge, no mask, no tracks
    rgb = vt.build_composite(er_disp, spt_disp, mask=None,
                             er_pct=(0, 100), er_gain=1.0,
                             spt_pct=(0, 100), spt_gain=1.0)
    ax = axes_row[2]
    ax.imshow(rgb, origin='upper')
    y0, y1, x0, x1 = crop if crop is not None else (0, er_disp.shape[0],
                                                    0, er_disp.shape[1])
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_xticks([]); ax.set_yticks([])
    vt._scale_bar(ax, x0, x1, y0, y1)
    if SHOW_TITLES:
        ax.set_title('merge (ER green / ribosome magenta)', fontsize=10)

    # 4 — mask + raw particles + uniform-magenta tracks (no classification)
    vt.mask_only_panel(
        axes_row[3], er_disp.shape, masks=[on_broad], mask_colors=[MASK_COLOR],
        crop=crop, spt_frame=spt_disp, spt_pct=(0, 100), spt_gain=1.0,
        spt_channel='magenta',
        traj_df=traj, frame_T=FRAME_T, color_by='fixed',
        fixed_color=TRACE_MAGENTA, tail_length=TAIL_LENGTH,
        show_tip_markers=True, show_localizations=show_rings,
        loc_color_by_zone=False, loc_marker_size=45,
        present_window=present_window,
        title=_title('ER mask + ribosome tracks'))

    # 5 — same, tracks coloured by zone classification, `ambiguous` omitted
    vt.mask_only_panel(
        axes_row[4], er_disp.shape, masks=[on_broad], mask_colors=[MASK_COLOR],
        crop=crop, spt_frame=spt_disp, spt_pct=(0, 100), spt_gain=1.0,
        spt_channel='magenta',
        traj_df=drop_zones(traj, zone_dict), frame_T=FRAME_T, color_by='zone',
        traj_zone_dict=zone_dict, tail_length=TAIL_LENGTH,
        show_tip_markers=True, show_localizations=show_rings,
        loc_color_by_zone=True, loc_marker_size=45,
        present_window=present_window,
        pool_on_near=True, pooled_color=ON_BROAD_COLOR,
        title=_title('zone classification (ambiguous omitted)'))

    if zoom_boxes:
        for j in range(5):
            for box, c in zip(zoom_boxes, ZOOM_BOX_COLORS):
                vt.draw_zoom_box(axes_row[j], box, color=c)


def add_zone_legend(ax):
    handles = [
        Line2D([], [], color=ON_BROAD_COLOR, lw=2.5, label='on ER (broad)'),
        Line2D([], [], color=OFF_ER_COLOR, lw=2.5, label='off ER'),
        Line2D([], [], color=vt.ZONE_COLORS['intermediate'], lw=2.5,
               label='intermediate'),
    ]
    ax.legend(handles=handles, loc='upper left', fontsize=8, framealpha=0.75,
              facecolor='black', edgecolor='none', labelcolor='white',
              title='ambiguous not shown', title_fontsize=7)
    ax.get_legend().get_title().set_color('#bbbbbb')


def build_one(stem, label):
    traj, er, on_broad, zone_dict, masks, spt_path = vt.load_fov(
        stem, paths, er_source='corrected', mask_name='on_ER_broad')

    spt_raw = vt.load_spt_frame(spt_path, FRAME_T, projection=SNAP_PROJ_N,
                                mode=SNAP_PROJ_MODE)

    # Normalize ONCE; every panel reuses these so the channels are identical
    # across the single-channel, merge and mask panels.
    er_disp = vt._norm(er, *ER_PCT) ** ER_GAMMA
    spt_disp = vt._norm(spt_raw, *SPT_PCT, gain=SPT_GAIN)

    crops = CROPS.get(stem) or pick_zoom_crops(stem, n=N_ZOOMS)
    n_rows = 1 + len(crops)

    fig, axes = plt.subplots(n_rows, 5, figsize=(25, 5.6 * n_rows))
    if n_rows == 1:
        axes = axes[None, :]

    panel_row(axes[0], traj, er_disp, spt_disp, on_broad, zone_dict,
              crop=None, show_rings=False, present_window=2,
              zoom_boxes=crops, row_label=f'{label} | full FOV')
    add_zone_legend(axes[0][4])

    for i, crop in enumerate(crops):
        panel_row(axes[1 + i], traj, er_disp, spt_disp, on_broad, zone_dict,
                  crop=crop, show_rings=True, present_window=1,
                  zoom_boxes=None, row_label=f'zoom {i + 1} {tuple(crop)}')

    fig.suptitle(
        f'{stem} | frame {FRAME_T} '
        f'({SNAP_PROJ_MODE}-proj {FRAME_T}-{FRAME_T + SNAP_PROJ_N}, '
        f'NOT instantaneous)', y=1.0, fontsize=11)
    fig.tight_layout()
    out_pdf = os.path.join(OUT, f'costain_{stem}_T{FRAME_T}.pdf')
    fig.savefig(out_pdf, bbox_inches='tight', dpi=200)
    plt.close(fig)
    print(f'wrote {out_pdf} ({n_rows} rows)')


def main():
    only = sys.argv[1:] or None
    for stem, label in FOVS:
        if only and stem not in only:
            continue
        build_one(stem, label)


if __name__ == '__main__':
    main()
