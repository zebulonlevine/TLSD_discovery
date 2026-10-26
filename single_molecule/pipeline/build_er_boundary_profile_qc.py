"""
build_er_boundary_profile_qc.py
=================================

ER mask QC (Fig S7H): boundary-crossing intensity profile.

The on_ER_broad mask is an adaptive LOCAL-brightness threshold on the
deconvolved ER image (run_filewise_er.py::_three_class_mask), so a raw
"in-mask brighter than off-mask" median comparison (build_er_brightness_qc.py)
is nearly circular. A boundary-crossing profile instead shows something a
trivial brightness cut would NOT guarantee: that the mask boundary sits on
a SHARP intensity transition and that the raw ER signal keeps FALLING AWAY
with distance outside the mask (the non-circular part -- a flat global-
level cut wouldn't produce a graded falloff at increasing distance from
the boundary; a mask tracing real ER structure does).

Per FOV (46 total, both degron backgrounds):
  1. raw ER = preprocessed_corrected/{stem}_er.npy
  2. mask = on_ER_broad = on_structure | near_ER (from {stem}_er_masks.npz)
  3. signed distance to the mask boundary via scipy.ndimage.distance_transform_edt:
       signed = edt(mask) - edt(~mask)   (+ inside, - outside), x pixel_size -> um
  4. z-score the raw ER image over ALL its pixels (per-FOV), so FOVs of very
     different absolute brightness (~150-820 a.u.) pool honestly
  5. mean z-scored intensity per signed-distance bin -> one profile per FOV

Plot: all per-FOV curves (faint) + across-FOV mean +/- SEM, split by degron,
dotted line at the boundary (distance 0). Per-FOV curves are the honest
independent unit (pixels within an FOV are spatially correlated), matching
build_er_brightness_qc.py's per-FOV-summary precedent.
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'er')))
import viz_tools as vt  # noqa: F401 (kept for ExpPaths/PIXEL_SIZE parity with sibling QC script)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
OUT_DIR = tlsd_paths.outp('PaperFigures_updated', 'figureS7', 'H')
SEG_DIR = os.path.join(tlsd_paths.spt('all_data', '2025_09_15_zeb_ER'), 'tracking_output_q=15p0', 'er_segmentation')
CORR_DIR = os.path.join(tlsd_paths.spt('all_data', '2025_09_15_zeb_ER'), 'preprocessed_corrected')

PIXEL_SIZE = 0.11  # um/px, matches viz_tools.PIXEL_SIZE

# Grey = noDegron, orange = Degron, using the same two shades as the Fig 5F
# zone posterior overlay's condition curves so the QC panel reads with the rest
# of the figure. (Was blue/red, which collided with the #3366cc/#cc3333 that
# encodes "CI excludes 0" in the region-effect panels.)
DEGRON_COLORS = {'noDegron': '#4d4d4d', 'Degron': '#e6550d'}

# FOV restriction: UNTREATED fields only. `HHT=000` is the no-drug condition;
# `HHT=005_t=000` is the drug series' own zero timepoint, i.e. imaged before
# run-off has had time to act, so it is legitimately untreated. Every later
# timepoint (t=0p5 ... t=10p5) is excluded.
#
# The trailing underscore in 't=000_' matters: without it the pattern also
# matches nothing useful but a bare 't=000' prefix test would wrongly admit
# 't=0p5' style stems if the naming ever changed to zero-padded halves.
#
# This is exactly the FOV set combinations.yaml pools as `noDegron` / `Degron`
# (both list HHT=000 + HHT=005_t=000), so it should reproduce the 10 / 14 FOV
# split the Fig 5F contrasts and the DiD use -- checked at runtime below.
UNTREATED_ONLY = True
EXPECTED_FOV = {'noDegron': 10, 'Degron': 14}


def is_untreated(stem):
    return ('HHT=000' in stem) or ('HHT=005_t=000_' in stem)

# signed-distance bin grid, in um. Bins are CENTERED on integer-pixel
# distances (0, +/-1px, +/-2px, ... x PIXEL_SIZE) with edges at half-pixel
# offsets. This matters: the Euclidean distance transform never yields
# |distance| < 1px (the nearest any pixel gets to the boundary is exactly
# 1px, either side), so pixels don't exist in the (-1px, +1px) sliver. If
# the edges instead sat ON 0 and integer-pixel multiples (as they did
# originally), np.digitize's left-closed [lo,hi) bins send the -1px pixels
# to the bin just below 0 but the +1px pixels to the bin at +1.5px, leaving
# an ASYMMETRIC one-sided empty bin at +0.5px. Centering on integer pixels
# instead puts +/-1px pixels in mirror-image bins and leaves a single
# SYMMETRIC empty bin at exactly 0 (honest: no pixel lies on the boundary).
DIST_MAX_UM = 2.2
_N_SIDE = int(round(DIST_MAX_UM / PIXEL_SIZE))
BIN_CENTERS_UM = np.arange(-_N_SIDE, _N_SIDE + 1) * PIXEL_SIZE
BIN_EDGES_UM = np.concatenate([BIN_CENTERS_UM - PIXEL_SIZE / 2,
                               [BIN_CENTERS_UM[-1] + PIXEL_SIZE / 2]])

# which mask boundaries to profile (label -> builder from the mask npz dict)
BOUNDARIES = {
    'on_ER_broad': lambda mm: mm['on_structure'].astype(bool) | mm['near_ER'].astype(bool),
    'on_structure': lambda mm: mm['on_structure'].astype(bool),
}


def discover_stems():
    stems = sorted(f[:-len('_er_masks.npz')] for f in os.listdir(SEG_DIR)
                   if f.endswith('_er_masks.npz'))
    if not UNTREATED_ONLY:
        return stems
    keep = [s for s in stems if is_untreated(s)]
    dropped = [s for s in stems if not is_untreated(s)]
    print(f'FOV filter: keeping {len(keep)} untreated, dropping '
          f'{len(dropped)} HHT-treated timepoints')
    for s in dropped:
        print(f'    dropped: {s}')
    return keep


def fov_profile(er_z, mask):
    """Mean z-scored intensity per signed-distance bin for one FOV/mask."""
    signed_px = distance_transform_edt(mask) - distance_transform_edt(~mask)
    signed_um = signed_px * PIXEL_SIZE
    idx = np.digitize(signed_um.ravel(), BIN_EDGES_UM) - 1
    flat_z = er_z.ravel()
    prof = np.full(len(BIN_CENTERS_UM), np.nan)
    for b in range(len(BIN_CENTERS_UM)):
        sel = idx == b
        if sel.any():
            prof[b] = flat_z[sel].mean()
    return prof


def main():
    stems = discover_stems()
    print(f'{len(stems)} FOVs with ER segmentation')

    profiles = {}   # (label, stem) -> profile array
    degron_of = {}
    rows = []
    for stem in stems:
        er_path = os.path.join(CORR_DIR, f'{stem}_er.npy')
        mask_path = os.path.join(SEG_DIR, f'{stem}_er_masks.npz')
        if not (os.path.exists(er_path) and os.path.exists(mask_path)):
            print(f'  SKIP {stem}: missing ER image or mask file')
            continue
        er = np.load(er_path).astype(np.float64)
        er_z = (er - er.mean()) / er.std()
        degron = 'Degron' if 'degron=001' in stem else 'noDegron'
        degron_of[stem] = degron
        with np.load(mask_path) as mm:
            for label, build in BOUNDARIES.items():
                mask = build(mm)
                prof = fov_profile(er_z, mask)
                profiles[(label, stem)] = prof
                for c, z in zip(BIN_CENTERS_UM, prof):
                    rows.append({'fov': stem, 'degron': degron, 'boundary': label,
                                 'distance_um': float(c), 'mean_z': float(z)})

    counts = {d: sum(1 for v in degron_of.values() if v == d)
              for d in ('noDegron', 'Degron')}
    print(f'FOVs used by condition: {counts}')
    if UNTREATED_ONLY and counts != EXPECTED_FOV:
        # Not fatal -- but it means this panel's FOV set has drifted from the
        # one combinations.yaml pools for the Fig 5F contrasts, so the QC would
        # no longer be describing the masks those contrasts actually used.
        print(f'  WARNING: expected {EXPECTED_FOV} (the combinations.yaml '
              f'noDegron/Degron pools). This panel no longer matches the FOV set '
              f'behind the Fig 5F contrasts.')

    df = pd.DataFrame(rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUT_DIR, 'er_boundary_profile_qc_per_fov.csv')
    df.to_csv(csv_path, index=False)
    n_fov = df['fov'].nunique()
    print(f'wrote {csv_path} ({n_fov} FOVs x {len(BOUNDARIES)} boundaries)')

    labels = list(BOUNDARIES)
    fig, axes = plt.subplots(1, len(labels), figsize=(7.5 * len(labels), 5.5),
                             sharey=True, squeeze=False)
    axes = axes[0]

    for ax, label in zip(axes, labels):
        stems_here = [s for (lab, s) in profiles if lab == label]
        for stem in stems_here:
            prof = profiles[(label, stem)]
            ax.plot(BIN_CENTERS_UM, prof, color=DEGRON_COLORS[degron_of[stem]],
                    lw=0.6, alpha=0.25, zorder=1)
        for degron, color in DEGRON_COLORS.items():
            mat = np.vstack([profiles[(label, s)] for s in stems_here
                             if degron_of[s] == degron])
            mean = np.nanmean(mat, axis=0)
            sem = np.nanstd(mat, axis=0, ddof=1) / np.sqrt(mat.shape[0])
            ax.plot(BIN_CENTERS_UM, mean, color=color, lw=2.4, zorder=3,
                    label=f'{degron} (n={mat.shape[0]})')
            ax.fill_between(BIN_CENTERS_UM, mean - sem, mean + sem,
                            color=color, alpha=0.25, lw=0, zorder=2)
        ax.axvline(0, color='k', ls=':', lw=1.2, zorder=4)
        ax.set_xlabel(f'signed distance to {label} boundary (µm)\n(− outside / + inside)')
        ax.set_title(f'{label} boundary')
        ax.legend(loc='upper left', fontsize=8)
    axes[0].set_ylabel('raw ER intensity (per-FOV z-score)')
    scope = ('untreated FOVs only (HHT=000 + HHT=005 t=0)' if UNTREATED_ONLY
             else 'all FOVs incl. HHT timepoints')
    fig.suptitle(f'ER mask boundary-crossing intensity profile '
                 f'(n={n_fov} FOVs; {scope})\n'
                 f'faint = per-FOV; bold = across-FOV mean ± SEM by degron',
                 y=1.02, fontsize=12)
    fig.tight_layout()
    out_pdf = os.path.join(OUT_DIR, 'er_boundary_profile_qc.pdf')
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f'wrote {out_pdf}')

    # --- quantify the non-circular signature, for the on_ER_broad boundary ---
    label = 'on_ER_broad'
    stems_here = [s for (lab, s) in profiles if lab == label]
    inside_bin = np.argmin(np.abs(BIN_CENTERS_UM - 1.0))     # ~ +1.0 um (plateau inside)
    far_out_bin = np.argmin(np.abs(BIN_CENTERS_UM - (-2.0)))  # ~ -2.0 um (far outside)
    near_out_bin = np.argmin(np.abs(BIN_CENTERS_UM - (-0.3)))  # just outside the edge
    n_monotonic = 0
    gaps = []
    falloffs = []
    for stem in stems_here:
        p = profiles[(label, stem)]
        gaps.append(p[inside_bin] - p[far_out_bin])
        falloffs.append(p[near_out_bin] - p[far_out_bin])  # >0 => signal falls away outward
        # monotonic rise across the boundary: inside-edge > just-outside > far-outside
        if p[inside_bin] > p[near_out_bin] > p[far_out_bin]:
            n_monotonic += 1
    gaps = np.array(gaps); falloffs = np.array(falloffs)
    print(f'\non_ER_broad boundary, across {len(stems_here)} FOVs:')
    print(f'  inside(+1µm) − far-outside(−2µm) z-gap: mean={gaps.mean():.2f} '
          f'median={np.median(gaps):.2f}  (all positive: {(gaps > 0).mean()*100:.0f}%)')
    print(f'  outside falloff [near-edge(−0.3µm) − far(−2µm)]: mean={falloffs.mean():.2f}  '
          f'(positive i.e. falls away outward: {(falloffs > 0).mean()*100:.0f}% of FOVs)')
    print(f'  monotonic rise inside>near-out>far-out: {n_monotonic}/{len(stems_here)} FOVs')


if __name__ == '__main__':
    main()
