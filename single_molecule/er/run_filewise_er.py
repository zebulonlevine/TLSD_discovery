#!/usr/bin/env python3
"""
run_filewise_er.py

Per-file processing for ER+SPT dual-channel experiments.

Steps per file:
  1. Read frame_interval and camera_bg from the source ND2 (NOT from the
     frame-stripped SPT TIF, which lacks ND2 metadata). Override via
     settings['override_frame_interval'].
  2. Track (quot) — spot detection + linking → _traj.csv
  3. ER segmentation — three-class adaptive thresholding → masks + processed image
  4. Trajectory classification — binomial enrichment test
     (on_structure / near_ER / off_ER)
  5. Export subset traj CSVs per --subsets (default: on_structure, on_ER_broad, off_ER)
  6. saSPT on ALL trajectories → posterior, MLE, rolling windows

Subset (zone) saSPT is no longer run per FOV. It is now performed by
run_combinationwise.py against per-condition pooled zone trajectories,
declared in combinations.yaml with a 'zone:' field (e.g. zone:
on_structure). Per-FOV n is too low per zone to fit a reliable saSPT
population prior; pooling per condition gives a far cleaner posterior.

ER segmentation uses adaptive local thresholding at two spatial scales:
  - Fine block (~25 px / ~2.75 µm at 0.11 µm/px): ER tubules/sheets
  - Broad block (~101 px / ~11.1 µm): ER-dense region
  - on_structure = fine AND broad (prevents spurious detections in dark regions)
  - Optional pct_above_strict: stricter threshold for isolated tubules
  - on_ER_broad = on_structure + near_ER combined

Trajectory classification uses binomial enrichment test:
  - H0: random placement gives prob = area_on of landing on structure
  - H1: trajectory is enriched on structure (one-sided binomtest, alpha=0.05)
  - Prevents 97% on-ER artifact when ER covers large image area
  - Auto-detects whether x/y coordinates are in pixels or microns

Usage (called by SLURM, not typically run directly):
    python run_filewise_er.py \\
        -i preprocessed/sample_spt.tif \\
        -s settings.pkl \\
        -c condition_name \\
        --er-image preprocessed_corrected/sample_er.npy
"""

from quot.core import track_file
import argparse
import pickle
import os
import re
import sys
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


def _read_spt_bg(spt_tif_path, n_frames=10):
    """
    Estimate camera background from the SPT TIF.

    The SPT TIF's pixel values are bit-identical to the corresponding
    ND2 frames (preprocess_er_spt_file.py extracts frames 100+ without
    rescaling), so reading bg from either is mathematically equivalent.
    We read from the TIF because pims-nd2reader returns nan-filled
    arrays for "gap frames" in multi-channel acquisitions, and that
    silently propagates to a nan camera_bg and blank fit columns. The
    TIF is a clean sequential dump with no such gaps.

    Reads the first n_frames of the SPT TIF and returns the median
    pixel value. Uses a small frame count for speed; camera offset is
    very stable across frames within a file.
    """
    import tifffile
    try:
        with tifffile.TiffFile(spt_tif_path) as tf:
            n_pages = len(tf.pages)
            sample_n = min(n_frames, n_pages)
            sample = np.stack([tf.pages[i].asarray() for i in range(sample_n)])
    except Exception as e:
        raise RuntimeError(
            f"Failed to read SPT TIF for camera_bg estimation "
            f"({spt_tif_path}): {e}.") from e
    bg = float(np.median(sample))
    if not np.isfinite(bg) or bg <= 0:
        raise RuntimeError(
            f"SPT TIF at {spt_tif_path} returned a non-finite or "
            f"non-positive median (bg={bg}). Refusing to proceed.")
    return bg


# ============================================================
# ER segmentation (inlined to avoid import dependency on er_analysis_lib
# in the SLURM environment)
# ============================================================

def _process_er(er_image, flatfield_sigma=100, psf_sigma=0.74,
                rl_iter=2, no_flatfield=False, no_deconv=False):
    """Residual flat-field + RL deconvolution."""
    processed = er_image.astype(np.float64).copy()

    if not no_flatfield:
        bg = gaussian_filter(processed, sigma=flatfield_sigma)
        bg = np.maximum(bg, 1.0)
        processed = (processed / bg) * processed.mean()

    if not no_deconv:
        try:
            from skimage.restoration import richardson_lucy
            psf_size = int(psf_sigma * 12) | 1
            y, x = np.mgrid[-psf_size//2:psf_size//2+1,
                             -psf_size//2:psf_size//2+1]
            psf = np.exp(-(x**2 + y**2) / (2 * psf_sigma**2))
            psf /= psf.sum()
            img_pos = np.maximum(processed - processed.min() + 1, 1)
            dec = richardson_lucy(img_pos, psf, num_iter=rl_iter, clip=False)
            processed = dec * (processed.mean() / dec.mean())
        except ImportError:
            print("  WARNING: skimage not available, skipping deconvolution")

    return processed


def _three_class_mask(image, pct_above=1.5, block_fine=25, block_broad=101,
                      min_size=500, pct_above_strict=None, pct_floor=50,
                      require_broad=True):
    """Three-class segmentation with optional strict threshold for isolated tubules."""
    try:
        from skimage.filters import threshold_local
        from skimage.morphology import remove_small_objects, closing, disk
    except ImportError:
        print("  WARNING: skimage not available, skipping ER segmentation")
        return None

    def _mask(img, pct, block, minsz):
        local_mean = threshold_local(img, block_size=block,
                                      method='gaussian', offset=0)
        m = img > local_mean * (1 + pct / 100.0)
        m = closing(m, disk(2))
        m = remove_small_objects(m, min_size=minsz)
        return m

    mask_fine = _mask(image, pct_above, block_fine, min_size)
    mask_broad = _mask(image, pct_above, block_broad, min_size)

    # Global percentile floor (per-FOV, on the processed image). A pixel
    # must exceed the pct_floor-th percentile of the whole image to count
    # as ER, on top of passing local thresholding. Suppresses
    # locally-bright-but-globally-dim detections in dark regions.
    # pct_floor=None disables (recovers pre-floor behavior).
    if pct_floor is not None:
        above_floor = image > np.percentile(image, pct_floor)
        mask_fine = mask_fine & above_floor
        mask_broad = mask_broad & above_floor

    # require_broad=True (default): on_structure must pass BOTH fine and
    # broad local thresholds ("on has to be near"). require_broad=False:
    # on_structure is just fine (+ floor), with no broad-density gate —
    # the percentile floor alone guards against dark-region false positives.
    if require_broad:
        on_in_broad = mask_fine & mask_broad
    else:
        on_in_broad = mask_fine.copy()

    if pct_above_strict is not None:
        mask_fine_strict = _mask(image, pct_above_strict, block_fine, min_size)
        if pct_floor is not None:
            mask_fine_strict = mask_fine_strict & above_floor
        on_outside_broad = mask_fine_strict & ~mask_broad
        on_structure = on_in_broad | on_outside_broad
    else:
        on_structure = on_in_broad

    return {
        'on_structure': on_structure,
        'near_ER': mask_broad & ~on_structure,
        'off_ER': ~mask_broad & ~on_structure,
    }


def _detect_coord_scale(traj_df, image_shape, pixel_size_um=0.11):
    """Auto-detect if coordinates are pixels (return 1.0) or microns (return pixel_size_um)."""
    max_x = traj_df['x'].max()
    max_y = traj_df['y'].max()
    img_h, img_w = image_shape[:2]
    px_match = max(max_x / img_w, max_y / img_h)
    um_match = max(max_x / (img_w * pixel_size_um), max_y / (img_h * pixel_size_um))
    if abs(px_match - 1.0) < abs(um_match - 1.0):
        print(f"  Coordinates in PIXELS (max_x={max_x:.0f}, img_w={img_w}) → scale=1.0")
        return 1.0
    else:
        print(f"  Coordinates in MICRONS (max_x={max_x:.1f}) → scale={pixel_size_um}")
        return pixel_size_um


def classify_trajectories(traj_df, er_processed, pixel_size_um,
                           pct_above=1.5, block_fine=25, block_broad=101,
                           min_size=500, fraction_threshold=0.95,
                           enrichment_alpha=0.05, min_frames=2,
                           pct_above_strict=None, pct_floor=50,
                           require_broad=True):
    """
    Classify trajectories using binomial enrichment test.

    A trajectory is on_structure only if its fraction on-structure is
    significantly higher than expected from random placement (given the
    ER area fraction). This prevents the 97% on-ER artifact when ER
    covers a large portion of the image.

    Auto-detects whether coordinates are in pixels or microns.

    Three classes:
      on_structure: enriched on fine ER structure
      near_ER: enriched on broad (on+near) but not fine
      off_ER: not enriched on either

    Returns
    -------
    traj_df : DataFrame with 'er_zone' and 'er_zone_label' columns added
    per_traj : DataFrame with trajectory, n_frames, frac_on, frac_near,
               enrichment_p, enrichment_p_near, er_zone_label
    masks : dict of boolean arrays
    """
    from scipy.stats import binomtest

    masks = _three_class_mask(er_processed, pct_above, block_fine,
                               block_broad, min_size, pct_above_strict,
                               pct_floor, require_broad)
    if masks is None:
        return traj_df, pd.DataFrame(), None

    df = traj_df.copy()

    # Auto-detect coordinate scale
    coord_scale = _detect_coord_scale(df, er_processed.shape, pixel_size_um)

    py = np.clip((df['y'].values / coord_scale).astype(int),
                  0, er_processed.shape[0] - 1)
    px = np.clip((df['x'].values / coord_scale).astype(int),
                  0, er_processed.shape[1] - 1)

    # Per-localization zone (three classes)
    zone = np.full(len(df), 'off_ER', dtype=object)
    zone[masks['near_ER'][py, px]] = 'near_ER'
    zone[masks['on_structure'][py, px]] = 'on_structure'
    df['er_zone'] = zone

    # Area fractions (null probability for binomial test)
    area_on = masks['on_structure'].mean()
    area_on_near = area_on + masks['near_ER'].mean()

    # --- Derived minimum trajectory lengths for testability ---------------
    # For a one-sided binomial at level `enrichment_alpha`, the strongest
    # possible evidence is k == n (every localization on the mask). Its
    # p-value is area**n. So a trajectory can only POSSIBLY reach
    # significance if area**n < alpha, i.e. n >= log(alpha)/log(area).
    # Below that length the test has no power and a trajectory cannot be
    # classified on that mask no matter where it sits — it is 'ambiguous',
    # not 'off_ER'. This makes the length cutoff a derived property of each
    # FOV's coverage rather than a guessed constant.
    def _n_min_binom(area, alpha):
        if not (0 < area < 1):
            return None                      # degenerate mask: untestable
        n = 1
        while area ** n >= alpha:
            n += 1
            if n > 10000:
                return None
        return n

    n_min_on   = _n_min_binom(area_on,      enrichment_alpha)
    n_min_near = _n_min_binom(area_on_near, enrichment_alpha)
    # `min_frames` defaults to 2 (the structural floor: saSPT needs >=2
    # localizations for >=1 displacement, and enrichment is trivial on 1
    # frame). By default it does nothing — the derived n_min rules the
    # classification cutoff. It is retained as a *sweepable robustness
    # lever*: raise it (5, 8, 12) to test whether excluding more short
    # trajectories changes the downstream diffusion conclusion. Effective
    # cutoff is the larger of the structural floor and the derived n_min.
    # None => that test is impossible (infinite cutoff).
    cutoff_on   = max(min_frames, n_min_on)   if n_min_on   is not None else np.inf
    cutoff_near = max(min_frames, n_min_near) if n_min_near is not None else np.inf

    # Per-trajectory stats
    per_traj = df.groupby('trajectory').agg(
        n_frames=('frame', 'count'),
        n_on=('er_zone', lambda s: (s == 'on_structure').sum()),
        n_near=('er_zone', lambda s: (s == 'near_ER').sum()),
        frac_on=('er_zone', lambda s: (s == 'on_structure').mean()),
        frac_near=('er_zone', lambda s: (s == 'near_ER').mean()),
    ).reset_index()

    # Binomial enrichment test for on_structure (FINE mask).
    enrichment_p = np.ones(len(per_traj))
    for i, row in per_traj.iterrows():
        n = int(row['n_frames'])
        k = int(row['n_on'])
        if n >= cutoff_on and 0 < area_on < 1:
            enrichment_p[i] = binomtest(k, n, area_on, alternative='greater').pvalue
    per_traj['enrichment_p'] = enrichment_p

    # Binomial enrichment test for near_ER (combined on+near, BROAD mask).
    enrichment_p_near = np.ones(len(per_traj))
    for i, row in per_traj.iterrows():
        n = int(row['n_frames'])
        k = int(row['n_on']) + int(row['n_near'])
        if n >= cutoff_near and 0 < area_on_near < 1:
            enrichment_p_near[i] = binomtest(k, n, area_on_near, alternative='greater').pvalue
    per_traj['enrichment_p_near'] = enrichment_p_near

    # Testability flags
    on_testable   = per_traj['n_frames'] >= cutoff_on
    near_testable = per_traj['n_frames'] >= cutoff_near

    # --- Classification ----------------------------------------------------
    # Cascade: default 'intermediate' (testable but spatially mixed) ->
    # ambiguous (too short to test) -> off_ER (>=95% off the broad mask) ->
    # near_ER -> on_structure.
    #
    # on_structure now uses the BROAD spatial-purity gate (frac_on+frac_near
    # >= threshold: the particle essentially never leaves the ER region) but
    # the FINE binomial enrichment (specifically concentrated on fine
    # structure beyond chance). This loosens the spatial requirement while
    # keeping the specificity test strict.
    # Default 'intermediate': testable but spatially mixed (broad fraction
    # in the 5-95% band -- neither >=95% on/near nor >=95% off). off_ER is a
    # POSITIVE label below, not the catch-all default.
    per_traj['er_zone_label'] = 'intermediate'

    # Ambiguous: too short for even the on-structure test (smallest n_min).
    ambiguous_mask = ~on_testable
    per_traj.loc[ambiguous_mask, 'er_zone_label'] = 'ambiguous'

    # off_ER (POSITIVE): >= fraction_threshold of localizations OUTSIDE the
    # broad mask (broad fraction <= 1 - fraction_threshold), and n >=
    # cutoff_on -- same floor as on_structure, so on-vs-off is length-matched.
    # No binomial: 'enriched off' is near-trivial given the large off area and
    # would force cutoff_off ~9, biasing off_ER toward long (slow) tracks.
    off_mask = on_testable & (
        (per_traj['frac_on'] + per_traj['frac_near']) <= (1.0 - fraction_threshold))
    per_traj.loc[off_mask, 'er_zone_label'] = 'off_ER'

    near_mask = (
        ((per_traj['frac_on'] + per_traj['frac_near']) >= fraction_threshold) &
        (per_traj['enrichment_p_near'] < enrichment_alpha) &
        near_testable
    )
    per_traj.loc[near_mask, 'er_zone_label'] = 'near_ER'

    on_mask = (
        ((per_traj['frac_on'] + per_traj['frac_near']) >= fraction_threshold) &
        (per_traj['enrichment_p'] < enrichment_alpha) &
        on_testable
    )
    per_traj.loc[on_mask, 'er_zone_label'] = 'on_structure'

    # Merge labels back onto df for export
    df = df.merge(per_traj[['trajectory', 'er_zone_label']],
                  on='trajectory', how='left')

    # Print classification summary
    n_total = len(per_traj)
    n_on = (per_traj['er_zone_label'] == 'on_structure').sum()
    n_near = (per_traj['er_zone_label'] == 'near_ER').sum()
    n_off = (per_traj['er_zone_label'] == 'off_ER').sum()
    n_amb = (per_traj['er_zone_label'] == 'ambiguous').sum()
    n_int = (per_traj['er_zone_label'] == 'intermediate').sum()
    # how many trajectories are too short for the (harder) near test but
    # long enough for the on test — these can be on/off but never near
    n_near_untestable = int((per_traj['n_frames'] < cutoff_near).sum())
    print(f"    Enrichment test (alpha={enrichment_alpha}, min_frames={min_frames}):")
    print(f"      Area on_structure: {area_on:.1%}   on_ER_broad: {area_on_near:.1%}")
    print(f"      n_min  on_structure: {n_min_on}  (cutoff {cutoff_on})   "
          f"near/broad: {n_min_near}  (cutoff {cutoff_near})")
    print(f"      on_structure: {n_on}/{n_total} ({n_on/n_total:.1%})")
    print(f"      near_ER:      {n_near}/{n_total} ({n_near/n_total:.1%})")
    print(f"      off_ER:       {n_off}/{n_total} ({n_off/n_total:.1%})")
    print(f"      ambiguous:    {n_amb}/{n_total} ({n_amb/n_total:.1%})  "
          f"[too short for on-structure test]")
    print(f"      intermediate: {n_int}/{n_total} ({n_int/n_total:.1%})  "
          f"[testable but mixed 5-95% broad fraction; excluded from pools]")
    print(f"      near-untestable: {n_near_untestable}/{n_total} "
          f"({n_near_untestable/n_total:.1%})  [too short for broad test; "
          f"high area_on_near => broad analysis power-limited]")

    return df, per_traj, masks


def export_subset_traj(traj_df, subset_name, output_path):
    """
    Export trajectories for a subset, re-indexing trajectory IDs.

    Supports combined subsets:
      'on_structure' -> on_structure only
      'near_ER'      -> near_ER only
      'on_ER_broad'  -> on_structure + near_ER combined
      'off_ER'       -> off_ER only
    """
    SUBSET_LABELS = {
        'on_structure': ['on_structure'],
        'near_ER':      ['near_ER'],
        'on_ER_broad':  ['on_structure', 'near_ER'],
        'off_ER':       ['off_ER'],
        'ambiguous':    ['ambiguous'],
    }
    labels = SUBSET_LABELS.get(subset_name, [subset_name])
    subset = traj_df[traj_df['er_zone_label'].isin(labels)].copy()
    if subset.empty:
        return 0

    # Drop ER columns, re-index trajectories from 0
    drop_cols = [c for c in subset.columns if c.startswith(('er_zone', 'er_intensity'))]
    orig_cols = [c for c in subset.columns if c not in drop_cols]
    subset = subset[orig_cols].copy()

    tmap = {old: new for new, old in enumerate(subset['trajectory'].unique())}
    subset['trajectory'] = subset['trajectory'].map(tmap)

    subset.to_csv(output_path, index=False)
    return len(tmap)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Per-file ER+SPT processing: track → classify → saSPT×3")
    parser.add_argument('-i', '--inputfile', required=True,
                        help='Input SPT image (TIF or ND2)')
    parser.add_argument('-s', '--settingsPKL', required=True,
                        help='Settings pickle path')
    parser.add_argument('-c', '--condition', required=False,
                        help='Condition name')
    parser.add_argument('--er-image', required=True,
                        help='Path to illumination-corrected ER .npy image')

    # ER segmentation parameters
    parser.add_argument('--pct-above', type=float, default=1.5)
    parser.add_argument('--pct-above-strict', type=float, default=None,
                        help='Stricter pct threshold for isolated tubules outside '
                             'broad ER region (e.g., 5.0). None = exclude entirely.')
    parser.add_argument('--block-fine', type=int, default=25)
    parser.add_argument('--block-broad', type=int, default=101)
    parser.add_argument('--min-object-size', type=int, default=500)
    parser.add_argument('--require-broad',
                        type=lambda s: str(s).strip().lower()
                            in ('true','1','yes','t','y'),
                        default=True,
                        help='true|false. false => on_structure is fine+'
                             'floor only, no broad-density gate.')
    parser.add_argument('--pct-floor', type=float, default=50,
                        help='Global percentile floor on the processed ER '
                             'image; pixels below this percentile cannot be '
                             'ER. Set to a negative value to disable.')
    parser.add_argument('--frac-threshold', type=float, default=0.95)
    parser.add_argument('--enrichment-alpha', type=float, default=0.05,
                        help='Significance level for binomial enrichment test')
    parser.add_argument('--min-frames', type=int, default=2,
                        help='Structural floor on trajectory length (default 2 '
                             '= saSPT minimum). Per-FOV derived n_min rules the '
                             'classification cutoff above this; raise --min-frames '
                             'only to sweep short-trajectory sensitivity.')
    parser.add_argument('--subsets', type=str, default='on_structure,on_ER_broad,off_ER',
                        help='Comma-separated subsets to run saSPT on. Options: '
                             'on_structure, near_ER, on_ER_broad (on+near), off_ER')
    parser.add_argument('--no-flatfield', action='store_true')
    parser.add_argument('--no-deconv', action='store_true')

    args = parser.parse_args()

    print("")
    print(f"Input SPT:  {args.inputfile}")
    print(f"ER image:   {args.er_image}")
    print(f"Settings:   {args.settingsPKL}")
    if args.condition:
        print(f"Condition:  {args.condition}")

    with open(args.settingsPKL, 'rb') as f:
        settings = pickle.load(f)

    sys.path.append(settings['io']['code_directory'])
    import trackingCodeFunctions

    input_basename = os.path.basename(args.inputfile)
    file_stem = re.sub(r'(_spt)?\.tif$|\.nd2$', '', input_basename)

    # Parse subset list
    subsets = [s.strip() for s in args.subsets.split(',')]

    # Output paths
    traj_all_path = os.path.join(
        settings['io']['traj_directory'],
        f'{file_stem}_traj.csv')

    # Create subset directories
    for subset_name in subsets:
        d = os.path.join(settings['io']['analysis_directory'],
                         f'trajectories_{subset_name}')
        os.makedirs(d, exist_ok=True)

    # Segmented ER output directory
    seg_dir = os.path.join(settings['io']['analysis_directory'],
                           'er_segmentation')
    os.makedirs(seg_dir, exist_ok=True)

    # ================================================================
    # Step 1: Update metadata and track
    # ================================================================
    # NOTE: do NOT call update_settings_with_image_metadata here. That
    # function probes settings['io']['data_directory'] for ND2 files; in
    # the ER pipeline data_directory is preprocessed/ (which contains
    # _spt.tif files with frames 1..99 removed and stripped of camera
    # metadata). It would silently fall into the TIFF branch and leave
    # frame_interval at whatever default was in settings (typically
    # 0.01s, while the real interval is ~0.015091s) — a silent ~30%
    # error in every D estimate.
    #
    # Instead, locate the source ND2 (under <exp>/raw_data/<basename>.nd2
    # per batch_er_workflow.py convention) and read frame_rate directly.
    # If settings has override_frame_interval=True (set by
    # batch_er_workflow when settings_override.yaml provides an
    # explicit frame_interval), respect that and skip the probe.
    if settings.get('override_frame_interval'):
        print(f"  frame_interval locked by override: "
              f"{settings['saspt']['frame_interval']} s")
        # Even with override, we still need camera_bg for quot's fitter.
        # Read it from the SPT TIF rather than relying on a default.
        bg_level = _read_spt_bg(args.inputfile)
        settings['quot']['localize']['camera_bg'] = bg_level
        print(f"  camera_bg from SPT TIF: {bg_level:.1f}")
    else:
        # We need two pieces of metadata:
        #   1. frame_interval  — from the ND2 (TIF has it stripped)
        #   2. camera_bg       — from the SPT TIF, NOT the ND2
        #
        # The SPT TIF is what quot will actually fit (it's a clean frame
        # 100+ extraction from the ND2's SPT channel). Its pixel values
        # are bit-identical to the corresponding ND2 frames.
        #
        # Why not just read camera_bg from the ND2 too? pims-nd2reader
        # represents "gap frames" (frames the acquisition skipped, common
        # in multi-channel interleaved acquisitions) as nan-filled arrays.
        # If the gap happens to be near frame 0, np.median(img[0])
        # returns nan; quot then gets camera_bg=nan and every Gaussian
        # fit produces nan output, with no error. The TIF, being a
        # sequential dump of real frames, has no such gaps.
        spt_dir = os.path.dirname(os.path.abspath(args.inputfile))
        exp_dir = os.path.dirname(spt_dir)
        nd2_basename = re.sub(r'_spt\.tif$|\.tif$', '', input_basename)
        nd2_path = os.path.join(exp_dir, 'raw_data', f'{nd2_basename}.nd2')
        if not os.path.exists(nd2_path):
            raise FileNotFoundError(
                f"Cannot find source ND2 at {nd2_path} (derived from "
                f"input {args.inputfile}). The ER pipeline requires the "
                f"original ND2 to read frame_interval, since the SPT TIF "
                f"lacks that metadata. Either:\n"
                f"  (a) Move/link the source ND2 to {nd2_path}, OR\n"
                f"  (b) Set override_frame_interval=True with an explicit "
                f"frame_interval in settings (or settings_override.yaml).")
        try:
            import pims
            img = pims.open(nd2_path)
            dt = 1.0 / img.frame_rate
        except Exception as e:
            raise RuntimeError(
                f"Failed to read frame_rate from {nd2_path}: {e}. "
                f"Either fix the ND2 or set override_frame_interval=True "
                f"with an explicit frame_interval in settings.") from e
        if not np.isfinite(dt) or dt <= 0:
            raise RuntimeError(
                f"ND2 at {nd2_path} returned a non-finite or non-positive "
                f"frame_interval (dt={dt}). Refusing to proceed.")
        bg_level = _read_spt_bg(args.inputfile)
        settings['quot']['localize']['camera_bg'] = bg_level
        settings['quot']['track']['frame_interval'] = dt
        settings['saspt']['frame_interval'] = dt
        print(f"  frame_interval from ND2 {nd2_path}: {dt:.6f}s")
        print(f"  camera_bg from SPT TIF: {bg_level:.1f}")

    pixel_size = settings['quot']['track']['pixel_size_um']

    print('\n=== Step 1: Tracking ===')
    track_file(args.inputfile, out_csv=traj_all_path, **settings['quot'])
    print(f'  Saved: {traj_all_path}')

    # ================================================================
    # Sanity check: did quot actually fit anything?
    # ================================================================
    # If camera_bg or other fit parameters are wildly off, quot can
    # produce a CSV with detection columns populated (frame, x_detect,
    # y_detect, loc_idx, trajectory) but FIT columns all blank
    # (x, y, I0, bg, *_err). saSPT and zone classification both need
    # the fit columns. We refuse to proceed and write zone subsets from
    # such data.
    _traj_df_check = pd.read_csv(traj_all_path, nrows=5000)
    n_total = len(_traj_df_check)
    if n_total == 0:
        raise RuntimeError(
            f"quot wrote an empty trajectory CSV at {traj_all_path}. "
            f"Tracking failed entirely.")
    n_with_fit = _traj_df_check['x'].notna().sum() if 'x' in _traj_df_check.columns else 0
    frac_with_fit = n_with_fit / n_total
    if frac_with_fit < 0.1:
        # Print enough context for debugging without dumping all of settings
        loc = settings['quot']['localize']
        raise RuntimeError(
            f"Quot's fit columns are blank for "
            f"{100*(1-frac_with_fit):.1f}% of detections "
            f"({n_total-n_with_fit}/{n_total} sampled). This usually means "
            f"camera_bg is wrong for the data, the localize method "
            f"choked, or the SPT TIF is degenerate. Aborting before "
            f"writing zone subsets from broken data.\n"
            f"  traj path:  {traj_all_path}\n"
            f"  camera_bg:  {loc.get('camera_bg')}\n"
            f"  method:     {loc.get('method')}\n"
            f"  frame_intv: {settings['quot']['track']['frame_interval']}\n"
            f"Inspect a few rows manually:  head -5 {traj_all_path}")
    del _traj_df_check

    # ================================================================
    # Step 2: ER segmentation
    # ================================================================
    print('\n=== Step 2: ER Segmentation ===')

    er_raw = np.load(args.er_image).astype(np.float64)
    er_processed = _process_er(er_raw,
                                no_flatfield=args.no_flatfield,
                                no_deconv=args.no_deconv)

    # Save processed ER and masks
    np.save(os.path.join(seg_dir, f'{file_stem}_er_processed.npy'),
            er_processed.astype(np.float32))

    masks = _three_class_mask(er_processed, args.pct_above,
                               args.block_fine, args.block_broad,
                               args.min_object_size,
                               args.pct_above_strict,
                               (None if args.pct_floor < 0 else args.pct_floor),
                               args.require_broad)

    if masks is not None:
        np.savez_compressed(
            os.path.join(seg_dir, f'{file_stem}_er_masks.npz'),
            on_structure=masks['on_structure'],
            near_ER=masks['near_ER'],
            off_ER=masks['off_ER'])

        print(f'  on_structure: {masks["on_structure"].mean():.1%} of pixels')
        print(f'  near_ER:      {masks["near_ER"].mean():.1%} of pixels')
        print(f'  off_ER:       {masks["off_ER"].mean():.1%} of pixels')
    else:
        print('  WARNING: Segmentation failed, running saSPT on all only')

    # ================================================================
    # Step 3: Classify trajectories
    # ================================================================
    print('\n=== Step 3: Classify Trajectories ===')

    traj_df = pd.read_csv(traj_all_path)

    if masks is not None:
        traj_df, per_traj, _ = classify_trajectories(
            traj_df, er_processed, pixel_size,
            pct_above=args.pct_above,
            block_fine=args.block_fine,
            block_broad=args.block_broad,
            min_size=args.min_object_size,
            pct_floor=(None if args.pct_floor < 0 else args.pct_floor),
            require_broad=args.require_broad,
            fraction_threshold=args.frac_threshold,
            enrichment_alpha=args.enrichment_alpha,
            min_frames=args.min_frames,
            pct_above_strict=args.pct_above_strict)

        # Save classification
        per_traj['source_file'] = input_basename
        per_traj.to_csv(
            os.path.join(seg_dir, f'{file_stem}_traj_classification.csv'),
            index=False)

        # Export subset traj CSVs
        for subset_name in subsets:
            subset_path = os.path.join(
                settings['io']['analysis_directory'],
                f'trajectories_{subset_name}',
                f'{file_stem}_traj.csv')
            n = export_subset_traj(traj_df, subset_name, subset_path)
            print(f'  Exported {subset_name}: {n} trajectories → {subset_path}')
    else:
        per_traj = None

    # ================================================================
    # Step 4: saSPT on ALL trajectories
    # ================================================================
    print('\n=== Step 4: saSPT (all) ===')
    trackingCodeFunctions.executeSASPT(traj_all_path, settings)

    print('\n  Rolling windows (all)')
    trackingCodeFunctions.analyze_rolling_windows(traj_all_path, settings)

    # ================================================================
    # Step 5: (removed) — subset saSPT now happens at combinationwise time
    # ================================================================
    # The previous per-FOV subset saSPT loop has been removed. Per-FOV
    # n is too low per zone to fit a reliable population prior in saSPT,
    # so the resulting posteriors were noisy and not directly comparable
    # across FOVs/conditions. Subset saSPT now runs in
    # run_combinationwise.py against properly pooled per-condition zone
    # data (declared in combinations.yaml with a 'zone:' field).
    #
    # What filewise still writes for downstream combinationwise to use:
    #   trajectories/{basename}_traj.csv            (all)
    #   trajectories_on_structure/{basename}_traj.csv
    #   trajectories_on_ER_broad/{basename}_traj.csv
    #   trajectories_off_ER/{basename}_traj.csv
    #   er_segmentation/{basename}_traj_classification.csv
    #   er_segmentation/{basename}_er_masks.npz
    #   er_segmentation/{basename}_er_processed.npy
    # And the per-FOV global saSPT outputs (posterior/, MLE/,
    # rolling_windows/) — those feed run_conditionwise.py.

    print(f'\n=== Done: {file_stem} ===')


if __name__ == '__main__':
    main()
