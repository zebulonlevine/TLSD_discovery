#!/usr/bin/env python3
"""
rerun_er_zoning.py

Re-run ONLY the ER segmentation + trajectory classification + zone-subset
export, reading existing per-FOV trajectories. Does NOT re-track (quot)
and does NOT run saSPT. Use this to iterate on segmentation parameters
(pct_floor, require_broad, pct_above, ...) without paying for a full
pipeline re-run.

For each FOV it reproduces filewise Steps 2-3 + the subset export:
  - load preprocessed_corrected/{stem}_er.npy  (or reuse existing
    er_segmentation/{stem}_er_processed.npy to skip re-deconvolution)
  - three_class_mask(...) with the new parameters
  - classify_trajectories(...) on existing trajectories/{stem}_traj.csv
  - OVERWRITE er_segmentation/{stem}_er_masks.npz
  - OVERWRITE er_segmentation/{stem}_er_processed.npy (only if recomputed)
  - OVERWRITE er_segmentation/{stem}_traj_classification.csv
  - OVERWRITE trajectories_{on_structure,on_ER_broad,off_ER}/{stem}_traj.csv

It reuses the EXACT classify_trajectories and export_subset_traj logic from
run_filewise_er.py (loaded without triggering quot import), so the zones it
writes are byte-identical to what a full filewise run would write with the
same parameters.

Usage (serial, all FOVs in an analysis dir):
    python rerun_er_zoning.py \\
        --analysis-dir /path/to/tracking_output_q=15p0 \\
        --er-dir /path/to/preprocessed_corrected \\
        --pct-floor 60 --require-broad false

Usage (SLURM array — one FOV per task):
    python rerun_er_zoning.py --analysis-dir ... --er-dir ... \\
        --pct-floor 60 --require-broad false --fov-index $SLURM_ARRAY_TASK_ID

Parameters default to filewise defaults except where you override.
"""

import argparse
import glob
import importlib.util
import os
import re
import sys

import numpy as np
import pandas as pd


ER_CODE = os.path.dirname(os.path.abspath(__file__))


def _load_filewise_helpers(filewise_path):
    """
    Import classify_trajectories, export_subset_traj, _three_class_mask,
    and _process_er from run_filewise_er.py WITHOUT executing its
    module-level `from quot.core import track_file` (which is slow and
    may be unavailable outside the tracking environment).

    Strategy: read the source, strip the top-level quot import line, exec
    the rest in a throwaway module namespace, and pull the functions out.
    The helpers we need don't use track_file, so removing that import is
    safe.
    """
    with open(filewise_path) as f:
        src = f.read()

    # Comment out the top-level quot import so module exec doesn't require it.
    src = re.sub(r'^(from quot\.core import track_file)',
                 r'# \1  # stripped by rerun_er_zoning', src,
                 flags=re.MULTILINE)
    # Guard: make sure the __main__ block doesn't execute on import.
    # run_filewise_er.py guards main() under if __name__ == '__main__',
    # so exec under a non-__main__ name is safe. We give it a private name.
    mod = type(sys)('._filewise_helpers')
    mod.__file__ = filewise_path
    mod.__name__ = '_filewise_helpers'
    code = compile(src, filewise_path, 'exec')
    exec(code, mod.__dict__)
    return mod


def parse_bool(tok):
    return str(tok).strip().lower() in ('true', '1', 'yes', 't', 'y')


def discover_fov_stems(seg_dir, traj_dir):
    """FOV stems = basenames of existing trajectories/{stem}_traj.csv."""
    traj_files = sorted(glob.glob(os.path.join(traj_dir, '*_traj.csv')))
    return [os.path.basename(p).replace('_traj.csv', '') for p in traj_files]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--analysis-dir', required=True,
                    help='tracking_output_q=*/ ; holds trajectories/, '
                         'er_segmentation/, trajectories_{zone}/')
    ap.add_argument('--er-dir', required=True,
                    help='preprocessed_corrected/ with {stem}_er.npy')
    ap.add_argument('--filewise', default=os.path.join(ER_CODE, 'run_filewise_er.py'),
                    help='Path to run_filewise_er.py (for helper reuse).')
    # Segmentation params (defaults mirror filewise argparse defaults)
    ap.add_argument('--pct-above', type=float, default=1.5)
    ap.add_argument('--pct-above-strict', type=float, default=None)
    ap.add_argument('--block-fine', type=int, default=25)
    ap.add_argument('--block-broad', type=int, default=101)
    ap.add_argument('--min-object-size', type=int, default=500)
    ap.add_argument('--pct-floor', type=float, default=50,
                    help='Global percentile floor; <0 disables.')
    ap.add_argument('--require-broad', type=parse_bool, default=True,
                    help="true|false. false => on_structure is fine+floor "
                         "only, no broad-density gate.")
    # Classification params
    ap.add_argument('--frac-threshold', type=float, default=0.95)
    ap.add_argument('--enrichment-alpha', type=float, default=0.05)
    ap.add_argument('--min-frames', type=int, default=2,
                    help='Structural floor (default 2). Derived per-FOV n_min '
                         'rules the classification cutoff; raise only to sweep '
                         'short-trajectory sensitivity.')
    ap.add_argument('--pixel-size', type=float, default=0.11)
    # Subsets to export
    ap.add_argument('--subsets', type=str,
                    default='on_structure,on_ER_broad,off_ER')
    # Behavior
    ap.add_argument('--reprocess-er', action='store_true',
                    help='Recompute _er_processed.npy from _er.npy via '
                         'flatfield+deconv. Default: reuse existing '
                         '_er_processed.npy if present (faster; valid '
                         'because only mask params changed, not processing).')
    ap.add_argument('--fov-index', type=int, default=None,
                    help='If set, process only the Nth FOV (for SLURM array).')
    args = ap.parse_args()

    pct_floor = None if args.pct_floor < 0 else args.pct_floor
    subsets = [s.strip() for s in args.subsets.split(',')]

    seg_dir = os.path.join(args.analysis_dir, 'er_segmentation')
    traj_dir = os.path.join(args.analysis_dir, 'trajectories')
    os.makedirs(seg_dir, exist_ok=True)
    for sub in subsets:
        os.makedirs(os.path.join(args.analysis_dir, f'trajectories_{sub}'),
                    exist_ok=True)

    # Load filewise helpers (without quot)
    fw = _load_filewise_helpers(args.filewise)
    classify_trajectories = fw.classify_trajectories
    export_subset_traj = fw.export_subset_traj
    _process_er = fw._process_er

    fov_stems = discover_fov_stems(seg_dir, traj_dir)
    if not fov_stems:
        raise SystemExit(f'No *_traj.csv in {traj_dir}')

    if args.fov_index is not None:
        if args.fov_index >= len(fov_stems):
            raise SystemExit(f'fov-index {args.fov_index} out of range '
                             f'(have {len(fov_stems)} FOVs)')
        fov_stems = [fov_stems[args.fov_index]]

    print(f'Re-zoning {len(fov_stems)} FOV(s)')
    print(f'  pct_floor={pct_floor}  require_broad={args.require_broad}  '
          f'pct_above={args.pct_above}  pct_above_strict={args.pct_above_strict}')
    print(f'  reprocess_er={args.reprocess_er}')

    for stem in fov_stems:
        print(f'\n=== {stem} ===')
        traj_path = os.path.join(traj_dir, f'{stem}_traj.csv')
        if not os.path.exists(traj_path):
            print(f'  skip: no traj CSV at {traj_path}')
            continue

        # --- processed ER: reuse or recompute ---
        proc_path = os.path.join(seg_dir, f'{stem}_er_processed.npy')
        if args.reprocess_er or not os.path.exists(proc_path):
            # Find the raw corrected ER. filewise names it {stem}_er.npy,
            # but stem may carry a trailing _spt; normalize.
            er_basename = re.sub(r'_spt$', '', stem) + '_er.npy'
            er_path = os.path.join(args.er_dir, er_basename)
            if not os.path.exists(er_path):
                # fall back to {stem}_er.npy as-is
                er_path = os.path.join(args.er_dir, f'{stem}_er.npy')
            if not os.path.exists(er_path):
                print(f'  skip: no ER npy ({er_basename} or {stem}_er.npy)')
                continue
            er_raw = np.load(er_path).astype(np.float64)
            er_processed = _process_er(er_raw)
            np.save(proc_path, er_processed.astype(np.float32))
            print(f'  processed ER from {os.path.basename(er_path)} (recomputed)')
        else:
            er_processed = np.load(proc_path).astype(np.float64)
            print(f'  reused existing {os.path.basename(proc_path)}')

        # --- classify (this internally calls _three_class_mask with our params) ---
        traj_df = pd.read_csv(traj_path)
        traj_df, per_traj, masks = classify_trajectories(
            traj_df, er_processed, args.pixel_size,
            pct_above=args.pct_above,
            block_fine=args.block_fine,
            block_broad=args.block_broad,
            min_size=args.min_object_size,
            pct_floor=pct_floor,
            require_broad=args.require_broad,
            fraction_threshold=args.frac_threshold,
            enrichment_alpha=args.enrichment_alpha,
            min_frames=args.min_frames,
            pct_above_strict=args.pct_above_strict)

        if masks is None:
            print('  WARNING: segmentation failed (skimage?); skipping FOV')
            continue

        # --- overwrite masks ---
        np.savez_compressed(
            os.path.join(seg_dir, f'{stem}_er_masks.npz'),
            on_structure=masks['on_structure'],
            near_ER=masks['near_ER'],
            off_ER=masks['off_ER'])
        a_on = masks['on_structure'].mean()
        a_near = masks['near_ER'].mean()
        a_off = masks['off_ER'].mean()
        print(f'  area on/near/off: {a_on:.1%} / {a_near:.1%} / {a_off:.1%}')

        # --- overwrite classification ---
        per_traj['source_file'] = f'{stem}_spt.tif'
        per_traj.to_csv(
            os.path.join(seg_dir, f'{stem}_traj_classification.csv'),
            index=False)
        nt = len(per_traj)
        if nt:
            n_on = (per_traj['er_zone_label'] == 'on_structure').sum()
            n_near = (per_traj['er_zone_label'] == 'near_ER').sum()
            n_off = (per_traj['er_zone_label'] == 'off_ER').sum()
            print(f'  traj on/near/off: {n_on}/{n_near}/{n_off} of {nt} '
                  f'({n_on/nt:.1%} / {n_near/nt:.1%} / {n_off/nt:.1%})')

        # --- overwrite subset exports ---
        for subset_name in subsets:
            subset_path = os.path.join(
                args.analysis_dir, f'trajectories_{subset_name}',
                f'{stem}_traj.csv')
            n = export_subset_traj(traj_df, subset_name, subset_path)
            print(f'  exported {subset_name}: {n} traj')

    print('\nRe-zoning complete.')
    print('NEXT: wipe stale combos, then re-run combinationwise:')
    print(f'  rm {args.analysis_dir}/traj_combinations/*pooled*')
    print(f'  rm {args.analysis_dir}/posterior_combinations/*pooled*')
    print(f'  rm {args.analysis_dir}/MLE_combinations/*pooled*')
    print(f'  rm {args.analysis_dir}/rolling_windows_combinations/*pooled*')


if __name__ == '__main__':
    main()
