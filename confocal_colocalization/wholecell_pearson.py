#!/usr/bin/env python
"""
wholecell_pearson.py — Per-cell whole-cell Pearson correlation (slim)
=====================================================================
Slim replacement for scripts/quantify_coloc.py: keeps only the whole-cell
Pearson calculation (no ER compartments, no Manders/ICQ/Spearman). Numbers
are identical to the `whole_cell` rows of quantify_coloc.py output.

Per image:
  1. Clip each channel to the [0.1, 99.9] percentile of its nonzero voxels
  2. Drop cells smaller than --min-cell-voxels
  3. Pearson r between every channel pair over all voxels in each 3D cell label

Inputs (per position, written by scripts/segment_cells.py):
  - {tag}_cell_labels.tif   in <input>/cell_masks/
  - {tag}_ch{N}.tif         in <input>/channels/

Output:
  - {tag}_wholecell_pearson.csv — one row per cell × channel pair

Usage:
  python wholecell_pearson.py \\
    --input ../output --output per_cell --index 0 \\
    --channel-names GFP-R3HDM4,Sec61B,L29 --pattern "*GFP-R3HDM4*"
"""

import argparse
import sys
import time
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import find_objects
from scipy.stats import pearsonr
from tifffile import imread

warnings.filterwarnings("ignore", category=RuntimeWarning)


def parse_args():
    p = argparse.ArgumentParser(description="Per-cell whole-cell Pearson r")
    p.add_argument("--input", required=True,
                   help="Base directory containing cell_masks/ and channels/")
    p.add_argument("--output", required=True, help="Output directory for per-image CSVs")
    p.add_argument("--index", type=int, default=None, help="Position index (SLURM array)")
    p.add_argument("--channels", default="0,1,2")
    p.add_argument("--channel-names", default=None,
                   help="Comma-separated names, same order as --channels")
    p.add_argument("--pattern", default="*",
                   help="Glob on the image tag to select positions")
    p.add_argument("--min-cell-voxels", type=int, default=500)
    p.add_argument("--outlier-clip", default="0.1,99.9",
                   help="Percentile range for per-channel clipping")
    return p.parse_args()


def find_positions(base_dir, ch_indices, pattern):
    """Match cell label masks with per-channel TIFs by tag."""
    mask_dir = Path(base_dir) / "cell_masks"
    ch_dir = Path(base_dir) / "channels"
    positions = []
    for mf in sorted(mask_dir.glob(f"{pattern}_cell_labels.tif")):
        tag = mf.stem.replace("_cell_labels", "")
        pos = {"tag": tag, "cell_labels": mf}
        ch_paths = [ch_dir / f"{tag}_ch{ci}.tif" for ci in ch_indices]
        missing = [p for p in ch_paths if not p.exists()]
        if missing:
            print(f"  WARNING: skipping {tag}, missing {missing}")
            continue
        for ci, cp in zip(ch_indices, ch_paths):
            pos[f"ch{ci}"] = cp
        positions.append(pos)
    return positions


def clip_channel(vol, clip_lo, clip_hi):
    """Clip to percentiles of nonzero voxels (as in quantify_coloc.py)."""
    vol = vol.astype(np.float32)
    nz = vol[vol > 0]
    lo, hi = np.percentile(nz, [clip_lo, clip_hi]) if len(nz) > 0 else (0, 1)
    return np.clip(vol, lo, hi)


def compute_pearson(a, b):
    """Pearson r and p, NaN for degenerate cases."""
    if len(a) < 10 or np.std(a) < 1e-10 or np.std(b) < 1e-10:
        return np.nan, np.nan
    try:
        return pearsonr(a, b)
    except Exception:
        return np.nan, np.nan


def quantify_position(cell_labels, clipped, ch_indices, ch_names, min_voxels):
    rows = []
    # Bounding boxes per label: avoids a full-volume comparison per cell
    for cid, sl in enumerate(find_objects(cell_labels), start=1):
        if sl is None:
            continue
        cell_mask = cell_labels[sl] == cid
        n_vox = int(cell_mask.sum())
        if n_vox < min_voxels:
            continue
        vals = {ci: clipped[ci][sl][cell_mask] for ci in ch_indices}
        for ci_a, ci_b in combinations(ch_indices, 2):
            r, p = compute_pearson(vals[ci_a], vals[ci_b])
            rows.append({
                "cell_id": cid,
                "n_voxels": n_vox,
                "ch_a": ch_names[ch_indices.index(ci_a)],
                "ch_b": ch_names[ch_indices.index(ci_b)],
                "pearson_r": r,
                "pearson_p": p,
            })
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    ch_indices = [int(c) for c in args.channels.split(",")]
    ch_names = (args.channel_names.split(",") if args.channel_names
                else [f"ch{i}" for i in ch_indices])
    assert len(ch_names) == len(ch_indices), "channel names/indices mismatch"
    clip_lo, clip_hi = [float(x) for x in args.outlier_clip.split(",")]

    positions = find_positions(args.input, ch_indices, args.pattern)
    if not positions:
        print(f"ERROR: no matching positions in {args.input}")
        sys.exit(1)

    if args.index is not None:
        if args.index >= len(positions):
            print(f"Index {args.index} out of range ({len(positions)} positions)")
            sys.exit(0)
        positions = [positions[args.index]]

    print(f"Channels: {list(zip(ch_indices, ch_names))}")
    for pos in positions:
        t0 = time.time()
        tag = pos["tag"]
        print(f"\n--- {tag} ---")
        cell_labels = imread(str(pos["cell_labels"]))
        clipped = {ci: clip_channel(imread(str(pos[f"ch{ci}"])), clip_lo, clip_hi)
                   for ci in ch_indices}

        df = quantify_position(cell_labels, clipped, ch_indices, ch_names,
                               args.min_cell_voxels)
        df.insert(0, "image", tag)
        df.to_csv(out_dir / f"{tag}_wholecell_pearson.csv", index=False)
        n_cells = df["cell_id"].nunique() if len(df) else 0
        print(f"  {n_cells} cells kept ({int(cell_labels.max())} labels), "
              f"{time.time() - t0:.0f} s")

    print("\nDone.")


if __name__ == "__main__":
    main()
