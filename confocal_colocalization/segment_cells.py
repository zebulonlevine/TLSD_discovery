#!/usr/bin/env python
"""
segment_cells.py — 3D Cell Segmentation with Cellpose-SAM (No Nuclear Stain)
==============================================================================
Reads ND2 files, sums all channels for a combined
cytoplasmic signal, and runs Cellpose-SAM for whole-cell 3D segmentation.

Designed for SLURM array parallelization: each job processes one file.

Inputs:
  - ND2 files (multi-channel z-stacks)

Outputs (per file, per position):
  - {tag}_cell_labels.tif       — 3D cell label mask (int32)
  - {tag}_qc_overlay.tif        — max-projection RGB with cell outlines
  - {tag}_qc_outlines_3d.tif    — 3D binary outline mask for Fiji overlay
  - {tag}_seg_meta.json         — metadata (n_cells, diameter, etc.)

Usage:
  python segment_cells.py \\
    --input /path/to/nd2_dir \\
    --output /path/to/output_dir \\
    --index 0 \\
    [--channels 0,1,2] \\
    [--cellpose-diam 120] \\
    [--anisotropy 2.46] \\
    [--min-cell-vol 500] \\
    [--flow-threshold 0.4] \\
    [--cellprob-threshold 0.0] \\
    [--no-gpu]
"""

import argparse
import json
import sys
import time
import warnings
from itertools import product as iterproduct
from pathlib import Path

import numpy as np
from nd2 import ND2File
from tifffile import imwrite
from skimage.segmentation import find_boundaries
from skimage.morphology import remove_small_objects

warnings.filterwarnings("ignore")


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description="3D cell segmentation with Cellpose-SAM on combined channels"
    )
    p.add_argument("--input", required=True, help="Directory with ND2 files")
    p.add_argument("--output", required=True, help="Output directory")
    p.add_argument("--index", type=int, default=None,
                   help="File index for SLURM array")
    p.add_argument("--channels", default="0,1,2",
                   help="Comma-separated channel indices to combine for segmentation")
    p.add_argument("--cellpose-diam", type=int, default=120,
                   help="Expected cell diameter in XY pixels (measure in Fiji)")
    p.add_argument("--anisotropy", type=float, default=None,
                   help="Z/XY voxel ratio (auto-detected from ND2 if not set)")
    p.add_argument("--min-cell-vol", type=int, default=500,
                   help="Remove cells smaller than this many voxels")
    p.add_argument("--flow-threshold", type=float, default=0.4,
                   help="Cellpose flow threshold (lower = stricter)")
    p.add_argument("--cellprob-threshold", type=float, default=0.0,
                   help="Cellpose cell probability threshold")
    p.add_argument("--stitch-threshold", type=float, default=0.0,
                   help="If >0, use 2D+stitch instead of true 3D (try 0.5 if 3D fragments)")
    p.add_argument("--no-gpu", action="store_true")
    p.add_argument("--pattern", default="*.nd2",
                   help="Glob pattern for input files")
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════════
# ND2 loading
# ═══════════════════════════════════════════════════════════════════════════

def load_nd2(path):
    """Load ND2 → (data, axes_dict, anisotropy)."""
    with ND2File(str(path)) as f:
        data = f.asarray()
        axes = f.sizes
        try:
            voxel = f.voxel_size()
            anisotropy = voxel.z / voxel.x if voxel.z else None
        except Exception:
            anisotropy = None
    return data, axes, anisotropy


def reshape_to_pzcyx(data, axes):
    """
    Reshape ND2 array to (P, Z, C, Y, X).
    Handles common ND2 dimension orderings.
    """
    # Build axis string from axes dict
    axis_order = list(axes.keys())

    # Normalize axis names
    rename = {"T": "P", "S": "P", "M": "P"}  # treat time/series/multipoint as P
    axis_order = [rename.get(a, a) for a in axis_order]

    # If no P dimension, add one
    if "P" not in axis_order:
        data = data[np.newaxis, ...]
        axis_order = ["P"] + axis_order

    # If no Z dimension, add one
    if "Z" not in axis_order:
        # Find where to insert Z (after P, before C)
        p_idx = axis_order.index("P")
        data = np.expand_dims(data, axis=p_idx + 1)
        axis_order.insert(p_idx + 1, "Z")

    # Target order: P Z C Y X
    target = ["P", "Z", "C", "Y", "X"]
    if axis_order != target:
        perm = [axis_order.index(a) for a in target]
        data = np.transpose(data, perm)

    return data


# ═══════════════════════════════════════════════════════════════════════════
# Segmentation
# ═══════════════════════════════════════════════════════════════════════════

def build_cpsam_input(volumes, ch_indices):
    """
    Build Cellpose-SAM input: normalize each channel to 0-1 independently,
    then take per-pixel max across channels.
    Returns (Z, Y, X) single-channel volume.
    """
    normalized = []
    for idx in ch_indices:
        ch = volumes[idx].astype(np.float32)
        cmax = np.percentile(ch, 99.9)
        if cmax > 0:
            ch = np.clip(ch / cmax, 0, 1)
        normalized.append(ch)

    return np.maximum.reduce(normalized)


def run_cellpose_sam(img, diameter, anisotropy, flow_thr, cellprob_thr,
                     stitch_thr, use_gpu, min_vol):
    """Run Cellpose-SAM and return cleaned label mask."""
    from cellpose import models

    model = models.CellposeModel(pretrained_model="cpsam", gpu=use_gpu)

    # Decide 3D vs 2D+stitch
    do_3d = stitch_thr <= 0
    stitch_threshold = stitch_thr if not do_3d else 0.0

    print(f"  Cellpose params: diam={diameter}, do_3D={do_3d}, "
          f"stitch={stitch_threshold}, cellprob={cellprob_thr}, flow={flow_thr}")

    masks, flows, _ = model.eval(
        img,
        diameter=diameter,
        z_axis=0,         # (Z, Y, X) → Z is axis 0
        do_3D=do_3d,
        anisotropy=anisotropy,
        flow_threshold=flow_thr,
        cellprob_threshold=cellprob_thr,
        stitch_threshold=stitch_threshold,
        min_size=min_vol,
    )

    n_before = int(masks.max())
    print(f"  Cellpose found {n_before} cells (before size filter)")

    # Remove small objects and relabel consecutively
    if min_vol > 0 and n_before > 0:
        masks = remove_small_objects(masks, min_size=min_vol)
        # Relabel consecutively
        unique_ids = np.unique(masks)
        unique_ids = unique_ids[unique_ids > 0]
        relabeled = np.zeros_like(masks)
        for new_id, old_id in enumerate(unique_ids, 1):
            relabeled[masks == old_id] = new_id
        masks = relabeled

    return masks


# ═══════════════════════════════════════════════════════════════════════════
# QC outputs
# ═══════════════════════════════════════════════════════════════════════════

def save_qc_overlay(signal_vol, masks, out_path):
    """Max-projection RGB with cell outlines in green."""
    mip_signal = signal_vol.max(axis=0).astype(np.float32)
    p1, p99 = np.percentile(mip_signal, [1, 99.5])
    mip_signal = np.clip((mip_signal - p1) / (p99 - p1 + 1e-8), 0, 1)

    rgb = np.stack([mip_signal] * 3, axis=-1)

    mip_masks = masks.max(axis=0)
    bounds = find_boundaries(mip_masks, mode="outer")
    rgb[bounds] = [0, 1, 0]  # green outlines

    imwrite(str(out_path), (rgb * 255).astype(np.uint8), compression="zlib")


def save_3d_outlines(masks, out_path):
    """Save 3D binary outline mask for Fiji overlay."""
    outlines = np.zeros(masks.shape, dtype=np.uint8)
    for z in range(masks.shape[0]):
        bounds = find_boundaries(masks[z], mode="outer")
        outlines[z] = bounds.astype(np.uint8) * 255
    imwrite(str(out_path), outlines, compression="zlib")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()
    input_dir = Path(args.input)
    out_dir = Path(args.output) / "cell_masks"
    qc_dir = Path(args.output) / "cell_qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    ch_indices = [int(c) for c in args.channels.split(",")]

    # Find ND2 files
    nd2_files = sorted(input_dir.glob(args.pattern))
    if not nd2_files:
        print(f"ERROR: No files matching {args.pattern} in {input_dir}")
        sys.exit(1)

    # SLURM array: select one file
    if args.index is not None:
        if args.index >= len(nd2_files):
            print(f"Index {args.index} out of range ({len(nd2_files)} files)")
            sys.exit(0)
        nd2_files = [nd2_files[args.index]]

    use_gpu = not args.no_gpu

    for nd2_path in nd2_files:
        t0 = time.time()
        stem = nd2_path.stem
        print(f"\n{'='*60}")
        print(f"Processing: {stem}")
        print(f"{'='*60}")

        # Load
        data, axes, auto_aniso = load_nd2(nd2_path)
        anisotropy = args.anisotropy or auto_aniso or 2.46  # fallback
        print(f"  Shape: {data.shape}, axes: {axes}")
        print(f"  Anisotropy: {anisotropy:.2f}")

        # Reshape to (P, Z, C, Y, X)
        data = reshape_to_pzcyx(data, axes)
        n_pos = data.shape[0]
        print(f"  Positions: {n_pos}, Z: {data.shape[1]}, "
              f"C: {data.shape[2]}, YX: {data.shape[3]}x{data.shape[4]}")

        for pos_idx in range(n_pos):
            tag = f"{stem}_p{pos_idx:03d}" if n_pos > 1 else stem
            print(f"\n  --- Position {pos_idx} ({tag}) ---")

            vol = data[pos_idx]  # (Z, C, Y, X)

            # Extract channels
            volumes = [vol[:, c, :, :] for c in range(vol.shape[1])]
            print(f"  Channels available: {len(volumes)}, "
                  f"using {ch_indices} for segmentation")

            # Build Cellpose-SAM input
            img_3ch = build_cpsam_input(volumes, ch_indices)
            print(f"  CPSAM input shape: {img_3ch.shape}")

            # Run segmentation
            masks = run_cellpose_sam(
                img_3ch,
                diameter=args.cellpose_diam,
                anisotropy=anisotropy,
                flow_thr=args.flow_threshold,
                cellprob_thr=args.cellprob_threshold,
                stitch_thr=args.stitch_threshold,
                use_gpu=use_gpu,
                min_vol=args.min_cell_vol,
            )

            n_cells = int(masks.max())
            print(f"  Segmented {n_cells} cells")

            # Save label mask
            imwrite(str(out_dir / f"{tag}_cell_labels.tif"),
                    masks.astype(np.int32), compression="zlib")

            # Save QC (same normalize-then-max as segmentation input)
            normed = []
            for i in ch_indices:
                ch = volumes[i].astype(np.float32)
                cmax = np.percentile(ch, 99.9)
                normed.append(ch / cmax if cmax > 0 else ch)
            combined_signal = np.maximum.reduce(normed)
            save_qc_overlay(combined_signal, masks,
                            qc_dir / f"{tag}_qc_overlay.tif")
            save_3d_outlines(masks, qc_dir / f"{tag}_qc_outlines_3d.tif")

            # Also save per-channel TIFs for quantification downstream
            ch_out_dir = Path(args.output) / "channels"
            ch_out_dir.mkdir(parents=True, exist_ok=True)
            for c_idx in range(vol.shape[1]):
                imwrite(
                    str(ch_out_dir / f"{tag}_ch{c_idx}.tif"),
                    volumes[c_idx], compression="zlib"
                )

            # Metadata
            meta = {
                "tag": tag,
                "source_file": str(nd2_path.name),
                "position_index": pos_idx,
                "n_cells": n_cells,
                "shape_zyx": list(masks.shape),
                "anisotropy": anisotropy,
                "cellpose_diameter": args.cellpose_diam,
                "flow_threshold": args.flow_threshold,
                "cellprob_threshold": args.cellprob_threshold,
                "stitch_threshold": args.stitch_threshold,
                "min_cell_vol": args.min_cell_vol,
                "channels_combined": ch_indices,
                "elapsed_sec": round(time.time() - t0, 1),
            }
            (out_dir / f"{tag}_seg_meta.json").write_text(
                json.dumps(meta, indent=2)
            )

    print("\nSegmentation complete.")


if __name__ == "__main__":
    main()
