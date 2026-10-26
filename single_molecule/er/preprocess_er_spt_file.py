#!/usr/bin/env python3
"""
preprocess_er_spt_file.py

Preprocess a single dual-channel ER + SPT ND2 file.
Called by SLURM for parallel processing.

Channel mapping (specific to your microscope setup):
  - Channel 0 = SPT (illuminated frames 100+)
  - Channel 1 = ER (illuminated frames 0-99)

Usage:
    python preprocess_er_spt_file.py -i input.nd2 -o output_dir/ [options]
"""

import os
import sys
import argparse
import numpy as np

import nd2reader
import tifffile
from PIL import Image


def process_file(input_path, output_dir, er_frames=100, spt_start=100, 
                 er_channel=1, spt_channel=0):
    """
    Extract ER reference image and SPT movie from dual-channel ND2.
    
    Parameters
    ----------
    input_path : str
        Path to input ND2 file
    output_dir : str
        Output directory
    er_frames : int
        Number of frames to use for ER projection (default: 100)
    spt_start : int
        Frame to start SPT extraction (default: 100)
    er_channel : int
        Channel index for ER (default: 1)
    spt_channel : int
        Channel index for SPT (default: 0)
    """
    
    basename = os.path.splitext(os.path.basename(input_path))[0]
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Processing: {input_path}")
    print(f"  ER channel: {er_channel}, SPT channel: {spt_channel}")
    print(f"  ER frames: 0-{er_frames-1}, SPT frames: {spt_start}+")
    
    with nd2reader.ND2Reader(input_path) as img:
        n_frames = img.sizes.get('t', len(img))
        n_channels = img.sizes.get('c', 1)
        
        print(f"  File: {n_channels} channels, {n_frames} frames")
        
        if n_channels < 2:
            raise ValueError(f"Expected 2 channels, found {n_channels}")
        
        img.bundle_axes = 'yx'
        
        # Extract ER channel - median projection of first er_frames
        print(f"  Extracting ER (channel {er_channel}, frames 0-{er_frames-1})...")
        er_stack = []
        for t in range(min(er_frames, n_frames)):
            frame = np.nan_to_num(img.get_frame_2D(c=er_channel, t=t))
            er_stack.append(frame)
        
        er_stack = np.array(er_stack)
        er_median = np.median(er_stack, axis=0).astype(np.float32)
        
        # Save ER as numpy array
        er_npy_path = os.path.join(output_dir, f"{basename}_er.npy")
        np.save(er_npy_path, er_median)
        print(f"    Saved: {er_npy_path}")
        
        # Save ER as PNG preview
        er_norm = (er_median - er_median.min()) / (er_median.max() - er_median.min() + 1e-8) * 255
        er_png_path = os.path.join(output_dir, f"{basename}_er.png")
        Image.fromarray(er_norm.astype(np.uint8)).save(er_png_path)
        print(f"    Saved: {er_png_path}")
        
        # Save ER as TIFF (for Fiji)
        er_tif_path = os.path.join(output_dir, f"{basename}_er.tif")
        tifffile.imwrite(er_tif_path, er_median)
        print(f"    Saved: {er_tif_path}")
        
        # Extract SPT channel - frames spt_start onwards
        print(f"  Extracting SPT (channel {spt_channel}, frames {spt_start}-{n_frames-1})...")
        spt_frames = []
        for t in range(spt_start, n_frames):
            frame = np.nan_to_num(img.get_frame_2D(c=spt_channel, t=t))
            spt_frames.append(frame)
        
        spt_stack = np.array(spt_frames).astype(np.uint16)
        
        # Save SPT as TIFF for quot
        spt_tif_path = os.path.join(output_dir, f"{basename}_spt.tif")
        tifffile.imwrite(spt_tif_path, spt_stack, imagej=True)
        print(f"    Saved: {spt_tif_path} ({spt_stack.shape[0]} frames)")
    
    print("  Done!")
    return er_npy_path, spt_tif_path


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess dual-channel ER+SPT ND2 file"
    )
    parser.add_argument("-i", "--input", required=True, help="Input ND2 file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--er-frames", type=int, default=100,
                        help="Number of ER frames (default: 100)")
    parser.add_argument("--spt-start", type=int, default=100,
                        help="SPT start frame (default: 100)")
    parser.add_argument("--er-channel", type=int, default=1,
                        help="ER channel index (default: 1)")
    parser.add_argument("--spt-channel", type=int, default=0,
                        help="SPT channel index (default: 0)")
    
    args = parser.parse_args()
    
    process_file(
        args.input,
        args.output,
        er_frames=args.er_frames,
        spt_start=args.spt_start,
        er_channel=args.er_channel,
        spt_channel=args.spt_channel
    )


if __name__ == "__main__":
    main()
