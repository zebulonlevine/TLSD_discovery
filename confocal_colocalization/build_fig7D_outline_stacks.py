"""
build_fig7D_outline_stacks.py -- Figure 7D: the 23 h GFP-R3HDM4 confocal
stacks as Fiji hyperstacks with the Cellpose cell boundaries added as an
extra channel.

ND2 can't be written back reliably, so the "modified ND2" is an ImageJ
hyperstack TIFF: the three original channels untouched (uint16, Nikon
batch-denoised, as in the source ND2) plus a 4th "cell outline" channel.

Masks are NOT re-segmented. They are the existing Cellpose-SAM labels in
<confocal data root>/output/cell_masks/ (segment_cells.py: cpsam,
diameter 260, 2D + stitch 0.5, channels 0-2 combined) -- the same labels
wholecell_pearson.py quantifies, so the outlines drawn are the
cells that were measured.

Per file, two outputs in <output root>/confocal_colocalization/fig7D_outline_stacks/:
  {tag}_outlines_hyperstack.tif  Z x C(4) x Y x X. Outline channel holds
      each slice's own cross-section boundary, so any single slice is
      correct. (A max projection of this channel overlays every slice's
      boundary and looks busy -- use the MIP file for a projection.)
  {tag}_outlines_MIP.tif  C(4) x Y x X. Channels 0-2 are the max
      projection; channel 3 is the boundary of the label max projection
      (the same construction as the pipeline's qc_overlay).

Outline pixels = 65535, else 0 (single dtype across a hyperstack).
Calibration: 0.1625 um/px, 0.4 um z-step (ND2 metadata).

Run:  python build_fig7D_outline_stacks.py   (nd2_coloc env; ~4 GB written, run under Slurm)
"""

import os
import sys

import nd2
import numpy as np
import tifffile
from skimage.segmentation import find_boundaries

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots

SRC = tlsd_paths.confocal('Denoised')
MASKS = tlsd_paths.confocal('output', 'cell_masks')
OUT = tlsd_paths.outp('confocal_colocalization', 'fig7D_outline_stacks')

TAGS = [
    'Batch Denoise_000_L29_GFP-R3HDM4_DMSO_23h_40x_live_immediate_003',
    'Batch Denoise_000_L29_GFP-R3HDM4_PT179_23h_40x_live_immediate_004',
]

OUTLINE = np.uint16(65535)

# Fiji LUTs (3 x 256): green, magenta, cyan for the data, white for outlines.
_ramp = np.arange(256, dtype=np.uint8)
_zero = np.zeros(256, np.uint8)
LUTS = [np.stack([_zero, _ramp, _zero]),
        np.stack([_ramp, _zero, _ramp]),
        np.stack([_zero, _ramp, _ramp]),
        np.stack([_ramp, _ramp, _ramp])]


def boundaries_per_slice(labels):
    """Inner boundary of every cell, computed within each z-slice."""
    out = np.zeros(labels.shape, dtype=np.uint16)
    for z in range(labels.shape[0]):
        out[z][find_boundaries(labels[z], mode='inner')] = OUTLINE
    return out


def write_hyperstack(path, data, axes, xy_um, z_um, names):
    meta = {'axes': axes, 'unit': 'um', 'mode': 'composite',
            'Labels': names, 'LUTs': LUTS}
    if z_um is not None:
        meta['spacing'] = z_um
    tifffile.imwrite(path, data, imagej=True,
                     resolution=(1 / xy_um, 1 / xy_um), metadata=meta)


def main():
    os.makedirs(OUT, exist_ok=True)
    for tag in TAGS:
        with nd2.ND2File(os.path.join(SRC, f'{tag}.nd2')) as f:
            assert list(f.sizes) == ['Z', 'C', 'Y', 'X'], f.sizes
            vox = f.voxel_size()
            names = [c.channel.name for c in f.metadata.channels]
            img = f.asarray()                                  # Z C Y X
        labels = tifffile.imread(os.path.join(MASKS, f'{tag}_cell_labels.tif'))
        assert labels.shape == (img.shape[0], *img.shape[2:]), \
            (labels.shape, img.shape)
        assert img.dtype == np.uint16, img.dtype
        names = names + ['cell outline']

        stack = np.concatenate([img, boundaries_per_slice(labels)[:, None]],
                               axis=1)
        p = os.path.join(OUT, f'{tag}_outlines_hyperstack.tif')
        write_hyperstack(p, stack, 'ZCYX', vox.x, vox.z, names)
        print(f'{tag}: {int(labels.max())} cells, {stack.shape} -> {p}')

        mip = img.max(axis=0)
        mip_outline = np.zeros(mip.shape[1:], np.uint16)
        mip_outline[find_boundaries(labels.max(axis=0), mode='inner')] = OUTLINE
        p = os.path.join(OUT, f'{tag}_outlines_MIP.tif')
        write_hyperstack(p, np.concatenate([mip, mip_outline[None]]), 'CYX',
                         vox.x, None, names)
        print(f'  MIP -> {p}  (xy {vox.x} um, z {vox.z} um, '
              f'channels {names})')


if __name__ == '__main__':
    main()
