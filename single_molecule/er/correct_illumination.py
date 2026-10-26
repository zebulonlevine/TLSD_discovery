#!/usr/bin/env python3
"""
correct_illumination.py

Correct uneven illumination in ER images before segmentation.

Changes from prior version:
  - REMOVED per-image 0-100 min-max rescaling that destroyed cross-file
    comparability. Corrected values now stay in natural units.
  - Added --global-normalize mode: compute global percentiles across all
    files first, then apply consistent rescaling so cross-file comparisons
    are valid.
  - Added --parallel flag for SLURM array compatibility (single-file mode).

Usage:
    # Single file preview
    python correct_illumination.py --er-file sample_er.npy --preview

    # Directory mode (natural units, no rescaling)
    python correct_illumination.py --er-dir preprocessed/ --output corrected/

    # Directory mode with global normalization (for cross-file comparison)
    python correct_illumination.py --er-dir preprocessed/ --output corrected/ --global-normalize

    # Single-file mode for SLURM array parallelization
    python correct_illumination.py --er-file preprocessed/sample_er.npy --output corrected/ --method gaussian
"""

import os
import glob
import argparse
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.ndimage import uniform_filter1d
import warnings


# ============================================================
# Correction methods
# ============================================================

def correct_rolling_ball(image, radius=50):
    """
    Rolling ball background subtraction (like ImageJ).
    Good for general uneven illumination.
    """
    from scipy.ndimage import grey_opening

    struct = np.ones((radius, radius))
    background = grey_opening(image, structure=struct)
    corrected = image - background

    return corrected, background


def correct_row_normalize(image, smooth_window=50):
    """
    Normalize each row by its smoothed mean.
    Good when illumination varies primarily along Y.
    """
    row_means = np.mean(image, axis=1)
    row_means_smooth = uniform_filter1d(row_means, size=smooth_window)
    row_means_smooth = np.maximum(row_means_smooth, 1)

    global_mean = np.mean(image)
    correction = global_mean / row_means_smooth
    corrected = image * correction[:, np.newaxis]

    return corrected, row_means_smooth


def correct_polynomial_surface(image, degree=2):
    """
    Fit polynomial surface to background and subtract.
    Good for smooth gradients.
    """
    ny, nx = image.shape
    y, x = np.mgrid[0:ny, 0:nx]

    x_flat = x.ravel()
    y_flat = y.ravel()
    z_flat = image.ravel()

    # Build design matrix for 2D polynomial
    terms = []
    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            terms.append((x_flat ** i) * (y_flat ** j))

    A = np.vstack(terms).T
    coeffs, _, _, _ = np.linalg.lstsq(A, z_flat, rcond=None)

    background = np.zeros_like(z_flat)
    for k, (i, j) in enumerate(
        [(i, j) for i in range(degree + 1) for j in range(degree + 1 - i)]
    ):
        background += coeffs[k] * (x_flat ** i) * (y_flat ** j)

    background = background.reshape(image.shape)
    corrected = image - background + np.mean(background)

    return corrected, background


def correct_gaussian_divide(image, sigma=100):
    """
    Divide by heavily blurred image (flat-field style).
    Good for multiplicative illumination variation.
    """
    blurred = ndimage.gaussian_filter(image.astype(float), sigma=sigma)
    blurred = np.maximum(blurred, 1)

    corrected = image / blurred
    corrected = corrected * np.mean(image)  # restore roughly original scale

    return corrected, blurred


def correct_percentile_rows(image, percentile=10, smooth_window=30):
    """
    Subtract row-wise background estimated from low percentile.
    Robust to cells — uses darkest pixels as background estimate.
    """
    row_background = np.percentile(image, percentile, axis=1)
    row_background_smooth = uniform_filter1d(row_background, size=smooth_window)

    corrected = image - row_background_smooth[:, np.newaxis]
    corrected = np.maximum(corrected, 0)

    return corrected, row_background_smooth


# ============================================================
# Apply correction (shared logic)
# ============================================================

METHODS = {
    'row_normalize': correct_row_normalize,
    'percentile_rows': correct_percentile_rows,
    'polynomial': lambda img: correct_polynomial_surface(img, degree=2),
    'gaussian': lambda img: correct_gaussian_divide(img, sigma=100),
    'rolling_ball': lambda img: correct_rolling_ball(img, radius=50),
}


def apply_correction(image, method='gaussian'):
    """Apply a correction method by name. Returns (corrected, profile)."""
    if method not in METHODS:
        raise ValueError(f"Unknown method: {method}. Choose from {list(METHODS.keys())}")
    return METHODS[method](image)


# ============================================================
# Preview mode
# ============================================================

def preview_correction(er_path, method='gaussian'):
    """Test correction on single file and save comparison."""

    er_image = np.load(er_path)

    print(f"Image shape: {er_image.shape}")
    print(f"Original range: {er_image.min():.1f} - {er_image.max():.1f}")

    corrected, profile = apply_correction(er_image, method)
    print(f"Corrected range: {corrected.min():.1f} - {corrected.max():.1f}")

    # Save comparison image (side by side)
    def normalize_for_display(img):
        img = img - np.percentile(img, 1)
        img = img / (np.percentile(img, 99) + 1e-8) * 255
        return np.clip(img, 0, 255).astype(np.uint8)

    orig_display = normalize_for_display(er_image)
    corr_display = normalize_for_display(corrected)
    comparison = np.hstack([orig_display, corr_display])

    output_dir = os.path.dirname(er_path)
    basename = os.path.basename(er_path).replace('_er.npy', '')

    comp_path = os.path.join(output_dir, f"{basename}_correction_compare_{method}.png")
    Image.fromarray(comparison).save(comp_path)
    print(f"Comparison saved: {comp_path}")

    corr_npy_path = os.path.join(output_dir, f"{basename}_er_corrected.npy")
    np.save(corr_npy_path, corrected.astype(np.float32))
    print(f"Corrected image saved: {corr_npy_path}")

    # Plot row profiles
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        axes[0].plot(np.mean(er_image, axis=1), label='Original', alpha=0.7)
        axes[0].plot(np.mean(corrected, axis=1), label='Corrected', alpha=0.7)
        axes[0].set_xlabel('Row (Y position)')
        axes[0].set_ylabel('Mean intensity')
        axes[0].legend()
        axes[0].set_title('Row-wise intensity profile')

        axes[1].hist(er_image.ravel(), bins=100, alpha=0.5, label='Original', density=True)
        axes[1].hist(corrected.ravel(), bins=100, alpha=0.5, label='Corrected', density=True)
        axes[1].set_xlabel('Intensity')
        axes[1].set_ylabel('Density')
        axes[1].legend()
        axes[1].set_title('Intensity distribution')

        plt.tight_layout()
        profile_path = os.path.join(output_dir, f"{basename}_correction_profile_{method}.png")
        plt.savefig(profile_path, dpi=150)
        plt.close()
        print(f"Profile plot saved: {profile_path}")
    except ImportError:
        print("Matplotlib not available for profile plot")


# ============================================================
# Directory processing
# ============================================================

def process_single_file(er_path, output_dir, method='gaussian'):
    """Correct a single ER file. Used for SLURM parallelization."""
    os.makedirs(output_dir, exist_ok=True)
    basename = os.path.basename(er_path)

    er_image = np.load(er_path)
    corrected, _ = apply_correction(er_image, method)

    # NO rescaling — keep natural corrected units
    out_path = os.path.join(output_dir, basename)
    np.save(out_path, corrected.astype(np.float32))
    print(f"  {basename}: range [{corrected.min():.1f}, {corrected.max():.1f}]")
    return out_path


def process_directory(er_dir, output_dir, method='gaussian', global_normalize=False):
    """Apply correction to all ER images in a directory."""

    os.makedirs(output_dir, exist_ok=True)

    er_files = sorted(glob.glob(os.path.join(er_dir, "*_er.npy")))
    if not er_files:
        print(f"No *_er.npy files found in {er_dir}")
        return

    print(f"Processing {len(er_files)} files with method: {method}")
    print(f"Global normalize: {global_normalize}")

    # First pass: correct all files
    corrected_images = {}
    for er_path in er_files:
        basename = os.path.basename(er_path)
        stem = basename.replace('_er.npy', '')
        print(f"  {stem}...", end=" ")

        er_image = np.load(er_path)
        corrected, _ = apply_correction(er_image, method)
        corrected_images[basename] = corrected
        print(f"range [{corrected.min():.1f}, {corrected.max():.1f}]")

    # Optional: global percentile normalization
    if global_normalize:
        print("\nComputing global percentiles for normalization...")
        all_values = np.concatenate([img.ravel() for img in corrected_images.values()])
        p1 = np.percentile(all_values, 1)
        p99 = np.percentile(all_values, 99)
        print(f"  Global 1st percentile: {p1:.2f}")
        print(f"  Global 99th percentile: {p99:.2f}")

        for basename in corrected_images:
            img = corrected_images[basename]
            img = (img - p1) / (p99 - p1 + 1e-8) * 100.0
            img = np.clip(img, 0, None)  # floor at 0, allow above 100
            corrected_images[basename] = img

    # Save all
    for basename, corrected in corrected_images.items():
        out_path = os.path.join(output_dir, basename)
        np.save(out_path, corrected.astype(np.float32))

    print(f"\nCorrected images saved to: {output_dir}")


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Correct uneven ER illumination")
    parser.add_argument("--er-file", help="Single ER .npy file")
    parser.add_argument("--er-dir", help="Directory of ER images")
    parser.add_argument("--output", "-o", help="Output directory")
    parser.add_argument("--method", default="gaussian",
                        choices=list(METHODS.keys()),
                        help="Correction method (default: gaussian)")
    parser.add_argument("--preview", action="store_true", help="Preview mode with comparison images")
    parser.add_argument("--global-normalize", action="store_true",
                        help="Apply global percentile normalization across all files "
                             "(ensures cross-file comparability)")

    args = parser.parse_args()

    if args.preview and args.er_file:
        preview_correction(args.er_file, args.method)
    elif args.er_file and args.output:
        # Single-file mode (for SLURM parallel)
        process_single_file(args.er_file, args.output, args.method)
    elif args.er_dir:
        output_dir = args.output or os.path.join(args.er_dir, 'corrected')
        process_directory(args.er_dir, output_dir, args.method,
                          global_normalize=args.global_normalize)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
