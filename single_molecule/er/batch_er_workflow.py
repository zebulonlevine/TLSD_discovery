#!/usr/bin/env python3
"""
batch_er_workflow.py

Complete SLURM workflow for ER + SPT dual-channel experiments.

Steps:
  1. Preprocess: Extract ER + SPT from ND2 (job array, parallelized)
  2. Correct illumination: Per-file SLURM jobs (parallelized)
  3. Track: quot/saspt filewise + conditionwise per condition
  4. ER Colocalization with MLE: Consolidated script with optional MLE merge
  5. Visualize: Generate plots including MLE × ER cross-analyses

Changes from prior version:
  - Preprocessing submitted as a single SLURM job array (one array ID to depend on)
  - Illumination correction parallelized: one job per file instead of serial
  - Uses consolidated run_er_colocalization.py (replaces 3 separate scripts)
  - Cleaner dependency management throughout

Directory structure expected:
    experiment/
        raw_data/
            ConditionA_001.nd2
            ConditionA_002.nd2
            ConditionB_001.nd2
            ...

Usage:
    python batch_er_workflow.py /path/to/experiment
    python batch_er_workflow.py --dry-run /path/to/experiment
    python batch_er_workflow.py --from-file experiments.txt
"""

import os
import sys
import glob
import re
import pickle
import argparse
import subprocess
import numpy as np
import yaml
from pathlib import Path


# =============================================================================
# CONFIGURATION - Edit these paths for your setup
# =============================================================================

PYTHON_PATH = sys.executable
ER_WORKFLOW_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIRECTORY = os.path.abspath(os.path.join(ER_WORKFLOW_DIR, '..', 'pipeline'))

# ER segmentation params forwarded to run_filewise_er.py.
# Edit these to change segmentation for a full pipeline run.
PCT_FLOOR = 60        # global percentile floor; <0 disables
REQUIRE_BROAD = False # on_structure = fine+floor only when False

# Cluster-specific; leave unset to use the Slurm defaults.
PARTITION = os.environ.get('TLSD_SLURM_PARTITION', '')
ACCOUNT = os.environ.get('TLSD_SLURM_ACCOUNT', '')

# ER extraction parameters
ER_FRAMES = 100
SPT_START_FRAME = 100
ER_CHANNEL = 1
SPT_CHANNEL = 0

# Illumination correction parameters
ILLUMINATION_METHOD = "gaussian"


# =============================================================================
sys.path.insert(0, CODE_DIRECTORY)
from init_settings import get_default_settings


def _sbatch_args():
    """--partition/--account only when configured."""
    return [a for a in (f'--partition={PARTITION}' if PARTITION else '',
                        f'--account={ACCOUNT}' if ACCOUNT else '') if a]


def _sbatch_header():
    return ''.join(f'#SBATCH {a}\n' for a in _sbatch_args())


def setup_directories(exp_dir, settings):
    """Create output directory structure."""
    t_val = settings['quot']['detect']['t']
    output_name = f"tracking_output_q={str(float(t_val)).replace('.', 'p')}"

    settings['io']['home_directory'] = exp_dir
    settings['io']['data_filepath'] = 'preprocessed'
    settings['io']['data_directory'] = os.path.join(exp_dir, 'preprocessed')
    settings['io']['analysis_directory'] = os.path.join(exp_dir, output_name)
    settings['io']['preprocessed_directory'] = os.path.join(exp_dir, 'preprocessed')
    settings['io']['corrected_er_directory'] = os.path.join(exp_dir, 'preprocessed_corrected')
    settings['io']['traj_directory'] = os.path.join(settings['io']['analysis_directory'], 'trajectories')
    settings['io']['post_directory'] = os.path.join(settings['io']['analysis_directory'], 'posterior')
    settings['io']['MLE_directory'] = os.path.join(settings['io']['analysis_directory'], 'MLE')
    settings['io']['plot_directory'] = os.path.join(settings['io']['analysis_directory'], 'plots')
    settings['io']['rolling_window_directory'] = os.path.join(settings['io']['analysis_directory'], 'rolling_windows')
    settings['io']['movie_directory'] = os.path.join(settings['io']['analysis_directory'], 'movies')
    settings['io']['subtrack_posterior_directory'] = os.path.join(settings['io']['analysis_directory'], 'subtrack_posteriors')
    settings['io']['code_directory'] = CODE_DIRECTORY

    for key, path in settings['io'].items():
        if key.endswith('_directory') and path:
            os.makedirs(path, exist_ok=True)

    logs_dir = os.path.join(settings['io']['analysis_directory'], 'logs')
    os.makedirs(logs_dir, exist_ok=True)

    return settings


def discover_files(raw_data_dir):
    """Find ND2 files and group by condition."""
    nd2_files = sorted(glob.glob(os.path.join(raw_data_dir, "*.nd2")))

    if not nd2_files:
        raise ValueError(f"No ND2 files found in {raw_data_dir}")

    conditions = {}
    for f in nd2_files:
        basename = os.path.basename(f)
        match = re.match(r"(.+?)_\d+\.nd2$", basename)
        cond = match.group(1) if match else basename.replace('.nd2', '')

        if cond not in conditions:
            conditions[cond] = []
        conditions[cond].append(f)

    return conditions


def sbatch(cmd, dry_run=False, dry_label=""):
    """Submit an sbatch command. Returns job ID string or dry-run label."""
    if dry_run:
        print(f"    [DRY] {dry_label}")
        return dry_label

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        job_id = result.stdout.strip()
        return job_id
    else:
        print(f"    ERROR: {result.stderr.strip()}")
        return None


def submit_experiment(exp_dir, dry_run=False):
    """Submit all SLURM jobs for one experiment."""

    exp_dir = os.path.abspath(exp_dir)
    exp_name = os.path.basename(exp_dir)
    raw_data_dir = os.path.join(exp_dir, 'raw_data')

    print(f"\n{'='*60}")
    print(f"Experiment: {exp_name}")
    print(f"{'='*60}")

    if not os.path.isdir(raw_data_dir):
        print(f"  ERROR: No raw_data/ directory found")
        return

    try:
        conditions = discover_files(raw_data_dir)
    except ValueError as e:
        print(f"  ERROR: {e}")
        return

    total_files = sum(len(f) for f in conditions.values())
    print(f"  Found {total_files} files in {len(conditions)} conditions:")
    for cond, files in conditions.items():
        print(f"    {cond}: {len(files)} files")

    # Initialize settings
    settings = get_default_settings()

    # Read settings_override.yaml if present
    override_path = os.path.join(exp_dir, 'settings_override.yaml')
    if os.path.exists(override_path):
        with open(override_path) as f:
            override_cfg = yaml.load(f, Loader=yaml.FullLoader)
        print(f"\n  Applying overrides from: {override_path}")

        # Import deep_update
        sys.path.insert(0, CODE_DIRECTORY)
        from trackingCodeFunctions import deep_update
        
        # Check if frame_interval is explicitly overridden
        has_fi_override = False
        if override_cfg:
            if override_cfg.get('saspt', {}).get('frame_interval') is not None:
                has_fi_override = True
            if override_cfg.get('quot', {}).get('track', {}).get('frame_interval') is not None:
                has_fi_override = True
        
        deep_update(settings, override_cfg)
        
        if has_fi_override:
            settings['override_frame_interval'] = True
            print(f"    frame_interval locked to {settings['saspt']['frame_interval']} s")
    else:
        print(f"\n  No settings_override.yaml in {exp_dir}")

    settings = setup_directories(exp_dir, settings)

    pkl_path = os.path.join(settings['io']['analysis_directory'], 'settings.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(settings, f)
    print(f"\n  Settings saved: {pkl_path}")

    logs_dir = os.path.join(settings['io']['analysis_directory'], 'logs')
    preprocessed_dir = settings['io']['preprocessed_directory']
    corrected_er_dir = settings['io']['corrected_er_directory']
    traj_dir = settings['io']['traj_directory']
    mle_dir = settings['io']['MLE_directory']

    # Flatten all ND2 files
    all_nd2 = []
    for files in conditions.values():
        all_nd2.extend(files)

    # =========================================================================
    # STEP 1: Preprocessing as a SLURM job array
    # =========================================================================
    print(f"\n  [Step 1] Submitting preprocessing job array ({len(all_nd2)} files)...")

    # Write file list for the array job
    filelist_path = os.path.join(logs_dir, 'preprocess_filelist.txt')
    with open(filelist_path, 'w') as f:
        for nd2 in all_nd2:
            f.write(nd2 + '\n')

    # Create a wrapper script for the array job
    array_script = os.path.join(logs_dir, 'preprocess_array.sh')
    with open(array_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --job-name=pp_{exp_name[:12]}\n")
        f.write(f"#SBATCH --output={logs_dir}/preprocess_%A_%a.log\n")
        f.write(_sbatch_header())
        f.write(f"#SBATCH --mem=16G\n")
        f.write(f"#SBATCH --cpus-per-task=1\n")
        f.write(f"#SBATCH --array=0-{len(all_nd2)-1}\n")
        f.write(f"\n")
        f.write(f"FILELIST={filelist_path}\n")
        f.write(f"ND2_PATH=$(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" $FILELIST)\n")
        f.write(f"\n")
        f.write(f"{PYTHON_PATH} {ER_WORKFLOW_DIR}/preprocess_er_spt_file.py \\\n")
        f.write(f"    -i \"$ND2_PATH\" -o {preprocessed_dir} \\\n")
        f.write(f"    --er-frames {ER_FRAMES} --spt-start {SPT_START_FRAME} \\\n")
        f.write(f"    --er-channel {ER_CHANNEL} --spt-channel {SPT_CHANNEL}\n")

    preprocess_job_id = sbatch(
        ["sbatch", "--parsable", array_script],
        dry_run=dry_run,
        dry_label="PP_ARRAY"
    )

    if not preprocess_job_id:
        print("  No preprocessing job submitted, stopping")
        return
    else:
        print(f"    Submitted array {preprocess_job_id}: {len(all_nd2)} preprocessing tasks")

    # =========================================================================
    # STEP 2: Illumination Correction (parallelized, one job per file)
    # =========================================================================
    print(f"\n  [Step 2] Submitting illumination correction jobs ({len(all_nd2)} files, parallelized)...")

    # Create array script for illumination correction too
    illum_filelist_path = os.path.join(logs_dir, 'illum_filelist.txt')
    with open(illum_filelist_path, 'w') as f:
        for nd2 in all_nd2:
            basename = os.path.splitext(os.path.basename(nd2))[0]
            er_npy = os.path.join(preprocessed_dir, f"{basename}_er.npy")
            f.write(er_npy + '\n')

    illum_array_script = os.path.join(logs_dir, 'illum_array.sh')
    with open(illum_array_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --job-name=illum_{exp_name[:10]}\n")
        f.write(f"#SBATCH --output={logs_dir}/illum_%A_%a.log\n")
        f.write(_sbatch_header())
        f.write(f"#SBATCH --mem=8G\n")
        f.write(f"#SBATCH --cpus-per-task=1\n")
        f.write(f"#SBATCH --array=0-{len(all_nd2)-1}\n")
        f.write(f"\n")
        f.write(f"FILELIST={illum_filelist_path}\n")
        f.write(f"ER_PATH=$(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" $FILELIST)\n")
        f.write(f"\n")
        f.write(f"{PYTHON_PATH} {ER_WORKFLOW_DIR}/correct_illumination.py \\\n")
        f.write(f"    --er-file \"$ER_PATH\" --output {corrected_er_dir} \\\n")
        f.write(f"    --method {ILLUMINATION_METHOD}\n")

    illum_job_id = sbatch(
        ["sbatch", "--parsable",
         f"--dependency=afterok:{preprocess_job_id}",
         illum_array_script],
        dry_run=dry_run,
        dry_label="ILLUM_ARRAY"
    )

    if illum_job_id:
        print(f"    Submitted array {illum_job_id}: {len(all_nd2)} illumination correction tasks")

    # =========================================================================
    # STEP 3: Filewise tracking + ER classification + saSPT (per file)
    # Depends on BOTH preprocessing AND illumination correction
    # =========================================================================
    print(f"\n  [Step 3] Submitting filewise ER+tracking jobs...")

    filewise_job_ids = []
    conditionwise_job_ids = []

    # Filewise needs both preprocess and illum to be done
    fw_dependency = f"{preprocess_job_id}"
    if illum_job_id:
        fw_dependency = f"{preprocess_job_id}:{illum_job_id}"

    # Write filewise array script — each task processes one file
    # Needs both the SPT tif and the corrected ER npy
    fw_filelist = os.path.join(logs_dir, 'filewise_er_filelist.txt')
    with open(fw_filelist, 'w') as f:
        for nd2 in all_nd2:
            basename = os.path.splitext(os.path.basename(nd2))[0]
            spt_tif = os.path.join(preprocessed_dir, f"{basename}_spt.tif")
            er_npy = os.path.join(corrected_er_dir, f"{basename}_er.npy")
            # One line per file: spt_path<TAB>er_path<TAB>condition
            cond_name = None
            for c, files in conditions.items():
                if nd2 in files:
                    cond_name = c
                    break
            f.write(f"{spt_tif}\t{er_npy}\t{cond_name}\n")

    fw_array_script = os.path.join(logs_dir, 'filewise_er_array.sh')
    with open(fw_array_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --job-name=fwER_{exp_name[:10]}\n")
        f.write(f"#SBATCH --output={logs_dir}/FILEWISE_ER_%A_%a.log\n")
        f.write(_sbatch_header())
        f.write(f"#SBATCH --mem=16G\n")
        f.write(f"#SBATCH --cpus-per-task=1\n")
        f.write(f"#SBATCH --array=0-{len(all_nd2)-1}\n")
        f.write(f"\n")
        f.write(f"FILELIST={fw_filelist}\n")
        f.write(f"LINE=$(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" $FILELIST)\n")
        f.write(f"SPT_PATH=$(echo \"$LINE\" | cut -f1)\n")
        f.write(f"ER_PATH=$(echo \"$LINE\" | cut -f2)\n")
        f.write(f"CONDITION=$(echo \"$LINE\" | cut -f3)\n")
        f.write(f"\n")
        f.write(f"{PYTHON_PATH} {ER_WORKFLOW_DIR}/run_filewise_er.py \\\n")
        f.write(f"    -i \"$SPT_PATH\" \\\n")
        f.write(f"    -s {pkl_path} \\\n")
        f.write(f"    -c \"$CONDITION\" \\\n")
        f.write(f"    --er-image \"$ER_PATH\" \\\n")
        f.write(f"    --pct-floor {PCT_FLOOR} \\\n")
        f.write(f"    --require-broad {str(REQUIRE_BROAD).lower()}\n")

    fw_array_job_id = sbatch(
        ["sbatch", "--parsable",
         f"--dependency=afterok:{fw_dependency}",
         fw_array_script],
        dry_run=dry_run,
        dry_label="FW_ER_ARRAY"
    )

    if fw_array_job_id:
        print(f"    Submitted array {fw_array_job_id}: {len(all_nd2)} filewise ER+tracking tasks")
        filewise_job_ids.append(fw_array_job_id)

    # Conditionwise aggregation (one per condition, depends on filewise array)
    for cond in conditions.keys():
        cw_job_id = sbatch(
            ["sbatch", "--parsable",
             f"--job-name=cw_{cond[:12]}",
             f"--output={logs_dir}/CONDITIONWISE_{cond}_%j.log",
             *_sbatch_args(),
             "--mem=10G",
             "--cpus-per-task=4",
             f"--dependency=afterok:{fw_array_job_id}" if fw_array_job_id else "",
             f"--export=PKL_PATH={pkl_path},CONDITION={cond},CODE_DIRECTORY={CODE_DIRECTORY}",
             os.path.join(CODE_DIRECTORY, "call_conditionwise_operations.slurm")],
            dry_run=dry_run,
            dry_label=f"CW_{cond}"
        )

        if cw_job_id:
            print(f"    Submitted {cw_job_id}: conditionwise {cond}")
            conditionwise_job_ids.append(cw_job_id)

    # =========================================================================
    # Summary
    # =========================================================================
    print(f"\n  All jobs submitted for {exp_name}")
    print(f"  Monitor with: squeue -u $USER")
    print(f"  Results will be in: {settings['io']['analysis_directory']}")
    print(f"\n  Output directories:")
    print(f"    Preprocessed ER/SPT: {preprocessed_dir}")
    print(f"    Corrected ER:        {corrected_er_dir}")
    print(f"    Trajectories:        {traj_dir}")
    print(f"    MLE:                 {mle_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Submit ER+SPT dual-channel analysis workflow to SLURM"
    )
    parser.add_argument("directories", nargs="*", help="Experiment directories")
    parser.add_argument("--from-file", "-f", help="Read directories from file")
    parser.add_argument("--dry-run", "-n", action="store_true",
                        help="Preview without submitting jobs")

    args = parser.parse_args()

    directories = list(args.directories)
    if args.from_file:
        with open(args.from_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    directories.append(line)

    if not directories:
        parser.print_help()
        sys.exit(1)

    print("ER + SPT Dual-Channel Workflow")
    print(f"Code: {CODE_DIRECTORY}")
    print(f"ER scripts: {ER_WORKFLOW_DIR}")
    print(f"\nWorkflow steps:")
    print(f"  1. Preprocess: Extract ER + SPT from ND2 (job array)")
    print(f"  2. Illumination correction: {ILLUMINATION_METHOD} (parallelized)")
    print(f"  3. Track: quot/saspt on SPT")

    if args.dry_run:
        print("\n*** DRY RUN - No jobs will be submitted ***")

    for exp_dir in directories:
        if not os.path.isdir(exp_dir):
            print(f"\nWarning: Directory not found: {exp_dir}")
            continue
        submit_experiment(exp_dir, dry_run=args.dry_run)

    print(f"\n{'='*60}")
    print("Done!")


if __name__ == "__main__":
    main()
