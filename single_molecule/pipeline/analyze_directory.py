#### take a filename and run quot on it
from quot.core import track_file
import init_settings
import trackingCodeFunctions
import plotFunctions

import os
import subprocess

import pickle

import pims
from collections import defaultdict
import re
import time
import yaml
import shutil

import argparse

import pandas as pd
import numpy as np

force = False

parser = argparse.ArgumentParser(description="Run full SASPT analysis")
parser.add_argument("--input_directory", help="Path to input directory")
parser.add_argument("--force", help="Force recalculation of spots and trajectories")

args = parser.parse_args()
directory_to_analyze = args.input_directory

print(f"Working from directory : {directory_to_analyze}")

if args.force != None:
    force = True
print(f"  force set to : {force}")



## identify conditions
    ## for each condition, generate would-be output filenames   x
    ## check if the output filenames exist. If they do NOT exist  x
        ## run file-wise SLURM job for each missing file, using a job array
            ## filewise script will run...
                ## QUOT
                ## rolling_windows
                ## SASPT
    ## continually check to see if any of the conditions' job arrays have finished. As they finish...
        ## run condition-wise SLURM job for newly finsihed condition
            ## condition-wise script will run...
                ## QC histograms / remove outliers
                ## plot histograms
                ## save df of histogram data
                ## plot barGraphs
                ## save df if barGraph data
                ## aggregate rolling_windows, run HMM
                ## print HMM output, plot empirical survival 1 plot per pop from HMM posterior
                ## save PLK of empirical survival curves from HMM posterior
                ## option to run the published HMM
                ## save example movies
## Aggregate plots
    ## aggregate df of histogram data, barGraph data, survival curves
    ## plot one file per data-stream (1 barGraph, 1 histogram, 1 survival)



## should only require 2 SLURM scripts... one file-wise and one condition-wise


########################
# setup block
####
## initialize settings and dump them to a pkl for use by sub-scripts
## settings file will live in the analysis directory, and will not be updated after it is created

settings = init_settings.get_default_settings()

# Default code_directory to the directory holding analyze_directory.py
# itself. A settings_override.yaml may overwrite this if the user wants
# to explicitly point elsewhere.
settings['io']['code_directory'] = os.path.dirname(
    os.path.abspath(__file__))

override_file = 'settings_override.yaml'
override_path = os.path.join(directory_to_analyze , override_file)
if os.path.exists(override_path):
    with open(override_path) as f:
        override_cfg = yaml.load(f, Loader=yaml.FullLoader)
    trackingCodeFunctions.deep_update(settings, override_cfg)
settings = init_settings.update_default_settings_for_analysis(settings , directory_to_analyze)
####

####
### dump settings file and print settings to stdout
pkl_path = os.path.join(settings['io']['analysis_directory'] , "settings.pkl")
with open(pkl_path, "wb") as f:
    pickle.dump(settings, f)

txt_path = os.path.join(settings['io']['analysis_directory'] , "settings.txt")
print('')
#print_nested_dict(settings)
trackingCodeFunctions.print_nested_dict(settings , file = txt_path)
print('')
####

####
data_dir = os.path.join(settings['io']['home_directory'] , settings['io']['data_filepath'])
#grouped_files = group_files_by_metadata(data_dir) ## NOTE FILE NAMING CONVENTIONS ARE NOT VERY FLEXIBLE

    
    
#### start to do the work
#### by default, this will *not* recalculate the spots and trajectories, but if the user provides --force then it will

if args.force != None:
    force = True

nd2_files = trackingCodeFunctions.identify_missing_filewise(settings , force = force)

filewise_script = os.path.join(settings['io']['code_directory'], 'call_filewise_operations.slurm')
conditionwise_script = os.path.join(settings['io']['code_directory'], 'call_conditionwise_operations.slurm')

# Dictionary to hold submitted job array IDs per condition
job_ids_by_condition = {}
for cond, fileList in nd2_files.items():
    fileList_string = ":".join(os.path.join(settings['io']['data_directory'], f) for f in fileList)
    
    # Submit the job array
    result = subprocess.run([
        "sbatch",
        f"--array=0-{len(fileList)-1}",
        '--parsable',  # Returns only the job ID
        '--export', f"ITEM_STRING={fileList_string},PKL_PATH={pkl_path},CONDITION={cond},CODE_DIRECTORY={settings['io']['code_directory']}",
        filewise_script
    ], capture_output=True, text=True)

    job_id = result.stdout.strip()  # SLURM job ID of the array
    job_ids_by_condition[cond] = job_id
    print(f"Submitted array for condition {cond}, Job ID: {job_id}")

print('')

## if there were no files previously run, go head and pull the conditions with empty dependencies associated
if not job_ids_by_condition:
    job_ids_by_condition = trackingCodeFunctions.identify_missing_filewise(settings , force = True)
    for key in job_ids_by_condition:
        job_ids_by_condition[key] = 0

# Now submit the "next step" jobs that depend on each array finishing
for cond, array_job_id in job_ids_by_condition.items():
    # The dependency string for SLURM: all jobs in array must finish
    dependency_string = f"afterok:{array_job_id}"
    
    if array_job_id == 0:
        job_string = [
            "sbatch",
            "--export", f"CONDITION={cond},PKL_PATH={pkl_path},CODE_DIRECTORY={settings['io']['code_directory']}",
            conditionwise_script
        ]
    else:
        job_string = [
            "sbatch",
            f"--dependency={dependency_string}",
            "--export", f"CONDITION={cond},PKL_PATH={pkl_path},CODE_DIRECTORY={settings['io']['code_directory']}",
            conditionwise_script
        ]
    
    result_next = subprocess.run(job_string, capture_output=True, text=True)

    if result_next.returncode != 0:
        print(f"Error submitting job for {cond}:")
        print(result_next.stderr)
    else:
        print(f"Submitted next step for {cond}, dependent on array {array_job_id}")
        print(result_next.stdout)


########################
# combination submission
####
## If combinations.yaml exists at the experiment root, submit one
## combinationwise SLURM job per combination, dependent on the file-wise
## arrays of its member conditions. Per-condition outputs continue
## independently; combinations run in parallel as soon as their member
## arrays finish.
import combination_functions

combinations_yaml = os.path.join(directory_to_analyze, 'combinations.yaml')
combinations = combination_functions.load_combinations(combinations_yaml)

combination_job_ids = {}  # combo_name -> SLURM job ID (megapool deps)

if combinations:
    combinationwise_script = os.path.join(
        settings['io']['code_directory'],
        'call_combinationwise_operations.slurm')

    print('')
    print(f"Found {len(combinations)} combination(s) in combinations.yaml.")

    combination_job_ids = {}  # combo_name -> SLURM job ID, for megapool deps
    force_env = "1" if force else "0"

    for combo_name, combo_cfg in combinations.items():
        members = combo_cfg['conditions']

        # Collect job IDs of all member conditions whose file-wise arrays
        # were submitted. Conditions not in job_ids_by_condition either
        # had no missing files (already done) or aren't in this dataset
        # — both fine; we just don't add them to the dependency.
        dep_ids = []
        unknown_members = []
        for m in members:
            if m in job_ids_by_condition and job_ids_by_condition[m] != 0:
                dep_ids.append(str(job_ids_by_condition[m]))
            elif m in job_ids_by_condition:
                # member is known but had no jobs (already complete)
                pass
            else:
                unknown_members.append(m)

        if unknown_members:
            print(f"  WARNING: combination '{combo_name}' references "
                  f"unknown conditions: {unknown_members}. They will be "
                  f"silently ignored unless their tracking outputs exist.")

        export_str = (
            f"COMBINATION={combo_name},"
            f"PKL_PATH={pkl_path},"
            f"CODE_DIRECTORY={settings['io']['code_directory']},"
            f"FORCE={force_env}"
        )

        if dep_ids:
            dependency_string = "afterok:" + ":".join(dep_ids)
            cmd = [
                "sbatch",
                f"--dependency={dependency_string}",
                "--export", export_str,
                combinationwise_script
            ]
        else:
            # No new file-wise jobs were submitted for any member condition,
            # so this combo can run immediately.
            cmd = [
                "sbatch",
                "--export", export_str,
                combinationwise_script
            ]

        result_combo = subprocess.run(cmd, capture_output=True, text=True)
        if result_combo.returncode != 0:
            print(f"Error submitting combination {combo_name}:")
            print(result_combo.stderr)
        else:
            dep_summary = (
                f"dependent on {len(dep_ids)} member array(s)"
                if dep_ids else "no dependencies (members already complete)")
            print(f"Submitted combination '{combo_name}' "
                  f"({dep_summary}): {result_combo.stdout.strip()}")
            # Capture the combo job ID for later megapool dependencies
            try:
                # sbatch stdout looks like "Submitted batch job 12345"
                combo_job_id = int(result_combo.stdout.strip().split()[-1])
                combination_job_ids[combo_name] = combo_job_id
            except Exception:
                pass


########################
# megapool submission
####
## If megapools.yaml exists, submit one megapool SLURM job per entry,
## with afterok dependencies on:
##   - the filewise array of each raw-condition member, AND
##   - the combinationwise job of each combination-nickname member.
## A megapool can only run after every contributing pool is materialized.
import megapool_functions

megapools_yaml = os.path.join(directory_to_analyze, 'megapools.yaml')
megapools = megapool_functions.load_megapools(megapools_yaml)

if megapools:
    megapool_script = os.path.join(
        settings['io']['code_directory'],
        'call_megapool_operations.slurm')

    print('')
    print(f"Found {len(megapools)} megapool(s) in megapools.yaml.")

    force_env = "1" if force else "0"

    for megapool_name, megapool_cfg in megapools.items():
        members = megapool_cfg['members']

        dep_ids = []
        unresolved = []
        for m in members:
            # A member can be a combination nickname or a raw condition
            if combinations and m in combinations:
                # Wait for its combination job
                if m in combination_job_ids:
                    dep_ids.append(str(combination_job_ids[m]))
                # else: combination already complete, no dep needed
            elif m in job_ids_by_condition and job_ids_by_condition[m] != 0:
                # Raw condition with a freshly-submitted filewise array
                dep_ids.append(str(job_ids_by_condition[m]))
            elif m in job_ids_by_condition:
                # Raw condition with no pending jobs (already complete)
                pass
            else:
                unresolved.append(m)

        if unresolved:
            print(f"  WARNING: megapool '{megapool_name}' references "
                  f"unknown members: {unresolved}. They will be silently "
                  f"ignored unless their tracking/combination outputs exist.")

        export_str = (
            f"MEGAPOOL={megapool_name},"
            f"PKL_PATH={pkl_path},"
            f"CODE_DIRECTORY={settings['io']['code_directory']},"
            f"FORCE={force_env}"
        )

        if dep_ids:
            dependency_string = "afterok:" + ":".join(dep_ids)
            cmd = [
                "sbatch",
                f"--dependency={dependency_string}",
                "--export", export_str,
                megapool_script
            ]
        else:
            cmd = [
                "sbatch",
                "--export", export_str,
                megapool_script
            ]

        result_mp = subprocess.run(cmd, capture_output=True, text=True)
        if result_mp.returncode != 0:
            print(f"Error submitting megapool {megapool_name}:")
            print(result_mp.stderr)
        else:
            dep_summary = (
                f"dependent on {len(dep_ids)} upstream job(s)"
                if dep_ids else "no dependencies")
            print(f"Submitted megapool '{megapool_name}' "
                  f"({dep_summary}): {result_mp.stdout.strip()}")


