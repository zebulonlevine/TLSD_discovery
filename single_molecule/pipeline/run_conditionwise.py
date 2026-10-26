## run_conditionwise.py
##
## Per-condition entry point. Thin CLI wrapper around
## conditionwise_core.process_condition(cond, settings).
##
## Called by SLURM as:
##   python run_conditionwise.py --condition COND --settingsPKL PKL


import argparse
import os
import pickle
import sys


# Put this script's own directory on sys.path so all sibling modules
# (conditionwise_core, plotFunctions, ...) are findable.
# Done before any project-module imports so the rest of the file can rely
# on it. Independent of whatever 'code_directory' the settings PKL has.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


print("")
print('parsing Conditionwise arguments')
parser = argparse.ArgumentParser(description="Per-condition downstream processing.")
parser.add_argument('-s', '--settingsPKL', required=True,
                    help='Input settings Pickle file path')
parser.add_argument('-c', '--condition', required=True,
                    help='Condition name to process')

args = parser.parse_args()

print(f"reading settings from: {args.settingsPKL}")
print(f"Working on condition: {args.condition}")

with open(args.settingsPKL, "rb") as f:
    settings = pickle.load(f)

import conditionwise_core

conditionwise_core.process_condition(args.condition, settings)
