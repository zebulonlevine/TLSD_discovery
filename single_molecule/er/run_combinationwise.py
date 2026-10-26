## run_combinationwise.py
##
## Entry point for pooled-saSPT analysis on a user-defined combination of
## conditions. Combinations are declared in `combinations.yaml` at the
## experiment root.
##
## Pipeline:
##   1. Look up the combination in combinations.yaml; collect its member
##      conditions and any optional cutoff overrides.
##   2. Skip-if-exists: if the three pooled saSPT outputs already exist
##      and --force is not set, jump straight to downstream processing.
##   3. Pool member conditions' trajectory CSVs (with traj-ID offsetting)
##      → traj_combinations/{combo}_pooled_traj.csv.
##   4. Run executeSASPT on the pooled CSV with redirected I/O →
##      posterior_combinations/, MLE_combinations/.
##   5. Run analyze_rolling_windows on the pooled CSV →
##      rolling_windows_combinations/.
##   6. Run conditionwise_core.process_condition with the combo name as
##      the "condition", honoring per-combo cutoff overrides via the
##      _combination_context settings injection.
##
## SLURM invocation (see call_combinationwise_operations.slurm):
##   python run_combinationwise.py \
##       --combination NAME --settingsPKL PKL [--force]


import argparse
import os
import pickle
import sys


# Put this script's own directory on sys.path so all sibling modules
# (combination_functions, conditionwise_core, cutoff_functions, ...) are
# findable. Done before any project-module imports. Independent of
# whatever 'code_directory' the settings PKL has.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# pipeline modules (trackingCodeFunctions, conditionwise_core) come after this
# directory, so the ER copy of combination_functions shadows the pipeline one,
# exactly as when this ran with PYTHONPATH=<pipeline> (er/combo_array.sh).
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pipeline')))


print("")
print('parsing combinationwise arguments')
parser = argparse.ArgumentParser(
    description="Pooled-saSPT analysis on a combination of conditions.")
parser.add_argument('-s', '--settingsPKL', required=True,
                    help='Settings PKL path')
parser.add_argument('-c', '--combination', required=True,
                    help='Combination name (must appear in combinations.yaml)')
parser.add_argument('--force', action='store_true',
                    help='Re-run pooling and saSPT even if outputs exist')

args = parser.parse_args()

print(f"reading settings from: {args.settingsPKL}")
print(f"Working on combination: {args.combination}")
print(f"force re-pool: {args.force}")

with open(args.settingsPKL, "rb") as f:
    settings = pickle.load(f)

import trackingCodeFunctions
import combination_functions
import conditionwise_core


# ============================================================
# 1. Resolve the combination from combinations.yaml
# ============================================================
combos_path = os.path.join(settings['io']['home_directory'], 'combinations.yaml')
combo_name, combo_cfg = combination_functions.find_combination(
    args.combination, combos_path)

if combo_name is None:
    raise SystemExit(
        f"Combination '{args.combination}' not found in {combos_path}. "
        f"Aborting.")

member_conditions = combo_cfg['conditions']
print(f"  member conditions: {member_conditions}")


# ============================================================
# 2. Apply combo-specific cutoff overrides + redirect I/O dirs
# ============================================================
combination_functions.apply_combo_cutoff_overrides(settings, combo_cfg)
redirected = combination_functions.redirect_settings_for_combo(settings)


# ============================================================
# 3. Skip-if-exists gate (saSPT pooling step only)
# ============================================================
outputs_exist = combination_functions.combination_outputs_exist(
    combo_name, redirected)

if outputs_exist and not args.force:
    print(f"\n  Pooled saSPT outputs already exist for '{combo_name}'. "
          f"Skipping pooling/saSPT/rolling. Use --force to re-run.")
else:
    if outputs_exist and args.force:
        print(f"\n  --force set; re-running pooling/saSPT/rolling for "
              f"'{combo_name}' even though outputs exist.")

    # 4a. Pool trajectories. zone is None for legacy combos (pulls
    # from the canonical trajectories/ dir); for zone-tagged combos it
    # points the pooler at trajectories_{zone}/ from run_filewise_er.py.
    print("\n  Pooling member-condition trajectory CSVs...")
    pooled_csv = combination_functions.pool_combination_trajectories(
        combo_name, member_conditions, redirected,
        zone=combo_cfg.get("zone"))

    # 4b. saSPT on pooled CSV
    print("\n  Running saSPT on pooled trajectories...")
    trackingCodeFunctions.executeSASPT(pooled_csv, redirected)

    # 4c. Rolling-window saSPT
    print("\n  Running rolling-window saSPT...")
    trackingCodeFunctions.analyze_rolling_windows(pooled_csv, redirected)


# ============================================================
# 4. Inject combination context and run downstream
# ============================================================
redirected['_combination_context'] = {
    'combo_name': combo_name,
    'member_conditions': member_conditions,
}

conditionwise_core.process_condition(combo_name, redirected)

print(f"\nCombinationwise analysis complete for '{combo_name}'")
