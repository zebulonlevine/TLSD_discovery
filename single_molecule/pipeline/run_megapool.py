## run_megapool.py
##
## Entry point for joint-saSPT analysis across a user-defined megapool
## of conditions and/or combinations.
##
## Pipeline:
##   1. Look up the megapool in megapools.yaml; collect its members and
##      any optional cutoff overrides.
##   2. Skip-if-exists: if the two pooled saSPT outputs already exist
##      and --force is not set, jump straight to partitioning + downstream.
##   3. Pool member trajectory CSVs (with offset trajectory IDs) into a
##      single CSV. Write a sidecar source map (trajectory_id → member).
##   4. Run executeSASPT on the pooled CSV with redirected I/O →
##      posterior_megapools/, MLE_megapools/.
##      (rolling-window analysis is intentionally skipped.)
##   5. Partition the resulting MLE rows back to source members via
##      the source map, joining on the input trajectory ID embedded in
##      orig_trajectory.
##   6. For each member: write a per-member MLE CSV and a duplicate
##      copy of the shared joint posterior, then run
##      conditionwise_core.process_condition on that member with a
##      _megapool_context flag so plots go to plots/plotting_pkls/
##      megapool_{NAME}/{member}.pkl.
##
## SLURM invocation:
##   python run_megapool.py --megapool NAME --settingsPKL PKL [--force]


import argparse
import os
import pickle
import sys


# Put this script's own directory on sys.path so sibling modules resolve
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


print("")
print('parsing megapool arguments')
parser = argparse.ArgumentParser(
    description="Joint-saSPT analysis across a megapool of members.")
parser.add_argument('-s', '--settingsPKL', required=True,
                    help='Settings PKL path')
parser.add_argument('-m', '--megapool', required=True,
                    help='Megapool name (must appear in megapools.yaml)')
parser.add_argument('--force', action='store_true',
                    help='Re-run pooling and saSPT even if outputs exist')

args = parser.parse_args()

print(f"reading settings from: {args.settingsPKL}")
print(f"Working on megapool: {args.megapool}")
print(f"force re-pool: {args.force}")

with open(args.settingsPKL, "rb") as f:
    settings = pickle.load(f)


import trackingCodeFunctions
import megapool_functions
import combination_functions   # to detect which members are combos
import conditionwise_core


# ============================================================
# 1. Resolve the megapool from megapools.yaml
# ============================================================
megapools_path = os.path.join(settings['io']['home_directory'],
                              'megapools.yaml')
megapool_name, megapool_cfg = megapool_functions.find_megapool(
    args.megapool, megapools_path)

if megapool_name is None:
    raise SystemExit(
        f"Megapool '{args.megapool}' not found in {megapools_path}. "
        f"Aborting.")

members = megapool_cfg['members']
max_per_member = megapool_cfg.get('max_trajectories_per_member', None)
max_total = megapool_cfg.get('max_trajectories_total', None)
balance_to_min_trajectories = bool(megapool_cfg.get('balance_to_min_trajectories', False))
random_seed = int(megapool_cfg.get('random_seed', 42))
print(f"  members: {members}")
print(f"  max_trajectories_per_member: {max_per_member}")
print(f"  max_trajectories_total: {max_total}")
print(f"  balance_to_min_trajectories: {balance_to_min_trajectories}")


# ============================================================
# 2. Load combinations.yaml so we know which members are combos
# ============================================================
combinations_path = os.path.join(settings['io']['home_directory'],
                                  'combinations.yaml')
combinations_cfg = combination_functions.load_combinations(combinations_path)


# ============================================================
# 3. Redirect I/O
# ============================================================
redirected = megapool_functions.redirect_settings_for_megapool(settings)


# ============================================================
# 4. Pool + run saSPT (skip-if-exists)
# ============================================================
outputs_exist = megapool_functions.megapool_outputs_exist(
    megapool_name, redirected)
source_map_csv = os.path.join(
    redirected['io']['traj_megapools_directory'],
    f"{megapool_name}_source_map.csv")
pooled_csv = os.path.join(
    redirected['io']['traj_megapools_directory'],
    f"{megapool_name}_pooled_traj.csv")

if outputs_exist and os.path.exists(source_map_csv) and not args.force:
    print(f"\n  Pooled saSPT outputs already exist for megapool "
          f"'{megapool_name}'. Skipping pooling/saSPT. Use --force to re-run.")
else:
    if outputs_exist and args.force:
        print(f"\n  --force set; re-running pooling/saSPT for "
              f"'{megapool_name}'.")

    print("\n  Pooling member trajectory CSVs...")
    pooled_csv, source_map_csv = megapool_functions.pool_megapool_trajectories(
        megapool_name, members, redirected,
        combinations_cfg=combinations_cfg,
        max_per_member=max_per_member,
        max_total=max_total,
        balance_to_min_trajectories=balance_to_min_trajectories,
        random_seed=random_seed)

    print("\n  Running joint saSPT on pooled trajectories...")
    print("  (rolling-window analysis is skipped for megapools)")
    trackingCodeFunctions.executeSASPT(pooled_csv, redirected)


# ============================================================
# 5. Partition MLE by source member
# ============================================================
joint_mle_csv = os.path.join(
    redirected['io']['MLE_megapools_directory'],
    f"{megapool_name}_pooled_MLE.csv")
joint_posterior_csv = os.path.join(
    redirected['io']['post_megapools_directory'],
    f"{megapool_name}_pooled_posterior.csv")

print(f"\n  Partitioning joint MLE by source member...")
mle_with_source, stats = megapool_functions.partition_mle_by_source(
    joint_mle_csv, source_map_csv)
print(f"    {stats['n_mle_rows']:,} MLE rows total")
print(f"    {stats['n_traj_id_extracted']:,} input trajectory IDs extracted")
print(f"    {stats['n_matched_to_source']:,} matched to source map")
print(f"    {stats['n_unmatched']:,} unmatched (saSPT subsampling artifact)")
print(f"    per-member subtrack counts:")
for m, c in sorted(stats['per_member_counts'].items()):
    print(f"      {m}: {c:,}")

joint_posterior = megapool_functions.partition_posterior_by_source(
    joint_posterior_csv, mle_with_source)


# ============================================================
# 6. Write per-member files and run process_condition per member
# ============================================================
print(f"\n  Running downstream processing per member...")

# Override plot_directory to a megapool subdirectory so per-member
# plotting PKLs land in plots/plotting_pkls/megapool_{NAME}/{member}.pkl
# rather than colliding with the per-condition or per-combo PKLs.
member_settings_base = redirected
megapool_plot_subdir = os.path.join(
    member_settings_base['io']['plot_directory'],
    f"megapool_{megapool_name}")
os.makedirs(megapool_plot_subdir, exist_ok=True)

# We'll keep the parent plot_directory for things like the GMM
# diagnostics file location, but the plotting_pkls write target gets
# redirected per-member below.

import copy

for member in sorted(set(members)):
    if member not in stats['per_member_counts']:
        print(f"    member '{member}': no subtracks; skipping")
        continue

    # Write member-tagged copies of MLE and posterior so
    # process_condition's existing glob '{cond}_pooled_*.csv' finds
    # them under a unique name.
    cond_name = f"{megapool_name}_{member}"
    mle_path, n_rows = megapool_functions.write_member_mle_csv(
        mle_with_source, megapool_name, member, member_settings_base)
    post_path = megapool_functions.write_member_posterior_csv(
        joint_posterior, megapool_name, member, member_settings_base,
        mle_df_with_source=mle_with_source)
    print(f"    member '{member}': {n_rows:,} subtracks → "
          f"{os.path.basename(mle_path)}, {os.path.basename(post_path)}")

    # Build per-member settings:
    #   - plotting_pkls write goes to plots/plotting_pkls/megapool_{name}/
    #     We achieve this by redirecting plot_directory inside settings
    #     to a fresh subdirectory whose plotting_pkls/ subdir we want.
    member_settings = copy.deepcopy(member_settings_base)
    member_plot_dir = os.path.join(
        member_settings_base['io']['plot_directory'],
        f"megapool_{megapool_name}")
    member_settings['io']['plot_directory'] = member_plot_dir
    os.makedirs(member_plot_dir, exist_ok=True)

    # Inject the megapool context so process_condition knows it's a
    # megapool partition (uses its own cutoff cache key, knows there's
    # no rolling-window data).
    member_settings['_megapool_context'] = {
        'megapool_name': megapool_name,
        'member': member,
        'member_list': members,
    }
    # For cutoff resolution, this looks like a combination_context too —
    # it's a pooled comparison-unit fit on shared data. Reuse that
    # branch in process_condition by also setting _combination_context.
    member_settings['_combination_context'] = {
        'combo_name': cond_name,
        'member_conditions': [member],
    }

    conditionwise_core.process_condition(cond_name, member_settings)

print(f"\nMegapool analysis complete for '{megapool_name}'")
print(f"Per-member plotting PKLs are at:")
print(f"  {megapool_plot_subdir}/plotting_pkls/")
