## combination_functions.py
##
## Helpers for pooled-saSPT analysis across user-defined combinations of
## conditions, declared in `combinations.yaml` at the experiment root.
##
## A combination is a "nickname" — a virtual condition name — that maps to
## a list of real condition names. Trajectories from the listed conditions
## are concatenated (with trajectory-ID offsetting to avoid collisions),
## written to a pooled CSV under traj_combinations/, and fed to saSPT and
## analyze_rolling_windows. Outputs land in posterior_combinations/,
## MLE_combinations/, and rolling_windows_combinations/. Downstream
## processing then runs through conditionwise_core.process_condition with
## the combination name in place of a real condition.


import copy
import glob
import os

import pandas as pd
import yaml


def load_combinations(yaml_path):
    """
    Load combination definitions from a YAML file.

    Expected format:
        combinations:
          all_drug:
            conditions: [drug_rep1, drug_rep2, drug_rep3]
            labels: [immobile, slow, fast]      # optional
            n_components: 3                     # optional
            initial_centers: [-1.3, -0.3, 0.5]  # optional

    Parameters
    ----------
    yaml_path : str
        Path to combinations.yaml.

    Returns
    -------
    dict : {combo_name: {'conditions': [...], ...}}
        Empty dict if the file does not exist.
    """
    if not os.path.exists(yaml_path):
        return {}

    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f) or {}

    combos = config.get('combinations', {})
    # Reject combos with empty/whitespace conditions or path-unsafe
    # nicknames. Combination names become file basenames and plot PKL
    # keys, so anything that would corrupt a path needs to fail fast.
    forbidden = set('/\\:*?"<>|\t\n\r ')
    for name, cfg in combos.items():
        if cfg is None:
            raise ValueError(
                f"Combination '{name}' has no body. Either remove the "
                f"entry or add a 'conditions:' list under it.")
        if 'conditions' not in cfg:
            raise ValueError(f"Combination '{name}' missing 'conditions' key")
        if not isinstance(cfg['conditions'], list) or len(cfg['conditions']) == 0:
            raise ValueError(
                f"Combination '{name}' must have a non-empty list of conditions")
        bad_chars = forbidden & set(name)
        if bad_chars:
            raise ValueError(
                f"Combination name '{name}' contains characters that are "
                f"not safe for filenames or paths: {sorted(bad_chars)}. "
                f"Replace with underscores (e.g., 'sgUBA5_sgR3HDM4').")
        # Optional ER zone filter. Validated at load time so a typo in
        # the YAML fails fast (before SLURM submission) rather than at
        # combinationwise time. None / missing / 'all' all mean "no
        # filter, pull from canonical traj_directory".
        if 'zone' in cfg and cfg['zone'] is not None:
            if not isinstance(cfg['zone'], str):
                raise ValueError(
                    f"Combination '{name}': 'zone' must be a string or "
                    f"omitted; got {type(cfg['zone']).__name__}.")
    return combos


def find_combination(combo_name, yaml_path):
    """Return (combo_name, combo_cfg) or (None, None)."""
    combos = load_combinations(yaml_path)
    if combo_name in combos:
        return combo_name, combos[combo_name]
    return None, None


def pool_combination_trajectories(combo_name, member_conditions, settings,
                                  nd2_by_cond=None, zone=None):
    """
    Concatenate per-FOV trajectory CSVs from a list of member conditions,
    offsetting trajectory IDs so they remain globally unique, and write
    the pooled CSV to traj_combinations/{combo_name}_pooled_traj.csv.

    Parameters
    ----------
    combo_name : str
        Combination nickname.
    member_conditions : list of str
        Real condition names whose tracking CSVs should be pooled.
    settings : dict
        Pipeline settings (must have I/O directories already initialized).
    nd2_by_cond : dict or None
        Optional preloaded {condition: [nd2_filenames]} mapping. If None,
        it is computed via trackingCodeFunctions.group_files_by_metadata.
    zone : str or None
        Optional ER zone filter. When set, the pooler reads per-FOV traj
        CSVs from <analysis_dir>/trajectories_{zone}/ instead of the
        canonical settings['io']['traj_directory']. The per-FOV CSVs in
        zone subdirectories are produced by run_filewise_er.py during
        the ER classification step and have re-indexed trajectory IDs
        (still 0-based within each FOV, so traj-ID offsetting works the
        same way). When None (default), behavior is unchanged from older
        callers. Valid values are anything that exists as a directory
        named trajectories_{zone}/; commonly: 'on_structure',
        'on_ER_broad', 'off_ER'. Use zone='all' as a no-op alias that
        explicitly pulls from the canonical traj_directory — useful for
        making an apples-to-apples joint-fit 'all' combo alongside the
        zone-filtered combos.

    Returns
    -------
    str
        Path to the pooled trajectory CSV.

    Notes
    -----
    Per-FOV traj CSVs use trajectory IDs starting at 0 within each file,
    so naive concatenation would collide IDs across files. We add a
    cumulative offset equal to the previous max ID + 1 before each new
    file's data. The 'trajectory' column is the only one saSPT groups by;
    other columns (frame, x, y, I0, etc.) need no rewriting.
    """
    import trackingCodeFunctions

    if nd2_by_cond is None:
        data_dir = os.path.join(settings['io']['home_directory'],
                                settings['io']['data_filepath'])
        # Prefer .nd2; fall back to .tif if no nd2 files exist.
        nd2_by_cond = trackingCodeFunctions.group_files_by_metadata(
            data_dir, extension='.nd2')
        if not nd2_by_cond:
            nd2_by_cond = trackingCodeFunctions.group_files_by_metadata(
                data_dir, extension='.tif')

    # Source directory for per-FOV traj CSVs: canonical or zone-filtered.
    # 'all' is treated as an explicit no-op alias for the canonical dir,
    # so combinations.yaml can carry zone: all entries for symmetry with
    # the three ER zone combos without changing pool behavior.
    if zone is None or zone == 'all':
        traj_dir = settings['io']['traj_directory']
    else:
        traj_dir = os.path.join(
            settings['io']['analysis_directory'],
            f'trajectories_{zone}')
        if not os.path.isdir(traj_dir):
            raise RuntimeError(
                f"Combination '{combo_name}' requests zone '{zone}' but "
                f"the per-FOV zone directory {traj_dir} does not exist. "
                f"Run run_filewise_er.py first to produce zone subsets.")

    out_dir = settings['io']['traj_combinations_directory']
    os.makedirs(out_dir, exist_ok=True)

    pieces = []
    fov_map_rows = []   # collect (trajectory_id, source_FOV_basename)
    offset = 0
    files_used = 0
    files_missing = 0

    for cond in member_conditions:
        if cond not in nd2_by_cond:
            print(f"  WARNING: combination '{combo_name}' references "
                  f"condition '{cond}' which has no source files.")
            continue

        for src_file in nd2_by_cond[cond]:
            base = os.path.splitext(src_file)[0]
            traj_path = os.path.join(traj_dir, base + '_traj.csv')
            if not os.path.exists(traj_path):
                files_missing += 1
                print(f"    missing trajectory CSV: {traj_path}")
                continue

            df = pd.read_csv(traj_path)
            if 'trajectory' not in df.columns:
                print(f"    skipping {traj_path}: no 'trajectory' column")
                continue

            df['trajectory'] = df['trajectory'] + offset
            offset = int(df['trajectory'].max()) + 1
            pieces.append(df)
            files_used += 1

            # Record (offset trajectory ID → source FOV) for each unique
            # trajectory in this file. The FOV name is the basename of the
            # per-FOV traj CSV (matches what process_condition would have
            # stamped as file_basename had this been per-FOV input).
            fov_basename = os.path.basename(traj_path)
            for tid in df['trajectory'].unique():
                fov_map_rows.append({
                    'trajectory': int(tid),
                    'condition': cond,
                    'fov_basename': fov_basename,
                })

    if not pieces:
        raise RuntimeError(
            f"Combination '{combo_name}': no trajectory CSVs found for "
            f"member conditions {member_conditions}")

    pooled = pd.concat(pieces, ignore_index=True)
    out_path = os.path.join(out_dir, f"{combo_name}_pooled_traj.csv")
    pooled.to_csv(out_path, index=False)

    # Write the sidecar map that lets downstream tooling recover per-FOV
    # provenance from a pooled trajectory ID. Used by spt_compare to do
    # per-FOV analysis on combinations.
    fov_map_path = os.path.join(out_dir, f"{combo_name}_fov_map.csv")
    pd.DataFrame(fov_map_rows).to_csv(fov_map_path, index=False)

    print(f"  Pooled {files_used} FOV(s) from {len(member_conditions)} "
          f"condition(s)"
          + (f" [zone={zone}]" if zone is not None else "")
          + f" into {out_path}")
    print(f"  Wrote FOV map sidecar: {fov_map_path}")
    print(f"  Total trajectories: {pooled['trajectory'].nunique()}; "
          f"total detections: {len(pooled)}")
    if files_missing:
        print(f"  WARNING: {files_missing} expected trajectory CSV(s) were "
              f"missing.")

    return out_path


def redirect_settings_for_combo(settings):
    """
    Return a deep-copied settings dict with post/MLE/rolling-window I/O
    keys pointing to the *_combinations directories.

    This lets executeSASPT and analyze_rolling_windows write to the
    combinations namespace without per-call argument plumbing, and lets
    process_condition aggregate from the same directories.

    Also forces saspt.sample_size to a huge value (1e9) so saSPT does
    not internally subsample. This preserves the 1:1 mapping between
    input trajectory IDs and orig_trajectory values in the MLE output,
    which is required for any downstream that needs to reconstruct
    parent trajectories from subtracks (lifetime analysis, per-FOV
    stratification on combos via the fov_map sidecar, etc.).

    Pool-size control is exercised at the combination_outputs level via
    the per-FOV cap inherited from filewise tracking. If you need to
    cap combo pool sizes further, do it by trimming the input traj CSVs
    or via a future combo-level cap setting.

    Parameters
    ----------
    settings : dict

    Returns
    -------
    dict
        Deep copy with three I/O keys redirected and saspt.sample_size
        forced to 1e9.
    """
    new = copy.deepcopy(settings)
    new['io']['post_directory'] = settings['io']['post_combinations_directory']
    new['io']['MLE_directory'] = settings['io']['MLE_combinations_directory']
    new['io']['rolling_window_directory'] = settings['io']['rolling_window_combinations_directory']
    new['saspt'] = dict(new['saspt'])
    new['saspt']['sample_size'] = int(1e9)
    return new


def combination_outputs_exist(combo_name, settings):
    """
    Return True iff all three saSPT outputs (posterior, MLE, rolling) for
    the combination already exist. Used to skip re-running saSPT on
    re-execution unless --force is passed.

    Parameters
    ----------
    combo_name : str
    settings : dict
        Should be the *redirected* settings (post/MLE/rolling pointing to
        the combinations dirs); but works with the original settings too,
        as it explicitly checks the *_combinations_directory keys.

    Returns
    -------
    bool
    """
    expected = [
        os.path.join(settings['io']['post_combinations_directory'],
                     f"{combo_name}_pooled_posterior.csv"),
        os.path.join(settings['io']['MLE_combinations_directory'],
                     f"{combo_name}_pooled_MLE.csv"),
        os.path.join(settings['io']['rolling_window_combinations_directory'],
                     f"{combo_name}_pooled_rollingMLE.csv"),
    ]
    return all(os.path.exists(p) for p in expected)


def apply_combo_cutoff_overrides(settings, combo_cfg):
    """
    Apply optional per-combination cutoff overrides (n_components, labels,
    initial_centers) to settings['cutoff'] in place.

    Mirrors the override mechanism in cutoff_functions for groups.yaml.
    """
    if 'n_components' in combo_cfg:
        n_comp = combo_cfg['n_components']
        settings['cutoff']['n_components'] = (
            'auto' if str(n_comp).lower() == 'auto' else int(n_comp))
    if 'labels' in combo_cfg:
        settings['cutoff']['labels'] = combo_cfg['labels']
    if 'initial_centers' in combo_cfg:
        settings['cutoff']['initial_centers_log10'] = combo_cfg['initial_centers']
