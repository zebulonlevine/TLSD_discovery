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
import re

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


def _discover_traj_by_cond(traj_dir):
    """
    Discover per-FOV trajectory CSVs directly from a traj_dir,
    grouping by condition. Mirrors trackingCodeFunctions.
    group_files_by_metadata but operates on `<cond>_NNN_traj.csv` files
    found on disk rather than `<cond>_NNN.<ext>` source files.

    Returns
    -------
    dict {condition: [traj_csv_filename, ...]}
        Empty dict if traj_dir does not exist.

    Notes
    -----
    Mostly an escape hatch for the ER pipeline, where
    settings['io']['data_filepath'] = 'preprocessed/' contains
    `_NNN_spt.tif` and `_NNN_er.tif` files that don't match the
    `(.*)_\\d{3}.tif$` regex group_files_by_metadata expects. The traj
    CSVs themselves are named canonically (`<cond>_NNN_traj.csv`) so we
    can use them as the source of truth instead.
    """
    from collections import defaultdict
    out = defaultdict(list)
    if not os.path.isdir(traj_dir):
        return dict(out)
    pat = re.compile(r"(.*)_\d{3}_traj\.csv$")
    for fname in os.listdir(traj_dir):
        if not fname.endswith('_traj.csv'):
            continue
        m = pat.match(fname)
        if m:
            out[m.group(1)].append(fname)
    return dict(out)


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

    # Build a fallback "discover traj CSVs directly from traj_dir" map.
    # This is needed when settings['io']['data_filepath'] points at a
    # directory that doesn't contain canonically-named source files —
    # e.g., the ER pipeline's preprocessed/ holds `<cond>_NNN_spt.tif`
    # and `<cond>_NNN_er.tif`, neither of which matches the
    # `(.*)_\d{3}.tif$` regex used by group_files_by_metadata. The
    # actual traj_dir always has the canonical `<cond>_NNN_traj.csv`
    # form because filewise wrote it that way. Discover from there.
    traj_by_cond = _discover_traj_by_cond(traj_dir)

    out_dir = settings['io']['traj_combinations_directory']
    os.makedirs(out_dir, exist_ok=True)

    pieces = []
    fov_map_rows = []   # collect (trajectory_id, source_FOV_basename)
    offset = 0
    files_used = 0
    files_missing = 0
    files_empty = 0   # exist but no trajectories in this zone

    for cond in member_conditions:
        # Resolve the list of per-FOV traj CSVs for this condition.
        # Prefer the on-disk discovery (authoritative — these are the
        # files we'll actually open). Fall back to nd2_by_cond only if
        # traj-dir discovery returned nothing, e.g. because the
        # condition uses a non-`_NNN_` naming scheme that the disk
        # discovery regex doesn't match.
        if cond in traj_by_cond:
            traj_filenames = traj_by_cond[cond]
        elif cond in nd2_by_cond:
            # Reconstruct traj filenames from source basenames; tolerate
            # ER-pipeline suffixes (_spt, _er) on the source side that
            # don't appear on the traj-CSV side.
            traj_filenames = []
            for src_file in nd2_by_cond[cond]:
                base = os.path.splitext(src_file)[0]
                base = re.sub(r'_(spt|er)$', '', base)
                traj_filenames.append(base + '_traj.csv')
        else:
            print(f"  WARNING: combination '{combo_name}' references "
                  f"condition '{cond}' which has no traj CSVs in "
                  f"{traj_dir} and no source files in data_filepath.")
            continue

        for traj_fname in traj_filenames:
            traj_path = os.path.join(traj_dir, traj_fname)
            if not os.path.exists(traj_path):
                files_missing += 1
                print(f"    missing trajectory CSV: {traj_path}")
                continue

            df = pd.read_csv(traj_path)
            if 'trajectory' not in df.columns:
                print(f"    skipping {traj_path}: no 'trajectory' column")
                continue
            if len(df) == 0:
                # Zone exports can legitimately be empty for a given FOV
                # (e.g. an FOV with zero on_structure trajectories).
                # Counted but not fatal.
                files_empty += 1
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
        # Build an informative diagnosis. The three failure modes are
        # distinguishable from the counters: a missing zone directory, a
        # filename-pattern mismatch (dir exists but no matching files),
        # or every-FOV-was-empty in the requested zone.
        diagnosis = []
        diagnosis.append(
            f"Combination '{combo_name}': no usable trajectory data.")
        diagnosis.append(f"  zone:     {zone if zone is not None else '(none, canonical)'}")
        diagnosis.append(f"  traj_dir: {traj_dir}")
        diagnosis.append(f"  exists:   {os.path.isdir(traj_dir)}")
        if os.path.isdir(traj_dir):
            n_in_dir = len([f for f in os.listdir(traj_dir)
                            if f.endswith('_traj.csv')])
            diagnosis.append(f"  *_traj.csv files in traj_dir: {n_in_dir}")
        diagnosis.append(f"  member_conditions: {member_conditions}")
        for cond in member_conditions:
            n_disk = len(traj_by_cond.get(cond, []))
            n_src = len(nd2_by_cond.get(cond, []))
            diagnosis.append(
                f"    {cond}: {n_disk} traj CSV(s) discovered on disk, "
                f"{n_src} source ND2/TIF(s) in data_filepath")
        diagnosis.append(
            f"  counts: missing={files_missing}, empty(zone)={files_empty}, "
            f"used={files_used}")
        if files_missing > 0 and files_empty == 0:
            diagnosis.append(
                "  -> Looks like per-FOV traj CSVs were not written for "
                "this zone. Re-run filewise (run_filewise_er.py) so it "
                f"populates trajectories_{zone}/.")
        elif files_empty > 0 and files_missing == 0:
            diagnosis.append(
                "  -> All FOVs returned zero trajectories in this zone. "
                "Either the zone classification is too stringent or this "
                "combo has genuinely empty data; widen "
                "fraction_threshold/enrichment_alpha, or drop this combo.")
        raise RuntimeError("\n".join(diagnosis))

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
