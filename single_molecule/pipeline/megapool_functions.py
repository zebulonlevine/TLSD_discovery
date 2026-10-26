## megapool_functions.py
##
## Helpers for joint-saSPT analysis across user-defined megapools of
## conditions and/or combinations, declared in `megapools.yaml` at the
## experiment root.
##
## A megapool concatenates trajectory CSVs from multiple member units
## (raw conditions or already-pooled combination nicknames), runs ONE
## joint saSPT fit on the combined detections, and then partitions the
## resulting per-subtrack MLE_D values back to source members. This
## yields per-member distributions all judged against the same joint
## occupations posterior — the statistically clean form of
## cross-condition comparison.
##
## Rolling-window analysis is intentionally NOT run for megapools.
## Bleaching curves and example movies remain available in the
## per-condition / per-combination plots that already exist upstream.
##
## Key invariant (verified against the saSPT runtime):
##   saSPT's `orig_trajectory` column in its MLE output equals the input
##   trajectory ID. When `sample_size` >= n_input_trajectories, every
##   input trajectory is preserved with a 1:1 mapping. We rely on this
##   for source-member recovery and set sample_size to a huge value to
##   force this regime.
##
## To partition by source member, we keep a sidecar source map written
## at pool time and join MLE_df against it via the input trajectory ID
## extracted from saSPT's orig_trajectory string.


import copy
import glob
import os
import re

import numpy as np
import pandas as pd
import yaml


# Regex matches the trailing ::<integer> that executeSASPT appends to
# orig_trajectory when it writes the MLE CSV.
TRAJ_ID_FROM_ORIG = re.compile(r'::(\d+)$')


# ============================================================
# YAML loading
# ============================================================

def load_megapools(yaml_path):
    """
    Load megapool definitions from a YAML file.

    Expected format:
        megapools:
          comparison_set_1:
            members:
              - sgNT                          # combo nickname OR
              - line=L29Avi_Genetic=...       # raw condition name
            max_trajectories_per_member: 100000   # optional
            max_trajectories_total: 1000000       # optional
            balance_to_min_trajectories: false    # optional; if true and
                                                   # max_trajectories_per_member
                                                   # is not set, caps every
                                                   # member to the SMALLEST
                                                   # member's post-filter
                                                   # trajectory count (whole
                                                   # trajectories, not
                                                   # subtracks), so the
                                                   # lowest-count condition
                                                   # isn't outweighed
            random_seed: 42                       # optional
            n_components: 4                       # optional GMM override
            labels: [immobile, slow, fast]        # optional GMM override
            initial_centers: [-1.3, -0.3, 0.5]    # optional GMM override

    Returns
    -------
    dict : {megapool_name: {'members': [...], ...}}
        Empty dict if the file does not exist.
    """
    if not os.path.exists(yaml_path):
        return {}

    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f) or {}

    pools = config.get('megapools', {})

    # Validation
    forbidden = set('/\\:*?"<>|\t\n\r ')
    for name, cfg in pools.items():
        if cfg is None:
            raise ValueError(
                f"Megapool '{name}' has no body. Either remove the "
                f"entry or add a 'members:' list under it.")
        if 'members' not in cfg:
            raise ValueError(f"Megapool '{name}' missing 'members' key")
        if not isinstance(cfg['members'], list) or len(cfg['members']) == 0:
            raise ValueError(
                f"Megapool '{name}' must have a non-empty list of members")
        bad = forbidden & set(name)
        if bad:
            raise ValueError(
                f"Megapool name '{name}' contains characters not safe "
                f"for filenames: {sorted(bad)}. Replace with underscores.")
    return pools


def find_megapool(megapool_name, yaml_path):
    """Return (name, cfg) or (None, None)."""
    pools = load_megapools(yaml_path)
    if megapool_name in pools:
        return megapool_name, pools[megapool_name]
    return None, None


# ============================================================
# Member resolution
# ============================================================

def resolve_member_to_traj_csvs(member, settings, combinations_cfg=None):
    """
    Given a member name (a combo nickname or a raw condition), return
    a list of trajectory CSV paths to read for that member.

    For a combination, returns the single pre-pooled CSV in
    traj_combinations/. For a raw condition, returns the per-FOV CSVs
    matching that condition's filename prefix in trajectories/.

    Returns
    -------
    (kind, paths) : (str, list of str)
    """
    if combinations_cfg and member in combinations_cfg:
        combo_csv = os.path.join(
            settings['io']['traj_combinations_directory'],
            f"{member}_pooled_traj.csv")
        if not os.path.exists(combo_csv):
            raise FileNotFoundError(
                f"Member '{member}' is a combination but its pooled "
                f"trajectory CSV is missing: {combo_csv}. "
                f"Run combinationwise processing first.")
        return ('combination', [combo_csv])

    traj_dir = settings['io']['traj_directory']
    # Per-FOV files follow {cond}_{NNN}_traj.csv
    pat = os.path.join(traj_dir, f"{member}_*_traj.csv")
    paths = sorted(glob.glob(pat))
    if not paths:
        raise FileNotFoundError(
            f"Member '{member}' is neither in combinations.yaml nor a "
            f"recognized condition: no files matching {pat}")
    return ('condition', paths)


# ============================================================
# Pooling
# ============================================================

def pre_filter_trajectories(traj_df, min_length=2,
                             dropna_cols=('x', 'y'),
                             verbose=False):
    """
    Drop trajectories that would produce no subtracks or have data
    integrity issues.

    Two criteria, both empirically motivated by saSPT's behavior
    (verified by direct test):

    1. Drop rows with NaN in `dropna_cols` (default 'x' and 'y').
       saSPT does not crash on NaN but issues a RuntimeWarning and
       produces degenerate (uninformative) subtracks for NaN-containing
       trajectories. We drop them for analysis cleanliness.

    2. Drop trajectories with length < `min_length` (default 2). saSPT
       requires at least 2 detections (1 jump) to produce a subtrack;
       length-1 trajectories yield zero subtracks and only waste
       random-sample slots when capping. Verified by direct test:
       length=1 → 0 subtracks; length=2..7 → 1 subtrack; etc.

    Parameters
    ----------
    traj_df : pd.DataFrame
        Must have 'trajectory' column. Other required columns depend
        on `dropna_cols`.
    min_length : int
        Minimum number of detections per trajectory (after NaN drop).
    dropna_cols : tuple of str
        Columns checked for NaN; rows with NaN in any of these are
        dropped. Default ('x', 'y') is the minimum for saSPT.
    verbose : bool
        Print filter stats.

    Returns
    -------
    (filtered_df, stats) : (pd.DataFrame, dict)
        filtered_df has the same column structure as traj_df.
        stats has keys:
          n_detections_before, n_detections_after_nan,
          n_detections_dropped_nan, n_trajectories_before,
          n_trajectories_after_nan, n_trajectories_kept,
          n_trajectories_dropped_short
    """
    n_det_before = len(traj_df)
    n_traj_before = traj_df['trajectory'].nunique()

    # Drop rows with NaN in critical columns
    dropna_cols = list(dropna_cols)
    missing = [c for c in dropna_cols if c not in traj_df.columns]
    if missing:
        raise ValueError(
            f"pre_filter_trajectories: dropna_cols {missing} not in "
            f"traj_df columns {list(traj_df.columns)}")
    filtered = traj_df.dropna(subset=dropna_cols)
    n_det_after_nan = len(filtered)
    n_traj_after_nan = filtered['trajectory'].nunique()

    # Compute post-NaN-drop length per trajectory and filter short ones
    lengths = filtered.groupby('trajectory').size()
    keep_ids = lengths[lengths >= min_length].index
    filtered = filtered[filtered['trajectory'].isin(keep_ids)].copy()
    n_traj_kept = filtered['trajectory'].nunique()

    stats = {
        'n_detections_before': int(n_det_before),
        'n_detections_after_nan': int(n_det_after_nan),
        'n_detections_dropped_nan': int(n_det_before - n_det_after_nan),
        'n_trajectories_before': int(n_traj_before),
        'n_trajectories_after_nan': int(n_traj_after_nan),
        'n_trajectories_kept': int(n_traj_kept),
        'n_trajectories_dropped_short': int(n_traj_after_nan - n_traj_kept),
        'min_length': int(min_length),
        'dropna_cols': list(dropna_cols),
    }

    if verbose:
        print(f"    pre_filter: {n_det_before:,} det / {n_traj_before:,} trajs "
              f"-> drop {stats['n_detections_dropped_nan']:,} NaN det "
              f"({100*stats['n_detections_dropped_nan']/max(n_det_before,1):.3f}%) "
              f"-> drop {stats['n_trajectories_dropped_short']:,} trajs with "
              f"length<{min_length} "
              f"({100*stats['n_trajectories_dropped_short']/max(n_traj_after_nan,1):.2f}%) "
              f"-> {n_traj_kept:,} trajs kept "
              f"({100*n_traj_kept/max(n_traj_before,1):.1f}% of input)")

    return filtered, stats


def _load_combo_fov_map(combo_csv_path):
    """
    Look for a sibling fov_map sidecar next to a combo's pooled CSV.
    Combo fov_maps are written by combination_functions.pool_combination_trajectories
    and have columns trajectory,fov_basename (and optionally a 'condition'
    column from the underlying conditionwise step — we ignore that here).

    Returns None if no sidecar exists.
    """
    parent = os.path.dirname(combo_csv_path)
    base = os.path.basename(combo_csv_path)
    # combo CSV is {combo}_pooled_traj.csv → sidecar is {combo}_fov_map.csv
    if base.endswith('_pooled_traj.csv'):
        combo = base[:-len('_pooled_traj.csv')]
        sidecar = os.path.join(parent, f"{combo}_fov_map.csv")
        if os.path.exists(sidecar):
            return pd.read_csv(sidecar)
    return None


def _count_member_trajectories(member, settings, combinations_cfg=None,
                               pre_filter=True, pre_filter_min_length=2,
                               pre_filter_dropna_cols=('x', 'y')):
    """Read + filter one member's trajectory CSVs the same way
    pool_megapool_trajectories does, but only to COUNT post-filter
    trajectories -- used by balance_to_min_trajectories to determine the
    common per-member cap before the real pooling pass. Duplicates some
    I/O against the main loop; acceptable since this runs once per
    megapool at pool-build time, not a hot path."""
    kind, paths = resolve_member_to_traj_csvs(
        member, settings, combinations_cfg=combinations_cfg)
    member_offset = 0
    pieces = []
    for p in paths:
        df = pd.read_csv(p)
        if 'trajectory' not in df.columns:
            continue
        df = df.copy()
        df['trajectory'] = df['trajectory'] + member_offset
        member_offset = int(df['trajectory'].max()) + 1
        pieces.append(df)
    if not pieces:
        return 0
    member_df = pd.concat(pieces, ignore_index=True)
    if pre_filter:
        member_df, _ = pre_filter_trajectories(
            member_df, min_length=pre_filter_min_length,
            dropna_cols=pre_filter_dropna_cols, verbose=False)
    return member_df['trajectory'].nunique()


def pool_megapool_trajectories(megapool_name, members, settings,
                                combinations_cfg=None,
                                max_per_member=None, max_total=None,
                                balance_to_min_trajectories=False,
                                random_seed=42,
                                pre_filter=True,
                                pre_filter_min_length=2,
                                pre_filter_dropna_cols=('x', 'y')):
    """
    Pool trajectory CSVs from a megapool's member units, with optional
    per-member and total trajectory caps. Returns paths to the pooled
    CSV and source map sidecar.

    Trajectory IDs are offset to keep them globally unique across the
    pool. saSPT's `orig_trajectory` preserves this offset-pool ID,
    letting us partition results by source member.

    Pipeline per member:
      1. Read trajectory CSVs and intra-member offset
      2. (Optional) pre-filter: drop NaN-containing rows, drop short
         trajectories. By default ON; min_length=2 (empirically
         determined as saSPT's minimum useful length).
      3. Per-member cap: if max_per_member is set and the FILTERED
         trajectory count exceeds it, subsample to max_per_member.
      4. Offset to global pool space.
      5. Record source mapping and (when possible) FOV mapping.

    balance_to_min_trajectories: if True and max_per_member is not
    explicitly set, a pre-pass counts each member's POST-FILTER
    trajectory count (whole trajectories, not subtracks -- subsampling
    always operates at the trajectory level here, before saSPT's own
    internal splitsize-based subtrack splitting) and uses the SMALLEST
    member's count as a uniform max_per_member, so the lowest-count
    condition isn't outweighed by larger ones in the joint fit. Ignored
    if max_per_member is explicitly given (explicit wins).

    After all members: optional global max_total cap.

    Sidecars written:
      - {megapool}_source_map.csv: trajectory → source_member
      - {megapool}_fov_map.csv:    trajectory → source_member, fov_basename

    Returns
    -------
    (pooled_csv_path, source_map_path) : (str, str)
    """
    rng = np.random.default_rng(random_seed)

    if balance_to_min_trajectories and max_per_member is None:
        counts = {}
        for member in members:
            counts[member] = _count_member_trajectories(
                member, settings, combinations_cfg=combinations_cfg,
                pre_filter=pre_filter,
                pre_filter_min_length=pre_filter_min_length,
                pre_filter_dropna_cols=pre_filter_dropna_cols)
        max_per_member = min(counts.values()) if counts else None
        print(f"  balance_to_min_trajectories: per-member post-filter counts "
              f"{counts} -> max_per_member={max_per_member}")

    out_dir = settings['io']['traj_megapools_directory']
    os.makedirs(out_dir, exist_ok=True)
    pooled_csv = os.path.join(out_dir, f"{megapool_name}_pooled_traj.csv")
    source_map_csv = os.path.join(out_dir, f"{megapool_name}_source_map.csv")
    fov_map_csv = os.path.join(out_dir, f"{megapool_name}_fov_map.csv")

    pieces = []
    source_map_rows = []
    fov_map_rows = []
    offset = 0
    member_stats = []

    for member in members:
        kind, paths = resolve_member_to_traj_csvs(
            member, settings, combinations_cfg=combinations_cfg)

        # Concatenate this member's FOV CSVs with intra-member offsetting.
        # Track per-trajectory FOV identity as we go.
        member_pieces = []
        # member_traj_id_to_fov: combo-local trajectory ID → FOV name
        member_traj_id_to_fov = {}
        member_offset = 0

        # Load combo-level fov_map if applicable
        combo_fov_map_df = None
        if kind == 'combination':
            combo_fov_map_df = _load_combo_fov_map(paths[0])
            if combo_fov_map_df is not None:
                # Build lookup: combo-traj-id → fov_basename
                member_traj_id_to_fov = dict(zip(
                    combo_fov_map_df['trajectory'].astype(int),
                    combo_fov_map_df['fov_basename']))

        for p in paths:
            df = pd.read_csv(p)
            if 'trajectory' not in df.columns:
                continue
            df = df.copy()
            df['trajectory'] = df['trajectory'] + member_offset

            if kind == 'condition':
                # For raw conditions, each path is one per-FOV CSV
                # → all its trajectories belong to this FOV
                fov_basename = os.path.basename(p)
                for tid in df['trajectory'].unique():
                    member_traj_id_to_fov[int(tid)] = fov_basename

            member_offset = int(df['trajectory'].max()) + 1
            member_pieces.append(df)

        if not member_pieces:
            print(f"  WARNING: member '{member}' contributed no detections")
            continue

        member_df = pd.concat(member_pieces, ignore_index=True)
        n_pre_filter = member_df['trajectory'].nunique()

        # Pre-filter: drop NaN-containing rows and trajectories with
        # length < min_length. Empirically determined to match saSPT's
        # minimum useful trajectory length.
        if pre_filter:
            member_df, filter_stats = pre_filter_trajectories(
                member_df,
                min_length=pre_filter_min_length,
                dropna_cols=pre_filter_dropna_cols,
                verbose=False)
            n_post_filter = filter_stats['n_trajectories_kept']
            if n_post_filter < n_pre_filter:
                print(f"  pre-filter '{member}': "
                      f"{n_pre_filter:,} -> {n_post_filter:,} trajectories "
                      f"(dropped {filter_stats['n_detections_dropped_nan']:,} "
                      f"NaN det, {filter_stats['n_trajectories_dropped_short']:,} "
                      f"short trajs)")
        else:
            filter_stats = None
            n_post_filter = n_pre_filter

        if len(member_df) == 0:
            print(f"  WARNING: member '{member}' empty after pre-filter")
            continue

        member_traj_ids = member_df['trajectory'].unique()
        n_member_trajs = len(member_traj_ids)

        # Per-member cap (random subsample of trajectories) — applied to
        # POST-filter trajectories so the cap means "N usable trajectories"
        if max_per_member is not None and n_member_trajs > max_per_member:
            keep = rng.choice(member_traj_ids, size=max_per_member,
                              replace=False)
            member_df = member_df[member_df['trajectory'].isin(keep)]
            n_kept = max_per_member
            print(f"  member '{member}' ({kind}): {n_member_trajs:,} "
                  f"trajectories → sampled to {n_kept:,}")
        else:
            n_kept = n_member_trajs
            print(f"  member '{member}' ({kind}): {n_member_trajs:,} "
                  f"trajectories (all kept)")

        # Offset to global pool space. The combo-local trajectory ID
        # of a row was already used to build member_traj_id_to_fov;
        # save that lookup keyed on the GLOBAL pool ID before applying
        # the offset shift.
        member_df = member_df.copy()
        # Snapshot pre-offset IDs so we can map to FOV info that uses
        # combo-local keys
        pre_offset_traj = member_df['trajectory'].values
        member_df['trajectory'] = member_df['trajectory'] + offset
        new_offset = int(member_df['trajectory'].max()) + 1
        pieces.append(member_df)

        # Record source mapping and FOV mapping
        if member_traj_id_to_fov:
            # combo-local → global pool ID via the offset
            for tid in member_df['trajectory'].unique():
                local_tid = int(tid) - offset
                fov = member_traj_id_to_fov.get(local_tid, '<unknown>')
                source_map_rows.append({'trajectory': int(tid),
                                        'source_member': member})
                fov_map_rows.append({'trajectory': int(tid),
                                      'source_member': member,
                                      'fov_basename': fov})
        else:
            # Combo with no fov_map sidecar — record source_member but
            # leave fov unknown
            for tid in member_df['trajectory'].unique():
                source_map_rows.append({'trajectory': int(tid),
                                        'source_member': member})
                fov_map_rows.append({'trajectory': int(tid),
                                      'source_member': member,
                                      'fov_basename': '<no_fov_map>'})

        member_stats.append({
            'member': member, 'kind': kind,
            'n_input_trajs': n_member_trajs,
            'n_kept': n_kept,
            'n_detections': len(member_df),
            'traj_id_range_lo': offset,
            'traj_id_range_hi': new_offset - 1,
            'has_fov_map': bool(member_traj_id_to_fov),
        })
        offset = new_offset

    if not pieces:
        raise RuntimeError(
            f"Megapool '{megapool_name}': no trajectories pooled from "
            f"members {members}")

    pooled = pd.concat(pieces, ignore_index=True)
    source_map = pd.DataFrame(source_map_rows)
    fov_map = pd.DataFrame(fov_map_rows)

    # Global cap (uniform random across pooled trajectories)
    if max_total is not None and len(source_map) > max_total:
        n_before = len(source_map)
        keep_ids = rng.choice(source_map['trajectory'].values,
                              size=max_total, replace=False)
        keep_set = set(int(x) for x in keep_ids)
        pooled = pooled[pooled['trajectory'].isin(keep_set)]
        source_map = source_map[source_map['trajectory'].isin(keep_set)]
        fov_map = fov_map[fov_map['trajectory'].isin(keep_set)]
        print(f"  global cap: {n_before:,} → {len(source_map):,} "
              f"trajectories")

    pooled.to_csv(pooled_csv, index=False)
    source_map.to_csv(source_map_csv, index=False)
    fov_map.to_csv(fov_map_csv, index=False)

    n_with_fov = int((fov_map['fov_basename'] != '<no_fov_map>').sum())
    print(f"\n  Pooled {len(pooled):,} detections from "
          f"{len(source_map):,} trajectories across {len(member_stats)} "
          f"members → {pooled_csv}")
    print(f"  Source map written to {source_map_csv}")
    print(f"  FOV map  written to {fov_map_csv} "
          f"({n_with_fov:,}/{len(fov_map):,} trajectories with FOV info)")
    return pooled_csv, source_map_csv


# ============================================================
# Settings redirection
# ============================================================

def redirect_settings_for_megapool(settings):
    """
    Return a deep-copied settings dict with post/MLE I/O keys pointing
    to megapool directories. Sets saspt.sample_size to a huge value so
    saSPT does not subsample internally — this preserves the 1:1
    orig_trajectory ↔ input_trajectory mapping we rely on for source
    partitioning.

    Does NOT redirect rolling_window_directory: rolling-window analysis
    is skipped for megapools entirely.
    """
    new = copy.deepcopy(settings)
    new['io']['post_directory'] = settings['io']['post_megapools_directory']
    new['io']['MLE_directory'] = settings['io']['MLE_megapools_directory']
    new['saspt'] = dict(new['saspt'])
    new['saspt']['sample_size'] = int(1e9)
    return new


# ============================================================
# Skip-if-exists
# ============================================================

def megapool_outputs_exist(megapool_name, settings):
    """True iff the two joint saSPT outputs already exist."""
    expected = [
        os.path.join(settings['io']['post_megapools_directory'],
                     f"{megapool_name}_pooled_posterior.csv"),
        os.path.join(settings['io']['MLE_megapools_directory'],
                     f"{megapool_name}_pooled_MLE.csv"),
    ]
    return all(os.path.exists(p) for p in expected)


# ============================================================
# Partitioning by source member
# ============================================================

def _extract_input_traj_id(orig_trajectory_value):
    """
    Extract the trailing integer (the input trajectory ID) from saSPT's
    orig_trajectory string of the form 'directory::file.csv::<INT>'.
    Returns int or None.
    """
    m = TRAJ_ID_FROM_ORIG.search(str(orig_trajectory_value))
    return int(m.group(1)) if m else None


def partition_mle_by_source(mle_csv, source_map_csv):
    """
    Load a megapool MLE CSV, extract input trajectory IDs from
    orig_trajectory, and join with the source map. Returns the
    augmented DataFrame and a diagnostic dict.

    Returns
    -------
    (mle_df, stats) : (pd.DataFrame, dict)
        mle_df has new columns '_input_traj_id' (int) and 'source_member'.
    """
    mle = pd.read_csv(mle_csv)
    sm = pd.read_csv(source_map_csv)

    mle['_input_traj_id'] = mle['orig_trajectory'].apply(_extract_input_traj_id)
    n_total = len(mle)
    n_extracted = int(mle['_input_traj_id'].notna().sum())

    src_dict = dict(zip(sm['trajectory'].astype(int).tolist(),
                        sm['source_member'].tolist()))
    mle['source_member'] = mle['_input_traj_id'].map(
        lambda x: src_dict.get(int(x)) if pd.notna(x) else None)

    n_matched = int(mle['source_member'].notna().sum())
    n_unmatched = n_total - n_matched

    stats = {
        'n_mle_rows': n_total,
        'n_traj_id_extracted': n_extracted,
        'n_matched_to_source': n_matched,
        'n_unmatched': n_unmatched,
        'per_member_counts': mle['source_member'].value_counts().to_dict(),
    }
    return mle, stats


def partition_posterior_by_source(posterior_csv, mle_df_with_source):
    """
    Load the megapool's JOINT posterior CSV (one shared P(D) from the
    joint fit across all members -- not any one member's own density).

    This function only loads it; the actual per-member synthesis (KDE
    on each member's own MLE_D values) happens downstream in
    write_member_posterior_csv when it's given mle_df_with_source. This
    function's `mle_df_with_source` argument is unused directly but kept
    for call-site symmetry with write_member_posterior_csv.

    Returns
    -------
    pd.DataFrame
        The joint posterior, unchanged.
    """
    post = pd.read_csv(posterior_csv)
    return post


def write_member_mle_csv(mle_df_with_source, megapool_name, member,
                          settings, fov_map_csv=None):
    """
    Write one member's MLE rows (extracted from the joint MLE) to a
    per-member CSV. Output filename matches the pattern
    process_condition's glob expects ({cond}_pooled_*.csv).

    Also writes a per-member fov_map sidecar to the plotting directory
    so SPTExperiment can find it via the standard fov_map lookup. The
    per-member sidecar has columns [trajectory, fov_basename] keyed on
    the saSPT input trajectory ID (i.e., the trailing integer of
    orig_trajectory).

    If `fov_map_csv` is not supplied, looks for the convention path
    {traj_megapools_directory}/{megapool_name}_fov_map.csv.
    """
    member_dir = os.path.join(settings['io']['MLE_megapools_directory'])
    os.makedirs(member_dir, exist_ok=True)
    out_path = os.path.join(member_dir,
                            f"{megapool_name}_{member}_pooled_MLE.csv")
    member_rows = mle_df_with_source[
        mle_df_with_source['source_member'] == member].copy()
    # Drop our scratch columns before writing so the file looks
    # identical to a normal saSPT MLE output (process_condition reads
    # source_file column).
    drop_cols = [c for c in ('_input_traj_id', 'source_member')
                 if c in member_rows.columns]
    member_rows_for_disk = member_rows.drop(columns=drop_cols)
    member_rows_for_disk.to_csv(out_path, index=False)

    # Auto-detect fov_map by convention if not passed
    if fov_map_csv is None:
        fov_map_csv = os.path.join(
            settings['io']['traj_megapools_directory'],
            f"{megapool_name}_fov_map.csv")

    # Write per-member fov_map sidecar if we have FOV info to propagate.
    if os.path.exists(fov_map_csv):
        try:
            full_fov_map = pd.read_csv(fov_map_csv)
            # Filter to this member's trajectories
            member_fov_map = full_fov_map[
                full_fov_map['source_member'] == member]
            # Keep only the columns SPTExperiment expects
            keep_cols = [c for c in ('trajectory', 'fov_basename')
                          if c in member_fov_map.columns]
            member_fov_map = member_fov_map[keep_cols]
            sidecar_dir = os.path.join(
                settings['io']['plot_directory'], 'plotting_pkls')
            os.makedirs(sidecar_dir, exist_ok=True)
            sidecar_path = os.path.join(
                sidecar_dir,
                f"{megapool_name}_{member}_fov_map.csv")
            member_fov_map.to_csv(sidecar_path, index=False)
        except Exception as e:
            print(f"  WARNING: could not write member fov_map for "
                  f"{megapool_name}/{member}: {e}")

    return out_path, len(member_rows)


def write_member_posterior_csv(posterior_df, megapool_name, member,
                                settings, mle_df_with_source=None, bw=0.1):
    """
    Write a member's posterior CSV.

    If `mle_df_with_source` is given, synthesizes a genuine per-member
    density via Gaussian KDE (scipy, same method/bandwidth as
    spt_compare.overlay_mle_kde) on that member's own MLE_D values,
    evaluated on the joint posterior's D grid (for direct comparability
    across members and against the joint fit). This is a KDE-smoothed
    proxy over hard per-subtrack MLE_D (argmax) assignments, NOT a true
    responsibility-weighted Bayesian posterior -- subtrack-level soft
    posteriors aren't available for megapool members (megapool never
    calls save_subtrack_posteriors; see docs/MEGAPOOL.md). Falls back to
    the (shared, non-member-specific) joint posterior with a warning if
    the member has too few finite MLE_D points to KDE meaningfully
    (<10, matching overlay_mle_kde's own minimum).

    If `mle_df_with_source` is not given, falls back to the old
    behavior: duplicate the joint posterior verbatim (all members would
    look identical -- kept only for backward compatibility, don't rely
    on it for real per-member comparisons).
    """
    member_dir = os.path.join(settings['io']['post_megapools_directory'])
    os.makedirs(member_dir, exist_ok=True)
    out_path = os.path.join(member_dir,
                            f"{megapool_name}_{member}_pooled_posterior.csv")

    if mle_df_with_source is None:
        posterior_df.to_csv(out_path, index=False)
        return out_path

    from scipy.stats import gaussian_kde
    member_mle = mle_df_with_source[mle_df_with_source['source_member'] == member]
    log_mle = np.log10(member_mle['MLE_D'].values)
    log_mle = log_mle[np.isfinite(log_mle)]

    if len(log_mle) < 10:
        print(f"    WARNING: member '{member}' has only {len(log_mle)} "
              f"finite MLE_D points -- too few to KDE a per-member "
              f"posterior (need >=10, matching overlay_mle_kde's own "
              f"minimum). Falling back to the joint (shared, non-member-"
              f"specific) posterior for this member.")
        out_df = posterior_df.copy()
    else:
        D_grid = np.sort(posterior_df['D'].unique())
        log_grid = np.log10(D_grid)
        kde = gaussian_kde(log_mle, bw_method=bw)
        density = kde(log_grid)
        density = density / density.sum()
        out_df = pd.DataFrame({'D': D_grid, 'density': density,
                               'file': f"{megapool_name}_{member}_pooled_traj.csv"})

    out_df.to_csv(out_path, index=False)
    return out_path
