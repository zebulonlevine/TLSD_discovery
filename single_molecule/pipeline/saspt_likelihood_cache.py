"""
saspt_likelihood_cache.py
==========================

Build ONE label-blind, per-subtrack saSPT likelihood cache from a pooled
trajectory CSV (a combinationwise or megapool pool), for reuse by every
future two-group occupancy comparison drawn from that pool
(occupancy_compare.py :: compare_groups / refit_permutation).

Gate-0 findings this module relies on (saspt 0.4.0; see
docs/differentiability_handoff.md section 2):
  * `StateArray.from_detections(df, **saspt_settings)` lazily computes
    `naive_assignment_probabilities` (shape (n_D, n_loc_err, n_subtracks),
    normalized per-subtrack over the full grid) via
    `likelihood.exp(likelihood(trajectories))` -- this is the label-blind,
    no-refit quantity we cache. It never touches VBEM/`infer_posterior`
    unless `.posterior_*` is explicitly accessed, so no `max_iter=0` is
    strictly required, but we set it anyway to make the "we do not refit
    here" intent explicit and cheap to grep for.
  * There is no supported API to inject a precomputed likelihood array
    into a fresh StateArray -- `StateArray.__init__` always (lazily) calls
    `likelihood(trajectories)` itself. We do not fabricate one; instead we
    compute the likelihood ONCE here and cache the resulting array
    directly, bypassing StateArray entirely for every later reuse.
  * Per-subtrack identity: `StateArray.trajectories.detections` carries
    both the post-split contiguous `trajectory` column (0..n_subtracks-1,
    the axis order used by `naive_assignment_probabilities`) and the
    pre-split `orig_trajectory` column (the raw `trajectory` id in the
    input pooled CSV, which is what `{pool}_fov_map.csv` is keyed on).
    This lets us join subtrack -> FOV directly, without the `orig_trajectory
    :: <int>` string round-trip the CSV-based pipeline stages use (that
    round-trip is for reading back *already-written* MLE/posterior CSVs;
    we never write one here).

Invariants enforced (see handoff section 1):
  * The saSPT settings block (likelihood_type, pixel_size_um,
    frame_interval, focal_depth, splitsize, diff_coefs, loc_errors) is
    inherited VERBATIM from the deployed settings.pkl. Only `sample_size`
    is overridden (forced to 1e9, matching the established
    combination/megapool convention -- see combination_functions.py:380,
    megapool_functions.py:557) so that every subtrack is retained
    deterministically.
  * Grid identity is checksummed on build and asserted on load.
"""

import copy
import hashlib
import json
import os
import pickle

import numpy as np
import pandas as pd

# The subset of settings['saspt'] that defines the likelihood grid itself.
# sample_size is deliberately excluded: we always override it to 1e9 here,
# so it is not part of "did the grid change" identity.
_GRID_KEYS = (
    'likelihood_type', 'pixel_size_um', 'frame_interval', 'focal_depth',
    'splitsize', 'diff_coefs', 'loc_errors',
)

FORCED_SAMPLE_SIZE = int(1e9)


def _grid_checksum(saspt_settings):
    """Stable checksum over the grid-defining subset of saspt settings."""
    payload = {}
    for k in _GRID_KEYS:
        v = saspt_settings.get(k)
        if isinstance(v, np.ndarray):
            v = v.tolist()
        payload[k] = v
    blob = json.dumps(payload, sort_keys=True).encode('utf-8')
    return hashlib.sha256(blob).hexdigest()


def _load_saspt_settings(settings_pkl):
    with open(settings_pkl, 'rb') as f:
        settings = pickle.load(f)
    sa = dict(settings['saspt'])
    # StateArray/StateArrayParameters does not accept this pipeline-internal key
    sa.pop('save_subtrack_posteriors', None)
    return sa


def _load_min_I0(settings_pkl):
    with open(settings_pkl, 'rb') as f:
        settings = pickle.load(f)
    return settings.get('quot', {}).get('track', {}).get('min_I0')


def build_likelihood_cache(pooled_traj_csv, settings_pkl, out_dir,
                            fov_map_csv=None, cache_name=None,
                            min_track_length=2, force=False, verbose=True):
    """
    Build (or load, if cached and not force) the label-blind per-subtrack
    likelihood cache for one pooled fit.

    Parameters
    ----------
    pooled_traj_csv : str
        Path to a pooled trajectory CSV (combinationwise `*_pooled_traj.csv`
        or megapool `*_pooled_traj.csv`). The 'trajectory' column here is
        the RAW (pre-split) id that `fov_map_csv` is keyed on.
    settings_pkl : str
        Path to the deployed pipeline settings.pkl. The saSPT grid
        (diff_coefs, loc_errors, frame_interval, pixel_size_um,
        focal_depth, splitsize, likelihood_type) is inherited verbatim
        from here. `sample_size` is always overridden to 1e9.
    out_dir : str
        Directory to write `{cache_name}_likelihood_cache.npz` +
        `{cache_name}_fov_metadata.csv` sidecar.
    fov_map_csv : str, optional
        Path to the `*_fov_map.csv` sidecar for `pooled_traj_csv`. Must
        have a 'trajectory' column (raw id) and 'fov_basename' column,
        plus arbitrary additional grouping columns (e.g. 'source_member',
        'condition', 'zone'). Defaults to
        `{pooled_traj_csv minus _pooled_traj.csv}_fov_map.csv`.
    cache_name : str, optional
        Basename for the cache files. Defaults to the pooled_traj_csv's
        basename with `_pooled_traj.csv` stripped.
    min_track_length : int
        Subtracks below this length are dropped before caching (matches
        saSPT's own >1-jump filter; kept as an explicit, overridable arg
        rather than silently relying on saSPT's internal default).
    force : bool
        Recompute even if a cache already exists.
    verbose : bool

    Returns
    -------
    str
        Path to the written `.npz` cache.
    """
    os.makedirs(out_dir, exist_ok=True)

    if cache_name is None:
        base = os.path.basename(pooled_traj_csv)
        cache_name = base[:-len('_pooled_traj.csv')] if base.endswith(
            '_pooled_traj.csv') else os.path.splitext(base)[0]

    if fov_map_csv is None:
        pool_dir = os.path.dirname(pooled_traj_csv)
        fov_map_csv = os.path.join(pool_dir, f"{cache_name}_fov_map.csv")

    npz_path = os.path.join(out_dir, f"{cache_name}_likelihood_cache.npz")
    meta_path = os.path.join(out_dir, f"{cache_name}_fov_metadata.csv")

    if os.path.exists(npz_path) and not force:
        if verbose:
            print(f"[{cache_name}] cached -> {npz_path}")
        return npz_path

    if not os.path.exists(fov_map_csv):
        raise FileNotFoundError(
            f"fov_map_csv not found: {fov_map_csv} (pass fov_map_csv= "
            f"explicitly if it doesn't follow the default naming)")

    from saspt import StateArray

    saspt_settings = _load_saspt_settings(settings_pkl)
    grid_checksum = _grid_checksum(saspt_settings)
    grid_snapshot = {k: (saspt_settings[k].tolist()
                          if isinstance(saspt_settings[k], np.ndarray)
                          else saspt_settings[k])
                      for k in _GRID_KEYS}

    run_settings = copy.deepcopy(saspt_settings)
    run_settings['sample_size'] = FORCED_SAMPLE_SIZE
    run_settings['max_iter'] = 0     # label-blind: never touch VBEM here
    run_settings['progress_bar'] = False

    if verbose:
        print(f"[{cache_name}] loading {pooled_traj_csv}")
    traj_df = pd.read_csv(pooled_traj_csv)

    # Reproduce trackingCodeFunctions.executeSASPT's pre-fit filter exactly
    # (the actual call path for combinationwise/megapool saSPT runs) -- the
    # `counts >= 0` step there is a no-op and is skipped. Without this, the
    # cache would contain low-intensity trajectories the real production
    # fit excluded, breaking "cached subtrack set == fitted set".
    min_I0 = _load_min_I0(settings_pkl)
    if min_I0 is not None and 'I0' in traj_df.columns:
        n_before = traj_df['trajectory'].nunique()
        traj_df = traj_df.groupby('trajectory').filter(
            lambda g: g['I0'].mean() > min_I0)
        n_after = traj_df['trajectory'].nunique()
        if verbose:
            print(f"[{cache_name}] min_I0={min_I0} filter (matches "
                  f"executeSASPT): {n_before:,} -> {n_after:,} raw "
                  f"trajectories")
    elif min_I0 is not None:
        print(f"WARN: [{cache_name}] settings has quot.track.min_I0="
              f"{min_I0} but pooled_traj_csv has no 'I0' column -- "
              f"filter NOT applied")

    if verbose:
        print(f"[{cache_name}] running StateArray.from_detections "
              f"(sample_size={FORCED_SAMPLE_SIZE:g}, label-blind only)")
    import saspt
    SA = StateArray.from_detections(traj_df, **run_settings)

    naive = SA.naive_assignment_probabilities  # (n_D, n_loc_err, n_subtracks)
    jumps_per_track = SA.jumps_per_track        # (n_subtracks,)
    diff_coefs = np.asarray(SA.likelihood.diff_coefs)
    loc_errors = np.asarray(SA.likelihood.loc_errors)
    n_subtracks = naive.shape[-1]

    dets = SA.trajectories.detections
    sub_info = (dets.groupby('trajectory')
                     .agg(orig_trajectory=('orig_trajectory', 'first'),
                          track_length=('track_length', 'first'))
                     .sort_index())
    if not np.array_equal(sub_info.index.to_numpy(),
                           np.arange(n_subtracks)):
        raise RuntimeError(
            "subtrack index is not contiguous 0..n_subtracks-1; saspt "
            "internals may have changed -- re-run Gate 0")
    orig_trajectory = sub_info['orig_trajectory'].to_numpy()
    track_length = sub_info['track_length'].to_numpy()

    keep = track_length >= min_track_length
    if not keep.all():
        n_drop = int((~keep).sum())
        if verbose:
            print(f"[{cache_name}] dropping {n_drop} subtracks with "
                  f"track_length < {min_track_length}")
        naive = naive[:, :, keep]
        jumps_per_track = jumps_per_track[keep]
        orig_trajectory = orig_trajectory[keep]
        track_length = track_length[keep]
        n_subtracks = naive.shape[-1]

    fov_map_raw = pd.read_csv(fov_map_csv)
    if 'trajectory' not in fov_map_raw.columns or \
            'fov_basename' not in fov_map_raw.columns:
        raise ValueError(
            f"{fov_map_csv} must have 'trajectory' and 'fov_basename' "
            f"columns; found {list(fov_map_raw.columns)}")

    raw_to_fov = dict(zip(fov_map_raw['trajectory'].astype(int),
                          fov_map_raw['fov_basename']))
    fov_of_subtrack = pd.Series(orig_trajectory).map(raw_to_fov)
    n_unmapped = int(fov_of_subtrack.isna().sum())
    if n_unmapped > 0:
        frac = n_unmapped / n_subtracks
        msg = (f"[{cache_name}] {n_unmapped}/{n_subtracks} "
               f"({frac:.2%}) subtracks have no matching raw trajectory "
               f"in {fov_map_csv}")
        if frac > 0.01:
            raise RuntimeError(msg + " -- exceeds 1% tolerance, refusing "
                                "to build a cache with unreliable FOV "
                                "identity. Check that fov_map_csv matches "
                                "pooled_traj_csv exactly.")
        print("WARN: " + msg + " (dropping)")
        keep2 = fov_of_subtrack.notna().to_numpy()
        naive = naive[:, :, keep2]
        jumps_per_track = jumps_per_track[keep2]
        orig_trajectory = orig_trajectory[keep2]
        track_length = track_length[keep2]
        fov_of_subtrack = fov_of_subtrack[keep2]
        n_subtracks = naive.shape[-1]
    fov_of_subtrack = fov_of_subtrack.to_numpy().astype('U')

    # fov_metadata: one row per unique FOV, with whatever grouping columns
    # the input fov_map carries beyond 'trajectory'/'fov_basename' (e.g.
    # 'source_member', 'condition', 'zone'). Assert each FOV maps to a
    # single, consistent value per grouping column (a FOV belongs to
    # exactly one condition/member).
    group_cols = [c for c in fov_map_raw.columns
                  if c not in ('trajectory', 'fov_basename')]
    fov_metadata = fov_map_raw[['fov_basename'] + group_cols].drop_duplicates()
    dup = fov_metadata['fov_basename'].duplicated()
    if dup.any():
        bad = fov_metadata.loc[dup, 'fov_basename'].tolist()
        raise RuntimeError(
            f"fov_map_csv has inconsistent grouping-column values within "
            f"a single FOV for: {bad[:5]}{'...' if len(bad) > 5 else ''}")
    fov_metadata = fov_metadata.reset_index(drop=True)

    # Store L subtrack-major: (n_subtracks, n_D, n_loc_err). This is a
    # transpose of saspt's native (n_D, n_loc_err, n_subtracks) axis order
    # (recorded below) -- subtrack-major is what every downstream
    # comparison actually slices/permutes on.
    L = np.moveaxis(naive, -1, 0).astype(np.float32)

    np.savez_compressed(
        npz_path,
        L=L,
        L_axis_order=np.array(['subtrack', 'diff_coef', 'loc_error']),
        native_axis_order=np.array(['diff_coef', 'loc_error', 'subtrack']),
        jumps_per_track=jumps_per_track,
        track_length=track_length,
        orig_trajectory=orig_trajectory,
        fov_of_subtrack=fov_of_subtrack,
        diff_coefs=diff_coefs,
        loc_errors=loc_errors,
        grid_checksum=grid_checksum,
        grid_snapshot=json.dumps(grid_snapshot),
        saspt_version=saspt.__version__,
        pooled_traj_csv=os.path.abspath(pooled_traj_csv),
        fov_map_csv=os.path.abspath(fov_map_csv),
        settings_pkl=os.path.abspath(settings_pkl),
        gate0_precomputed_likelihood_injection='not supported in saspt '
            '0.4.0; StateArray always computes likelihood(trajectories) '
            'lazily -- see module docstring',
        gate0_naive_vs_posterior='naive_assignment_probabilities = '
            'likelihood.exp(log_L), per-subtrack normalized over the full '
            'grid, no VBEM; posterior_occs/posterior_assignment_'
            'probabilities require infer_posterior() (VBEM) and are not '
            'used by compare_groups',
    )
    fov_metadata.to_csv(meta_path, index=False)

    if verbose:
        print(f"[{cache_name}] wrote {npz_path}: "
              f"n_subtracks={n_subtracks:,}, n_FOV={len(fov_metadata)}, "
              f"grid={L.shape[1]}x{L.shape[2]}, checksum={grid_checksum[:12]}")
        print(f"[{cache_name}] wrote {meta_path}")

    return npz_path


def merge_likelihood_caches(member_caches, out_dir, cache_name,
                             group_column='source_member', member_metadata=None,
                             force=False, verbose=True):
    """
    Merge N independently-built per-condition caches (each from
    `build_likelihood_cache` on a SINGLE-condition pool) into one combined
    cache with a synthetic grouping column -- the alternative to building
    a brand-new JOINT saSPT megapool just to get two/more conditions into
    one comparator cache.

    Why this is valid: the naive per-subtrack likelihood (what these
    caches store) is provably a function of that subtrack's own jump
    vector and the shared grid/settings ONLY -- saspt's RBMELikelihood.
    __call__ computes each track's covariance-based likelihood
    independently (verified directly: it loops over trajectories/length
    classes and never references any OTHER trajectory). Concatenating N
    independently-fit caches therefore gives byte-identical `L` values to
    what ONE joint saSPT call over the union of their detections would
    have produced, as long as every member was built from the SAME
    deployed settings.pkl (checksummed below -- refuses to merge on any
    mismatch). No new saSPT execution needed; this is a pure numpy
    concatenation of already-cached arrays.

    Parameters
    ----------
    member_caches : dict
        {label: npz_path}. Each npz_path is a cache built by
        `build_likelihood_cache` from a single-condition (or any) pool.
        `label` becomes that member's value in the new `group_column`.
    out_dir, cache_name : where to write `{cache_name}_likelihood_cache.npz`
        + `{cache_name}_fov_metadata.csv`.
    group_column : name of the new grouping column added to fov_metadata
        (default 'source_member', matching the megapool convention this
        mirrors).
    member_metadata : dict, optional
        {label: {col: val, ...}} -- extra metadata columns beyond
        `group_column` (e.g. {'construct': '101', 'drug': 'HHT'}), so one
        merged cache spanning many conditions can serve many contrasts via
        compare_groups' filter-dict group resolution, instead of building
        a separate merged cache per pair. `group_column` is still set to
        `label` for every member regardless.
    force : bool
        Recompute even if a cache already exists at the target path.

    Returns
    -------
    str
        Path to the written merged `.npz` cache.
    """
    os.makedirs(out_dir, exist_ok=True)
    npz_path = os.path.join(out_dir, f"{cache_name}_likelihood_cache.npz")
    meta_path = os.path.join(out_dir, f"{cache_name}_fov_metadata.csv")
    if os.path.exists(npz_path) and not force:
        if verbose:
            print(f"[{cache_name}] cached -> {npz_path}")
        return npz_path

    if len(member_caches) < 2:
        raise ValueError("merge_likelihood_caches needs >= 2 member caches "
                          "(use build_likelihood_cache directly for a "
                          "single pool)")

    members = {}
    for label, path in member_caches.items():
        members[label] = load_likelihood_cache(path, assert_checksum=False)

    checksums = {label: c['grid_checksum'] for label, c in members.items()}
    unique_checksums = set(checksums.values())
    if len(unique_checksums) > 1:
        raise RuntimeError(
            f"grid checksum mismatch across member caches -- refusing to "
            f"merge (their per-subtrack likelihoods are NOT guaranteed "
            f"comparable): {checksums}")
    grid_checksum = unique_checksums.pop()

    L_parts, jumps_parts, tl_parts, orig_parts, fov_parts = [], [], [], [], []
    fov_metadata_parts = []
    orig_offset = 0
    for label, c in members.items():
        n_sub = c['L'].shape[0]
        L_parts.append(c['L'])
        jumps_parts.append(c['jumps_per_track'])
        tl_parts.append(c['track_length'])
        orig_parts.append(c['orig_trajectory'].astype(np.int64) + orig_offset)
        fov_parts.append(c['fov_of_subtrack'])
        orig_offset += int(c['orig_trajectory'].max()) + 1 if n_sub > 0 else 0

        meta = c['fov_metadata'][['fov_basename']].drop_duplicates().copy()
        meta[group_column] = label
        for col, val in (member_metadata or {}).get(label, {}).items():
            meta[col] = val
        fov_metadata_parts.append(meta)

        if verbose:
            print(f"  [{label}] {n_sub:,} subtracks, "
                  f"{c['fov_metadata']['fov_basename'].nunique()} FOVs")

    fov_metadata = pd.concat(fov_metadata_parts, ignore_index=True)
    dup = fov_metadata['fov_basename'].duplicated()
    if dup.any():
        bad = fov_metadata.loc[dup, 'fov_basename'].tolist()
        raise RuntimeError(
            f"the same FOV basename appears in more than one member cache "
            f"-- real collision, not expected across distinct conditions: "
            f"{bad[:5]}{'...' if len(bad) > 5 else ''}")

    diff_coefs = next(iter(members.values()))['diff_coefs']
    loc_errors = next(iter(members.values()))['loc_errors']
    grid_snapshot = next(iter(members.values()))['grid_snapshot']

    L = np.concatenate(L_parts, axis=0)
    jumps_per_track = np.concatenate(jumps_parts, axis=0)
    track_length = np.concatenate(tl_parts, axis=0)
    orig_trajectory = np.concatenate(orig_parts, axis=0)
    fov_of_subtrack = np.concatenate(fov_parts, axis=0).astype('U')

    np.savez_compressed(
        npz_path,
        L=L,
        L_axis_order=np.array(['subtrack', 'diff_coef', 'loc_error']),
        native_axis_order=np.array(['diff_coef', 'loc_error', 'subtrack']),
        jumps_per_track=jumps_per_track,
        track_length=track_length,
        orig_trajectory=orig_trajectory,
        fov_of_subtrack=fov_of_subtrack,
        diff_coefs=diff_coefs,
        loc_errors=loc_errors,
        grid_checksum=grid_checksum,
        grid_snapshot=grid_snapshot,
        merged_from=json.dumps({label: os.path.abspath(p)
                                 for label, p in member_caches.items()}),
    )
    fov_metadata.to_csv(meta_path, index=False)

    if verbose:
        print(f"[{cache_name}] wrote {npz_path}: "
              f"n_subtracks={L.shape[0]:,}, n_FOV={len(fov_metadata)}, "
              f"grid={L.shape[1]}x{L.shape[2]}, checksum={grid_checksum[:12]} "
              f"(merged from {len(member_caches)} member caches, no new "
              f"saSPT fit)")

    return npz_path


def load_likelihood_cache(npz_path, settings_pkl=None, assert_checksum=True):
    """
    Load a likelihood cache built by `build_likelihood_cache`.

    Parameters
    ----------
    npz_path : str
    settings_pkl : str, optional
        If given, re-derive the grid checksum from this deployed
        settings.pkl and assert it matches the cache's stored checksum --
        refuses to load on mismatch (handoff invariant: "assert match on
        load").
    assert_checksum : bool

    Returns
    -------
    dict with keys: L (n_subtracks, n_D, n_loc_err), jumps_per_track,
        track_length, orig_trajectory, fov_of_subtrack, diff_coefs,
        loc_errors, grid_checksum, fov_metadata (pd.DataFrame),
        plus the raw npz for anything else.
    """
    meta_path = npz_path.replace('_likelihood_cache.npz', '_fov_metadata.csv')
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"missing fov_metadata sidecar: {meta_path}")

    npz = np.load(npz_path, allow_pickle=True)  # trusted, self-produced cache
    out = {k: npz[k] for k in npz.files}
    out['grid_checksum'] = str(out['grid_checksum'])
    out['fov_metadata'] = pd.read_csv(meta_path)

    if settings_pkl is not None and assert_checksum:
        live_settings = _load_saspt_settings(settings_pkl)
        live_checksum = _grid_checksum(live_settings)
        if live_checksum != out['grid_checksum']:
            raise RuntimeError(
                f"Grid checksum mismatch: cache={out['grid_checksum'][:12]} "
                f"vs deployed settings.pkl={live_checksum[:12]}. The "
                f"likelihood cache was built on a different saSPT grid "
                f"than the currently deployed settings -- rebuild it "
                f"(force=True) before trusting any comparison drawn from "
                f"it.")

    return out
