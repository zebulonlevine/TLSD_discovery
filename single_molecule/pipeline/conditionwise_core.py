## conditionwise_core.py
##
## Shared per-condition (and per-combination) downstream processing.
##
## After file-wise saSPT outputs (posterior, MLE, rollingMLE) are written to
## settings['io']['post_directory'], settings['io']['MLE_directory'], and
## settings['io']['rolling_window_directory'], this function aggregates them
## across files for a single condition (or pooled combination), determines
## diffusion-coefficient cutoffs, generates plots, and produces movies and
## bleaching curves.
##
## HMM and ExTrack post-processing have been removed from this pipeline.
##
## Callers:
##   - run_conditionwise.py        (per-condition; dirs = standard)
##   - run_combinationwise.py      (planned; redirects io dirs to *_combinations/)


import os
import re
import pickle
import numpy as np
import pandas as pd

import trackingCodeFunctions
import plotFunctions


def process_condition(cond, settings):
    """
    Run all per-condition downstream analysis for a given condition.

    Aggregates posterior, MLE, and rollingMLE CSVs from the directories
    specified in settings['io'], removes outlier files, determines bar-graph
    cutoffs (group-aware GMM with caching), generates plots, saves the
    plotting PKL, computes bleaching curves, and generates example movies.

    For a combination (pooled) run, the caller should deepcopy settings,
    redirect the three input directories to the *_combinations/ variants,
    and pass the combination name as cond. Plot output naturally lands in
    settings['io']['plot_directory'] labeled with cond.

    Parameters
    ----------
    cond : str
        Condition name (or combination name when called from
        run_combinationwise).
    settings : dict
        Loaded settings dictionary. Should already have any directory
        redirection applied by the caller. This function does not mutate
        settings beyond writing back the resolved barGraphBreaks /
        barGraphLabels at the end of cutoff determination.
    """

    # ============================================================
    # Aggregate per-file outputs for this condition
    # ============================================================
    # For combinations, the pooled CSVs are named exactly
    # {combo}_pooled_{posterior,MLE,rollingMLE}.csv. A bare {cond}*
    # prefix glob would over-match when one combo name is a prefix of
    # another (e.g., 'oeR3HDM4' globs both oeR3HDM4_pooled_* AND
    # oeR3HDM4_sgUBA5_pooled_*). Use the explicit '_pooled' infix for
    # combos to make matches unambiguous.
    combo_ctx_check = settings.get('_combination_context')
    if combo_ctx_check is not None:
        agg_pat = f"{cond}_pooled_*.csv"
    else:
        # Per-FOV per-condition files follow {cond}_{NNN}_{type}.csv.
        # Anchor on the underscore to avoid the same prefix-collision
        # issue between conditions (e.g., 'cond' vs 'cond_v2').
        agg_pat = f"{cond}_*.csv"

    posterior_df = trackingCodeFunctions.aggregate_csv(
        settings['io']['post_directory'], agg_pat)
    MLE_df = trackingCodeFunctions.aggregate_csv(
        settings['io']['MLE_directory'], agg_pat)
    # rolling_window_directory may not exist in megapool settings (rolling
    # is skipped). The aggregate may also find zero matching files, which
    # makes aggregate_csv raise on its internal pd.concat([]). Treat both
    # as "no rolling data" and proceed.
    rolling_dir = settings['io'].get('rolling_window_directory')
    if rolling_dir is None or not os.path.isdir(rolling_dir):
        rollingMLE = pd.DataFrame()
    else:
        try:
            rollingMLE = trackingCodeFunctions.aggregate_csv(
                rolling_dir, agg_pat)
        except ValueError:
            # 'No objects to concatenate' — empty glob
            rollingMLE = pd.DataFrame()

    # When called from run_megapool, rolling-window analysis is skipped
    # entirely and the rollingMLE directory contains no matching files.
    # Detect that here and proceed without rolling-window-derived output.
    has_rolling = (rollingMLE is not None
                   and len(rollingMLE) > 0
                   and 'ur_trajectory' in rollingMLE.columns)

    # Stamp condition column directly. The aggregate_csv glob filtered on
    # `{cond}*.csv`, so every row belongs to this condition (or to the
    # pooled combination whose pooled CSV starts with cond). Skipping the
    # filename-regex extraction makes this robust to combination files
    # that don't follow the per-FOV `_NNN_traj.csv` naming convention.
    posterior_df['condition'] = cond

    if 'source_file' in MLE_df.columns:
        MLE_df['file_basename'] = [os.path.basename(f) for f in MLE_df['source_file']]
    elif 'orig_trajectory' in MLE_df.columns:
        MLE_df['file_basename'] = MLE_df['orig_trajectory'].str.split('::').str[1]
    else:
        MLE_df['file_basename'] = ''
    MLE_df['condition'] = cond

    if has_rolling:
        rollingMLE['file_basename'] = [
            re.sub("::[0-9]*$", '.csv', traj) for traj in rollingMLE['ur_trajectory']]
        rollingMLE['condition'] = cond


    # ============================================================
    # Outlier removal (per-FOV; no-op when only one file is present,
    # which is the case for pooled combinations)
    # ============================================================
    outliers = plotFunctions.find_outlier_filenames(
        posterior_df, mult_on_sd=settings['plot']['mult_on_sd'])

    print('')
    print(f"allowed mult_on_sd is {settings['plot']['mult_on_sd']}")
    if outliers:
        print('Removing outliers:')
        for ol in outliers:
            print(f"    {ol}")
    else:
        print('No outliers removed.')
    print('')

    posterior_df = posterior_df[~posterior_df['file'].isin(outliers)]
    MLE_df = MLE_df[~MLE_df['file_basename'].isin(outliers)]
    if has_rolling:
        rollingMLE = rollingMLE[~rollingMLE['file_basename'].isin(outliers)]


    # ============================================================
    # Cutoffs: manual only (auto_gmm/hierarchical_gmm modes and
    # cutoff_functions.py removed 2026-07-12 -- no scientific conclusion
    # depended on them; the actual cutoff methodology in use is
    # cutoff_floor_analysis.py's agreed_cutoffs(), a separate notebook-
    # side, cutoff-free/permutation-floor-based method). barGraphBreaks/
    # barGraphLabels pass through unchanged from settings['plot']
    # (defaults, or per-condition/combo/megapool overrides already
    # applied upstream).
    # ============================================================
    megapool_ctx = settings.get('_megapool_context')
    combo_ctx = settings.get('_combination_context')
    if megapool_ctx is not None:
        print(f"\nProcessing megapool '{megapool_ctx['megapool_name']}' "
              f"member '{megapool_ctx['member']}'")
    elif combo_ctx is not None:
        print(f"\nProcessing combination '{combo_ctx['combo_name']}'")

    gmm_result = None
    print(f"Using barGraphBreaks: {settings['plot']['barGraphBreaks']}")
    print(f"Using barGraphLabels: {settings['plot']['barGraphLabels']}")


    # ============================================================
    # Plots
    # ============================================================
    posterior_plot_df = plotFunctions.generate_posterior_plotting_df(posterior_df)
    # Override condition label so combination plots are labeled by combo name
    posterior_plot_df['condition'] = cond

    plotFunctions.generate_plots(
        MLE_df, posterior_plot_df, posterior_df, settings, showPlot=False)


    # ============================================================
    # Save plotting PKL (no HMM-derived columns)
    # ============================================================
    plotting_pkl_dir = f"{settings['io']['plot_directory']}/plotting_pkls"
    os.makedirs(plotting_pkl_dir, exist_ok=True)

    plotting_pkl = f"{plotting_pkl_dir}/{cond}.pkl"
    with open(plotting_pkl, "wb") as f:
        pickle.dump({
            'MLE_df': MLE_df,
            'posterior_plot_df': posterior_plot_df,
            'rollingMLE': rollingMLE,
            'barGraphBreaks': settings['plot']['barGraphBreaks'],
            'barGraphLabels': settings['plot']['barGraphLabels'],
            'gmm_result': gmm_result,
        }, f)


    # ============================================================
    # Bleaching curves and example movies
    # (skipped for megapools, which don't run rolling-window analysis)
    # ============================================================
    if has_rolling:
        bc = plotFunctions.aggregate_bleaching_data(
            rollingMLE, settings, n_particles=250)
        plotFunctions.plot_bleaching_curves(settings, bc)
        plotFunctions.generate_movies(rollingMLE, settings)
    else:
        print('  No rollingMLE data; skipping bleaching curves and movies.')

    print(f'\nConditionwise analysis complete for {cond}')
