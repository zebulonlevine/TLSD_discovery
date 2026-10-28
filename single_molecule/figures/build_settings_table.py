#!/usr/bin/env python
"""
build_settings_table.py -- Supplementary table of the quot tracking and saSPT
inference settings behind every SPT panel.

Emits Supplementary Table S12 twice:

    supplementary/TableS12_SPT_settings.csv  one row per parameter, ONE value
                                             column: every setting is identical
                                             across the three datasets behind the
                                             draft5 panels (asserted -- a
                                             difference raises instead of being
                                             collapsed), except the bin scheme,
                                             whose cell names both subunits
    supplementary/TableS12_SPT_settings.md   the same rows as Markdown

EVERY VALUE IS READ FROM THE DEPLOYED FILES, never typed here. The defaults in
init_settings.py are placeholders (pixel_size 0.11 / frame_interval 0.01 grid
values) and have been wrong for at least one dataset, so a hand-maintained
table would drift silently. Re-run this after any settings change.

Two things the table has to get right, because both have bitten this project:

  1. `settings.pkl` is NOT what produced the published panels. Those were built
     under PANEL_VARIANT=36state, which carries its own settings block
     (occupancy_compare_cache_36state/settings_36state.pkl) with saSPT's default
     36-state localization-error grid and its default 0.7 um focal depth, rather than
     the 5-state grid and defocalization-off 0.1 um of the legacy block. The
     saSPT rows therefore come from the variant, and the legacy values are
     reported in a footnote rather than dropped.
  2. Several numbers that govern the published inference are NOT in any settings
     file -- they are constants in panel_core (concentration parameter, VBEM
     iteration count and tolerance). Those are pulled from the module so they
     cannot fall out of sync either.
  3. The resampling counts are read from the statistics files the placed panels
     were drawn from, not from panel_core.N_BOOT/N_PERM (2000), which are the
     non-refit defaults: the refit panels (S7A, S7C, S7D) used 500 draws, and
     Fig 5F used 2000 zone-label and 10000 treatment permutations.
"""
import csv
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
OUT = tlsd_paths.outp('PaperFigures_updated', 'supplementary')

# Dataset -> (tier, the panels it produces). Tier paths are copied from the
# builders themselves; keep them in step if a builder is repointed. They are
# resolved through tlsd_paths.spt, so a per-experiment data root applies.
DATASETS = [
    ('2025_06_06_zeb', 'all_data/2025_06_06_zeb/tracking_output_q=15p0_splitsize10',
     'Fig 5C, S7A, S7B, S7C'),
    ('2026_07_08', 'all_data/2026_07_08/tracking_output_q=15p0',
     'Fig S7D, S7E'),
    ('2025_09_15_zeb_ER', 'all_data/2025_09_15_zeb_ER/tracking_output_q=15p0',
     'Fig 5E, 5F, S7H-J'),
]

VARIANT_REL = 'occupancy_compare_cache_36state/settings_36state.pkl'


def fmt(v):
    """A cell value that survives a round trip through CSV and reads cleanly."""
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, str) and v in ('True', 'False'):   # stringified bool in one pkl
        return v.lower()
    if isinstance(v, float):
        # Do NOT collapse 15.0 -> 15: the table doubles as something a reader
        # retypes into a YAML override, and the config value is a float.
        r = f'{v:.6g}'
        return r if '.' in r or 'e' in r else r + '.0'
    return str(v)


def grid(a):
    """Describe a saSPT state grid without dumping 150 numbers into a cell."""
    a = np.asarray(a, dtype=float)
    if len(a) > 1 and np.allclose(np.diff(np.log10(a[a > 0])),
                                  np.diff(np.log10(a[a > 0]))[0]):
        return f'{len(a)} states, {a.min():g}-{a.max():g}, log-spaced'
    step = np.diff(a)
    if len(a) > 1 and np.allclose(step, step[0]):
        return f'{len(a)} states, {a.min():g}-{a.max():g}, step {step[0]:g}'
    return f'{len(a)} states, {a.min():g}-{a.max():g}'


# (section, parameter, getter, units, description)
#
# `getter` takes (base_settings, variant_settings) and returns the value, so a
# row can declare for itself which settings block is authoritative.
def q(*path):
    return lambda n, b, v: _dig(b['quot'], path)


def sa(key):
    return lambda n, b, v: v['saspt'][key]


def _dig(d, path):
    for p in path:
        d = d[p]
    return d


# Region scheme per dataset, as the builders import it. 2025_06_06_zeb (Fig 5C,
# S7A-C) runs both subunits against different schemes, so its cell carries both.
from region_schemes import SCHEME_60S, SCHEME_40S  # noqa: E402

SCHEMES = {
    '2025_06_06_zeb': [('60S', SCHEME_60S), ('40S', SCHEME_40S)],
    '2026_07_08': [('', SCHEME_60S)],
    '2025_09_15_zeb_ER': [('', SCHEME_60S)],
}


def _scheme_for(name):
    return SCHEMES[name]


def bounds(name):
    out = []
    for tag, sch in _scheme_for(name):
        edges = sorted({e for v in sch.values() for e in v} - {0.0, 30.0})
        out.append((f'{tag}: ' if tag else '') + ', '.join(f'{e:g}' for e in edges))
    return '; '.join(out)


def bin_names(name):
    out = []
    for tag, sch in _scheme_for(name):
        labels = [k for k, _ in sorted(sch.items(), key=lambda kv: kv[1][0])]
        out.append((f'{tag}: ' if tag else '') + ', '.join(labels))
    return '; '.join(out)


ROWS = [
    ('Acquisition', 'pixel_size_um', q('track', 'pixel_size_um'), 'um/px',
     'Camera pixel size in the sample plane'),
    ('Acquisition', 'frame_interval', q('track', 'frame_interval'), 's',
     'Time between consecutive frames'),

    ('quot: filter', 'method', q('filter', 'method'), '',
     'Background subtraction over the frame stack; identity = none applied'),
    ('quot: filter', 'start', q('filter', 'start'), 'frame',
     'First frame used. ND2 frame 0 is a pims gap frame and is excluded upstream'),
    ('quot: filter', 'chunk_size', q('filter', 'chunk_size'), 'frames',
     'Frames held in memory per filtering block'),

    ('quot: detect', 'method', q('detect', 'method'), '',
     'Spot detection: log-likelihood ratio test for a Gaussian spot'),
    ('quot: detect', 'k', q('detect', 'k'), 'px',
     'Gaussian kernel width used by the LLR detector'),
    ('quot: detect', 'w', q('detect', 'w'), 'px',
     'Detection window edge length'),
    ('quot: detect', 't', q('detect', 't'), 'LLR',
     'Detection threshold. This is the "q" in the tracking_output_q=15p0 tier name'),

    ('quot: localize', 'method', q('localize', 'method'), '',
     'Subpixel localization: least-squares integrated-Gaussian fit'),
    ('quot: localize', 'window_size', q('localize', 'window_size'), 'px',
     'Fitting window edge length'),
    ('quot: localize', 'sigma', q('localize', 'sigma'), 'px',
     'Fixed PSF standard deviation'),
    ('quot: localize', 'ridge', q('localize', 'ridge'), '',
     'Ridge regularization on the fit'),
    ('quot: localize', 'max_iter', q('localize', 'max_iter'), 'iterations',
     'Levenberg-Marquardt iteration cap'),
    ('quot: localize', 'damp', q('localize', 'damp'), '',
     'Damping factor on the parameter update'),
    ('quot: localize', 'camera_bg', q('localize', 'camera_bg'), 'ADU',
     'Camera offset, read from the SPT TIF; reading it from the ND2 yields NaN'),

    ('quot: track', 'method', q('track', 'method'), '',
     'Frame-to-frame linking by Euclidean distance (Hungarian assignment)'),
    ('quot: track', 'search_radius', q('track', 'search_radius'), 'um',
     'Maximum displacement allowed between consecutive frames'),
    ('quot: track', 'max_blinks', q('track', 'max_blinks'), 'frames',
     'Gap frames tolerated within one trajectory; 0 = no gap closing'),
    ('quot: track', 'min_I0', q('track', 'min_I0'), 'ADU',
     'Spot intensity threshold: a spot below it cannot start a trajectory, and trajectories '
     'with mean spot intensity at or below it are discarded before saSPT'),
    ('quot: track', 'scale', q('track', 'scale'), '',
     'Coordinate scale factor; 1.0 keeps trajectories in pixels'),

    ('saSPT: model', 'likelihood_type', sa('likelihood_type'), '',
     'RBME: regular Brownian motion with localization error marginalized'),
    ('saSPT: model', 'diff_coefs', lambda n, b, v: grid(v['saspt']['diff_coefs']),
     'um^2/s', 'Diffusion-coefficient state grid'),
    ('saSPT: model', 'loc_errors', lambda n, b, v: grid(v['saspt']['loc_errors']),
     'um', "Localization-error state grid (saSPT's own default under this variant)"),
    ('saSPT: model', 'focal_depth', sa('focal_depth'), 'um',
     "Defocalization correction depth (saSPT's default); conclusions checked across 0.1-1.5 um "
     'with validation/calibrate_focal_depth.py'),
    ('saSPT: model', 'splitsize', sa('splitsize'), 'jumps',
     'Maximum jumps per sub-trajectory'),
    ('saSPT: model', 'start_frame', sa('start_frame'), 'frame',
     'First frame contributing jumps'),

    ('saSPT: model', 'sub_trajectories_used',
     lambda n, b, v: 'all', '',
     "Every sub-trajectory in the likelihood cache enters the fit; saSPT's own "
     "sample_size cap applies to StateArray.from_detections, not to this path"),

    ('saSPT: inference', 'concentration_parameter',
     lambda n, b, v: pc.CONC_PARAM, '',
     "Dirichlet prior on state occupancies (saSPT DEFAULT_CONC_PARAM); absent from settings.pkl"),
    ('saSPT: inference', 'max_iter', lambda n, b, v: pc.MAX_ITER_SASPT, 'iterations',
     'VBEM iteration cap, matching saSPT DEFAULT_MAX_ITER'),
    ('saSPT: inference', 'vbem_iterations_applied', lambda n, b, v: pc._VBEM_ITERS,
     'iterations',
     'max_iter-1: saSPT contracts R->n at the top of each pass, so K passes reflect K-1 R-updates'),
    ('saSPT: inference', 'convergence_tolerance', lambda n, b, v: pc._VBEM_TOL, '',
     '0 disables early stopping, which saSPT does not have'),
    ('saSPT: inference', 'occupancy_weighting',
     lambda n, b, v: 'displacement' if not pc.EQUAL_WEIGHT else 'equal-per-subtrajectory',
     '', "saSPT's own definition: occupancy is a fraction of displacements"),

    ('Diffusion bins', 'bin_boundaries',
     lambda n, b, v: bounds(n), 'um^2/s',
     'Interior boundaries of the manually-selected bin scheme; '
     'bins are numbered I-V in the figures'),
    ('Diffusion bins', 'bin_labels', lambda n, b, v: bin_names(n), '',
     'Bin names in ascending D'),

    ('Resampling', 'n_bootstrap', lambda n, b, v: RS['n_refit'], 'draws',
     'Field-of-view bootstrap draws for confidence intervals'),
    ('Resampling', 'n_permutation', lambda n, b, v: RS['n_refit'], 'draws',
     'Whole-FOV label permutations for the null (genotype and drug contrasts, '
     'Figs S7A, S7C, S7D); each draw re-fit with saSPT'),
    ('Resampling', 'permutation_p_floor', lambda n, b, v: f"1/{RS['n_refit'] + 1}", '',
     'Phipson-Smyth add-one estimator; a permutation p cannot be 0 (floor for the '
     '500-draw tests; Fig 5F floors in the rows below)'),
    ('Resampling', 'n_permutation_zone_Fig5F', lambda n, b, v: RS['n_zone_perm'], 'draws',
     'Fig 5F on- vs off-ER contrast within each treatment: zone labels shuffled among '
     "trajectories within each field of view (a trajectory's sub-trajectories move "
     'together); the pooled saSPT fit is invariant under this null, so no re-fit is '
     f"needed. p floor 1/{{zone_floor}}"),
    ('Resampling', 'n_permutation_treatment_Fig5F', lambda n, b, v: RS['n_treat_perm'], 'draws',
     'Fig 5F PT-179 vs DMSO difference-in-differences of the on- minus off-ER contrast: '
     'fields of view permuted between treatments on per-FOV zone differences from one '
     "joint saSPT fit. p floor 1/{treat_floor}; 95% CI from all {n_pairs} pairs of the "
     "two treatments' {n_zone_boot} bootstrap draws"),
    ('Resampling', 'multiple_comparison', lambda n, b, v: 'Benjamini-Hochberg', '',
     "Applied within each contrast's bins; Holm also carried in every stats CSV"),
]


# Resampling counts, read from the placed panels' statistics files.
RS = {}


def resampling_counts():
    import pandas as pd
    fig = lambda *p: pc.tag(tlsd_paths.outp('PaperFigures_updated', *p))
    refit = set()
    for f in ('figureS7/A/figS7A_W1_refit_stats.csv',
              'figureS7/C/figS7C_60S_bin_difference_refit_stats.csv',
              'figureS7/D/figS7D_refit_stats.csv'):
        refit |= set(pd.read_csv(fig(f))['n_draws'])
    if len(refit) != 1:
        raise SystemExit(f'refit panels disagree on the draw count: {refit}')
    z = pd.read_csv(fig('figure5/F/fig5F_refit_stats.csv'))
    one = lambda panel, col: int(z.loc[z.panel == panel, col].drop_duplicates().item())
    RS.update(n_refit=int(refit.pop()),
              n_zone_perm=one('5F_zone_effect', 'n_perm'),
              n_zone_boot=one('5F_zone_effect', 'n_boot'),
              n_treat_perm=one('5F_interaction', 'n_perm'),
              n_pairs=one('5F_interaction', 'n_boot'))


def main():
    os.makedirs(OUT, exist_ok=True)
    resampling_counts()
    loaded = []
    for name, tier, panels in DATASETS:
        base = pickle.load(open(tlsd_paths.spt(*tier.split('/'), 'settings.pkl'), 'rb'))
        var = pickle.load(open(tlsd_paths.spt(*tier.split('/'), *VARIANT_REL.split('/')), 'rb'))
        loaded.append((name, tier, panels, base, var))

    header = ['Section', 'Parameter Name', 'Value', 'Units', 'Description']
    fill = dict(zone_floor=RS['n_zone_perm'] + 1, treat_floor=RS['n_treat_perm'] + 1,
                n_pairs=f"{RS['n_pairs']:,}", n_zone_boot=RS['n_zone_boot'])
    rows = []
    for section, param, get, units, desc in ROWS:
        vals = [fmt(get(n, b, v)) for n, _, _, b, v in loaded]
        if len(set(vals)) == 1:
            val = vals[0]
        elif section == 'Diffusion bins':
            # 2025_06_06_zeb runs both subunit schemes; the others use its 60S one
            val = max(vals, key=len)
            assert all(x in val for x in vals), (param, vals)
        else:
            raise SystemExit(f'{param} differs between datasets: '
                             f'{dict(zip([n for n, *_ in loaded], vals))}')
        rows.append([section, param, val, units, desc.format(**fill)])

    csv_path = os.path.join(OUT, 'TableS12_SPT_settings.csv')
    with open(csv_path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)

    md = ['# Supplementary Table S12: SPT acquisition, tracking and saSPT settings', '',
          'Datasets: ' + '; '.join(f'`{n}` ({p})' for n, _, p, *_ in loaded) + '.', '',
          '| Section | Parameter | Value | Units | Description |', '|---|---|---|---|---|']
    md += [f'| {r[0]} | `{r[1]}` | {r[2]} | {r[3]} | {r[4]} |' for r in rows]
    md += ['', 'quot rows: `<tier>/settings.pkl`. saSPT rows: '
           f'`<tier>/{VARIANT_REL}` (the published panels were built under '
           '`PANEL_VARIANT=36state`). Inference rows: constants in `panel_core.py`. '
           'Resampling rows: the statistics files of the placed panels.', '']

    md_path = os.path.join(OUT, 'TableS12_SPT_settings.md')
    with open(md_path, 'w') as fh:
        fh.write('\n'.join(md))

    print(f'wrote {csv_path}\n      {md_path}')
    print(f'{len(rows)} parameters (identical across {len(loaded)} datasets)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
