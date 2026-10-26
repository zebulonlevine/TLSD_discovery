"""
panel_core.py -- the shared numeric + plotting core for the updated SPT panels.

Every panel in `PaperFigures_updated/` is built through this module so that all
of them share one fit configuration, one resampling scheme, one multiplicity
correction, and one visual style. Nothing here re-implements saSPT: the fit is
`occupancy_compare._vbem_infer` (the verified port of
`saspt.sa.StateArray.infer_posterior`) run on likelihoods produced by a real
`saspt.StateArray`, and the defocalization correction is saSPT's own
`Likelihood.correct_for_defocalization`.

Fit configuration (the "exactly as saSPT" settings)
---------------------------------------------------
`EQUAL_WEIGHT = False` -- each sub-trajectory is weighted by its number of
displacements, which is saSPT's own definition of occupancy
(`saspt/sa.py:356`: `n = (R * jumps_per_track).sum(axis=-1)`). The earlier
figure build used equal weighting; that is not what saSPT does.

Defocalization is applied exactly as saSPT applies it. Note a structural
consequence, verified bitwise by `test_saspt_equivalence.py` T2(iii): saSPT
applies the correction AFTER the VBEM loop, to the Dirichlet sufficient
statistic only, so the converged `exp_log_tau` -- and therefore every
per-sub-trajectory posterior built from it -- is independent of `focal_depth`.
Per-FOV soft occupancy, and every bin panel derived from it, is therefore
*structurally immune* to the focal-depth choice. Only the population posterior
curve and W1 move with `focal_depth`, which is why those are the quantities the
dz sweep in `calibrate_focal_depth.py` covers.

Statistics
----------
- Per-FOV soft occupancy is the unit of analysis; the FOV is the exchangeable
  unit (sub-trajectory-level resampling is anticonservative because
  sub-trajectories cluster within a cell/FOV).
- Uncertainty: FOV cluster bootstrap, 95% percentile interval.
- Null: whole-FOV relabelling between the two groups, preserving group sizes.
- Permutation p uses the add-one estimator (1+k)/(1+B).
- Multiplicity: Benjamini-Hochberg across the bins of a panel, matching the
  convention used throughout the rest of the paper (the previous build used
  Holm; Holm is also reported in the stats CSV as the conservative alternative).
"""

import json
import os
import re
import sys

import numpy as np
import pandas as pd

_WORKFLOW = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pipeline'))
if _WORKFLOW not in sys.path:
    sys.path.insert(0, _WORKFLOW)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from occupancy_compare import (_make_likelihood_from_cache, _resolve_group,
                               _subband_summary, _w1_logD, _vbem_infer)
from saspt_likelihood_cache import load_likelihood_cache, _load_saspt_settings
from softpost_region_panels import group_soft_per_fov
from multiplicity import bh, holm, perm_p_addone
from saspt.defoc import f_remain_rbm

import paperfig_style as ps

# saSPT's own weighting: occupancy is a fraction of DISPLACEMENTS. Set
# PANEL_EQUAL_WEIGHT=1 in the environment to rebuild everything under the
# equal-per-sub-trajectory weighting the previous figure build used, which is
# what `compare_weighting.py` diffs against. Outputs are suffixed so the two
# builds never overwrite each other.
EQUAL_WEIGHT = os.environ.get('PANEL_EQUAL_WEIGHT', '') == '1'

# PANEL_VARIANT selects an alternative likelihood-cache layout and settings
# block, so a re-fit under different inference parameters lands beside the
# committed panels instead of replacing them. Currently one variant:
#
#   PANEL_VARIANT=36state   saSPT's default 36-state localization-error grid
#                           (0-0.070 um) and focal_depth 0.7 um, with caches
#                           built PER CONDITION rather than merged -- at 36
#                           states a merged claim-2 cache would be ~25 GB
#                           stored and ~50 GB during saSPT's float64 build.
#
# The variant also carries its own settings.pkl, written next to its caches, so
# load_likelihood_cache's grid-checksum assertion has the matching reference
# and cannot silently pair a cache with the wrong grid.
VARIANT = os.environ.get('PANEL_VARIANT', '')
CACHE_SUBDIR = ('occupancy_compare_cache'
                + (f'_{VARIANT}' if VARIANT else ''))
SETTINGS_REL = (os.path.join(CACHE_SUBDIR, f'settings_{VARIANT}.pkl')
                if VARIANT else 'settings.pkl')

SUFFIX = ('_equalweight' if EQUAL_WEIGHT else '') + (f'_{VARIANT}' if VARIANT else '')
CONC_PARAM = 1.0        # saspt DEFAULT_CONC_PARAM; settings.pkl omits it
MAX_ITER_SASPT = 200    # saspt DEFAULT_MAX_ITER; settings.pkl omits it

# saspt contracts R -> n at the TOP of each pass and updates R at the bottom,
# so the `n` it returns after K passes reflects K-1 R-updates. `_vbem_infer`
# folds both into a single recurrence, so after K passes its `n` reflects K
# R-updates. Calling it with K-1 therefore reproduces saspt(max_iter=K) to
# machine precision -- measured at 3.0e-16 by test_saspt_equivalence.py T1,
# versus 1.2e-5 if both are called with the same K. `tol=0` disables the
# convergence early-stop saspt does not have.
_VBEM_ITERS = MAX_ITER_SASPT - 1
_VBEM_TOL = 0.0

N_BOOT = 2000
N_PERM = 2000


# ---------------------------------------------------------------------------
# loading / fitting
# ---------------------------------------------------------------------------

def settings_path(tier):
    """The settings.pkl that matches the active variant's caches."""
    return os.path.join(tier, SETTINGS_REL)


def load(tier, cache_rel):
    """Load a likelihood cache with its grid checksum asserted against the
    settings block that matches it, plus the likelihood object and the D grid.

    `cache_rel` may be given relative to the tier with any cache subdirectory;
    under an active PANEL_VARIANT the subdirectory is rewritten to the
    variant's, so callers can keep passing one path.
    """
    spkl = settings_path(tier)
    base = os.path.basename(cache_rel)
    cache = load_likelihood_cache(os.path.join(tier, CACHE_SUBDIR, base),
                                  settings_pkl=spkl)
    return dict(cache=cache,
                likelihood=_make_likelihood_from_cache(cache),
                settings=_load_saspt_settings(spkl),
                diff_coefs=np.asarray(cache['diff_coefs']))


def condition_ctx(tier, condition, merged_cache):
    """Return (ctx, group) for one condition, hiding the cache-layout
    difference between variants.

    Default layout: one merged multi-condition cache, and the condition is
    selected with a `source_member` metadata filter.
    Variant layout: one cache per condition, so the group is simply every FOV
    in that cache.

    Returning the group alongside the context keeps the builders on a single
    code path instead of branching on the layout.
    """
    if VARIANT:
        ctx = load(tier, f'{condition}_likelihood_cache.npz')
        return ctx, list(ctx['cache']['fov_metadata']['fov_basename'])
    ctx = load(tier, merged_cache)
    return ctx, {'source_member': condition}


def _fit_idx(ctx, sub_idx):
    """One VBEM fit over explicit sub-trajectory ROW INDICES, at saSPT's exact
    iteration semantics.

    Indices rather than a boolean mask so a resample can repeat a field of
    view: `np.isin` collapses duplicates, which would silently turn a
    with-replacement bootstrap into subsampling without replacement.
    """
    cache = ctx['cache']
    L_sub = cache['L'][sub_idx]
    jumps = np.asarray(cache['jumps_per_track'])[sub_idx].astype(np.float32)
    prior = np.full(L_sub.shape[1:], CONC_PARAM, dtype=np.float32)
    w = np.ones_like(jumps) if EQUAL_WEIGHT else None
    n_raw, _, tau = _vbem_infer(L_sub, jumps, prior, ctx['likelihood'],
                                _VBEM_ITERS, tol=_VBEM_TOL, weights=w,
                                apply_defoc=False)
    return L_sub, np.asarray(n_raw, np.float64), tau


def _fit(ctx, fovs):
    """One VBEM fit over a group's sub-trajectories. Returns
    (L_sub, fov_of_subtrack, n_raw, exp_log_tau)."""
    cache = ctx['cache']
    sub_idx = np.flatnonzero(np.isin(cache['fov_of_subtrack'], fovs))
    L_sub, n_raw, tau = _fit_idx(ctx, sub_idx)
    return L_sub, cache['fov_of_subtrack'][sub_idx], n_raw, tau


def _marginals(L_sub, tau):
    """Each sub-trajectory's posterior marginal over D, from a converged fit.

    One einsum pass instead of a 36-step loop over the non-contiguous slices
    `L_sub[:, :, j]`, and the normaliser taken as the row sum of the
    unnormalised result rather than by a second full pass over L -- they are
    the same quantity, since sum_d sum_l L[n,d,l] tau[d,l] is exactly that row
    sum. Measured 12x faster (5.74 s -> 0.48 s on 83k sub-trajectories) and
    agrees with the previous implementation to 2.6e-07 relative, i.e. float32
    round-off (diagnostics/bench_marginals.py).
    """
    marginal = np.einsum('ndl,dl->nd', L_sub, tau, optimize=True).astype(np.float64)
    norm = marginal.sum(axis=1)
    norm[norm == 0] = 1.0
    marginal /= norm[:, None]
    return marginal


def occupancy_from_draw(ctx, draw_fovs):
    """Per-FOV occupancy for a possibly-REPEATING list of fields of view.

    One row per draw SLOT, not per distinct FOV, so a bootstrap resample that
    draws the same field of view twice contributes it twice -- which is the
    whole point of resampling with replacement.

    Identical to `per_fov_occupancy` when `draw_fovs` has no repeats, and it
    is the same fit and the same marginals underneath, so a null or interval
    built from this is computed by the same code as the observed statistic
    rather than by a parallel implementation that might drift from it.
    """
    fov_of_subtrack = ctx['cache']['fov_of_subtrack']
    rows_by_fov = {f: np.flatnonzero(fov_of_subtrack == f)
                   for f in np.unique(np.asarray(draw_fovs))}
    blocks = [rows_by_fov[f] for f in draw_fovs]
    sub_idx = np.concatenate(blocks)
    L_sub, _, tau = _fit_idx(ctx, sub_idx)
    marg = _marginals(L_sub, tau)

    occ = np.full((len(draw_fovs), marg.shape[1]), np.nan)
    off = 0
    for i, blk in enumerate(blocks):
        if len(blk):
            occ[i] = marg[off:off + len(blk)].mean(axis=0)
        off += len(blk)
    return occ


def per_fov_occupancy(ctx, group):
    """Per-FOV soft posterior occupancy (n_fov, n_D), rows summing to 1.

    One VBEM fit over the group's sub-trajectories, then the per-FOV mean of
    each sub-trajectory's full posterior over D. This reproduces
    `softpost_region_panels.group_soft_per_fov`, but pins the iteration count
    and tolerance to saSPT's exact semantics rather than relying on the early
    stop not having fired. Defocalization-invariant by construction (see the
    module docstring).
    """
    fovs = _resolve_group(ctx['cache'], group)
    L_sub, fov_local, _, tau = _fit(ctx, fovs)
    marginal = _marginals(L_sub, tau)
    n_D = marginal.shape[1]

    occ = np.full((len(fovs), n_D), np.nan)
    n_per_fov = np.zeros(len(fovs), dtype=int)
    for i, f in enumerate(fovs):
        rows = fov_local == f
        n_per_fov[i] = int(rows.sum())
        if rows.any():
            occ[i] = marginal[rows].mean(axis=0)
    return occ, fovs, n_per_fov


def population_pmf(ctx, group, dz=None):
    """saSPT's own population posterior over D for a group, at an arbitrary
    focal depth.

    Runs the fit once with the correction off to obtain the raw Dirichlet
    statistic, then applies saSPT's correction analytically. This is exactly
    equivalent to refitting at each dz (T2(iii)) and makes a dz sweep free.
    `dz=None` uses the deployed settings value.
    """
    s = ctx['settings']
    dz = float(s['focal_depth']) if dz is None else dz
    dt = float(s['frame_interval'])
    fovs = _resolve_group(ctx['cache'], group)
    L_sub, _, n_raw, _ = _fit(ctx, fovs)
    if not np.isinf(dz):
        fr = np.array([f_remain_rbm(d, 1, dt, dz)[0] for d in ctx['diff_coefs']])
        n_raw = (n_raw.T / fr).T
    p = n_raw.sum(axis=1)
    return p / p.sum(), int(L_sub.shape[0])


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def bin_effect(ctx, occ_a, occ_b, scheme, name_a, name_b, seed=0):
    """Per-bin difference (b - a) in per-FOV soft occupancy, with a FOV
    bootstrap 95% CI, a whole-FOV relabelling null, and BH-corrected
    permutation p-values across the bins of this panel."""
    rows = _subband_summary(ctx['diff_coefs'], occ_a, occ_b, scheme,
                            None, None, seed=seed, n_boot=N_BOOT, n_perm=N_PERM)
    df = pd.DataFrame(rows)
    df['p_perm'] = [perm_p_addone(p, N_PERM) for p in df['p_value_two_sided']]
    df['p_bh'] = bh(df['p_perm'].to_numpy())
    df['p_holm'] = holm(df['p_perm'].to_numpy())     # conservative alternative
    df['stars_bh'] = [ps.stars(p) for p in df['p_bh']]
    df['group_a'], df['group_b'] = name_a, name_b
    df['n_fov_a'], df['n_fov_b'] = occ_a.shape[0], occ_b.shape[0]
    df['equal_weight'] = EQUAL_WEIGHT
    return df


def plot_w1_vs_null(series, out_pdf, size=None):
    """W1 against its whole permutation null, one row per contrast.

    Shows the null DISTRIBUTION rather than only its 95th percentile: a
    threshold line cannot distinguish an effect sitting just past the shoulder
    from one past every draw, and those mean different things.

    `series` is a list of dicts: {label, w1, boot, null}.
    """
    ps.apply()
    import matplotlib.pyplot as plt
    n = len(series)
    fig, axes = plt.subplots(n, 1, figsize=(3.6, 1.30 * n), sharex=True)
    axes = np.atleast_1d(axes)
    fig.subplots_adjust(left=0.27, right=0.97, bottom=0.17 if n > 1 else 0.24,
                        top=0.80, hspace=0.40)
    for ax, s in zip(axes, series):
        null, boot, w1 = np.asarray(s['null']), np.asarray(s['boot']), s['w1']
        p95 = float(np.percentile(null, 95))
        first = ax is axes[0]
        ax.hist(null, bins=60, color='0.78', edgecolor='none',
                label='FOV-relabelling null' if first else None)
        ax.axvline(p95, color='#CC2222', lw=0.9, ls='--',
                   label='null 95th pct' if first else None)
        y = ax.get_ylim()[1]
        ax.plot([np.percentile(boot, 2.5), np.percentile(boot, 97.5)],
                [y * 0.55] * 2, color=ps.PERTURBED, lw=1.4, solid_capstyle='butt')
        ax.plot([w1], [y * 0.55], 'o', color=ps.PERTURBED, ms=4.5, zorder=5,
                label='observed $W_1$ (95% CI)' if first else None)
        ax.set_ylim(0, y * 1.05)
        ax.set_yticks([])
        ax.set_ylabel(s['label'], rotation=0, ha='right', va='center',
                      fontsize=6.5, labelpad=6)
        note = ('p < 0.001, beyond all null draws' if null.max() < w1
                else 'p < 0.001' if w1 > p95 else 'n.s.')
        ax.text(0.985, 0.86, f"$W_1$={w1:.3f}   {note}", transform=ax.transAxes,
                ha='right', va='top', fontsize=6,
                color=ps.PERTURBED if w1 > p95 else '0.35')
    axes[-1].set_xlabel(r'$W_1$ on the log$_{10}$D axis')
    axes[0].legend(loc='lower center', bbox_to_anchor=(0.42, 1.10), ncol=3,
                   fontsize=5.4, handlelength=1.2, columnspacing=1.0)
    ps.save(fig, tag(out_pdf))
    return tag(out_pdf)


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def density(pmf, diff_coefs):
    """Occupancy mass -> probability density per unit log10 D, so the y-axis is
    a density that integrates to 1 over the plotted log axis and is comparable
    between panels with different grid spacings. See
    PaperFigures/METHODS_posterior_density_note.md."""
    dlog = float(np.diff(np.log10(diff_coefs)).mean())
    return np.asarray(pmf) / dlog


# Density panels are capped at this value rather than auto-scaled. The spike at
# the grid floor (everything slower than the localization-error limit piles onto
# the lowest state) is several-fold taller than any resolved peak, so
# auto-scaling flattens the part of the curve the figure is actually about.
# saSPT's own plot_posterior has a `truncate_y_axis` flag for the same reason.
# Curves exceeding the cap are clipped and `plot_overlay` prints how much.
#
# 1.0 is chosen so it trims only the floor-state pile-up: in the published
# (displacement-weighted) panels every clipped peak sits at D = 0.001, the
# lowest grid state, while the tallest resolved feature -- the free-protein peak
# of Fig 5C at D = 7.2 -- reaches 0.96 and is left intact. The equal-weight
# robustness rebuild does clip that real peak (it reaches 2.04 there); the two
# are kept on the same axis so they stay comparable, and the clip is reported.
YMAX_DENSITY = 1.0


# A safe CSV column key from a legend label: "sgNT (n=83,158)" -> "sgNT".
_LABEL_N = re.compile(r'\s*\(n\s*=\s*[\d,]+\)\s*$')


def _curve_key(label, taken):
    k = re.sub(r'[^0-9A-Za-z]+', '_', _LABEL_N.sub('', str(label))).strip('_')
    k = k or 'curve'
    if k in taken:                      # two curves with the same legend text
        i = 2
        while f'{k}_{i}' in taken:
            i += 1
        k = f'{k}_{i}'
    taken.add(k)
    return k


def write_overlay_curves(curves, diff_coefs, scheme, out_pdf, **spec):
    """Persist an overlay's plotted curves next to its PDF, so the panel can be
    restyled later without refitting.

    Two files, each with one job:

      <stem>_curves.csv    D plus one `pmf_<curve>` column per curve. This is
                           the OLD pipeline's format verbatim
                           (run_equalweight_queue.py:111-112), so the tables are
                           interchangeable with the ones under PaperFigures/ and
                           drop straight into Prism or Illustrator. Only the pmf
                           is stored, never the density -- `density()` derives it
                           from the grid, and storing both invites them to
                           disagree.
      <stem>_plotspec.json Everything else needed to redraw byte-identically:
                           curve order, colour, linestyle, fill state, the region
                           scheme, axis labels and limits. Consumed by
                           replot_overlays.py.

    Why this exists: every OTHER panel type already persisted its plotted
    numbers -- write_stats() for the bin panels, the refit statistics CSVs for
    the W1 panels -- because those numbers are cited in the draft. The overlays cite
    none, so plot_overlay was written as a pure renderer and the sidecar the old
    pipeline had was not carried into the rewrite. The cost of that landed in
    Sep 2026: changing the FILL on Fig 5C needed a 19-minute VBEM rebuild,
    because the curve 5C plots (population_pmf, defocalization-corrected) was
    not written down anywhere. The refit caches under refit_resampling/ look
    like they cover it and do not -- those hold per-FOV MEAN occupancy, a
    different estimand, and only for contrasts that had a refit.
    """
    tagged = tag(out_pdf)
    stem = os.path.splitext(tagged)[0]
    os.makedirs(os.path.dirname(stem) or '.', exist_ok=True)
    taken, cols, meta = set(), {'D': np.asarray(diff_coefs)}, []
    for c in curves:
        k = _curve_key(c['label'], taken)
        cols[f'pmf_{k}'] = np.asarray(c['pmf'])
        meta.append(dict(key=k, label=c['label'], color=c['color'],
                         ls=c.get('ls', '-'), lw=c.get('lw', 1.0),
                         fill=bool(c.get('fill')), alpha=c.get('alpha')))
    pd.DataFrame(cols).to_csv(f'{stem}_curves.csv', index=False)
    with open(f'{stem}_plotspec.json', 'w') as fh:
        json.dump(dict(
            curves_csv=os.path.basename(f'{stem}_curves.csv'),
            # No output path is stored, deliberately. The spec sits beside its
            # own PDF, so a reader derives the path from where the spec IS --
            # which keeps the sidecar valid if the tree moves, and keeps
            # absolute lab paths out of the tlsd-spt-analysis export (whose own
            # check only greps src/*.py, so a path in here would slip through).
            # `suffix` is what lets a reader strip the variant tag back off to
            # recover the untagged path plot_overlay expects.
            suffix=SUFFIX,
            scheme={k: list(v) for k, v in scheme.items()},
            curves=meta, **spec), fh, indent=2, default=float)
    return f'{stem}_curves.csv'


def plot_overlay(curves, diff_coefs, scheme, out_pdf, size='overlay',
                 xlabel=None, ylabel=None, xlim=(1e-3, 3e1), legend_loc='upper left',
                 ymax=YMAX_DENSITY, fill_alpha=None):
    """Density-vs-D overlay with dotted bin boundaries and roman-numeral bin
    headers (Fig 5C, 5F-left, S7B, S7E, S7J).

    `curves` is a list of dicts: {pmf, label, color, ls, fill, alpha}.
    `ymax` caps the y-axis (None to auto-scale); clipping is reported on stdout
    so a truncated peak is never silent.

    Filling is opt-in per curve (`fill=True`) so a panel that wants stroke-only
    curves -- 5F-left, where zone is encoded by dash pattern -- gets exactly
    that. Where fills ARE
    asked for, every filled curve is tinted in its OWN colour, because the y-axis
    is a density meant to be read by area (see paperfig_style.LEGEND_DENSITY_NOTE)
    and filling only one of two curves invites reading by peak height instead.
    This matches the older PaperFigures overlays (occupancy_compare.py:1498-1500).

    Alpha is chosen from how many curves are filled: two translucent fills at
    0.18 stay legible where they overlap, but three at 0.18 stack into a tint
    darker than any single curve's, so 3+ drops to 0.12. `fill_alpha` overrides
    the panel default and a per-curve `alpha` key overrides both.
    """
    fig, ax = ps.figure(size)
    n_fill = sum(1 for c in curves if c.get('fill'))
    a_default = fill_alpha if fill_alpha is not None else (0.18 if n_fill <= 2
                                                          else 0.12)
    clipped = []
    for c in curves:
        y = density(c['pmf'], diff_coefs)
        if ymax is not None and y.max() > ymax:
            d_at = float(diff_coefs[int(np.argmax(y))])
            clipped.append((c['label'], float(y.max()), d_at))
        ax.plot(diff_coefs, y, color=c['color'], ls=c.get('ls', '-'),
                lw=c.get('lw', 1.0), label=c['label'], zorder=3)
        if c.get('fill'):
            ax.fill_between(diff_coefs, 0, y, color=c['color'],
                            alpha=c.get('alpha', a_default), lw=0, zorder=2)
    ax.set_xscale('log')
    ax.set_xlim(*xlim)
    ax.set_ylim(0, ymax) if ymax is not None else ax.set_ylim(bottom=0)
    if clipped:
        # tag(), not out_pdf: under PANEL_EQUAL_WEIGHT the file written is the
        # _equalweight variant, and naming the untagged file here made a clip in
        # the robustness rebuild look like a clip in the published panel.
        print(f'  [{os.path.basename(tag(out_pdf))}] y-axis capped at {ymax}; '
              f'clipped: ' + '; '.join(
                  f'{lab} peaks {pk:.2f} at D={d:.3g}' for lab, pk, d in clipped))
    ax.set_xlabel(xlabel or ps.XLABEL_D)
    ax.set_ylabel(ylabel or ps.YLABEL_DENSITY)
    ax.legend(loc=legend_loc)
    ps.annotate_bins_log(ax, scheme)
    ps.save(fig, tag(out_pdf))
    # Unconditional, not opt-in: a sidecar a call site can forget to request is
    # exactly how the last one went missing.
    write_overlay_curves(curves, diff_coefs, scheme, out_pdf, size=size,
                         xlabel=xlabel, ylabel=ylabel, xlim=list(xlim),
                         legend_loc=legend_loc, ymax=ymax,
                         fill_alpha=fill_alpha)
    return tag(out_pdf)


def plot_bin_effect(series, out_pdf, size='effect', ylabel=None, brackets=None,
                    title=None, zero_line=True, legend_ncol=None):
    """Per-bin effect panel: point estimate + 95% bootstrap CI, one or two
    series per bin (Fig 5F-right, S7C, S7D).

    `series` is a list of dicts: {df, label, color, marker}. Each `df` is a
    `bin_effect` frame; all must share bin order.
    `brackets` optionally maps bin index -> significance string, drawn as a
    bracket spanning the two series of that bin (Fig 5F-right).
    """
    fig, ax = ps.figure(size)
    n_bins = len(series[0]['df'])
    k = len(series)
    offs = np.linspace(-0.16, 0.16, k) if k > 1 else [0.0]

    for si, s in enumerate(series):
        df = s['df']
        x = np.arange(n_bins) + offs[si]
        y = df['observed_diff'].to_numpy()
        lo = y - df['ci_lo'].to_numpy()
        hi = df['ci_hi'].to_numpy() - y
        ax.errorbar(x, y, yerr=[lo, hi], fmt=s.get('marker', 'o'),
                    ms=s.get('ms', 3.2), color=s['color'], ecolor=s['color'],
                    elinewidth=0.8, capsize=0, lw=0, label=s['label'], zorder=3)
        if k == 1:
            for xi, yi, st in zip(x, df['ci_hi'].to_numpy(), df['stars_bh']):
                if st:
                    ax.text(xi, yi, st, ha='center', va='bottom', fontsize=7)

    if zero_line:
        ax.axhline(0, color='0.35', lw=0.6, zorder=1)
    # headroom so stars/brackets drawn at ci_hi are never clipped
    ax.margins(y=0.18)
    ax.set_xlim(-0.5, n_bins - 0.5)
    ax.set_xticks([])
    ax.set_ylabel(ylabel or ps.YLABEL_DELTA)
    ps.annotate_bins_categorical(ax, n_bins)
    if brackets:
        span = abs(offs[-1] - offs[0]) / 2 if k > 1 else 0.2
        top = max(s['df']['ci_hi'].max() for s in series)
        rng_y = ax.get_ylim()[1] - ax.get_ylim()[0]
        for b, txt in brackets.items():
            if txt:
                ps.sig_bracket(ax, b - span, b + span,
                               top + 0.04 * rng_y, txt, dy=0.015 * rng_y)
    if len(series) > 1:
        # below the axes, never inside: these panels are dense enough that any
        # in-axes legend lands on a data point or gets clipped
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.02),
                  ncol=legend_ncol or min(len(series), 2), handletextpad=0.4,
                  columnspacing=1.0, borderpad=0.2)
    if title:
        ax.set_title(title, pad=14)
    ps.save(fig, tag(out_pdf))
    return tag(out_pdf)


def write_stats(df, out_csv, **extra):
    for k, v in extra.items():
        df[k] = v
    df.to_csv(tag(out_csv), index=False)
    return tag(out_csv)


def tag(path):
    """Insert the weighting suffix before the extension, so an equal-weight
    rebuild sits beside the saSPT-weighted panel rather than replacing it."""
    if not SUFFIX:
        return path
    root, ext = os.path.splitext(path)
    return f'{root}{SUFFIX}{ext}'
