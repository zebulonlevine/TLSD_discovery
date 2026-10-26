"""
occupancy_compare.py
=====================

Reusable two-group occupancy comparator off one cached saSPT likelihood
(saspt_likelihood_cache.py). See docs/differentiability_handoff.md for the
full design rationale; this module implements sections 4 (core: naive,
label-blind, no-refit) and 5-B (optional escalation: conditionwise VBEM-
refit permutation).

Fixed framing (decided in the handoff, do not relitigate here):
  * FOV is the exchangeable permutation unit. Never permute subtracks.
  * The core comparator is defended on naive (prior-free) occupancy, not
    on joint-fit or per-condition posterior slices.
  * No translation/delta statistic anywhere.
  * Regional (sub-band) tests are opt-in per comparison (`subbands=`).

Explicit, flagged decisions (handoff section 7 -- "flag, don't silently
choose"):
  1. Defocalization correction is applied to `naive_occs`-style POPULATION
     aggregates in the deployed pipeline, but is NOT applied here to the
     per-subtrack / per-FOV naive occupancy vectors that `compare_groups`
     is built on. It's a fixed multiplicative function of D alone (same
     for both groups), so it cancels in a *signed difference* comparison
     up to a renormalization constant -- omitting it doesn't change the
     conclusion (parsimony), but flagging it since it does change the
     absolute shape of plot 1's curves relative to `naive_occs`. Pass
     `correct_defocalization=True` to apply it if an apples-to-apples
     comparison against `naive_occs` plots is ever needed.
  2. Immobile-spike handling: full grid range by default (`d_range=None`),
     not mobile-range-only -- matches `naive_occs`' own full-grid
     normalization. Restrict with `d_range=(d_lo, d_hi)` if needed.
  3. Band level: 95% max-T by default (`band_level=0.95`), exposed as an
     arg.
  4. Ribbon: SEM by default (`ribbon='sem'`), matches the pipeline's
     existing per-FOV-SEM convention (spt_compare.bar_graph); `'iqr'` or
     `'both'` also supported.
"""

import json
import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from saspt_likelihood_cache import load_likelihood_cache

_EPS = 1e-20


# ---------------------------------------------------------------------------
# group resolution / validation
# ---------------------------------------------------------------------------

def _load_cache_if_path(cache):
    if isinstance(cache, (str, os.PathLike)):
        return load_likelihood_cache(cache)
    return cache


def _resolve_group(cache, spec):
    """spec: dict of {metadata_column: value_or_list} filter, or an
    explicit list/array of fov_basenames."""
    meta = cache['fov_metadata']
    if isinstance(spec, dict):
        mask = np.ones(len(meta), dtype=bool)
        for col, val in spec.items():
            vals = val if isinstance(val, (list, tuple, set, np.ndarray)) else [val]
            mask &= meta[col].isin(vals).to_numpy()
        fovs = meta.loc[mask, 'fov_basename'].tolist()
    else:
        fovs = list(spec)
    return fovs


def _validate_groups(cache, fovs_a, fovs_b, negative_control):
    set_a, set_b = set(fovs_a), set(fovs_b)
    if not negative_control:
        overlap = set_a & set_b
        if overlap:
            raise ValueError(
                f"group_a and group_b overlap in {len(overlap)} FOV(s): "
                f"{sorted(overlap)[:5]}{'...' if len(overlap) > 5 else ''} "
                f"-- pass negative_control=True if this is intentional")
    for name, fovs in (('group_a', fovs_a), ('group_b', fovs_b)):
        if len(fovs) < 4:
            raise ValueError(
                f"{name} has only {len(fovs)} FOV(s); need >= 4 for a "
                f"meaningful FOV-level permutation floor")
    if len(fovs_a) > 0 and len(fovs_b) > 0:
        ratio = max(len(fovs_a), len(fovs_b)) / min(len(fovs_a), len(fovs_b))
        if ratio >= 2.0:
            print(f"WARN: large FOV-count imbalance (n_a={len(fovs_a)}, "
                  f"n_b={len(fovs_b)}, ratio={ratio:.1f}x) -- widens the "
                  f"permutation band")


def _subtrack_counts(cache, fovs):
    mask = np.isin(cache['fov_of_subtrack'], fovs)
    return int(mask.sum())


# ---------------------------------------------------------------------------
# core naive-occupancy machinery
# ---------------------------------------------------------------------------

def _marginal_pmf_per_subtrack(cache, representation='naive'):
    """Sum the loc_error axis; already per-subtrack normalized (L sums to
    1 over the full (D, loc_error) grid per subtrack).

    representation='naive' (default): the soft per-subtrack probability
    curve itself -- what compare_groups has always used.

    representation='mle_d': collapse each subtrack to a hard point
    estimate (the D bin with the highest naive likelihood, matching how
    the rest of this pipeline defines MLE_D -- trackingCodeFunctions.
    executeSASPT's `diff_coefs[marginal_D.argmax(axis=0)]` -- except here
    from the naive, no-refit likelihood rather than the VBEM posterior,
    since this is the label-blind core comparator). One-hot at the argmax
    bin, so a per-FOV MEAN of these vectors (the same downstream machinery
    naive occupancy uses) is exactly that FOV's MLE_D histogram --
    directly comparable to this pipeline's existing MLE_D histogram plots
    (e.g. histogram_mle_*.pdf), not just an abstract soft-probability
    curve.
    """
    pmf = cache['L'].sum(axis=2).astype(np.float64)  # (n_subtracks, n_D)
    if representation == 'naive':
        return pmf
    if representation == 'mle_d':
        idx = np.argmax(pmf, axis=1)
        one_hot = np.zeros_like(pmf)
        one_hot[np.arange(pmf.shape[0]), idx] = 1.0
        return one_hot
    raise ValueError(f"unknown representation: {representation!r} "
                      f"(expected 'naive' or 'mle_d')")


def _per_fov_occupancy(pmf_D, fov_of_subtrack, fovs):
    """Returns (occ, n_subtracks_per_fov): occ has shape (len(fovs), n_D),
    row i = mean over FOV i's subtracks of the per-subtrack P(D)."""
    n_D = pmf_D.shape[1]
    occ = np.full((len(fovs), n_D), np.nan)
    n_sub = np.zeros(len(fovs), dtype=int)
    for i, fov in enumerate(fovs):
        rows = fov_of_subtrack == fov
        n = int(rows.sum())
        n_sub[i] = n
        if n > 0:
            occ[i] = pmf_D[rows].mean(axis=0)
    return occ, n_sub


def _studentize(occ_a, occ_b):
    """Delta(d) = mean_B - mean_A, studentized by SE of the difference
    (pooled across-FOV variance in each group, standard two-sample SE)."""
    na, nb = occ_a.shape[0], occ_b.shape[0]
    mean_a, mean_b = np.nanmean(occ_a, axis=0), np.nanmean(occ_b, axis=0)
    var_a = np.nanvar(occ_a, axis=0, ddof=1) if na > 1 else np.zeros_like(mean_a)
    var_b = np.nanvar(occ_b, axis=0, ddof=1) if nb > 1 else np.zeros_like(mean_b)
    se = np.sqrt(var_a / na + var_b / nb)
    delta = mean_b - mean_a
    delta_star = delta / (se + _EPS)
    return delta, se, delta_star


def _fov_permutation_null(occ_all, n_a, n_perm, seed):
    """FOV-block permutation: relabel the pooled FOV set into pseudo-A/B of
    the real sizes; whole FOVs move, subtracks never do. `occ_all` is the
    ALREADY-COMPUTED per-FOV occupancy matrix (n_total_fov, n_D) for the
    concatenation fovs_a+fovs_b -- a FOV's own occupancy vector is fixed
    regardless of which pseudo-group it's relabeled into, so permutation
    only needs to reshuffle rows of this small matrix, never touch
    subtrack-level data again. Returns the null array of
    M_b = max_d |Delta*_b(d)| (each Delta*_b studentized by ITS OWN
    permutation's SE, i.e. a proper max-T self-normalizing statistic)."""
    rng = np.random.default_rng(seed)
    n_total = occ_all.shape[0]
    M = np.empty(n_perm)
    for b in range(n_perm):
        perm = rng.permutation(n_total)
        occ_a = occ_all[perm[:n_a]]
        occ_b = occ_all[perm[n_a:]]
        _, _, delta_star_b = _studentize(occ_a, occ_b)
        M[b] = np.nanmax(np.abs(delta_star_b))
    return M


def _apply_d_range(diff_coefs, d_range):
    """Inclusive [d_lo, d_hi] analysis window. Deliberately CLOSED, unlike the
    region-band masks (which are half-open [lo, hi) because they partition the
    axis and must not double-count an on-grid boundary point). This is a single
    window, not a partition, so including both endpoints is the intent."""
    if d_range is None:
        return np.ones(len(diff_coefs), dtype=bool)
    d_lo, d_hi = d_range
    return (diff_coefs >= d_lo) & (diff_coefs <= d_hi)


# ---------------------------------------------------------------------------
# subbands
# ---------------------------------------------------------------------------

def _subband_summary(diff_coefs, occ_a, occ_b, subbands, n_perm_null_occ_a,
                      n_perm_null_occ_b, seed, n_boot=2000, n_perm=2000):
    """Per-subband summed-occupancy contrast: a FOV-block relabeling null
    (the "reshuffle" -- is the observed group difference bigger than whole-
    FOV label noise?) reported alongside a FOV bootstrap 95% CI on the
    observed difference itself (the "resample" -- how uncertain is the
    observed effect, holding group labels fixed?). A region reads as real
    when its bootstrap CI clears the null's 95% band; see
    plot_region_effect_vs_null.

    n_perm was previously hard-coded at 2000 and not recorded on the output
    rows, so downstream code could not tell what resolution floor the
    permutation p carried. It is now a parameter (DEFAULT UNCHANGED at 2000 --
    do not alter it casually, the RNG stream and therefore every committed
    number depends on the draw count) and is emitted on each row."""
    out = []
    rng = np.random.default_rng(seed)
    for label, (d_lo, d_hi) in subbands.items():
        band_mask = (diff_coefs >= d_lo) & (diff_coefs < d_hi)
        sum_a = occ_a[:, band_mask].sum(axis=1)
        sum_b = occ_b[:, band_mask].sum(axis=1)
        obs = float(np.nanmean(sum_b) - np.nanmean(sum_a))

        # FOV-block relabeling null (reshuffle): pool both groups' per-FOV
        # sums and re-split into pseudo-groups of the real sizes.
        pooled = np.concatenate([sum_a, sum_b])
        na = len(sum_a)
        null = np.empty(n_perm)
        for i in range(n_perm):
            perm = rng.permutation(len(pooled))
            null[i] = pooled[perm[na:]].mean() - pooled[perm[:na]].mean()
        p_one_sided = float((np.abs(null) >= abs(obs)).mean())

        # FOV bootstrap (resample): resample FOVs with replacement WITHIN
        # each group (group labels held fixed, unlike the null above).
        boot = np.empty(n_boot)
        na_, nb_ = len(sum_a), len(sum_b)
        for i in range(n_boot):
            ra = rng.integers(0, na_, na_)
            rb = rng.integers(0, nb_, nb_)
            boot[i] = sum_b[rb].mean() - sum_a[ra].mean()
        ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])

        out.append({
            'label': label, 'd_lo': d_lo, 'd_hi': d_hi,
            'observed_diff': obs, 'p_value_two_sided': p_one_sided,
            'null_mean': float(null.mean()), 'null_sd': float(null.std()),
            'null_p2p5': float(np.percentile(null, 2.5)),
            'null_p97p5': float(np.percentile(null, 97.5)),
            'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi),
            'n_boot': n_boot, 'n_perm': n_perm,
        })
    return out


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def _ribbon(occ, ribbon):
    n = occ.shape[0]
    mean = np.nanmean(occ, axis=0)
    lo_sem = mean - np.nanstd(occ, axis=0, ddof=1) / np.sqrt(max(n, 1))
    hi_sem = mean + np.nanstd(occ, axis=0, ddof=1) / np.sqrt(max(n, 1))
    lo_iqr = np.nanpercentile(occ, 25, axis=0)
    hi_iqr = np.nanpercentile(occ, 75, axis=0)
    if ribbon == 'sem':
        return [('SEM', lo_sem, hi_sem, 0.3)]
    if ribbon == 'iqr':
        return [('IQR', lo_iqr, hi_iqr, 0.15)]
    return [('SEM', lo_sem, hi_sem, 0.35), ('IQR', lo_iqr, hi_iqr, 0.12)]


def _style_logD_axis(ax):
    ax.set_xscale('log')
    ax.set_xlabel(r'Diffusion coefficient D ($\mu m^2$/s)')


def plot_naive_occupancy(diff_coefs, occ_a, occ_b, n_sub_a, n_sub_b,
                          name_a, name_b, ribbon='sem', out_pdf=None,
                          color_a='#3366cc', color_b='#cc3333',
                          y_label='mean naive occupation (per FOV)'):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for occ, n_sub, name, color in (
            (occ_a, n_sub_a, name_a, color_a), (occ_b, n_sub_b, name_b, color_b)):
        mean = np.nanmean(occ, axis=0)
        label = f"{name} (n_FOV={occ.shape[0]}, n_subtracks={int(n_sub.sum()):,})"
        ax.plot(diff_coefs, mean, color=color, lw=1.8, label=label)
        for rname, lo, hi, alpha in _ribbon(occ, ribbon):
            ax.fill_between(diff_coefs, lo, hi, color=color, alpha=alpha, lw=0)
    _style_logD_axis(ax)
    ax.set_ylabel(y_label)
    ax.legend(loc='upper right', fontsize=8)
    ax.set_title(f'naive occupancy: {name_a} vs {name_b}')
    fig.tight_layout()
    if out_pdf:
        fig.savefig(out_pdf)
    return fig


def plot_occupancy_difference(diff_coefs, delta, name_a, name_b,
                               diff_sign='b_minus_a', out_pdf=None):
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.axhline(0, color='k', lw=0.8)
    ax.plot(diff_coefs, delta, color='#333333', lw=1.5)
    _style_logD_axis(ax)
    ax.set_ylabel(r'$\Delta$(occupation)')
    ax.set_title(f'occupancy difference ({diff_sign}: {name_b} - {name_a})')
    fig.tight_layout()
    if out_pdf:
        fig.savefig(out_pdf)
    return fig


def plot_permutation_threshold(diff_coefs, delta, se, c95, name_a, name_b,
                                n_perm, seed, out_pdf=None):
    band = c95 * se
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axhline(0, color='k', lw=0.8)
    ax.fill_between(diff_coefs, -band, band, color='gray', alpha=0.3,
                     label=f'permutation noise floor ({int(c95*100)/100:.2f}$\\sigma$, '
                           f'n_perm={n_perm}, seed={seed})')
    exceeds = np.abs(delta) > band
    ax.plot(diff_coefs, delta, color='#333333', lw=1.3, label=r'$\Delta$(D)')
    if exceeds.any():
        ax.plot(diff_coefs[exceeds], delta[exceeds], 'o', color='#cc3333',
                ms=3, label='exceeds floor', zorder=5)
    _style_logD_axis(ax)
    ax.set_ylabel(r'$\Delta$(occupation)')
    ax.set_title(f'permutation threshold: {name_a} vs {name_b}')
    ax.legend(loc='upper right', fontsize=7)
    fig.tight_layout()
    if out_pdf:
        fig.savefig(out_pdf)
    return fig


def plot_combined(diff_coefs, occ_a, occ_b, n_sub_a, n_sub_b, delta, se, c95,
                   name_a, name_b, n_perm, seed, ribbon='sem',
                   diff_sign='b_minus_a', out_pdf=None,
                   y_label='mean naive occupation'):
    fig, axes = plt.subplots(3, 1, figsize=(7, 11), sharex=True)
    for occ, n_sub, name, color in (
            (occ_a, n_sub_a, name_a, '#3366cc'), (occ_b, n_sub_b, name_b, '#cc3333')):
        mean = np.nanmean(occ, axis=0)
        label = f"{name} (n_FOV={occ.shape[0]}, n_subtracks={int(n_sub.sum()):,})"
        axes[0].plot(diff_coefs, mean, color=color, lw=1.8, label=label)
        for rname, lo, hi, alpha in _ribbon(occ, ribbon):
            axes[0].fill_between(diff_coefs, lo, hi, color=color, alpha=alpha, lw=0)
    axes[0].set_ylabel(y_label)
    axes[0].legend(loc='upper right', fontsize=7)
    axes[0].set_title(f'{name_a} vs {name_b}')

    axes[1].axhline(0, color='k', lw=0.8)
    axes[1].plot(diff_coefs, delta, color='#333333', lw=1.5)
    axes[1].set_ylabel(r'$\Delta$(occupation)')
    axes[1].set_title(f'difference ({diff_sign})')

    band = c95 * se
    axes[2].axhline(0, color='k', lw=0.8)
    axes[2].fill_between(diff_coefs, -band, band, color='gray', alpha=0.3,
                          label=f'floor (n_perm={n_perm}, seed={seed})')
    exceeds = np.abs(delta) > band
    axes[2].plot(diff_coefs, delta, color='#333333', lw=1.3)
    if exceeds.any():
        axes[2].plot(diff_coefs[exceeds], delta[exceeds], 'o', color='#cc3333', ms=3)
    axes[2].set_ylabel(r'$\Delta$(occupation)')
    axes[2].legend(loc='upper right', fontsize=7)
    _style_logD_axis(axes[2])
    fig.tight_layout()
    if out_pdf:
        fig.savefig(out_pdf)
    return fig


def plot_region_effect_vs_null(subband_rows, name_a, name_b,
                                diff_sign='b_minus_a', out_pdf=None,
                                positive_label=None, negative_label=None,
                                style='estimation'):
    """Per-named-region observed occupancy-mass difference (point + FOV-
    bootstrap 95% CI, the "resample") plotted against the FOV-reshuffle
    permutation null's 95% band (the "reshuffle"). Replaces the per-FOV
    strip + Welch-t-test panel: a region reads as real when its CI clears
    the null band, not when a per-region t-test crosses p<0.05.

    `name_a`/`name_b` MUST correspond to whatever produced `observed_diff`
    as mean(group/effect B) - mean(group/effect A) in `subband_rows` --
    get this backwards and every sign in the plot reads opposite to the
    legend. `positive_label`/`negative_label` state in plain language what
    each direction means (default: "more in {name_b}" / "more in
    {name_a}"); ALWAYS pass these explicitly for a difference-in-
    differences panel, where "positive" means an interaction/relative
    comparison, not a simple "more in group B"."""
    import matplotlib.patches as mpatches

    if positive_label is None:
        positive_label = f'more in {name_b}'
    if negative_label is None:
        negative_label = f'more in {name_a}'

    labels = [r['label'] for r in subband_rows]
    x = np.arange(len(labels))
    obs = np.array([r['observed_diff'] for r in subband_rows])
    ci_lo = np.array([r['ci_lo'] for r in subband_rows])
    ci_hi = np.array([r['ci_hi'] for r in subband_rows])
    null_lo = np.array([r['null_p2p5'] for r in subband_rows])
    null_hi = np.array([r['null_p97p5'] for r in subband_rows])
    clears = (ci_lo > null_hi) | (ci_hi < null_lo)
    excl_zero = (ci_lo > 0) | (ci_hi < 0)
    estimation = (style == 'estimation')

    fig, ax = plt.subplots(figsize=(1.1 * len(labels) + 2, 4.5))
    width = 0.6
    if not estimation:
        # legacy view: gray FOV-reshuffle null band behind each region
        for xi, lo, hi in zip(x, null_lo, null_hi):
            ax.add_patch(mpatches.Rectangle((xi - width / 2, lo), width, hi - lo,
                                             facecolor='gray', alpha=0.35, lw=0, zorder=1))
    ax.axhline(0, color='k', lw=0.8, zorder=2)

    # estimation view (paper default): colour by whether the 95% CI excludes 0;
    # legacy view: colour by whether the CI clears the null band.
    flag = excl_zero if estimation else clears
    colors = ['#cc3333' if c else '#3366cc' for c in flag]
    yerr = np.vstack([obs - ci_lo, ci_hi - obs])
    ax.errorbar(x, obs, yerr=yerr, fmt='none', ecolor='k', elinewidth=1.2,
                capsize=3, zorder=3)
    ax.scatter(x, obs, c=colors, s=55, zorder=4, edgecolor='k', linewidth=0.6)

    # Holm-adjusted significance markers, only if a caller (annotate_multiplicity)
    # has attached p_holm. Absent the key this is a no-op, so every existing
    # caller is unaffected. Star sits outside the whisker, on the side the
    # effect points.
    if subband_rows and 'p_holm' in subband_rows[0]:
        from multiplicity import stars as _stars
        span = float(np.nanmax(ci_hi) - np.nanmin(ci_lo)) or 1.0
        pad = 0.04 * span
        n_star = 0
        for xi, o, lo, hi, r in zip(x, obs, ci_lo, ci_hi, subband_rows):
            mark = _stars(r.get('p_holm'))
            if not mark:
                continue
            n_star += 1
            if o >= 0:
                ax.text(xi, hi + pad, mark, ha='center', va='bottom',
                        fontsize=11, fontweight='bold', zorder=5)
            else:
                ax.text(xi, lo - pad, mark, ha='center', va='top',
                        fontsize=11, fontweight='bold', zorder=5)
        fam = subband_rows[0].get('family', '')
        m = subband_rows[0].get('family_size', len(subband_rows))
        ax.set_title(ax.get_title() + (f'\n* Holm-adjusted p<0.05 within family '
                                       f'{fam} (m={m}); {n_star}/{len(subband_rows)} regions'),
                     fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right')
    if estimation:
        ax.set_ylabel('Δ occupancy (fraction of population)\n'
                      '(point = observed; whisker = FOV-bootstrap 95% CI)')
        ax.set_title(f'{name_a} vs {name_b} — region effect (95% CI)')
    else:
        ax.set_ylabel('Δ(region occupancy mass)\n'
                      '(point = observed; whisker = FOV-bootstrap 95% CI)')
        ax.set_title(f'{name_a} vs {name_b} — region effect vs FOV-reshuffle null')

    # Explicit, in-figure statement of what each sign means -- this is the
    # part that was ambiguous before, especially for DiD panels.
    ax.annotate(f'▲ positive = {positive_label}', xy=(0.02, 0.98),
                xycoords='axes fraction', ha='left', va='top', fontsize=8,
                style='italic', color='#333333')
    ax.annotate(f'▼ negative = {negative_label}', xy=(0.02, 0.02),
                xycoords='axes fraction', ha='left', va='bottom', fontsize=8,
                style='italic', color='#333333')

    if estimation:
        handles = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#cc3333',
                       markeredgecolor='k', label='95% CI excludes 0'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#3366cc',
                       markeredgecolor='k', label='95% CI includes 0'),
        ]
    else:
        handles = [
            mpatches.Patch(facecolor='gray', alpha=0.35, label='FOV-reshuffle null (95%)'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#cc3333',
                       markeredgecolor='k', label='clears null'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#3366cc',
                       markeredgecolor='k', label='within null'),
        ]
    ax.legend(handles=handles, fontsize=7, loc='center right')
    fig.tight_layout()
    if out_pdf:
        fig.savefig(out_pdf)
        # Always persist the underlying stats table next to the panel so the
        # per-region observed diff / bootstrap CI / null band / permutation
        # p-value are recorded on disk, not only in the figure/stdout.
        # NOTE: the two effect_* columns are DiD-only and were previously
        # dropped here, so the difference-in-differences sidecars silently lost
        # the two per-group effects the DiD is built from. The multiplicity
        # columns are added by annotate_multiplicity.py and are absent on a
        # fresh write; the comprehension below filters by presence, so listing
        # them here is safe for every caller.
        cols = ['label', 'd_lo', 'd_hi', 'observed_diff',
                'effect_a_drug_minus_none', 'effect_b_drug_minus_none',
                'ci_lo', 'ci_hi', 'null_p2p5', 'null_p97p5', 'null_mean',
                'null_sd', 'p_value_two_sided', 'p_perm_addone', 'p_holm',
                'p_bh', 'holm_reject', 'family', 'family_size',
                'n_perm', 'n_boot']
        df = pd.DataFrame(subband_rows)
        df['clears_null'] = (df['ci_lo'] > df['null_p97p5']) | (df['ci_hi'] < df['null_p2p5'])
        df['ci_excludes_zero'] = (df['ci_lo'] > 0) | (df['ci_hi'] < 0)
        df['name_a'] = name_a
        df['name_b'] = name_b
        df[[c for c in cols if c in df.columns] + ['clears_null', 'ci_excludes_zero', 'name_a', 'name_b']].to_csv(
            out_pdf[:-4] + '_stats.csv' if out_pdf.endswith('.pdf') else out_pdf + '_stats.csv',
            index=False)
    return fig


# ---------------------------------------------------------------------------
# difference-in-differences (drug x genotype interaction)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# core deliverable
# ---------------------------------------------------------------------------

def compare_groups(cache, group_a, group_b, name_a="A", name_b="B",
                    out_dir=".", n_perm=2000, seed=0, diff_sign="b_minus_a",
                    subbands=None, negative_control=False, band_level=0.95,
                    ribbon='sem', d_range=None, correct_defocalization=False,
                    representation='naive', verbose=True):
    """
    Compare two pre-specified groups of FOVs on prior-free (naive)
    occupancy, with an FOV-block permutation noise threshold. See module
    docstring / docs/differentiability_handoff.md section 4.

    Parameters
    ----------
    cache : dict or str
        A loaded cache (saspt_likelihood_cache.load_likelihood_cache) or a
        path to one.
    group_a, group_b : dict or list
        Each a metadata filter dict (matched against cache['fov_metadata']
        columns) OR an explicit list of fov_basenames.
    negative_control : bool
        If True, group_b is ignored; group_a's FOVs are split into two
        random (seeded) halves and compared against each other. Run this
        before trusting any real result on a new cache.
    representation : 'naive' or 'mle_d'
        'naive' (default): per-subtrack soft probability curve, as
        always. 'mle_d': collapse each subtrack to its single most-likely
        D bin first (matching this pipeline's own MLE_D convention, just
        from the naive/no-refit likelihood -- see
        `_marginal_pmf_per_subtrack`), so the comparison is against a
        per-FOV MLE_D histogram instead of a soft occupancy curve. Same
        FOV-permutation machinery either way; output files get a
        `_mle_d` tag so they don't collide with the naive run.

    Returns
    -------
    dict
        results_dict; also written to `{out_dir}/results_{name_a}_vs_
        {name_b}[_negctrl].json` and `.csv` (subband rows only, if any).
    """
    if correct_defocalization:
        raise NotImplementedError(
            "correct_defocalization=True not implemented -- flagged as an "
            "explicit open decision in the module docstring; the core "
            "comparator intentionally omits it (see decision 1)")

    cache = _load_cache_if_path(cache)
    os.makedirs(out_dir, exist_ok=True)
    diff_coefs = cache['diff_coefs']
    fov_of_subtrack = cache['fov_of_subtrack']

    fovs_a = _resolve_group(cache, group_a)
    suffix = ''
    if negative_control:
        rng = np.random.default_rng(seed)
        perm = rng.permutation(fovs_a)
        half = len(perm) // 2
        fovs_b = list(perm[half:])
        fovs_a = list(perm[:half])
        name_b = f"{name_a}_half2"
        name_a = f"{name_a}_half1"
        suffix = '_negctrl'
    else:
        fovs_b = _resolve_group(cache, group_b)

    _validate_groups(cache, fovs_a, fovs_b, negative_control)

    n_sub_a_total = _subtrack_counts(cache, fovs_a)
    n_sub_b_total = _subtrack_counts(cache, fovs_b)
    if verbose:
        print(f"[{name_a} vs {name_b}] n_FOV_a={len(fovs_a)}, "
              f"n_FOV_b={len(fovs_b)}, n_subtracks_a={n_sub_a_total:,}, "
              f"n_subtracks_b={n_sub_b_total:,}, n_perm={n_perm}, seed={seed}")

    pmf_D = _marginal_pmf_per_subtrack(cache, representation=representation)
    d_mask = _apply_d_range(diff_coefs, d_range)
    y_label = ('mean naive occupation (per FOV)' if representation == 'naive'
               else 'MLE_D fraction (per FOV)')

    # One pass over subtracks for the whole fovs_a+fovs_b union: a FOV's
    # occupancy vector never changes across permutations, so this is the
    # only place subtrack-level data is touched.
    all_fovs = list(fovs_a) + list(fovs_b)
    n_a = len(fovs_a)
    occ_all, n_sub_all = _per_fov_occupancy(pmf_D, fov_of_subtrack, all_fovs)
    occ_a, n_sub_a = occ_all[:n_a], n_sub_all[:n_a]
    occ_b, n_sub_b = occ_all[n_a:], n_sub_all[n_a:]
    delta, se, delta_star = _studentize(occ_a, occ_b)

    M_null = _fov_permutation_null(occ_all, n_a, n_perm=n_perm, seed=seed)
    c95 = float(np.quantile(M_null, band_level))
    band = c95 * se

    exceeds = np.abs(delta) > band
    if d_mask is not None and not d_mask.all():
        exceeds = exceeds & d_mask
    excursions = _excursion_intervals(diff_coefs, exceeds, delta)

    rep_suffix = '' if representation == 'naive' else f'_{representation}'
    tag = f"{name_a}_vs_{name_b}{suffix}{rep_suffix}"
    plot_naive_occupancy(diff_coefs, occ_a, occ_b, n_sub_a, n_sub_b,
                          name_a, name_b, ribbon=ribbon, y_label=y_label,
                          out_pdf=os.path.join(out_dir, f"naive_occupancy_{tag}.pdf"))
    plot_occupancy_difference(diff_coefs, delta, name_a, name_b, diff_sign,
                               out_pdf=os.path.join(out_dir, f"occupancy_difference_{tag}.pdf"))
    plot_permutation_threshold(diff_coefs, delta, se, c95, name_a, name_b,
                                n_perm, seed,
                                out_pdf=os.path.join(out_dir, f"permutation_threshold_{tag}.pdf"))
    plot_combined(diff_coefs, occ_a, occ_b, n_sub_a, n_sub_b, delta, se, c95,
                  name_a, name_b, n_perm, seed, ribbon=ribbon, diff_sign=diff_sign,
                  y_label=y_label,
                  out_pdf=os.path.join(out_dir, f"combined_{tag}.pdf"))
    plt.close('all')

    subband_rows = None
    if subbands:
        subband_rows = _subband_summary(diff_coefs, occ_a, occ_b, subbands,
                                         None, None, seed)
        plot_region_effect_vs_null(
            subband_rows, name_a, name_b, diff_sign=diff_sign,
            out_pdf=os.path.join(out_dir, f"region_effect_vs_null_{tag}.pdf"))
        plt.close('all')

    results = {
        'name_a': name_a, 'name_b': name_b, 'diff_sign': diff_sign,
        'negative_control': negative_control,
        'n_fov_a': len(fovs_a), 'n_fov_b': len(fovs_b),
        'n_subtracks_a': n_sub_a_total, 'n_subtracks_b': n_sub_b_total,
        'n_perm': n_perm, 'seed': seed, 'band_level': band_level, 'c95': c95,
        'excursion_intervals': excursions,
        'fovs_a': fovs_a, 'fovs_b': fovs_b,
        'grid_checksum': cache['grid_checksum'],
        'subbands': subband_rows,
    }
    with open(os.path.join(out_dir, f"results_{tag}.json"), 'w') as f:
        json.dump(results, f, indent=2, default=_json_default)

    per_d = pd.DataFrame({
        'D': diff_coefs, 'mean_occ_a': np.nanmean(occ_a, axis=0),
        'mean_occ_b': np.nanmean(occ_b, axis=0), 'delta': delta, 'se': se,
        'delta_star': delta_star, 'band': band, 'exceeds_band': exceeds,
    })
    per_d.to_csv(os.path.join(out_dir, f"results_{tag}.csv"), index=False)

    if verbose:
        verdict = ("no excursion (Delta stays inside the permutation floor "
                    "everywhere)" if not excursions else
                    f"{len(excursions)} excursion interval(s) exceed the floor")
        print(f"[{name_a} vs {name_b}] c95={c95:.3f}, {verdict}")

    results['per_d'] = per_d
    results['occ_a'] = occ_a
    results['occ_b'] = occ_b
    results['diff_coefs'] = diff_coefs
    return results


def _excursion_intervals(diff_coefs, exceeds, delta):
    intervals = []
    in_run = False
    for i, e in enumerate(exceeds):
        if e and not in_run:
            start_i = i
            in_run = True
        if in_run and (not e or i == len(exceeds) - 1):
            end_i = i if e else i - 1
            sign = 'positive' if np.mean(delta[start_i:end_i + 1]) > 0 else 'negative'
            intervals.append({
                'd_lo': float(diff_coefs[start_i]), 'd_hi': float(diff_coefs[end_i]),
                'sign': sign,
            })
            in_run = False
    return intervals


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(f"not JSON serializable: {type(o)}")


# ---------------------------------------------------------------------------
# Fork B -- conditionwise VBEM-refit permutation (escalation, opt-in)
# ---------------------------------------------------------------------------

def _vbem_infer(L_sub, jumps_sub, prior_dirichlet_param, likelihood, max_iter,
                 tol=1e-5, weights=None, apply_defoc=True):
    """Algebraically-reformulated port of saspt.sa.StateArray.
    infer_posterior's VBEM loop (saspt/sa.py:325-379) -- see Gate 0 notes
    in this module's docstring and saspt_likelihood_cache.py: no supported
    API exists to hand a precomputed likelihood array to a StateArray, so
    we reproduce the inference loop directly against the cached array
    rather than reconstructing detections/TrajectoryGroup objects. This IS
    the escalation Fork-B "Tier 1" path from the handoff (slice L rows,
    run VBEM loop only, no IO / no likelihood recompute).

    saspt's own loop materializes the full per-subtrack responsibility
    matrix R (n_sub, n_D, n_loc) every iteration just to immediately
    contract it back down to a (n_D, n_loc) sufficient statistic. That
    R was originally measured at ~0.7s/iteration for a ~150k-subtrack
    group (dominated by allocating/writing/reducing the full array, not
    FLOPs) -- 150s/fit at max_iter=200, which is why the first working
    version of Fork B failed to complete even a 20-permutation dry run
    inside a 2-hour SLURM job. Substituting the recurrence

        R[s,i,j] = L[s,i,j] * exp_log_tau[i,j] / norm[s]
        norm[s]  = sum_ij L[s,i,j] * exp_log_tau[i,j]
        n[i,j]   = sum_s jumps[s] * R[s,i,j]
                 = exp_log_tau[i,j] * sum_s (jumps[s]/norm[s]) * L[s,i,j]

    replaces that with two (n_sub x n_D*n_loc) matrix-vector products per
    iteration and never materializes R at all -- ~11x faster on a
    correctness-checked subset (identical to the original port to
    float32 precision, ~1e-5 max relative difference), and the ratio
    only improves at the full subtrack counts this is actually run at.
    Initializing `exp_log_tau` to all-ones reproduces the original loop's
    `R = L_sub.copy()` starting point exactly, since L is already
    per-subtrack normalized (norm[s] = 1 on the first iteration).

    Also adds a convergence check on top of saspt's fixed-`max_iter` loop
    (which always runs the full 200 regardless of convergence): stopping
    once the per-state occupation estimate stabilizes (max abs change <
    `tol` between iterations, relative to a floor of 1.0) computes the
    SAME fixed point with fewer wasted iterations. `max_iter` remains a
    hard cap.
    """
    from scipy.special import digamma
    n_sub, n_D, n_loc = L_sub.shape
    L_flat = L_sub.reshape(n_sub, n_D * n_loc)      # view, not a copy
    prior_flat = prior_dirichlet_param.reshape(-1)
    exp_log_tau = np.ones(n_D * n_loc, dtype=L_sub.dtype)
    # `weights` replaces the per-subtrack jump counts in the sufficient statistic.
    # weights=None keeps saSPT's own definition (occupancy = fraction of
    # DISPLACEMENTS, so a 10-jump subtrack counts 10x a 1-jump one). Passing
    # ones() makes it fraction of SUBTRACKS instead -- equal weighting.
    w_vec = jumps_sub if weights is None else np.asarray(
        weights, dtype=jumps_sub.dtype)
    norm = L_flat @ exp_log_tau                      # iteration-0 norm == 1
    norm[norm == 0] = 1.0
    n = (w_vec / norm) @ L_flat                       # matches max_iter=0 fallback
    n_prev = None
    for _ in range(max_iter):
        if n_prev is not None:
            denom = np.maximum(np.abs(n_prev), 1.0)
            if np.max(np.abs(n - n_prev) / denom) < tol:
                break
        n_prev = n
        m = n + prior_flat
        exp_log_tau = np.exp(digamma(m))
        norm = L_flat @ exp_log_tau
        norm[norm == 0] = 1.0
        w = w_vec / norm
        n_raw = w @ L_flat
        n = exp_log_tau * n_raw
    exp_log_tau_2d = exp_log_tau.reshape(n_D, n_loc)
    n = n.reshape(n_D, n_loc)
    # apply_defoc=False leaves the occupancy uncorrected for the fact that fast
    # states exit the focal slab and are under-observed (the correction divides by
    # f_remain_rbm(D, 1 frame, dt, focal_depth) -- up to ~x23 at D=30 when
    # focal_depth=0.1). Only ever affects post_mean/pmf, never exp_log_tau, so the
    # per-subtrack marginals callers build downstream are unaffected either way.
    if apply_defoc:
        n = likelihood.correct_for_defocalization(n, normalize=False)
    total = n.sum()
    post_mean = n / total if total > 0 else prior_dirichlet_param / prior_dirichlet_param.sum()
    return n, post_mean, exp_log_tau_2d


def _make_likelihood_from_cache(cache):
    import saspt
    grid_snapshot = json.loads(str(cache['grid_snapshot']))
    kwargs = dict(grid_snapshot)
    kwargs['diff_coefs'] = np.asarray(kwargs['diff_coefs'])
    kwargs['loc_errors'] = np.asarray(kwargs['loc_errors'])
    likelihood_type = kwargs.pop('likelihood_type')
    kwargs.pop('splitsize', None)
    return saspt.make_likelihood(likelihood_type, **kwargs)


def _w1_logD(log_grid, pmf_a, pmf_b):
    Fa = np.cumsum(pmf_a); Fa = Fa / Fa[-1]
    Fb = np.cumsum(pmf_b); Fb = Fb / Fb[-1]
    trapz = getattr(np, 'trapezoid', getattr(np, 'trapz'))
    return float(trapz(np.abs(Fa - Fb), log_grid))


def _group_posterior_full(cache, likelihood, fovs, conc_param, max_iter,
                          equal_weight=False, apply_defoc=True):
    """One VBEM fit per group/pseudo-group, returning everything derived
    from it: the population-level soft posterior curve (pmf, what Fork B
    has always used for the global W1 test) AND a per-FOV MLE_D-under-
    POSTERIOR histogram (occ, n_sub) -- the posterior-level analog of
    compare_groups' representation='mle_d'.

    The per-FOV breakdown needs the population fit's converged
    exp_log_tau (unlike the naive/no-refit MLE_D, a FOV's posterior-MLE_D
    histogram is NOT fixed independent of which pool it was fit in -- it
    depends on the whole group's shared posterior "tilt"). Getting it is
    still cheap: one extra pass over L_sub per fit (a handful of
    elementwise (n_sub, n_D) ops, not materializing the full (n_sub, n_D,
    n_loc) array), not a second VBEM fit -- so this doesn't change Fork
    B's fit count (still 2 observed + 2*n_perm permutation fits), only
    adds a small constant per fit.
    """
    # float32 throughout: occupancy fractions don't need float64 precision,
    # and scipy.special.digamma accepts float32 directly.
    mask = np.isin(cache['fov_of_subtrack'], fovs)
    L_sub = cache['L'][mask]  # already float32 in the cache
    jumps_sub = cache['jumps_per_track'][mask].astype(np.float32)
    prior = np.full(L_sub.shape[1:], conc_param, dtype=np.float32)
    # equal_weight -> occupancy = fraction of SUBTRACKS, not of displacements;
    # apply_defoc=False -> no focal-slab correction. Both only affect post_mean
    # (i.e. pmf/W1); the per-FOV MLE_D histogram built below is unaffected either
    # way, since it comes from the marginals, not from post_mean.
    _w = np.ones_like(jumps_sub) if equal_weight else None
    _, post_mean, exp_log_tau_2d = _vbem_infer(
        L_sub, jumps_sub, prior, likelihood, max_iter, weights=_w, apply_defoc=apply_defoc)
    pmf = post_mean.sum(axis=1)   # marginalize loc_error -> pmf over D

    n_sub_local, n_D, n_loc = L_sub.shape
    L_flat = L_sub.reshape(n_sub_local, n_D * n_loc)
    norm = L_flat @ exp_log_tau_2d.reshape(-1)
    norm[norm == 0] = 1.0
    marginal = np.zeros((n_sub_local, n_D), dtype=np.float32)
    for j in range(n_loc):
        marginal += L_sub[:, :, j] * exp_log_tau_2d[None, :, j]
    marginal /= norm[:, None]
    mle_d_idx = np.argmax(marginal, axis=1)
    one_hot = np.zeros((n_sub_local, n_D), dtype=np.float64)
    one_hot[np.arange(n_sub_local), mle_d_idx] = 1.0

    fov_local = cache['fov_of_subtrack'][mask]
    occ = np.full((len(fovs), n_D), np.nan)
    n_sub = np.zeros(len(fovs), dtype=int)
    for i, fov in enumerate(fovs):
        rows = fov_local == fov
        n = int(rows.sum())
        n_sub[i] = n
        if n > 0:
            occ[i] = one_hot[rows].mean(axis=0)
    return pmf, occ, n_sub


def _resolve_permutation_split(all_fovs, n_a, seed, n_perm, k):
    """The FOV relabeling for permutation index k, out of n_perm total,
    under `seed` -- deterministic and INDEPENDENT of execution order or
    of any other index, via numpy's SeedSequence.spawn (keyed only by
    (seed, k)). This is what makes refit_permutation's null loop
    embarrassingly parallel: a worker handling only k in [50, 60) gets
    byte-identical draws to running the full n_perm serially with the
    same seed -- batching changes wall-clock only, never the statistics.
    `n_perm` must match across all workers computing splits for the same
    `seed` (it's part of what SeedSequence.spawn's key derivation uses),
    which is why every batch/array worker is always told the FULL n_perm,
    not just its own batch size.
    """
    child = np.random.SeedSequence(seed).spawn(n_perm)[k]
    rng_k = np.random.default_rng(child)
    perm = rng_k.permutation(len(all_fovs))
    return all_fovs[perm[:n_a]], all_fovs[perm[n_a:]]


def refit_permutation_null_batch(cache, group_a, group_b, k_indices, seed,
                                  n_perm, conc_param=1.0, max_iter=200,
                                  out_path=None, verbose=True,
                                  equal_weight=False, apply_defoc=True):
    """
    Compute the null contribution for a SUBSET of permutation indices
    (`k_indices`) -- e.g. one batch out of many for a SLURM array. Pure
    function of (cache, group_a, group_b, seed, n_perm, k) per index, so
    any worker can run any subset of indices independently and the
    results concatenate into exactly what a serial `refit_permutation`
    call with the same seed/n_perm would have produced (see
    `_resolve_permutation_split`).

    This is the expensive part of Fork B (2 VBEM fits per index) -- the
    thing to actually distribute across an array. `refit_permutation`
    itself calls this once with k_indices=range(n_perm) for the
    interactive/serial case; for a parallel run, call it once per batch
    (each with a disjoint slice of range(n_perm)) and combine with
    `refit_permutation_finalize`.

    Returns
    -------
    dict with keys: k_indices, null (W1 per index), null_curves
    (per-D difference per index), M_mled_perm (posterior-MLE_D per-FOV
    max-T stat per index). Also written to `out_path` (.npz) if given.
    """
    import time
    cache = _load_cache_if_path(cache)
    likelihood = _make_likelihood_from_cache(cache)
    diff_coefs = cache['diff_coefs']
    log_grid = np.log10(diff_coefs)
    n_D = len(diff_coefs)

    fovs_a = _resolve_group(cache, group_a)
    fovs_b = _resolve_group(cache, group_b)
    all_fovs = np.array(list(fovs_a) + list(fovs_b))
    n_a = len(fovs_a)

    k_indices = np.asarray(list(k_indices))
    null = np.empty(len(k_indices))
    null_curves = np.empty((len(k_indices), n_D))
    M_mled_perm = np.empty(len(k_indices))
    t0 = time.time()
    for i, k in enumerate(k_indices):
        pseudo_a, pseudo_b = _resolve_permutation_split(all_fovs, n_a, seed, n_perm, int(k))
        pa, occ_pa_mled, _ = _group_posterior_full(
            cache, likelihood, pseudo_a, conc_param, max_iter,
            equal_weight=equal_weight, apply_defoc=apply_defoc)
        pb, occ_pb_mled, _ = _group_posterior_full(
            cache, likelihood, pseudo_b, conc_param, max_iter,
            equal_weight=equal_weight, apply_defoc=apply_defoc)
        null[i] = _w1_logD(log_grid, pa, pb)
        null_curves[i] = pb - pa
        _, _, delta_star_mled_perm = _studentize(occ_pa_mled, occ_pb_mled)
        M_mled_perm[i] = np.nanmax(np.abs(delta_star_mled_perm))
    elapsed = time.time() - t0

    if verbose:
        print(f"[refit_permutation_null_batch] {len(k_indices)} indices "
              f"({k_indices.min()}..{k_indices.max()}) took {elapsed:.1f}s "
              f"({elapsed/max(len(k_indices),1):.2f}s/perm)")

    out = {'k_indices': k_indices, 'null': null, 'null_curves': null_curves,
           'M_mled_perm': M_mled_perm, 'elapsed_sec': elapsed}
    if out_path:
        np.savez_compressed(out_path, **out)
        if verbose:
            print(f"wrote {out_path}")
    return out


# ----------------------------------------------------------------------------
# Display normalization: posterior-occupation PMF -> probability density
# ----------------------------------------------------------------------------
# The VBEM-refit posterior is a probability *mass* function over the discrete,
# log-spaced diffusion-coefficient grid (sums to 1 across the ~150 grid points).
# For SHAPE overlays it reads more naturally as a probability *density* per unit
# log10(D): f = p / Δlog10D, so that  ∫ f d(log10 D) = 1.  This is a single
# constant rescaling (by 1/Δlog10D ≈ 30 for this grid) applied only at plot
# time -- W1, region-occupancy masses, and the permutation nulls are computed
# from the PMF elsewhere and are unaffected (W1 renormalizes; a constant scale
# cancels in exceedance tests).
def _logD_bin_width(diff_coefs):
    """Uniform log10(D) grid spacing (decade width represented per grid point)."""
    logd = np.log10(np.asarray(diff_coefs, dtype=float))
    return float(np.mean(np.diff(logd)))


def _to_logD_density(pmf, diff_coefs):
    """Convert a posterior-occupation PMF (Σ=1 over the grid) to a probability
    density per unit log10(D): f = p / Δlog10D, with ∫ f d(log10 D) = 1."""
    return np.asarray(pmf, dtype=float) / _logD_bin_width(diff_coefs)


_DENSITY_YLABEL = r'VBEM-refit posterior density (per unit log$_{10}$D)'
_DENSITY_DELTA_YLABEL = r'$\Delta$ posterior density (per unit log$_{10}$D)'


def plot_refit_posterior_difference(diff_coefs, pmf_a, pmf_b, name_a, name_b,
                                     w1_obs, out_pdf=None):
    """Fork-B analog of plot_occupancy_difference: the two groups'
    independently VBEM-refit posterior curves (not naive), shown as a
    probability density per unit log10(D), and their difference vs D."""
    pmf_a = _to_logD_density(pmf_a, diff_coefs)
    pmf_b = _to_logD_density(pmf_b, diff_coefs)
    fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    axes[0].plot(diff_coefs, pmf_a, color='#3366cc', lw=1.8, label=name_a)
    axes[0].plot(diff_coefs, pmf_b, color='#cc3333', lw=1.8, label=name_b)
    axes[0].set_ylabel(_DENSITY_YLABEL)
    axes[0].legend(loc='upper right', fontsize=8)
    axes[0].set_title(f'refit posteriors: {name_a} vs {name_b} (W1={w1_obs:.4f})')

    axes[1].axhline(0, color='k', lw=0.8)
    axes[1].plot(diff_coefs, pmf_b - pmf_a, color='#333333', lw=1.5)
    axes[1].set_ylabel(_DENSITY_DELTA_YLABEL)
    axes[1].set_title(f'difference (b_minus_a: {name_b} - {name_a})')
    _style_logD_axis(axes[1])
    fig.tight_layout()
    if out_pdf:
        fig.savefig(out_pdf)
    return fig


# Shade ramps for single-hue condition series, light -> dark. Same hues as
# plot_genotype_drug_overlay (grey = control/NT background, orange = sgUBA5
# background) so a shaded series reads consistently across figures.
# Wider than plot_genotype_drug_overlay's ramps: with three curves in ONE
# panel (rather than three per panel across two hues) the light and mid tones
# have to separate on their own, without a hue cue to fall back on.
SHADES_GREY = ('#bdbdbd', '#737373', '#1a1a1a')
SHADES_ORANGE = ('#fdae6b', '#e6550d', '#7f2704')


def plot_refit_null_distribution(null, w1_obs, p95, name_a, name_b, n_perm,
                                  seed, out_pdf=None):
    """Histogram of the FOV-relabeling null W1 distribution with the
    observed W1 and its 95th-percentile threshold marked."""
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(null, bins=30, color='gray', alpha=0.6,
            label=f'relabeling null (n_perm={n_perm}, seed={seed})')
    ax.axvline(p95, color='k', ls='--', lw=1.2, label=f'null p95={p95:.4f}')
    ax.axvline(w1_obs, color='#cc3333', lw=2,
               label=f'observed W1={w1_obs:.4f}')
    ax.set_xlabel(r'W1(log$_{10}$ D)')
    ax.set_ylabel('count')
    ax.set_title(f'Fork B relabeling null: {name_a} vs {name_b}')
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout()
    if out_pdf:
        fig.savefig(out_pdf)
    return fig


def plot_refit_permutation_threshold(diff_coefs, delta_obs, s_null, c95,
                                      name_a, name_b, n_perm, seed,
                                      out_pdf=None):
    """Per-D analog of plot_permutation_threshold, but for Fork B: the
    observed VBEM-refit posterior difference curve against a band built
    from the SAME n_perm FOV-relabeling null curves already computed for
    the global W1 test (pointwise null SD per D bin, studentized max-T
    threshold across D) -- i.e. per-diffusion-state significance, not
    just the single aggregate W1 number."""
    band = c95 * s_null
    # display in probability-density units (per unit log10 D); the constant
    # rescale cancels in the exceedance test, so significance is unchanged.
    du = _logD_bin_width(diff_coefs)
    delta_obs = np.asarray(delta_obs, dtype=float) / du
    band = np.asarray(band, dtype=float) / du
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axhline(0, color='k', lw=0.8)
    ax.fill_between(diff_coefs, -band, band, color='gray', alpha=0.3,
                     label=f'relabeling-null floor (n_perm={n_perm}, seed={seed})')
    exceeds = np.abs(delta_obs) > band
    ax.plot(diff_coefs, delta_obs, color='#333333', lw=1.3, label=r'$\Delta$(D)')
    if exceeds.any():
        ax.plot(diff_coefs[exceeds], delta_obs[exceeds], 'o', color='#cc3333',
                ms=3, label='exceeds floor', zorder=5)
    _style_logD_axis(ax)
    ax.set_ylabel(_DENSITY_DELTA_YLABEL)
    ax.set_title(f'Fork B per-D permutation threshold: {name_a} vs {name_b}')
    ax.legend(loc='upper right', fontsize=7)
    fig.tight_layout()
    if out_pdf:
        fig.savefig(out_pdf)
    return fig


def _compute_observed_fits(cache, likelihood, fovs_a, fovs_b, conc_param, max_iter,
                            equal_weight=False, apply_defoc=True):
    """The two REAL-group VBEM fits (not permutation nulls) -- shared by
    `refit_permutation` (serial) and the array aggregator."""
    import time
    t0 = time.time()
    pmf_a, occ_a_mled, n_sub_a_mled = _group_posterior_full(
        cache, likelihood, fovs_a, conc_param, max_iter,
        equal_weight=equal_weight, apply_defoc=apply_defoc)
    pmf_b, occ_b_mled, n_sub_b_mled = _group_posterior_full(
        cache, likelihood, fovs_b, conc_param, max_iter,
        equal_weight=equal_weight, apply_defoc=apply_defoc)
    per_fit_sec = (time.time() - t0) / 2.0
    delta_mled_obs, se_mled_obs, _ = _studentize(occ_a_mled, occ_b_mled)
    return {
        'pmf_a': pmf_a, 'pmf_b': pmf_b,
        'occ_a_mled': occ_a_mled, 'occ_b_mled': occ_b_mled,
        'n_sub_a_mled': n_sub_a_mled, 'n_sub_b_mled': n_sub_b_mled,
        'delta_mled_obs': delta_mled_obs, 'se_mled_obs': se_mled_obs,
        'per_fit_sec': per_fit_sec,
    }


def refit_permutation_finalize(cache, name_a, name_b, out_dir, n_perm, seed,
                                observed, null, null_curves, M_mled_perm,
                                per_fit_sec=None, total_perm_sec=None,
                                conc_param=None, max_iter=None, verbose=True):
    """
    Statistics + plots + save, given the observed fits (from
    `_compute_observed_fits`) and the FULL null arrays (concatenated
    across however many batches computed them -- see
    `refit_permutation_null_batch`). This is the shared tail between the
    serial `refit_permutation` and the array-aggregation path; splitting
    the expensive part (null fits) from this cheap part (pure numpy,
    seconds) is what makes the array path possible.

    Two per-D statistics computed here, both against the SAME null:
      1. Soft posterior: pmf_b(d)-pmf_a(d) against a max-T band
         studentized by the null curves' own pointwise SD (no per-FOV
         replicate at the pooled-curve level).
      2. Posterior MLE_D (per FOV): studentized by ACROSS-FOV variance,
         like compare_groups' own tests -- real per-FOV replicates exist
         here. Generally the more conservative of the two (see
         docs/HANDOFF or project memory for a case where it found much
         less signal than #1 despite visually striking mean curves).
    """
    diff_coefs = cache['diff_coefs']
    log_grid = np.log10(diff_coefs)
    pmf_a, pmf_b = observed['pmf_a'], observed['pmf_b']
    occ_a_mled, occ_b_mled = observed['occ_a_mled'], observed['occ_b_mled']
    n_sub_a_mled, n_sub_b_mled = observed['n_sub_a_mled'], observed['n_sub_b_mled']
    delta_mled_obs, se_mled_obs = observed['delta_mled_obs'], observed['se_mled_obs']
    per_fit_sec = per_fit_sec if per_fit_sec is not None else observed.get('per_fit_sec')

    os.makedirs(out_dir, exist_ok=True)
    w1_obs = _w1_logD(log_grid, pmf_a, pmf_b)
    delta_obs = pmf_b - pmf_a

    p95 = float(np.percentile(null, 95))
    p_value = float((null >= w1_obs).mean())
    verdict = ('within relabeling null' if w1_obs <= p95 else
               'separation IS above the relabeling null')

    # Per-D max-T band (soft posterior): studentize each null curve by
    # the null's OWN pointwise SD (there's no per-FOV replicate to draw a
    # variance estimate from at the pooled-curve level).
    s_null = null_curves.std(axis=0, ddof=1)
    T_null = null_curves / (s_null[None, :] + _EPS)
    M_perm = np.nanmax(np.abs(T_null), axis=1)
    c95_perD = float(np.quantile(M_perm, 0.95))
    band = c95_perD * s_null
    exceeds = np.abs(delta_obs) > band
    excursion_intervals = _excursion_intervals(diff_coefs, exceeds, delta_obs)

    # Per-D max-T band (posterior MLE_D, per FOV): c95 from the
    # permutation max-stats (each already studentized by ITS OWN
    # within-permutation across-FOV SD); band uses the REAL data's own
    # se for the final threshold -- same convention as compare_groups.
    c95_mled = float(np.quantile(M_mled_perm, 0.95))
    band_mled = c95_mled * se_mled_obs
    exceeds_mled = np.abs(delta_mled_obs) > band_mled
    excursion_intervals_mled = _excursion_intervals(diff_coefs, exceeds_mled, delta_mled_obs)

    result = {
        'name_a': name_a, 'name_b': name_b, 'w1_observed': w1_obs,
        'null_p95': p95, 'null_mean': float(null.mean()), 'p_value': p_value,
        'verdict': verdict, 'n_perm': n_perm, 'seed': seed,
        'conc_param': conc_param, 'max_iter': max_iter,
        'timing_per_fit_sec': per_fit_sec,
        'timing_total_permutation_sec': total_perm_sec,
        'n_fov_a': int(occ_a_mled.shape[0]), 'n_fov_b': int(occ_b_mled.shape[0]),
        'per_D_c95': c95_perD, 'excursion_intervals': excursion_intervals,
        'mle_d_per_D_c95': c95_mled,
        'mle_d_excursion_intervals': excursion_intervals_mled,
    }
    tag = f"{name_a}_vs_{name_b}"
    with open(os.path.join(out_dir, f"refit_permutation_{tag}.json"), 'w') as f:
        json.dump({k: v for k, v in result.items()}, f, indent=2, default=_json_default)

    # Save the raw curves/null too -- these are what plotting needs, and
    # the null array in particular is expensive (the whole point of this
    # function) so it should never have to be silently discarded.
    np.savez_compressed(
        os.path.join(out_dir, f"refit_permutation_{tag}.npz"),
        diff_coefs=diff_coefs, pmf_a=pmf_a, pmf_b=pmf_b, null=null,
        null_curves=null_curves, delta_obs=delta_obs, band=band, s_null=s_null,
        occ_a_mled=occ_a_mled, occ_b_mled=occ_b_mled,
        n_sub_a_mled=n_sub_a_mled, n_sub_b_mled=n_sub_b_mled,
        delta_mled_obs=delta_mled_obs, se_mled_obs=se_mled_obs,
        band_mled=band_mled, M_mled_perm=M_mled_perm,
    )

    plot_refit_permutation_threshold(
        diff_coefs, delta_obs, s_null, c95_perD, name_a, name_b, n_perm, seed,
        out_pdf=os.path.join(out_dir, f"refit_permutation_threshold_{tag}.pdf"))
    plot_refit_posterior_difference(
        diff_coefs, pmf_a, pmf_b, name_a, name_b, w1_obs,
        out_pdf=os.path.join(out_dir, f"refit_posterior_difference_{tag}.pdf"))
    plot_refit_null_distribution(
        null, w1_obs, p95, name_a, name_b, n_perm, seed,
        out_pdf=os.path.join(out_dir, f"refit_null_distribution_{tag}.pdf"))

    mled_y_label = 'MLE_D fraction (posterior, per FOV)'
    plot_naive_occupancy(
        diff_coefs, occ_a_mled, occ_b_mled, n_sub_a_mled, n_sub_b_mled,
        name_a, name_b, y_label=mled_y_label,
        out_pdf=os.path.join(out_dir, f"posterior_mle_d_occupancy_{tag}.pdf"))
    plot_occupancy_difference(
        diff_coefs, delta_mled_obs, name_a, name_b,
        out_pdf=os.path.join(out_dir, f"posterior_mle_d_difference_{tag}.pdf"))
    plot_permutation_threshold(
        diff_coefs, delta_mled_obs, se_mled_obs, c95_mled, name_a, name_b,
        n_perm, seed,
        out_pdf=os.path.join(out_dir, f"posterior_mle_d_threshold_{tag}.pdf"))
    plot_combined(
        diff_coefs, occ_a_mled, occ_b_mled, n_sub_a_mled, n_sub_b_mled,
        delta_mled_obs, se_mled_obs, c95_mled, name_a, name_b, n_perm, seed,
        y_label=mled_y_label,
        out_pdf=os.path.join(out_dir, f"posterior_mle_d_combined_{tag}.pdf"))
    plt.close('all')

    if verbose:
        print(f"[refit_permutation {name_a} vs {name_b}] observed W1="
              f"{w1_obs:.4f}, null p95={p95:.4f}, p={p_value:.4f} -> {verdict}")
        if per_fit_sec is not None and total_perm_sec is not None:
            print(f"  ~{per_fit_sec:.2f}s per VBEM fit; {n_perm} perms took "
                  f"{total_perm_sec:.1f}s ({total_perm_sec/max(n_perm,1):.2f}s/perm)")
        perD_verdict = ("no excursion (per-D soft-posterior difference stays "
                         "inside the relabeling floor everywhere)" if not
                         excursion_intervals else f"{len(excursion_intervals)} "
                         f"per-D excursion interval(s) exceed the floor")
        print(f"  per-D (soft posterior): {perD_verdict}")
        mled_verdict = ("no excursion (posterior-MLE_D per-FOV difference "
                         "stays inside the floor everywhere)" if not
                         excursion_intervals_mled else
                         f"{len(excursion_intervals_mled)} per-D excursion "
                         f"interval(s) exceed the floor")
        print(f"  per-D (posterior MLE_D, per FOV): {mled_verdict}")

    result['pmf_a'] = pmf_a
    result['pmf_b'] = pmf_b
    result['null'] = null
    result['diff_coefs'] = diff_coefs
    result['occ_a_mled'] = occ_a_mled
    result['occ_b_mled'] = occ_b_mled
    return result


def refit_permutation(cache, group_a, group_b, name_a="A", name_b="B",
                       out_dir=".", n_perm=200, seed=0, conc_param=1.0,
                       max_iter=200, verbose=True):
    """
    Fork B -- tests whether per-condition-fit sharpness (a VBEM posterior,
    not the naive prior-free curve) is real vs prior-manufactured. Only
    run this when a specific comparison needs to defend that claim -- not
    routinely (handoff section 5).

    Serial convenience wrapper: computes the 2 observed fits, then ALL
    n_perm permutation fits in one process (`refit_permutation_null_batch`
    with k_indices=range(n_perm)), then `refit_permutation_finalize`. For
    n_perm large enough that this is slow, run the null fits as a SLURM
    array instead (`refit_permutation_null_batch` per batch, one array
    task each) and call `refit_permutation_finalize` once at the end --
    see run_refit_permutation_batch.py / run_refit_permutation_aggregate.py
    / call_refit_permutation_array.slurm. Both paths give byte-identical
    results for the same (seed, n_perm): permutation index k's FOV split
    depends only on (seed, k), never on execution order (see
    `_resolve_permutation_split`), and this function always tells the
    batch worker the FULL n_perm even though it's requesting every index
    -- exactly what a parallel caller would do per-batch.

    Returns
    -------
    See `refit_permutation_finalize`.
    """
    cache = _load_cache_if_path(cache)
    likelihood = _make_likelihood_from_cache(cache)
    fovs_a = _resolve_group(cache, group_a)
    fovs_b = _resolve_group(cache, group_b)
    _validate_groups(cache, fovs_a, fovs_b, negative_control=False)

    observed = _compute_observed_fits(cache, likelihood, fovs_a, fovs_b,
                                       conc_param, max_iter)
    null_result = refit_permutation_null_batch(
        cache, group_a, group_b, range(n_perm), seed, n_perm,
        conc_param=conc_param, max_iter=max_iter, verbose=False)

    return refit_permutation_finalize(
        cache, name_a, name_b, out_dir, n_perm, seed, observed,
        null_result['null'], null_result['null_curves'],
        null_result['M_mled_perm'], per_fit_sec=observed['per_fit_sec'],
        total_perm_sec=null_result['elapsed_sec'], conc_param=conc_param,
        max_iter=max_iter, verbose=verbose)
