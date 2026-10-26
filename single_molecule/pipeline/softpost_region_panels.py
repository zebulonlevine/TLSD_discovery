"""
softpost_region_panels.py
===========================

SOFT-posterior region-effect panels -- the posterior-occupancy analog of
forkb_region_panels.region_effect_from_forkb (which uses the HARD per-
trajectory posterior-MLE_D argmax). Here each trajectory contributes its
full VBEM posterior over D (soft), averaged per FOV, so the per-region
value is a genuine difference in posterior occupancy MASS (not a fraction
of argmax-classified trajectories).

Kept SEPARATE from the MLE_D panels (distinct `_softpost` filenames) so
both representations coexist:
  - region_effect_forkb_*.pdf     : hard posterior-MLE_D (argmax) fractions
  - region_effect_softpost_*.pdf  : soft posterior occupancy mass

The saved refit_permutation npz store only the population soft curve
(pmf_a/pmf_b) and the per-FOV HARD histograms (occ_*_mled) -- not per-FOV
soft occupancy -- so this recomputes the observed VBEM fit per group (one
fit each; NOT the permutation refits). The region-effect null is the same
FOV-block relabel + bootstrap of the observed per-FOV vectors as the
MLE_D panels (occupancy_compare._subband_summary), applied to soft vectors.
"""

import numpy as np

from occupancy_compare import (
    _vbem_infer, _make_likelihood_from_cache, _resolve_group, _subband_summary,
    _w1_logD, plot_region_effect_vs_null, load_likelihood_cache,
)


def group_soft_per_fov(cache, likelihood, fovs, conc_param=1.0, max_iter=200, equal_weight=False):
    """Per-FOV SOFT posterior occupancy (n_fov, n_D) for one group: one VBEM
    fit over the group's subtracks, then per-FOV mean of the per-subtrack
    soft posterior marginal over D (each row sums to 1). The soft analog of
    occupancy_compare._group_posterior_full's one-hot occ.
    equal_weight=True -> fit weights each subtrack equally (fraction of subtracks)
    instead of by displacement count; only the shared exp_log_tau fit changes
    (per-FOV aggregation is already unweighted).

    TWO PROPERTIES WORTH NOT RE-DERIVING:

    (1) INDEPENDENT OF DEFOCALIZATION. _vbem_infer applies the defocalization
        correction after the VBEM loop, to `n` only, and this function uses the
        third return value (exp_log_tau_2d), discarding n/post_mean. So the
        per-FOV soft occupancy returned here is identical whether apply_defoc is
        True or False -- which is why it can be gated against an apply_defoc=False
        permutation null without a mismatch. There is no apply_defoc parameter
        here because one would be provably dead code.

    (2) equal_weight ONLY RETILTS THE SHARED FIT. Aggregation is already
        unweighted at both levels: each subtrack contributes its own marginal,
        and each FOV is a plain mean over its subtracks. So equal_weight=True
        makes the fit's weighting agree with the aggregation already in use,
        which is why the FOV-mean estimator and the pooled-population estimator
        converge under it."""
    mask = np.isin(cache['fov_of_subtrack'], fovs)
    L_sub = cache['L'][mask]
    jumps = cache['jumps_per_track'][mask].astype(np.float32)
    prior = np.full(L_sub.shape[1:], conc_param, dtype=np.float32)
    _w = np.ones_like(jumps) if equal_weight else None
    _, _, exp_log_tau_2d = _vbem_infer(L_sub, jumps, prior, likelihood, max_iter, weights=_w)

    n_sub, n_D, n_loc = L_sub.shape
    L_flat = L_sub.reshape(n_sub, n_D * n_loc)
    norm = L_flat @ exp_log_tau_2d.reshape(-1)
    norm[norm == 0] = 1.0
    marginal = np.zeros((n_sub, n_D), dtype=np.float64)
    for j in range(n_loc):
        marginal += L_sub[:, :, j] * exp_log_tau_2d[None, :, j]
    marginal /= norm[:, None]                      # per-subtrack SOFT posterior over D

    fov_local = cache['fov_of_subtrack'][mask]
    occ = np.full((len(fovs), n_D), np.nan)
    n_sub_fov = np.zeros(len(fovs), dtype=int)
    for i, f in enumerate(fovs):
        rows = fov_local == f
        n_sub_fov[i] = int(rows.sum())
        if rows.any():
            occ[i] = marginal[rows].mean(axis=0)
    return occ, n_sub_fov


def did_from_occ(diff, oa_n, oa_d, ob_n, ob_d, regions, n_perm=2000,
                 n_boot=2000, seed=0):
    """Per-region difference-in-differences [a_drug-a_none]-[b_drug-b_none]
    from four PRECOMPUTED per-FOV soft-occupancy matrices (no fitting) --
    lets a driver fit each group's soft occ once and reuse it across many
    DiDs. Mirrors forkb_region_panels.genotype_drug_interaction_from_forkb."""
    rng = np.random.default_rng(seed)
    out = []
    for label, (lo, hi) in regions.items():
        m = (diff >= lo) & (diff < hi)
        sa_d = oa_d[:, m].sum(1); sa_n = oa_n[:, m].sum(1)
        sb_d = ob_d[:, m].sum(1); sb_n = ob_n[:, m].sum(1)
        eff_a = float(np.nanmean(sa_d) - np.nanmean(sa_n))
        eff_b = float(np.nanmean(sb_d) - np.nanmean(sb_n))
        obs = eff_a - eff_b
        na_d, na_n, nb_d, nb_n = len(sa_d), len(sa_n), len(sb_d), len(sb_n)
        boot = np.empty(n_boot)
        for i in range(n_boot):
            ea = sa_d[rng.integers(0, na_d, na_d)].mean() - sa_n[rng.integers(0, na_n, na_n)].mean()
            eb = sb_d[rng.integers(0, nb_d, nb_d)].mean() - sb_n[rng.integers(0, nb_n, nb_n)].mean()
            boot[i] = ea - eb
        ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
        pooled_d = np.concatenate([sa_d, sb_d]); pooled_n = np.concatenate([sa_n, sb_n])
        null = np.empty(n_perm)
        for i in range(n_perm):
            pd = rng.permutation(len(pooled_d)); pn = rng.permutation(len(pooled_n))
            null[i] = ((pooled_d[pd[:na_d]].mean() - pooled_n[pn[:na_n]].mean())
                       - (pooled_d[pd[na_d:]].mean() - pooled_n[pn[na_n:]].mean()))
        out.append({'label': label, 'd_lo': lo, 'd_hi': hi, 'observed_diff': obs,
                    'effect_a_drug_minus_none': eff_a, 'effect_b_drug_minus_none': eff_b,
                    'ci_lo': float(ci_lo), 'ci_hi': float(ci_hi), 'n_boot': n_boot,
                    'null_mean': float(null.mean()), 'null_sd': float(null.std()),
                    'null_p2p5': float(np.percentile(null, 2.5)),
                    'null_p97p5': float(np.percentile(null, 97.5)),
                    'p_value_two_sided': float((np.abs(null) >= abs(obs)).mean()), 'n_perm': n_perm})
    return out
