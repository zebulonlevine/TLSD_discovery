"""
multiplicity.py
===============

Shared multiple-comparison helpers for the SPT panels.

Exists because seven near-identical `holm()` implementations had accumulated
across the HMM/rate and topology scripts. The canonical one here is lifted
from spt_compare._holm -- the only one of the seven that is NaN-safe, which
matters because region-effect sidecars can carry a NaN p for a region with
no mass.

Deliberately numpy-only. It does NOT import occupancy_compare (matplotlib +
saspt_likelihood_cache), so the rate/topology scripts can use it without
pulling in the occupancy engine.

PERMUTATION P-VALUES: use perm_p_addone() before correcting. A permutation
test cannot produce a p of exactly 0 -- `(np.abs(null) >= abs(obs)).mean()`
returning 0.0 means "no draw exceeded the observation", which bounds p below
1/n_perm, it does not establish p=0. Holm propagates a literal 0 as an
adjusted 0, which overclaims. The Phipson & Smyth (2010) add-one estimator
(1 + k) / (1 + n_perm) is the unbiased fix and floors p at 1/(1+n_perm).
"""

import numpy as np


def holm(p_values):
    """Holm-Bonferroni step-down adjusted p-values, aligned with input order.

    NaN-in produces NaN-out, and NaNs are excluded from the family size
    (a region with no testable mass should not inflate the penalty on the
    regions that do)."""
    p_values = np.asarray(p_values, dtype=float)
    out = np.full(len(p_values), np.nan)
    valid = ~np.isnan(p_values)
    if not valid.any():
        return out
    idx_valid = np.where(valid)[0]
    sorted_local = np.argsort(p_values[idx_valid])
    n = len(idx_valid)
    running_max = 0.0
    for rank, local_i in enumerate(sorted_local):
        i = idx_valid[local_i]
        adj = min(p_values[i] * (n - rank), 1.0)
        running_max = max(running_max, adj)
        out[i] = running_max
    return out


def bh(p_values):
    """Benjamini-Hochberg FDR-adjusted p-values, aligned with input order.
    Reported alongside Holm for reference; Holm is the default because FWER
    is the conservative choice at the small family sizes used here (m ~ 5-6)
    and because the failure mode to guard against in an interaction test is a
    false positive, not a missed discovery."""
    p_values = np.asarray(p_values, dtype=float)
    out = np.full(len(p_values), np.nan)
    valid = ~np.isnan(p_values)
    if not valid.any():
        return out
    idx_valid = np.where(valid)[0]
    p = p_values[idx_valid]
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    running_min = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        running_min = min(running_min, p[i] * n / (rank + 1))
        adj[i] = min(running_min, 1.0)
    out[idx_valid] = adj
    return out


def perm_p_addone(p_raw, n_perm):
    """Phipson-Smyth add-one permutation p from a raw exceedance fraction.

    p_raw is the fraction of permutation draws at least as extreme as the
    observation (what _subband_summary / did_from_occ report). Recovers the
    integer exceedance count k = round(p_raw * n_perm) and returns
    (1 + k) / (1 + n_perm), which floors at 1/(1 + n_perm) instead of 0."""
    p_raw = np.asarray(p_raw, dtype=float)
    k = np.rint(p_raw * n_perm)
    out = (1.0 + k) / (1.0 + n_perm)
    return np.where(np.isnan(p_raw), np.nan, out)


def stars(p, alpha=0.05):
    """Significance marker for an already-adjusted p. '' if not significant
    at alpha (or NaN)."""
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ''
    if p < 1e-3:
        return '***'
    if p < 1e-2:
        return '**'
    if p < alpha:
        return '*'
    return ''
