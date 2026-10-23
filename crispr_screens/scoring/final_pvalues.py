#!/usr/bin/env python3
"""Per-sgRNA p-values and gene-level roll-up for the weighted synthesis/degradation score.

Final design and why
--------------------
The statistic tested is the pipeline's own score,
    score = parallel*|parallel| / (|parallel| + perp),
nonlinear in the two LFCs -- so this is not a test on `parallel` alone.

Because the score is homogeneous of degree 1 in (lfc_x, lfc_y) and sequencing
depth is shared between the two comparisons for a given guide, per-guide noise
factors out: score / a_i is comparable across guides, where a_i is the guide's
depth-driven noise scale (from a log-log trend of pooled within-group replicate
residuals -- signal-free, 4 df per guide before shrinkage).

The null for score/a_i is then taken EMPIRICALLY from the non-targeting guides.
Three parametric nulls were tried and rejected:

  null                                 NTC p<0.05 (target 0.05)
  analytic h(u,v), iid equal-variance  0.100 synth / 0.022-0.029 deg
  + single scale fitted on NTCs        0.073-0.084 synth / 0.051 deg
  + measured 2x2 LFC covariance        0.135-0.149 synth / 0.039-0.059 deg

Degradation can be made to calibrate; synthesis cannot, because the NTC null is
genuinely heavy-tailed and skewed (S3/Synthesis kurtosis 14.1, skew 1.85 vs 1.59
and 0.0 for the Gaussian-derived null) -- a handful of outlier non-targeting
guides set the tail. No Gaussian model reproduces that, so the empirical null is
the honest choice.

Cost of the empirical null: the p-value floor is 1/(n_NTC+1) ~ 0.002 with 490
NTC guides. That is fine for the gene-level roll-up, which aggregates guides,
but it means a single sgRNA can never be called more significant than ~2e-3.
"""
import os
import numpy as np
import pandas as pd
import sys
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sgrna_score_pvalues as M
from common import paths  # noqa: E402  (path set up by sgrna_score_pvalues)


def sgrna_table(screen, score, cond='wt'):
    """Per-sgRNA p-values for the PIPELINE's published line score.

    The score tested is the one the pipeline already publishes
    (<scores>/<dir>/*_sgrna_scores.csv), i.e. built from MAGeCK's LFCs. Only the
    noise scale comes from the variance model here.

    An earlier version used a score recomputed from raw counts with its own
    median-of-ratios + log2(x+32) LFCs. That correlated r=0.996-0.998 with the
    pipeline score, but it silently bundled an LFC redefinition into what was
    meant to be a change of gene-level selection only -- and the two effects on
    gene epsilon were comparable in size (selection alone r=0.992-0.996, LFC
    recomputation alone r=0.995-0.998). Using the pipeline score keeps
    Table S4's LFC columns coherent with its score column and isolates the
    change to the p-ranking and matched null.
    """
    import glob as _glob
    t = M.build(screen, score, cond=cond)
    t['a'] = np.sqrt(t['s2_fit'].values)          # depth-driven noise scale
    # parallel/perp from the recomputed LFCs are retained for reference only
    par, perp = t['parallel'].values, t['perp'].values
    t['score_recomputed'] = par * np.abs(par) / (np.abs(par) + perp + 1e-10)

    sdir = f'{screen}_{score}' + ('' if cond == 'wt' else '_sgUBA5')
    # Tolerant glob: the sgUBA5 directories have inconsistent filenames --
    # L29_Degradation_sgUBA5/ omits the _sgUBA5 tag and S3_Synthesis_sgUBA5/ has
    # a doubled suffix ('..._sgrna_scores.csvsgrna_scores.csv'), so a
    # '*_sgrna_scores.csv' pattern silently misses it.
    d = paths.inp('Sublibrary_results/scores', sdir)
    cand = sorted(f for f in _glob.glob(os.path.join(d, '*'))
                  if 'sgrna_scores' in os.path.basename(f) and f.endswith('.csv'))
    if len(cand) != 1:
        raise FileNotFoundError(
            f'expected exactly one sgrna_scores csv in scores/{sdir}, found {cand}')
    pipe = pd.read_csv(cand[0])[['sgrna', 'score']]
    n_before = len(t)
    t = t.merge(pipe, on='sgrna', how='inner')
    if len(t) != n_before:
        print(f'  {sdir}: {n_before} -> {len(t)} guides after joining the published score')
    t['t_stat'] = t['score'] / t['a']

    ntc_abs = np.sort(np.abs(t.loc[t['is_ntc'], 't_stat'].values))
    rank = np.searchsorted(ntc_abs, np.abs(t['t_stat'].values), side='left')
    t['pvalue'] = (len(ntc_abs) - rank + 1) / (len(ntc_abs) + 1)
    t['p_floor'] = 1.0 / (len(ntc_abs) + 1)
    return t


def gene_table(t, top_n=3, fdr=0.05, seed=42):
    """Gene scores from the per-sgRNA p-values, MAGeCK-iNC structure.

    Differences from the current sgrna_line_score_utest():
      1. guides are ranked by p-value and the Mann-Whitney test compares
         p-values to the NTC p-value pool -- ranking statistic == tested
         statistic. Ranking by |score| and then testing signed scores mixes
         directions inside a gene and cancels in the reported mean.
      2. each pseudogene draws as many NTC guides as its matched real gene has,
         not a fixed 5. Guides per gene here is modally 6 (up to 19), so a
         fixed-5 null makes the real top-3 selection more extreme than the
         null's, which is anti-conservative.
    """
    rng = np.random.default_rng(seed)
    ntc, targ = t[t['is_ntc']], t[~t['is_ntc']]
    ntc_p = ntc['pvalue'].values
    ntc_score = ntc['score'].values

    rows = []
    for gene, sub in targ.groupby('gene', sort=False):
        # Rank by p-value, breaking ties on |score| (largest first).
        # The empirical NTC null has a floor of 1/(n_NTC+1) ~ 2e-3, so for a
        # strong gene every guide can sit AT that floor; nsmallest() would then
        # break the tie on row order and epsilon would average three arbitrary
        # guides instead of the three largest-magnitude ones. 19-23% of genes
        # hit a tie at the cut, and up to 116 per screen have all three at the
        # floor.
        s = sub.assign(_absscore=sub['score'].abs()) \
               .sort_values(['pvalue', '_absscore'], ascending=[True, False]) \
               .head(top_n)
        _, pv = stats.mannwhitneyu(s['pvalue'].values, ntc_p, alternative='two-sided')
        rows.append((gene, s['score'].mean(), pv, len(sub)))
    genes = pd.DataFrame(rows, columns=['gene_id', 'epsilon', 'pvalue', 'n_sgrnas'])

    nulls = []
    for i, n in enumerate(genes['n_sgrnas'].values):
        k = min(int(n), len(ntc_p))
        idx = rng.choice(len(ntc_p), size=k, replace=False)
        # Select the SAME guides for epsilon and for the test, exactly as a real
        # gene does. Taking the 3 lowest-p values for the test but the mean of 3
        # arbitrary drawn guides for epsilon under-disperses the null: p-selection
        # favours large |score|, so real genes would get systematically more
        # extreme epsilon than the null, pushing the product threshold too low.
        # identical tiebreak in the null, otherwise the pseudogenes' selection is
        # less extreme than the real genes' and the FDR threshold drifts low
        chosen = idx[np.lexsort((-np.abs(ntc_score[idx]), ntc_p[idx]))[:top_n]]
        _, pv = stats.mannwhitneyu(ntc_p[chosen], ntc_p, alternative='two-sided')
        nulls.append((f'NTC_{i}', ntc_score[chosen].mean(), pv, n))
    genes = pd.concat([genes, pd.DataFrame(nulls, columns=genes.columns)], ignore_index=True)

    is_null = genes['gene_id'].str.startswith('NTC_')
    genes['epsilon'] -= genes.loc[is_null, 'epsilon'].median()
    genes['product'] = genes['epsilon'] * -np.log10(genes['pvalue'].clip(lower=1e-300))

    thr = float('inf')
    for cut in np.arange(0, np.abs(genes['product']).max(), 0.02):
        above = np.abs(genes['product']) > cut
        if above.sum() and (is_null & above).sum() / above.sum() < fdr:
            thr = cut
            break
    genes['is_hit'] = (np.abs(genes['product']) > thr) & ~is_null
    genes.attrs['threshold'] = thr
    return genes
