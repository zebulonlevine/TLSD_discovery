"""
build_figS7DE_panels.py -- Supplementary Figures S7D and S7E.

Dataset `2026_07_08`, tier `q=15p0`, 15.091 ms / 0.11 um. Endogenously tagged
RPL29 (60S); conditions `60S_{NT,UBA5}_{none,HT10,eme360}`. The raw acquisition
labels of this experiment were swapped; the names here are the true identity
(see "Subunit labels" in single_molecule/README.md).

Drugs, added 1 h before imaging: harringtonine (`HT10`, 10 ug/mL, initiation
inhibitor -> run-off to free subunits) and emetine (`eme360`, 360 uM, traps
elongating ribosomes -> polysomes stabilised, free subunits depleted). The
opposite expected directions are what let each diffusion bin be assigned as
translating or free.

Panels
------
  S7D  per-bin occupancy shift on drug treatment, within each genotype
  S7E  the same numbers as a translation-state assignment table. The
       assignment is STATED (see ASSIGNMENT below), identical for both
       genotypes; the two drug responses are reported as a test of it, not as
       its source.

Per the draft legend, comparisons are made only drug-vs-untreated within a
genotype, and genotype-vs-genotype at matched treatment. BH correction is
applied across the bins of each panel.

Run:  sbatch call_build_panels.slurm build_figS7DE_panels.py
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
import paperfig_style as ps
from occupancy_compare import _subband_summary
from multiplicity import bh, holm, perm_p_addone
from region_schemes import SCHEME_60S

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
TIER = tlsd_paths.spt('all_data', '2026_07_08', 'tracking_output_q=15p0')
# merged multi-condition cache of the default (non-36state) layout; with
# PANEL_VARIANT=36state each condition has its own cache and this is unused
CACHE = 'occupancy_compare_cache/60S_did_merged_likelihood_cache.npz'
OUT = tlsd_paths.outp('PaperFigures_updated')
SCHEME = SCHEME_60S

DRUG_LABEL = {'none': 'untreated', 'HT10': 'harringtonine',
              'eme360': 'emetine'}
# Distinct hues for the two drugs because their expected effects are OPPOSITE
# (harringtonine runs ribosomes off, emetine traps them), so the two curves
# should not read as shades of one perturbation.
DRUG_COLOR = {'none': ps.BLACK, 'HT10': ps.PERTURBED, 'eme360': ps.BLUE}
GENO_LABEL = {'NT': 'sgNT', 'UBA5': 'sgUBA5'}
GENO_COLOR = {'NT': ps.CONTROL, 'UBA5': ps.PERTURBED}
MEMBER = '60S_{geno}_{drug}'


def contrast(ctx, occ, a, b, kind, seed=0):
    """Per-bin (b - a) occupancy difference with FOV bootstrap CI, whole-FOV
    relabelling null, and BH across this panel's bins."""
    rows = _subband_summary(ctx['diff_coefs'], occ[a], occ[b], SCHEME,
                            None, None, seed=seed, n_boot=pc.N_BOOT,
                            n_perm=pc.N_PERM)
    df = pd.DataFrame(rows)
    df['p_perm'] = [perm_p_addone(p, pc.N_PERM) for p in df['p_value_two_sided']]
    df['p_bh'] = bh(df['p_perm'].to_numpy())
    df['p_holm'] = holm(df['p_perm'].to_numpy())
    df['stars_bh'] = [ps.stars(p) for p in df['p_bh']]
    df['reference'], df['condition'], df['contrast_kind'] = a, b, kind
    df['n_fov_ref'], df['n_fov_cond'] = occ[a].shape[0], occ[b].shape[0]
    df['equal_weight'] = pc.EQUAL_WEIGHT
    return df


# Translation-state assignment, by bin index, stated in advance and applied
# IDENTICALLY to both genotypes (Zeb, 2026-09-16).
#
# Deriving the call per genotype from the drug responses -- the earlier
# approach -- made the two genotypes disagree about which bin was the
# polysome (sgNT said bin II, sgUBA5 said bin III), which is an artefact of
# reading a label off the same data the label is then used to interpret. The
# assignment is a property of the diffusion regime, so it is fixed here and
# the drug responses are reported as corroboration.
#
# Bin V (> 3 um2/s) is too fast for a 60S subunit; it is the free-protein
# bin and is named as such in SCHEME_60S. No translation state is claimed.
ASSIGNMENT = {
    0: ('polysome', 'Y'),        # I    <= 0.004
    1: ('polysome', 'Y'),        # II   0.004 - 0.025
    2: ('ambiguous', '?'),       # III  0.025 - 0.35
    3: ('free subunit', 'N'),    # IV   0.35 - 3.0
    4: ('free protein', '-'),    # V    > 3.0
}

# Expected sign of (harringtonine, emetine) for each assigned state.
# Harringtonine blocks initiation and lets ribosomes run off, so a
# translating species DECREASES and free subunits ACCUMULATE. Emetine
# freezes elongating ribosomes, so a translating species is STABILISED while
# free subunits are depleted. States with no directional prediction are not
# tested.
PREDICTED_SIGN = {'polysome': (-1, +1), 'free subunit': (+1, -1)}


def check_assignment(bin_index, hht, eme):
    """Do this bin's two drug responses agree with its stated state?

    Returns (prediction, verdict). Bins whose state makes no directional
    prediction return two empty strings and are not counted either way.
    """
    state = ASSIGNMENT[bin_index][0]
    want = PREDICTED_SIGN.get(state)
    if want is None:
        return '', ''
    pred = (f'HHT {"down" if want[0] < 0 else "up"}, '
            f'eme {"up" if want[1] > 0 else "down"}')
    ok = (np.sign(hht) == want[0]) and (np.sign(eme) == want[1])
    return pred, 'consistent' if ok else 'INCONSISTENT'


def main():
    ps.apply()
    occ, n_sub, n_fov = {}, {}, {}
    for geno in ('NT', 'UBA5'):
        for drug in ('none', 'HT10', 'eme360'):
            m = MEMBER.format(geno=geno, drug=drug)
            print(f'fitting {m} ...', flush=True)
            c, grp = pc.condition_ctx(TIER, m, CACHE)
            o, fovs, per_fov = pc.per_fov_occupancy(c, grp)
            occ[m], n_fov[m], n_sub[m] = o, len(fovs), int(per_fov.sum())
            # Keep only a lightweight context. Retaining one ctx per condition
            # would hold every per-condition cache resident -- the six
            # 2026_07_08 L arrays are 25.4 GB between them -- and nothing downstream needs
            # more than the D grid and the settings block.
            ctx = dict(diff_coefs=c['diff_coefs'], settings=c['settings'])
            del c

    print('\nFOV counts (draft S7D legend says "n=20 FOV per condition"):')
    for m in occ:
        flag = '' if n_fov[m] == 20 else '   <-- not 20'
        print(f'  {m:<18} n_fov={n_fov[m]:>3}  n_sub={n_sub[m]:,}{flag}')

    # --- S7D: drug-vs-untreated within genotype -------------------------
    drug_dfs = {}
    for geno in ('NT', 'UBA5'):
        series = []
        for drug, marker in (('HT10', 'o'), ('eme360', 's')):
            a = MEMBER.format(geno=geno, drug='none')
            b = MEMBER.format(geno=geno, drug=drug)
            df = contrast(ctx, occ, a, b, f'{drug}_vs_none')
            drug_dfs[(geno, drug)] = df
            series.append(dict(df=df, color=GENO_COLOR[geno], marker=marker,
                               label=DRUG_LABEL[drug]))
        # the two drugs are drawn with different markers in one colour so the
        # genotype stays readable while the drug direction is comparable
        series[1]['color'] = ps.CONTROL_LIGHT if geno == 'NT' else ps.PERTURBED_LIGHT
        pc.plot_bin_effect(
            series, f'{OUT}/figureS7/D/figS7D_drug_shift_{GENO_LABEL[geno]}.pdf',
            ylabel=ps.YLABEL_DELTA + '\n(drug - untreated)',
            title=f'drug response, {GENO_LABEL[geno]}')

    # cross-genotype at matched treatment (legend's second allowed family)
    geno_dfs = {}
    for drug in ('none', 'HT10', 'eme360'):
        a = MEMBER.format(geno='NT', drug=drug)
        b = MEMBER.format(geno='UBA5', drug=drug)
        geno_dfs[drug] = contrast(ctx, occ, a, b, f'UBA5_vs_NT_at_{drug}')
    pc.plot_bin_effect(
        [dict(df=geno_dfs[d], color=c, marker=m,
              label=f'{DRUG_LABEL[d]}')
         for d, c, m in (('none', ps.CONTROL, 'o'),
                         ('HT10', ps.PERTURBED, '^'),
                         ('eme360', ps.PERTURBED_LIGHT, 's'))],
        f'{OUT}/figureS7/D/figS7D_genotype_effect_by_drug.pdf',
        ylabel=ps.YLABEL_DELTA + '\n(sgUBA5 - sgNT)',
        title='genotype effect at matched treatment')

    all_stats = pd.concat(list(drug_dfs.values()) + list(geno_dfs.values()))
    all_stats.to_csv(pc.tag(f'{OUT}/figureS7/D/figS7D_stats.csv'), index=False)

    # --- S7E overlays: the distributions the table summarises -----------
    # One panel per genotype, the three drug levels overlaid. Plotted from the
    # per-FOV mean soft occupancy, i.e. the SAME quantity the S7D bin panels
    # integrate, so a bin difference in S7D is the area between these curves
    # over that bin rather than a separate estimand.
    for geno in ('NT', 'UBA5'):
        curves = []
        for drug in ('none', 'HT10', 'eme360'):
            m = MEMBER.format(geno=geno, drug=drug)
            curves.append(dict(pmf=np.nanmean(occ[m], axis=0),
                               color=DRUG_COLOR[drug],
                               label=f'{DRUG_LABEL[drug]} (n={n_sub[m]:,})',
                               fill=True))
        pc.plot_overlay(
            curves, ctx['diff_coefs'], SCHEME,
            f'{OUT}/figureS7/E/figS7E_overlay_{GENO_LABEL[geno]}.pdf',
            xlabel=r'Estimated Diffusion Coefficient 60S (µm$^2$/sec)')

    # --- S7E: translation-state assignment table ------------------------
    rows = []
    for geno in ('NT', 'UBA5'):
        h, e = drug_dfs[(geno, 'HT10')], drug_dfs[(geno, 'eme360')]
        for i, label in enumerate(h['label']):
            call, code = ASSIGNMENT[i]
            hd, ed = h['observed_diff'].iloc[i], e['observed_diff'].iloc[i]
            pred, verdict = check_assignment(i, hd, ed)
            rows.append(dict(
                genotype=GENO_LABEL[geno], bin=ps.ROMAN[i], state=label,
                d_lo=h['d_lo'].iloc[i], d_hi=h['d_hi'].iloc[i],
                assignment=call, translating=code,
                harringtonine_delta=hd, harringtonine_p_bh=h['p_bh'].iloc[i],
                emetine_delta=ed, emetine_p_bh=e['p_bh'].iloc[i],
                predicted=pred, drug_consistent=verdict))
    tab = pd.DataFrame(rows)
    tab.to_csv(pc.tag(f'{OUT}/figureS7/E/figS7E_translation_state_assignment.csv'),
               index=False)
    print('\nS7E translation-state assignment (stated; drugs are the test):')
    print(tab[['genotype', 'bin', 'assignment', 'harringtonine_delta',
               'emetine_delta', 'predicted', 'drug_consistent']]
          .to_string(index=False))
    tested = tab[tab.drug_consistent != '']
    n_ok = int((tested.drug_consistent == 'consistent').sum())
    print(f'  {n_ok}/{len(tested)} testable bin-genotype pairs agree with '
          f'the stated assignment')
    return 0


if __name__ == '__main__':
    sys.exit(main())
