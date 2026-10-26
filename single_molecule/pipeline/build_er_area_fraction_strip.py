"""
build_er_area_fraction_strip.py
===============================

Per-FOV ER AREA fraction (% of image pixels inside the ER mask) for Degron vs
noDegron, as a strip plot with a Welch t-test.

Why this matters beyond description: the zone contrasts classify trajectories by
which mask they sit in, so a systematic difference in mask COVERAGE between the
two cell lines is a confound. It changes how much of each field is available to
be called on-ER, and it directly changes the per-FOV binomial-test power
(CLAUDE.md: n_min_binom is computed per FOV from area, and high-coverage FOVs make
most trajectories untestable). So this is a control on the zone results, not just
a morphology readout.

`% not-ER` is NOT plotted separately: on_structure / near_ER / off_ER partition
the image exactly (verified at load -- the three coverages sum to 1.0 per FOV), so
not-ER is 100 - on_ER_broad and would be a mirror-image panel carrying no extra
information. Stated here so its absence isn't mistaken for an omission.

PRE-SPECIFIED, to avoid picking the mask definition after seeing the p-values
(CLAUDE.md sec 9 forbids that):
  PRIMARY   on_ER_broad  (= on_structure | near_ER) -- the definition the headline
                          zone contrast uses
  SECONDARY on_structure  -- the strict definition; nested inside the primary, so
                          the two tests are strongly correlated and NOT
                          independent comparisons to correct across.
Welch (unequal-variance) t-test is the headline; Mann-Whitney U is reported
alongside because n = 10 / 14 is small enough that normality is not testable.

FOV set: untreated only (HHT=000 + HHT=005_t=000), the same 24 fields as
build_er_boundary_profile_qc.py and the Fig 5F contrasts. Asserted at runtime.

============ WHAT THIS PANEL DOES AND DOES NOT SHOW -- READ THIS ============
run_filewise_er._three_class_mask applies pct_floor=50: a pixel must exceed the
FIELD'S OWN MEDIAN to be eligible for any ER class (`above_floor = image >
percentile(image, 50)`), and both mask_broad and on_structure are intersected with
it. Since on_ER_broad = on_structure | mask_broad, BOTH definitions plotted here
are HARD-CAPPED AT 50% of the field by construction. On top of that the local
criterion is `img > local_mean * (1 + 1.5/100)`, i.e. "brighter than your own
neighbourhood", which on a smoothly varying background selects a fairly stable
pixel fraction whatever the absolute ER content.

Consequences:
  VALID    as a CONFOUND CONTROL for the zone contrasts -- mask coverage is
           matched between the lines, so on-vs-off-ER comparisons are not driven
           by one line simply having more of its field labelled on-ER, and the
           per-FOV binomial power (n_min_binom, which is computed from area) is
           comparable across lines.
  NOT VALID as a measurement of ER ABUNDANCE or morphology. An adaptive local
           threshold with a median floor would return ~30% coverage in both lines
           even if their true ER content differed substantially, so this test has
           little power to detect a real difference. A null here is close to
           circular for that question -- the same near-circularity
           build_er_boundary_profile_qc.py exists to work around for brightness.
Do not report this as "the two lines have equivalent ER".
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
SEG_DIR = os.path.join(tlsd_paths.spt('all_data', '2025_09_15_zeb_ER'), 'tracking_output_q=15p0', 'er_segmentation')
OUT_DIR = tlsd_paths.outp('PaperFigures_updated', 'figureS7', 'I')

CONDITIONS = ['noDegron', 'Degron']
COLORS = {'noDegron': '#4d4d4d', 'Degron': '#e6550d'}
MARKERS = {'noDegron': 'o', 'Degron': 's'}
MASKS = [('on_ER_broad', 'PRIMARY'), ('on_structure', 'SECONDARY')]
EXPECTED_FOV = {'noDegron': 10, 'Degron': 14}
SEED = 0


def is_untreated(stem):
    return ('HHT=000' in stem) or ('HHT=005_t=000_' in stem)


def collect():
    rows = []
    for f in sorted(os.listdir(SEG_DIR)):
        if not f.endswith('_er_masks.npz'):
            continue
        stem = f[:-len('_er_masks.npz')]
        if not is_untreated(stem):
            continue
        with np.load(os.path.join(SEG_DIR, f)) as mm:
            on = mm['on_structure'].astype(bool)
            near = mm['near_ER'].astype(bool)
            off = mm['off_ER'].astype(bool)
        total = on.mean() + near.mean() + off.mean()
        if abs(total - 1.0) > 1e-6:
            raise SystemExit(f'{stem}: on_structure+near_ER+off_ER covers '
                             f'{total:.6f} of the image, not 1.0 -- the three '
                             f'labels are not the partition this script assumes, '
                             f'so "not-ER = 100 - on_ER_broad" would be wrong')
        rows.append(dict(
            fov=stem,
            degron='Degron' if 'degron=001' in stem else 'noDegron',
            on_structure=100 * on.mean(),
            near_ER=100 * near.mean(),
            on_ER_broad=100 * (on | near).mean(),
            off_ER=100 * off.mean()))
    return pd.DataFrame(rows)


def annotate_sig(ax, x1, x2, y, p, h=None):
    stars = ('n.s.' if p >= 0.05 else '*' if p >= 0.01
             else '**' if p >= 0.001 else '***')
    h = h if h is not None else 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color='k', lw=1.1)
    ax.text((x1 + x2) / 2, y + h, f'{stars}\np={p:.3g}', ha='center',
            va='bottom', fontsize=8.5)


def main():
    df = collect()
    counts = df['degron'].value_counts().to_dict()
    print(f'FOVs: {counts}  (untreated only)')
    if counts != EXPECTED_FOV:
        raise SystemExit(f'FOV counts {counts} != expected {EXPECTED_FOV}')
    os.makedirs(OUT_DIR, exist_ok=True)
    df.to_csv(os.path.join(OUT_DIR, 'er_area_fraction_per_fov.csv'), index=False)

    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(1, len(MASKS), figsize=(4.4 * len(MASKS), 4.8))
    stat_rows = []

    for ax, (mask, role) in zip(np.atleast_1d(axes), MASKS):
        vals = {c: df.loc[df['degron'] == c, mask].to_numpy() for c in CONDITIONS}
        t, p_t = stats.ttest_ind(vals['noDegron'], vals['Degron'],
                                 equal_var=False)          # Welch
        u, p_u = stats.mannwhitneyu(vals['noDegron'], vals['Degron'],
                                    alternative='two-sided')
        # Hedges-corrected Cohen's d, pooled SD
        n1, n2 = len(vals['noDegron']), len(vals['Degron'])
        sp = np.sqrt(((n1 - 1) * vals['noDegron'].var(ddof=1)
                      + (n2 - 1) * vals['Degron'].var(ddof=1)) / (n1 + n2 - 2))
        d = (vals['Degron'].mean() - vals['noDegron'].mean()) / sp if sp else np.nan
        g = d * (1 - 3 / (4 * (n1 + n2) - 9))

        for i, c in enumerate(CONDITIONS):
            v = vals[c]
            x = i + rng.uniform(-0.10, 0.10, len(v))
            ax.scatter(x, v, s=46, facecolor=COLORS[c], marker=MARKERS[c],
                       edgecolor='k', linewidth=0.6, alpha=0.9, zorder=3)
            m, sem = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
            ax.plot([i - 0.24, i + 0.24], [m, m], color='k', lw=2.0, zorder=4)
            ax.errorbar([i], [m], yerr=[sem], fmt='none', ecolor='k',
                        elinewidth=1.3, capsize=5, zorder=4)
            print(f'  {mask:13s} {c:9s} n={len(v):2d}  mean={m:5.2f}%  '
                  f'SEM={sem:4.2f}  range {v.min():5.2f}-{v.max():5.2f}')

        top = max(v.max() for v in vals.values())
        ax.set_ylim(0, top * 1.34)
        annotate_sig(ax, 0, 1, top * 1.06, p_t)
        ax.set_xticks(range(len(CONDITIONS)))
        ax.set_xticklabels(CONDITIONS)
        ax.set_xlim(-0.6, len(CONDITIONS) - 0.4)
        ax.set_title(f'{mask}  [{role}]', fontsize=10)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        print(f'  {mask}: Welch t={t:+.3f} p={p_t:.4g} | MWU U={u:.0f} '
              f'p={p_u:.4g} | Hedges g={g:+.2f}')
        stat_rows.append(dict(mask=mask, role=role,
                              n_noDegron=n1, n_Degron=n2,
                              mean_noDegron=vals['noDegron'].mean(),
                              mean_Degron=vals['Degron'].mean(),
                              sem_noDegron=vals['noDegron'].std(ddof=1)/np.sqrt(n1),
                              sem_Degron=vals['Degron'].std(ddof=1)/np.sqrt(n2),
                              welch_t=t, welch_p=p_t, mwu_U=u, mwu_p=p_u,
                              hedges_g=g))

    np.atleast_1d(axes)[0].set_ylabel('ER mask area (% of field)')
    fig.suptitle('Per-FOV ER mask area fraction by degron '
                 '(untreated FOVs: HHT=000 + HHT=005 t=0)\n'
                 'dot = one FOV; bar = mean ± SEM; Welch two-sample t-test',
                 y=1.02, fontsize=10.5)
    fig.tight_layout()
    pdf = os.path.join(OUT_DIR, 'er_area_fraction_strip.pdf')
    fig.savefig(pdf, bbox_inches='tight')
    plt.close(fig)
    sp = os.path.join(OUT_DIR, 'er_area_fraction_stats.csv')
    pd.DataFrame(stat_rows).to_csv(sp, index=False)
    print(f'\nwrote {pdf}\nwrote {sp}')


if __name__ == '__main__':
    main()
