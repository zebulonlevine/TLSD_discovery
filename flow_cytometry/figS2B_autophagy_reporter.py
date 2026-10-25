"""Fig S2B: GFP-LC3-RFP-LC3ΔG autophagic-flux reporter in the RPS3 (S3) and RPL29 (L29) screening lines.

GFP-LC3 is degraded in the lysosome during autophagy while RFP-LC3ΔG is not, so blocking autophagy
raises the GFP/RFP ratio. Each well is one value: the median per-cell GFP/RFP ratio of single cells,
exported from FlowJo (column `Ratio`). 3 wells per clone x medium x perturbation, one experiment
(2024-06-25). Media: 10% FBS (growing), overnight EBSS (amino-acid starvation), 1% FBS (quiescent).

  1. Each knockdown (RB1CC1 or ATG7) is compared with the non-targeting control arm in the same
     clone and medium (vehicle wells).
  2. Effect = 100 x (mean_test / mean_control - 1), delta-method SE of the ratio of two independent
     means; effective df by Welch-Satterthwaite.
  3. Clone-pooled effect = unweighted mean of the S3 and L29 effects, SE in quadrature.
  4. Brackets: (EBSS - 10% FBS) and (1% FBS - 10% FBS) contrasts of the pooled effect, in percentage
     points; t = contrast / SE at the Welch-Satterthwaite df, two-sided, not adjusted for multiple
     comparisons. * p < 0.05, ** p < 0.01, *** p < 0.001.
The BafA1 and EBSS + STLC wells in the export are not analysed.

Outputs (output root, flow_cytometry/):
  figS2B_autophagy_flux.pdf
  figS2B_effects.csv     per clone and pooled: arm means, absolute and percent effects, SE, df, p
  figS2B_contrasts.csv   starvation - 10% FBS contrasts of the percent effect, pooled and per clone

    python flow_cytometry/figS2B_autophagy_reporter.py
"""
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from scipy import stats as sp_stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import paths  # noqa: E402

INPUT = ('flow_cytometry', 'S2B', 'autophagyData.csv')     # FlowJo per-well export
OUT_DIR = 'flow_cytometry'

BASES = ['10%', 'EBSS', '1%']
CLONES = ['S3', 'L29']
FACTORS = ['RB1CC1', 'ATG7']
FACTOR_LABELS = {'RB1CC1': 'RB1CC1 KD', 'ATG7': 'ATG7 KD'}
# export Treatment (vehicle wells) -> medium; BafA1 and EBSS + STLC wells are not analysed
TREATMENTS = {'10p_noBaf': '10%', 'EBSS ON': 'EBSS', '1p': '1%'}

COLORS = ['#534AB7', '#BA7517']
CTRL_COL = '#888880'
S3_COL, L29_COL = '#3266ad', '#c94040'


def out(name):
    return paths.outp(OUT_DIR, name, mkdir=True)


# ---------------------------------------------------------------------------
# estimates with propagated error
# ---------------------------------------------------------------------------

class Est:
    """A scalar estimate carrying a standard error and effective df."""
    __slots__ = ("val", "se", "df")

    def __init__(self, val, se, df):
        self.val, self.se, self.df = float(val), float(se), float(df)

    def __repr__(self):
        return f"Est({self.val:+.4g} +/- {self.se:.4g}, df={self.df:.2f})"


def mean_est(vals):
    """Group mean with its SEM; df = n - 1."""
    a = np.asarray(vals, dtype=float)
    n = a.size
    return Est(a.mean(), a.std(ddof=1) / math.sqrt(n), n - 1)


def lincomb(terms):
    """Linear combination sum(c * Est) with quadrature error propagation
    and Welch-Satterthwaite effective df.  Components must be independent
    (they are here: every term comes from a disjoint set of wells)."""
    val   = sum(c * e.val for c, e in terms)
    parts = [(c * e.se) ** 2 for c, e in terms]
    var   = sum(parts)
    denom = sum(p ** 2 / e.df for p, (_, e) in zip(parts, terms) if p > 0)
    df    = var ** 2 / denom if denom > 0 else float("inf")
    return Est(val, math.sqrt(var), df)


def pct_change(test, ctrl):
    """100 * (mean_test / mean_ctrl - 1), delta-method SE for a ratio of
    two independent means, Welch-Satterthwaite df on the relative terms."""
    R  = test.val / ctrl.val
    va = (test.se / test.val) ** 2
    vb = (ctrl.se / ctrl.val) ** 2
    se = abs(R) * math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / test.df + vb ** 2 / ctrl.df)
    return Est(100.0 * (R - 1.0), 100.0 * se, df)


def pval(e):
    """Two-sided Welch's t: estimate / propagated SE at its effective df."""
    if e.se == 0:
        return float("nan")
    return float(2 * sp_stats.t.sf(abs(e.val / e.se), e.df))


def sig(p):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------

def load_wells():
    """(clone, medium, knockdown) -> per-well median GFP/RFP ratios."""
    wells = defaultdict(list)
    with open(paths.inp(*INPUT), "r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            d = {k.strip(): v.strip() for k, v in row.items()}
            base = TREATMENTS.get(d["Treatment"])
            if base:
                wells[(d["CellLine"], base, d["KD"])].append(float(d["Ratio"]))
    return wells


def arms(wells, factor, clone, base):
    """(knockdown wells, non-targeting wells) for a factor in one clone and medium."""
    return wells[(clone, base, factor)], wells[(clone, base, "None")]


def effects(wells):
    """Per-clone effects, clone-pooled effects, and starvation - 10% FBS contrasts."""
    eff = {}
    for f in FACTORS:
        for cl in CLONES:
            for b in BASES:
                t_vals, c_vals = arms(wells, f, cl, b)
                t, c = mean_est(t_vals), mean_est(c_vals)
                eff[(f, cl, b)] = {
                    "test": t, "ctrl": c,
                    "abs": lincomb([(1.0, t), (-1.0, c)]), "pct": pct_change(t, c),
                    # Welch's two-sample t on the wells (same as the propagated test on 'abs')
                    "p_welch": float(sp_stats.ttest_ind(t_vals, c_vals, equal_var=False)[1]),
                    "n_test": len(t_vals), "n_ctrl": len(c_vals),
                }
    pooled, pooled_abs = {}, {}
    for f in FACTORS:
        for b in BASES:
            pooled[(f, b)] = lincomb([(0.5, eff[(f, cl, b)]["pct"]) for cl in CLONES])
            pooled_abs[(f, b)] = lincomb([(0.5, eff[(f, cl, b)]["abs"]) for cl in CLONES])
    contrast, contrast_clone = {}, {}
    for f in FACTORS:
        for b in BASES[1:]:
            contrast[(f, b)] = lincomb([(1.0, pooled[(f, b)]), (-1.0, pooled[(f, "10%")])])
            for cl in CLONES:
                contrast_clone[(f, cl, b)] = lincomb(
                    [(1.0, eff[(f, cl, b)]["pct"]), (-1.0, eff[(f, cl, "10%")]["pct"])])
    return eff, pooled, pooled_abs, contrast, contrast_clone


def effects_table(eff, pooled, pooled_abs):
    rows = []
    for f in FACTORS:
        for b in BASES:
            for cl in CLONES:
                e = eff[(f, cl, b)]
                rows.append(dict(readout=FACTOR_LABELS[f], clone=cl, medium=b,
                                 n_test=e["n_test"], n_ctrl=e["n_ctrl"],
                                 ctrl_mean=e["ctrl"].val, ctrl_se=e["ctrl"].se,
                                 test_mean=e["test"].val, test_se=e["test"].se,
                                 abs_diff=e["abs"].val, abs_se=e["abs"].se, abs_df=e["abs"].df,
                                 pct_change=e["pct"].val, pct_se=e["pct"].se, pct_df=e["pct"].df,
                                 p_welch=e["p_welch"]))
            p, pa = pooled[(f, b)], pooled_abs[(f, b)]
            rows.append(dict(readout=FACTOR_LABELS[f], clone='pooled', medium=b,
                             abs_diff=pa.val, abs_se=pa.se, abs_df=pa.df,
                             pct_change=p.val, pct_se=p.se, pct_df=p.df, p_pct_vs_0=pval(p)))
    return pd.DataFrame(rows)


def contrasts_table(contrast, contrast_clone):
    rows = []
    for f in FACTORS:
        for b in BASES[1:]:
            for cl, c in [('pooled', contrast[(f, b)])] + [(cl, contrast_clone[(f, cl, b)]) for cl in CLONES]:
                rows.append(dict(readout=FACTOR_LABELS[f], clone=cl, contrast=f'{b} - 10%',
                                 diff_pp=c.val, se=c.se, df=c.df, t=c.val / c.se, p=pval(c),
                                 stars=sig(pval(c)), plotted=cl == 'pooled'))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------

def plot(eff, pooled, contrast, path):
    """Percent effect per medium, clones pooled (bars) with each clone's estimate (open markers)."""
    fig = plt.figure(figsize=(6.0, 3.6))
    gs = fig.add_gridspec(1, len(FACTORS), wspace=0.08,
                          left=0.156, right=0.958, bottom=0.17, top=0.85)

    axes_top = [fig.add_subplot(gs[0, i]) for i in range(len(FACTORS))]
    for ax in axes_top[1:]:
        ax.sharey(axes_top[0])

    for ax_i, (ax, f, col) in enumerate(zip(axes_top, FACTORS, COLORS)):
        x = np.arange(len(BASES))
        means = [pooled[(f, b)].val for b in BASES]
        sems  = [pooled[(f, b)].se  for b in BASES]
        colors = [CTRL_COL if i == 0 else col for i in range(len(BASES))]

        ax.bar(x, means, width=0.55, color=colors, alpha=0.75,
               edgecolor="none", zorder=2)
        ax.errorbar(x, means, yerr=sems, fmt="none", ecolor="black",
                    elinewidth=1.2, capsize=4, capthick=1.2, zorder=3)

        # per-clone estimates, each with its own propagated error
        for cl, dx, mk, ccol in [("S3", -0.17, "o", S3_COL), ("L29", 0.17, "D", L29_COL)]:
            v  = [eff[(f, cl, b)]["pct"].val for b in BASES]
            sv = [eff[(f, cl, b)]["pct"].se  for b in BASES]
            ax.errorbar(x + dx, v, yerr=sv, fmt=mk, markersize=4,
                        markerfacecolor="white", markeredgecolor=ccol,
                        markeredgewidth=1.0, ecolor=ccol, elinewidth=0.8,
                        capsize=2, capthick=0.8, linestyle="none", zorder=4)

        # brackets: EBSS vs 10%, 1% vs 10%  (Welch's t on the pp contrast)
        tops = [max(pooled[(f, b)].val + pooled[(f, b)].se,
                    max(eff[(f, cl, b)]["pct"].val + eff[(f, cl, b)]["pct"].se
                        for cl in CLONES)) for b in BASES]
        span = max(tops) - min(0.0, min(pooled[(f, b)].val - pooled[(f, b)].se for b in BASES))
        step = 0.11 * max(span, 1.0)
        tier = 0
        for i, b in enumerate(BASES[1:], 1):
            c = contrast[(f, b)]
            p = pval(c)
            y_bar = max(tops[: i + 1]) + 0.35 * step + tier * step
            ax.plot([0, 0, i, i],
                    [y_bar - 0.18 * step, y_bar, y_bar, y_bar - 0.18 * step],
                    color="0.4", linewidth=0.7, zorder=5)
            ax.text(i / 2, y_bar + 0.05 * step, sig(p), ha="center", va="bottom",
                    fontsize=8.5, fontweight="bold" if p < 0.05 else "normal",
                    color="#A32D2D" if p < 0.05 else "0.5", zorder=5)
            tier += 1

        ax.set_title(FACTOR_LABELS[f], fontsize=10, fontweight="medium", pad=8)
        ax.set_xticks(x)
        ax.set_xticklabels(["10%\nFBS", "EBSS", "1%\nFBS"], fontsize=9)
        ax.axhline(0, color="0.7", linewidth=0.6, zorder=1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.6)
        ax.spines["left"].set_linewidth(0.6)
        ax.tick_params(axis="both", length=3, width=0.6)
        if ax_i == 0:
            ax.set_ylabel("Change in GFP/RFP ratio\nvs matched control (%)", fontsize=9)
        else:
            ax.tick_params(labelleft=False)

    ymin = min(ax.get_ylim()[0] for ax in axes_top)
    ymax = max(ax.get_ylim()[1] for ax in axes_top)
    pad = 0.04 * (ymax - ymin)
    for ax in axes_top:
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:+.0f}"))

    handles = [
        Line2D([0], [0], marker="o", color="w", markeredgecolor=S3_COL,
               markerfacecolor="w", markersize=5, linewidth=0, label="S3"),
        Line2D([0], [0], marker="D", color="w", markeredgecolor=L29_COL,
               markerfacecolor="w", markersize=4, linewidth=0, label="L29"),
    ]
    axes_top[0].legend(handles=handles, loc="upper left", frameon=False,
                       fontsize=8, handletextpad=0.3, borderpad=0.2)

    fig.text(0.557, 0.965,
             "Knockdown effect by medium (clones pooled; brackets = t vs 10% FBS)",
             ha="center", fontsize=10, fontweight="medium")

    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def main():
    eff, pooled, pooled_abs, contrast, contrast_clone = effects(load_wells())
    effects_table(eff, pooled, pooled_abs).to_csv(out('figS2B_effects.csv'), index=False)
    contrasts_table(contrast, contrast_clone).to_csv(out('figS2B_contrasts.csv'), index=False)
    plot(eff, pooled, contrast, out('figS2B_autophagy_flux.pdf'))
    for f in FACTORS:
        print(f'{FACTOR_LABELS[f]:>10s}: ' + '  '.join(
            f'{b} {pooled[(f, b)].val:+.1f}% {sig(pval(contrast[(f, b)])) if b != "10%" else ""}'.rstrip()
            for b in BASES))
    print('wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
