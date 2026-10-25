#!/usr/bin/env python3
"""Figure S4A (40S vs 60S degradation, and 40S vs 60S synthesis) coloured with
the bioBySubunit scheme used by Fig 2C/D (Fig2.ipynb cells 47/48).

Fig2 cells 23/29 plot these two comparisons with the flat 3-category
GenomeWideCats_threeOptions palette. This script keeps their geometry (diagonal,
axis labels, point size, per-axis top-5 labelling) but swaps in the
bioBySubunit category table and three-shades-of-green palette, so S4A is
consistent with 2C/D.

Category construction is copied verbatim from cell 47 so the two panels can't
diverge: biogenesis split by subunit from GenomeWideCats_synth_bySubunit.csv,
RNA Pol I forced into 'Both subunits', non-biogenesis carried over from
GenomeWideCats_threeOptions.csv, and RNASEK/ATG14 hand-assigned to
Lysosome/Autophagy.

Writes NEW filenames; the cell 23/29 outputs are untouched.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, REPO)
LIB = os.path.join(REPO, 'screen_figures', 'lib')
sys.path.insert(0, LIB)
from common import paths  # noqa: E402
sys.modules.pop('screen_viz', None)
import screen_viz  # noqa: E402
assert os.path.dirname(os.path.abspath(screen_viz.__file__)) == LIB

CATS = paths.repo('screen_figures', 'categories')

RNA_POL_I = ['POLR1A', 'POLR1B', 'POLR1C', 'POLR1D',
             'POLR1E', 'POLR1F', 'POLR1G', 'POLR1H']
LYSO_MANUAL = ['RNASEK', 'ATG14']

GREEN_COLORS = {
    'Both subunits':           '#00441b',   # dark green (includes RNA Pol I)
    '60S Biogenesis/Proteins': '#238b45',   # medium green
    '40S Biogenesis/Proteins': '#74c476',   # light green
    'Lysosome/Autophagy':      '#c51b7d',
    'UFMylation':              '#E69F00',
}
BIO_ORDER = ['Both subunits', '60S Biogenesis/Proteins', '40S Biogenesis/Proteins',
             'Lysosome/Autophagy', 'UFMylation']


def bysubunit_categories():
    """Exactly the table Fig2 cell 47 builds."""
    bio = screen_viz.load_categories(
        os.path.join(CATS, 'GenomeWideCats_synth_bySubunit.csv'))[['gene', 'category']].copy()
    bio['category'] = bio['category'].replace({'Both Biogenesis': 'Both subunits'})
    bio.loc[bio['gene'].isin(RNA_POL_I), 'category'] = 'Both subunits'

    base = screen_viz.load_categories(
        os.path.join(CATS, 'GenomeWideCats_threeOptions.csv'))[['gene', 'category']]
    non_bio = base.loc[base['category'] != 'Biogenesis']

    cats = pd.concat([bio, non_bio], ignore_index=True)
    cats = cats[~cats['gene'].isin(LYSO_MANUAL)]
    cats = pd.concat([cats, pd.DataFrame({'gene': LYSO_MANUAL,
                                          'category': ['Lysosome/Autophagy'] * len(LYSO_MANUAL)})],
                     ignore_index=True)
    return cats


def top_each_direction(cx, cy, n=5):
    """Label the n most extreme genes in EACH direction on EACH axis.

    Copied from Fig2 cell 47. label_top_n_x/label_top_n_y (what cells 23/29 use)
    rank by |epsilon|, which picks 5 genes total per axis and in practice lands
    them all at one end; this returns both tails of both axes, and combined with
    adjust_labels=True it stops the pile-ups.
    """
    dx = screen_viz.deduplicate_by_max_abs_epsilon(cx.results, 'gene_id', 'epsilon')
    dy = screen_viz.deduplicate_by_max_abs_epsilon(cy.results, 'gene_id', 'epsilon')
    m = dx[['gene_id', 'gene', 'epsilon']].merge(
        dy[['gene_id', 'epsilon']], on='gene_id', suffixes=('_x', '_y'))
    genes = set()
    for col in ['epsilon_x', 'epsilon_y']:
        genes |= set(m.nlargest(n, col)['gene'])
        genes |= set(m.nsmallest(n, col)['gene'])
    return sorted(genes)


def load(score_dir, label):
    return screen_viz.load_comparison(
        os.path.join(score_dir, 'gene_scores_all_genes.csv'),
        os.path.join(score_dir, 'gene_scores_hits_fdr0.05.csv'), label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--published', action='store_true',
                    help='use the July line-score gene tables (Sublibrary_results/scores/)')
    ap.add_argument('--scores-root', default=None,
                    help='any gene-score tree laid out like fixed_NotI/scores/ (default)')
    ap.add_argument('--suffix', default=None)
    ap.add_argument('--outdir', default=paths.outp('fixed_NotI', 'plots'))
    args = ap.parse_args()

    root = args.scores_root or (paths.inp('Sublibrary_results', 'scores') if args.published
                                else paths.inp('fixed_NotI', 'scores'))
    sfx = args.suffix if args.suffix is not None else (
        '_published' if args.published else '_fixedNotI')
    os.makedirs(args.outdir, exist_ok=True)
    cats = bysubunit_categories()
    print(f"  category table: {len(cats)} genes, "
          f"{cats.category.nunique()} categories -> {sorted(cats.category.unique())}")

    jobs = [
        ('Degradation', 'L29_Degradation', 'S3_Degradation',
         '60S Degradation Score', '40S Degradation Score',
         f'compare_40S60SDegscores_bioBySubunit{sfx}.pdf'),
        ('Synthesis', 'L29_Synthesis', 'S3_Synthesis',
         '60S Synthesis Score', '40S Synthesis Score',
         f'compare_40S60SSynthscores_bioBySubunit{sfx}.pdf'),
    ]
    for what, xdir, ydir, xlab, ylab, out in jobs:
        c1 = load(os.path.join(root, xdir), xlab)
        c2 = load(os.path.join(root, ydir), ylab)
        # Axes geometry kept from cells 23/29 (diagonal + axis labels, which are
        # meaningful for a cross-subunit comparison); colouring AND labelling
        # taken from the bioBySubunit style of cell 47.
        common = dict(
            categories=cats,
            title=None,
            label_genes=top_each_direction(c1, c2, 5),
            label_top_n_x=0, label_top_n_y=0, adjust_labels=True,
            show_diagonal=True, show_correlation=False,
            figsize=(5.5, 3.2), legend_outside=True,
            size_by_hits=False, point_size=60, fontsize=8,
            colors=GREEN_COLORS, category_order=BIO_ORDER,
            dpi=600,
        )
        screen_viz.plot_comparison_scatter(c1, c2, os.path.join(args.outdir, out), **common)
        screen_viz.plot_comparison_scatter(
            c1, c2, os.path.join(args.outdir, out.replace('.pdf', '_noLegend.pdf')),
            show_legend=False, **common)
        n = cats['gene'].isin(set(c1.results['gene'].astype(str))).sum()
        print(f"  {what}: {n} categorised genes present in the sublibrary")


if __name__ == '__main__':
    main()
