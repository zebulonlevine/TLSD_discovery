#!/usr/bin/env python3
"""
Regenerate Figure S3B -- genes selected into the sublibrary, shown on the
genome-wide new/old ratio axes.

Equivalent to Fig2.ipynb cells 0, 1, 3, 4. Lifted into a script so it can be
re-run after a genome-wide rescore without executing the other 50 cells of
Fig2.ipynb, which are sublibrary analyses that are deliberately NOT being
re-run.

The coordinates are genome-wide (L29/S3 _TopvBottom); the sublibrary results
supply only the list of genes that were selected into the sublibrary, so this
panel does move when the genome-wide scores change.

Adds the legend-less PDF that never existed -- previously only a PNG was
written, which is un-editable in Illustrator.
"""

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

GW = paths.inp('GenomeWide_results')
SUB = paths.inp('Sublibrary_results')
OUTDIR = paths.outp('Sublibrary_results', 'plots', 'GeneScatter')
os.makedirs(OUTDIR, exist_ok=True)

SUBLIB_COMPARISONS = [
    'S3_unsorted_v_early_wt', 'L29_unsorted_v_early_wt',
    'S3_old_top_v_bottom_wt', 'L29_old_top_v_bottom_wt',
    'S3_new_top_v_bottom_wt', 'L29_new_top_v_bottom_wt',
    'S3_ratioON_top_v_bottom_wt', 'L29_ratioON_top_v_bottom_wt',
]


def load(base, name, label):
    return screen_viz.load_comparison(
        results_file=os.path.join(base, 'results', name, name + '_all_genes.csv'),
        hits_file=os.path.join(base, 'results', name, name + '_hits_fdr0.05.csv'),
        name=label)


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    L29 = load(GW, 'L29_TopvBottom', 'L29 Ratio Genomewide')
    S3 = load(GW, 'S3_TopvBottom', 'S3 Ratio Genomewide')

    # Union of genes scored in any sublibrary comparison = the sublibrary
    # content. 'NTC' here is the pseudogene label, excluded by exact match.
    genes = []
    seen = set()
    for name in SUBLIB_COMPARISONS:
        c = load(SUB, name, name)
        for g in c.results['gene'].astype(str):
            if g != 'NTC' and g not in seen:
                seen.add(g)
                genes.append(g)
    cat_sublib = pd.DataFrame({'gene': genes, 'category': 'sublibrary'})
    print(f"  sublibrary genes: {len(cat_sublib)}")

    present = set(L29.results['gene'].astype(str)) | set(S3.results['gene'].astype(str))
    print(f"  of those, present in the genome-wide results: {len(seen & present)}")

    common = dict(
        categories=cat_sublib,
        title='Sublibrary Selection from Genome-Wide',
        figsize=(6, 4),
        show_correlation=False,
        ref_lines=None,
        show_diagonal=False,
        show_axis_lines=False,
        colors={'sublibrary': '#0072B2'},
    )
    stem = os.path.join(OUTDIR, 'Genomewide_includedGenes_sublib')
    # PDF added here; previously only the PNG and the legend PDF existed, so the
    # clean panel could not be edited as vector art.
    screen_viz.plot_comparison_scatter(L29, S3, output_file=stem + '.pdf',
                                       show_legend=False, **common)
    screen_viz.plot_comparison_scatter(L29, S3, output_file=stem + '.png',
                                       show_legend=False, **common)
    screen_viz.plot_comparison_scatter(L29, S3, output_file=stem + '_legend.pdf',
                                       show_legend=True, legend_outside=True, **common)
    print('\nDone.')


if __name__ == '__main__':
    main()
