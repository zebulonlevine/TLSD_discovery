#!/usr/bin/env python3
"""
Regenerate the screen-derived panels of Figure 1 / S1 / S2A / S6B from the current
genome-wide scores.

Equivalent to running Fig1.ipynb cells 0, 2, 4, 5, 6, 7, 11, 12, 25 -- lifted
into a script because several other cells in that notebook cannot be executed
(cell 26 is a bare `y = `, cell 36 references an undefined name, cells 29-35
make live network calls, and cell 40 rewrites decoupler inside site-packages).
Inputs are read from the data root and outputs written under the output root (common/paths).

Adds the two panels that had no code: S1D (correlation heatmap) and 1F (GO
enrichment).

    python fig1_figS1_figS2A_figS6B.py [--panels 1E,S1D,S1E,S1F,S2A,S6B,1F]
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')

# Stale-import guard. There are ~10 copies of screen_viz.py in this tree and
# only claude_tools/ is live; importlib.reload cannot switch between them, so a
# stale import fails later with a confusing TypeError about an unexpected
# keyword argument, having written no file.
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, REPO)
LIB = os.path.join(REPO, 'screen_figures', 'lib')
sys.path.insert(0, LIB)
from common import paths  # noqa: E402
sys.modules.pop('screen_viz', None)
import screen_viz  # noqa: E402
assert os.path.dirname(os.path.abspath(screen_viz.__file__)) == LIB

CATS = paths.repo('screen_figures', 'categories')
RES = paths.inp('GenomeWide_results', 'results')
PLOTS_FIG1 = paths.outp('GenomeWide_results', 'plots', 'Fig1')
PLOTS_ESS = paths.outp('GenomeWide_results', 'plots', 'EssentialsScatter')
os.makedirs(PLOTS_FIG1, exist_ok=True)
os.makedirs(PLOTS_ESS, exist_ok=True)

L29NewOld, S3NewOld = 'L29_TopvBottom', 'S3_TopvBottom'
L29Growth, S3Growth = 'L29_UnsortedvEarly', 'S3_UnsortedvEarly'

# Fig 1E palette, reused for the 1F label colours so the two panels agree.
FIG1E_COLORS = {
    'Biogenesis': '#4e9342',
    'Lysosome/Autophagy': '#c51b7d',
    'UFMylation': '#E69F00',
}

# Pinned from the default CATEGORY_PALETTE assignment that produced the draft
# panel (sorted category order -> palette index). Pinning stops the colours
# shuffling if a category is ever added to Cats_Hippo.csv.
HIPPO_COLORS = {
    'core kinase cascade': '#0072B2',
    'negative regulator': '#E69F00',
    'positive regulator': '#009E73',
    'transcriptional module': '#CC79A7',
}
HIPPO_ORDER = ['positive regulator', 'negative regulator',
               'core kinase cascade', 'transcriptional module']

# 'RQC' was in the notebook's category_order but is absent from manyCats2.csv,
# so it silently drew nothing. Dropped rather than left to reappear.
MANYCATS_ORDER = ['Both Biogenesis', '40S Biogenesis/Proteins', '60S Biogenesis/Proteins',
                  'Autophagy', 'Lysosome', 'Proteasome', 'E3', 'UFMylation']

TERM_LABELS = {
    'Protein Ufmylation (GO:0071569)': 'Protein UFMylation',
    'Regulation of Intracellular Steroid Hormone Receptor Signaling Pathway (GO:0033143)':
        'Reg. of Steroid Hormone Signalling',
    'RNA N6-methyladenosine Methyltransferase Complex (GO:0036396)':
        'RNA m6A Methyltransferase Complex',
    'Positive Regulation of TORC1 Signaling (GO:1904263)': 'Pos. Reg. of TORC1 Signaling',
}


def comp(name, label):
    return screen_viz.load_comparison(
        results_file=os.path.join(RES, name, name + '_all_genes.csv'),
        hits_file=os.path.join(RES, name, name + '_hits_fdr0.05.csv'),
        name=label)


def cat(filename):
    return screen_viz.load_categories(os.path.join(CATS, filename))


def three_ways(fn, outdir, stem, **kw):
    """PDF / PNG / legend-PDF, the pattern every Fig1 cell uses."""
    fn(output_file=os.path.join(outdir, stem + '.pdf'), show_legend=False, **kw)
    fn(output_file=os.path.join(outdir, stem + '.png'), show_legend=False, **kw)
    fn(output_file=os.path.join(outdir, stem + '_withLegend.pdf'),
       show_legend=True, legend_outside=True, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panels', default='1E,S1D,S1E,S1F,S2A,S6B,1F')
    args = ap.parse_args()
    want = {p.strip().upper() for p in args.panels.split(',')}

    os.makedirs(PLOTS_FIG1, exist_ok=True)
    os.makedirs(PLOTS_ESS, exist_ok=True)

    from functools import partial

    need_scatter = want & {'1E', 'S1E', 'S1F', 'S6B'}
    if need_scatter:
        L29NewOld_comp = comp(L29NewOld, 'L29 NewOld')
        S3NewOld_comp = comp(S3NewOld, 'S3 NewOld')
        L29Growth_comp = comp(L29Growth, 'L29 Growth')
        S3Growth_comp = comp(S3Growth, 'S3 Growth')

    if 'S6B' in want:
        print('\n=== Fig S6B: growth scores, Hippo pathway ===')
        three_ways(
            partial(screen_viz.plot_comparison_scatter,
                    L29Growth_comp, S3Growth_comp,
                    categories=cat('Cats_Hippo.csv'),
                    colors=HIPPO_COLORS, category_order=HIPPO_ORDER,
                    title='Growth Comparison Quiescent Screen Hippo',
                    figsize=(6, 4), show_correlation=False),
            PLOTS_ESS, 'GrowthL29_S3comparisonHippo')

    if 'S1E' in want:
        print('\n=== Fig S1E: growth scores, common essentials ===')
        three_ways(
            partial(screen_viz.plot_comparison_scatter,
                    L29Growth_comp, S3Growth_comp,
                    categories=cat('CommonEssentials.csv'),
                    title='Growth Comparison Quiescent Screen Common Essentials',
                    figsize=(6, 4), show_correlation=False),
            PLOTS_ESS, 'GrowthL29_S3comparisonEssentials')

    if 'S1F' in want:
        print('\n=== Fig S1F: new/old ratio, common essentials ===')
        three_ways(
            partial(screen_viz.plot_comparison_scatter,
                    L29NewOld_comp, S3NewOld_comp,
                    categories=cat('CommonEssentials.csv'),
                    title='NewOld Comparison Quiescent Screen Common Essentials',
                    figsize=(6, 4), show_correlation=False),
            PLOTS_ESS, 'NewOldL29_S3comparisonEssentials')

    if '1E' in want:
        print('\n=== Fig 1E: new/old ratio, three categories ===')
        cats3 = cat('GenomeWideCats_threeOptions.csv')
        three_ways(
            partial(screen_viz.plot_comparison_scatter,
                    L29NewOld_comp, S3NewOld_comp,
                    categories=cats3, colors=FIG1E_COLORS,
                    title='NewOld Comparison Quiescent Screen',
                    figsize=(5, 5), show_correlation=False),
            PLOTS_FIG1, 'NewOldL29_S3comparison3cats')
        print('  Fig 1E category sizes:')
        for c, n in cats3['category'].value_counts().items():
            print(f'    {c:24s} {n}')

    if 'S2A' in want:
        print('\n=== Fig S2A: category violins ===')
        catsMany = cat('manyCats2.csv')
        grey = {x: 'darkgrey' for x in catsMany.category.unique()}
        for screen, stem, title in ((S3NewOld, 'S3boxplot_many', '40S degradation pathways'),
                                    (L29NewOld, 'L29boxplot_many', '60S degradation pathways')):
            for ext in ('.png', '.pdf'):
                screen_viz.plot_category_enrichment(
                    results_file=os.path.join(RES, screen, screen + '_all_genes.csv'),
                    output_file=os.path.join(PLOTS_FIG1, stem + ext),
                    categories=catsMany,
                    hits_file=os.path.join(RES, screen, screen + '_hits_fdr0.05.csv'),
                    title=title, figsize=(13, 3.5), orientation='vertical',
                    show_legend=False, colors=grey,
                    category_order=MANYCATS_ORDER, reference_violins=True,
                    show_counts=True)
        # The per-category counts are now drawn on the panel itself
        # (show_counts=True), so the Illustrator labels no longer need pasting.
        # They are still printed -- for both screens -- so the numbers can be
        # quoted in the legend and cross-checked.
        import pandas as pd
        # Must mirror plot_category_enrichment: a gene listed under two
        # categories is counted in BOTH violins, so use the multi-map. A plain
        # dict(zip(...)) here silently under-counts every shared gene (it kept
        # only the last row).
        g2cs = screen_viz._gene_to_cats_map(catsMany)
        for screen, lab in ((L29NewOld, '60S (L29)'), (S3NewOld, '40S (S3)')):
            allg = pd.read_csv(os.path.join(RES, screen, screen + '_all_genes.csv'))
            is_p = screen_viz.is_ntc(allg['gene_id'])
            real = allg[~is_p].copy()
            real['_cats'] = real['gene'].astype(str).map(lambda g: g2cs.get(g, ()))
            print(f'  S2A category counts, {lab} screen:')
            print(f'    {"category":28s} {"symbols":>8s} {"pts drawn":>10s}')
            print(f'    {"All genes":28s} {real.gene.nunique():8d} {len(real):10d}')
            print(f'    {"NTC":28s} {int(is_p.sum()):8d} {int(is_p.sum()):10d}')
            for c in MANYCATS_ORDER:
                sub = real[real._cats.map(lambda cs: c in cs)]
                print(f'    {c:28s} {sub.gene.nunique():8d} {len(sub):10d}')
            print()

    if 'S1D' in want:
        print('\n=== Fig S1D: fold-change correlation heatmap ===')
        order = [('L29_TopvBottom', '60S New/Old (quiescent)'),
                 ('S3_TopvBottom', '40S New/Old (quiescent)'),
                 ('L29_TopvBottom_growing', '60S New/Old (growing)'),
                 ('S3_TopvBottom_growing', '40S New/Old (growing)'),
                 ('L29_UnsortedvEarly', '60S Growth Phenotype'),
                 ('S3_UnsortedvEarly', '40S Growth Phenotype')]
        screen_viz.plot_foldchange_heatmap(
            results_files=[os.path.join(RES, c, c + '_all_genes.csv') for c, _ in order],
            comparison_names=[n for _, n in order],
            output_file=os.path.join(PLOTS_FIG1, 'S1D_foldchange_correlation.pdf'),
            cluster=False, annot=False, mask_upper=True,
            cmap='RdBu_r', vmin=0.0, vmax=1.0,
            title='Pearson Correlation of $\\log_2$ Fold Change',
            figsize=(4.5, 3.6))

    if '1F' in want:
        print('\n=== Fig 1F: GO enrichment ===')
        out = os.path.join(PLOTS_FIG1, 'GO_enrichment_1F.pdf')
        df = screen_viz.plot_go_enrichment_diverging(
            screen_results=[('40S Screen', paths.inp('GSEA_Fig1/S3/final_resultsS3Genomewide.csv')),
                            ('60S Screen', paths.inp('GSEA_Fig1/L29/decoupler_results60S.csv'))],
            output_file=out, top_n=4,
            categories=cat('GenomeWideCats_threeOptions.csv'),
            colors=FIG1E_COLORS, term_labels=TERM_LABELS,
            break_range=(160, 60), figsize=(7.2, 5.0))
        df.to_csv(os.path.join(PLOTS_FIG1, 'GO_enrichment_1F_terms.csv'), index=False)

    print('\nDone.')


if __name__ == '__main__':
    main()
