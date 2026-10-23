"""
Fig5 -- growth-screen epistasis (sgUBA5 vs sgNT), using "unsorted v early" as the
growth phenotype.

Generates, for each reporter line (L29 = 60S, S3 = 40S):
  * growth(sgNT) vs growth(sgUBA5) scatter          -> growth_epistasis_{60S,40S}_scatter.pdf
  * waterfall of epistasis (sgUBA5 - sgNT)          -> growth_epistasis_{60S,40S}_waterfall.pdf
and for the phenotype averaged between the two reporter lines:
  * averaged per-condition tables                   -> growth_avgLines_{sgNT,sgUBA5}_all_genes.csv
  * averaged growth phenotype waterfall (sgNT)      -> growth_avgLines_waterfall_phenotype.pdf
  * averaged scatter + epistasis waterfall          -> growth_avgLines_{scatter,waterfall_epistasis}.pdf
plus a between-line QC scatter that motivates the averaging.

Fig S7E = growth_avgLines_scatter*.pdf, Fig S7F = growth_avgLines_waterfall_epistasis.pdf
(written to fixed_NotI/plots/growth/ under the output root).
Run:  python figS7EF_growth_epistasis.py      (envs/crispri-screen.yml)
"""
import os
import sys

import pandas as pd

if 'ipykernel' not in sys.modules:            # headless when run from the CLI
    import matplotlib
    matplotlib.use('Agg')
from matplotlib import pyplot as plt

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, REPO)
LIB = os.path.join(REPO, 'screen_figures', 'lib')
sys.path.insert(0, LIB)
from common import paths  # noqa: E402
sys.modules.pop('screen_viz', None)
import screen_viz  # noqa: E402
assert os.path.dirname(os.path.abspath(screen_viz.__file__)) == LIB
print('screen_viz:', screen_viz.__file__)

# Growth gene tables from the re-scored tree (53 NotI/BamHI guides excluded).
RES = paths.inp('fixed_NotI', 'results')
print('growth results:', RES)
CATS = paths.repo('screen_figures', 'categories')
# Every output below is a bare filename; write them all to one folder.
os.makedirs(paths.outp('fixed_NotI', 'plots', 'growth'), exist_ok=True)
os.chdir(paths.outp('fixed_NotI', 'plots', 'growth'))

# =========================================================================
# Categories: biogenesis split by subunit (3 greens) + Lysosome/Autophagy +
# UFMylation + EIF3 + R3HDM4. Identical construction to the compare_*_bioBySubunit
# cells. Row order matters: plot_waterfall takes its category order (and therefore
# its draw order) from order of first appearance, so R3HDM4 goes last = on top.
# =========================================================================
rna_pol_I = ['POLR1A', 'POLR1B', 'POLR1C', 'POLR1D', 'POLR1E', 'POLR1F', 'POLR1G', 'POLR1H']

bio_bysub = screen_viz.load_categories(os.path.join(CATS, 'GenomeWideCats_synth_bySubunit.csv'))[['gene', 'category']].copy()
bio_bysub['category'] = bio_bysub['category'].replace({'Both Biogenesis': 'Both subunits'})
bio_bysub.loc[bio_bysub['gene'].isin(rna_pol_I), 'category'] = 'Both subunits'

_base = screen_viz.load_categories(os.path.join(CATS, 'GenomeWideCats_threeOptions.csv'))[['gene', 'category']]
_non_bio = _base.loc[_base['category'] != 'Biogenesis']
categories_bysub = pd.concat([bio_bysub, _non_bio], ignore_index=True)

_lyso_manual = ['RNASEK', 'ATG14']
categories_bysub = categories_bysub[~categories_bysub['gene'].isin(_lyso_manual)]
categories_bysub = pd.concat(
    [categories_bysub, pd.DataFrame({'gene': _lyso_manual,
                                     'category': ['Lysosome/Autophagy'] * len(_lyso_manual)})],
    ignore_index=True)

_eif3 = ['EIF3A', 'EIF3B', 'EIF3C', 'EIF3D', 'EIF3E', 'EIF3F', 'EIF3G',
         'EIF3H', 'EIF3I', 'EIF3J', 'EIF3K', 'EIF3L', 'EIF3M']   # only D/F/H/M are in the sublibrary
_extra = _eif3 + ['R3HDM4']
categories_bysub = categories_bysub[~categories_bysub['gene'].isin(_extra)]
categories_bysub = pd.concat(
    [categories_bysub, pd.DataFrame({'gene': _extra,
                                     'category': ['EIF3'] * len(_eif3) + ['R3HDM4']})],
    ignore_index=True)

# plot_waterfall takes its category (legend + draw) order from order of first
# appearance in this table, so sort the rows into bio_order to match the scatters.
categories_bysub = categories_bysub.sort_values(
    'category', key=lambda c: c.map({k: i for i, k in enumerate(
        ['Both subunits', '60S Biogenesis/Proteins', '40S Biogenesis/Proteins',
         'Lysosome/Autophagy', 'UFMylation', 'EIF3', 'R3HDM4'])}),
    kind='stable').reset_index(drop=True)

green_colors = {
    'Both subunits'           : '#00441b',  # dark green (includes RNA Pol I)
    '60S Biogenesis/Proteins' : '#238b45',  # medium green
    '40S Biogenesis/Proteins' : '#74c476',  # light green
    'Lysosome/Autophagy'      : '#c51b7d',
    'UFMylation'              : '#E69F00',
    'EIF3'                    : '#5CB1E3',
    'R3HDM4'                  : '#231F20',  # black, matches Fig5 cell 10
}

# R3HDM4 drawn fully opaque; every other category keeps the default 0.7.
cat_alpha = {'R3HDM4': 1.0}
bio_order = ['Both subunits', '60S Biogenesis/Proteins', '40S Biogenesis/Proteins',
             'Lysosome/Autophagy', 'UFMylation', 'EIF3', 'R3HDM4']


def _top_each_direction(cx, cy, n=5):
    """The n most extreme genes in each direction on each axis."""
    dx = screen_viz.deduplicate_by_max_abs_epsilon(cx.results, 'gene_id', 'epsilon')
    dy = screen_viz.deduplicate_by_max_abs_epsilon(cy.results, 'gene_id', 'epsilon')
    m = dx[['gene_id', 'gene', 'epsilon']].merge(dy[['gene_id', 'epsilon']], on='gene_id', suffixes=('_x', '_y'))
    m = m[~m['gene_id'].astype(str).str.contains('NTC', case=False)]
    genes = set()
    for col in ['epsilon_x', 'epsilon_y']:
        genes |= set(m.nlargest(n, col)['gene'])
        genes |= set(m.nsmallest(n, col)['gene'])
    return sorted(genes)


def _growth_files(prefix, cond):
    """(all_genes.csv, hits_fdr0.05.csv) for one line/condition."""
    d = os.path.join(RES, f'{prefix}_unsorted_v_early_{cond}')
    stem = f'{prefix}_unsorted_v_early_{cond}'
    return os.path.join(d, f'{stem}_all_genes.csv'), os.path.join(d, f'{stem}_hits_fdr0.05.csv')


SCATTER_KW = dict(
    categories=categories_bysub, title=None,
    show_diagonal=True, show_correlation=False, figsize=(4.6, 3.0),
    legend_outside=True, size_by_hits=False, point_size=30, fontsize=8,
    dpi=600, colors=green_colors, category_order=bio_order,
    category_alpha=cat_alpha,
)
WATERFALL_KW = dict(
    categories=categories_bysub, figsize=(5.5, 3.0), fontsize=8,
    legend_outside=True, colors=green_colors, category_alpha=cat_alpha,
)

# =========================================================================
# 1. Per reporter line: growth vs growth scatter + epistasis waterfall
# =========================================================================
LINES = [('60S', 'L29'), ('40S', 'S3')]

for label, prefix in LINES:
    wt_all, wt_hits = _growth_files(prefix, 'wt')
    ub_all, ub_hits = _growth_files(prefix, 'sgUBA5')

    comp_wt = screen_viz.load_comparison(wt_all, wt_hits, f'Growth, sgNT ({label} line)')
    comp_ub = screen_viz.load_comparison(ub_all, ub_hits, f'Growth, sgUBA5 ({label} line)')
    label_list = sorted(set(_top_each_direction(comp_wt, comp_ub, 5)) | {'R3HDM4'})

    for fname, leg in [(f'growth_epistasis_{label}_scatter.pdf', True),
                       (f'growth_epistasis_{label}_scatter_noLegend.pdf', False)]:
        screen_viz.plot_comparison_scatter(comp_wt, comp_ub, fname,
                                           label_genes=label_list, adjust_labels=True,
                                           show_legend=leg, **SCATTER_KW)
        plt.close('all')

    # Epistasis = sgUBA5 - sgNT, ranked low to high. NTCs (715 of 1430 rows) are
    # kept in black as the null distribution for the difference.
    screen_viz.plot_waterfall_diff(
        ub_all, wt_all, f'growth_epistasis_{label}_waterfall.pdf',
        value_col='epsilon', name_a='sgUBA5', name_b='sgNT',
        ylabel=f'Δ growth ε, sgUBA5 − sgNT ({label} line)',
        label_top_n=8, label_genes=['R3HDM4'],
        diff_csv=f'growth_epistasis_{label}.csv', **WATERFALL_KW)
    plt.close('all')

# =========================================================================
# 2. Growth phenotype averaged between the two reporter lines
#    epsilon = mean(60S, 40S); pvalue = max of the two (conservative, so the
#    column still means "worst case"); hits = called in BOTH lines.
# =========================================================================
KEEP = ['gene_id', 'epsilon', 'pvalue', 'gene', 'product']


def average_lines(cond, out_stem):
    a_all, a_hits = _growth_files('L29', cond)
    b_all, b_hits = _growth_files('S3', cond)
    a = pd.read_csv(a_all)[KEEP]
    b = pd.read_csv(b_all)[KEEP]
    m = a.merge(b, on=['gene_id', 'gene'], suffixes=('_60S', '_40S'))
    m['epsilon'] = m[['epsilon_60S', 'epsilon_40S']].mean(axis=1)
    m['product'] = m[['product_60S', 'product_40S']].mean(axis=1)
    m['pvalue'] = m[['pvalue_60S', 'pvalue_40S']].max(axis=1)
    m = m[['gene_id', 'epsilon', 'pvalue', 'gene', 'product',
           'epsilon_60S', 'epsilon_40S']].sort_values('epsilon').reset_index(drop=True)

    hits = set(pd.read_csv(a_hits)['gene_id']) & set(pd.read_csv(b_hits)['gene_id'])
    h = m[m['gene_id'].isin(hits)]

    all_csv, hits_csv = f'{out_stem}_all_genes.csv', f'{out_stem}_hits_fdr0.05.csv'
    m.to_csv(all_csv, index=False)
    h.to_csv(hits_csv, index=False)
    print(f'Saved: {all_csv} ({len(m)} genes), {hits_csv} ({len(h)} hits in both lines)')
    return all_csv, hits_csv


avg_wt_all, avg_wt_hits = average_lines('wt', 'growth_avgLines_sgNT')
avg_ub_all, avg_ub_hits = average_lines('sgUBA5', 'growth_avgLines_sgUBA5')

comp_avg_wt = screen_viz.load_comparison(avg_wt_all, avg_wt_hits, 'Growth, sgNT (mean of lines)')
comp_avg_ub = screen_viz.load_comparison(avg_ub_all, avg_ub_hits, 'Growth, sgUBA5 (mean of lines)')
avg_labels = sorted(set(_top_each_direction(comp_avg_wt, comp_avg_ub, 5)) | {'R3HDM4'})

for fname, leg in [('growth_avgLines_scatter.pdf', True),
                   ('growth_avgLines_scatter_noLegend.pdf', False)]:
    screen_viz.plot_comparison_scatter(comp_avg_wt, comp_avg_ub, fname,
                                       label_genes=avg_labels, adjust_labels=True,
                                       show_legend=leg, **SCATTER_KW)
    plt.close('all')

# The averaged growth phenotype itself (sgNT), ranked. log2(1.5) guides kept here
# since this is an LFC-like quantity, unlike the difference plots.
screen_viz.plot_waterfall(
    avg_wt_all, 'growth_avgLines_waterfall_phenotype.pdf', hits_file=avg_wt_hits,
    ylabel='Growth ε, mean of 60S & 40S lines (sgNT)',
    label_top_n=8, label_genes=['R3HDM4'], **WATERFALL_KW)
plt.close('all')

screen_viz.plot_waterfall_diff(
    avg_ub_all, avg_wt_all, 'growth_avgLines_waterfall_epistasis.pdf',
    value_col='epsilon', name_a='sgUBA5', name_b='sgNT',
    ylabel='Δ growth ε, sgUBA5 − sgNT (mean of lines)',
    label_top_n=8, label_genes=['R3HDM4'],
    diff_csv='growth_epistasis_avgLines.csv', **WATERFALL_KW)
plt.close('all')

# =========================================================================
# 3. QC: the two reporter lines agree on growth, which is what justifies (2).
# =========================================================================
qc_wt_60 = screen_viz.load_comparison(*_growth_files('L29', 'wt'), 'Growth sgNT (60S line)')
qc_wt_40 = screen_viz.load_comparison(*_growth_files('S3', 'wt'), 'Growth sgNT (40S line)')
screen_viz.plot_comparison_scatter(
    qc_wt_60, qc_wt_40, 'growth_QC_60Svs40S_sgNT_scatter.pdf',
    label_genes=sorted(set(_top_each_direction(qc_wt_60, qc_wt_40, 5)) | {'R3HDM4'}),
    adjust_labels=True, show_legend=True,
    **{**SCATTER_KW, 'show_correlation': True})
plt.close('all')

print('done')
