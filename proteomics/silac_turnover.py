"""SILAC pulse-chase proteomics in quiescent RPS3-HaloTag cells (sgNT vs sgDDRGK1).

Cells were fully heavy-labelled, chased into light medium, and collected 3 and 8 days after the
shift (3 replicates each), then multiplexed by TMT. Heavy = old protein, light = new protein.

  1. Each sample's heavy and light signals are normalized by that sample's total (H+L) signal.
  2. Replicate means and SEMs; day-8/day-3 ratios (old: heavy, new: light) with propagated
     error, log2.
  3. sgDDRGK1 - sgNT difference of the log2 ratios, t = difference / SEM, two-sided p from t with
     4 df.
  4. Localization categories from the organelle-IP table (nucleolus counted as nucleus; stress
     granule, P-body and translation counted as cytosol), ribosomal proteins by gene list.

Outputs (output root, proteomics/):
  TableS8_SILAC_annotated.csv / .xlsx
  fig3B_delta_old_by_category.pdf             + fig3B_delta_old_mwu.csv
  figS4B_intensity_violins.pdf
  figS4D_old_turnover_by_category.pdf,        + figS4DE_turnover_{old,new}_{NT,sgDDRGK1}_mwu.csv
  figS4E_new_turnover_by_category.pdf
  figS4F_er_gene_sets_strip.pdf               + figS4F_values.csv
  fig7C_new_protein_volcano.pdf               + fig7C_new_protein_stats.csv

    python proteomics/silac_turnover.py
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402
import silac_proteomics as sp  # noqa: E402
from settings import (  # noqa: E402
    SILAC_TSV, LOCALIZATION_XLSX, LOCALIZATION_SHEET, OUT_DIR,
    BAIT, CHANNELS, SAMPLE_NAMES, LABELS, RAW_LABEL_TOKENS, RAW_TEMPLATE, SAMPLE_TEMPLATE,
    TURNOVER_COMPARISONS, DELTA_COMPARISONS, T_TEST_DF,
    KEEP_CATEGORIES, CATEGORY_OVERRIDES, DESCRIPTION_RULES, CYTO_CATEGORIES, NUC_CATEGORIES, ADJ_CATEGORY,
    MWU, DELTA_SETS, BACKGROUND, SET_COLORS, FOCUS_GENE,
)

DAYS = ['3d', '8d']
REPS = [1, 2, 3]
ABUNDANCE_COL = 'rel_abundance'
SIZE_LABEL = 'abundance\n(total-normalized, a.u.)'
NEW_PROTEIN_FC = 'light_8d_sgDDRGK1_vs_NT_log2FC'


def out(name):
    return paths.outp(OUT_DIR, name, mkdir=True)


def load_localization():
    return pd.read_excel(paths.inp(LOCALIZATION_XLSX), sheet_name=LOCALIZATION_SHEET)


def adjusted_categories(df, loc_df, gene_col='Gene'):
    """Localization categories used in the figures and tables (see settings.ADJ_CATEGORY)."""
    df = sp.annotate_localization(df, loc_df, gene_col=gene_col, out_col=ADJ_CATEGORY)
    df = sp.assign_category_overrides(df, CATEGORY_OVERRIDES, gene_col=gene_col, out_col=ADJ_CATEGORY)
    merges = {**{c: 'cytosol' for c in CYTO_CATEGORIES}, **{c: 'nucleus' for c in NUC_CATEGORIES}}
    df[ADJ_CATEGORY] = df[ADJ_CATEGORY].replace(merges)
    return sp.collapse_categories(df, keep=KEEP_CATEGORIES, source_col=ADJ_CATEGORY, out_col=ADJ_CATEGORY)


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------

def build_silac(loc_df=None):
    """Protein table with normalized signals, turnover ratios, deltas and categories."""
    df = sp.drop_contaminants_and_decoys(sp.load_protein_quant(paths.inp(SILAC_TSV)))
    channel_map = sp.build_channel_map(raw_template=RAW_TEMPLATE, sample_template=SAMPLE_TEMPLATE,
                                       channels=CHANNELS, sample_names=SAMPLE_NAMES, labels=LABELS,
                                       raw_label_tokens=RAW_LABEL_TOKENS)
    df = sp.rename_channels(df, channel_map)
    sample_cols = list(channel_map.values())

    # each sample's H and L share one total-signal denominator
    df = sp.normalize_by_total(df, sp.pair_silac_columns(sample_cols, labels=LABELS), suffix='_norm', rescale=True)
    norm_cols = [c + '_norm' for c in sample_cols]
    groups = sp.replicate_groups_from_names(norm_cols, replicate_field=-3)   # RPS3_NT_3d_<rep>_H_norm
    df = sp.group_mean_sem(df, groups, mean_suffix='_avg', sem_suffix='_sem')

    for name, condition, late, early, label in TURNOVER_COMPARISONS:
        df = sp.ratio_with_error(df, f'{BAIT}_{condition}_{late}_{label}_norm_avg',
                                 f'{BAIT}_{condition}_{early}_{label}_norm_avg',
                                 num_sem=f'{BAIT}_{condition}_{late}_{label}_norm_sem',
                                 den_sem=f'{BAIT}_{condition}_{early}_{label}_norm_sem', out=name, log2=True)
    for name, a, b in DELTA_COMPARISONS:
        df = sp.difference_with_error(df, a, b, a_sem=a + '_sem', b_sem=b + '_sem', out=name)
        df = sp.two_sided_t_test(df, name, name + '_sem', df_resid=T_TEST_DF)

    # sgDDRGK1 / sgNT per day and label (used in the Fig S4C correlations)
    for day in DAYS:
        for label in LABELS:
            df = sp.ratio_with_error(df, f'{BAIT}_sgDDRGK1_{day}_{label}_norm_avg', f'{BAIT}_NT_{day}_{label}_norm_avg',
                                     num_sem=f'{BAIT}_sgDDRGK1_{day}_{label}_norm_sem',
                                     den_sem=f'{BAIT}_NT_{day}_{label}_norm_sem',
                                     out=f'log2_sgDDRGK1_vs_NT_{day}_{label}', log2=True)

    loc_df = load_localization() if loc_df is None else loc_df
    df = sp.annotate_localization(df, loc_df)
    df = sp.assign_category_overrides(df, CATEGORY_OVERRIDES)
    df = sp.annotate_by_description(df, DESCRIPTION_RULES)
    df = sp.collapse_categories(df, keep=KEEP_CATEGORIES, source_col='category', out_col='category_simplified')
    df = df.rename(columns={'Gene Symbol': 'Gene', 'Protein Id': 'Reference'})

    # Fig S4F point size: total-normalized H+L signal averaged over both days, one global scale
    abundance_cols = [f'{BAIT}_{c}_{t}_{lab}_norm_avg' for c in ('NT', 'sgDDRGK1') for t in DAYS for lab in LABELS]
    df[ABUNDANCE_COL] = (df[abundance_cols] / df[abundance_cols].stack().median()).mean(axis=1)

    return adjusted_categories(df, loc_df)


def new_protein_tests(df, days=DAYS):
    """Fig 7C: per protein and day, light (new) signal sgDDRGK1 vs sgNT on log2 replicate values,
    3 vs 3 Student t-test, BH across proteins quantified in all six channels."""
    light = {(g, d, r): f'{BAIT}_{g}_{d}_{r}_L_norm' for g in ('NT', 'sgDDRGK1') for d in DAYS for r in REPS}
    light_l2 = np.log2(df[list(light.values())] / df[list(light.values())].stack().median())
    rows = []
    for day in days:
        A = light_l2[[light[('NT', day, r)] for r in REPS]].to_numpy()
        B = light_l2[[light[('sgDDRGK1', day, r)] for r in REPS]].to_numpy()
        ok = np.isfinite(A).all(axis=1) & np.isfinite(B).all(axis=1)
        t, p = stats.ttest_ind(B, A, axis=1)
        padj = np.full(len(p), np.nan)
        padj[ok] = sp.benjamini_hochberg(p[ok])
        se = np.sqrt((A.var(axis=1, ddof=1) + B.var(axis=1, ddof=1)) / len(REPS))   # pooled SE, equal n
        rows.append(pd.DataFrame({'Gene': df.Gene, 'day': day, 'log2FC': B.mean(axis=1) - A.mean(axis=1),
                                  'se': se, 't': t, 'p': p, 'padj': padj})[ok])
    return pd.concat(rows)            # index = row of df


def table_S8(df, tests):
    """Table S8: identifiers, raw and normalized signals, turnover, deltas, Fig 7C test, annotations."""
    keep = [c for c in df.columns
            if c in ('Reference', 'Gene', 'Description', 'Group ID', 'heavy Peptides', 'light Peptides',
                     'Collapsed', 'Parsimony')
            or (c.startswith(BAIT + '_') and c.endswith(('_H', '_L', '_norm', '_norm_avg', '_norm_sem')))
            or c.startswith(('turnover_', 'delta_'))]
    t = df[keep].copy()
    for name, _, _ in DELTA_COMPARISONS:
        cols = list(t.columns)
        t.insert(cols.index(name + '_p') + 1, name + '_p_BH', sp.benjamini_hochberg(t[name + '_p']))
    d8 = tests[tests.day == '8d']
    for suffix, src in (('', 'log2FC'), ('_sem', 'se'), ('_t', 't'), ('_p', 'p'), ('_p_BH', 'padj')):
        t[NEW_PROTEIN_FC + suffix] = d8[src].reindex(t.index)
    t['category_localization'] = df[ADJ_CATEGORY]
    gene_to_set = {}
    for name, genes in DELTA_SETS.items():
        for g in genes:
            gene_to_set.setdefault(g, name)
    t['ER_GeneSets_delta'] = df['Gene'].map(gene_to_set)
    return t.sort_values('Gene', key=lambda g: g.astype(str).str.lower(), kind='stable')   # gene order


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def fig3B(df):
    """Change in old-protein turnover (sgDDRGK1 - sgNT) by localization."""
    sp.category_violin(df, category_col=ADJ_CATEGORY, value_col='delta_old_log2', categories=KEEP_CATEGORIES,
                       title=f'{BAIT} delta_old_log2 by localization (nucleolus in nucleus; '
                             'SG / P-body / translation in cytosol; no collagen label)',
                       ylabel='delta_old_log2', ylim=(-1, 1), save=out('fig3B_delta_old_by_category.pdf'))
    plt.close('all')
    sp.mann_whitney_subsampled(df, category_col=ADJ_CATEGORY, value_col='delta_old_log2',
                               categories=KEEP_CATEGORIES, **MWU).to_csv(out('fig3B_delta_old_mwu.csv'), index=False)


def figS4B(df):
    """Heavy (old) and light (new) intensity distributions, day 3 vs day 8."""
    plot_df = pd.melt(df, value_vars=[f'{BAIT}_{c}_{t}_{lab}_norm_avg' for lab in LABELS
                                      for c in ('NT', 'sgDDRGK1') for t in DAYS],
                      var_name='condition', value_name='intensity')
    plot_df['log_intensity'] = np.log10(plot_df['intensity'])
    plot_df['timepoint'] = plot_df['condition'].str.split('_').str[2]
    plot_df['Genetics_x_mass'] = plot_df['condition'].apply(lambda x: '_'.join([x.split('_')[1], x.split('_')[3]]))
    plt.figure()
    sns.violinplot(data=plot_df, x='Genetics_x_mass', y='log_intensity', hue='timepoint', split=True,
                   palette=['lightgrey', 'orange'], inner='quartile')
    plt.savefig(out('figS4B_intensity_violins.pdf'), dpi=300, bbox_inches='tight')
    plt.close('all')


def figS4DE(df):
    """Old (S4D) and new (S4E) protein day-8/day-3 ratios, sgNT vs sgDDRGK1, by localization."""
    for age, panel in (('old', 'figS4D'), ('new', 'figS4E')):
        nt, ko = f'turnover_{age}_NT_log2', f'turnover_{age}_sgDDRGK1_log2'
        path = out(f'{panel}_{age}_turnover_by_category.pdf')
        sp.category_violin(df, category_col=ADJ_CATEGORY, value_cols=[nt, ko], value_labels=('NT', 'sgDDRGK1'),
                           categories=KEEP_CATEGORIES, title=f'{BAIT} {age}-protein turnover: NT vs sgDDRGK1',
                           ylabel='log2(d8/d3)', ylim=(min(df[nt].min(), df[ko].min()), max(df[nt].max(), df[ko].max())),
                           save=path, split_palette=('lightgrey', 'orange'), median_line=False, figsize=(6, 3))
        plt.hlines(df[nt].median(), xmin=-0.5, xmax=9.5, linestyle='--', color='dimgrey', linewidth=1)
        plt.hlines(df[ko].median(), xmin=-0.5, xmax=9.5, linestyle='--', color='darkorange', linewidth=1)
        plt.savefig(path, dpi=300, bbox_inches='tight')
        plt.close('all')
        for cond, col in (('NT', nt), ('sgDDRGK1', ko)):
            sp.mann_whitney_subsampled(df, category_col=ADJ_CATEGORY, value_col=col, categories=KEEP_CATEGORIES,
                                       **MWU).to_csv(out(f'figS4DE_turnover_{age}_{cond}_mwu.csv'), index=False)


def figS4F(df):
    """Change in old (top) and new (bottom) protein for translocon complexes and ribosomal subunits."""
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.4), sharex=True, sharey=True)
    tidy_all = []
    for ax, (col, age) in zip(axes, [('delta_old_log2', 'old (heavy)'), ('delta_new_log2', 'new (light)')]):
        _, tidy = sp.gene_set_strip(df, value_col=col, gene_sets=DELTA_SETS, background=BACKGROUND,
                                    point_colors=SET_COLORS, size_col=ABUNDANCE_COL, size_label=SIZE_LABEL,
                                    size_legend=ax is axes[0], separator_after=['ER', 'EMC'], mean_bar=True,
                                    label_mode='outliers', title=f'$\\Delta$ {age}',
                                    ylabel='$\\Delta$log2(d8/d3)  sgDDRGK1 - sgNT', ylim=(-1.4, 1.2), ax=ax)
        tidy_all.append(tidy.assign(comparison=col))
    plt.savefig(out('figS4F_er_gene_sets_strip.pdf'), dpi=300, bbox_inches='tight')
    plt.close('all')
    values = pd.concat(tidy_all, ignore_index=True)
    values[~values.is_background].to_csv(out('figS4F_values.csv'), index=False)


def fig7C(tests, days=('8d',), path=None, fdr=0.05):
    """Volcano of new (light) protein, sgDDRGK1 vs sgNT, with TLSD1 (R3HDM4) marked."""
    fig, axes = plt.subplots(1, len(days), figsize=(2.4 * len(days), 3.0), sharex=True, sharey=True, squeeze=False)
    axes = axes[0]
    for ax, day in zip(axes, days):
        d = tests[tests.day == day]
        ax.scatter(d.log2FC, -np.log10(d.padj), s=3, color='#B3B3B3', alpha=0.35, linewidth=0, rasterized=True)
        ax.axhline(-np.log10(fdr), linestyle='--', color='#909090', linewidth=0.6)
        ax.axvline(0, color='#909090', linewidth=0.6)
        r = d[d.Gene == FOCUS_GENE].iloc[0]
        ax.scatter([r.log2FC], [-np.log10(r.padj)], s=26, facecolor='#E8871A', edgecolor='white', linewidth=0.5, zorder=5)
        ax.annotate(f'{FOCUS_GENE}\n{r.log2FC:+.2f}, padj {r.padj:.0e}', xy=(r.log2FC, -np.log10(r.padj)),
                    xytext=(-5, 12), textcoords='offset points', ha='right', fontsize=5, color='#404040',
                    arrowprops=dict(arrowstyle='-', color='#909090', linewidth=0.4))
        ax.set_title(f'{day}, new protein', fontsize=7.5)
        ax.set_xlabel('log2 FC (sgDDRGK1 / sgNT)', fontsize=7)
        ax.set_ylabel('-log10(BH-adjusted p)' if ax is axes[0] else '', fontsize=7)
        ax.tick_params(labelsize=6.5)
        sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(path or out('fig7C_new_protein_volcano.pdf'), dpi=300, bbox_inches='tight')
    plt.close('all')


def write_table(t, stem, sheet=None):
    t.to_csv(out(stem + '.csv'), index=False)
    if sheet:
        try:
            t.to_excel(out(stem + '.xlsx'), sheet_name=sheet, index=False)
        except ImportError:
            print('openpyxl not installed -- skipped', stem + '.xlsx')


def main():
    warnings.filterwarnings('ignore', message='After omitting NaNs')   # SEM of proteins seen in < 2 replicates = NaN
    df = build_silac()
    tests = new_protein_tests(df)
    tests[['Gene', 'day', 'log2FC', 'p', 'padj']].to_csv(out('fig7C_new_protein_stats.csv'), index=False)
    write_table(table_S8(df, tests), 'TableS8_SILAC_annotated', sheet='SupplementaryTable8_SILAC')
    fig3B(df)
    figS4B(df)
    figS4DE(df)
    figS4F(df)
    fig7C(tests)
    r = tests[(tests.day == '8d') & (tests.Gene == FOCUS_GENE)].iloc[0]
    print(f'{len(df)} proteins; {FOCUS_GENE} day-8 new protein: log2FC {r.log2FC:+.3f}, BH p {r.padj:.2g}')
    print('wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
