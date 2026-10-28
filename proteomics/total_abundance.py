"""Total-protein TMT proteomics: sgDDRGK1 vs sgNT in growing and quiescent cells.

The input table holds log2 TMT abundances, one column per sample (3 replicates each of RPS3-HaloTag
sgNT / sgDDRGK1, growing / quiescent, and RPL29-HaloTag sgDDRGK1 quiescent).

  1. Each sample is median-centred to the mean of the sample medians.
  2. Replicate means and SEMs; log2FC = difference of the means, SEM by root-sum-square,
     t = log2FC / SEM, two-sided p from t with 4 df, BH across proteins.
  3. Localization categories (as in silac_turnover.py) and membrane topology from
     annotation/uniprot_topology_pipeline.py (UniProt features + MitoCarta 3.0).

Outputs (output root, proteomics/):
  TableS9_nonSILAC_proteomics.csv / .xlsx
  fig4A_figS6A_total_by_category.pdf          + fig4A_figS6A_{growing,quiescent}_mwu.csv
  fig4E_figS6C_total_by_topology.pdf          + fig4E_figS6C_{growing,quiescent}_mwu.csv
  figS6D_volcano_quiescent.pdf, figS6E_volcano_growing.pdf

    python proteomics/total_abundance.py
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

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
from common import paths  # noqa: E402
import silac_proteomics as sp  # noqa: E402
from silac_turnover import adjusted_categories, load_localization, write_table  # noqa: E402
from settings import (  # noqa: E402
    TOTAL_CSV, TOTAL_TOPOLOGY, OUT_DIR, TOTAL_COLUMNS, TOTAL_REPLICATES, TOTAL_COMPARISONS, T_TEST_DF,
    QUIESCENT_FC, GROWING_FC, KEEP_CATEGORIES, CATEGORY_OVERRIDES, DESCRIPTION_RULES, ADJ_CATEGORY,
    MWU, TOPOLOGY_CATEGORIES, TOPOLOGY_MAP, TRANSLOCON, UPR_ER,
)

CONDITIONS = {'quiescent': QUIESCENT_FC, 'growing': GROWING_FC}
GENE_SET_COL = 'isTransloconOrUPR'


def out(name):
    return paths.outp(OUT_DIR, name, mkdir=True)


def uniprot_accession(reference):
    """'sp|P12931-2|SRC_HUMAN' -> 'P12931-2'."""
    return reference.str.split('|').apply(lambda x: x[1] if len(x) > 1 else x[0])


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------

def build_total(loc_df=None, topology=True):
    """Protein table with normalized abundances, comparisons, categories and topology."""
    df = pd.read_csv(paths.repo(TOTAL_CSV))
    df.columns = TOTAL_COLUMNS

    loc_df = load_localization() if loc_df is None else loc_df
    df = sp.annotate_localization(df, loc_df, gene_col='Gene')
    df = sp.assign_category_overrides(df, CATEGORY_OVERRIDES, gene_col='Gene')
    df = sp.annotate_by_description(df, DESCRIPTION_RULES, description_col='Annotation')
    df = sp.collapse_categories(df, keep=KEEP_CATEGORIES, source_col='category', out_col='category_simplified')

    # values are already log2: centre each sample's median on the mean of the medians
    mean_median = df[TOTAL_REPLICATES].median().mean()
    for c in TOTAL_REPLICATES:
        df[c + '_norm'] = df[c] - df[c].median() + mean_median
    groups = sp.replicate_groups_from_names([c + '_norm' for c in TOTAL_REPLICATES], replicate_field=-2)
    df = sp.group_mean_sem(df, groups, mean_suffix='_avg', sem_suffix='_sem')

    for a, b in TOTAL_COMPARISONS:
        fc = f'{a}_vs_{b}_log2FC'
        df = sp.difference_with_error(df, a=a + '_norm_avg', b=b + '_norm_avg', out=fc,
                                      a_sem=a + '_norm_sem', b_sem=b + '_norm_sem')
        df = sp.two_sided_t_test(df, estimate=fc, sem=fc + '_sem', df_resid=T_TEST_DF,
                                 out_t=fc + '_t', out_p=fc + '_p')
        df[fc + '_p_BH'] = sp.benjamini_hochberg(df[fc + '_p'])

    if topology:
        df = add_topology(df)
    df[GENE_SET_COL] = 'Other'
    df.loc[df.Gene.isin(TRANSLOCON), GENE_SET_COL] = 'Translocon'
    df.loc[df.Gene.isin(UPR_ER), GENE_SET_COL] = 'UPR'
    return adjusted_categories(df, loc_df)


def add_topology(df):
    """Membrane topology class per protein; '_mito' appended for MitoCarta 3.0 proteins."""
    annot = pd.read_csv(paths.repo(TOTAL_TOPOLOGY), sep='\t')
    annot.loc[annot.in_mitocarta3, 'class'] = annot['class'][annot.in_mitocarta3] + '_mito'
    df['accession'] = uniprot_accession(df.Reference)
    df = pd.merge(df, annot[['input_id', 'class']].rename(columns={'class': 'topology'}),
                  left_on='accession', right_on='input_id', how='left')
    df['simplifiedTopology'] = df['topology'].map(TOPOLOGY_MAP).fillna('')
    return df.drop(columns=['accession', 'input_id'])


def table_S9(df):
    """Table S9: identifiers, raw and normalized abundances, comparisons, annotations."""
    keep = (['Reference', 'Gene', 'Annotation'] + TOTAL_REPLICATES + [c + '_norm' for c in TOTAL_REPLICATES]
            + [c for c in df.columns if c.endswith(('_norm_avg', '_norm_sem'))]
            + [c for c in df.columns if '_vs_' in c]
            + ['simplifiedTopology', GENE_SET_COL, ADJ_CATEGORY])
    t = df[keep].rename(columns={'Gene': 'Gene Symbol'})
    return t.sort_values(QUIESCENT_FC + '_t', kind='stable')          # quiescent t, most decreased first


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def split_box(df, group_col, order, title, path):
    """sgDDRGK1 vs sgNT log2FC per group, quiescent (grey) next to growing (blue);
    dashed lines = proteome median per condition."""
    cols, labels = list(CONDITIONS.values()), list(CONDITIONS)
    plot_df = pd.melt(df[df[group_col].isin(order)], id_vars=group_col, value_vars=cols,
                      var_name='condition', value_name='log2fc')
    plot_df['condition'] = plot_df['condition'].map(dict(zip(cols, labels)))
    fig, ax = plt.subplots(figsize=(6, 3))
    sns.boxplot(data=plot_df, x=group_col, y='log2fc', hue='condition', hue_order=labels, order=order,
                palette=['lightgrey', 'lightblue'], width=0.7, fliersize=2, linewidth=0.8, ax=ax)
    ax.set_title(title)
    ax.set_ylabel('log2 fold change sgDDRGK1 vs sgNT')
    ax.set_ylim(-2, 2)
    ax.axhline(df[cols[0]].median(), linestyle='--', color='dimgrey', linewidth=1)
    ax.axhline(df[cols[1]].median(), linestyle='--', color='lightblue', linewidth=1)
    plt.setp(ax.get_xticklabels(), rotation=90, ha='right')
    sns.despine(ax=ax)
    ax.legend(frameon=False)
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close('all')


def mwu_by_group(df, group_col, groups, stem):
    """Each group vs the rest of the proteome, per condition (subsampled MWU, BH)."""
    for cond, col in CONDITIONS.items():
        sp.mann_whitney_subsampled(df, category_col=group_col, value_col=col, categories=groups,
                                   **MWU).to_csv(out(f'{stem}_{cond}_mwu.csv'), index=False)


def fig4A_S6A(df, category_col=ADJ_CATEGORY):
    """Total protein change by localization; Fig 4A shows the RPL / RPS groups."""
    split_box(df, category_col, KEEP_CATEGORIES, 'changes quiescent vs growing by localization',
              out('fig4A_figS6A_total_by_category.pdf'))
    mwu_by_group(df, category_col, KEEP_CATEGORIES, 'fig4A_figS6A')


def fig4E_S6C(df):
    """Total protein change by membrane topology."""
    split_box(df, 'simplifiedTopology', TOPOLOGY_CATEGORIES, 'changes quiescent vs growing by membrane topology',
              out('fig4E_figS6C_total_by_topology.pdf'))
    mwu_by_group(df, 'simplifiedTopology', TOPOLOGY_CATEGORIES, 'fig4E_figS6C')


def figS6DE(df):
    """Volcanoes, quiescent (S6D) and growing (S6E): Sec61 translocon and ER UPR genes highlighted;
    dashed lines at 1.5-fold and p = 0.05."""
    for panel, cond in (('figS6D', 'quiescent'), ('figS6E', 'growing')):
        col = CONDITIONS[cond]
        sp.volcano(df, effect_col=col, pvalue_col=col + '_p', title=f'S3 {cond} log2FC sgDDRGK1 vs NT',
                   category_col=GENE_SET_COL, highlight={'Other': 'lightgrey', 'Translocon': 'blue', 'UPR': 'orange'},
                   base_size=27, highlight_size=40.5)          # 1.5x the default marker diameter (area x 2.25)
        plt.vlines([np.log2(1.5), np.log2(1 / 1.5)], ymin=0, ymax=7, linestyles='dashed')
        plt.hlines([-np.log10(0.05)], xmin=-6, xmax=2.5, linestyles='dashed')
        plt.xlim(-6, 2.5)
        # Points are rasterized (sp.volcano); dpi sets their resolution (default 100 is too coarse).
        # TrueType text (fonttype 42) stays editable in Illustrator; matplotlib's default Type 3 does not.
        with matplotlib.rc_context({'pdf.fonttype': 42, 'ps.fonttype': 42}):
            plt.savefig(out(f'{panel}_volcano_{cond}.pdf'), dpi=300)
        plt.close('all')


def main():
    warnings.filterwarnings('ignore', message='After omitting NaNs')   # SEM of proteins seen in < 2 replicates = NaN
    df = build_total()
    write_table(table_S9(df), 'TableS9_nonSILAC_proteomics', sheet='SupplementaryTable9')   # xlsx keeps MARCH5 etc. as text
    fig4A_S6A(df)
    fig4E_S6C(df)
    figS6DE(df)
    print(f'{len(df)} proteins; wrote', paths.outp(OUT_DIR))


if __name__ == '__main__':
    main()
