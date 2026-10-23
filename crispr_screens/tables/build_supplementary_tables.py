#!/usr/bin/env python3
"""
Supplementary table builders for the genome-wide ribosome-turnover screens.

Three merges over the pipeline's own output, one per table:

    Table S1  sgRNA level  -- counts / normalized means / LFC / p per comparison
    Table S2  gene level   -- epsilon / pvalue / product / hit flag per comparison,
                              plus the category annotations used by the figures
    Table S3  GO level     -- decoupler ULM results for both screens, side by side

The merges are the reusable part; the table layout is just column naming on top
of them. Written as CSV so nothing depends on openpyxl (absent from the
crispri-screen env).

Usage:
    python build_supplementary_tables.py --config config_genomewide.yaml
    python build_supplementary_tables.py --config config_genomewide.yaml --tables s2
"""

import argparse
import os
import sys

import pandas as pd
import numpy as np

# Reuse the pipeline's NTC definition rather than re-deriving it here.
_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, 'crispr_screens', 'pipeline'))
from common import paths  # noqa: E402
from screen_pipeline import load_config, is_ntc, is_ntc_pseudogene  # noqa: E402
import mageck_counts
import gene_aliases  # noqa: E402


CATEGORY_DIR = paths.repo('screen_figures', 'categories')

# The 8 comparisons that appear in the paper. The other 8 the pipeline scores
# (TopvUnsorted, UnsortedvBottom, +/- _growing) are intermediates and stay out.
#
# 'label' is the Table S1 naming ("Unsorted" is called "End" in the paper),
# 'prefix' is the Table S2 naming.
COMPARISONS = [
    dict(comp='L29_TopvBottom',             label='L29_TopvBottom',          prefix='L29_RatioNewOld_Quiescent'),
    dict(comp='S3_TopvBottom',              label='S3_TopvBottom',           prefix='S3_RatioNewOld_Quiescent'),
    dict(comp='L29_UnsortedvEarly',         label='L29_EndvEarly',           prefix='L29_Growth_Quiescent'),
    dict(comp='S3_UnsortedvEarly',          label='S3_EndvEarly',            prefix='S3_Growth_Quiescent'),
    dict(comp='L29_TopvBottom_growing',     label='L29_TopvBottom_growing',  prefix='L29_RatioNewOld_growing'),
    dict(comp='S3_TopvBottom_growing',      label='S3_TopvBottom_growing',   prefix='S3_RatioNewOld_growing'),
    dict(comp='L29_UnsortedvEarly_growing', label='L29_EndvEarly_growing',   prefix='L29_Growth_growing'),
    dict(comp='S3_UnsortedvEarly_growing',  label='S3_EndvEarly_growing',    prefix='S3_Growth_growing'),
]

# Table S1 column order, matching the existing workbook exactly so the CSV is a
# drop-in replacement. Derived from draft_paperFigure_actual/Table_S1_sgRNAs.xlsx.
# Sub-pool labels. The raw annotation carries the provider's internal name;
# these are what the guides are called in the literature.
LIBRARY_LABELS = {'Weissman': 'CRISPRi v2', 'Doench': 'Dolcetto'}

S1_COLUMN_ORDER = [
    'Gene', 'sgRNA', 'sgRNA Sequence', 'Library',
    'S3_Bottom_growing_count', 'S3_Early_growing_count', 'S3_Bottom_count', 'S3_Early_count',
    'L29_Bottom_growing_count', 'L29_Early_growing_count', 'L29_Bottom_count', 'L29_Early_count',
    'S3_Top_growing_count', 'S3_End_growing_count', 'S3_Top_count', 'S3_End_count',
    'L29_Top_growing_count', 'L29_End_growing_count', 'L29_Top_count', 'L29_End_count',
    'S3_Bottom_growing_mean', 'S3_Early_growing_mean', 'S3_Bottom_mean', 'S3_Early_mean',
    'L29_Bottom_growing_mean', 'L29_Early_growing_mean', 'L29_Bottom_mean', 'L29_Early_mean',
    'S3_Top_growing_mean', 'S3_End_growing_mean', 'S3_Top_mean', 'S3_End_mean',
    'L29_Top_growing_mean', 'L29_End_growing_mean', 'L29_Top_mean', 'L29_End_mean',
    'S3_TopvBottom_growing_LFC', 'S3_EndvEarly_growing_LFC', 'S3_TopvBottom_LFC', 'S3_EndvEarly_LFC',
    'L29_TopvBottom_growing_LFC', 'L29_EndvEarly_growing_LFC', 'L29_TopvBottom_LFC', 'L29_EndvEarly_LFC',
    'S3_TopvBottom_growing_p.twosided', 'S3_EndvEarly_growing_p.twosided',
    'S3_TopvBottom_p.twosided', 'S3_EndvEarly_p.twosided',
    'L29_TopvBottom_growing_p.twosided', 'L29_EndvEarly_growing_p.twosided',
    'L29_TopvBottom_p.twosided', 'L29_EndvEarly_p.twosided',
]

# sample-name suffix -> condition name used in the Table S1 column headers
CONDITION_NAMES = {
    'bottom': 'Bottom',
    'top': 'Top',
    'early': 'Early',
    'unsorted': 'End',
}

# Table S2 annotation columns -> category file
ANNOTATION_FILES = {
    'CommonEssential':        'CommonEssentials.csv',
    'Fig1FCategories':        'GenomeWideCats_threeOptions.csv',
    'FigS1FCategories':       'manyCats2.csv',
    'HippoPathway_FigureS5A': 'Cats_Hippo.csv',
}


# =============================================================================
# Helpers
# =============================================================================

def _condition(sample_names: str) -> str:
    """Map a comparison's treatment/control sample list to a condition name.

    'L29rep1_bottom,L29rep2_bottom' -> 'Bottom'. Replicates must agree.
    """
    suffixes = {s.strip().rsplit('_', 1)[-1].lower()
                for s in str(sample_names).split(',') if s.strip()}
    if len(suffixes) != 1:
        raise ValueError(f"mixed conditions in sample list: {sample_names!r}")
    suffix = suffixes.pop()
    if suffix not in CONDITION_NAMES:
        raise ValueError(f"unknown condition suffix {suffix!r} in {sample_names!r}")
    return CONDITION_NAMES[suffix]


def _line_and_state(comp: str) -> tuple:
    """'L29_TopvBottom_growing' -> ('L29', '_growing')."""
    line = comp.split('_')[0]
    state = '_growing' if comp.endswith('_growing') else ''
    return line, state


# Columns whose figure plots a DISTRIBUTION (violins), where a gene listed
# under two categories legitimately contributes to both. Everything else backs
# a scatter, where the point carries one colour.
MULTI_CATEGORY_COLUMNS = {'FigS1FCategories'}


def _load_category(filename: str, multi: bool = False) -> pd.Series:
    """gene -> category, from a two-column Categories CSV.

    A handful of genes are listed under more than one category (ATG16L1 is both
    Autophagy and E3; RACK1 is both Biogenesis and Lysosome/Autophagy). The
    resolution mirrors screen_viz:

    multi=False  one category, autophagy preferred, else last row wins
                 (screen_viz._gene_to_cat_map -- what the scatters colour by)
    multi=True   every category, '; '-joined
                 (screen_viz._gene_to_cats_map -- what the violins group by)
    """
    path = os.path.join(CATEGORY_DIR, filename)
    df = pd.read_csv(path)
    gene_col = 'gene' if 'gene' in df.columns else df.columns[0]
    cat_col = 'category' if 'category' in df.columns else df.columns[-1]
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df[cat_col] = df[cat_col].astype(str).str.strip()
    df = df[df[gene_col] != '']

    if multi:
        return (df.drop_duplicates(subset=[gene_col, cat_col])
                  .groupby(gene_col)[cat_col].apply('; '.join))

    out = {}
    for gene, cat in zip(df[gene_col], df[cat_col]):
        held = out.get(gene)
        if held is None or 'autophagy' not in held.lower():
            out[gene] = cat
    return pd.Series(out)


# =============================================================================
# Table S1 -- sgRNA level
# =============================================================================

def build_sgrna_table(results_dir: str, comparisons_table: str) -> pd.DataFrame:
    """Full sgRNA universe, left-joined with per-comparison MAGeCK statistics.

    Keyed on counts_for_mageck.txt (every annotated guide), NOT on the
    sgrna_summary files, which are already count-thresholded. Guides dropped by
    thresholding keep a row with blank statistics.
    """
    counts_file = os.path.join(results_dir, 'counts_for_mageck.txt')
    ann_file = os.path.join(results_dir, 'sgrna_annotation.csv')
    agg_file = os.path.join(results_dir, 'aggregated_counts.csv')

    counts = pd.read_table(counts_file, usecols=['sgRNA', 'Gene'])
    print(f"  sgRNA universe: {len(counts)} guides from counts_for_mageck.txt")

    if os.path.exists(ann_file):
        ann = pd.read_csv(ann_file)[['sgRNA', 'sgRNA Sequence']]
        table = counts.merge(ann, on='sgRNA', how='left')
        n_missing = table['sgRNA Sequence'].isna().sum()
        if n_missing:
            print(f"  WARNING: {n_missing} guides have no sequence in {ann_file}")
    else:
        # Pre-dates the sgrna_annotation.csv export; emit without sequence
        # rather than guessing via row order.
        print(f"  WARNING: {ann_file} not found -- 'sgRNA Sequence' will be blank. "
              f"Re-run the aggregate stage to generate it.")
        table = counts.copy()
        table['sgRNA Sequence'] = pd.NA

    # Which sub-pool each guide came from. aggregated_counts.csv is row-aligned
    # with counts_for_mageck.txt (both are written by the same aggregate stage,
    # before any guide is dropped), so this is a positional read, not a join on
    # sequence -- 14 spacers occur in BOTH pools and a sequence join would
    # collapse them onto one label.
    if os.path.exists(agg_file):
        agg = pd.read_csv(agg_file, usecols=['Library', 'sgRNA Sequence'], low_memory=False)
        assert len(agg) == len(table), (len(agg), len(table))
        assert (agg['sgRNA Sequence'].astype(str).str.upper().values
                == table['sgRNA Sequence'].astype(str).str.upper().values).all(), \
            'aggregated_counts.csv is not row-aligned with counts_for_mageck.txt'
        unknown = set(agg['Library']) - set(LIBRARY_LABELS)
        assert not unknown, f'unlabelled sub-pools: {sorted(unknown)}'
        table['Library'] = agg['Library'].map(LIBRARY_LABELS).values
        print('  sub-pools: ' + ', '.join(
            f'{k} = {v}' for k, v in table['Library'].value_counts().items()))
    else:
        print(f"  WARNING: {agg_file} not found -- 'Library' will be blank")
        table['Library'] = pd.NA

    table = table[['Gene', 'sgRNA', 'sgRNA Sequence', 'Library']]

    comps = pd.read_csv(comparisons_table).set_index('comparison_name')

    count_cols, mean_cols, stat_cols = [], [], []
    backfill_spec = {}
    for spec in COMPARISONS:
        comp, label = spec['comp'], spec['label']
        path = os.path.join(results_dir, comp, f'{comp}.sgrna_summary.txt')
        if not os.path.exists(path):
            print(f"  WARNING: missing {path} -- columns for {label} will be blank")
            continue

        line, state = _line_and_state(comp)
        treat_cond = _condition(comps.loc[comp, 'treatment'])
        ctrl_cond = _condition(comps.loc[comp, 'control'])

        sg = pd.read_table(path, usecols=['sgrna', 'control_count', 'treatment_count',
                                          'control_mean', 'treat_mean', 'LFC', 'p.twosided'])
        rename = {
            'treatment_count': f'{line}_{treat_cond}{state}_count',
            'control_count':   f'{line}_{ctrl_cond}{state}_count',
            'treat_mean':      f'{line}_{treat_cond}{state}_mean',
            'control_mean':    f'{line}_{ctrl_cond}{state}_mean',
            'LFC':             f'{label}_LFC',
            'p.twosided':      f'{label}_p.twosided',
        }
        sg = sg.rename(columns=rename).rename(columns={'sgrna': 'sgRNA'})

        count_cols += [rename['treatment_count'], rename['control_count']]
        split = lambda v: [x.strip() for x in str(v).split(',') if x.strip()]
        backfill_spec[rename['treatment_count']] = split(comps.loc[comp, 'treatment'])
        backfill_spec[rename['control_count']] = split(comps.loc[comp, 'control'])
        mean_cols += [rename['treat_mean'], rename['control_mean']]
        stat_cols += [rename['LFC'], rename['p.twosided']]

        # Count/mean columns repeat across comparisons that share a sample
        # (none do among these 8, but guard rather than silently duplicating).
        new = [c for c in sg.columns if c == 'sgRNA' or c not in table.columns]
        table = table.merge(sg[new], on='sgRNA', how='left')
        print(f"  {comp:30s} -> {label}")

    # Counts exist for every guide whether or not a comparison scored it, so
    # fill the cells mageck never emitted. Everything downstream of counting
    # (mean/LFC/p) stays NaN -- that is what marks a guide as unscored.
    raw = pd.read_table(counts_file).set_index('sgRNA').drop(columns='Gene').astype(float)
    filled = 0
    for col, samples in backfill_spec.items():
        if col not in table.columns:
            continue
        before = int(table[col].isna().sum())
        series = mageck_counts.backfill_counts(table.set_index('sgRNA')[col], samples, raw)
        table[col] = table[col].fillna(table['sgRNA'].map(series))
        assert table[col].notna().all(), f'{col} still has gaps'
        filled += before
    print(f'  backfilled {filled} empty count cells across '
          f'{len(backfill_spec)} columns ({len(table)} guides)')

    built = set(count_cols) | set(mean_cols) | set(stat_cols) | {
        'Gene', 'sgRNA', 'sgRNA Sequence', 'Library'}
    missing = [c for c in S1_COLUMN_ORDER if c not in table.columns]
    extra = sorted(built - set(S1_COLUMN_ORDER))
    if missing:
        print(f"  WARNING: expected columns not produced: {missing}")
    if extra:
        print(f"  WARNING: produced columns not in S1_COLUMN_ORDER: {extra}")
    ordered = [c for c in S1_COLUMN_ORDER if c in table.columns] + extra
    return table[ordered]


# =============================================================================
# Table S2 -- gene level
# =============================================================================

def build_gene_table(results_dir: str) -> pd.DataFrame:
    """Outer join of the 8 comparisons' gene scores, plus figure categories."""
    table = None
    for spec in COMPARISONS:
        comp, prefix = spec['comp'], spec['prefix']
        all_path = os.path.join(results_dir, comp, f'{comp}_all_genes.csv')
        hits_path = os.path.join(results_dir, comp, f'{comp}_hits_fdr0.05.csv')
        if not os.path.exists(all_path):
            print(f"  WARNING: missing {all_path} -- skipping {prefix}")
            continue

        df = pd.read_csv(all_path)[['gene_id', 'gene', 'epsilon', 'pvalue', 'product']]
        hit_ids = set()
        if os.path.exists(hits_path):
            hits = pd.read_csv(hits_path)
            hit_ids = set(hits['gene_id'])
        else:
            print(f"  WARNING: missing {hits_path} -- {prefix}_hit will be all False")

        df[f'{prefix}_hit'] = df['gene_id'].isin(hit_ids)
        df = df.rename(columns={'epsilon': f'{prefix}_epsilon',
                                'pvalue': f'{prefix}_pvalue',
                                'product': f'{prefix}_product'})
        table = df if table is None else table.merge(
            df, on=['gene_id', 'gene'], how='outer')
        print(f"  {comp:30s} -> {prefix:32s} ({len(df)} rows, {len(hit_ids)} hits)")

    if table is None:
        raise RuntimeError("no gene-score files found")

    # A gene absent from one comparison has no hit call there; False (not NaN)
    # is the honest reading -- it was not called a hit.
    for spec in COMPARISONS:
        col = f"{spec['prefix']}_hit"
        if col in table.columns:
            table[col] = table[col].fillna(False).astype(bool)

    # Figure category annotations, from the same files the figures read.
    for out_col, filename in ANNOTATION_FILES.items():
        try:
            mapping = _load_category(filename, multi=out_col in MULTI_CATEGORY_COLUMNS)
        except FileNotFoundError:
            print(f"  WARNING: {filename} not found -- {out_col} left blank")
            table[out_col] = pd.NA
            continue
        table[out_col] = table['gene'].astype(str).str.strip().map(mapping)
        print(f"  {out_col:24s} <- {filename:38s} ({table[out_col].notna().sum()} genes annotated)")

    ordered = ['gene_id', 'gene']
    for spec in COMPARISONS:
        p = spec['prefix']
        ordered += [c for c in (f'{p}_epsilon', f'{p}_pvalue', f'{p}_product', f'{p}_hit')
                    if c in table.columns]
    ordered += list(ANNOTATION_FILES.keys())
    return table[ordered]


# =============================================================================
# Table S3 -- GO / decoupler ULM
# =============================================================================

def build_go_table(s3_file: str, l29_file: str) -> pd.DataFrame:
    """Side-by-side decoupler ULM results for the 40S (RPS3) and 60S (RPL29) screens.

    Both source columns already exist upstream -- FDR_adj_p_value and
    neglog10_FDR -- so this is a rename, not a derivation. Verified against the
    existing workbook: all four numeric columns match to <1e-13.
    """
    keep = {
        't_statistic': 't_statistic',
        'FDR_adj_p_value': 'FDR_adj_p_value',
        'neglog10_FDR': 'neg_log10_pvalue',
        'n_overlap': 'n_overlap',
        'genes_in_set_and_screen': 'genes_in_set_and_screen',
    }
    out = None
    for path, suffix in ((s3_file, 'RPS3'), (l29_file, 'RPL29')):
        if not os.path.exists(path):
            print(f"  WARNING: missing {path} -- {suffix} columns omitted")
            continue
        df = pd.read_csv(path)
        df = df[['Gene_Set'] + list(keep)].rename(
            columns={src: f'{dst}_{suffix}' for src, dst in keep.items()})
        out = df if out is None else out.merge(df, on='Gene_Set', how='outer')
        print(f"  {os.path.basename(path):34s} -> *_{suffix} ({len(df)} gene sets)")

    if out is None:
        raise RuntimeError("no decoupler result files found")

    ordered = ['Gene_Set']
    for suffix in ('RPS3', 'RPL29'):
        ordered += [f'{dst}_{suffix}' for dst in keep.values()
                    if f'{dst}_{suffix}' in out.columns]
    return out[ordered]


# =============================================================================
# CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True, help='pipeline config YAML')
    ap.add_argument('--results-dir', default=paths.inp('GenomeWide_results/results'),
                    help='genome-wide pipeline results (aggregated counts, MAGeCK and gene tables)')
    ap.add_argument('--outdir', default=paths.outp('GenomeWide_results/supplementary_tables'),
                    help='default: <results parent>/supplementary_tables')
    ap.add_argument('--tables', default='s1,s2,s3',
                    help='comma-separated subset: s1,s2,s3')
    ap.add_argument('--gene-labels', choices=['aliased', 'plain'], default='aliased',
                    help="'aliased' applies the combined display labels "
                         "(TLSD1/R3HDM4 ...) to Table S1's Gene column; 'plain' "
                         "leaves the raw annotation symbol. The manuscript draft "
                         "copy uses plain, the pipeline copy uses aliased -- "
                         "see claude_docs/METHODS_DRAFT_data_analysis.md 6b.")
    ap.add_argument('--go-s3', default=paths.inp('GSEA_Fig1/S3/final_resultsS3Genomewide.csv'))
    ap.add_argument('--go-l29', default=paths.inp('GSEA_Fig1/L29/decoupler_results60S.csv'))
    args = ap.parse_args()

    config = load_config(args.config)
    results_dir = args.results_dir
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    wanted = {t.strip().lower() for t in args.tables.split(',')}

    if 's1' in wanted:
        print("\n=== Table S1: sgRNA level ===")
        t1 = build_sgrna_table(results_dir, config.comparisons_table)
        if args.gene_labels == 'aliased':
            t1 = gene_aliases.apply_to(t1, suffixed_cols=('Gene',))
        else:
            print('  Gene labels: plain (combined display labels NOT applied)')
        p1 = os.path.join(outdir, 'Table_S1_sgRNAs.csv')
        t1.to_csv(p1, index=False)
        print(f"  wrote {p1}  ({len(t1)} rows x {len(t1.columns)} cols)")

    if 's2' in wanted:
        print("\n=== Table S2: gene level ===")
        t2 = build_gene_table(results_dir)
        # label column only; gene_id stays unaliased as the join key
        t2 = gene_aliases.apply_to(t2, plain_cols=('gene',))
        p2 = os.path.join(outdir, 'TableS2_GeneLevel_genomewide.csv')
        t2.to_csv(p2, index=False)
        print(f"  wrote {p2}  ({len(t2)} rows x {len(t2.columns)} cols)")
        n_pseudo = int(is_ntc_pseudogene(t2['gene_id']).sum())
        print(f"  NTC pseudogenes: {n_pseudo}; real genes: {len(t2) - n_pseudo}")
        for g in ('KNTC1', 'NEGR1'):
            hit = t2[t2['gene'].astype(str) == g]
            print(f"  {g}: {len(hit)} row(s)"
                  + (f" -- present as a real gene" if len(hit) else " -- ABSENT"))

    if 's3' in wanted:
        print("\n=== Table S3: GO / decoupler ULM ===")
        t3 = build_go_table(args.go_s3, args.go_l29)
        p3 = os.path.join(outdir, 'TableS3_GO_Genomewide_screens.csv')
        t3.to_csv(p3, index=False)
        print(f"  wrote {p3}  ({len(t3)} rows x {len(t3.columns)} cols)")

    print("\nDone.")


if __name__ == '__main__':
    main()
