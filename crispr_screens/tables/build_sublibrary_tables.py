#!/usr/bin/env python3
"""Supplementary Tables S4, S5, S10, S11 for the sublibrary screens.

  S4   sgRNA level, sgNT arm    -- counts / means / LFC / p per comparison,
                                   plus the four line scores AND (new) their
                                   per-sgRNA p-values
  S5   gene level, sgNT arm     -- epsilon / pvalue / product / hit per
                                   comparison and per line score
  S10  sgRNA level, sgUBA5 arm  -- S4 with every column suffixed _sgUBA5
  S11  gene level, sgUBA5 arm   -- S5 for sgUBA5, plus sgUBA5-minus-sgNT
                                   epistasis columns for the four line scores

Provenance (all verified against the previous hand-built tables):
  *_Top/Bottom_count      <- mageck treatment_count / control_count (kept as the
                             slash-joined per-replicate strings)
  *_Top/Bottom_mean       <- treat_mean / control_mean
  *_LFC, *_pvalue         <- LFC, p.twosided
  *_score                 <- line score from <scores>/<dir>/*_sgrna_scores.csv
  *_score_pvalue          <- NEW: per-sgRNA p-value for the score
  *_lfc_{epsilon,...}     <- results/<comp>/<comp>_all_genes.csv (+ hits file)
  *_score_{epsilon,...}   <- <scores>/<dir>/gene_scores_all_genes.csv (+ hits)

RatioON column labelling
------------------------
`*_ratioON_top_v_bottom_*` is defined with treatment=Bottom, control=Top --
reversed relative to the new/old comparisons despite its name. Traced to the raw
counts to be sure: mageck's treatment_count for that comparison matches the raw
ratioON_Bottom samples to within normalisation (~1.00x) while Top is 3-5x off.

The previous S4/S10 put treatment_count under `RatioON_Top_count`, so those
columns held BOTTOM counts and the `RatioON_TopvBottom_LFC` name described the
opposite of what it contained. Default here is `--ratioon-labels fixed`, which
puts Top counts under Top and names the fold change `RatioON_BottomvTop_LFC`.
No number changes -- only labels. Use `--ratioon-labels legacy` to reproduce the
old (incorrect) headers exactly.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, 'crispr_screens', 'pipeline'))
from common import paths  # noqa: E402
import gene_aliases  # noqa: E402
import mageck_counts  # noqa: E402
import restriction_sites  # noqa: E402

# Published tables are built from the re-scored tree (53 NotI/BamHI guides excluded).
RES = paths.inp('fixed_NotI/results')
COMPARISONS = paths.repo('crispr_screens/pipeline/comparisons_sublibrary.csv')
MAGECK_COUNTS = restriction_sites.MAGECK_COUNTS

# line-score dir -> (subunit label for gene tables, line label for sgRNA tables)
SCORES = [('L29_Degradation', '60S_Degradation', 'L29_Degradation'),
          ('L29_Synthesis',   '60S_Synthesis',   'L29_Synthesis'),
          ('S3_Degradation',  '40S_Degradation', 'S3_Degradation'),
          ('S3_Synthesis',    '40S_Synthesis',   'S3_Synthesis')]

# per-comparison blocks: (line, stem, sgRNA-table prefix, gene-table prefix)
BLOCKS = [('L29', 'new_top_v_bottom',     'L29_New',      '60S_New_TopvBottom'),
          ('L29', 'old_top_v_bottom',     'L29_Old',      '60S_Old_TopvBottom'),
          ('L29', 'ratioON_top_v_bottom', 'L29_RatioON',  '60S_RatioON_TopvBottom'),
          ('L29', 'unsorted_v_early',     'L29',          '60S_UnsortedvEarly'),
          ('S3',  'new_top_v_bottom',     'S3_New',       '40S_New_TopvBottom'),
          ('S3',  'old_top_v_bottom',     'S3_Old',       '40S_Old_TopvBottom'),
          ('S3',  'ratioON_top_v_bottom', 'S3_RatioON',   '40S_RatioON_TopvBottom'),
          ('S3',  'unsorted_v_early',     'S3',           '40S_UnsortedvEarly')]


def arms(comp):
    tab = pd.read_csv(COMPARISONS).set_index('comparison_name')
    row = tab.loc[comp]
    split = lambda v: [x.strip() for x in str(v).split(',') if x.strip()]
    return split(row['treatment']), split(row['control'])


def treatment_is_bottom(comp):
    t, _ = arms(comp)
    return all('_Bottom_' in s or s.endswith('_Bottom') for s in t)


def library_frame():
    """The guides that belong in S4/S10, in sublibrary.csv row order.

    Columns: sgrna, gene, spacer_sequence. Guides whose ordered oligo carries a
    NotI/BamHI site are dropped, so both tables are built on one identical row
    set. See restriction_sites.library_annotation() for the name-reconstruction
    rule and why re-running the pipeline to get these names is not an option.

    spacer_sequence is the DESIGNED 20-nt protospacer (2026-10-23; previously the
    counted 20-mer, which for the 1,632 +1G Dolcetto guides is G + 19 nt), so the
    sublibrary tables match Table S1 and Table S7.
    """
    lib = restriction_sites.library_annotation()
    lib = lib[~lib['has_site']].drop(columns='has_site').reset_index(drop=True)
    lib['spacer_sequence'] = lib.pop('designed_spacer')
    assert len(lib) == restriction_sites.N_LIBRARY - restriction_sites.N_SITE_GUIDES
    return lib


def raw_counts():
    """Unfiltered per-sample counts for all library guides, indexed by sgRNA."""
    return pd.read_table(MAGECK_COUNTS).set_index('sgRNA').drop(columns='Gene').astype(float)


# =============================================================================
# sgRNA level (S4 / S10)
# =============================================================================

def build_sgrna_table(cond, scores_root, ratioon_labels='fixed', lib=None):
    sfx = '' if cond == 'wt' else '_sgUBA5'
    ct = sfx                      # condition tag inserted before quantity markers
    out = None
    gene_maps = []
    raw = raw_counts()
    backfill = {}
    for line, stem, pre, _ in BLOCKS:
        comp = f'{line}_{stem}_{cond}'
        path = os.path.join(RES, comp, f'{comp}.sgrna_summary.txt')
        if not os.path.exists(path):
            print(f'  WARNING: missing {path}'); continue
        sg = pd.read_table(path, usecols=['sgrna', 'Gene', 'control_count', 'treatment_count',
                                          'control_mean', 'treat_mean', 'LFC', 'p.twosided'])
        if stem == 'unsorted_v_early':
            a, b, lfc_name = 'Unsorted', 'Early', f'{pre}_UnsortedvEarly'
        elif stem == 'ratioON_top_v_bottom' and treatment_is_bottom(comp):
            if ratioon_labels == 'fixed':
                # treatment is Bottom: put it under Bottom, and name the fold
                # change for the direction it actually measures.
                a, b, lfc_name = 'Bottom', 'Top', f'{pre}_BottomvTop'
            else:
                a, b, lfc_name = 'Top', 'Bottom', f'{pre}_TopvBottom'   # legacy, mislabelled
        else:
            a, b, lfc_name = 'Top', 'Bottom', f'{pre}_TopvBottom'

        # Condition tag sits immediately after the sample/comparison identifier
        # and BEFORE the quantity marker (_count/_mean/_LFC/_pvalue), so the
        # header reads L29_New_Top_sgUBA5_count rather than
        # L29_New_Top_count_sgUBA5.
        ren = {'treatment_count': f'{pre}_{a}{ct}_count', 'control_count': f'{pre}_{b}{ct}_count',
               'treat_mean': f'{pre}_{a}{ct}_mean', 'control_mean': f'{pre}_{b}{ct}_mean',
               'LFC': f'{lfc_name}{ct}_LFC', 'p.twosided': f'{lfc_name}{ct}_pvalue'}
        # Collect this block's Gene labels; they are coalesced across ALL blocks
        # below. Taking the label from only one source (the score files, or a
        # single comparison) leaves NaN for the handful of guides that source
        # happens to drop.
        sg = sg.rename(columns=ren)
        # the column named for `a` always holds the treatment arm (see `ren`)
        t_samples, c_samples = arms(comp)
        for qcol, reps in ((f'{pre}_{a}{ct}_count', t_samples),
                           (f'{pre}_{b}{ct}_count', c_samples)):
            backfill[qcol] = mageck_counts.backfill_counts(
                sg.set_index('sgrna')[qcol], reps, raw)
        gene_maps.append(sg[['sgrna', 'Gene']].dropna())
        out = sg.drop(columns=['Gene']) if out is None else \
            out.merge(sg.drop(columns=['Gene']), on='sgrna', how='outer')

    # line scores + the NEW per-sgRNA score p-values
    for sdir, _, slab in SCORES:
        d = os.path.join(scores_root, sdir + sfx if cond != 'wt' else sdir)
        f = glob.glob(os.path.join(d, '*_sgrna_scores.csv'))
        if not f:
            print(f'  WARNING: no sgrna_scores.csv in {d}'); continue
        sc = pd.read_csv(f[0])
        cols = {'score': f'{slab}{ct}_score'}
        if 'pvalue' in sc.columns:
            cols['pvalue'] = f'{slab}{ct}_score_pvalue'
        else:
            print(f'  NOTE: {os.path.basename(f[0])} has no per-sgRNA pvalue column')
        keep = ['sgrna'] + [c for c in cols if c in sc.columns]
        sc = sc[keep].rename(columns=cols)
        out = out.merge(sc, on='sgrna', how='outer')

    # coalesce the gene label across every comparison that saw the guide
    if gene_maps:
        gmap = pd.concat(gene_maps, ignore_index=True).drop_duplicates('sgrna')
        out = out.merge(gmap.rename(columns={'Gene': '_gene'}), on='sgrna', how='left')
        out['gene'] = out['_gene'] if 'gene' not in out.columns else \
            out['gene'].fillna(out['_gene'])
        out = out.drop(columns='_gene')

    # Reindex onto the library. This restores the guides no comparison could
    # score -- they failed the 50-read per-replicate-pair threshold everywhere,
    # which is a sampling fact, not a library fact -- and drops the
    # restriction-site guides, leaving S4 and S10 on one identical row set.
    if lib is None:
        lib = library_frame()
    removed = sorted(set(out['sgrna']) - set(lib['sgrna']))
    print(f'  restored {len(lib) - len(out) + len(removed)} under-sampled guides, '
          f'removed {len(removed)} restriction-site guides: {removed}')
    out = lib.merge(out, on='sgrna', how='left', suffixes=('_lib', ''))
    out['gene'] = out['gene'].fillna(out['gene_lib'])
    out = out.drop(columns='gene_lib')

    # Counts exist for every guide whether or not it was scored, so fill the
    # cells mageck never emitted. Everything downstream of counting stays NaN.
    for col, series in backfill.items():
        out[col] = out[col].fillna(out['sgrna'].map(series))
        assert out[col].notna().all(), f'{col} still has gaps'

    n_null = out['gene'].isna().sum()
    if n_null:
        print(f'  WARNING: {n_null} guides still have no gene label')

    key = ['sgrna', 'gene', 'spacer_sequence']
    return out[key + [c for c in out.columns if c not in key]]


# =============================================================================
# gene level (S5 / S11)
# =============================================================================

def _gene_block(path_all, prefix):
    if not os.path.exists(path_all):
        print(f'  WARNING: missing {path_all}'); return None
    d = pd.read_csv(path_all)
    hits_path = path_all.replace('_all_genes.csv', '_hits_fdr0.05.csv')
    hits = set(pd.read_csv(hits_path)['gene_id']) if os.path.exists(hits_path) else set()
    d = d[['gene_id', 'gene', 'epsilon', 'pvalue', 'product']].copy()
    d[f'{prefix}_hit'] = d['gene_id'].isin(hits)
    return d.rename(columns={'epsilon': f'{prefix}_epsilon', 'pvalue': f'{prefix}_pvalue',
                             'product': f'{prefix}_product'})


def build_gene_table(cond, scores_root):
    # Condition tag before the _lfc / _score marker, e.g.
    # 60S_New_TopvBottom_sgUBA5_lfc_epsilon, 60S_Degradation_sgUBA5_score_hit.
    ct = '' if cond == 'wt' else '_sgUBA5'
    out = None
    for line, stem, _, gpre in BLOCKS:
        comp = f'{line}_{stem}_{cond}'
        blk = _gene_block(os.path.join(RES, comp, f'{comp}_all_genes.csv'), f'{gpre}{ct}_lfc')
        if blk is None: continue
        out = blk if out is None else out.merge(blk, on=['gene_id', 'gene'], how='outer')
    for sdir, glab, _ in SCORES:
        d = os.path.join(scores_root, sdir if cond == 'wt' else sdir + '_sgUBA5')
        blk = _gene_block(os.path.join(d, 'gene_scores_all_genes.csv'), f'{glab}{ct}_score')
        if blk is None: continue
        out = out.merge(blk, on=['gene_id', 'gene'], how='outer')
    return out


def add_epistasis(s11, s5):
    """sgUBA5 - sgNT epsilon difference for each line score.

    Left NaN for the NTC_{j} pseudogene rows: each condition draws its own
    pseudogenes independently, so NTC_17 in the sgUBA5 table and NTC_17 in the
    sgNT table are different random guide sets and their difference carries no
    meaning. Half the rows (715 of 1430) are pseudogenes, so including them
    would also distort any summary statistic taken over the column.
    """
    is_pseudo = s11['gene_id'].astype(str).str.startswith('NTC_')
    for _, glab, _ in SCORES:
        ub_col = f'{glab}_sgUBA5_score_epsilon'     # tagged, in S11
        nt_col = f'{glab}_score_epsilon'            # untagged, in S5
        m = s5[['gene_id', nt_col]].rename(columns={nt_col: '_nt'})
        s11 = s11.merge(m, on='gene_id', how='left')
        diff = s11[ub_col] - s11['_nt']
        s11[f'{glab}_score_epistasis_sgUBA5_minus_sgNT'] = diff.where(~is_pseudo)
        s11 = s11.drop(columns='_nt')
    return s11


def main():
    global RES
    ap = argparse.ArgumentParser()
    ap.add_argument('--scores-root',
                    default=paths.inp('fixed_NotI/scores'))
    ap.add_argument('--outdir',
                    default=paths.outp('fixed_NotI/supplementary_tables'))
    ap.add_argument('--results-root', default=RES,
                    help='per-comparison mageck/U-test results; point this at a '
                         're-scored tree (e.g. fixed_NotI/results) to rebuild '
                         'S5/S11 from it')
    ap.add_argument('--ratioon-labels', choices=['fixed', 'legacy'], default='fixed')
    ap.add_argument('--tag', default='')
    args = ap.parse_args()
    RES = args.results_root
    os.makedirs(args.outdir, exist_ok=True)
    t = args.tag

    # one shared row set, so S4 and S10 cannot drift apart
    lib = library_frame()
    print(f'=== library: {len(lib)} guides after dropping {restriction_sites.N_SITE_GUIDES} '
          f'NotI/BamHI oligos ===')

    print('=== Table S4: sgRNA level, sgNT ===')
    s4 = build_sgrna_table('wt', args.scores_root, args.ratioon_labels, lib=lib)
    s4 = gene_aliases.apply_to(s4, plain_cols=('gene',))
    s4.to_csv(os.path.join(args.outdir, f'TableS4_sgRNA_merged{t}.csv'), index=False)
    print(f'  {s4.shape[0]} rows x {s4.shape[1]} cols')

    print('=== Table S5: gene level, sgNT ===')
    s5 = build_gene_table('wt', args.scores_root)
    s5 = gene_aliases.apply_to(s5, plain_cols=('gene',))
    s5.to_csv(os.path.join(args.outdir, f'TableS5_GeneLevel_sublib{t}.csv'), index=False)
    print(f'  {s5.shape[0]} rows x {s5.shape[1]} cols')

    print('=== Table S10: sgRNA level, sgUBA5 ===')
    s10 = build_sgrna_table('sgUBA5', args.scores_root, args.ratioon_labels, lib=lib)
    s10 = gene_aliases.apply_to(s10, plain_cols=('gene',))
    s10.to_csv(os.path.join(args.outdir, f'TableS10_sgRNA_sgUBA5_merged{t}.csv'), index=False)
    print(f'  {s10.shape[0]} rows x {s10.shape[1]} cols')

    print('=== Table S11: gene level, sgUBA5 + epistasis ===')
    s11 = build_gene_table('sgUBA5', args.scores_root)
    s11 = add_epistasis(s11, s5)
    s11 = gene_aliases.apply_to(s11, plain_cols=('gene',))
    s11.to_csv(os.path.join(args.outdir, f'TableS11_GeneLevel_sublib_sgUBA5{t}.csv'), index=False)
    print(f'  {s11.shape[0]} rows x {s11.shape[1]} cols')
    ep = [c for c in s11.columns if 'epistasis' in c]
    print(f'  epistasis columns: {ep}')
    for c in ep:
        print(f'    {c:52s} n={s11[c].notna().sum():5d} '
              f'mean={s11[c].mean():+.4f} sd={s11[c].std():.4f}')


if __name__ == '__main__':
    main()
