#!/usr/bin/env python3
"""Table S6: GO-term enrichment (decoupler ULM) for the four sublibrary line
scores, pooled into one sheet.

One row per GO term that passes the size filter in the sublibrary (1,074).
`n_overlap` / `genes_in_set_and_screen` are identical for all four scores (they
share the 715-gene universe), so they are taken once. Rows are ordered by the
smallest FDR across the four scores, ties (BH produces ~355) by term name.

This reproduces the layout of the hand-assembled draft workbook
`TableS6_GO_pooled_sublibraryScores.xlsx`. Rebuilt from the July GO tables
(`--suffix '' --check <draft.xlsx>`) every term's values match the draft; only
the order *within* tied FDRs differs, because the draft was sorted by hand.

    python build_tableS6.py        # published: re-scored GO tables (*_fixedNotI.csv)
    python build_tableS6.py --go-dir <dir> --suffix '' --check <draft.xlsx>

Writes <out>.csv always and <out>.xlsx when openpyxl is importable (the
crispri-screen env lacks it; the mageck env has it).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from common import paths  # noqa: E402

GO_DIR = paths.inp('fixed_NotI/GO_decoupler')
SCORES = ['40S_Synthesis', '40S_Degradation', '60S_Synthesis', '60S_Degradation']
SHEET = 'GO enrichment'


def build(go_dir, suffix):
    tabs = {s: pd.read_csv(os.path.join(go_dir, f'{s}_GO_decoupler{suffix}.csv'))
               .set_index('Gene_Set') for s in SCORES}
    base = tabs[SCORES[0]]
    for s, t in tabs.items():
        if not t.index.sort_values().equals(base.index.sort_values()):
            sys.exit(f'{s}: GO term set differs from {SCORES[0]}')
        for c in ('n_overlap', 'genes_in_set_and_screen'):
            if not (t.loc[base.index, c] == base[c]).all():
                sys.exit(f'{s}: {c} differs from {SCORES[0]} -- gene universes disagree')
    out = base[['n_overlap', 'genes_in_set_and_screen']].copy()
    for s in SCORES:
        lab = s.replace('_', ' ')
        out[f't_statistic ({lab})'] = tabs[s]['t_statistic']
        out[f'FDR_adj_p_value ({lab})'] = tabs[s]['FDR_adj_p_value']
    fdr = out[[c for c in out.columns if c.startswith('FDR_adj_p_value')]].min(axis=1)
    out = (out.assign(_k=fdr.values).reset_index()
              .sort_values(['_k', 'Gene_Set'], kind='mergesort').drop(columns='_k'))
    return out.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--go-dir', default=GO_DIR)
    ap.add_argument('--suffix', default='_fixedNotI',
                    help="GO csv suffix: '_fixedNotI' (published), '_pRanked' (superseded), '' (July)")
    ap.add_argument('--out', default=paths.outp('fixed_NotI/supplementary_tables/TableS6_GO_pooled_sublibraryScores'),
        help='output path prefix (.csv / .xlsx appended)')
    ap.add_argument('--check', metavar='XLSX',
                    help='compare the build against an existing workbook instead of writing')
    args = ap.parse_args()

    t = build(args.go_dir, args.suffix)
    if args.check:
        ref = pd.read_excel(args.check, sheet_name=SHEET)
        same_cols = list(ref.columns) == list(t.columns)
        same_terms = set(ref['Gene_Set']) == set(t['Gene_Set']) and len(ref) == len(t)
        r, b = ref.set_index('Gene_Set'), t.set_index('Gene_Set').loc[ref['Gene_Set']]
        num = [c for c in t.columns if c.startswith(('t_statistic', 'FDR'))]
        dmax = float(np.nanmax(np.abs(r[num].values - b[num].values)))
        txt = all((r[c] == b[c]).all() for c in ('n_overlap', 'genes_in_set_and_screen'))
        fdr = r[[c for c in num if c.startswith('FDR')]].min(axis=1).values
        ref_sorted = bool(np.all(np.diff(fdr) >= 0))
        print(f'columns identical: {same_cols} | same {len(t)} terms: {same_terms} | '
              f'n_overlap/genes identical: {txt} | max |numeric diff|: {dmax:.3g} | '
              f'reference also ordered by min FDR: {ref_sorted}')
        return 0 if (same_cols and same_terms and txt and dmax < 1e-12 and ref_sorted) else 1

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    t.to_csv(args.out + '.csv', index=False)
    print(f'wrote {args.out}.csv  ({len(t)} terms)')
    try:
        import openpyxl  # noqa: F401
    except ImportError:
        print('openpyxl not importable -- .xlsx skipped (run with the mageck env for it)')
        return 0
    t.to_excel(args.out + '.xlsx', sheet_name=SHEET, index=False)
    print(f'wrote {args.out}.xlsx')
    return 0


if __name__ == '__main__':
    sys.exit(main())
