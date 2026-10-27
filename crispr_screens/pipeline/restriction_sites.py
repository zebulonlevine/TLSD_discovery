#!/usr/bin/env python3
"""Sublibrary guide annotation and the NotI/BamHI exclusion list.

Single source of truth shared by build_sublibrary_tables.py, the fixed_NotI
re-scoring driver and the figure regenerators, so the exclusion cannot drift
between the tables and the scores computed from them.

A guide is excluded if `GCGGCCGC` (NotI) or `GGATCC` (BamHI) appears anywhere in
its full ordered oligo (`toOrder`), which includes sites formed at the junction
with the constant flanks. Both sites are palindromes, so one strand suffices.
That is 53 guides across 39 genes; 46 of them never passed the 50-read
threshold in any comparison and 16 have no reads at all, against a library
median of 127,792 reads per guide.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from common import paths  # noqa: E402

LIBINFO = paths.repo('crispr_screens', 'library')   # shipped with the code
MAGECK_COUNTS = paths.inp('Sublibrary_results/results/counts_for_mageck.txt')

RESTRICTION_SITES = {'NotI': 'GCGGCCGC', 'BamHI': 'GGATCC'}   # both palindromic
N_LIBRARY = 4928
N_SITE_GUIDES = 53


def library_annotation():
    """One row per library guide in sublibrary.csv order.

    Columns: sgrna, gene, spacer_sequence (the 20-mer count_spacers.py matched
    reads against), designed_spacer (the 20-nt protospacer as designed, column
    'sgRNA Sequence' of sublibrary.csv), has_site.

    The two differ for the 1,632 Dolcetto guides that were cloned with a 5' G
    prepended: their counted 20-mer is G + the first 19 nt of the design. The
    supplementary tables report designed_spacer (as Table S1 and S7 do).

    The sgrna names are reconstructed, not read back from the pipeline:
    prepare_mageck_input() names guides positionally, Gene_Single + '_sg' + a
    0-based counter within each Gene_Single group, in library row order
    (screen_pipeline.py:841-845). Re-running that stage today is NOT an option --
    a later NTC-pseudogene renaming block (lines 807-818) post-dates the
    sublibrary's 2026-01-15 run and would renumber every non-targeting guide,
    invalidating all current downstream results. So the rule is reproduced here
    and asserted against counts_for_mageck.txt.
    """
    ann = pd.read_csv(os.path.join(LIBINFO, 'sublibrary.csv'), low_memory=False)
    with open(os.path.join(LIBINFO, 'sublibrary_spacers.csv')) as fh:
        spacers = [ln.strip().upper() for ln in fh if ln.strip()]
    assert len(ann) == len(spacers) == N_LIBRARY, (len(ann), len(spacers))
    assert {len(s) for s in spacers} == {20}

    gene = ann['Gene_Single'].astype(str)
    assert not gene.isin(['nan', '']).any(), 'library has rows with no Gene_Single'

    oligo = ann['toOrder'].astype(str).str.upper()
    site = pd.Series(False, index=ann.index)
    for s in RESTRICTION_SITES.values():
        site |= oligo.str.contains(s, regex=False)
    assert site.sum() == N_SITE_GUIDES, f'{site.sum()} site guides, expected {N_SITE_GUIDES}'

    designed = ann['sgRNA Sequence'].astype(str).str.upper()
    assert set(designed.str.len()) == {20}
    # the counted 20-mer is the design itself, or G + its first 19 nt (+1G guides)
    sp = pd.Series(spacers)
    assert ((sp == designed) | (sp == 'G' + designed.str[:19])).all(), 'counted vs designed spacer drift'

    lib = pd.DataFrame({
        'sgrna': gene + '_sg' + ann.groupby('Gene_Single').cumcount().astype(str),
        'gene': gene,
        'spacer_sequence': spacers,
        'designed_spacer': designed.values,
        'has_site': site.values,
    })
    assert not lib['sgrna'].duplicated().any()
    names = pd.read_table(MAGECK_COUNTS, usecols=['sgRNA'])['sgRNA']
    assert (names.values == lib['sgrna'].values).all(), 'sgRNA name reconstruction drifted'
    return lib


def site_guides():
    """The 53 sgRNA names to exclude."""
    lib = library_annotation()
    return frozenset(lib.loc[lib['has_site'], 'sgrna'])


if __name__ == '__main__':
    lib = library_annotation()
    bad = lib[lib['has_site']]
    print(f'{len(lib)} library guides, {len(bad)} carry a NotI/BamHI site '
          f'across {bad.gene.nunique()} genes')
    print(bad[['sgrna', 'gene', 'spacer_sequence']].to_string(index=False))
