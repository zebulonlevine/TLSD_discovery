"""Display aliases for gene symbols in the supplementary tables.

Applied to the human-readable gene LABEL columns only, at write time, after all
joins. `gene_id` is deliberately left unaliased: it is the join key used against
the score files, the hits files and the category tables, and the audit's identity
checks merge on it.

Genome-wide labels carry an isoform suffix (ACAD11_P1P2), so the alias is applied
to the leading token and the suffix preserved.
"""

GENE_ALIASES = {
    'R3HDM4': 'TLSD1/R3HDM4',
    'ACAD11': 'UBA5/ACAD11',
    'ANKZF1': 'ATG9A/ANKZF1',
    # This gene is annotated under DIFFERENT symbols in the two libraries --
    # NEPRO in the genome-wide annotation, C3orf17 in the sublibrary -- and both
    # are retired symbols for HGNC:24496, whose current symbol is RMP64 (the
    # label the figures use). Both map to RMP64 so the gene reads identically in
    # every table and figure; gene_id keeps the library's own symbol.
    # (2026-10-23: was 'C3orf17/NEPRO'.)
    'C3orf17': 'RMP64',
    'NEPRO':   'RMP64',
}


def alias_symbol(s):
    """'ACAD11' -> 'UBA5/ACAD11'; anything else unchanged."""
    return GENE_ALIASES.get(str(s), s)


def alias_suffixed(s):
    """'ACAD11_P1P2' -> 'UBA5/ACAD11_P1P2'; leading token aliased, suffix kept."""
    s = str(s)
    head, sep, tail = s.partition('_')
    if head in GENE_ALIASES:
        return GENE_ALIASES[head] + sep + tail
    return s


def apply_to(df, plain_cols=(), suffixed_cols=()):
    """Return df with the named label columns aliased. Does not touch gene_id."""
    out = df.copy()
    for c in plain_cols:
        if c in out.columns:
            out[c] = out[c].map(alias_symbol)
    for c in suffixed_cols:
        if c in out.columns:
            out[c] = out[c].map(alias_suffixed)
    return out


def report(df, cols):
    n = {}
    for c in cols:
        if c in df.columns:
            n[c] = int(df[c].astype(str).str.contains('/', regex=False).sum())
    return n
