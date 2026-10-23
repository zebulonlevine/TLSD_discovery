#!/usr/bin/env python3
"""Put sub-threshold guides' read counts back on MAGeCK's normalised scale.

The `*_count` columns of Tables S1/S4/S10 hold MAGeCK's median-normalised
counts as slash-joined per-replicate strings. Guides that failed a comparison's
minimum-count filter never reach MAGeCK, so those cells are empty even though
the guide was sequenced.

MAGeCK normalises by a single scalar per sample, so `normalised / raw` is
constant across guides -- measured spread over the scored guides of
`L29_new_top_v_bottom_wt` is <=1e-4 relative, which is just its
5-significant-figure output. Recovering that scalar lets the missing cells be
filled on the same scale as the rest of the column, rather than mixing raw and
normalised numbers in one column.

Round-trip check on the sublibrary: over 306,912 guide x sample cells, 0.74%
of reconstructions differ from MAGeCK's own string and none by more than one
unit in the 5th significant figure.
"""
import numpy as np
import pandas as pd

MAGECK_SIGFIGS = 5


def size_factor(norm, raw, max_spread=1e-3):
    """Median of normalised/raw over the guides MAGeCK scored.

    `norm` and `raw` are Series on the same index; guides MAGeCK dropped are NaN
    in `norm`. A spread above `max_spread` means the sample mapping or the
    replicate order is wrong, which must fail loudly rather than silently
    producing plausible-looking counts.
    """
    ok = (raw > 0) & norm.notna()
    assert ok.sum() > 100, f'only {int(ok.sum())} guides to fit a size factor'
    r = norm[ok] / raw[ok]
    f = float(r.median())
    spread = float((r.max() - r.min()) / f)
    assert spread < max_spread, f'size factor not constant (relative spread {spread:.2e})'
    return f


def backfill_counts(norm_col, samples, raw):
    """Normalised counts for every guide in `raw`, in mageck's own format.

    `norm_col`   Series of mageck's values, indexed by guide name (NaN where the
                 guide was thresholded out). Slash-joined strings when the
                 comparison has several replicates, plain floats when it has one.
    `samples`    raw-count column names, in the order they appear in the string.
    `raw`        DataFrame of raw counts indexed by guide name.

    Returns a Series over `raw.index` with the SAME dtype as `norm_col` -- string
    for multi-replicate columns, float for single-replicate ones -- so filling
    does not turn a numeric column into object. Intended for `.fillna()`;
    mageck's own values should be left exactly as they are.

    The single-replicate branch matters: the genome-wide `growing` comparisons
    have one replicate each, so those eight columns are float64 and have no
    slash to split on.
    """
    s = norm_col.reindex(raw.index)
    single = len(samples) == 1
    if single:
        parts = pd.DataFrame({0: pd.to_numeric(s, errors='coerce')})
    else:
        parts = s.astype('string').str.split('/', expand=True).astype(float)
        assert parts.shape[1] == len(samples), \
            f'{parts.shape[1]} replicates in the string, {len(samples)} samples given'

    fmt = ('{:.%dg}' % MAGECK_SIGFIGS).format
    out = None
    for i, name in enumerate(samples):
        f = size_factor(parts[i], raw[name])
        scaled = raw[name] * f
        if single:
            # keep the column numeric, but at mageck's precision
            return scaled.map(lambda v: float(fmt(v)))
        c = scaled.map(fmt)
        out = c if out is None else out + '/' + c
    return out
