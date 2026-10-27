"""Helper functions for the SILAC + TMT proteomics analyses (Figs 3B, 4A, 4E, 7C; S5B-F; S6A, S6C-E;
Tables S8, S9).

Loading and channel naming, total-signal normalization, replicate means with SEM, ratios and
differences with propagated error, two-sided t-tests, Benjamini-Hochberg correction, the
subsampled Mann-Whitney test used for every category comparison, localization categories, and
the violin / strip / volcano plots. Layout-agnostic: the TMT channel map and replicate groups are
passed in (see `silac_turnover.py` and `total_abundance.py`).

Sample-name convention: <bait>_<condition>_<timepoint>_<replicate>_<H|L>, e.g. RPS3_NT_3d_1_H.
"""
from __future__ import annotations

from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib.pyplot as plt
import seaborn as sns


# ---------------------------------------------------------------------------
# Loading, channel naming and normalization
# ---------------------------------------------------------------------------


def load_protein_quant(path: str, sep: str = "\t") -> pd.DataFrame:
    """Read a protein-quant TSV from the search engine output.

    Thin wrapper around `pd.read_csv` so the rest of the pipeline has a
    single, named entry point for raw input.
    """
    return pd.read_csv(path, sep=sep)


def drop_contaminants_and_decoys(
    df: pd.DataFrame,
    *,
    gene_col: str = "Gene Symbol",
    protein_id_col: str = "Protein Id",
    decoy_marker: str = "##",
    contaminant_marker: str = "_contaminant",
) -> pd.DataFrame:
    """Remove decoy hits (gene symbol contains `##`) and contaminants.

    Returns a new DataFrame with the index reset.
    """
    out = df.copy()
    out[gene_col] = out[gene_col].astype(str)
    decoy_mask = out[gene_col].str.contains(decoy_marker, regex=False, na=False)
    contam_mask = out[protein_id_col].astype(str).str.contains(
        contaminant_marker, regex=False, na=False
    )
    return out.loc[~(decoy_mask | contam_mask)].reset_index(drop=True)


def build_channel_map(
    raw_template: str,
    sample_template: str,
    channels: Sequence[str],
    sample_names: Sequence[str],
    *,
    labels: Sequence[str] = ("H", "L"),
    raw_label_tokens: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Build a `{raw_col: sample_name}` dict from templates.

    Example
    -------
    >>> channels = ['126', '127n', '127c', '128n']
    >>> sample_names = ['NT_3d_1', 'NT_3d_2', 'NT_3d_3', 'sgUFBP1_3d_1']
    >>> build_channel_map(
    ...     raw_template="{label_token}~rq_{channel}_sn sum",
    ...     sample_template="RPS3_{sample}_{label}",
    ...     channels=channels,
    ...     sample_names=sample_names,
    ...     labels=("H", "L"),
    ...     raw_label_tokens={"H": "heavy", "L": "light"},
    ... )

    The two templates use the placeholders:
      raw_template:    {label_token}, {channel}
      sample_template: {label}, {sample}
    """
    if len(channels) != len(sample_names):
        raise ValueError("channels and sample_names must be the same length")
    raw_label_tokens = raw_label_tokens or {label: label for label in labels}

    out: dict[str, str] = {}
    for label in labels:
        token = raw_label_tokens[label]
        for channel, sample in zip(channels, sample_names):
            raw = raw_template.format(label_token=token, channel=channel)
            new = sample_template.format(sample=sample, label=label)
            out[raw] = new
    return out


def rename_channels(df: pd.DataFrame, channel_map: Mapping[str, str]) -> pd.DataFrame:
    """Rename TMT channel columns to descriptive sample names.

    Parameters
    ----------
    df : DataFrame from `load_protein_quant`
    channel_map : dict
        Mapping from raw channel column name (e.g. ``'heavy~rq_126_sn sum'``)
        to sample name (e.g. ``'RPS3_NT_3d_1_H'``).

    Returns
    -------
    DataFrame with renamed columns (in place rename, returns same df).
    """
    return df.rename(columns=dict(channel_map))


def pair_silac_columns(
    sample_cols: Sequence[str],
    *,
    labels: Sequence[str] = ("H", "L"),
    sep: str = "_",
) -> dict[str, list[str]]:
    """Group H/L (or other label-suffixed) columns by their shared sample stem.

    Convenience for building the dict that `normalize_by_total` expects when
    you want to share a denominator across SILAC channels.

    Example
    -------
    >>> pair_silac_columns(['RPS3_NT_3d_1_H', 'RPS3_NT_3d_1_L',
    ...                     'RPS3_NT_3d_2_H', 'RPS3_NT_3d_2_L'])
    {'RPS3_NT_3d_1': ['RPS3_NT_3d_1_H', 'RPS3_NT_3d_1_L'],
     'RPS3_NT_3d_2': ['RPS3_NT_3d_2_H', 'RPS3_NT_3d_2_L']}

    Columns whose final token is not in `labels` are skipped.
    """
    label_set = set(labels)
    groups: dict[str, list[str]] = {}
    for col in sample_cols:
        parts = col.rsplit(sep, 1)
        if len(parts) != 2 or parts[1] not in label_set:
            continue
        stem = parts[0]
        groups.setdefault(stem, []).append(col)
    return groups


def normalize_by_total(
    df: pd.DataFrame,
    sample_cols: Sequence[str] | Mapping[str, Sequence[str]],
    *,
    suffix: str = "_norm",
    rescale: bool = True,
) -> pd.DataFrame:
    """Total-protein normalization, with optional shared denominators.

    Two ways to call this:

    1. **Independent columns** (`sample_cols` is a flat list).
       Each column is divided by its own column sum::

           normalize_by_total(df, ['s1_H', 's1_L', 's2_H', 's2_L'])

    2. **Grouped columns** (`sample_cols` is a dict).
       All columns in a group share a single denominator equal to the sum of
       *all* values across *all* columns in that group. This is the SILAC
       use case: H and L of the same sample must be divided by H+L combined,
       so their ratio is preserved::

           normalize_by_total(df, {
               'RPS3_NT_3d_1': ['RPS3_NT_3d_1_H', 'RPS3_NT_3d_1_L'],
               'RPS3_NT_3d_2': ['RPS3_NT_3d_2_H', 'RPS3_NT_3d_2_L'],
               ...
           })

       For SILAC data, **always use the dict form** with H and L paired —
       per-channel normalization (the flat-list version) discards the
       heavy/light ratio that turnover analysis depends on.

    If `rescale=True`, multiply through by the mean of the per-group sums so
    normalized values stay on roughly the input scale (useful before log
    transforms). Zero values are replaced with NaN.

    Parameters
    ----------
    df : DataFrame
    sample_cols : list of column names, or dict {group_label: [col1, col2, ...]}
    suffix : appended to each output column name (e.g. 's1_H' -> 's1_H_norm')
    rescale : if True, multiply by mean(group_sums) so values stay on input scale

    Returns
    -------
    DataFrame (same object) with normalized columns added.
    """
    if isinstance(sample_cols, Mapping):
        groups = {k: list(v) for k, v in sample_cols.items()}
    else:
        # Flat list: each column is its own group (independent normalization)
        groups = {col: [col] for col in sample_cols}

    # Per-group denominator = sum of *all* values across *all* columns in the group
    group_sums = {g: df[cols].to_numpy().sum() for g, cols in groups.items()}
    scale = (np.mean(list(group_sums.values())) if rescale else 1.0)

    for group, cols in groups.items():
        denom = group_sums[group]
        for col in cols:
            new_col = col + suffix
            if denom == 0 or not np.isfinite(denom):
                df[new_col] = np.nan
            else:
                df[new_col] = df[col] / denom * scale
                df.loc[df[new_col] == 0, new_col] = np.nan
    return df


def replicate_groups_from_names(
    sample_names: Sequence[str],
    *,
    sep: str = "_",
    replicate_field: int = -1,
    name_template: str = "{prefix}",
) -> dict[str, list[str]]:
    """Group sample names by stripping the replicate field.

    Example
    -------
    >>> replicate_groups_from_names(
    ...     ['RPS3_NT_3d_1_H', 'RPS3_NT_3d_2_H', 'RPS3_NT_3d_3_H',
    ...      'RPS3_NT_8d_1_H', 'RPS3_NT_8d_2_H', 'RPS3_NT_8d_3_H'])
    {'RPS3_NT_3d_H': ['RPS3_NT_3d_1_H', ...], 'RPS3_NT_8d_H': ['RPS3_NT_8d_1_H', ...]}

    `replicate_field` is the index (default -1, i.e. last token) of the
    field to strip. The remaining tokens are joined back with `sep`. Pass
    `replicate_field=-2` if your replicate is the second-to-last field
    (e.g. `..._1_H`, `..._2_H`).
    """
    groups: dict[str, list[str]] = {}
    for name in sample_names:
        tokens = name.split(sep)
        idx = replicate_field if replicate_field >= 0 else len(tokens) + replicate_field
        prefix = sep.join(tokens[:idx] + tokens[idx + 1:])
        key = name_template.format(prefix=prefix)
        groups.setdefault(key, []).append(name)
    return groups


def group_mean_sem(
    df: pd.DataFrame,
    replicate_groups: Mapping[str, Sequence[str]],
    *,
    mean_suffix: str = "_avg",
    sem_suffix: str = "_sem",
) -> pd.DataFrame:
    """Compute mean and SEM across replicates for each named group.

    Parameters
    ----------
    df : DataFrame
    replicate_groups : dict
        ``{'d3_NT_H': ['RPS3_NT_3d_1_H', 'RPS3_NT_3d_2_H', 'RPS3_NT_3d_3_H'], ...}``
        Keys are arbitrary group labels; values are the replicate columns.

    Returns
    -------
    DataFrame with `<group><mean_suffix>` and `<group><sem_suffix>` columns added.
    """
    for group, cols in replicate_groups.items():
        df[group + mean_suffix] = df[list(cols)].mean(axis=1)
        df[group + sem_suffix] = stats.sem(df[list(cols)].to_numpy(), axis=1, nan_policy="omit")
    return df


# ---------------------------------------------------------------------------
# Fold changes, error propagation and tests
# ---------------------------------------------------------------------------


def ratio_with_error(
    df: pd.DataFrame,
    numerator: str,
    denominator: str,
    *,
    out: str,
    num_sem: str | None = None,
    den_sem: str | None = None,
    log2: bool = False,
) -> pd.DataFrame:
    """Compute `numerator / denominator` with propagated SEM.

    Optionally also add the log2-transformed ratio and its propagated SEM.

    Standard relative-error propagation:
        SEM(num/den) = (num/den) * sqrt((SEM_num/num)^2 + (SEM_den/den)^2)
    Log2 propagation: SEM(log2(x)) = SEM(x) / x / ln(2)

    Parameters
    ----------
    df : DataFrame
    numerator, denominator : column names of the means
    out : base name for output columns. Creates `<out>` (the ratio) and,
        if SEM columns are provided, `<out>_sem`. If `log2=True` also creates
        `<out>_log2` and `<out>_log2_sem`.
    num_sem, den_sem : optional column names of the SEMs of the means
    log2 : also add log2 ratio and its propagated SEM
    """
    ratio = df[numerator] / df[denominator]
    ratio = ratio.replace([np.inf, -np.inf], np.nan)
    df[out] = ratio

    if num_sem is not None and den_sem is not None:
        rel_err = np.sqrt((df[num_sem] / df[numerator]) ** 2 + (df[den_sem] / df[denominator]) ** 2)
        df[out + "_sem"] = rel_err * ratio

    if log2:
        df[out + "_log2"] = np.log2(ratio)
        if num_sem is not None and den_sem is not None:
            # SEM(log2(r)) = SEM(r) / r / ln(2)
            df[out + "_log2_sem"] = df[out + "_sem"] / ratio / np.log(2)
    return df


def difference_with_error(
    df: pd.DataFrame,
    a: str,
    b: str,
    *,
    out: str,
    a_sem: str | None = None,
    b_sem: str | None = None,
) -> pd.DataFrame:
    """Compute `a - b` with propagated SEM (sqrt(SEM_a^2 + SEM_b^2))."""
    df[out] = df[a] - df[b]
    if a_sem is not None and b_sem is not None:
        df[out + "_sem"] = np.sqrt(df[a_sem] ** 2 + df[b_sem] ** 2)
    return df


def two_sided_t_test(
    df: pd.DataFrame,
    estimate: str,
    sem: str,
    *,
    df_resid: int,
    out_t: str | None = None,
    out_p: str | None = None,
) -> pd.DataFrame:
    """One-sample two-sided t-test that `estimate` differs from zero.

    Treats `estimate / sem` as a t-statistic with `df_resid` degrees of
    freedom. For a comparison of two groups of n replicates each (no
    pooled variance assumption beyond what's already in the SEMs),
    `df_resid = 2*n - 2` is the usual choice.

    Parameters
    ----------
    df : DataFrame
    estimate : column with the point estimate (e.g. log2 fold change)
    sem : column with its propagated SEM
    df_resid : degrees of freedom for the t-distribution
    out_t, out_p : output column names. Default: `<estimate>_t`, `<estimate>_p`.
    """
    out_t = out_t or estimate + "_t"
    out_p = out_p or estimate + "_p"
    t = df[estimate] / df[sem]
    df[out_t] = t
    df[out_p] = stats.t.sf(np.abs(t), df_resid) * 2
    return df


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted p-values.

    NaNs in `p_values` propagate to NaNs in the output (they are excluded
    from the rank/comparison count, then re-inserted at their original
    positions). Adjusted values are clipped to <=1.
    """
    p = np.asarray(p_values, dtype=float)
    out = np.full_like(p, np.nan)
    valid = ~np.isnan(p)
    pv = p[valid]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    # adjusted = min over k>=i of (ranked[k] * n / (k+1)), clipped to 1
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    # restore original order
    inverse = np.empty_like(order)
    inverse[order] = np.arange(n)
    out[valid] = adj[inverse]
    return out


def mann_whitney_subsampled(
    df: pd.DataFrame,
    *,
    category_col: str,
    value_col: str,
    categories: Sequence[str] | None = None,
    alternative: str = "two-sided",
    max_n: int = 50,
    n_resamples: int = 1000,
    summary: str = "median",
    random_state: int | None = 0,
    fdr: bool = True,
    exclude_categories: Sequence[str] = ("",),
) -> pd.DataFrame:
    """Mann-Whitney U with both groups subsampled to equal n.

    Why this exists: the regular `mann_whitney_by_category` returns vanishingly
    small p-values for tiny shifts in large categories, because rank-sum tests
    are "consistent" -- power -> 1 with n for any non-zero shift. Subsampling
    both groups to a common, capped n removes the n-dominance: the test asks
    "if I had only `max_n` proteins per group, would the shift still look
    significant?" Repeating the draw and reporting the median p-value gives a
    stable summary.

    Suggested figure-legend wording (one liner):
        "Mann-Whitney U on 1,000 random subsamples of n=50 per group; median p
         reported. *p<0.05, **p<0.01, ***p<0.001."

    Parameters
    ----------
    df : DataFrame
    category_col, value_col : column names
    categories : optional list to restrict / order
    alternative : ``'two-sided'``, ``'less'``, ``'greater'``
    max_n : per-group sample size for each draw. Each draw uses
        ``n_draw = min(n_in, n_out, max_n)``. Use 30-50 for large studies.
    n_resamples : number of subsampled MWU tests per category. The
        reported p-value is the median across draws (`summary='median'`)
        or the mean (`summary='mean'`).
    summary : ``'median'`` (default) or ``'mean'`` for combining the
        per-draw p-values.
    random_state : seed.
    fdr : add Benjamini-Hochberg adjusted p-values.
    exclude_categories : never reported as a row (still kept in the rest pool).

    Returns
    -------
    DataFrame with one row per category, columns:
        category, n_in, n_out, n_draw, median_in, median_out, median_diff,
        p_value [, p_adj]

    Notes
    -----
    `n_draw` is the per-group size actually used for each subsample (which
    may be smaller than `max_n` for sparse categories). The `n_in` / `n_out`
    columns report the full available counts for context.
    """
    if alternative not in {"two-sided", "less", "greater"}:
        raise ValueError("alternative must be 'two-sided', 'less', or 'greater'")
    if summary not in {"median", "mean"}:
        raise ValueError("summary must be 'median' or 'mean'")

    rng = np.random.default_rng(random_state)

    sub = df[[category_col, value_col]].dropna(subset=[value_col]).copy()
    sub[category_col] = sub[category_col].fillna("")

    if categories is None:
        categories = [
            c for c in sub[category_col].unique()
            if c not in set(exclude_categories)
        ]

    rows = []
    for cat in categories:
        in_v = sub.loc[sub[category_col] == cat, value_col].to_numpy()
        out_v = sub.loc[sub[category_col] != cat, value_col].to_numpy()
        n_in, n_out = in_v.size, out_v.size
        n_draw = int(min(n_in, n_out, max_n))

        if n_draw < 2:
            rows.append({
                "category": cat, "n_in": n_in, "n_out": n_out, "n_draw": n_draw,
                "median_in": (np.median(in_v) if n_in else np.nan),
                "median_out": (np.median(out_v) if n_out else np.nan),
                "median_diff": np.nan,
                "p_value": np.nan,
            })
            continue

        ps = np.empty(n_resamples, dtype=float)
        for k in range(n_resamples):
            in_b = rng.choice(in_v, size=n_draw, replace=False)
            out_b = rng.choice(out_v, size=n_draw, replace=False)
            ps[k] = stats.mannwhitneyu(in_b, out_b, alternative=alternative).pvalue

        agg = float(np.median(ps)) if summary == "median" else float(np.mean(ps))

        rows.append({
            "category": cat, "n_in": n_in, "n_out": n_out, "n_draw": n_draw,
            "median_in": float(np.median(in_v)),
            "median_out": float(np.median(out_v)),
            "median_diff": float(np.median(in_v) - np.median(out_v)),
            "p_value": agg,
        })

    result = pd.DataFrame(rows)
    if fdr and len(result):
        result["p_adj"] = benjamini_hochberg(result["p_value"].to_numpy())
    return result.sort_values("p_value", na_position="last").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Localization categories
# ---------------------------------------------------------------------------


def annotate_localization(
    df: pd.DataFrame,
    localization_df: pd.DataFrame,
    *,
    gene_col: str = "Gene Symbol",
    loc_gene_col: str = "Gene_name_canonical",
    loc_label_col: str = "graph_localization_annotation",
    out_col: str = "category",
    default: str = "",
) -> pd.DataFrame:
    """Map each protein in `df` to a localization label from `localization_df`.

    The localization table is typically a published OpenCell-style annotation
    (e.g. supplementary table with `Gene_name_canonical` and a localization
    column). Genes not present in the localization table get `default`.

    Returns
    -------
    DataFrame (same object) with `out_col` added/overwritten.
    """
    df[out_col] = default
    for label in localization_df[loc_label_col].dropna().unique():
        genes = localization_df.loc[
            localization_df[loc_label_col] == label, loc_gene_col
        ]
        df.loc[df[gene_col].isin(genes), out_col] = label
    return df


def assign_category_overrides(
    df: pd.DataFrame,
    overrides: Mapping[str, Iterable[str]],
    *,
    gene_col: str = "Gene Symbol",
    out_col: str = "category",
) -> pd.DataFrame:
    """Manually assign categories to specific gene lists (e.g. RPL, RPS).

    Parameters
    ----------
    overrides : dict ``{category_label: [gene1, gene2, ...]}``

    Later entries overwrite earlier ones; useful for forcing ribosomal
    proteins into a single bin even if they're scattered in the localization
    annotation.
    """
    for label, genes in overrides.items():
        df.loc[df[gene_col].isin(list(genes)), out_col] = label
    return df


def annotate_by_description(
    df: pd.DataFrame,
    rules: Mapping[str, str],
    *,
    description_col: str = "Description",
    out_col: str = "category",
    case: bool = False,
) -> pd.DataFrame:
    """Assign categories based on substring matches in a description column.

    Parameters
    ----------
    rules : dict ``{category_label: substring}``
        Each entry sets `out_col = label` for rows whose description contains
        `substring`. Later rules overwrite earlier ones.
    case : passed to ``str.contains`` (default case-insensitive).
    """
    desc = df[description_col].astype(str)
    for label, needle in rules.items():
        df.loc[desc.str.contains(needle, case=case, na=False), out_col] = label
    return df


def collapse_categories(
    df: pd.DataFrame,
    keep: Iterable[str],
    *,
    source_col: str = "category",
    out_col: str = "category_simplified",
    other_label: str = "",
) -> pd.DataFrame:
    """Make a simplified category column that keeps only listed labels.

    Anything not in `keep` becomes `other_label`. Useful for plotting when the
    full localization annotation has too many categories.
    """
    keep = set(keep)
    df[out_col] = df[source_col].where(df[source_col].isin(keep), other_label)
    return df


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def _despine_and_rotate(ax: plt.Axes, rotation: int = 90) -> None:
    plt.setp(ax.get_xticklabels(), rotation=rotation, ha="right")
    sns.despine(ax=ax)


def category_violin(
    df: pd.DataFrame,
    *,
    category_col: str,
    value_col: str | None = None,
    value_cols: Sequence[str] | None = None,
    value_labels: Sequence[str] | None = None,
    categories: Sequence[str] | None = None,
    title: str = "",
    ylabel: str | None = None,
    ylim: tuple[float, float] | None = None,
    median_line: bool = True,
    color: str = "lightgrey",
    split_palette: Sequence[str] = ("#4C72B0", "#DD8452"),
    figsize: tuple[float, float] = (15, 10),
    save: str | None = None,
    ax: plt.Axes | None = None,
    **violin_kwargs,
) -> plt.Axes:
    """Violin plot of one or two value columns grouped by `category_col`.

    Two modes:

    1. **Single** (pass `value_col`): one violin per category, filled grey,
       with an internal box plot showing quartiles + median + whiskers.
    2. **Split** (pass `value_cols=[col_a, col_b]`): a two-color split violin
       per category — left half is `col_a`, right half is `col_b`. Internal
       quartile lines are drawn instead of a box (seaborn's box-inside
       doesn't render on split halves). Use `value_labels=('NT','sgUFBP1')`
       to control the legend.

    Parameters
    ----------
    df : DataFrame
    category_col : column to group by (x-axis)
    value_col : single y-axis column (single mode)
    value_cols : two y-axis columns (split mode); mutually exclusive with `value_col`
    value_labels : optional ``(label_a, label_b)`` for the legend in split mode.
        Defaults to the column names themselves.
    categories : optional list of categories to keep AND order on the x-axis.
        If `None`, all categories present are shown in alphabetical order.
    median_line : if True, draws a dashed horizontal line at the overall
        median of the value(s) being plotted.
    color : fill color for single-mode violins (default light grey)
    split_palette : two colors for the split halves (default seaborn-ish blue/orange)
    **violin_kwargs : forwarded to `sns.violinplot` (e.g. `cut=0`, `bw_adjust=0.7`).

    Returns
    -------
    The matplotlib Axes.
    """
    # --- validate the mode ---
    if (value_col is None) == (value_cols is None):
        raise ValueError("Pass exactly one of `value_col` or `value_cols`.")
    split_mode = value_cols is not None
    if split_mode and len(value_cols) != 2:
        raise ValueError("`value_cols` must contain exactly two column names.")

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    # --- subset / order categories ---
    if categories is not None:
        plot_df = df[df[category_col].isin(list(categories))].copy()
        order = list(categories)
    else:
        plot_df = df.copy()
        order = sorted(plot_df[category_col].dropna().unique().tolist())

    # --- single mode ---
    if not split_mode:
        sns.violinplot(
            data=plot_df,
            x=category_col,
            y=value_col,
            order=order,
            color=color,
            inner="box",
            ax=ax,
            **violin_kwargs,
        )
        if median_line:
            ax.axhline(df[value_col].median(), linestyle="--", color="black", linewidth=1)
        ax.set_ylabel(ylabel or value_col)

    # --- split mode ---
    else:
        col_a, col_b = value_cols
        labels = list(value_labels) if value_labels is not None else [col_a, col_b]
        if len(labels) != 2:
            raise ValueError("`value_labels` must have length 2.")

        # Melt to long format: one row per (category, hue, value)
        long_df = pd.melt(
            plot_df,
            id_vars=[category_col],
            value_vars=[col_a, col_b],
            var_name="_series",
            value_name="_value",
        )
        long_df["_series"] = long_df["_series"].map({col_a: labels[0], col_b: labels[1]})

        sns.violinplot(
            data=long_df,
            x=category_col,
            y="_value",
            hue="_series",
            hue_order=labels,
            order=order,
            split=True,
            inner="quart",
            palette=list(split_palette),
            ax=ax,
            **violin_kwargs,
        )
        if median_line:
            both = pd.concat([df[col_a], df[col_b]])
            ax.axhline(both.median(), linestyle="--", color="black", linewidth=1)
        ax.set_ylabel(ylabel or f"{col_a} / {col_b}")
        ax.legend(title=None, frameon=False)

    ax.set_title(title)
    if ylim:
        ax.set_ylim(*ylim)
    _despine_and_rotate(ax)

    if save:
        ax.figure.savefig(save, bbox_inches="tight")
    return ax


def gene_set_long(
    df: pd.DataFrame,
    *,
    gene_sets: Mapping[str, Sequence[str]],
    value_cols: Sequence[str],
    value_labels: Sequence[str] | None = None,
    background: Mapping[str, object] | None = None,
    extra_cols: Sequence[str] | None = None,
    gene_col: str = "Gene",
    category_col: str = "category_simplified",
    group_col: str = "gene_set",
    series_col: str = "series",
    value_name: str = "value",
) -> pd.DataFrame:
    """Tidy one row per (protein, group, series) for gene-set strip plots.

    `gene_sets` maps a display label -> list of gene symbols. `background`
    maps a display label -> one of:
      * ``None`` / ``'ALL'``  -> every row in `df` (proteome-wide)
      * a string or list      -> rows whose `category_col` matches
      * a boolean Series      -> used directly as a mask
    A protein can appear in several groups (e.g. SEC61A1 is both in 'Sec61'
    and in the 'ER' background); each membership is its own row.

    `extra_cols` are per-protein columns (e.g. an abundance measure) carried
    through unmelted, so they can be used to style the points.

    The returned frame carries a boolean `is_background` column so plotting
    code can draw backgrounds as violins and gene sets as points.
    """
    labels = list(value_labels) if value_labels is not None else list(value_cols)
    if len(labels) != len(value_cols):
        raise ValueError("`value_labels` must be the same length as `value_cols`")
    id_cols = [gene_col] + [c for c in (extra_cols or []) if c != gene_col]

    pieces: list[pd.DataFrame] = []

    def _add(label: str, sub: pd.DataFrame, is_background: bool) -> None:
        if sub.empty:
            return
        long = pd.melt(
            sub[id_cols + list(value_cols)],
            id_vars=id_cols,
            value_vars=list(value_cols),
            var_name=series_col,
            value_name=value_name,
        )
        long[series_col] = long[series_col].map(dict(zip(value_cols, labels)))
        long[group_col] = label
        long["is_background"] = is_background
        pieces.append(long.dropna(subset=[value_name]))

    for label, spec in (background or {}).items():
        if spec is None or (isinstance(spec, str) and spec == "ALL"):
            mask = pd.Series(True, index=df.index)
        elif isinstance(spec, pd.Series):
            mask = spec.reindex(df.index).fillna(False).astype(bool)
        else:
            wanted = [spec] if isinstance(spec, str) else list(spec)
            mask = df[category_col].isin(wanted)
        _add(label, df.loc[mask], True)

    for label, genes in gene_sets.items():
        _add(label, df.loc[df[gene_col].isin(list(genes))], False)

    out = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(
        columns=id_cols + [series_col, value_name, group_col, "is_background"]
    )
    order = list(background or {}) + list(gene_sets)
    out[group_col] = pd.Categorical(out[group_col], categories=order, ordered=True)
    out[series_col] = pd.Categorical(out[series_col], categories=labels, ordered=True)
    return out.sort_values([group_col, series_col, gene_col]).reset_index(drop=True)


def gene_set_strip(
    df: pd.DataFrame,
    *,
    gene_sets: Mapping[str, Sequence[str]],
    value_col: str | None = None,
    value_cols: Sequence[str] | None = None,
    value_labels: Sequence[str] | None = None,
    background: Mapping[str, object] | None = None,
    gene_col: str = "Gene",
    category_col: str = "category_simplified",
    title: str = "",
    ylabel: str | None = None,
    ylim: tuple[float, float] | None = None,
    reference_line: float | None = 0.0,
    median_line: bool = True,
    point_color: str = "#E8871A",
    point_colors: Mapping[str, str] | None = None,
    split_palette: Sequence[str] = ("#B3B3B3", "#E8871A"),
    background_color: str = "#B3B3B3",
    point_size: float = 26.0,
    size_col: str | None = None,
    size_range: tuple[float, float] = (14.0, 90.0),
    size_log: bool = True,
    size_clip: tuple[float, float] = (0.02, 0.98),
    size_label: str | None = None,
    size_legend: bool = True,
    size_legend_values: Sequence[float] | None = None,
    label_genes: bool = True,
    label_mode: str = "all",
    outlier_k: float = 1.5,
    label_all_below: int = 4,
    label_max_n: int = 8,
    mean_bar: bool = False,
    connect: bool = True,
    pair_offset: float = 0.07,
    separator_after: Sequence[str] | None = None,
    figsize: tuple[float, float] = (7.5, 3.2),
    save: str | None = None,
    ax: plt.Axes | None = None,
    seed: int = 0,
) -> tuple[plt.Axes, pd.DataFrame]:
    """Strip plot of named gene sets against violins of background populations.

    Two modes, mirroring `category_violin`:

    1. **Single** (`value_col`): one dot per protein per gene set.
    2. **Split** (`value_cols=[col_a, col_b]`): the two columns are drawn at
       offset x positions within each group (e.g. NT vs sgDDRGK1), with a
       split violin for the background groups. Sets with
       ``n <= label_max_n`` are drawn as one tight pair per protein
       (+/- `pair_offset`) spread across the slot and joined by a thin line
       (`connect`); larger sets are dodged jittered clouds that line up with
       the violin halves.

    `separator_after` takes group labels after which a faint vertical rule is
    drawn, to set the background block apart from the gene sets.

    `label_mode` is ``'all'`` (label every gene in a set small enough to take
    labels) or ``'outliers'`` (label only genes outside the Tukey fence
    ``Q1 - outlier_k*IQR .. Q3 + outlier_k*IQR`` of their own set, which also
    works for large sets; sets with fewer than `label_all_below` members are
    labelled in full, since a fence is meaningless there). `mean_bar` draws a
    short horizontal bar at each set's mean.

    `size_col` scales the marker area by a per-protein measure such as total
    abundance. Marker *area* is linear in log10(`size_col`) by default
    (`size_log`), because abundance spans decades; the mapping is anchored on
    the `size_clip` quantiles of the whole `df`, so sizes are comparable
    between figures drawn from the same table. A size legend is drawn
    whenever `size_col` is set - label it with `size_label`.

    Background groups (`background`, e.g.
    ``{'All proteins': None, 'Cytosol': 'cytosol', 'ER': 'ER'}``) are drawn
    as violins at the left of the axes; gene sets are drawn as points, in the
    order given. Points are spread deterministically for sets with
    ``n <= label_max_n`` (so gene labels line up) and jittered for larger sets.

    Returns
    -------
    ``(ax, plot_df)`` - the Axes and the tidy frame that was plotted, so the
    exact values behind each dot can be written out alongside the figure.
    """
    if (value_col is None) == (value_cols is None):
        raise ValueError("Pass exactly one of `value_col` or `value_cols`.")
    cols = [value_col] if value_col is not None else list(value_cols)
    if value_cols is not None and len(cols) != 2:
        raise ValueError("`value_cols` must contain exactly two column names.")
    split_mode = value_cols is not None

    plot_df = gene_set_long(
        df,
        gene_sets=gene_sets,
        value_cols=cols,
        value_labels=value_labels,
        background=background,
        extra_cols=[size_col] if size_col else None,
        gene_col=gene_col,
        category_col=category_col,
    )
    labels = list(plot_df["series"].cat.categories)
    order = list(plot_df["gene_set"].cat.categories)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    bg = plot_df[plot_df["is_background"]]
    sets = plot_df[~plot_df["is_background"]]

    # --- background populations as violins ---
    if not bg.empty:
        violin_kws = dict(
            data=bg, x="gene_set", y="value", order=order,
            cut=0, linewidth=0.6, width=0.8, ax=ax,
        )
        if split_mode:
            sns.violinplot(hue="series", hue_order=labels, split=True, inner="quart",
                           palette=list(split_palette), **violin_kws)
        else:
            sns.violinplot(color=background_color, inner="box", **violin_kws)

    # --- gene sets as points ---
    # marker area <- log10(abundance), anchored on the whole table so that
    # the same abundance gives the same dot size in every panel
    if size_col is not None:
        anchor = pd.to_numeric(df[size_col], errors="coerce")
        if size_log:
            anchor = np.log10(anchor.where(anchor > 0))
        size_lo, size_hi = anchor.quantile(list(size_clip))
        if not np.isfinite([size_lo, size_hi]).all() or size_hi <= size_lo:
            size_lo, size_hi = np.nanmin(anchor), np.nanmax(anchor)

    def _sizes(values: pd.Series) -> np.ndarray:
        if size_col is None:
            return np.full(len(values), point_size)
        v = pd.to_numeric(values, errors="coerce")
        if size_log:
            v = np.log10(v.where(v > 0))
        frac = np.clip((v - size_lo) / (size_hi - size_lo), 0.0, 1.0)
        out = size_range[0] + frac * (size_range[1] - size_range[0])
        return out.fillna(point_size).to_numpy()

    rng = np.random.default_rng(seed)
    series_colors = (
        dict(zip(labels, list(split_palette))) if split_mode
        else {labels[0]: point_color}
    )

    for name, group in sets.groupby("gene_set", observed=True, sort=False):
        i = order.index(name)
        small = group[gene_col].nunique() <= label_max_n
        color_override = (point_colors or {}).get(name)

        if small:
            # Order small sets by value so x tracks y: the dots (and their
            # labels) ladder upwards instead of piling into one column. In
            # split mode each gene becomes a tight pair centred on its own
            # sub-slot, which keeps the join lines short and untangled.
            genes = (group.groupby(gene_col, observed=True)["value"].mean()
                     .sort_values().index.tolist())
            slot = 0.72
            centres = (np.linspace(-0.5, 0.5, len(genes)) * slot
                       if len(genes) > 1 else np.zeros(1))
            gene_centre = dict(zip(genes, i + centres))
            # shrink the pair offset if the sub-slots are close together, so a
            # dense set's pairs stay inside their own sub-slot
            gap = slot / (len(genes) - 1) if len(genes) > 1 else slot
            offset = min(pair_offset, 0.38 * gap)
            offsets = ({labels[0]: -offset, labels[1]: offset}
                       if split_mode else {labels[0]: 0.0})
        else:
            # Large sets stay as jittered clouds dodged to match the split
            # violins; a shared per-gene jitter keeps the two series aligned.
            genes = sorted(group[gene_col].unique())
            slot = 0.3 if split_mode else 0.62
            gene_centre = dict(zip(genes,
                                   i + rng.uniform(-slot / 2, slot / 2, len(genes))))
            offsets = ({labels[0]: -0.2, labels[1]: 0.2}
                       if split_mode else {labels[0]: 0.0})

        def _x(sub: pd.DataFrame, series: str) -> np.ndarray:
            return sub[gene_col].map(gene_centre).to_numpy() + offsets[series]

        if split_mode and connect and small:
            wide = group.pivot_table(index=gene_col, columns="series",
                                     values="value", observed=True)
            for gene, row in wide.iterrows():
                if row.reindex(labels).isna().any():
                    continue
                ax.plot(
                    [gene_centre[gene] + offsets[labels[0]],
                     gene_centre[gene] + offsets[labels[1]]],
                    [row[labels[0]], row[labels[1]]],
                    color="dimgrey", linewidth=0.5, zorder=2,
                )

        for series in labels:
            sub = group[group["series"] == series]
            if sub.empty:
                continue
            ax.scatter(
                _x(sub, series), sub["value"].to_numpy(),
                s=_sizes(sub[size_col]) if size_col else point_size,
                facecolor=color_override or series_colors[series],
                edgecolor="white", linewidth=0.5,
                zorder=3, clip_on=True,
            )

        if mean_bar:
            # short bar at the set mean, one per series
            for series in labels:
                sub = group[group["series"] == series]
                if sub.empty:
                    continue
                half = (0.16 if split_mode else 0.34)
                centre = i + offsets[series]
                ax.plot([centre - half, centre + half],
                        [sub["value"].mean()] * 2,
                        color="#404040", linewidth=1.4, solid_capstyle="butt",
                        zorder=4)

        if label_genes:
            if label_mode == "outliers":
                # Tukey fence on this set, so only genes that stand out get named
                q1, q3 = group["value"].quantile([0.25, 0.75])
                iqr = q3 - q1
                flagged = group.loc[
                    (group["value"] < q1 - outlier_k * iqr)
                    | (group["value"] > q3 + outlier_k * iqr), gene_col
                ].unique()
                to_label = [g for g in genes if g in set(flagged)]
                # a 2- or 3-member set has no outliers to speak of, and no
                # clutter either, so name all of them
                if len(genes) < label_all_below:
                    to_label = list(genes)
            elif small:
                to_label = genes
            else:
                to_label = []

            # left-to-right, so the alternating label rows separate neighbours
            to_label = sorted(to_label, key=lambda g: gene_centre[g])
            # labels alternate over two to four rows so crowded sets stay readable
            rows_used = 2 if len(to_label) <= 5 else 3 if len(to_label) <= 7 else 4
            set_median = group["value"].median()
            for rank, gene in enumerate(to_label):
                vals = group.loc[group[gene_col] == gene, "value"]
                # low outliers get labelled downwards, so the text does not land
                # on top of the rest of the set
                below = vals.max() < set_median
                pad = 5 + 13 * (rank % rows_used)
                # fan the labels out sideways too - in a jittered set several
                # outliers can sit at almost the same x, and the leader line
                # keeps each one tied to its own dot
                dx = ((rank - (len(to_label) - 1) / 2) * 6.5
                      if label_mode == "outliers" and len(to_label) > 1 else 0.0)
                ax.annotate(
                    gene,
                    xy=(gene_centre[gene], vals.min() if below else vals.max()),
                    xytext=(dx, -pad if below else pad),
                    textcoords="offset points",
                    rotation=90, ha="center", va="top" if below else "bottom",
                    fontsize=5, color="#404040", annotation_clip=True,
                    arrowprops=dict(arrowstyle="-", color="#C0C0C0",
                                    linewidth=0.4, shrinkA=0, shrinkB=2),
                )

    # --- reference lines ---
    if reference_line is not None:
        ax.axhline(reference_line, linestyle="-", color="#909090", linewidth=0.6,
                   zorder=1)
    if median_line:
        for series, color in zip(labels, list(split_palette) if split_mode
                                 else [background_color]):
            med = df[dict(zip(labels, cols))[series]].median()
            ax.axhline(med, linestyle="--", color=color, linewidth=1, zorder=1)

    # --- axes cosmetics ---
    if ylim is None:
        lo_bg, hi_bg = (bg["value"].quantile([0.005, 0.995]) if not bg.empty
                        else (plot_df["value"].min(), plot_df["value"].max()))
        lo = min(lo_bg, sets["value"].min()) if not sets.empty else lo_bg
        hi = max(hi_bg, sets["value"].max()) if not sets.empty else hi_bg
        pad = 0.18 * (hi - lo)
        ylim = (lo - 0.05 * (hi - lo), hi + pad)   # headroom for gene labels
    ax.set_ylim(*ylim)
    ax.set_xlim(-0.6, len(order) - 0.4)
    for label in separator_after or []:
        ax.axvline(order.index(label) + 0.5, color="#DDDDDD", linewidth=0.6,
                   zorder=0)
    # a gene set with nothing quantified gets flagged rather than left blank
    plotted = set(sets["gene_set"].unique())
    for label in gene_sets:
        if label not in plotted:
            ax.text(order.index(label), ylim[0] + 0.04 * (ylim[1] - ylim[0]),
                    "n.d.", ha="center", va="bottom", fontsize=5,
                    color="#909090")
    ax.set_xlabel("")
    ax.set_ylabel(ylabel or (cols[0] if not split_mode else " / ".join(cols)))
    ax.set_title(title)

    hue_legend = None
    if split_mode:
        hue_legend = ax.legend(title=None, frameon=False, loc="upper left",
                               fontsize=7)
    if size_col is not None and size_legend:
        ticks = size_legend_values
        if ticks is None:
            # whole decades inside the mapped range, else its two endpoints
            if size_log:
                # tolerance so a range of e.g. -0.99..1.02 still yields 0.1/1/10
                decades = np.arange(np.ceil(size_lo - 0.05),
                                    np.floor(size_hi + 0.05) + 1)
                ticks = (10.0 ** decades if len(decades) >= 2
                         else 10.0 ** np.array([size_lo, size_hi]))
            else:
                ticks = np.linspace(size_lo, size_hi, 3)
        handles = [
            ax.scatter([], [], s=_sizes(pd.Series([v]))[0], facecolor="none",
                       edgecolor="#606060", linewidth=0.6)
            for v in ticks
        ]
        legend_of_sizes = ax.legend(
            handles, [f"{v:g}" for v in ticks],
            title=size_label or size_col, loc="upper left",
            bbox_to_anchor=(1.01, 1.0),          # outside, clear of the data
            frameon=False, fontsize=6, title_fontsize=6,
            labelspacing=1.2, handletextpad=0.9, borderpad=0.4,
        )
        if hue_legend is not None:
            ax.add_artist(hue_legend)
    _despine_and_rotate(ax, rotation=45)

    if save:
        ax.figure.savefig(save, dpi=300, bbox_inches="tight")
    return ax, plot_df


def volcano(
    df: pd.DataFrame,
    *,
    effect_col: str,
    pvalue_col: str,
    highlight: Mapping[str, Sequence[str]] | None = None,
    category_col: str = "category",
    title: str = "",
    xlabel: str | None = None,
    xlim: tuple[float, float] | None = None,
    base_color: str = "lightgrey",
    base_alpha: float = 0.3,
    figsize: tuple[float, float] = (8, 7),
    save: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Volcano plot: effect size vs -log10(p).

    `highlight` lets you overlay specific categories in distinct colors:
        ``{'RPL': 'red', 'RPS': 'orange'}``  -- value can be color string OR
        ``{'RPL': ['RPL10A', 'RPL3', ...], 'RPS': [...]}``  -- list of genes
    Only the dict-of-color form uses `category_col`; the dict-of-list form
    matches gene symbols directly via the index of the DataFrame's `Gene
    Symbol` column (assumed present).
    """
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    x = df[effect_col]
    y = -np.log10(df[pvalue_col])
    ax.scatter(x, y, alpha=base_alpha, color=base_color, s=12, label="all", rasterized=True)

    if highlight:
        for label, spec in highlight.items():
            if isinstance(spec, str):
                # spec is a color, filter on category_col
                mask = df[category_col] == label
                color = spec
            else:
                # spec is a list of gene symbols
                if "Gene Symbol" not in df.columns:
                    raise KeyError("highlight by gene list requires a 'Gene Symbol' column")
                mask = df["Gene Symbol"].isin(list(spec))
                color = None
            ax.scatter(x[mask], y[mask], alpha=0.7, s=18, label=label, color=color, rasterized=True)

    ax.set_xlabel(xlabel or effect_col)
    ax.set_ylabel(f"-log10({pvalue_col})")
    ax.set_title(title)
    if xlim:
        ax.set_xlim(*xlim)
    if highlight:
        ax.legend(frameon=False, loc="best")
    sns.despine(ax=ax)
    if save:
        ax.figure.savefig(save, bbox_inches="tight")
    return ax
