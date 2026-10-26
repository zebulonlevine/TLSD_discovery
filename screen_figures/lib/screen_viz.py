#!/usr/bin/env python3
"""
CRISPRi Screen Visualization Toolkit
====================================

Tools for visualizing and comparing CRISPRi screen results.

Features:
    - Replicate correlation scatterplots
    - Clustered correlation heatmaps (counts or fold-changes)
    - Cross-comparison scatter plots with hit status sizing
    - Gene category coloring
    - Customizable labeling

Usage:
    # Replicate correlation scatter
    python screen_viz.py scatter-replicates --counts counts.csv --x Sample1 --y Sample2 -o plot.pdf

    # Correlation heatmap
    python screen_viz.py heatmap --counts counts.csv --samples S1,S2,S3,S4 -o heatmap.pdf

    # Cross-comparison scatter
    python screen_viz.py compare --results comp1_all_genes.csv comp2_all_genes.csv \
        --hits comp1_hits.csv comp2_hits.csv --categories categories.csv -o compare.pdf
"""

import argparse
import re
import sys
import textwrap
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass, field

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram


# =============================================================================
# Figure output settings (must be set before any figure is created)
# =============================================================================

# matplotlib defaults to pdf.fonttype=3 (Type 3), whose glyphs import into
# Illustrator as un-editable outline art. 42 = TrueType keeps text as live,
# re-typable text. svg.fonttype='none' does the equivalent for SVG.
# These only take effect on re-run, so existing PDFs are unchanged until a
# figure is deliberately regenerated.
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Nimbus Sans', 'DejaVu Sans']
# Keep mathtext in Arial too, otherwise $...$ silently falls back to DejaVu.
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'


# =============================================================================
# Negative-control (NTC) identification
# =============================================================================

NTC_LABELS = frozenset({
    'negative_control', 'negative control', 'non-targeting', 'nontargeting',
    'ntc', 'neg_control', 'negative', 'safe_harbor',
})


def is_ntc(names: pd.Series) -> pd.Series:
    """Exact-label NTC detection on the leading token of a gene id.

    Substring matching false-flags real genes -- KNTC1 contains "NTC" and
    NEGR1 contains "NEG" -- which would draw them as light-grey control points
    and fold them into the NTC reference violin. See CLAUDE.md section 3.

    Handles both levels this module plots at, because both are 'SYMBOL_suffix':
    gene-level ids from *_all_genes.csv ('KNTC1_P1P2', 'NTC_412') and
    sgRNA-level gene labels ('non-targeting_7_negative_control').

    Mirrors screen_pipeline.is_ntc(); the two modules are independently
    importable, so the definition is intentionally duplicated rather than
    creating an import dependency between them.
    """
    head = names.astype(str).str.strip().str.split('_').str[0].str.lower()
    return head.isin(NTC_LABELS)



# Display-only gene labels. Mirrors Claude_ScreenAnalysis_v4/gene_aliases.py, but
# deliberately holds a SMALLER set: only genes whose plotted point would otherwise
# be misread. ANKZF1's sublibrary guides overlap the ATG9A locus, and it is
# categorised as Lysosome/Autophagy for that reason, so the point must carry both
# symbols. Applied at draw time only -- never to category keys, join keys or
# hit-membership, so nothing about which point is plotted or coloured changes.
DISPLAY_ALIASES = {
    'ANKZF1': 'ATG9A/ANKZF1',
}


def display_label(g):
    """'ANKZF1' -> 'ATG9A/ANKZF1'; 'ANKZF1_P1P2' -> 'ATG9A/ANKZF1_P1P2'."""
    head, sep, tail = str(g).partition('_')
    return DISPLAY_ALIASES.get(head, head) + sep + tail


# ---------------------------------------------------------------------------
# Category membership
#
# The category CSVs legitimately list some genes under MORE THAN ONE category
# (ATG16L1 is both Autophagy and E3; RACK1 is both Biogenesis and
# Lysosome/Autophagy). The two plot families need different resolutions:
#
#   violins  -- a gene can contribute to every category it belongs to, since a
#               distribution can legitimately count it twice. _gene_to_cats_map.
#   scatters -- a point carries exactly one colour, so the conflict must be
#               resolved. Autophagy wins; otherwise the last row wins, which is
#               the dict(zip(...)) behaviour every figure used previously.
# ---------------------------------------------------------------------------
PREFER_CATEGORY = ('autophagy',)


def _gene_to_cat_map(categories, prefer=PREFER_CATEGORY):
    """gene -> ONE category, resolving multi-category genes in favour of `prefer`."""
    out = {}
    for gene, cat in zip(categories.iloc[:, 0].astype(str),
                         categories.iloc[:, 1].astype(str)):
        held = out.get(gene)
        if held is None or not any(p in held.lower() for p in prefer):
            out[gene] = cat
    return out


def _gene_to_cats_map(categories):
    """gene -> LIST of every category it is listed under, in file order."""
    out = {}
    for gene, cat in zip(categories.iloc[:, 0].astype(str),
                         categories.iloc[:, 1].astype(str)):
        lst = out.setdefault(gene, [])
        if cat not in lst:
            lst.append(cat)
    return out

# =============================================================================
# Color Palettes
# =============================================================================

# COLORS = {
#     'orange': '#E69F00',
#     'sky_blue': '#56B4E9',
#     'green': '#009E73',
#     'yellow': '#F0E442',
#     'blue': '#0072B2',
#     'vermillion': '#D55E00',
#     'purple': '#CC79A7',
#     'black': '#000000',
#     'grey': '#999999',
#     'light_grey': '#CCCCCC',
#     'dark_grey': '#4D4D4D',
#     'red': '#D55E00',
#     'cyan': '#56B4E9',
#     'pink': '#CC79A7',
# }

COLORS = {
    'orange': '#E69F00',
    'sky_blue': '#56B4E9',
    'green': '#009E73',
    'yellow': '#F0E442',
    'blue': '#0072B2',
    'vermillion': '#D55E00',
    'purple': '#CC79A7',
    'black': '#000000',
    'grey': '#676566',#was '#4D4D4D'
    'light_grey': '#999999',
    'dark_grey': '#676566', #was '#4D4D4D'
    'red': '#D55E00',
    'cyan': '#56B4E9',
    'pink': '#CC79A7',
}

CATEGORY_PALETTE = [
    '#0072B2', '#E69F00', '#009E73', '#CC79A7', '#D55E00',
    '#56B4E9', '#F0E442', '#999999', '#882255', '#44AA99',
    '#332288', '#DDCC77', '#117733', '#AA4499', '#88CCEE',
]


# =============================================================================
# Data Loading
# =============================================================================

def load_counts(filepath: str) -> pd.DataFrame:
    if filepath.endswith('.txt') or filepath.endswith('.tsv'):
        return pd.read_csv(filepath, sep='\t')
    return pd.read_csv(filepath)


def load_results(filepath: str) -> pd.DataFrame:
    return pd.read_csv(filepath)


def load_categories(filepath: str) -> pd.DataFrame:
    return pd.read_csv(filepath)


def load_colors(filepath: str) -> Dict[str, str]:
    """
    Load a custom color mapping from a CSV file.

    Supported formats:

    1. Two-column (category → color):
        category,color
        UFMylation,#E69F00
        Tumor suppressor,#0072B2

    2. One-column (ordered palette, replaces CATEGORY_PALETTE order):
        color
        #E69F00
        #0072B2

    Returns
    -------
    dict
        {category_name: hex_color} for two-column input, or
        {'__palette__': [color, ...]} for one-column input.
    """
    df = pd.read_csv(filepath)
    df.columns = [c.strip() for c in df.columns]
    if df.shape[1] >= 2:
        col_cat = df.columns[0]
        col_col = df.columns[1]
        return dict(zip(df[col_cat].astype(str), df[col_col].astype(str)))
    else:
        col = df.columns[0]
        return {'__palette__': df[col].tolist()}


def _resolve_cat_colors(
    cat_names: List[str],
    colors: Optional[Dict] = None
) -> Dict[str, str]:
    """
    Build a category → color dict.

    Priority:
      1. Explicit category→color entries from `colors`
      2. Ordered palette from `colors` (if __palette__ key present)
      3. Default CATEGORY_PALETTE
    """
    if colors and '__palette__' in colors:
        palette = colors['__palette__']
    else:
        palette = CATEGORY_PALETTE

    cat_colors = {cat: palette[i % len(palette)] for i, cat in enumerate(cat_names)}

    if colors:
        for cat in cat_names:
            if cat in colors:
                cat_colors[cat] = colors[cat]

    return cat_colors


def deduplicate_by_max_abs_epsilon(
    df: pd.DataFrame,
    key_col: str,
    epsilon_col: str = 'epsilon'
) -> pd.DataFrame:
    if key_col not in df.columns:
        return df
    if not df[key_col].duplicated().any():
        return df
    n_before = len(df)
    df = df.copy()
    df['_abs_epsilon'] = df[epsilon_col].abs()
    df = df.sort_values('_abs_epsilon', ascending=False)
    df = df.drop_duplicates(subset=key_col, keep='first')
    df = df.drop(columns=['_abs_epsilon'])
    n_after = len(df)
    if n_before != n_after:
        print(f"  Deduplicated '{key_col}': {n_before} → {n_after} rows (kept max |{epsilon_col}|)")
    return df.reset_index(drop=True)


# =============================================================================
# Correlation Analysis
# =============================================================================

def compute_correlation_matrix(
    df: pd.DataFrame,
    samples: List[str],
    method: str = 'pearson',
    log_transform: bool = True,
    pseudocount: float = 1.0
) -> pd.DataFrame:
    data = df[samples].copy()
    if log_transform:
        data = np.log10(data + pseudocount)
    if method == 'pearson':
        corr = data.corr(method='pearson')
    else:
        corr = data.corr(method='spearman')
    return corr


# =============================================================================
# Replicate Scatter Plots
# =============================================================================

def plot_replicate_scatter(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    output_file: str,
    gene_col: Optional[str] = None,
    highlight_genes: Optional[List[str]] = None,
    categories: Optional[pd.DataFrame] = None,
    log_scale: bool = True,
    pseudocount: float = 1.0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 8),
    point_size: int = 10,
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    ref_lines: Optional[List[Tuple[float, float]]] = None,
    show_correlation: bool = True,
    show_legend: bool = True,
    colors: Optional[Dict[str, str]] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_diagonal: bool = True,
    show_axis_lines: bool = False,
    diagonal_style: str = '--',
    axis_line_style: str = '--'
) -> None:
    fig, ax = plt.subplots(figsize=figsize)

    x = df[x_col].values.astype(float)
    y = df[y_col].values.astype(float)

    if log_scale:
        x = np.log10(x + pseudocount)
        y = np.log10(y + pseudocount)

    mask = np.isfinite(x) & np.isfinite(y)
    r_pearson, _ = pearsonr(x[mask], y[mask])
    r_spearman, _ = spearmanr(x[mask], y[mask])

    gene_to_cat = {}
    cat_names = []
    cat_colors = {}
    categorized_genes = set()

    if categories is not None and gene_col is not None:
        gene_to_cat = _gene_to_cat_map(categories)
        cat_names = list(categories.iloc[:, 1].unique())
        cat_colors = _resolve_cat_colors(cat_names, colors)
        categorized_genes = set(gene_to_cat.keys())

    highlight_set = set()
    if highlight_genes:
        highlight_set = set(g.upper() for g in highlight_genes)

    general_mask = np.ones(len(df), dtype=bool)
    if gene_col is not None:
        genes = df[gene_col].astype(str)
        for i, gene in enumerate(genes):
            if gene in categorized_genes or gene.upper() in highlight_set:
                general_mask[i] = False

    ax.scatter(x[general_mask], y[general_mask],
               c=COLORS['grey'], s=point_size, alpha=0.4,
               edgecolors='none', zorder=1)

    if gene_col is not None:
        genes = df[gene_col].astype(str)
        ntc_mask = is_ntc(genes)
        if ntc_mask.any():
            ax.scatter(x[ntc_mask], y[ntc_mask],
                       c=COLORS['black'], s=point_size, alpha=0.4,
                       edgecolors='none', zorder=2)

    if categories is not None and gene_col is not None:
        genes = df[gene_col].astype(str)
        for cat in cat_names:
            cat_mask = np.array([gene_to_cat.get(g, None) == cat for g in genes])
            if cat_mask.any():
                ax.scatter(x[cat_mask], y[cat_mask],
                           c=cat_colors[cat], s=point_size * 1.5, alpha=0.7,
                           edgecolors='none', zorder=3)

    if highlight_genes and gene_col:
        for gene in highlight_genes:
            idx = df[df[gene_col].astype(str).str.upper() == gene.upper()].index
            if len(idx) > 0:
                i = idx[0]
                xi, yi = x[i], y[i]
                ax.scatter(xi, yi, c=COLORS['red'], s=50, alpha=0.7,
                           zorder=5, edgecolors='black', linewidths=0.5)
                ax.annotate(display_label(gene), (xi, yi), fontsize=fontsize or 8,
                           xytext=(5, 5), textcoords='offset points')

    # Apply axis limits
    ax.set_xlim(xlim if xlim is not None else ax.get_xlim())
    ax.set_ylim(ylim if ylim is not None else ax.get_ylim())

    # Diagonal reference line
    if show_diagonal:
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        lims = [min(cur_xlim[0], cur_ylim[0]), max(cur_xlim[1], cur_ylim[1])]
        ax.plot(lims, lims, diagonal_style, c=COLORS['grey'], linewidth=1, zorder=0)
        ax.set_xlim(cur_xlim)
        ax.set_ylim(cur_ylim)

    # Extra reference lines (y = slope*x + intercept)
    if ref_lines:
        cur_xlim = ax.get_xlim()
        x_range = np.array(cur_xlim)
        for slope, intercept in ref_lines:
            y_line = slope * x_range + intercept
            ax.plot(x_range, y_line, diagonal_style, c=COLORS['grey'], linewidth=0.8, zorder=0)

    # Horizontal / vertical axis lines
    if show_axis_lines:
        ax.axhline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)
        ax.axvline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    legend_fontsize = (fontsize - 2) if fontsize else 8
    annot_fontsize = fontsize or 10

    if show_legend:
        handles = []
        if categories is not None:
            for cat in cat_names:
                handles.append(mpatches.Patch(color=cat_colors[cat], label=cat))
        handles.append(mpatches.Patch(color=COLORS['grey'], label='Other'))
        if gene_col is not None:
            handles.append(mpatches.Patch(color=COLORS['black'], label='NTC'))
        if handles:
            if legend_outside:
                ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=legend_fontsize)
            else:
                ax.legend(handles=handles, loc='lower right', fontsize=legend_fontsize)

    xlabel = f'$\\log_{{10}}$({x_col} + {pseudocount})' if log_scale else x_col
    ylabel = f'$\\log_{{10}}$({y_col} + {pseudocount})' if log_scale else y_col
    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)

    if title:
        ax.set_title(title, fontsize=title_fontsize)

    if show_correlation:
        ax.text(0.05, 0.95, f'Pearson r = {r_pearson:.3f}\nSpearman ρ = {r_spearman:.3f}',
                transform=ax.transAxes, fontsize=annot_fontsize, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {output_file}")


# =============================================================================
# Correlation Heatmaps
# =============================================================================

def plot_correlation_heatmap(
    df: pd.DataFrame,
    samples: List[str],
    output_file: str,
    method: str = 'pearson',
    log_transform: bool = True,
    cluster: bool = True,
    sample_order: Optional[List[str]] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 8),
    cmap: str = 'RdBu_r',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    annot: bool = True,
    mask_upper: bool = False
) -> None:
    """Correlation heatmap.

    mask_upper=True hides the strict upper triangle (diagonal kept), for the
    triangle-only layout used in Fig S1D. Ignored when cluster=True, since a
    clustermap reorders rows and columns and a fixed triangle would be
    meaningless.
    """
    corr = compute_correlation_matrix(df, samples, method, log_transform)

    if sample_order is not None:
        missing = set(sample_order) - set(corr.columns)
        if missing:
            raise ValueError(f"Samples not found in data: {missing}")
        corr = corr.loc[sample_order, sample_order]
        cluster = False

    if vmin is None or vmax is None:
        mask = ~np.eye(len(corr), dtype=bool)
        off_diag = corr.values[mask]
        data_min = off_diag.min()
        data_max = off_diag.max()
        if data_min >= 0:
            vmin, vmax = 0.0, 1.0
        elif data_max <= 0:
            vmin, vmax = -1.0, 0.0
        else:
            vmin, vmax = -1.0, 1.0

    n_samples = len(samples)
    if n_samples <= 8:       annot_fontsize = 10
    elif n_samples <= 12:    annot_fontsize = 8
    elif n_samples <= 20:    annot_fontsize = 6
    elif n_samples <= 30:    annot_fontsize = 5
    else:                    annot_fontsize = 4

    if n_samples > 15:
        scale = n_samples / 10
        figsize = (figsize[0] * scale, figsize[1] * scale)

    annot_kws = {'size': annot_fontsize}

    if cluster and mask_upper:
        print("  NOTE: mask_upper ignored because cluster=True (clustering reorders axes)")

    if cluster:
        g = sns.clustermap(
            corr, figsize=figsize, cmap=cmap, vmin=vmin, vmax=vmax,
            annot=annot, fmt='.2f', linewidths=0.5, method='average',
            annot_kws=annot_kws
        )
        tick_fontsize = min(10, max(6, 12 - n_samples // 4))
        g.ax_heatmap.tick_params(axis='both', labelsize=tick_fontsize)
        if title:
            g.fig.suptitle(title, y=1.02, fontsize=14)
        g.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
    else:
        fig, ax = plt.subplots(figsize=figsize)
        tri_mask = np.triu(np.ones_like(corr, dtype=bool), k=1) if mask_upper else None
        sns.heatmap(
            corr, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
            annot=annot, fmt='.2f', linewidths=0.5, square=True,
            annot_kws=annot_kws, mask=tri_mask
        )
        tick_fontsize = min(10, max(6, 12 - n_samples // 4))
        ax.tick_params(axis='both', labelsize=tick_fontsize)
        if title:
            ax.set_title(title, fontsize=14)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()

    plt.close()
    print(f"Saved: {output_file}")


def plot_foldchange_heatmap(
    results_files: List[str],
    output_file: str,
    comparison_names: Optional[List[str]] = None,
    method: str = 'pearson',
    cluster: bool = True,
    sample_order: Optional[List[str]] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 8),
    epsilon_col: str = 'epsilon',
    merge_on: Optional[str] = None,
    cmap: str = 'RdBu_r',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    annot: bool = True,
    mask_upper: bool = False
) -> None:
    all_fc = {}
    for i, f in enumerate(results_files):
        df = load_results(f)
        name = comparison_names[i] if comparison_names else Path(f).stem
        if merge_on:
            gene_col = merge_on
        elif 'gene_id' in df.columns:
            gene_col = 'gene_id'
        elif 'gene' in df.columns:
            gene_col = 'gene'
        else:
            gene_col = df.columns[0]
        df = deduplicate_by_max_abs_epsilon(df, gene_col, epsilon_col)
        all_fc[name] = df.set_index(gene_col)[epsilon_col]

    fc_df = pd.DataFrame(all_fc).dropna()
    print(f"Fold-change heatmap: {len(fc_df)} genes in common (merged on '{gene_col}')")

    plot_correlation_heatmap(
        fc_df, samples=list(all_fc.keys()), output_file=output_file,
        method=method, log_transform=False, cluster=cluster,
        sample_order=sample_order, title=title or 'Fold-change Correlation',
        figsize=figsize, cmap=cmap, vmin=vmin, vmax=vmax,
        annot=annot, mask_upper=mask_upper
    )


# =============================================================================
# Cross-Comparison Scatter Plots
# =============================================================================

@dataclass
class ComparisonData:
    name: str
    results: pd.DataFrame
    hits: pd.DataFrame

    @property
    def hit_genes(self) -> Set[str]:
        if 'gene_id' in self.hits.columns:
            return set(self.hits['gene_id'])
        elif 'gene' in self.hits.columns:
            return set(self.hits['gene'])
        return set(self.hits.iloc[:, 0])


def load_comparison(
    results_file: str,
    hits_file: str,
    name: Optional[str] = None
) -> ComparisonData:
    results = load_results(results_file)
    hits = load_results(hits_file)
    if name is None:
        name = Path(results_file).stem.replace('_all_genes', '')
    return ComparisonData(name=name, results=results, hits=hits)


def _line_score(parallel_dist, perp_dist, mode: str = 'perp_weighted', eps: float = 1e-10):
    """Per-sgRNA line score from signed along-line (parallel) and off-line
    (perpendicular) distances to the reference line.

    mode='perp_weighted' (default, original behavior):
        parallel * |parallel| / (|parallel| + perp). The signed projection is
        scaled by the fraction |parallel| / (|parallel| + perp), so guides that
        lie on the line approach their full projection while guides far off the
        line shrink toward zero.
    mode='parallel':
        the raw signed distance along the line, with no perpendicular
        normalization.
    """
    parallel_dist = np.asarray(parallel_dist, dtype=float)
    perp_dist = np.asarray(perp_dist, dtype=float)
    if mode == 'parallel':
        return parallel_dist
    if mode == 'perp_weighted':
        abs_parallel = np.abs(parallel_dist)
        return parallel_dist * abs_parallel / (abs_parallel + perp_dist + eps)
    raise ValueError(f"score_mode must be 'perp_weighted' or 'parallel', got {mode!r}")


def calculate_line_score(
    comp1: ComparisonData,
    comp2: ComparisonData,
    slope: float,
    intercept: float,
    epsilon_col: str = 'epsilon',
    merge_on: Optional[str] = None,
    score_mode: str = 'perp_weighted'
) -> pd.DataFrame:
    if merge_on:
        gene_id_col1 = merge_on
        gene_id_col2 = merge_on
    else:
        gene_id_col1 = 'gene_id' if 'gene_id' in comp1.results.columns else 'gene'
        gene_id_col2 = 'gene_id' if 'gene_id' in comp2.results.columns else 'gene'
    gene_col1 = 'gene' if 'gene' in comp1.results.columns else gene_id_col1
    gene_col2 = 'gene' if 'gene' in comp2.results.columns else gene_id_col2

    results1 = deduplicate_by_max_abs_epsilon(comp1.results, gene_id_col1, epsilon_col)
    results2 = deduplicate_by_max_abs_epsilon(comp2.results, gene_id_col2, epsilon_col)

    df1 = results1[[gene_id_col1, epsilon_col]].copy()
    df1.columns = ['gene_id', 'x_epsilon']
    df1['gene'] = results1[gene_col1].values

    df2 = results2[[gene_id_col2, epsilon_col]].copy()
    df2.columns = ['gene_id', 'y_epsilon']

    merged = df1.merge(df2, on='gene_id', how='inner')

    x = merged['x_epsilon'].values
    y = merged['y_epsilon'].values
    m, b = slope, intercept
    norm = np.sqrt(1 + m**2)

    merged['perp_dist'] = np.abs(y - m*x - b) / norm
    merged['parallel_dist'] = (x + m * (y - b)) / norm

    parallel = merged['parallel_dist'].values
    perp = merged['perp_dist'].values
    merged['score'] = _line_score(parallel, perp, mode=score_mode)

    merged['hit_in_' + comp1.name] = merged['gene_id'].isin(comp1.hit_genes)
    merged['hit_in_' + comp2.name] = merged['gene_id'].isin(comp2.hit_genes)

    return merged.sort_values('score', ascending=False).reset_index(drop=True)


def sgrna_line_score_utest(
    sgrna_file1: str,
    sgrna_file2: str,
    slope: float,
    intercept: float,
    output_prefix: str,
    comp1_name: Optional[str] = None,
    comp2_name: Optional[str] = None,
    top_n: int = 3,
    fdr: float = 0.05,
    min_count: int = 20,
    score_mode: str = 'perp_weighted'
) -> Tuple[pd.DataFrame, pd.DataFrame, float]:
    from scipy.stats import mannwhitneyu
    import random
    from random import shuffle

    df1 = pd.read_table(sgrna_file1)
    df2 = pd.read_table(sgrna_file2)
    df1 = df1[df1['treat_mean'] > min_count].copy()
    df2 = df2[df2['treat_mean'] > min_count].copy()

    if comp1_name is None: comp1_name = Path(sgrna_file1).parent.name
    if comp2_name is None: comp2_name = Path(sgrna_file2).parent.name

    print(f"sgRNA line score U-test analysis:")
    print(f"  Comparison 1: {comp1_name} ({len(df1)} sgRNAs)")
    print(f"  Comparison 2: {comp2_name} ({len(df2)} sgRNAs)")
    print(f"  Reference line: y = {slope}x + {intercept}")
    print(f"  Top N sgRNAs per gene: {top_n}")

    merged = df1[['sgrna', 'Gene', 'LFC']].merge(
        df2[['sgrna', 'LFC']], on='sgrna', how='inner', suffixes=('_1', '_2')
    )
    merged.columns = ['sgrna', 'gene', 'lfc1', 'lfc2']
    print(f"  Merged sgRNAs: {len(merged)}")

    x, y = merged['lfc1'].values, merged['lfc2'].values
    m, b = slope, intercept
    norm = np.sqrt(1 + m**2)
    perp_dist = np.abs(y - m*x - b) / norm
    parallel_dist = (x + m * (y - b)) / norm
    merged['score'] = _line_score(parallel_dist, perp_dist, mode=score_mode)
    merged['parallel_dist'] = parallel_dist
    merged['perp_dist'] = perp_dist

    ntc_mask = is_ntc(merged['gene'])
    df_ntc = merged[ntc_mask]
    df_targeting = merged[~ntc_mask]
    print(f"  NTC sgRNAs: {len(df_ntc)}")
    print(f"  Targeting sgRNAs: {len(df_targeting)}")

    ntc_scores = df_ntc['score'].tolist()
    genes = df_targeting['gene'].unique()
    gene_results = {}

    for gene in genes:
        df_gene = df_targeting[df_targeting['gene'] == gene].copy()
        df_gene['abs_score'] = df_gene['score'].abs()
        df_gene = df_gene.sort_values('abs_score', ascending=False)
        top_sgrnas = df_gene.head(top_n)
        mean_score = top_sgrnas['score'].mean()
        if len(ntc_scores) > 0 and len(top_sgrnas) > 0:
            _, pvalue = mannwhitneyu(top_sgrnas['score'].tolist(), ntc_scores, alternative='two-sided')
        else:
            pvalue = 1.0
        gene_results[gene] = [mean_score, pvalue, len(df_gene)]

    random.seed(42)
    ntc_scores_for_null = list(zip(df_ntc['score'], df_ntc['gene']))
    for j in range(len(genes)):
        if len(ntc_scores_for_null) >= 5:
            shuffle(ntc_scores_for_null)
            selected = ntc_scores_for_null[:5]
            selected_scores = [v[0] for v in selected]
            selected_sorted = sorted(selected_scores, key=lambda v: abs(v), reverse=True)[:top_n]
            mean_score = np.mean(selected_sorted)
            _, pval = mannwhitneyu(selected_sorted, ntc_scores, alternative='two-sided')
            gene_results[f'NTC_{j}'] = [mean_score, pval, 5]

    results = pd.DataFrame(gene_results).T
    results.columns = ['epsilon', 'pvalue', 'n_sgrnas']
    results = results.reset_index().rename(columns={'index': 'gene_id'})
    results['gene'] = results['gene_id'].apply(lambda v: v.split('_')[0] if not v.startswith('NTC_') else 'NTC')

    ntc_med = results[results['gene'] == 'NTC']['epsilon'].median()
    if pd.notna(ntc_med):
        results['epsilon'] -= ntc_med
        print(f"  NTC median score: {ntc_med:.4f} (subtracted)")

    results['product'] = results['epsilon'] * (-np.log10(results['pvalue'].clip(lower=1e-300)))

    threshold = 0
    for t in np.arange(0, abs(results['product']).max(), 0.1):
        above = results[abs(results['product']) > t]
        if len(above) > 0:
            ntc_frac = above['gene_id'].str.startswith('NTC_').sum() / len(above)
            if ntc_frac < fdr:
                threshold = t
                break

    hits = results[abs(results['product']) > threshold].copy()
    hits = hits[~hits['gene_id'].str.startswith('NTC_')]
    print(f"  Product threshold (FDR={fdr}): {threshold:.2f}")
    print(f"  Hits: {len(hits)} genes")

    results.sort_values('product', ascending=False).to_csv(f'{output_prefix}_all_genes.csv', index=False)
    hits.sort_values('product', ascending=False).to_csv(f'{output_prefix}_hits_fdr{fdr}.csv', index=False)
    merged.to_csv(f'{output_prefix}_sgrna_scores.csv', index=False)
    print(f"  Saved: {output_prefix}_all_genes.csv")
    print(f"  Saved: {output_prefix}_hits_fdr{fdr}.csv")
    print(f"  Saved: {output_prefix}_sgrna_scores.csv")

    return results, hits, threshold


def sgrna_scores_to_genes(
    sgrna_scores_csv: str,
    output_prefix: str,
    top_n: int = 3,
    fdr: float = 0.05
) -> Tuple[pd.DataFrame, pd.DataFrame, float]:
    from scipy.stats import mannwhitneyu
    import random
    from random import shuffle

    merged = pd.read_csv(sgrna_scores_csv)
    print(f"sgRNA scores to gene-level analysis:")
    print(f"  Input: {sgrna_scores_csv}")
    print(f"  sgRNAs: {len(merged)}")
    print(f"  Top N per gene: {top_n}")

    ntc_mask = is_ntc(merged['gene'])
    df_ntc = merged[ntc_mask]
    df_targeting = merged[~ntc_mask]
    print(f"  NTC sgRNAs: {len(df_ntc)}")
    print(f"  Targeting sgRNAs: {len(df_targeting)}")

    ntc_scores = df_ntc['score'].tolist()
    genes = df_targeting['gene'].unique()
    gene_results = {}

    for gene in genes:
        df_gene = df_targeting[df_targeting['gene'] == gene].copy()
        df_gene['abs_score'] = df_gene['score'].abs()
        df_gene = df_gene.sort_values('abs_score', ascending=False)
        top_sgrnas = df_gene.head(top_n)
        mean_score = top_sgrnas['score'].mean()
        if len(ntc_scores) > 0 and len(top_sgrnas) > 0:
            _, pvalue = mannwhitneyu(top_sgrnas['score'].tolist(), ntc_scores, alternative='two-sided')
        else:
            pvalue = 1.0
        gene_results[gene] = [mean_score, pvalue, len(df_gene)]

    random.seed(42)
    ntc_scores_for_null = list(zip(df_ntc['score'], df_ntc['gene']))
    for j in range(len(genes)):
        if len(ntc_scores_for_null) >= 5:
            shuffle(ntc_scores_for_null)
            selected = ntc_scores_for_null[:5]
            selected_scores = [v[0] for v in selected]
            selected_sorted = sorted(selected_scores, key=lambda v: abs(v), reverse=True)[:top_n]
            mean_score = np.mean(selected_sorted)
            _, pval = mannwhitneyu(selected_sorted, ntc_scores, alternative='two-sided')
            gene_results[f'NTC_{j}'] = [mean_score, pval, 5]

    results = pd.DataFrame(gene_results).T
    results.columns = ['epsilon', 'pvalue', 'n_sgrnas']
    results = results.reset_index().rename(columns={'index': 'gene_id'})
    results['gene'] = results['gene_id'].apply(lambda v: v.split('_')[0] if not v.startswith('NTC_') else 'NTC')

    ntc_med = results[results['gene'] == 'NTC']['epsilon'].median()
    if pd.notna(ntc_med):
        results['epsilon'] -= ntc_med
        print(f"  NTC median score: {ntc_med:.4f} (subtracted)")

    results['product'] = results['epsilon'] * (-np.log10(results['pvalue'].clip(lower=1e-300)))

    threshold = 0
    for t in np.arange(0, abs(results['product']).max(), 0.1):
        above = results[abs(results['product']) > t]
        if len(above) > 0:
            ntc_frac = above['gene_id'].str.startswith('NTC_').sum() / len(above)
            if ntc_frac < fdr:
                threshold = t
                break

    hits = results[abs(results['product']) > threshold].copy()
    hits = hits[~hits['gene_id'].str.startswith('NTC_')]
    print(f"  Product threshold (FDR={fdr}): {threshold:.2f}")
    print(f"  Hits: {len(hits)} genes")

    results.sort_values('product', ascending=False).to_csv(f'{output_prefix}_all_genes.csv', index=False)
    hits.sort_values('product', ascending=False).to_csv(f'{output_prefix}_hits_fdr{fdr}.csv', index=False)
    print(f"  Saved: {output_prefix}_all_genes.csv")
    print(f"  Saved: {output_prefix}_hits_fdr{fdr}.csv")

    return results, hits, threshold


def plot_gene_line_scores(
    gene_scores_csv: str,
    comp1_results: str,
    comp2_results: str,
    output_file: str,
    comp1_name: Optional[str] = None,
    comp2_name: Optional[str] = None,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 0,
    label_bottom_n: int = 0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 10),
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    cmap: str = 'coolwarm',
    merge_on: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_legend: bool = True,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_axis_lines: bool = True,
    axis_line_style: str = '-'
) -> None:
    gene_scores = pd.read_csv(gene_scores_csv)
    df1 = pd.read_csv(comp1_results)
    df2 = pd.read_csv(comp2_results)

    if comp1_name is None: comp1_name = Path(comp1_results).parent.name
    if comp2_name is None: comp2_name = Path(comp2_results).parent.name

    if merge_on:
        gene_col1 = merge_on
        gene_col2 = merge_on
    else:
        gene_col1 = 'gene_id' if 'gene_id' in df1.columns else 'gene'
        gene_col2 = 'gene_id' if 'gene_id' in df2.columns else 'gene'

    df1 = deduplicate_by_max_abs_epsilon(df1, gene_col1, 'epsilon')
    df2 = deduplicate_by_max_abs_epsilon(df2, gene_col2, 'epsilon')

    df1_sub = df1[[gene_col1, 'epsilon']].copy()
    df1_sub.columns = ['gene_id', 'x_epsilon']
    df2_sub = df2[[gene_col2, 'epsilon']].copy()
    df2_sub.columns = ['gene_id', 'y_epsilon']

    merged = df1_sub.merge(df2_sub, on='gene_id', how='inner')

    if 'epsilon' in gene_scores.columns:
        gene_scores_sub = gene_scores[['gene_id', 'epsilon', 'product']].copy()
        gene_scores_sub.columns = ['gene_id', 'line_score', 'product']
        merged = merged.merge(gene_scores_sub, on='gene_id', how='left')
        merged['line_score'] = merged['line_score'].fillna(0)

    ntc_mask = is_ntc(merged['gene_id'])
    df_ntc = merged[ntc_mask]
    df_targeting = merged[~ntc_mask]

    fig, ax = plt.subplots(figsize=figsize)
    scores = df_targeting['line_score'].values

    if vmin is None or vmax is None:
        auto_vmax = max(abs(scores.min()), abs(scores.max())) if len(scores) > 0 else 1
        if vmax is None: vmax = auto_vmax
        if vmin is None: vmin = -vmax

    if len(df_ntc) > 0:
        ax.scatter(df_ntc['x_epsilon'], df_ntc['y_epsilon'],
                   c=COLORS['grey'], s=8, alpha=0.3, zorder=1, label=f'NTC (n={len(df_ntc)})')

    scatter = ax.scatter(df_targeting['x_epsilon'], df_targeting['y_epsilon'],
                         c=scores, cmap=cmap, s=20, alpha=0.8,
                         vmin=vmin, vmax=vmax, edgecolors='none', zorder=2)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label('Line Score (ε)', fontsize=fontsize or 10)

    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)

    if show_axis_lines:
        ax.axhline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)
        ax.axvline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    legend_fontsize = (fontsize - 2) if fontsize else 8
    annot_fontsize = (fontsize - 2) if fontsize else 8

    genes_to_label = set()
    if label_genes: genes_to_label.update(label_genes)
    if label_top_n > 0:
        genes_to_label.update(df_targeting.nlargest(label_top_n, 'line_score')['gene_id'].tolist())
    if label_bottom_n > 0:
        genes_to_label.update(df_targeting.nsmallest(label_bottom_n, 'line_score')['gene_id'].tolist())

    for gene in genes_to_label:
        row = df_targeting[df_targeting['gene_id'] == gene]
        if len(row) == 0:
            row = df_targeting[df_targeting['gene_id'].str.upper() == gene.upper()]
        if len(row) > 0:
            r = row.iloc[0]
            ax.scatter(r['x_epsilon'], r['y_epsilon'],
                      c='none', s=80, edgecolors='black', linewidths=1.5, zorder=10)
            ax.annotate(display_label(gene), (r['x_epsilon'], r['y_epsilon']),
                       fontsize=annot_fontsize, xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel(f'{comp1_name} ε', fontsize=label_fontsize)
    ax.set_ylabel(f'{comp2_name} ε', fontsize=label_fontsize)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize)
    else:
        ax.set_title(f'Gene Scores ({len(df_targeting)} genes)', fontsize=title_fontsize)

    if show_legend:
        if legend_outside:
            ax.legend(fontsize=legend_fontsize, loc='center left', bbox_to_anchor=(1.15, 0.5))
        else:
            ax.legend(fontsize=legend_fontsize, loc='best')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved gene score plot: {output_file}")


def line_score_full_pipeline(
    sgrna_file1: str,
    sgrna_file2: str,
    comp1_results: str,
    comp2_results: str,
    slope: float,
    intercept: float,
    output_dir: str,
    comp1_name: Optional[str] = None,
    comp2_name: Optional[str] = None,
    top_n: int = 3,
    fdr: float = 0.05,
    min_count: int = 20,
    label_top_n: int = 10,
    label_bottom_n: int = 5,
    figsize: Tuple[float, float] = (10, 10),
    fontsize: Optional[int] = None,
    cmap: str = 'coolwarm',
    also_png: bool = False,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    import os
    os.makedirs(output_dir, exist_ok=True)

    if comp1_name is None: comp1_name = Path(sgrna_file1).parent.name
    if comp2_name is None: comp2_name = Path(sgrna_file2).parent.name

    print(f"\n{'='*60}")
    print(f"LINE SCORE FULL PIPELINE")
    print(f"{'='*60}")
    print(f"Comparison 1: {comp1_name}")
    print(f"Comparison 2: {comp2_name}")
    print(f"Reference line: y = {slope}x + {intercept}")
    print(f"Output: {output_dir}")

    print(f"\n--- Step 1: sgRNA-level scores ---")
    sgrna_csv = os.path.join(output_dir, 'sgrna_scores.csv')
    sgrna_plot = os.path.join(output_dir, 'sgrna_scatter.pdf')

    sgrna_scores = plot_sgrna_score_scatter(
        sgrna_file1, sgrna_file2, slope, intercept, sgrna_plot,
        comp1_name=comp1_name, comp2_name=comp2_name, min_count=min_count,
        label_top_n=label_top_n, label_bottom_n=label_bottom_n,
        figsize=figsize, fontsize=fontsize, cmap=cmap,
        score_csv=sgrna_csv, vmin=vmin, vmax=vmax
    )

    if also_png:
        sgrna_plot_png = os.path.join(output_dir, 'sgrna_scatter.png')
        plot_sgrna_score_scatter(
            sgrna_file1, sgrna_file2, slope, intercept, sgrna_plot_png,
            comp1_name=comp1_name, comp2_name=comp2_name, min_count=min_count,
            label_top_n=label_top_n, label_bottom_n=label_bottom_n,
            figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax
        )

    print(f"\n--- Step 2: Gene-level aggregation ---")
    gene_prefix = os.path.join(output_dir, 'gene_scores')
    gene_results, gene_hits, threshold = sgrna_scores_to_genes(
        sgrna_csv, gene_prefix, top_n=top_n, fdr=fdr
    )

    print(f"\n--- Step 3: Gene-level plot ---")
    gene_plot = os.path.join(output_dir, 'gene_scatter.pdf')
    gene_csv = f"{gene_prefix}_all_genes.csv"

    plot_gene_line_scores(
        gene_csv, comp1_results, comp2_results, gene_plot,
        comp1_name=comp1_name, comp2_name=comp2_name,
        label_top_n=label_top_n, label_bottom_n=label_bottom_n,
        figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax
    )

    if also_png:
        gene_plot_png = os.path.join(output_dir, 'gene_scatter.png')
        plot_gene_line_scores(
            gene_csv, comp1_results, comp2_results, gene_plot_png,
            comp1_name=comp1_name, comp2_name=comp2_name,
            label_top_n=label_top_n, label_bottom_n=label_bottom_n,
            figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax
        )

    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"Outputs:")
    print(f"  sgRNA scores: {sgrna_csv}")
    print(f"  sgRNA plot:   {sgrna_plot}")
    print(f"  Gene scores:  {gene_csv}")
    print(f"  Gene hits:    {gene_prefix}_hits_fdr{fdr}.csv")
    print(f"  Gene plot:    {gene_plot}")
    print(f"\nHits: {len(gene_hits)} genes (FDR={fdr})")

    return sgrna_scores, gene_results, gene_hits


def plot_score_comparison(
    score_files: List[str],
    score_names: List[str],
    x_score: str,
    y_score: str,
    output_file: str,
    gene_col: str = 'gene_id',
    score_col: str = 'epsilon',
    color_by: Optional[str] = None,
    color_score_col: str = 'product',
    categories: Optional[pd.DataFrame] = None,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 10),
    fontsize: Optional[int] = None,
    cmap: str = 'coolwarm',
    show_correlation: bool = True,
    show_legend: bool = True,
    colors: Optional[Dict[str, str]] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_diagonal: bool = True,
    show_axis_lines: bool = True,
    diagonal_style: str = '--',
    axis_line_style: str = '-'
) -> pd.DataFrame:
    dfs = {}
    for filepath, name in zip(score_files, score_names):
        df = pd.read_csv(filepath)
        ntc_mask = is_ntc(df[gene_col])
        df = df[~ntc_mask].copy()
        dfs[name] = df

    merged = dfs[score_names[0]][[gene_col, score_col]].copy()
    merged.columns = [gene_col, score_names[0]]
    for name in score_names[1:]:
        df = dfs[name][[gene_col, score_col]].copy()
        df.columns = [gene_col, name]
        merged = merged.merge(df, on=gene_col, how='inner')

    if color_by and color_by in dfs:
        color_df = dfs[color_by][[gene_col, color_score_col]].copy()
        color_df.columns = [gene_col, 'color_value']
        merged = merged.merge(color_df, on=gene_col, how='left')
        merged['color_value'] = merged['color_value'].fillna(0)

    print(f"Score comparison:")
    print(f"  Scores: {score_names}")
    print(f"  X-axis: {x_score}")
    print(f"  Y-axis: {y_score}")
    print(f"  Genes: {len(merged)}")

    gene_to_cat = {}
    cat_names = []
    cat_colors = {}
    categorized_genes = set()

    if categories is not None:
        gene_to_cat = _gene_to_cat_map(categories)
        cat_names = category_order if category_order is not None else sorted(categories.iloc[:, 1].unique())
        cat_colors = _resolve_cat_colors(cat_names, colors)
        categorized_genes = set(gene_to_cat.keys())

    genes_to_label = set()
    if label_genes: genes_to_label.update(g.upper() for g in label_genes)
    if label_top_n > 0:
        merged['combined_effect'] = np.abs(merged['x_epsilon']) + np.abs(merged['y_epsilon'])
        genes_to_label.update(g.upper() for g in merged.nlargest(label_top_n, 'combined_effect')['gene'].tolist())(categories.iloc[:, 1].unique())
        cat_colors = _resolve_cat_colors(cat_names, colors)
        categorized_genes = set(gene_to_cat.keys())
        merged['category'] = merged[gene_col].astype(str).map(lambda g: gene_to_cat.get(g, None))

    fig, ax = plt.subplots(figsize=figsize)
    x = merged[x_score].values
    y = merged[y_score].values

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    annot_fontsize = (fontsize - 2) if fontsize else 8
    legend_fontsize = (fontsize - 2) if fontsize else 8

    if categories is not None:
        non_cat = merged[merged['category'].isna()]
        ax.scatter(non_cat[x_score], non_cat[y_score],
                   c=COLORS['grey'], s=15, alpha=0.4, zorder=1, label=f'Other (n={len(non_cat)})')
        for cat in cat_names:
            cat_df = merged[merged['category'] == cat]
            ax.scatter(cat_df[x_score], cat_df[y_score],
                       c=cat_colors[cat], s=30, alpha=0.8, label=f'{cat} (n={len(cat_df)})',
                       edgecolors='black', linewidths=0.3, zorder=2)
    elif color_by:
        colors_vals = merged['color_value'].values
        cv_max = max(abs(colors_vals.min()), abs(colors_vals.max())) if len(colors_vals) > 0 else 1
        scatter = ax.scatter(x, y, c=colors_vals, cmap=cmap, s=20, alpha=0.7,
                             vmin=-cv_max, vmax=cv_max, edgecolors='none')
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
        cbar.set_label(f'{color_by} {color_score_col}', fontsize=fontsize or 10)
    else:
        ax.scatter(x, y, c=COLORS['blue'], s=20, alpha=0.7)

    # Apply axis limits
    ax.set_xlim(xlim if xlim is not None else ax.get_xlim())
    ax.set_ylim(ylim if ylim is not None else ax.get_ylim())

    # Diagonal reference line (y=x)
    if show_diagonal:
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        lims = [min(cur_xlim[0], cur_ylim[0]), max(cur_xlim[1], cur_ylim[1])]
        ax.plot(lims, lims, diagonal_style, c='grey', linewidth=1, zorder=0, label='y=x')
        ax.set_xlim(cur_xlim)
        ax.set_ylim(cur_ylim)

    if show_axis_lines:
        ax.axhline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)
        ax.axvline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)

    genes_to_label = set()
    if label_genes: genes_to_label.update(label_genes)
    if label_top_n > 0:
        merged['diff'] = (merged[x_score] - merged[y_score]).abs()
        genes_to_label.update(merged.nlargest(label_top_n, 'diff')[gene_col].tolist())

    for gene in genes_to_label:
        row = merged[merged[gene_col] == gene]
        if len(row) == 0:
            row = merged[merged[gene_col].str.upper() == gene.upper()]
        if len(row) > 0:
            r = row.iloc[0]
            ax.scatter(r[x_score], r[y_score],
                      c='none', s=80, edgecolors='black', linewidths=1.5, zorder=10)
            ax.annotate(display_label(gene), (r[x_score], r[y_score]),
                       fontsize=annot_fontsize, xytext=(5, 5), textcoords='offset points')

    from scipy.stats import pearsonr, spearmanr
    r_pearson, _ = pearsonr(x, y)
    r_spearman, _ = spearmanr(x, y)

    ax.set_xlabel(f'{x_score} {score_col}', fontsize=label_fontsize)
    ax.set_ylabel(f'{y_score} {score_col}', fontsize=label_fontsize)

    if title:
        ax.set_title(title, fontsize=title_fontsize)
    elif show_correlation:
        ax.set_title(f'Score Comparison (r={r_pearson:.3f}, ρ={r_spearman:.3f})', fontsize=title_fontsize)
    else:
        ax.set_title('Score Comparison', fontsize=title_fontsize)

    if show_legend:
        ax.legend(fontsize=legend_fontsize)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()

    print(f"  Pearson r: {r_pearson:.3f}")
    print(f"  Spearman ρ: {r_spearman:.3f}")
    print(f"  Saved: {output_file}")

    return merged


def plot_sgrna_strip(
    sgrna_data: pd.DataFrame,
    genes: List[str],
    output_file: str,
    score_col: str = 'score',
    gene_col: str = 'gene',
    sgrna_col: str = 'sgrna',
    categories: Optional[pd.DataFrame] = None,
    title: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    fontsize: Optional[int] = None,
    show_median: bool = True,
    show_ntc: bool = True,
    jitter: float = 0.2,
    show_legend: bool = True,
    colors: Optional[Dict[str, str]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_axis_lines: bool = True,
    axis_line_style: str = '--'
) -> None:
    gene_to_cat = {}
    cat_colors = {}

    if categories is not None:
        gene_to_cat = _gene_to_cat_map(categories)
        cat_names = sorted(categories.iloc[:, 1].unique())
        cat_colors = _resolve_cat_colors(cat_names, colors)

    genes_upper = [g.upper() for g in genes]
    sgrna_data = sgrna_data.copy()
    sgrna_data['_gene_upper'] = sgrna_data[gene_col].astype(str).str.upper()

    ntc_mask = is_ntc(sgrna_data['_gene_upper'])
    ntc_data = sgrna_data[ntc_mask]
    gene_data = sgrna_data[sgrna_data['_gene_upper'].isin(genes_upper)]

    n_genes = len(genes) + (1 if show_ntc and len(ntc_data) > 0 else 0)
    if figsize is None:
        figsize = (max(6, 1.2 * n_genes), 6)

    fig, ax = plt.subplots(figsize=figsize)

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    tick_fontsize = (fontsize - 2) if fontsize else 10

    x_positions = []
    x_labels = []
    x_idx = 0

    if show_ntc and len(ntc_data) > 0:
        y_vals = ntc_data[score_col].values
        x_vals = np.random.uniform(-jitter, jitter, len(y_vals)) + x_idx
        ax.scatter(x_vals, y_vals, c=COLORS['grey'], s=25, alpha=0.5, edgecolors='none', zorder=2)
        if show_median:
            ax.hlines(np.median(y_vals), x_idx - 0.3, x_idx + 0.3, colors='black', linewidths=2, zorder=3)
        x_positions.append(x_idx)
        x_labels.append(f'NTC\n(n={len(y_vals)})')
        x_idx += 1

    for gene in genes:
        gene_upper = gene.upper()
        gene_subset = gene_data[gene_data['_gene_upper'] == gene_upper]
        if len(gene_subset) == 0:
            print(f"  Warning: No sgRNAs found for {gene}")
            continue
        y_vals = gene_subset[score_col].values
        x_vals = np.random.uniform(-jitter, jitter, len(y_vals)) + x_idx
        color = cat_colors.get(gene_to_cat[gene], COLORS['blue']) if gene in gene_to_cat else COLORS['blue']
        ax.scatter(x_vals, y_vals, c=color, s=40, alpha=0.7, edgecolors='black', linewidths=0.3, zorder=2)
        if show_median:
            ax.hlines(np.median(y_vals), x_idx - 0.3, x_idx + 0.3, colors='black', linewidths=2, zorder=3)
        x_positions.append(x_idx)
        x_labels.append(f'{gene}\n(n={len(y_vals)})')
        x_idx += 1

    if show_axis_lines:
        ax.axhline(0, c=COLORS['grey'], linewidth=1, linestyle=axis_line_style, zorder=1)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=tick_fontsize)
    ax.set_xlim(-0.5, x_idx - 0.5)
    ax.set_ylabel('Score', fontsize=label_fontsize)

    if ylim is not None:
        ax.set_ylim(ylim)

    if title:
        ax.set_title(title, fontsize=title_fontsize)

    if show_legend and categories is not None and len(cat_colors) > 0:
        legend_elements = [
            plt.scatter([], [], c=cat_colors[cat], s=40, edgecolors='black', linewidths=0.3, label=cat)
            for cat in sorted(cat_colors.keys())
            if any(gene_to_cat.get(g, None) == cat for g in genes)
        ]
        if legend_elements:
            ax.legend(handles=legend_elements, loc='best', fontsize=tick_fontsize)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved strip plot: {output_file}")


def load_sgrna_scores_for_strip(
    source: str,
    source2: Optional[str] = None,
    slope: Optional[float] = None,
    intercept: Optional[float] = None,
    min_count: int = 20,
    score_col: str = 'LFC',
    score_mode: str = 'perp_weighted'
) -> pd.DataFrame:
    if source.endswith('.csv'):
        df = pd.read_csv(source)
        if 'Gene' in df.columns:
            df = df.rename(columns={'Gene': 'gene'})
        return df

    df1 = pd.read_table(source)
    df1 = df1[df1['treat_mean'] > min_count].copy()

    if source2 is None:
        result = df1[['sgrna', 'Gene', score_col]].copy()
        result.columns = ['sgrna', 'gene', 'score']
        return result

    if slope is None:
        raise ValueError("Must provide slope and intercept for line score computation")

    df2 = pd.read_table(source2)
    df2 = df2[df2['treat_mean'] > min_count].copy()

    merged = df1[['sgrna', 'Gene', 'LFC']].merge(
        df2[['sgrna', 'LFC']], on='sgrna', how='inner', suffixes=('_1', '_2')
    )
    merged.columns = ['sgrna', 'gene', 'lfc1', 'lfc2']

    x, y = merged['lfc1'].values, merged['lfc2'].values
    m, b = slope, intercept
    norm = np.sqrt(1 + m**2)
    perp_dist = np.abs(y - m*x - b) / norm
    parallel_dist = (x + m * (y - b)) / norm
    merged['score'] = _line_score(parallel_dist, perp_dist, mode=score_mode)

    return merged


def plot_sgrna_score_scatter(
    sgrna_file1: str,
    sgrna_file2: str,
    slope: float,
    intercept: float,
    output_file: str,
    comp1_name: Optional[str] = None,
    comp2_name: Optional[str] = None,
    min_count: int = 20,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 0,
    label_bottom_n: int = 0,
    highlight_genes: Optional[List[str]] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 10),
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    cmap: str = 'coolwarm',
    score_csv: Optional[str] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_legend: bool = True,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_ref_line: bool = True,
    show_axis_lines: bool = True,
    ref_line_style: str = '--',
    axis_line_style: str = '-',
    plain: bool = False,
    ntc_color: str = '#DDDDDD',
    score_mode: str = 'perp_weighted',
    dpi: int = 300,
    cbar_separate: bool = False,
    score_sign: float = 1
) -> pd.DataFrame:
    df1 = pd.read_table(sgrna_file1)
    df2 = pd.read_table(sgrna_file2)
    df1 = df1[df1['treat_mean'] > min_count].copy()
    df2 = df2[df2['treat_mean'] > min_count].copy()

    if comp1_name is None: comp1_name = Path(sgrna_file1).parent.name
    if comp2_name is None: comp2_name = Path(sgrna_file2).parent.name

    merged = df1[['sgrna', 'Gene', 'LFC']].merge(
        df2[['sgrna', 'LFC']], on='sgrna', how='inner', suffixes=('_1', '_2')
    )
    merged.columns = ['sgrna', 'gene', 'lfc1', 'lfc2']

    x, y = merged['lfc1'].values, merged['lfc2'].values
    m, b = slope, intercept
    norm = np.sqrt(1 + m**2)
    perp_dist = np.abs(y - m*x - b) / norm
    parallel_dist = (x + m * (y - b)) / norm
    # score_sign=-1 reverses the along-line direction (e.g. to match a score defined with x and y swapped)
    merged['score'] = score_sign * _line_score(parallel_dist, perp_dist, mode=score_mode)
    merged['parallel_dist'] = parallel_dist
    merged['perp_dist'] = perp_dist

    ntc_mask = is_ntc(merged['gene'])

    if score_csv:
        merged.to_csv(score_csv, index=False)
        print(f"Saved sgRNA scores: {score_csv} ({len(merged)} sgRNAs)")

    fig, ax = plt.subplots(figsize=figsize)
    scores = merged['score'].values

    if vmin is None or vmax is None:
        auto_vmax = max(abs(scores.min()), abs(scores.max()))
        if vmax is None: vmax = auto_vmax
        if vmin is None: vmin = -vmax

    df_ntc = merged[ntc_mask]
    df_targeting = merged[~ntc_mask]

    if len(df_ntc) > 0:
        ntc_zorder = 2 if plain else 1
        ax.scatter(df_ntc['lfc1'], df_ntc['lfc2'], c=ntc_color, s=5, alpha=0.3,
                   zorder=ntc_zorder, label=f'NTC (n={len(df_ntc)})', rasterized=True)

    if plain:
        ax.scatter(df_targeting['lfc1'], df_targeting['lfc2'],
                   c=COLORS['grey'], s=12, alpha=0.5, edgecolors='none', zorder=1,
                   rasterized=True)
        cbar = None
    else:
        scatter = ax.scatter(df_targeting['lfc1'], df_targeting['lfc2'],
                             c=df_targeting['score'], cmap=cmap, s=12, alpha=0.7,
                             vmin=vmin, vmax=vmax, edgecolors='none', zorder=2,
                             rasterized=True)
        if cbar_separate:
            # Colorbar created AFTER tight_layout (below) in its own axes anchored
            # outside the main axes, so the plot region is identical to the
            # plain/flat version. Keep a handle to the scatter for later.
            cbar = None
        else:
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label('Line Score', fontsize=fontsize or 10)

    # Apply axis limits before drawing reference lines
    ax.set_xlim(xlim if xlim is not None else ax.get_xlim())
    ax.set_ylim(ylim if ylim is not None else ax.get_ylim())

    if show_ref_line:
        cur_xlim = ax.get_xlim()
        x_line = np.linspace(cur_xlim[0], cur_xlim[1], 100)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, ref_line_style, c='black', linewidth=1.5, zorder=1,
                label=f'y = {slope:.2f}x + {intercept:.2f}')

    if show_axis_lines:
        ax.axhline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)
        ax.axvline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    legend_fontsize = (fontsize - 2) if fontsize else 8
    annot_fontsize = (fontsize - 2) if fontsize else 7

    if highlight_genes:
        highlight_colors = plt.cm.tab10(np.linspace(0, 1, len(highlight_genes)))
        for i, gene in enumerate(highlight_genes):
            gene_mask = merged['gene'].str.upper() == gene.upper()
            gene_df = merged[gene_mask]
            if len(gene_df) > 0:
                ax.scatter(gene_df['lfc1'], gene_df['lfc2'],
                          c='none', s=60, edgecolors=highlight_colors[i],
                          linewidths=1.5, zorder=10, label=f'{gene} (n={len(gene_df)})')

    sgrnas_to_label = set()
    if label_genes:
        for gene in label_genes:
            sgrnas_to_label.update(merged[merged['gene'].str.upper() == gene.upper()]['sgrna'].tolist())
    if label_top_n > 0:
        sgrnas_to_label.update(df_targeting.nlargest(label_top_n, 'score')['sgrna'].tolist())
    if label_bottom_n > 0:
        sgrnas_to_label.update(df_targeting.nsmallest(label_bottom_n, 'score')['sgrna'].tolist())

    for sgrna in sgrnas_to_label:
        row = merged[merged['sgrna'] == sgrna]
        if len(row) > 0:
            r = row.iloc[0]
            ax.annotate(display_label(r['gene']), (r['lfc1'], r['lfc2']),
                       fontsize=annot_fontsize, xytext=(3, 3), textcoords='offset points', alpha=0.8)

    ax.set_xlabel(f'{comp1_name} LFC', fontsize=label_fontsize)
    ax.set_ylabel(f'{comp2_name} LFC', fontsize=label_fontsize)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize)
    else:
        ax.set_title(f'sgRNA Line Scores ({len(df_targeting)} targeting sgRNAs)', fontsize=title_fontsize)

    if show_legend:
        if legend_outside:
            ax.legend(fontsize=legend_fontsize, loc='center left', bbox_to_anchor=(1.15, 0.5))
        else:
            ax.legend(fontsize=legend_fontsize, loc='best')

    plt.tight_layout()

    # Separate colorbar: keep the main axes at exactly the size/position it has in
    # the plain version, then widen the figure and add the colorbar in its own
    # fixed axes to the right. Using fig.add_axes (not an anchored inset) avoids a
    # matplotlib bug where rasterized artists + bbox_inches='tight' + an anchored
    # inset_axes rasterize at the wrong (tiny) resolution.
    if not plain and cbar_separate:
        fig.canvas.draw()
        pos = ax.get_position()
        fw, fh = fig.get_size_inches()
        al, ab, aw, ah = pos.x0 * fw, pos.y0 * fh, pos.width * fw, pos.height * fh
        pad_in, cbar_w_in = 0.08, 0.12
        new_fw = al + aw + pad_in + cbar_w_in + 0.6   # spare width; tight-bbox trims it
        fig.set_size_inches(new_fw, fh)
        ax.set_position([al / new_fw, ab / fh, aw / new_fw, ah / fh])
        cax = fig.add_axes([(al + aw + pad_in) / new_fw, ab / fh,
                            cbar_w_in / new_fw, ah / fh])
        cbar = fig.colorbar(scatter, cax=cax)
        cbar.set_label('Line Score', fontsize=fontsize or 10)

    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.show()
    plt.close()

    print(f"Saved sgRNA scatter plot: {output_file}")
    print(f"  {len(df_targeting)} targeting sgRNAs, {len(df_ntc)} NTCs")

    return merged


def plot_score_scatter(
    score_df: pd.DataFrame,
    output_file: str,
    comp1_name: str,
    comp2_name: str,
    slope: float,
    intercept: float,
    categories: Optional[pd.DataFrame] = None,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 0,
    label_bottom_n: int = 0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 10),
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    cmap: str = 'coolwarm',
    show_legend: bool = True,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_ref_line: bool = True,
    show_axis_lines: bool = True,
    ref_line_style: str = '--',
    axis_line_style: str = '-'
) -> None:
    fig, ax = plt.subplots(figsize=figsize)

    x = score_df['x_epsilon'].values
    y = score_df['y_epsilon'].values
    scores = score_df['score'].values

    vmax = max(abs(scores.min()), abs(scores.max()))

    scatter = ax.scatter(x, y, c=scores, cmap=cmap, s=15, alpha=0.7,
                         vmin=-vmax, vmax=vmax, edgecolors='none', zorder=2)
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label('Score', fontsize=fontsize or 10)

    # Apply axis limits
    ax.set_xlim(xlim if xlim is not None else ax.get_xlim())
    ax.set_ylim(ylim if ylim is not None else ax.get_ylim())

    if show_ref_line:
        cur_xlim = ax.get_xlim()
        x_line = np.array(cur_xlim)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, ref_line_style, c='black', linewidth=1.5, zorder=1,
                label=f'y = {slope:.2f}x + {intercept:.2f}')

    if show_axis_lines:
        ax.axhline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)
        ax.axvline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)

    gene_to_cat = {}
    if categories is not None:
        gene_to_cat = _gene_to_cat_map(categories)

    genes_to_label = set()
    if label_genes: genes_to_label.update(g.upper() for g in label_genes)
    if label_top_n > 0:
        genes_to_label.update(g.upper() for g in score_df.nlargest(label_top_n, 'score')['gene'].tolist())
    if label_bottom_n > 0:
        genes_to_label.update(g.upper() for g in score_df.nsmallest(label_bottom_n, 'score')['gene'].tolist())

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    legend_fontsize = (fontsize - 2) if fontsize else 8
    annot_fontsize = fontsize or 8

    for gene in genes_to_label:
        row = score_df[score_df['gene'].str.upper() == gene.upper()]
        if len(row) > 0:
            r = row.iloc[0]
            ax.scatter(r['x_epsilon'], r['y_epsilon'],
                      c='none', s=80, edgecolors='black', linewidths=1.5, zorder=10)
            ax.annotate(display_label(gene), (r['x_epsilon'], r['y_epsilon']),
                       fontsize=annot_fontsize, xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel(f'{comp1_name} (ε)', fontsize=label_fontsize)
    ax.set_ylabel(f'{comp2_name} (ε)', fontsize=label_fontsize)

    if title:
        ax.set_title(title, fontsize=title_fontsize)

    if show_legend:
        if legend_outside:
            ax.legend(loc='center left', bbox_to_anchor=(1.15, 0.5), fontsize=legend_fontsize)
        else:
            ax.legend(loc='upper left', fontsize=legend_fontsize)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {output_file}")


def plot_comparison_scatter(
    comp1: ComparisonData,
    comp2: ComparisonData,
    output_file: str,
    categories: Optional[pd.DataFrame] = None,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 0,
    label_top_n_x: int = 0,
    label_top_n_y: int = 0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 10),
    epsilon_col: str = 'epsilon',
    pvalue_col: str = 'pvalue',
    product_col: str = 'product',
    size_by_hits: bool = True,
    show_ntc: bool = True,
    hits_csv: Optional[str] = None,
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    ref_lines: Optional[List[Tuple[float, float]]] = None,
    merge_on: Optional[str] = None,
    show_correlation: bool = True,
    show_legend: bool = True,
    colors: Optional[Dict[str, str]] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_diagonal: bool = True,
    show_axis_lines: bool = True,
    diagonal_style: str = '--',
    axis_line_style: str = '-',
    point_size: int = 20,
    other_color: str = COLORS['grey'],   # dark grey #4D4D4D
    ntc_color: str = '#DDDDDD',          # light grey
    category_order: Optional[List[str]] = None,
    dpi: int = 300,
    show_axis_labels: bool = True,
    adjust_labels: bool = False,
    category_alpha: Optional[Dict[str, float]] = None,
    rasterize_points: bool = True
) -> None:
    if merge_on:
        gene_id_col1 = merge_on
        gene_id_col2 = merge_on
        gene_id_col1 = merge_on
        gene_id_col2 = merge_on
    else:
        gene_id_col1 = 'gene_id' if 'gene_id' in comp1.results.columns else 'gene'
        gene_id_col2 = 'gene_id' if 'gene_id' in comp2.results.columns else 'gene'

    gene_col1 = 'gene' if 'gene' in comp1.results.columns else gene_id_col1
    gene_col2 = 'gene' if 'gene' in comp2.results.columns else gene_id_col2

    results1 = deduplicate_by_max_abs_epsilon(comp1.results, gene_id_col1, epsilon_col)
    results2 = deduplicate_by_max_abs_epsilon(comp2.results, gene_id_col2, epsilon_col)

    cols1 = [gene_id_col1, epsilon_col, pvalue_col]
    if product_col in results1.columns: cols1.append(product_col)
    df1 = results1[cols1].copy()
    col_names1 = ['gene_id', 'x_epsilon', 'x_pvalue']
    if product_col in results1.columns: col_names1.append('x_product')
    df1.columns = col_names1
    df1['gene'] = results1[gene_col1].values if gene_col1 != gene_id_col1 else df1['gene_id']

    cols2 = [gene_id_col2, epsilon_col, pvalue_col]
    if product_col in results2.columns: cols2.append(product_col)
    df2 = results2[cols2].copy()
    col_names2 = ['gene_id', 'y_epsilon', 'y_pvalue']
    if product_col in results2.columns: col_names2.append('y_product')
    df2.columns = col_names2
    df2['gene'] = results2[gene_col2].values if gene_col2 != gene_id_col2 else df2['gene_id']

    merged = df1.merge(df2, on='gene_id', how='inner', suffixes=('', '_y'))
    if 'gene_y' in merged.columns:
        merged = merged.drop(columns=['gene_y'])

    hits1 = comp1.hit_genes
    hits2 = comp2.hit_genes
    merged['hit_status'] = 'neither'
    merged.loc[merged['gene_id'].isin(hits1) | merged['gene_id'].isin(hits2), 'hit_status'] = 'one'
    merged.loc[merged['gene_id'].isin(hits1) & merged['gene_id'].isin(hits2), 'hit_status'] = 'both'
    merged['is_ntc'] = is_ntc(merged['gene_id'])

    fig, ax = plt.subplots(figsize=figsize)

    size_map = {'neither': point_size, 'one': point_size * 2.5, 'both': point_size * 5} if size_by_hits else {'neither': point_size, 'one': point_size, 'both': point_size}

    gene_to_cat = {}
    cat_names = []
    cat_colors = {}
    categorized_genes = set()

    if categories is not None:
        gene_to_cat = _gene_to_cat_map(categories)
        cat_names = category_order if category_order is not None else sorted(categories.iloc[:, 1].unique())
        cat_colors = _resolve_cat_colors(cat_names, colors)
        categorized_genes = set(gene_to_cat.keys())

    genes_to_label = set()
    if label_genes: genes_to_label.update(g.upper() for g in label_genes)
    if label_top_n > 0:
        merged['combined_effect'] = np.abs(merged['x_epsilon']) + np.abs(merged['y_epsilon'])
        genes_to_label.update(g.upper() for g in merged.nlargest(label_top_n, 'combined_effect')['gene'].tolist())
    if label_top_n_x > 0:
        top_x_idx = merged['x_epsilon'].abs().nlargest(label_top_n_x).index
        genes_to_label.update(g.upper() for g in merged.loc[top_x_idx, 'gene'].tolist())
    if label_top_n_y > 0:
        top_y_idx = merged['y_epsilon'].abs().nlargest(label_top_n_y).index
        genes_to_label.update(g.upper() for g in merged.loc[top_y_idx, 'gene'].tolist())

    merged['is_categorized'] = merged['gene'].astype(str).isin(categorized_genes)

    if hits_csv:
        hits_df = merged[merged['hit_status'].isin(['one', 'both'])].copy()
        hits_df['category'] = hits_df['gene'].astype(str).map(
            lambda g: gene_to_cat.get(g, '') if gene_to_cat else '')
        hits_df[f'hit_in_{comp1.name}'] = hits_df['gene_id'].isin(hits1)
        hits_df[f'hit_in_{comp2.name}'] = hits_df['gene_id'].isin(hits2)
        if 'x_product' in hits_df.columns and 'y_product' in hits_df.columns:
            hits_df['combined_product'] = hits_df['x_product'].abs() + hits_df['y_product'].abs()
            hits_df = hits_df.sort_values('combined_product', ascending=False)
        elif 'combined_effect' not in hits_df.columns:
            hits_df['combined_effect'] = np.abs(hits_df['x_epsilon']) + np.abs(hits_df['y_epsilon'])
            hits_df = hits_df.sort_values('combined_effect', ascending=False)
        else:
            hits_df = hits_df.sort_values('combined_effect', ascending=False)
        out_cols = ['gene_id', 'gene', 'category', 'hit_status',
                    f'hit_in_{comp1.name}', f'hit_in_{comp2.name}']
        rename_map = {
            'x_epsilon': f'{comp1.name}_epsilon', 'x_pvalue': f'{comp1.name}_pvalue',
            'y_epsilon': f'{comp2.name}_epsilon', 'y_pvalue': f'{comp2.name}_pvalue',
        }
        if 'x_product' in hits_df.columns: rename_map['x_product'] = f'{comp1.name}_product'
        if 'y_product' in hits_df.columns: rename_map['y_product'] = f'{comp2.name}_product'
        for old in rename_map:
            if old in hits_df.columns: out_cols.append(old)
        hits_df[out_cols].rename(columns=rename_map).to_csv(hits_csv, index=False)
        print(f"Saved hits CSV: {hits_csv} ({len(hits_df)} genes)")

    # Layer 1: general genes
    general = merged[~merged['is_ntc'] & ~merged['is_categorized']]
    for status in ['neither', 'one', 'both']:
        subset = general[general['hit_status'] == status]
        if len(subset) > 0:
            ax.scatter(subset['x_epsilon'], subset['y_epsilon'],
                      c=other_color, s=size_map[status], alpha=0.4, edgecolors='none', zorder=1,
                      rasterized=rasterize_points)

    # Layer 2: NTCs
    if show_ntc:
        ntc = merged[merged['is_ntc']]
        for status in ['neither', 'one', 'both']:
            subset = ntc[ntc['hit_status'] == status]
            if len(subset) > 0:
                ax.scatter(subset['x_epsilon'], subset['y_epsilon'],
                          c=ntc_color, s=size_map[status], alpha=0.4, edgecolors='none', zorder=2,
                          rasterized=rasterize_points)

    # Layer 3: categorized genes
    if categories is not None:
        for cat in cat_names:
            cat_mask = merged['gene'].astype(str).map(lambda g: gene_to_cat.get(g, None) == cat)
            cat_data = merged[cat_mask & ~merged['is_ntc']]
            for status in ['neither', 'one', 'both']:
                subset = cat_data[cat_data['hit_status'] == status]
                if len(subset) > 0:
                    ax.scatter(subset['x_epsilon'], subset['y_epsilon'],
                              c=cat_colors[cat], s=size_map[status] * 1.2,
                              alpha=(category_alpha or {}).get(cat, 0.7),
                              edgecolors='none', zorder=3,
                              rasterized=rasterize_points)

    # Apply axis limits
    ax.set_xlim(xlim if xlim is not None else ax.get_xlim())
    ax.set_ylim(ylim if ylim is not None else ax.get_ylim())

    # Diagonal reference line (y=x)
    if show_diagonal:
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        lims = [min(cur_xlim[0], cur_ylim[0]), max(cur_xlim[1], cur_ylim[1])]
        ax.plot(lims, lims, diagonal_style, c=COLORS['grey'], linewidth=0.8, zorder=0)
        ax.set_xlim(cur_xlim)
        ax.set_ylim(cur_ylim)

    # Horizontal / vertical axis lines
    if show_axis_lines:
        ax.axhline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)
        ax.axvline(0, c=COLORS['grey'], linewidth=0.5, linestyle=axis_line_style, zorder=0)

    # Extra reference lines
    if ref_lines:
        cur_xlim = ax.get_xlim()
        x_range = np.array(cur_xlim)
        for rslope, rintercept in ref_lines:
            y_line = rslope * x_range + rintercept
            ax.plot(x_range, y_line, diagonal_style, c=COLORS['grey'], linewidth=0.8, zorder=0)

    # Layer 4: labeled genes
    _label_texts = []
    for gene in genes_to_label:
        row = merged[merged['gene'].str.upper() == gene.upper()]
        if len(row) > 0:
            r = row.iloc[0]
            ax.scatter(r['x_epsilon'], r['y_epsilon'],
                      c='none', s=100, edgecolors=COLORS['black'],
                      linewidths=1.5, alpha=0.7, zorder=10)
            if adjust_labels:
                _label_texts.append(ax.text(r['x_epsilon'], r['y_epsilon'], display_label(gene),
                                            fontsize=fontsize or 8, zorder=11))
            else:
                ax.annotate(display_label(gene), (r['x_epsilon'], r['y_epsilon']),
                           fontsize=fontsize or 8, xytext=(5, 5), textcoords='offset points')

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    legend_fontsize = (fontsize - 2) if fontsize else 8
    annot_fontsize = fontsize or 10

    mask = np.isfinite(merged['x_epsilon']) & np.isfinite(merged['y_epsilon'])
    r, _ = pearsonr(merged.loc[mask, 'x_epsilon'], merged.loc[mask, 'y_epsilon'])

    if show_axis_labels:
        ax.set_xlabel(f'{comp1.name} (ε)', fontsize=label_fontsize)
        ax.set_ylabel(f'{comp2.name} (ε)', fontsize=label_fontsize)

    if title:
        ax.set_title(title, fontsize=title_fontsize)

    if show_correlation:
        ax.text(0.05, 0.95, f'r = {r:.3f}\nn = {len(merged)}',
                transform=ax.transAxes, fontsize=annot_fontsize, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    if show_legend:
        handles = []
        if size_by_hits:
            handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['blue'],
                                 markersize=4, label='Hit in neither'))
            handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['blue'],
                                 markersize=7, label='Hit in one'))
            handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['blue'],
                                 markersize=10, label='Hit in both'))
        if categories is not None:
            handles.append(Line2D([0], [0], marker='', color='w', label=''))
            for cat in cat_names:
                handles.append(mpatches.Patch(color=cat_colors[cat], label=cat))
            handles.append(mpatches.Patch(color=other_color, label='Other'))
        if show_ntc:
            handles.append(mpatches.Patch(color=ntc_color, label='NTC'))
        if handles:
            if legend_outside:
                ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=legend_fontsize)
            else:
                ax.legend(handles=handles, loc='lower right', fontsize=legend_fontsize)

    # Repel gene labels so they don't overlap each other / the points
    if adjust_labels and _label_texts:
        try:
            from adjustText import adjust_text
            adjust_text(
                _label_texts, ax=ax,
                expand=(1.3, 1.5),
                arrowprops=dict(arrowstyle='-', color='0.5', lw=0.4),
            )
        except Exception as _e:
            print(f"adjust_text skipped ({_e}); labels left un-repelled")

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {output_file}")


def plot_multi_comparison_matrix(
    comparisons: List[ComparisonData],
    output_file: str,
    categories: Optional[pd.DataFrame] = None,
    figsize_per_plot: float = 4,
    epsilon_col: str = 'epsilon'
) -> None:
    n = len(comparisons)
    fig, axes = plt.subplots(n, n, figsize=(figsize_per_plot * n, figsize_per_plot * n))

    for i, comp1 in enumerate(comparisons):
        for j, comp2 in enumerate(comparisons):
            ax = axes[i, j] if n > 1 else axes

            if i == j:
                gene_col = 'gene_id' if 'gene_id' in comp1.results.columns else 'gene'
                ax.hist(comp1.results[epsilon_col], bins=50, color=COLORS['blue'], alpha=0.7)
                ax.set_xlabel(f'{comp1.name} ε')
                ax.set_ylabel('Count')
            elif i > j:
                gene_col1 = 'gene_id' if 'gene_id' in comp1.results.columns else 'gene'
                gene_col2 = 'gene_id' if 'gene_id' in comp2.results.columns else 'gene'
                df1 = comp1.results[[gene_col1, epsilon_col]].copy()
                df1.columns = ['gene', 'x']
                df2 = comp2.results[[gene_col2, epsilon_col]].copy()
                df2.columns = ['gene', 'y']
                merged = df1.merge(df2, on='gene')
                ax.scatter(merged['x'], merged['y'], s=5, alpha=0.7, c=COLORS['blue'])
                r, _ = pearsonr(merged['x'], merged['y'])
                ax.text(0.05, 0.95, f'r={r:.2f}', transform=ax.transAxes, fontsize=8)
                ax.set_xlabel(f'{comp2.name}')
                ax.set_ylabel(f'{comp1.name}')
            else:
                gene_col1 = 'gene_id' if 'gene_id' in comp1.results.columns else 'gene'
                gene_col2 = 'gene_id' if 'gene_id' in comp2.results.columns else 'gene'
                df1 = comp1.results[[gene_col1, epsilon_col]].copy()
                df1.columns = ['gene', 'x']
                df2 = comp2.results[[gene_col2, epsilon_col]].copy()
                df2.columns = ['gene', 'y']
                merged = df1.merge(df2, on='gene')
                r, _ = pearsonr(merged['x'], merged['y'])
                ax.text(0.5, 0.5, f'{r:.3f}', transform=ax.transAxes,
                       fontsize=20, ha='center', va='center',
                       color=COLORS['red'] if r > 0.7 else COLORS['blue'] if r > 0.4 else COLORS['grey'])
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {output_file}")


# =============================================================================
# Category Enrichment / Violin Plot
# =============================================================================

def plot_category_enrichment(
    results_file: str,
    output_file: str,
    categories: pd.DataFrame,
    hits_file: Optional[str] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 4),
    epsilon_col: str = 'epsilon',
    orientation: str = 'vertical',
    category_order: Optional[List[str]] = None,
    colors: Optional[Dict[str, str]] = None,
    show_legend: bool = True,
    fontsize: Optional[int] = None,
    show_ntc: bool = False,
    violin: bool = True,
    strip: bool = True,
    jitter: float = 0.15,
    reference_violins: bool = False,
    show_counts: bool = False
) -> None:
    """
    Violin/strip plot of epsilon scores split by gene category.

    Parameters
    ----------
    results_file : str
        Path to _all_genes.csv
    output_file : str
        Output path (extension determines format)
    categories : DataFrame
        Gene categories (gene, category columns)
    hits_file : str, optional
        Path to hits CSV — hit genes are marked with larger points
    title : str, optional
        Plot title
    figsize : tuple
        Figure size
    epsilon_col : str
        Column to use as score (default: 'epsilon')
    orientation : str
        'vertical' (categories on x axis) or 'horizontal' (categories on y axis)
    category_order : list, optional
        Order of categories. Defaults to sorted unique values.
    colors : dict, optional
        {category: color} mapping — output of load_colors() or hand-built dict.
        If None, uses CATEGORY_PALETTE.
    show_legend : bool
        Show category color legend
    show_ntc : bool
        Include NTC distribution as a reference category
    violin : bool
        Show violin plot layer
    strip : bool
        Show individual point strip layer
    jitter : float
        Horizontal jitter for strip points (0–0.4)
    reference_violins : bool
        Prepend violin-only background distributions — 'All genes' (all
        targeting genes) and 'NTC' — as reference columns to compare the
        category violins against. These are drawn without the strip layer.
        When True, it supersedes show_ntc's separate NTC column to avoid
        drawing NTC twice.
    """
    df = load_results(results_file)
    gene_id_col = 'gene_id' if 'gene_id' in df.columns else 'gene'
    gene_col = 'gene' if 'gene' in df.columns else gene_id_col

    # Load hits
    hit_genes = set()
    if hits_file:
        hits = load_results(hits_file)
        hit_genes = set(hits[gene_id_col] if gene_id_col in hits.columns else hits.iloc[:, 0])

    # Build category mapping
    gene_to_cats = _gene_to_cats_map(categories)
    all_cats = list(categories.iloc[:, 1].unique())
    cat_names = category_order if category_order is not None else sorted(all_cats)
    cat_colors = _resolve_cat_colors(cat_names, colors)

    # Assign categories to genes
    df['_cats'] = df[gene_col].astype(str).map(lambda g: gene_to_cats.get(g, ()))
    df['_is_hit'] = df[gene_id_col].isin(hit_genes)
    df['_is_ntc'] = is_ntc(df[gene_id_col])

    # Build list of (label, scores, color, hit_mask, violin_only) per group
    plot_groups = []

    # Violin-only background distributions for comparison.
    # Greys match plot_comparison_scatter: 'other'/all genes dark, NTC light.
    if reference_violins:
        all_vals = df.loc[~df['_is_ntc'], epsilon_col].dropna().values
        if len(all_vals) > 0:
            n_all = df.loc[~df['_is_ntc'] & df[epsilon_col].notna(), gene_col].nunique()
            plot_groups.append(('All genes', all_vals, COLORS['grey'],
                                np.zeros(len(all_vals), dtype=bool), True, n_all))
        ref_ntc_vals = df.loc[df['_is_ntc'], epsilon_col].dropna().values
        if len(ref_ntc_vals) > 0:
            n_ntc = df.loc[df['_is_ntc'] & df[epsilon_col].notna(), gene_id_col].nunique()
            plot_groups.append(('NTC', ref_ntc_vals, '#DDDDDD',
                                np.zeros(len(ref_ntc_vals), dtype=bool), True, n_ntc))

    # NTC as a full violin+strip category — skipped if already a reference violin
    if show_ntc and not reference_violins:
        ntc_vals = df.loc[df['_is_ntc'], epsilon_col].dropna().values
        if len(ntc_vals) > 0:
            n_ntc = df.loc[df['_is_ntc'] & df[epsilon_col].notna(), gene_id_col].nunique()
            plot_groups.append(('NTC', ntc_vals, COLORS['grey'],
                                np.zeros(len(ntc_vals), dtype=bool), False, n_ntc))

    for cat in cat_names:
        cat_df = df[df['_cats'].map(lambda cs: cat in cs)]
        vals = cat_df[epsilon_col].dropna().values
        hits_mask = cat_df.loc[cat_df[epsilon_col].notna(), '_is_hit'].values
        if len(vals) > 0:
            n_sym = cat_df.loc[cat_df[epsilon_col].notna(), gene_col].nunique()
            plot_groups.append((cat, vals, cat_colors[cat], hits_mask, False, n_sym))

    label_fontsize = fontsize or 11
    title_fontsize = (fontsize + 2) if fontsize else 13
    tick_fontsize = (fontsize - 1) if fontsize else 9
    legend_fontsize = (fontsize - 2) if fontsize else 8

    fig, ax = plt.subplots(figsize=figsize)

    vert = orientation == 'vertical'

    for idx, (cat, vals, color, hits_mask, violin_only, _n) in enumerate(plot_groups):
        if violin and len(vals) >= 3:
            vp = ax.violinplot([vals], positions=[idx], vert=vert,
                               showmedians=True, showextrema=False)
            for pc in vp['bodies']:
                pc.set_facecolor(color)
                pc.set_alpha(0.5)
                pc.set_edgecolor('none')
            vp['cmedians'].set_color('black')
            vp['cmedians'].set_linewidth(1.5)

        if strip and not violin_only:
            np.random.seed(idx)
            jit = np.random.uniform(-jitter, jitter, len(vals))
            # non-hit points
            non_hit = ~hits_mask
            if non_hit.any():
                pos = (idx + jit[non_hit]) if vert else None
                # Strip points are left VECTOR so each gene stays an individually
                # selectable object in Illustrator. Rasterizing them was a wrong
                # call earlier: it assumed the 'All genes' / 'NTC' reference
                # columns contributed strip points, but those are violin-only
                # (reference_violins), so the real point count is only ~1350 for
                # the genome-wide panels and ~143 for the sublibrary ones --
                # nothing like the ~42,000 in the comparison scatters, which do
                # need rasterizing. Note the PDF backend also merges contiguous
                # rasterized artists into ONE image, so rasterizing here would
                # collapse every category into a single un-separable object.
                if vert:
                    ax.scatter(idx + jit[non_hit], vals[non_hit],
                               c=color, s=15, alpha=0.6, edgecolors='none', zorder=3)
                else:
                    ax.scatter(vals[non_hit], idx + jit[non_hit],
                               c=color, s=15, alpha=0.6, edgecolors='none', zorder=3)
            # hit points (larger, outlined)
            if hits_mask.any():
                if vert:
                    ax.scatter(idx + jit[hits_mask], vals[hits_mask],
                               c=color, s=40, alpha=0.9,
                               edgecolors='black', linewidths=0.5, zorder=4)
                else:
                    ax.scatter(vals[hits_mask], idx + jit[hits_mask],
                               c=color, s=40, alpha=0.9,
                               edgecolors='black', linewidths=0.5, zorder=4)

    # Axis zero line
    if vert:
        ax.axhline(0, c='black', linewidth=0.8, linestyle='--', zorder=0)
    else:
        ax.axvline(0, c='black', linewidth=0.8, linestyle='--', zorder=0)

    # Tick labels
    # Distinct gene symbols, not points: a gene with two TSS entries (_P1/_P2)
    # draws two points but is one gene, and the pasted draft labels use genes.
    if show_counts:
        _sep = '\n' if vert else ' '
        labels = [f'{g[0]}{_sep}({g[5]:,})' for g in plot_groups]
    else:
        labels = [g[0] for g in plot_groups]
    positions = list(range(len(plot_groups)))

    if vert:
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=tick_fontsize)
        ax.set_ylabel('Score (ε)', fontsize=label_fontsize)
    else:
        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=tick_fontsize)
        ax.set_xlabel('Score (ε)', fontsize=label_fontsize)

    if title:
        ax.set_title(title, fontsize=title_fontsize)

    if show_legend:
        handles = [mpatches.Patch(color=cat_colors[c], label=c) for c in cat_names if c in cat_colors]
        if handles:
            ax.legend(handles=handles, fontsize=legend_fontsize,
                      loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {output_file}")




def make_diff_table(
    file_a: str,
    file_b: str,
    value_col: str = 'score',
    name_a: str = 'A',
    name_b: str = 'B',
    merge_on: Optional[str] = None,
    diff_col: str = 'diff'
) -> pd.DataFrame:
    """
    Merge two per-gene score/effect-size tables and compute a per-gene
    difference (A minus B) suitable for feeding to plot_waterfall.

    Both inputs must contain `value_col` and a gene identifier. Genes are
    matched on `merge_on` if given, else on 'gene_id' (falling back to 'gene').
    Only genes present in both tables are kept.

    Returns a DataFrame with columns:
        gene_id, gene, {value_col}_{name_a}, {value_col}_{name_b}, {diff_col}
    where {diff_col} = A - B, sorted ascending by {diff_col}.
    """
    a = load_results(file_a)
    b = load_results(file_b)

    def _ids(df):
        gid = merge_on or ('gene_id' if 'gene_id' in df.columns else 'gene')
        gname = 'gene' if 'gene' in df.columns else gid
        return gid, gname

    gid_a, gname_a = _ids(a)
    gid_b, _ = _ids(b)

    if value_col not in a.columns or value_col not in b.columns:
        raise ValueError(f"'{value_col}' must be present in both tables; "
                         f"found {list(a.columns)} and {list(b.columns)}")

    da = a[[gid_a, value_col]].rename(columns={gid_a: 'gene_id', value_col: f'{value_col}_{name_a}'})
    da['gene'] = a[gname_a].values
    db = b[[gid_b, value_col]].rename(columns={gid_b: 'gene_id', value_col: f'{value_col}_{name_b}'})

    merged = da.merge(db, on='gene_id', how='inner')
    merged[diff_col] = merged[f'{value_col}_{name_a}'] - merged[f'{value_col}_{name_b}']
    return merged.sort_values(diff_col, ascending=True).reset_index(drop=True)


def _waterfall_from_df(
    df: pd.DataFrame,
    output_file: str,
    value_col: str = 'epsilon',
    hit_genes: Optional[Set[str]] = None,
    categories: Optional[pd.DataFrame] = None,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 6),
    ylabel: Optional[str] = None,
    ref_lines: Optional[Tuple[float, ...]] = (0.585, -0.585),
    ascending: bool = True,
    show_ntc: bool = True,
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    show_legend: bool = True,
    colors: Optional[Dict[str, str]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    category_alpha: Optional[Dict[str, float]] = None
) -> None:
    """Core waterfall renderer operating on an in-memory DataFrame.

    Ranks every gene by `value_col` and plots rank (x) vs value (y).
    `value_col` can be an effect size (epsilon), a composite score
    (product, line score), or a between-condition difference (see
    make_diff_table). `ref_lines` are dashed horizontal guides; pass None
    or () to omit them (sensible for difference/score plots where the
    LFC = log2(1.5) guide is meaningless).
    """
    df = df.copy()
    hit_genes = hit_genes or set()
    gene_id_col = 'gene_id' if 'gene_id' in df.columns else 'gene'
    gene_col = 'gene' if 'gene' in df.columns else gene_id_col

    df['is_ntc'] = is_ntc(df[gene_id_col])
    if not show_ntc:
        df = df[~df['is_ntc']].copy()

    df = df.sort_values(value_col, ascending=ascending).reset_index(drop=True)
    df['rank'] = range(len(df))

    gene_to_cat = {}
    cat_names = []
    cat_colors = {}
    categorized_genes = set()

    if categories is not None:
        gene_to_cat = _gene_to_cat_map(categories)
        cat_names = list(categories.iloc[:, 1].unique())
        cat_colors = _resolve_cat_colors(cat_names, colors)
        categorized_genes = set(gene_to_cat.keys())

    genes_to_label = set()
    if label_genes: genes_to_label.update(g.upper() for g in label_genes)
    if label_top_n > 0:
        genes_to_label.update(g.upper() for g in df.nlargest(label_top_n, value_col)[gene_col].tolist())
        genes_to_label.update(g.upper() for g in df.nsmallest(label_top_n, value_col)[gene_col].tolist())

    df['is_hit'] = df[gene_id_col].isin(hit_genes)
    df['is_categorized'] = df[gene_col].astype(str).isin(categorized_genes)
    df['size'] = df['is_hit'].map({True: 30, False: 8})

    fig, ax = plt.subplots(figsize=figsize)

    general = df[~df['is_ntc'] & ~df['is_categorized']]
    ax.scatter(general['rank'], general[value_col],
               c=COLORS['grey'], s=general['size'], alpha=0.4, edgecolors='none', zorder=1)

    if show_ntc:
        ntc = df[df['is_ntc']]
        ax.scatter(ntc['rank'], ntc[value_col],
                   c=COLORS['black'], s=ntc['size'], alpha=0.4, edgecolors='none', zorder=2)

    if categories is not None:
        for cat in cat_names:
            cat_mask = df[gene_col].astype(str).map(lambda g: gene_to_cat.get(g, None) == cat)
            cat_df = df[cat_mask]
            ax.scatter(cat_df['rank'], cat_df[value_col],
                       c=cat_colors[cat], s=cat_df['size'] * 1.5,
                       alpha=(category_alpha or {}).get(cat, 0.7),
                       edgecolors='none', zorder=3)

    ax.axhline(0, c='black', linewidth=1, zorder=0)
    for _r in (ref_lines or ()):
        ax.axhline(_r, c='black', linewidth=0.8, linestyle='--', zorder=0)

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    legend_fontsize = (fontsize - 2) if fontsize else 8
    annot_fontsize = (fontsize - 3) if fontsize else 7

    for gene in genes_to_label:
        row = df[df[gene_col].str.upper() == gene.upper()]
        if len(row) > 0:
            r = row.iloc[0]
            ax.annotate(display_label(gene), (r['rank'], r[value_col]),
                        fontsize=annot_fontsize, xytext=(3, 3), textcoords='offset points',
                        rotation=45, ha='left', va='bottom')

    if ylabel is None:
        ylabel = '$\\log_2$ Fold Change (ε)' if value_col == 'epsilon' else value_col
    ax.set_xlabel('Rank', fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize)
    if ylim is not None:
        ax.set_ylim(ylim)

    if show_legend:
        handles = []
        if categories is not None:
            for cat in cat_names:
                handles.append(mpatches.Patch(color=cat_colors[cat], label=cat))
        handles.append(mpatches.Patch(color=COLORS['grey'], label='Other'))
        if show_ntc:
            handles.append(mpatches.Patch(color=COLORS['black'], label='NTC'))
        if hit_genes:
            handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['grey'],
                                  markersize=6, label='Hit'))
            handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['grey'],
                                  markersize=3, label='Non-hit'))
        if handles:
            if legend_outside:
                ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=legend_fontsize)
            else:
                ax.legend(handles=handles, loc='upper left', fontsize=legend_fontsize)

    ax.set_xlim(-len(df) * 0.02, len(df) * 1.02)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {output_file}")


def plot_waterfall(
    results_file: str,
    output_file: str,
    hits_file: Optional[str] = None,
    categories: Optional[pd.DataFrame] = None,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 6),
    epsilon_col: str = 'epsilon',
    value_col: Optional[str] = None,
    ylabel: Optional[str] = None,
    ref_lines: Optional[Tuple[float, ...]] = (0.585, -0.585),
    ascending: bool = True,
    show_ntc: bool = True,
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    show_legend: bool = True,
    colors: Optional[Dict[str, str]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    category_alpha: Optional[Dict[str, float]] = None
) -> None:
    """Ranked waterfall of a per-gene results table.

    By default ranks/plots `epsilon` (effect size) — unchanged behavior.
    Set `value_col` to rank by any other column instead (e.g. 'product',
    'score', a composite). `epsilon_col` is retained for backward
    compatibility and is used when `value_col` is None. For non-LFC
    quantities, pass ref_lines=None (or ()) to drop the log2(1.5) guides
    and `ylabel` to relabel the axis.
    """
    vcol = value_col or epsilon_col
    df = load_results(results_file)
    gene_id_col = 'gene_id' if 'gene_id' in df.columns else 'gene'
    hit_genes: Set[str] = set()
    if hits_file:
        hits = load_results(hits_file)
        hit_genes = set(hits[gene_id_col] if gene_id_col in hits.columns else hits.iloc[:, 0])

    _waterfall_from_df(
        df, output_file, value_col=vcol, hit_genes=hit_genes,
        categories=categories, label_genes=label_genes, label_top_n=label_top_n,
        title=title, figsize=figsize, ylabel=ylabel, ref_lines=ref_lines,
        ascending=ascending, show_ntc=show_ntc, fontsize=fontsize,
        legend_outside=legend_outside, show_legend=show_legend,
        colors=colors, ylim=ylim, category_alpha=category_alpha
    )


def plot_waterfall_diff(
    file_a: str,
    file_b: str,
    output_file: str,
    value_col: str = 'score',
    name_a: str = 'A',
    name_b: str = 'B',
    merge_on: Optional[str] = None,
    categories: Optional[pd.DataFrame] = None,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 6),
    ylabel: Optional[str] = None,
    ref_lines: Optional[Tuple[float, ...]] = None,
    show_ntc: bool = True,
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    show_legend: bool = True,
    colors: Optional[Dict[str, str]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    diff_csv: Optional[str] = None,
    category_alpha: Optional[Dict[str, float]] = None
) -> None:
    """Waterfall of the per-gene DIFFERENCE in `value_col` between two tables.

    Each bar is (A - B) for one gene; genes are ranked low-to-high so the
    extremes are condition-specific. Useful for differential sensitivity /
    differential dependency between two screens or scoring schemes. By
    default no log2(1.5) guides are drawn (ref_lines=None) since a
    difference has no LFC threshold.
    """
    diff = make_diff_table(file_a, file_b, value_col=value_col,
                           name_a=name_a, name_b=name_b, merge_on=merge_on)
    if diff_csv:
        diff.to_csv(diff_csv, index=False)
        print(f"Saved diff table: {diff_csv} ({len(diff)} genes)")
    if ylabel is None:
        ylabel = f'Δ {value_col} ({name_a} − {name_b})'
    _waterfall_from_df(
        diff, output_file, value_col='diff', hit_genes=set(),
        categories=categories, label_genes=label_genes, label_top_n=label_top_n,
        title=title, figsize=figsize, ylabel=ylabel, ref_lines=ref_lines,
        ascending=True, show_ntc=show_ntc, fontsize=fontsize,
        legend_outside=legend_outside, show_legend=show_legend,
        colors=colors, ylim=ylim, category_alpha=category_alpha
    )

def plot_volcano(
    results_file: str,
    hits_file: str,
    output_file: str,
    categories: Optional[pd.DataFrame] = None,
    label_genes: Optional[List[str]] = None,
    label_top_n: int = 10,
    title: Optional[str] = None,
    fdr: float = 0.05,
    figsize: Tuple[float, float] = (10, 8),
    epsilon_col: str = 'epsilon',
    pvalue_col: str = 'pvalue',
    product_col: str = 'product',
    fontsize: Optional[int] = None,
    legend_outside: bool = False,
    show_legend: bool = True,
    colors: Optional[Dict[str, str]] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None
) -> None:
    df = load_results(results_file)
    hits = load_results(hits_file)

    gene_id_col = 'gene_id' if 'gene_id' in df.columns else 'gene'
    gene_col = 'gene' if 'gene' in df.columns else gene_id_col
    hit_genes = set(hits[gene_id_col] if gene_id_col in hits.columns else hits.iloc[:, 0])

    fig, ax = plt.subplots(figsize=figsize)

    df['is_ntc'] = is_ntc(df[gene_id_col])
    df['is_hit'] = df[gene_id_col].isin(hit_genes)
    df['is_positive'] = df['is_hit'] & (df[epsilon_col] > 0)
    df['is_negative'] = df['is_hit'] & (df[epsilon_col] < 0)

    categorized_genes = set()
    cat_names = []
    cat_colors = {}
    gene_to_cat = {}
    if categories is not None:
        gene_to_cat = _gene_to_cat_map(categories)
        cat_names = list(categories.iloc[:, 1].unique())
        cat_colors = _resolve_cat_colors(cat_names, colors)
        categorized_genes = set(gene_to_cat.keys())

    genes_to_label = set()
    if label_genes: genes_to_label.update(g.upper() for g in label_genes)
    if label_top_n > 0:
        top = hits.nlargest(label_top_n, product_col) if product_col in hits.columns else hits.head(label_top_n)
        label_col = gene_col if gene_col in top.columns else (gene_id_col if gene_id_col in top.columns else top.columns[0])
        genes_to_label.update(str(g).upper() for g in top[label_col])

    df['is_categorized'] = df[gene_col].astype(str).isin(categorized_genes)

    other = df[~df['is_hit'] & ~df['is_ntc'] & ~df['is_categorized']]
    ax.scatter(other[epsilon_col], -np.log10(other[pvalue_col]),
              c=COLORS['grey'], s=5, alpha=0.4, zorder=1)

    ntc = df[df['is_ntc']]
    ax.scatter(ntc[epsilon_col], -np.log10(ntc[pvalue_col]),
              c=COLORS['black'], s=3, alpha=0.4, zorder=2)

    pos = df[df['is_positive'] & ~df['is_categorized']]
    ax.scatter(pos[epsilon_col], -np.log10(pos[pvalue_col]),
              c=COLORS['red'], s=15, alpha=0.7, zorder=3,
              label=f'Positive (n={len(df[df["is_positive"]])})')

    neg = df[df['is_negative'] & ~df['is_categorized']]
    ax.scatter(neg[epsilon_col], -np.log10(neg[pvalue_col]),
              c=COLORS['blue'], s=15, alpha=0.7, zorder=3,
              label=f'Negative (n={len(df[df["is_negative"]])})')

    if categories is not None:
        for cat in cat_names:
            cat_mask = df[gene_col].astype(str).map(lambda g: gene_to_cat.get(g, None) == cat)
            cat_df = df[cat_mask]
            ax.scatter(cat_df[epsilon_col], -np.log10(cat_df[pvalue_col]),
                      c=cat_colors[cat], s=30, alpha=0.7, label=cat,
                      edgecolors='black', linewidths=0.5, zorder=4)

    if product_col in df.columns and len(hits) > 0:
        threshold = abs(hits[product_col]).min()
        x_curve = np.linspace(0.01, 10, 100)
        y_curve = threshold / x_curve
        ax.plot(x_curve, y_curve, '--k', linewidth=0.8)
        ax.plot(-x_curve, y_curve, '--k', linewidth=0.8, label=f'FDR={fdr}')

    label_fontsize = fontsize or 12
    title_fontsize = (fontsize + 2) if fontsize else 14
    legend_fontsize = (fontsize - 2) if fontsize else 8
    annot_fontsize = fontsize or 8

    for gene in genes_to_label:
        row = df[df[gene_col].str.upper() == gene.upper()]
        if len(row) > 0:
            r = row.iloc[0]
            ax.annotate(display_label(gene), (r[epsilon_col], -np.log10(r[pvalue_col])),
                       fontsize=annot_fontsize, xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel('Phenotype (ε)', fontsize=label_fontsize)
    ax.set_ylabel(r'$-\log_{10}(p)$', fontsize=label_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize)

    # Auto xlim if not provided
    if xlim is None:
        xlim_val = max(abs(df[epsilon_col].min()), df[epsilon_col].max()) + 0.5
        ax.set_xlim(-xlim_val, xlim_val)
    else:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if show_legend:
        if legend_outside:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=legend_fontsize)
        else:
            ax.legend(loc='upper right', fontsize=legend_fontsize)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {output_file}")


# =============================================================================
# CLI
# =============================================================================

def plot_go_enrichment_diverging(
    screen_results: List[Tuple[str, str]],
    output_file: str,
    top_n: int = 4,
    term_col: str = 'Gene_Set',
    stat_col: str = 't_statistic',
    neglog_col: str = 'neglog10_FDR',
    genes_col: str = 'genes_in_set_and_screen',
    categories: Optional[pd.DataFrame] = None,
    colors: Optional[Dict[str, str]] = None,
    term_colors: Optional[Dict[str, str]] = None,
    default_label_color: str = '#000000',
    term_labels: Optional[Dict[str, str]] = None,
    min_overlap_frac: float = 0.25,
    min_dominance: float = 2.0,
    bar_color: str = '#D3D3D3',
    break_range: Optional[Tuple[float, float]] = None,
    tick_step: float = 20.0,
    figsize: Tuple[float, float] = (7.0, 5.0),
    title: Optional[str] = None,
    xlabel: str = 'Adjusted p-value',
    left_header: str = 'Decreased New/Old Ratio',
    right_header: str = 'Increased New/Old Ratio',
    fontsize: Optional[int] = None,
    dpi: int = 300,
    strip_go_id: bool = True,
    wrap_width: int = 22,
) -> pd.DataFrame:
    """Diverging horizontal bar chart of GO enrichment for two screens (Fig 1F).

    Each screen contributes `top_n` terms per direction. A row pairs the i-th
    decreased term (bar growing left from the centre) with the i-th increased
    term (growing right), both ranked by `neglog_col`. Direction comes from the
    sign of `stat_col`.

    The x axis is signed -log10(adjusted p), so distance from the centre is
    significance in either direction; ticks are relabelled as powers of ten and
    the centre reads 1.0.

    screen_results
        [(screen_label, path_to_decoupler_csv), ...] -- drawn top to bottom.
    break_range
        (lo, hi) in -log10 units, e.g. (160, 60): hides the left-hand span
        between 10**-160 and 10**-60, which is otherwise mostly empty. The
        break is drawn on the axis and across any bar crossing it.
    categories
        Optional (gene, category) table -- the same one Fig 1E uses. Each term's
        gene list is scored against every category and the label takes that
        category's colour when the overlap fraction reaches `min_overlap_frac`.
    term_colors
        Explicit {term: colour} overrides, applied last.

    Returns the plotted table, so the selection can be inspected or exported.
    """
    label_fontsize = fontsize or 9
    header_fontsize = (fontsize or 9) + 1

    # ---- assemble the rows -------------------------------------------------
    rows = []
    for screen_label, path in screen_results:
        df = pd.read_csv(path)
        for col in (term_col, stat_col, neglog_col):
            if col not in df.columns:
                raise ValueError(f"{path}: missing column '{col}'")
        df = df.dropna(subset=[stat_col, neglog_col])
        dec = (df[df[stat_col] < 0].sort_values(neglog_col, ascending=False)
               .head(top_n).reset_index(drop=True))
        inc = (df[df[stat_col] > 0].sort_values(neglog_col, ascending=False)
               .head(top_n).reset_index(drop=True))
        for i in range(max(len(dec), len(inc))):
            entry = {'screen': screen_label, 'rank': i}
            for side, sub, sign in (('left', dec, -1.0), ('right', inc, 1.0)):
                if i < len(sub):
                    r = sub.iloc[i]
                    entry[f'{side}_term'] = r[term_col]
                    entry[f'{side}_x'] = sign * float(r[neglog_col])
                    entry[f'{side}_genes'] = r.get(genes_col, '')
            rows.append(entry)
    plotted = pd.DataFrame(rows)
    if plotted.empty:
        raise ValueError("no terms selected")

    # ---- label colours -----------------------------------------------------
    cat_genes: Dict[str, set] = {}
    if categories is not None and len(categories):
        gcol = 'gene' if 'gene' in categories.columns else categories.columns[0]
        ccol = 'category' if 'category' in categories.columns else categories.columns[-1]
        for cat, sub in categories.groupby(ccol):
            cat_genes[str(cat)] = set(sub[gcol].astype(str).str.strip())

    def _label_color(term, genes_field):
        if term_colors and term in term_colors:
            return term_colors[term]
        if not cat_genes or not colors:
            return default_label_color
        genes = {g.strip() for g in str(genes_field).split(',') if g.strip()}
        if not genes:
            return default_label_color
        fracs = sorted(
            ((len(genes & members) / len(genes), cat)
             for cat, members in cat_genes.items() if cat in colors),
            reverse=True)
        if not fracs:
            return default_label_color
        best_frac, best = fracs[0]
        runner_up = fracs[1][0] if len(fracs) > 1 else 0.0
        # Require a dominant category, not merely the largest one. Terms that
        # straddle two categories (e.g. Reticulophagy sits ~0.26 UFMylation /
        # ~0.21 Lysosome-Autophagy) are genuinely ambiguous and stay black.
        if best_frac < min_overlap_frac or best_frac < runner_up * min_dominance:
            return default_label_color
        return colors[best]

    def _clean(term):
        if term_labels and term in term_labels:
            t = term_labels[term]
            return '\n'.join(textwrap.wrap(t, wrap_width)) if wrap_width else t
        t = str(term)
        if strip_go_id:
            t = re.sub(r'\s*\(GO:\d+\)\s*$', '', t)
        return '\n'.join(textwrap.wrap(t, wrap_width)) if wrap_width else t

    # ---- geometry ----------------------------------------------------------
    xs = [v for c in ('left_x', 'right_x') if c in plotted
          for v in plotted[c].dropna()]
    xmin, xmax = min(xs), max(xs)
    pad_l, pad_r = abs(xmin) * 0.02 + 1, abs(xmax) * 0.02 + 1

    screens = list(dict.fromkeys(plotted['screen']))
    gap = 0.8
    ypos, yticks_blocks = [], {}
    y = 0.0
    for s in screens:
        n = int((plotted['screen'] == s).sum())
        block = []
        for _ in range(n):
            ypos.append(y)
            block.append(y)
            y += 1.0
        yticks_blocks[s] = block
        y += gap
    plotted['_y'] = ypos

    fig = plt.figure(figsize=figsize)
    if break_range:
        lo, hi = break_range           # e.g. (160, 60): hide -160 .. -60
        left_span = (-xmin - pad_l) - lo
        main_span = (hi + xmax + pad_r)
        gs = fig.add_gridspec(1, 2, width_ratios=[max(left_span, 1e-6), main_span],
                              wspace=0.04)
        ax_l = fig.add_subplot(gs[0, 0])
        ax_m = fig.add_subplot(gs[0, 1], sharey=ax_l)
        # Nudge both limits past the break so the 10^-lo and 10^-hi tick
        # labels don't collide on the seam.
        seam = tick_step * 0.35
        ax_l.set_xlim(xmin - pad_l, -lo + seam)
        ax_m.set_xlim(-hi - seam, xmax + pad_r)
        axes = [ax_l, ax_m]
    else:
        ax_m = fig.add_subplot(1, 1, 1)
        ax_m.set_xlim(xmin - pad_l, xmax + pad_r)
        ax_l = None
        axes = [ax_m]

    # ---- bars (drawn on every axis; each clips to its own xlim) ------------
    for ax in axes:
        for _, r in plotted.iterrows():
            for side in ('left', 'right'):
                xv = r.get(f'{side}_x')
                if pd.isna(xv):
                    continue
                ax.barh(r['_y'], xv, left=0.0, height=0.62,
                        color=bar_color, edgecolor='none', zorder=2)
        ax.axvline(0, color='black', linewidth=1.2, zorder=4)
        ax.invert_yaxis()

    # ---- ticks, gridlines --------------------------------------------------
    def _ticks(lo_x, hi_x):
        first = np.ceil(lo_x / tick_step) * tick_step
        return np.arange(first, hi_x + 1e-9, tick_step)

    for ax in axes:
        lo_x, hi_x = ax.get_xlim()
        ticks = _ticks(lo_x, hi_x)
        ax.set_xticks(ticks)
        ax.set_xticklabels(
            ['1.0' if abs(t) < 1e-9 else f'$10^{{-{int(abs(t))}}}$' for t in ticks],
            fontsize=label_fontsize - 1)
        for t in ticks:
            if abs(t) > 1e-9:
                ax.axvline(t, color='#BBBBBB', linewidth=0.5,
                           linestyle='--', zorder=0)
        ax.set_yticks([])
        for spine in ('top', 'right', 'left'):
            ax.spines[spine].set_visible(False)
        ax.spines['bottom'].set_linewidth(1.0)

    # ---- term labels -------------------------------------------------------
    # Placed just outside the bar tip when there is room on that axis,
    # otherwise inside against the centre line (long bars would run off-figure).
    for _, r in plotted.iterrows():
        for side, ha_out, ha_in, sign in (('left', 'right', 'right', -1.0),
                                          ('right', 'left', 'left', 1.0)):
            xv = r.get(f'{side}_x')
            if pd.isna(xv):
                continue
            ax = ax_m
            if break_range and xv < -break_range[1]:
                ax = ax_l if xv < -break_range[0] else ax_m
            lo_x, hi_x = ax.get_xlim()
            pad = (hi_x - lo_x) * 0.015
            color = _label_color(r[f'{side}_term'], r.get(f'{side}_genes', ''))
            outside_x = xv - sign * pad
            fits = (outside_x > lo_x) if sign < 0 else (outside_x < hi_x)
            if fits and (ax is ax_m or break_range is None):
                ax.text(outside_x, r['_y'], _clean(r[f'{side}_term']),
                        ha=ha_out, va='center', fontsize=label_fontsize,
                        color=color, zorder=5)
            else:
                ax_m.text(-sign * pad, r['_y'], _clean(r[f'{side}_term']),
                          ha=ha_in, va='center', fontsize=label_fontsize,
                          color=color, zorder=5)

    # ---- break markers -----------------------------------------------------
    if break_range:
        d = 0.012
        for ax, at in ((ax_l, 1.0), (ax_m, 0.0)):
            ax.spines['bottom'].set_visible(True)
            kw = dict(transform=ax.transAxes, color='black',
                      clip_on=False, linewidth=1.0)
            ax.plot((at - d, at + d), (-d, +d), **kw)
            ax.plot((at - d, at + d), (-d * 3, d * 3), **kw)
        ax_l.spines['right'].set_visible(False)
        ax_m.spines['left'].set_visible(False)
        # slashes across any bar that crosses the break
        for _, r in plotted.iterrows():
            xv = r.get('left_x')
            if pd.isna(xv) or xv > -break_range[0]:
                continue
            for ax, frac in ((ax_l, 0.97), (ax_m, 0.015)):
                lo_x, hi_x = ax.get_xlim()
                xc = lo_x + (hi_x - lo_x) * frac
                w = (hi_x - lo_x) * 0.02
                for off in (-w, w):
                    ax.plot([xc + off - w * 0.6, xc + off + w * 0.6],
                            [r['_y'] + 0.30, r['_y'] - 0.30],
                            color='#333333', linewidth=1.0, zorder=3)

    # ---- headers, screen brackets, labels ---------------------------------
    lo_x, hi_x = ax_m.get_xlim()
    ax_m.text(lo_x + (0 - lo_x) / 2, -0.9, left_header, ha='center', va='bottom',
              fontsize=header_fontsize)
    ax_m.text((hi_x - 0) / 2, -0.9, right_header, ha='center', va='bottom',
              fontsize=header_fontsize)

    bracket_ax = ax_l if break_range else ax_m
    bl, _ = bracket_ax.get_xlim()
    for s, block in yticks_blocks.items():
        y0, y1 = min(block) - 0.45, max(block) + 0.45
        bracket_ax.plot([bl, bl], [y0, y1], color='black', linewidth=1.2,
                        clip_on=False, zorder=5)
        bracket_ax.text(bl - (0 - bl) * 0.02, (y0 + y1) / 2, s, ha='right',
                        va='center', fontsize=header_fontsize, clip_on=False)

    ax_m.set_xlabel(xlabel, fontsize=label_fontsize + 1)
    if title:
        fig.suptitle(title, fontsize=header_fontsize + 2)

    ax_m.set_ylim(max(ypos) + 0.7, -1.3)

    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    png = str(Path(output_file).with_suffix('.png'))
    if png != str(output_file):
        plt.savefig(png, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")
    return plotted


def plot_go_score_panels(
    panels: List[Dict],
    output_file: str,
    top_n: int = 6,
    term_col: str = 'Gene_Set',
    stat_col: str = 't_statistic',
    neglog_col: str = 'neglog10_FDR',
    genes_col: str = 'genes_in_set_and_screen',
    count_col: str = 'n_overlap',
    show_counts: bool = True,
    categories: Optional[pd.DataFrame] = None,
    colors: Optional[Dict[str, str]] = None,
    term_colors: Optional[Dict[str, str]] = None,
    term_labels: Optional[Dict[str, str]] = None,
    default_label_color: str = '#000000',
    min_overlap_frac: float = 0.25,
    min_dominance: float = 2.0,
    bar_color: str = '#D3D3D3',
    tick_step: float = 5.0,
    figsize: Tuple[float, float] = (12.0, 6.0),
    wspace: float = 0.30,
    xlabel: str = 'corrected p-value enrichment',
    title: Optional[str] = None,
    fontsize: Optional[int] = None,
    dpi: int = 300,
    strip_go_id: bool = True,
    wrap_width: int = 26,
) -> pd.DataFrame:
    """Multi-panel diverging GO bar chart for score-level enrichment (Fig S3E).

    Same visual grammar as plot_go_enrichment_diverging (Fig 1F) but laid out as
    several side-by-side panels that share a row grid, which is what S3E needs:
    one panel per score type (Degradation, Synthesis), each with its own x axis,
    and row groups (60S-enriched, 40S-enriched) shared down the figure.

    Kept as a separate function rather than generalising the Fig 1F one: this
    repo's convention is parallel duplicated blocks over shared plotting helpers,
    so each figure's output stays independently reproducible, and the Fig 1F code
    path is not disturbed.

    panels
        [{'title': 'Degradation Score', 'noun': 'degradation',
          'groups': [('60S-enriched GO Terms', path), ('40S-enriched GO Terms', path)]},
         ...]
        Every panel must have the same number of groups, and groups are matched
        row-for-row across panels.
    show_counts
        Append the gene-overlap count to each term label, e.g. "Protein
        UFMylation (5)" -- the draft labels carry these.
    categories / colors / min_dominance
        Term-label colouring by gene overlap, identical rule to Fig 1F: a term
        takes a category's colour only if that category's overlap fraction is at
        least `min_overlap_frac` AND at least `min_dominance` times the runner-up,
        so terms straddling two categories stay black.

    Returns the plotted table.
    """
    label_fontsize = fontsize or 8
    header_fs = (fontsize or 8) + 1

    n_groups = {len(p['groups']) for p in panels}
    if len(n_groups) != 1:
        raise ValueError('every panel must have the same number of groups')

    # ---- assemble ----------------------------------------------------------
    rows = []
    for pi, panel in enumerate(panels):
        for gi, (glabel, path) in enumerate(panel['groups']):
            df = pd.read_csv(path)
            for c in (term_col, stat_col, neglog_col):
                if c not in df.columns:
                    raise ValueError(f"{path}: missing column '{c}'")
            df = df.dropna(subset=[stat_col, neglog_col])
            dec = (df[df[stat_col] < 0].sort_values(neglog_col, ascending=False)
                   .head(top_n).reset_index(drop=True))
            inc = (df[df[stat_col] > 0].sort_values(neglog_col, ascending=False)
                   .head(top_n).reset_index(drop=True))
            for i in range(top_n):
                e = {'panel': pi, 'panel_title': panel['title'],
                     'group': gi, 'group_label': glabel, 'rank': i}
                for side, sub, sign in (('left', dec, -1.0), ('right', inc, 1.0)):
                    if i < len(sub):
                        r = sub.iloc[i]
                        e[f'{side}_term'] = r[term_col]
                        e[f'{side}_x'] = sign * float(r[neglog_col])
                        e[f'{side}_genes'] = r.get(genes_col, '')
                        e[f'{side}_n'] = r.get(count_col, None)
                rows.append(e)
    plotted = pd.DataFrame(rows)
    if plotted.empty:
        raise ValueError('no terms selected')

    # ---- label colour / text ----------------------------------------------
    cat_genes: Dict[str, set] = {}
    if categories is not None and len(categories):
        gcol = 'gene' if 'gene' in categories.columns else categories.columns[0]
        ccol = 'category' if 'category' in categories.columns else categories.columns[-1]
        for cat, sub in categories.groupby(ccol):
            cat_genes[str(cat)] = set(sub[gcol].astype(str).str.strip())

    def _colour(term, genes_field):
        if term_colors and term in term_colors:
            return term_colors[term]
        if not cat_genes or not colors:
            return default_label_color
        genes = {g.strip() for g in str(genes_field).split(',') if g.strip()}
        if not genes:
            return default_label_color
        fr = sorted(((len(genes & m) / len(genes), c)
                     for c, m in cat_genes.items() if c in colors), reverse=True)
        if not fr:
            return default_label_color
        best, cat = fr[0]
        runner = fr[1][0] if len(fr) > 1 else 0.0
        if best < min_overlap_frac or best < runner * min_dominance:
            return default_label_color
        return colors[cat]

    def _text(term, n):
        t = term_labels[term] if (term_labels and term in term_labels) else str(term)
        if term_labels is None or term not in term_labels:
            if strip_go_id:
                t = re.sub(r'\s*\(GO:\d+\)\s*$', '', t)
        if show_counts and n is not None and not pd.isna(n):
            t = f'{t} ({int(n)})'
        return '\n'.join(textwrap.wrap(t, wrap_width)) if wrap_width else t

    # ---- row grid (shared across panels) -----------------------------------
    ng = n_groups.pop()
    gap = 1.0
    ypos, blocks = {}, {}
    y = 0.0
    for gi in range(ng):
        block = []
        for i in range(top_n):
            ypos[(gi, i)] = y
            block.append(y)
            y += 1.0
        blocks[gi] = block
        y += gap
    plotted['_y'] = [ypos[(r.group, r.rank)] for r in plotted.itertuples()]

    fig, axes = plt.subplots(1, len(panels), figsize=figsize, sharey=True)
    if len(panels) == 1:
        axes = [axes]
    fig.subplots_adjust(wspace=wspace)

    for pi, panel in enumerate(panels):
        ax = axes[pi]
        sub = plotted[plotted['panel'] == pi]
        xs = [v for c in ('left_x', 'right_x') if c in sub for v in sub[c].dropna()]
        xmin, xmax = (min(xs), max(xs)) if xs else (-1, 1)
        pad = max(abs(xmin), abs(xmax)) * 0.04 + 0.4
        ax.set_xlim(xmin - pad, xmax + pad)

        for _, r in sub.iterrows():
            for side in ('left', 'right'):
                xv = r.get(f'{side}_x')
                if pd.isna(xv):
                    continue
                ax.barh(r['_y'], xv, left=0.0, height=0.62,
                        color=bar_color, edgecolor='none', zorder=2)
        ax.axvline(0, color='black', linewidth=1.2, zorder=4)

        lo_x, hi_x = ax.get_xlim()
        ticks = np.arange(np.ceil(lo_x / tick_step) * tick_step, hi_x + 1e-9, tick_step)
        ax.set_xticks(ticks)
        ax.set_xticklabels(['1.0' if abs(t) < 1e-9 else f'$10^{{-{int(abs(t))}}}$'
                            for t in ticks], fontsize=label_fontsize - 0.5)
        for t in ticks:
            if abs(t) > 1e-9:
                ax.axvline(t, color='#BBBBBB', linewidth=0.5, linestyle='--', zorder=0)
        ax.set_yticks([])
        for sp in ('top', 'right', 'left'):
            ax.spines[sp].set_visible(False)
        ax.set_xlabel(xlabel, fontsize=label_fontsize)

        # term labels: just outside the bar tip when it fits, else inside at the centre
        for _, r in sub.iterrows():
            for side, sign, ha_out, ha_in in (('left', -1.0, 'right', 'right'),
                                              ('right', 1.0, 'left', 'left')):
                xv = r.get(f'{side}_x')
                if pd.isna(xv):
                    continue
                pad_t = (hi_x - lo_x) * 0.015
                col = _colour(r[f'{side}_term'], r.get(f'{side}_genes', ''))
                txt = _text(r[f'{side}_term'], r.get(f'{side}_n'))
                out_x = xv - sign * pad_t
                fits = (out_x > lo_x) if sign < 0 else (out_x < hi_x)
                if fits:
                    ax.text(out_x, r['_y'], txt, ha=ha_out, va='center',
                            fontsize=label_fontsize, color=col, zorder=5,
                            clip_on=False)
                else:
                    # Bar runs to the panel edge: sit the label inside it,
                    # anchored just off the centre line. clip_on keeps it from
                    # spilling into the next panel.
                    ax.text(-sign * pad_t, r['_y'], txt, ha=ha_in, va='center',
                            fontsize=label_fontsize, color=col, zorder=5,
                            clip_on=True)

        noun = panel.get('noun', 'score')
        ax.text(lo_x + (0 - lo_x) / 2, -1.05, f'decreased {noun}\nupon knockdown',
                ha='center', va='bottom', fontsize=label_fontsize)
        ax.text((hi_x - 0) / 2, -1.05, f'increased {noun}\nupon knockdown',
                ha='center', va='bottom', fontsize=label_fontsize)
        ax.text(0, -2.05, panel['title'], ha='center', va='bottom',
                fontsize=header_fs + 1)

    # group brackets + labels on the leftmost panel only (rows are shared)
    ax0 = axes[0]
    bl, _ = ax0.get_xlim()
    for gi in range(ng):
        y0, y1 = min(blocks[gi]) - 0.45, max(blocks[gi]) + 0.45
        ax0.plot([bl, bl], [y0, y1], color='black', linewidth=1.2,
                 clip_on=False, zorder=5)
        glabel = plotted.loc[plotted['group'] == gi, 'group_label'].iloc[0]
        ax0.text(bl - (0 - bl) * 0.03, (y0 + y1) / 2, glabel, ha='right',
                 va='center', fontsize=header_fs, clip_on=False)

    axes[0].set_ylim(max(ypos.values()) + 0.7, -2.6)
    if title:
        fig.suptitle(title, fontsize=header_fs + 3)
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    png = str(Path(output_file).with_suffix('.png'))
    if png != str(output_file):
        plt.savefig(png, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")
    return plotted


def parse_figsize(figsize_str: Optional[str], default) -> Optional[Tuple[float, float]]:
    if figsize_str is None:
        return default
    try:
        parts = figsize_str.split(',')
        return (float(parts[0]), float(parts[1]))
    except (ValueError, IndexError):
        print(f"Warning: Invalid figsize '{figsize_str}', using default {default}")
        return default


def parse_lim(lim_str: Optional[str]) -> Optional[Tuple[float, float]]:
    """Parse axis limit string like '-3,3' into tuple."""
    if lim_str is None:
        return None
    try:
        parts = lim_str.split(',')
        return (float(parts[0]), float(parts[1]))
    except (ValueError, IndexError):
        print(f"Warning: Invalid limit '{lim_str}', ignoring")
        return None


def parse_ref_lines(ref_lines_str: Optional[str]) -> Optional[List[Tuple[float, float]]]:
    if ref_lines_str is None:
        return None
    try:
        lines = []
        for line in ref_lines_str.split(';'):
            parts = line.strip().split(',')
            lines.append((float(parts[0]), float(parts[1])))
        return lines
    except (ValueError, IndexError) as e:
        print(f"Warning: Invalid ref-lines '{ref_lines_str}': {e}")
        return None


def _add_colors_arg(p):
    p.add_argument(
        '--colors',
        help=(
            'CSV file for custom category colors. '
            'Two-column (category,color): maps category names to hex. '
            'One-column (color): ordered palette replacing defaults.'
        )
    )


def _add_scatter_ref_args(p, has_diagonal=True, has_ref_line=False):
    """Add shared axis/reference line arguments to a scatter subparser."""
    p.add_argument('--xlim', help='X-axis limits as "min,max" (e.g., "-3,3")')
    p.add_argument('--ylim', help='Y-axis limits as "min,max" (e.g., "-3,3")')
    if has_diagonal:
        p.add_argument('--no-diagonal', action='store_true', help='Hide diagonal y=x reference line')
        p.add_argument('--diagonal-style', default='--', help='Linestyle for diagonal (default: --)')
    if has_ref_line:
        p.add_argument('--no-ref-line', action='store_true', help='Hide slope/intercept reference line')
        p.add_argument('--ref-line-style', default='--', help='Linestyle for reference line (default: --)')
    p.add_argument('--no-axis-lines', action='store_true', help='Hide horizontal/vertical lines at 0')
    p.add_argument('--axis-line-style', default='-', help='Linestyle for axis lines (default: -)')


def main():
    parser = argparse.ArgumentParser(
        description='CRISPRi Screen Visualization Toolkit',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', help='Visualization type')

    # scatter-replicates
    scatter_p = subparsers.add_parser('scatter-replicates')
    scatter_p.add_argument('--counts', required=True)
    scatter_p.add_argument('--x', required=True)
    scatter_p.add_argument('--y', required=True)
    scatter_p.add_argument('-o', '--output', required=True)
    scatter_p.add_argument('--gene-col')
    scatter_p.add_argument('--highlight')
    scatter_p.add_argument('--categories')
    scatter_p.add_argument('--no-log', action='store_true')
    scatter_p.add_argument('--title')
    scatter_p.add_argument('--figsize')
    scatter_p.add_argument('--fontsize', type=int)
    scatter_p.add_argument('--legend-outside', action='store_true')
    scatter_p.add_argument('--ref-lines')
    scatter_p.add_argument('--no-correlation', action='store_true')
    scatter_p.add_argument('--no-legend', action='store_true')
    _add_scatter_ref_args(scatter_p, has_diagonal=True, has_ref_line=False)
    _add_colors_arg(scatter_p)

    # heatmap
    heatmap_p = subparsers.add_parser('heatmap')
    heatmap_p.add_argument('--counts', required=True)
    heatmap_p.add_argument('--samples', required=True)
    heatmap_p.add_argument('-o', '--output', required=True)
    heatmap_p.add_argument('--method', choices=['pearson', 'spearman'], default='pearson')
    heatmap_p.add_argument('--no-log', action='store_true')
    heatmap_p.add_argument('--no-cluster', action='store_true')
    heatmap_p.add_argument('--title')
    heatmap_p.add_argument('--order')
    heatmap_p.add_argument('--figsize')

    # heatmap-fc
    heatmap_fc_p = subparsers.add_parser('heatmap-fc')
    heatmap_fc_p.add_argument('--results', nargs='+', required=True)
    heatmap_fc_p.add_argument('--names')
    heatmap_fc_p.add_argument('-o', '--output', required=True)
    heatmap_fc_p.add_argument('--method', choices=['pearson', 'spearman'], default='pearson')
    heatmap_fc_p.add_argument('--no-cluster', action='store_true')
    heatmap_fc_p.add_argument('--title')
    heatmap_fc_p.add_argument('--order')
    heatmap_fc_p.add_argument('--figsize')
    heatmap_fc_p.add_argument('--merge-on', choices=['gene', 'gene_id'])

    # compare
    compare_p = subparsers.add_parser('compare')
    compare_p.add_argument('--results', nargs=2, required=True)
    compare_p.add_argument('--hits', nargs=2, required=True)
    compare_p.add_argument('--names')
    compare_p.add_argument('-o', '--output', required=True)
    compare_p.add_argument('--categories')
    compare_p.add_argument('--label')
    compare_p.add_argument('--label-top', type=int, default=0)
    compare_p.add_argument('--no-size', action='store_true')
    compare_p.add_argument('--no-ntc', action='store_true')
    compare_p.add_argument('--title')
    compare_p.add_argument('--hits-csv')
    compare_p.add_argument('--figsize')
    compare_p.add_argument('--fontsize', type=int)
    compare_p.add_argument('--legend-outside', action='store_true')
    compare_p.add_argument('--ref-lines')
    compare_p.add_argument('--merge-on', choices=['gene', 'gene_id'])
    compare_p.add_argument('--no-correlation', action='store_true')
    compare_p.add_argument('--no-legend', action='store_true')
    _add_scatter_ref_args(compare_p, has_diagonal=True, has_ref_line=False)
    _add_colors_arg(compare_p)

    # compare-matrix
    matrix_p = subparsers.add_parser('compare-matrix')
    matrix_p.add_argument('--results', nargs='+', required=True)
    matrix_p.add_argument('--hits', nargs='+', required=True)
    matrix_p.add_argument('--names')
    matrix_p.add_argument('-o', '--output', required=True)

    # volcano
    volcano_p = subparsers.add_parser('volcano')
    volcano_p.add_argument('--results', required=True)
    volcano_p.add_argument('--hits', required=True)
    volcano_p.add_argument('-o', '--output', required=True)
    volcano_p.add_argument('--categories')
    volcano_p.add_argument('--label')
    volcano_p.add_argument('--label-top', type=int, default=10)
    volcano_p.add_argument('--title')
    volcano_p.add_argument('--fdr', type=float, default=0.05)
    volcano_p.add_argument('--figsize')
    volcano_p.add_argument('--fontsize', type=int)
    volcano_p.add_argument('--legend-outside', action='store_true')
    volcano_p.add_argument('--no-legend', action='store_true')
    volcano_p.add_argument('--xlim', help='X-axis limits as "min,max"')
    volcano_p.add_argument('--ylim', help='Y-axis limits as "min,max"')
    _add_colors_arg(volcano_p)

    # waterfall
    waterfall_p = subparsers.add_parser('waterfall')
    waterfall_p.add_argument('--results', required=True)
    waterfall_p.add_argument('--hits')
    waterfall_p.add_argument('-o', '--output', required=True)
    waterfall_p.add_argument('--categories')
    waterfall_p.add_argument('--label')
    waterfall_p.add_argument('--label-top', type=int, default=0)
    waterfall_p.add_argument('--title')
    waterfall_p.add_argument('--no-ntc', action='store_true')
    waterfall_p.add_argument('--figsize')
    waterfall_p.add_argument('--fontsize', type=int)
    waterfall_p.add_argument('--legend-outside', action='store_true')
    waterfall_p.add_argument('--no-legend', action='store_true')
    waterfall_p.add_argument('--ylim', help='Y-axis limits as "min,max"')
    waterfall_p.add_argument('--value-col', help='Column to rank/plot (default: epsilon)')
    waterfall_p.add_argument('--ylabel', help='Y-axis label override')
    waterfall_p.add_argument('--ref-lines', default='0.585,-0.585',
                             help='Comma-separated dashed guide y-values, or "none" to disable')
    waterfall_p.add_argument('--descending', action='store_true',
                             help='Rank high-to-low instead of low-to-high')
    _add_colors_arg(waterfall_p)

    # waterfall-diff (per-gene difference between two score/effect tables)
    wdiff_p = subparsers.add_parser('waterfall-diff')
    wdiff_p.add_argument('--results', nargs=2, required=True,
                         metavar=('TABLE_A', 'TABLE_B'),
                         help='Two per-gene tables; plotted as A - B')
    wdiff_p.add_argument('-o', '--output', required=True)
    wdiff_p.add_argument('--value-col', default='score',
                         help='Column present in both tables (default: score)')
    wdiff_p.add_argument('--names', help='Comma-separated A,B labels')
    wdiff_p.add_argument('--merge-on', help='Gene id column to merge on')
    wdiff_p.add_argument('--diff-csv', help='Optional path to save the merged diff table')
    wdiff_p.add_argument('--categories')
    wdiff_p.add_argument('--label')
    wdiff_p.add_argument('--label-top', type=int, default=0)
    wdiff_p.add_argument('--title')
    wdiff_p.add_argument('--ylabel')
    wdiff_p.add_argument('--ref-lines', default='none',
                         help='Comma-separated dashed guide y-values (default: none)')
    wdiff_p.add_argument('--no-ntc', action='store_true')
    wdiff_p.add_argument('--figsize')
    wdiff_p.add_argument('--fontsize', type=int)
    wdiff_p.add_argument('--legend-outside', action='store_true')
    wdiff_p.add_argument('--no-legend', action='store_true')
    wdiff_p.add_argument('--ylim', help='Y-axis limits as "min,max"')
    _add_colors_arg(wdiff_p)

    # score
    score_p = subparsers.add_parser('score')
    score_p.add_argument('--results', nargs=2, required=True)
    score_p.add_argument('--hits', nargs=2, required=True)
    score_p.add_argument('--names')
    score_p.add_argument('--line', required=True)
    score_p.add_argument('-o', '--output', required=True)
    score_p.add_argument('--categories')

    # score-plot
    score_plot_p = subparsers.add_parser('score-plot')
    score_plot_p.add_argument('--results', nargs=2, required=True)
    score_plot_p.add_argument('--hits', nargs=2, required=True)
    score_plot_p.add_argument('--names')
    score_plot_p.add_argument('--line', required=True)
    score_plot_p.add_argument('-o', '--output', required=True)
    score_plot_p.add_argument('--score-csv')
    score_plot_p.add_argument('--categories')
    score_plot_p.add_argument('--label')
    score_plot_p.add_argument('--label-top', type=int, default=0)
    score_plot_p.add_argument('--label-bottom', type=int, default=0)
    score_plot_p.add_argument('--title')
    score_plot_p.add_argument('--figsize')
    score_plot_p.add_argument('--fontsize', type=int)
    score_plot_p.add_argument('--legend-outside', action='store_true')
    score_plot_p.add_argument('--cmap', default='coolwarm')
    score_plot_p.add_argument('--no-legend', action='store_true')
    _add_scatter_ref_args(score_plot_p, has_diagonal=False, has_ref_line=True)

    # sgrna-score-utest
    sgrna_score_p = subparsers.add_parser('sgrna-score-utest')
    sgrna_score_p.add_argument('--sgrna1', required=True)
    sgrna_score_p.add_argument('--sgrna2', required=True)
    sgrna_score_p.add_argument('--names')
    sgrna_score_p.add_argument('--line', required=True)
    sgrna_score_p.add_argument('-o', '--output', required=True)
    sgrna_score_p.add_argument('--top-n', type=int, default=3)
    sgrna_score_p.add_argument('--fdr', type=float, default=0.05)
    sgrna_score_p.add_argument('--min-count', type=int, default=20)
    sgrna_score_p.add_argument('--score-mode', choices=['perp_weighted', 'parallel'],
                               default='perp_weighted',
                               help="'perp_weighted' (default) downweights off-line guides; "
                                    "'parallel' uses raw distance along the line")

    # sgrna-score-plot
    sgrna_plot_p = subparsers.add_parser('sgrna-score-plot')
    sgrna_plot_p.add_argument('--sgrna1', required=True)
    sgrna_plot_p.add_argument('--sgrna2', required=True)
    sgrna_plot_p.add_argument('--names')
    sgrna_plot_p.add_argument('--line', required=True)
    sgrna_plot_p.add_argument('-o', '--output', required=True)
    sgrna_plot_p.add_argument('--min-count', type=int, default=20)
    sgrna_plot_p.add_argument('--label')
    sgrna_plot_p.add_argument('--label-top', type=int, default=0)
    sgrna_plot_p.add_argument('--label-bottom', type=int, default=0)
    sgrna_plot_p.add_argument('--highlight')
    sgrna_plot_p.add_argument('--title')
    sgrna_plot_p.add_argument('--figsize')
    sgrna_plot_p.add_argument('--fontsize', type=int)
    sgrna_plot_p.add_argument('--legend-outside', action='store_true')
    sgrna_plot_p.add_argument('--cmap', default='coolwarm')
    sgrna_plot_p.add_argument('--score-csv')
    sgrna_plot_p.add_argument('--vmin', type=float)
    sgrna_plot_p.add_argument('--vmax', type=float)
    sgrna_plot_p.add_argument('--no-legend', action='store_true')
    sgrna_plot_p.add_argument('--plain', action='store_true',
                               help='Plot all targeting sgRNAs in grey without score coloring')
    sgrna_plot_p.add_argument('--score-mode', choices=['perp_weighted', 'parallel'],
                               default='perp_weighted',
                               help="'perp_weighted' (default) downweights off-line guides; "
                                    "'parallel' uses raw distance along the line")
    sgrna_plot_p.add_argument('--ntc-color', default='#DDDDDD',
                               help='Color for NTC points (default: light grey #DDDDDD)')
    _add_scatter_ref_args(sgrna_plot_p, has_diagonal=False, has_ref_line=True)

    # sgrna-to-genes
    sgrna_to_genes_p = subparsers.add_parser('sgrna-to-genes')
    sgrna_to_genes_p.add_argument('--scores', required=True)
    sgrna_to_genes_p.add_argument('-o', '--output', required=True)
    sgrna_to_genes_p.add_argument('--top-n', type=int, default=3)
    sgrna_to_genes_p.add_argument('--fdr', type=float, default=0.05)

    # line-score-pipeline
    pipeline_p = subparsers.add_parser('line-score-pipeline')
    pipeline_p.add_argument('--sgrna1', required=True)
    pipeline_p.add_argument('--sgrna2', required=True)
    pipeline_p.add_argument('--results1', required=True)
    pipeline_p.add_argument('--results2', required=True)
    pipeline_p.add_argument('--names')
    pipeline_p.add_argument('--line', required=True)
    pipeline_p.add_argument('-o', '--output', required=True)
    pipeline_p.add_argument('--top-n', type=int, default=3)
    pipeline_p.add_argument('--fdr', type=float, default=0.05)
    pipeline_p.add_argument('--min-count', type=int, default=20)
    pipeline_p.add_argument('--label-top', type=int, default=10)
    pipeline_p.add_argument('--label-bottom', type=int, default=5)
    pipeline_p.add_argument('--figsize')
    pipeline_p.add_argument('--fontsize', type=int)
    pipeline_p.add_argument('--cmap', default='coolwarm')
    pipeline_p.add_argument('--png', action='store_true')
    pipeline_p.add_argument('--vmin', type=float)
    pipeline_p.add_argument('--vmax', type=float)

    # gene-score-plot
    gene_plot_p = subparsers.add_parser('gene-score-plot')
    gene_plot_p.add_argument('--scores', required=True)
    gene_plot_p.add_argument('--results1', required=True)
    gene_plot_p.add_argument('--results2', required=True)
    gene_plot_p.add_argument('--names')
    gene_plot_p.add_argument('-o', '--output', required=True)
    gene_plot_p.add_argument('--label')
    gene_plot_p.add_argument('--label-top', type=int, default=0)
    gene_plot_p.add_argument('--label-bottom', type=int, default=0)
    gene_plot_p.add_argument('--title')
    gene_plot_p.add_argument('--figsize')
    gene_plot_p.add_argument('--fontsize', type=int)
    gene_plot_p.add_argument('--legend-outside', action='store_true')
    gene_plot_p.add_argument('--cmap', default='coolwarm')
    gene_plot_p.add_argument('--merge-on', choices=['gene', 'gene_id'])
    gene_plot_p.add_argument('--vmin', type=float)
    gene_plot_p.add_argument('--vmax', type=float)
    gene_plot_p.add_argument('--no-legend', action='store_true')
    gene_plot_p.add_argument('--xlim', help='X-axis limits as "min,max"')
    gene_plot_p.add_argument('--ylim', help='Y-axis limits as "min,max"')
    gene_plot_p.add_argument('--no-axis-lines', action='store_true')
    gene_plot_p.add_argument('--axis-line-style', default='-')

    # compare-scores
    compare_scores_p = subparsers.add_parser('compare-scores')
    compare_scores_p.add_argument('--scores', nargs='+', required=True)
    compare_scores_p.add_argument('--names', required=True)
    compare_scores_p.add_argument('--x', required=True)
    compare_scores_p.add_argument('--y', required=True)
    compare_scores_p.add_argument('-o', '--output', required=True)
    compare_scores_p.add_argument('--gene-col', default='gene_id')
    compare_scores_p.add_argument('--score-col', default='epsilon')
    compare_scores_p.add_argument('--color-by')
    compare_scores_p.add_argument('--color-col', default='product')
    compare_scores_p.add_argument('--categories')
    compare_scores_p.add_argument('--label')
    compare_scores_p.add_argument('--label-top', type=int, default=0)
    compare_scores_p.add_argument('--title')
    compare_scores_p.add_argument('--figsize')
    compare_scores_p.add_argument('--fontsize', type=int)
    compare_scores_p.add_argument('--cmap', default='coolwarm')
    compare_scores_p.add_argument('--no-correlation', action='store_true')
    compare_scores_p.add_argument('--no-legend', action='store_true')
    _add_scatter_ref_args(compare_scores_p, has_diagonal=True, has_ref_line=False)
    _add_colors_arg(compare_scores_p)

    # strip-plot
    strip_p = subparsers.add_parser('strip-plot')
    strip_p.add_argument('--scores')
    strip_p.add_argument('--sgrna1')
    strip_p.add_argument('--sgrna2')
    strip_p.add_argument('--line')
    strip_p.add_argument('--genes', required=True)
    strip_p.add_argument('-o', '--output', required=True)
    strip_p.add_argument('--score-col', default='score')
    strip_p.add_argument('--gene-col', default='gene')
    strip_p.add_argument('--categories')
    strip_p.add_argument('--title')
    strip_p.add_argument('--figsize')
    strip_p.add_argument('--fontsize', type=int)
    strip_p.add_argument('--no-median', action='store_true')
    strip_p.add_argument('--no-ntc', action='store_true')
    strip_p.add_argument('--min-count', type=int, default=20)
    strip_p.add_argument('--no-legend', action='store_true')
    strip_p.add_argument('--ylim', help='Y-axis limits as "min,max"')
    strip_p.add_argument('--no-axis-lines', action='store_true')
    strip_p.add_argument('--axis-line-style', default='--')
    _add_colors_arg(strip_p)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return

    def get_colors(a):
        return load_colors(a.colors) if getattr(a, 'colors', None) else None

    if args.command == 'scatter-replicates':
        df = load_counts(args.counts)
        categories = load_categories(args.categories) if args.categories else None
        figsize = parse_figsize(args.figsize, (8, 8))
        plot_replicate_scatter(
            df, args.x, args.y, args.output,
            gene_col=args.gene_col,
            highlight_genes=args.highlight.split(',') if args.highlight else None,
            categories=categories, log_scale=not args.no_log,
            title=args.title, figsize=figsize, fontsize=args.fontsize,
            legend_outside=args.legend_outside,
            ref_lines=parse_ref_lines(getattr(args, 'ref_lines', None)),
            show_correlation=not args.no_correlation,
            show_legend=not args.no_legend,
            colors=get_colors(args),
            xlim=parse_lim(args.xlim),
            ylim=parse_lim(args.ylim),
            show_diagonal=not args.no_diagonal,
            show_axis_lines=not args.no_axis_lines,
            diagonal_style=args.diagonal_style,
            axis_line_style=args.axis_line_style
        )

    elif args.command == 'heatmap':
        df = load_counts(args.counts)
        samples = [s.strip() for s in args.samples.split(',')]
        sample_order = [s.strip() for s in args.order.split(',')] if args.order else None
        figsize = parse_figsize(args.figsize, (10, 8))
        plot_correlation_heatmap(
            df, samples, args.output, method=args.method,
            log_transform=not args.no_log, cluster=not args.no_cluster,
            sample_order=sample_order, title=args.title, figsize=figsize
        )

    elif args.command == 'heatmap-fc':
        names = [s.strip() for s in args.names.split(',')] if args.names else None
        sample_order = [s.strip() for s in args.order.split(',')] if args.order else None
        figsize = parse_figsize(args.figsize, (10, 8))
        plot_foldchange_heatmap(
            args.results, args.output, comparison_names=names,
            method=args.method, cluster=not args.no_cluster,
            sample_order=sample_order, title=args.title, figsize=figsize,
            merge_on=getattr(args, 'merge_on', None)
        )

    elif args.command == 'compare':
        names = [s.strip() for s in args.names.split(',')] if args.names else [None, None]
        comp1 = load_comparison(args.results[0], args.hits[0], names[0])
        comp2 = load_comparison(args.results[1], args.hits[1], names[1] if len(names) > 1 else None)
        categories = load_categories(args.categories) if args.categories else None
        label_genes = [g.strip() for g in args.label.split(',')] if args.label else None
        figsize = parse_figsize(args.figsize, (10, 10))
        plot_comparison_scatter(
            comp1, comp2, args.output, categories=categories,
            label_genes=label_genes, label_top_n=args.label_top,
            size_by_hits=not args.no_size, show_ntc=not args.no_ntc,
            title=args.title, hits_csv=args.hits_csv, figsize=figsize,
            fontsize=args.fontsize, legend_outside=args.legend_outside,
            ref_lines=parse_ref_lines(getattr(args, 'ref_lines', None)),
            merge_on=getattr(args, 'merge_on', None),
            show_correlation=not args.no_correlation,
            show_legend=not args.no_legend,
            colors=get_colors(args),
            xlim=parse_lim(args.xlim),
            ylim=parse_lim(args.ylim),
            show_diagonal=not args.no_diagonal,
            show_axis_lines=not args.no_axis_lines,
            diagonal_style=args.diagonal_style,
            axis_line_style=args.axis_line_style
        )

    elif args.command == 'compare-matrix':
        names = [s.strip() for s in args.names.split(',')] if args.names else [None] * len(args.results)
        comparisons = [load_comparison(r, h, names[i] if i < len(names) else None)
                       for i, (r, h) in enumerate(zip(args.results, args.hits))]
        plot_multi_comparison_matrix(comparisons, args.output)

    elif args.command == 'volcano':
        categories = load_categories(args.categories) if args.categories else None
        label_genes = [g.strip() for g in args.label.split(',')] if args.label else None
        figsize = parse_figsize(args.figsize, (10, 8))
        plot_volcano(
            args.results, args.hits, args.output, categories=categories,
            label_genes=label_genes, label_top_n=args.label_top,
            title=args.title, fdr=args.fdr, figsize=figsize,
            fontsize=args.fontsize, legend_outside=args.legend_outside,
            show_legend=not args.no_legend,
            colors=get_colors(args),
            xlim=parse_lim(args.xlim),
            ylim=parse_lim(args.ylim)
        )

    elif args.command == 'waterfall':
        categories = load_categories(args.categories) if args.categories else None
        label_genes = [g.strip() for g in args.label.split(',')] if args.label else None
        figsize = parse_figsize(args.figsize, (12, 6))
        rl = args.ref_lines
        ref_lines = None if (rl is None or rl.strip().lower() in ('none', '')) \
            else tuple(float(v) for v in rl.split(','))
        plot_waterfall(
            args.results, args.output, hits_file=args.hits,
            categories=categories, label_genes=label_genes,
            label_top_n=args.label_top, title=args.title,
            value_col=args.value_col, ylabel=args.ylabel,
            ref_lines=ref_lines, ascending=not args.descending,
            show_ntc=not args.no_ntc, figsize=figsize,
            fontsize=args.fontsize, legend_outside=args.legend_outside,
            show_legend=not args.no_legend,
            colors=get_colors(args),
            ylim=parse_lim(args.ylim)
        )

    elif args.command == 'waterfall-diff':
        categories = load_categories(args.categories) if args.categories else None
        label_genes = [g.strip() for g in args.label.split(',')] if args.label else None
        figsize = parse_figsize(args.figsize, (12, 6))
        names = [s.strip() for s in args.names.split(',')] if args.names else ['A', 'B']
        rl = args.ref_lines
        ref_lines = None if (rl is None or rl.strip().lower() in ('none', '')) \
            else tuple(float(v) for v in rl.split(','))
        plot_waterfall_diff(
            args.results[0], args.results[1], args.output,
            value_col=args.value_col, name_a=names[0],
            name_b=names[1] if len(names) > 1 else 'B',
            merge_on=args.merge_on, categories=categories,
            label_genes=label_genes, label_top_n=args.label_top,
            title=args.title, ylabel=args.ylabel, ref_lines=ref_lines,
            show_ntc=not args.no_ntc, figsize=figsize,
            fontsize=args.fontsize, legend_outside=args.legend_outside,
            show_legend=not args.no_legend, colors=get_colors(args),
            ylim=parse_lim(args.ylim), diff_csv=args.diff_csv
        )

    elif args.command == 'score':
        names = [s.strip() for s in args.names.split(',')] if args.names else [None, None]
        comp1 = load_comparison(args.results[0], args.hits[0], names[0])
        comp2 = load_comparison(args.results[1], args.hits[1], names[1] if len(names) > 1 else None)
        slope, intercept = (float(v) for v in args.line.split(','))
        score_df = calculate_line_score(comp1, comp2, slope, intercept)
        if args.categories:
            cats = load_categories(args.categories)
            g2c = dict(zip(cats.iloc[:, 0].astype(str), cats.iloc[:, 1].astype(str)))
            score_df['category'] = score_df['gene'].astype(str).map(lambda g: g2c.get(g, ''))
        score_df.to_csv(args.output, index=False)
        print(f"Saved scores: {args.output} ({len(score_df)} genes)")
        print(f"Score range: {score_df['score'].min():.3f} to {score_df['score'].max():.3f}")

    elif args.command == 'score-plot':
        names = [s.strip() for s in args.names.split(',')] if args.names else [None, None]
        comp1 = load_comparison(args.results[0], args.hits[0], names[0])
        comp2 = load_comparison(args.results[1], args.hits[1], names[1] if len(names) > 1 else None)
        slope, intercept = (float(v) for v in args.line.split(','))
        score_df = calculate_line_score(comp1, comp2, slope, intercept)
        categories = None
        if args.categories:
            categories = load_categories(args.categories)
            g2c = dict(zip(categories.iloc[:, 0].astype(str), categories.iloc[:, 1].astype(str)))
            score_df['category'] = score_df['gene'].astype(str).map(lambda g: g2c.get(g, ''))
        if args.score_csv:
            score_df.to_csv(args.score_csv, index=False)
            print(f"Saved scores: {args.score_csv} ({len(score_df)} genes)")
        label_genes = [g.strip() for g in args.label.split(',')] if args.label else None
        figsize = parse_figsize(args.figsize, (10, 10))
        plot_score_scatter(
            score_df, args.output, comp1_name=comp1.name, comp2_name=comp2.name,
            slope=slope, intercept=intercept, categories=categories,
            label_genes=label_genes, label_top_n=args.label_top,
            label_bottom_n=args.label_bottom, title=args.title,
            figsize=figsize, fontsize=args.fontsize,
            legend_outside=args.legend_outside, cmap=args.cmap,
            show_legend=not args.no_legend,
            plain=getattr(args, 'plain', False),
            xlim=parse_lim(args.xlim),
            ylim=parse_lim(args.ylim),
            show_ref_line=not args.no_ref_line,
            show_axis_lines=not args.no_axis_lines,
            ref_line_style=args.ref_line_style,
            axis_line_style=args.axis_line_style
        )

    elif args.command == 'sgrna-score-utest':
        names = [s.strip() for s in args.names.split(',')] if args.names else [None, None]
        slope, intercept = (float(v) for v in args.line.split(','))
        sgrna_line_score_utest(
            args.sgrna1, args.sgrna2, slope, intercept, args.output,
            comp1_name=names[0], comp2_name=names[1] if len(names) > 1 else None,
            top_n=args.top_n, fdr=args.fdr, min_count=args.min_count,
            score_mode=args.score_mode
        )

    elif args.command == 'sgrna-score-plot':
        names = [s.strip() for s in args.names.split(',')] if args.names else [None, None]
        slope, intercept = (float(v) for v in args.line.split(','))
        label_genes = [g.strip() for g in args.label.split(',')] if args.label else None
        highlight_genes = [g.strip() for g in args.highlight.split(',')] if args.highlight else None
        figsize = parse_figsize(args.figsize, (10, 10))
        plot_sgrna_score_scatter(
            args.sgrna1, args.sgrna2, slope, intercept, args.output,
            comp1_name=names[0], comp2_name=names[1] if len(names) > 1 else None,
            min_count=args.min_count, label_genes=label_genes,
            label_top_n=args.label_top, label_bottom_n=args.label_bottom,
            highlight_genes=highlight_genes, title=args.title,
            figsize=figsize, fontsize=args.fontsize,
            legend_outside=args.legend_outside, cmap=args.cmap,
            score_csv=args.score_csv, vmin=args.vmin, vmax=args.vmax,
            show_legend=not args.no_legend,
            xlim=parse_lim(args.xlim),
            ylim=parse_lim(args.ylim),
            show_ref_line=not args.no_ref_line,
            show_axis_lines=not args.no_axis_lines,
            ref_line_style=args.ref_line_style,
            axis_line_style=args.axis_line_style,
            score_mode=args.score_mode,
            ntc_color=args.ntc_color
        )

    elif args.command == 'sgrna-to-genes':
        sgrna_scores_to_genes(args.scores, args.output, top_n=args.top_n, fdr=args.fdr)

    elif args.command == 'line-score-pipeline':
        names = [s.strip() for s in args.names.split(',')] if args.names else [None, None]
        slope, intercept = (float(v) for v in args.line.split(','))
        figsize = parse_figsize(args.figsize, (10, 10))
        line_score_full_pipeline(
            args.sgrna1, args.sgrna2, args.results1, args.results2,
            slope, intercept, args.output,
            comp1_name=names[0], comp2_name=names[1] if len(names) > 1 else None,
            top_n=args.top_n, fdr=args.fdr, min_count=args.min_count,
            label_top_n=args.label_top, label_bottom_n=args.label_bottom,
            figsize=figsize, fontsize=args.fontsize, cmap=args.cmap,
            also_png=args.png, vmin=args.vmin, vmax=args.vmax
        )

    elif args.command == 'gene-score-plot':
        names = [s.strip() for s in args.names.split(',')] if args.names else [None, None]
        label_genes = [g.strip() for g in args.label.split(',')] if args.label else None
        figsize = parse_figsize(args.figsize, (10, 10))
        plot_gene_line_scores(
            args.scores, args.results1, args.results2, args.output,
            comp1_name=names[0], comp2_name=names[1] if len(names) > 1 else None,
            label_genes=label_genes, label_top_n=args.label_top,
            label_bottom_n=args.label_bottom, title=args.title,
            figsize=figsize, fontsize=args.fontsize,
            legend_outside=args.legend_outside, cmap=args.cmap,
            merge_on=getattr(args, 'merge_on', None),
            vmin=args.vmin, vmax=args.vmax,
            show_legend=not args.no_legend,
            xlim=parse_lim(args.xlim),
            ylim=parse_lim(args.ylim),
            show_axis_lines=not args.no_axis_lines,
            axis_line_style=args.axis_line_style
        )

    elif args.command == 'compare-scores':
        score_names = [s.strip() for s in args.names.split(',')]
        label_genes = [g.strip() for g in args.label.split(',')] if args.label else None
        figsize = parse_figsize(args.figsize, (10, 10))
        categories = load_categories(args.categories) if args.categories else None
        plot_score_comparison(
            args.scores, score_names, args.x, args.y, args.output,
            gene_col=args.gene_col, score_col=args.score_col,
            color_by=args.color_by, color_score_col=args.color_col,
            categories=categories, label_genes=label_genes,
            label_top_n=args.label_top, title=args.title,
            figsize=figsize, fontsize=args.fontsize, cmap=args.cmap,
            show_correlation=not args.no_correlation,
            show_legend=not args.no_legend,
            colors=get_colors(args),
            xlim=parse_lim(args.xlim),
            ylim=parse_lim(args.ylim),
            show_diagonal=not args.no_diagonal,
            show_axis_lines=not args.no_axis_lines,
            diagonal_style=args.diagonal_style,
            axis_line_style=args.axis_line_style
        )

    elif args.command == 'strip-plot':
        genes = [g.strip() for g in args.genes.split(',')]
        categories = load_categories(args.categories) if args.categories else None
        figsize = parse_figsize(args.figsize, None)
        if args.scores:
            sgrna_data = load_sgrna_scores_for_strip(args.scores)
        elif args.sgrna1 and args.sgrna2 and args.line:
            slope, intercept = (float(v) for v in args.line.split(','))
            sgrna_data = load_sgrna_scores_for_strip(
                args.sgrna1, args.sgrna2, slope=slope, intercept=intercept,
                min_count=args.min_count
            )
        else:
            print("Error: Must provide either --scores or (--sgrna1, --sgrna2, --line)")
            return
        plot_sgrna_strip(
            sgrna_data, genes, args.output,
            score_col=args.score_col, gene_col=args.gene_col,
            categories=categories, title=args.title,
            figsize=figsize, fontsize=args.fontsize,
            show_median=not args.no_median, show_ntc=not args.no_ntc,
            show_legend=not args.no_legend,
            colors=get_colors(args),
            ylim=parse_lim(args.ylim),
            show_axis_lines=not args.no_axis_lines,
            axis_line_style=args.axis_line_style
        )


if __name__ == '__main__':
    main()
