"""Waterfall plotting library for MAGeCK+MWU screen results.

This module houses the plotting code so it can be used from both the CLI
(03_waterfall.py) and notebooks. Import from a notebook via:

    import sys; sys.path.insert(0, '/path/to/ricin_pipeline')
    import waterfall

    fig, ax = waterfall.plot_waterfall(
        gene_scores='results/K562_vs_Ricin_all_genes.csv',
        rank_by='epsilon',
        background_known=True,
        highlights=[
            dict(label='R3HDM4', color='#ffcc00', size=450, marker='*',
                 genes=['R3HDM4'], annotate=True),
            dict(label='60S RPs', color='#9467bd', size=45,
                 genes='@rpl_genes.txt'),
        ],
    )

See the module docstring in 03_waterfall.py for the highlight-spec DSL used
by the CLI. From Python, pass highlights as a list of dicts with the same
keys.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------------
# Literature-curated background gene sets (for --background-known)
# ---------------------------------------------------------------------------
KNOWN_SENSITIZERS = {
    # mTORC1 negative regulators
    "NPRL2", "NPRL3", "DEPDC5", "TSC1", "TSC2", "SZT2", "KPTN", "ITFG2", "C12orf66",
    # WDR11/FAM91A1 trafficking
    "WDR11", "FAM91A1", "C17orf75",
    # AP-1 (TGN clathrin adaptor)
    "AP1G1", "AP1B1", "AP1M1", "AP1S1",
    # Cohesin / chromatin
    "STAG2", "SMC1A", "SMC3", "RAD21", "ADNP", "POGZ", "HDAC8", "CBX3",
    # COPI coatomer
    "COPA", "COPB1", "COPB2", "COPG1", "ARCN1",
    # TM9SF / SLC39
    "TM9SF1", "TM9SF4", "SLC39A9", "SLC39A3",
    # Validated Bassik 2013 / Tian 2018 hits
    "SPI1", "TIPRL", "ASTE1", "IPO5", "UXS1", "MGAT3", "ST3GAL2",
}

KNOWN_RESISTANCE = {
    # N-glycosylation
    "MGAT1", "MGAT2", "SLC35A2", "GMDS", "TSTA3", "GALE", "SLC35C1", "UGP2", "HSP90B1",
    "ALG8", "ALG5", "ALG6", "PMM2",
    # ER export / COPII
    "SEC23B", "SEC24A", "SEC24C", "SEC31A", "SEC23A",
    # TRAPP
    "TRAPPC4", "TRAPPC5", "TRAPPC8", "TRAPPC11", "TRAPPC12", "TRAPPC13", "TRAPPC2L",
    # GARP
    "VPS53", "VPS54",
    # M6P
    "IGF2R", "NAGPA",
    # OST
    "STT3A", "OST4",
    # AP-3
    "AP3B1", "AP3D1",
    # Other
    "PIK3C3", "RABGEF1", "DHRSX", "ERP44", "TARDBP",
}


# ---------------------------------------------------------------------------
# Highlight parsing (used by both CLI DSL and notebook dicts)
# ---------------------------------------------------------------------------
_HIGHLIGHT_DEFAULTS = {
    "label": None,
    "color": "#ff7f0e",
    "size": 30.0,
    "marker": "o",
    "alpha": 0.9,
    "edgecolor": "white",
    "edgewidth": 0.4,
    "zorder": None,
    "genes": None,
    "annotate": False,
}


def parse_gene_list(spec: Union[str, Sequence[str]]) -> List[str]:
    """Normalize a gene spec to a list of gene symbols.

    Accepts:
      - a list/tuple of strings: returned as-is
      - "g1,g2,g3": comma-separated inline
      - "@path/to/file.txt": one gene per line (# comments allowed)
    """
    if isinstance(spec, (list, tuple, set)):
        return [str(g) for g in spec]
    if not isinstance(spec, str):
        raise TypeError(f"gene spec must be str or sequence, got {type(spec).__name__}")
    if spec.startswith("@"):
        path = spec[1:]
        if not os.path.isfile(path):
            raise FileNotFoundError(f"gene file not found: {path}")
        genes: List[str] = []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                for tok in line.replace(",", " ").split():
                    if tok:
                        genes.append(tok)
        return genes
    return [g.strip() for g in spec.split(",") if g.strip()]


def parse_highlight_dsl(spec: str) -> Dict:
    """Parse one CLI-style highlight string, e.g.
        'label=R3HDM4;color=#ffcc00;size=400;genes=R3HDM4;annotate=true'
    into the internal dict form used by make_waterfall.
    """
    opts = dict(_HIGHLIGHT_DEFAULTS)
    for piece in spec.split(";"):
        piece = piece.strip()
        if not piece:
            continue
        if "=" not in piece:
            raise ValueError(f"highlight spec '{spec}' has token without '=': '{piece}'")
        k, v = piece.split("=", 1)
        k = k.strip()
        v = v.strip()
        if k not in opts:
            raise ValueError(f"unknown highlight key '{k}' (valid: {sorted(opts)})")
        if k in ("size", "edgewidth", "alpha"):
            opts[k] = float(v)
        elif k == "zorder":
            opts[k] = int(v)
        elif k == "annotate":
            opts[k] = v.lower() in ("true", "yes", "1", "t", "y")
        elif k == "genes":
            opts[k] = parse_gene_list(v)
        else:
            opts[k] = v
    if opts["label"] is None:
        raise ValueError(f"highlight spec '{spec}' missing required key 'label'")
    if not opts["genes"]:
        raise ValueError(f"highlight spec '{spec}' has no genes")
    return opts


def _normalize_highlight(h: Union[Dict, "Highlight"]) -> Dict:
    """Coerce a Highlight dataclass or partial dict into a filled-in dict."""
    if isinstance(h, Highlight):
        h = h.as_dict()
    if not isinstance(h, dict):
        raise TypeError(f"highlight must be dict or Highlight, got {type(h).__name__}")
    if "label" not in h or h["label"] is None:
        raise ValueError(f"highlight missing 'label': {h}")
    if "genes" not in h or not h["genes"]:
        raise ValueError(f"highlight '{h.get('label')}' has no genes")
    out = dict(_HIGHLIGHT_DEFAULTS)
    out.update(h)
    out["genes"] = parse_gene_list(out["genes"])
    return out


@dataclass
class Highlight:
    """Convenience dataclass for building highlight groups from Python.

    Gives you autocomplete on the field names in editors/notebooks. Pass
    a list of Highlight instances (or plain dicts) to plot_waterfall via
    the ``highlights=`` argument.
    """
    label: str
    genes: Union[str, Sequence[str]]
    color: str = "#ff7f0e"
    size: float = 30.0
    marker: str = "o"
    alpha: float = 0.9
    edgecolor: str = "white"
    edgewidth: float = 0.4
    zorder: Optional[int] = None
    annotate: bool = False

    def as_dict(self) -> Dict:
        return dict(
            label=self.label, genes=self.genes, color=self.color, size=self.size,
            marker=self.marker, alpha=self.alpha, edgecolor=self.edgecolor,
            edgewidth=self.edgewidth, zorder=self.zorder, annotate=self.annotate,
        )


# ---------------------------------------------------------------------------
# Core plotting function
# ---------------------------------------------------------------------------
def make_waterfall(
    df: pd.DataFrame,
    rank_col: str = "epsilon",
    highlights: Optional[Sequence[Union[Dict, Highlight]]] = None,
    output: Optional[str] = None,
    background_known: bool = False,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (13.5, 7.5),
    annotate_fontsize: float = 8.0,
) -> Tuple[Figure, Axes]:
    """Draw the waterfall plot and return the matplotlib figure & axes.

    Parameters
    ----------
    df :
        Gene-level score table (e.g. from 02_gene_scores.py). Must contain
        columns 'gene' plus whatever is given as rank_col.
        NTC pseudo-genes (gene starting with 'NTC_') are dropped
        automatically.
    rank_col :
        Column name to rank genes by (default 'epsilon').
    highlights :
        List of highlight groups. Each may be a Highlight dataclass or a
        dict with the same keys.
    output :
        If given, save PDF (at this path) and PNG (same path with .png ext).
        If None, no files are written.
    background_known :
        If True, colour all literature-curated ricin genes as a backdrop
        layer (blue = sensitizer set, red = resistance set) using
        KNOWN_SENSITIZERS and KNOWN_RESISTANCE.
    title :
        Plot title. If None, an auto title is generated.
    figsize :
        Figure size in inches.
    annotate_fontsize :
        Font size for gene labels next to highlighted points (when
        annotate=True on a highlight).

    Returns
    -------
    (fig, ax) so you can further customize before saving/displaying.
    """
    if rank_col not in df.columns:
        raise ValueError(f"rank column '{rank_col}' not in df; have: {list(df.columns)}")

    df = df[~df["gene"].astype(str).str.startswith("NTC_")].copy()
    df = df.sort_values(rank_col, ascending=True).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    n = len(df)

    highlights_norm: List[Dict] = [_normalize_highlight(h) for h in (highlights or [])]

    fig, ax = plt.subplots(figsize=figsize)

    # Background layer
    if background_known:
        colors_bg = []
        for g in df["gene"]:
            if g in KNOWN_SENSITIZERS:
                colors_bg.append("#1f4f73")
            elif g in KNOWN_RESISTANCE:
                colors_bg.append("#8a1f1f")
            else:
                colors_bg.append("#bbbbbb")
        ax.scatter(df["rank"], df[rank_col],
                   s=3, c=colors_bg, alpha=0.55, linewidth=0, zorder=2)
        mask_known = df["gene"].isin(KNOWN_SENSITIZERS | KNOWN_RESISTANCE)
        ax.scatter(df.loc[mask_known, "rank"], df.loc[mask_known, rank_col],
                   s=18,
                   c=[c for c, m in zip(colors_bg, mask_known) if m],
                   alpha=0.9, edgecolor="white", linewidth=0.3, zorder=3)
    else:
        ax.scatter(df["rank"], df[rank_col],
                   s=3, c="#bbbbbb", alpha=0.55, linewidth=0, zorder=2)

    ax.axhline(0, color="black", linewidth=0.6, alpha=0.5)

    # Highlight layers
    default_zorder = 10
    legend_handles: List[Line2D] = []
    if background_known:
        legend_handles += [
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="#1f4f73", markersize=8,
                   label=f"Known ricin sensitizers (n={sum(g in KNOWN_SENSITIZERS for g in df['gene'])})"),
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="#8a1f1f", markersize=8,
                   label=f"Known ricin resistance (n={sum(g in KNOWN_RESISTANCE for g in df['gene'])})"),
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="#bbbbbb", markersize=5,
                   label="Other genes"),
        ]

    for i, h in enumerate(highlights_norm):
        sub = df[df["gene"].isin(h["genes"])]
        missing = set(h["genes"]) - set(df["gene"])
        if missing:
            print(f"[waterfall] '{h['label']}': {len(missing)} gene(s) not in score table, "
                  f"skipped: {sorted(missing)[:5]}{'...' if len(missing) > 5 else ''}",
                  file=sys.stderr)
        if len(sub) == 0:
            print(f"[waterfall] '{h['label']}': no matching genes, skipping", file=sys.stderr)
            continue
        z = h["zorder"] if h["zorder"] is not None else (default_zorder + i)
        ax.scatter(sub["rank"], sub[rank_col],
                   s=h["size"], c=h["color"], marker=h["marker"],
                   alpha=h["alpha"],
                   edgecolor=h["edgecolor"], linewidth=h["edgewidth"],
                   zorder=z)
        legend_handles.append(
            Line2D([0], [0], marker=h["marker"], color="w",
                   markerfacecolor=h["color"], markeredgecolor=h["edgecolor"],
                   markersize=max(8, np.sqrt(h["size"]) * 0.7),
                   label=f"{h['label']} (n={len(sub)})")
        )
        if h["annotate"]:
            for _, row in sub.iterrows():
                ax.annotate(
                    row["gene"],
                    xy=(row["rank"], row[rank_col]),
                    xytext=(row["rank"] + n * 0.012, row[rank_col]),
                    fontsize=annotate_fontsize, color=h["color"], fontweight="bold",
                    ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", color=h["color"], lw=0.5, alpha=0.6),
                    zorder=z + 1,
                )

    # Direction labels
    ax.text(0.015, 0.97, "← Sensitizers\n(depleted in ricin)",
            transform=ax.transAxes, fontsize=11, color="#1f4f73",
            ha="left", va="top", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#1f4f73", alpha=0.9))
    ax.text(0.985, 0.03, "Resistance →\n(enriched in ricin)",
            transform=ax.transAxes, fontsize=11, color="#8a1f1f",
            ha="right", va="bottom", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#8a1f1f", alpha=0.9))

    ax.set_xlim(-n * 0.025, n * 1.025)
    ax.set_xlabel(f"Gene rank by {rank_col}  (1 = strongest sensitizer  →  {n:,} = strongest resistance)",
                  fontsize=12)
    ylabel = {"epsilon": r"ε (mean LFC of top-N sgRNAs, NTC-centered)",
              "product": r"Product score (ε × −log$_{10}$ p)"}.get(rank_col, rank_col)
    ax.set_ylabel(ylabel, fontsize=12)
    if title is not None:
        ax.set_title(title, fontsize=13)
    else:
        ax.set_title(f"Waterfall plot — genes ranked by {rank_col}", fontsize=13)
    ax.grid(True, alpha=0.2)
    if legend_handles:
        ax.legend(handles=legend_handles, loc="center", fontsize=9, framealpha=0.95)

    fig.tight_layout()

    if output is not None:
        base, _ = os.path.splitext(output)
        fig.savefig(base + ".pdf", bbox_inches="tight")
        fig.savefig(base + ".png", dpi=160, bbox_inches="tight")
        print(f"[waterfall] wrote {base}.pdf and {base}.png", file=sys.stderr)

    return fig, ax


def plot_waterfall(
    gene_scores: Union[str, pd.DataFrame],
    rank_by: str = "epsilon",
    highlights: Optional[Sequence[Union[Dict, Highlight]]] = None,
    output: Optional[str] = None,
    background_known: bool = False,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (13.5, 7.5),
    annotate_fontsize: float = 8.0,
) -> Tuple[Figure, Axes]:
    """Convenience wrapper: accepts a CSV path OR a DataFrame.

    Same as make_waterfall but reads the gene-scores CSV from disk if
    given a string path.
    """
    if isinstance(gene_scores, str):
        df = pd.read_csv(gene_scores)
    elif isinstance(gene_scores, pd.DataFrame):
        df = gene_scores
    else:
        raise TypeError(f"gene_scores must be path or DataFrame, got {type(gene_scores).__name__}")
    return make_waterfall(
        df=df, rank_col=rank_by, highlights=highlights, output=output,
        background_known=background_known, title=title,
        figsize=figsize, annotate_fontsize=annotate_fontsize,
    )
