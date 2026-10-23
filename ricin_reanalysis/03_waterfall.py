#!/usr/bin/env python3
"""Waterfall plot CLI. Thin wrapper around waterfall.py.

For programmatic use (e.g. from a Jupyter notebook), import waterfall.py
directly instead of shelling out to this script:

    import waterfall
    fig, ax = waterfall.plot_waterfall(
        gene_scores='results/K562_vs_Ricin_all_genes.csv',
        rank_by='epsilon',
        highlights=[dict(label='R3HDM4', color='#ffcc00', size=450,
                         marker='*', genes=['R3HDM4'], annotate=True)],
        background_known=True,
    )

CLI usage
---------
  --gene-scores PATH      CSV from 02_gene_scores.py
  --rank-by COL           'epsilon' or 'product' (default: epsilon)
  --output PATH           Output file (PDF + PNG saved)
  --highlight SPEC        One or more times. See DSL below.
  --background-known      Add literature-curated ricin-gene backdrop
  --title TEXT            Plot title
  --figsize WxH           Figure size in inches (default: 13.5x7.5)

Highlight DSL: one string of semicolon-separated key=value pairs:
    label=NAME           legend label (required)
    color=#HEX           matplotlib color (default #ff7f0e)
    size=N               marker size (default 30)
    marker=M             matplotlib marker (default 'o')
    alpha=F              opacity (default 0.9)
    edgecolor=COL        marker edge color (default white)
    edgewidth=N          marker edge width (default 0.4)
    zorder=N             stacking order (default by declaration order)
    genes=g1,g2,...      inline gene list
    genes=@PATH          file with one gene per line (# comments ok)
    annotate=true|false  draw gene labels next to points

Example:
    python 03_waterfall.py \
      --gene-scores results/K562_vs_Ricin_all_genes.csv \
      --rank-by epsilon \
      --output results/waterfall.pdf \
      --background-known \
      --highlight 'label=R3HDM4;color=#ffcc00;size=450;marker=*;genes=R3HDM4;annotate=true' \
      --highlight 'label=60S RPs;color=#9467bd;size=45;genes=@rpl_genes.txt'
"""

import argparse
import os
import sys

# Ensure the sibling waterfall.py is importable when this script is run directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import waterfall  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--gene-scores", required=True, help="CSV from 02_gene_scores.py")
    ap.add_argument("--rank-by", default="epsilon", choices=["epsilon", "product"])
    ap.add_argument("--output", required=True, help="Output file path (PDF+PNG saved)")
    ap.add_argument("--highlight", action="append", default=[],
                    help="Highlight spec (see header). May be given multiple times.")
    ap.add_argument("--background-known", action="store_true")
    ap.add_argument("--title", default=None)
    ap.add_argument("--figsize", default="13.5x7.5")
    ap.add_argument("--annotate-fontsize", type=float, default=8.0)
    args = ap.parse_args()

    fw, fh = (float(x) for x in args.figsize.lower().split("x"))

    highlights = []
    for spec in args.highlight:
        try:
            highlights.append(waterfall.parse_highlight_dsl(spec))
        except (ValueError, FileNotFoundError) as e:
            sys.exit(f"ERROR: {e}")

    print(f"[waterfall] {len(highlights)} highlight group(s)", file=sys.stderr)
    waterfall.plot_waterfall(
        gene_scores=args.gene_scores,
        rank_by=args.rank_by,
        highlights=highlights,
        output=args.output,
        background_known=args.background_known,
        title=args.title,
        figsize=(fw, fh),
        annotate_fontsize=args.annotate_fontsize,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
