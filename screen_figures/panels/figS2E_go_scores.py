#!/usr/bin/env python3
"""Figure S3C -- GO enrichment of the sublibrary synthesis / degradation scores.

Backed by the same tables as Supplementary Table S6
(Sublibrary_results/GO_decoupler/*_GO_decoupler*.csv, produced by
Sublibrary_decoupler_GO.py). No code previously plotted this panel; it was
assembled by hand in Illustrator.

Writes the p-ranked version by default; --published reproduces the panel from
the published scores for side-by-side comparison.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, REPO)
LIB = os.path.join(REPO, 'screen_figures', 'lib')
sys.path.insert(0, LIB)
from common import paths  # noqa: E402
sys.modules.pop('screen_viz', None)
import screen_viz  # noqa: E402
assert os.path.dirname(os.path.abspath(screen_viz.__file__)) == LIB

GO = paths.inp('fixed_NotI', 'GO_decoupler')          # re-scored GO (go_sublibrary.py)
GO_JULY = paths.inp('Sublibrary_results', 'GO_decoupler')
CATS = paths.repo('screen_figures', 'categories')

# Draft palette. Lysosome+Autophagy and the three biogenesis categories are
# merged before colouring: the dominance rule compares categories against each
# other, so leaving them split would make a lysosomal term look ambiguous
# (Lysosome 0.3 vs Autophagy 0.25) and fall through to black.
MERGE = {
    'Lysosome': 'Lysosome/Autophagy',
    'Autophagy': 'Lysosome/Autophagy',
    'Both Biogenesis': 'Biogenesis',
    '40S Biogenesis/Proteins': 'Biogenesis',
    '60S Biogenesis/Proteins': 'Biogenesis',
}
COLORS = {
    'Biogenesis': '#4e9342',
    'Lysosome/Autophagy': '#c51b7d',
    'UFMylation': '#E69F00',
}
# Abbreviations used in the draft labels.
TERM_LABELS = {
    'Proton Transmembrane Transport (GO:1902600)': 'H+ Transmembrane Transport',
    'Vacuolar Proton-Transporting V-type ATPase Complex (GO:0016471)':
        'Vacuolar H+ Transporting V-ATPase',
    'Cation-Transporting ATPase Complex (GO:0090533)': 'Cation Transporting ATPase Complex',
    'Protein Ufmylation (GO:0071569)': 'Protein UFMylation',
    'Regulation of Intracellular Steroid Hormone Receptor Signaling Pathway (GO:0033143)':
        'Reg. of Intracellular Steroid Hormone Signaling',
    'Regulation of Intracellular Estrogen Receptor Signaling Pathway (GO:0033146)':
        'Reg. of Intracellular Estrogen Receptor Signaling',
    'Regulation of Proteolysis Involved in Protein Catabolic Process (GO:1903050)':
        'Reg. of Proteolysis Involved in Catabolic Process',
    'Regulation of Proteasomal Ubiquitin-Dependent Protein Catabolic Process (GO:0032434)':
        'Reg. of Proteasomal Ubiquitin-Dep. Protein Catabolic Process',
}


def cats():
    df = screen_viz.load_categories(os.path.join(CATS, 'manyCats.csv'))
    df['category'] = df['category'].astype(str).str.strip().replace(MERGE)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--published', action='store_true',
                    help='use the July GO tables (Sublibrary_results/GO_decoupler/) instead')
    ap.add_argument('--top-n', type=int, default=6)
    ap.add_argument('--go-dir', default=None, help='directory of *_GO_decoupler*.csv')
    ap.add_argument('--suffix', default=None)
    ap.add_argument('--outdir', default=paths.outp('fixed_NotI', 'plots'))
    args = ap.parse_args()

    sfx = args.suffix if args.suffix is not None else ('' if args.published else '_fixedNotI')
    GO_DIR = args.go_dir or (GO_JULY if args.published else GO)
    out = os.path.join(args.outdir, f'FigS3C_GO_scores{sfx or "_published"}.pdf')
    os.makedirs(os.path.dirname(out), exist_ok=True)

    def p(name):
        return os.path.join(GO_DIR, f'{name}_GO_decoupler{sfx}.csv')

    panels = [
        dict(title='Degradation Score', noun='degradation',
             groups=[('60S-enriched\nGO Terms', p('60S_Degradation')),
                     ('40S-enriched\nGO Terms', p('40S_Degradation'))]),
        dict(title='Synthesis Score', noun='synthesis',
             groups=[('60S-enriched\nGO Terms', p('60S_Synthesis')),
                     ('40S-enriched\nGO Terms', p('40S_Synthesis'))]),
    ]
    df = screen_viz.plot_go_score_panels(
        panels, out, top_n=args.top_n,
        categories=cats(), colors=COLORS, term_labels=TERM_LABELS,
        tick_step=5.0, figsize=(15.0, 6.6), wrap_width=24, wspace=0.34,
        fontsize=7)
    df.to_csv(out.replace('.pdf', '_terms.csv'), index=False)
    print(f"  {len(df)} rows; terms table -> {os.path.basename(out).replace('.pdf','_terms.csv')}")


if __name__ == '__main__':
    main()
