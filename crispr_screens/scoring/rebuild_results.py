#!/usr/bin/env python3
"""Re-run the per-comparison gene-level U-test with the NotI/BamHI guides removed.

Produces the `<comp>_all_genes.csv` / `<comp>_hits_fdr0.05.csv` pair that feeds
the `*_lfc_*` blocks of Tables S5/S11, into fixed_NotI/results/<comp>/.

MAGeCK itself is NOT re-run: its sgRNA-level output is per-guide, so dropping
guides from `<comp>.sgrna_summary.txt` and re-running the gene roll-up is
equivalent to having excluded them from the gene test, without touching the
LFCs of the guides that remain. The filtered summary is written alongside so the
input to each gene table is inspectable.

Only the 16 comparisons the supplementary tables actually use are rebuilt
(BLOCKS x {wt, sgUBA5}); the other 28 feed nothing downstream.

Verification: with --no-filter this reproduces the committed gene tables
byte-for-byte, which is what pins fdr/min_treat_mean to the values the original
run used.
"""
import argparse
import os
import sys

import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'crispr_screens', 'pipeline'))
from common import paths                      # noqa: E402
import restriction_sites                      # noqa: E402
import screen_pipeline as SP                  # noqa: E402

RES = paths.inp('Sublibrary_results/results')
OUT = paths.outp('fixed_NotI/results')

# from build_sublibrary_tables.BLOCKS -- the comparisons S5/S11 read
STEMS = ['new_top_v_bottom', 'old_top_v_bottom', 'ratioON_top_v_bottom', 'unsorted_v_early']
FDR = 0.05
MIN_TREAT_MEAN = 5          # config_sublibrary.yaml


def comparisons():
    for cond in ('wt', 'sgUBA5'):
        for line in ('L29', 'S3'):
            for stem in STEMS:
                yield f'{line}_{stem}_{cond}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default=OUT)
    ap.add_argument('--no-filter', action='store_true',
                    help='keep every guide; used to prove the roll-up reproduces '
                         'the committed tables')
    args = ap.parse_args()

    drop = frozenset() if args.no_filter else restriction_sites.site_guides()
    print(f'excluding {len(drop)} NotI/BamHI guides\n')
    print(f"{'comparison':34s} {'guides':>14s} {'genes':>6s} {'hits':>5s} {'thr':>6s}")
    for comp in comparisons():
        src = os.path.join(RES, comp, f'{comp}.sgrna_summary.txt')
        if not os.path.exists(src):
            print(f'  WARNING: missing {src}'); continue
        d = os.path.join(args.outdir, comp)
        os.makedirs(d, exist_ok=True)
        # Filter as TEXT, not through pandas. Round-tripping the summary through
        # read_table/to_csv perturbs ~10 of the p-values by one ulp (they run down
        # to 1e-36, where the C parser and shortest-repr disagree), which then
        # shows up as a spurious S4/S10 diff against the published tables.
        with open(src) as fh:
            lines = fh.readlines()
        header, body = lines[0], lines[1:]
        n0 = len(body)
        body = [ln for ln in body if ln.split('\t', 1)[0] not in drop]
        filt = os.path.join(d, f'{comp}.sgrna_summary.txt')
        with open(filt, 'w') as fh:
            fh.write(header)
            fh.writelines(body)
        sg = body
        _, hits, thr = SP.mann_whitney_analysis(
            filt, os.path.join(d, comp), fdr=FDR, min_treat_mean=MIN_TREAT_MEAN)
        allg = pd.read_csv(os.path.join(d, f'{comp}_all_genes.csv'))
        print(f'{comp:34s} {n0:6d} -> {len(sg):5d} {len(allg):6d} {len(hits):5d} {thr:6.1f}')


if __name__ == '__main__':
    main()
