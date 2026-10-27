#!/usr/bin/env python3
"""Write the p-ranked / matched-null gene scores (synthesis and degradation, sgNT
and sgUBA5 arms) into fixed_NotI/scores/<screen>_<score>[_sgUBA5]/ under the
output root. These are the published gene-level scores (Tables S5/S11, Fig 2,
Fig 6, Fig S2/S3/S7).

Inputs: the per-sgRNA line scores from run_line_scores.sh
(Sublibrary_results/scores/) and counts_for_mageck.txt, under the data root, and
the sublibrary design files shipped in crispr_screens/library/.
The 53 NotI/BamHI-site guides are excluded by default; --keep-restriction-sites
reproduces the superseded pre-2026-10-21 scores.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import final_pvalues as F
import sgrna_score_pvalues as M

DEFAULT_OUT = M.paths.outp('fixed_NotI/scores')
COLS = ['gene_id', 'epsilon', 'pvalue', 'n_sgrnas', 'gene', 'product']


def main(out_root=None, exclude_restriction_sites=True):
    """out_root lets a caller (e.g. audit_pipeline.py) write to a scratch
    directory instead of overwriting the committed scores -- an audit must not
    mutate the artifacts it is auditing.

    exclude_restriction_sites drops the 53 guides whose ordered oligo carries a
    NotI or BamHI site before the gene roll-up. They could not be cloned, but
    six of them were sequenced deeply enough to be scored, so they are in the
    published gene scores. Dropping them changes epsilon for 6 genes and 2 hit
    calls directly; the rest of the churn is the resampled FDR threshold moving.
    On by default (adopted 2026-10-21); off reproduces the superseded scores."""
    OUT = out_root or DEFAULT_OUT
    drop = frozenset()
    if exclude_restriction_sites:
        sys.path.insert(0, M.paths.repo('crispr_screens', 'pipeline'))
        import restriction_sites
        drop = restriction_sites.site_guides()
        print(f'excluding {len(drop)} NotI/BamHI guides from the gene roll-up')
    print(f"{'dir':28s} {'all_genes':>10s} {'hits':>6s} {'thr':>6s}")
    for screen in ['L29', 'S3']:
        for score in ['Synthesis', 'Degradation']:
            for cond in ['wt', 'sgUBA5']:
                name = f'{screen}_{score}' + ('' if cond == 'wt' else '_sgUBA5')
                d = os.path.join(OUT, name)
                os.makedirs(d, exist_ok=True)

                t = F.sgrna_table(screen, score, cond)
                if drop:
                    n0 = len(t)
                    t = t[~t['sgrna'].isin(drop)].reset_index(drop=True)
                    print(f'  {name}: {n0} -> {len(t)} guides after the site filter')
                g = F.gene_table(t)
                is_null = g['gene_id'].str.startswith('NTC_')
                g['gene'] = np.where(is_null, 'NTC', g['gene_id'])
                g['n_sgrnas'] = g['n_sgrnas'].astype(float)

                allg = g[COLS].sort_values('product', ascending=False)
                hits = g.loc[g['is_hit'], COLS].sort_values('product', ascending=False)
                allg.to_csv(os.path.join(d, 'gene_scores_all_genes.csv'), index=False)
                hits.to_csv(os.path.join(d, 'gene_scores_hits_fdr0.05.csv'), index=False)
                # sgRNA-level table, with the new per-guide p-values
                t.to_csv(os.path.join(d, f'{name}_sgrna_scores.csv'), index=False)
                print(f"{name:28s} {len(allg):10d} {len(hits):6d} {g.attrs['threshold']:6.2f}")


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default=None)
    ap.add_argument('--exclude-restriction-sites', action='store_true', default=True,
                    help='drop the 53 NotI/BamHI guides before the gene roll-up (default)')
    ap.add_argument('--keep-restriction-sites', dest='exclude_restriction_sites',
                    action='store_false',
                    help='keep them: reproduces the superseded pre-2026-10-21 scores')
    a = ap.parse_args()
    main(a.outdir, a.exclude_restriction_sites)
