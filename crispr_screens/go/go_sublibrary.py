#!/usr/bin/env python
"""Sublibrary-background GO enrichment (decoupler v2 ULM) of the synthesis and
degradation gene scores for both subunits (60S = L29, 40S = S3): Fig S2E and
Table S6.

Same method and GMT as go_genomewide.py. The difference is that the input is a
sublibrary score file (715 genes), so the ULM background is exactly the
sublibrary. Merges Sublibrary_decoupler_GO.py with its re-run wrapper
run_go_pranked.py.

Run with envs/dc.yml:
    python go_sublibrary.py            # published: re-scored tree, *_fixedNotI.csv
Reads   ${TLSD_DATA}/fixed_NotI/scores/<dir>/gene_scores_all_genes.csv  (export_scores.py)
        ${TLSD_DATA}/GSEA_Fig1/S3/gene_sets.gmt
Writes  ${TLSD_OUT}/fixed_NotI/GO_decoupler/{60S,40S}_{Synthesis,Degradation}_GO_decoupler_fixedNotI.csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
import decoupler as dc

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from common import paths  # noqa: E402

GMT = paths.inp("GSEA_Fig1", "S3", "gene_sets.gmt")

# (label, score-dir, output filename)
RUNS = [
    ("60S Synthesis",   "L29_Synthesis",   "60S_Synthesis_GO_decoupler.csv"),
    ("60S Degradation", "L29_Degradation", "60S_Degradation_GO_decoupler.csv"),
    ("40S Synthesis",   "S3_Synthesis",    "40S_Synthesis_GO_decoupler.csv"),
    ("40S Degradation", "S3_Degradation",  "40S_Degradation_GO_decoupler.csv"),
]


def max_abs(s):
    """Keep the value with the largest |epsilon| (preserves sign). Positional -> safe
    with duplicate index labels."""
    return s.iloc[s.abs().values.argmax()]


def load_go_network(path):
    """Parse an Enrichr-style GMT into decoupler v2 net format (source, target)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")   # TAB, not whitespace
            if len(fields) < 3:
                continue
            term, genes = fields[0], fields[2:]       # field 1 = library name, skip
            rows.extend((term, g) for g in genes if g)
    return pd.DataFrame(rows, columns=["source", "target"]).drop_duplicates()


def bh_log(logp):
    """Global Benjamini-Hochberg in log space (avoids underflow to 0)."""
    order = np.argsort(logp)
    m = len(logp)
    ranks = np.arange(1, m + 1)
    adj = logp[order] + np.log(m) - np.log(ranks)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 0.0)
    return out


def run_one(label, scores_root, score_dir, out_path, go_network):
    src = os.path.join(scores_root, score_dir, "gene_scores_all_genes.csv")
    df = pd.read_csv(src)
    df = df.loc[df.gene != "NTC", :]

    # 1 observation ("epsilon") x N genes, background = these genes (the sublibrary)
    crispr_matrix = df.set_index("gene")[["epsilon"]].T
    crispr_matrix = crispr_matrix.T.groupby(level=0).agg(max_abs).T
    n_genes = crispr_matrix.shape[1]

    scores, padj = dc.mt.ulm(data=crispr_matrix, net=go_network, tmin=3)

    res = (
        scores.stack().rename("t_statistic").to_frame()
              .join(padj.stack().rename("FDR_adj_p_value"))
              .reset_index()
    )
    res.columns = ["Score", "Gene_Set", "t_statistic", "FDR_adj_p_value"]
    res["Score"] = label

    # manual global log-space BH (decoupler's padj underflows to 0 for top terms)
    ddf = n_genes - 2
    t = res["t_statistic"].abs().to_numpy()
    logp_raw = np.log(2) + stats.t.logsf(t, ddf)
    res["log_padj"] = bh_log(logp_raw)
    res["neglog10_FDR"] = -res["log_padj"] / np.log(10)

    # overlap annotation vs the sublibrary background
    screen_genes = set(crispr_matrix.columns)
    overlap = (go_network.groupby("source")["target"]
               .apply(lambda tt: sorted(screen_genes & set(tt))))
    res["genes_in_set_and_screen"] = res["Gene_Set"].map(overlap.map(", ".join))
    res["n_overlap"] = res["Gene_Set"].map(overlap.map(len))

    res = res.sort_values("t_statistic").reset_index(drop=True)
    res.to_csv(out_path, index=False)
    print(f"[{label}] genes={n_genes}  GO terms tested={len(res)}  -> {out_path}")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--scores', default=paths.inp('fixed_NotI', 'scores'),
                    help='gene-score tree laid out like fixed_NotI/scores/')
    ap.add_argument('--outdir', default=paths.outp('fixed_NotI', 'GO_decoupler'))
    ap.add_argument('--suffix', default='_fixedNotI', help='tag appended to each output name')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    go_network = load_go_network(GMT)
    print(f"GMT: {go_network['source'].nunique()} terms, {len(go_network)} edges")
    for label, score_dir, out_name in RUNS:
        out = os.path.join(args.outdir, out_name.replace('.csv', f'{args.suffix}.csv'))
        run_one(label, args.scores, score_dir, out, go_network)


if __name__ == "__main__":
    main()
