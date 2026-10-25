#!/usr/bin/env python3
"""GO-term enrichment of the genome-wide new/old-ratio gene scores (Fig 1F, Table S3).

Univariate linear model (decoupler v2 `dc.mt.ulm`, tmin = 3) of gene epsilon on
membership in each GO term (Enrichr 2025 GO Biological Process + Molecular
Function + Cellular Component, `GSEA_Fig1/S3/gene_sets.gmt`). p-values are
recomputed from the t-statistics in log space and Benjamini-Hochberg-corrected
across all terms, because decoupler's own padj underflows to 0 for the top terms.

Extracted from Fig1Decoupler.ipynb cells 0-21. The two screens are processed in
the notebook's order: 40S (S3) adds the log-space p columns before the overlap
columns, 60S (L29) adds them after. The published files' column order follows
that, so it is kept.

Run with envs/dc.yml:
    python go_genomewide.py
Reads   ${TLSD_DATA}/GenomeWide_results/results/{S3,L29}_TopvBottom/*_all_genes.csv
        ${TLSD_DATA}/GSEA_Fig1/S3/gene_sets.gmt
Writes  ${TLSD_OUT}/GSEA_Fig1/S3/final_resultsS3Genomewide.csv
        ${TLSD_OUT}/GSEA_Fig1/L29/decoupler_results60S.csv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
import decoupler as dc

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from common import paths  # noqa: E402

GMT = paths.inp('GSEA_Fig1/S3/gene_sets.gmt')


def read_gmt(path):
    """Enrichr-style GMT -> decoupler net (source, target)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")   # TAB, not whitespace
            if len(fields) < 3:
                continue
            term, genes = fields[0], fields[2:]       # field 1 is the library name, skip it
            rows.extend((term, g) for g in genes if g)
    return pd.DataFrame(rows, columns=["source", "target"]).drop_duplicates()


def max_abs(s):
    return s.iloc[s.abs().values.argmax()]   # positional -> safe with duplicate index labels


def bh_log(logp):
    """Global Benjamini-Hochberg in log space."""
    order = np.argsort(logp)
    m = len(logp); ranks = np.arange(1, m + 1)
    adj = logp[order] + np.log(m) - np.log(ranks)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m); out[order] = np.minimum(adj, 0.0)
    return out


def gene_matrix(results):
    """1 x genes epsilon matrix (NTC pseudogenes dropped; duplicate symbols -> max |epsilon|)."""
    m = results.loc[results.gene != 'NTC', :].set_index("gene")[["epsilon"]].T
    return m.T.groupby(level=0).agg(max_abs).T


def run_ulm(crispr_matrix, go_network):
    scores, padj = dc.mt.ulm(data=crispr_matrix, net=go_network, tmin=3)
    final_results = (
        scores.stack().rename("t_statistic").to_frame()
              .join(padj.stack().rename("FDR_adj_p_value"))
              .reset_index()
    )
    final_results.columns = ["CRISPR_Condition", "Gene_Set", "t_statistic", "FDR_adj_p_value"]
    return final_results.sort_values("t_statistic")


def add_log_padj(final_results, crispr_matrix):
    df = crispr_matrix.shape[1] - 2     # genes in the regression, minus 2 -- same for all sets
    t = final_results["t_statistic"].abs().to_numpy()
    logp_raw = np.log(2) + stats.t.logsf(t, df)
    final_results["log_padj"] = bh_log(logp_raw)          # global BH, one family
    final_results["neglog10_FDR"] = -final_results["log_padj"] / np.log(10)


def add_overlap(final_results, crispr_matrix, go_network):
    screen_genes = set(crispr_matrix.columns)
    overlap = (go_network.groupby("source")["target"]
                         .apply(lambda t: sorted(screen_genes & set(t))))
    final_results["genes_in_set_and_screen"] = final_results["Gene_Set"].map(overlap.map(", ".join))
    final_results["n_overlap"] = final_results["Gene_Set"].map(overlap.map(len))


def main():
    res = lambda c: paths.inp(f'GenomeWide_results/results/{c}_TopvBottom/{c}_TopvBottom_all_genes.csv')
    go_network = read_gmt(GMT)

    # 40S (RPS3 reporter): ULM -> log-space BH -> overlap
    m = gene_matrix(pd.read_csv(res('S3')))
    fr = run_ulm(m, go_network)
    add_log_padj(fr, m)
    add_overlap(fr, m, go_network)
    fr.to_csv(paths.outp('GSEA_Fig1/S3/final_resultsS3Genomewide.csv', mkdir=True))

    # 60S (RPL29 reporter): ULM -> overlap -> log-space BH
    m = gene_matrix(pd.read_csv(res('L29')))
    fr = run_ulm(m, go_network)
    add_overlap(fr, m, go_network)
    add_log_padj(fr, m)
    fr.to_csv(paths.outp('GSEA_Fig1/L29/decoupler_results60S.csv', mkdir=True))
    print(f'GO terms tested: {len(fr)}')


if __name__ == '__main__':
    main()
