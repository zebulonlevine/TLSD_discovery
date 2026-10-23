#!/usr/bin/env python3
"""Gene-level scoring via Mann-Whitney U test of top-N sgRNAs vs NTC pool.

Replicates the existing screen_pipeline.py `line-score-utest` workflow as a
standalone, self-contained script.

Inputs
------
  --sgrna-summary    MAGeCK sgrna_summary.txt
                     (must contain columns: sgrna, Gene, LFC, control_mean, treat_mean)

Outputs
-------
  PREFIX_sgrna_scores.csv    per-sgRNA table: sgRNA, gene, ctrl_mean, trt_mean, lfc, |lfc|
  PREFIX_all_genes.csv       per-gene table: gene, epsilon, pvalue, product, n_sgrnas, n_top
                             (includes NTC_* pseudo-genes for FDR null estimation)
  PREFIX_hits_fdr{FDR}.csv   gene-level hits passing the empirical FDR threshold

Method
------
  1. Median-normalized sgRNA LFCs come from MAGeCK (already in sgrna_summary.txt).
  2. Per gene, take the top-N sgRNAs by |LFC|.
  3. MWU two-sided test: top-N LFCs vs subsampled NTC pool (default 2000).
  4. NTC pseudo-genes built by drawing 5 random NTC sgRNAs and applying the
     same top-N selection + MWU. Used both for re-centering epsilon and for
     the empirical FDR on product score.
  5. Score outputs:
       epsilon  = mean LFC of top-N sgRNAs, minus NTC pseudo median epsilon
       product  = epsilon * -log10(pvalue)
"""

import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--sgrna-summary", required=True,
                    help="MAGeCK sgrna_summary.txt")
    ap.add_argument("--output-prefix", required=True,
                    help="Output prefix (incl. directory)")
    ap.add_argument("--ntc-gene-name", default="non-targeting",
                    help="Gene label for non-targeting controls in sgrna_summary.txt (default: non-targeting)")
    ap.add_argument("--top-n", type=int, default=3,
                    help="Top-N sgRNAs per gene by |LFC| for the U-test (default: 3)")
    ap.add_argument("--ntc-pool-size", type=int, default=2000,
                    help="Subsampled NTC pool size for U-test (default: 2000; lower = faster, less power)")
    ap.add_argument("--n-pseudo-genes", type=int, default=3000,
                    help="Number of NTC pseudo-genes for the empirical null (default: 3000)")
    ap.add_argument("--min-control-count", type=float, default=30.0,
                    help="Drop sgRNAs with control-mean below this threshold (default: 30)")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="Target empirical FDR (default: 0.05)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    # ------------------------------------------------------------------
    # 1) Load MAGeCK sgrna_summary.txt
    # ------------------------------------------------------------------
    print(f"[gene_scores] loading {args.sgrna_summary}", file=sys.stderr)
    sg = pd.read_csv(args.sgrna_summary, sep="\t")
    # MAGeCK column names: sgrna, Gene, control_count, treatment_count,
    # control_mean, treat_mean, LFC, control_var, adj_var, score, p.low, p.high, p.twosided, FDR, high_in_treatment
    expected = {"sgrna", "Gene", "LFC", "control_mean", "treat_mean"}
    missing = expected - set(sg.columns)
    if missing:
        sys.exit(f"ERROR: MAGeCK sgrna_summary missing columns: {missing}")

    sg = sg.rename(columns={"sgrna": "sgRNA", "Gene": "gene",
                            "control_mean": "ctrl_mean",
                            "treat_mean": "trt_mean",
                            "LFC": "lfc"})

    # Filter by control count to avoid noisy low-count sgRNAs
    n_before = len(sg)
    kept = sg[sg["ctrl_mean"] >= args.min_control_count].copy()
    print(f"[gene_scores] kept {len(kept)} / {n_before} sgRNAs (ctrl_mean >= {args.min_control_count})", file=sys.stderr)

    ntc_mask = kept["gene"] == args.ntc_gene_name
    ntc_lfcs = kept.loc[ntc_mask, "lfc"].to_numpy()
    print(f"[gene_scores] NTC pool: {len(ntc_lfcs)} sgRNAs (median LFC = {np.median(ntc_lfcs):.4f})",
          file=sys.stderr)
    if len(ntc_lfcs) < 10:
        sys.exit("ERROR: NTC pool too small to run U-test")

    # Save per-sgRNA table
    kept["abs_lfc"] = kept["lfc"].abs()
    sgrna_out = kept[["sgRNA", "gene", "ctrl_mean", "trt_mean", "lfc", "abs_lfc"]]
    sgrna_out.to_csv(f"{args.output_prefix}_sgrna_scores.csv", index=False)

    # ------------------------------------------------------------------
    # 2) Subsample NTC pool for U-test speed
    # ------------------------------------------------------------------
    if len(ntc_lfcs) > args.ntc_pool_size:
        ntc_pool = rng.choice(ntc_lfcs, size=args.ntc_pool_size, replace=False)
        print(f"[gene_scores] subsampled NTC pool to {args.ntc_pool_size} for U-test", file=sys.stderr)
    else:
        ntc_pool = ntc_lfcs

    # ------------------------------------------------------------------
    # 3) Per-gene U-test on top-N sgRNAs
    # ------------------------------------------------------------------
    kept_sorted = kept.sort_values(["gene", "abs_lfc"], ascending=[True, False])
    top = kept_sorted.groupby("gene", sort=False).head(args.top_n)
    n_total_per_gene = kept.groupby("gene").size()

    real_genes = [g for g in kept["gene"].unique() if g != args.ntc_gene_name]
    print(f"[gene_scores] running MWU for {len(real_genes)} genes...", file=sys.stderr)

    top_by_gene = {g: v.to_numpy() for g, v in top.groupby("gene", sort=False)["lfc"]}
    results = []
    for i, gene in enumerate(real_genes):
        arr = top_by_gene.get(gene)
        if arr is None or len(arr) == 0:
            continue
        try:
            _, p = mannwhitneyu(arr, ntc_pool, alternative="two-sided")
        except ValueError:
            p = 1.0
        results.append((gene, arr.mean(), p, int(n_total_per_gene.get(gene, 0)), len(arr)))
        if (i + 1) % 2000 == 0:
            print(f"  ... {i + 1}/{len(real_genes)}", file=sys.stderr)

    # ------------------------------------------------------------------
    # 4) NTC pseudo-genes for empirical null
    # ------------------------------------------------------------------
    print(f"[gene_scores] generating {args.n_pseudo_genes} NTC pseudo-genes", file=sys.stderr)
    ntc_list = list(ntc_lfcs)
    pseudo = []
    for j in range(args.n_pseudo_genes):
        sel = rng.choice(ntc_list, size=5, replace=False)
        sel_top = sel[np.argsort(-np.abs(sel))[:args.top_n]]
        try:
            _, p = mannwhitneyu(sel_top, ntc_pool, alternative="two-sided")
        except ValueError:
            p = 1.0
        pseudo.append((f"NTC_{j}", float(np.mean(sel_top)), float(p), 5, args.top_n))

    res = pd.DataFrame(results, columns=["gene", "epsilon", "pvalue", "n_sgrnas", "n_top"])
    res_pseudo = pd.DataFrame(pseudo, columns=["gene", "epsilon", "pvalue", "n_sgrnas", "n_top"])
    res = pd.concat([res, res_pseudo], ignore_index=True)

    # ------------------------------------------------------------------
    # 5) Re-center epsilon on NTC pseudo median; compute product score
    # ------------------------------------------------------------------
    ntc_med = res.loc[res["gene"].str.startswith("NTC_"), "epsilon"].median()
    print(f"[gene_scores] subtracting NTC pseudo median epsilon = {ntc_med:.4f}", file=sys.stderr)
    res["epsilon"] = res["epsilon"] - ntc_med
    res["product"] = res["epsilon"] * -np.log10(res["pvalue"].clip(lower=1e-300))

    # ------------------------------------------------------------------
    # 6) Empirical FDR on product score
    # ------------------------------------------------------------------
    abs_max = float(res["product"].abs().max())
    threshold = 0.0
    for cand in np.arange(0.0, abs_max + 0.05, 0.05):
        above = res[res["product"].abs() > cand]
        if len(above) == 0:
            continue
        frac_ntc = above["gene"].str.startswith("NTC_").mean()
        if frac_ntc < args.fdr:
            threshold = float(cand)
            break

    hits = res[
        (res["product"].abs() > threshold) & (~res["gene"].str.startswith("NTC_"))
    ].copy()
    print(f"[gene_scores] product-score threshold at FDR={args.fdr}: {threshold:.3f}", file=sys.stderr)
    print(f"[gene_scores] hits: {len(hits)}", file=sys.stderr)

    res.sort_values("product").to_csv(f"{args.output_prefix}_all_genes.csv", index=False)
    hits.sort_values("product").to_csv(f"{args.output_prefix}_hits_fdr{args.fdr}.csv", index=False)

    print(f"[gene_scores] wrote:", file=sys.stderr)
    print(f"  {args.output_prefix}_sgrna_scores.csv", file=sys.stderr)
    print(f"  {args.output_prefix}_all_genes.csv", file=sys.stderr)
    print(f"  {args.output_prefix}_hits_fdr{args.fdr}.csv", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
