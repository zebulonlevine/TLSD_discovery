#!/usr/bin/env python
"""
plot_pearson_swarm.py — Per-cell Pearson swarm plot, control vs treatment
=========================================================================
Rebuilds output/pearson_L29_swarm_DMSO_vs_PT179_23h.pdf.

One panel per channel pair; one point per cell (whole-cell Pearson r);
horizontal bar = median; p = two-sided Mann-Whitney U test.

Images are selected by filename: the tag must contain `_{line}_` and
`_{condition}_{timepoint}_` (so DMSO_23h does not match DMSOGNE_23h).

Input is either the per_cell/ directory written by wholecell_pearson.py, or
the original merged all_pairwise.csv from quantify_coloc.py (whole_cell rows
are used, and --channel-map renames ch0/ch1/ch2).

Usage:
  python plot_pearson_swarm.py --input per_cell --output figures \\
    --line L29 --timepoint 23h --control DMSO --treatment PT179
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

matplotlib.rcParams["pdf.fonttype"] = 42  # editable text in Illustrator
matplotlib.rcParams["ps.fonttype"] = 42

# Point color / median-bar color per group (control grey, treatment orange)
CONTROL_COLORS = ("#bdbdbd", "#4b5563")
TREATMENT_COLORS = ("#f28e4c", "#c2410c")


def parse_args():
    p = argparse.ArgumentParser(description="Per-cell Pearson swarm plot")
    p.add_argument("--input", required=True,
                   help="per_cell/ directory or a merged all_pairwise.csv")
    p.add_argument("--output", required=True, help="Output directory")
    p.add_argument("--line", default="L29")
    p.add_argument("--timepoint", default="23h")
    p.add_argument("--control", default="DMSO")
    p.add_argument("--treatment", default="PT179")
    p.add_argument("--channel-map", default=None,
                   help="Rename channels, e.g. ch0=GFP-R3HDM4,ch1=Sec61B,ch2=L29")
    p.add_argument("--pairs", default=None,
                   help="Semicolon-separated pairs to plot, e.g. 'GFP-R3HDM4:Sec61B'")
    p.add_argument("--name", default=None, help="Output file stem")
    p.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"])
    return p.parse_args()


def load_pairwise(input_path):
    input_path = Path(input_path)
    if input_path.is_dir():
        files = sorted(input_path.glob("*_wholecell_pearson.csv"))
        if not files:
            raise SystemExit(f"No *_wholecell_pearson.csv in {input_path}")
        return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = pd.read_csv(input_path)
    if "compartment" in df.columns:
        df = df[df["compartment"] == "whole_cell"]
    return df


def select_group(df, line, condition, timepoint):
    pat = re.compile(rf"_{re.escape(line)}_.*_{re.escape(condition)}_{re.escape(timepoint)}_")
    images = sorted(i for i in df["image"].unique() if pat.search(i))
    return df[df["image"].isin(images)], images


def main():
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_pairwise(args.input)
    if args.channel_map:
        cmap = dict(kv.split("=") for kv in args.channel_map.split(","))
        df = df.assign(ch_a=df["ch_a"].replace(cmap), ch_b=df["ch_b"].replace(cmap))
    df = df.assign(pair=df["ch_a"] + " vs " + df["ch_b"])

    ctrl, ctrl_imgs = select_group(df, args.line, args.control, args.timepoint)
    trt, trt_imgs = select_group(df, args.line, args.treatment, args.timepoint)
    if not ctrl_imgs or not trt_imgs:
        raise SystemExit(f"No images matched: control={ctrl_imgs}, treatment={trt_imgs}")
    print(f"{args.control}: {ctrl_imgs}")
    print(f"{args.treatment}: {trt_imgs}")

    data = pd.concat([ctrl.assign(group=args.control), trt.assign(group=args.treatment)])
    data = data.dropna(subset=["pearson_r"])

    if args.pairs:
        pairs = [p.replace(":", " vs ") for p in args.pairs.split(";")]
    else:
        pairs = list(dict.fromkeys(df["pair"]))  # preserve file order

    order = [args.control, args.treatment]
    colors = {args.control: CONTROL_COLORS, args.treatment: TREATMENT_COLORS}

    fig, axes = plt.subplots(1, len(pairs), figsize=(5.2 * len(pairs), 5.2),
                             sharey=True, squeeze=False)
    stats_rows = []
    for ax, pair in zip(axes[0], pairs):
        sub = data[data["pair"] == pair]
        sns.swarmplot(data=sub, x="group", y="pearson_r", order=order,
                      hue="group", hue_order=order,
                      palette={g: c[0] for g, c in colors.items()},
                      size=3, linewidth=0.3, edgecolor="white", legend=False, ax=ax)

        x = sub.loc[sub["group"] == args.control, "pearson_r"].to_numpy()
        y = sub.loc[sub["group"] == args.treatment, "pearson_r"].to_numpy()
        u, p = mannwhitneyu(x, y, alternative="two-sided")

        for i, (g, vals) in enumerate(zip(order, (x, y))):
            med = np.median(vals)
            ax.hlines(med, i - 0.3, i + 0.3, color=colors[g][1], lw=2.5, zorder=5)
            ax.annotate(f"n={len(vals)}", (i, 0), xycoords=("data", "axes fraction"),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=9, color="#555555")
            stats_rows.append({"pair": pair, "group": g, "n_cells": len(vals),
                               "median_r": med, "mean_r": vals.mean()})
        stats_rows[-1].update({"mwu_U": u, "mwu_p_two_sided": p})

        ax.set_title(pair, fontsize=12, fontweight="bold")
        ax.text(0.5, 0.97, f"p={p:.1e}", transform=ax.transAxes,
                ha="center", va="top", fontsize=9, color="#374151")
        ax.set_xlabel("")
        ax.set_ylabel("Pearson r" if ax is axes[0][0] else "")
        ax.grid(axis="y", color="#e5e5e5", lw=0.6)
        ax.set_axisbelow(True)
        sns.despine(ax=ax)

    fig.suptitle(f"{args.line} — Pearson correlation at {args.timepoint} (whole cell)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()

    name = args.name or (f"pearson_{args.line}_swarm_{args.control}_vs_"
                         f"{args.treatment}_{args.timepoint}")
    fig.savefig(out_dir / f"{name}.{args.format}", dpi=300, bbox_inches="tight")
    stats = pd.DataFrame(stats_rows)
    stats.to_csv(out_dir / f"{name}_stats.csv", index=False)
    print(stats.to_string(index=False))
    print(f"\nSaved {out_dir / f'{name}.{args.format}'}")


if __name__ == "__main__":
    main()
