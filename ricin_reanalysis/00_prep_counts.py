#!/usr/bin/env python3
"""Convert the Morgens 2017 ncomms15178-s4.xlsx (or any similarly-structured xlsx)
into a MAGeCK count table.

The xlsx has one sheet per replicate, with no header row, two columns
(sgRNA, count).  Gene names are parsed from the sgRNA name using the convention:

  0None_none_*    -> non-targeting
  0Safe_safe_*    -> non-targeting
  ENSG#####_GSYM_* -> gene GSYM

Usage:
  python 00_prep_counts.py \\
      --xlsx ncomms15178-s4.xlsx \\
      --ctrl-sheets "K562 replicate 1" "K562 replicate 2" \\
      --trt-sheets  "Ricin replicate 1" "Ricin replicate 2" \\
      --ctrl-names  K562_rep1 K562_rep2 \\
      --trt-names   Ricin_rep1 Ricin_rep2 \\
      --output      counts.txt
"""

import argparse
import sys
import pandas as pd


def gene_from_sgrna(s: str) -> str:
    """Assign gene symbol from sgRNA name per Morgens library convention."""
    if s.startswith("0None_") or s.startswith("0Safe_"):
        return "non-targeting"
    parts = s.split("_")
    if parts[0].startswith("ENSG") and len(parts) >= 2:
        return parts[1]
    return parts[0]


def load_sheet(xlsx: str, sheet: str, sample_name: str) -> pd.DataFrame:
    df = pd.read_excel(xlsx, sheet_name=sheet, header=None, names=["sgRNA", sample_name])
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True, help="Input Excel file")
    ap.add_argument("--ctrl-sheets", nargs="+", required=True, help="Sheet names for control samples")
    ap.add_argument("--trt-sheets", nargs="+", required=True, help="Sheet names for treated samples")
    ap.add_argument("--ctrl-names", nargs="+", required=True, help="Column names to use for ctrl samples in output")
    ap.add_argument("--trt-names", nargs="+", required=True, help="Column names to use for trt samples in output")
    ap.add_argument("--output", required=True, help="Output count file (MAGeCK format)")
    args = ap.parse_args()

    if len(args.ctrl_sheets) != len(args.ctrl_names):
        sys.exit("--ctrl-sheets and --ctrl-names must be same length")
    if len(args.trt_sheets) != len(args.trt_names):
        sys.exit("--trt-sheets and --trt-names must be same length")

    print(f"[prep_counts] reading {args.xlsx}", file=sys.stderr)
    merged = None
    pairs = list(zip(args.ctrl_sheets + args.trt_sheets,
                     args.ctrl_names + args.trt_names))
    for sheet, name in pairs:
        df = load_sheet(args.xlsx, sheet, name)
        print(f"  sheet '{sheet}' -> column '{name}': {len(df)} rows", file=sys.stderr)
        merged = df if merged is None else merged.merge(df, on="sgRNA", how="outer")

    all_cols = args.ctrl_names + args.trt_names
    n_na = merged[all_cols].isna().sum().sum()
    if n_na:
        print(f"  filling {int(n_na)} missing counts with 0 (outer-merged)", file=sys.stderr)
    for c in all_cols:
        merged[c] = merged[c].fillna(0).astype(int)

    merged["gene"] = merged["sgRNA"].apply(gene_from_sgrna)

    out = merged[["sgRNA", "gene"] + all_cols]
    out.to_csv(args.output, sep="\t", index=False)

    n_ntc = (out["gene"] == "non-targeting").sum()
    n_targeting = (out["gene"] != "non-targeting").sum()
    n_genes = out.loc[out["gene"] != "non-targeting", "gene"].nunique()
    print(f"[prep_counts] wrote {len(out)} sgRNAs to {args.output}", file=sys.stderr)
    print(f"  non-targeting: {n_ntc}", file=sys.stderr)
    print(f"  gene-targeting: {n_targeting} ({n_genes} unique genes)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
