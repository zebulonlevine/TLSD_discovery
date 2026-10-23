# Ricin screen pipeline (Morgens 2017 ncomms15178)

Modular pipeline that takes the Morgens 2017 raw count xlsx, runs MAGeCK for
median normalization + sgRNA fold-changes, applies the existing
Mann-Whitney U-test gene scoring (top-N sgRNAs vs NTC pool), and produces
a waterfall plot with user-controllable highlight sets.

## Files

```
ricin_pipeline/
├── run_all.sh                 orchestrator (edit CONDA_ENV at top)
├── 00_prep_counts.py          xlsx → MAGeCK count.txt
├── 01_run_mageck.sh           conda activate + mageck test
├── 02_gene_scores.py          MWU on top-N sgRNAs vs NTC pool
├── waterfall.py               importable library: plot_waterfall(), Highlight
├── 03_waterfall.py            thin CLI wrapper around waterfall.py
├── fig7F_ricin_waterfall.py   the manuscript panel (Fig 7F)
├── rpl_genes.txt              60S ribosomal proteins
└── ricin_trafficking.txt      genes known in ricin trafficking
```

## Setup

```bash
# 1. Download Supplementary Data 4 of Morgens et al. 2017 (Nat Commun 8:15178),
#    ncomms15178-s4.xlsx, into ${TLSD_DATA}/ricin_data_reanalysis/.
# 2. Create the env: conda env create -f ../envs/mageck.yml  (MAGeCK 0.5.9.5 + pandas/openpyxl)
# 3. With that env active:
bash run_all.sh                      # counts -> MAGeCK -> gene scores -> waterfall
python fig7F_ricin_waterfall.py      # the manuscript panel
```
Independent validation of the hits uses Bassik et al. 2013 (Cell, PMID 23394947), not
Morgens 2017 itself.

Outputs land in `${TLSD_OUT}/ricin_data_reanalysis/results/`:
- `counts.txt` — MAGeCK count table
- `K562_vs_Ricin.sgrna_summary.txt` — MAGeCK per-sgRNA stats (with LFCs)
- `K562_vs_Ricin_sgrna_scores.csv` — identically-processed per-sgRNA scores
- `K562_vs_Ricin_all_genes.csv` — gene-level table (ε, p, product, n_sgrnas)
- `K562_vs_Ricin_hits_fdr0.05.csv` — gene hits at FDR 0.05
- `K562_vs_Ricin_waterfall.pdf` / `.png` — waterfall plot

## Iterating on highlights

The waterfall step is independent of MAGeCK — once you have
`*_all_genes.csv`, you can re-render the plot with new highlights in
seconds.

### From Python

`fig7F_ricin_waterfall.py` is the manuscript example. In a notebook,
`import waterfall` and type `Highlight(` / `plot_waterfall(` for the fields.

Basic call:

```python
import waterfall
from waterfall import plot_waterfall, Highlight

fig, ax = plot_waterfall(
    gene_scores="results/K562_vs_Ricin_all_genes.csv",
    rank_by="epsilon",
    background_known=True,
    highlights=[
        Highlight(label="R3HDM4", color="#ffcc00", size=450, marker="*",
                  genes=["R3HDM4"], annotate=True),
        Highlight(label="60S RPs", color="#9467bd", size=45,
                  genes="@rpl_genes.txt"),
    ],
)
```

`plot_waterfall()` returns `(fig, ax)`, so you can further customize the
plot before saving. Pass `output="path.pdf"` to save both PDF and PNG.

### From the command line

```bash
python 03_waterfall.py \
    --gene-scores results/K562_vs_Ricin_all_genes.csv \
    --rank-by epsilon \
    --output results/waterfall_custom.pdf \
    --background-known \
    --highlight 'label=R3HDM4;color=#ffcc00;size=450;marker=*;genes=R3HDM4;annotate=true' \
    --highlight 'label=60S RPs;color=#9467bd;size=45;genes=@rpl_genes.txt' \
    --highlight 'label=Ricin trafficking;color=#e41a1c;size=35;genes=@ricin_trafficking.txt'
```

### Highlight spec syntax

Each `--highlight` value is a single string of semicolon-separated `key=value`
pairs:

| key         | meaning                                  | default       |
|-------------|------------------------------------------|---------------|
| `label`     | legend label (required)                  | —             |
| `color`     | matplotlib color (hex or name)           | `#ff7f0e`     |
| `size`      | marker size                              | `30`          |
| `marker`    | matplotlib marker (`o`, `*`, `s`, `D`)   | `o`           |
| `alpha`     | opacity 0–1                              | `0.9`         |
| `edgecolor` | marker edge color                        | `white`       |
| `edgewidth` | marker edge width                        | `0.4`         |
| `zorder`    | plot stacking order (higher = on top)    | auto by order |
| `genes`     | inline list `g1,g2,...` or `@file.txt`   | required      |
| `annotate`  | draw gene labels next to points          | `false`       |

Highlights drawn in order; later highlights overlay earlier ones. To force a
specific layering, set `zorder` explicitly.

### Switching the rank axis

By default genes are ranked by ε (effect size). To rank by product score
(epsilon × −log10 p) instead:

```bash
python 03_waterfall.py ... --rank-by product
```

## Notes

- MAGeCK normalization uses non-targeting sgRNAs as the size-factor anchors.
- The MWU pipeline uses `min-control-count=30` by default to drop noisy
  low-coverage sgRNAs; tune via `--min-control-count` in `02_gene_scores.py`.
- For speed, the NTC pool is subsampled to 2000 per U-test. The empirical
  FDR threshold reflects this; if you want stricter FDR, increase
  `--ntc-pool-size` (slower but more powerful).
- The pseudo-gene null re-centers ε so that the median NTC pseudo-gene sits at
  ε = 0; this is consistent with the existing screen_pipeline.py behavior.

## Parallel-pipeline notes

This re-runs MAGeCK from scratch each invocation. If you want to feed
MAGeCK's output into your existing `screen_pipeline.py line-score-utest`
subcommand instead of `02_gene_scores.py`, the `sgrna_summary.txt` file
produced by step 01 is the correct input.
