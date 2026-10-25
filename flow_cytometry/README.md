# Flow cytometry (Fig 1C, 3H, 4B, 4D, 5G, 6C, 7A; Fig S2B, S5I, S7G, S8D, S9B–D)

Gating and per-well median fluorescence were done in FlowJo. Fig 5G is the exception: its gating was
done in Python (`fcs_gating.py`). The scripts here start from the per-well exports, normalize, run
the statistics, and draw the panels. The exports and FCS files are not part of this repository
(data deposit TBD). Place them under `<data>/flow_cytometry/<folder>/` as listed below.

Panels drawn entirely in FlowJo (no code): Fig S1A right (EdU), the per-cell histograms in Fig 1C,
and Fig S1C.

| Panel | Script | Input (`<data>/flow_cytometry/…`) | Statistics |
|---|---|---|---|
| 1C | `fig1C_ribohalo_turnover.py` | `1C/day{1,3,…,15}.csv`, `1C/PlateMap_Feb.xlsx` | A·exp(−kt), weighted fit; t½ = ln2/k with 95% CI |
| S2B | `figS2B_autophagy_reporter.py` | `S2B/autophagyData.csv` | % change vs matched control; delta-method SE; t at Welch–Satterthwaite df, unadjusted |
| 3H, S5I | `fig3H_figS5I_riok3_uba5.py` | `3H_S5I/RIOK3.xlsx` (sheet `guide2`, sgRIOK3-g3) | Welch's t, BH within each panel |
| 4B, 4D | `fig4BD_serum_confluence_gne.py` | `4B_4D/Confluence_GNE_5d.xlsx` | Welch's t sgUBA5 vs sgNT per condition, BH across the 4 conditions |
| 5G | `fig5G_lov_birA.py` (+ `fcs_gating.py`) | `5G/channel_statistics_dark.csv` (or FCS files: `--fcs-dir`, `--metadata`) | t with propagated error, Welch–Satterthwaite df |
| 6C, S8D | `fig6C_figS8D_mkeima.py` | `6C_S8D/DataWorkup_L29_mKeima_CJZ18.xlsx` | Welch's t; BH across the 6 BafA1 tests and across the 15 between-group pairs |
| 7A, S9B–D | `fig7A_figS9BD_uba5_tlsd1.py` | `7A_S9B-D/R3HDM4_L29_knockdown_cover_3xHA.xlsx` | Welch's t, BH within each panel |
| S7G | `figS7G_degron_timecourse.py` | `S7G/flow_timecourse_2603.xlsx` | Welch's t, DMSO vs PT-179, BH across the experiment's 96 tests |

`strip_plots.py` holds the shared normalization, strip plot and pairwise-test code for 3H, 4B, 4D,
7A, S5I and S9B–D. Each script's docstring gives the processing steps and its output files. Every
script also writes its plotted values and statistics as CSV next to the PDFs in `<out>/flow_cytometry/`.

## Run
```
python flow_cytometry/<script>.py          # seconds each; envs/flow_cytometry.yml
```
Fig 5G uses the Arial font (matplotlib substitutes a default font if Arial is not installed).
To rebuild the Fig 5G per-well table from FCS files, run
`python flow_cytometry/fig5G_lov_birA.py --fcs-dir <dir> --metadata <csv>`.

## Verification
Each script was checked against the original notebook or script run on the same exports
(2026-10-24, pandas 2.3.3). All statistics and plotted values are identical (max |Δ| = 0), and all
panels are pixel-identical at 100 dpi, with these exceptions:
- **Strip-plot point jitter (3H, 4B, 4D, 7A, S5I, S9B–D).** The originals jittered points with an
  unseeded random generator; the scripts seed it. The panels were compared with jitter off.
- **Fig 5G.** The table's row order differs from the original metadata, so points scatter slightly
  differently. Values are identical.
- **Fig S7G Total.** The placed panel was drawn without the notebook's final per-day adjustment; the
  script reproduces the placed panel. No star changes either way.
- **Fig S2B.** Only the panel shown in the paper is drawn: the RB1CC1 and ATG7 columns of the
  original figure's top row. The BafA1 column and its table rows were dropped on 2026-10-25; the
  remaining values are unchanged.
- **Fig 6C, S8D (changed 2026-10-25).** The original corrected with Holm over all 30 ordered pairs, so
  each pair was counted twice. The script uses BH over the 15 unique pairs, as the legends and
  Methods state. Raw p-values are unchanged.

## Removed from the uploaded notebooks (not in the paper)
- Old/new ratio, CV and raw-median plots.
- The 3xHA-TLSD1 rescue lines.
- The second RIOK3 sgRNA experiment (sheet `guide1`, sgRIOK3-2), which is not shown in 3H/S5I.
- The S2B BafA1 column, the basal-normalized row, and v1 of the S2B script (`make_figure.py`), which
  the placed panel does not use.
- The 40S confluence plots.
- Earlier trials in the 1C notebook (RPS25 line, RPE1 background, viability, unweighted fits).
- The overview and genotype-comparison plots from S7G, plus an earlier draft notebook and an
  exploratory script.
- Gate QC plots in `fcs_gating.py`, and the rest of the original 2,451-line FCS library.
