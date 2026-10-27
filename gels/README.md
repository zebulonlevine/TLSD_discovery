# Gel and immunoblot quantification (Fig 3F, 5D, 7B; Fig S5H, S7F, S9F, S9G, S9I)

Bands were quantified in ImageStudio (LI-COR) with local background. The scripts here start from the
exported band intensities, normalize them, run the statistics, and draw the panels. The exports and
the gel scans are not part of this repository (Zenodo data record, DOI TBD). Place the exports under
`<data>/gels/<folder>/` as listed below.

Gel and blot images in the figures (Fig 3E, S5G, S7F right, S9E) are the scans themselves, with no code.

| Panel | Script | Input (`<data>/gels/…`) | Statistics |
|---|---|---|---|
| 3F, S5H | `fig3F_figS5H_rpl26_pulldown.py` | `3F_S5H/QuiescentPulldownFluorGels.xlsx` (sheet `data`) | Welch's t on log2 new/old, elute vs flowthrough within each construct, unadjusted |
| 5D (graph), S7F | `fig5D_figS7F_pt179_ufm1.py` | `5D_S7F/blotData_fullLane_UFM1.xlsx` | none (replicates and mean) |
| 7B, S9F, S9G, S9I | `fig7B_figS9FGI_lov_birA.py` | `7B_S9FGI/BirAgels.xlsx` (sheet `gelRedo`) | Welch's t, BH within each family of comparisons (see the script docstring); 95% CI from the t distribution |

Each script's docstring gives the processing steps and its output files. Every script also writes
its plotted values and statistics as CSV next to the PDFs in `<out>/gels/`. Significance brackets are
not drawn; the stars are in the statistics tables.

The LOV-BirA sheet holds a repeat run of one sgNT and one sgUBA5 replicate (gel 8, `extra` = Y).
As in the updated analysis notebook, these lanes count toward the gel-8 normalization and the sgNT
references, but are excluded from the plotted values and the tests, so every group has n = 6.
`--include-repeat-lanes` keeps them (n = 7 for sgNT and sgUBA5) and writes to
`<out>/gels/include_repeat_lanes/`.

## Run
```
python gels/<script>.py          # seconds each; envs/flow_cytometry.yml
```

## Verification
Each script was checked against the original notebook run unmodified on the same inputs
(2026-10-26, `flow` environment, pandas 2.3.3):
- **Values.** All plotted values, means, SEMs and 95% CIs are identical (max |Δ| = 0). The CSV files
  carry pandas' default float formatting, which rounds the last digit (≤ 4e-16).
- **Fig 7B, S9F, S9G, S9I drawing.** The original plotting function, given the same data, axis labels
  and limits, gives renders that are pixel-identical at 100 dpi, including the seeded point jitter.
- **Fig 3F, S5H, 5D, S7F drawing.** The original notebooks do not draw these layouts: the pulldown
  notebook plots all four fractions as bars, and the PT-179 notebook only loads the table. The scripts
  draw the published layouts, checked point by point against the figures.
- **Statistics.**
  - Fig 3F follows the legend (Welch's t). The notebook used Student's t: the 5K5R p is 0.044 instead
    of 0.041, and both are *.
  - The LOV-BirA panels follow the analysis notebook as updated on 2026-10-26, with one change: the
    notebook's families (Benjamini-Hochberg within each) and its exclusion of the repeat lanes after
    normalization are kept, but the test is Welch's t, as in the legends and Methods, instead of the
    notebook's Student's t. With Student's t, the tests equal the notebook's exactly. Plotted values
    are identical either way.

## Removed from the uploaded notebooks (not in the paper)
- Pulldown: old/new and percent-tagged plots, the per-fraction 549/646 plots, and the t-tests between
  fractions other than elute vs flowthrough.
- LOV-BirA: per-gel QC strip plots of every band, the HaloTag-normalized ratios, the within-genotype
  (day-0-normalized) quantities, the combined HaloTag/JF549 total, and the unused sheet `Sheet1`.
