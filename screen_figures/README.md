# Screen figure panels

Plotting library: `lib/screen_viz.py`. Category tables: `categories/`. Environment:
`envs/crispri-screen.yml`. Render everything with `bash run_all_figures.sh` (under Slurm,
`--ntasks=1`). Outputs go to the output root, at the same relative paths as the published files.

| Panel (draft5) | Script (`panels/`) | Output |
|---|---|---|
| 1F, 1G, S1C, S1D, S1E, S1F, S5A | `fig1_figS1_figS5A.py [--panels ...]` | `GenomeWide_results/plots/{Fig1,EssentialsScatter}/` |
| S2B | `figS2B_sublibrary_selection.py` | `Sublibrary_results/plots/GeneScatter/Genomewide_includedGenes_sublib*.pdf` |
| 2C | `fig2C_synth_vs_deg.py` | `fixed_NotI/plots/compare_{60S,40S}scores_bioBySubunit[_noLegend]_fixedNotI.pdf` |
| 2D, S3B | `fig2D_figS3B_category_violins.py` | `fixed_NotI/plots/{60S,40S}{deg,synth}_violinPlots_fixedNotI.pdf` |
| S2C | `figS2C_score_construction.py` | `fixed_NotI/plots/{60S,40S}_{synth,degrad}_score_dev[_flat]_square_fixedNotI.pdf` |
| S2D, S7A | `figS2D_figS7A_lfc_correlation.py` | `fixed_NotI/plots/heatmap_correlation_{wt,sgUBA5}_fixedNotI.pdf` |
| S2E | `figS2E_go_scores.py` | `fixed_NotI/plots/FigS3C_GO_scores_fixedNotI.pdf` (+ `_terms.csv`) |
| S2F | `figS2F_growth_vs_scores.py` | `fixed_NotI/plots/heatmap_growth_vs_synthdeg_wt_fixedNotI.pdf` |
| S3A | `figS3A_subunit_scatter.py` | `fixed_NotI/plots/compare_40S60S{Synth,Deg}scores_bioBySubunit_fixedNotI[_noLegend].pdf` |
| 6B, S7B, S7C | `fig6B_figS7BC_suppressor_scatter.py` | `fixed_NotI/plots/compare_{60S,40S}{Deg,Synth}scores_sgUBA5_R3HDM4cats_fixedNotI.pdf` |
| S7E, S7F | `figS7EF_growth_epistasis.py` | `fixed_NotI/plots/growth/growth_avgLines_{scatter*,waterfall_epistasis}.pdf` |

## Notes
- Sublibrary panels use the re-scored tree, `fixed_NotI/`: the 53 NotI/BamHI-site guides are excluded
  (see `../crispr_screens/scoring/`). `_fixedNotI` in the output names marks that tree.
- Some output names keep pre-draft5 numbering (`FigS3C_GO_scores` is draft5 Fig S2E).
- Six scripts are notebook cells (`Fig2.ipynb`, `Fig5.ipynb`) lifted verbatim, with path substitutions
  only; each header names its source cells.
- Gene labels in scatters are placed by adjustText, which moves them slightly from run to run.
  Labelled genes, data and axes are identical.
- Colours: non-targeting controls `#DDDDDD`, other genes `#4D4D4D`.
- Fig S2C, degradation: the panels pass x = old, y = ratioON with slope −1. The plotted colour is
  therefore (old − ratio)/√2, the opposite sign of the published degradation score (ratio − old)/√2.
  Kept as drawn in the manuscript; flagged for the authors.
- `categories/make_common_essentials.py` rebuilds `CommonEssentials.csv` from the DepMap Achilles
  common-essential list (download separately).
