# Screen figure panels

Plotting library: `lib/screen_viz.py`. Category tables: `categories/`. Environment:
`envs/crispri-screen.yml`. Render everything with `bash run_all_figures.sh` (under Slurm,
`--ntasks=1`). Outputs go to the output root, at the same relative paths as the published files.

| Panel | Script (`panels/`) | Output |
|---|---|---|
| 1E, 1F, S1D, S1E, S1F, S2A, S6A | `fig1_figS1_figS2A_figS6A.py [--panels ...]` | `GenomeWide_results/plots/{Fig1,EssentialsScatter}/` |
| S3B | `figS3B_sublibrary_selection.py` | `Sublibrary_results/plots/GeneScatter/Genomewide_includedGenes_sublib*.pdf` |
| 2C | `fig2C_synth_vs_deg.py` | `fixed_NotI/plots/compare_{60S,40S}scores_bioBySubunit[_noLegend]_fixedNotI.pdf` |
| 2D, S4B | `fig2D_figS4B_category_violins.py` | `fixed_NotI/plots/{60S,40S}{deg,synth}_violinPlots_fixedNotI.pdf` |
| S3C | `figS3C_score_construction.py` | `fixed_NotI/plots/{60S,40S}_{synth,degrad}_score_dev[_flat]_square_fixedNotI.pdf` |
| S3D, S8A | `figS3D_figS8A_lfc_correlation.py` | `fixed_NotI/plots/heatmap_correlation_{wt,sgUBA5}_fixedNotI.pdf` |
| S3E | `figS3E_go_scores.py` | `fixed_NotI/plots/FigS3E_GO_scores_fixedNotI.pdf` (+ `_terms.csv`) |
| S3F | `figS3F_growth_vs_scores.py` | `fixed_NotI/plots/heatmap_growth_vs_synthdeg_wt_fixedNotI.pdf` |
| S4A | `figS4A_subunit_scatter.py` | `fixed_NotI/plots/compare_40S60S{Synth,Deg}scores_bioBySubunit_fixedNotI[_noLegend].pdf` |
| 6B, S8B, S8C | `fig6B_figS8BC_suppressor_scatter.py` | `fixed_NotI/plots/compare_{60S,40S}{Deg,Synth}scores_sgUBA5_R3HDM4cats_fixedNotI.pdf` |
| S8E, S8F | `figS8EF_growth_epistasis.py` | `fixed_NotI/plots/growth/growth_avgLines_{scatter*,waterfall_epistasis}.pdf` |

## Notes
- Sublibrary panels use the re-scored tree, `fixed_NotI/`: the 53 NotI/BamHI-site guides are excluded
  (see `../crispr_screens/scoring/`). `_fixedNotI` in the output names marks that tree.
- Six scripts are notebook cells (`Fig2.ipynb`, `Fig5.ipynb`) lifted verbatim, with path substitutions
  only; each header names its source cells.
- Gene labels in scatters are placed by adjustText, which moves them slightly from run to run.
  Labelled genes, data and axes are identical.
- Colours: non-targeting controls `#DDDDDD`, other genes `#4D4D4D`.
- Fig S3C, degradation: the panels pass x = old, y = ratioON with slope −1. The plotted colour is
  therefore (old − ratio)/√2, the opposite sign of the published degradation score (ratio − old)/√2.
  Kept as drawn in the manuscript; flagged for the authors.
- `categories/make_common_essentials.py` rebuilds `CommonEssentials.csv` from the DepMap Achilles
  common-essential list (download separately).
