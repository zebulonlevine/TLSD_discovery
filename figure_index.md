# Figure index

This index lists every panel of Figures 1–7 and Supplementary Figures S1–S10, and every Supplementary
Table S1–S12, with the script that makes it and the source of its inputs. Panel letters follow the
manuscript figure legends. Script paths are relative to the repository root. Panels made without code
give the reason instead of a script. Run order, environments, output file names, expected numbers and
verification notes are in the module READMEs: `crispr_screens/pipeline/README.md`,
`crispr_screens/scoring/README.md`, `crispr_screens/go/README.md`, `crispr_screens/tables/README.md`,
`screen_figures/README.md`, `ricin_reanalysis/README.md`, `proteomics/README.md`,
`flow_cytometry/README.md`, `gels/README.md`, `single_molecule/README.md`,
`confocal_colocalization/README.md` and `simulation/README.md`. Set the data and output roots as in
the root `README.md` (Setup).

## Input sources
| Short name | What and where |
|---|---|
| GEO genome-wide counts | FASTQs and `genomewide_aggregated_counts.csv` at GEO (accession TBD), with the library design files in `crispr_screens/library/` (ship with the code). Processed by `crispr_screens/pipeline/`; GO by `crispr_screens/go/go_genomewide.py`. The results are Supplementary Tables S1–S3. |
| GEO sublibrary counts | FASTQs and `sublibrary_aggregated_counts.csv` at GEO (accession TBD), with the library design files in `crispr_screens/library/`. Processed by `crispr_screens/pipeline/`, then `crispr_screens/scoring/`; GO by `crispr_screens/go/go_sublibrary.py`. The results are Supplementary Tables S4–S6, S10 and S11. |
| `screen_figures/categories/` | Gene category tables; ship with the code. |
| `proteomics/input/` | Processed TMT quantification tables and the topology annotation; ship with the code. Raw files and searches: PRIDE (accession PXD084612). |
| Zenodo `<name>` | `<name>.zip` in the Zenodo data record (DOI [10.5281/zenodo.23688082](https://doi.org/10.5281/zenodo.23688082)). The record's `MANIFEST_public.csv` lists every file with its panels and the script that reads it. |
| Morgens 2017 SD4 | Morgens *et al.* 2017 Supplementary Data 4 (`ncomms15178-s4.xlsx`); third-party download. |
| Enrichr GO 2025 | Enrichr `GO_Biological_Process_2025`, `GO_Molecular_Function_2025` and `GO_Cellular_Component_2025`, concatenated into one GMT; third-party download. |
| OpenCell mmc4 | Organelle-IP localization table (Hein *et al.* 2025, Cell, supplementary file mmc4); third-party download. Read by every proteomics script. |
| DepMap | DepMap Achilles common essentials; third-party download. Already converted to `screen_figures/categories/CommonEssentials.csv` by `screen_figures/categories/make_common_essentials.py`. |
| MitoCarta 3.0 | Third-party download. Already merged into `proteomics/input/total_topology.tsv` by `proteomics/annotation/uniprot_topology_pipeline.py`. |

The single-molecule panels also need the 36-state likelihood caches, and Fig S7A, S7C and S7D the
joint caches. Both are rebuilt from the deposited trajectories (`single_molecule/README.md`, run order
steps 2 and 3).

## Figure 1
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| 1A | HaloTag7 knock-in at RPS3 and RPL29 | No code: structure render (ChimeraX, PDB 3J7R) | — |
| 1B | RiboHalo pulse-chase labeling assay | No code: schematic | — |
| 1C | Turnover curves and per-cell new/old histograms | `flow_cytometry/fig1C_ribohalo_turnover.py` (curves); histograms drawn in FlowJo, no code | Zenodo `flow_RiboHalo-turnover-timecourse` |
| 1D | Genome-wide new/old ratio screen | No code: schematic | — |
| 1E | 40S vs 60S new/old gene scores | `screen_figures/panels/fig1_figS1_figS2A_figS6B.py --panels 1E` | GEO genome-wide counts; `screen_figures/categories/` |
| 1F | GO terms enriched for new/old ratio | `screen_figures/panels/fig1_figS1_figS2A_figS6B.py --panels 1F` | GEO genome-wide counts; Enrichr GO 2025 |

## Figure 2
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| 2A | Synthesis vs degradation effects on ratio, abundance | No code: schematic | — |
| 2B | Multi-sort screen and score construction | No code: schematic | — |
| 2C | Gene synthesis vs degradation scores, 60S, 40S | `screen_figures/panels/fig2C_synth_vs_deg.py` | GEO sublibrary counts; `screen_figures/categories/` |
| 2D | Degradation scores by gene category | `screen_figures/panels/fig2D_figS4B_category_violins.py` | GEO sublibrary counts; `screen_figures/categories/` |

## Figure 3
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| 3A | SILAC-TMT pulse-chase design | No code: schematic | — |
| 3B | Old-protein turnover change by category | `proteomics/silac_turnover.py` | `proteomics/input/protein_quant_30972.tsv`; OpenCell mmc4 |
| 3C | RPL26 UFMylation cascade at the ER | No code: schematic | — |
| 3D | Non-UFMylatable RPL26 pulldown assay | No code: schematic | — |
| 3E | UFM1 immunoblot of RPL26 variants | No code: representative immunoblot | Raw scans: Zenodo `gel_RPL26-3xFLAG_pulldown` |
| 3F | New/old ratio, tagged vs untagged ribosomes | `gels/fig3F_figS5H_rpl26_pulldown.py` | Zenodo `gel_RPL26-3xFLAG_pulldown` |
| 3G | RNF10-RIOK3 40S degradation pathway | No code: schematic | — |
| 3H | Old 40S, 60S after RIOK3, UBA5 knockdown | `flow_cytometry/fig3H_figS5I_riok3_uba5.py` | Zenodo `flow_RIOK3-UBA5-knockdown` |

## Figure 4
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| 4A | RPL, RPS protein levels after DDRGK1 knockdown | `proteomics/total_abundance.py` | `proteomics/input/all_ttest_results.csv`; OpenCell mmc4 |
| 4B | Old 60S by serum and confluence | `flow_cytometry/fig4BD_serum_confluence_gne.py` | Zenodo `flow_confluence-GNE-7883` |
| 4C | Hippo pathway and confluence | No code: schematic | — |
| 4D | Old 60S with GNE-7883 at low confluence | `flow_cytometry/fig4BD_serum_confluence_gne.py` | Zenodo `flow_confluence-GNE-7883` |
| 4E | Translocon-dependent protein levels after DDRGK1 knockdown | `proteomics/total_abundance.py` | `proteomics/input/all_ttest_results.csv`, `proteomics/input/total_topology.tsv`; OpenCell mmc4 |

## Figure 5
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| 5A | Sucrose gradient, sgUBA5 vs sgNT | No code: sucrose gradient traces analysed manually | — |
| 5B | Single-molecule imaging schematic and raw snapshots | `single_molecule/figures/build_fig5B_snapshot.py --render` (snapshots; frames in `single_molecule/figures/fig5B_snapshot_spec.json`); schematic drawn by hand | Zenodo `SPT_60S_UBA5-knockdown_2025-06-06` |
| 5C | 60S diffusion distribution, sgUBA5 vs sgNT | `single_molecule/figures/build_fig5C_figS7B_panels.py` | Zenodo `SPT_60S_UBA5-knockdown_2025-06-06` |
| 5D | PT-179 degron; RPL26-UFM1 loss | Left: no code, schematic. Right: `gels/fig5D_figS7F_pt179_ufm1.py` | Zenodo `gel_PT-179_UFM1-immunoblot` |
| 5E | ER co-imaging and trajectory classification | `single_molecule/pipeline/build_er_costain_figure.py line=L29_degron=000_HHT=000_001`; schematic drawn by hand | Zenodo `SPT_60S_ER-costain_UBA5-degron_2025-09-15` |
| 5F | On-ER vs off-ER 60S diffusion, PT-179 | Left: `single_molecule/figures/build_fig5F_figS7J_panels.py`. Right: `single_molecule/figures/run_zone_refit_batch.py` → `single_molecule/figures/aggregate_zone_refit.py` → `single_molecule/figures/replot_refit_panels.py` | Zenodo `SPT_60S_ER-costain_UBA5-degron_2025-09-15` |
| 5G | LOV-BirA proximity labeling, flow cytometry | `flow_cytometry/fig5G_lov_birA.py` (gating from FCS files: `flow_cytometry/fcs_gating.py`); schematic drawn by hand | Zenodo `flow_LOV-BirA-ER-labeling` |

## Figure 6
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| 6A | Dual-guide suppressor screen and epistasis logic | No code: schematic | — |
| 6B | Degradation scores, sgUBA5 vs sgNT screens | `screen_figures/panels/fig6B_figS8BC_suppressor_scatter.py` | GEO sublibrary counts; `screen_figures/categories/` |
| 6C | BafA1-dependent acidic RPL29-mKeima signal | `flow_cytometry/fig6C_figS8D_mkeima.py`; schematic drawn by hand | Zenodo `flow_RPL29-mKeima-BafA1` |

## Figure 7
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| 7A | Old 60S after UBA5, TLSD1 knockdown | `flow_cytometry/fig7A_figS9BD_uba5_tlsd1.py` | Zenodo `flow_TLSD1-UBA5-knockdown` |
| 7B | ER-labeled ribosome fraction, days 0 and 5 | `gels/fig7B_figS9FGI_lov_birA.py` | Zenodo `gel_LOV-BirA_ER-labeling` |
| 7C | Day-8 new-protein volcano, TLSD1 highlighted | `proteomics/silac_turnover.py` | `proteomics/input/protein_quant_30972.tsv`; OpenCell mmc4 |
| 7D | GFP-TLSD1 and Sec61β colocalization | Right: `confocal_colocalization/segment_cells.py` → `confocal_colocalization/wholecell_pearson.py` → `confocal_colocalization/plot_pearson_swarm.py` (driver `confocal_colocalization/submit_pearson_repro.sh`). Left (inset): `confocal_colocalization/build_fig7D_figS10_outline_stacks.py`, assembled in Fiji | Zenodo `confocal_GFP-TLSD1_ER_PT-179` |
| 7E | Sucrose gradients, TLSD1 and UBA5 knockdown | No code: sucrose gradient traces analysed manually | — |
| 7F | Ricin screen waterfall, TLSD1 highlighted | Top: no code, schematic. Bottom: `ricin_reanalysis/run_all.sh`, then `ricin_reanalysis/fig7F_ricin_waterfall.py` | Morgens 2017 SD4; the `run_all.sh` output `K562_vs_Ricin_all_genes.csv` is also in Zenodo `ricin-screen-reanalysis_Morgens2017` |
| 7G | Model of TLSD | No code: schematic | — |

## Figure S1
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S1A | Quiescence induction; EdU flow cytometry | No code: schematic; EdU plots drawn in FlowJo | — |
| S1B | Immunoblots of HaloTag knock-in clones | No code: representative immunoblots | — |
| S1C | New/old ratio histograms, screen populations | No code: drawn in FlowJo | — |
| S1D | Gene-score correlations across screens | `screen_figures/panels/fig1_figS1_figS2A_figS6B.py --panels S1D` | GEO genome-wide counts |
| S1E | Growth scores, common essentials highlighted | `screen_figures/panels/fig1_figS1_figS2A_figS6B.py --panels S1E` | GEO genome-wide counts; DepMap |
| S1F | New/old scores, common essentials highlighted | `screen_figures/panels/fig1_figS1_figS2A_figS6B.py --panels S1F` | GEO genome-wide counts; DepMap |

## Figure S2
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S2A | New/old scores by gene category | `screen_figures/panels/fig1_figS1_figS2A_figS6B.py --panels S2A` | GEO genome-wide counts; `screen_figures/categories/` |
| S2B | Autophagy reporter flux after ATG7, RB1CC1 knockdown | `flow_cytometry/figS2B_autophagy_reporter.py`; schematic drawn by hand | Zenodo `flow_autophagy-reporter` |

## Figure S3
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S3A | Simulated synthesis vs degradation changes | `simulation/figS3A_turnover_simulation.py` (model: `simulation/biomolecule_turnover.py`) | None (simulation) |
| S3B | Genes selected for the sublibrary | `screen_figures/panels/figS3B_sublibrary_selection.py` | GEO genome-wide counts; GEO sublibrary counts |
| S3C | sgRNA-level score construction | `screen_figures/panels/figS3C_score_construction.py` | GEO sublibrary counts |
| S3D | Genome-wide vs sublibrary fold-change correlation | `screen_figures/panels/figS3D_figS8A_lfc_correlation.py` | GEO genome-wide counts; GEO sublibrary counts |
| S3E | GO terms for synthesis, degradation scores | `screen_figures/panels/figS3E_go_scores.py` | GEO sublibrary counts; Enrichr GO 2025; `screen_figures/categories/` |
| S3F | Growth vs synthesis, degradation score correlation | `screen_figures/panels/figS3F_growth_vs_scores.py` | GEO sublibrary counts |

## Figure S4
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S4A | 40S vs 60S synthesis, degradation scores | `screen_figures/panels/figS4A_subunit_scatter.py` | GEO sublibrary counts; `screen_figures/categories/` |
| S4B | Synthesis scores by gene category | `screen_figures/panels/fig2D_figS4B_category_violins.py` | GEO sublibrary counts; `screen_figures/categories/` |

## Figure S5
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S5A | Proteomics datasets collected | No code: schematic | — |
| S5B | Old and new protein intensity by genotype | `proteomics/silac_turnover.py` | `proteomics/input/protein_quant_30972.tsv`; OpenCell mmc4 |
| S5C | Correlation of proteome fold changes | `proteomics/figS5C_correlation.py` | `proteomics/input/protein_quant_30972.tsv`, `proteomics/input/all_ttest_results.csv`; OpenCell mmc4 |
| S5D | Old-protein change by compartment | `proteomics/silac_turnover.py` | `proteomics/input/protein_quant_30972.tsv`; OpenCell mmc4 |
| S5E | New-protein change by compartment | `proteomics/silac_turnover.py` | `proteomics/input/protein_quant_30972.tsv`; OpenCell mmc4 |
| S5F | New and old ER protein categories | `proteomics/silac_turnover.py` | `proteomics/input/protein_quant_30972.tsv`; OpenCell mmc4 |
| S5G | RPL26 immunoblots of pulldown fractions | No code: representative immunoblots | Raw scans: Zenodo `gel_RPL26-3xFLAG_pulldown` |
| S5H | Ribosome pellet new/old ratio before pulldown | `gels/fig3F_figS5H_rpl26_pulldown.py` | Zenodo `gel_RPL26-3xFLAG_pulldown` |
| S5I | New 40S, 60S after RIOK3, UBA5 knockdown | `flow_cytometry/fig3H_figS5I_riok3_uba5.py` | Zenodo `flow_RIOK3-UBA5-knockdown` |

## Figure S6
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S6A | Protein categories after DDRGK1 knockdown | `proteomics/total_abundance.py` (`fig4A_figS6A_*`) | `proteomics/input/all_ttest_results.csv`; OpenCell mmc4 |
| S6B | Quiescent growth scores, Hippo pathway highlighted | `screen_figures/panels/fig1_figS1_figS2A_figS6B.py --panels S6B` | GEO genome-wide counts; `screen_figures/categories/` |
| S6C | Membrane-protein categories after DDRGK1 knockdown | `proteomics/total_abundance.py` | `proteomics/input/all_ttest_results.csv`, `proteomics/input/total_topology.tsv`; OpenCell mmc4 |
| S6D | Volcano, sgDDRGK1 vs sgNT, quiescent | `proteomics/total_abundance.py` | `proteomics/input/all_ttest_results.csv`; OpenCell mmc4 |
| S6E | Volcano, sgDDRGK1 vs sgNT, growing | `proteomics/total_abundance.py` | `proteomics/input/all_ttest_results.csv`; OpenCell mmc4 |

## Figure S7
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S7A | Wasserstein distance, sgUBA5 vs sgNT, 60S, 40S | `single_molecule/figures/run_panel_refit_batch.py` → `single_molecule/figures/aggregate_panel_refit.py` → `single_molecule/figures/replot_refit_panels.py` | Zenodo `SPT_60S_UBA5-knockdown_2025-06-06`, `SPT_40S_UBA5-knockdown_2025-06-06` |
| S7B | 40S diffusion distribution, sgUBA5 vs sgNT | `single_molecule/figures/build_fig5C_figS7B_panels.py` | Zenodo `SPT_40S_UBA5-knockdown_2025-06-06` |
| S7C | 60S diffusion-bin differences | `single_molecule/figures/run_panel_refit_batch.py` → `single_molecule/figures/aggregate_panel_refit.py` → `single_molecule/figures/replot_refit_panels.py` | Zenodo `SPT_60S_UBA5-knockdown_2025-06-06` |
| S7D | Diffusion bins after harringtonine or emetine | `single_molecule/figures/run_panel_refit_batch.py` → `single_molecule/figures/aggregate_panel_refit.py` → `single_molecule/figures/replot_refit_panels.py`; schematic drawn by hand | Zenodo `SPT_60S_harringtonine-emetine_2026-07-08` |
| S7E | Translation-state assignment and drug overlays | `single_molecule/figures/build_figS7DE_panels.py` (table); `single_molecule/figures/replot_refit_panels.py` (overlays) | Zenodo `SPT_60S_harringtonine-emetine_2026-07-08` |
| S7F | SD40-UBA5 degradation efficiency | `gels/fig5D_figS7F_pt179_ufm1.py` (graph); blots are representative scans, no code | Zenodo `gel_PT-179_UFM1-immunoblot` |
| S7G | Ribosome levels during PT-179 time course | `flow_cytometry/figS7G_degron_timecourse.py` | Zenodo `flow_UBA5-degron_PT-179-timecourse` |
| S7H | ER signal vs distance to mask boundary | `single_molecule/pipeline/build_er_boundary_profile_qc.py` | Zenodo `SPT_60S_ER-costain_UBA5-degron_2025-09-15` |
| S7I | ER mask area per field of view | `single_molecule/pipeline/build_er_area_fraction_strip.py` | Zenodo `SPT_60S_ER-costain_UBA5-degron_2025-09-15` |
| S7J | Total 60S diffusion, ignoring ER mask | `single_molecule/figures/build_fig5F_figS7J_panels.py` | Zenodo `SPT_60S_ER-costain_UBA5-degron_2025-09-15` |

## Figure S8
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S8A | Fold-change correlation, sgUBA5 vs sgNT screens | `screen_figures/panels/figS3D_figS8A_lfc_correlation.py` | GEO sublibrary counts |
| S8B | 60S synthesis scores, sgUBA5 vs sgNT | `screen_figures/panels/fig6B_figS8BC_suppressor_scatter.py` | GEO sublibrary counts; `screen_figures/categories/` |
| S8C | 40S synthesis scores, sgUBA5 vs sgNT | `screen_figures/panels/fig6B_figS8BC_suppressor_scatter.py` | GEO sublibrary counts; `screen_figures/categories/` |
| S8D | Neutral RPL29-mKeima signal with BafA1 | `flow_cytometry/fig6C_figS8D_mkeima.py` | Zenodo `flow_RPL29-mKeima-BafA1` |
| S8E | Growth scores, sgUBA5 vs sgNT screens | `screen_figures/panels/figS8EF_growth_epistasis.py` | GEO sublibrary counts; `screen_figures/categories/` |
| S8F | Growth epistasis waterfall | `screen_figures/panels/figS8EF_growth_epistasis.py` | GEO sublibrary counts; `screen_figures/categories/` |

## Figure S9
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S9A | AlphaFold model of TLSD1 | No code: AlphaFold model render (ChimeraX) | — |
| S9B | Old 40S after UBA5, TLSD1 knockdown | `flow_cytometry/fig7A_figS9BD_uba5_tlsd1.py` | Zenodo `flow_TLSD1-UBA5-knockdown` |
| S9C | New 60S after UBA5, TLSD1 knockdown | `flow_cytometry/fig7A_figS9BD_uba5_tlsd1.py` | Zenodo `flow_TLSD1-UBA5-knockdown` |
| S9D | New 40S after UBA5, TLSD1 knockdown | `flow_cytometry/fig7A_figS9BD_uba5_tlsd1.py` | Zenodo `flow_TLSD1-UBA5-knockdown` |
| S9E | LOV-BirA gels and immunoblots | No code: representative gel and immunoblot images | Raw scans: Zenodo `gel_LOV-BirA_ER-labeling` |
| S9F | Total JF549-labeled RPL29-HaloTag per tubulin | `gels/fig7B_figS9FGI_lov_birA.py` | Zenodo `gel_LOV-BirA_ER-labeling` |
| S9G | Biotinylated RPL29-HaloTag-Avi per tubulin | `gels/fig7B_figS9FGI_lov_birA.py` | Zenodo `gel_LOV-BirA_ER-labeling` |
| S9H | HaloTag lysosomal delivery assay | No code: schematic | — |
| S9I | Lysosomal HaloTag band ratio (60S ribophagy) | `gels/fig7B_figS9FGI_lov_birA.py` | Zenodo `gel_LOV-BirA_ER-labeling` |

## Figure S10
| Panel | What it shows | Script | Inputs |
|---|---|---|---|
| S10 (single panel) | GFP-TLSD1 and Sec61β, full confocal view | `confocal_colocalization/build_fig7D_figS10_outline_stacks.py` (outline stacks), assembled in Fiji | Zenodo `confocal_GFP-TLSD1_ER_PT-179` |

## Supplementary tables
| Table | What it shows | Script | Inputs |
|---|---|---|---|
| S1 | Genome-wide screens, sgRNA-level data | `crispr_screens/tables/build_supplementary_tables.py` (genome-wide config) | GEO genome-wide counts |
| S2 | Genome-wide screens, gene-level data | `crispr_screens/tables/build_supplementary_tables.py` (genome-wide config) | GEO genome-wide counts; `screen_figures/categories/` (includes DepMap) |
| S3 | GO enrichment, genome-wide screens | `crispr_screens/tables/build_supplementary_tables.py` (genome-wide config) | GEO genome-wide counts; Enrichr GO 2025 |
| S4 | Sublibrary sgRNA-level data and scores | `crispr_screens/tables/build_sublibrary_tables.py` | GEO sublibrary counts |
| S5 | Sublibrary gene-level data and scores | `crispr_screens/tables/build_sublibrary_tables.py` | GEO sublibrary counts |
| S6 | GO enrichment, synthesis and degradation scores | `crispr_screens/go/build_tableS6.py` | GEO sublibrary counts; Enrichr GO 2025 |
| S7 | sgRNA spacers and qPCR validation | No code: reagent table | — |
| S8 | SILAC-TMT new/old protein data | `proteomics/silac_turnover.py` | `proteomics/input/protein_quant_30972.tsv`; OpenCell mmc4 |
| S9 | Total-abundance TMT data with annotations | `proteomics/total_abundance.py` | `proteomics/input/all_ttest_results.csv`, `proteomics/input/total_topology.tsv`; OpenCell mmc4; MitoCarta 3.0 |
| S10 | sgUBA5 sublibrary sgRNA-level data and scores | `crispr_screens/tables/build_sublibrary_tables.py` | GEO sublibrary counts |
| S11 | sgUBA5 sublibrary gene-level data and scores | `crispr_screens/tables/build_sublibrary_tables.py` | GEO sublibrary counts |
| S12 | Single-particle tracking analysis settings | `single_molecule/figures/build_settings_table.py` | Zenodo `SPT_60S_UBA5-knockdown_2025-06-06`, `SPT_40S_UBA5-knockdown_2025-06-06`, `SPT_60S_harringtonine-emetine_2026-07-08`, `SPT_60S_ER-costain_UBA5-degron_2025-09-15`; refit statistics of Fig S7A, S7C, S7D and 5F |
