# Code for: Failure to release from the ER translocon induces 60S ribosome degradation

Levine *et al.* (in preparation). Analysis and figure code.

> Work in progress: this repository is being assembled. See `docs/CONVENTIONS.md`.

## Repository map
| Folder | Contents |
|---|---|
| `crispr_screens/` | sgRNA counting → MAGeCK → gene scores → synthesis/degradation line scores → per-sgRNA p-values → GO → supplementary tables |
| `screen_figures/` | Plotting library and one script per screen-derived figure panel |
| `ricin_reanalysis/` | Reanalysis of Morgens *et al.* 2017 ricin screen (Fig 7F) |
| `single_molecule/` | Single-particle tracking (quot + saSPT), ER masking, resampling statistics (Fig 5, S7, Table S12) |
| `confocal_colocalization/` | Cellpose segmentation and per-cell colocalization (Fig 7D, S10) |
| `gels/` | In-gel fluorescence and immunoblot quantification (Fig 3F, 5D, 7B; S5H, S7F, S9F, S9G, S9I) |
| `proteomics/`, `flow_cytometry/`, `simulation/` | Remaining quantitative panels |
| `common/` | Path configuration shared by all scripts |
| `envs/` | Conda environment files |

The sucrose gradient traces (Fig 5A, 7E) were analysed manually; no code accompanies them.

## Setup
1. Create the conda environments in `envs/`.
2. Copy `config/paths.example.yaml` to `config/paths.yaml`, and point `data` and `out` at the data
   directory (see Data availability).
3. Place the inputs under `data:` at the relative paths the scripts read (the module READMEs and the Zenodo
   README list them).

## Data availability
- **Code:** this repository, archived on Zenodo as a Software record (DOI TBD).
- **CRISPRi screens: GEO (accession TBD).** The FASTQs and the per-sample count tables
  (`genomewide_aggregated_counts.csv`, `sublibrary_aggregated_counts.csv`). The screen results and gene
  scores are Supplementary Tables S1–S6, S10 and S11. `crispr_screens/pipeline/README.md` describes how
  to place the count tables.
- **Proteomics: ProteomeXchange/PRIDE (TBD).** Raw files and the processed quantification tables; the
  per-protein results are Tables S8 and S9 (see `proteomics/README.md`).
- **Zenodo (TBD).** Single-particle tracking inputs and resampling draws, the confocal stacks and masks,
  flow-cytometry exports, gel quantification sheets and raw scans, and the ricin reanalysis output.
- **Third-party downloads** (not redistributed):
  - Enrichr `GO_Biological_Process_2025`, `GO_Molecular_Function_2025` and `GO_Cellular_Component_2025`;
  - the OpenCell organelle-IP table (ref. 170, supplementary file mmc4);
  - Morgens *et al.* 2017 Supplementary Data 4;
  - the DepMap Achilles common essentials;
  - MitoCarta 3.0.

## License
MIT; see `LICENSE`. Copyright (c) 2026 Whitehead Institute for Biomedical Research. Third-party tools
(MAGeCK, saSPT, quot, Cellpose and others) are installed separately under their own licenses and are not
included here.

## AI disclosure
Claude Opus models 4–5 were used in the generation of analysis and plotting code during this work. All
workflows were generated as directed by humans, and the code was verified for accuracy.
