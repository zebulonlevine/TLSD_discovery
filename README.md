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
| `proteomics/`, `flow_cytometry/`, `gels/`, `sucrose_gradients/`, `simulation/` | Remaining quantitative panels |
| `common/` | Path configuration shared by all scripts |
| `envs/` | Conda environment files |

## Setup
1. Create the conda environments in `envs/`.
2. Copy `config/paths.example.yaml` to `config/paths.yaml`, and point `data` and `out` at the data
   directory (see Data availability).

## Data availability
- CRISPRi screen sequencing: GEO (accession TBD).
- Proteomics: ProteomeXchange/PRIDE (TBD).
- Single-particle trajectories: Zenodo (TBD).

## AI disclosure
Claude Opus models 4–5 were used in the generation of analysis and plotting code during this work. All
workflows were generated as directed by humans, and the code was verified for accuracy.
