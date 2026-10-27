# Single-particle tracking (Fig 5B, 5C, 5E, 5F; Supplementary Fig S7A–E, S7H–J; Table S12)

Code for the single-molecule panels of *Failure to release from the ER translocon induces 60S ribosome
degradation*.

**Pipeline:**
1. quot tracking;
2. saSPT state-array likelihoods on a 36-state grid;
3. per-field-of-view (FOV) VBEM occupancies;
4. panels and resampling statistics.

**Further documentation:**
- Methods detail: `docs/METHODS.md`.
- Stage diagram: `docs/SPT_ANALYSIS_FLOWCHART.md`.
- Megapools: `docs/MEGAPOOL.md`.

**Origin.** Imported from the lab analysis tree at commit `8c70d70`. The tracking/saSPT workflow in
`pipeline/` builds on G. Schlissel's `quot_saspt_workflow` (gschliss/quot_saspt_workflow).

**Dependencies.** quot (A. Heckert, commit `1b9051e`) and saSPT 0.4.0. Both are installed from upstream
(`envs/SPT_env.yml`), not vendored.

## Data
The inputs are on Zenodo: DOI [10.5281/zenodo.23688082](https://doi.org/10.5281/zenodo.23688082).
That deposit has:
- trajectories, pooled per condition (`traj_combinations/`);
- ER masks;
- the settings of the 36-state likelihood caches (`occupancy_compare_cache_36state/settings_36state.pkl`).
  The caches and their FOV tables are not deposited; step 2 of the run order rebuilds them;
- the refit resampling draws.

Point `TLSD_SPT_DATA` at the unpacked tree, which holds `all_data/<experiment>/...`.

**Subunit labels (2026_07_08).** The 2026_07_08 movies were acquired with the two tag labels swapped:
raw `halo=S3` is RPL29 (60S) and raw `halo=L29` is RPS3 (40S), confirmed from the constructs
(2026-08-03). The deposited data use the true identity. The 60S arm is renamed `halo=60S_...` (conditions)
and `60S_*` (combinations, caches), and the raw-filename typo `UBA5drug=eme360` is normalised to
`UBA5_drug=eme360`. The 40S arm is used by no panel and is not deposited.

Three experiments feed the panels:

| Experiment | Tier | Panels |
|---|---|---|
| `2025_06_06_zeb` | `tracking_output_q=15p0_splitsize10` | 5B, 5C, S7A–C |
| `2026_07_08` | `tracking_output_q=15p0` | S7D, S7E |
| `2025_09_15_zeb_ER` | `tracking_output_q=15p0` | 5E, 5F, S7H–J |

## Layout
The three folders mirror the three code trees the analysis ran from, so every script imports the same
module versions it ran with:

| Folder | Origin | Contents |
|---|---|---|
| `pipeline/` | quot/saSPT workflow | tracking (`run_filewise.py` …), condition/combination pooling, likelihood caches, ER zone caches, 36-state cache builder, ER panel builders |
| `er/` | ER-SPT workflow | ER preprocessing, illumination correction, ER zoning, ER-specific combinationwise, `viz_tools.py`. Its `combination_functions.py` and `run_combinationwise.py` add per-combination cutoff overrides and shadow the `pipeline/` copies. |
| `figures/` | figure builders | panel builders, refit (resampling) chain, Table S12 builder, shared plotting style |
| `validation/` | — | focal-depth calibration and saSPT-equivalence checks cited in `docs/METHODS.md` |

## Panels
Outputs go under `$TLSD_OUT/PaperFigures_updated/`, one directory per panel. Panel files from the
36-state fits carry the `_36state` suffix.

| Panel | Script | Output |
|---|---|---|
| 5B | `figures/build_fig5B_snapshot.py --render` (FOVs/frames/crops: `figures/fig5B_snapshot_spec.json`) ¹ | `figure5/B/fig5B_snapshot.pdf` |
| 5C | `figures/build_fig5C_figS7B_panels.py` | `figure5/C/fig5C_60S_posterior_overlay_36state.pdf` |
| 5E | `pipeline/build_er_costain_figure.py line=L29_degron=000_HHT=000_001` ² | `figure5/E/costain_line=L29_degron=000_HHT=000_001_T300.pdf` |
| 5F left | `figures/build_fig5F_figS7J_panels.py` | `figure5/F/fig5F_left_zone_posterior_overlay_36state.pdf` |
| 5F right | refit chain: `run_zone_refit_batch.py` → `aggregate_zone_refit.py` → `replot_refit_panels.py` | `figure5/F/fig5F_right_zone_bin_effect_refit_36state.pdf`, `fig5F_refit_stats_36state.csv` |
| S7A | refit chain: `run_panel_refit_batch.py` → `aggregate_panel_refit.py` → `replot_refit_panels.py` | `figureS7/A/figS7A_W1_60S_40S_refit_36state.pdf`, `figS7A_W1_refit_stats_36state.csv` |
| S7B | `figures/build_fig5C_figS7B_panels.py` | `figureS7/B/figS7B_40S_posterior_overlay_36state.pdf` |
| S7C | refit chain (as S7A) | `figureS7/C/figS7C_60S_bin_difference_refit_36state.pdf` and its stats CSV |
| S7D | refit chain (as S7A) | `figureS7/D/figS7D_drug_shift_both_genotypes_refit_36state.pdf`, `figS7D_refit_stats_36state.csv` |
| S7E | `figures/replot_refit_panels.py` (overlays) | `figureS7/E/figS7E_overlay_{sgNT,sgUBA5}_refit_36state.pdf` |
| S7E table | `figures/build_figS7DE_panels.py` | `figureS7/E/figS7E_translation_state_assignment_36state.csv` |
| S7H | `pipeline/build_er_boundary_profile_qc.py` | `figureS7/H/er_boundary_profile_qc.pdf` |
| S7I | `pipeline/build_er_area_fraction_strip.py` | `figureS7/I/er_area_fraction_strip.pdf`, `er_area_fraction_stats.csv` |
| S7J | `figures/build_fig5F_figS7J_panels.py` | `figureS7/J/figS7J_total_posterior_overlay_36state.pdf` |
| Table S12 | `figures/build_settings_table.py` | `supplementary/TableS12_SPT_settings.{csv,md}` |

**Fits.** Every panel fits each condition independently.

**Refit chain (S7A, S7C, S7D, 5F right).**
- The refit runs saSPT again on every resampling draw. It uses 500 FOV-bootstrap and 500 FOV-permutation
  draws, so the smallest permutation p is 1/501 (≈ 0.002).
- These panels therefore carry at most **, and the stats CSVs flag p-values at that floor (`p_at_floor`).
- The 5F zone-label permutation needs no refit; it uses 2000 draws, and the treatment permutation uses
  10 000.

**No-refit companion outputs.** `build_figS7DE_panels.py` and `build_fig5F_figS7J_panels.py` also write no-refit
versions of the S7D, S7E and 5F-right statistics, using 2000 resamples without refitting. These are not
the placed panels.

**Overlay sidecars.** Every density overlay writes `_curves.csv` and `_plotspec.json` sidecars, and
`figures/replot_overlays.py` redraws an overlay from them without refitting.

¹ **Fig 5B.** `fig5B_snapshot_spec.json` selects the two raw frames placed in 5B, both NT-guide (sgNT):
- left: FOV 007, frame 450;
- right: FOV 004, frame 110.

Each is an 80-px crop, chosen from the candidates `--scout` ranks. The panel was drafted as sgNT vs sgUBA5,
but the right-hand crop is an NT-guide FOV, so the placed panel carries no condition labels (2026-10-26).
Tracks are drawn in the magenta of the placed panel (`#EC008C`). The panel was first rendered in yellow
(`#ffe600`) and recoloured by hand. The following were also done by hand:
- the legend ("past trajectory", "current position");
- the track line weight;
- the scale-bar label size.

² **Fig 5E.** The placed panels come from the FOV 001 co-stain figure:
- the full-FOV merge (row 0, column 3);
- the zoom-1 merge, mask + tracks and classification panels (row 1, columns 3–5).

The classification colours are the published ones: on ER `#23AAE1`, off ER `#ED1C24`, undetermined
(intermediate) `#662D91`. `ambiguous` is omitted. These colours were applied by hand to the original render,
so they are set in this builder only; `er/viz_tools.py` keeps its shared palette (lab commit `a3d40f7`).

Other edits were made by hand in the figure and are not in the code:
- the titles, white marker outlines and cyan zoom-2 box are removed;
- the classification strokes are heavier;
- the scale bars are relabelled.

## Environment
Use `envs/SPT_env.yml`. Set these variables, or the matching keys in `config/paths.yaml`:

| Variable | Meaning |
|---|---|
| `TLSD_REPO` | repository root (used by the Slurm wrappers) |
| `TLSD_SPT_DATA` | SPT data root: holds `all_data/<experiment>/...` (movies, tracking tiers, likelihood caches) |
| `TLSD_SPT_DATA_<experiment>` | optional: a separate root (holding `all_data/<experiment>/...`) for one experiment; config key `spt_data_<experiment>` |
| `TLSD_OUT` | output root |
| `PANEL_VARIANT` | `36state` for every panel |
| `TLSD_PYTHON` | interpreter used by the Slurm wrappers (default `python3`) |
| `TLSD_SLURM_PARTITION`, `TLSD_SLURM_ACCOUNT` | optional, for the jobs `er/batch_er_workflow.py` submits |

**Slurm.** The wrappers carry no partition or account; pass them to `sbatch`
(e.g. `sbatch -p <partition> ...`). The panel builders take 5–20 min and up to ~40 GB RAM each.

**Writes into the data tree.** Steps 3–6 only read the data tree and write under `TLSD_OUT`. The joint caches
of step 3 go to `$TLSD_OUT/spt_joint_caches/<experiment>/<tier>/`. Steps 1 and 2 regenerate inputs from the
raw movies and write into the experiment's own tree (tracking outputs, per-condition caches).

## Run order
With the deposit, run steps 2, 3, 4 and 6.
- Step 2 builds the 12 per-condition likelihood caches from the deposited pooled trajectories.
- Step 3 concatenates them into the joint caches.
- Neither set of caches is deposited. Both rebuilds are bit-identical to the caches behind the published
  panels: every array, the dtypes and the FOV tables; only the stored input paths differ.
- Steps 1 and 5 regenerate the trajectories and the resampling draws from the raw movies.
- The deposit holds trajectories pooled per condition (`traj_combinations/`), not per field of view.

1. **Tracking.**
   - Run `pipeline/call_filewise_operations.slurm` (`run_filewise.py`). For the ER data, run
     `er/batch_er_workflow.py`, which runs `run_filewise_er.py`.
   - Then run `call_conditionwise_operations.slurm` and `call_combinationwise_operations.slurm`.
     - ER zones: `er/generate_er_combinations.py` + `er/combo_array.sh`.
     - Re-zoning: `er/rezone_array.sh --pct-floor 60 --require-broad false --min-frames 2`.
   - Then run `call_megapool_operations.slurm` (2026_07_08).
   - The 2026_07_08 raw movies carry the 60S names (see Subunit labels); that experiment's
     `combinations.yaml` maps `60S_*` to the `halo=60S_...` conditions.
   - **Background level and frame interval.** `run_filewise.py` takes `camera_bg` and `dt` for a condition
     from one raw movie.
     - The published runs used the first movie `glob` returned, i.e. directory order.
     - `pipeline/image_metadata_sources.csv` records that movie and its values for each experiment and
       condition, so re-tracking reproduces the published trajectories on any copy, even without that
       movie. Verified byte-identical: 6 FOVs of 2026_07_08 (one per condition) and 5 of 2025_06_06_zeb.
     - A condition not in the table uses the first movie in sorted order.
     - The ER experiment (`run_filewise_er.py`) reads `camera_bg` from each FOV's own SPT TIF and does not
       use the table.
   - **Pooling order.** Pooling still takes the per-FOV files in `os.listdir` order, so a copied tree can
     pool the FOVs in a different order. The trajectories of each FOV are unchanged.
2. **36-state likelihood caches.**
   - Run `pipeline/submit_36state_caches.py`. It submits one Slurm job per condition, with memory sized
     per condition: `call_variant_cache.slurm` → `build_variant_cache_one.py`.
   - Or run one condition directly:
     `python pipeline/build_variant_cache_one.py --dataset <experiment> --condition <condition> --variant 36state`.
   - Conditions:
     - `2025_06_06_zeb`: `sgNT_60S`, `sgUBA5_60S`, `sgNT_40S`, `sgUBA5_40S`;
     - `2026_07_08`: `60S_{NT,UBA5}_{none,HT10,eme360}`;
     - `2025_09_15_zeb_ER`: `noDegron`, `Degron`.
   - Measured on 4 CPUs: 3–15 min and 11–36 GB peak per condition. Run in sequence, that is about
     20 min for 2025_06_06_zeb, 60 min for 2026_07_08 and 10 min for 2025_09_15_zeb_ER.
   - Output: `<condition>_likelihood_cache.npz` and `<condition>_fov_metadata.csv` in the tier's
     `occupancy_compare_cache_36state/` (36 GB in total).
   - The Fig 5F builder also builds the two 2025_09_15_zeb_ER caches itself if they are missing.
3. **Joint caches for the permutation nulls.** For each contrast in `figures/submit_panel_refit.py`, run
   `figures/merge_contrast_cache.py --tier <tier> --cond-a A --cond-b B`. Contrast keys name their panel
   (`figS7AC_60S`, `figS7A_40S`, `figS7D_*`).
4. **Panels.**
   - `sbatch pipeline/call_build_panels.slurm build_fig5C_figS7B_panels.py`, and the same for
     `build_figS7DE_panels.py` and `build_fig5F_figS7J_panels.py`.
   - `pipeline/build_er_boundary_profile_qc.py`, `pipeline/build_er_area_fraction_strip.py`,
     `pipeline/build_er_costain_figure.py line=L29_degron=000_HHT=000_001`.
   - `sbatch pipeline/call_build_panels.slurm build_fig5B_snapshot.py --render`.
5. **Resampling draws.**
   - `figures/submit_panel_refit_all.py` (→ `run_panel_refit_batch.py`).
   - For Fig 5F, `run_zone_refit_batch.py` for each condition and scheme.
   - `figures/fill_refit_gaps.py` resubmits any missing batch.
6. **Statistics and refit panels.**
   - `figures/aggregate_panel_refit.py --key <key>` for each contrast in `submit_panel_refit.CONTRASTS`.
   - Then `figures/write_refit_sidecars.py`, `figures/aggregate_zone_refit.py` and
     `figures/replot_refit_panels.py --all`.
   - Last, `figures/build_settings_table.py`, which builds Table S12 from the refit statistics.
