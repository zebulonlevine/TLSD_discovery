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
The inputs are on Zenodo: DOI **TBD**. That deposit has:
- trajectories;
- ER masks;
- the 36-state likelihood caches (`occupancy_compare_cache_36state/`);
- the refit resampling draws.

Point `TLSD_SPT_DATA` at the unpacked tree, which holds `all_data/<experiment>/...`.

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
Outputs go under `$TLSD_OUT`. Panel files carry the `_36state` suffix. The ER builders write to
`PaperFigures/figure4_claim3/`; everything else writes to `PaperFigures_updated/`.

| Panel | Script | Output |
|---|---|---|
| 5B | `figures/build_fig5B_snapshot.py --render` (FOVs/frames/crops: `figures/fig5B_snapshot_spec.json`) ¹ | `figure5/B/fig5B_snapshot.pdf` |
| 5C | `figures/build_claim1_panels.py` | `figure5/C/fig5C_60S_posterior_overlay_36state.pdf` |
| 5E | `pipeline/build_er_costain_figure.py line=L29_degron=000_HHT=000_001` ² | `costain/costain_line=L29_degron=000_HHT=000_001_T300.pdf` |
| 5F left | `figures/build_claim3_panels.py` | `figure5/F/fig5F_left_zone_posterior_overlay_36state.pdf` |
| 5F right | refit chain: `run_zone_refit_batch.py` → `aggregate_zone_refit.py` → `replot_refit_panels.py` | `figure5/F/fig5F_right_zone_bin_effect_refit_36state.pdf`, `fig5F_refit_stats_36state.csv` |
| S7A | refit chain: `run_panel_refit_batch.py` → `aggregate_panel_refit.py` → `replot_refit_panels.py` | `figureS7/A/figS7A_W1_60S_40S_refit_36state.pdf`, `figS7A_W1_refit_stats_36state.csv` |
| S7B | `figures/build_claim1_panels.py` | `figureS7/B/figS7B_40S_posterior_overlay_36state.pdf` |
| S7C | refit chain (as S7A) | `figureS7/C/figS7C_60S_bin_difference_refit_36state.pdf` and its stats CSV |
| S7D | refit chain (as S7A) | `figureS7/D/figS7D_drug_shift_both_genotypes_refit_36state.pdf`, `figS7D_refit_stats_36state.csv` |
| S7E | `figures/replot_refit_panels.py` (overlays) | `figureS7/E/figS7E_overlay_{sgNT,sgUBA5}_refit_36state.pdf` |
| S7E table | `figures/build_claim2_panels.py` | `figureS7/E/figS7E_translation_state_assignment_36state.csv` |
| S7H | `pipeline/build_er_boundary_profile_qc.py` | `er_boundary_profile_qc.pdf` |
| S7I | `pipeline/build_er_area_fraction_strip.py` | `er_area_fraction_strip.pdf`, `er_area_fraction_stats.csv` |
| S7J | `figures/build_claim3_panels.py` | `figureS7/J/figS7J_total_posterior_overlay_36state.pdf` |
| Table S12 | `figures/build_settings_table.py` | `supplementary/TableS12_SPT_settings.{csv,md}` |

**Fits.** Every panel fits each condition independently.

**Refit chain (S7A, S7C, S7D, 5F right).**
- The refit runs saSPT again on every resampling draw. It uses 500 FOV-bootstrap and 500 FOV-permutation
  draws, so the smallest permutation p is 1/501 (≈ 0.002).
- These panels therefore carry at most **, and the stats CSVs flag p-values at that floor (`p_at_floor`).
- The 5F zone-label permutation needs no refit; it uses 2000 draws, and the treatment permutation uses
  10 000.

**No-refit companion outputs.** `build_claim2_panels.py` and `build_claim3_panels.py` also write no-refit
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
| `TLSD_OUT` | output root |
| `PANEL_VARIANT` | `36state` for every panel |
| `TLSD_PYTHON` | interpreter used by the Slurm wrappers (default `python3`) |
| `TLSD_SLURM_PARTITION`, `TLSD_SLURM_ACCOUNT` | optional, for the jobs `er/batch_er_workflow.py` submits |

**Slurm.** The wrappers carry no partition or account; pass them to `sbatch`
(e.g. `sbatch -p <partition> ...`). The panel builders take 5–20 min and up to ~40 GB RAM each.

**Writes into the data tree.** Only `merge_contrast_cache.py` writes there, adding the joint cache to the
tier's `occupancy_compare_cache_36state/`. Every other script only reads the data tree and writes under
`TLSD_OUT`.

## Run order
With the deposited caches and refit draws, only steps 3, 4 and 6 are needed. Step 3 rebuilds the joint caches,
which are exact concatenations of the deposited per-condition caches and are not deposited. Steps 1, 2 and 5
regenerate the other inputs from the raw movies. The deposit holds trajectories pooled per condition
(`traj_combinations/`), not per field of view.

1. **Tracking.**
   - Run `pipeline/call_filewise_operations.slurm` (`run_filewise.py`). For the ER data, run
     `er/batch_er_workflow.py`, which runs `run_filewise_er.py`.
   - Then run `call_conditionwise_operations.slurm` and `call_combinationwise_operations.slurm`.
     - ER zones: `er/generate_er_combinations.py` + `er/combo_array.sh`.
     - Re-zoning: `er/rezone_array.sh --pct-floor 60 --require-broad false --min-frames 2`.
   - Then run `call_megapool_operations.slurm` (2026_07_08).
2. **36-state likelihood caches.** Run `pipeline/submit_36state_caches.py`, which calls
   `call_variant_cache.slurm` → `build_variant_cache_one.py`.
3. **Joint caches for the permutation nulls.** For each contrast in `figures/submit_panel_refit.py`, run
   `figures/merge_contrast_cache.py --tier <tier> --cond-a A --cond-b B`.
4. **Panels.**
   - `sbatch pipeline/call_build_panels.slurm build_claim1_panels.py`, and the same for
     `build_claim2_panels.py` and `build_claim3_panels.py`.
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
