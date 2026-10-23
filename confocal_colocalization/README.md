# Confocal colocalization (Fig 7D)

Whole-cell Pearson correlation between GFP-R3HDM4 (TLSD1), Sec61β (ER, JF552) and L29 (60S, JFX650)
in live L29 cells, DMSO vs PT-179 (23 h files). Methods, acquisition settings and the channel map are
in `METHODS.md`.

| Step | Script | Output (under `${TLSD_OUT}/confocal_colocalization/`) |
|---|---|---|
| 1 | `segment_cells.py`: Cellpose-SAM (`cpsam`), 2D per plane + IoU stitching | `<seg-out>/cell_masks/`, `channels/`, `cell_qc/` |
| 2 | `wholecell_pearson.py`: per-cell Pearson r for each channel pair | `per_cell/*_wholecell_pearson.csv` |
| 3 | `plot_pearson_swarm.py`: swarm plot, medians, two-sided Mann-Whitney U | `figures/pearson_L29_swarm_DMSO_vs_PT179_23h.{pdf,_stats.csv}` |
| inset | `build_fig7D_outline_stacks.py`: Fiji hyperstacks with the Cellpose outlines as a 4th channel | `fig7D_outline_stacks/*.tif` |

`build_fig7D_outline_stacks.py` was written after draft5 and is likely the source of the 7D image inset.
It does not re-segment; it reads the step-1 masks.

## Environment
- Conda env: `envs/nd2_coloc.yml`. Step 1 needs a GPU.
- `TLSD_CONFOCAL_DATA`: experiment folder holding `Denoised/` (NIS Batch Denoise ND2s) and `output/cell_masks/`
  (the original masks; read only). It can also be set as `confocal_data` in `config/paths.yaml` for the
  Python scripts.
- `TLSD_OUT`: output root.
- `TLSD_PYTHON`: interpreter of the env (default `python3`).
- `TLSD_GPU_PARTITION` (step 1), `TLSD_SLURM_PARTITION` (steps 2–3; optional), `TLSD_SLURM_ACCOUNT` (optional),
  `TLSD_JOB_PREFIX` (default `${USER}_tlsd_confocal`).

## Run
```bash
bash submit_pearson_repro.sh                      # steps 2-3 on the original masks
bash submit_pearson_repro.sh --start-step 1 \
     --seg-out /path/to/new_segmentation          # full rerun incl. Cellpose-SAM (GPU)
python build_fig7D_outline_stacks.py              # inset stacks; ~4 GB written, run under Slurm
```

## Expected numbers (`*_stats.csv`)
| Pair | n (DMSO / PT-179) | median r (DMSO / PT-179) | U | p (two-sided) |
|---|---|---|---|---|
| GFP-R3HDM4 vs Sec61B | 276 / 192 | 0.376 / 0.799 | 10750 | 7.3 × 10⁻²⁸ |
| GFP-R3HDM4 vs L29 | 276 / 192 | 0.753 / 0.755 | 26308 | 0.90 |
| Sec61B vs L29 | 276 / 192 | 0.440 / 0.680 | 14063 | 5.7 × 10⁻¹⁸ |

Not shipped: the superseded 4-step pipeline (`quantify_coloc.py`, `make_er_mask.py`, `merge_and_plot.py`,
`submit_coloc_pipeline.sh`) and the exploratory notebook. Whole-cell Pearson does not use the ER mask.
