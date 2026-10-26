# Sublibrary synthesis / degradation scores and per-sgRNA p-values

Run in this order after the pipeline's MAGeCK stage (`../pipeline/`). Environment
`envs/crispri-screen.yml`. Run every step under Slurm (`srun --ntasks=1`).

Inputs: the MAGeCK outputs under `${TLSD_DATA}/Sublibrary_results/results/` (built from the GEO count
table; see `../pipeline/README.md`) and the sublibrary library files. The published scores are
Supplementary Tables S4, S5, S10 and S11.

1. `bash run_line_scores.sh`
   - Per-sgRNA line scores for the 8 sublibrary screens (60S/40S × synthesis/degradation × sgNT/sgUBA5),
     written to `Sublibrary_results/scores/<dir>/<dir>_sgrna_scores.csv`.
   - x is the ratioON comparison. y is the new (synthesis, line y = x) or old (degradation,
     line y = −x) comparison.
   - score = parallel·|parallel| / (|parallel| + perp).
2. `python export_scores.py`
   - Per-sgRNA p-values: depth-dependent noise model (`sgrna_score_pvalues.build`) with an empirical
     non-targeting null (`final_pvalues.sgrna_table`).
   - Gene-level roll-up: top 3 guides by p (ties broken by |score|), Mann–Whitney U against the NTC
     p-values, a guide-count-matched pseudogene null, and 5% FDR on epsilon × −log10 p
     (`final_pvalues.gene_table`).
   - Output: `fixed_NotI/scores/`. These are the published gene scores.
   - **The 53 NotI/BamHI-site guides** (`../pipeline/restriction_sites.py`) are excluded before the
     gene roll-up. `--keep-restriction-sites` reproduces the superseded pre-2026-10-21 scores.
3. `python rebuild_results.py`
   - Re-runs the per-comparison gene-level U-test with the same 53 guides removed (text-filtered
     `sgrna_summary.txt`; MAGeCK itself is not re-run).
   - Output: `fixed_NotI/results/`, used for the log2FC and growth blocks of Tables S5/S11 and for
     Figs S3C/D/F and S8A/E/F.

All draws are seeded (seed 42); re-running reproduces the published files byte-for-byte.
