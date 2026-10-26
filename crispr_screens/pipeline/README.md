# CRISPRi screen pipeline: FASTQ → counts → MAGeCK → gene scores

Used for the genome-wide screens (`config_genomewide.yaml`) and the sublibrary screens
(`config_sublibrary.yaml`). Environment: `envs/crispri-screen.yml`. Jobs run under Slurm.

## Stages
1. **Demultiplex** with bcl-convert (outside this repo).
2. **Count**: `python screen_pipeline.py count --config <config>`. This writes one Slurm job per
   FASTQ × library file, each calling `count_spacers.py` (mU6 preset, key `CACCTTGTTG`, `-no-g`). Each
   sample is counted only against its own library.
3. **Aggregate**: `python screen_pipeline.py aggregate --config <config>` →
   `results/{aggregated_counts.csv, counts_for_mageck.txt}`.
   > **Do not re-run this stage for the sublibrary.** sgRNA names are assigned positionally, and the
   > NTC-pseudogene renaming block in `screen_pipeline.py` post-dates the sublibrary's 2026-01-15 run.
   > Re-running would renumber every non-targeting guide and invalidate all downstream sublibrary
   > results. Use the deposited `counts_for_mageck.txt` instead.
4. **MAGeCK + gene scores**: `python screen_pipeline.py mageck --config <config>`. This writes one job
   per comparison (`comparisons_*.csv`), each running three steps:
   - `threshold`: keeps a guide only if every replicate pair has ≥50 summed reads;
   - `mageck test -k -t -c -n` (MAGeCK 0.5.9.5, median normalisation);
   - `utest`: MAGeCK-iNC-style gene scores. These are the mean LFC of the 3 guides with the smallest
     MAGeCK p, with a two-sided Mann–Whitney U test against the non-targeting p-values, against NTC
     pseudogenes (seed 42). Settings: `min_treat_mean` 5, FDR 0.05.
5. **Summary**: `python screen_pipeline.py results --config <config>`.

`python screen_pipeline.py run --config <config> --wait` chains stages 2–5.

## Orientation warning
The `ratioON` comparisons use the **bottom** sorted fraction as treatment and the **top** as control,
so their fold change is log2(bottom/top). This is deliberate: it makes the slope +1/−1 geometry of the
synthesis/degradation scores correct (see `../scoring/`). Check any "Top"/"Bottom" table column for
these comparisons against `comparisons_*.csv`, not against the name.

## Files
| File | Purpose |
|---|---|
| `screen_pipeline.py` | Stage driver and gene-level U-test |
| `count_spacers.py` | Exact-match spacer counting after the key sequence (+1G aware) |
| `restriction_sites.py` | Sublibrary guide annotation and the 53 NotI/BamHI-site guides excluded from all analyses |
| `gene_aliases.py` | Display labels for the supplementary tables (TLSD1/R3HDM4, UBA5/ACAD11, ATG9A/ANKZF1, RMP64) |
| `mageck_counts.py` | Puts sub-threshold guides' counts on MAGeCK's normalised scale for the tables |
| `config_*.yaml`, `sample_table_*.csv`, `library_table_*.csv`, `comparisons_*.csv` | Screen definitions |

## Inputs not in this repository
- **FASTQs (GEO).** Place them under `${TLSD_DATA}/fastq/{genomewide,sublibrary}/`, at the paths given in
  `sample_table_*.csv`.
- **Count tables (GEO).** To start at the MAGeCK stage (stage 4), build `counts_for_mageck.txt` from the
  GEO count tables. Rows are in the same order and the counts are the ones the published results used.
  - Genome-wide, `${TLSD_DATA}/GenomeWide_results/results/`:
    - `counts_for_mageck.txt` from `genomewide_aggregated_counts.csv`: tab-separated; columns `sgRNA`,
      `Gene` (the sgRNA name without its `_sg<n>` suffix), then the 22 sample columns.
    - For Table S1, save the GEO table itself as `aggregated_counts.csv`; it supplies `Library` and
      `sgRNA Sequence`, row-aligned.
    - Also write `sgrna_annotation.csv` with the columns `sgRNA`, `Gene`, `sgRNA Sequence`.
  - Sublibrary, `${TLSD_DATA}/Sublibrary_results/results/`:
    - `counts_for_mageck.txt` from `sublibrary_aggregated_counts.csv`: rename `sgrna` → `sgRNA` and
      `gene` → `Gene`, then add the sample columns.
    - GEO omits the eight `outgrown` samples; no figure or table uses them.
  - The published per-comparison results and gene scores are Supplementary Tables S1–S6, S10 and S11.
- **Library design files.** These are not deposited:
  - `${TLSD_DATA}/GenomeWide_results/libraryInfo/`: `221010_Weissman_Lib1A_targets.csv`,
    `221111_lib1_targets.csv`, `230215_extra_isoforms_oligos_single.csv` (counting from FASTQs only);
  - `${TLSD_DATA}/Sublibrary_results/libraryInfo/`: `sublibrary.csv`, `sublibrary_spacers.csv` (read by
    `restriction_sites.py`).

  Guide names, genes and spacers are in the GEO count tables, which also flag the 53 NotI/BamHI-site
  sublibrary guides, and in Tables S1, S4 and S10.

## Notes
- For comparisons with one replicate per arm (the growing-condition genome-wide screens), MAGeCK fits
  its mean–variance model to the pooled treatment and control samples. That makes those p-values
  conservative.
- `slurm.partition` in the configs is the Whitehead cluster's queue name; change it for other clusters.

Credits: sgRNA counting based on code by Julia Joung; MAGeCK-iNC method from the Kampmann lab.
