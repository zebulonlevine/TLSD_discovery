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

## Ratio sorts: sorter names vs new/old names
The sorter gated the ratio sorts on the **old/new** ratio (JF549/JF646), in both screens. The sample
names here keep the sorter's gate names, as do the FASTQ files: genome-wide `top`/`bottom` and
sublibrary `ratioON_Top`/`ratioON_Bottom` are the top and bottom 25% of **old/new**.

The paper, the supplementary tables and GEO report these samples by **new/old** ratio. The sorter's
**bottom** gate is the top 25% of new/old, and its **top** gate is the bottom 25%.

The ratio comparisons (`*_TopvBottom*` genome-wide, `*_ratioON_top_v_bottom_*` sublibrary) take the
sorter's bottom gate as treatment. Their fold change is therefore log2(new/old top / new/old bottom),
the same direction as the new and old sorts. The synthesis/degradation scores rely on this (see
`../scoring/`). Every value in the tables follows this orientation; only the column names say new/old:
- Table S1: `RatioNewOld_Top` = sorter `bottom`.
- Tables S4/S10: `RatioNewOld_Top` = sorter `ratioON_Bottom`.
- S1, S4/S10 and S5/S11 all name the fold change `RatioNewOld_TopvBottom`.

`geo_sample_names.csv` maps each GEO sample name (series GSE352053) to the sample name used here.

## Files
| File | Purpose |
|---|---|
| `screen_pipeline.py` | Stage driver and gene-level U-test |
| `count_spacers.py` | Exact-match spacer counting after the key sequence (+1G aware) |
| `restriction_sites.py` | Sublibrary guide annotation and the 53 NotI/BamHI-site guides excluded from all analyses |
| `gene_aliases.py` | Display labels for the supplementary tables (TLSD1/R3HDM4, UBA5/ACAD11, ATG9A/ANKZF1, RMP64) |
| `geo_sample_names.csv` | GSM accession and GEO sample name → sample name used here, with each sample's fraction. The ratio sorts are new/old at GEO and old/new (sorter) here |
| `mageck_counts.py` | Puts sub-threshold guides' counts on MAGeCK's normalised scale for the tables |
| `config_*.yaml`, `sample_table_*.csv`, `library_table_*.csv`, `comparisons_*.csv` | Screen definitions |

## Inputs not in this repository
- **FASTQs (GEO, GSE352053).** Place them under `${TLSD_DATA}/fastq/{genomewide,sublibrary}/`, at the paths given in
  `sample_table_*.csv`.
- **Count tables (GEO, GSE352053).** To start at the MAGeCK stage (stage 4), build `counts_for_mageck.txt` from the
  GEO count tables. Rows are in the same order and the counts are the ones the published results used.
  GEO names the ratio-sorted samples by new/old ratio (e.g. `L29rep1_ratioNewOld_top`). Rename the
  sample columns to the names used here with `geo_sample_names.csv` (`geo_name` → `sample_name`; e.g.
  `L29rep1_ratioNewOld_top` → `L29rep1_bottom`). The other samples keep their names.
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
- **Library design files** ship with the code in `../library/`:
  - genome-wide: `221010_Weissman_Lib1A_targets.csv`, `221111_lib1_targets.csv` and
    `230215_extra_isoforms_oligos_single.csv`, read through `library_table_genomewide.csv` for counting
    from FASTQs;
  - sublibrary: `sublibrary.csv` and `sublibrary_spacers.csv`, read through `library_table_sublibrary.csv`
    for counting and by `restriction_sites.py` (the NotI/BamHI filter used by scoring and Tables S4, S5,
    S10, S11).

  Guide names, genes and spacers are in the GEO count tables, which also flag the 53 NotI/BamHI-site
  sublibrary guides, and in Tables S1, S4 and S10.

## Notes
- For comparisons with one replicate per arm (the growing-condition genome-wide screens), MAGeCK fits
  its mean–variance model to the pooled treatment and control samples. That makes those p-values
  conservative.
- `slurm.partition` in the configs is the Whitehead cluster's queue name; change it for other clusters.

Credits: sgRNA counting based on code by Julia Joung; MAGeCK-iNC method from the Kampmann lab.
