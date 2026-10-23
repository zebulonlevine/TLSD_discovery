# Supplementary tables S1–S6, S10, S11

Environment `envs/crispri-screen.yml`. Run after `../pipeline`, `../scoring` and `../go`.

| Table | Command | Output (under the output root) |
|---|---|---|
| S1, S2, S3 | `python build_supplementary_tables.py --config ../pipeline/config_genomewide.yaml` | `GenomeWide_results/supplementary_tables/{Table_S1_sgRNAs, TableS2_GeneLevel_genomewide, TableS3_GO_Genomewide_screens}.csv` |
| S4, S5, S10, S11 | `python build_sublibrary_tables.py` | `fixed_NotI/supplementary_tables/TableS{4,5,10,11}_*.csv` |
| S6 | `python ../go/build_tableS6.py` | `fixed_NotI/supplementary_tables/TableS6_GO_pooled_sublibraryScores.{csv,xlsx}` |

**Notes**
- `--tables s1,s2,s3` selects a subset.
- Gene labels: the combined display labels (TLSD1/R3HDM4, UBA5/ACAD11, ATG9A/ANKZF1, RMP64) are used
  in every table. `--gene-labels plain` affects only S1's `Gene` column. `gene_id` is never aliased.
- S4/S10 list the 4,875 cloned guides (the 53 NotI/BamHI-site guides are dropped). The spacer column
  holds the designed 20-nt spacer.
- Count columns hold MAGeCK median-normalised counts. Sub-threshold guides are filled on the same scale
  (`../pipeline/mageck_counts.py`).
- `build_sublibrary_tables.py --scores-root <SUB>/scores_pRanked --results-root <SUB>/results` rebuilds
  the superseded pre-re-scoring tables.
- Manuscript file names: S1 is shipped as `TableS1_sgRNAs_genomewide.csv`; the other tables keep the
  names above. Tables S7–S9 and S12 come from `../../proteomics/`, `../../single_molecule/`, or are
  reagent tables.
