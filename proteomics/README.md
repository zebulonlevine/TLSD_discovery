# Proteomics (Fig 3B, 4A, 4E, 7C; Fig S5B–F, S6B–E; Tables S8, S9)

Two TMT experiments.
- **SILAC pulse-chase.** Quiescent RPS3-HaloTag cells, sgNT or sgDDRGK1, were fully heavy-labelled,
  chased into light medium, and collected at day 3 and day 8, with 3 replicates each. Heavy is old
  protein; light is new protein.
- **Total protein.** RPS3-HaloTag sgNT or sgDDRGK1 cells, growing or quiescent, plus quiescent
  RPL29-HaloTag sgDDRGK1, with 3 replicates each.

| File | Role |
|---|---|
| `silac_turnover.py` | SILAC: normalization, day-8/day-3 turnover, sgDDRGK1 − sgNT deltas, Fig 7C test. Writes Table S8 and Fig 3B, S5B, S5D–F and 7C |
| `total_abundance.py` | Total protein: median normalization and sgDDRGK1 vs sgNT tests. Writes Table S9 and Fig 4A/S6B, 4E/S6C, S6D/E |
| `figS5C_correlation.py` | Correlation between the SILAC and total-protein comparisons (Fig S5C) |
| `silac_proteomics.py` | Helper library: loading, normalization, error propagation, tests, categories, plots |
| `settings.py` | Input paths, plex layout, comparisons, category and gene-set definitions |
| `annotation/export_accessions.py` | Writes the UniProt accession list for each table |
| `annotation/uniprot_topology_pipeline.py` | Membrane topology from UniProt features, plus MitoCarta 3.0 flags (Fig 4E, S6C) |

## Statistics
- **SILAC normalization.** Each sample's heavy and light signals are divided by that sample's total
  (H + L) signal.
- **Turnover.** day-8 / day-3 ratio of the replicate means. SEM is propagated, then the ratio is
  log2-transformed.
- **Delta.** sgDDRGK1 − sgNT on the log2 turnover. Its SEM is the root-sum-square of the two SEMs;
  t = delta / SEM; two-sided p from t with 4 df.
- **Total protein.** The log2 abundances of each sample are median-centred. log2FC is the difference
  of the replicate means, tested the same way (t with 4 df). p is BH-adjusted across proteins.
- **Fig 7C.** Per protein, the light (new) signal at day 8 is log2-transformed for each replicate.
  sgDDRGK1 and sgNT are compared with a Student t-test (3 vs 3), and p is BH-adjusted across proteins
  quantified in all six channels.
- **Category and topology tests.** Each group is compared with the rest of the proteome by
  Mann-Whitney U. Groups are subsampled to ≤ 50 proteins; the reported p is the median over 1,000
  seeded draws, BH-adjusted across groups.
- **Localization categories.** Taken from the organelle-IP table, with these adjustments:
  - nucleolus counts as nucleus;
  - stress granule, P-body and "translation" count as cytosol;
  - ribosomal proteins are assigned from the RPL/RPS gene lists.

  This is the `category_locAdj` column, used in Fig 3B, 4A, S5D/E and S6B.

## Outputs (`<out>/proteomics/`)
| Panel / table | File(s) |
|---|---|
| Fig 3B | `fig3B_delta_old_by_category.pdf`, `fig3B_delta_old_mwu.csv` |
| Fig S5B | `figS5B_intensity_violins.pdf` |
| Fig S5C | `figS5C_correlation_heatmap.pdf`, `figS5C_correlation.csv` |
| Fig S5D, S5E | `figS5D_old_turnover_by_category.pdf`, `figS5E_new_turnover_by_category.pdf`, `figS5DE_turnover_{old,new}_{NT,sgDDRGK1}_mwu.csv` |
| Fig S5F | `figS5F_er_gene_sets_strip.pdf`, `figS5F_values.csv` |
| Fig 7C | `fig7C_new_protein_volcano.pdf`, `fig7C_new_protein_stats.csv` (days 3 and 8) |
| Fig 4A, S6B | `fig4A_figS6B_total_by_category.pdf` (Fig 4A shows the RPL / RPS groups), `fig4A_figS6B_{growing,quiescent}_mwu.csv` |
| Fig 4E, S6C | `fig4E_figS6C_total_by_topology.pdf`, `fig4E_figS6C_{growing,quiescent}_mwu.csv` |
| Fig S6D, S6E | `figS6D_volcano_quiescent.pdf`, `figS6E_volcano_growing.pdf` |
| Table S8 | `TableS8_SILAC_annotated.csv`, `.xlsx` |
| Table S9 | `TableS9_nonSILAC_proteomics.csv`, `.xlsx` |

Fig 7C labels TLSD1 by its gene symbol, R3HDM4.

## Inputs
| File | Content | Location |
|---|---|---|
| `protein_quant_30972.tsv` | SILAC-TMT protein quantification (heavy and light signal-to-noise sums for each of 12 channels) | ships with the code: `input/` |
| `all_ttest_results.csv` | Total-protein TMT table (log2 abundances for 15 samples) | ships with the code: `input/` |
| `silac_topology.tsv`, `total_topology.tsv` | Membrane topology and MitoCarta flags | ships with the code: `input/`; written by `annotation/uniprot_topology_pipeline.py` (UniProt release 2026_01; see Run) |
| `mmc4_published_OpenCell.xlsx` | Organelle-IP localization table (sheet `subcell_annotation_summary`) | download (ref. 170, supplementary file mmc4) to `<data>/proteomics/` |

The two quantification tables are the processed outputs of the searches deposited at
ProteomeXchange/PRIDE (accession TBD), which holds the raw files. `all_ttest_results.csv` (9.1 MB) is
the one file in the repository over 5 MB. The per-protein results are Supplementary Tables S8 and S9.

## Run
```
python proteomics/silac_turnover.py        # ~1 min each; proteomics env
python proteomics/total_abundance.py
python proteomics/figS5C_correlation.py
```
Expected console output: `6075 proteins; R3HDM4 day-8 new protein: log2FC +2.665, BH p 0.00079`,
then `8023 proteins`, then `3600 proteins in both experiments`.

Topology annotation (`membrane_analysis` env; needs network access to rest.uniprot.org). MitoCarta 3.0
is the `Human.MitoCarta3.0.xls` download from the Broad Institute.
```
python proteomics/annotation/export_accessions.py --outdir topo/
python proteomics/annotation/uniprot_topology_pipeline.py --input topo/total_uniprotIDS.csv \
    --mitocarta Human.MitoCarta3.0.xls --output total_topology.tsv --workdir topo/work
# and the same for topo/silac_uniprotIDs.csv -> silac_topology.tsv
```
A fresh run queries the current UniProt release and records it in `<output>.meta.txt`. The tables in
`input/` come from release 2026_01 (run 2026-05-11); to use a fresh annotation instead, write it to
`input/` under the same names.

## Verification
All runs used the original analysis notebook (cells 0–123) as the reference, on 2026-10-23, with the
notebook's ribosomal gene lists. The list fix below was applied afterwards (2026-10-24), and only the
RPL/RPS-dependent outputs changed.
- **Tables and statistics.** Every statistic and every shared table column is identical
  (max |Δ| = 0).
- **Figures.** Every figure PDF renders pixel-identical at 100 dpi. The one exception is the
  Fig 4E/S6C title, which was corrected (see below).
- **Topology.** Re-running the topology pipeline from the cached UniProt entries reproduces both
  topology tables byte-for-byte.

## Changes from the original notebook
- **Fig 4A/S6B** are grouped by `category_locAdj`, like Fig 3B and S5D/E. The notebook grouped
  them by the unadjusted categories. The significance stars are unchanged.
- **Fig 4E/S6C title** now reads "by membrane topology" (it said "by localization").
- **Ribosomal gene lists** (`settings.py`), 2026-10-24:
  - Current symbols: RPLP0 (was `RPP0`) and RPS4X (was `RPS4`). Under the old names, neither
    protein was assigned to RPL/RPS.
  - RPL7L1 removed from RPL. It is a 60S assembly factor, not a ribosomal protein.
  - Affected: Fig 3B, 4A, S5D–F, S6B and the category columns of Tables S8/S9.
  - Text changes in the placed panels: only the RPL7L1 label in S5F is gone.
  - The subsampled MWU shares one random stream across categories, so changing a group size
    shifts every later category's draws. Calls near a threshold can move for that reason alone.
- **Table S8:**
  - The Fig 7C columns (`light_8d_sgDDRGK1_vs_NT_log2FC*`) carry the test plotted in Fig 7C;
    `_sem` is the pooled standard error.
  - `category_localization` is `category_locAdj`.
  - `ER_GeneSets_delta` gives each protein's Fig S5F gene set.
- **Table S9:** `Annotation` keeps the protein description. The notebook had overwritten it with the
  UniProt accession.
- **Removed (not in the paper):**
  - RPL29 SILAC replicate experiment;
  - amino-acid composition;
  - Adamson 2016 UPR-branch, ER-phagy, ERAD and RBP subsets;
  - direct sgDDRGK1 vs sgNT difference tests per day and label;
  - TLSD1-vs-group panels;
  - replicate correlation maps;
  - single-condition violins and box plots;
  - exploratory volcanoes;
  - clipboard exports.
