# GO-term enrichment (decoupler ULM)

Environment `envs/dc.yml` (decoupler 2.1.6), except `build_tableS6.py`, which needs pandas + openpyxl
(`envs/crispri-screen.yml`).

| Script | Output | Manuscript |
|---|---|---|
| `go_genomewide.py` | `GSEA_Fig1/{S3,L29}/…csv` | Fig 1G, Table S3 (via `../tables/`) |
| `go_sublibrary.py` | `fixed_NotI/GO_decoupler/*_GO_decoupler_fixedNotI.csv` | Fig S2E, Table S6 |
| `build_tableS6.py` | `fixed_NotI/supplementary_tables/TableS6_GO_pooled_sublibraryScores.{csv,xlsx}` | Table S6 |

**Method.**
- A univariate linear model of gene epsilon on GO-term membership (`dc.mt.ulm`, tmin = 3).
- p-values are recomputed from t in log space and Benjamini–Hochberg-corrected across all terms of one
  screen.
- The background is the genes in the screen: genome-wide for Fig 1G, the 715-gene sublibrary for
  Fig S2E and Table S6.

**Gene sets** (not committed). Enrichr 2025 GO Biological Process, GO Molecular Function and GO
Cellular Component, concatenated into `${TLSD_DATA}/GSEA_Fig1/S3/gene_sets.gmt` (6,922 terms; 1.9 MB).
Download from https://maayanlab.cloud/Enrichr/#libraries (`GO_Biological_Process_2025`,
`GO_Molecular_Function_2025`, `GO_Cellular_Component_2025`).

**Reproducibility.**
- `go_sublibrary.py` and `build_tableS6.py` reproduce the published files byte-for-byte.
- `go_genomewide.py` reproduces the published genome-wide tables to floating-point precision
  (|Δt| ≤ 2e-13, identical FDRs, gene lists and overlaps). It is deterministic across thread counts;
  the residual comes from decoupler's compiled ULM on different CPUs. Near-tied terms (|Δt| < 2e-14)
  can swap row order.
