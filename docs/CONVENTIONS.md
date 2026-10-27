# Conventions for code in this repository

- **Paths**: never hardcode absolute paths. Use `common/paths.py`:
  - `paths.inp(...)` for reads under the data root;
  - `paths.outp(...)` for writes under the output root;
  - `paths.repo(...)` for files shipped in this repo (categories, configs, library design tables, the small
    inputs in `crispr_screens/library/` and `proteomics/input/`).
  - YAML configs may use `${TLSD_DATA}`, `${TLSD_OUT}` and `${TLSD_REPO}`.
- **Layout**: both roots mirror the original analysis tree (see `config/paths.example.yaml`). Scripts keep
  the original relative sub-paths so outputs can be compared file-for-file with the published ones.
- **Verification**: each stage is checked by reading the published tree (`TLSD_DATA`) and writing to scratch
  (`TLSD_OUT`). A copied script must reproduce the published output exactly. Figures are compared on label
  sets, axis ticks and data, because adjustText re-places labels on every run. Any numeric difference is
  reported, never "fixed".
- **Heavy steps** (MAGeCK, scoring, figure rebuilds) run under Slurm (`srun --ntasks=1`); see each
  module's README.
- **Colours** (fixed): non-targeting controls `#DDDDDD`, other genes `#4D4D4D`.
- **Panel scripts** are named by manuscript panel (draft5 numbering).
- **Code brought in from outside the analysis tree** (laptop notebooks/scripts) is confirmed, cleaned and
  made presentable:
  1. **Confirm**: re-run it and show the cleaned version reproduces the original outputs exactly
     (arrays / tables bit-identical; figures text- and pixel-identical), and that it matches the
     manuscript panel.
  2. **Clean**: remove dead code, demo blocks, unused imports, sandbox or absolute paths, and stale
     comments.
  3. **Present**: module docstring stating the panel and method, functions plus `main()`, outputs via
     `common/paths`, a per-folder README with the run command and expected numbers, and source-data
     CSVs where cheap.
