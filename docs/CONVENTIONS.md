# Conventions for code in this repository

- **Paths**: never hardcode absolute paths. Use `common/paths.py`:
  - `paths.inp(...)` for reads under the data root;
  - `paths.outp(...)` for writes under the output root;
  - `paths.repo(...)` for files shipped in this repo (categories, configs, library design tables).
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
