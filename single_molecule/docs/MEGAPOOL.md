# Megapool: Joint-saSPT Across Conditions

## Purpose

When comparing diffusion distributions across conditions (or across
combinations), running saSPT separately on each yields posteriors that
are each fit with their own variational dynamics. Small differences
between per-condition fits can introduce noise that confounds
inter-condition comparisons.

A **megapool** runs ONE joint saSPT fit across all conditions you want
to compare, then partitions the per-subtrack MLE values back to source
conditions. All members are judged against the same joint occupations
posterior and the same diffusion-state cutoffs. This is the
statistically clean form of cross-condition comparison.

## How saSPT preserves the source mapping

Empirically verified against the saSPT runtime (Heckert 2022):

1. saSPT's `TrajectoryGroup.preprocess` removes singlets, splits long
   trajectories at `splitsize` jumps, and renumbers trajectories
   contiguously 0..N-1 internally.
2. saSPT preserves the input trajectory ID as `orig_trajectory` in
   `SA.trajectories.detections`.
3. `executeSASPT` writes that to the MLE CSV as
   `orig_trajectory = "trajectories::pooled_file.csv::<INPUT_TRAJ_ID>"`.
4. **When `saspt.sample_size` >= pool size, this is a 1:1 mapping**:
   every input trajectory has its MLE rows tagged with its input ID.
5. When saSPT *subsamples* (`sample_size < pool size`), it renumbers
   ~7% of trajectories, breaking the 1:1 mapping for those.

The megapool pipeline forces `sample_size = 1e9` so saSPT never
internally subsamples. Pool-size control happens via
`max_trajectories_per_member` and `max_trajectories_total` in
`megapools.yaml` — these subsample at *pool construction time*, before
saSPT sees the data.

## Architecture

```
megapools.yaml                              (user-provided, experiment root)
    │
    ▼  analyze_directory.py
    │   submits one megapool SLURM job per entry,
    │   with afterok deps on each member's combo/filewise job
    │
    ▼  run_megapool.py
    │
    │  1. pool_megapool_trajectories
    │       for each member: read traj CSV(s), offset trajectory IDs
    │       write traj_megapools/{name}_pooled_traj.csv
    │       write traj_megapools/{name}_source_map.csv
    │                                  (trajectory_id → source_member)
    │
    │  2. executeSASPT on the pooled CSV
    │       writes posterior_megapools/{name}_pooled_posterior.csv
    │       writes MLE_megapools/{name}_pooled_MLE.csv
    │       (rolling-window analysis is SKIPPED)
    │
    │  3. partition_mle_by_source
    │       extract input traj ID from orig_trajectory via regex
    │       join with source map → source_member per subtrack
    │       write {name}_{member}_pooled_MLE.csv per member
    │       duplicate joint posterior as {name}_{member}_pooled_posterior.csv
    │
    │  4. process_condition per member
    │       reads per-member MLE + (joint) posterior
    │       determines/loads cutoffs (shared across all members)
    │       plots histograms, bar graphs
    │       saves plotting PKL
    │       does NOT do bleaching curves or movies (rolling skipped)
    │
    ▼  plots/megapool_{name}/plotting_pkls/{name}_{member}.pkl
       one PKL per member, all with the same posterior_plot_df
       (the shared joint posterior) and per-member MLE_df partitions
```

## Files

| file | role |
|---|---|
| `megapool_functions.py` | YAML loader, pool builder, partition logic, settings redirector |
| `run_megapool.py` | SLURM-callable entry point |
| `call_megapool_operations.slurm` | SLURM wrapper |
| `megapools.example.yaml` | annotated example to copy into experiment dir |

## `megapools.yaml` schema

Place at experiment root. Picked up automatically by `analyze_directory.py`.

```yaml
megapools:
  knockdown_comparison:
    members:                              # required, list of:
      - sgNT                              #   - combination nicknames (from combinations.yaml), OR
      - sgUBA5                            #   - raw condition names
      - sgR3HDM4
      - sgUBA5_sgR3HDM4
      - oeR3HDM4
      - oeR3HDM4_sgUBA5
    max_trajectories_per_member: 100000   # optional; balances members of different size
    max_trajectories_total: null          # optional; cap after per-member balancing
    random_seed: 42                       # optional; for reproducible subsampling
    n_components: 4                       # optional GMM cutoff override
    labels: [immobile, slow, fast, very_fast]  # optional
    initial_centers: [-1.3, -0.5, 0.3, 1.0]    # optional
```

**Naming rule:** megapool names must not contain `/ \ : * ? " < > |` or
whitespace. Validation rejects these.

**Member rule:** members can be combo nicknames or raw conditions. Mix
freely. Each member's pre-pooled trajectory data is what gets pooled
again into the megapool.

## Running

### Standard flow (via `analyze_directory.py`)

When `megapools.yaml` exists at the experiment root, the orchestrator
submits one megapool SLURM job per entry, with `afterok` dependencies
on each member's contributing pipeline step:
- Combo members wait for their `run_combinationwise.py` job
- Raw-condition members wait for their filewise array

The megapool job runs after every member's tracking + saSPT output
exists.

```bash
python analyze_directory.py --input_directory /path/to/experiment
```

### Re-running just one megapool

```bash
EXP=/path/to/experiment
CODE=/path/to/code/dir
PKL=$EXP/tracking_output_q=*/settings.pkl

sbatch --export=MEGAPOOL=knockdown_comparison,PKL_PATH=$PKL,CODE_DIRECTORY=$CODE,FORCE=1 \
       $CODE/call_megapool_operations.slurm
```

`FORCE=0` (default): if the joint saSPT outputs already exist, skip
pooling and saSPT and re-run only partitioning + downstream processing.
Useful for trying different cutoffs.

`FORCE=1`: re-pool, re-run saSPT, re-run downstream.

## Output layout

Under `tracking_output_q={t}/`:

```
traj_megapools/
  {name}_pooled_traj.csv              # the pooled detection CSV fed to saSPT
  {name}_source_map.csv               # trajectory_id → source_member

posterior_megapools/
  {name}_pooled_posterior.csv         # the JOINT posterior (shared across members)
  {name}_{member}_pooled_posterior.csv  # duplicate per member (for process_condition glob)

MLE_megapools/
  {name}_pooled_MLE.csv               # the JOINT MLE (all subtracks)
  {name}_{member}_pooled_MLE.csv      # member-partitioned MLE rows

plots/megapool_{name}/
  plotting_pkls/
    {name}_{member}.pkl               # one per member, ready for spt_compare
    global_cutoffs_megapool_{name}.pkl  # shared cutoff cache (if auto_gmm)
  <plot files>                        # histograms, bar graphs per member
```

## Visualizing with `spt_compare.py`

Point `SPTExperiment` at the megapool subdirectory:

```python
from spt_compare import SPTExperiment

exp = SPTExperiment(
    '/path/to/.../plots/megapool_knockdown_comparison/plotting_pkls')

# All members are now in exp.conditions, labeled
# 'knockdown_comparison_sgNT', 'knockdown_comparison_sgUBA5', etc.
exp.conditions

# Optionally rename for cleaner plots
display_names = {
    'knockdown_comparison_sgNT':            'control',
    'knockdown_comparison_sgUBA5':          'sgUBA5',
    'knockdown_comparison_sgR3HDM4':        'sgR3HDM4',
    'knockdown_comparison_sgUBA5_sgR3HDM4': 'sgUBA5+sgR3HDM4',
    'knockdown_comparison_oeR3HDM4':        'oeR3HDM4',
    'knockdown_comparison_oeR3HDM4_sgUBA5': 'oeR3HDM4+sgUBA5',
}
exp = SPTExperiment(
    '/path/to/.../plots/megapool_knockdown_comparison/plotting_pkls',
    display_names=display_names)

# Overlay MLE KDEs (per-member subtrack distributions, joint-fit)
exp.overlay_mle_kde(exp.conditions, xlim=(-3, 2), ylim=(0, 1.5))

# Bar graph of state fractions (all using the shared cutoffs)
exp.bar_graph(exp.conditions)
```

## Two important caveats

### 1. Posterior is shared, not per-member

Every member's `posterior_plot_df` is the **same joint posterior** P(D)
from the joint fit. It's not a per-member posterior. The per-member
content is in:
- `MLE_df` — partitioned subtracks, one set per member
- Subsequent KDE/bar-graph operations on `MLE_df`

If you plot "P(D) overlay" with `overlay_posterior`, all members will
look identical (they share the same posterior). Use `overlay_mle_kde`
or `overlay_histogram` to see actual per-member differences in
diffusion-state assignments.

### 2. Memory and time scaling

From the benchmark (sgUBA5 pool, 334k trajectories):

| sample_size | peak RSS | wall time |
|---|---|---|
| 100,000 | 2.4 GB | 5 min |
| 200,000 | 4.1 GB | 7.5 min |
| 334,000 | 5.1 GB | 9 min |

For a typical 6-combo megapool (~1.5M trajectories pool size after
balancing), expect ~25 GB peak RSS and ~40 min wall time. The SLURM
wrapper requests 64 GB to leave comfortable headroom. Adjust the
`#SBATCH --mem=` line if your data is much larger.

## What `analyze_directory.py` flags do

- `--input_directory PATH` (required) — experiment root
- `--force` — re-run filewise tracking for ALL files, AND re-run combos
  AND megapools. Propagated as `FORCE=1` to SLURM jobs.

There is currently no flag for "re-megapool only without re-tracking";
use the manual `sbatch` invocation above for that.

## When NOT to use a megapool

- **Comparing very dissimilar experiments** (different probe, different
  cell line, different imaging condition). The joint fit will average
  characteristics across them. Use per-condition fits with shared
  cutoff cache instead (the existing `groups.yaml` mechanism).
- **You only have one condition to plot.** Just use the existing
  per-condition or per-combo run.
- **You need per-condition rolling-window data** (bleaching curves,
  movies). Megapools skip these. Use the per-condition / per-combo
  plots for those.

## Troubleshooting

**"Megapool 'X' not found in megapools.yaml"**
Check the spelling matches exactly. Names are case-sensitive.

**"Member 'X' is a combination but its pooled trajectory CSV is missing"**
The combination's `run_combinationwise.py` hasn't run yet, or the
`--force` flag was passed but combinationwise outputs were deleted.
Run combinationwise first, then megapool.

**"Member 'X' is neither in combinations.yaml nor a recognized condition"**
The member name doesn't match either a combo nickname or a raw
condition's file prefix. Check both `combinations.yaml` and your
`.nd2`/`.tif` file naming convention.

**Lots of `n_unmatched` in the partition stats**
The 1:1 mapping was broken — likely because `saspt.sample_size`
wasn't set high enough. The pipeline forces it to 1e9, so this
shouldn't happen unless something in `executeSASPT` changed.

**Memory failures**
Bump `--mem` in `call_megapool_operations.slurm`, or lower
`max_trajectories_per_member` in `megapools.yaml`.
