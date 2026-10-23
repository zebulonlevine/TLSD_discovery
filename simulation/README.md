# Turnover simulation (Fig S2A)

Cohort simulation of zero-order synthesis and first-order degradation. Starting from steady state
(t½ = 7 d, synthesis 100 molecules/day, time step 0.05 d), the run either raises synthesis 3× or raises
degradation by the factor (≈1.46×) that gives the same new/old ratio 14 days after the change.
"New" means no older than 7 days, the initial half-life.

| File | Role |
|---|---|
| `biomolecule_turnover.py` | Simulation (cohort bookkeeping; `python biomolecule_turnover.py` runs a steady-state self-check) |
| `figS2A_turnover_simulation.py` | Runs both scenarios and writes the outputs below to `<out>/simulation/` |

Outputs:
- `figS2A_ratio_and_abundance.pdf`: the two Fig S2A panels.
- `figS2A_curves.csv`: source data.
- `increased_ratio_summary.pdf`: a diagnostic view of total, old and new molecules per scenario.

```
python simulation/figS2A_turnover_simulation.py      # ~10 s; crispri-screen env (numpy, scipy, pandas, matplotlib)
```
Expected console output: degradation ×1.459, with new/old = 1.515 at day 14 in both scenarios.
