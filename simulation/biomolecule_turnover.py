"""Cohort simulation of biomolecule turnover (Fig S2A).

Zero-order synthesis and first-order degradation,

    dN/dt = k_syn - k_deg * N,

integrated with explicit time steps. Every time step adds a cohort of
k_syn * dt molecules, and every cohort decays by exp(-k_deg * dt) per step.
Molecules are called "new" if they are no older than the INITIAL half-life
(ln 2 / k_deg_initial) and "old" otherwise, so at steady state with unchanged
rates new/old = 1. The cohort bookkeeping lets the age split stay exact when
the rates change part-way through a run.
"""
from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class SimulationParams:
    """Rates, initial abundance and time grid for one simulation."""
    k_syn_initial: float                  # synthesis rate (molecules / day)
    k_deg_initial: float                  # degradation rate constant (1 / day)
    N0: float                             # abundance at t = 0
    t_max: float                          # simulated duration (days)
    dt: float = 0.01                      # time step (days)
    t_change: Optional[float] = None      # time at which the rates change (None = never)
    k_syn_final: Optional[float] = None   # synthesis rate after t_change (None = unchanged)
    k_deg_final: Optional[float] = None   # degradation rate after t_change (None = unchanged)


@dataclass
class SimulationResults:
    time: np.ndarray
    total: np.ndarray
    old: np.ndarray                 # molecules older than the initial half-life
    new: np.ndarray                 # molecules no older than the initial half-life
    ratio_new_to_old: np.ndarray
    params: SimulationParams
    t_half_initial: float           # the age cutoff (initial half-life)


def get_rates(t: float, params: SimulationParams) -> tuple:
    """(k_syn, k_deg) in force at time t."""
    if params.t_change is not None and t >= params.t_change:
        k_syn = params.k_syn_final if params.k_syn_final is not None else params.k_syn_initial
        k_deg = params.k_deg_final if params.k_deg_final is not None else params.k_deg_initial
    else:
        k_syn = params.k_syn_initial
        k_deg = params.k_deg_initial
    return k_syn, k_deg


def calculate_steady_state(k_syn: float, k_deg: float) -> float:
    """Steady-state abundance, N_ss = k_syn / k_deg."""
    if k_deg == 0:
        return np.inf
    return k_syn / k_deg


def calculate_half_life(k_deg: float) -> float:
    """Half-life, t_1/2 = ln 2 / k_deg."""
    if k_deg == 0:
        return np.inf
    return np.log(2) / k_deg


def simulate_turnover(params: SimulationParams) -> SimulationResults:
    """Run the cohort simulation and split abundance into new and old molecules."""
    n_steps = int(params.t_max / params.dt) + 1
    time = np.linspace(0, params.t_max, n_steps)
    t_half_initial = calculate_half_life(params.k_deg_initial)

    # Start from the steady-state age distribution: cohorts going back ~5
    # half-lives, each the survivors of one time step's synthesis, rescaled so
    # the total equals N0.
    cohort_birth_times = []
    cohort_amounts = []
    n_initial_cohorts = int(5 * t_half_initial / params.dt)
    for i in range(n_initial_cohorts):
        age = i * params.dt
        surviving = params.k_syn_initial * params.dt * np.exp(-params.k_deg_initial * age)
        if surviving > 1e-10:
            cohort_birth_times.append(-age)
            cohort_amounts.append(surviving)
    total_initial = sum(cohort_amounts)
    if total_initial > 0:
        scale = params.N0 / total_initial
        cohort_amounts = [a * scale for a in cohort_amounts]
    cohort_birth_times = np.array(cohort_birth_times)
    cohort_amounts = np.array(cohort_amounts)

    total = np.zeros(n_steps)
    old = np.zeros(n_steps)
    new = np.zeros(n_steps)

    for i, t in enumerate(time):
        k_syn, k_deg = get_rates(t, params)

        if i > 0:                                   # synthesis during this step
            cohort_birth_times = np.append(cohort_birth_times, t)
            cohort_amounts = np.append(cohort_amounts, k_syn * params.dt)

        cohort_amounts = cohort_amounts * np.exp(-k_deg * params.dt)   # degradation

        ages = t - cohort_birth_times
        total[i] = np.sum(cohort_amounts)
        new[i] = np.sum(cohort_amounts[ages <= t_half_initial])
        old[i] = np.sum(cohort_amounts[ages > t_half_initial])

        # drop cohorts that have decayed below 1e-6 molecules
        keep = cohort_amounts > 1e-6
        cohort_birth_times = cohort_birth_times[keep]
        cohort_amounts = cohort_amounts[keep]

    old = np.maximum(old, 0)
    new = np.maximum(new, 0)
    total = old + new
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(old > 1e-10, new / old, np.inf)

    return SimulationResults(time=time, total=total, old=old, new=new,
                             ratio_new_to_old=ratio, params=params,
                             t_half_initial=t_half_initial)


if __name__ == '__main__':
    # Self-check: with unchanged rates and N0 at steady state, new/old sits at ~1
    # (1.07 at t=0 from discretising the initial age distribution, 1.004 by t=25).
    p = SimulationParams(k_syn_initial=100, k_deg_initial=0.1, N0=1000, t_max=50, dt=0.05)
    r = simulate_turnover(p)
    print(f'age cutoff (t1/2) = {r.t_half_initial:.2f} d')
    for t in (0, 25, 50):
        print(f'new/old at t={t:>2}: {r.ratio_new_to_old[int(round(t / p.dt))]:.3f}')
