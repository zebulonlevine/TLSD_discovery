#!/usr/bin/env python
"""
run_batch_2025_06_06.py
=========================

60S (RPL29) and 40S (RPS24) sgNT vs sgUBA5-LT at 15ms, from
all_data/2025_06_06_zeb/. No joint megapool exists for either pair, but
each single condition's own combinationwise pool already does (combos
sgNT_60S/sgUBA5_60S/sgNT_40S/sgUBA5_40S under tracking_output_q=15p0_
splitsize10/traj_combinations/) -- build a standalone cache per condition
and merge_likelihood_caches instead of running a new saSPT joint fit
(see saspt_likelihood_cache.merge_likelihood_caches docstring for why
this is exactly equivalent).

Runs Gate 3 only (naive + mle_d representations) -- Fork B is a separate,
opt-in escalation.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pipeline')))  # the modules under test live in ../pipeline
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
from saspt_likelihood_cache import (
    build_likelihood_cache, load_likelihood_cache, merge_likelihood_caches,
)
from occupancy_compare import compare_groups

DATASET_ROOT = tlsd_paths.spt('all_data', '2025_06_06_zeb')
TO = os.path.join(DATASET_ROOT, "tracking_output_q=15p0_splitsize10")
SETTINGS_PKL = os.path.join(TO, "settings.pkl")
COMBO_DIR = os.path.join(TO, "traj_combinations")
CACHE_DIR = os.path.join(TO, "occupancy_compare_cache")
RUNS_DIR = os.path.join(TO, "occupancy_compare_runs", "batch_2025_06_06")

# (cache_name, [(member_label, combo_name), ...], (name_a, name_b))
JOBS = [
    ("60S_bothgeno", [("sgNT_60S", "sgNT_60S"), ("sgUBA5_60S", "sgUBA5_60S")],
     ("sgNT_60S", "sgUBA5_60S")),
    ("40S_bothgeno", [("sgNT_40S", "sgNT_40S"), ("sgUBA5_40S", "sgUBA5_40S")],
     ("sgNT_40S", "sgUBA5_40S")),
]


def main():
    os.makedirs(RUNS_DIR, exist_ok=True)

    for cache_name, members, (name_a, name_b) in JOBS:
        print("=" * 70, flush=True)
        print(f"cache: {cache_name}", flush=True)
        print("=" * 70, flush=True)

        member_caches = {}
        for label, combo in members:
            npz_path = build_likelihood_cache(
                pooled_traj_csv=os.path.join(COMBO_DIR, f"{combo}_pooled_traj.csv"),
                settings_pkl=SETTINGS_PKL, out_dir=CACHE_DIR,
                fov_map_csv=os.path.join(COMBO_DIR, f"{combo}_fov_map.csv"),
                cache_name=f"_standalone_{combo}")
            member_caches[label] = npz_path

        merged_path = merge_likelihood_caches(
            member_caches, out_dir=CACHE_DIR, cache_name=cache_name)
        cache = load_likelihood_cache(merged_path, settings_pkl=SETTINGS_PKL)

        for representation in ("naive", "mle_d"):
            tag = f"{name_a}_vs_{name_b}"
            out_dir = os.path.join(RUNS_DIR, tag)
            t0 = time.time()
            result = compare_groups(
                cache, group_a={"source_member": name_a},
                group_b={"source_member": name_b},
                name_a=name_a, name_b=name_b, out_dir=out_dir,
                n_perm=2000, seed=0, representation=representation)
            print(f"  [{tag}] representation={representation}: "
                  f"{len(result['excursion_intervals'])} excursion(s), "
                  f"{time.time()-t0:.1f}s", flush=True)

    print("ALL DONE", flush=True)


if __name__ == '__main__':
    main()
