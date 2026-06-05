"""Frame A driver: within-animal learning trajectory of the communication subspace.

Per learner, slide a window over trials; in each window fit pCCA on the
continuous running bins (control all other recorded areas), read the dominant
held-out canonical correlation, and record it against the window-centre trial
fraction (0 = first engaged trial, 1 = last). Per (animal, pair) fit a slope of
CC1 vs trial fraction; across animals test whether the slope is consistently
signed (the powered alternative to the 3-epoch contrast).

Writes:
* results/trajectory_windows.csv  -- one row per (animal, pair, window)
* results/trajectory_slopes.csv   -- one row per (animal, pair): slope, r
and prints the per-pair across-animal sign test.

Continuous regime (validated by prototype_continuous_pcca.py): 25 ms bins,
PCA->K, whole-trial CV. No depth/ISI/waveforms; FS via idx_fs only.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import (config, dataio, paired_stats, partial,  # noqa: E402
                     trajectory)

K = 30
N_FOLDS = 5
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=25)
    p.add_argument("--window", type=int, default=30, help="trials per window")
    p.add_argument("--step", type=int, default=15, help="trial step between windows")
    p.add_argument("--min-window", type=int, default=20)
    return p.parse_args()


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms)
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    learners = [a for a in animals if a.animal_id in entries]
    thr = cfg.velocity_thresh_cm_s

    win_rows = []
    slope_rows = []
    print(f"Frame A: bin={args.bin_ms}ms, window={args.window} trials, "
          f"step={args.step}, PCA K={K}, pCCA, {N_FOLDS}-fold CV\n")

    for a in learners:
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids = streams.trial_idx_50ms[run]
        uniq = np.unique(trial_ids)
        if uniq.size < args.min_window:
            continue
        t0, t1 = float(uniq.min()), float(uniq.max())
        span = (t1 - t0) if t1 > t0 else 1.0
        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run]
        windows = trajectory.sliding_windows(uniq, args.window, args.step,
                                             args.min_window)
        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            centers, ccs = [], []
            for w in windows:
                mask = np.isin(trial_ids, w)
                if mask.sum() < K * 50:                  # keep the sample regime
                    continue
                Xi, Yi = X[mask], Y[mask]
                if Z is not None:
                    Zi = Z[mask]
                    Xi = partial.partial_out(Xi, Zi)
                    Yi = partial.partial_out(Yi, Zi)
                cc_te, _ = trajectory.heldout_cc1(
                    Xi, Yi, trial_ids[mask], K, N_FOLDS)
                frac = (float(np.mean(w)) - t0) / span
                centers.append(frac)
                ccs.append(cc_te)
                win_rows.append({"animal": a.animal_id, "pair": f"{ax}-{ay}",
                                 "trial_frac": round(frac, 4),
                                 "n_bins": int(mask.sum()),
                                 "cc1": round(cc_te, 4) if np.isfinite(cc_te) else ""})
            if sum(np.isfinite(ccs)) >= 4:
                slope, r = trajectory.linear_slope(centers, ccs)
                slope_rows.append({"animal": a.animal_id, "pair": f"{ax}-{ay}",
                                   "n_windows": int(sum(np.isfinite(ccs))),
                                   "slope": slope, "r": r,
                                   "cc1_first": next((c for c in ccs if np.isfinite(c)), np.nan),
                                   "cc1_last": next((c for c in reversed(ccs) if np.isfinite(c)), np.nan)})
        print(f"  animal {a.animal_id}: {len([s for s in slope_rows if s['animal']==a.animal_id])} pair-trajectories")

    rdir = config.RESULTS_DIR
    with open(rdir / "trajectory_windows.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["animal", "pair", "trial_frac", "n_bins", "cc1"],
                           lineterminator="\n")
        w.writeheader(); w.writerows(win_rows)
    with open(rdir / "trajectory_slopes.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["animal", "pair", "n_windows", "slope", "r",
                                          "cc1_first", "cc1_last"], lineterminator="\n")
        w.writeheader(); w.writerows(slope_rows)
    print(f"\nwrote {rdir/'trajectory_windows.csv'} ({len(win_rows)} rows)")
    print(f"wrote {rdir/'trajectory_slopes.csv'} ({len(slope_rows)} rows)\n")

    # Across-animal sign test on per-animal slopes, per pair.
    by_pair = defaultdict(list)
    for s in slope_rows:
        if np.isfinite(s["slope"]):
            by_pair[s["pair"]].append(s["slope"])
    print("Across-animal trajectory test (does CA1/... communication track learning?):")
    print(f"  {'pair':9s} {'n_animals':>9s} {'median_slope':>13s} {'n_up':>5s} {'p_signrank':>11s}")
    for pair in [f"{x}-{y}" for x, y in PAIRS]:
        sl = np.array(by_pair.get(pair, []))
        if sl.size < 3:
            continue
        n_up = int(np.sum(sl > 0))
        _, med, _, p = paired_stats.wilcoxon_signed(sl.tolist())
        print(f"  {pair:9s} {sl.size:>9d} {med:>+13.4f} {n_up:>3d}/{sl.size:<2d} {p:>11.4g}")


if __name__ == "__main__":
    main()
