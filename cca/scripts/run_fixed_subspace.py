"""FIXED-subspace lag curves per epoch — meeting 2026-07-28 items 5 and 7.

    "Go back to also doing one CC and lag it across time."
    "Check with fixed subspaces, identified across all trials, and compare naive/exp."

Identify the communication subspace ONCE per (animal, pair) on a trial set balanced
across the three learning epochs, freeze it, then project each epoch through the
identical weights and lag one canonical component across time. Because the weights never
change, a naive-vs-expert difference in the lag curve is a difference in the ACTIVITY,
not in the fit — which is the confound in every refit-per-epoch contrast the project has
run so far.

The correlations here are IN-SAMPLE by construction (the frozen fit saw every epoch) and
must not be quoted as held-out coupling strengths; `lag_subspaces` is the leak-free arm.
What this buys is a clean epoch CONTRAST, since both epochs pass through one set of
weights and the balanced fit means neither dominates it.

Writes one row per (animal, pair, epoch, dim, lag):
    results/fixed_subspace_bin{10}{,_fsincl}.csv

Usage: PYTHONPATH=src python scripts/run_fixed_subspace.py --bin-ms 10 --smooth-ms 2.5
       [--include-fs]
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, dataio, fixed_subspace, partial  # noqa: E402

K = 30
MAX_LAG = 25                 # +/-250 ms at 10 ms bins, matching run_lag_curves
N_DIMS = 3
MAX_BINS_PER_TRIAL = 600
EPOCHS = ["naive", "intermediate", "expert"]
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]
FIELDS = ["animal", "learner", "pair", "epoch", "dim", "bin_ms", "lag_bins", "lag_ms",
          "r", "n_bins", "n_fit_trials", "r_fit"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=10)
    p.add_argument("--include-fs", action="store_true")
    p.add_argument("--smooth-ms", type=float, default=2.5)
    p.add_argument("--max-lag", type=int, default=MAX_LAG)
    p.add_argument("--out", default="")
    return p.parse_args()


def _cap_bins(trial_ids, max_per_trial=MAX_BINS_PER_TRIAL):
    """Keep a CONSECUTIVE leading block of each trial — preserves bin adjacency for the
    lagged correlation (an even stride would destroy it; see GOTCHAS.md)."""
    keep = [np.where(trial_ids == t)[0][:max_per_trial]
            for t in np.unique(trial_ids)]
    return np.sort(np.concatenate(keep)) if keep else np.array([], dtype=int)


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs,
                              gaussian_sd_ms=args.smooth_ms)
    suffix = "_fsincl" if args.include_fs else ""
    stem = args.out or f"fixed_subspace_bin{args.bin_ms}{suffix}"
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    lags = np.arange(-args.max_lag, args.max_lag + 1)
    print(f"FIXED SUBSPACE | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} | "
          f"lags +/-{args.max_lag} bins | FS={'incl' if args.include_fs else 'excl'}\n")

    out = config.RESULTS_DIR / f"{stem}.csv"
    rows = []

    def _flush():
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    for a in animals:
        if a.animal_id not in entries:            # epoch split needs a learning point
            continue
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids_full = streams.trial_idx_50ms[run]
        uniq = np.unique(trial_ids_full)
        if uniq.size < 6:
            continue
        ew = dataio.epoch_windows(int(entries[a.animal_id].lp), uniq.size, cfg)
        if ew is None:
            print(f"  animal {a.animal_id}: epoch_windows None (skip)"); continue
        cap = _cap_bins(trial_ids_full)
        trial_ids = trial_ids_full[cap]
        epoch_of_trial = {}
        for epoch in EPOCHS:
            pos = ew[epoch]
            pos = pos[(pos >= 0) & (pos < uniq.size)]
            for t in uniq[pos]:
                epoch_of_trial[int(t)] = epoch
        fit_trials = fixed_subspace.balanced_trials(epoch_of_trial, EPOCHS, seed=0)
        if not fit_trials:
            print(f"  animal {a.animal_id}: no balanced trial set (skip)"); continue

        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run][cap]

        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            fit = fixed_subspace.fit_fixed(X, Y, trial_ids, Z=Z, k=K,
                                           trials=fit_trials)
            if fit is None:
                print(f"  fit fail {a.animal_id} {ax}-{ay}"); continue
            # project through the FROZEN weights; residualise once, in-sample, so every
            # epoch sees the identical map (the whole point of this arm)
            Xr = partial.partial_out(X, Z) if Z is not None else X
            Yr = partial.partial_out(Y, Z) if Z is not None else Y
            n_dims = int(min(N_DIMS, fit.wx.shape[1], fit.wy.shape[1]))
            for epoch in EPOCHS:
                trials = [t for t, e in epoch_of_trial.items() if e == epoch]
                m = np.isin(trial_ids, trials)
                if m.sum() < 3 * args.max_lag:
                    continue
                for dim in range(n_dims):
                    u, v = fixed_subspace.project(Xr[m], Yr[m], fit, dim=dim)
                    r = fixed_subspace.variate_lag_curve(u, v, trial_ids[m], lags)
                    for lag, val in zip(lags, r):
                        rows.append({
                            "animal": a.animal_id, "learner": 1,
                            "pair": f"{ax}-{ay}", "epoch": epoch, "dim": dim + 1,
                            "bin_ms": args.bin_ms, "lag_bins": int(lag),
                            "lag_ms": int(lag) * args.bin_ms,
                            "r": round(float(val), 5) if np.isfinite(val) else "",
                            "n_bins": int(m.sum()),
                            "n_fit_trials": fit.n_trials,
                            "r_fit": round(float(fit.r_fit[dim]), 4)
                                     if dim < fit.r_fit.size else "",
                        })
        _flush()
        print(f"  animal {a.animal_id}: {len(rows)} rows total")
    _flush()
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
