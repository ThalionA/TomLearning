"""Fixed 5-trial-bin subspace readout over the first 30 trials (Tom's B/D).

Unlike run_trajectory's sliding 15-trial windows, this fits the full pCCA readout in
NON-overlapping 5-trial ordinal bins — trials [1-5],[6-10],[11-15],[16-20],[21-25],[26-30]
— per (animal, pair, bin), for learners AND non-learners. Gives the metrics as a function
of early trial NUMBER (not trial fraction). Exports metrics, per-dim rows (for mincc), and
per-neuron contrib (for the per-area Gini).

Usage: PYTHONPATH=src python scripts/run_trajectory_bins.py --bin-ms 10 --max-lag 5 \
       --smooth-ms 2.5 --tag _bin10 [--include-fs]
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, dataio, subspace_window  # noqa: E402

K, N_FOLDS, N_SHUFFLES, MAX_SAMPLES = 30, 5, config.SURROGATE_SHUFFLES, 8000
BIN_SIZE, N_BINS = 5, 6                     # 5-trial bins, first 30 trials
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]
FIELDS = ["animal", "learner", "pair", "bin", "trial_lo", "trial_hi", "n_bins",
          "cc1", "n_sig", "mi_sig", "ifi", "optimal_lag", "gini_x", "gini_y"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=10)
    p.add_argument("--max-lag", type=int, default=5)
    p.add_argument("--smooth-ms", type=float, default=0.0)
    p.add_argument("--tag", default="")
    p.add_argument("--include-fs", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs,
                              gaussian_sd_ms=args.smooth_ms)
    suffix = f"{args.tag}{'_fsincl' if args.include_fs else ''}"
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    print(f"5-trial bins (first {BIN_SIZE*N_BINS} trials) | bin={args.bin_ms}ms "
          f"smooth={args.smooth_ms}ms K={K} lag+-{args.max_lag} | "
          f"FS {'incl' if args.include_fs else 'excl'}\n")

    rows, dim_rows, wt_rows = [], [], []
    for a in animals:
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids = streams.trial_idx_50ms[run]
        uniq = np.unique(trial_ids)
        if uniq.size < BIN_SIZE:
            continue
        is_learner = a.animal_id in entries
        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run]
        bins = [uniq[i * BIN_SIZE:(i + 1) * BIN_SIZE] for i in range(N_BINS)]
        n_rows = 0
        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            for bi, tb in enumerate(bins):
                if tb.size < BIN_SIZE:
                    continue
                idx = np.where(np.isin(trial_ids, tb))[0]
                if idx.size > MAX_SAMPLES:        # even stride, not first-N: keep ALL trials
                    idx = idx[np.linspace(0, idx.size - 1, MAX_SAMPLES).astype(int)]
                if idx.size < K * 50 or np.unique(trial_ids[idx]).size < N_FOLDS:
                    continue
                groups = trial_ids[idx]
                Zi = Z[idx] if Z is not None else None
                ws = subspace_window.window_subspace(
                    X[idx], Y[idx], groups, Z=Zi, k=K, max_lag=args.max_lag,
                    n_shuffles=N_SHUFFLES, n_folds=N_FOLDS)
                cc1 = float(ws.cc[0]) if ws.cc.size else float("nan")
                if not np.isfinite(cc1):
                    continue
                rows.append({"animal": a.animal_id, "learner": int(is_learner),
                             "pair": f"{ax}-{ay}", "bin": bi + 1,
                             "trial_lo": int(bi * BIN_SIZE + 1), "trial_hi": int((bi + 1) * BIN_SIZE),
                             "n_bins": int(idx.size), "cc1": round(cc1, 4),
                             "n_sig": ws.n_sig, "mi_sig": round(ws.mi_sig, 4),
                             "ifi": round(ws.ifi, 4), "optimal_lag": ws.optimal_lag,
                             "gini_x": round(ws.gini_x, 4), "gini_y": round(ws.gini_y, 4)})
                n_rows += 1
                for d in range(int(ws.cc.size)):
                    dim_rows.append({"animal": a.animal_id, "learner": int(is_learner),
                                     "pair": f"{ax}-{ay}", "bin": bi + 1, "dim": d + 1,
                                     "cc": round(float(ws.cc[d]), 4),
                                     "sig": int(bool(ws.sig_mask[d])) if d < ws.sig_mask.size else 0,
                                     "ifi": round(float(ws.ifi_per_dim[d]), 4)
                                     if d < ws.ifi_per_dim.size else ""})
                for area, W in ((ax, ws.weights_x), (ay, ws.weights_y)):
                    for u, c in enumerate(np.linalg.norm(np.atleast_2d(W), axis=1)):
                        wt_rows.append({"animal": a.animal_id, "learner": int(is_learner),
                                        "pair": f"{ax}-{ay}", "bin": bi + 1, "area": area,
                                        "unit": u, "contrib": round(float(c), 6)})
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): {n_rows} rows")

    R = config.RESULTS_DIR
    for name, data, cols in (
            (f"trajectory_bins{suffix}", rows, FIELDS),
            (f"trajectory_bins_dims{suffix}", dim_rows,
             ["animal", "learner", "pair", "bin", "dim", "cc", "sig", "ifi"]),
            (f"trajectory_bins_weights{suffix}", wt_rows,
             ["animal", "learner", "pair", "bin", "area", "unit", "contrib"])):
        with open(R / f"{name}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
            w.writeheader(); w.writerows(data)
        print(f"wrote {name}.csv ({len(data)} rows)")


if __name__ == "__main__":
    main()
