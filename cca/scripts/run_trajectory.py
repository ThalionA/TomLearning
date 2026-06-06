"""Frame A driver: FULL subspace readout per window over each animal's trials.

For every learner (and non-learner, for control), slide a window over trials and
in each window fit pCCA on the continuous running bins (control all other areas),
then read the full subspace suite via ``subspace_window``:
held-out CC1, n_sig, mi_sig, IFI, optimal lag, Gini (x/y), plus cross-window
principal-angle ROTATION and Jaccard membership OVERLAP relative to the previous
window. Each row is also tagged with three learning axes (trial fraction,
behavioural performance, LP-relative trial position) and the learner flag.

This driver only WRITES the rich per-window table (results/trajectory_windows.csv)
-- the expensive part. Slopes / across-animal sign tests / learner splits are
computed cheaply and re-runnably by scripts/analyze_trajectory.py from that CSV.

Continuous regime (validated by prototype_continuous_pcca.py): 25 ms, PCA->K,
whole-trial CV, drop saturated (CC>=0.99) windows. No depth/ISI/waveforms.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import (config, dataio, membership, partial,  # noqa: E402
                     subspace, subspace_window, trajectory)

K = 30
N_FOLDS = 5
MAX_LAG = 10                  # +-10 bins = +-250 ms at 25 ms
N_SHUFFLES = 20
SAT = 0.99
MAX_SAMPLES = 6000            # cap bins/window: ~50x>K for a stable fit, keeps
                              # the contiguous lag structure, ~7x faster fits
ROT_DIMS = 3                  # dims used for the rotation (dominant subspace)
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=25)
    p.add_argument("--window", type=int, default=30)
    p.add_argument("--step", type=int, default=15)
    p.add_argument("--min-window", type=int, default=20)
    p.add_argument("--include-fs", action="store_true",
                   help="keep fast-spiking units (default: FS-excluded in "
                        "V1/RSC/CA1/CA3 per Tom convention); writes a _fsincl CSV")
    return p.parse_args()


def _performance(animal) -> np.ndarray:
    try:
        with h5py.File(animal.streams_path, "r") as f:
            return np.asarray(f["analysis_behaviour"]["lick_ratio"]).ravel()
    except Exception:
        return np.array([])


def _rotation(prev_w, cur_w, d_use=ROT_DIMS) -> float:
    """Max principal angle (deg) between the dominant subspaces of two windows."""
    if prev_w is None:
        return float("nan")
    d = int(min(d_use, prev_w.shape[1], cur_w.shape[1]))
    if d < 1:
        return float("nan")
    qa, _ = np.linalg.qr(prev_w[:, :d])
    qb, _ = np.linalg.qr(cur_w[:, :d])
    ang = subspace.principal_angles(qa, qb)
    return float(np.degrees(np.nanmax(ang)))


FIELDS = ["animal", "learner", "pair", "center_trial", "trial_frac",
          "performance", "lp_rel", "n_bins", "cc1", "n_sig", "mi_sig", "ifi",
          "optimal_lag", "gini_x", "gini_y", "rot_x", "rot_y", "jac_x", "jac_y",
          "sh_x", "sh_y"]


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs)
    suffix = "_fsincl" if args.include_fs else ""
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    n_learn = sum(1 for a in animals if a.animal_id in entries)
    print(f"Frame A FULL suite | ALL {len(animals)} animals ({n_learn} learners) "
          f"| bin={args.bin_ms}ms window={args.window} step={args.step} "
          f"K={K} pCCA lag+-{MAX_LAG} shuffles={N_SHUFFLES}\n")

    out = config.RESULTS_DIR / f"trajectory_windows{suffix}.csv"

    def _write(rows):
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    rows = []
    for a in animals:
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
        is_learner = a.animal_id in entries
        lp = float(entries[a.animal_id].lp) if is_learner else float("nan")
        perf = _performance(a)
        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run]
        windows = trajectory.sliding_windows(uniq, args.window, args.step,
                                             args.min_window)
        n_rows_animal = 0
        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            pwx = pwy = pmx = pmy = None
            for w in windows:
                idx_all = np.where(np.isin(trial_ids, w))[0]
                if idx_all.size < K * 50:
                    continue
                # contiguous cap (keeps lag structure), but grow until it spans
                # enough whole trials for the held-out CV folds.
                cap = MAX_SAMPLES
                idx = idx_all[:cap]
                while (np.unique(trial_ids[idx]).size < N_FOLDS + 1
                       and cap < idx_all.size):
                    cap += 2000
                    idx = idx_all[:cap]
                groups = trial_ids[idx]
                Xi, Yi = X[idx], Y[idx]
                if Z is not None:
                    Zi = Z[idx]
                    Xi = partial.partial_out(Xi, Zi)
                    Yi = partial.partial_out(Yi, Zi)
                ws = subspace_window.window_subspace(
                    Xi, Yi, groups, k=K, max_lag=MAX_LAG,
                    n_shuffles=N_SHUFFLES, n_folds=N_FOLDS)
                cc1 = float(ws.cc[0]) if ws.cc.size else float("nan")
                if not np.isfinite(cc1) or cc1 >= SAT:
                    continue
                center = float(np.mean(w))
                rows.append({
                    "animal": a.animal_id, "learner": int(is_learner),
                    "pair": f"{ax}-{ay}", "center_trial": round(center, 1),
                    "trial_frac": round((center - t0) / span, 4),
                    "performance": round(float(np.nanmean(perf[(w - 1).astype(int)])), 4)
                    if perf.size else "",
                    "lp_rel": round(center - lp, 1) if np.isfinite(lp) else "",
                    "n_bins": int(idx.size), "cc1": round(cc1, 4),
                    "n_sig": ws.n_sig, "mi_sig": round(ws.mi_sig, 4),
                    "ifi": round(ws.ifi, 4), "optimal_lag": ws.optimal_lag,
                    "gini_x": round(ws.gini_x, 4), "gini_y": round(ws.gini_y, 4),
                    "rot_x": round(_rotation(pwx, ws.weights_x), 2),
                    "rot_y": round(_rotation(pwy, ws.weights_y), 2),
                    "jac_x": round(membership.jaccard(pmx, ws.member_x), 4)
                    if pmx is not None else "",
                    "jac_y": round(membership.jaccard(pmy, ws.member_y), 4)
                    if pmy is not None else "",
                    "sh_x": round(ws.split_half_x, 2) if np.isfinite(ws.split_half_x) else "",
                    "sh_y": round(ws.split_half_y, 2) if np.isfinite(ws.split_half_y) else "",
                })
                n_rows_animal += 1
                pwx, pwy = ws.weights_x, ws.weights_y
                pmx, pmy = ws.member_x, ws.member_y
        _write(rows)                                  # incremental: survive a kill
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): "
              f"{n_rows_animal} window-rows ({len(rows)} total)", flush=True)

    _write(rows)
    print(f"\nwrote {config.RESULTS_DIR/'trajectory_windows.csv'} ({len(rows)} rows). "
          "Analyse with scripts/analyze_trajectory.py", flush=True)


if __name__ == "__main__":
    main()
