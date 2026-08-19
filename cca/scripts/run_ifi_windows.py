"""Directionality with a SWEPT lag-integration window — which window is cleanest?

For each animal x area-pair, fit the session subspace (residualise the third area,
PCA->K on all running bins), then compute the **held-out, segment-aware** lag curve
(`lagged.heldout_lag_curve_flat` — CCA refit per fold at each lag; the lag never
crosses a trial boundary, so neither the in-sample bias nor the running-bin
concatenation contaminates it). From that curve the IFI is integrated over every
window |lag| <= w (`lagged.ifi_by_window`). Writes one row per (animal, pair, w):

    results/ifi_windows_bin{25,50}{,_fsincl}.csv

so analyze_ifi_windows.py can ask, per pair, which integration window gives the most
consistent across-animal directionality (narrow misses an off-zero peak; wide dilutes
with far-lag noise). Run at both bin widths. Moderately expensive (refits CCA at every
lag x fold on the whole session) — launch after the trajectory/KCCA batch.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, core, dataio, lagged, partial  # noqa: E402

K = 30
N_FOLDS = 5
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]
FIELDS = ["animal", "learner", "pair", "bin_ms", "n_bins", "window_bins",
          "window_ms", "ifi", "cc_pos", "cc_neg", "peak_lag_bins", "peak_lag_ms"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=25)
    p.add_argument("--include-fs", action="store_true")
    p.add_argument("--max-lag", type=int, default=12,
                   help="widest lag (bins) of the sweep; windows run w=1..max-lag")
    p.add_argument("--out", default="", help="output CSV stem override")
    p.add_argument("--smooth-ms", type=float, default=0.0,
                   help="Gaussian s.d. (ms) for spike-train smoothing (Buzsáki = 2.5)")
    return p.parse_args()



def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs,
                              gaussian_sd_ms=args.smooth_ms)
    suffix = "_fsincl" if args.include_fs else ""
    stem = args.out or f"ifi_windows_bin{args.bin_ms}{suffix}"
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    print(f"IFI-WINDOW sweep | bin={args.bin_ms}ms K={K} max_lag={args.max_lag} "
          f"({args.max_lag * args.bin_ms} ms) | FS={'incl' if args.include_fs else 'excl'}\n")

    out = config.RESULTS_DIR / f"{stem}.csv"
    rows = []

    def _flush():
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    for a in animals:
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids = streams.trial_idx_50ms[run]
        if np.unique(trial_ids).size < N_FOLDS + 1:
            continue
        is_learner = a.animal_id in entries
        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run]

        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            # full-session residualise + PCA -> scores (PCA fit once; the CCA in the
            # lag curve is held-out per fold, so the directionality CC is unbiased)
            Xr = partial.partial_out(X, Z) if Z is not None else X
            Yr = partial.partial_out(Y, Z) if Z is not None else Y
            Sx, _ = core.pca_scores(Xr, min(K, Xr.shape[1]))
            Sy, _ = core.pca_scores(Yr, min(K, Yr.shape[1]))
            try:
                lags, cc = lagged.heldout_lag_curve_flat(
                    Sx, Sy, trial_ids, max_lag=args.max_lag, n_folds=N_FOLDS)
            except Exception as e:
                print(f"  fail {a.animal_id} {ax}-{ay}: {e}"); continue
            ifi_w = lagged.ifi_by_window(lags, cc)          # w = 1..max_lag
            peak_lag = int(lags[np.nanargmax(cc)]) if np.any(np.isfinite(cc)) else 0
            for w in range(1, ifi_w.size + 1):
                mask = np.abs(lags) <= w
                cc_pos, cc_neg = lagged.ifi_sides(lags[mask], cc[mask])
                rows.append({
                    "animal": a.animal_id, "learner": int(is_learner),
                    "pair": f"{ax}-{ay}", "bin_ms": args.bin_ms,
                    "n_bins": int(trial_ids.size), "window_bins": w,
                    "window_ms": w * args.bin_ms,
                    "ifi": round(float(ifi_w[w - 1]), 4) if np.isfinite(ifi_w[w - 1]) else "",
                    "cc_pos": round(cc_pos, 4) if np.isfinite(cc_pos) else "",
                    "cc_neg": round(cc_neg, 4) if np.isfinite(cc_neg) else "",
                    "peak_lag_bins": peak_lag, "peak_lag_ms": peak_lag * args.bin_ms})
        _flush()
        print(f"  animal {a.animal_id}: done ({len(rows)} rows)")
    _flush()
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
