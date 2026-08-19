"""Directionality with a SWEPT lag-integration window — which window is cleanest?

For each animal x area-pair, fit the session subspace (residualise the third area,
PCA->K on all running bins), then compute the **held-out, segment-aware** lag curve
(`lagged.heldout_lag_curve_flat` — CCA refit per fold at each lag; the lag never
crosses a trial boundary, so neither the in-sample bias nor the running-bin
concatenation contaminates it). From that curve the IFI is integrated over every
window |lag| <= w (`lagged.ifi_by_window`). Writes one row per (animal, pair, w):

    results/ifi_windows_bin{10,25,50}{,_fsincl}.csv

so analyze_ifi_windows.py can ask, per pair, which integration window gives the most
consistent across-animal directionality (narrow misses an off-zero peak; wide dilutes
with far-lag noise). Run at both bin widths. Moderately expensive (refits CCA at every
lag x fold on the whole session) — launch after the trajectory/KCCA batch.
"""

from __future__ import annotations

import csv

import numpy as np

from _common import PAIRS, TEMPORAL, animals_filter, cfg_from_args, config, fs_suffix, temporal_parser
from tom_cca import dataio, lagged, preprocess

K = TEMPORAL.k
N_FOLDS = TEMPORAL.n_folds
FIELDS = ["animal", "learner", "pair", "bin_ms", "n_bins", "window_bins",
          "window_ms", "ifi", "cc_pos", "cc_neg", "peak_lag_bins", "peak_lag_ms"]


def parse_args():
    p = temporal_parser(__doc__.splitlines()[0], shuffles=False)
    p.set_defaults(max_lag=12)          # windows run w = 1..max_lag (historical default)
    return p.parse_args()


def main():
    args = parse_args()
    cfg = cfg_from_args(args)
    stem = args.out or f"ifi_windows_bin{args.bin_ms}{fs_suffix(args.include_fs)}"
    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(animals, cfg,
                                     behaviour_lookup=dataio.load_learning_points())
    print(f"IFI-WINDOW sweep | bin={args.bin_ms}ms K={K} max_lag={args.max_lag} "
          f"({args.max_lag * args.bin_ms} ms) | FS={'incl' if args.include_fs else 'excl'}\n")
    only = animals_filter(args.animals)

    out = config.RESULTS_DIR / f"{stem}.csv"
    rows = []

    def _flush():
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    for a in animals:
        if only is not None and int(a.animal_id) not in only:
            continue
        sess = preprocess.load_running_session(
            a, cfg, entries, max_samples=args.max_samples, max_per_trial=args.max_per_trial,
            min_trials=N_FOLDS + 1)
        if isinstance(sess, str):
            print(f"  animal {a.animal_id}: SKIP — {sess}"); continue
        is_learner = sess.learner
        trial_ids = sess.trial_ids

        for ax, ay in PAIRS:
            pair = sess.pair(ax, ay)
            if pair is None:
                continue
            X, Y, Z = pair
            # full-session residualise + PCA -> scores (PCA fit once; the CCA in the
            # lag curve is held-out per fold, so the directionality CC is unbiased)
            Sx, Sy = preprocess.residual_pca_scores(X, Y, Z, K)
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
                # a side with NO finite CC is written blank (as before), not 0.0
                has_pos = np.any(np.isfinite(cc[mask & (lags > 0)]))
                has_neg = np.any(np.isfinite(cc[mask & (lags < 0)]))
                rows.append({
                    "animal": a.animal_id, "learner": int(is_learner),
                    "pair": f"{ax}-{ay}", "bin_ms": args.bin_ms,
                    "n_bins": int(trial_ids.size), "window_bins": w,
                    "window_ms": w * args.bin_ms,
                    "ifi": round(float(ifi_w[w - 1]), 4) if np.isfinite(ifi_w[w - 1]) else "",
                    "cc_pos": round(cc_pos, 4) if has_pos else "",
                    "cc_neg": round(cc_neg, 4) if has_neg else "",
                    "peak_lag_bins": peak_lag, "peak_lag_ms": peak_lag * args.bin_ms})
        _flush()
        print(f"  animal {a.animal_id}: done ({len(rows)} rows)")
    _flush()
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
