"""Lagged communication SUBSPACES — meeting 2026-07-28 items 2, 3 and 4.

`run_lag_curves.py` exports canonical *correlations* at each lag; this exports the
*subspace* at each lag, so two questions become answerable:

  item 3  How stable/similar are the CCs across time lags?
          -> principal angle between the lag-0 subspace and each lagged subspace,
             against the split-half noise floor measured at lag 0. An angle at or below
             the floor means "the same subspace, read at a different delay"; an angle
             above it means the coupling at that delay is carried by different neurons.

  items 2/4  Separate the subspace into feedforward and feedback.
          -> FF = the fit at +TAU_BINS (X leads), FB = the fit at -TAU_BINS (Y leads),
             TAU = the report's headline IFI window (+/-50 ms). Defining FF/FB by the FIT
             LAG rather than by the sign of a per-dim IFI is deliberate: `perdim_ifi`
             showed the significant tail is dominated by dims sitting at the held-out CC
             floor, where the IFI sign is near a coin flip, so a per-dim tag would be
             noise-splitting followed by a circular test.

Gini is the CONNECTION-SPECIFIC definition (`subspace_contribution_connection`), not the
partner-invariant one -- see GOTCHAS.md (2026-07-28).

Mirrors run_lag_curves.py's loading exactly: engaged running bins, consecutive
within-trial cap, partial against all other recorded areas, PCA K=30, 5-fold by-trial CV.

Writes one row per (animal, pair, lag):
    results/lag_subspaces_bin{10}{,_fsincl}.csv

Usage: PYTHONPATH=src python scripts/run_lag_subspaces.py --bin-ms 10 --smooth-ms 2.5
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
from tom_cca import (config, dataio, lag_subspace, membership,  # noqa: E402
                     subspace_window)

K = 30
N_FOLDS = 5
MAX_SAMPLES = 12000
MAX_BINS_PER_TRIAL = 600
TAU_BINS = 5                  # +/-50 ms at 10 ms bins — the report's headline IFI window
ANGLE_DIMS = 3                # subspace dimensionality for the principal-angle readout
LAGS = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 7, -7, 10, -10, 15, -15, 20, -20,
        25, -25]
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]
EPOCHS = ["naive", "intermediate", "expert"]
FIELDS = ["animal", "learner", "pair", "epoch", "bin_ms", "lag_bins", "lag_ms",
          "n_samples",
          "cc1", "cc_mean3", "n_sig", "angle_x", "angle_y", "angle_x_cc1",
          "angle_y_cc1", "floor_x", "floor_y", "floor_cc1", "gini_x_conn",
          "gini_y_conn", "angle_ff_fb_x", "angle_ff_fb_y", "angle_ff_fb_cc1"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=10)
    p.add_argument("--include-fs", action="store_true")
    p.add_argument("--smooth-ms", type=float, default=2.5)
    p.add_argument("--epochs", action="store_true",
                   help="resolve every readout by learning epoch (naive/intermediate/"
                        "expert) instead of pooling the session. Learners only — an "
                        "epoch split needs a learning point. Serves meeting items 4 "
                        "and 6 (lagged curves and integration windows, naive vs exp) "
                        "and the FF/FB evolution in item 2.")
    p.add_argument("--out", default="")
    return p.parse_args()


def _capped_index(trial_ids, max_samples=MAX_SAMPLES,
                  max_per_trial=MAX_BINS_PER_TRIAL):
    """CONSECUTIVE within-trial blocks of whole trials until ~max_samples — preserves
    the bin adjacency the segment-aware lag pairing needs (see GOTCHAS.md)."""
    keep, total = [], 0
    for t in np.unique(trial_ids):
        pos = np.where(trial_ids == t)[0][:max_per_trial]
        keep.append(pos)
        total += pos.size
        if total >= max_samples:
            break
    return np.sort(np.concatenate(keep))


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs,
                              gaussian_sd_ms=args.smooth_ms)
    suffix = "_fsincl" if args.include_fs else ""
    ep = "_epochs" if args.epochs else ""
    stem = args.out or f"lag_subspaces_bin{args.bin_ms}{ep}{suffix}"
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    print(f"LAG SUBSPACES | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} | "
          f"{len(LAGS)} lags, tau={TAU_BINS} bins | "
          f"FS={'incl' if args.include_fs else 'excl'}\n")

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
        trial_ids_full = streams.trial_idx_50ms[run]
        if np.unique(trial_ids_full).size < N_FOLDS + 1:
            continue
        is_learner = a.animal_id in entries
        segments = _segments(a, trial_ids_full, entries, cfg, args.epochs)
        if not segments:
            continue
        area_full = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                area_full[area] = m[run]

        for epoch, sel in segments:
            trial_ids = trial_ids_full[sel]
            if np.unique(trial_ids).size < N_FOLDS + 1:
                continue
            present = {k: v[sel] for k, v in area_full.items()}
            _do_pairs(rows, a, is_learner, epoch, present, trial_ids, args)
        _flush()
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): "
              f"{len(rows)} rows total")
    _flush()
    print(f"\n-> {out}")


def _segments(a, trial_ids_full, entries, cfg, by_epoch: bool):
    """``[(label, index_into_running_bins)]`` — one entry per analysis segment.

    Session mode gives a single capped segment. Epoch mode gives one per learning
    epoch, and the sample cap is applied WITHIN each epoch rather than to the session:
    capping first would spend the whole budget on the earliest trials and leave the
    expert epoch empty.
    """
    if not by_epoch:
        return [("session", _capped_index(trial_ids_full))]
    if a.animal_id not in entries:                 # an epoch split needs a learning pt
        return []
    uniq = np.unique(trial_ids_full)
    ew = dataio.epoch_windows(int(entries[a.animal_id].lp), uniq.size, cfg)
    if ew is None:
        return []
    segs = []
    for epoch in EPOCHS:
        pos = ew[epoch]
        pos = pos[(pos >= 0) & (pos < uniq.size)]
        idx = np.where(np.isin(trial_ids_full, uniq[pos]))[0]
        if idx.size == 0:
            continue
        segs.append((epoch, idx[_capped_index(trial_ids_full[idx])]))
    return segs


def _do_pairs(rows, a, is_learner, epoch, present, trial_ids, args):
    """Sweep every area pair for one animal-segment, appending rows in place."""
    for ax, ay in PAIRS:
        if ax not in present or ay not in present:
            continue
        X, Y = present[ax], present[ay]
        others = [present[z] for z in present if z not in (ax, ay)]
        Z = np.concatenate(others, axis=1) if others else None
        tag = f"{a.animal_id} {ax}-{ay} [{epoch}]"
        try:
            fits = lag_subspace.lag_sweep(X, Y, trial_ids, Z=Z, lags=LAGS, k=K,
                                          n_folds=N_FOLDS)
        except Exception as e:
            print(f"  sweep fail {tag}: {e}"); continue
        if 0 not in fits:
            print(f"  no lag-0 fit {tag}"); continue
        ref = fits[0]
        try:
            ws = subspace_window.window_subspace(
                X, Y, trial_ids, Z=Z, k=K, max_lag=5,
                n_shuffles=config.SURROGATE_SHUFFLES, n_folds=N_FOLDS)
            n_sig = int(np.sum(np.asarray(ws.sig_mask, dtype=bool)))
        except Exception:
            n_sig = 0
        floor_x, floor_y = lag_subspace.split_half_floor(
            X, Y, trial_ids, Z=Z, lag=0, k=K, d_use=ANGLE_DIMS)
        floor1_x, _ = lag_subspace.split_half_floor(
            X, Y, trial_ids, Z=Z, lag=0, k=K, d_use=1)
        # items 2/4: is the FEEDFORWARD subspace (+tau, X leads) the same set of
        # neurons as the FEEDBACK one (-tau)? Compared against the same lag-0
        # split-half floor — an FF/FB angle at the floor means one subspace read at
        # two delays, not two separable directions of flow.
        ff, fb = fits.get(TAU_BINS), fits.get(-TAU_BINS)
        if ff is not None and fb is not None:
            a_ffx = lag_subspace.subspace_angle(ff.wx, fb.wx, ANGLE_DIMS)
            a_ffy = lag_subspace.subspace_angle(ff.wy, fb.wy, ANGLE_DIMS)
            a_ff1 = lag_subspace.subspace_angle(ff.wx, fb.wx, 1)
        else:
            a_ffx = a_ffy = a_ff1 = float("nan")
        for lag in LAGS:
            f = fits.get(lag)
            if f is None:
                continue
            cc = np.asarray(f.cc, dtype=float)
            gx = membership.gini(membership.subspace_contribution_connection(
                f.wx, cc))
            gy = membership.gini(membership.subspace_contribution_connection(
                f.wy, cc))
            rows.append({
                "animal": a.animal_id, "learner": int(is_learner),
                "pair": f"{ax}-{ay}", "epoch": epoch, "bin_ms": args.bin_ms,
                "lag_bins": int(lag), "lag_ms": int(lag) * args.bin_ms,
                "n_samples": f.n_samples,
                "cc1": round(float(cc[0]), 4) if cc.size else "",
                "cc_mean3": round(float(np.nanmean(cc[:3])), 4) if cc.size else "",
                "n_sig": n_sig,
                "angle_x": round(lag_subspace.subspace_angle(ref.wx, f.wx,
                                                             ANGLE_DIMS), 2),
                "angle_y": round(lag_subspace.subspace_angle(ref.wy, f.wy,
                                                             ANGLE_DIMS), 2),
                "angle_x_cc1": round(lag_subspace.subspace_angle(ref.wx, f.wx, 1), 2),
                "angle_y_cc1": round(lag_subspace.subspace_angle(ref.wy, f.wy, 1), 2),
                "floor_x": round(floor_x, 2) if np.isfinite(floor_x) else "",
                "floor_y": round(floor_y, 2) if np.isfinite(floor_y) else "",
                "floor_cc1": round(floor1_x, 2) if np.isfinite(floor1_x) else "",
                "gini_x_conn": round(float(gx), 4),
                "gini_y_conn": round(float(gy), 4),
                "angle_ff_fb_x": round(a_ffx, 2) if np.isfinite(a_ffx) else "",
                "angle_ff_fb_y": round(a_ffy, 2) if np.isfinite(a_ffy) else "",
                "angle_ff_fb_cc1": round(a_ff1, 2) if np.isfinite(a_ff1) else "",
            })


if __name__ == "__main__":
    main()
