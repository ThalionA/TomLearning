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

Usage: python scripts/run_lag_subspaces.py [--epochs] [--include-fs]
       (defaults = config.TEMPORAL; ⚠ this driver's historical cap is --max-samples 12000
       --max-per-trial 600 = the first ~20 trials (per epoch in --epochs mode) — kept as
       the default so numbers reproduce; pass --max-samples 0 --max-per-trial 0 to uncap.)
"""

from __future__ import annotations

import csv

import numpy as np

from _common import (PAIRS, TEMPORAL, animals_filter, cap_desc, cfg_from_args, config,
                     fs_suffix, temporal_parser)
from tom_cca import dataio, lag_subspace, membership, preprocess, subspace_window

K = TEMPORAL.k
N_FOLDS = TEMPORAL.n_folds
TAU_BINS = TEMPORAL.label_w_bins  # +/-50 ms at 10 ms bins — the report's headline IFI window
ANGLE_DIMS = 3                # subspace dimensionality for the principal-angle readout
LAGS = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 7, -7, 10, -10, 15, -15, 20, -20,
        25, -25]
EPOCHS = list(config.EPOCH_NAMES)
FIELDS = ["animal", "learner", "pair", "epoch", "bin_ms", "lag_bins", "lag_ms",
          "n_samples",
          "cc1", "cc_mean3", "n_sig", "angle_x", "angle_y", "angle_x_cc1",
          "angle_y_cc1", "floor_x", "floor_y", "floor_x_cc1", "floor_y_cc1",
          "gini_x_conn",
          "gini_y_conn", "angle_ff_fb_x", "angle_ff_fb_y", "angle_ff_fb_cc1"]


def parse_args():
    p = temporal_parser(__doc__.splitlines()[0], shuffles=False, max_lag=False)
    p.set_defaults(max_samples=12000, max_per_trial=600)   # historical cap (see docstring)
    p.add_argument("--epochs", action="store_true",
                   help="resolve every readout by learning epoch (naive/intermediate/"
                        "expert) instead of pooling the session. Learners only — an "
                        "epoch split needs a learning point. Serves meeting items 4 "
                        "and 6 (lagged curves and integration windows, naive vs exp) "
                        "and the FF/FB evolution in item 2.")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = cfg_from_args(args)
    ep = "_epochs" if args.epochs else ""
    stem = args.out or f"lag_subspaces_bin{args.bin_ms}{ep}{fs_suffix(args.include_fs)}"
    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(animals, cfg,
                                     behaviour_lookup=dataio.load_learning_points())
    print(f"LAG SUBSPACES | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} | "
          f"{len(LAGS)} lags, tau={TAU_BINS} bins | "
          f"FS={'incl' if args.include_fs else 'excl'} | cap: {cap_desc(args)}"
          f"{' (applied WITHIN each epoch)' if args.epochs else ''}\n")
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
        # load UNCAPPED; the cap is applied per analysis segment in _segments
        full = preprocess.load_running_session(a, cfg, entries, min_trials=N_FOLDS + 1)
        if isinstance(full, str):
            print(f"  animal {a.animal_id}: SKIP — {full}"); continue
        is_learner = full.learner
        segments = _segments(full, cfg, args)
        if not segments:
            continue

        for epoch, sel in segments:
            sess = full.subset(sel)
            trial_ids = sess.trial_ids
            if np.unique(trial_ids).size < N_FOLDS + 1:
                continue
            n_tr = int(np.unique(trial_ids).size)
            print(f"    {a.animal_id} [{epoch}]: {sel.size} bins from {n_tr} trials "
                  f"(ids {int(trial_ids.min())}-{int(trial_ids.max())})")
            _do_pairs(rows, a, is_learner, epoch, sess, args)
        _flush()
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): "
              f"{len(rows)} rows total")
    _flush()
    print(f"\n-> {out}")


def _segments(full, cfg, args):
    """``[(label, index_into_running_bins)]`` — one entry per analysis segment.

    ⚠ "session" IS A MISNOMER under the historical cap (12 000 bins at <= 600 per
    trial): that "session" segment is the FIRST ~20 TRIALS, not the whole session — for
    an animal whose learning point is at trial 69 that is entirely pre-learning. The
    driver prints the retained trial count so this is visible; pass
    ``--max-samples 0 --max-per-trial 0`` for the true session.

    Epoch mode gives one segment per learning epoch, with the cap applied WITHIN each
    epoch rather than to the session — capping first would spend the whole budget on the
    earliest trials and leave the expert epoch empty.
    """
    tid = full.trial_ids
    if not args.epochs:
        return [("session", preprocess.cap_running_bins(tid, args.max_samples, args.max_per_trial))]
    if not full.learner or full.lp is None:        # an epoch split needs a learning pt
        return []
    eot = preprocess.epoch_of_trial(full.lp, full.trials_full, cfg, EPOCHS)
    if eot is None:
        return []
    segs = []
    for epoch in EPOCHS:
        trials = [t for t, e in eot.items() if e == epoch]
        idx = np.flatnonzero(np.isin(tid, trials))
        if idx.size == 0:
            continue
        segs.append((epoch, idx[preprocess.cap_running_bins(tid[idx], args.max_samples,
                                                             args.max_per_trial)]))
    return segs


def _do_pairs(rows, a, is_learner, epoch, sess, args):
    """Sweep every area pair for one animal-segment, appending rows in place."""
    trial_ids = sess.trial_ids
    for ax, ay in PAIRS:
        pair = sess.pair(ax, ay)
        if pair is None:
            continue
        X, Y, Z = pair
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
        # BOTH d=1 floors. Exporting only the X floor and subtracting it from both
        # areas' angles (as this driver did until 2026-08-03) reintroduces exactly the
        # shared-floor flaw that split_half_floor returns two values to avoid: the floor
        # rises as an area has fewer units, and these pairs are lopsided (CA1 vs SUB).
        floor1_x, floor1_y = lag_subspace.split_half_floor(
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
                "floor_x_cc1": round(floor1_x, 2) if np.isfinite(floor1_x) else "",
                "floor_y_cc1": round(floor1_y, 2) if np.isfinite(floor1_y) else "",
                "gini_x_conn": round(float(gx), 4),
                "gini_y_conn": round(float(gy), 4),
                "angle_ff_fb_x": round(a_ffx, 2) if np.isfinite(a_ffx) else "",
                "angle_ff_fb_y": round(a_ffy, 2) if np.isfinite(a_ffy) else "",
                "angle_ff_fb_cc1": round(a_ff1, 2) if np.isfinite(a_ff1) else "",
            })


if __name__ == "__main__":
    main()
