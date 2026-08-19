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

Usage: python scripts/run_fixed_subspace.py [--include-fs]
       (defaults = config.TEMPORAL; ⚠ this driver's historical cap is --max-per-trial 600,
       i.e. the first 6 s of EVERY trial — kept as the default so numbers reproduce; pass
       --max-per-trial 0 to uncap. See HANDOFF.md §2.)
"""

from __future__ import annotations

import csv

import numpy as np

from _common import (PAIRS, TEMPORAL, animals_filter, cap_desc, cfg_from_args, config,
                     fs_suffix, temporal_parser)
from tom_cca import dataio, fixed_subspace, partial, preprocess

K = TEMPORAL.k
N_DIMS = 3
EPOCHS = list(config.EPOCH_NAMES)
FIELDS = ["animal", "learner", "pair", "epoch", "dim", "bin_ms", "lag_bins", "lag_ms",
          "r", "n_bins", "n_fit_trials", "r_fit"]


def parse_args():
    p = temporal_parser(__doc__.splitlines()[0], shuffles=False)
    p.set_defaults(max_per_trial=600)     # historical cap of THIS driver (see docstring)
    return p.parse_args()


def main():
    args = parse_args()
    cfg = cfg_from_args(args)
    stem = args.out or f"fixed_subspace_bin{args.bin_ms}{fs_suffix(args.include_fs)}"
    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(animals, cfg,
                                     behaviour_lookup=dataio.load_learning_points())
    lags = np.arange(-args.max_lag, args.max_lag + 1)
    print(f"FIXED SUBSPACE | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} | "
          f"lags +/-{args.max_lag} bins | FS={'incl' if args.include_fs else 'excl'} | "
          f"cap: {cap_desc(args)}\n")
    only = animals_filter(args.animals)

    out = config.RESULTS_DIR / f"{stem}.csv"
    rows = []

    def _flush():
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    for a in animals:
        if a.animal_id not in entries:            # epoch split needs a learning point
            continue
        if only is not None and int(a.animal_id) not in only:
            continue
        sess = preprocess.load_running_session(
            a, cfg, entries, max_samples=args.max_samples, max_per_trial=args.max_per_trial)
        if isinstance(sess, str):
            print(f"  animal {a.animal_id}: SKIP — {sess}"); continue
        epoch_of_trial = preprocess.epoch_of_trial(sess.lp, sess.trials_full, cfg, EPOCHS)
        if epoch_of_trial is None:
            print(f"  animal {a.animal_id}: epoch_windows None (skip)"); continue
        fit_trials = fixed_subspace.balanced_trials(epoch_of_trial, EPOCHS, seed=0)
        if not fit_trials:
            print(f"  animal {a.animal_id}: no balanced trial set (skip)"); continue
        trial_ids = sess.trial_ids

        for ax, ay in PAIRS:
            pair = sess.pair(ax, ay)
            if pair is None:
                continue
            X, Y, Z = pair
            # Residualise ONCE, with coefficients estimated on the balanced fit trials
            # only, and apply them to every row. Then fit and project in that SAME space
            # (Z=None below).
            #
            # Until 2026-08-03 this was inconsistent: fit_fixed residualised internally
            # using fit-trial coefficients while the projection used all-trial
            # coefficients, so the frozen weights were applied in a different residual
            # space from the one they were fitted in. The residualisation is part of the
            # frozen map and has to travel with it — otherwise "identical weights across
            # epochs" is not actually an identical transform.
            if Z is not None:
                train = np.isin(trial_ids, fit_trials)
                Xr = partial.partial_out_cv(X, Z, train)
                Yr = partial.partial_out_cv(Y, Z, train)
            else:
                Xr, Yr = X, Y
            fit = fixed_subspace.fit_fixed(Xr, Yr, trial_ids, Z=None, k=K,
                                           trials=fit_trials)
            if fit is None:
                print(f"  fit fail {a.animal_id} {ax}-{ay}"); continue
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
