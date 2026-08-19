"""Item 2 — label each CC feedforward/feedback ONCE on all trials, then follow those
same CCs across the learning epochs.

    "Separate subspaces into FF/FB meant: find CCs which have positive or negative IFIs
     and check their changes with learning specifically."

Why the axes must be frozen. A per-epoch refit reorders canonical dimensions, so "CC2 in
naive" and "CC2 in expert" are different components and *tracking* a CC is undefined. The
subspace is therefore identified once, frozen, and every epoch projected through the
identical map. Labelling once also avoids the trap item 1 exposed: at +/-50 ms the sign
of a CC's IFI is close to a coin flip, so re-labelling per epoch would manufacture
apparent FF->FB "switching" out of noise.

Pipeline per (animal, pair):

1. FIT SET = **every running bin of every trial** — no sample cap. The other drivers cap
   at ~12k bins because they refit CCA 255 times per pair; this one fits once. The cap
   was not merely unnecessary here, it was harmful: it kept only the first 0.4-1.8 s of
   each ~30 s trial, and running speed rises twice as much naive->expert on those
   trial-onset bins (+12.0 cm/s) as it does over whole trials (+6.6 cm/s).
2. Residualise the third areas ONCE with fit-set coefficients, then fit and project in
   that same residual space (the frozen map includes its residualisation -- getting this
   wrong was the August fit/projection bug).
3. Freeze the CCA. Per-dimension significance from `fixed_subspace.frozen_perm_null` —
   a circular-shift permutation test on THE FROZEN FIT'S OWN canonical correlations,
   rank for rank, so observed and null are the same estimator and no cross-fit
   correspondence exists to mis-attribute.
4. LABEL each dim from the whole fit set: project -> lagged correlation -> IFI at
   +/-LABEL_W bins -> sign. FF = +IFI (first-named area leads), FB = -IFI.
5. TRACK: project each epoch's bins through the SAME axes and export the lagged
   correlation curve per dim. A leave-epoch-out label (`label_loo`) is exported too: the
   all-trials label shares data with every epoch, so an epoch's own IFI sign matches it
   ~60% of the time by construction and cannot be tested against 50%.

Positive lag => the FIRST-named area leads.

⚠ The correlations are IN-SAMPLE by construction (the frozen fit saw every epoch). They
are contrast statistics, not coupling strengths; the leak-free numbers live in the
`lag_curves` arm. What the freeze buys is that an epoch difference is a difference in the
ACTIVITY rather than in the fit.

Writes one row per (animal, pair, dim, epoch, lag):
    results/cc_label_track_bin10{,_fsincl}.csv

Usage: python scripts/run_cc_label_track.py [--include-fs]   (defaults = config.TEMPORAL; no cap)
"""

from __future__ import annotations

import csv

import numpy as np

from _common import (PAIRS, TEMPORAL, animals_filter, cfg_from_args, config, fs_suffix,
                     temporal_parser)
from tom_cca import dataio, fixed_subspace, lagged, partial, preprocess

K = TEMPORAL.k
N_FOLDS = TEMPORAL.n_folds
LABEL_W = TEMPORAL.label_w_bins   # bins = +/-50 ms — the window the FF/FB label is taken at
N_DIMS = 10                       # dims carried through (FDR family is the same 10)
FDR_DIMS = TEMPORAL.fdr_dims
# NO SAMPLE CAP. The other drivers cap at ~12k bins because they refit CCA 255 times
# per pair (51 lags x 5 folds); this one fits ONCE, so the cap bought nothing and cost a
# great deal. A trial is ~2900 bins (~30 s of running); the previous cap kept 44-184 of
# them, i.e. the first 0.4-1.8 s, which is the trial-onset phase -- and running speed
# rises +12.0 cm/s naive->expert on exactly those bins (11/12 animals, p = 0.0015)
# against +6.6 cm/s over whole trials, so the cap DOUBLED the speed confound on the very
# contrast being measured. Measured cost of the full session (370k bins): residualise
# 4.3 s, PCA 1.1 s, one CCA fit 0.9 s.
EPOCHS = list(config.EPOCH_NAMES)
FIELDS = ["animal", "learner", "pair", "dim", "epoch", "bin_ms", "lag_bins", "lag_ms",
          "r", "label", "ifi_fit", "label_loo", "ifi_loo", "label_xv", "ifi_xv",
          "n_trials_xv", "sig", "p_perdim", "cc_heldout", "n_bins", "n_fit_trials"]


def parse_args():
    return temporal_parser(__doc__.splitlines()[0], cap=False).parse_args()


def main():
    args = parse_args()
    cfg = cfg_from_args(args)
    stem = args.out or f"cc_label_track_bin{args.bin_ms}{fs_suffix(args.include_fs)}"
    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(animals, cfg,
                                     behaviour_lookup=dataio.load_learning_points())
    lags = np.arange(-args.max_lag, args.max_lag + 1)
    # lags are in BINS, so the label window is a bin comparison (LABEL_W = 5 bins =
    # +/-50 ms at 10 ms). information_flow_index only uses the sign of the lag, but
    # keeping the units consistent here matters because the same array indexes r.
    m_lab = np.abs(lags) <= LABEL_W
    print(f"CC LABEL+TRACK | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} | "
          f"label at ±{LABEL_W * args.bin_ms} ms | {args.n_shuffles} shuffles | "
          f"FS={'incl' if args.include_fs else 'excl'}\n")
    only = animals_filter(args.animals)

    out = config.RESULTS_DIR / f"{stem}.csv"
    rows = []

    def _flush():
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    for a in animals:
        if a.animal_id not in entries:              # epochs need a learning point
            continue
        if only is not None and int(a.animal_id) not in only:
            continue
        sess = preprocess.load_running_session(a, cfg, entries)     # every bin, every trial
        if isinstance(sess, str):
            print(f"  animal {a.animal_id}: SKIP — {sess}"); continue
        epoch_of_trial = preprocess.epoch_of_trial(sess.lp, sess.trials_full, cfg, EPOCHS)
        if epoch_of_trial is None:
            print(f"  animal {a.animal_id}: epoch_windows None (skip)"); continue
        trial_ids = sess.trial_ids

        for ax, ay in PAIRS:
            pair = sess.pair(ax, ay)
            if pair is None:
                continue
            X, Y, Z = pair
            # Residualise once, on the fit rows (= all rows here), and fit/project in
            # that same space. Z=None below so fit_fixed does not residualise again.
            if Z is not None:
                allrows = np.ones(X.shape[0], dtype=bool)
                Xr = partial.partial_out_cv(X, Z, allrows)
                Yr = partial.partial_out_cv(Y, Z, allrows)
            else:
                Xr, Yr = X, Y
            fit = fixed_subspace.fit_fixed(Xr, Yr, trial_ids, Z=None, k=K)
            if fit is None or fit.cx is None:
                print(f"  fit fail {a.animal_id} {ax}-{ay}"); continue
            d = int(min(N_DIMS, fit.wx.shape[1], fit.wy.shape[1]))

            # Significance of the FROZEN fit's own canonical correlations, rank for
            # rank. The previous version took them from heldout_lag_curve_flat_perdim,
            # which refits CCA per fold and averages BY RANK, then attached them to the
            # frozen fit's dim by bare index -- the same mis-alignment class as d38a833.
            # Its signature was visible in the output: cc_heldout rose with rank in 38/38
            # cells and went negative in 27/38. frozen_perm_null cannot have that failure
            # because observed and null are the same estimator on the same fit.
            Sx = (Xr - fit.x_mean) @ fit.cx
            Sy = (Yr - fit.y_mean) @ fit.cy
            try:
                res = fixed_subspace.frozen_perm_null(
                    Sx, Sy, n_shuffles=args.n_shuffles, seed=0, correct="fdr",
                    fdr_dims=FDR_DIMS)
                sig_mask, pvals, cc_ho = res.mask, res.p, res.r_obs
            except Exception as e:
                print(f"  sig fail {a.animal_id} {ax}-{ay}: {e}")
                sig_mask = np.zeros(d, dtype=bool)
                pvals, cc_ho = np.ones(d), np.full(d, np.nan)

            for dim in range(d):
                # ONE pairing pass per dim. Every curve below — whole-session label,
                # each epoch, and each leave-that-epoch-out label — is an aggregation of
                # these per-(trial, lag) sufficient statistics, because they are additive
                # over trials. Recomputing the curve per subset instead made the
                # uncapped run intractable (it timed out at 41 min on 3 shuffles).
                u, v = fixed_subspace.project(Xr, Yr, fit, dim=dim)
                trials_m, M = fixed_subspace.trial_lag_moments(u, v, trial_ids, lags)

                r_fit = fixed_subspace.curve_from_moments(M)
                ifi_fit = lagged.information_flow_index(lags[m_lab], r_fit[m_lab])
                label = "FF" if ifi_fit > 0 else ("FB" if ifi_fit < 0 else "none")
                # Cross-validated label (Tom, 2026-08-18): the sign of the IFI on the
                # trials OUTSIDE both plotted epochs (everything but naive and expert),
                # so a naive-vs-expert cross-correlogram split by this label never saw
                # the data it is plotted on. ⚠ The frozen WEIGHTS still come from all
                # trials; only the FF/FB assignment is held out.
                plot_trials = [t for t, e in epoch_of_trial.items()
                               if e in ("naive", "expert")]
                sel_xv = ~np.isin(trials_m, plot_trials)
                r_xv = fixed_subspace.curve_from_moments(M, sel_xv)
                ifi_xv = (lagged.information_flow_index(lags[m_lab], r_xv[m_lab])
                          if np.any(np.isfinite(r_xv)) else float("nan"))
                label_xv = ("FF" if ifi_xv > 0 else ("FB" if ifi_xv < 0 else "none")
                            if np.isfinite(ifi_xv) else "none")
                n_trials_xv = int(sel_xv.sum())

                for epoch in EPOCHS:
                    ep_trials = [t for t, e in epoch_of_trial.items() if e == epoch]
                    sel_t = np.isin(trials_m, ep_trials)
                    if not sel_t.any():
                        continue
                    r_ep = fixed_subspace.curve_from_moments(M, sel_t)
                    if not np.any(np.isfinite(r_ep)):
                        continue
                    # label recomputed WITHOUT this epoch, so persistence has a genuine
                    # 50% chance level (the all-trials label shares data with every
                    # epoch and matches it ~60% of the time by construction)
                    r_out = fixed_subspace.curve_from_moments(M, ~sel_t)
                    ifi_out = lagged.information_flow_index(lags[m_lab], r_out[m_lab])
                    label_out = ("FF" if ifi_out > 0
                                 else ("FB" if ifi_out < 0 else "none"))
                    n_ep = int(M[sel_t, lags.size // 2, 0].sum())
                    for lag, val in zip(lags, r_ep):
                        rows.append({
                            "animal": a.animal_id, "learner": 1,
                            "pair": f"{ax}-{ay}", "dim": dim + 1, "epoch": epoch,
                            "bin_ms": args.bin_ms, "lag_bins": int(lag),
                            "lag_ms": int(lag) * args.bin_ms,
                            "r": round(float(val), 5) if np.isfinite(val) else "",
                            "label": label,
                            "ifi_fit": round(float(ifi_fit), 5),
                            "label_loo": label_out,
                            "ifi_loo": round(float(ifi_out), 5)
                                       if np.isfinite(ifi_out) else "",
                            "label_xv": label_xv,
                            "ifi_xv": round(float(ifi_xv), 5)
                                      if np.isfinite(ifi_xv) else "",
                            "n_trials_xv": n_trials_xv,
                            "sig": int(sig_mask[dim]) if dim < sig_mask.size else 0,
                            "p_perdim": round(float(pvals[dim]), 5)
                                        if dim < pvals.size else "",
                            "cc_heldout": round(float(cc_ho[dim]), 5)
                                          if dim < cc_ho.size
                                          and np.isfinite(cc_ho[dim]) else "",
                            "n_bins": n_ep,
                            "n_fit_trials": fit.n_trials,
                        })
        _flush()
        print(f"  animal {a.animal_id}: {len(rows)} rows total")
    _flush()
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
