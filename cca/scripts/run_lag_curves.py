"""Held-out PER-DIMENSION lagged canonical-correlation curves (the cross-correlation
profile behind the IFI), for the lag-CC figure (report R2).

For each animal x area-pair this mirrors ``run_ifi_windows`` exactly — residualise the
third area, PCA->K on all running bins, then the **held-out, segment-aware** lag curve
(CCA refit per fold at each lag; the lag never crosses a trial boundary) — but exports
the full per-dimension curve (``lagged.heldout_lag_curve_flat_perdim``) rather than the
IFI summary. Significance of each canonical dim is computed on THESE scores from THIS curve's own
lag-0 held-out CC (``lagged.perdim_significance``, the same circular-shift dominant-dim
surrogate the rest of the report uses), so the flag and the curve describe the same
canonical dimension. The per-dim curves can then be pooled two ways in the figure:

  * n = ANIMALS    -- the dominant-dim (CC1) curve, one per animal, averaged across animals
  * n = SUBSPACES  -- every SIGNIFICANT canonical dim's curve, pooled across animals & dims

Writes one row per (animal, pair, lag, dim):
    results/lag_curves_bin{10,25,50}{,_fsincl}.csv

Usage: python scripts/run_lag_curves.py [--include-fs] [--animals 36 --max-samples 3000]
       (defaults = config.TEMPORAL: 10 ms, sigma 2.5 ms, K 30, +/-25 lags, 200 shuffles, uncapped)
Significance uses a PER-DIM null computed HELD-OUT: dimension j's observed held-out CC
is compared to the shuffled distribution of dimension j, evaluated through the same
fold structure. The earlier dominant-dim null compared a held-out observed value to an
IN-SAMPLE shuffled one and gated every dim on the top dim's null -- doubly conservative,
leaving 0.7 significant dims per cell.

Expensive: the null costs n_shuffles x n_folds CCA fits per pair on top of the 51-lag
curve (~70 min per FS condition at 200 shuffles). Launch detached.
"""

from __future__ import annotations

import csv

import numpy as np

from _common import (PAIRS, TEMPORAL, animals_filter, cap_desc, cfg_from_args, config,
                     fs_suffix, temporal_parser)
from tom_cca import dataio, lagged, preprocess

K = TEMPORAL.k
N_FOLDS = TEMPORAL.n_folds
FDR_DIMS = TEMPORAL.fdr_dims       # BH family = leading dims; see lagged.perdim_significance
FIELDS = ["animal", "learner", "pair", "bin_ms", "lag_bins", "lag_ms",
          "dim", "cc", "sig", "sig_uncorr", "p_perdim", "n_sig"]


def parse_args():
    return temporal_parser(__doc__.splitlines()[0]).parse_args()


def main():
    args = parse_args()
    cfg = cfg_from_args(args)
    stem = args.out or f"lag_curves_bin{args.bin_ms}{fs_suffix(args.include_fs)}"
    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(animals, cfg,
                                     behaviour_lookup=dataio.load_learning_points())
    print(f"LAG-CC curves | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} "
          f"max_lag={args.max_lag} ({args.max_lag * args.bin_ms} ms) | "
          f"FS={'incl' if args.include_fs else 'excl'} | cap: {cap_desc(args)}\n")
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
        tr = sess.trials
        print(f"  animal {a.animal_id}: {sess.n_bins} of {sess.n_running_total} running bins, "
              f"trials {int(tr.min())}..{int(tr.max())} ({tr.size} trials)")
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
            d = int(min(Sx.shape[1], Sy.shape[1]))
            try:
                lags, cc = lagged.heldout_lag_curve_flat_perdim(
                    Sx, Sy, trial_ids, max_lag=args.max_lag, n_dims=d, n_folds=N_FOLDS)
            except Exception as e:
                print(f"  curve fail {a.animal_id} {ax}-{ay}: {e}"); continue
            # Significance of each canonical dim, computed on THESE SCORES and THIS
            # curve's own lag-0 held-out CC. It used to come from a separate
            # window_subspace fit (own residualisation, own PCA, own folds) and was
            # attached by bare index, so dim k of that fit was matched to dim k of this
            # one: 19% of flagged dims had a NEGATIVE held-out CC and the largest CC in
            # the dataset was flagged non-significant. Same circular-shift
            # dominant-dim null as before — only the alignment is fixed.
            zero = int(np.flatnonzero(lags == 0)[0])
            try:
                res = lagged.perdim_significance(
                    Sx, Sy, cc[zero, :], groups=trial_ids,
                    n_shuffles=args.n_shuffles, n_folds=N_FOLDS, seed=0,
                    null_mode="perdim", correct="fdr", fdr_dims=FDR_DIMS)
                sig_mask, pvals = res.mask, res.p
                n_sig = int(sig_mask.sum())
            except Exception as e:
                print(f"  sig fail {a.animal_id} {ax}-{ay}: {e}")
                sig_mask = np.zeros(d, dtype=bool)
                pvals, n_sig = np.ones(d), 0
            for li, lag in enumerate(lags):
                for dim in range(d):
                    c = cc[li, dim]
                    rows.append({
                        "animal": a.animal_id, "learner": int(is_learner),
                        "pair": f"{ax}-{ay}", "bin_ms": args.bin_ms,
                        "lag_bins": int(lag), "lag_ms": int(lag) * args.bin_ms,
                        "dim": dim + 1,
                        "cc": round(float(c), 4) if np.isfinite(c) else "",
                        "sig": int(sig_mask[dim]) if dim < sig_mask.size else 0,
                        "sig_uncorr": int(pvals[dim] < 0.05)
                                      if dim < pvals.size else 0,
                        "p_perdim": round(float(pvals[dim]), 5)
                                    if dim < pvals.size else "",
                        "n_sig": n_sig})
        _flush()
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): done "
              f"({len(rows)} rows)")
    _flush()
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
