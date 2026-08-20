"""Trial 1 vs trial 2 through a FROZEN reference subspace — per-trial lag curves.

Does communication change between the very first and second running trial of the
session? A single trial (~500-2900 running bins) cannot support a 30-dim CCA refit
(see `early_trials.py`), so the subspace is fit ONCE per (animal, pair) on trial
ordinals > N_HELD (the sample-rich later trials, `early_trials.reference_fit`) and
each of the first N_HELD=10 running trials is PROJECTED through it, leak-free —
ordinals 3..10 are equally out-of-fit and form the adjacent-trial control band.

This driver exports per-trial lagged correlations of the frozen canonical variates
(`fixed_subspace.variate_lag_curve`), one row per (animal, pair, dim, ordinal,
matched, lag). Every IFI (all integration windows, degenerate flags) is derived
downstream by `analyze_trial12.py` through `perdim_ifi.windows_table`, so the
degenerate-zero convention cannot leak into means here.

MEASURED CONFOUND this design must carry: trial 1 has systematically MORE running
bins than later trials (vs trial 4: median +815 bins ~ +8 s of running, paired-t
p = 0.003, 12/15 animals — `early_trials_projected_bin10.csv`). r0 dispersion and
IFI clipping are both n-dependent, so every metric is exported three ways
(`matched` column):

  0  raw — all of the trial's running bins
  1  pairwise-matched — ordinals 1 and 2 only, cut to min(n1, n2) LEADING bins
  2  common-matched — every ordinal cut to the min bin count over the gated
     ordinals 1..N_HELD (the arm the adjacent-trial control is standardised on)

Leading CONSECUTIVE blocks, never strided (segment-aware lag pairing needs bin
adjacency — GOTCHAS.md). "Ordinal" means the o-th RUNNING trial: a trial with no
supra-threshold running bin never enters `trial_ids`, so ordinal 1 is the first
trial the animal actually ran in; `trial_id` is exported so the mapping stays
auditable. Per-dim significance = `fixed_subspace.frozen_perm_null` on the TRAIN
rows only (the fit's own canonical correlations — readout trials stay out; the
all-rows variant in `run_cc_label_track.py` has no held-out set to protect);
`r_frozen` = those in-sample canonical correlations, while the per-trial `r` IS
held out. Both FS conditions: run once per `--include-fs` setting.

Writes results/trial12_curves_bin10{,_fsincl}.csv   (large, regenerable — gitignored)
       results/trial12_trials_bin10{,_fsincl}.csv   (per-trial behaviour: n_bins,
                                                     n_gaps, speed — dim-invariant
                                                     columns repeated per dim)

Usage: python scripts/run_trial12.py [--include-fs] [--animals 36]   (defaults = config.TEMPORAL)
"""

from __future__ import annotations

import csv

import numpy as np

from _common import (PAIRS, TEMPORAL, animals_filter, cfg_from_args, results_path,
                     temporal_parser)
from tom_cca import config, dataio, early_trials, fixed_subspace, preprocess

K = TEMPORAL.k
N_DIMS = 10                       # dims carried through (FDR family is the same 10)
FDR_DIMS = TEMPORAL.fdr_dims
N_HELD = 10                       # ordinals 1..N_HELD read out; fit on the rest
MIN_TRIALS = 2 * N_HELD           # animal gate: >= 10 held-out + >= 10 fit trials
CURVE_FIELDS = ["animal", "learner", "pair", "dim", "ordinal", "trial_id", "matched",
                "bin_ms", "lag_bins", "lag_ms", "r", "sig", "p_perdim", "r_frozen",
                "n_bins", "n_fit_trials"]
TRIAL_FIELDS = ["animal", "learner", "pair", "dim", "ordinal", "trial_id", "matched",
                "r0", "sig", "p_perdim", "r_frozen", "n_bins", "n_gaps",
                "vel_mean", "vel_median", "n_fit_trials"]


def parse_args():
    p = temporal_parser(__doc__.splitlines()[0], cap=False)
    return p.parse_args()


def pair_rows(u, v, trial_ids, uniq, lags, d, sig_mask, pvals, r_frozen,
              velocity, stream_index, min_bins, n_held=N_HELD, bin_ms=TEMPORAL.bin_ms):
    """Curve + trial rows for ordinals 1..n_held of ONE (animal, pair).

    ``u``/``v`` are the frozen canonical variates on ALL rows, (n_bins, >= d).
    Returns ``(curve_rows, trial_rows)`` without the animal/learner/pair keys (the
    caller adds them). A trial below ``min_bins`` (after the matched cut) writes
    trial rows with a blank ``r0`` and no curve rows — gated, never silently
    dropped.
    """
    bins_by_o = {o: np.flatnonzero(trial_ids == uniq[o - 1])
                 for o in range(1, min(n_held, uniq.size) + 1)}
    n_full = {o: b.size for o, b in bins_by_o.items()}
    passing = [o for o, nb in n_full.items() if nb >= min_bins]
    n12 = min(n_full.get(1, 0), n_full.get(2, 0))
    ncom = min(n_full[o] for o in passing) if passing else 0
    arms = {0: None, 1: n12, 2: ncom}
    lags = np.asarray(lags)
    i0 = int(np.flatnonzero(lags == 0)[0])
    curve_rows, trial_rows = [], []
    for o, bins_full in sorted(bins_by_o.items()):
        for matched, ncap in arms.items():
            if matched == 1 and o not in (1, 2):
                continue
            bins = bins_full if ncap is None else bins_full[:ncap]
            gated = bins.size < min_bins
            n_gaps = (int(np.sum(np.diff(stream_index[bins]) > 1))
                      if stream_index is not None and bins.size else 0)
            vel = velocity[bins] if velocity is not None and bins.size else None
            for dim in range(d):
                meta = {
                    "dim": dim + 1, "ordinal": o, "trial_id": int(uniq[o - 1]),
                    "matched": matched,
                    "sig": int(sig_mask[dim]) if dim < sig_mask.size else 0,
                    "p_perdim": (round(float(pvals[dim]), 5)
                                 if dim < pvals.size else ""),
                    "r_frozen": (round(float(r_frozen[dim]), 5)
                                 if dim < r_frozen.size and np.isfinite(r_frozen[dim])
                                 else ""),
                    "n_bins": int(bins.size),
                }
                behaviour = {
                    "n_gaps": n_gaps,
                    "vel_mean": round(float(np.mean(vel)), 3) if vel is not None else "",
                    "vel_median": round(float(np.median(vel)), 3) if vel is not None else "",
                }
                if gated:
                    trial_rows.append({**meta, **behaviour, "r0": ""})
                    continue
                r = fixed_subspace.variate_lag_curve(u[bins, dim], v[bins, dim],
                                                     trial_ids[bins], lags)
                trial_rows.append({
                    **meta, **behaviour,
                    "r0": round(float(r[i0]), 5) if np.isfinite(r[i0]) else ""})
                for lag, val in zip(lags, r):
                    curve_rows.append({
                        **meta, "bin_ms": bin_ms, "lag_bins": int(lag),
                        "lag_ms": int(lag) * bin_ms,
                        "r": round(float(val), 5) if np.isfinite(val) else ""})
    return curve_rows, trial_rows


def main():
    args = parse_args()
    cfg = cfg_from_args(args)
    min_bins = max(3 * args.max_lag, 200)
    lags = np.arange(-args.max_lag, args.max_lag + 1)
    curves_out = results_path(args.out or "trial12_curves", fs=args.include_fs,
                              bin_ms=args.bin_ms)
    trials_out = results_path("trial12_trials", fs=args.include_fs, bin_ms=args.bin_ms)
    print(f"TRIAL 1 vs 2 | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} | "
          f"fit on ordinals >{N_HELD}, read ordinals 1..{N_HELD} | "
          f"gates: >={MIN_TRIALS} running trials, >={min_bins} bins/trial | "
          f"{args.n_shuffles} shuffles | FS={'incl' if args.include_fs else 'excl'}\n")

    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(animals, cfg,
                                     behaviour_lookup=dataio.load_learning_points())
    only = animals_filter(args.animals)
    curve_rows, trial_rows = [], []

    def _flush():
        for path, fields, rows in ((curves_out, CURVE_FIELDS, curve_rows),
                                   (trials_out, TRIAL_FIELDS, trial_rows)):
            with open(path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
                w.writeheader(); w.writerows(rows)

    n_pairs_done = 0
    for a in animals:
        if only is not None and int(a.animal_id) not in only:
            continue
        sess = preprocess.load_running_session(a, cfg, entries, min_trials=MIN_TRIALS)
        if isinstance(sess, str):
            print(f"  animal {a.animal_id}: SKIP — {sess}"); continue
        uniq = sess.trials
        train_ids = uniq[N_HELD:]
        train_mask = np.isin(sess.trial_ids, train_ids)

        for ax, ay in PAIRS:
            pair = sess.pair(ax, ay)
            if pair is None:
                continue
            X, Y, Z = pair
            try:
                fit = early_trials.reference_fit(X, Y, train_mask, Z=Z, k=K)
            except Exception as e:                              # noqa: BLE001
                print(f"  fit fail {a.animal_id} {ax}-{ay}: {e}"); continue
            d = int(min(N_DIMS, fit.model.A.shape[1], fit.model.B.shape[1]))
            try:
                res = fixed_subspace.frozen_perm_null(
                    fit.Sx[train_mask], fit.Sy[train_mask],
                    n_shuffles=args.n_shuffles, seed=0, correct="fdr",
                    fdr_dims=FDR_DIMS)
                sig_mask, pvals, r_frozen = res.mask, res.p, res.r_obs
                if not np.allclose(r_frozen[:d], fit.model.r[:d], atol=1e-6):
                    print(f"  ⚠ {a.animal_id} {ax}-{ay}: frozen_perm_null r_obs != "
                          "model.r — the null is not testing the fit it claims to")
            except Exception as e:                              # noqa: BLE001
                print(f"  sig fail {a.animal_id} {ax}-{ay}: {e}")
                sig_mask = np.zeros(d, dtype=bool)
                pvals, r_frozen = np.ones(d), np.full(d, np.nan)

            u_all, v_all = early_trials.variates(fit, slice(None))
            crows, trows = pair_rows(
                u_all, v_all, sess.trial_ids, uniq, lags, d, sig_mask, pvals,
                r_frozen, sess.velocity, sess.stream_index, min_bins,
                bin_ms=args.bin_ms)
            base = {"animal": a.animal_id, "learner": int(sess.learner),
                    "pair": f"{ax}-{ay}", "n_fit_trials": int(train_ids.size)}
            curve_rows += [{**base, **r} for r in crows]
            trial_rows += [{**base, **r} for r in trows]
            n_pairs_done += 1
        _flush()
        print(f"  animal {a.animal_id}: {uniq.size} running trials, "
              f"{len(trial_rows)} trial rows total")
    _flush()
    print(f"\n{n_pairs_done} animal-pairs -> {curves_out}\n{' ' * 16}-> {trials_out}")


if __name__ == "__main__":
    main()
