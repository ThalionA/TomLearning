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

Usage: PYTHONPATH=src python scripts/run_lag_curves.py --bin-ms 10 --max-lag 25 \
       --smooth-ms 2.5 [--include-fs]
Significance uses a PER-DIM null computed HELD-OUT: dimension j's observed held-out CC
is compared to the shuffled distribution of dimension j, evaluated through the same
fold structure. The earlier dominant-dim null compared a held-out observed value to an
IN-SAMPLE shuffled one and gated every dim on the top dim's null -- doubly conservative,
leaving 0.7 significant dims per cell.

Expensive: the null costs n_shuffles x n_folds CCA fits per pair on top of the 51-lag
curve (~70 min per FS condition at 200 shuffles). Launch detached.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, dataio, lagged, partial  # noqa: E402

K = 30
N_FOLDS = 5
MAX_SAMPLES = 12000          # cap per pair — these sessions run to ~370k engaged bins
MAX_BINS_PER_TRIAL = 600     # <= 6 s/trial kept; consecutive, so lag adjacency survives
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]
FDR_DIMS = 10                # BH family = leading dims; see lagged.perdim_significance
FIELDS = ["animal", "learner", "pair", "bin_ms", "lag_bins", "lag_ms",
          "dim", "cc", "sig", "sig_uncorr", "p_perdim", "n_sig"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=10)
    p.add_argument("--include-fs", action="store_true")
    p.add_argument("--max-lag", type=int, default=25,
                   help="widest lag (bins); curve spans -max_lag..+max_lag")
    p.add_argument("--smooth-ms", type=float, default=0.0,
                   help="Gaussian s.d. (ms) for spike-train smoothing (Buzsáki = 2.5)")
    p.add_argument("--n-shuffles", type=int, default=200,
                   help="circular-shift surrogates for the per-dim null. The empirical "
                        "p floor is 1/(n+1); BH over FDR_DIMS needs it below "
                        "alpha/FDR_DIMS, i.e. >199 for 10 dims at alpha=0.05.")
    p.add_argument("--out", default="", help="output CSV stem override")
    return p.parse_args()


def _pca_scores(M, k):
    mu = M.mean(0)
    _, _, vt = np.linalg.svd(M - mu, full_matrices=False)
    return (M - mu) @ vt[:k].T


def _capped_index(trial_ids, max_samples=MAX_SAMPLES, max_per_trial=MAX_BINS_PER_TRIAL):
    """Indices keeping a CONSECUTIVE within-trial block (<= max_per_trial bins) from
    each whole trial in turn until ~max_samples — bounds cost on these long sessions
    while preserving the bin adjacency the segment-aware lag pairing needs (an even
    stride would break it)."""
    keep, total = [], 0
    for t in np.unique(trial_ids):
        pos = np.where(trial_ids == t)[0][:max_per_trial]   # temporal order within trial
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
    stem = args.out or f"lag_curves_bin{args.bin_ms}{suffix}"
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    print(f"LAG-CC curves | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} "
          f"max_lag={args.max_lag} ({args.max_lag * args.bin_ms} ms) | "
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
        cap = _capped_index(trial_ids_full)
        trial_ids = trial_ids_full[cap]
        is_learner = a.animal_id in entries
        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run][cap]

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
            Sx = _pca_scores(Xr, min(K, Xr.shape[1]))
            Sy = _pca_scores(Yr, min(K, Yr.shape[1]))
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
