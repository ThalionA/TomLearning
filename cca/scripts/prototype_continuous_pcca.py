"""Prototype: pCCA on CONTINUOUS running data at the paper-correct sample regime.

Demonstrates the fix for the overfitting that plagued the epoch/landmark sweep
(`STATE.md` §4: high-k configs saturate held-out CC -> 0.999). Both reference
papers fit (p)CCA on a long continuous span with ~>=50 samples per variable,
PCA-reduced to a fixed modest number of components. Here we:

* load each learner's whole-session activity at 25 ms (reusing the temporal
  streams loader -- no new I/O, no depth/ISI/waveforms),
* keep running, in-trial bins only,
* residualise CA1 and the partner area on ALL other recorded areas (pCCA),
* PCA-reduce each to K components, fit CCA with whole-trial 5-fold CV,
* report n_samples, samples/component, and the held-out dominant CC.

Success = tens of thousands of samples, samples/K >> 50, and a STABLE,
non-saturated held-out CC -- unlike the ~2,000-sample epoch regime.

Read-only, exploratory. Not wired into the sweep.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, core, dataio, partial  # noqa: E402

K = 30                       # fixed PCA components per area (paper regime)
N_FOLDS = 5
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG")]


def _pca_project(train, test, k):
    mu = train.mean(0)
    _, _, vt = np.linalg.svd(train - mu, full_matrices=False)
    comp = vt[:k].T
    return (train - mu) @ comp, (test - mu) @ comp


def _heldout_cc1(X, Y, trial_ids, k, n_folds):
    """Whole-trial CV dominant held-out canonical correlation."""
    trials = np.unique(trial_ids)
    if trials.size < n_folds:
        return np.nan, np.nan
    rng = np.random.default_rng(0)
    folds = np.array_split(rng.permutation(trials), n_folds)
    cc1_test, cc1_train = [], []
    for test_tr in folds:
        te = np.isin(trial_ids, test_tr)
        tr = ~te
        kk = int(min(k, X.shape[1], Y.shape[1], np.sum(tr) - 1))
        if kk < 2:
            continue
        sx_tr, sx_te = _pca_project(X[tr], X[te], kk)
        sy_tr, sy_te = _pca_project(Y[tr], Y[te], kk)
        model = core.cca_fit(sx_tr, sy_tr)
        r_te = core.cca_score(sx_te, sy_te, model)
        r_tr = core.cca_score(sx_tr, sy_tr, model)
        if r_te.size:
            cc1_test.append(float(r_te[0]))
            cc1_train.append(float(r_tr[0]))
    if not cc1_test:
        return np.nan, np.nan
    return float(np.nanmean(cc1_test)), float(np.nanmean(cc1_train))


def main():
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=25)
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    learners = [a for a in animals if a.animal_id in entries]
    thr = cfg.velocity_thresh_cm_s

    print(f"Regime: bin={cfg.temporal_bin_ms} ms | PCA K={K} | pCCA (control all "
          f"other recorded areas) | whole-trial {N_FOLDS}-fold CV")
    print(f"running threshold = {thr} cm/s\n")
    print(f"{'animal':>6} {'pair':<9} {'n_run_bins':>10} {'samp/K':>7} "
          f"{'cc1_test':>9} {'cc1_train':>9} {'gap':>6}")
    print("-" * 62)

    rows = []
    for a in learners:
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids = streams.trial_idx_50ms[run]
        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run]                       # (n_run, n_units)
        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            n = X.shape[0]
            if others:
                Z = np.concatenate(others, axis=1)
                Xr = partial.partial_out(X, Z)
                Yr = partial.partial_out(Y, Z)
            else:
                Xr, Yr = X, Y
            cc_te, cc_tr = _heldout_cc1(Xr, Yr, trial_ids, K, N_FOLDS)
            gap = (cc_tr - cc_te) if np.isfinite(cc_tr) else np.nan
            print(f"{a.animal_id:>6} {ax+'-'+ay:<9} {n:>10} {n/K:>7.0f} "
                  f"{cc_te:>9.3f} {cc_tr:>9.3f} {gap:>6.3f}")
            rows.append((ax + "-" + ay, n, cc_te, cc_tr))

    if rows:
        ns = [r[1] for r in rows]
        ccs = [r[2] for r in rows if np.isfinite(r[2])]
        sat = sum(1 for c in ccs if c >= 0.99)
        print("\nSummary:")
        print(f"  fits: {len(rows)} | median n_run_bins: {int(np.median(ns))} "
              f"(samp/K ~ {np.median(ns)/K:.0f})")
        print(f"  held-out cc1: median {np.median(ccs):.3f}, "
              f"range [{min(ccs):.3f}, {max(ccs):.3f}]")
        print(f"  saturated (cc1>=0.99): {sat}/{len(ccs)}  "
              "<- should be ~0, vs the epoch regime where high-k saturated")


if __name__ == "__main__":
    main()
