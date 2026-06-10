"""Prototype: regularised KERNEL CCA on the continuous running data.

The Step-1 superiority/sanity test for going nonlinear (see the methods note). For
each learner and area pair, on the SAME PCA->30 inputs the linear pipeline uses,
fit regularised RBF-KCCA with whole-trial train/test and ask the two questions
that gate any further investment:

  (1) SATURATION — does the held-out KCCA correlation blow up / does in-sample
      saturate at OUR per-window sample count? Reported as the in-sample vs
      held-out gap across a ridge grid (a large gap at small reg = overfitting).
  (2) SURROGATE — does the real held-out KCCA correlation clear a circular-shift
      null (roll one area, preserving its autocorrelation, destroying cross-area
      timing; refit; 95th percentile of the held-out null)?

Plus the head-to-head: KCCA held-out CC vs LINEAR CCA held-out CC on identical
folds (does nonlinearity buy out-of-sample predictive power?).

KCCA is O(n^3), so each fit is sub-sampled to N_CAP running bins. Read-only,
exploratory; not wired into the sweep.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, core, dataio, kernel_cca, partial  # noqa: E402

K = 30                       # PCA components per area (same as the linear pipeline)
N_CAP = 1500                 # sub-sample cap per fit (KCCA is O(n^3))
TEST_FRAC = 0.3              # whole-trial held-out fraction
REGS = (0.1, 1.0, 10.0)      # ridge grid for the saturation sweep
REG_MAIN = 1.0               # ridge used for the surrogate + head-to-head
N_SHUF = 20                  # circular-shift surrogate draws
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA3", "DG"),
         ("CA1", "V1")]


def _pca_project(train, test, k):
    mu = train.mean(0)
    _, _, vt = np.linalg.svd(train - mu, full_matrices=False)
    comp = vt[:k].T
    return (train - mu) @ comp, (test - mu) @ comp


def _trial_split(trial_ids, frac, seed):
    trials = np.unique(trial_ids)
    if trials.size < 4:
        return None, None
    rng = np.random.default_rng(seed)
    perm = rng.permutation(trials)
    n_test = max(1, int(round(frac * trials.size)))
    test_tr = perm[:n_test]
    te = np.isin(trial_ids, test_tr)
    return ~te, te


def _subsample(n, cap, seed=0):
    if n <= cap:
        return np.ones(n, dtype=bool)
    idx = np.sort(np.random.default_rng(seed).choice(n, cap, replace=False))
    m = np.zeros(n, dtype=bool); m[idx] = True
    return m


def _heldout(Xr, Yr, tr, te, reg):
    """(linear cc1_test, kcca cc1_train(in-sample), kcca cc1_test) on one split."""
    kk = int(min(K, Xr.shape[1], Yr.shape[1], np.sum(tr) - 1))
    sx_tr, sx_te = _pca_project(Xr[tr], Xr[te], kk)
    sy_tr, sy_te = _pca_project(Yr[tr], Yr[te], kk)
    lin = core.cca_fit(sx_tr, sy_tr)
    lin_te = float(core.cca_score(sx_te, sy_te, lin)[0])
    m = kernel_cca.kcca_fit(sx_tr, sy_tr, reg=reg, n_components=1)
    k_tr = float(m.r[0])
    k_te = float(kernel_cca.kcca_score(sx_te, sy_te, m)[0])
    return lin_te, k_tr, k_te


def _surrogate(Xr, Yr, trial_ids, reg, seed):
    """95th-pct held-out KCCA CC under circular shift of Y (autocorr-preserving)."""
    n = Xr.shape[0]
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(N_SHUF):
        shift = int(rng.integers(1, n))
        Ys = np.roll(Yr, shift, axis=0)
        tr, te = _trial_split(trial_ids, TEST_FRAC, seed=int(rng.integers(1e6)))
        if tr is None:
            continue
        try:
            null.append(_heldout(Xr, Ys, tr, te, reg)[2])
        except Exception:
            continue
    return float(np.nanquantile(null, 0.95)) if null else np.nan


def main():
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=25)
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    learners = [a for a in animals if a.animal_id in entries]
    thr = cfg.velocity_thresh_cm_s

    print(f"KCCA prototype | bin=25ms PCA K={K} | cap={N_CAP} bins/fit | "
          f"reg grid {REGS} | surrogate {N_SHUF} shuffles @ reg={REG_MAIN}")
    print(f"saturation = in-sample - held-out gap (large = overfit); "
          f"sig = held-out > 95th-pct shift null\n")
    hdr = (f"{'animal':>6} {'pair':<9} {'n':>5} | {'lin_te':>7} | "
           + " ".join(f"k_te@{r}".rjust(9) for r in REGS) + " | "
           + f"{'k_tr@1':>7} {'gap@1':>6} | {'null95':>7} {'sig?':>5}")
    print(hdr); print("-" * len(hdr))

    rows = []
    for a in learners:
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids_all = streams.trial_idx_50ms[run]
        present = {}
        for area in config.AREAS:
            mm, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = mm[run]
        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            Xr = partial.partial_out(X, Z) if Z is not None else X
            Yr = partial.partial_out(Y, Z) if Z is not None else Y
            sub = _subsample(Xr.shape[0], N_CAP, seed=a.animal_id)
            Xr, Yr, tids = Xr[sub], Yr[sub], trial_ids_all[sub]
            tr, te = _trial_split(tids, TEST_FRAC, seed=a.animal_id)
            if tr is None:
                continue
            try:
                k_te = {}
                for reg in REGS:
                    lin_te, k_tr, kt = _heldout(Xr, Yr, tr, te, reg)
                    k_te[reg] = kt
                    if reg == REG_MAIN:
                        gap = k_tr - kt; ktr_main = k_tr
                null95 = _surrogate(Xr, Yr, tids, REG_MAIN, seed=a.animal_id + 7)
            except Exception as exc:                 # noqa: BLE001
                print(f"{a.animal_id:>6} {ax+'-'+ay:<9}  fit failed: {type(exc).__name__}")
                continue
            sig = (np.isfinite(null95) and k_te[REG_MAIN] > null95)
            print(f"{a.animal_id:>6} {ax+'-'+ay:<9} {Xr.shape[0]:>5} | "
                  f"{lin_te:>7.3f} | " + " ".join(f"{k_te[r]:>9.3f}" for r in REGS)
                  + f" | {ktr_main:>7.3f} {gap:>6.3f} | {null95:>7.3f} "
                  f"{'YES' if sig else 'no':>5}", flush=True)
            rows.append((ax + "-" + ay, lin_te, k_te[REG_MAIN], ktr_main, gap,
                         null95, sig))

    if rows:
        lin = np.array([r[1] for r in rows]); kcca = np.array([r[2] for r in rows])
        gaps = np.array([r[4] for r in rows]); sigs = sum(r[6] for r in rows)
        sat = sum(1 for r in rows if r[2] >= 0.99)
        print("\nSUMMARY")
        print(f"  fits: {len(rows)} | KCCA held-out cc1>=0.99 (hard saturation): "
              f"{sat}/{len(rows)}")
        print(f"  in-sample vs held-out gap @reg={REG_MAIN}: "
              f"median {np.nanmedian(gaps):.3f} (large => overfit at this n)")
        print(f"  surrogate: significant in {sigs}/{len(rows)} fits")
        print(f"  KCCA vs LINEAR held-out cc1: median KCCA {np.nanmedian(kcca):.3f} "
              f"vs linear {np.nanmedian(lin):.3f} "
              f"(Δ={np.nanmedian(kcca - lin):+.3f}); "
              f"KCCA>linear in {int(np.sum(kcca > lin))}/{len(rows)}")


if __name__ == "__main__":
    main()
