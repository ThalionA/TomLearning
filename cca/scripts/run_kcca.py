"""Driver: NONLINEAR (kernel) CCA full suite, per (animal, pair, scope).

The kernel analogue of run_epochs.py. For EVERY animal (learners and
non-learners) and area pair it fits the full `kcca_window` suite on:
  * scope="session" — the whole running session (the all-animals structural /
    superiority readout; the only scope available for non-learners), and
  * scope in {naive, intermediate, expert} — the learning epochs (learners only).
Each cell reports held-out KCCA CC, the held-out LINEAR CC on the same folds
(superiority), significance (shift surrogate), IFI / optimal lag (direction), and
structure-coefficient Gini (membership). FS-excluded by default; --include-fs for
the FS-included variant. KCCA is O(n^3) so each cell is sub-sampled to CAP bins.

Writes results/kcca_metrics.csv (or _fsincl). Use --animals a,b,c to run a chunk
(for cross-process parallelism); chunks share the same CSV name per FS variant, so
run chunks to distinct --out names and merge, or run the whole cohort in one go.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, dataio, kcca_window  # noqa: E402

K = 20
REG = 10.0                    # fixed ridge (prototype sweet spot; fair vs linear)
CAP = 900                     # sub-sample cap per cell (KCCA is O(n^3))
MAX_LAG = 8                   # upgraded 4->8 (review 2026-06-11): wider lag window
N_SHUFFLES = 30              # upgraded 10->30 (review 2026-06-11): stabler surrogate
N_FOLDS = 5
D = 3
EPOCHS = ["naive", "intermediate", "expert"]
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]
FIELDS = ["animal", "learner", "pair", "scope", "n_bins", "cc_kcca", "cc_linear",
          "n_sig", "mi_sig", "ifi", "optimal_lag", "gini_x", "gini_y"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=25)
    p.add_argument("--include-fs", action="store_true")
    p.add_argument("--animals", default="", help="comma-separated subset (default all)")
    p.add_argument("--out", default="", help="override output CSV stem")
    return p.parse_args()


def _subsample(n, cap, seed):
    if n <= cap:
        return np.ones(n, dtype=bool)
    idx = np.sort(np.random.default_rng(seed).choice(n, cap, replace=False))
    m = np.zeros(n, dtype=bool); m[idx] = True
    return m


def _cell(X, Y, Z, groups, seed):
    r = kcca_window.kcca_window(X, Y, groups, Z=Z, reg=REG, k=K, max_lag=MAX_LAG,
                                n_shuffles=N_SHUFFLES, n_folds=N_FOLDS, d=D, seed=seed)
    cck = float(r.cc_kcca[0]) if r.cc_kcca.size else float("nan")
    ccl = float(r.cc_linear[0]) if r.cc_linear.size else float("nan")
    return r, cck, ccl


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs)
    suffix = "_fsincl" if args.include_fs else ""
    stem = args.out or f"kcca_metrics{suffix}"
    only = {int(x) for x in args.animals.split(",") if x.strip()} if args.animals else None
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    print(f"KCCA suite | bin={args.bin_ms}ms K={K} reg={REG} cap={CAP} | "
          f"FS {'incl' if args.include_fs else 'excl'} | lags+-{MAX_LAG} "
          f"shuffles={N_SHUFFLES}\n", flush=True)

    out = config.RESULTS_DIR / f"{stem}.csv"
    rows = []

    def _write():
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    for a in animals:
        if only is not None and a.animal_id not in only:
            continue
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids = streams.trial_idx_50ms[run]
        uniq = np.unique(trial_ids)
        is_learner = a.animal_id in entries
        ew = (dataio.epoch_windows(int(entries[a.animal_id].lp), uniq.size, cfg)
              if is_learner else None)
        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run]
        n_rows = 0
        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            # scope rows: whole session (all animals) + epochs (learners)
            scopes = [("session", np.ones(trial_ids.size, dtype=bool))]
            if ew is not None:
                for ep in EPOCHS:
                    pos = ew[ep]; pos = pos[(pos >= 0) & (pos < uniq.size)]
                    scopes.append((ep, np.isin(trial_ids, uniq[pos])))
            for scope, mask in scopes:
                if mask.sum() < K * 30:
                    continue
                sub = _subsample(int(mask.sum()), CAP, seed=a.animal_id)
                idx = np.where(mask)[0][sub]
                Xi, Yi = X[idx], Y[idx]
                Zi = Z[idx] if Z is not None else None
                try:
                    r, cck, ccl = _cell(Xi, Yi, Zi, trial_ids[idx],
                                        seed=a.animal_id)
                except Exception as exc:             # noqa: BLE001
                    print(f"  {a.animal_id} {ax}-{ay} {scope}: "
                          f"fail {type(exc).__name__}", flush=True)
                    continue
                rows.append({
                    "animal": a.animal_id, "learner": int(is_learner),
                    "pair": f"{ax}-{ay}", "scope": scope, "n_bins": int(idx.size),
                    "cc_kcca": round(cck, 4), "cc_linear": round(ccl, 4),
                    "n_sig": r.n_sig, "mi_sig": round(r.mi_sig, 4),
                    "ifi": round(r.ifi, 4), "optimal_lag": r.optimal_lag,
                    "gini_x": round(r.gini_x, 4), "gini_y": round(r.gini_y, 4)})
                n_rows += 1
        _write()
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): "
              f"{n_rows} rows ({len(rows)} total)", flush=True)

    _write()
    print(f"\nwrote {out} ({len(rows)} rows). Analyse with scripts/analyze_kcca.py",
          flush=True)


if __name__ == "__main__":
    main()
