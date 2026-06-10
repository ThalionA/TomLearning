"""Driver: NONLINEAR (kernel) CCA on the uncued->cued condition contrast.

The kernel analogue of run_transition.py. For each animal with both conditions,
sample-matched, fits the `kcca_window` suite separately on the uncued (world 4)
and cued (world 3) running bins and reports per-condition held-out KCCA CC,
linear CC (superiority), significance, IFI (direction), Gini (membership), plus
the cued-uncued delta. FS-excluded by default; --include-fs for the variant.
Writes results/kcca_transition.csv (or _fsincl).
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, dataio, kcca_window  # noqa: E402

K = 15
REG = 10.0
CAP = 900
MAX_LAG = 4
N_SHUFFLES = 10
N_FOLDS = 4
D = 3
WORLD_UNCUED, WORLD_CUED = 4, 3
METRICS = ["cc_kcca", "cc_linear", "n_sig", "mi_sig", "ifi", "gini_x"]
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--include-fs", action="store_true")
    return p.parse_args()


def _world_50ms(animal, cfg, n_bins):
    bm = int(cfg.temporal_bin_ms)
    with h5py.File(animal.streams_path, "r") as f:
        w1 = np.asarray(f["data_behaviour"]["world_binned"]).ravel()
    c = w1[bm // 2:: bm][:n_bins]
    if c.size < n_bins:
        c = np.concatenate([c, np.full(n_bins - c.size, np.nan)])
    return c


def _metrics(X, Y, Z, groups, seed):
    r = kcca_window.kcca_window(X, Y, groups, Z=Z, reg=REG, k=K, max_lag=MAX_LAG,
                                n_shuffles=N_SHUFFLES, n_folds=N_FOLDS, d=D, seed=seed)
    return {"cc_kcca": float(r.cc_kcca[0]) if r.cc_kcca.size else np.nan,
            "cc_linear": float(r.cc_linear[0]) if r.cc_linear.size else np.nan,
            "n_sig": r.n_sig, "mi_sig": r.mi_sig, "ifi": r.ifi, "gini_x": r.gini_x}


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=25,
                              exclude_fast_spiking=not args.include_fs)
    suffix = "_fsincl" if args.include_fs else ""
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    rng = np.random.default_rng(0)
    print(f"KCCA uncued->cued | K={K} reg={REG} cap={CAP} | "
          f"FS {'incl' if args.include_fs else 'excl'}\n", flush=True)

    rows = []
    out = config.RESULTS_DIR / f"kcca_transition{suffix}.csv"
    cols = (["animal", "learner", "pair", "n_matched"]
            + [f"{p}_{mt}" for mt in METRICS for p in ("unc", "cue", "d")])

    def _write():
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    for a in animals:
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        n_bins = streams.vel_50ms.size
        world = _world_50ms(a, cfg, n_bins)
        run = streams.vel_50ms >= thr
        unc = run & (world == WORLD_UNCUED)
        cue = run & (world == WORLD_CUED) & (~np.isnan(streams.trial_idx_50ms))
        if unc.sum() < K * 30 or cue.sum() < K * 30:
            continue
        is_learner = a.animal_id in entries
        cue_idx = np.where(cue)[0]
        m = int(min(unc.sum(), cue_idx.size, CAP))
        cue_sel = np.zeros(n_bins, bool); cue_sel[cue_idx[:m]] = True
        unc_idx = np.where(unc)[0]
        if unc_idx.size > m:
            unc_idx = rng.choice(unc_idx, m, replace=False)
        unc_sel = np.zeros(n_bins, bool); unc_sel[np.sort(unc_idx)] = True
        # Position-based CV folds for BOTH conditions: the sample-matched blocks
        # are short (matched to the few-trial uncued phase), so they often span
        # < N_FOLDS whole trials -> whole-trial CV is impossible. Contiguous
        # position folds are used identically for cued and uncued, so the
        # cued-uncued delta stays fair (mild within-trial CV leakage, both sides).
        cue_groups = np.floor(np.linspace(0, N_FOLDS, int(cue_sel.sum()),
                                          endpoint=False)).astype(int)
        unc_groups = np.floor(np.linspace(0, N_FOLDS, int(unc_sel.sum()),
                                          endpoint=False)).astype(int)
        present = {}
        for area in config.AREAS:
            mm, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = mm
        n_pairs = 0
        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            try:
                u = _metrics(X[unc_sel], Y[unc_sel],
                             Z[unc_sel] if Z is not None else None, unc_groups,
                             seed=a.animal_id)
                c = _metrics(X[cue_sel], Y[cue_sel],
                             Z[cue_sel] if Z is not None else None, cue_groups,
                             seed=a.animal_id + 1)
            except Exception as exc:                 # noqa: BLE001
                print(f"  {a.animal_id} {ax}-{ay}: fail {type(exc).__name__}",
                      flush=True)
                continue
            row = {"animal": a.animal_id, "learner": int(is_learner),
                   "pair": f"{ax}-{ay}", "n_matched": m}
            for mt in METRICS:
                row[f"unc_{mt}"] = round(float(u[mt]), 4)
                row[f"cue_{mt}"] = round(float(c[mt]), 4)
                row[f"d_{mt}"] = round(float(c[mt]) - float(u[mt]), 4)
            rows.append(row); n_pairs += 1
        _write()
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): "
              f"{n_pairs} pairs ({len(rows)} total)", flush=True)

    _write()
    print(f"\nwrote {out} ({len(rows)} rows). Analyse with scripts/analyze_kcca.py",
          flush=True)


if __name__ == "__main__":
    main()
