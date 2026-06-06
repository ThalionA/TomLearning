"""Full subspace readout on the naive / intermediate / expert learning epochs.

For each learner and area pair, fit pCCA (the full `subspace_window` suite) on
the continuous running bins of each of the three project epochs
(`dataio.epoch_windows` from the animal's learning point), giving per
(animal, pair, epoch) metrics for early/intermediate/post-LP contrasts. Writes
results/epoch_metrics.csv (or _fsincl). Learners only (epochs need an LP).
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, dataio, partial, subspace_window  # noqa: E402

K = 20
N_FOLDS = 5
MAX_LAG = 8
N_SHUFFLES = 100              # circular-shift surrogate for dim significance (95th pct)
SAT = 0.99
EPOCHS = ["naive", "intermediate", "expert"]
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]
FIELDS = ["animal", "learner", "pair", "epoch", "n_bins", "cc1", "n_sig",
          "mi_sig", "ifi", "optimal_lag", "gini_x", "gini_y"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=25)
    p.add_argument("--include-fs", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs)
    suffix = "_fsincl" if args.include_fs else ""
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    print(f"Epochs (naive/intermediate/expert) | bin={args.bin_ms}ms K={K} pCCA "
          f"| FS {'included' if args.include_fs else 'excluded'}\n")

    rows = []
    dim_rows = []                                     # one row per significant canonical dim
    for a in animals:
        if a.animal_id not in entries:               # epochs require an LP
            continue
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids = streams.trial_idx_50ms[run]
        uniq = np.unique(trial_ids)
        ew = dataio.epoch_windows(int(entries[a.animal_id].lp), uniq.size, cfg)
        if ew is None:
            print(f"  animal {a.animal_id}: epoch_windows None (skip)")
            continue
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
            for epoch in EPOCHS:
                pos = ew[epoch]
                pos = pos[(pos >= 0) & (pos < uniq.size)]
                trials = uniq[pos]
                mask = np.isin(trial_ids, trials)
                if mask.sum() < K * 30:
                    continue
                Xi, Yi = X[mask], Y[mask]
                if Z is not None:
                    Zi = Z[mask]
                    Xi = partial.partial_out(Xi, Zi)
                    Yi = partial.partial_out(Yi, Zi)
                ws = subspace_window.window_subspace(
                    Xi, Yi, trial_ids[mask], k=K, max_lag=MAX_LAG,
                    n_shuffles=N_SHUFFLES, n_folds=N_FOLDS)
                cc1 = float(ws.cc[0]) if ws.cc.size else float("nan")
                if not np.isfinite(cc1) or cc1 >= SAT:
                    continue
                rows.append({
                    "animal": a.animal_id, "learner": 1, "pair": f"{ax}-{ay}",
                    "epoch": epoch, "n_bins": int(mask.sum()), "cc1": round(cc1, 4),
                    "n_sig": ws.n_sig, "mi_sig": round(ws.mi_sig, 4),
                    "ifi": round(ws.ifi, 4), "optimal_lag": ws.optimal_lag,
                    "gini_x": round(ws.gini_x, 4), "gini_y": round(ws.gini_y, 4)})
                n_rows += 1
                for d in range(int(ws.cc.size)):      # per-dimension (dims-as-n)
                    dim_rows.append({
                        "animal": a.animal_id, "pair": f"{ax}-{ay}", "epoch": epoch,
                        "dim": d + 1, "peak_cc": round(float(ws.cc[d]), 4),
                        "sig": int(bool(ws.sig_mask[d])) if d < ws.sig_mask.size else 0,
                        "ifi": round(float(ws.ifi_per_dim[d]), 4)
                        if d < ws.ifi_per_dim.size else "",
                        "lag": int(ws.lag_per_dim[d]) if d < ws.lag_per_dim.size else ""})
        print(f"  animal {a.animal_id}: lp={entries[a.animal_id].lp} {n_rows} rows")

    out = config.RESULTS_DIR / f"epoch_metrics{suffix}.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    out_d = config.RESULTS_DIR / f"epoch_dims{suffix}.csv"
    with open(out_d, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["animal", "pair", "epoch", "dim",
                                          "peak_cc", "sig", "ifi", "lag"],
                           lineterminator="\n")
        w.writeheader(); w.writerows(dim_rows)
    print(f"\nwrote {out} ({len(rows)} rows) + {out_d} ({len(dim_rows)} dim-rows). "
          "Analyse with scripts/analyze_epochs.py / analyze_dims_as_n.py")


if __name__ == "__main__":
    main()
