"""Very-early-trial (pre-10) readout driver — does the communication subspace
change in the FIRST few trials and then plateau?

Two complementary readouts per animal × area-pair (see `src/tom_cca/early_trials.py`
for why the split is forced by sample size — a single trial is ~580 running bins,
far too few to re-fit a 30-component (p)CCA):

  (A) PROJECTED per-trial readout  -> results/early_trials_projected{,_fsincl}.csv
      Fit the subspace ONCE on the animal's later trials (ordinal >10, sample-rich)
      and project trials 1 / 4 / 7 / 10 through it, leak-free. Gives a true
      per-trial CC, IFI / optimal lag, and participation-Gini.

  (B) BLOCK refit readout          -> results/early_trials_blocks{,_fsincl}.csv
      Re-fit the full suite on cumulative early blocks [1-5] / [1-7] / [1-10] and a
      late reference (last 10 trials) for the metrics that NEED a fit: held-out CC,
      #sig dims, MI, weight-Gini, IFI, split-half floor, the rotation angle vs the
      late subspace, and (unless --no-kcca) kernel CCA. Smallest unit = 5 trials.

The analysis (`analyze_early_trials.py`) does trial-1-vs-4/7/10 (A) and
block-vs-late (B) paired Wilcoxon across animals + the plateau check.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import (config, dataio, early_trials, kcca_window,  # noqa: E402
                     subspace_window)

K = 30
N_FOLDS = 5
MAX_LAG = 10
N_SHUFFLES = config.SURROGATE_SHUFFLES
SAT = 0.99
ROT_DIMS = 3
EARLY_ORDINALS = [1, 2, 4, 7, 10]     # the trials the user asked about (2 added 2026-08-20)
BLOCK_SIZES = [5, 7, 10]              # cumulative early blocks [1..n]
LATE_N = 10                           # late reference = last LATE_N trials
KCCA_CAP, KCCA_REG, KCCA_LAG, KCCA_D, KCCA_SHUF = 900, 10.0, 4, 3, 10
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]

PROJ_FIELDS = ["animal", "learner", "pair", "ordinal", "trial_id", "n_bins",
               "cc1", "ifi", "optimal_lag", "gini_part_x", "gini_part_y"]
BLOCK_FIELDS = ["animal", "learner", "pair", "block", "n_trials", "n_bins",
                "cc1", "n_sig", "mi_sig", "ifi", "optimal_lag", "gini_x", "gini_y",
                "sh_x", "sh_y", "angle_x", "angle_y",
                "cc_kcca", "gini_kcca_x", "gini_kcca_y"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--include-fs", action="store_true",
                   help="include fast-spiking units (default: exclude)")
    p.add_argument("--no-kcca", action="store_true",
                   help="skip the (slow) kernel-CCA block metrics")
    p.add_argument("--bin-ms", type=int, default=config.DEFAULT.temporal_bin_ms)
    p.add_argument("--smooth-ms", type=float, default=0.0,
                   help="Gaussian s.d. (ms) for spike-train smoothing (Buzsáki = 2.5)")
    return p.parse_args()


def _subsample(n, cap, seed):
    if n <= cap:
        return np.ones(n, dtype=bool)
    idx = np.sort(np.random.default_rng(seed).choice(n, cap, replace=False))
    m = np.zeros(n, dtype=bool); m[idx] = True
    return m


def _rotation(ref_w, cur_w, d_use=ROT_DIMS) -> float:
    """Max principal angle (deg) between two dominant subspaces; NaN if no ref."""
    if ref_w is None or cur_w is None:
        return float("nan")
    d = int(min(d_use, ref_w.shape[1], cur_w.shape[1]))
    if d < 1:
        return float("nan")
    qa, _ = np.linalg.qr(ref_w[:, :d])
    qb, _ = np.linalg.qr(cur_w[:, :d])
    from tom_cca import subspace
    return float(np.degrees(np.nanmax(subspace.principal_angles(qa, qb))))


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs,
                              gaussian_sd_ms=args.smooth_ms)
    bintag = f"_bin{args.bin_ms}" if args.bin_ms != 25 else ""
    suffix = bintag + ("_fsincl" if args.include_fs else "")
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    n_learn = sum(1 for a in animals if a.animal_id in entries)
    print(f"EARLY-TRIALS | ALL {len(animals)} animals ({n_learn} learners) | "
          f"bin={args.bin_ms}ms K={K} ordinals={EARLY_ORDINALS} blocks={BLOCK_SIZES}+late "
          f"| kcca={'off' if args.no_kcca else 'on'}\n")

    proj_out = config.RESULTS_DIR / f"early_trials_projected{suffix}.csv"
    block_out = config.RESULTS_DIR / f"early_trials_blocks{suffix}.csv"
    proj_rows, block_rows = [], []

    def _flush():
        with open(proj_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=PROJ_FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(proj_rows)
        with open(block_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=BLOCK_FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(block_rows)

    for a in animals:
        try:
            streams = dataio._load_temporal_streams(a, cfg)
        except Exception:
            continue
        run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= thr)
        trial_ids = streams.trial_idx_50ms[run]
        uniq = np.unique(trial_ids)
        if uniq.size < max(BLOCK_SIZES) + LATE_N:        # need early blocks + a disjoint late ref
            continue
        is_learner = a.animal_id in entries
        present = {}
        for area in config.AREAS:
            m, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = m[run]

        early_ids = [uniq[o - 1] for o in EARLY_ORDINALS]   # 1st,4th,7th,10th trial
        train_ids = uniq[max(EARLY_ORDINALS):]              # reference = ordinal >10
        late_ids = uniq[-LATE_N:]
        blocks = [(f"t1-{n}", uniq[:n]) for n in BLOCK_SIZES] + [("late", late_ids)]

        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            pair = f"{ax}-{ay}"

            # ---- (A) projected per-trial readout (fit on later trials) ----
            try:
                metr = early_trials.early_trial_metrics(
                    X, Y, trial_ids, early_ids, train_ids, Z=Z, k=K, max_lag=MAX_LAG)
                for o, r in zip(EARLY_ORDINALS, metr):
                    proj_rows.append({
                        "animal": a.animal_id, "learner": int(is_learner), "pair": pair,
                        "ordinal": o, "trial_id": int(r["trial"]), "n_bins": r["n_bins"],
                        "cc1": round(r["cc1"], 4) if np.isfinite(r["cc1"]) else "",
                        "ifi": round(r["ifi"], 4) if np.isfinite(r["ifi"]) else "",
                        "optimal_lag": r["optimal_lag"],
                        "gini_part_x": round(r["gini_part_x"], 4) if np.isfinite(r["gini_part_x"]) else "",
                        "gini_part_y": round(r["gini_part_y"], 4) if np.isfinite(r["gini_part_y"]) else ""})
            except Exception as e:
                print(f"  proj fail {a.animal_id} {pair}: {e}")

            # ---- (B) block refit readout (fit late first for the angle ref) ----
            ws_by_block, late_w = {}, (None, None)
            for name, btr in blocks:
                idx = np.where(np.isin(trial_ids, btr))[0]
                groups = trial_ids[idx]
                if np.unique(groups).size < N_FOLDS:
                    continue
                ws = subspace_window.window_subspace(
                    X[idx], Y[idx], groups, Z=Z[idx] if Z is not None else None,
                    k=K, max_lag=MAX_LAG, n_shuffles=N_SHUFFLES, n_folds=N_FOLDS)
                ws_by_block[name] = (ws, idx.size, int(np.unique(groups).size))
                if name == "late":
                    late_w = (ws.weights_x, ws.weights_y)
            for name, btr in blocks:
                if name not in ws_by_block:
                    continue
                ws, nb, ntr = ws_by_block[name]
                cc1 = float(ws.cc[0]) if ws.cc.size else float("nan")
                row = {"animal": a.animal_id, "learner": int(is_learner), "pair": pair,
                       "block": name, "n_trials": ntr, "n_bins": nb,
                       "cc1": round(cc1, 4) if np.isfinite(cc1) else "",
                       "n_sig": ws.n_sig, "mi_sig": round(ws.mi_sig, 4),
                       "ifi": round(ws.ifi, 4), "optimal_lag": ws.optimal_lag,
                       "gini_x": round(ws.gini_x, 4), "gini_y": round(ws.gini_y, 4),
                       "sh_x": round(ws.split_half_x, 2) if np.isfinite(ws.split_half_x) else "",
                       "sh_y": round(ws.split_half_y, 2) if np.isfinite(ws.split_half_y) else "",
                       "angle_x": "", "angle_y": "",
                       "cc_kcca": "", "gini_kcca_x": "", "gini_kcca_y": ""}
                if name != "late":
                    ax_ang = _rotation(late_w[0], ws.weights_x)
                    ay_ang = _rotation(late_w[1], ws.weights_y)
                    row["angle_x"] = round(ax_ang, 2) if np.isfinite(ax_ang) else ""
                    row["angle_y"] = round(ay_ang, 2) if np.isfinite(ay_ang) else ""
                if not args.no_kcca:
                    idx = np.where(np.isin(trial_ids, btr))[0]
                    sub = _subsample(idx.size, KCCA_CAP, seed=a.animal_id)
                    gi = idx[sub]; grp = trial_ids[gi]
                    if np.unique(grp).size >= N_FOLDS:
                        try:
                            kr = kcca_window.kcca_window(
                                X[gi], Y[gi], grp, Z=Z[gi] if Z is not None else None,
                                reg=KCCA_REG, k=K, max_lag=KCCA_LAG,
                                n_shuffles=KCCA_SHUF, n_folds=N_FOLDS, d=KCCA_D,
                                seed=a.animal_id)
                            row["cc_kcca"] = round(float(kr.cc_kcca[0]), 4) if kr.cc_kcca.size else ""
                            row["gini_kcca_x"] = round(kr.gini_x, 4)
                            row["gini_kcca_y"] = round(kr.gini_y, 4)
                        except Exception as e:
                            print(f"  kcca fail {a.animal_id} {pair} {name}: {e}")
                block_rows.append(row)
        _flush()
        print(f"  animal {a.animal_id}: {len(uniq)} trials done "
              f"(proj rows {len(proj_rows)}, block rows {len(block_rows)})")
    _flush()
    print(f"\nprojected -> {proj_out}\nblocks    -> {block_out}")


if __name__ == "__main__":
    main()
