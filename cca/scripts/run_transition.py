"""Uncued -> cued transition: does the communication subspace change at task onset?

Most animals run several trials in the *uncued* corridor (world 4) before the
*cued* task (world 3). For the animals that have it, we fit the full pCCA
subspace readout (`subspace_window`) on the uncued running bins and on the
matched first-cued block (sample-count-matched to remove the N->correlation
confound), and compare the whole suite -- CC1, n_sig, mi_sig, IFI, optimal lag,
Gini -- plus the principal-angle ROTATION and member Jaccard BETWEEN the uncued
and cued subspaces (does the subspace reorient at task onset?). Across animals we
sign-test each delta, split by learner status.

Continuous regime; K is smaller (the uncued phase is short). No depth/ISI/
waveforms. Reuses tested subspace_window / membership / subspace.
"""

from __future__ import annotations

import csv
import dataclasses
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import (config, dataio, membership, paired_stats,  # noqa: E402
                     partial, subspace, subspace_window)

K = 15
N_FOLDS = 4
N_SHUFFLES = config.SURROGATE_SHUFFLES   # centralised (see config)
MAX_LAG = 8
ROT_DIMS = 3
MAX_SAMPLES = 6000            # cap matched block size (speed; stays >> 50*K)
WORLD_UNCUED, WORLD_CUED = 4, 3
DELTAS = ["cc1", "n_sig", "mi_sig", "ifi", "optimal_lag", "gini_x", "gini_y"]
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]


def _world_50ms(animal, cfg, n_bins):
    bm = int(cfg.temporal_bin_ms)
    with h5py.File(animal.streams_path, "r") as f:
        w1 = np.asarray(f["data_behaviour"]["world_binned"]).ravel()
    c = w1[bm // 2:: bm][:n_bins]
    if c.size < n_bins:
        c = np.concatenate([c, np.full(n_bins - c.size, np.nan)])
    return c


def _fit(X, Y, Z, groups):
    if Z is not None:
        X = partial.partial_out(X, Z)
        Y = partial.partial_out(Y, Z)
    return subspace_window.window_subspace(
        X, Y, groups, k=K, max_lag=MAX_LAG, n_shuffles=N_SHUFFLES, n_folds=N_FOLDS)


def _angle(wa, wb, d_use=ROT_DIMS):
    d = int(min(d_use, wa.shape[1], wb.shape[1]))
    if d < 1:
        return float("nan")
    qa, _ = np.linalg.qr(wa[:, :d])
    qb, _ = np.linalg.qr(wb[:, :d])
    return float(np.degrees(np.nanmax(subspace.principal_angles(qa, qb))))


def main():
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=25)
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    rng = np.random.default_rng(0)

    rows = []
    print(f"Uncued->cued transition (full suite) | K={K} pCCA, sample-matched\n")
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
        m = int(min(unc.sum(), cue_idx.size, MAX_SAMPLES))
        cue_early = np.zeros(n_bins, bool); cue_early[cue_idx[:m]] = True
        unc_idx = np.where(unc)[0]
        if unc_idx.size > m:
            unc_idx = rng.choice(unc_idx, m, replace=False)
        unc_sel = np.zeros(n_bins, bool); unc_sel[np.sort(unc_idx)] = True
        cue_groups = streams.trial_idx_50ms[cue_early]
        unc_groups = np.floor(np.linspace(0, N_FOLDS, unc_sel.sum(),
                                          endpoint=False)).astype(int)
        present = {}
        for area in config.AREAS:
            mm, idx = dataio.area_activity_50ms(a, area, cfg)
            if len(idx) >= cfg.min_units:
                present[area] = mm
        for ax, ay in PAIRS:
            if ax not in present or ay not in present:
                continue
            X, Y = present[ax], present[ay]
            others = [present[z] for z in present if z not in (ax, ay)]
            Z = np.concatenate(others, axis=1) if others else None
            wu = _fit(X[unc_sel], Y[unc_sel],
                      Z[unc_sel] if Z is not None else None, unc_groups)
            wc = _fit(X[cue_early], Y[cue_early],
                      Z[cue_early] if Z is not None else None, cue_groups)
            cc_u = float(wu.cc[0]) if wu.cc.size else np.nan
            cc_c = float(wc.cc[0]) if wc.cc.size else np.nan
            if not (np.isfinite(cc_u) and np.isfinite(cc_c)):
                continue
            row = {"animal": a.animal_id, "learner": int(is_learner),
                   "pair": f"{ax}-{ay}", "n_matched_bins": m,
                   "angle_x": round(_angle(wu.weights_x, wc.weights_x), 2),
                   "angle_y": round(_angle(wu.weights_y, wc.weights_y), 2),
                   "jac_x": round(membership.jaccard(wu.member_x, wc.member_x), 4),
                   "jac_y": round(membership.jaccard(wu.member_y, wc.member_y), 4)}
            for mt in DELTAS:
                u = wu.cc[0] if mt == "cc1" else getattr(wu, mt)
                c = wc.cc[0] if mt == "cc1" else getattr(wc, mt)
                row[f"unc_{mt}"] = round(float(u), 4)
                row[f"cue_{mt}"] = round(float(c), 4)
                row[f"d_{mt}"] = round(float(c) - float(u), 4)
            rows.append(row)
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): "
              f"uncued_run={int(unc.sum())} bins, "
              f"{len([r for r in rows if r['animal']==a.animal_id])} pairs")

    rdir = config.RESULTS_DIR
    cols = (["animal", "learner", "pair", "n_matched_bins", "angle_x", "angle_y",
             "jac_x", "jac_y"]
            + [f"{p}_{mt}" for mt in DELTAS for p in ("unc", "cue", "d")])
    with open(rdir / "transition_uncued_cued.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {rdir/'transition_uncued_cued.csv'} ({len(rows)} rows)\n")

    def tables(subset, title):
        print(f"[{title}]")
        for pair in [f"{x}-{y}" for x, y in PAIRS]:
            sub = [r for r in subset if r["pair"] == pair]
            if len(sub) < 3:
                continue
            line = f"  {pair:9s} n={len(sub):<2d}"
            ang = np.array([r["angle_x"] for r in sub if np.isfinite(r["angle_x"])])
            line += f" rot_x={np.nanmedian(ang):>5.1f}deg"
            for mt in ["cc1", "n_sig", "mi_sig", "ifi"]:
                d = np.array([r[f"d_{mt}"] for r in sub])
                _, med, _, p = paired_stats.wilcoxon_signed(d.tolist())
                s = "*" if (np.isfinite(p) and p < 0.05) else " "
                line += f" | d_{mt}={med:>+6.3f}p={p:.2g}{s}"
            print(line)

    print("Transition: delta = cued - uncued (sign test); rot = uncued->cued subspace angle")
    tables(rows, "all")
    tables([r for r in rows if r["learner"]], "learners")
    tables([r for r in rows if not r["learner"]], "non-learners")


if __name__ == "__main__":
    main()
