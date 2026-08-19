"""Cosine similarity of each canonical component to ITSELF across lags, with CCA
refit at every lag.

Asked before meeting item 4: *"look at cosine similarity of the significant CCs to
themselves across lags (when they're refit at each lag obviously)"*.

The lag curves already tell us how the canonical CORRELATION varies with lag. This asks
whether it is the same COMPONENT being measured at each lag, or a different one. The PCA
is fit once on the whole session and CCA is refit at every lag on those shared scores, so
the canonical weight vectors live in one common PC basis and are directly comparable
across lags.

Two readouts per (dim, lag), for each area separately:

  cos_same  |cos| between dim k's weight vector at lag 0 and dim k's at lag L.
  cos_best  the largest |cos| between dim k at lag 0 and ANY dim at lag L, plus which
            dim achieved it (`best_dim`).

Reporting both matters. If `cos_best` is high while `cos_same` is low, the component has
not disappeared — the canonical ORDER has swapped, and any analysis that matches
dimensions across fits by bare rank is silently comparing different components. That
failure has bitten this project twice (d38a833, 1bae90e), so it is measured here rather
than assumed.

**Absolute cosine** throughout: CCA determines each canonical vector only up to sign, so
a sign flip between lags is not a change in the component.

Significance is computed IN THIS DRIVER from its own lag-0 fit
(`fixed_subspace.frozen_perm_null`) rather than imported from `lag_curves`. Importing a
mask keyed to another fit's dimension order is precisely the bug this script exists to
quantify.

Data prep mirrors `run_lag_curves.py` (same 12k consecutive-block cap): this is a
question about subspace geometry across lag, not across epochs, so the cap's
trial-onset bias is not in play.

Writes one row per (animal, pair, dim, lag):
    results/lag_cosine_bin10{,_fsincl}.csv

Usage: PYTHONPATH=src python scripts/run_lag_cosine.py --bin-ms 10 --smooth-ms 2.5
       [--include-fs]
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import (config, core, dataio, fixed_subspace, lag_subspace,  # noqa: E402
                     partial)

K = 30
MAX_SAMPLES = 12000
MAX_BINS_PER_TRIAL = 600
N_DIMS = 10
FDR_DIMS = 10
PAIRS = [("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
         ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC")]
FIELDS = ["animal", "learner", "pair", "dim", "bin_ms", "lag_bins", "lag_ms",
          "cos_same_x", "cos_best_x", "best_dim_x",
          "cos_same_y", "cos_best_y", "best_dim_y",
          "cos_split_x", "cos_split_y", "split_best_x", "split_best_y",
          "r_lag", "sig", "p_perdim", "r_zero", "n_pairs"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin-ms", type=int, default=10)
    p.add_argument("--include-fs", action="store_true")
    p.add_argument("--smooth-ms", type=float, default=2.5)
    p.add_argument("--max-lag", type=int, default=25)
    p.add_argument("--n-shuffles", type=int, default=200)
    p.add_argument("--out", default="")
    return p.parse_args()



def _capped_index(trial_ids, max_samples=MAX_SAMPLES,
                  max_per_trial=MAX_BINS_PER_TRIAL):
    keep, total = [], 0
    for t in np.unique(trial_ids):
        pos = np.where(trial_ids == t)[0][:max_per_trial]
        keep.append(pos)
        total += pos.size
        if total >= max_samples:
            break
    return np.sort(np.concatenate(keep))


def _abs_cos(a, B):
    """|cos| between vector ``a`` and every column of ``B``."""
    na = np.linalg.norm(a)
    nb = np.linalg.norm(B, axis=0)
    ok = (na > 0) & (nb > 0)
    out = np.full(B.shape[1], np.nan)
    if na > 0:
        out[ok] = np.abs(a @ B[:, ok]) / (na * nb[ok])
    return out


def main():
    args = parse_args()
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=args.bin_ms,
                              exclude_fast_spiking=not args.include_fs,
                              gaussian_sd_ms=args.smooth_ms)
    suffix = "_fsincl" if args.include_fs else ""
    stem = args.out or f"lag_cosine_bin{args.bin_ms}{suffix}"
    animals = dataio.load_animals(config.DATA_DIR)
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, cfg, behaviour_lookup=behaviour)
    thr = cfg.velocity_thresh_cm_s
    lags = np.arange(-args.max_lag, args.max_lag + 1)
    print(f"LAG COSINE | bin={args.bin_ms}ms smooth={args.smooth_ms}ms K={K} | "
          f"lags ±{args.max_lag} bins | FS={'incl' if args.include_fs else 'excl'}\n")

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
        tid_full = streams.trial_idx_50ms[run]
        if np.unique(tid_full).size < 6:
            continue
        cap = _capped_index(tid_full)
        trial_ids = tid_full[cap]
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
            Xr = partial.partial_out(X, Z) if Z is not None else X
            Yr = partial.partial_out(Y, Z) if Z is not None else Y
            Sx, _ = core.pca_scores(Xr, min(K, Xr.shape[1]))
            Sy, _ = core.pca_scores(Yr, min(K, Yr.shape[1]))

            # --- SPLIT-HALF control ------------------------------------------
            # A cross-lag cosine is uninterpretable without knowing how much two CCA
            # fits disagree from sampling alone. Both terms are therefore HALF-data
            # fits: the reference is half 1 at lag 0, and every lag (including 0) is
            # fitted on half 2. cos(H1@0, H2@0) is then the pure-noise floor and
            # cos(H1@0, H2@L) adds only the lag, with the split noise held constant.
            # Comparing a full-data cross-lag cosine to a half-data floor would mix
            # two different noise levels — the mistake that made the 3-dim subspace
            # angle unmeasurable in item 3.
            uniq_t = np.unique(trial_ids)
            rng = np.random.default_rng(0)
            perm = rng.permutation(uniq_t)
            h1 = np.isin(trial_ids, perm[: perm.size // 2])
            h2 = ~h1
            A_h1 = B_h1 = None
            base_h1 = lag_subspace.segment_lag(Sx[h1], Sy[h1], None,
                                               trial_ids[h1], 0)
            if base_h1 is not None:
                try:
                    m_h1 = core.cca_fit(base_h1[0], base_h1[1])
                    A_h1, B_h1 = m_h1.A, m_h1.B
                except Exception:
                    A_h1 = B_h1 = None

            # --- reference fit at lag 0, and its own significance -------------
            base = lag_subspace.segment_lag(Sx, Sy, None, trial_ids, 0)
            if base is None:
                continue
            bx, by, _, _, _ = base
            try:
                m0 = core.cca_fit(bx, by)
                sig = fixed_subspace.frozen_perm_null(
                    bx, by, n_shuffles=args.n_shuffles, seed=0, correct="fdr",
                    fdr_dims=FDR_DIMS)
            except Exception as e:
                print(f"  base fail {a.animal_id} {ax}-{ay}: {e}"); continue
            A0, B0 = m0.A, m0.B
            d = int(min(N_DIMS, A0.shape[1], B0.shape[1]))
            r0 = np.asarray(m0.r, dtype=float)

            for lag in lags:
                packed = lag_subspace.segment_lag(Sx, Sy, None, trial_ids, int(lag))
                if packed is None:
                    continue
                xp, yp, _, _, _ = packed
                try:
                    mL = core.cca_fit(xp, yp)
                except Exception:
                    continue
                AL, BL = mL.A, mL.B
                rL = np.asarray(mL.r, dtype=float)
                # half-2 fit at this lag, for the split-half-matched comparison
                A2 = B2 = None
                if A_h1 is not None:
                    p2 = lag_subspace.segment_lag(Sx[h2], Sy[h2], None,
                                                  trial_ids[h2], int(lag))
                    if p2 is not None:
                        try:
                            m2 = core.cca_fit(p2[0], p2[1])
                            A2, B2 = m2.A, m2.B
                        except Exception:
                            A2 = B2 = None
                for dim in range(d):
                    cx = _abs_cos(A0[:, dim], AL)
                    cy = _abs_cos(B0[:, dim], BL)
                    jx = int(np.nanargmax(cx)) if np.any(np.isfinite(cx)) else -1
                    jy = int(np.nanargmax(cy)) if np.any(np.isfinite(cy)) else -1
                    if A_h1 is not None and A2 is not None and dim < A_h1.shape[1]:
                        sx_ = _abs_cos(A_h1[:, dim], A2)
                        sy_ = _abs_cos(B_h1[:, dim], B2)
                        kx_ = int(np.nanargmax(sx_)) if np.any(np.isfinite(sx_)) else -1
                        ky_ = int(np.nanargmax(sy_)) if np.any(np.isfinite(sy_)) else -1
                        cs_x = sx_[dim] if dim < sx_.size else np.nan
                        cs_y = sy_[dim] if dim < sy_.size else np.nan
                        sb_x = sx_[kx_] if kx_ >= 0 else np.nan
                        sb_y = sy_[ky_] if ky_ >= 0 else np.nan
                    else:
                        cs_x = cs_y = sb_x = sb_y = np.nan
                    rows.append({
                        "animal": a.animal_id, "learner": int(is_learner),
                        "pair": f"{ax}-{ay}", "dim": dim + 1, "bin_ms": args.bin_ms,
                        "lag_bins": int(lag), "lag_ms": int(lag) * args.bin_ms,
                        "cos_same_x": round(float(cx[dim]), 5)
                                      if dim < cx.size and np.isfinite(cx[dim]) else "",
                        "cos_best_x": round(float(cx[jx]), 5) if jx >= 0 else "",
                        "best_dim_x": jx + 1 if jx >= 0 else "",
                        "cos_same_y": round(float(cy[dim]), 5)
                                      if dim < cy.size and np.isfinite(cy[dim]) else "",
                        "cos_best_y": round(float(cy[jy]), 5) if jy >= 0 else "",
                        "best_dim_y": jy + 1 if jy >= 0 else "",
                        "cos_split_x": round(float(cs_x), 5)
                                       if np.isfinite(cs_x) else "",
                        "cos_split_y": round(float(cs_y), 5)
                                       if np.isfinite(cs_y) else "",
                        "split_best_x": round(float(sb_x), 5)
                                        if np.isfinite(sb_x) else "",
                        "split_best_y": round(float(sb_y), 5)
                                        if np.isfinite(sb_y) else "",
                        "r_lag": round(float(rL[dim]), 5) if dim < rL.size else "",
                        "sig": int(sig.mask[dim]) if dim < sig.mask.size else 0,
                        "p_perdim": round(float(sig.p[dim]), 5)
                                    if dim < sig.p.size else "",
                        "r_zero": round(float(r0[dim]), 5) if dim < r0.size else "",
                        "n_pairs": int(xp.shape[0]),
                    })
        _flush()
        print(f"  animal {a.animal_id} ({'L' if is_learner else 'n'}): "
              f"{len(rows)} rows total")
    _flush()
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
