"""Trajectory dims-as-n: do the per-dimension directionality (IFI) and strength (CC)
change with learning when the unit is the canonical DIMENSION (Gonzalez/Buzsáki),
and does that agree with the animal-level conclusion?

Reads results/trajectory_w15_bin{25,50}{,_fsincl}_dims.csv (one row per canonical
dim per window: animal, pair, trial_frac/lp_rel/performance, dim, cc, ifi, lag, sig).
For each pair × metric (ifi, cc) × learning axis, reports THREE views:
  * ANIMALS-as-n  — per-animal OLS slope (over that animal's significant-dim windows),
                    signed-rank across animals (the inferential unit);
  * dims-as-n LMM — random-intercept LMM over all sig-dim windows (cluster-robust:
                    the powerful-but-honest dims unit);
  * dims-as-n OLS — naive pooled OLS over all sig-dim windows (the *inflated* literal
                    Buzsáki unit, shown to expose pseudoreplication).
Where OLS is significant but the animal level / LMM is not, that gap is the inflation.

Run with no args for 25 ms FS-excluded; pass e.g. ``bin50`` / ``fsincl``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, mixed_effects, paired_stats, trajectory  # noqa: E402

PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
AXES = ["trial_frac", "lp_rel", "performance"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _ols_p(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 5 or np.std(x) == 0:
        return np.nan, np.nan
    n = x.size
    b1 = np.polyfit(x, y, 1)[0]
    yhat = np.polyval([b1, np.mean(y) - b1 * np.mean(x)], x)
    sse = np.sum((y - yhat) ** 2)
    se = np.sqrt(sse / (n - 2)) / (np.std(x) * np.sqrt(n)) if n > 2 else np.nan
    from scipy import stats
    p = float(2 * stats.t.sf(abs(b1 / se), n - 2)) if se and np.isfinite(se) and se > 0 else np.nan
    return float(b1), p


def _per_animal_slopes(g, metric, axis):
    out = []
    for _, ga in g.groupby("animal"):
        s, _ = trajectory.linear_slope(_num(ga[axis]).to_numpy(), _num(ga[metric]).to_numpy())
        if np.isfinite(s):
            out.append(s)
    return out


def load(arg=""):
    """Read the per-dim trajectory CSV named by ``arg`` (``""``=25 ms FS-excl)."""
    binms = "50" if "50" in arg else "25"
    suffix = "_fsincl" if "fsincl" in arg else ""
    # driver writes "{out}_dims{suffix}.csv" → _dims BEFORE _fsincl (run_trajectory.py)
    path = config.RESULTS_DIR / f"trajectory_w15_bin{binms}_dims{suffix}.csv"
    if not path.is_file():
        return None, path
    return pd.read_csv(path), path


def compute_table(df, min_animals=3, min_dimwin=8):
    """Per pair × metric (ifi, cc) × axis, the three views as tidy records.

    Returns a list of dicts with keys: pair, metric, axis, n_an, n_dimwin,
    anim_slope/anim_p (signed-rank over per-animal slopes), lmm_b/lmm_p
    (cluster-robust random-intercept LMM), ols_b/ols_p (naive pooled OLS).
    Only learner + SIGNIFICANT dims are used (the dims-as-n convention)."""
    learn = df[(df["learner"] == 1) & (df["sig"] == 1)].copy()
    out = []
    for metric in ["ifi", "cc"]:
        for axis in AXES:
            for pair in PAIRS:
                g = learn[learn["pair"] == pair].dropna(subset=[metric, axis])
                if g["animal"].nunique() < min_animals or len(g) < min_dimwin:
                    continue
                sl = _per_animal_slopes(g, metric, axis)
                a_med, a_p = (np.nan, np.nan)
                if len(sl) >= 3:
                    _, a_med, _, a_p = paired_stats.wilcoxon_signed(sl)
                recs = [{"animal_id": int(r["animal"]), "axis": float(r[axis]),
                         "value": float(r[metric])} for _, r in g.iterrows()]
                lmm = mixed_effects.lmm_slope(recs, value="value", axis="axis")
                ob, op = _ols_p(g[axis], g[metric])
                out.append(dict(pair=pair, metric=metric, axis=axis,
                                n_an=len(sl), n_dimwin=len(g),
                                anim_slope=a_med, anim_p=a_p,
                                lmm_b=lmm.get("estimate", np.nan), lmm_p=lmm.get("p", np.nan),
                                ols_b=ob, ols_p=op))
    return out


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else ""
    df, path = load(arg)
    if df is None:
        print(f"missing {path.name} — run run_trajectory (per-dim) first"); return
    learn = df[(df["learner"] == 1) & (df["sig"] == 1)]
    print(f"{path.name}: {len(df)} dim-rows, {df['animal'].nunique()} animals; "
          f"{len(learn)} significant learner dim-windows\n")
    tab = compute_table(df)
    star = lambda p: "*" if (np.isfinite(p) and p < 0.05) else " "  # noqa: E731
    for metric in ["ifi", "cc"]:
        for axis in AXES:
            print("=" * 92)
            print(f"[{metric} vs {axis}]  ANIMALS-as-n (signed-rank) | dims-as-n LMM (cluster) "
                  "| dims-as-n OLS (inflated)")
            print("=" * 92)
            for r in [r for r in tab if r["metric"] == metric and r["axis"] == axis]:
                print(f"  {r['pair']:8s} n_an={r['n_an']:<2d} n_dimwin={r['n_dimwin']:<4d} | "
                      f"ANIM med={r['anim_slope']:+.3f} p={r['anim_p']:.2g}{star(r['anim_p'])} | "
                      f"LMM b={r['lmm_b']:+.3f} p={r['lmm_p']:.2g}{star(r['lmm_p'])} | "
                      f"OLS b={r['ols_b']:+.3f} p={r['ols_p']:.2g}{star(r['ols_p'])}")
            print()


if __name__ == "__main__":
    main()
