"""Fast analysis of the full-suite trajectory CSV (no re-fitting).

Reads results/trajectory_windows.csv (written by run_trajectory.py) and reports,
for every subspace metric:

* LEVELS  -- per-pair mean over learner windows, with a one-sample sign test
  vs 0 for the directionality metrics (IFI, optimal lag): is communication
  consistently directed?
* SLOPES  -- per (animal, pair) slope of the metric on each learning axis
  (trial fraction / performance / LP-relative), then an across-animal sign test
  per pair, split learners vs non-learners: does the metric track learning?

Re-runnable in seconds; all the expensive pCCA fitting lives in run_trajectory.py.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, paired_stats, trajectory  # noqa: E402

AXES = ["trial_frac", "performance", "lp_rel"]
SLOPE_METRICS = ["cc1", "n_sig", "mi_sig", "ifi", "gini_x", "rot_x", "jac_x"]
LEVEL_METRICS = ["cc1", "n_sig", "mi_sig", "ifi", "optimal_lag", "gini_x",
                 "rot_x", "jac_x"]
PAIR_ORDER = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG",
              "CA1-SUB", "RSC-SUB", "V1-RSC"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def per_animal_slopes(df, metric, axis):
    out = []
    for (an, pr), g in df.groupby(["animal", "pair"]):
        x = _num(g[axis]).to_numpy()
        y = _num(g[metric]).to_numpy()
        if np.sum(np.isfinite(x) & np.isfinite(y)) >= 4:
            s, _ = trajectory.linear_slope(x, y)
            out.append((pr, s))
    return out


def sign_table(pairs_slopes, indent="    "):
    by_pair = defaultdict(list)
    for pr, s in pairs_slopes:
        if np.isfinite(s):
            by_pair[pr].append(s)
    for pair in PAIR_ORDER:
        sl = np.array(by_pair.get(pair, []))
        if sl.size < 3:
            continue
        _, med, _, p = paired_stats.wilcoxon_signed(sl.tolist())
        star = "  *" if (np.isfinite(p) and p < 0.05) else ""
        print(f"{indent}{pair:9s} n={sl.size:<2d} med={med:>+9.4f} "
              f"up={int(np.sum(sl > 0))}/{sl.size:<2d} p={p:.3g}{star}")


def main():
    path = config.RESULTS_DIR / "trajectory_windows.csv"
    df = pd.read_csv(path)
    learn = df[df["learner"] == 1]
    print(f"{path.name}: {len(df)} window-rows, "
          f"{df['animal'].nunique()} animals ({learn['animal'].nunique()} learners)\n")

    print("=" * 70)
    print("LEVELS (mean over learner windows; IFI/lag one-sample sign test vs 0)")
    print("=" * 70)
    hdr = "  ".join(f"{m:>9}" for m in LEVEL_METRICS)
    print(f"{'pair':9s}  {hdr}")
    for pair in PAIR_ORDER:
        g = learn[learn["pair"] == pair]
        if g.empty:
            continue
        cells = []
        for m in LEVEL_METRICS:
            v = _num(g[m]).to_numpy()
            v = v[np.isfinite(v)]
            cells.append(f"{np.mean(v):>9.3f}" if v.size else f"{'-':>9}")
        print(f"{pair:9s}  " + "  ".join(cells))
    # directionality: is IFI consistently non-zero per pair (per-animal mean)?
    print("\n  IFI directionality (per-animal mean IFI, sign test vs 0; "
          ">0 = X leads Y):")
    for pair in PAIR_ORDER:
        g = learn[learn["pair"] == pair]
        if g.empty:
            continue
        per_an = [_num(sub["ifi"]).mean() for _, sub in g.groupby("animal")]
        per_an = [x for x in per_an if np.isfinite(x)]
        if len(per_an) < 3:
            continue
        _, med, _, p = paired_stats.wilcoxon_signed(per_an)
        star = "  *" if (np.isfinite(p) and p < 0.05) else ""
        print(f"    {pair:9s} n={len(per_an):<2d} mean_IFI={med:>+7.4f} "
              f"p={p:.3g}{star}")

    print("\n" + "=" * 70)
    print("SLOPES vs learning axes (across-animal sign test of per-animal slope)")
    print("=" * 70)
    for metric in SLOPE_METRICS:
        print(f"\n[{metric}]")
        for axis in AXES:
            sl_learn = per_animal_slopes(learn, metric, axis)
            print(f"  axis={axis} (learners):")
            sign_table(sl_learn)
            if axis != "lp_rel":
                nonl = df[df["learner"] == 0]
                sl_non = per_animal_slopes(nonl, metric, axis)
                if any(np.isfinite(s) for _, s in sl_non):
                    print(f"  axis={axis} (non-learners):")
                    sign_table(sl_non)


if __name__ == "__main__":
    main()
