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
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, mixed_effects, paired_stats, trajectory  # noqa: E402

AXES = ["trial_frac", "performance", "lp_rel"]
# gini_pearson_x = CCA-independent control for the gini_x de-sparsification headline.
# It is interpreted via its within-animal SLOPE only — its absolute level carries a
# partner-unit-count noise floor, so it is excluded from the cross-pair LEVELS table.
# gini_*_conn / gini_*_sig = the CONNECTION-SPECIFIC contribution Ginis
# (membership.subspace_contribution_connection). The plain gini_x/gini_y are
# provably partner-invariant (area-intrinsic, PROJECT_LOG 2026-07-28), so the
# "participation broadens" headline must be re-tested on these. ⚠ gini_*_sig is
# NaN whenever n_sig = 0, so its slopes are conditioned on windows that had a
# significant dim — read it alongside the n_sig slope, never alone.
SLOPE_METRICS = ["cc1", "n_sig", "mi_sig", "ifi", "gini_x", "gini_pearson_x",
                 "gini_x_conn", "gini_y_conn", "gini_x_sig", "gini_y_sig",
                 "rot_x", "jac_x"]
LEVEL_METRICS = ["cc1", "n_sig", "mi_sig", "ifi", "optimal_lag", "gini_x",
                 "gini_x_conn", "gini_x_sig",
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
    name = sys.argv[1] if len(sys.argv) > 1 else "trajectory_windows.csv"
    path = config.RESULTS_DIR / name
    df = pd.read_csv(path)
    # Drop any battery metric absent from this CSV (e.g. gini_pearson_x in a
    # pre-control run) so the script still runs on older tables.
    present = set(df.columns)
    slope_metrics = [m for m in SLOPE_METRICS if m in present]
    level_metrics = [m for m in LEVEL_METRICS if m in present]
    learn = df[df["learner"] == 1]
    print(f"{path.name}: {len(df)} window-rows, "
          f"{df['animal'].nunique()} animals ({learn['animal'].nunique()} learners)\n")

    # Rotation vs the split-half noise floor (per-animal median; sign test rot-floor>0).
    # Run for the top-3 subspace AND, when the columns are present (post-#14 runs),
    # the CC1-only dominant dim — whose floor is far tighter, so a CC1-only rotation
    # is a more sensitive reorientation test than the top-3 angle (§3.4).
    def _rot_floor(rot_col, sh_col, label):
        if rot_col not in df.columns or sh_col not in df.columns:
            return
        print("=" * 70)
        print(f"ROTATION vs NOISE FLOOR — {label} (per-animal median; sign test rot - floor > 0)")
        print("=" * 70)
        for pair in PAIR_ORDER:
            g = learn[learn["pair"] == pair]
            diffs = []
            for _, sub in g.groupby("animal"):
                r = _num(sub[rot_col]).median(); f = _num(sub[sh_col]).median()
                if np.isfinite(r) and np.isfinite(f):
                    diffs.append(r - f)
            if len(diffs) < 3:
                continue
            _, med, _, p = paired_stats.wilcoxon_signed(diffs)
            star = "  *" if (np.isfinite(p) and p < 0.05) else ""
            mr = _num(g[rot_col]).median(); mf = _num(g[sh_col]).median()
            print(f"    {pair:9s} n={len(diffs):<2d} rot={mr:>5.1f}° floor={mf:>5.1f}° "
                  f"Δ(med)={med:>+6.1f}° p={p:.3g}{star}")
        print()

    _rot_floor("rot_x", "sh_x", "top-3 subspace")
    _rot_floor("rot_x_cc1", "sh_x_cc1", "CC1 only (dominant dim)")

    print("=" * 70)
    print("LEVELS (mean over learner windows; IFI/lag one-sample sign test vs 0)")
    print("=" * 70)
    hdr = "  ".join(f"{m:>14}" for m in level_metrics)
    print(f"{'pair':9s}  {hdr}")
    for pair in PAIR_ORDER:
        g = learn[learn["pair"] == pair]
        if g.empty:
            continue
        cells = []
        for m in level_metrics:
            v = _num(g[m]).to_numpy()
            v = v[np.isfinite(v)]
            cells.append(f"{np.mean(v):>14.3f}" if v.size else f"{'-':>14}")
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
    for metric in slope_metrics:
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

    # Parametric / powerful tests for the key metrics: per-animal-slope t-test
    # AND a random-slope LMM population slope (pools all windows).
    print("\n" + "=" * 70)
    print("PARAMETRIC slope tests (learners): Wilcoxon | t-test(slopes) | LMM(all windows)")
    print("=" * 70)
    param_metrics = [m for m in ["gini_x", "gini_x_conn", "gini_x_sig",
                                 "gini_y_conn", "gini_y_sig",
                                 "gini_pearson_x", "ifi", "mi_sig",
                                 "cc1", "n_sig"] if m in present]
    for metric in param_metrics:
        print(f"\n[{metric}]")
        for axis in AXES:
            print(f"  axis={axis}")
            sl = per_animal_slopes(learn, metric, axis)
            for pair in PAIR_ORDER:
                vals = [s for (pr, s) in sl if pr == pair and np.isfinite(s)]
                if len(vals) < 3:
                    continue
                _, _, wp = (None, None, paired_stats.wilcoxon_signed(vals)[3])
                t_p = float(stats.ttest_1samp(vals, 0.0).pvalue)
                g = learn[learn["pair"] == pair][["animal", axis, metric]].rename(
                    columns={"animal": "animal_id", axis: "axis", metric: "value"})
                lmm = mixed_effects.lmm_slope(g.to_dict("records"))
                lp_p = lmm["p"] if lmm["ok"] else float("nan")
                mark = lambda p: "*" if (np.isfinite(p) and p < 0.05) else " "
                print(f"    {pair:9s} n={len(vals):<2d} med={np.median(vals):>+8.4f} "
                      f"| W p={wp:.3g}{mark(wp)} | t p={t_p:.3g}{mark(t_p)} "
                      f"| LMM β={lmm['estimate']:+.4f} p={lp_p:.3g}{mark(lp_p)}")


if __name__ == "__main__":
    main()
