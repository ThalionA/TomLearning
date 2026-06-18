"""Trajectory dims-as-n figure: per-dimension directionality (IFI) and strength (CC)
change with learning, shown under THREE units side by side so the pseudoreplication
inflation is visible at a glance.

For each pair (y-axis) and each learning axis (columns: trial_frac, lp_rel,
performance), three markers give the slope estimate under:
  * ANIM — animals-as-n, per-animal slope + signed-rank (the honest inferential unit);
  * LMM  — dims-as-n, cluster-robust random-intercept mixed model (powerful-but-honest);
  * OLS  — dims-as-n, naive pooled OLS (the inflated literal Buzsáki unit).
Filled marker = p < 0.05; open = n.s. Where only OLS is filled, that is manufactured
significance; where LMM is filled but ANIM is at its small-n floor, the dims unit is
out-powering a real, same-signed animal-level trend (read case by case).

Rows: IFI (top), CC (bottom). Usage: ``python figs_trajectory_dims.py [bin50|fsincl]``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402
from analyze_trajectory_dims import AXES, PAIRS, compute_table, load  # noqa: E402

figstyle.apply()
ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"

UNITS = [  # (key prefix, label, marker, colour, y-offset)
    ("anim", "animals-as-n (paired t)", "o", "#23231f", +0.24),
    ("lmm", "dims-as-n LMM (cluster-robust)", "s", "#1f78b4", 0.0),
    ("ols", "dims-as-n OLS (inflated)", "^", "#c0392b", -0.24),
]
AXIS_LABEL = {"trial_frac": "trial fraction", "lp_rel": "trials rel. learning point",
              "performance": "performance"}


def _slope_key(unit):
    return "anim_slope" if unit == "anim" else f"{unit}_b"


def _panel(ax, tab, metric, axis, order, ylabels):
    rows = {r["pair"]: r for r in tab if r["metric"] == metric and r["axis"] == axis}
    ypos = {p: len(order) - 1 - i for i, p in enumerate(order)}  # first pair on top
    for unit, _lab, mk, col, dy in UNITS:
        xs, ys, fill = [], [], []
        for p in order:
            r = rows.get(p)
            if r is None:
                continue
            b, pv = r[_slope_key(unit)], r[f"{unit}_p"]
            if not np.isfinite(b):
                continue
            xs.append(b); ys.append(ypos[p] + dy)
            fill.append(np.isfinite(pv) and pv < 0.05)
        xs, ys, fill = np.array(xs), np.array(ys), np.array(fill, bool)
        if xs.size:
            if fill.any():
                ax.scatter(xs[fill], ys[fill], marker=mk, s=52, c=col,
                           edgecolors=col, zorder=5)
            if (~fill).any():
                ax.scatter(xs[~fill], ys[~fill], marker=mk, s=52, facecolors="none",
                           edgecolors=col, linewidths=1.4, zorder=4)
    ax.axvline(0, color="#9a9990", lw=1, ls="--", zorder=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(list(reversed(order)) if ylabels else [])
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_title(f"{metric.upper()} vs {AXIS_LABEL[axis]}", fontsize=12)
    ax.tick_params(labelsize=10)


def main():
    arg = " ".join(sys.argv[1:])
    df, path = load(arg)
    if df is None:
        print(f"missing {path.name} — run run_trajectory (per-dim) first"); return
    pool = "pool" in arg
    cohort = df if pool else df[df["learner"] == 1]
    n_an = cohort["animal"].nunique()
    tab = compute_table(df, pool=pool)
    tag = ("fsincl" if "fsincl" in arg else "fsexcl") + ("_pooled" if pool else "")
    binms = "10" if "10" in arg else ("50" if "50" in arg else "25")
    order = [p for p in PAIRS if any(r["pair"] == p for r in tab)]  # fixed rows, all panels

    fig, axes = plt.subplots(2, len(AXES), figsize=(4.4 * len(AXES), 8.4))
    for i, metric in enumerate(["ifi", "cc"]):
        for j, axis in enumerate(AXES):
            _panel(axes[i, j], tab, metric, axis, order, ylabels=(j == 0))
    for j, axis in enumerate(AXES):
        axes[1, j].set_xlabel("slope (per-dim metric vs learning axis)")

    handles = [plt.Line2D([0], [0], marker=mk, color="none", markerfacecolor=col,
                          markeredgecolor=col, markersize=9, label=lab)
               for _, lab, mk, col, _ in UNITS]
    handles.append(plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="none",
                              markeredgecolor="#444", markersize=9, label="open marker = n.s. (p≥0.05)"))
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=10.5,
               frameon=False)
    cohort_lbl = f"all {n_an} animals (learners + non)" if pool else f"{n_an} learners"
    fig.suptitle(f"Directionality (IFI) & strength (CC) vs learning — three units of analysis "
                 f"\n({binms} ms, {tag}, {cohort_lbl}; filled = p<0.05)", fontsize=12.5)
    out = ATT / f"HCV1_CCA_{tag}_trajdims_bin{binms}.png"
    figstyle.save(fig, out)
    print(f"wrote {out.name}  ({len(tab)} pair×metric×axis cells, {n_an} animals)")


if __name__ == "__main__":
    main()
