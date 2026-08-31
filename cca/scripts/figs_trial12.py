"""Figures for the trial-1-vs-2 frozen-subspace arm (built 2026-08-20).

Renders the already-computed verdict tables (`trial12_delta_tests_bin10*.csv`,
`trial12_control_bin10*.csv` from analyze_trial12.py) — no recomputation. The
subspace is fit once per (animal, pair) on running-trial ordinals 11+ and
ordinals 1..10 are projected leak-free; deltas are trial 1 − trial 2 within
animal, tested animals-as-n. Summary bars are mean ± SEM across animals
(per-animal deltas live only in the 50 MB curves CSV, so no dot cloud here);
star = Wilcoxon p < 0.05 (uncorrected), filled square = BH survivor (none).

Outputs per FS condition:
  * HCV1_trial12_deltas_<fs>_bin10  — Δ strength | Δ IFI | Δ bins | Δ speed
  * HCV1_trial12_control_<fs>_bin10 — the 1→2 step z-scored against the 9
    adjacent-trial steps (ordinals 3..10 control band)

Usage: PYTHONPATH=src python scripts/figs_trial12.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402

figstyle.apply()
RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-CA3", "CA1-DG", "CA3-DG", "CA1-SUB", "CA1-RSC", "CA1-V1",
         "RSC-SUB", "V1-RSC"]
BAR = "#3b6fb6"
HILITE = "#c8552c"          # V1-RSC, the pre-registered candidate pair


def summary_axis(ax, df, metric, matched, title, xlab):
    sub = df[(df["metric"] == metric) & (df["matched"] == matched)]
    for y, pair in enumerate(PAIRS):
        r = sub[sub["pair"] == pair]
        if r.empty:
            continue
        r = r.iloc[0]
        col = HILITE if pair == "V1-RSC" else BAR
        ax.errorbar(r["mean"], y, xerr=r["sem"], fmt="s" if r["bh_pass"] else "o",
                    ms=6, color=col, mec="black", mew=0.6, capsize=3, zorder=4)
        if np.isfinite(r["w_p"]) and r["w_p"] < 0.05:
            ax.text(0.93, y, "*", transform=ax.get_yaxis_transform(),
                    ha="center", va="center", fontsize=15, fontweight="bold",
                    color=figstyle.STAR_COLOR, zorder=20)
        ax.text(1.02, y, f"n={int(r['n'])}", transform=ax.get_yaxis_transform(),
                va="center", fontsize=7, color="0.45")
    ax.axvline(0, color="0.6", lw=0.8, zorder=1)
    ax.set_yticks(range(len(PAIRS)), PAIRS, fontsize=8)
    ax.set_ylim(len(PAIRS) - 0.5, -0.9)
    ax.set_xlabel(xlab, fontsize=9)
    ax.set_title(title, fontsize=9)


def make_deltas(fs, tag):
    df = pd.read_csv(RES / f"trial12_delta_tests_bin10{tag}.csv")
    fig, axes = plt.subplots(1, 4, figsize=(12.6, 3.4), sharey=True,
                             constrained_layout=True)
    summary_axis(axes[0], df, "r0_allsig_w", 1,
                 "Δ strength (r0, all sig CCs,\nbin-matched)", "Δ r (trial 1 − 2)")
    summary_axis(axes[1], df, "ifi_allsig_w", 1,
                 "Δ direction (IFI ±50 ms,\nall sig CCs, bin-matched)",
                 "Δ IFI (trial 1 − 2)")
    summary_axis(axes[2], df, "d_n_bins", 0,
                 "Δ behaviour: running bins\n(raw)", "Δ bins (trial 1 − 2)")
    summary_axis(axes[3], df, "d_vel_mean", 0,
                 "Δ behaviour: speed\n(raw)", "Δ cm/s (trial 1 − 2)")
    fig.suptitle(f"Trial 1 vs trial 2 through a frozen subspace — communication "
                 f"null, behaviour robust ({'FS excluded' if fs == 'fsexcl' else 'FS included'})",
                 fontsize=11)
    figstyle.save(fig, f"HCV1_trial12_deltas_{fs}_bin10")


def make_control(fs, tag):
    df = pd.read_csv(RES / f"trial12_control_bin10{tag}.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4), sharey=True,
                             constrained_layout=True)
    for ax, metric, title in [(axes[0], "r0_cc1", "strength (CC1 r0)"),
                              (axes[1], "ifi_cc1", "direction (CC1 IFI)")]:
        sub = df[(df["metric"] == metric) & (df["stat"] == "z")]
        for y, pair in enumerate(PAIRS):
            r = sub[sub["pair"] == pair]
            if r.empty:
                continue
            r = r.iloc[0]
            col = HILITE if pair == "V1-RSC" else BAR
            ax.errorbar(r["mean"], y, xerr=r["sem"], fmt="o", ms=6, color=col,
                        mec="black", mew=0.6, capsize=3, zorder=4)
            if np.isfinite(r["w_p"]) and r["w_p"] < 0.05:
                ax.text(0.93, y, "*", transform=ax.get_yaxis_transform(),
                        ha="center", va="center", fontsize=15,
                        fontweight="bold", color=figstyle.STAR_COLOR, zorder=20)
            ax.text(1.02, y, f"n={int(r['n'])}", transform=ax.get_yaxis_transform(),
                    va="center", fontsize=7, color="0.45")
        ax.axvline(0, color="0.6", lw=0.8)
        ax.set_yticks(range(len(PAIRS)), PAIRS, fontsize=8)
        ax.set_ylim(len(PAIRS) - 0.5, -0.9)
        ax.set_xlabel("z of 1→2 step vs adjacent steps", fontsize=9)
        ax.set_title(title, fontsize=10)
    fig.suptitle(f"Is the first-experience step special? 1→2 vs the ordinal "
                 f"3..10 control band ({'FS excluded' if fs == 'fsexcl' else 'FS included'})",
                 fontsize=11)
    figstyle.save(fig, f"HCV1_trial12_control_{fs}_bin10")


def main():
    for fs, tag in [("fsexcl", ""), ("fsincl", "_fsincl")]:
        make_deltas(fs, tag)
        make_control(fs, tag)
        print(f"{fs}: 2 figures written")


if __name__ == "__main__":
    main()
