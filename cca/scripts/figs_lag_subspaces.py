"""Lagged-subspace figures — meeting 2026-07-28 items 2, 3, 4.

Reads results/lag_subspaces_bin10{,_fsincl}.csv.

Figure 1  HCV1_lagsubspace_stability_<fs>_bin10
          Item 3. Principal angle between the lag-0 subspace and the lagged one, vs lag,
          per pair. One line per area; the grey band is that pair's split-half noise
          floor. Inside the band = the same subspace read at a different delay.

Figure 2  HCV1_lagsubspace_fffb_<fs>_bin10
          Items 2/4. Feedforward (+50 ms) vs feedback (-50 ms), one line per animal:
          held-out CC1 and connection-specific Gini. The right-hand panel shows the
          FF-vs-FB subspace angle against the same floor — the gate on whether FF and FB
          are separable at all.

Usage: python scripts/figs_lag_subspaces.py [fsincl]
"""
from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from _common import RES, TEMPORAL, config, fs_from_argv
import figstyle

figstyle.apply()

PAIRS = list(config.PAIR_NAMES)
C_X, C_Y, C_FLOOR = "#c0392b", "#2c6fbb", "#95a5a6"
TAU_MS = TEMPORAL.label_w_bins * TEMPORAL.bin_ms   # 50 ms


def _mean_sem(piv, lags):
    M = piv.reindex(columns=lags).to_numpy(float)
    n = np.sum(np.isfinite(M), axis=0)
    mean = np.nanmean(M, axis=0)
    sd = np.nanstd(M, axis=0, ddof=1)
    return mean, np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)


def fig_stability(df: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[df["pair"] == pair]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        ax_name, ay_name = pair.split("-")
        lags = sorted(sub["lag_ms"].unique())
        # CC1 is primary: the 3-dim split-half floor is ~79 deg (two halves of the SAME
        # data are near-orthogonal), so the 3-dim curve has no headroom and is drawn
        # faintly only to show that. See analyze_lag_subspaces.stability.
        for col, colour, label in [("angle_x_cc1", C_X, ax_name),
                                   ("angle_y_cc1", C_Y, ay_name)]:
            piv = sub.pivot_table(index="animal", columns="lag_ms", values=col)
            mean, sem = _mean_sem(piv, lags)
            ax.plot(lags, mean, "-o", color=colour, lw=1.6, ms=3, zorder=3,
                    label=f"{label} CC₁ (n={piv.shape[0]})")
            ax.fill_between(lags, mean - sem, mean + sem, color=colour, alpha=0.18,
                            zorder=2)
        for col, colour in [("angle_x", C_X), ("angle_y", C_Y)]:
            piv = sub.pivot_table(index="animal", columns="lag_ms", values=col)
            mean, _ = _mean_sem(piv, lags)
            ax.plot(lags, mean, ":", color=colour, lw=1.0, alpha=0.5, zorder=2)
        fl = float(pd.concat([sub["floor_x_cc1"], sub["floor_y_cc1"]]).mean())
        fl_sd = float(pd.concat([sub["floor_x_cc1"],
                                 sub["floor_y_cc1"]]).std(ddof=1))
        ax.axhspan(fl - fl_sd, fl + fl_sd, color=C_FLOOR, alpha=0.30, zorder=0)
        ax.axhline(fl, color=C_FLOOR, lw=1.2, ls="--", zorder=1)
        fl3 = np.nanmean([sub["floor_x"].mean(), sub["floor_y"].mean()])
        ax.axhline(fl3, color=C_FLOOR, lw=1.0, ls=":", zorder=1)
        ax.set_ylim(0, 92)
        ax.axvline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.set_title(pair, fontsize=10)
        ax.set_xlabel(f"lag (ms)   +: {ax_name} leads", fontsize=8)
        ax.set_ylabel("angle vs lag-0 subspace (°)", fontsize=8)
        ax.legend(fontsize=6.5, loc="lower right", framealpha=0.6)
    fig.suptitle(f"Subspace stability across lag — {fs}, 10 ms smoothed | solid = CC₁, "
                 "dotted = 3-dim (not estimable) | grey band = CC₁ split-half floor",
                 fontsize=12)
    figstyle.save(fig, f"HCV1_lagsubspace_stability_{fs}_bin10.png")
    print(f"wrote HCV1_lagsubspace_stability_{fs}_bin10.png")


def fig_fffb(df: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[df["pair"] == pair]
        ff = sub[sub["lag_ms"] == TAU_MS].set_index("animal")
        fb = sub[sub["lag_ms"] == -TAU_MS].set_index("animal")
        common = [a for a in ff.index if a in fb.index]
        if not common:
            ax.set_title(f"{pair}\n(no ±{TAU_MS} ms fit)", fontsize=9)
            ax.axis("off"); continue
        ax_name = pair.split("-")[0]
        for a in common:
            ax.plot([0, 1], [ff.loc[a, "cc1"], fb.loc[a, "cc1"]], "-o",
                    color=C_FLOOR, lw=0.9, ms=4, alpha=0.75, zorder=2)
        m = [np.nanmean([ff.loc[a, "cc1"] for a in common]),
             np.nanmean([fb.loc[a, "cc1"] for a in common])]
        ax.plot([0, 1], m, "-o", color=C_X, lw=2.4, ms=6, zorder=4)
        ax.set_xlim(-0.35, 1.35)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"FF (+{TAU_MS} ms)\n{ax_name} leads",
                            f"FB (−{TAU_MS} ms)"], fontsize=7)
        ax.set_ylabel("held-out CC₁", fontsize=8)
        ang = float(ff.loc[common, "angle_ff_fb_cc1"].mean())
        flr = float(ff.loc[common, "floor_x_cc1"].mean())
        ax.set_title(f"{pair}  (n={len(common)})\n"
                     f"FF/FB angle {ang:.0f}° vs floor {flr:.0f}°", fontsize=8.5)
    fig.suptitle(f"Feedforward vs feedback subspace — {fs}, 10 ms smoothed | "
                 "one line per animal; title compares the FF/FB angle to the noise "
                 "floor", fontsize=12)
    figstyle.save(fig, f"HCV1_lagsubspace_fffb_{fs}_bin10.png")
    print(f"wrote HCV1_lagsubspace_fffb_{fs}_bin10.png")


def main():
    suf = fs_from_argv()                       # "" (FS-excluded) or "_fsincl"
    fs = "fsincl" if suf else "fsexcl"
    src = RES / f"lag_subspaces_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/run_lag_subspaces.py first")
    df = pd.read_csv(src)
    for c in ("cc1", "angle_x", "angle_y", "angle_x_cc1", "angle_y_cc1", "floor_x",
              "floor_y", "floor_x_cc1", "floor_y_cc1", "angle_ff_fb_x",
              "angle_ff_fb_y",
              "angle_ff_fb_cc1", "gini_x_conn", "gini_y_conn"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    fig_stability(df, fs)
    fig_fffb(df, fs)


if __name__ == "__main__":
    main()
