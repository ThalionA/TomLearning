"""Item 1 figures — IFI of each CC across integration windows, and the sign structure.

Reads results/cc_ifi_windows_bin10{,_fsincl}.csv (analyze_cc_ifi_signs.py).

Figure 1  HCV1_cc_ifi_windows_<fs>_bin10
          "Calculate IFI at different lags for each CC and look at them." Per pair,
          IFI vs integration window (|lag| <= w), one line per canonical dimension,
          mean +/- SEM across animals. The zero line is the question: do different CCs
          sit on opposite sides of it?

Figure 2  HCV1_cc_ifi_signmap_<fs>_bin10
          The same thing per animal, at the headline +/-50 ms window: a diverging
          heatmap of IFI, rows = animals, columns = canonical dimension. Red = the
          first-named area leads, blue = it lags. A row with both colours is an animal
          whose CCs disagree in direction.

Positive IFI => the FIRST-named area leads.

Usage: PYTHONPATH=src python scripts/figs_cc_ifi_signs.py [fsincl]
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
ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"
RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
SHOW_DIMS = 6          # leading canonical dims to draw
HEADLINE_W = 5         # bins = +/-50 ms
MAP_DIMS = 8


def fig_windows(df: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    cmap = plt.get_cmap("viridis")
    for ax, pair in zip(axes, PAIRS):
        # ONLY CCs that beat the circular-shift null. Ungated, CC3-CC8 clear it in just
        # 14-17% of cells with held-out CC at the floor, so their IFI sign is noise and
        # plotting them manufactures apparent sign heterogeneity.
        sub = df[(df["pair"] == pair) & (df["dim"] <= SHOW_DIMS) &
                 (df["degenerate"] == 0) & (df["sig"] == 1)]
        allsub = df[(df["pair"] == pair) & (df["dim"] <= SHOW_DIMS)]
        if sub.empty:
            ax.set_title(f"{pair}\n(no significant CCs)", fontsize=9)
            ax.axis("off"); continue
        x_lead = pair.split("-")[0]
        n_an = sub["animal"].nunique()
        for d in range(1, SHOW_DIMS + 1):
            g = sub[sub["dim"] == d]
            if g["animal"].nunique() < 3:      # need a few animals to average
                continue
            rate = (allsub[allsub["dim"] == d].groupby("animal")["sig"].max().mean()
                    if len(allsub[allsub["dim"] == d]) else np.nan)
            piv = g.pivot_table(index="animal", columns="window_ms", values="ifi")
            wins = np.array(sorted(piv.columns), dtype=float)
            M = piv[sorted(piv.columns)].to_numpy(float)
            n = np.sum(np.isfinite(M), axis=0)
            mean = np.nanmean(M, axis=0)
            sd = np.nanstd(M, axis=0, ddof=1)
            sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
            col = cmap((d - 1) / max(1, SHOW_DIMS - 1))
            ax.plot(wins, mean, "-", color=col, lw=1.7, zorder=3,
                    label=f"CC{d} ({rate:.0%} sig, n={piv.shape[0]})")
            ax.fill_between(wins, mean - sem, mean + sem, color=col, alpha=0.13,
                            zorder=2)
        ax.axhline(0.0, color="k", lw=1.0, ls="--", zorder=1)
        ax.axvline(HEADLINE_W * 10, color="#888", lw=0.8, ls=":", zorder=0)
        ax.set_title(f"{pair}  (n={n_an})", fontsize=10)
        ax.set_xlabel("integration window ±w (ms)", fontsize=8)
        ax.set_ylabel(f"IFI   +: {x_lead} leads", fontsize=8)
        ax.legend(fontsize=6, ncol=2, loc="best", framealpha=0.6)
    fig.suptitle(f"IFI of each canonical dimension vs integration window — {fs}, "
                 "10 ms smoothed | CCA refit at every lag | ONLY CCs beating the "
                 "circular-shift null | +IFI = first area leads", fontsize=12)
    figstyle.save(fig, ATT / f"HCV1_cc_ifi_windows_{fs}_bin10.png")
    print(f"wrote HCV1_cc_ifi_windows_{fs}_bin10.png")


def fig_signmap(df: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    sub_all = df[(df["window_bins"] == HEADLINE_W) & (df["dim"] <= MAP_DIMS) &
                 (df["sig"] == 1)]
    vmax = float(np.nanpercentile(np.abs(sub_all["ifi"].to_numpy(float)), 95)) or 1.0
    im = None
    for ax, pair in zip(axes, PAIRS):
        sub = sub_all[sub_all["pair"] == pair]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        piv = sub.pivot_table(index="animal", columns="dim", values="ifi")
        piv = piv.reindex(columns=range(1, MAP_DIMS + 1))
        M = piv.to_numpy(float)
        # NaN = that CC did not beat the null for that animal. In RdBu_r a NaN and a
        # genuine IFI ~ 0 would both render near-white, so non-significant cells get an
        # explicit grey instead — "no reliable coupling" must not look like "balanced".
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#dcdcdc")
        im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-vmax, vmax=vmax,
                       aspect="auto")
        # mark animals whose CCs disagree in sign
        n_mixed = 0
        for r in range(M.shape[0]):
            row = M[r][np.isfinite(M[r]) & (M[r] != 0)]
            if row.size and (row > 0).any() and (row < 0).any():
                n_mixed += 1
                ax.text(-0.9, r, "◀", va="center", ha="center", fontsize=7,
                        color="#c0392b")
        x_lead = pair.split("-")[0]
        ax.set_title(f"{pair}   {n_mixed}/{M.shape[0]} animals mixed", fontsize=9)
        ax.set_xticks(range(MAP_DIMS))
        ax.set_xticklabels([f"CC{d}" for d in range(1, MAP_DIMS + 1)], fontsize=6)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels([f"{a} ({int(np.isfinite(M[i]).sum())})"
                            for i, a in enumerate(piv.index)], fontsize=6)
        ax.set_xlabel(f"red: {x_lead} leads", fontsize=7)
    if im is not None:
        cb = fig.colorbar(im, ax=list(axes), fraction=0.015, pad=0.01)
        cb.set_label("IFI (red = first area leads)", fontsize=8)
    fig.suptitle(f"Sign of each CC's IFI per animal, ±{HEADLINE_W * 10} ms window — "
                 f"{fs} | ONLY CCs beating the null (grey = not significant) | "
                 "◀ = that animal's significant CCs disagree in direction", fontsize=12)
    figstyle.save(fig, ATT / f"HCV1_cc_ifi_signmap_{fs}_bin10.png")
    print(f"wrote HCV1_cc_ifi_signmap_{fs}_bin10.png")


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    src = RES / f"cc_ifi_windows_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/analyze_cc_ifi_signs.py first")
    df = pd.read_csv(src)
    df["ifi"] = pd.to_numeric(df["ifi"], errors="coerce")
    fig_windows(df, fs)
    fig_signmap(df, fs)


if __name__ == "__main__":
    main()
