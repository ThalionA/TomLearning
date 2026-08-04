"""Fixed-subspace figures — meeting 2026-07-28 items 5, 6, 7.

Reads results/fixed_subspace_bin10{,_fsincl}.csv and the reduced
fixed_subspace_epoch_bin10*.csv written by analyze_fixed_subspace.py.

Figure 1  HCV1_fixedsubspace_lagcurves_<fs>_bin10
          Items 5/6/7. One canonical component, frozen weights, lagged across time —
          one curve per learning epoch, mean ± shaded SEM across animals. This is the
          "lagged curves for naive vs exp" plot, with the refit confound removed.

Figure 2  HCV1_fixedsubspace_windows_<fs>_bin10
          Item 6. The integration window (half-max width of the same curve) and the
          curve's asymmetry (IFI), naive vs expert, one line per animal.

Positive lag = the FIRST-named area leads.

Usage: PYTHONPATH=src python scripts/figs_fixed_subspace.py [fsincl]
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
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
EPOCH_COLOUR = {"naive": "#2c6fbb", "intermediate": "#95a5a6", "expert": "#c0392b"}
EPOCHS = ["naive", "intermediate", "expert"]


def fig_curves(df: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[(df["pair"] == pair) & (df["dim"] == 1)]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        x_lead = pair.split("-")[0]
        for epoch in EPOCHS:
            e = sub[sub["epoch"] == epoch]
            if e.empty:
                continue
            piv = e.pivot_table(index="animal", columns="lag_ms", values="r")
            lags = np.array(sorted(piv.columns), dtype=float)
            M = piv[sorted(piv.columns)].to_numpy(float)
            n = np.sum(np.isfinite(M), axis=0)
            mean = np.nanmean(M, axis=0)
            sd = np.nanstd(M, axis=0, ddof=1)
            sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
            ax.plot(lags, mean, "-", color=EPOCH_COLOUR[epoch], lw=1.8, zorder=3,
                    label=f"{epoch} (n={piv.shape[0]})")
            ax.fill_between(lags, mean - sem, mean + sem,
                            color=EPOCH_COLOUR[epoch], alpha=0.16, zorder=2)
        ax.axvline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axhline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.set_title(pair, fontsize=10)
        ax.set_xlabel(f"lag (ms)   +: {x_lead} leads", fontsize=8)
        ax.set_ylabel("r (frozen CC₁ variates)", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper right", framealpha=0.6)
    fig.suptitle(f"One canonical component lagged across time, FIXED subspace — {fs}, "
                 "10 ms smoothed | weights identical across epochs (in-sample r)",
                 fontsize=12)
    figstyle.save(fig, ATT / f"HCV1_fixedsubspace_lagcurves_{fs}_bin10.png")
    print(f"wrote HCV1_fixedsubspace_lagcurves_{fs}_bin10.png")


def fig_windows(tab: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = tab[tab["pair"] == pair]
        nv = sub[sub["epoch"] == "naive"].set_index("animal")
        ex = sub[sub["epoch"] == "expert"].set_index("animal")
        common = [a for a in nv.index if a in ex.index]
        if not common:
            ax.set_title(f"{pair}\n(no paired epochs)", fontsize=9)
            ax.axis("off"); continue
        for a in common:
            ax.plot([0, 1], [nv.loc[a, "width_ms"], ex.loc[a, "width_ms"]], "-o",
                    color="#95a5a6", lw=0.9, ms=4, alpha=0.75, zorder=2)
        m = [np.nanmean([nv.loc[a, "width_ms"] for a in common]),
             np.nanmean([ex.loc[a, "width_ms"] for a in common])]
        ax.plot([0, 1], m, "-o", color="#c0392b", lw=2.4, ms=6, zorder=4)
        ax.set_xlim(-0.35, 1.35)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["naive", "expert"], fontsize=8)
        ax.set_ylabel("integration window (ms, half-max)", fontsize=8)
        # Flag ringing pairs: on a theta-oscillating curve this width is a half-period,
        # NOT an integration window (see analyze_fixed_subspace + GOTCHAS 2026-07-29).
        ringing = sub.dropna(subset=["side_peak_ms"])
        ringing = ringing[ringing["side_peak_ratio"] > 0.5]
        frac = len(ringing) / len(sub) if len(sub) else 0.0
        if frac >= 0.5:
            ax.set_title(f"{pair}  (n={len(common)})", fontsize=10, color="#c0392b")
            ax.text(0.5, 0.94, "θ-ringing: width = ½-period, NOT integration",
                    transform=ax.transAxes, ha="center", va="top", fontsize=6.5,
                    color="#c0392b",
                    bbox=dict(boxstyle="round,pad=0.2", fc="#fdecea", ec="#c0392b",
                              lw=0.6))
        else:
            ax.set_title(f"{pair}  (n={len(common)})  [readable window]", fontsize=9)
    fig.suptitle(f"Integration window of the frozen component, naive vs expert — {fs}, "
                 "10 ms smoothed | red = theta-ringing (width is a half-period, not an "
                 "integration window)", fontsize=11)
    figstyle.save(fig, ATT / f"HCV1_fixedsubspace_windows_{fs}_bin10.png")
    print(f"wrote HCV1_fixedsubspace_windows_{fs}_bin10.png")


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    src = RES / f"fixed_subspace_bin10{suf}.csv"
    red = RES / f"fixed_subspace_epoch_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/run_fixed_subspace.py first")
    df = pd.read_csv(src)
    df["r"] = pd.to_numeric(df["r"], errors="coerce")
    fig_curves(df, fs)
    if red.exists():
        fig_windows(pd.read_csv(red), fs)
    else:
        print(f"skip windows figure — {red.name} not found "
              "(run scripts/analyze_fixed_subspace.py)")


if __name__ == "__main__":
    main()
