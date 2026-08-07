"""Item 4 figures — integration window vs IFI, naive vs experienced, frozen axes.

Reads results/ifi_windows_epochs_bin10{,_fsincl}.csv
(analyze_ifi_windows_epochs.py).

Figure 1  HCV1_ifi_windows_epochs_<fs>_bin10
          THE ASKED-FOR PLOT. Per pair, CC1's IFI as a function of the integration
          window (|lag| <= w, 10..250 ms), one line per learning epoch, mean +/- SEM
          across animals. Zero = no directional asymmetry; the dotted vertical is the
          report's headline +/-50 ms window.

Figure 2  HCV1_ifi_windows_epochs_fffb_<fs>_bin10
          The same, but pooling all significant CCs within their FF/FB label instead of
          taking CC1. Signed IFI is never averaged across opposite labels — that would
          cancel by construction.

Positive IFI => the FIRST-named area leads. Correlations are in-sample (frozen fit),
so these are contrast statistics.

Usage: PYTHONPATH=src python scripts/figs_ifi_windows_epochs.py [fsincl]
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
EPOCHS = ["naive", "intermediate", "expert"]
EPOCH_COLOUR = {"naive": "#2c6fbb", "intermediate": "#95a5a6", "expert": "#c0392b"}
HEADLINE_MS = 50


def _curve(ax, sub, pair):
    """Mean +/- SEM IFI vs window, one line per epoch. Animals-as-n: an animal's CCs are
    averaged before the across-animal mean, so an animal with many CCs cannot dominate."""
    drawn = False
    for epoch in EPOCHS:
        g = sub[sub["epoch"] == epoch]
        if g.empty:
            continue
        piv = g.pivot_table(index="animal", columns="window_ms", values="ifi",
                            aggfunc="mean")
        if piv.empty:
            continue
        wins = np.array(sorted(piv.columns), dtype=float)
        M = piv[sorted(piv.columns)].to_numpy(float)
        n = np.sum(np.isfinite(M), axis=0)
        mean = np.nanmean(M, axis=0)
        sd = np.nanstd(M, axis=0, ddof=1)
        sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
        ax.plot(wins, mean, "-", color=EPOCH_COLOUR[epoch], lw=2.0, zorder=3,
                label=f"{epoch} (n={piv.shape[0]})")
        ax.fill_between(wins, mean - sem, mean + sem, color=EPOCH_COLOUR[epoch],
                        alpha=0.16, zorder=2)
        drawn = True
    if drawn:
        ax.axhline(0.0, color="k", lw=1.0, ls="--", zorder=1)
        ax.axvline(HEADLINE_MS, color="#666", lw=0.9, ls=":", zorder=1)
        x_lead = pair.split("-")[0]
        ax.set_xlabel("integration window ±w (ms)", fontsize=8)
        ax.set_ylabel(f"IFI   +: {x_lead} leads", fontsize=8)
        ax.legend(fontsize=6.5, loc="best", framealpha=0.6)
    return drawn


def fig_cc1(df, fs):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[(df["pair"] == pair) & (df["dim"] == 1) & (df["sig"] == 1) &
                 (df["degenerate"] == 0)]
        if sub.empty or not _curve(ax, sub, pair):
            ax.set_title(f"{pair}\n(no significant CC1)", fontsize=9)
            ax.axis("off"); continue
        ax.set_title(pair, fontsize=10)
    fig.suptitle(f"Integration window vs IFI, by learning epoch — {fs} | CC1, FROZEN "
                 f"axes so the same component is compared across epochs", fontsize=11)
    figstyle.save(fig, ATT / f"HCV1_ifi_windows_epochs_{fs}_bin10.png")
    print(f"wrote HCV1_ifi_windows_epochs_{fs}_bin10.png")


def fig_fffb(df, fs):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        base = df[(df["pair"] == pair) & (df["sig"] == 1) & (df["degenerate"] == 0)]
        if base.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        x_lead = pair.split("-")[0]
        any_drawn = False
        for label, ls in [("FF", "-"), ("FB", "--")]:
            sub = base[base["label"] == label]
            if sub.empty:
                continue
            for epoch in EPOCHS:
                g = sub[sub["epoch"] == epoch]
                if g.empty:
                    continue
                piv = g.pivot_table(index="animal", columns="window_ms", values="ifi",
                                    aggfunc="mean")
                if piv.empty:
                    continue
                wins = np.array(sorted(piv.columns), dtype=float)
                M = piv[sorted(piv.columns)].to_numpy(float)
                ax.plot(wins, np.nanmean(M, axis=0), ls,
                        color=EPOCH_COLOUR[epoch], lw=1.6, zorder=3,
                        label=f"{label} {epoch}")
                any_drawn = True
        if not any_drawn:
            ax.set_title(f"{pair}\n(no labelled CCs)", fontsize=9)
            ax.axis("off"); continue
        ax.axhline(0.0, color="k", lw=1.0, ls="--", zorder=1)
        ax.axvline(HEADLINE_MS, color="#666", lw=0.9, ls=":", zorder=1)
        ax.set_title(pair, fontsize=10)
        ax.set_xlabel("integration window ±w (ms)", fontsize=8)
        ax.set_ylabel(f"IFI   +: {x_lead} leads", fontsize=8)
        ax.legend(fontsize=5.5, ncol=2, loc="best", framealpha=0.6)
    fig.suptitle(f"Integration window vs IFI by epoch, split by FF/FB label — {fs} | "
                 f"solid = FF, dashed = FB; signed IFI never averaged across labels",
                 fontsize=10.5)
    figstyle.save(fig, ATT / f"HCV1_ifi_windows_epochs_fffb_{fs}_bin10.png")
    print(f"wrote HCV1_ifi_windows_epochs_fffb_{fs}_bin10.png")


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    src = RES / f"ifi_windows_epochs_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/analyze_ifi_windows_epochs.py first")
    df = pd.read_csv(src)
    df["ifi"] = pd.to_numeric(df["ifi"], errors="coerce")
    fig_cc1(df, fs)
    fig_fffb(df, fs)


if __name__ == "__main__":
    main()
