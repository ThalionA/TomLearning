"""Figure for the IFI lag-integration-window sweep — which window is cleanest?

From results/ifi_windows_bin{25,50}{,_fsincl}.csv: one panel per pair, across-animal
mean IFI +/- SEM vs the integration window |lag| <= w (faint per-animal lines), with
the cleanest window (max |mean/SEM|) marked. Renders into the vault attachments.
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
from tom_cca import config  # noqa: E402

ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _fig(df, bin_ms, tag):
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharex=True)
    ws = np.array(sorted(df["window_bins"].unique()))
    xs = ws * bin_ms
    for ax, pair in zip(axes.ravel(), PAIRS):
        g = df[df["pair"] == pair]
        if g.empty:
            ax.set_visible(False); continue
        per_an = {}
        for an, ga in g.groupby("animal"):
            per_an[an] = {int(r): v for r, v in zip(ga["window_bins"], _num(ga["ifi"]))}
        for an, m in per_an.items():
            ax.plot(xs, [m.get(int(w), np.nan) for w in ws], color="#bdc3c7",
                    lw=0.6, alpha=0.6, zorder=1)
        means = np.array([_num(g[g["window_bins"] == w]["ifi"]).mean() for w in ws])
        sems = np.array([_num(g[g["window_bins"] == w]["ifi"]).std(ddof=1)
                         / np.sqrt(max(1, _num(g[g["window_bins"] == w]["ifi"]).notna().sum()))
                         for w in ws])
        ax.errorbar(xs, means, yerr=sems, color="#c0392b", lw=2, marker="o",
                    capsize=3, zorder=3)
        tstat = np.abs(means / np.where(sems > 0, sems, np.nan))
        if np.any(np.isfinite(tstat)):
            wbest = xs[int(np.nanargmax(tstat))]
            ax.axvline(wbest, color="#2980b9", ls="--", lw=1.2,
                       label=f"cleanest = {int(wbest)} ms")
            ax.legend(fontsize=8, loc="best")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(pair, fontsize=10)
    for ax in axes[1, :]:
        ax.set_xlabel("integration window |lag| (ms)")
    for ax in axes[:, 0]:
        ax.set_ylabel("IFI (>0: first area leads)")
    fig.suptitle(f"Directionality vs lag-integration window — {tag} "
                 "(held-out, segment-aware lag curve; faint = animals, red = mean ± SEM; "
                 "dashed = cleanest, max |mean/SEM|)", fontsize=12)
    fig.tight_layout()
    fig.savefig(ATT / f"HCV1_ifi_windows_{tag}.png", dpi=140)
    plt.close(fig)


def main():
    ATT.mkdir(parents=True, exist_ok=True)
    any_done = False
    for binms in ("25", "50"):
        for suffix, tag in (("", f"bin{binms}_fsexcl"), ("_fsincl", f"bin{binms}_fsincl")):
            path = config.RESULTS_DIR / f"ifi_windows_bin{binms}{suffix}.csv"
            if not path.is_file():
                continue
            df = pd.read_csv(path)
            _fig(df, int(df["bin_ms"].iloc[0]), tag)
            print(f"  wrote HCV1_ifi_windows_{tag}.png")
            any_done = True
    if not any_done:
        print("  (no ifi_windows_*.csv yet — run run_ifi_windows.py first)")
    print("figures ->", ATT)


if __name__ == "__main__":
    main()
