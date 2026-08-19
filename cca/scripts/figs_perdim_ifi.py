"""Per-dimension IFI figures — meeting 2026-07-28 item 1.

Reads results/perdim_ifi_bin10{,_fsincl}.csv (written by analyze_perdim_ifi.py).

Figure 1  HCV1_perdim_ifi_<fs>_bin10   — IFI vs canonical RANK, 8-pair grid.
          Thin grey lines = individual animals; heavy line = mean +/- SEM across
          animals; significant dims only. y = 0 is "no directional asymmetry".
Figure 2  HCV1_perdim_ifi_cc1vstail_<fs>_bin10 — per animal, CC1's IFI against the
          mean IFI of its significant tail, paired (the test in the tables).

Positive IFI = the FIRST-named area leads.

Usage: python scripts/figs_perdim_ifi.py [fsincl]
"""
from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from _common import RES, config, fs_from_argv
import figstyle
from tom_cca import perdim_ifi

figstyle.apply()

PAIRS = list(config.PAIR_NAMES)
C_MEAN, C_ANIM, C_CC1 = "#c0392b", "#95a5a6", "#2c6fbb"
MIN_FRAC_ANIMALS = 0.5      # a rank joins the connected mean only if >=half contribute


def _consecutive_runs(vals):
    """Split a sorted integer array into runs of consecutive values."""
    runs, cur = [], [vals[0]]
    for v in vals[1:]:
        if v == cur[-1] + 1:
            cur.append(v)
        else:
            runs.append(cur); cur = [v]
    runs.append(cur)
    return runs


def fig_rank(df: pd.DataFrame, fs: str):
    """IFI as a function of canonical rank, per pair.

    Significance is sparse and non-contiguous at high rank (a rank-19 dim may be
    significant in 2 of 9 animals), so a connected line across all ranks would imply
    a smooth function of rank that the data do not support — and would average over a
    *selected* subset of animals. The heavy line is therefore drawn only over
    consecutive ranks carried by at least half the pair's animals; sparser ranks appear
    as unconnected points sized by how many animals contribute.
    """
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[(df["pair"] == pair) & (df["sig"] == 1)]
        if sub.empty:
            ax.set_title(f"{pair}\n(no significant dims)", fontsize=9)
            ax.axis("off"); continue
        x_lead = pair.split("-")[0]
        n_animals_pair = df[df["pair"] == pair]["animal"].nunique()
        piv = sub.pivot_table(index="animal", columns="dim", values="ifi")
        dims = np.array(sorted(piv.columns), dtype=int)
        M = piv[sorted(piv.columns)].to_numpy(float)
        n = np.sum(np.isfinite(M), axis=0)
        mean = np.nanmean(M, axis=0)
        sd = np.nanstd(M, axis=0, ddof=1)
        sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
        # every animal's per-rank IFI as unconnected dots (no line across rank gaps)
        for row in M:
            ok = np.isfinite(row)
            ax.plot(dims[ok], row[ok], ".", color=C_ANIM, ms=3.5, alpha=0.65,
                    zorder=1)
        # sparse ranks: point + errorbar only, sized by contributing animals
        ax.errorbar(dims, mean, yerr=sem, fmt="o", color=C_MEAN, ms=3, lw=1.0,
                    alpha=0.45, capsize=2, zorder=2)
        dense = dims[n >= max(2, np.ceil(MIN_FRAC_ANIMALS * n_animals_pair))]
        if dense.size:
            for run in _consecutive_runs(dense):
                if len(run) < 2:
                    continue
                m = np.isin(dims, run)
                ax.plot(dims[m], mean[m], "-o", color=C_MEAN, lw=2.0, ms=5, zorder=4)
                ax.fill_between(dims[m], (mean - sem)[m], (mean + sem)[m],
                                color=C_MEAN, alpha=0.2, zorder=3)
        ax.axhline(0.0, color="k", lw=0.8, ls=":", zorder=0)
        ax.set_title(f"{pair}  ({n_animals_pair} animals)", fontsize=10)
        ax.set_xlabel("canonical rank (dim)", fontsize=8)
        ax.set_ylabel(f"IFI   +: {x_lead} leads", fontsize=8)
        ax.set_xlim(0.3, dims.max() + 0.7)
        ax.plot([], [], "-o", color=C_MEAN, lw=2.0, ms=5,
                label=f"mean ± SEM (≥{int(np.ceil(MIN_FRAC_ANIMALS * n_animals_pair))}"
                      " animals)")
        ax.plot([], [], "o", color=C_MEAN, ms=3, alpha=0.45, label="sparser ranks")
        ax.legend(fontsize=6, loc="upper right", framealpha=0.6)
    fig.suptitle(f"Information Flow Index by canonical rank — {fs}, 10 ms smoothed | "
                 "significant dims only; +IFI = first area leads", fontsize=12)
    figstyle.save(fig, f"HCV1_perdim_ifi_{fs}_bin10.png")
    print(f"wrote HCV1_perdim_ifi_{fs}_bin10.png")


def fig_cc1_vs_tail(df: pd.DataFrame, fs: str):
    """Paired per-animal CC1 IFI vs the mean IFI of its significant tail."""
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[df["pair"] == pair].sort_values(["animal", "dim"])
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        x_lead = pair.split("-")[0]
        pts = []
        for animal, grp in sub.groupby("animal", sort=True):
            grp = grp.sort_values("dim")
            cc1, rest = perdim_ifi.rank_split(grp["ifi"].to_numpy(float),
                                              grp["sig"].to_numpy(int))
            if np.isfinite(cc1) and np.isfinite(rest):
                pts.append((cc1, rest))
        if not pts:
            ax.set_title(f"{pair}\n(no significant tail)", fontsize=9)
            ax.axis("off"); continue
        for cc1, rest in pts:
            ax.plot([0, 1], [cc1, rest], "-", color=C_ANIM, lw=0.9, alpha=0.7,
                    zorder=1)
            ax.plot(0, cc1, "o", color=C_CC1, ms=5, zorder=3)
            ax.plot(1, rest, "o", color=C_MEAN, ms=5, zorder=3)
        arr = np.array(pts)
        ax.plot([0, 1], arr.mean(axis=0), "-", color="k", lw=2.2, zorder=4)
        ax.axhline(0.0, color="k", lw=0.8, ls=":", zorder=0)
        ax.set_xlim(-0.35, 1.35)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["CC₁", "sig tail"], fontsize=8)
        ax.set_title(f"{pair}  (n={len(pts)})", fontsize=10)
        ax.set_ylabel(f"IFI   +: {x_lead} leads", fontsize=8)
    fig.suptitle(f"CC₁ vs the rest of the significant spectrum — {fs}, 10 ms smoothed | "
                 "one line per animal", fontsize=12)
    figstyle.save(fig, f"HCV1_perdim_ifi_cc1vstail_{fs}_bin10.png")
    print(f"wrote HCV1_perdim_ifi_cc1vstail_{fs}_bin10.png")


def main():
    suf = fs_from_argv()                       # "" (FS-excluded) or "_fsincl"
    fs = "fsincl" if suf else "fsexcl"
    src = RES / f"perdim_ifi_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/analyze_perdim_ifi.py first")
    df = pd.read_csv(src)
    fig_rank(df, fs)
    fig_cc1_vs_tail(df, fs)


if __name__ == "__main__":
    main()
