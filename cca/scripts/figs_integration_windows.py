"""Integration windows and why they must not be read naively — meeting item 6.

Reads results/fixed_subspace_bin10{,_fsincl}.csv and the reduced
results/fixed_subspace_epoch_bin10{,_fsincl}.csv.

The raw half-max widths look like a clean hierarchy of integration timescales
(intra-hippocampal ~20-30 ms vs cortical ~180-260 ms). That reading is mostly wrong: the
frozen canonical component's lag curve RINGS at theta, so the half-max region is truncated
by the first trough and the "width" is a half-period of the rhythm. These two figures make
that visible instead of leaving it in a table.

Figure 1  HCV1_integration_window_anatomy_<fs>_bin10
          Per pair, the across-animal mean lag curve with (a) the half-max level, (b) the
          CONTIGUOUS half-max region shaded — this is what curve_half_width measures — and
          (c) the detected secondary peak marked with its implied frequency. Where a side
          peak exists, the shaded region is bounded by the ring, not by a decay.

Figure 2  HCV1_integration_window_summary_<fs>_bin10
          One panel: mean half-max width against the fraction of curves that ring, per pair.
          Pairs at high ringing fraction (upper region) have widths that are half-periods;
          only the low-ringing pairs have widths readable as integration windows.

Usage: PYTHONPATH=src python scripts/figs_integration_windows.py [fsincl]
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
from tom_cca import fixed_subspace  # noqa: E402

figstyle.apply()
ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"
RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
C_CURVE, C_HALF, C_RING = "#c0392b", "#2c6fbb", "#8e44ad"
RING_RATIO = 0.5          # same gate analyze_fixed_subspace uses
MAJORITY = 0.5


def _half_max_span(lags, r):
    """(lo, hi, half) of the contiguous >= half-max region around the peak — mirrors
    fixed_subspace.curve_half_width so the shading shows exactly what it measures."""
    ok = np.isfinite(lags) & np.isfinite(r)
    lags, r = lags[ok], r[ok]
    order = np.argsort(lags)
    lags, r = lags[order], r[order]
    pk = int(np.argmax(r))
    if r[pk] <= 0:
        return None
    half = r[pk] / 2.0
    lo = hi = pk
    while lo - 1 >= 0 and r[lo - 1] >= half:
        lo -= 1
    while hi + 1 < r.size and r[hi + 1] >= half:
        hi += 1
    return float(lags[lo]), float(lags[hi]), float(half)


def fig_anatomy(df: pd.DataFrame, red: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[(df["pair"] == pair) & (df["dim"] == 1)]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        x_lead = pair.split("-")[0]
        piv = sub.pivot_table(index=["animal", "epoch"], columns="lag_ms", values="r")
        lags = np.array(sorted(piv.columns), dtype=float)
        M = piv[sorted(piv.columns)].to_numpy(float)
        n = np.sum(np.isfinite(M), axis=0)
        mean = np.nanmean(M, axis=0)
        sd = np.nanstd(M, axis=0, ddof=1)
        sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
        ax.plot(lags, mean, "-", color=C_CURVE, lw=1.8, zorder=4)
        ax.fill_between(lags, mean - sem, mean + sem, color=C_CURVE, alpha=0.16,
                        zorder=3)
        span = _half_max_span(lags, mean)
        if span is not None:
            lo, hi, half = span
            ax.axhline(half, color=C_HALF, lw=1.0, ls="--", zorder=2)
            ax.axvspan(lo, hi, color=C_HALF, alpha=0.16, zorder=1)
        side_ms, side_ratio = fixed_subspace.side_peak(lags, mean)
        if np.isfinite(side_ms) and side_ms > 0:
            pk_lag = lags[int(np.nanargmax(mean))]
            for sgn in (+1, -1):
                xpos = pk_lag + sgn * side_ms
                if lags.min() <= xpos <= lags.max():
                    ax.axvline(xpos, color=C_RING, lw=1.1, ls=":", zorder=2)
            freq = 1000.0 / side_ms
            ring_txt = f"ring {side_ms:.0f} ms = {freq:.1f} Hz"
        else:
            ring_txt = "no ring"
        g = red[red["pair"] == pair]
        w = g["width_ms"].mean() if len(g) else np.nan
        frac = ((g["side_peak_ratio"] > RING_RATIO) & g["side_peak_ms"].notna()).mean() \
            if len(g) else np.nan
        verdict = "half-period" if (np.isfinite(frac) and frac >= MAJORITY) else "window"
        shaded = f"{span[1] - span[0]:.0f}" if span is not None else "—"
        # Two different widths, deliberately distinguished: the shaded band is the
        # half-max span of THIS displayed mean curve, whereas `width` is the mean of the
        # per-animal widths that the statistics use. Averaging curves with different
        # peak lags narrows the mean curve, so they need not agree.
        ax.set_title(f"{pair}   per-animal width {w:.0f} ms → {verdict}\n"
                     f"{ring_txt} · rings in {frac:.0%} of curves · "
                     f"band on mean curve {shaded} ms", fontsize=8)
        ax.set_xlabel(f"lag (ms)   +: {x_lead} leads", fontsize=8)
        ax.set_ylabel("r (frozen CC₁)", fontsize=8)
        ax.axvline(0, color="k", lw=0.6, ls=":", zorder=0)
    fig.suptitle(
        f"What the integration window actually measures — {fs}, 10 ms smoothed | blue = "
        "half-max level and the contiguous region `curve_half_width` returns; purple = "
        "detected ring", fontsize=11.5)
    figstyle.save(fig, ATT / f"HCV1_integration_window_anatomy_{fs}_bin10.png")
    print(f"wrote HCV1_integration_window_anatomy_{fs}_bin10.png")


def fig_summary(red: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(1, ncols=1, panel=(7.2, 5.2))
    ax = axes[0]
    for pair in PAIRS:
        g = red[red["pair"] == pair]
        if g.empty:
            continue
        frac = ((g["side_peak_ratio"] > RING_RATIO) & g["side_peak_ms"].notna()).mean()
        w = g["width_ms"].mean()
        w_sem = g["width_ms"].std(ddof=1) / np.sqrt(g["width_ms"].notna().sum())
        intra = pair in ("CA1-CA3", "CA1-DG", "CA3-DG")
        colour = "#c0392b" if intra else "#2c6fbb"
        ax.errorbar(frac, w, yerr=w_sem, fmt="o", color=colour, ms=9, capsize=3,
                    zorder=3)
        ax.annotate(pair, (frac, w), xytext=(6, 6), textcoords="offset points",
                    fontsize=8.5, zorder=4)
    ax.axvline(MAJORITY, color="k", lw=1.0, ls="--", alpha=0.7, zorder=1)
    ax.text(MAJORITY + 0.02, ax.get_ylim()[1] * 0.96,
            "majority of curves ring →\nwidth is a theta half-period",
            fontsize=8.5, va="top", color="#555555")
    ax.text(MAJORITY - 0.02, ax.get_ylim()[1] * 0.96,
            "← width readable as\nan integration window",
            fontsize=8.5, va="top", ha="right", color="#555555")
    ax.plot([], [], "o", color="#c0392b", ms=8, label="intra-hippocampal")
    ax.plot([], [], "o", color="#2c6fbb", ms=8, label="hippocampal–cortical / cortical")
    ax.set_xlabel("fraction of curves with a secondary peak (ringing)", fontsize=9)
    ax.set_ylabel("mean half-max width (ms)", fontsize=9)
    # Deliberately NOT claiming a correlation. Across pairs Spearman rho = -0.55,
    # p = 0.154 (n = 8) FS-excluded and -0.50, p = 0.207 FS-included — not significant,
    # and pairs share areas and animals so a cross-pair correlation is pseudoreplicated
    # regardless. What the panel shows is WHICH SIDE of the majority line each pair falls
    # on, which is a per-pair classification, not a trend.
    ax.set_title(f"Integration window vs ringing — {fs}, 10 ms smoothed\n"
                 "only pairs LEFT of the line have a width readable as an "
                 "integration window", fontsize=10.5)
    ax.legend(fontsize=8, loc="center right", framealpha=0.7)
    figstyle.save(fig, ATT / f"HCV1_integration_window_summary_{fs}_bin10.png")
    print(f"wrote HCV1_integration_window_summary_{fs}_bin10.png")


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    src = RES / f"fixed_subspace_bin10{suf}.csv"
    red_src = RES / f"fixed_subspace_epoch_bin10{suf}.csv"
    if not src.exists() or not red_src.exists():
        sys.exit("run run_fixed_subspace.py and analyze_fixed_subspace.py first")
    df = pd.read_csv(src)
    df["r"] = pd.to_numeric(df["r"], errors="coerce")
    red = pd.read_csv(red_src)
    fig_anatomy(df, red, fs)
    fig_summary(red, fs)


if __name__ == "__main__":
    main()
