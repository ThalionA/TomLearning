"""Held-out lagged canonical-correlation curves per pair — the cross-correlation
profile behind the IFI (report R2). Reads results/lag_curves_bin10{,_fsincl}.csv.

Per area-pair (8-pair grid) two curves are overlaid, both mean ± shaded SEM:
  * solid  — CC1 (dominant dim), one curve per ANIMAL, averaged across animals (n=animals)
  * dashed — every SIGNIFICANT canonical dim's curve, pooled across animals & dims
             (n=subspaces, the Buzsáki-style dims-as-n unit)
x = bin lag in ms; POSITIVE lag = the FIRST-named area leads (X leads Y). The vertical
line marks lag 0; the horizontal line marks CC = 0.

Usage: python scripts/figs_lag_curves.py [fsincl]
       (default FS-excl; pass 'fsincl' for the FS-included variant)
"""
from __future__ import annotations


import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from _common import RES, config, fs_from_argv
import figstyle

figstyle.apply()

PAIRS = list(config.PAIR_NAMES)
C_ANIM, C_DIM = "#c0392b", "#2c6fbb"


def _mean_sem(piv):
    """piv: rows = curves, columns = lag_ms. Returns (lags, mean, sem) over rows,
    NaN-aware, with SEM from the count of finite curves at each lag."""
    lags = np.array(sorted(piv.columns))
    M = piv[lags].to_numpy(float)
    n = np.sum(np.isfinite(M), axis=0)
    mean = np.nanmean(M, axis=0)
    sd = np.nanstd(M, axis=0, ddof=1)
    sem = np.where(n > 1, sd / np.sqrt(n), 0.0)
    return lags, mean, sem


def main():
    suf = fs_from_argv()                       # "" (FS-excluded) or "_fsincl"
    fs = "fsincl" if suf else "fsexcl"
    df = pd.read_csv(RES / f"lag_curves_bin10{suf}.csv")
    df["cc"] = pd.to_numeric(df["cc"], errors="coerce")

    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[df["pair"] == pair]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        x_lead, y_lead = pair.split("-")
        ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="#b4b2a9", lw=0.8, ls=":")
        # n = ANIMALS: dominant-dim (CC1) curve, one per animal
        a1 = sub[sub["dim"] == 1].pivot_table(index="animal", columns="lag_ms", values="cc")
        if not a1.empty:
            lags, m, se = _mean_sem(a1)
            ax.plot(lags, m, "-", color=C_ANIM, lw=1.8, zorder=3,
                    label=f"CC₁ · n={a1.shape[0]} animals")
            ax.fill_between(lags, m - se, m + se, color=C_ANIM, alpha=0.18, zorder=1)
        # n = SUBSPACES: every significant canonical dim, pooled across animals & dims
        sd = sub[sub["sig"] == 1].pivot_table(index=["animal", "dim"],
                                              columns="lag_ms", values="cc")
        if not sd.empty:
            lags, m, se = _mean_sem(sd)
            ax.plot(lags, m, "--", color=C_DIM, lw=1.6, zorder=3,
                    label=f"sig dims · n={sd.shape[0]}")
            ax.fill_between(lags, m - se, m + se, color=C_DIM, alpha=0.15, zorder=1)
        ax.set_title(pair, fontsize=10)
        ax.set_xlabel("lag (ms)   +: %s leads" % x_lead, fontsize=8)
        ax.set_ylabel("held-out CC", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper right", framealpha=0.6)
    fig.suptitle(f"Lagged canonical correlation (held-out, segment-aware) — {fs}, "
                 f"10 ms smoothed | +lag = first area leads", fontsize=12)
    figstyle.save(fig, f"HCV1_lagcc_{fs}_bin10.png")
    print(f"wrote HCV1_lagcc_{fs}_bin10.png")


if __name__ == "__main__":
    main()
