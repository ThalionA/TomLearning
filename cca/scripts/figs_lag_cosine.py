"""Cosine similarity of each canonical component to itself across lags.

Reads results/lag_cosine_bin10{,_fsincl}.csv (run_lag_cosine.py).

Figure 1  HCV1_lag_cosine_<fs>_bin10
          Per pair: |cos| between dim k's canonical weight vector at lag 0 and at lag L,
          significant CCs only, mean +/- SEM across animals, one line per dim.
          SOLID  = same rank at both lags (`cos_same`) — "is CC k still CC k?"
          DASHED = best match to ANY dim at that lag (`cos_best`) — "is the component
                   still there, possibly at another rank?"
          A dashed line well above a solid one is canonical-order SWAPPING, not loss of
          the component. Both sides (X and Y area) are averaged.

Figure 2  HCV1_lag_cosine_swap_<fs>_bin10
          How often the best match is a DIFFERENT dim, as a function of |lag| — the
          rank-swap rate. This is the quantity that silently breaks any analysis
          matching dimensions across fits by bare rank.

Usage: python scripts/figs_lag_cosine.py [fsincl]
"""
from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import RES, config, fs_from_argv
import figstyle

figstyle.apply()

PAIRS = list(config.PAIR_NAMES)
SHOW_DIMS = 4


RANDOM_COS = 0.146     # E|cos| between random unit vectors in R^30 = sqrt(2/(pi*30))


def _prep(df):
    for c in ("cos_same_x", "cos_best_x", "cos_same_y", "cos_best_y",
              "cos_split_x", "cos_split_y", "split_best_x", "split_best_y",
              "best_dim_x", "best_dim_y", "r_zero"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    # average the two areas — the question is about the component, not the side
    df["cos_same"] = df[["cos_same_x", "cos_same_y"]].mean(axis=1)
    df["cos_best"] = df[["cos_best_x", "cos_best_y"]].mean(axis=1)
    # matched split-half: reference is half 1 at lag 0, every lag fitted on half 2, so
    # sampling noise is present at EVERY point including lag 0 and only the lag varies
    df["cos_split"] = df[["cos_split_x", "cos_split_y"]].mean(axis=1)
    df["split_best"] = df[["split_best_x", "split_best_y"]].mean(axis=1)
    df["swapped"] = (((df["best_dim_x"] != df["dim"]).astype(float) +
                      (df["best_dim_y"] != df["dim"]).astype(float)) / 2.0)
    return df


def fig_cosine(df, fs):
    """Matched split-half cosine vs lag, with the no-lag noise floor drawn on.

    The naive version (full-data lag 0 vs full-data lag L) starts at 1.0 by construction
    — lag 0 is compared to itself — so its decline conflates "a different fit" with "a
    different lag". Here BOTH terms are half-data fits, so the lag-0 point already
    carries the sampling noise and the curve isolates the lag.
    """
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    cmap = plt.get_cmap("viridis")
    for ax, pair in zip(axes, PAIRS):
        sub = df[(df["pair"] == pair) & (df["sig"] == 1) & (df["dim"] <= SHOW_DIMS)]
        if sub.empty:
            ax.set_title(f"{pair}\n(no significant CCs)", fontsize=9)
            ax.axis("off"); continue
        for d in range(1, SHOW_DIMS + 1):
            g = sub[sub["dim"] == d]
            if g["animal"].nunique() < 3:
                continue
            col = cmap((d - 1) / max(1, SHOW_DIMS - 1))
            piv = g.pivot_table(index="animal", columns="lag_ms", values="cos_split")
            lags = np.array(sorted(piv.columns), dtype=float)
            M = piv[sorted(piv.columns)].to_numpy(float)
            n = np.sum(np.isfinite(M), axis=0)
            mean = np.nanmean(M, axis=0)
            sd = np.nanstd(M, axis=0, ddof=1)
            sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
            ax.plot(lags, mean, "-", color=col, lw=1.9, zorder=3,
                    label=f"CC{d} (n={piv.shape[0]})")
            ax.fill_between(lags, mean - sem, mean + sem, color=col, alpha=0.13,
                            zorder=2)
            # the no-lag floor for this dim: same measure evaluated at lag 0
            z = g[g["lag_ms"] == 0]["cos_split"].mean()
            if np.isfinite(z):
                ax.axhline(z, color=col, lw=0.8, ls=":", alpha=0.75, zorder=1)
        ax.axhline(RANDOM_COS, color="k", lw=1.1, ls="--", zorder=1)
        ax.text(0.02, RANDOM_COS + 0.02, "random vectors", transform=
                ax.get_yaxis_transform(), fontsize=6, color="k")
        ax.set_ylim(0, 1.02)
        ax.axvline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.set_title(pair, fontsize=10)
        ax.set_xlabel("lag (ms)", fontsize=8)
        ax.set_ylabel("|cos|  (split-half matched)", fontsize=8)
        ax.legend(fontsize=6, loc="upper right", framealpha=0.6)
    fig.suptitle(f"Is it the same canonical component at every lag? — {fs} | BOTH terms "
                 f"are half-data fits, so the dotted line at lag 0 is pure sampling "
                 f"noise, not a lag effect", fontsize=10.5)
    figstyle.save(fig, f"HCV1_lag_cosine_{fs}_bin10.png")
    print(f"wrote HCV1_lag_cosine_{fs}_bin10.png")


def fig_swap(df, fs):
    fig, axes = figstyle.grid(2, ncols=2, panel=(6.4, 4.2))
    sub = df[(df["sig"] == 1) & (df["dim"] <= SHOW_DIMS)].copy()
    sub["abs_lag"] = sub["lag_ms"].abs()

    ax = axes[0]
    for d in range(1, SHOW_DIMS + 1):
        g = sub[sub["dim"] == d]
        if g.empty:
            continue
        m = g.groupby("abs_lag")["swapped"].mean()
        ax.plot(m.index, m.to_numpy(), "-o", ms=3, label=f"CC{d}")
    ax.set_xlabel("|lag| (ms)", fontsize=9)
    ax.set_ylabel("fraction where the best match is a DIFFERENT dim", fontsize=9)
    ax.set_ylim(0, 1)
    ax.set_title("rank-swap rate vs lag", fontsize=10)
    ax.legend(fontsize=7, framealpha=0.6)

    ax = axes[1]
    m = sub.groupby("pair")[["cos_same", "cos_best"]].mean().reindex(
        [p for p in PAIRS if p in sub["pair"].unique()])
    y = np.arange(len(m))
    ax.barh(y - 0.2, m["cos_same"].to_numpy(), height=0.4, color="#c0392b",
            label="same rank")
    ax.barh(y + 0.2, m["cos_best"].to_numpy(), height=0.4, color="#95a5a6",
            label="best match")
    ax.set_yticks(y); ax.set_yticklabels(m.index, fontsize=8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("mean |cos| across all lags", fontsize=9)
    ax.set_title("same-rank vs best-match, per pair", fontsize=10)
    ax.legend(fontsize=7, framealpha=0.6)
    fig.suptitle(f"Canonical-order swapping across lag — {fs} | the gap between the two "
                 f"bars is what bare-rank matching gets wrong", fontsize=11)
    figstyle.save(fig, f"HCV1_lag_cosine_swap_{fs}_bin10.png")
    print(f"wrote HCV1_lag_cosine_swap_{fs}_bin10.png")


def main():
    suf = fs_from_argv()                       # "" (FS-excluded) or "_fsincl"
    fs = "fsincl" if suf else "fsexcl"
    src = RES / f"lag_cosine_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/run_lag_cosine.py first")
    df = _prep(pd.read_csv(src))
    fig_cosine(df, fs)
    fig_swap(df, fs)


if __name__ == "__main__":
    main()
