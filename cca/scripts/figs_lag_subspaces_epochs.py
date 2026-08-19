"""FF/FB evolution across learning — meeting 2026-08-07, item 4 (and item 2's second half).

Reads results/lag_subspaces_bin10_epochs{,_fsincl}.csv (run_lag_subspaces.py --epochs).

Figure 1  HCV1_lagsubspace_evolution_<fs>_bin10
          Per pair, the held-out CC₁ of the FEEDFORWARD fit (+50 ms, first area leads) and the
          FEEDBACK fit (−50 ms) across naive → intermediate → expert, mean ± SEM across animals.
          The gap between the two lines is the directional asymmetry; the question item 4 asks is
          whether that gap MOVES. Sub-title carries the expert−naive asymmetry change and its p.

Figure 2  HCV1_lagsubspace_evolution_curves_<fs>_bin10
          The same data as full lag curves — held-out CC₁ vs lag, one curve per epoch. This is the
          "separate lagged curves ... also for naive/exp" plot. The FF and FB readouts above are
          the ±50 ms points of these curves.

⚠ Read both against the session-level gate: the FF and FB subspaces are NOT separable (their
principal angle never clears the split-half noise floor, 0/8 pairs, both FS). So the gap plotted
here is an asymmetry in coupling strength WITHIN one subspace, not two streams.

Usage: python scripts/figs_lag_subspaces_epochs.py [fsincl]
"""
from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from _common import EPOCHS, RES, TEMPORAL, config, fs_from_argv
import figstyle

figstyle.apply()

PAIRS = list(config.PAIR_NAMES)
EPOCH_COLOUR = {"naive": "#2c6fbb", "intermediate": "#95a5a6", "expert": "#c0392b"}
C_FF, C_FB = "#c0392b", "#2c6fbb"
TAU_MS = TEMPORAL.label_w_bins * TEMPORAL.bin_ms   # 50 ms


def _mean_sem(vals):
    v = np.asarray([x for x in vals if np.isfinite(x)], dtype=float)
    if v.size == 0:
        return np.nan, np.nan
    return float(v.mean()), (float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1
                             else np.nan)


def fig_evolution(df: pd.DataFrame, evo: pd.DataFrame, fs: str):
    """FF and FB coupling strength across the three learning epochs."""
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    x = np.arange(len(EPOCHS))
    for ax, pair in zip(axes, PAIRS):
        sub = df[df["pair"] == pair]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        x_lead, y_lead = pair.split("-")
        n_seen = 0
        for lag, colour, label in [(TAU_MS, C_FF, f"FF (+{TAU_MS} ms) {x_lead}→{y_lead}"),
                                   (-TAU_MS, C_FB, f"FB (−{TAU_MS} ms) {y_lead}→{x_lead}")]:
            means, sems = [], []
            for epoch in EPOCHS:
                v = sub[(sub["epoch"] == epoch) & (sub["lag_ms"] == lag)]["cc1"]
                m, s = _mean_sem(v.to_numpy(float))
                means.append(m); sems.append(s)
                n_seen = max(n_seen, int(v.notna().sum()))
            ax.errorbar(x, means, yerr=sems, fmt="-o", color=colour, lw=1.9, ms=5,
                        capsize=3, zorder=3, label=label)
        r = evo[evo["pair"] == pair]
        if len(r):
            r = r.iloc[0]
            sub_t = (f"Δasym (exp−naive) = {r['d_asym']:+.3f}, p = {r['p_asym']:.2f}"
                     f"  (n={r['n_asym']:.0f})")
        else:
            sub_t = ""
        ax.set_title(f"{pair}\n{sub_t}", fontsize=8.5)
        ax.set_xticks(x); ax.set_xticklabels(EPOCHS, fontsize=8)
        ax.set_ylabel("held-out CC₁", fontsize=8)
        ax.legend(fontsize=6, loc="best", framealpha=0.6)
    fig.suptitle(
        f"Feedforward vs feedback coupling across learning — {fs}, 10 ms smoothed | "
        "the GAP is the directional asymmetry; item 4 asks whether it moves "
        "(it does not: 0/8 pairs)", fontsize=11.5)
    figstyle.save(fig, f"HCV1_lagsubspace_evolution_{fs}_bin10.png")
    print(f"wrote HCV1_lagsubspace_evolution_{fs}_bin10.png")


def fig_epoch_curves(df: pd.DataFrame, fs: str):
    """Full held-out CC₁ lag curves, one per learning epoch."""
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = df[df["pair"] == pair]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        x_lead = pair.split("-")[0]
        for epoch in EPOCHS:
            e = sub[sub["epoch"] == epoch]
            if e.empty:
                continue
            piv = e.pivot_table(index="animal", columns="lag_ms", values="cc1")
            lags = np.array(sorted(piv.columns), dtype=float)
            M = piv[sorted(piv.columns)].to_numpy(float)
            n = np.sum(np.isfinite(M), axis=0)
            mean = np.nanmean(M, axis=0)
            sd = np.nanstd(M, axis=0, ddof=1)
            sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
            ax.plot(lags, mean, "-", color=EPOCH_COLOUR[epoch], lw=1.7, zorder=3,
                    label=f"{epoch} (n={piv.shape[0]})")
            ax.fill_between(lags, mean - sem, mean + sem,
                            color=EPOCH_COLOUR[epoch], alpha=0.16, zorder=2)
        for v, ls in [(0, ":"), (TAU_MS, "--"), (-TAU_MS, "--")]:
            ax.axvline(v, color="k", lw=0.6, ls=ls, alpha=0.5, zorder=0)
        ax.set_title(pair, fontsize=10)
        ax.set_xlabel(f"lag (ms)   +: {x_lead} leads", fontsize=8)
        ax.set_ylabel("held-out CC₁", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper right", framealpha=0.6)
    fig.suptitle(f"Held-out lagged CC₁ curves by learning epoch — {fs}, 10 ms smoothed | "
                 f"dashed lines mark the ±{TAU_MS} ms FF/FB readouts", fontsize=12)
    figstyle.save(fig, f"HCV1_lagsubspace_evolution_curves_{fs}_bin10.png")
    print(f"wrote HCV1_lagsubspace_evolution_curves_{fs}_bin10.png")


def main():
    suf = fs_from_argv()                       # "" (FS-excluded) or "_fsincl"
    fs = "fsincl" if suf else "fsexcl"
    src = RES / f"lag_subspaces_bin10_epochs{suf}.csv"
    evo_src = RES / f"lag_subspaces_evolution_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run run_lag_subspaces.py --epochs first")
    df = pd.read_csv(src)
    df["cc1"] = pd.to_numeric(df["cc1"], errors="coerce")
    evo = pd.read_csv(evo_src) if evo_src.exists() else pd.DataFrame(columns=["pair"])
    fig_evolution(df, evo, fs)
    fig_epoch_curves(df, fs)


if __name__ == "__main__":
    main()
