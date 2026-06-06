"""Animals-as-n vs significant-dimensions-as-n figures, with stats, per metric.

For the per-DIMENSION metrics (held-out CC, IFI, optimal lag) the unit of
analysis can be the animal or the significant canonical dimension (the
Gonzalez/Buzsáki unit). This renders, for each such metric, a figure with TWO
rows — ANIMALS-as-n (top) and DIMS-as-n (bottom) — × 8 pairs, showing the metric
across the naive/intermediate/expert epochs with per-unit points + mean ± SEM,
a one-sample test vs 0 at each epoch (stars), and the naive→expert contrast in
each panel title (animals: paired t; dims: rank-sum). Side by side, the figure
makes the pseudoreplication of the dims unit visually obvious (tight SEMs, many
points, smaller p) against the honest animal unit.

(Population-level metrics — Gini, n_sig, mi_sig, rotation — have no per-dimension
analogue and are animals-as-n only; see HCV1_CCA_*_epochs.png.)

Reads results/epoch_metrics.csv (per-animal) + results/epoch_dims.csv (per-dim).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config  # noqa: E402

ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")
EPOCHS = ["naive", "intermediate", "expert"]
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
# (file label, animal column, dim column, y-label)
METRICS = [("cc", "cc1", "peak_cc", "held-out CC"),
           ("ifi", "ifi", "ifi", "IFI (>0: X leads)"),
           ("lag", "optimal_lag", "lag", "optimal lag (bins)")]


CAP_DIMS = 3                  # cap significant dims per animal per epoch (top-N by CC)


def _f(s):
    v = pd.to_numeric(s, errors="coerce").to_numpy()
    return v[np.isfinite(v)]


def _cap(g, cap=CAP_DIMS):
    """Keep the top-`cap` significant dims per (animal, epoch) by held-out CC,
    bounding each animal's contribution and curbing the dims-as-n imbalance."""
    g = g.copy()
    g["peak_cc"] = pd.to_numeric(g["peak_cc"], errors="coerce")
    return g.sort_values("peak_cc", ascending=False).groupby(
        ["animal", "epoch"]).head(cap)


def _star(p):
    return "*" if (np.isfinite(p) and p < 0.05) else ""


def _panel(ax, by_ep, hline0):
    xs = np.arange(3)
    means, sems = [], []
    for i, e in enumerate(EPOCHS):
        v = by_ep[e]
        jit = (np.random.default_rng(i).random(v.size) - 0.5) * 0.25
        ax.scatter(np.full(v.size, i) + jit, v, s=8, color="#34495e", alpha=0.6, zorder=2)
        means.append(np.mean(v) if v.size else np.nan)
        sems.append(np.std(v, ddof=1) / np.sqrt(v.size) if v.size > 1 else 0.0)
        # vs-0 star
        if v.size >= 3:
            p0 = float(stats.ttest_1samp(v, 0.0).pvalue)
            if p0 < 0.05:
                ax.text(i, np.max(v), "*", ha="center", va="bottom", fontsize=11)
    ax.errorbar(xs, means, yerr=sems, color="#c0392b", lw=1.8, capsize=3, zorder=3)
    if hline0:
        ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(["nv", "int", "exp"], fontsize=7)


def fig_metric(em, ed, label, acol, dcol, ylab):
    hline0 = (label != "cc")
    fig, axes = plt.subplots(2, len(PAIRS), figsize=(2.1 * len(PAIRS), 6),
                             squeeze=False)
    for j, pair in enumerate(PAIRS):
        # ANIMALS-as-n (top)
        g = em[em["pair"] == pair]
        by_ep = {e: _f(g[g["epoch"] == e][acol]) for e in EPOCHS}
        ax = axes[0][j]
        _panel(ax, by_ep, hline0)
        piv = g.pivot_table(index="animal", columns="epoch", values=acol, aggfunc="mean")
        ne = ""
        if {"naive", "expert"} <= set(piv.columns):
            d = piv[["naive", "expert"]].dropna()
            if len(d) >= 3:
                p = float(stats.ttest_rel(d["expert"], d["naive"]).pvalue)
                ne = f" e–n p={p:.2g}{_star(p)}"
        ax.set_title(f"{pair}{ne}", fontsize=7)
        if j == 0:
            ax.set_ylabel(f"ANIMALS-as-n\n{ylab}", fontsize=8)
        # DIMS-as-n (bottom) — capped to top-CAP_DIMS sig dims per animal/epoch
        gd = _cap(ed[(ed["pair"] == pair) & (ed["sig"] == 1)])
        by_ep = {e: _f(gd[gd["epoch"] == e][dcol]) for e in EPOCHS}
        ax = axes[1][j]
        _panel(ax, by_ep, hline0)
        ne = ""
        if by_ep["naive"].size >= 3 and by_ep["expert"].size >= 3:
            p = float(stats.mannwhitneyu(by_ep["expert"], by_ep["naive"],
                                         alternative="two-sided").pvalue)
            ne = f" e–n U p={p:.2g}{_star(p)}"
        nd = sum(by_ep[e].size for e in EPOCHS)
        ax.set_title(f"{pair} (nd={nd}){ne}", fontsize=7)
        if j == 0:
            ax.set_ylabel(f"DIMS-as-n\n{ylab}", fontsize=8)
    fig.suptitle(f"{ylab}: animals-as-n (top) vs significant-dims-as-n "
                 f"(bottom, ≤{CAP_DIMS} dims/animal) across epochs — points=units, "
                 "red=mean±SEM, *=vs0 p<0.05; title=naive→expert p", fontsize=11)
    fig.tight_layout()
    ATT.mkdir(parents=True, exist_ok=True)
    fig.savefig(ATT / f"HCV1_CCA_fsexcl_units_{label}.png", dpi=150)
    plt.close(fig)


def main():
    em = pd.read_csv(config.RESULTS_DIR / "epoch_metrics.csv")
    ed = pd.read_csv(config.RESULTS_DIR / "epoch_dims.csv")
    for label, acol, dcol, ylab in METRICS:
        fig_metric(em, ed, label, acol, dcol, ylab)
        print(f"  wrote HCV1_CCA_fsexcl_units_{label}.png")
    print("units-comparison figures ->", ATT)


if __name__ == "__main__":
    main()
