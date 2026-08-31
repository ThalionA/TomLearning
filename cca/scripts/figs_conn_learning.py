"""Figures: naive-epoch communication vs learning speed (exploratory, null).

Reads results/conn_learning_bin10*.csv and epoch_metrics_bin10*.csv (written
by analyze_conn_learning.py / run_epochs.py). Outcome = learning point LP
(trials to expert criterion; LOWER = faster). Predictors measured in the
naive epoch (first 10 running trials) — before learning.

Output per FS: HCV1_connlearning_<fs>_bin10 —
  A  Spearman rho forest per pair: naive cc1 (strength) and naive IFI
     (direction) vs LP; grey = |n| < 5 pairs excluded upstream
  B  scatter: global naive coupling (within-pair z of cc1, animal mean) vs LP
  C  scatter: CA1-DG naive cc1 vs LP — the one repeatable (n.s.) lean

Usage: PYTHONPATH=src python scripts/figs_conn_learning.py
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
from analyze_conn_learning import PAIRS  # noqa: E402
from figs_contrib_reliability import COL_X, COL_Y  # noqa: E402
from scipy import stats  # noqa: E402
from tom_cca import config, dataio  # noqa: E402
import dataclasses  # noqa: E402

figstyle.apply()
RES = Path(__file__).resolve().parents[1] / "results"


def naive_with_lp(tag):
    em = pd.read_csv(RES / f"epoch_metrics_bin10{tag}.csv")
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=10,
                              gaussian_sd_ms=2.5)
    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(
        animals, cfg, behaviour_lookup=dataio.load_learning_points())
    lp = {a: e.lp for a, e in entries.items()
          if e.role == "learner" and e.lp is not None}
    naive = em[(em["epoch"] == "naive") & (em["learner"] == 1)].copy()
    naive["lp"] = naive["animal"].map(lp)
    return naive[np.isfinite(naive["lp"])]


def rho_forest(ax, cl):
    y = 0
    ticks, labels = [], []
    for metric, col, mlab in (("cc1", COL_X, "strength (cc1)"),
                              ("ifi", COL_Y, "direction (IFI)")):
        for pair in PAIRS:
            r = cl[(cl["metric"] == metric) & (cl["pair"] == pair)]
            if r.empty:
                continue
            r = r.iloc[0]
            ticks.append(y); labels.append(f"{mlab[:4]} · {pair}")
            ax.plot(r["rho"], y, "o", ms=6, color=col, mec="black", mew=0.6)
            if r["p"] < 0.05:
                ax.text(0.93, y, "*", transform=ax.get_yaxis_transform(),
                        ha="center", va="center", fontsize=15,
                        fontweight="bold", color=figstyle.STAR_COLOR)
            ax.text(1.02, y, f"n={int(r['n'])}",
                    transform=ax.get_yaxis_transform(), va="center",
                    fontsize=7, color="0.45")
            y += 1
        y += 0.6
    ax.axvline(0, color="0.6", lw=0.8)
    ax.set_yticks(ticks, labels, fontsize=7)
    ax.set_ylim(y - 1.1, -0.9)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("Spearman ρ vs LP\n(−ve = stronger → faster)", fontsize=8)
    ax.set_title("A  naive-epoch predictors of learning speed", fontsize=9)


def scatter(ax, x, y_lp, ids, xlab, title):
    ax.scatter(x, y_lp, s=28, color=COL_X, alpha=0.8, zorder=3)
    for xi, yi, an in zip(x, y_lp, ids):
        ax.annotate(str(int(an)), (xi, yi), xytext=(3, 3),
                    textcoords="offset points", fontsize=6, color="0.4")
    r = stats.spearmanr(x, y_lp)
    ax.set_xlabel(xlab, fontsize=8)
    ax.set_ylabel("learning point (trials; lower = faster)", fontsize=8)
    ax.set_title(f"{title}\nρ = {r.statistic:+.2f}, p = {r.pvalue:.3f}, "
                 f"n = {len(x)}", fontsize=8)


def main():
    for fs, tag in [("fsexcl", ""), ("fsincl", "_fsincl")]:
        cl = pd.read_csv(RES / f"conn_learning_bin10{tag}.csv")
        naive = naive_with_lp(tag)
        fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8),
                                 constrained_layout=True)
        rho_forest(axes[0], cl)
        z = naive.copy()
        z["cc1_z"] = z.groupby("pair")["cc1"].transform(
            lambda v: (v - v.mean()) / v.std(ddof=1))
        ga = z.groupby("animal").agg(cc1_z=("cc1_z", "mean"), lp=("lp", "first"))
        scatter(axes[1], ga["cc1_z"], ga["lp"], ga.index,
                "global naive coupling (mean within-pair z of cc1)",
                "B  global coupling vs learning speed")
        g = naive[naive["pair"] == "CA1-DG"][["animal", "cc1", "lp"]].dropna()
        scatter(axes[2], g["cc1"], g["lp"], g["animal"],
                "CA1-DG naive cc1 (held-out)",
                "C  the one repeatable lean (n.s.)")
        fig.suptitle("Does early communication predict how fast animals learn? "
                     f"— no survivor ({'FS excluded' if fs == 'fsexcl' else 'FS included'})",
                     fontsize=11)
        figstyle.save(fig, f"HCV1_connlearning_{fs}_bin10")
        print(f"{fs}: 1 figure written")


if __name__ == "__main__":
    main()
