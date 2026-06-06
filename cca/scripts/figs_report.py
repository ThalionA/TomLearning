"""Generate report figures from the committed full-suite CSVs into the vault.

Reads results/trajectory_windows.csv and results/transition_uncued_cued.csv
(no re-fitting) and writes labelled PNGs to the ResearchVault attachments folder
for embedding in the Hippocampus-V1 report. Every figure has axis labels with
units, a title, and a legend where multiple series are shown.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, paired_stats, trajectory  # noqa: E402

ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
RES = config.RESULTS_DIR


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def per_animal_slopes(df, metric, axis):
    out = defaultdict(list)
    for (an, pr), g in df.groupby(["animal", "pair"]):
        x = _num(g[axis]).to_numpy(); y = _num(g[metric]).to_numpy()
        if np.sum(np.isfinite(x) & np.isfinite(y)) >= 4:
            s, _ = trajectory.linear_slope(x, y)
            if np.isfinite(s):
                out[pr].append(s)
    return out


def fig_levels(df):
    learn = df[df["learner"] == 1]
    metrics = [("cc1", "held-out CC$_1$"), ("n_sig", "# significant dims"),
               ("mi_sig", "MI over sig. dims (nats)")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    for ax, (m, lab) in zip(axes, metrics):
        vals = [np.nanmean(_num(learn[learn["pair"] == p][m])) for p in PAIRS]
        order = np.argsort(vals)[::-1]
        ax.bar([PAIRS[i] for i in order], [vals[i] for i in order],
               color="#34495e")
        ax.set_ylabel(lab); ax.set_xticklabels([PAIRS[i] for i in order],
                                               rotation=45, ha="right", fontsize=8)
    fig.suptitle("Communication-subspace levels by area pair "
                 "(continuous running, pCCA, learners)")
    fig.tight_layout(); fig.savefig(ATT / "HCV1_CCA_levels.png", dpi=150)
    plt.close(fig)


def fig_gini(df):
    learn = df[df["learner"] == 1]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    # left: CA1-RSC Gini trajectories per animal
    ax = axes[0]
    for an, g in learn[learn["pair"] == "CA1-RSC"].groupby("animal"):
        x = _num(g["trial_frac"]); y = _num(g["gini_x"])
        ax.plot(x, y, "-o", ms=3, alpha=0.5, color="#c0392b")
    ax.set_xlabel("trial fraction (0 = first, 1 = last)")
    ax.set_ylabel("Gini (CA1 subspace-weight sparsity)")
    ax.set_title("CA1–RSC: subspace de-sparsifies with learning")
    # right: per-pair median Gini slope vs trial_frac, sign-test stars
    ax = axes[1]
    sl = per_animal_slopes(learn, "gini_x", "trial_frac")
    pairs = [p for p in PAIRS if len(sl.get(p, [])) >= 3]
    meds = [np.median(sl[p]) for p in pairs]
    ps = [paired_stats.wilcoxon_signed(sl[p])[3] for p in pairs]
    colours = ["#c0392b" if (np.isfinite(pp) and pp < 0.05) else "#95a5a6"
               for pp in ps]
    ax.barh(pairs, meds, color=colours)
    ax.axvline(0, color="k", lw=0.6)
    for i, (p, pp) in enumerate(zip(pairs, ps)):
        if np.isfinite(pp) and pp < 0.05:
            ax.text(meds[i], i, " *", va="center", fontsize=14)
    ax.set_xlabel("median per-animal slope  d(Gini)/d(trial fraction)")
    ax.set_title("Gini slope with learning (red = sign-test p<0.05)")
    fig.tight_layout(); fig.savefig(ATT / "HCV1_CCA_gini.png", dpi=150)
    plt.close(fig)


def fig_direction(df):
    learn = df[df["learner"] == 1]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    # left: per-pair mean IFI (per-animal mean then average), sign test vs 0
    ax = axes[0]
    meds, ps = [], []
    for p in PAIRS:
        per_an = [_num(sub["ifi"]).mean()
                  for _, sub in learn[learn["pair"] == p].groupby("animal")]
        per_an = [v for v in per_an if np.isfinite(v)]
        meds.append(np.median(per_an) if per_an else np.nan)
        ps.append(paired_stats.wilcoxon_signed(per_an)[3] if len(per_an) >= 3 else np.nan)
    colours = ["#2c3e50" if (np.isfinite(pp) and pp < 0.05) else "#95a5a6" for pp in ps]
    ax.bar(PAIRS, meds, color=colours)
    ax.axhline(0, color="k", lw=0.6)
    for i, pp in enumerate(ps):
        if np.isfinite(pp) and pp < 0.05:
            ax.text(i, meds[i], "*", ha="center", fontsize=14)
    ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("median IFI   (>0: first area leads)")
    ax.set_title("Directionality by pair (dark = p<0.05 vs 0)")
    # right: CA1-DG IFI trajectories per animal
    ax = axes[1]
    for an, g in learn[learn["pair"] == "CA1-DG"].groupby("animal"):
        ax.plot(_num(g["trial_frac"]), _num(g["ifi"]), "-o", ms=3, alpha=0.5,
                color="#2980b9")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("trial fraction"); ax.set_ylabel("IFI (CA1→DG)")
    ax.set_title("CA1→DG directionality rises with learning")
    fig.tight_layout(); fig.savefig(ATT / "HCV1_CCA_direction.png", dpi=150)
    plt.close(fig)


def fig_transition():
    path = RES / "transition_uncued_cued.csv"
    if not path.is_file():
        return
    t = pd.read_csv(path)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    data = [_num(t[t["pair"] == p]["angle_x"]).dropna().to_numpy() for p in PAIRS]
    data = [d for d in data if d.size]
    labs = [p for p, d in zip(PAIRS, [_num(t[t["pair"] == p]["angle_x"]).dropna().to_numpy()
            for p in PAIRS]) if d.size]
    ax.boxplot(data, labels=labs, vert=True)
    ax.axhline(90, color="grey", ls=":", label="orthogonal (90°)")
    ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("uncued→cued subspace angle (°)")
    ax.set_title("Subspace rotation at task onset (CA3–DG most stable)")
    ax.legend(fontsize=8)
    ax = axes[1]
    g = t[(t["pair"] == "CA1-DG") & (t["learner"] == 1)]
    d = _num(g["d_ifi"]).dropna().to_numpy()
    ax.plot(np.zeros_like(d), d, "o", color="#2980b9", alpha=0.7)
    ax.boxplot([d], positions=[0], widths=0.3)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks([]); ax.set_ylabel("Δ IFI  (cued − uncued), CA1→DG")
    ax.set_title("CA1→DG flow higher in cued task (p=0.016)")
    fig.tight_layout(); fig.savefig(ATT / "HCV1_CCA_transition.png", dpi=150)
    plt.close(fig)


def main():
    ATT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(RES / "trajectory_windows.csv")
    fig_levels(df); fig_gini(df); fig_direction(df); fig_transition()
    print("wrote figures to", ATT)
    for p in ["HCV1_CCA_levels.png", "HCV1_CCA_gini.png",
              "HCV1_CCA_direction.png", "HCV1_CCA_transition.png"]:
        print("  ", p, "ok" if (ATT / p).is_file() else "MISSING")


if __name__ == "__main__":
    main()
