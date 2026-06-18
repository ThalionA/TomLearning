"""B/D 5-trial-bin trajectories (Tom): each metric across the 6 ordinal trial bins
[1-5],[6-10],…,[26-30], with LEARNERS vs NON-LEARNERS as two mean ± shaded-SEM lines,
per area-pair (8-pair grid). The trial-NUMBER analogue of the sliding-window trajectory.

Reads results/trajectory_bins_bin10{,_fsincl}.csv (+ _dims for the sig-dim strength).
Metrics: cc1, n_sig, ifi, gini_x (per-pair), mincc = mean held-out CC over a bin's
significant canonical dims (per-bin strength-at-significant-dims proxy).

Usage: PYTHONPATH=src python scripts/figs_trajectory_bins.py [fsincl]
       (default FS-excl; pass 'fsincl' for the FS-included variant)
Note: Gini here is the per-PAIR weight Gini; a per-AREA-across-partners Gini over bins
(Tom's directive #1) would use trajectory_bins_weights — add as a variant if wanted.
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
METRICS = [("cc1", "CC₁"), ("n_sig", "# sig dims"), ("ifi", "IFI (>0: first leads)"),
           ("gini_x", "weight Gini (X)"), ("mincc", "mean CC over sig dims")]
COL = {1: "#c0392b", 0: "#2c6fbb"}
LAB = {1: "learners", 0: "non-learners"}


def _line(sub, metric):
    """mean ± SEM across animals per bin for one (pair, learner) slice."""
    xs, m, se = [], [], []
    for b in sorted(sub["bin"].unique()):
        v = pd.to_numeric(sub[sub["bin"] == b][metric], errors="coerce").to_numpy()
        v = v[np.isfinite(v)]
        if v.size:
            xs.append(int(b)); m.append(v.mean())
            se.append(v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else 0.0)
    return np.array(xs), np.array(m), np.array(se)


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    df = pd.read_csv(RES / f"trajectory_bins_bin10{suf}.csv")
    dd = pd.read_csv(RES / f"trajectory_bins_dims_bin10{suf}.csv")
    dd["cc"] = pd.to_numeric(dd["cc"], errors="coerce")
    mincc = (dd[dd["sig"] == 1].groupby(["animal", "learner", "pair", "bin"])["cc"]
             .mean().reset_index(name="mincc"))
    df = df.merge(mincc, on=["animal", "learner", "pair", "bin"], how="left")
    blab = df.groupby("bin").agg(lo=("trial_lo", "min"), hi=("trial_hi", "max"))
    xticklab = {int(b): f"{int(r.lo)}-{int(r.hi)}" for b, r in blab.iterrows()}

    for metric, ylab in METRICS:
        fig, axes = figstyle.grid(len(PAIRS), ncols=4)
        for ax, pair in zip(axes, PAIRS):
            sp = df[df["pair"] == pair]
            if sp.empty:
                ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
            for lrn in (1, 0):
                s = sp[sp["learner"] == lrn]
                if s.empty:
                    continue
                xs, m, se = _line(s, metric)
                if xs.size:
                    na = int(s.groupby("bin")["animal"].nunique().max())
                    ax.plot(xs, m, "-o", color=COL[lrn], ms=4, lw=1.6,
                            label=f"{LAB[lrn]} (n≤{na})")
                    ax.fill_between(xs, m - se, m + se, color=COL[lrn], alpha=0.15)
            if metric == "ifi":
                ax.axhline(0, color="k", lw=0.5)
            ax.set_title(pair, fontsize=10)
            ks = sorted(xticklab)
            ax.set_xticks(ks); ax.set_xticklabels([xticklab[b] for b in ks],
                                                  fontsize=6, rotation=45)
            ax.set_xlabel("trial bin", fontsize=8); ax.set_ylabel(ylab, fontsize=8)
            ax.legend(fontsize=6.5, loc="best", framealpha=0.6)
        fig.suptitle(f"{ylab} across 5-trial bins — learners vs non-learners — "
                     f"{fs}, 10 ms smoothed", fontsize=12)
        figstyle.save(fig, ATT / f"HCV1_trajbins_{metric}_{fs}_bin10.png")
        print(f"wrote HCV1_trajbins_{metric}_{fs}_bin10.png")


if __name__ == "__main__":
    main()
