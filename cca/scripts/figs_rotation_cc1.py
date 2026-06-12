"""Figure: subspace rotation vs the split-half noise floor — top-3 vs CC1-only.

The §3.4 point made visual. Two panels (top-3 subspace | CC1 dominant dim). Per pair:
cross-window rotation and the within-window split-half floor, each as mean ± SEM over
per-animal medians with per-animal dots; a star marks pairs where rotation is
significantly BELOW the floor (signed-rank on per-animal rot−floor) — i.e. the
direction is *more* stable than estimation noise, so it does not reorient. The CC1
floor is far tighter than the top-3's near-orthogonal ~80°, making the test meaningful.

Needs a trajectory CSV with the CC1 columns (rot_x_cc1 / sh_x_cc1) — i.e. a post-#14
run (`trajectory_w15_bin50` currently). Usage: python figs_rotation_cc1.py <stem>.
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
from tom_cca import config, paired_stats  # noqa: E402
import figstyle  # noqa: E402

figstyle.apply()
ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _per_animal_median(g, col):
    return [float(_num(s[col]).median()) for _, s in g.groupby("animal")
            if np.isfinite(_num(s[col]).median())]


def _panel(ax, learn, rot_col, sh_col, title):
    xs = np.arange(len(PAIRS)); w = 0.38
    for off, col, lab, c in [(-w / 2, rot_col, "rotation", "#8e44ad"),
                             (+w / 2, sh_col, "split-half floor", "#95a5a6")]:
        means, sems = [], []
        for i, p in enumerate(PAIRS):
            v = _per_animal_median(learn[learn["pair"] == p], col)
            m = float(np.mean(v)) if v else np.nan
            s = float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0
            means.append(m); sems.append(s)
            jit = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.25
            ax.scatter(np.full(len(v), i) + off + jit, v, s=10, color="#2c3e50",
                       alpha=0.7, zorder=3)
        ax.bar(xs + off, means, w, yerr=sems, capsize=3, color=c, alpha=0.7,
               label=lab, zorder=1)
    # star where rotation is significantly BELOW the floor (stable direction)
    for i, p in enumerate(PAIRS):
        g = learn[learn["pair"] == p]
        diffs = [rm - fm for rm, fm in zip(
            (_num(s[rot_col]).median() for _, s in g.groupby("animal")),
            (_num(s[sh_col]).median() for _, s in g.groupby("animal")))
            if np.isfinite(rm) and np.isfinite(fm)]
        if len(diffs) >= 3:
            pv = paired_stats.wilcoxon_signed(diffs)[3]
            if np.isfinite(pv) and pv < 0.05 and np.median(diffs) < 0:
                top = max(np.nanmean(_per_animal_median(g, rot_col) or [0]),
                          np.nanmean(_per_animal_median(g, sh_col) or [0]))
                figstyle.star(ax, i, top, always=True)
    ax.set_xticks(xs); ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("max principal angle (deg)"); ax.set_ylim(0, 95)
    ax.axhline(90, color="k", lw=0.4, ls=":")
    ax.set_title(title, fontsize=10); ax.legend(fontsize=8, loc="lower right")


def main():
    args = " ".join(sys.argv[1:])
    stem = next((a for a in sys.argv[1:] if a.startswith("trajectory")), "trajectory_w15_bin50")
    pool = "pool" in args
    df = pd.read_csv(config.RESULTS_DIR / f"{stem}.csv")
    if "rot_x_cc1" not in df.columns:
        print(f"  {stem}.csv has no CC1-rotation columns (needs a post-#14 run)"); return
    data = df if pool else df[df["learner"] == 1]
    n = data["animal"].nunique()
    cohort = f"all {n} animals (learners + non)" if pool else f"{n} learners"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    _panel(a1, data, "rot_x", "sh_x", "Top-3 subspace (floor ≈ orthogonal — weak test)")
    _panel(a2, data, "rot_x_cc1", "sh_x_cc1", "CC1 dominant dim (tight floor — meaningful)")
    fig.suptitle(f"Subspace rotation stays at/below the split-half noise floor — no reorientation — "
                 f"{cohort}\n(bars = mean ± SEM over animals; dots = animals; "
                 "* = rotation sig. BELOW floor, signed-rank)", fontsize=11)
    out = f"HCV1_CCA_fsexcl{'_pooled' if pool else ''}_rotation_cc1"
    figstyle.save(fig, ATT / out)
    print(f"wrote {out}.{{png,svg}} ->", ATT)


if __name__ == "__main__":
    main()
