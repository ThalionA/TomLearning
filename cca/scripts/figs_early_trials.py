"""Figures for the very-early-trial (pre-10) analysis.

From results/early_trials_projected{,_fsincl}.csv and early_trials_blocks{,_fsincl}.csv,
renders into the ResearchVault attachments, in the report house style (faint
per-animal lines + across-animal mean ± SEM, dots = animals, * = paired Wilcoxon):

  * HCV1_early_proj_<metric>   -- per-trial CC / IFI / participation-Gini at trial
                                  1/4/7/10, one panel per pair; * = trial-1-vs-k p<0.05
  * HCV1_early_block_<metric>  -- held-out CC / weight-Gini / angle-vs-late / KCCA at
                                  blocks t1-5/t1-7/t1-10/late; * = block-vs-late p<0.05
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


def _sem(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if v.size == 0:
        return np.nan, np.nan
    return float(v.mean()), (float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else 0.0)


def _series(df, pair, metric, level_col, level):
    g = df[(df["pair"] == pair) & (df[level_col] == level)]
    return {int(r): float(v) for r, v in zip(g["animal"], _num(g[metric])) if np.isfinite(v)}


def _star(base_map, other_map):
    common = sorted(set(base_map) & set(other_map))
    d = [other_map[a] - base_map[a] for a in common]
    return paired_stats.wilcoxon_signed(d)[3] if len(d) >= 3 else np.nan


def _pair_grid(df, metric, mlab, level_col, levels, base, xlabels, tag, fname,
               numeric_x=True):
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), sharex=True)
    xs = (np.array(levels, float) if numeric_x else np.arange(len(levels)))
    base_idx = levels.index(base)
    for ax, pair in zip(axes.ravel(), PAIRS):
        maps = {lv: _series(df, pair, metric, level_col, lv) for lv in levels}
        if all(len(m) == 0 for m in maps.values()):
            ax.set_visible(False); continue
        animals = sorted(set().union(*[set(m) for m in maps.values()]))
        for an in animals:                                   # faint per-animal lines
            ys = [maps[lv].get(an, np.nan) for lv in levels]
            ax.plot(xs, ys, color="#bdc3c7", lw=0.6, alpha=0.6, zorder=1)
        means = [_sem(list(maps[lv].values()))[0] for lv in levels]
        sems = [_sem(list(maps[lv].values()))[1] for lv in levels]
        ax.errorbar(xs, means, yerr=sems, color="#c0392b", lw=2, marker="o",
                    capsize=3, zorder=3, label="mean ± SEM")
        base_map = maps[base]
        for i, lv in enumerate(levels):                      # stars vs base
            if lv == base:
                continue
            p = _star(base_map, maps[lv])
            if np.isfinite(p) and p < 0.05:
                top = (means[i] or 0) + (sems[i] or 0)
                figstyle.star(ax, xs[i], top, always=True)
        ax.set_title(pair, fontsize=10)
        ax.set_xticks(xs); ax.set_xticklabels(xlabels, fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel(mlab)
    fig.suptitle(f"{mlab} over very early trials — {tag} "
                 f"(faint = animals, red = mean ± SEM; * = paired Wilcoxon vs {base})",
                 fontsize=12)
    figstyle.save(fig, ATT / fname)


def main():
    ATT.mkdir(parents=True, exist_ok=True)
    for tag in ("fsexcl", "fsincl"):
        suf = "" if tag == "fsexcl" else "_fsincl"
        pp = config.RESULTS_DIR / f"early_trials_projected{suf}.csv"
        bp = config.RESULTS_DIR / f"early_trials_blocks{suf}.csv"
        if not pp.is_file() or not bp.is_file():
            print(f"  (missing {pp.name}/{bp.name} — skip {tag})"); continue
        proj = pd.read_csv(pp); block = pd.read_csv(bp)

        for metric, mlab in [("cc1", "per-trial CC$_1$"), ("ifi", "per-trial IFI"),
                             ("gini_part_x", "participation-Gini (X)"),
                             ("gini_part_y", "participation-Gini (Y)")]:
            _pair_grid(proj, metric, mlab, "ordinal", [1, 4, 7, 10], 1,
                       ["1", "4", "7", "10"], tag,
                       f"HCV1_early_proj_{metric}_{tag}.png", numeric_x=True)

        blevels = ["t1-5", "t1-7", "t1-10", "late"]
        for metric, mlab in [("cc1", "held-out CC$_1$"), ("gini_x", "weight-Gini (X)"),
                             ("n_sig", "# sig dims"), ("ifi", "IFI"),
                             ("angle_x", "angle vs late (deg, X)"),
                             ("cc_kcca", "KCCA CC$_1$"),
                             ("gini_kcca_x", "KCCA struct-coeff Gini (X)")]:
            _pair_grid(block, metric, mlab, "block", blevels, "late",
                       blevels, tag, f"HCV1_early_block_{metric}_{tag}.png",
                       numeric_x=False)
        print(f"  {tag}: early-trial figures written")
    print("figures ->", ATT)


if __name__ == "__main__":
    main()
