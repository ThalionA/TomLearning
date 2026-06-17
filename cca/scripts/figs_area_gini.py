"""Per-AREA Gini (across all partners) + weight cumulative distribution (Fig-4b style),
from the exported per-neuron contributions. Tom's directive: Gini per area, pooling every
pair the area participates in.

Aggregation: for area A in condition C, each partner-pair gives A's per-neuron contribution
vector; normalise each to unit sum (so every communication channel weighs equally), average
across partners → per-neuron mean participation barc; report Gini(barc) and the Lorenz curve
(cumulative share of participation vs neuron fraction — the curve whose Gini is barc's).

Outputs (smoothed 10 ms, FS-incl):
  * HCV1_areagini_<mode>_fsincl_bin10.png      — paired per-area Gini (cond A vs B)
  * HCV1_weightcdf_<mode>_fsincl_bin10.png     — per-area Lorenz curves, mean±SEM, A vs B

Usage: PYTHONPATH=src python scripts/figs_area_gini.py epochs|transition
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402
from tom_cca import membership, paired_stats  # noqa: E402

figstyle.apply()
ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"
RES = Path(__file__).resolve().parents[1] / "results"
AREAS = ["CA1", "CA3", "DG", "SUB", "RSC", "V1"]
GRID = np.linspace(0, 1, 101)          # common neuron-fraction grid for averaging Lorenz curves
SIG = "#1a7f4b"


def _agg_area(sub):
    """sub: rows (pair, unit, contrib) for one (animal, area, condition).
    Returns (gini, lorenz_on_GRID) of the partner-averaged per-neuron participation."""
    piv = sub.pivot_table(index="unit", columns="pair", values="contrib", aggfunc="mean")
    cols = []
    for p in piv.columns:
        c = piv[p].to_numpy(float)
        s = np.nansum(c)
        if s > 0:
            cols.append(np.nan_to_num(c) / s)        # normalise each partner to unit sum
    if not cols:
        return np.nan, None
    barc = np.mean(np.vstack(cols), axis=0)           # per-neuron mean participation
    barc = barc[np.isfinite(barc)]
    if barc.size < 2 or barc.sum() == 0:
        return np.nan, None
    s = np.sort(barc)
    lor = np.concatenate([[0], np.cumsum(s) / s.sum()])      # Lorenz, n+1 points
    x = np.linspace(0, 1, lor.size)
    return float(membership.gini(barc)), np.interp(GRID, x, lor)


def main():
    mode = "transition" if "transition" in sys.argv[1:] else "epochs"
    if mode == "epochs":
        wt = pd.read_csv(RES / "epoch_weights_bin10_fsincl.csv")
        condcol, A, B, labels = "epoch", "naive", "expert", ("naïve", "trained")
    else:
        wt = pd.read_csv(RES / "transition_weights_bin10_fsincl.csv")
        condcol, A, B, labels = "cond", "uncued", "cued", ("uncued", "cued")

    # per (animal, area, condition): Gini + Lorenz
    gini = {}                       # (area, cond) -> {animal: gini}
    lor = {}                        # (area, cond) -> list of Lorenz arrays
    for (an, area, cond), sub in wt.groupby(["animal", "area", condcol]):
        if cond not in (A, B):
            continue
        g, lc = _agg_area(sub)
        if np.isfinite(g):
            gini.setdefault((area, cond), {})[an] = g
            if lc is not None:
                lor.setdefault((area, cond), []).append(lc)

    # ---- paired per-area Gini ----
    fig, axes = figstyle.grid(len(AREAS), ncols=3)
    for ax, area in zip(axes, AREAS):
        ga, gb = gini.get((area, A), {}), gini.get((area, B), {})
        ans = sorted(set(ga) & set(gb))
        va = np.array([ga[a] for a in ans]); vb = np.array([gb[a] for a in ans])
        if va.size < 1:
            ax.set_title(f"{area}\n(no data)", fontsize=10); ax.axis("off"); continue
        for a, b in zip(va, vb):
            ax.plot([0, 1], [a, b], "-", color="#b4b2a9", lw=0.8, alpha=0.7)
        ax.scatter(np.zeros(va.size), va, s=22, color="#2c3e50", zorder=3)
        ax.scatter(np.ones(vb.size), vb, s=22, color="#2c3e50", zorder=3)
        for x, v in [(0, va), (1, vb)]:
            ax.plot([x - .18, x + .18], [v.mean(), v.mean()], color="#c0392b", lw=2.4)
        p_t = st.ttest_rel(va, vb)[1] if va.size >= 2 else np.nan
        _, _, _, p_w = paired_stats.wilcoxon_signed((vb - va).tolist())
        col = SIG if (np.isfinite(p_t) and p_t < .05) else "#3a3a36"
        ax.set_title(f"{area}\nt p={p_t:.3f} · W p={p_w:.2g} · n={va.size}", fontsize=10,
                     color=col, fontweight=("bold" if col == SIG else "normal"))
        ax.set_xlim(-.5, 1.5); ax.set_xticks([0, 1]); ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("per-area Gini")
    fig.suptitle(f"Per-area Gini (across all partners) — {labels[0]} vs {labels[1]} — "
                 "FS-incl, 10 ms smoothed", fontsize=12)
    figstyle.save(fig, ATT / f"HCV1_areagini_{mode}_fsincl_bin10.png")
    print(f"wrote HCV1_areagini_{mode}_fsincl_bin10.png")

    # ---- weight cumulative distribution (Lorenz), mean±SEM, A vs B ----
    fig, axes = figstyle.grid(len(AREAS), ncols=3)
    for ax, area in zip(axes, AREAS):
        ax.plot([0, 1], [0, 1], ":", color="#b4b2a9", lw=1)          # equality line
        for cond, c in ((A, "#5f5e5a"), (B, "#c0392b")):
            L = lor.get((area, cond))
            if not L:
                continue
            M = np.vstack(L); m = M.mean(0); se = M.std(0, ddof=1) / np.sqrt(len(L))
            ax.plot(GRID, m, color=c, lw=2, label=f"{labels[0] if cond==A else labels[1]} (n={len(L)})")
            ax.fill_between(GRID, m - se, m + se, color=c, alpha=0.2)
        ax.set_title(area, fontsize=11); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("fraction of neurons"); ax.set_ylabel("cum. participation")
        ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(f"Weight cumulative distribution (Lorenz; Gini = area between curve & diagonal) — "
                 f"{labels[0]} vs {labels[1]} — FS-incl, 10 ms smoothed", fontsize=12)
    figstyle.save(fig, ATT / f"HCV1_weightcdf_{mode}_fsincl_bin10.png")
    print(f"wrote HCV1_weightcdf_{mode}_fsincl_bin10.png")


if __name__ == "__main__":
    main()
