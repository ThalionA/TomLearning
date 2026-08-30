"""Figures for pooled contribution × reliability and cross-partner member overlap.

Reads results/contrib_pooled_bin10*.csv and contrib_overlap_bin10*.csv (written
by analyze_contrib_pooled_overlap.py). Per-animal value everywhere = Fisher-z
mean over epochs (and partner pairs, for the overlap); stars = Wilcoxon
signed-rank vs 0, animals-as-n per area.

Outputs per FS condition (<fs> = fsexcl | fsincl):
  * HCV1_contribpool_forest_<fs>_bin10  — pooled (all-partner) contribution vs
    reliability, raw vs rate-partialled, per area
  * HCV1_contriboverlap_<fs>_bin10      — cross-partner overlap: observed rho,
    intrinsic-geometry ceiling, partner-specific residual, top-quartile Jaccard
    excess over chance
  * HCV1_contriboverlap_ca1matrix_<fs>_bin10 — CA1 partner × partner matrices
    (observed | residual)

Usage: PYTHONPATH=src python scripts/figs_contrib_overlap.py
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
from figs_contrib_reliability import COL_X, COL_Y, forest_axis  # noqa: E402

figstyle.apply()
RES = Path(__file__).resolve().parents[1] / "results"
AREAS = ["CA1", "CA3", "DG", "SUB", "RSC", "V1"]      # HC formation, then cortex
AREA_COL = {a: (COL_X if a in ("CA1", "CA3", "DG", "SUB") else COL_Y) for a in AREAS}
CA1_PARTNERS = ["CA3", "DG", "SUB", "RSC", "V1"]


def per_animal_by_area(df, col, fisher=True):
    """{area: per-animal array}, Fisher-z mean over epochs (and partner pairs)."""
    v = df[col].clip(-0.999, 0.999)
    z = np.arctanh(v) if fisher else v
    m = z.groupby([df["area"], df["animal"]]).mean()
    if fisher:
        m = np.tanh(m)
    return {area: g.to_numpy() for area, g in m.groupby(level=0)}


def area_rows(vals):
    return [(a, vals.get(a, []), AREA_COL[a]) for a in AREAS]


def make_pooled(pooled, fs):
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), sharey=True,
                             constrained_layout=True)
    for ax, col, title in [(axes[0], "rho_pooled", "raw"),
                           (axes[1], "rho_pooled_ratepart", "rate-partialled")]:
        forest_axis(ax, area_rows(per_animal_by_area(pooled, col)), title)
    h = [plt.Line2D([], [], marker="o", ls="", color=c, label=l)
         for c, l in [(COL_X, "hippocampal formation"), (COL_Y, "cortex")]]
    axes[0].legend(handles=h, loc="lower left", fontsize=7, framealpha=0.9)
    fig.suptitle(f"Pooled contribution (mean over ALL partners) vs spatial reliability "
                 f"({'FS excluded' if fs == 'fsexcl' else 'FS included'})", fontsize=11)
    figstyle.save(fig, f"HCV1_contribpool_forest_{fs}_bin10")


def make_overlap(overlap, fs):
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.4), sharey=True,
                             constrained_layout=True)
    panels = [("rho_conn", "observed: contrib_conn\ncross-partner ρ", True),
              ("rho_intr", "ceiling: intrinsic contrib\n(shared geometry)", True),
              ("rho_resid", "partner-specific residual\n(geometry partialled out)", True),
              ("jaccard_excess", "top-quartile members:\nJaccard − chance", False)]
    for ax, (col, title, fisher) in zip(axes, panels):
        xlim = (-1, 1) if fisher else (-0.5, 0.75)
        forest_axis(ax, area_rows(per_animal_by_area(overlap, col, fisher)),
                    title, xlim=xlim)
        ax.set_xlabel("Spearman ρ" if fisher else "Jaccard − chance")
    fig.suptitle(f"Do the same units serve several communication subspaces? "
                 f"({'FS excluded' if fs == 'fsexcl' else 'FS included'})", fontsize=11)
    figstyle.save(fig, f"HCV1_contriboverlap_{fs}_bin10")


def make_ca1_matrix(overlap, fs):
    sub = overlap[overlap["area"] == "CA1"]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), constrained_layout=True)
    for ax, col, title in [(axes[0], "rho_conn", "observed contrib_conn ρ"),
                           (axes[1], "rho_resid", "partner-specific residual ρ")]:
        M = np.full((5, 5), np.nan)
        for i, p1 in enumerate(CA1_PARTNERS):
            for j, p2 in enumerate(CA1_PARTNERS):
                if i >= j:
                    continue
                pairset = {f"CA1-{p1}", f"CA1-{p2}"}
                g = sub[(sub["partner1"].isin(pairset)) & (sub["partner2"].isin(pairset))]
                if g.empty:
                    continue
                z = np.arctanh(g[col].clip(-0.999, 0.999))
                r = float(np.tanh(z.groupby(g["animal"]).mean().mean()))
                M[i, j] = M[j, i] = r
        im = ax.imshow(M, vmin=-1, vmax=1, cmap="RdBu_r")
        ax.set_xticks(range(5), CA1_PARTNERS, fontsize=8)
        ax.set_yticks(range(5), CA1_PARTNERS, fontsize=8)
        for i in range(5):
            for j in range(5):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center",
                            fontsize=8,
                            color="white" if abs(M[i, j]) > 0.6 else "black")
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.85, label="Spearman ρ (animal mean)")
    fig.suptitle(f"CA1 units: contribution similarity across partner subspaces "
                 f"({'FS excluded' if fs == 'fsexcl' else 'FS included'})", fontsize=11)
    figstyle.save(fig, f"HCV1_contriboverlap_ca1matrix_{fs}_bin10")


def main():
    for fs, tag in [("fsexcl", ""), ("fsincl", "_fsincl")]:
        pooled = pd.read_csv(RES / f"contrib_pooled_bin10{tag}.csv")
        overlap = pd.read_csv(RES / f"contrib_overlap_bin10{tag}.csv")
        make_pooled(pooled, fs)
        make_overlap(overlap, fs)
        make_ca1_matrix(overlap, fs)
        print(f"{fs}: 3 figures written")


if __name__ == "__main__":
    main()
