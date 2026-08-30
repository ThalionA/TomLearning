"""Figures for contribution × spatial reliability (Fig-5 leg: are the units that
carry a communication subspace spatially special?).

Reads results/contrib_reliability_bin10*.csv (per animal × pair × epoch × area
Spearman rhos, written by analyze_contrib_reliability.py) and the per-unit join
results/contrib_reliability_units_bin10*.csv. Per-animal value everywhere =
Fisher-z mean of the per-epoch rho, shown back in rho units; stars = Wilcoxon
signed-rank of the per-animal values vs 0 (animals-as-n, per (pair, area) cell).

Outputs per FS condition (<fs> = fsexcl | fsincl):
  * HCV1_contribrel_forest_<fs>_bin10   — 16-cell forest: raw vs rate-partialled rho
  * HCV1_contribrel_controls_<fs>_bin10 — intrinsic-contrib control vs rate–reliability confound
  * HCV1_contribrel_scatter_<fs>_bin10  — per-unit scatters, median animal, robust cells
  * HCV1_contribrel_epochs_<fs>_bin10   — rate-partialled rho across epochs, per pair

Usage: PYTHONPATH=src python scripts/figs_contrib_reliability.py
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
from tom_cca import paired_stats  # noqa: E402

figstyle.apply()
RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-CA3", "CA1-DG", "CA3-DG", "CA1-SUB", "CA1-RSC", "CA1-V1",
         "RSC-SUB", "V1-RSC"]          # intra-HC, HC->out, cortico-cortical
EPOCHS = ["naive", "intermediate", "expert"]
COL_X, COL_Y = "#3b6fb6", "#c8552c"    # first-named (X) / second-named (Y) area
ROBUST_CELLS = [("CA1-RSC", "RSC"), ("CA1-V1", "V1"), ("V1-RSC", "RSC")]


def per_animal(df, col):
    """Fisher-z mean of per-epoch rho within each (animal, pair, area) -> rho units."""
    z = np.arctanh(df[col].clip(-0.999, 0.999))
    m = z.groupby([df["animal"], df["pair"], df["area"]]).mean()
    return np.tanh(m).rename(col).reset_index()


def cells(pair):
    return [(pair, pair.split("-")[0], COL_X), (pair, pair.split("-")[1], COL_Y)]


def forest_axis(ax, rows, title, xlim=(-1, 1)):
    """Generic per-animal forest. ``rows``: list of None (gap) or
    (label, per-animal-values array, colour). Dots + mean±SEM + Wilcoxon star."""
    ypos, labels = [], []
    y = 0.0
    for row in rows:
        if row is None:
            y += 0.6                               # gap between groups
            continue
        label, v, col = row
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        ypos.append(y); labels.append(label)
        if v.size:
            ax.scatter(v, np.full(v.size, y), s=14, color=col, alpha=0.55, zorder=3)
            ax.errorbar(v.mean(), y, xerr=v.std(ddof=1) / np.sqrt(v.size), fmt="o",
                        ms=6, color=col, mec="black", mew=0.6, capsize=3, zorder=4)
            _, _, _, p = paired_stats.wilcoxon_signed(v.tolist())
            # figstyle.star() offsets upward for bar charts; on this inverted
            # categorical axis that lands between rows, so anchor va=center.
            # x in AXES fraction (get_yaxis_transform) so any xlim works.
            if np.isfinite(p) and p < 0.05:
                ax.text(0.93, y, "*", transform=ax.get_yaxis_transform(),
                        ha="center", va="center", fontsize=15, fontweight="bold",
                        color=figstyle.STAR_COLOR, zorder=20)
            ax.text(1.02, y, f"n={v.size}", transform=ax.get_yaxis_transform(),
                    va="center", fontsize=7, color="0.45")
        y += 1
    ax.axvline(0, color="0.6", lw=0.8, zorder=1)
    ax.set_yticks(ypos, labels, fontsize=8)
    ax.set_ylim(y - 0.5, -0.9)
    ax.set_xlim(*xlim)
    ax.set_xlabel("Spearman ρ (per-animal Fisher-z mean)")
    ax.set_title(title, fontsize=10)


def pair_rows(vals):
    """rows for the 16 (pair, area) cells with a gap between pairs."""
    rows = []
    for pair in PAIRS:
        for _, area, col in cells(pair):
            rows.append((f"{pair}  {area}", vals.get((pair, area), []), col))
        rows.append(None)
    return rows[:-1]


def make_forest(fit, fs, stem, cols_titles, suptitle="Unit contribution to the "
                "communication subspace vs spatial reliability"):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 6.4), sharey=True,
                             constrained_layout=True)
    for ax, (col, title) in zip(axes, cols_titles):
        pa = per_animal(fit, col)
        vals = {(p, a): g[col].to_numpy()
                for (p, a), g in pa.groupby(["pair", "area"])}
        forest_axis(ax, pair_rows(vals), title)
    h = [plt.Line2D([], [], marker="o", ls="", color=c, label=l)
         for c, l in [(COL_X, "first-named area (X)"), (COL_Y, "second-named area (Y)")]]
    axes[0].legend(handles=h, loc="lower left", fontsize=7, framealpha=0.9)
    fig.suptitle(f"{suptitle} ({'FS excluded' if fs == 'fsexcl' else 'FS included'})",
                 fontsize=11)
    figstyle.save(fig, f"{stem}_{fs}_bin10")


def make_scatter(units, fit, fs):
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9), constrained_layout=True)
    for ax, (pair, area) in zip(axes, ROBUST_CELLS):
        pa = per_animal(fit[(fit["pair"] == pair) & (fit["area"] == area)],
                        "rho_contrib_conn_ratepart")
        pa = pa[np.isfinite(pa["rho_contrib_conn_ratepart"])]
        an = pa.sort_values("rho_contrib_conn_ratepart")["animal"].iloc[len(pa) // 2]
        sub = units[(units["pair"] == pair) & (units["area"] == area)
                    & (units["animal"] == an)]
        for epoch, mk in zip(EPOCHS, ["o", "s", "^"]):
            e = sub[sub["epoch"] == epoch]
            ax.scatter(e["contrib_conn"], e["reliability"], s=16, marker=mk,
                       alpha=0.65, label=epoch)
        rho = float(pa.loc[pa["animal"] == an, "rho_contrib_conn_ratepart"].iloc[0])
        ax.set_title(f"{pair} · {area} units — animal {an} (median)\n"
                     f"rate-partialled ρ̄ = {rho:+.2f}", fontsize=9)
        ax.set_xlabel("contribution to subspace (a.u.)")
        ax.set_ylabel("spatial reliability (mean ±2-trial map r)")
        ax.legend(fontsize=7, title="epoch", title_fontsize=7)
    fig.suptitle(f"Per-unit contribution vs reliability — median animal per robust cell "
                 f"({'FS excluded' if fs == 'fsexcl' else 'FS included'})", fontsize=11)
    figstyle.save(fig, f"HCV1_contribrel_scatter_{fs}_bin10")


def make_epochs(fit, fs):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4, panel=(2.6, 2.3))
    for ax, pair in zip(axes, PAIRS):
        for _, area, col in cells(pair):
            g = fit[(fit["pair"] == pair) & (fit["area"] == area)]
            m, s = [], []
            for epoch in EPOCHS:
                v = g.loc[g["epoch"] == epoch, "rho_contrib_conn_ratepart"].to_numpy()
                v = v[np.isfinite(v)]
                m.append(v.mean() if v.size else np.nan)
                s.append(v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else np.nan)
            ax.errorbar(range(3), m, yerr=s, marker="o", ms=4, capsize=3,
                        color=col, label=area)
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_xticks(range(3), ["naive", "interm.", "expert"], fontsize=7)
        ax.set_ylim(-0.65, 0.85)
        ax.set_title(pair, fontsize=9)
        ax.set_ylabel("rate-partialled ρ", fontsize=8)
        ax.legend(fontsize=6, loc="upper left")
    fig.suptitle(f"Contribution–reliability link across learning epochs, animals mean±SEM "
                 f"({'FS excluded' if fs == 'fsexcl' else 'FS included'})", fontsize=11)
    figstyle.save(fig, f"HCV1_contribrel_epochs_{fs}_bin10")


def main():
    for fs, tag in [("fsexcl", ""), ("fsincl", "_fsincl")]:
        fit = pd.read_csv(RES / f"contrib_reliability_bin10{tag}.csv")
        units = pd.read_csv(RES / f"contrib_reliability_units_bin10{tag}.csv")
        make_forest(fit, fs, "HCV1_contribrel_forest",
                    [("rho_contrib_conn", "raw"),
                     ("rho_contrib_conn_ratepart", "rate-partialled")])
        make_forest(fit, fs, "HCV1_contribrel_controls",
                    [("rho_contrib", "area-intrinsic contribution (control)"),
                     ("rho_rate_rel", "mean rate vs reliability (confound)")])
        if "rho_contrib_conn_tune" in fit.columns:
            make_forest(fit, fs, "HCV1_contribtune_forest",
                        [("rho_contrib_conn_tune", "raw"),
                         ("rho_contrib_conn_tune_ratepart", "rate-partialled")],
                        suptitle="Unit contribution to the communication subspace "
                                 "vs spatial tuning (z vs shuffle null)")
        make_scatter(units, fit, fs)
        make_epochs(fit, fs)
        print(f"{fs}: 5 figures written")


if __name__ == "__main__":
    main()
