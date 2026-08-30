"""Figures: do the contribution × spatial-property links change with experience?

Reads the fit-level CSVs of the contribution thread and renders the epoch
story (naive = first 10 running trials, intermediate = 10 ending at LP,
expert = 10 after LP). Delta forests show WITHIN-animal epoch contrasts in
Fisher-z units (jaccard_excess raw), Wilcoxon stars, animals-as-n — the
computation is imported from analyze_contrib_epochs (no duplication).

Outputs per FS condition (<fs> = fsexcl | fsincl):
  * HCV1_contribepochs_pooled_traj_<fs>_bin10  — pooled rel/tune (rate-part)
    trajectories per area
  * HCV1_contribepochs_delta_pooled_<fs>_bin10 — pooled delta forest
  * HCV1_contribepochs_delta_overlap_<fs>_bin10— overlap delta forest
  * HCV1_contribepochs_delta_pairs_<fs>_bin10  — per-(pair, area) delta forest,
    contrib_conn vs reliability rate-partialled

Usage: PYTHONPATH=src python scripts/figs_contrib_epochs.py
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
from analyze_contrib_epochs import CONTRASTS, EPOCHS, deltas, epoch_values  # noqa: E402
from figs_contrib_overlap import AREA_COL, AREAS  # noqa: E402
from figs_contrib_reliability import COL_X, COL_Y, forest_axis, pair_rows  # noqa: E402

figstyle.apply()
RES = Path(__file__).resolve().parents[1] / "results"
CTITLES = {("expert", "naive"): "expert − naive",
           ("intermediate", "naive"): "intermediate − naive",
           ("expert", "intermediate"): "expert − intermediate"}


def delta_forest(fig_axes, ev_by_metric, row_specs, stem, fs, suptitle,
                 xlab="Δ Fisher-z(ρ)", xlim=(-1.5, 1.5)):
    """row_specs: list of None | (label, metric, cell, colour)."""
    fig, axes = fig_axes
    for ax, contrast in zip(axes, CONTRASTS):
        rows = []
        for spec in row_specs:
            if spec is None:
                rows.append(None)
                continue
            label, metric, cell, col = spec
            rows.append((label, deltas(ev_by_metric[metric], cell, *contrast), col))
        forest_axis(ax, rows, CTITLES[contrast], xlim=xlim)
        ax.set_xlabel(xlab)
    fig.suptitle(f"{suptitle} ({'FS excluded' if fs == 'fsexcl' else 'FS included'})",
                 fontsize=11)
    figstyle.save(fig, f"{stem}_{fs}_bin10")


def make_pooled_traj(pooled, fs):
    fig, axes = figstyle.grid(len(AREAS), ncols=3, panel=(2.7, 2.3))
    evs = {m: epoch_values(pooled, m, ["area"])
           for m in ("rho_pooled_ratepart", "rho_pooled_tune_ratepart")}
    for ax, area in zip(axes, AREAS):
        for metric, col, label in [("rho_pooled_ratepart", COL_X, "reliability"),
                                   ("rho_pooled_tune_ratepart", COL_Y, "tuning")]:
            m, s = [], []
            for epoch in EPOCHS:
                try:
                    v = evs[metric].loc[(area, slice(None), epoch)].to_numpy(float)
                except KeyError:
                    v = np.array([])
                v = np.tanh(v[np.isfinite(v)])
                m.append(v.mean() if v.size else np.nan)
                s.append(v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else np.nan)
            ax.errorbar(range(3), m, yerr=s, marker="o", ms=4, capsize=3,
                        color=col, label=label)
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_xticks(range(3), ["naive\n(first 10)", "interm.\n(10 at LP)",
                                 "expert\n(10 after)"], fontsize=7)
        ax.set_ylim(-0.5, 0.75)
        ax.set_title(area, fontsize=9)
        ax.set_ylabel("rate-partialled ρ", fontsize=8)
        ax.legend(fontsize=6, loc="upper left")
    fig.suptitle(f"Pooled contribution links across learning epochs, animals "
                 f"mean±SEM ({'FS excluded' if fs == 'fsexcl' else 'FS included'})",
                 fontsize=11)
    figstyle.save(fig, f"HCV1_contribepochs_pooled_traj_{fs}_bin10")


def main():
    for fs, tag in [("fsexcl", ""), ("fsincl", "_fsincl")]:
        rel = pd.read_csv(RES / f"contrib_reliability_bin10{tag}.csv")
        pooled = pd.read_csv(RES / f"contrib_pooled_bin10{tag}.csv")
        overlap = pd.read_csv(RES / f"contrib_overlap_bin10{tag}.csv")

        make_pooled_traj(pooled, fs)

        ev = {m: epoch_values(pooled, m, ["area"])
              for m in ("rho_pooled_ratepart", "rho_pooled_tune_ratepart")}
        specs = ([(f"rel · {a}", "rho_pooled_ratepart", a, AREA_COL[a])
                  for a in AREAS] + [None]
                 + [(f"tune · {a}", "rho_pooled_tune_ratepart", a, AREA_COL[a])
                    for a in AREAS])
        delta_forest(plt.subplots(1, 3, figsize=(12.0, 4.6), sharey=True,
                                  constrained_layout=True),
                     ev, specs, "HCV1_contribepochs_delta_pooled", fs,
                     "Epoch change — pooled contribution vs reliability / tuning "
                     "(rate-partialled)")

        ev = {m: epoch_values(overlap, m, ["area"])
              for m in ("rho_conn", "rho_resid", "jaccard_excess")}
        specs = ([(f"ρ_conn · {a}", "rho_conn", a, AREA_COL[a]) for a in AREAS]
                 + [None]
                 + [(f"ρ_resid · {a}", "rho_resid", a, AREA_COL[a]) for a in AREAS]
                 + [None]
                 + [(f"Jacc · {a}", "jaccard_excess", a, AREA_COL[a])
                    for a in AREAS])
        delta_forest(plt.subplots(1, 3, figsize=(12.0, 6.2), sharey=True,
                                  constrained_layout=True),
                     ev, specs, "HCV1_contribepochs_delta_overlap", fs,
                     "Epoch change — cross-partner overlap (Δz; Jaccard raw Δ)")

        ev = {"rel": epoch_values(rel, "rho_contrib_conn_ratepart",
                                  ["pair", "area"])}
        vals = {}   # pair_rows wants {(pair, area): rows}; build spec-style via forest
        specs = []
        for row in pair_rows({}):
            if row is None:
                specs.append(None)
            else:
                label = row[0]
                pair, area = label.split()
                specs.append((label, "rel", (pair, area),
                              COL_X if area == pair.split("-")[0] else COL_Y))
        delta_forest(plt.subplots(1, 3, figsize=(12.0, 6.6), sharey=True,
                                  constrained_layout=True),
                     ev, specs, "HCV1_contribepochs_delta_pairs", fs,
                     "Epoch change — per-pair contribution vs reliability "
                     "(rate-partialled)")
        print(f"{fs}: 4 figures written")


if __name__ == "__main__":
    main()
