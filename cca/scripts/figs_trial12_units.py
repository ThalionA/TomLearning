"""Figure: per-unit carrying of the frozen subspace, trial 1 vs 2.

Renders analyze_trial12_units.py's three questions as one forest figure per FS
(computation imported from the analyzer — no duplication). Primary arm
matched=1 (common bin count). Rows = 16 (pair, side) cells + a GLOBAL row
(one value per animal, mean over its cells); dots = animals, mean ± SEM,
Wilcoxon stars (per-cell stars uncorrected; none survives BH — the global row
carries the inference).

Output: HCV1_trial12_units_deltas_<fs>_bin10

Usage: PYTHONPATH=src python scripts/figs_trial12_units.py
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
from analyze_trial12_units import collect_deltas, per_animal_cells  # noqa: E402
from figs_contrib_reliability import COL_X, COL_Y, forest_axis  # noqa: E402

figstyle.apply()
RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-CA3", "CA1-DG", "CA3-DG", "CA1-SUB", "CA1-RSC", "CA1-V1",
         "RSC-SUB", "V1-RSC"]
GLOBAL_COL = "#444444"


def rows_for(d):
    rows = []
    for pair in PAIRS:
        for side, col in (("x", COL_X), ("y", COL_Y)):
            area = pair.split("-")[0 if side == "x" else 1]
            rows.append((f"{pair}  {area}", d.get((pair, side), []), col))
        rows.append(None)
    glob = [float(np.mean(v)) for v in per_animal_cells(d).values()]
    rows.append(("GLOBAL (per animal)", glob, GLOBAL_COL))
    return rows


def main():
    for fs, tag in [("fsexcl", ""), ("fsincl", "_fsincl")]:
        df = pd.read_csv(RES / f"trial12_units_bin10{tag}.csv")
        d_strength, _, d_sim12, d_conv = collect_deltas(df, matched=1)
        fig, axes = plt.subplots(1, 3, figsize=(12.6, 6.6), sharey=True,
                                 constrained_layout=True)
        panels = [
            (d_strength, "Δ carrying strength\nmean |r| (trial 1 − 2)",
             "Δ Fisher-z", (-1.0, 1.0)),
            (d_sim12, "membership stability\nsim(1,2) − adjacent-step sim",
             "excess (Fisher-z)", (-1.0, 1.0)),
            (d_conv, "convergence to trained profile\nΔ sim-to-reference (1 − 2)",
             "Δ Fisher-z", (-1.0, 1.0)),
        ]
        for ax, (d, title, xlab, xlim) in zip(axes, panels):
            forest_axis(ax, rows_for(d), title, xlim=xlim)
            ax.set_xlabel(xlab)
        fig.suptitle(
            "Which units carry the channel on trials 1 and 2? Frozen subspace, "
            "common-bin arm — strength stable, membership reorders 1→2 "
            f"({'FS excluded' if fs == 'fsexcl' else 'FS included'})", fontsize=11)
        figstyle.save(fig, f"HCV1_trial12_units_deltas_{fs}_bin10")
        print(f"{fs}: 1 figure written")


if __name__ == "__main__":
    main()
