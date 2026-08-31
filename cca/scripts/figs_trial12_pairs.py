"""Per-PAIR trial-1-vs-2 figures — the V1-RSC deep-dive layout for all 8 pairs.

One figure per (pair, FS): the first 10 running trials through the frozen
subspace (fit on ordinals 11+, leak-free; analyze_trial12's exact chain via
figs_trial12_v1rsc.load_tables — degenerate-drop, dim-matching, cc_peak
weighting).

Panels (2 × 3):
  A  frozen CC1 lag curve: ordinal 1 vs 2 vs the 3..10 band (common-bin arm)
  B  strength r0 across ordinals — all-sig CCs, r_frozen-weighted
  C  strength r0 across ordinals — CC1 only
  D  IFI ±50 ms across ordinals — all-sig CCs, cc_peak-weighted
  E  paired trial 1 vs 2, r0 (all-sig weighted, bin-matched arm) — paired t
  F  paired trial 1 vs 2, IFI ±50 ms (dim-matched, bin-matched arm) — paired t

The paired panels report the paired t-test (identical to the one-sample t on
within-animal differences) with Wilcoxon alongside; ⚠ per-pair family,
uncorrected — the 2026-08-20 verdict (0 BH survivors) stands.

Output: HCV1_trial12_pair_<PAIR>_<fs>_bin10  (16 figures)
Usage: PYTHONPATH=src python scripts/figs_trial12_pairs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402
from analyze_trial12 import HEADLINE_W  # noqa: E402
from figs_trial12_v1rsc import (load_tables, panel_lagcurve,  # noqa: E402
                                panel_ordinal, panel_paired)

figstyle.apply()
PAIRS = ["CA1-CA3", "CA1-DG", "CA3-DG", "CA1-SUB", "CA1-RSC", "CA1-V1",
         "RSC-SUB", "V1-RSC"]


def main():
    for fs, tag in [("fsexcl", ""), ("fsincl", "_fsincl")]:
        for pair in PAIRS:
            curves, allsig12, allsig_all, trials_allsig, trials_cc1 = \
                load_tables(tag, pair=pair)
            if curves.empty:
                print(f"  {pair} ({fs}): no data, skipped")
                continue
            lead = pair.split("-")[0]
            ifi12 = allsig12[allsig12["window_bins"] == HEADLINE_W]
            ifi_traj = allsig_all[allsig_all["window_bins"] == HEADLINE_W]
            fig, axes = plt.subplots(2, 3, figsize=(11.8, 7.0),
                                     constrained_layout=True)
            panel_lagcurve(axes[0, 0], curves, lead=lead)
            panel_ordinal(axes[0, 1], trials_allsig, "wmean",
                          "r0 (all-sig weighted)",
                          "B  strength across ordinals (all sig CCs)")
            panel_ordinal(axes[0, 2], trials_cc1, "r0", "r0 (CC1 only)",
                          "C  strength across ordinals (CC1 only)")
            panel_ordinal(axes[1, 0], ifi_traj, "wmean",
                          f"IFI ±50 ms (+ve = {lead} leads)",
                          "D  direction across ordinals (all sig CCs)")
            panel_paired(axes[1, 1], trials_allsig, "wmean",
                         "r0 (all-sig weighted)",
                         "E  paired 1 vs 2 — strength")
            panel_paired(axes[1, 2], ifi12, "wmean", "IFI ±50 ms",
                         "F  paired 1 vs 2 — direction")
            fig.suptitle(f"{pair} through the frozen subspace, first 10 running "
                         f"trials ({'FS excluded' if fs == 'fsexcl' else 'FS included'})",
                         fontsize=12)
            figstyle.save(fig, f"HCV1_trial12_pair_{pair}_{fs}_bin10")
        print(f"{fs}: done")


if __name__ == "__main__":
    main()
