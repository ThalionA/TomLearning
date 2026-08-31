"""Deep-dive: the V1-RSC trial-1-vs-2 candidate, in full.

The 2026-08-20 arm left ONE candidate: V1-RSC's 1→2 strength step is the most
negative of the adjacent-step band in both FS (uncorrected), and its all-CC
IFI delta is the only direction delta with p < 0.05 (uncorrected, no BH).
This script renders what that actually looks like and stress-tests it —
reusing analyze_trial12's exact table chain (degenerate-drop, dim-matching,
cc_peak weighting), no re-derivation.

Per FS condition, one figure HCV1_trial12_v1rsc_<fs>_bin10:
  A  CC1 frozen lag curve r(lag): ordinal 1 vs 2 vs the 3..10 band (matched=2)
  B  r0 (all-sig weighted) across ordinals 1..10, per-animal traces + mean
  C  IFI ±50 ms (all-sig weighted) across ordinals, same layout
  D  per-animal paired 1→2 lines for r0 and IFI (the n=9 the tests rest on)

Printed + written (results/trial12_v1rsc_sides_bin10*.csv): the IFI delta
decomposed into its sides — does V1→RSC (positive-lag mass) fall on trial 1,
or does RSC→V1 (negative-lag mass) rise? — plus arm (0/1/2) and window
sensitivity of the V1-RSC deltas.

Usage: PYTHONPATH=src python scripts/figs_trial12_v1rsc.py
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
from analyze_trial12 import (BIN_MS, HEADLINE_W, degenerate_trials,  # noqa: E402
                             dim_matched, drop_degenerate_trials,
                             per_trial_windows)
from tom_cca import cc_aggregate, lagged, paired_stats  # noqa: E402
from scipy import stats  # noqa: E402

figstyle.apply()
RES = Path(__file__).resolve().parents[1] / "results"
PAIR = "V1-RSC"
C1, C2, CB = "#c8552c", "#3b6fb6", "0.55"      # trial 1 / trial 2 / band 3..10


def load_tables(tag, pair=PAIR):
    curves = pd.read_csv(RES / f"trial12_curves_bin{BIN_MS}{tag}.csv")
    trials = pd.read_csv(RES / f"trial12_trials_bin{BIN_MS}{tag}.csv")
    flags = degenerate_trials(trials)
    curves = drop_degenerate_trials(curves, flags)
    trials = drop_degenerate_trials(trials, flags)
    curves = curves[curves["pair"] == pair]
    trials = trials[trials["pair"] == pair]
    win = per_trial_windows(curves)
    # contrast frame: dims matched between ordinals 1 and 2 (the tested delta)
    allsig12 = cc_aggregate.per_animal_mean(
        dim_matched(win, "ifi"), value="ifi",
        by=["ordinal", "matched", "window_bins"], weight="cc_peak",
        sig_only=False, drop_degenerate=False)
    # descriptive frame: ALL ordinals (per-ordinal sig/degenerate gates) — for
    # the trajectory panels only; dim_matched keeps just ordinals 1-2 by design
    allsig_all = cc_aggregate.per_animal_mean(
        win, value="ifi", by=["ordinal", "matched", "window_bins"],
        weight="cc_peak")
    trials_allsig = cc_aggregate.per_animal_mean(
        trials, value="r0", by=["ordinal", "matched"], weight="r_frozen")
    trials_cc1 = trials[trials["dim"] == 1][
        ["animal", "ordinal", "matched", "r0"]].copy()
    trials_cc1["r0"] = pd.to_numeric(trials_cc1["r0"], errors="coerce")
    return curves, allsig12, allsig_all, trials_allsig, trials_cc1


def panel_lagcurve(ax, curves, lead="V1"):
    cc1 = curves[(curves["dim"] == 1) & (curves["matched"] == 2)]
    for sel, col, label in [(cc1["ordinal"] == 1, C1, "trial 1"),
                            (cc1["ordinal"] == 2, C2, "trial 2"),
                            (cc1["ordinal"] >= 3, CB, "trials 3–10")]:
        g = cc1[sel]
        # mean per (animal, lag) first so the band is not dominated by 8 ordinals
        pa = g.groupby(["animal", "lag_ms"])["r"].mean().reset_index()
        m = pa.groupby("lag_ms")["r"].agg(["mean", "sem"])
        ax.plot(m.index, m["mean"], color=col, lw=1.6, label=label)
        ax.fill_between(m.index, m["mean"] - m["sem"], m["mean"] + m["sem"],
                        color=col, alpha=0.25, lw=0)
    ax.axvline(0, color="0.7", lw=0.7)
    ax.axhline(0, color="0.7", lw=0.7)
    ax.set_xlabel(f"lag (ms; +ve = {lead} leads)")
    ax.set_ylabel("frozen CC1 r (held out)")
    ax.set_title("A  CC1 lag curve (common-bin arm)", fontsize=9)
    ax.legend(fontsize=7)


def panel_ordinal(ax, pa, value, ylab, title):
    sub = pa[pa["matched"] == 2]
    for an, g in sub.groupby("animal"):
        g = g.sort_values("ordinal")
        ax.plot(g["ordinal"], g[value], color="0.75", lw=0.7, zorder=1)
    m = sub.groupby("ordinal")[value].agg(["mean", "sem"])
    ax.errorbar(m.index, m["mean"], yerr=m["sem"], color="black", marker="o",
                ms=4, capsize=3, zorder=3)
    for o, col in ((1, C1), (2, C2)):
        if o in m.index:
            ax.errorbar([o], [m.loc[o, "mean"]], yerr=[m.loc[o, "sem"]],
                        color=col, marker="o", ms=6, capsize=3, zorder=4)
    ax.axhline(0, color="0.7", lw=0.7)
    ax.set_xticks(range(1, 11))
    ax.set_xlabel("running-trial ordinal")
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=9)


def panel_paired(ax, pa, value, ylab, title):
    sub = pa[(pa["matched"] == 1) & (pa["ordinal"].isin([1, 2]))]
    wide = sub.pivot_table(index="animal", columns="ordinal", values=value)
    wide = wide.dropna()
    for an, row in wide.iterrows():
        ax.plot([0, 1], [row[1], row[2]], color="0.6", lw=0.9, marker="o", ms=3)
    ax.plot([0, 1], [wide[1].mean(), wide[2].mean()], color="black", lw=2.2,
            marker="o", ms=6)
    _, med, _, w_p = paired_stats.wilcoxon_signed((wide[1] - wide[2]).tolist())
    t_p = float(stats.ttest_1samp((wide[1] - wide[2]).to_numpy(), 0.0).pvalue)
    ax.set_xticks([0, 1], ["trial 1", "trial 2"])
    ax.set_xlim(-0.35, 1.35)
    ax.axhline(0, color="0.7", lw=0.7)
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nn={len(wide)}, Δmed={med:+.3f}, "
                 f"paired t p={t_p:.3g}, W p={w_p:.3g}", fontsize=8)


def sides_and_sensitivity(curves, allsig12, trials_allsig, tag, fs):
    """IFI side decomposition + arm/window sensitivity for V1-RSC. Returns rows."""
    print(f"\n=== {fs}: V1-RSC sensitivity ===")
    rows = []
    # sides of the ±HEADLINE_W window, CC1, per (animal, ordinal), matched=2
    cc1 = curves[(curves["dim"] == 1) & (curves["matched"] == 2)
                 & (np.abs(curves["lag_bins"]) <= HEADLINE_W)]
    for an, g in cc1.groupby("animal"):
        for o, go in g.groupby("ordinal"):
            lg, c = go["lag_bins"].to_numpy(float), go["r"].to_numpy(float)
            m = np.isfinite(c)
            if not (np.any(lg[m] > 0) and np.any(lg[m] < 0)):
                continue
            pm, nm = lagged.ifi_sides(lg[m], c[m])
            rows.append(dict(animal=an, ordinal=int(o), pos_mass=pm, neg_mass=nm))
    sides = pd.DataFrame(rows)
    for col, lab in (("pos_mass", "V1→RSC (positive-lag) mass"),
                     ("neg_mass", "RSC→V1 (negative-lag) mass")):
        w = sides[sides["ordinal"].isin([1, 2])].pivot_table(
            index="animal", columns="ordinal", values=col).dropna()
        d = (w[1] - w[2]).to_numpy()
        if d.size >= 4:
            _, med, _, w_p = paired_stats.wilcoxon_signed(d.tolist())
            t_p = float(stats.ttest_1samp(d, 0.0).pvalue)
            print(f"  Δ(1−2) {lab}: n={d.size} med={med:+.4f} "
                  f"t p={t_p:.3g} W p={w_p:.3g}")
    # arm sensitivity of the headline deltas
    for label, pa, val in (("IFI±50", allsig12[allsig12["window_bins"] == HEADLINE_W],
                            "wmean"), ("r0", trials_allsig, "wmean")):
        for arm in (0, 1, 2):
            w = pa[(pa["matched"] == arm) & (pa["ordinal"].isin([1, 2]))]
            w = w.pivot_table(index="animal", columns="ordinal", values=val).dropna()
            if len(w) < 4:
                continue
            d = (w[1] - w[2]).to_numpy()
            _, med, _, w_p = paired_stats.wilcoxon_signed(d.tolist())
            t_p = float(stats.ttest_1samp(d, 0.0).pvalue)
            print(f"  {label} arm={arm}: n={len(w)} Δmed={med:+.4f} "
                  f"t p={t_p:.3g} W p={w_p:.3g}")
    sides.to_csv(RES / f"trial12_v1rsc_sides_bin{BIN_MS}{tag}.csv", index=False)
    return sides


def main():
    for fs, tag in [("fsexcl", ""), ("fsincl", "_fsincl")]:
        curves, allsig12, allsig_all, trials_allsig, _ = load_tables(tag)
        ifi12 = allsig12[allsig12["window_bins"] == HEADLINE_W]
        ifi_traj = allsig_all[allsig_all["window_bins"] == HEADLINE_W]
        fig, axes = plt.subplots(1, 4, figsize=(13.6, 3.6),
                                 constrained_layout=True)
        panel_lagcurve(axes[0], curves)
        panel_ordinal(axes[1], trials_allsig, "wmean", "r0 (all-sig weighted)",
                      "B  strength across ordinals")
        panel_ordinal(axes[2], ifi_traj, "wmean", "IFI ±50 ms (all-sig weighted)",
                      "C  direction across ordinals")
        panel_paired(axes[3], ifi12, "wmean", "IFI ±50 ms",
                     "D  paired 1 vs 2 (bin-matched)")
        fig.suptitle(f"V1-RSC through the frozen subspace, first 10 running "
                     f"trials ({'FS excluded' if fs == 'fsexcl' else 'FS included'})",
                     fontsize=11)
        figstyle.save(fig, f"HCV1_trial12_v1rsc_{fs}_bin10")
        sides_and_sensitivity(curves, allsig12, trials_allsig, tag, fs)
        print(f"{fs}: figure written")


if __name__ == "__main__":
    main()
