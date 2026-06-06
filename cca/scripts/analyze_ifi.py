"""Directionality (IFI) battery across learning epochs — animals-as-n AND dims-as-n.

For each area pair, the per-dimension information-flow index (IFI; >0 = first area
leads) is analysed across the naive / intermediate / expert epochs with the full
test set the user requested:

  ANIMALS-as-n (unit = animal; per-animal dominant-dim IFI from epoch_metrics):
    * IFI vs 0 at each epoch         -- one-sample t-test + Wilcoxon
    * naive vs expert                -- paired t-test (+ paired Wilcoxon)
    * across 3 epochs                -- RM-ANOVA (+ Friedman) with Holm-corrected
                                        pairwise post-hocs

  DIMS-as-n (unit = significant canonical dimension, pooled across animals;
  the Gonzalez/Buzsáki unit — pseudoreplicated, shown for comparison):
    * IFI vs 0 at each epoch         -- one-sample t-test + Wilcoxon
    * naive vs expert                -- rank-sum (Mann–Whitney U)
    * across 3 epochs                -- Kruskal–Wallis

Reads results/epoch_metrics.csv (per-animal) and results/epoch_dims.csv (per-dim).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config  # noqa: E402

EPOCHS = ["naive", "intermediate", "expert"]
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
# Set in main() from --metric: directionality readout to test.
#   ifi -> animal col 'ifi' (dom-dim IFI), dim col 'ifi'
#   lag -> animal col 'optimal_lag' (dom-dim lag, bins), dim col 'lag'
ANIMAL_COL = "ifi"
DIM_COL = "ifi"


def _f(x):
    x = np.asarray(x, float)
    return x[np.isfinite(x)]


def _vs0(vals):
    v = _f(vals)
    if v.size < 3:
        return v.size, np.nan, np.nan, np.nan
    t_p = float(stats.ttest_1samp(v, 0.0).pvalue)
    try:
        w_p = float(stats.wilcoxon(v).pvalue)
    except ValueError:
        w_p = np.nan
    return v.size, float(np.mean(v)), t_p, w_p


def _star(p):
    return "*" if (np.isfinite(p) and p < 0.05) else " "


def animals_block(em, pair):
    g = em[em["pair"] == pair]
    by_ep = {e: _f(pd.to_numeric(g[g["epoch"] == e][ANIMAL_COL], errors="coerce"))
             for e in EPOCHS}
    print(f"  ANIMALS-as-n:")
    for e in EPOCHS:
        n, m, tp, wp = _vs0(by_ep[e])
        if n >= 3:
            print(f"    dir vs0 {e:12s} n={n} mean={m:+.4f} "
                  f"t={tp:.3g}{_star(tp)} W={wp:.3g}{_star(wp)}")
    # naive vs expert (paired, by animal)
    gp = g.pivot_table(index="animal", columns="epoch", values=ANIMAL_COL, aggfunc="mean")
    if {"naive", "expert"} <= set(gp.columns):
        pair_df = gp[["naive", "expert"]].dropna()
        if len(pair_df) >= 3:
            tp = float(stats.ttest_rel(pair_df["expert"], pair_df["naive"]).pvalue)
            try:
                wp = float(stats.wilcoxon(pair_df["expert"], pair_df["naive"]).pvalue)
            except ValueError:
                wp = np.nan
            print(f"    nai→exp paired   n={len(pair_df)} "
                  f"Δ={float((pair_df['expert']-pair_df['naive']).mean()):+.4f} "
                  f"t={tp:.3g}{_star(tp)} W={wp:.3g}{_star(wp)}")
    # RM-ANOVA + Friedman + Holm post-hoc on complete subjects
    comp = gp.dropna(subset=EPOCHS) if set(EPOCHS) <= set(gp.columns) else pd.DataFrame()
    if len(comp) >= 3:
        long = comp[EPOCHS].reset_index().melt(id_vars="animal", var_name="epoch",
                                               value_name="ifi")
        try:
            aov = AnovaRM(long, "ifi", "animal", within=["epoch"]).fit()
            fp = float(aov.anova_table["Pr > F"].iloc[0])
        except Exception:
            fp = np.nan
        try:
            frp = float(stats.friedmanchisquare(*[comp[e] for e in EPOCHS]).pvalue)
        except ValueError:
            frp = np.nan
        # Holm-corrected pairwise paired t
        raw, labs = [], []
        for a, b in [("expert", "naive"), ("expert", "intermediate"),
                     ("intermediate", "naive")]:
            raw.append(float(stats.ttest_rel(comp[a], comp[b]).pvalue)); labs.append(f"{a[:3]}-{b[:3]}")
        holm = multipletests(raw, method="holm")[1]
        ph = "  ".join(f"{l}={p:.3g}{_star(p)}" for l, p in zip(labs, holm))
        print(f"    RM-ANOVA n={len(comp)} F-p={fp:.3g}{_star(fp)} "
              f"Friedman={frp:.3g}{_star(frp)} | Holm post-hoc: {ph}")


def dims_block(ed, pair):
    g = ed[(ed["pair"] == pair) & (ed["sig"] == 1)].copy()
    g["peak_cc"] = pd.to_numeric(g["peak_cc"], errors="coerce")
    g = g.sort_values("peak_cc", ascending=False).groupby(
        ["animal", "epoch"]).head(3)                     # cap top-3 dims/animal/epoch
    g[DIM_COL] = pd.to_numeric(g[DIM_COL], errors="coerce")
    by_ep = {e: _f(g[g["epoch"] == e][DIM_COL]) for e in EPOCHS}
    if all(v.size < 3 for v in by_ep.values()):
        return
    print(f"  DIMS-as-n (sig dims pooled):")
    for e in EPOCHS:
        n, m, tp, wp = _vs0(by_ep[e])
        if n >= 3:
            print(f"    dir vs0 {e:12s} nd={n} mean={m:+.4f} "
                  f"t={tp:.3g}{_star(tp)} W={wp:.3g}{_star(wp)}")
    if by_ep["naive"].size >= 3 and by_ep["expert"].size >= 3:
        up = float(stats.mannwhitneyu(by_ep["expert"], by_ep["naive"],
                                      alternative="two-sided").pvalue)
        print(f"    nai vs exp rank-sum  nd={by_ep['naive'].size}/{by_ep['expert'].size} "
              f"U-p={up:.3g}{_star(up)}")
    present = [by_ep[e] for e in EPOCHS if by_ep[e].size >= 3]
    if len(present) == 3:
        try:
            kp = float(stats.kruskal(*present).pvalue)
            print(f"    3-epoch Kruskal–Wallis  p={kp:.3g}{_star(kp)}")
        except ValueError:
            print("    3-epoch Kruskal–Wallis  n/a (ties)")


def main():
    global ANIMAL_COL, DIM_COL
    metric = sys.argv[1] if len(sys.argv) > 1 else "ifi"
    if metric not in ("ifi", "lag"):
        raise SystemExit("usage: analyze_ifi.py [ifi|lag]")
    ANIMAL_COL = "ifi" if metric == "ifi" else "optimal_lag"
    DIM_COL = "ifi" if metric == "ifi" else "lag"
    unit = ">0 = first area leads" if metric == "ifi" else "bins; >0 = first area leads (×25 ms)"
    em = pd.read_csv(config.RESULTS_DIR / "epoch_metrics.csv")
    ed = pd.read_csv(config.RESULTS_DIR / "epoch_dims.csv")
    print(f"DIRECTIONALITY battery — metric = {metric.upper()} ({unit}) | "
          f"epoch_metrics {len(em)} rows, epoch_dims {len(ed)} dim-rows\n")
    for pair in PAIRS:
        if em[em["pair"] == pair].empty:
            continue
        print(f"\n### {pair}")
        animals_block(em, pair)
        dims_block(ed, pair)


if __name__ == "__main__":
    main()
