"""Pooled contribution vs reliability + cross-partner overlap of subspace members.

Two questions on top of analyze_contrib_reliability.py's per-unit join
(results/contrib_reliability_units_bin10*.csv — per animal × pair × epoch ×
area × unit rows: contrib, contrib_conn, reliability, mean_rate):

1. POOLED: average each unit's connection-specific contribution over ALL the
   area's partners (each partner's vector normalised to unit sum first, so every
   communication channel weighs equally — same aggregation as figs_area_gini),
   then Spearman vs reliability per (animal, area, epoch); Fisher-z over epochs;
   animals-as-n per AREA. Raw + log-rate rank-partialled.

2. OVERLAP: for each (animal, area, epoch) and each unordered partner pair
   (P1, P2): Spearman of the two contrib_conn vectors across the shared units.
   ⚠ The area-intrinsic contribution is ~partner-invariant (median r 0.981,
   PROJECT_LOG 2026-07-28), so a high cross-partner correlation is expected from
   shared loading geometry alone. Three numbers are therefore reported per cell:
     rho_conn      observed contrib_conn overlap,
     rho_intr      the intrinsic-contribution overlap (the geometry ceiling),
     rho_resid     Spearman of the two vectors after rank-partialling the mean
                   intrinsic contribution out of both (the partner-SPECIFIC part).
   Plus member-set overlap: Jaccard of the top-quartile sets
   (membership.member_mask) minus the value expected for independent same-size
   draws (E[|A∩B|] = k1 k2 / n).

Writes results/contrib_pooled_bin10{,_fsincl}.csv  (per animal × area × epoch)
       results/contrib_overlap_bin10{,_fsincl}.csv (per animal × area × epoch × partner-pair)
and prints animals-as-n tables for both.

Usage: PYTHONPATH=src python scripts/analyze_contrib_pooled_overlap.py [--include-fs]
"""
from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, membership, paired_stats  # noqa: E402

EPOCHS = ["naive", "intermediate", "expert"]
AREAS = ["CA1", "CA3", "DG", "SUB", "RSC", "V1"]
MIN_UNITS = 5


def rank_partial_rho(x, y, z):
    """Spearman of x vs y with z rank-partialled out of both (pairwise-complete)."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < MIN_UNITS:
        return np.nan
    rx, ry, rz = (stats.rankdata(v[ok]) for v in (x, y, z))
    res = []
    for r in (rx, ry):
        b = np.polyfit(rz, r, 1)
        res.append(r - np.polyval(b, rz))
    return float(stats.pearsonr(res[0], res[1])[0])


def spear(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < MIN_UNITS:
        return np.nan
    return float(stats.spearmanr(x[ok], y[ok])[0])


def pooled_contrib(sub):
    """sub: one (animal, area, epoch), rows per (pair, raw_unit). Returns a
    per-unit table with the partner-averaged normalised contrib_conn."""
    piv = sub.pivot_table(index="raw_unit", columns="pair", values="contrib_conn")
    cols = []
    for p in piv.columns:
        c = piv[p].to_numpy(float)
        s = np.nansum(c)
        if s > 0:
            cols.append(c / s)
    if not cols:
        return None
    pooled = pd.Series(np.nanmean(np.vstack(cols), axis=0), index=piv.index,
                       name="pooled_conn")
    per_unit = sub.groupby("raw_unit")[["reliability", "tuning_z",
                                        "reliability_tom_z", "mean_rate"]].mean()
    return per_unit.join(pooled)


def expected_jaccard(mask_a, mask_b):
    n = mask_a.size
    k1, k2 = int(mask_a.sum()), int(mask_b.sum())
    inter = k1 * k2 / n
    return inter / (k1 + k2 - inter) if (k1 + k2 - inter) > 0 else np.nan


def animals_as_n(df, value_col, group_cols, label):
    print("=" * 78)
    print(f"{label} — per-animal Fisher-z mean, one-sample t + Wilcoxon vs 0")
    print("=" * 78)
    z = np.arctanh(df[value_col].clip(-0.999, 0.999))
    per_an = np.tanh(z.groupby([df[c] for c in group_cols + ["animal"]]).mean())
    for key, g in per_an.groupby(level=list(range(len(group_cols)))):
        v = g.to_numpy()
        v = v[np.isfinite(v)]
        if v.size < 3:
            continue
        _, _, _, wp = paired_stats.wilcoxon_signed(v.tolist())
        tp = float(stats.ttest_1samp(v, 0.0).pvalue)
        mark = lambda p: "*" if (np.isfinite(p) and p < 0.05) else " "
        name = key if isinstance(key, str) else " ".join(map(str, key))
        print(f"  {name:16s} n={v.size:<2d} mean={np.mean(v):>+7.3f} "
              f"med={np.median(v):>+7.3f} | t p={tp:.3g}{mark(tp)} | W p={wp:.3g}{mark(wp)}")
    print()
    return per_an


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-fs", action="store_true")
    args = ap.parse_args()
    tag = "_fsincl" if args.include_fs else ""
    units = pd.read_csv(config.RESULTS_DIR / f"contrib_reliability_units_bin10{tag}.csv")

    pooled_rows, overlap_rows = [], []
    for (an, area, epoch), sub in units.groupby(["animal", "area", "epoch"]):
        piv = sub.pivot_table(index="raw_unit", columns="pair", values="contrib_conn")
        partners = list(piv.columns)

        # -- 1. pooled contribution vs reliability ---------------------------
        per_unit = pooled_contrib(sub)
        if per_unit is not None and len(per_unit) >= MIN_UNITS:
            logr = np.log(per_unit["mean_rate"].to_numpy() + 1e-6)
            pc = per_unit["pooled_conn"].to_numpy()
            pooled_rows.append(dict(
                animal=an, area=area, epoch=epoch, n_units=len(per_unit),
                n_partners=len(partners),
                rho_pooled=spear(pc, per_unit["reliability"].to_numpy()),
                rho_pooled_ratepart=rank_partial_rho(
                    pc, per_unit["reliability"].to_numpy(), logr),
                rho_pooled_tune=spear(pc, per_unit["tuning_z"].to_numpy()),
                rho_pooled_tune_ratepart=rank_partial_rho(
                    pc, per_unit["tuning_z"].to_numpy(), logr),
                rho_pooled_tomrel=spear(pc, per_unit["reliability_tom_z"].to_numpy()),
                rho_pooled_tomrel_ratepart=rank_partial_rho(
                    pc, per_unit["reliability_tom_z"].to_numpy(), logr)))

        # -- 2. cross-partner overlap ---------------------------------------
        if len(partners) < 2:
            continue
        intr = sub.pivot_table(index="raw_unit", columns="pair", values="contrib")
        for p1, p2 in combinations(partners, 2):
            both = piv[[p1, p2]].dropna()
            if len(both) < MIN_UNITS:
                continue
            a = both[p1].to_numpy()
            b = both[p2].to_numpy()
            mi = intr.loc[both.index, [p1, p2]].mean(axis=1).to_numpy()
            ma = membership.member_mask(a)
            mb = membership.member_mask(b)
            overlap_rows.append(dict(
                animal=an, area=area, epoch=epoch, partner1=p1, partner2=p2,
                n_units=len(both),
                rho_conn=spear(a, b),
                rho_intr=spear(intr.loc[both.index, p1].to_numpy(),
                               intr.loc[both.index, p2].to_numpy()),
                rho_resid=rank_partial_rho(a, b, mi),
                jaccard=membership.jaccard(ma, mb),
                jaccard_excess=membership.jaccard(ma, mb) - expected_jaccard(ma, mb)))

    pooled = pd.DataFrame(pooled_rows)
    overlap = pd.DataFrame(overlap_rows)
    pooled.to_csv(config.RESULTS_DIR / f"contrib_pooled_bin10{tag}.csv", index=False)
    overlap.to_csv(config.RESULTS_DIR / f"contrib_overlap_bin10{tag}.csv", index=False)
    print(f"pooled: {len(pooled)} rows | overlap: {len(overlap)} rows "
          f"({'FS incl' if args.include_fs else 'FS excl'})\n")

    for col, label in [("rho_pooled", "POOLED contrib (all partners) vs reliability, raw"),
                       ("rho_pooled_ratepart", "POOLED contrib vs reliability, rate-partialled"),
                       ("rho_pooled_tune", "POOLED contrib vs tuning z, raw"),
                       ("rho_pooled_tune_ratepart", "POOLED contrib vs tuning z, rate-partialled"),
                       ("rho_pooled_tomrel", "POOLED contrib vs Tom's reliability z, raw"),
                       ("rho_pooled_tomrel_ratepart", "POOLED contrib vs Tom's reliability z, rate-partialled")]:
        animals_as_n(pooled, col, ["area"], label)
    for col, label in [("rho_conn", "OVERLAP: contrib_conn cross-partner Spearman"),
                       ("rho_intr", "OVERLAP ceiling: intrinsic contrib cross-partner"),
                       ("rho_resid", "OVERLAP partner-specific: residual cross-partner")]:
        animals_as_n(overlap, col, ["area"], label)
    # Jaccard is not a correlation — plain per-animal mean, test excess vs 0.
    print("=" * 78)
    print("MEMBER-SET overlap: top-quartile Jaccard, excess over independent-draw chance")
    print("=" * 78)
    per_an = overlap.groupby(["area", "animal"])[["jaccard", "jaccard_excess"]].mean()
    for area, g in per_an.groupby(level=0):
        v = g["jaccard_excess"].to_numpy()
        v = v[np.isfinite(v)]
        if v.size < 3:
            continue
        _, _, _, wp = paired_stats.wilcoxon_signed(v.tolist())
        tp = float(stats.ttest_1samp(v, 0.0).pvalue)
        mark = lambda p: "*" if (np.isfinite(p) and p < 0.05) else " "
        print(f"  {area:5s} n={v.size:<2d} J={g['jaccard'].mean():.3f} "
              f"excess={np.mean(v):>+6.3f} med={np.median(v):>+6.3f} "
              f"| t p={tp:.3g}{mark(tp)} | W p={wp:.3g}{mark(wp)}")


if __name__ == "__main__":
    main()
