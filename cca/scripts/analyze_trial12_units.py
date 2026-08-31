"""Trial 1 vs 2, per-unit: does WHO carries the communication change?

Reads results/trial12_units_bin10*.csv (run_trial12_units.py: carry_r = each
unit's correlation with the partner area's frozen CC1, per early ordinal;
ordinal 0 = the reference profile from the fit trials). Three questions, all
animals-as-n (paired t + Wilcoxon; BH across (pair, side) cells per metric):

1. STRENGTH of carrying: per (animal, pair, side, ordinal) mean |carry_r| over
   units; Δ(trial 1 − trial 2), plus the 1→2 step z-scored against the
   adjacent steps (3→4 .. 9→10 — steps free of trial-1 behaviour).
2. PROFILE stability: Spearman between the unit-profiles of ordinals (1, 2),
   vs the adjacent-pair profile correlations (3,4)..(9,10). Excess =
   z(sim_1,2) − mean z(sim_adjacent): does the very first experience reorder
   who carries the channel more than any later step does?
3. CONVERGENCE to the trained membership: sim(ordinal, reference profile);
   Δ(sim_1 − sim_2) — is trial 1's carrying pattern further from the eventual
   one than trial 2's?

Primary arm: `matched` = 1 (all ordinals cut to a common bin count — the
correlation-vs-n confound is removed). Raw arm printed for sensitivity.

Writes results/trial12_units_tests_bin10{,_fsincl}.csv and prints the tables.
Usage: PYTHONPATH=src python scripts/analyze_trial12_units.py [--include-fs]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, paired_stats  # noqa: E402

PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
MIN_UNITS = 5
ADJ = [(o, o + 1) for o in range(3, 10)]      # control steps, trial 1 excluded


def _z(r):
    return np.arctanh(np.clip(r, -0.999, 0.999))


def profile_sim(g, o1, o2):
    """Spearman between the (unit-aligned) profiles of two ordinals."""
    p1 = g[g["ordinal"] == o1].set_index("unit")["carry_r"]
    p2 = g[g["ordinal"] == o2].set_index("unit")["carry_r"]
    j = p1.index.intersection(p2.index)
    if len(j) < MIN_UNITS:
        return np.nan
    return float(stats.spearmanr(p1.loc[j], p2.loc[j]).statistic)


def bh(pvals):
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    out = np.zeros(p.size, dtype=bool)
    if ok.sum() == 0:
        return out
    q = stats.false_discovery_control(p[ok])
    out[ok] = q < 0.05
    return out


def report(per_cell, label, rows_out, metric):
    """per_cell: {(pair, side): per-animal array}. Wilcoxon + t + BH row."""
    print(f"\n[{metric}] {label}")
    keys, pv = [], []
    for key in sorted(k for k in per_cell if k[0] != "_an"):
        v = np.asarray(per_cell[key], float)
        v = v[np.isfinite(v)]
        if v.size < 3:
            continue
        t_p = float(stats.ttest_1samp(v, 0.0).pvalue)
        _, med, _, w_p = paired_stats.wilcoxon_signed(v.tolist())
        keys.append((key, v, med, t_p, w_p)); pv.append(w_p)
    passed = bh(pv)
    for (key, v, med, t_p, w_p), b in zip(keys, passed):
        pair, side = key
        mark = "*" if (np.isfinite(w_p) and w_p < 0.05) else " "
        print(f"  {pair:9s} {side} n={v.size:<2d} med={med:>+7.3f} "
              f"mean={v.mean():>+7.3f} | t p={t_p:.3g} | W p={w_p:.3g}{mark}"
              f"{'  BH' if b else ''}")
        rows_out.append(dict(metric=metric, pair=pair, side=side, n=v.size,
                             med=round(float(med), 4), mean=round(float(v.mean()), 4),
                             t_p=round(t_p, 5), w_p=round(w_p, 5), bh_pass=bool(b)))


def collect_deltas(df, matched=1):
    """The four per-cell delta dicts, keyed (pair, side) with parallel
    ("_an", pair, side) animal-id lists (see _push). Import point for figures."""
    early = df[(df["matched"] == matched) & (df["ordinal"] > 0)]
    ref = df[df["ordinal"] == 0]
    d_strength, d_strength_z = {}, {}
    d_sim12, d_conv = {}, {}

    def _push(d, key, an, val):
        """Append val under key AND its animal id under a parallel ("_an",)+key
        list, so per-animal aggregation can recover who is who."""
        d.setdefault(key, []).append(val)
        d.setdefault(("_an",) + key, []).append(an)

    for (an, pair, side), g in early.groupby(["animal", "pair", "side"]):
        key = (pair, side)
        # 1. mean |carry| per ordinal
        s = g.groupby("ordinal")["carry_r"].apply(
            lambda c: np.nanmean(np.abs(c)))
        if 1 in s.index and 2 in s.index:
            _push(d_strength, key, an, _z(s[1]) - _z(s[2]))
            adj = [_z(s[a]) - _z(s[b]) for a, b in ADJ
                   if a in s.index and b in s.index]
            if len(adj) >= 4 and np.std(adj) > 0:
                _push(d_strength_z, key, an,
                      ((_z(s[1]) - _z(s[2])) - np.mean(adj)) / np.std(adj, ddof=1))
        # 2. profile similarity 1-2 vs adjacent
        s12 = profile_sim(g, 1, 2)
        adj_sims = [profile_sim(g, a, b) for a, b in ADJ]
        adj_sims = [x for x in adj_sims if np.isfinite(x)]
        if np.isfinite(s12) and len(adj_sims) >= 4:
            _push(d_sim12, key, an,
                  _z(s12) - float(np.mean(_z(np.array(adj_sims)))))
        # 3. convergence to the reference profile
        gr = ref[(ref["animal"] == an) & (ref["pair"] == pair)
                 & (ref["side"] == side)]
        if not gr.empty:
            gg = pd.concat([g, gr])
            c1 = profile_sim(gg, 1, 0)
            c2 = profile_sim(gg, 2, 0)
            if np.isfinite(c1) and np.isfinite(c2):
                _push(d_conv, key, an, _z(c1) - _z(c2))

    return d_strength, d_strength_z, d_sim12, d_conv


def per_animal_cells(d):
    out = {}
    for key, vals in d.items():
        if key[0] == "_an":
            continue
        for an, val in zip(d.get(("_an",) + key, []), vals):
            if np.isfinite(val):
                out.setdefault(an, []).append(val)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-fs", action="store_true")
    ap.add_argument("--matched", type=int, default=1)
    args = ap.parse_args()
    tag = "_fsincl" if args.include_fs else ""
    df = pd.read_csv(config.RESULTS_DIR / f"trial12_units_bin10{tag}.csv")
    print(f"trial12 units ({'FS incl' if args.include_fs else 'FS excl'}), "
          f"arm matched={args.matched}: {len(df)} rows, "
          f"{df['animal'].nunique()} animals")
    d_strength, d_strength_z, d_sim12, d_conv = collect_deltas(df, args.matched)

    rows_out = []
    report(d_strength, "Δ mean |carry_r| (trial 1 − 2), Fisher-z", rows_out,
           "carry_strength_d12")
    report(d_strength_z, "1→2 strength step z-scored vs adjacent steps",
           rows_out, "carry_strength_zadj")
    report(d_sim12, "profile sim(1,2) − mean adjacent sim, Fisher-z", rows_out,
           "profile_sim12_excess")
    report(d_conv, "Δ sim-to-reference (trial 1 − 2), Fisher-z", rows_out,
           "profile_convergence_d12")

    # GLOBAL per-animal tests: the per-cell medians lean one way (e.g. sim(1,2)
    # below the adjacent band in most cells) without any cell surviving BH —
    # the honest summary is ONE value per animal (mean over its cells; cells
    # share animals so they are not independent tests), then animals-as-n.
    for d, metric in ((d_sim12, "profile_sim12_excess"),
                      (d_conv, "profile_convergence_d12"),
                      (d_strength, "carry_strength_d12")):
        per_an = {}
        for (an_key, vals) in per_animal_cells(d).items():
            per_an[an_key] = float(np.mean(vals))
        v = np.asarray(list(per_an.values()), float)
        v = v[np.isfinite(v)]
        if v.size < 3:
            continue
        t_p = float(stats.ttest_1samp(v, 0.0).pvalue)
        _, med, _, w_p = paired_stats.wilcoxon_signed(v.tolist())
        mark = "*" if (np.isfinite(w_p) and w_p < 0.05) else " "
        print(f"\n[GLOBAL {metric}] one value per animal (mean over its cells): "
              f"n={v.size} med={med:+.4f} mean={v.mean():+.4f} "
              f"| t p={t_p:.3g} | W p={w_p:.3g}{mark}")
        rows_out.append(dict(metric=f"GLOBAL_{metric}", pair="ALL", side="-",
                             n=v.size, med=round(float(med), 4),
                             mean=round(float(v.mean()), 4),
                             t_p=round(t_p, 5), w_p=round(w_p, 5),
                             bh_pass=bool(np.isfinite(w_p) and w_p < 0.05)))

    out = pd.DataFrame(rows_out)
    path = config.RESULTS_DIR / f"trial12_units_tests_bin10{tag}.csv"
    out.to_csv(path, index=False)
    n_star = int((out["w_p"] < 0.05).sum())
    print(f"\nwrote {path.name} ({len(out)} cells; {n_star} W-starred, "
          f"{int(out['bh_pass'].sum())} BH survivors)")


if __name__ == "__main__":
    main()
