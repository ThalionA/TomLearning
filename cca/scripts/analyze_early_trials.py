"""Analyse the very-early-trial readouts — does any metric change in the first
~10 trials and then plateau?

Reads results/early_trials_projected{,_fsincl}.csv (per-trial, ordinals 1/4/7/10)
and results/early_trials_blocks{,_fsincl}.csv (cumulative blocks t1-5/t1-7/t1-10
vs a late reference). For every (pair, metric) it runs paired Wilcoxon across
animals:

  * PROJECTED: trial 1 vs 4, 1 vs 7, 1 vs 10 (the contrasts the user asked for).
  * BLOCK:     t1-5, t1-7, t1-10 each vs the late reference.

and flags an **early-effect-then-plateau**: an early contrast is significant
(p<0.05) AND the latest contrast (7-vs-10, or t1-10-vs-late) is not — i.e. the
metric moves early and has stopped moving by trial 10. Per-pair family, no
cross-pair correction (as elsewhere in the report); credibility rests on
consistency across pairs/metrics, not isolated p-values. Writes
results/early_trials_summary{,_fsincl}.csv.

Run with no args for FS-excluded; pass ``fsincl`` for the FS-included tables.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, paired_stats  # noqa: E402

PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
PROJ_METRICS = [("cc1", "CC"), ("ifi", "IFI"), ("optimal_lag", "lag"),
                ("gini_part_x", "partGini_x"), ("gini_part_y", "partGini_y")]
BLOCK_METRICS = [("cc1", "heldoutCC"), ("n_sig", "nsig"), ("mi_sig", "MIsig"),
                 ("ifi", "IFI"), ("gini_x", "weightGini_x"), ("gini_y", "weightGini_y"),
                 ("cc_kcca", "KCCA_cc"), ("gini_kcca_x", "KCCA_gini_x")]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _per_animal_at(df, pair, metric, level_col, level):
    """{animal -> value} for one pair at one trial-ordinal / block level."""
    g = df[(df["pair"] == pair) & (df[level_col] == level)]
    return {int(r["animal"]): v for r, v in
            zip(g.to_dict("records"), _num(g[metric])) if np.isfinite(v)}


def _paired(a_map, b_map):
    """Paired Wilcoxon on the per-animal (b-a) deltas over shared animals."""
    common = sorted(set(a_map) & set(b_map))
    d = [b_map[an] - a_map[an] for an in common]
    if len(d) < 3:
        return np.nan, np.nan, len(d)
    _, med, _, p = paired_stats.wilcoxon_signed(d)
    return med, p, len(d)


def _contrasts(df, level_col, base, others, metrics, label):
    print("=" * 92)
    print(f"{label} — paired Wilcoxon vs {level_col}={base} "
          "(Δ = later − base; * p<0.05; PLATEAU = early sig & last n.s.)")
    print("=" * 92)
    rows = []
    for pair in PAIRS:
        if df[df["pair"] == pair].empty:
            continue
        for metric, mlab in metrics:
            base_map = _per_animal_at(df, pair, metric, level_col, base)
            cells, res = [], {}
            for lv in others:
                m, p, n = _paired(base_map, _per_animal_at(df, pair, metric, level_col, lv))
                res[lv] = (m, p, n)
                star = "*" if (np.isfinite(p) and p < 0.05) else " "
                cells.append(f"{lv}:Δ{m:+.3f} p={p:.2g}{star}" if np.isfinite(m) else f"{lv}: -")
            early_sig = any(np.isfinite(res[lv][1]) and res[lv][1] < 0.05 for lv in others[:-1])
            last = res[others[-1]]
            last_ns = np.isfinite(last[1]) and last[1] >= 0.05
            plateau = early_sig and last_ns
            n = max((res[lv][2] for lv in others), default=0)
            tag = "  <-- EARLY→PLATEAU" if plateau else ""
            print(f"  {pair:8s} {mlab:13s} n={n:<2d} | " + " | ".join(cells) + tag)
            row = {"readout": label, "pair": pair, "metric": mlab, "n": n, "plateau": int(plateau)}
            for lv in others:
                row[f"d_{lv}"], row[f"p_{lv}"] = round(res[lv][0], 4) if np.isfinite(res[lv][0]) else "", \
                    round(res[lv][1], 4) if np.isfinite(res[lv][1]) else ""
            rows.append(row)
    return rows


def angle_vs_floor(df, blocks=("t1-5", "t1-7", "t1-10")):
    """Is the early block's subspace rotated from the late subspace BEYOND the
    within-block split-half noise floor? (Δ = angle − floor; signed-rank vs 0;
    >0 & significant = genuine early reorientation, not estimation noise — cf. §3.4.)"""
    print("=" * 92)
    print("BLOCK subspace ROTATION from late, above the split-half noise floor "
          "(Δ = angle − floor; signed-rank vs 0)")
    print("=" * 92)
    rows = []
    for pair in PAIRS:
        if df[df["pair"] == pair].empty:
            continue
        cells = []
        for blk in blocks:
            g = df[(df["pair"] == pair) & (df["block"] == blk)]
            d = [a - s for a, s in zip(_num(g["angle_x"]), _num(g["sh_x"]))
                 if np.isfinite(a) and np.isfinite(s)]
            if len(d) < 3:
                cells.append(f"{blk}: -"); continue
            _, med, _, p = paired_stats.wilcoxon_signed(d)
            star = "*" if (np.isfinite(p) and p < 0.05) else " "
            cells.append(f"{blk}:Δ{med:+.1f}° p={p:.2g}{star}")
            rows.append({"readout": "BLOCK angle", "pair": pair,
                         "metric": f"angle_minus_floor_{blk}", "n": len(d),
                         "plateau": 0, f"d_{blk}": round(med, 2),
                         f"p_{blk}": round(p, 4) if np.isfinite(p) else ""})
        print(f"  {pair:8s} angle−floor (X) | " + " | ".join(cells))
    return rows


def main():
    tag = "_fsincl" if (len(sys.argv) > 1 and "fsincl" in sys.argv[1]) else ""
    proj_p = config.RESULTS_DIR / f"early_trials_projected{tag}.csv"
    block_p = config.RESULTS_DIR / f"early_trials_blocks{tag}.csv"
    if not proj_p.is_file() or not block_p.is_file():
        print(f"missing {proj_p.name} / {block_p.name} — run run_early_trials.py first")
        return
    proj = pd.read_csv(proj_p)
    block = pd.read_csv(block_p)
    print(f"{proj_p.name}: {len(proj)} rows, {proj['animal'].nunique()} animals; "
          f"{block_p.name}: {len(block)} rows\n")

    rows = []
    rows += _contrasts(proj, "ordinal", 1, [4, 7, 10], PROJ_METRICS,
                       "PROJECTED per-trial (trial 1 vs 4/7/10)")
    print()
    rows += _contrasts(block, "block", "late", ["t1-5", "t1-7", "t1-10"], BLOCK_METRICS,
                       "BLOCK refit (early block vs late reference)")
    print()
    rows += angle_vs_floor(block)

    out = config.RESULTS_DIR / f"early_trials_summary{tag}.csv"
    pd.DataFrame(rows).to_csv(out, index=False, lineterminator="\n")
    plateaus = [r for r in rows if r["plateau"]]
    print(f"\n{len(plateaus)} (pair,metric) cells flagged EARLY→PLATEAU. summary -> {out}")
    for r in plateaus:
        print(f"   {r['readout'][:9]:9s} {r['pair']:8s} {r['metric']}")


if __name__ == "__main__":
    main()
