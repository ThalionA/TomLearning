"""Mixed-effects naive->expert test on the landmark-arm cohort result.

The statistically correct counterpart to ``learning_changes.py``: instead of
pooling the 6 landmarks per animal into an independent-sample Wilcoxon (which
pseudoreplicates -- see ``STATE.md`` §3), this fits a per-pair linear mixed
model with a per-animal random slope on epoch, so the learning effect is tested
against BETWEEN-animal variation (n = animals, the honest unit).

For each pair and each stat (``mi_sig`` strength, ``ifi_weighted`` direction):
    value ~ C(epoch) + C(landmark),  random slope ~C(epoch) per animal.
P-values are Wald tests on the epoch contrasts, then Benjamini-Hochberg
corrected WITHIN each pair across its 3 contrasts (no cross-pair correction --
each area pair is its own pre-specified family).

Writes ``results/learning_changes_mixed_<tag>.csv``. Read-only over an existing
landmark pkl -- no re-run needed.
"""

from __future__ import annotations

import argparse
import csv
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import (config, mixed_effects, paired_stats,  # noqa: E402
                     subspace_stats)

STATS = ("mi_sig", "ifi_weighted")
CONTRASTS = mixed_effects.DEFAULT_CONTRASTS


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tag", default="landmark50_res_samp15")
    p.add_argument("--alpha", type=float, default=0.05,
                   help="per-dim significance alpha for the cell stats")
    p.add_argument("--fdr-q", type=float, default=0.05)
    return p.parse_args()


def _load(tag: str):
    pkl = config.RESULTS_DIR / f"{tag}.pkl"
    if not pkl.is_file():
        raise SystemExit(f"missing {pkl}")
    with open(pkl, "rb") as f:
        raw = pickle.load(f)
    return [r for r in raw["results"] if hasattr(r, "cells")]


def _records_by_pair(fits, alpha):
    """(area_x, area_y) -> list of {animal_id, landmark, epoch, mi_sig, ifi_weighted}."""
    by_pair: dict = defaultdict(list)
    for r in fits:
        for (epoch, lm), cell in r.cells.items():
            s = subspace_stats.cell_subspace_stats(cell, alpha=alpha)
            if s is None:
                continue
            by_pair[(r.area_x, r.area_y)].append({
                "animal_id": r.animal_id,
                "landmark": int(lm),
                "epoch": epoch,
                "mi_sig": s.mi_sig,
                "ifi_weighted": s.ifi_weighted,
            })
    return by_pair


def main():
    args = parse_args()
    fits = _load(args.tag)
    by_pair = _records_by_pair(fits, args.alpha)
    pairs = [p for p in config.PAIRS if p in by_pair]

    # Two honest estimators: the per-animal collapse (robust, always estimable,
    # underpowered) and the random-slope LMM (uses landmark structure, only
    # where there are enough animals). Report both; lead with the collapse.
    methods = {"collapsed": mixed_effects.collapsed_epoch_contrasts,
               "lmm_randslope": mixed_effects.lmm_epoch_contrasts}
    ordered = [f"{a}-{b}" for a, b in CONTRASTS]
    rows = []
    for pair in pairs:
        recs = by_pair[pair]
        for method, fn in methods.items():
            for stat in STATS:
                kw = {"value": stat, "contrasts": CONTRASTS}
                if method == "lmm_randslope":
                    kw["landmark"] = "landmark"
                res = fn(recs, **kw)
                pvals = np.array([res[c]["p"] for c in ordered])
                passed = paired_stats.fdr_bh(pvals, q=args.fdr_q)
                for c, ok in zip(ordered, passed):
                    d = res[c]
                    rows.append({
                        "method": method, "pair": f"{pair[0]}-{pair[1]}",
                        "stat": stat, "contrast": c,
                        "n_animals": d["n_animals"], "n_obs": d["n_obs"],
                        "estimate": d["estimate"], "p": d["p"],
                        "p_within_pair_fdr_pass": bool(ok),
                        "ok": d["ok"], "reason": d["reason"],
                    })

    out_csv = config.RESULTS_DIR / f"learning_changes_mixed_{args.tag}.csv"
    fields = ["method", "pair", "stat", "contrast", "n_animals", "n_obs",
              "estimate", "p", "p_within_pair_fdr_pass", "ok", "reason"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out_csv}  ({len(rows)} rows)\n")

    for method in methods:
        surv = [r for r in rows
                if r["method"] == method and r["p_within_pair_fdr_pass"]]
        label = ("PER-ANIMAL COLLAPSE (robust, n=animals)" if method == "collapsed"
                 else "RANDOM-SLOPE LMM (where N permits)")
        print(f"[{label}]")
        if surv:
            for r in sorted(surv, key=lambda x: x["p"]):
                print(f"  {r['pair']:9s} {r['stat']:12s} {r['contrast']:20s} "
                      f"est={r['estimate']:+.4f}  p={r['p']:.4g}  "
                      f"(n_animals={r['n_animals']})")
        else:
            print("  no within-pair-FDR-significant learning effect.")
        if method == "lmm_randslope":
            nbad = sum(1 for r in rows if r["method"] == method and not r["ok"])
            print(f"  ({nbad}/{len([r for r in rows if r['method']==method])} "
                  "contrasts unestimable — small-N pairs can't fit a random slope.)")
        print()


if __name__ == "__main__":
    main()
