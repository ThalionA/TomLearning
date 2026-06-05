"""Write the landmark config-pruning CSVs.

Scans ``results/`` for every done landmark config, pools the per-(config, pair,
epoch) significant-subspace-dimension values with the *same* aggregation the
sweep figures use (``summarise_landmark_sweep._per_pair_stats``), and emits two
flat tables via ``tom_cca.prune_table``:

  results/landmark_prune_table.csv     one row per (config, pair); wide over
                                       epochs; held-out CC / IFI / dimensionality
                                       metrics, overfitting diagnostics
                                       (max_cc, frac_cc_ge_099), the six headline
                                       p-values (CC & IFI: naive/expert vs 0,
                                       naive vs expert), and an overfit_flag.
  results/landmark_prune_summary.csv   one row per config; rollup over pairs:
                                       how many pairs overfit, worst held-out CC
                                       anywhere, typical expert strength, and how
                                       many pairs show significant expert CC /
                                       a significant naive->expert increase.

No re-fitting: operates purely on existing pkls. The p-values reproduce exactly
what ``figures/landmark_sweep/sweep_pvalues_summary.png`` draws.

Run:  PYTHONPATH=src python scripts/build_prune_table.py
"""

from __future__ import annotations

import argparse
import csv
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, prune_table  # noqa: E402

# Reuse the sweep's tested aggregation rather than re-deriving it.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarise_landmark_sweep import (  # noqa: E402
    _find_done_configs, _parse_tag, _per_pair_stats,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--alpha", type=float, default=0.05,
                   help="per-dim significance level for sig-dim selection AND "
                        "the keep/drop p-value threshold in the rollup")
    return p.parse_args()


def main():
    args = parse_args()
    configs = _find_done_configs()
    if not configs:
        raise SystemExit("no landmark*.done files in results/")
    print(f"found {len(configs)} done configs")

    per_cpe: dict = {}
    for tag, pkl in configs:
        if _parse_tag(tag) is None:
            print(f"  skip (unparseable tag): {tag}")
            continue
        with open(pkl, "rb") as f:
            raw = pickle.load(f)
        fits = [r for r in raw["results"] if hasattr(r, "cells")]
        by_pair = _per_pair_stats(fits, alpha=args.alpha)
        for pair, by_epoch in by_pair.items():
            for epoch in prune_table.EPOCHS:
                per_cpe[(tag, pair, epoch)] = by_epoch[epoch]
        print(f"  {tag}: {len(fits)} fits, {len(by_pair)} pairs")

    rows = prune_table.build_prune_rows(per_cpe)
    roll = prune_table.rollup_rows(rows, alpha=args.alpha)

    out_detail = config.RESULTS_DIR / "landmark_prune_table.csv"
    with open(out_detail, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=prune_table.prune_fieldnames(),
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_detail}  ({len(rows)} rows)")

    out_summary = config.RESULTS_DIR / "landmark_prune_summary.csv"
    with open(out_summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=prune_table.rollup_fieldnames(),
                           lineterminator="\n")
        w.writeheader()
        w.writerows(roll)
    print(f"wrote {out_summary}  ({len(roll)} rows)")

    n_overfit = sum(1 for r in roll if r["frac_overfit_pairs"] > 0)
    print(f"\nconfigs with >=1 overfit pair: {n_overfit}/{len(roll)}")
    flagged = sorted(roll, key=lambda r: -r["frac_overfit_pairs"])
    for r in flagged:
        if r["frac_overfit_pairs"] > 0:
            print(f"  {r['tag']:24s} overfit_pairs={r['n_overfit_pairs']}/"
                  f"{r['n_pairs']}  max_cc={r['max_cc_any']:.3f}")


if __name__ == "__main__":
    main()
