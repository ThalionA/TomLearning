"""Does early communication predict how fast an animal learns?

EXPLORATORY (pre-registered expectation: null at this n). Predictor = the
NAIVE-epoch (first 10 running trials) communication readouts from
results/epoch_metrics_bin10*.csv — measured BEFORE learning, so temporal
precedence is clean:
  * cc1  — held-out CC1 (coupling strength)
  * ifi  — information-flow index (direction; +ve = first-named area leads)
Outcome = the learning point LP (trial index at which the animal reaches
expert criterion; lower = faster learner), learners only (n = 12; per pair
4-12 depending on coverage).

Tests: Spearman rho(predictor, LP) per (pair, metric); BH within metric
family. GLOBAL strength predictor: cc1 z-scored WITHIN pair across animals
(pairs have very different coupling levels), then averaged per animal — one
value per animal, one test. Direction has no meaningful cross-pair pooling
(signs are pair-specific conventions) — per-pair only.

Confound: an animal that runs more in the naive epoch may both couple better
and learn faster — naive n_bins is tested against LP alongside.

Writes results/conn_learning_bin10{,_fsincl}.csv and prints the tables.
Usage: PYTHONPATH=src python scripts/analyze_conn_learning.py [--include-fs]
"""
from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, dataio  # noqa: E402

PAIRS = ["CA1-CA3", "CA1-DG", "CA3-DG", "CA1-SUB", "CA1-RSC", "CA1-V1",
         "RSC-SUB", "V1-RSC"]
MIN_N = 5


def bh(pvals):
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    out = np.zeros(p.size, dtype=bool)
    if ok.sum():
        out[ok] = stats.false_discovery_control(p[ok]) < 0.05
    return out


def spear(x, y):
    r = stats.spearmanr(x, y)
    return float(r.statistic), float(r.pvalue)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-fs", action="store_true")
    args = ap.parse_args()
    tag = "_fsincl" if args.include_fs else ""
    em = pd.read_csv(config.RESULTS_DIR / f"epoch_metrics_bin10{tag}.csv")

    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=10,
                              gaussian_sd_ms=2.5)
    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(
        animals, cfg, behaviour_lookup=dataio.load_learning_points())
    lp = {a: e.lp for a, e in entries.items()
          if e.role == "learner" and e.lp is not None}

    naive = em[(em["epoch"] == "naive") & (em["learner"] == 1)].copy()
    naive["lp"] = naive["animal"].map(lp)
    naive = naive[np.isfinite(naive["lp"])]
    print(f"naive-epoch predictors vs LP "
          f"({'FS incl' if args.include_fs else 'FS excl'}): "
          f"{naive['animal'].nunique()} learners, LP range "
          f"{int(naive['lp'].min())}-{int(naive['lp'].max())}")

    rows = []
    for metric in ("cc1", "ifi"):
        print(f"\n[naive {metric} vs LP]  (negative rho = stronger/more-directed "
              f"-> FASTER learning)")
        cells, pv = [], []
        for pair in PAIRS:
            g = naive[naive["pair"] == pair][["animal", metric, "lp"]].dropna()
            if len(g) < MIN_N:
                continue
            rho, p = spear(g[metric], g["lp"])
            cells.append((pair, len(g), rho, p)); pv.append(p)
        for (pair, n, rho, p), b in zip(cells, bh(pv)):
            mark = "*" if p < 0.05 else " "
            print(f"  {pair:9s} n={n:<2d} rho={rho:+.3f} p={p:.3g}{mark}"
                  f"{'  BH' if b else ''}")
            rows.append(dict(metric=metric, pair=pair, n=n,
                             rho=round(rho, 4), p=round(p, 5), bh_pass=bool(b)))

    # global strength: z-score cc1 within pair, mean per animal
    z = naive.copy()
    z["cc1_z"] = z.groupby("pair")["cc1"].transform(
        lambda v: (v - v.mean()) / v.std(ddof=1))
    ga = z.groupby("animal").agg(cc1_z=("cc1_z", "mean"), lp=("lp", "first"),
                                 n_pairs=("pair", "nunique"))
    rho, p = spear(ga["cc1_z"], ga["lp"])
    print(f"\n[GLOBAL naive coupling (within-pair z of cc1, mean per animal) vs LP] "
          f"n={len(ga)} rho={rho:+.3f} p={p:.3g}"
          f"{' *' if p < 0.05 else ''}")
    rows.append(dict(metric="GLOBAL_cc1_z", pair="ALL", n=len(ga),
                     rho=round(rho, 4), p=round(p, 5), bh_pass=p < 0.05))

    # behaviour confound: naive running amount vs LP
    beh = naive.groupby("animal").agg(n_bins=("n_bins", "first"),
                                      lp=("lp", "first"))
    rho, p = spear(beh["n_bins"], beh["lp"])
    print(f"[CONFOUND naive n_bins (running amount) vs LP] n={len(beh)} "
          f"rho={rho:+.3f} p={p:.3g}")
    rows.append(dict(metric="CONFOUND_n_bins", pair="ALL", n=len(beh),
                     rho=round(rho, 4), p=round(p, 5), bh_pass=False))

    out = pd.DataFrame(rows)
    path = config.RESULTS_DIR / f"conn_learning_bin10{tag}.csv"
    out.to_csv(path, index=False)
    print(f"\nwrote {path.name} ({len(out)} rows; "
          f"{int(out['bh_pass'].sum())} BH/global survivors)")


if __name__ == "__main__":
    main()
