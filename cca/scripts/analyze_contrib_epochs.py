"""Do the contribution × spatial-property links change with experience?

Epoch contrasts (naive = first 10 running trials, intermediate = 10 trials
ending at the learning point, expert = 10 after it — `dataio.epoch_windows`,
`trials_per_epoch = 10`) for every metric in the contribution thread:

  * per (pair, area): contrib_conn vs reliability / tuning (raw + rate-part)
    from results/contrib_reliability_bin10*.csv
  * per area, pooled over partners: same four metrics
    from results/contrib_pooled_bin10*.csv
  * per area, cross-partner overlap: rho_conn, rho_resid, jaccard_excess
    from results/contrib_overlap_bin10*.csv (Fisher-z mean over partner pairs
    within epoch; jaccard_excess plain mean)

Per animal and cell the epoch value is Fisher-z(rho) (jaccard_excess raw);
contrasts are WITHIN-animal deltas (expert-naive, intermediate-naive,
expert-intermediate), tested animals-as-n with paired t + Wilcoxon vs 0.

⚠ Power: n = 4-10 animals per cell and the naive CCA fit rests on the fewest
running bins (trial 1 is longest/slowest — PROJECT_LOG 2026-08-20), so naive
rhos are the noisiest. Nulls here are weak evidence of stability.

Writes results/contrib_epochs_deltas_bin10{,_fsincl}.csv and prints the tables.

Usage: PYTHONPATH=src python scripts/analyze_contrib_epochs.py [--include-fs]
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

EPOCHS = ["naive", "intermediate", "expert"]
CONTRASTS = [("expert", "naive"), ("intermediate", "naive"),
             ("expert", "intermediate")]
PAIR_METRICS = ["rho_contrib_conn", "rho_contrib_conn_ratepart",
                "rho_contrib_conn_tune", "rho_contrib_conn_tune_ratepart"]
POOLED_METRICS = ["rho_pooled", "rho_pooled_ratepart",
                  "rho_pooled_tune", "rho_pooled_tune_ratepart"]
OVERLAP_METRICS = ["rho_conn", "rho_resid", "jaccard_excess"]


def epoch_values(df, metric, cell_cols):
    """(cell..., animal, epoch) -> value; Fisher-z for rho metrics, raw for
    jaccard. Overlap rows repeat per partner-pair -> mean within epoch."""
    v = df[metric].astype(float)
    if metric != "jaccard_excess":
        v = np.arctanh(v.clip(-0.999, 0.999))
    keys = [df[c] for c in cell_cols] + [df["animal"], df["epoch"]]
    return v.groupby(keys).mean()


def deltas(ev, cell, hi, lo):
    """Per-animal hi-epoch minus lo-epoch for one cell; both epochs required."""
    try:
        sub = ev.loc[cell]
    except KeyError:
        return np.array([])
    piv = sub.unstack("epoch")
    if hi not in piv.columns or lo not in piv.columns:
        return np.array([])
    d = (piv[hi] - piv[lo]).to_numpy(float)
    return d[np.isfinite(d)]


def report(ev, cells, metric, rows_out, fam):
    print(f"\n[{fam}: {metric}]  (Δ in Fisher-z units"
          f"{' — raw Δ' if metric == 'jaccard_excess' else ''})")
    for cell in cells:
        name = cell if isinstance(cell, str) else " ".join(cell)
        line = f"  {name:14s}"
        printable = False
        for hi, lo in CONTRASTS:
            d = deltas(ev, cell, hi, lo)
            if d.size < 3:
                line += f" | {hi[:3]}-{lo[:3]}: {'n<3':>22s}"
                continue
            printable = True
            t_p = float(stats.ttest_1samp(d, 0.0).pvalue)
            _, med, _, w_p = paired_stats.wilcoxon_signed(d.tolist())
            mark = "*" if (np.isfinite(w_p) and w_p < 0.05) else " "
            tmark = "*" if (np.isfinite(t_p) and t_p < 0.05) else " "
            line += (f" | {hi[:3]}-{lo[:3]}: n={d.size} med={med:+.3f} "
                     f"t{tmark}{t_p:.3g} W{mark}{w_p:.3g}")
            rows_out.append(dict(family=fam, metric=metric, cell=name,
                                 contrast=f"{hi}-{lo}", n=d.size,
                                 med_delta=round(float(med), 4),
                                 t_p=round(t_p, 5), w_p=round(w_p, 5)))
        if printable:
            print(line)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-fs", action="store_true")
    args = ap.parse_args()
    tag = "_fsincl" if args.include_fs else ""
    rel = pd.read_csv(config.RESULTS_DIR / f"contrib_reliability_bin10{tag}.csv")
    pooled = pd.read_csv(config.RESULTS_DIR / f"contrib_pooled_bin10{tag}.csv")
    overlap = pd.read_csv(config.RESULTS_DIR / f"contrib_overlap_bin10{tag}.csv")
    print(f"epoch contrasts ({'FS incl' if args.include_fs else 'FS excl'}); "
          f"naive = first 10 running trials, intermediate = 10 ending at LP, "
          f"expert = 10 after LP")

    rows_out = []
    pair_cells = [(p, a) for p in sorted(rel["pair"].unique())
                  for a in p.split("-")]
    for metric in PAIR_METRICS:
        if metric in rel.columns:
            ev = epoch_values(rel, metric, ["pair", "area"])
            report(ev, pair_cells, metric, rows_out, "pair")
    areas = [a for a in ["CA1", "CA3", "DG", "SUB", "RSC", "V1"]
             if a in set(pooled["area"])]
    for metric in POOLED_METRICS:
        if metric in pooled.columns:
            report(epoch_values(pooled, metric, ["area"]), areas, metric,
                   rows_out, "pooled")
    for metric in OVERLAP_METRICS:
        report(epoch_values(overlap, metric, ["area"]), areas, metric,
               rows_out, "overlap")

    out = pd.DataFrame(rows_out)
    path = config.RESULTS_DIR / f"contrib_epochs_deltas_bin10{tag}.csv"
    out.to_csv(path, index=False)
    n_sig = int((out["w_p"] < 0.05).sum())
    print(f"\nwrote {path.name} ({len(out)} contrast-rows; "
          f"{n_sig} with Wilcoxon p<0.05 of {len(out)} tests, uncorrected)")


if __name__ == "__main__":
    main()
