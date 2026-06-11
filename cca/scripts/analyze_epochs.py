"""Early / intermediate / post-LP epoch analysis of the full subspace suite.

Reads results/epoch_metrics.csv (from run_epochs.py) and, per area pair and
metric, contrasts the three learning epochs with BOTH nonparametric and
parametric tests:
  * paired Wilcoxon signed-rank (per-animal)
  * paired t-test (per-animal, parametric)
  * random-slope linear mixed model across animals (`mixed_effects.lmm_epoch_contrasts`)
for every pair (all relationships). Also renders an epoch figure: each metric x
pair across the 3 epochs, per-animal points + mean ± SEM + paired lines.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, mixed_effects, paired_stats  # noqa: E402
import figstyle  # noqa: E402

figstyle.apply()

ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")
EPOCHS = ["naive", "intermediate", "expert"]
CONTRASTS = [("expert", "naive"), ("expert", "intermediate"), ("intermediate", "naive")]
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
METRICS = ["cc1", "n_sig", "mi_sig", "ifi", "gini_x"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _per_animal(df, pair, metric):
    """animal -> {epoch: value} for one pair/metric."""
    out = defaultdict(dict)
    g = df[df["pair"] == pair]
    for _, r in g.iterrows():
        v = pd.to_numeric(r[metric], errors="coerce")
        if np.isfinite(v):
            out[r["animal"]][r["epoch"]] = float(v)
    return out


def _paired(table, a, b):
    pa, pb = [], []
    for an, d in table.items():
        if a in d and b in d:
            pa.append(d[a]); pb.append(d[b])
    return np.array(pa), np.array(pb)


def analyse(df, tag):
    print(f"\n{'='*78}\nEPOCH CONTRASTS — {tag} (W=Wilcoxon, t=paired t, LMM=random-slope)\n{'='*78}")
    for metric in METRICS:
        print(f"\n[{metric}]")
        for pair in PAIRS:
            table = _per_animal(df, pair, metric)
            if len(table) < 3:
                continue
            # LMM across epochs (records: one row per animal-epoch)
            recs = [{"animal_id": an, "epoch": e, "value": v}
                    for an, d in table.items() for e, v in d.items()]
            lmm = mixed_effects.lmm_epoch_contrasts(recs, landmark=None)
            line = f"  {pair:9s}"
            for a, b in CONTRASTS:
                pa, pb = _paired(table, a, b)
                if pa.size < 3:
                    line += f" | {a[:3]}-{b[:3]}: n<3"
                    continue
                d = pa - pb
                wp = paired_stats.wilcoxon_signed(d.tolist())[3]
                tp = float(stats.ttest_rel(pa, pb).pvalue)
                lp = lmm.get(f"{a}-{b}", {}).get("p", float("nan"))
                mk = lambda p: "*" if (np.isfinite(p) and p < 0.05) else ""
                line += (f" | {a[:3]}-{b[:3]} Δ={np.median(d):+.3f} "
                         f"W={wp:.2g}{mk(wp)} t={tp:.2g}{mk(tp)} L={lp:.2g}{mk(lp)}")
            print(line)


def fig_epochs(df, tag):
    fig, axes = plt.subplots(len(METRICS), len(PAIRS),
                             figsize=(2.0 * len(PAIRS), 2.0 * len(METRICS)),
                             squeeze=False)
    xs = np.arange(len(EPOCHS))
    for i, metric in enumerate(METRICS):
        for j, pair in enumerate(PAIRS):
            ax = axes[i][j]
            table = _per_animal(df, pair, metric)
            for an, d in table.items():
                ys = [d.get(e, np.nan) for e in EPOCHS]
                ax.plot(xs, ys, "-", color="#bdc3c7", lw=0.6, alpha=0.6, zorder=1)
            means, sems = [], []
            for e in EPOCHS:
                v = np.array([d[e] for d in table.values() if e in d])
                means.append(np.mean(v) if v.size else np.nan)
                sems.append(np.std(v, ddof=1) / np.sqrt(v.size) if v.size > 1 else 0.0)
            ax.errorbar(xs, means, yerr=sems, color="#c0392b", lw=2, capsize=3, zorder=3)
            ax.set_xticks(xs)
            ax.set_xticklabels(["nv", "int", "exp"] if i == len(METRICS) - 1 else [],
                               fontsize=7)
            if j == 0:
                ax.set_ylabel(metric, fontsize=8)
            if i == 0:
                ax.set_title(pair, fontsize=8)
    fig.suptitle(f"Subspace metrics across learning epochs — {tag} "
                 "(red = mean ± SEM, faint = animals)", fontsize=12)
    ATT.mkdir(parents=True, exist_ok=True)
    figstyle.save(fig, ATT / f"HCV1_CCA_{tag}_epochs.png")


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "epoch_metrics.csv"
    tag = "fsincl" if "fsincl" in name else "fsexcl"
    df = pd.read_csv(config.RESULTS_DIR / name)
    print(f"{name}: {len(df)} rows, {df['animal'].nunique()} learners")
    analyse(df, tag)
    fig_epochs(df, tag)
    print(f"\nwrote epoch figure -> {ATT}/HCV1_CCA_{tag}_epochs.png")


if __name__ == "__main__":
    main()
