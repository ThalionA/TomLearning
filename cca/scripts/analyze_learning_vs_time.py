"""Is the Gini↓ (de-sparsification) effect LEARNING-specific or just TIME-on-task?

The one confound on the headline result: Gini also drifts down weakly in the 4
non-learners. Three within-data tests, per pair, on the de-sparsification metric
(`gini_x` vs the trial-fraction axis), from results/trajectory_windows.csv:

  1. Learners vs non-learners — per-animal slope distributions + Mann–Whitney
     between groups (is the drop steeper in learners?).
  2. Interaction LMM — gini_x ~ trial_frac * learner + (trial_frac | animal)
     across all 16 animals; the trial_frac×learner term tests whether the slope
     differs by group (the direct learning-vs-time test, pooling all windows).
  3. Post-LP plateau (learners only) — Gini slope on PRE-LP (lp_rel<0) vs
     POST-LP (lp_rel≥0) windows. Time-on-task ⇒ keeps dropping post-LP;
     learning ⇒ saturates at the LP (post-LP slope ≈ 0). This needs no
     non-learners, so it sidesteps the n=4 problem.
"""

from __future__ import annotations

import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, paired_stats, trajectory  # noqa: E402

PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
METRIC = "gini_x"


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _slopes(df, axis):
    out = []
    for an, g in df.groupby("animal"):
        x = _num(g[axis]).to_numpy(); y = _num(g[METRIC]).to_numpy()
        if np.sum(np.isfinite(x) & np.isfinite(y)) >= 4:
            s, _ = trajectory.linear_slope(x, y)
            if np.isfinite(s):
                out.append(s)
    return out


def _interaction_p(dfp):
    """trial_frac×learner interaction from a random-slope LMM over all animals."""
    d = dfp.dropna(subset=[METRIC, "trial_frac", "learner", "animal"]).copy()
    if d["animal"].nunique() < 6 or d["learner"].nunique() < 2:
        return np.nan, np.nan
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = smf.mixedlm(f"{METRIC} ~ trial_frac * learner", d,
                            groups=d["animal"], re_formula="~trial_frac")
            f = m.fit(reml=False, method="lbfgs")
    except Exception:
        return np.nan, np.nan
    term = next((t for t in f.params.index if "trial_frac:" in t), None)
    if term is None:
        return np.nan, np.nan
    return float(f.params[term]), float(f.pvalues[term])


def main():
    df = pd.read_csv(config.RESULTS_DIR / "trajectory_windows.csv")
    learn = df[df["learner"] == 1]
    nonl = df[df["learner"] == 0]
    print(f"Learning vs time-on-task — metric={METRIC}, "
          f"{learn['animal'].nunique()} learners, {nonl['animal'].nunique()} non-learners\n")
    for pair in PAIRS:
        gl = learn[learn["pair"] == pair]
        gn = nonl[nonl["pair"] == pair]
        slL = _slopes(gl, "trial_frac")
        slN = _slopes(gn, "trial_frac")
        if len(slL) < 3:
            continue
        # 1. learner vs non-learner slopes
        line = f"### {pair}\n  slopes: learners n={len(slL)} med={np.median(slL):+.4f}"
        if len(slN) >= 1:
            line += f" | non-learners n={len(slN)} med={np.median(slN):+.4f}"
        if len(slL) >= 3 and len(slN) >= 3:
            u = float(stats.mannwhitneyu(slL, slN, alternative="two-sided").pvalue)
            line += f" | between-group U-p={u:.3g}"
        print(line)
        # 2. interaction LMM
        b, p = _interaction_p(df[df["pair"] == pair])
        print(f"  interaction LMM (trial_frac×learner): β={b:+.4f} p={p:.3g}"
              + ("  *" if (np.isfinite(p) and p < 0.05) else ""))
        # 3. post-LP plateau (learners): pre-LP vs post-LP slope
        pre = _slopes(gl[_num(gl["lp_rel"]) < 0], "lp_rel")
        post = _slopes(gl[_num(gl["lp_rel"]) >= 0], "lp_rel")
        def _sgn(sl):
            if len(sl) < 3:
                return f"n={len(sl)} (too few)"
            _, med, _, pp = paired_stats.wilcoxon_signed(sl)
            return f"n={len(sl)} med={med:+.5f} p={pp:.3g}" + ("*" if pp < 0.05 else "")
        print(f"  plateau: pre-LP slope {_sgn(pre)} | post-LP slope {_sgn(post)}")
        print()


if __name__ == "__main__":
    main()
