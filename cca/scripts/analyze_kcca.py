"""Analyse the nonlinear (kernel) CCA outputs — strength, direction, membership,
and the KCCA-vs-linear superiority test, with honest session-as-unit stats.

Reads results/kcca_metrics{,_fsincl}.csv (scope = session | naive | intermediate |
expert) and results/kcca_transition{,_fsincl}.csv (cued vs uncued). Reports:

  1. SESSION levels (all animals): per pair held-out KCCA vs LINEAR CC (paired
     Wilcoxon — the superiority test), significance rate, IFI direction, Gini.
  2. SUPERIORITY summary across every session cell.
  3. EPOCH learning contrasts (learners): naive/intermediate/expert deltas in
     KCCA CC, IFI, Gini — Wilcoxon + random-slope LMM (vs the linear epoch result).
  4. LEARNER vs NON-LEARNER session levels.
  5. CUED vs UNCUED deltas.

Re-runnable in seconds (no re-fitting). Run with no args for FS-excluded; pass
``fsincl`` to read the FS-included tables.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, mixed_effects, paired_stats  # noqa: E402

PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
EPOCHS = ["naive", "intermediate", "expert"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _per_animal(df, pair, col, scope=None):
    g = df[df["pair"] == pair]
    if scope is not None:
        g = g[g["scope"] == scope]
    return [float(_num(s[col]).mean()) for _, s in g.groupby("animal")]


def _w(vals):
    vals = [v for v in vals if np.isfinite(v)]
    if len(vals) < 3:
        return np.nan, np.nan, len(vals)
    _, med, _, p = paired_stats.wilcoxon_signed(vals)
    return med, p, len(vals)


def session_levels(df):
    print("=" * 78)
    print("1. SESSION LEVELS (all animals) — KCCA vs LINEAR held-out CC1 "
          "(paired Wilcoxon)")
    print("=" * 78)
    print(f"{'pair':9s} {'n':>3} {'kcca':>7} {'linear':>7} {'Δ(k-l)':>7} "
          f"{'p_sup':>7} {'n_sig%':>7} {'IFI':>7} {'p_ifi':>6} {'gini':>6}")
    sess = df[df["scope"] == "session"]
    for pair in PAIRS:
        g = sess[sess["pair"] == pair]
        if g.empty:
            continue
        k = _per_animal(sess, pair, "cc_kcca"); l = _per_animal(sess, pair, "cc_linear")
        paired = [ki - li for ki, li in zip(k, l) if np.isfinite(ki) and np.isfinite(li)]
        _, p_sup, n = _w(paired)
        ifi = _per_animal(sess, pair, "ifi"); _, p_ifi, _ = _w(ifi)
        nsig_rate = float(np.mean(_num(g["n_sig"]) > 0))
        gini = np.nanmean(_per_animal(sess, pair, "gini_x"))
        print(f"{pair:9s} {n:>3} {np.nanmean(k):>7.3f} {np.nanmean(l):>7.3f} "
              f"{np.nanmean(paired):>+7.3f} {p_sup:>7.3g} {nsig_rate:>7.0%} "
              f"{np.nanmean(ifi):>+7.3f} {p_ifi:>6.2g} {gini:>6.3f}")


def superiority(df):
    sess = df[df["scope"] == "session"].copy()
    d = (_num(sess["cc_kcca"]) - _num(sess["cc_linear"])).dropna()
    win = int((d > 0).sum())
    print("\n" + "=" * 78)
    print("2. SUPERIORITY (every session cell): does KCCA beat LINEAR out-of-sample?")
    print("=" * 78)
    print(f"  cells: {len(d)} | median Δ(KCCA−linear) = {d.median():+.3f} "
          f"| KCCA>linear in {win}/{len(d)} ({win/len(d):.0%})")
    print(f"  median KCCA cc1 {_num(sess['cc_kcca']).median():.3f} vs "
          f"linear {_num(sess['cc_linear']).median():.3f}")


def epoch_contrasts(df):
    print("\n" + "=" * 78)
    print("3. EPOCH LEARNING CONTRASTS (learners) — Δ across naive→expert "
          "(Wilcoxon | LMM)")
    print("=" * 78)
    ep = df[df["scope"].isin(EPOCHS)].copy()
    for metric in ["cc_kcca", "ifi", "gini_x"]:
        print(f"\n[{metric}]")
        for pair in PAIRS:
            g = ep[ep["pair"] == pair]
            if g["animal"].nunique() < 3:
                continue
            recs = [{"animal_id": r["animal"], "epoch": r["scope"],
                     "value": float(_num(pd.Series([r[metric]]))[0])}
                    for _, r in g.iterrows()]
            coll = mixed_effects.collapsed_epoch_contrasts(recs)
            lmm = mixed_effects.lmm_epoch_contrasts(recs, landmark=None)
            cells = []
            for c in ("expert-naive", "intermediate-naive", "expert-intermediate"):
                cv, lv = coll.get(c, {}), lmm.get(c, {})
                est = cv.get("estimate", np.nan); wp = cv.get("p", np.nan)
                lp = lv.get("p", np.nan) if lv.get("ok") else np.nan
                star = "*" if (np.isfinite(wp) and wp < 0.05) else ""
                cells.append(f"{c.split('-')[0][:3]}-{c.split('-')[1][:3]} "
                             f"Δ={est:+.3f} W={wp:.2g}{star} L={lp:.2g}")
            print(f"  {pair:9s} n={g['animal'].nunique():<2d} | " + " | ".join(cells))


def learner_split(df):
    print("\n" + "=" * 78)
    print("4. LEARNER vs NON-LEARNER (session held-out KCCA CC1)")
    print("=" * 78)
    sess = df[df["scope"] == "session"]
    for pair in PAIRS:
        out = []
        for lv, lab in ((1, "L"), (0, "n")):
            sub = sess[(sess["pair"] == pair) & (sess["learner"] == lv)]
            v = _num(sub["cc_kcca"]).dropna()
            out.append(f"{lab}: {v.mean():.3f} (n={len(v)})" if len(v) else f"{lab}: -")
        print(f"  {pair:9s} " + "  ".join(out))


def transition(tag):
    path = config.RESULTS_DIR / f"kcca_transition{tag}.csv"
    if not path.is_file():
        print(f"\n(no {path.name})")
        return
    t = pd.read_csv(path)
    print("\n" + "=" * 78)
    print(f"5. CUED vs UNCUED (Δ = cued − uncued; learners) — {path.name}")
    print("=" * 78)
    learn = t[t["learner"] == 1]
    for pair in PAIRS:
        g = learn[learn["pair"] == pair]
        if len(g) < 3:
            continue
        cells = []
        for mt in ["cc_kcca", "ifi", "gini_x"]:
            d = _num(g[f"d_{mt}"]).dropna().tolist()
            _, med, _, p = paired_stats.wilcoxon_signed(d) if len(d) >= 3 else (0, np.nan, 0, np.nan)
            star = "*" if (np.isfinite(p) and p < 0.05) else ""
            cells.append(f"d_{mt}={med:+.3f} p={p:.2g}{star}")
        print(f"  {pair:9s} n={len(g):<2d} | " + " | ".join(cells))


def main():
    tag = "_fsincl" if (len(sys.argv) > 1 and "fsincl" in sys.argv[1]) else ""
    path = config.RESULTS_DIR / f"kcca_metrics{tag}.csv"
    if not path.is_file():
        print(f"missing {path}"); return
    df = pd.read_csv(path)
    print(f"{path.name}: {len(df)} rows, {df['animal'].nunique()} animals "
          f"({df[df['learner']==1]['animal'].nunique()} learners), "
          f"scopes={sorted(df['scope'].unique())}\n")
    session_levels(df)
    superiority(df)
    epoch_contrasts(df)
    learner_split(df)
    transition(tag)


if __name__ == "__main__":
    main()
