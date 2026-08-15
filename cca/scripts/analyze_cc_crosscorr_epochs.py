"""2026-08-07 meeting ask 2 — cross-correlograms of naive and expert, whole and split
by positive/negative-grouped (FF/FB) CCs.

"Cross-correlogram" here = the lagged correlation curve r(lag) of the canonical
variates, −250..+250 ms, on FROZEN axes: one CCA fit per (animal, pair) on every
running bin of every trial, then each epoch projected through identical weights
(`run_cc_label_track.py`, `results/cc_label_track_bin10{,_fsincl}.csv`). Frozen axes
are the only family in which a canonical dimension is the SAME component across
epochs (HANDOFF.md §0), which is what "naive vs expert" needs.

Two aggregations, both per-ANIMAL-first (an animal with ten CCs counts once):
    group = "all"   every CC that beats the frozen fit's per-dim null
    group = "FF"/"FB"  the same CCs split by their whole-session label (sign of the
                    IFI at ±50 ms on the all-trials curve — identical to the label on
                    slide 10 of CCA_update_20260807.pptx; `--loo` switches to the
                    leave-epoch-out label)

⚠ `r` is IN-SAMPLE for the whole session (the frozen fit saw every epoch). It is a
contrast statistic across epochs, never a coupling strength — do not compare its
level with the held-out numbers elsewhere.

⚠ In the uncapped frozen fit almost every leading dim is "significant" (470/475
cells FS-excl): the gate is nearly a no-op here, unlike in the refit-per-lag arm.

⚠ LEVEL contrasts against the naive epoch are biased on this fit. Found while
building this (2026-08-15): mean peak r is naive 0.052 → intermediate 0.067 →
expert 0.068 — intermediate (pre-LP) already equals expert, and expert beats
intermediate in only 56 % of animal-pairs. The first 10 trials are uniquely low
because an all-trials fit represents the rest of the session better than the
session's atypical opening; the BALANCED-trial frozen fit (`fixed_subspace_stats`)
shows no such rise (Δpeak r mixed-sign, all n.s.). So: trust the IFI (a shape
ratio, null here) for naive-vs-expert; read peak-r / curve-height differences as
"first trials are different", not learning. The stats file carries
intermediate − naive next to expert − naive as the built-in diagnostic.

Writes:
    results/cc_crosscorr_epochs_bin10{,_fsincl}.csv        per (animal, pair, epoch, group, lag)
    results/cc_crosscorr_epochs_stats_bin10{,_fsincl}.csv  per pair: expert−naive, group=all
    results/cc_crosscorr_epochs_tables.md

Usage: PYTHONPATH=src python scripts/analyze_cc_crosscorr_epochs.py [--loo]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import cc_aggregate, paired_stats  # noqa: E402

RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
EPOCHS = ["naive", "intermediate", "expert"]
GROUPS = ["all", "FF", "FB"]


def build_curves(df: pd.DataFrame, label_col: str = "label",
                 sig_only: bool = True) -> pd.DataFrame:
    """One row per (animal, pair, epoch, group, lag_ms): mean r over that animal's CCs.

    ``group`` is "all" (every significant CC) or the label value ("FF" / "FB") from
    ``label_col``. Cells with no CC in a group are absent, not NaN.
    """
    df = df.copy()
    df["r"] = pd.to_numeric(df["r"], errors="coerce")
    parts = []
    for group in GROUPS:
        sub = df if group == "all" else df[df[label_col] == group]
        if sub.empty:
            continue
        pa = cc_aggregate.per_animal_mean(sub, value="r", by=["epoch", "lag_ms"],
                                          sig_only=sig_only, drop_degenerate=False)
        pa.insert(2, "group", group)
        parts.append(pa)
    out = pd.concat(parts, ignore_index=True).rename(columns={"mean": "r_mean"})
    return out[["animal", "pair", "group", "epoch", "lag_ms", "r_mean", "n_ccs"]]


def cross_animal(curves: pd.DataFrame, pair: str, epoch: str, group: str):
    """Mean ± SEM across animals of the per-animal curves for one panel line.

    Returns ``(lags_ms, mean, sem, n_animals)``; empty arrays and n = 0 if no data.
    """
    sub = curves[(curves["pair"] == pair) & (curves["epoch"] == epoch) &
                 (curves["group"] == group)]
    if sub.empty:
        return np.array([]), np.array([]), np.array([]), 0
    piv = sub.pivot_table(index="animal", columns="lag_ms", values="r_mean")
    lags = np.array(sorted(piv.columns), dtype=float)
    M = piv[sorted(piv.columns)].to_numpy(float)
    n = np.sum(np.isfinite(M), axis=0)
    mean = np.nanmean(M, axis=0)
    sd = np.nanstd(M, axis=0, ddof=1)
    sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
    return lags, mean, sem, int(piv.shape[0])


def epoch_contrast(reduced: pd.DataFrame, metric: str, pairs=PAIRS,
                   a: str = "naive", b: str = "expert", min_n: int = 3) -> pd.DataFrame:
    """Per pair: paired *t* of (``b`` − ``a``) on the per-animal mean of ``metric``
    over all significant CCs. ``reduced`` is one row per (animal, pair, dim, epoch)
    — `cc_label_track_epoch_bin10*.csv` from analyze_cc_label_track.py — so the
    per-CC reductions (IFI at ±50 ms, peak r) are the ones already on disk.

    Animals lacking either epoch are dropped. Eight pairs = eight families
    (project policy); ``bh_pass`` across pairs is a sensitivity column.
    """
    pa = cc_aggregate.per_animal_mean(reduced, value=metric, by=["epoch"],
                                      drop_degenerate=False)
    rows = []
    for pair in pairs:
        w = pa[pa["pair"] == pair].pivot_table(index="animal", columns="epoch",
                                               values="mean")
        if a not in w.columns or b not in w.columns:
            continue
        d = (w[b] - w[a]).to_numpy(float)
        d = d[np.isfinite(d)]
        if d.size < min_n:
            continue
        n, _, t, p = paired_stats.paired_t(d)
        rows.append({"pair": pair, "metric": metric, "contrast": f"{b}-{a}", "n": int(n),
                     "mean_a": float(np.nanmean(w[a])), "mean_b": float(np.nanmean(w[b])),
                     "delta": float(np.mean(d)), "t": t, "p": p})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["bh_pass"] = paired_stats.fdr_bh(out["p"].to_numpy(float))
    return out


def _md(stats: pd.DataFrame, curves: pd.DataFrame, fs: str, label_col: str) -> str:
    lines = [f"### {fs} — cross-correlograms by epoch, all significant CCs "
             f"(frozen axes, per-animal-first)", "",
             f"Label column for the FF/FB split: `{label_col}`. Paired *t* across "
             f"animals on each animal's mean over its significant CCs of the per-CC "
             f"reduction already in `cc_label_track_epoch_*` (IFI at ±50 ms; peak r). "
             f"Positive IFI ⇒ first-named area leads. Per-pair families; BH across the "
             f"8 pairs is a sensitivity column.", "",
             "> ⚠ **Read peak-r (level) contrasts against naive with the caveat in the "
             "script docstring:** on this all-trials fit the naive epoch is uniquely "
             "low and intermediate (pre-LP) already equals expert — `intermediate-naive` "
             "is printed next to `expert-naive` so that is visible. The balanced-trial "
             "frozen fit shows no peak-r rise. IFI (a shape ratio) is the trustworthy "
             "naive-vs-expert statistic here.", "",
             "| pair | n | metric | contrast | a | b | Δ (b − a) | t | p | BH |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in stats.iterrows():
        star = "**" if r["p"] < 0.05 else ""
        lines.append(f"| {r['pair']} | {int(r['n'])} | {r['metric']} | {r['contrast']} | "
                     f"{r['mean_a']:+.3f} | {r['mean_b']:+.3f} | "
                     f"{star}{r['delta']:+.3f}{star} | {r['t']:+.2f} | "
                     f"{star}{r['p']:.3g}{star} | {'yes' if bool(r['bh_pass']) else 'no'} |")
    # how many CCs / animals sit behind each panel line
    lines += ["", "**Animals per (pair, epoch) — group = all / FF / FB:**", "",
              "| pair | " + " | ".join(EPOCHS) + " |", "|---|" + "---|" * len(EPOCHS)]
    for pair in PAIRS:
        cells = []
        for ep in EPOCHS:
            ns = []
            for grp in GROUPS:
                s = curves[(curves["pair"] == pair) & (curves["epoch"] == ep) &
                           (curves["group"] == grp)]
                ns.append(str(s["animal"].nunique()))
            cells.append("/".join(ns))
        lines.append(f"| {pair} | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def main():
    loo = "--loo" in sys.argv[1:]
    label_col = "label_loo" if loo else "label"
    tag = "_loo" if loo else ""
    md = ["# 2026-08-07 ask 2 — cross-correlograms of naive and expert, whole and by "
          "FF/FB label", "",
          "Frozen-axes lag curves (`cc_label_track_bin10*.csv`), averaged over each "
          "animal's significant CCs first, then across animals. `r` is in-sample for "
          "the whole session: a contrast across epochs, not a coupling strength.", ""]
    for suf, fs in [("", "FS-excluded"), ("_fsincl", "FS-included")]:
        src = RES / f"cc_label_track_bin10{suf}.csv"
        red_src = RES / f"cc_label_track_epoch_bin10{suf}.csv"
        if not src.exists() or not red_src.exists():
            print(f"skip {fs}: {src.name} or {red_src.name} missing"); continue
        df = pd.read_csv(src)
        curves = build_curves(df, label_col=label_col)
        curves.to_csv(RES / f"cc_crosscorr_epochs{tag}_bin10{suf}.csv", index=False,
                      lineterminator="\n")
        red = pd.read_csv(red_src)
        # the epoch file has no `sig` column: it was written sig-gated already
        stats = pd.concat(
            [epoch_contrast(red, "ifi"), epoch_contrast(red, "peak_r"),
             epoch_contrast(red, "peak_r", a="naive", b="intermediate")],
            ignore_index=True)
        stats.to_csv(RES / f"cc_crosscorr_epochs_stats{tag}_bin10{suf}.csv",
                     index=False, lineterminator="\n")
        md.append(_md(stats, curves, fs, label_col))
        n_cells = df.groupby(["animal", "pair", "dim"]).ngroups
        n_sig = df[df["sig"] == 1].groupby(["animal", "pair", "dim"]).ngroups
        print(f"\n[{fs}] {n_sig}/{n_cells} (animal, pair, dim) cells significant; "
              f"{len(curves)} curve rows")
        print("  contrasts on the all-CC per-animal mean (paired t):")
        for _, r in stats.iterrows():
            flag = " *" if r["p"] < 0.05 else ""
            print(f"    {r['pair']:8s} {r['metric']:7s} {r['contrast']:19s} "
                  f"n={int(r['n']):<2d} Δ={r['delta']:+.3f} t={r['t']:+.2f} "
                  f"p={r['p']:.3g}{flag}")
        lvl = red.groupby(["animal", "pair", "epoch"])["peak_r"].mean().unstack("epoch")
        print(f"  ⚠ level check — cohort mean peak r: naive {lvl['naive'].mean():.3f} → "
              f"intermediate {lvl['intermediate'].mean():.3f} → expert "
              f"{lvl['expert'].mean():.3f}; expert>intermediate in "
              f"{(lvl['expert'] > lvl['intermediate']).mean():.0%} of animal-pairs")
    (RES / f"cc_crosscorr_epochs{tag}_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / f'cc_crosscorr_epochs{tag}_tables.md'}")


if __name__ == "__main__":
    main()
