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

⚠ LEVEL contrasts against the naive epoch are NOT a learning readout on this fit.
Found while building this (2026-08-15) and adversarially verified: mean peak r over
animal-pairs is naive 0.052 → intermediate 0.067 → expert 0.067 (FS-incl 0.064 →
0.078 → 0.079); intermediate > naive in 88 % of animal-pairs, expert > intermediate
in 56 % / 47 % (chance). What the evidence supports:
  * the naive deficit is ~60 % a lag-INDEPENDENT offset — the whole curve, baseline
    at |lag| ≥ 200 ms included, sits lower in naive for the cortical pairs (CA1-RSC,
    CA1-V1, V1-RSC, RSC-SUB); only CA3-DG (n = 4) shows a peak-specific change. A
    flat offset on frozen axes is slow co-modulation, and running speed (+6.6 cm/s
    naive→expert, HANDOFF §6.1) is the obvious candidate. `curve_metrics` therefore
    reports peak r MINUS the far-lag baseline as the coupling-specific statistic;
  * the naive→intermediate rise is independent of the learning point (per-animal
    Δ(int − naive) vs LP: Spearman ρ ≈ −0.01, verifier check on
    `trajectory_w15_bin10.csv`), and the HELD-OUT per-epoch REFIT arms show no
    naive deficit at all — `epoch_metrics_bin10` CC₁ 0.149 / 0.139 / 0.145
    (int − naive p = 0.84, exp − naive p = 0.91), `lag_subspaces_bin10_epochs`
    lag-0 CC₁ 0.149 / 0.170 / 0.147 (p = 0.97), FS-incl likewise. That is the
    decisive evidence: coupling strength does not change; what changes is how
    well ONE whole-session set of axes fits the session's opening;
  * intermediate ≈ expert on its own is NOT decisive — they are adjacent 10-trial
    blocks, so gradual or pre-LP change would look the same; and the balanced-trial
    frozen fit (`fixed_subspace_stats`, Δpeak r mixed-sign, all n.s.) is CC₁-only on
    600 bins/trial and a 30-trial fit set, so it is corroboration, not a like-for-
    like control.
So: IFI (a shape ratio; insensitive to proportional scaling, though not to an
additive baseline shift) is the naive-vs-expert statistic and is null; read curve
height as "the first ~10 trials are different, mostly as an offset", not as
learning. The stats file carries intermediate − naive next to expert − naive, and
peak_minus_far / far_r next to peak_r, as the built-in diagnostics.

Writes:
    results/cc_crosscorr_epochs_bin10{,_fsincl}.csv        per (animal, pair, epoch, group, lag)
    results/cc_crosscorr_epochs_stats_bin10{,_fsincl}.csv  per pair: expert−naive, group=all
    results/cc_crosscorr_epochs_tables.md

Usage: PYTHONPATH=src python scripts/analyze_cc_crosscorr_epochs.py [--label-col label|label_loo|label_int|label_xv]
    label      whole-session label (default; NOT cross-validated)
    label_int  intermediate-epoch IFI sign — disjoint from naive/expert (Tom)
    label_xv   IFI sign on all trials outside naive ∪ expert (needs the 2026-08-18 driver)
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


def curve_metrics(curves: pd.DataFrame, group: str = "all", far_ms: int = 200) -> pd.DataFrame:
    """Per (animal, pair, epoch): peak r, the far-lag baseline (mean r at
    |lag| >= ``far_ms``) and their difference — from the per-animal mean curves.

    ``peak_minus_far`` is the coupling-specific statistic: a lag-INDEPENDENT offset
    of the whole cross-correlogram (slow co-modulation on the frozen axes — running
    speed is the obvious candidate) moves ``peak_r`` and ``far_r`` together and leaves
    this unchanged; a change in lagged coupling moves it. Verifiers found (2026-08-15)
    that ~60 % of the naive-vs-later peak-r difference is such an offset in the
    cortical pairs, so peak r alone overstates what changed.
    """
    sub = curves[curves["group"] == group]
    rows = []
    for (an, pair, ep), g in sub.groupby(["animal", "pair", "epoch"]):
        r = pd.to_numeric(g["r_mean"], errors="coerce").to_numpy(float)
        lag = g["lag_ms"].to_numpy(float)
        ok = np.isfinite(r)
        if not ok.any():
            continue
        far = ok & (np.abs(lag) >= far_ms)
        far_r = float(np.mean(r[far])) if far.any() else np.nan
        peak = float(np.max(r[ok]))
        rows.append({"animal": an, "pair": pair, "epoch": ep, "peak_r": peak,
                     "far_r": far_r, "peak_minus_far": peak - far_r})
    return pd.DataFrame(rows)


def epoch_contrast(reduced: pd.DataFrame, metric: str, pairs=PAIRS,
                   a: str = "naive", b: str = "expert", min_n: int = 3) -> pd.DataFrame:
    """Per pair: paired *t* of (``b`` − ``a``) on the per-animal mean of ``metric``
    over all significant CCs. ``reduced`` is one row per (animal, pair, dim, epoch)
    — `cc_label_track_epoch_bin10*.csv` from analyze_cc_label_track.py — so the
    per-CC reductions (IFI at ±50 ms, peak r) are the ones already on disk.

    Animals lacking either epoch are dropped. Eight pairs = eight families
    (project policy); ``bh_pass`` across pairs is a sensitivity column.
    """
    # `cc_label_track_epoch_*` was written sig-gated by analyze_cc_label_track
    # (reduce_curves, sig_only=True) and carries no `sig` column — say so explicitly.
    pa = cc_aggregate.per_animal_mean(reduced, value=metric, by=["epoch"],
                                      sig_only="sig" in reduced.columns,
                                      drop_degenerate=False)
    rows = []
    for pair in pairs:
        w = pa[pa["pair"] == pair].pivot_table(index="animal", columns="epoch",
                                               values="mean")
        if a not in w.columns or b not in w.columns:
            continue
        both = np.isfinite(w[a].to_numpy(float)) & np.isfinite(w[b].to_numpy(float))
        if both.sum() < min_n:
            continue
        va, vb = w[a].to_numpy(float)[both], w[b].to_numpy(float)[both]
        n, _, t, p = paired_stats.paired_t(vb - va)
        rows.append({"pair": pair, "metric": metric, "contrast": f"{b}-{a}",
                     "unit": "animals", "n": int(n),
                     "mean_a": float(va.mean()), "mean_b": float(vb.mean()),
                     "delta": float((vb - va).mean()), "t": t, "p": p})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["bh_pass"] = paired_stats.fdr_bh(out["p"].to_numpy(float))
    return out


def epoch_contrast_ccs(reduced: pd.DataFrame, metric: str, pairs=PAIRS,
                       a: str = "naive", b: str = "expert", min_n: int = 3) -> pd.DataFrame:
    """The same contrast with **CCs as n**: every significant (animal, dim) is one
    paired sample, pooled across animals within a pair. Frozen axes make the pairing
    exact — dim k of animal A in the naive epoch IS dim k of animal A in expert.

    ⚠ This is the field convention (Gonzalez & Buzsáki) and the more powerful unit,
    but CCs are nested in animals — an animal with ten CCs weighs ten times one with
    one, and the ten are not independent — so the project reports it as a POWER
    CHECK next to animals-as-n, never instead of it (STATE.md §3.0 units policy).
    """
    if "sig" in reduced.columns:
        reduced = reduced[pd.to_numeric(reduced["sig"], errors="coerce") == 1]
    rows = []
    for pair in pairs:
        w = reduced[reduced["pair"] == pair].pivot_table(index=["animal", "dim"],
                                                         columns="epoch", values=metric)
        if a not in w.columns or b not in w.columns:
            continue
        both = np.isfinite(w[a].to_numpy(float)) & np.isfinite(w[b].to_numpy(float))
        if both.sum() < min_n:
            continue
        va, vb = w[a].to_numpy(float)[both], w[b].to_numpy(float)[both]
        n, _, t, p = paired_stats.paired_t(vb - va)
        rows.append({"pair": pair, "metric": metric, "contrast": f"{b}-{a}",
                     "unit": "ccs", "n": int(n),
                     "mean_a": float(va.mean()), "mean_b": float(vb.mean()),
                     "delta": float((vb - va).mean()), "t": t, "p": p})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["bh_pass"] = paired_stats.fdr_bh(out["p"].to_numpy(float))
    return out


def per_cc_curve_metrics(df: pd.DataFrame, far_ms: int = 200,
                         sig_only: bool = True) -> pd.DataFrame:
    """Per (animal, pair, dim, epoch) from the RAW label-track curves: peak r, the
    far-lag baseline (|lag| >= ``far_ms``) and peak − baseline — the per-CC version
    of :func:`curve_metrics`, so the CCs-as-n contrast can be run on the same three
    quantities as the animals-as-n one.
    """
    if sig_only and "sig" in df.columns:
        df = df[pd.to_numeric(df["sig"], errors="coerce") == 1]
    rows = []
    for (an, pair, dim, ep), g in df.groupby(["animal", "pair", "dim", "epoch"]):
        r = pd.to_numeric(g["r"], errors="coerce").to_numpy(float)
        lag = pd.to_numeric(g["lag_ms"], errors="coerce").to_numpy(float)
        ok = np.isfinite(r) & np.isfinite(lag)
        if not ok.any():
            continue
        far = ok & (np.abs(lag) >= far_ms)
        far_r = float(np.mean(r[far])) if far.any() else np.nan
        peak = float(np.max(r[ok]))
        rows.append({"animal": an, "pair": pair, "dim": int(dim), "epoch": ep,
                     "peak_r": peak, "far_r": far_r, "peak_minus_far": peak - far_r})
    return pd.DataFrame(rows)


LABEL_TAGS = {"label": "", "label_loo": "_loo", "label_int": "_labint", "label_xv": "_labxv"}
LABEL_DESC = {
    "label": "whole-session fit at ±50 ms (as on slide 10) — NOT cross-validated: the "
             "plotted epochs contributed to the label",
    "label_loo": "leave-the-scored-epoch-out (each epoch's curve labelled without its own "
                 "trials; the FF set can differ between the naive and expert panels)",
    "label_int": "sign of the INTERMEDIATE epoch's own IFI at ±50 ms — 10 trials disjoint "
                 "from both plotted epochs (Tom's suggestion; noisy, few trials)",
    "label_xv": "sign of the IFI at ±50 ms on ALL trials outside naive ∪ expert (~100+ "
                "trials, disjoint from both plotted epochs) — the cross-validated label",
}


def attach_label(df: pd.DataFrame, red: pd.DataFrame, label_col: str) -> pd.DataFrame:
    """Make ``label_col`` available on the raw curve table.

    ``label`` / ``label_loo`` / ``label_xv`` are columns the driver writes;
    ``label_int`` is derived here from the intermediate epoch's per-CC IFI in the
    reduced epoch file (`cc_label_track_epoch_*`): FF if > 0, FB if < 0.
    """
    if label_col == "label_int":
        mid = red[red["epoch"] == "intermediate"][["animal", "pair", "dim", "ifi"]].copy()
        mid["label_int"] = np.where(mid["ifi"] > 0, "FF",
                                    np.where(mid["ifi"] < 0, "FB", "none"))
        return df.merge(mid[["animal", "pair", "dim", "label_int"]],
                        on=["animal", "pair", "dim"], how="left")
    if label_col not in df.columns:
        raise KeyError(f"{label_col} not in the label-track file — re-run "
                       f"run_cc_label_track.py (label_xv was added 2026-08-18)")
    return df


def label_contrast(red: pd.DataFrame, labels: pd.DataFrame, label_col: str,
                   metric: str = "ifi", a: str = "naive", b: str = "expert") -> pd.DataFrame:
    """Per pair and per label group (FF / FB): the ``b − a`` paired contrast on the
    per-animal mean of ``metric`` over that animal's CCs carrying that label — the
    per-label version of :func:`epoch_contrast`, for whichever label column is in
    use, so the figure titles never quote a test made with a different label.
    ``labels`` = one row per (animal, pair, dim) with ``label_col``."""
    lab = labels[["animal", "pair", "dim", label_col]].drop_duplicates()
    r = red.drop(columns=[c for c in (label_col,) if c in red.columns]).merge(
        lab, on=["animal", "pair", "dim"], how="left")
    parts = []
    for grp in ("FF", "FB"):
        sub = r[r[label_col] == grp]
        if sub.empty:
            continue
        s = epoch_contrast(sub, metric, a=a, b=b)
        if not s.empty:
            s.insert(1, "group", grp)
            parts.append(s)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _md(stats: pd.DataFrame, curves: pd.DataFrame, fs: str, label_col: str) -> str:
    lines = [f"### {fs} — cross-correlograms by epoch, all significant CCs "
             f"(frozen axes, per-animal-first)", "",
             f"Label column for the FF/FB split: `{label_col}` — {LABEL_DESC[label_col]}. "
             f"Paired *t* across "
             f"animals on each animal's mean over its significant CCs of the per-CC "
             f"reduction already in `cc_label_track_epoch_*` (IFI at ±50 ms; peak r). "
             f"Positive IFI ⇒ first-named area leads. Per-pair families; BH across the "
             f"8 pairs is a sensitivity column.", "",
             "> ⚠ **Curve-height (peak r) contrasts against naive are not a learning "
             "readout on this all-trials fit** — see the script docstring: naive is "
             "uniquely low, ~60 % of that is a lag-independent OFFSET (`far_r`, the "
             "|lag| ≥ 200 ms baseline, moves with it — slow co-modulation, speed the "
             "obvious candidate), the rise is LP-independent, and intermediate (pre-LP) "
             "already equals expert. `peak_minus_far` is the coupling-specific "
             "statistic; `intermediate-naive` sits next to `expert-naive`. IFI (shape) "
             "is the naive-vs-expert statistic and is null.", "",
             "> Rows with n ≤ 4 animals (3 df) are descriptive and are not bolded.", "",
             "> **Two units per row.** *animals-as-n* (inferential unit): each animal's "
             "significant CCs are averaged first, then a paired *t* across animals. "
             "*CCs-as-n* (power check, field convention): every significant (animal, "
             "dim) is one paired sample, pooled across animals — CCs are nested in "
             "animals, so this over-counts and is never the inferential statement.", "",
             "| pair | metric | contrast | a → b (animals) | Δ | animals: n, t, p, BH | "
             "CCs: n, Δ, t, p, BH |",
             "|---|---|---|---|---|---|---|"]
    g_all = stats[stats["group"] == "all"] if "group" in stats.columns else stats
    an = g_all[g_all["unit"] == "animals"].set_index(["pair", "metric", "contrast"])
    cc = g_all[g_all["unit"] == "ccs"].set_index(["pair", "metric", "contrast"])
    for key, r in an.iterrows():
        star = "**" if (r["p"] < 0.05 and r["n"] > 4) else ""
        c = cc.loc[key] if key in cc.index else None
        cs = (f"{int(c['n'])}, {c['delta']:+.3f}, {c['t']:+.2f}, "
              f"{'**' if c['p'] < 0.05 else ''}{c['p']:.3g}{'**' if c['p'] < 0.05 else ''}, "
              f"{'yes' if bool(c['bh_pass']) else 'no'}" if c is not None else "—")
        lines.append(f"| {key[0]} | {key[1]} | {key[2]} | "
                     f"{r['mean_a']:+.3f} → {r['mean_b']:+.3f} | {star}{r['delta']:+.3f}{star} | "
                     f"{int(r['n'])}, {r['t']:+.2f}, {star}{r['p']:.3g}{star}, "
                     f"{'yes' if bool(r['bh_pass']) else 'no'} | {cs} |")
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
    if "group" in stats.columns and (stats["group"] != "all").any():
        lines += ["", f"**Per-label expert − naive (label = `{label_col}`), animals-as-n:**",
                  "", "| pair | label | metric | n | naive → expert | Δ | t | p |",
                  "|---|---|---|---|---|---|---|---|"]
        for _, r in stats[stats["group"] != "all"].iterrows():
            star = "**" if (r["p"] < 0.05 and r["n"] > 4) else ""
            lines.append(f"| {r['pair']} | {r['group']} | {r['metric']} | {int(r['n'])} | "
                         f"{r['mean_a']:+.3f} → {r['mean_b']:+.3f} | {star}{r['delta']:+.3f}{star} | "
                         f"{r['t']:+.2f} | {star}{r['p']:.3g}{star} |")
    return "\n".join(lines) + "\n"


def _label_col_from_argv(argv) -> str:
    if "--loo" in argv:
        return "label_loo"
    if "--label-col" in argv:
        v = argv[argv.index("--label-col") + 1]
        if v not in LABEL_TAGS:
            sys.exit(f"--label-col must be one of {sorted(LABEL_TAGS)}")
        return v
    return "label"


def main():
    label_col = _label_col_from_argv(sys.argv[1:])
    tag = LABEL_TAGS[label_col]
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
        red = pd.read_csv(red_src)
        df = attach_label(df, red, label_col)
        curves = build_curves(df, label_col=label_col)
        curves.to_csv(RES / f"cc_crosscorr_epochs{tag}_bin10{suf}.csv", index=False,
                      lineterminator="\n")
        # the epoch file has no `sig` column: it was written sig-gated already
        cm = curve_metrics(curves)                 # per animal (mean curve)
        pcm = per_cc_curve_metrics(df)             # per CC (raw curves)
        contrasts = [("ifi", "naive", "expert"), ("peak_r", "naive", "expert"),
                     ("peak_r", "naive", "intermediate"),
                     ("far_r", "naive", "expert"), ("peak_minus_far", "naive", "expert"),
                     ("peak_minus_far", "naive", "intermediate")]
        parts = []
        for metric, ea, eb in contrasts:
            src_an = red if metric in ("ifi", "peak_r") else cm
            src_cc = red if metric in ("ifi", "peak_r") else pcm
            parts.append(epoch_contrast(src_an, metric, a=ea, b=eb))
            parts.append(epoch_contrast_ccs(src_cc, metric, a=ea, b=eb))
        stats = pd.concat(parts, ignore_index=True)
        stats.insert(1, "group", "all")
        labels = df[["animal", "pair", "dim", label_col]].drop_duplicates()
        stats = pd.concat([stats, label_contrast(red, labels, label_col, "ifi"),
                           label_contrast(red, labels, label_col, "peak_r")],
                          ignore_index=True)
        stats.to_csv(RES / f"cc_crosscorr_epochs_stats{tag}_bin10{suf}.csv",
                     index=False, lineterminator="\n")
        md.append(_md(stats, curves, fs, label_col))
        n_cells = df.groupby(["animal", "pair", "dim"]).ngroups
        n_sig = df[df["sig"] == 1].groupby(["animal", "pair", "dim"]).ngroups
        print(f"\n[{fs}] {n_sig}/{n_cells} (animal, pair, dim) cells significant; "
              f"{len(curves)} curve rows")
        print("  contrasts (paired t) — animals-as-n | CCs-as-n:")
        g_all = stats[stats["group"] == "all"]
        an = g_all[g_all["unit"] == "animals"].set_index(["pair", "metric", "contrast"])
        cc = g_all[g_all["unit"] == "ccs"].set_index(["pair", "metric", "contrast"])
        for key, r in an.iterrows():
            flag = " *" if r["p"] < 0.05 else ""
            c = cc.loc[key] if key in cc.index else None
            cs = (f"| CCs n={int(c['n']):<3d} Δ={c['delta']:+.3f} p={c['p']:.3g}"
                  f"{' *' if c['p'] < 0.05 else ''}" if c is not None else "")
            print(f"    {key[0]:8s} {key[1]:14s} {key[2]:19s} n={int(r['n']):<2d} "
                  f"Δ={r['delta']:+.3f} p={r['p']:.3g}{flag:3s} {cs}")
        lvl = red.groupby(["animal", "pair", "epoch"])["peak_r"].mean().unstack("epoch")
        print(f"  ⚠ level check (mean over animal-pairs) — peak r: naive "
              f"{lvl['naive'].mean():.3f} → intermediate {lvl['intermediate'].mean():.3f} "
              f"→ expert {lvl['expert'].mean():.3f}; expert>intermediate in "
              f"{(lvl['expert'] > lvl['intermediate']).mean():.0%} of animal-pairs")
        off = cm.pivot_table(index=["animal", "pair"], columns="epoch",
                             values=["far_r", "peak_minus_far"])
        d_far = (off["far_r"]["intermediate"] - off["far_r"]["naive"]).mean()
        d_pmf = (off["peak_minus_far"]["intermediate"] - off["peak_minus_far"]["naive"]).mean()
        print(f"  ⚠ offset check — int − naive: far-lag baseline Δ{d_far:+.4f}, "
              f"peak-minus-baseline Δ{d_pmf:+.4f} (mean over animal-pairs)")
    (RES / f"cc_crosscorr_epochs{tag}_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / f'cc_crosscorr_epochs{tag}_tables.md'}")


if __name__ == "__main__":
    main()
