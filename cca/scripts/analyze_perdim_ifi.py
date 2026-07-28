"""Per-canonical-dimension IFI — meeting 2026-07-28 item 1: do different CCs have
different IFI?

Pure re-analysis of `results/lag_curves_bin10{,_fsincl}.csv` (held-out, segment-aware
per-dim lag curves, dims 1-30 x lags +/-25 bins). No refits.

Sign convention: POSITIVE IFI => the FIRST-named area leads.

Per (animal, pair) it reduces each canonical dimension's lag curve to one IFI, then asks
three questions per area pair, animals-as-n (the project's inferential unit, STATE.md
S3.0):

  1. Does the dominant dim lead differently from the rest of the spectrum?
     paired t of (IFI_CC1 - mean IFI over the SIGNIFICANT tail), across animals.
  2. Is the direction consistent across canonical rank?
     per-animal fraction of significant tail dims agreeing in sign with CC1, one-sample
     t against the 0.5 chance level.
  3. Does IFI drift systematically with rank?
     random-slope LMM `ifi ~ dim + (1+dim|animal)` over significant dims.

Writes:
    results/perdim_ifi_bin10{,_fsincl}.csv   one row per (animal, pair, dim)
    results/perdim_ifi_tables.md             per-pair stats, both FS conditions

Usage: PYTHONPATH=src python scripts/analyze_perdim_ifi.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import mixed_effects, paired_stats, perdim_ifi  # noqa: E402

RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
FIELDS = ["animal", "learner", "pair", "dim", "ifi", "peak_lag_ms", "cc_max",
          "cc_at_zero", "sig", "n_sig"]


def build_table(df: pd.DataFrame) -> pd.DataFrame:
    """Long lag-curve rows -> one row per (animal, pair, dim)."""
    df = df.copy()
    df["cc"] = pd.to_numeric(df["cc"], errors="coerce")
    rows = []
    for (animal, learner, pair), grp in df.groupby(["animal", "learner", "pair"],
                                                   sort=True):
        dims, ifi, peak = perdim_ifi.ifi_by_dim(grp["lag_ms"].to_numpy(),
                                                grp["cc"].to_numpy(),
                                                grp["dim"].to_numpy())
        for d, i, p in zip(dims, ifi, peak):
            sub = grp[grp["dim"] == d]
            at0 = sub.loc[sub["lag_ms"] == 0, "cc"]
            rows.append({
                "animal": animal, "learner": learner, "pair": pair, "dim": int(d),
                "ifi": i, "peak_lag_ms": p,
                "cc_max": float(np.nanmax(sub["cc"])) if sub["cc"].notna().any()
                          else np.nan,
                "cc_at_zero": float(at0.iloc[0]) if len(at0) else np.nan,
                "sig": int(sub["sig"].max()),
                "n_sig": int(sub["n_sig"].max()),
            })
    return pd.DataFrame(rows, columns=FIELDS)


def per_pair_stats(tab: pd.DataFrame) -> pd.DataFrame:
    """Animals-as-n tests, one row per area pair."""
    out = []
    for pair in PAIRS:
        sub = tab[tab["pair"] == pair].sort_values(["animal", "dim"])
        if sub.empty:
            continue
        deltas, agrees, cc1s, rests, lmm_records = [], [], [], [], []
        for animal, grp in sub.groupby("animal", sort=True):
            grp = grp.sort_values("dim")
            ifi = grp["ifi"].to_numpy(float)
            sig = grp["sig"].to_numpy(int)
            cc1, rest = perdim_ifi.rank_split(ifi, sig)
            # sign agreement uses CC1 as reference against the significant tail only
            ref_tail = np.concatenate([ifi[:1], ifi[1:][sig[1:] == 1]])
            agree = perdim_ifi.sign_agreement(ref_tail)
            if np.isfinite(cc1) and np.isfinite(rest):
                deltas.append(cc1 - rest); cc1s.append(cc1); rests.append(rest)
            if np.isfinite(agree):
                agrees.append(agree)
            for d, v, s in zip(grp["dim"].to_numpy(int), ifi, sig):
                if s == 1 and np.isfinite(v):
                    lmm_records.append({"animal_id": animal, "axis": float(d),
                                        "value": float(v)})
        n, _, t, p = paired_stats.paired_t(deltas) if deltas else (0, np.nan, np.nan,
                                                                  np.nan)
        if agrees:
            n_ag, _, t_ag, p_ag = paired_stats.paired_t(np.asarray(agrees) - 0.5)
        else:
            n_ag, t_ag, p_ag = 0, np.nan, np.nan
        lmm = mixed_effects.lmm_slope(lmm_records) if lmm_records else {}
        out.append({
            "pair": pair,
            "n_animals": n,
            "ifi_cc1": float(np.mean(cc1s)) if cc1s else np.nan,
            "ifi_tail": float(np.mean(rests)) if rests else np.nan,
            "delta": float(np.mean(deltas)) if deltas else np.nan,
            "t": t, "p": p,
            "n_agree": n_ag,
            "sign_agree": float(np.mean(agrees)) if agrees else np.nan,
            "t_agree": t_ag, "p_agree": p_ag,
            "lmm_slope": lmm.get("estimate", np.nan),
            "lmm_p": lmm.get("p", np.nan),
            "lmm_n_obs": lmm.get("n_obs", np.nan),
            "lmm_reason": lmm.get("reason", "not_fit"),
        })
    return pd.DataFrame(out)


def rank_diagnostics(tab: pd.DataFrame) -> dict:
    """Is canonical RANK a meaningful ordering out of sample, and is the significant
    'tail' made of real dims or floor-level ones?

    Motivated by the rank figure: significance turns out to be scattered across all 30
    ranks rather than concentrated at the top. The circular-shift threshold is the
    (1-alpha) quantile of the *dominant*-dim null, which sounds conservative — but in a
    weakly-coupled cell that null is itself tiny, so the bar is tiny and dims sitting at
    the held-out CC floor clear it.
    """
    sig = tab[tab["sig"] == 1]
    by_rank = tab.groupby("dim")["cc_max"].mean()
    ok = tot = 0
    for _, g in tab.groupby(["animal", "pair"]):
        d = np.sort(g.loc[g["sig"] == 1, "dim"].to_numpy())
        if d.size == 0:
            continue
        tot += 1
        ok += int(np.array_equal(d, np.arange(1, d.size + 1)))
    return {
        "n_sig_dims": int(len(sig)),
        "frac_rank_le_5": float((sig["dim"] <= 5).mean()) if len(sig) else np.nan,
        "frac_rank_gt_20": float((sig["dim"] > 20).mean()) if len(sig) else np.nan,
        "cc_rank1": float(by_rank.get(1, np.nan)),
        "cc_rank5": float(by_rank.get(5, np.nan)),
        "cc_rank20": float(by_rank.get(20, np.nan)),
        "spearman_rank_cc": float(tab[["dim", "cc_max"]].corr(
            method="spearman").iloc[0, 1]),
        "frac_contiguous_cells": float(ok / tot) if tot else np.nan,
        "n_cells": tot,
    }


def _fmt_diag(d: dict, fs: str) -> str:
    return (
        f"**{fs}.** {d['n_sig_dims']} significant dims. Held-out CC collapses with rank "
        f"(rank 1 = {d['cc_rank1']:.3f}, rank 5 = {d['cc_rank5']:.3f}, rank 20 = "
        f"{d['cc_rank20']:.3f}; Spearman rank vs CC = {d['spearman_rank_cc']:+.3f}), so "
        f"rank *is* meaningful on average. But only {d['frac_rank_le_5']:.0%} of "
        f"significant dims sit at rank ≤ 5 and {d['frac_rank_gt_20']:.0%} sit beyond "
        f"rank 20, where the mean held-out CC is at the floor; the significant set is a "
        f"contiguous leading run in only {d['frac_contiguous_cells']:.0%} of "
        f"{d['n_cells']} cells. **The significant tail is therefore dominated by "
        f"floor-level dims** — gate on CC magnitude, not the `sig` flag alone, before "
        f"assigning any per-dim label (e.g. FF/FB).\n")


def _fmt(stats: pd.DataFrame, fs: str) -> str:
    lines = [f"### {fs}", "",
             "| pair | n | IFI CC₁ | IFI sig-tail | Δ (CC₁−tail) | t | p | "
             "sign-agree | p vs 0.5 | LMM slope/dim | LMM p |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in stats.iterrows():
        lines.append(
            f"| {r['pair']} | {r['n_animals']:.0f} | {r['ifi_cc1']:+.3f} | "
            f"{r['ifi_tail']:+.3f} | {r['delta']:+.3f} | {r['t']:.2f} | "
            f"{r['p']:.3g} | {r['sign_agree']:.2f} | {r['p_agree']:.3g} | "
            f"{r['lmm_slope']:+.4f} | {r['lmm_p']:.3g} |")
    return "\n".join(lines) + "\n"


def main():
    md = ["# Per-dimension IFI — does directionality depend on canonical rank?",
          "",
          "Meeting 2026-07-28 item 1. Re-analysis of the held-out segment-aware per-dim",
          "lag curves (`lag_curves_bin10*.csv`); no refits. **Positive IFI = the",
          "first-named area leads.** Animals-as-n throughout; 8 pairs, **no cross-pair",
          "correction** (per-pair family, STATE.md §3.0 policy).",
          "",
          "`sign-agree` = per-animal fraction of significant tail dims whose IFI shares",
          "CC₁'s sign; 0.5 = chance. LMM = `ifi ~ dim + (1+dim|animal)` over significant",
          "dims only.", "",
          "> **Caveat.** CCA is refit at every lag, so `dim` is canonical *rank*, not a",
          "> tracked component. These are statements about rank, not about a persistent",
          "> subspace — see meeting item 3.", ""]
    diag = ["", "## Is canonical rank meaningful out of sample?", ""]
    for suf, fs in [("", "FS-excluded"), ("_fsincl", "FS-included")]:
        src = RES / f"lag_curves_bin10{suf}.csv"
        if not src.exists():
            print(f"skip {fs}: {src.name} not found"); continue
        tab = build_table(pd.read_csv(src))
        out_csv = RES / f"perdim_ifi_bin10{suf}.csv"
        tab.to_csv(out_csv, index=False, lineterminator="\n")
        stats = per_pair_stats(tab)
        stats.to_csv(RES / f"perdim_ifi_stats_bin10{suf}.csv", index=False,
                     lineterminator="\n")
        md.append(_fmt(stats, fs))
        diag.append(_fmt_diag(rank_diagnostics(tab), fs))
        n_sig_dims = int((tab["sig"] == 1).sum())
        print(f"{fs}: {len(tab)} (animal,pair,dim) rows, {n_sig_dims} significant "
              f"dims, {tab['animal'].nunique()} animals -> {out_csv.name}")
        hits = stats[stats["p"] < 0.05]
        if len(hits):
            for _, r in hits.iterrows():
                print(f"   CC1-vs-tail: {r['pair']} Δ={r['delta']:+.3f} "
                      f"t={r['t']:.2f} p={r['p']:.4g} (n={r['n_animals']:.0f})")
        else:
            print("   CC1-vs-tail: no pair at p<0.05")
    (RES / "perdim_ifi_tables.md").write_text("\n".join(md + diag))
    print(f"wrote {RES / 'perdim_ifi_tables.md'}")


if __name__ == "__main__":
    main()
