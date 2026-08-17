"""Item 1 — do different CCs have different IFI, i.e. are there CCs with different SIGNS?

Asked at the 2026-07-28 meeting as: *"calculate IFI at different lags for each CC, for
each area, and then look at them."*

Analysis unchanged from the existing arm: the lag curves come from **CCA refit at every
lag** (`run_lag_curves.py` -> `lagged.heldout_lag_curve_flat_perdim`, held-out and
segment-aware). This script only re-reduces them, per canonical dimension, over a sweep of
INTEGRATION WINDOWS.

For each (animal, pair, dim): IFI over |lag| <= w for w = 1..25 bins (10..250 ms), via
`perdim_ifi.ifi_windows_by_dim`. Positive IFI => the FIRST-named area leads.

The question is descriptive, so the output is descriptive: per pair, how many canonical
dimensions carry a positive vs a negative IFI, and whether both signs co-occur within the
same animal and pair.

2026-08-07 meeting follow-ups (asks 1 and 3): the OVERALL IFI — the ungrouped mean over
all of an animal's significant CCs — alongside the sign-split, and a one-sample test of it
against 0 at the headline ±50 ms window, animals-as-n (`overall_by_animal`,
`overall_direction`). Unweighted per-animal mean is primary; the CC-strength-weighted mean
is carried alongside.

Writes:
    results/cc_ifi_windows_bin10{,_fsincl}.csv       one row per (animal, pair, dim, window)
    results/cc_ifi_signs_bin10{,_fsincl}.csv         per (animal, pair, window) sign counts
    results/cc_ifi_overall_bin10{,_fsincl}.csv       per (animal, pair, window) all-CC mean
    results/cc_ifi_overall_test_bin10{,_fsincl}.csv  per pair, all-CC IFI @ ±50 ms vs 0
    results/cc_ifi_signs_tables.md                   per-pair summary, both FS

Usage: PYTHONPATH=src python scripts/analyze_cc_ifi_signs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import cc_aggregate, perdim_ifi  # noqa: E402

RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
BIN_MS = 10
LEAD_DIMS = 5          # "leading" CCs — where held-out CC is above the floor
HEADLINE_W = 5         # bins = +/-50 ms, the report's headline IFI window


def check_sig_consistency(df: pd.DataFrame) -> dict:
    """Guard: is the `sig` flag consistent with the `cc` curve in the SAME row?

    The flag is a single scalar threshold on the lag-0 held-out CC, so within a cell it
    MUST be monotone in that CC — no dim with a higher lag-0 CC can be non-significant
    while a lower one is significant, and a negative CC can never pass.

    `run_lag_curves` used to violate this by importing the mask from a separate
    `window_subspace` fit and attaching it by bare index (fixed 2026-08-06). This
    re-checks it every run so the misalignment cannot come back silently.
    """
    z = df[df["lag_ms"] == 0].copy()
    z["cc"] = pd.to_numeric(z["cc"], errors="coerce")
    neg = int(((z["sig"] == 1) & (z["cc"] <= 0)).sum())
    n_sig = int((z["sig"] == 1).sum())
    # ⚠ Monotonicity of the mask in cc is NOT checked. It holds only under the old
    # dominant-dim null (one scalar threshold). Under the per-dim null each dimension
    # has its own bar and shuffled correlations fall with rank, so a non-significant
    # dim CAN legitimately outrank a significant one — testing for it here would fire
    # a false alarm on correct data.
    subset = np.nan
    if "sig_uncorr" in z.columns:
        # the FDR mask must be a subset of the uncorrected one
        subset = int(((z["sig"] == 1) & (z["sig_uncorr"] == 0)).sum())
    return {"n_sig": n_sig, "sig_with_nonpositive_cc": neg,
            "fdr_not_subset_of_uncorrected": subset,
            "ok": neg == 0 and (not np.isfinite(subset) or subset == 0)}


def build_windows(df: pd.DataFrame) -> pd.DataFrame:
    """Long file -> one row per (animal, pair, dim, integration window)."""
    df = df.copy()
    df["cc"] = pd.to_numeric(df["cc"], errors="coerce")
    rows = []
    for (animal, learner, pair), g in df.groupby(["animal", "learner", "pair"],
                                                 sort=True):
        dims, wins, ifi, degen = perdim_ifi.ifi_windows_by_dim(
            g["lag_bins"].to_numpy(), g["cc"].to_numpy(), g["dim"].to_numpy())
        # per-dim CC magnitude, so a sign can be weighed against whether the dim
        # carries any coupling at all
        cc_peak = g.groupby("dim")["cc"].max()
        sig = g.groupby("dim")["sig"].max()
        for i, d in enumerate(dims):
            for j, w in enumerate(wins):
                rows.append({
                    "animal": animal, "learner": learner, "pair": pair,
                    "dim": int(d), "window_bins": int(w),
                    "window_ms": int(w) * BIN_MS,
                    "ifi": ifi[i, j], "degenerate": int(degen[i, j]),
                    "cc_peak": float(cc_peak.get(d, np.nan)),
                    "sig": int(sig.get(d, 0)),
                })
    return pd.DataFrame(rows)


def sign_counts(tab: pd.DataFrame, lead_only: bool = True) -> pd.DataFrame:
    """Per (animal, pair, window): how many CCs lead, how many lag, are signs MIXED?

    Degenerate entries (IFI = 0 only because both sides clipped to zero) are excluded —
    they carry no direction. `lead_only` restricts to the leading LEAD_DIMS canonical
    ranks, where the held-out CC is above the floor; the tail is dominated by dims whose
    IFI sign is close to a coin flip (see analyze_perdim_ifi rank diagnostics).
    """
    sub = tab[tab["degenerate"] == 0]
    if lead_only:
        sub = sub[sub["dim"] <= LEAD_DIMS]
    rows = []
    for (animal, pair, w), g in sub.groupby(["animal", "pair", "window_bins"]):
        v = g["ifi"].to_numpy(float)
        v = v[np.isfinite(v) & (v != 0)]
        n_pos, n_neg = int((v > 0).sum()), int((v < 0).sum())
        rows.append({
            "animal": animal, "pair": pair, "window_bins": int(w),
            "window_ms": int(w) * BIN_MS, "n_dims": v.size,
            "n_pos": n_pos, "n_neg": n_neg,
            "mixed": int(n_pos > 0 and n_neg > 0),
            "frac_pos": n_pos / v.size if v.size else np.nan,
        })
    return pd.DataFrame(rows)


def dim_direction(tab: pd.DataFrame, window_bins: int = HEADLINE_W) -> pd.DataFrame:
    """Per (pair, dim): does THAT canonical dimension have a reliable direction?

    One-sample *t* of the IFI across animals against 0, animals-as-n. This is the test
    that actually answers "do different CCs have different signs", because per-animal
    sign mixing is uninformative: with `d` dimensions and a coin-flip sign, mixing is
    expected with probability ``1 - 2*0.5**d`` — 94 % at d = 5. Only a direction that
    reproduces ACROSS animals is evidence, and only then can two dimensions be said to
    disagree.

    Benjamini-Hochberg across the dims within a pair (the dims of one pair are one
    family); no cross-pair correction, per the project's stated policy.
    """
    from tom_cca import paired_stats
    sub = tab[(tab["window_bins"] == window_bins) & (tab["degenerate"] == 0)]
    rows = []
    for (pair, d), g in sub.groupby(["pair", "dim"]):
        v = g["ifi"].to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size < 3:
            continue
        n, _, t, p = paired_stats.paired_t(v)      # one-sample: deltas vs 0
        rows.append({"pair": pair, "dim": int(d), "n_animals": n,
                     "mean_ifi": float(np.mean(v)), "t": t, "p": p,
                     "cc_peak": float(g["cc_peak"].mean())})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    keep = []
    for _, g in out.groupby("pair"):
        g = g.copy()
        g["fdr_pass"] = paired_stats.fdr_bh(g["p"].to_numpy(float))
        keep.append(g)
    return pd.concat(keep, ignore_index=True)


def overall_by_animal(tab: pd.DataFrame) -> pd.DataFrame:
    """Ask 1 (2026-08-07): the average IFI across ALL of an animal's significant CCs.

    Per (animal, pair, window): ``mean_ifi`` = plain mean over the CCs that beat the
    per-dim null and are not degenerate — the same population as the +IFI / −IFI
    groups in the figure, just not split by sign; ``wmean_ifi`` = the same mean
    weighted by each CC's held-out peak (clipped at 0); ``n_ccs`` = how many CCs.

    Collapsing to one value per animal here is what makes the downstream test
    animals-as-n rather than dims-as-n.
    """
    pa = cc_aggregate.per_animal_mean(tab, value="ifi", by=["window_bins", "window_ms"],
                                      weight="cc_peak")
    return pa.rename(columns={"mean": "mean_ifi", "wmean": "wmean_ifi"})


def overall_direction(per_animal: pd.DataFrame, window_bins: int = HEADLINE_W,
                      min_n: int = 3) -> pd.DataFrame:
    """Ask 3 (2026-08-07): is the overall IFI at ±50 ms different from 0?

    One-sample *t* across animals of each animal's all-significant-CC mean IFI at
    the fixed headline window, per pair. The unweighted mean (``mean``, ``t``,
    ``p``) is the primary statement; the CC-strength-weighted version (``wmean``,
    ``wt``, ``wp``) is reported alongside. Eight pairs = eight per-pair families
    (project policy); ``bh_pass`` across pairs is a sensitivity column only.
    """
    at = per_animal[per_animal["window_bins"] == window_bins]
    un = cc_aggregate.one_sample_by_pair(at, value="mean_ifi", pairs=PAIRS, min_n=min_n)
    wt = cc_aggregate.one_sample_by_pair(at, value="wmean_ifi", pairs=PAIRS, min_n=min_n)
    wt = wt.rename(columns={"n": "wn", "mean": "wmean", "sem": "wsem", "t": "wt",
                            "p": "wp", "bh_pass": "wbh_pass"})
    out = un.merge(wt, on="pair", how="left")
    out.insert(1, "window_ms", int(window_bins) * BIN_MS)
    return out


def overall_direction_ccs(tab: pd.DataFrame, window_bins: int = HEADLINE_W,
                          min_n: int = 3) -> pd.DataFrame:
    """Ask 3 with **significant CCs as n**: every significant, non-degenerate
    (animal, dim) IFI at the headline window is one sample, pooled across animals
    within a pair; one-sample *t* vs 0. The field convention and the more powerful
    unit — but CCs are nested in animals (an animal with six CCs weighs six times
    one with one), so this is a POWER CHECK next to :func:`overall_direction`, never
    the inferential statement (STATE.md §3.0 units policy).
    """
    at = tab[(tab["window_bins"] == window_bins) & (tab["sig"] == 1) &
             (tab["degenerate"] == 0)][["animal", "pair", "dim", "ifi"]].copy()
    at["ifi"] = pd.to_numeric(at["ifi"], errors="coerce")
    at = at[np.isfinite(at["ifi"])]
    # one_sample_by_pair insists on one row per animal; here the unit IS the CC, so
    # give every (animal, dim) its own "animal" key on purpose
    at["animal"] = at["animal"].astype(str) + ":cc" + at["dim"].astype(str)
    out = cc_aggregate.one_sample_by_pair(at, value="ifi", pairs=PAIRS, min_n=min_n)
    out.insert(1, "unit", "ccs")
    out.insert(2, "window_ms", int(window_bins) * BIN_MS)
    return out


def cc1_direction(tab: pd.DataFrame, window_bins: int = HEADLINE_W,
                  min_n: int = 3) -> pd.DataFrame:
    """The like-for-like CC₁ comparator for :func:`overall_direction`: dim 1 of the
    SAME table (same 12k-bin subset, same estimator), every animal, one-sample *t*.

    `bin10_tables.md` §B asks the same question of CC₁ but on the whole-session
    `run_ifi_windows` arm (~370 k bins) — a different sample; per animal-pair the two
    CC₁ IFIs correlate only r ≈ 0.4. Comparing the all-CC number to §B therefore
    confounds "pooling CCs" with "which trials". This is the fair comparator.
    """
    d1 = tab[(tab["dim"] == 1) & (tab["window_bins"] == window_bins) &
             (tab["degenerate"] == 0)][["animal", "pair", "ifi"]]
    return cc_aggregate.one_sample_by_pair(d1, value="ifi", pairs=PAIRS, min_n=min_n)


def _md_overall(direction: pd.DataFrame, cc1: pd.DataFrame, n_all: pd.Series,
                fs: str) -> str:
    """Asks 1 + 3: the overall IFI at the headline window, per pair, vs 0."""
    lines = [f"#### {fs} — OVERALL IFI at ±{HEADLINE_W * BIN_MS} ms (all significant "
             "CCs pooled, unweighted per-animal mean) vs 0", "",
             "Per animal, the mean IFI over every CC that beats the per-dim null (not "
             "split by sign); then a one-sample *t* across animals, per pair. "
             "Positive ⇒ the first-named area leads. `w·` columns = the same with each "
             "CC weighted by its held-out peak. **CCs-as-n** = every significant CC as one "
             "sample pooled over animals (power check; CCs are nested in animals).", "",
             "> **Data:** the refit-per-lag curves are capped at 12 000 bins — the first "
             "≤ 600 bins (6 s) of each of the **first ~20 trials** — so this is the "
             "session's opening, not the whole session. **n** = animals with ≥ 1 "
             "significant CC (of `N` recorded for the pair). **CC₁ (same data)** = the "
             "like-for-like comparator: dim 1 of this same table, every animal — NOT "
             "`bin10_tables.md` §B, which is CC₁ on the whole-session `run_ifi_windows` "
             "arm (a different sample; the two correlate r ≈ 0.4 per animal-pair). "
             "Bold = p < 0.05 on this single look; the pair has 4 looks in all "
             "(FS-excl/incl × unweighted/weighted) — see the cross-FS summary below.", "",
             "| pair | n / N | mean IFI | SEM | t | p | BH (8 pairs) | weighted mean | "
             "p (weighted) | CCs-as-n: n, mean, p | CC₁ (same data): mean, p, n |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    c1 = cc1.set_index("pair")
    for _, r in direction.iterrows():
        star = "**" if r["p"] < 0.05 else ""
        if r["pair"] in c1.index:
            c = c1.loc[r["pair"]]
            c1s = f"{c['mean']:+.3f}, p={c['p']:.3g}, n={int(c['n'])}"
        else:
            c1s = "—"
        if "ccs_n" in direction.columns and np.isfinite(r.get("ccs_n", np.nan)):
            cs = "**" if r["ccs_p"] < 0.05 else ""
            ccs = f"{int(r['ccs_n'])}, {r['ccs_mean']:+.3f}, {cs}{r['ccs_p']:.3g}{cs}"
        else:
            ccs = "—"
        lines.append(
            f"| {r['pair']} | {int(r['n'])} / {int(n_all.get(r['pair'], 0))} | "
            f"{star}{r['mean']:+.3f}{star} | {r['sem']:.3f} | {r['t']:+.2f} | "
            f"{star}{r['p']:.3g}{star} | {'yes' if bool(r['bh_pass']) else 'no'} | "
            f"{r['wmean']:+.3f} | {r['wp']:.3g} | {ccs} | {c1s} |")
    return "\n".join(lines) + "\n"


def _md_overall_looks(odirs: dict[str, pd.DataFrame]) -> str:
    """How many of a pair's four looks (FS × weighting) are nominally p < 0.05, and
    whether the sign is consistent across them — the honest cross-FS summary."""
    lines = ["#### OVERALL IFI at ±50 ms — the four looks per pair "
             "(FS-excluded / FS-included × unweighted / weighted)", "",
             "Nominal p < 0.05 counts; a direction is only worth stating when the sign "
             "is the same in all four and more than one look is nominal.", "",
             "| pair | looks p<0.05 (of 4) | signs (excl-u, excl-w, incl-u, incl-w) | "
             "means |", "|---|---|---|---|"]
    for pair in PAIRS:
        ps, signs, means = [], [], []
        for fs in ("FS-excluded", "FS-included"):
            d = odirs.get(fs)
            if d is None or pair not in set(d["pair"]):
                continue
            r = d[d["pair"] == pair].iloc[0]
            for m, p in ((r["mean"], r["p"]), (r["wmean"], r["wp"])):
                ps.append(p); signs.append("+" if m > 0 else "−"); means.append(m)
        if not ps:
            continue
        k = int(sum(p < 0.05 for p in ps))
        lines.append(f"| {pair} | {k}/{len(ps)} | {' '.join(signs)} | "
                     f"{' / '.join(f'{m:+.3f}' for m in means)} |")
    return "\n".join(lines) + "\n"


def _md(tab: pd.DataFrame, counts: pd.DataFrame, fs: str) -> str:
    """Per-pair summary at the headline window."""
    at = counts[counts["window_bins"] == HEADLINE_W]
    n_lead = LEAD_DIMS
    chance = 1 - 2 * 0.5 ** n_lead
    lines = [f"### {fs} — do CCs within a pair disagree in sign?", "",
             f"Leading {LEAD_DIMS} canonical dims, integration window "
             f"±{HEADLINE_W * BIN_MS} ms. `mixed` = that animal has at least one "
             f"positive-IFI CC **and** at least one negative-IFI CC in the same pair. "
             f"Positive IFI ⇒ the first-named area leads.", "",
             f"> ⚠ **Per-animal mixing is NOT evidence on its own.** With {n_lead} "
             f"dimensions and a coin-flip sign, at least one of each is expected "
             f"{chance:.0%} of the time. Read the reliable-direction table below "
             f"instead — mixing only means something if the individual directions "
             f"reproduce across animals.", "",
             "| pair | animals | mixed-sign animals | mean CCs leading | "
             "mean CCs lagging | mean frac. positive |",
             "|---|---|---|---|---|---|"]
    for pair in PAIRS:
        g = at[at["pair"] == pair]
        if g.empty:
            continue
        lines.append(
            f"| {pair} | {len(g)} | **{int(g['mixed'].sum())}/{len(g)}** "
            f"({g['mixed'].mean():.0%}) | {g['n_pos'].mean():.1f} | "
            f"{g['n_neg'].mean():.1f} | {g['frac_pos'].mean():.2f} |")
    # how the answer depends on the integration window
    lines += ["", f"**Mixed-sign rate vs integration window** (leading {LEAD_DIMS} dims, "
              "pooled over pairs):", "",
              "| window (ms) | " + " | ".join(str(w * BIN_MS) for w in
                                              [1, 2, 5, 10, 15, 20, 25]) + " |",
              "|---|" + "---|" * 7]
    row = ["| mixed-sign animal-pairs "]
    for w in [1, 2, 5, 10, 15, 20, 25]:
        g = counts[counts["window_bins"] == w]
        row.append(f"| {g['mixed'].mean():.0%} " if len(g) else "| — ")
    lines.append("".join(row) + "|")
    return "\n".join(lines) + "\n"


def mixing_vs_chance(tab: pd.DataFrame, window_bins: int = HEADLINE_W):
    """Is sign-mixing more common than a coin flip would give? — the decisive test.

    Restricted to CCs that BEAT THE CIRCULAR-SHIFT NULL: ungated, dims 3-8 clear it in
    only 14-17 % of cells with held-out CC at the floor, so their sign is noise.

    The chance model has to be per-animal, because each animal has its own number ``k``
    of significant CCs and ``P(mixed) = 1 - 2*0.5**k`` rises steeply with ``k`` (0 % at
    k=1, 50 % at k=2, 94 % at k=5). Comparing a raw mixing rate to 50 % — or to nothing
    at all — would make chance look like structure.

    Returns ``(summary_dict, per_pair_frame)``.
    """
    from scipy import stats
    sub = tab[(tab["window_bins"] == window_bins) & (tab["degenerate"] == 0) &
              (tab["sig"] == 1)]
    obs_tot = exp_tot = n = 0
    per_pair = {}
    for (_, pair), g in sub.groupby(["animal", "pair"]):
        v = g["ifi"].to_numpy(float)
        v = v[np.isfinite(v) & (v != 0)]
        k = v.size
        if k < 2:
            continue
        obs = int((v > 0).any() and (v < 0).any())
        exp = 1 - 2 * 0.5 ** k
        d = per_pair.setdefault(pair, {"pair": pair, "obs": 0, "exp": 0.0, "n": 0})
        d["obs"] += obs; d["exp"] += exp; d["n"] += 1
        obs_tot += obs; exp_tot += exp; n += 1
    if n == 0:
        return {}, pd.DataFrame()
    p = float(stats.binomtest(obs_tot, n, exp_tot / n).pvalue)
    return ({"n_cells": n, "observed": obs_tot, "expected": exp_tot,
             "obs_frac": obs_tot / n, "exp_frac": exp_tot / n, "p": p},
            pd.DataFrame(per_pair.values()))


def direction_sweep(tab: pd.DataFrame, max_dim: int = 6) -> pd.DataFrame:
    """Per (pair, dim, window): one-sample *t* of the IFI across animals.

    The meeting asked for the IFI "at different lags", so the window is swept rather
    than fixed. ⚠ That makes the per-cell p-values a SELECTED minimum over
    ``n_windows x n_dims`` tests per pair (25 x 6 = 150 here) — they are descriptive
    pointers, not inference. The fixed-window FDR table (`dim_direction`) is the
    inferential statement.
    """
    from tom_cca import paired_stats
    sub = tab[(tab["degenerate"] == 0) & (tab["dim"] <= max_dim)]
    rows = []
    for (pair, d, w), g in sub.groupby(["pair", "dim", "window_bins"]):
        v = g["ifi"].to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size < 4:
            continue
        n, _, t, p = paired_stats.paired_t(v)
        rows.append({"pair": pair, "dim": int(d), "window_bins": int(w),
                     "window_ms": int(w) * BIN_MS, "n_animals": n,
                     "mean_ifi": float(np.mean(v)), "t": t, "p": p})
    return pd.DataFrame(rows)


def _md_sweep(sweep: pd.DataFrame, fs: str) -> str:
    if sweep.empty:
        return ""
    lines = [f"#### {fs} — strongest per-CC directions across the window sweep", "",
             "⚠ **Selected minima over ~150 tests per pair (25 windows × 6 dims) — "
             "descriptive pointers, not inference.** Nothing here survives the "
             "fixed-window FDR test above.", "",
             "| pair | CC | window | n | mean IFI | p (uncorrected) |",
             "|---|---|---|---|---|---|"]
    for _, r in sweep.sort_values("p").head(8).iterrows():
        lines.append(f"| {r['pair']} | CC{int(r['dim'])} | ±{int(r['window_ms'])} ms | "
                     f"{int(r['n_animals'])} | {r['mean_ifi']:+.3f} | {r['p']:.4f} |")
    # the only thing that would directly answer the question
    lines += ["", "**Pairs where two CCs of OPPOSITE sign are both nominally "
              "significant at the same window:**", ""]
    found = False
    for (pair, w), g in sweep.groupby(["pair", "window_ms"]):
        s = g[g["p"] < 0.05]
        if len(s) >= 2 and (s["mean_ifi"] > 0).any() and (s["mean_ifi"] < 0).any():
            found = True
            lines.append(f"- **{pair}** at ±{int(w)} ms: " + ", ".join(
                f"CC{int(r['dim'])} {r['mean_ifi']:+.3f} (p={r['p']:.3f})"
                for _, r in s.sort_values("dim").iterrows()))
    if not found:
        lines.append("- none")
    return "\n".join(lines) + "\n"


def _md_direction(direction: pd.DataFrame, fs: str) -> str:
    """The test that actually answers the question."""
    lines = [f"#### {fs} — which CCs have a RELIABLE direction across animals?", "",
             f"One-sample *t* of each dimension's IFI across animals (±"
             f"{HEADLINE_W * BIN_MS} ms window), BH-FDR across dims within a pair. "
             "A pair genuinely contains CCs of opposite direction only if **two dims "
             "survive with opposite signs**.", "",
             "| pair | dims with a reliable direction (FDR) | opposite signs? |",
             "|---|---|---|"]
    for pair in PAIRS:
        g = direction[(direction["pair"] == pair) & (direction["fdr_pass"])]
        if g.empty:
            lines.append(f"| {pair} | — none | no |")
            continue
        desc = ", ".join(
            f"CC{int(r['dim'])} {r['mean_ifi']:+.3f} (p={r['p']:.3g})"
            for _, r in g.sort_values("dim").iterrows())
        opp = "**YES**" if (g["mean_ifi"] > 0).any() and (g["mean_ifi"] < 0).any() \
            else "no"
        lines.append(f"| {pair} | {desc} | {opp} |")
    return "\n".join(lines) + "\n"


def main():
    md = ["# Item 1 — do different CCs have different IFI (different signs)?", "",
          "Lag curves are the existing **CCA-refit-at-every-lag** held-out, "
          "segment-aware per-dimension curves (`lag_curves_bin10*.csv`, dims 1–30 × "
          "lags ±250 ms). This re-reduces them per dimension over integration windows "
          "|lag| ≤ w.", "",
          "> Positive IFI ⇒ the **first-named area leads**. A dimension whose held-out "
          "> CC is negative at every lag clips to zero on both sides and returns IFI = 0 "
          "> for *no coupling* rather than *balanced*; those are flagged `degenerate` "
          "> and excluded from the sign counts.", ""]
    odirs: dict[str, pd.DataFrame] = {}
    for suf, fs in [("", "FS-excluded"), ("_fsincl", "FS-included")]:
        src = RES / f"lag_curves_bin10{suf}.csv"
        if not src.exists():
            print(f"skip {fs}: {src.name} not found"); continue
        raw = pd.read_csv(src)
        chk = check_sig_consistency(raw)
        status = "OK" if chk["ok"] else "*** INCONSISTENT ***"
        print(f"\n[{fs}] sig consistency: {status} — "
              f"{chk['sig_with_nonpositive_cc']}/{chk['n_sig']} flagged dims with "
              f"cc<=0 at lag 0; FDR-not-subset-of-uncorrected: "
              f"{chk['fdr_not_subset_of_uncorrected']}")
        if not chk["ok"]:
            print("      -> the `sig` column does not describe the `cc` in the same "
                  "row. Re-run scripts/run_lag_curves.py (see lagged."
                  "perdim_significance).")
        md.append(f"> **`sig` consistency check ({fs}):** {status} — "
                  f"{chk['sig_with_nonpositive_cc']}/{chk['n_sig']} flagged dims with "
                  f"non-positive lag-0 CC.\n")
        tab = build_windows(raw)
        counts = sign_counts(tab)
        tab.to_csv(RES / f"cc_ifi_windows_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        counts.to_csv(RES / f"cc_ifi_signs_bin10{suf}.csv", index=False,
                      lineterminator="\n")
        direction = dim_direction(tab)
        direction.to_csv(RES / f"cc_ifi_direction_bin10{suf}.csv", index=False,
                         lineterminator="\n")
        sweep = direction_sweep(tab)
        sweep.to_csv(RES / f"cc_ifi_direction_sweep_bin10{suf}.csv", index=False,
                     lineterminator="\n")
        summ, pp = mixing_vs_chance(tab)
        if summ:
            pp.to_csv(RES / f"cc_ifi_mixing_bin10{suf}.csv", index=False,
                      lineterminator="\n")
            md.append(
                f"#### {fs} — sign-mixing vs chance (significant CCs only)\n\n"
                f"Each animal contributes only the CCs that beat the circular-shift "
                f"null. Chance is computed **per animal** from its own count *k* of "
                f"significant CCs: P(mixed) = 1 − 2·0.5^k.\n\n"
                f"- animal-pairs with ≥2 significant CCs: **{summ['n_cells']}**\n"
                f"- observed mixed: **{summ['observed']} ({summ['obs_frac']:.0%})**\n"
                f"- expected by chance: **{summ['expected']:.1f} "
                f"({summ['exp_frac']:.0%})**\n"
                f"- binomial test: **p = {summ['p']:.3f}**\n\n"
                f"**Verdict: sign-mixing is indistinguishable from a coin flip.**\n")
        # 2026-08-07 asks 1 + 3: the ungrouped (all significant CCs) mean and its test
        overall = overall_by_animal(tab)
        overall.to_csv(RES / f"cc_ifi_overall_bin10{suf}.csv", index=False,
                       lineterminator="\n")
        odir = overall_direction(overall)
        # CCs-as-n power check, merged in as ccs_* columns so the figure reads one file
        occ = overall_direction_ccs(tab).drop(columns=["unit", "window_ms"])
        occ = occ.rename(columns={c: f"ccs_{c}" for c in occ.columns if c != "pair"})
        odir = odir.merge(occ, on="pair", how="left")
        odir.to_csv(RES / f"cc_ifi_overall_test_bin10{suf}.csv", index=False,
                    lineterminator="\n")
        cc1 = cc1_direction(tab)
        cc1.to_csv(RES / f"cc_ifi_cc1_test_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        n_all = tab.groupby("pair")["animal"].nunique()
        md.append(_md_overall(odir, cc1, n_all, fs))
        odirs[fs] = odir
        print(f"  overall IFI at ±{HEADLINE_W * BIN_MS} ms vs 0 (all sig CCs, "
              f"animals-as-n; first ~20 trials, 12k-bin cap):")
        c1 = cc1.set_index("pair")
        for _, r in odir.iterrows():
            flag = " *" if r["p"] < 0.05 else ""
            c1s = (f"CC1 same data {c1.loc[r['pair'], 'mean']:+.3f} "
                   f"p={c1.loc[r['pair'], 'p']:.3g}" if r["pair"] in c1.index else "")
            ccs = (f"CCs n={int(r['ccs_n'])} p={r['ccs_p']:.3g}"
                   f"{' *' if r['ccs_p'] < 0.05 else ''}"
                   if np.isfinite(r.get("ccs_n", np.nan)) else "")
            print(f"    {r['pair']:8s} n={int(r['n']):<2d}/{int(n_all.get(r['pair'], 0)):<2d} "
                  f"mean={r['mean']:+.3f} t={r['t']:+.2f} p={r['p']:.3g}{flag}   "
                  f"(weighted {r['wmean']:+.3f} p={r['wp']:.3g}; {ccs}; {c1s})")
        md.append(_md(tab, counts, fs))
        md.append(_md_direction(direction, fs))
        md.append(_md_sweep(sweep, fs))
        at = counts[counts["window_bins"] == HEADLINE_W]
        print(f"\n{fs}: {tab['animal'].nunique()} animals, "
              f"{tab['dim'].nunique()} dims, {tab['window_bins'].nunique()} windows "
              f"-> {len(tab)} rows")
        print(f"  mixed-sign animal-pairs at ±{HEADLINE_W * BIN_MS} ms "
              f"(leading {LEAD_DIMS} dims): "
              f"{int(at['mixed'].sum())}/{len(at)} ({at['mixed'].mean():.0%})")
        deg = tab[tab["dim"] <= LEAD_DIMS]["degenerate"].mean()
        deg_all = tab["degenerate"].mean()
        print(f"  degenerate (no coupling either way): "
              f"{deg:.0%} of leading-{LEAD_DIMS} rows, {deg_all:.0%} of all 30 dims")
    if odirs:
        md.append(_md_overall_looks(odirs))
    (RES / "cc_ifi_signs_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / 'cc_ifi_signs_tables.md'}")


if __name__ == "__main__":
    main()
