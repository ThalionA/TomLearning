"""Item 2, factorial version — epoch x FF/FB as a 2-way design, not a single contrast.

Replaces the ΔFF − ΔFB difference-of-differences used in `analyze_cc_label_track.py`.
That contrast has two defects Theo flagged:

  * it discards the INTERMEDIATE epoch entirely, using only expert − naive;
  * it collapses main effects and interaction into one number, so "no divergence" and
    "no change at all" are indistinguishable.

Here: a **2-way repeated-measures ANOVA per pair**, epoch (naive / intermediate /
expert) x label (FF / FB), within animal, reporting all three terms —

    effect of EPOCH        does the readout change with learning at all?
    effect of LABEL        do FF- and FB-labelled components differ?
    EPOCH x LABEL          do they change DIFFERENTLY with learning? (item 2's question)

Animals-as-n: an animal's several components in a cell are averaged first
(`aggregate_func="mean"`), so a component-level test cannot pseudoreplicate. Animals
missing any of the 6 cells are dropped — 47 of 50 animal-pairs are complete, so this
costs little.

⚠ **A LABEL main effect on SIGNED IFI is circular.** Components are labelled FF/FB by the
sign of their whole-session IFI (FF mean +0.075, FB -0.090), so that term is guaranteed
and means nothing. It is reported for completeness and marked. The non-circular readouts
are `peak_r`, `abs_ifi` and `peak_lag_ms` — the label was not defined using any of them.

A pooled linear mixed model across pairs is also fitted, since a per-pair ANOVA rests on
4-10 animals: `value ~ epoch * label` with a random intercept per animal.

Writes results/cc_label_anova_bin10{,_fsincl}.csv + cc_label_anova_tables.md

Usage: PYTHONPATH=src python scripts/analyze_cc_label_anova.py
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd


from _common import EPOCHS, FS_CONDITIONS, RES, config
PAIRS = list(config.PAIR_NAMES)
# (column, label, is the LABEL main effect circular for this DV?)
DVS = [("peak_r", "peak r", False),
       ("abs_ifi", "|IFI|", False),
       ("peak_lag_ms", "peak lag (ms)", False),
       ("ifi", "IFI (signed)", True)]
TERMS = ["epoch", "label", "epoch:label"]


def _prep(tab: pd.DataFrame) -> pd.DataFrame:
    t = tab[tab["label"].isin(["FF", "FB"])].copy()
    t["abs_ifi"] = t["ifi"].abs()
    return t


def _complete(sub: pd.DataFrame, dv: str) -> pd.DataFrame:
    """Animals having all 6 epoch x label cells with a finite DV."""
    g = sub.dropna(subset=[dv])
    keep = [a for a, ga in g.groupby("animal")
            if len(set(zip(ga["epoch"], ga["label"]))) == len(EPOCHS) * 2]
    return g[g["animal"].isin(keep)]


def anova_per_pair(tab: pd.DataFrame) -> pd.DataFrame:
    from statsmodels.stats.anova import AnovaRM
    rows = []
    for pair in PAIRS:
        sub = tab[tab["pair"] == pair]
        if sub.empty:
            continue
        for dv, _, circ in DVS:
            g = _complete(sub, dv)
            n = g["animal"].nunique()
            row = {"pair": pair, "dv": dv, "n_animals": n,
                   "label_effect_circular": circ}
            if n >= 3:
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        fit = AnovaRM(data=g, depvar=dv, subject="animal",
                                      within=["epoch", "label"],
                                      aggregate_func="mean").fit()
                    tb = fit.anova_table
                    for term in TERMS:
                        if term in tb.index:
                            row[f"F_{term}"] = float(tb.loc[term, "F Value"])
                            row[f"p_{term}"] = float(tb.loc[term, "Pr > F"])
                except Exception as exc:                    # noqa: BLE001
                    row["error"] = type(exc).__name__
            rows.append(row)
    return pd.DataFrame(rows)


def lmm_pooled(tab: pd.DataFrame) -> pd.DataFrame:
    """Pooled across pairs: value ~ epoch * label, random intercept per animal.

    A per-pair ANOVA rests on 4-10 animals. Pooling buys power but the pairs share
    animals and areas, so `pair` enters as a fixed covariate and the animal random
    intercept absorbs the repeated measurement.
    """
    import statsmodels.formula.api as smf
    rows = []
    for dv, _, circ in DVS:
        g = tab.dropna(subset=[dv]).copy()
        g = g.groupby(["animal", "pair", "epoch", "label"], as_index=False)[dv].mean()
        if g["animal"].nunique() < 4:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = smf.mixedlm(f"{dv} ~ C(epoch) * C(label) + C(pair)", g,
                                groups=g["animal"]).fit(reml=False, method="lbfgs")
            # Wald joint test per term via the fitted parameter names
            for term, patt in [("epoch", "C(epoch)"), ("label", "C(label)"),
                               ("epoch:label", "C(epoch)[T.")]:
                if term == "epoch:label":
                    names = [p for p in m.params.index
                             if "C(epoch)" in p and "C(label)" in p]
                else:
                    names = [p for p in m.params.index
                             if patt in p and not ("C(epoch)" in p
                                                   and "C(label)" in p)]
                if not names:
                    continue
                R = np.zeros((len(names), len(m.params)))
                for i, nm in enumerate(names):
                    R[i, list(m.params.index).index(nm)] = 1.0
                wt = m.wald_test(R, scalar=True)
                rows.append({"dv": dv, "term": term,
                             "chi2": float(np.squeeze(wt.statistic)),
                             "p": float(np.squeeze(wt.pvalue)),
                             "n_animals": int(g["animal"].nunique()),
                             "n_obs": int(len(g)),
                             "label_effect_circular": circ})
        except Exception as exc:                            # noqa: BLE001
            rows.append({"dv": dv, "term": "fit_failed",
                         "chi2": np.nan, "p": np.nan,
                         "n_animals": int(g["animal"].nunique()),
                         "n_obs": int(len(g)), "label_effect_circular": circ,
                         "error": type(exc).__name__})
    return pd.DataFrame(rows)


def _md(av: pd.DataFrame, lmm: pd.DataFrame, fs: str) -> str:
    lines = [f"### {fs} — 2-way RM-ANOVA per pair: epoch × FF/FB label", "",
             "Animals-as-n (an animal's components in a cell are averaged first); animals "
             "missing any of the 6 cells are dropped. **The label main effect on signed "
             "IFI is circular by construction** — components are labelled by the sign of "
             "their whole-session IFI — and is marked ⚠.", ""]
    for dv, name, circ in DVS:
        sub = av[av["dv"] == dv]
        if sub.empty:
            continue
        lines += [f"**{name}**", "",
                  "| pair | n | F epoch | p | F label | p | F epoch×label | p |",
                  "|---|---|---|---|---|---|---|---|"]
        for _, r in sub.iterrows():
            def g(k):
                v = r.get(k, np.nan)
                return f"{v:.3g}" if isinstance(v, float) and np.isfinite(v) else "—"
            mark = " ⚠" if circ else ""
            lines.append(
                f"| {r['pair']} | {r['n_animals']:.0f} | {g('F_epoch')} | "
                f"{g('p_epoch')} | {g('F_label')}{mark} | {g('p_label')}{mark} | "
                f"{g('F_epoch:label')} | {g('p_epoch:label')} |")
        lines.append("")
    if not lmm.empty:
        lines += ["", "#### Pooled across pairs — LMM `value ~ epoch * label + pair`, "
                  "random intercept per animal", "",
                  "| DV | term | χ² | p | animals | obs |", "|---|---|---|---|---|---|"]
        for _, r in lmm.iterrows():
            mark = " ⚠" if (r["label_effect_circular"] and r["term"] == "label") else ""
            lines.append(f"| {r['dv']} | {r['term']}{mark} | {r['chi2']:.3g} | "
                         f"{r['p']:.3g} | {r['n_animals']:.0f} | {r['n_obs']:.0f} |")
    return "\n".join(lines) + "\n"


def main():
    md = ["# Item 2 — factorial: effect of epoch, effect of FF/FB, and their interaction",
          "",
          "Replaces the ΔFF − ΔFB contrast, which discarded the intermediate epoch and "
          "could not separate main effects from the interaction.", ""]
    for suf, fs in FS_CONDITIONS:
        src = RES / f"cc_label_track_epoch_bin10{suf}.csv"
        if not src.exists():
            print(f"skip {fs}: {src.name} not found"); continue
        tab = _prep(pd.read_csv(src))
        av = anova_per_pair(tab)
        lmm = lmm_pooled(tab)
        av.to_csv(RES / f"cc_label_anova_bin10{suf}.csv", index=False,
                  lineterminator="\n")
        lmm.to_csv(RES / f"cc_label_anova_pooled_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        md.append(_md(av, lmm, fs))
        print(f"\n=== {fs} ===")
        for dv, name, circ in DVS:
            sub = av[av["dv"] == dv]
            if sub.empty:
                continue
            bits = []
            for term in TERMS:
                col = f"p_{term}"
                if col not in sub:
                    continue
                hits = sub[sub[col] < 0.05]
                tag = " (CIRCULAR)" if (circ and term == "label") else ""
                bits.append(f"{term}{tag}: {len(hits)}/{sub[col].notna().sum()}")
            print(f"  {name:14s} pairs at p<0.05 — " + " | ".join(bits))
        if not lmm.empty:
            print("  pooled LMM:")
            for _, r in lmm.iterrows():
                tag = " (CIRCULAR)" if (r["label_effect_circular"]
                                        and r["term"] == "label") else ""
                print(f"     {r['dv']:12s} {r['term']:12s} p={r['p']:.3g}{tag}")
    (RES / "cc_label_anova_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / 'cc_label_anova_tables.md'}")


if __name__ == "__main__":
    main()
