"""Item 4 — integration window vs IFI, naive vs experienced, on FROZEN axes.

    "Plot lagged curves and integration windows for naive vs exp."
    "Where are the integration window vs ifi curves, for naive vs experienced?"

Reads results/cc_label_track_bin10{,_fsincl}.csv (run_cc_label_track.py): one CCA fit per
(animal, pair) on every running bin of every trial, frozen, then each epoch projected
through the identical weights and lagged across +/-250 ms.

Frozen axes rather than refit-per-lag, deliberately. Refitting gives a different
component at every lag AND every epoch: measured on this data, the best-matching
dimension at another lag is a DIFFERENT dimension 82% of the time, and split-half
reproducibility of a canonical vector at a FIXED lag is only |cos| = 0.59 for CC1 and
~0.2 (near the 0.146 random baseline) by CC4. A naive-vs-expert contrast built on
refit-per-lag IFI would therefore compare different components across the very
conditions being contrasted. With frozen axes the component is the same by construction.

For each (animal, pair, dim, epoch): IFI over |lag| <= w for w = 1..25 bins (10..250 ms)
via `perdim_ifi.ifi_windows_by_dim`. Positive IFI => the FIRST-named area leads.

Two readouts:
  * CC1 only — the dominant component, and the only one with real identity in this data.
  * all significant CCs grouped by their FF/FB label (item 2 showed the label persists
    across epochs at 62-66% against a genuine 50% chance). Signed IFI is NOT averaged
    across CCs of opposite label, which would cancel by construction.

Unit of analysis is the ANIMAL: an animal's CCs are averaged before any test.

⚠ Correlations are IN-SAMPLE (the frozen fit saw every epoch). They are contrast
statistics, not coupling strengths.

Writes results/ifi_windows_epochs_bin10{,_fsincl}.csv + ifi_windows_epochs_tables.md

Usage: PYTHONPATH=src python scripts/analyze_ifi_windows_epochs.py
"""
from __future__ import annotations


import numpy as np
import pandas as pd


from _common import FS_CONDITIONS, RES, TEMPORAL, config
from tom_cca import paired_stats, perdim_ifi
PAIRS = list(config.PAIR_NAMES)
BIN_MS = TEMPORAL.bin_ms
HEADLINE_W = TEMPORAL.label_w_bins          # bins = +/-50 ms
EXPECTED_ANIMALS = config.EXPECTED_LEARNERS


def build_windows(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (animal, pair, dim, epoch, integration window)."""
    df = df.copy()
    df["r"] = pd.to_numeric(df["r"], errors="coerce")
    rows = []
    for (animal, pair, dim, epoch), g in df.groupby(["animal", "pair", "dim", "epoch"]):
        g = g.sort_values("lag_bins")
        dims, wins, ifi, degen = perdim_ifi.ifi_windows_by_dim(
            g["lag_bins"].to_numpy(), g["r"].to_numpy(),
            np.full(len(g), int(dim)))
        if dims.size == 0:
            continue
        for j, w in enumerate(wins):
            rows.append({
                "animal": animal, "pair": pair, "dim": int(dim), "epoch": epoch,
                "window_bins": int(w), "window_ms": int(w) * BIN_MS,
                "ifi": ifi[0, j], "degenerate": int(degen[0, j]),
                "label": g["label"].iloc[0],
                "sig": int(g["sig"].iloc[0]),
            })
    return pd.DataFrame(rows)


def epoch_contrast(tab: pd.DataFrame, dim_filter=None, label=None) -> pd.DataFrame:
    """Expert - naive at each integration window, animals-as-n paired t."""
    sub = tab[(tab["sig"] == 1) & (tab["degenerate"] == 0)]
    if dim_filter is not None:
        sub = sub[sub["dim"] == dim_filter]
    if label is not None:
        sub = sub[sub["label"] == label]
    out = []
    for pair in PAIRS:
        g = sub[sub["pair"] == pair]
        if g.empty:
            continue
        for w in sorted(g["window_bins"].unique()):
            gw = g[g["window_bins"] == w]
            # collapse an animal's CCs FIRST — otherwise a paired test over CCs is
            # pseudoreplicated (STATE.md §3.1)
            piv = gw.pivot_table(index="animal", columns="epoch", values="ifi",
                                 aggfunc="mean")
            if not {"naive", "expert"} <= set(piv.columns):
                continue
            d = (piv["expert"] - piv["naive"]).to_numpy(float)
            d = d[np.isfinite(d)]
            n, mean, t, p = (paired_stats.paired_t(d) if d.size >= 3
                             else (d.size, np.nan, np.nan, np.nan))
            out.append({
                "pair": pair, "window_bins": int(w), "window_ms": int(w) * BIN_MS,
                "n_animals": n, "ifi_naive": float(piv["naive"].mean()),
                "ifi_expert": float(piv["expert"].mean()),
                "delta": float(np.mean(d)) if d.size else np.nan, "t": t, "p": p,
            })
    return pd.DataFrame(out)


def _md(cc1: pd.DataFrame, fs: str) -> str:
    at = cc1[cc1["window_bins"] == HEADLINE_W]
    lines = [f"### {fs} — CC1 IFI, expert vs naive, at ±{HEADLINE_W * BIN_MS} ms", "",
             "Frozen axes; animals-as-n. Positive IFI ⇒ the first-named area leads.", "",
             "| pair | n | IFI naive | IFI expert | Δ | p |",
             "|---|---|---|---|---|---|"]
    for _, r in at.iterrows():
        lines.append(f"| {r['pair']} | {r['n_animals']:.0f} | {r['ifi_naive']:+.3f} | "
                     f"{r['ifi_expert']:+.3f} | {r['delta']:+.3f} | {r['p']:.3g} |")
    # does any window show a difference? report the best per pair, flagged as selected
    lines += ["", "**Swept across all 25 windows** (the smallest p per pair — a SELECTED "
              "minimum over 25 nested, highly correlated windows, so it is a pointer, "
              "not a test):", "",
              "| pair | best window | Δ | p (uncorrected) |", "|---|---|---|---|"]
    for pair, g in cc1.groupby("pair"):
        g = g.dropna(subset=["p"])
        if g.empty:
            continue
        b = g.loc[g["p"].idxmin()]
        lines.append(f"| {pair} | ±{b['window_ms']:.0f} ms | {b['delta']:+.3f} | "
                     f"{b['p']:.3g} |")
    return "\n".join(lines) + "\n"


def main():
    md = ["# Item 4 — integration window vs IFI, naive vs experienced", "",
          "Frozen axes (one CCA per animal-pair on all trials, each epoch projected "
          "through identical weights), so the same component is compared across epochs "
          "and windows. IFI integrated over |lag| ≤ w for w = 10…250 ms.", "",
          "> Refit-per-lag was rejected for this contrast: the best-matching dimension "
          "> at another lag is a different dimension 82 % of the time, and split-half "
          "> |cos| of a canonical vector is 0.59 (CC1) falling to ~0.2 by CC4 against a "
          "> 0.146 random baseline — so a refit contrast would compare different "
          "> components across the conditions being contrasted.", "",
          "> Correlations are IN-SAMPLE by construction; these are contrast statistics.",
          ""]
    for suf, fs in FS_CONDITIONS:
        src = RES / f"cc_label_track_bin10{suf}.csv"
        if not src.exists():
            print(f"skip {fs}: {src.name} not found"); continue
        raw = pd.read_csv(src)
        n_an = raw["animal"].nunique()
        if n_an < EXPECTED_ANIMALS:
            print(f"  ⚠ {fs}: only {n_an}/{EXPECTED_ANIMALS} animals — PARTIAL")
        tab = build_windows(raw)
        tab.to_csv(RES / f"ifi_windows_epochs_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        cc1 = epoch_contrast(tab, dim_filter=1)
        cc1.to_csv(RES / f"ifi_windows_epochs_cc1_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        md.append(_md(cc1, fs))
        print(f"\n{fs}: {tab['animal'].nunique()} animals, "
              f"{tab['window_bins'].nunique()} windows")
        at = cc1[cc1["window_bins"] == HEADLINE_W]
        hits = at[at["p"] < 0.05]
        print(f"  CC1 expert-vs-naive at ±{HEADLINE_W * BIN_MS} ms: "
              f"{len(hits)}/{len(at)} pairs at p<0.05"
              + ("  — " + ", ".join(f"{r['pair']} {r['delta']:+.3f} (p={r['p']:.3g})"
                                    for _, r in hits.iterrows()) if len(hits) else ""))
        anyw = cc1[cc1["p"] < 0.05]
        print(f"  any window: {anyw['pair'].nunique()}/{cc1['pair'].nunique()} pairs "
              f"have at least one window at p<0.05 "
              f"(expected ~{0.05 * cc1['window_bins'].nunique():.1f} windows per pair "
              f"by chance)")
        for label in ("FF", "FB"):
            lc = epoch_contrast(tab, label=label)
            lat = lc[lc["window_bins"] == HEADLINE_W]
            lh = lat[lat["p"] < 0.05]
            print(f"  {label}-labelled CCs at ±{HEADLINE_W * BIN_MS} ms: "
                  f"{len(lh)}/{len(lat)} pairs at p<0.05")
            lc.to_csv(RES / f"ifi_windows_epochs_{label}_bin10{suf}.csv", index=False,
                      lineterminator="\n")
    (RES / "ifi_windows_epochs_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / 'ifi_windows_epochs_tables.md'}")


if __name__ == "__main__":
    main()
