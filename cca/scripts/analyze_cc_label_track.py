"""Item 2 — what do the FF- and FB-labelled CCs do across learning?

Reads results/cc_label_track_bin10{,_fsincl}.csv (run_cc_label_track.py): each canonical
dimension labelled ONCE from the whole-session fit, then followed through
naive/intermediate/expert on frozen axes.

Per (animal, pair, dim, epoch) the lagged-correlation curve is reduced to:
    peak_r      peak correlation of that CC in that epoch  (IN-SAMPLE — a contrast
                statistic, not a coupling strength)
    ifi         directional asymmetry of the curve, +/-LABEL_W window
    peak_lag_ms where the curve peaks   (+ = the first-named area leads)

Then three questions, animals-as-n:

  Q1  Do the FF CCs and the FB CCs change DIFFERENTLY with learning?
      The interaction: (expert-naive) for FF CCs vs (expert-naive) for FB CCs.
      This is the question item 2 actually asks.
  Q2  Does either group change at all? expert-naive within each label.
  Q3  Do the labelled CCs KEEP their direction across epochs, or does the sign flip?
      A label assigned once should persist if it reflects something real.

⚠ Each animal contributes several CCs, so a per-CC test would be pseudoreplicated.
Everything below therefore collapses to ONE value per animal per label group first, and
n = animals.

⚠ Only CCs that beat the per-dim held-out null are used (item 1: ungated, the extreme
IFIs sit on the weakest CCs).

Writes results/cc_label_track_epoch_bin10{,_fsincl}.csv + cc_label_track_tables.md

Usage: PYTHONPATH=src python scripts/analyze_cc_label_track.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import lagged, paired_stats, perdim_ifi  # noqa: E402

RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
EPOCHS = ["naive", "intermediate", "expert"]
LABEL_W = 5           # bins; must match run_cc_label_track.LABEL_W
METRICS = [("peak_r", "peak r"), ("ifi", "IFI"), ("peak_lag_ms", "peak lag (ms)")]


def reduce_curves(df: pd.DataFrame, sig_only: bool = True) -> pd.DataFrame:
    """One row per (animal, pair, dim, epoch)."""
    if sig_only:
        df = df[df["sig"] == 1]
    rows = []
    for (animal, pair, dim, epoch), g in df.groupby(["animal", "pair", "dim", "epoch"]):
        g = g.sort_values("lag_ms")
        lag = g["lag_ms"].to_numpy(float)
        r = pd.to_numeric(g["r"], errors="coerce").to_numpy(float)
        if not np.any(np.isfinite(r)):
            continue
        m = np.abs(lag) <= LABEL_W * int(g["bin_ms"].iloc[0])
        rows.append({
            "animal": animal, "pair": pair, "dim": int(dim), "epoch": epoch,
            "label": g["label"].iloc[0],
            "ifi_fit": float(g["ifi_fit"].iloc[0]),
            "cc_heldout": pd.to_numeric(g["cc_heldout"], errors="coerce").iloc[0],
            "peak_r": float(np.nanmax(r)),
            "ifi": float(lagged.information_flow_index(lag[m], r[m])),
            "peak_lag_ms": perdim_ifi.curve_peak_lag(lag, r),
        })
    return pd.DataFrame(rows)


def _per_animal(tab: pd.DataFrame, metric: str, label: str) -> pd.DataFrame:
    """Collapse an animal's CCs of one label to a single value per epoch.

    Without this each animal would contribute up to 10 CCs and a paired test would be
    pseudoreplicated — the trap that produced the landmark arm's inflated p-values
    (STATE.md §3.1).
    """
    sub = tab[tab["label"] == label]
    return sub.pivot_table(index=["animal", "pair"], columns="epoch", values=metric,
                           aggfunc="mean")


def epoch_tests(tab: pd.DataFrame) -> pd.DataFrame:
    """Q1/Q2 per pair and metric: does each label group change, and do they differ?"""
    out = []
    for pair in PAIRS:
        sub = tab[tab["pair"] == pair]
        if sub.empty:
            continue
        row = {"pair": pair}
        for metric, _ in METRICS:
            deltas = {}
            for label in ("FF", "FB"):
                piv = _per_animal(sub, metric, label)
                if piv.empty or not {"naive", "expert"} <= set(piv.columns):
                    deltas[label] = np.array([])
                    row |= {f"n_{label}_{metric}": 0, f"d_{label}_{metric}": np.nan,
                            f"p_{label}_{metric}": np.nan}
                    continue
                d = (piv["expert"] - piv["naive"]).to_numpy(float)
                d = d[np.isfinite(d)]
                deltas[label] = d
                n, mean, t, p = (paired_stats.paired_t(d) if d.size >= 2
                                 else (d.size, np.nan, np.nan, np.nan))
                row |= {f"n_{label}_{metric}": n,
                        f"d_{label}_{metric}": float(np.mean(d)) if d.size else np.nan,
                        f"p_{label}_{metric}": p}
            # Q1 interaction: paired within animal-pair where BOTH labels are present
            pf = _per_animal(sub, metric, "FF")
            pb = _per_animal(sub, metric, "FB")
            inter = np.array([])
            if (not pf.empty and not pb.empty
                    and {"naive", "expert"} <= set(pf.columns)
                    and {"naive", "expert"} <= set(pb.columns)):
                common = pf.index.intersection(pb.index)
                if len(common):
                    a = (pf.loc[common, "expert"] - pf.loc[common, "naive"]).to_numpy(
                        float)
                    b = (pb.loc[common, "expert"] - pb.loc[common, "naive"]).to_numpy(
                        float)
                    inter = (a - b)[np.isfinite(a - b)]
            n, mean, t, p = (paired_stats.paired_t(inter) if inter.size >= 2
                             else (inter.size, np.nan, np.nan, np.nan))
            row |= {f"n_int_{metric}": n,
                    f"d_int_{metric}": float(np.mean(inter)) if inter.size else np.nan,
                    f"p_int_{metric}": p}
        out.append(row)
    return pd.DataFrame(out)


def sign_persistence(tab: pd.DataFrame) -> dict:
    """Q3: does a CC keep the direction it was labelled with?

    For each labelled CC, the fraction of epochs whose own IFI sign matches the label.
    A label assigned once from the whole session should persist if it reflects something
    stable; chance is 0.5.
    """
    keep = []
    for _, g in tab.groupby(["animal", "pair", "dim"]):
        lab = g["label"].iloc[0]
        if lab not in ("FF", "FB"):
            continue
        want = 1 if lab == "FF" else -1
        v = g["ifi"].to_numpy(float)
        v = v[np.isfinite(v) & (v != 0)]
        if v.size == 0:
            continue
        keep.append(float(np.mean(np.sign(v) == want)))
    if not keep:
        return {}
    arr = np.array(keep)
    n, mean, t, p = paired_stats.paired_t(arr - 0.5)
    return {"n_ccs": arr.size, "mean_persistence": float(arr.mean()),
            "p_vs_chance": p}


def _md(stats: pd.DataFrame, persist: dict, tab: pd.DataFrame, fs: str) -> str:
    n_ff = (tab.drop_duplicates(["animal", "pair", "dim"])["label"] == "FF").sum()
    n_fb = (tab.drop_duplicates(["animal", "pair", "dim"])["label"] == "FB").sum()
    lines = [f"### {fs} — do FF- and FB-labelled CCs change differently with learning?",
             "",
             f"{n_ff} FF and {n_fb} FB canonical dimensions (significant CCs only), "
             "labelled once on the whole-session fit and followed through the epochs on "
             "frozen axes. Animals-as-n: an animal's CCs of one label are averaged "
             "before testing. `Δ` is expert − naive; `interaction` is ΔFF − ΔFB, which "
             "is the question item 2 asks.", "",
             "| pair | metric | n | Δ FF | p | Δ FB | p | interaction | p |",
             "|---|---|---|---|---|---|---|---|---|"]
    for _, r in stats.iterrows():
        for metric, lab in METRICS:
            lines.append(
                f"| {r['pair']} | {lab} | {r[f'n_int_{metric}']:.0f} | "
                f"{r[f'd_FF_{metric}']:+.3g} | {r[f'p_FF_{metric}']:.3g} | "
                f"{r[f'd_FB_{metric}']:+.3g} | {r[f'p_FB_{metric}']:.3g} | "
                f"{r[f'd_int_{metric}']:+.3g} | {r[f'p_int_{metric}']:.3g} |")
    if persist:
        lines += ["", f"**Label persistence.** Across {persist['n_ccs']} labelled CCs, "
                  f"an epoch's own IFI sign matches the CC's label "
                  f"{persist['mean_persistence']:.0%} of the time "
                  f"(chance 50 %, p = {persist['p_vs_chance']:.3g}). A label that does "
                  f"not persist is not a property of the CC.", ""]
    return "\n".join(lines) + "\n"


def main():
    md = ["# Item 2 — FF/FB CCs labelled once, tracked across learning", "",
          "Canonical dimensions labelled by the sign of their IFI on the **whole-session "
          "fit**, then followed through naive/intermediate/expert through **frozen "
          "axes** — so a change is a change in the activity, not in the fit, and the "
          "label cannot be re-drawn by noise each epoch.", "",
          "> Correlations are IN-SAMPLE by construction (the frozen fit saw every "
          "> epoch). They are contrast statistics, not coupling strengths.", ""]
    for suf, fs in [("", "FS-excluded"), ("_fsincl", "FS-included")]:
        src = RES / f"cc_label_track_bin10{suf}.csv"
        if not src.exists():
            print(f"skip {fs}: {src.name} not found"); continue
        raw = pd.read_csv(src)
        tab = reduce_curves(raw)
        if tab.empty:
            print(f"skip {fs}: no significant CCs"); continue
        tab.to_csv(RES / f"cc_label_track_epoch_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        stats = epoch_tests(tab)
        persist = sign_persistence(tab)
        stats.to_csv(RES / f"cc_label_track_stats_bin10{suf}.csv", index=False,
                     lineterminator="\n")
        md.append(_md(stats, persist, tab, fs))
        u = tab.drop_duplicates(["animal", "pair", "dim"])
        print(f"\n{fs}: {tab['animal'].nunique()} animals, {len(u)} labelled CCs "
              f"({(u['label'] == 'FF').sum()} FF / {(u['label'] == 'FB').sum()} FB)")
        if persist:
            print(f"  label persists across epochs "
                  f"{persist['mean_persistence']:.0%} of the time "
                  f"(chance 50%, p={persist['p_vs_chance']:.3g})")
        for metric, lab in METRICS:
            hits = stats[stats[f"p_int_{metric}"] < 0.05]
            print(f"  interaction ΔFF−ΔFB on {lab}: {len(hits)}/{len(stats)} pairs"
                  + ("  — " + ", ".join(
                      f"{r['pair']} {r[f'd_int_{metric}']:+.3g} "
                      f"(p={r[f'p_int_{metric}']:.3g})"
                      for _, r in hits.iterrows()) if len(hits) else ""))
    (RES / "cc_label_track_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / 'cc_label_track_tables.md'}")


if __name__ == "__main__":
    main()
