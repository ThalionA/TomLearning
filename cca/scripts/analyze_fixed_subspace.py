"""Fixed-subspace epoch contrasts — meeting 2026-07-28 items 5, 6 and 7.

Reads results/fixed_subspace_bin10{,_fsincl}.csv (written by run_fixed_subspace.py):
the lagged correlation of ONE canonical component, projected through a subspace that was
identified once across balanced trials and then frozen.

Per (animal, pair, epoch) the CC1 lag curve is reduced to four readouts:
    peak_r        peak correlation (IN-SAMPLE — a contrast statistic, not a coupling
                  strength; see the module docstring of tom_cca.fixed_subspace)
    peak_lag_ms   where it peaks       (+ = the first-named area leads)
    ifi           asymmetry of the curve
    width_ms      half-max integration window (meeting item 6)

Then naive vs expert, animals-as-n paired t, per pair. Because the subspace is frozen and
balanced across epochs, a difference here is a difference in the ACTIVITY rather than in
the fit — the confound in every refit-per-epoch contrast run so far.

Writes results/fixed_subspace_epoch_bin10{,_fsincl}.csv + fixed_subspace_tables.md

Usage: PYTHONPATH=src python scripts/analyze_fixed_subspace.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import fixed_subspace, paired_stats, perdim_ifi  # noqa: E402

RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
EPOCHS = ["naive", "intermediate", "expert"]
METRICS = [("peak_r", "peak r"), ("ifi", "IFI"), ("peak_lag_ms", "peak lag (ms)"),
           ("width_ms", "half-max width (ms)")]


def reduce_curves(df: pd.DataFrame, dim: int = 1) -> pd.DataFrame:
    """One row per (animal, pair, epoch) from the CC`dim` lag curve."""
    rows = []
    sub = df[df["dim"] == dim]
    for (animal, pair, epoch), g in sub.groupby(["animal", "pair", "epoch"]):
        g = g.sort_values("lag_ms")
        lag = g["lag_ms"].to_numpy(float)
        r = g["r"].to_numpy(float)
        if not np.any(np.isfinite(r)):
            continue
        side_ms, side_ratio = fixed_subspace.side_peak(lag, r)
        rows.append({
            "animal": animal, "pair": pair, "epoch": epoch,
            "peak_r": float(np.nanmax(r)),
            "peak_lag_ms": perdim_ifi.curve_peak_lag(lag, r),
            "ifi": perdim_ifi.curve_ifi(lag, r),
            "width_ms": fixed_subspace.curve_half_width(lag, r),
            "side_peak_ms": side_ms,
            "side_peak_ratio": side_ratio,
            "n_bins": int(g["n_bins"].iloc[0]),
            "r_fit": float(g["r_fit"].iloc[0]) if pd.notna(g["r_fit"].iloc[0])
                     else np.nan,
        })
    return pd.DataFrame(rows)


def epoch_contrasts(tab: pd.DataFrame) -> pd.DataFrame:
    """Expert vs naive, animals-as-n paired t, per pair and metric."""
    out = []
    for pair in PAIRS:
        sub = tab[tab["pair"] == pair]
        if sub.empty:
            continue
        piv = {e: sub[sub["epoch"] == e].set_index("animal") for e in EPOCHS}
        common = sorted(set(piv["naive"].index) & set(piv["expert"].index))
        row = {"pair": pair, "n_animals": len(common)}
        for metric, _ in METRICS:
            if not common:
                row |= {f"{metric}_naive": np.nan, f"{metric}_expert": np.nan,
                        f"d_{metric}": np.nan, f"p_{metric}": np.nan}
                continue
            nv = piv["naive"].loc[common, metric].to_numpy(float)
            ex = piv["expert"].loc[common, metric].to_numpy(float)
            d = ex - nv
            d = d[np.isfinite(d)]
            if d.size >= 2:
                _, _, t, p = paired_stats.paired_t(d)
            else:
                t = p = np.nan
            row |= {f"{metric}_naive": float(np.nanmean(nv)),
                    f"{metric}_expert": float(np.nanmean(ex)),
                    f"d_{metric}": float(np.nanmean(d)) if d.size else np.nan,
                    f"t_{metric}": t, f"p_{metric}": p}
        out.append(row)
    return pd.DataFrame(out)


def _md(stats: pd.DataFrame, fs: str) -> str:
    lines = [f"### {fs} — expert vs naive through a FROZEN subspace", "",
             "Subspace identified once on trials balanced across epochs, then both "
             "epochs projected through the identical weights; animals-as-n paired *t*. "
             "`peak r` is in-sample by construction — read the Δ, not the level.", "",
             "| pair | n | " + " | ".join(
                 f"{lab} naive | {lab} exp | Δ | p" for _, lab in METRICS) + " |",
             "|---|---|" + "---|" * (4 * len(METRICS))]
    for _, r in stats.iterrows():
        cells = [f"{r['pair']}", f"{r['n_animals']:.0f}"]
        for metric, _ in METRICS:
            cells += [f"{r[f'{metric}_naive']:.3g}", f"{r[f'{metric}_expert']:.3g}",
                      f"{r[f'd_{metric}']:+.3g}", f"{r[f'p_{metric}']:.3g}"]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def _md_ringing(tab: pd.DataFrame, fs: str) -> str:
    """How to read the integration window, per pair — decaying or ringing?"""
    lines = [f"#### {fs} — is the lag curve decaying or RINGING?", "",
             "A secondary peak above half the central peak means the curve oscillates "
             "rather than decays. On a ringing curve the half-max width is the "
             "**half-period of the rhythm, not a temporal integration window**, and the "
             "IFI compares two lobes of that rhythm rather than a lead.", "",
             "| pair | curves with a side peak | side-peak lag | implied freq | "
             "mean width | read the width as |",
             "|---|---|---|---|---|---|"]
    for pair, g in tab.groupby("pair"):
        ring = g[(g["side_peak_ratio"] > 0.5) & g["side_peak_ms"].notna()]
        frac = len(ring) / len(g) if len(g) else np.nan
        off = float(ring["side_peak_ms"].mean()) if len(ring) else np.nan
        freq = f"{1000 / off:.1f} Hz" if np.isfinite(off) and off > 0 else "—"
        # The verdict has to follow the MAJORITY of curves. A pair where a third of
        # animals ring is not a pair whose width is a half-period.
        if not np.isfinite(frac) or frac == 0:
            verdict = "integration window"
        elif frac >= 0.5:
            verdict = "**half-period of a rhythm**"
        else:
            verdict = "mixed — mostly integration window"
        lines.append(
            f"| {pair} | {len(ring)}/{len(g)} ({frac:.0%}) | "
            f"{'%.0f ms' % off if np.isfinite(off) else '—'} | {freq} | "
            f"{g['width_ms'].mean():.0f} ms | {verdict} |")
    return "\n".join(lines) + "\n"


def main():
    md = ["# Fixed-subspace lag curves by epoch — meeting items 5, 6, 7", "",
          "One canonical component, identified once across balanced trials and frozen, "
          "then lagged across time within each learning epoch. Because the weights are "
          "identical across epochs, an epoch difference is a difference in the "
          "**activity**, not in the fit.", "",
          "> **Not a coupling strength.** The frozen fit saw every epoch, so `peak r` is "
          "> in-sample and optimistic. It is a contrast statistic only; the leak-free "
          "> numbers live in `lag_subspaces_tables.md`.", ""]
    for suf, fs in [("", "FS-excluded"), ("_fsincl", "FS-included")]:
        src = RES / f"fixed_subspace_bin10{suf}.csv"
        if not src.exists():
            print(f"skip {fs}: {src.name} not found"); continue
        df = pd.read_csv(src)
        df["r"] = pd.to_numeric(df["r"], errors="coerce")
        tab = reduce_curves(df)
        stats = epoch_contrasts(tab)
        tab.to_csv(RES / f"fixed_subspace_epoch_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        stats.to_csv(RES / f"fixed_subspace_stats_bin10{suf}.csv", index=False,
                     lineterminator="\n")
        md.append(_md(stats, fs))
        md.append(_md_ringing(tab, fs))
        print(f"\n{fs}: {tab['animal'].nunique()} animals, {len(tab)} "
              f"(animal,pair,epoch) curves")
        ring = tab.dropna(subset=["side_peak_ms"])
        ring = ring[ring["side_peak_ratio"] > 0.5]
        if len(ring):
            per_pair = ring.groupby("pair")["side_peak_ms"].mean()
            print("  RINGING curves (width = half-period, NOT an integration window): "
                  + ", ".join(f"{p} ~{v:.0f} ms ({1000 / v:.1f} Hz)"
                              for p, v in per_pair.items()))
        for metric, lab in METRICS:
            hits = stats[stats[f"p_{metric}"] < 0.05]
            print(f"  {lab}: {len(hits)}/{len(stats)} pairs at p<0.05" +
                  ("  — " + ", ".join(f"{r['pair']} Δ={r[f'd_{metric}']:+.3g} "
                                      f"p={r[f'p_{metric}']:.3g}"
                                      for _, r in hits.iterrows()) if len(hits) else ""))
    (RES / "fixed_subspace_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / 'fixed_subspace_tables.md'}")


if __name__ == "__main__":
    main()
