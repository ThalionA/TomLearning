"""Significant-canonical-dimensions-as-n analysis (the Gonzalez/Buzsáki unit).

Instead of the animal as the experimental unit, pool every SIGNIFICANT canonical
dimension across animals and treat each dim as an independent sample — the unit
used in Gonzalez & Buzsáki 2026 (e.g. directionality, n = canonical dimensions).
We then compare the per-dimension held-out CC between learning epochs with a
**rank-sum (Mann–Whitney U)** test, per pair, and print it ALONGSIDE the honest
per-animal n, so the inflation is explicit.

CAVEAT (stated loudly): dims are nested in cells nested in animals — this is the
MOST pseudoreplicated unit here (see repo STATE.md / earlier n-ladder). It is
included for comparison with the literature, NOT as our inferential test.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config  # noqa: E402

PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
CONTRASTS = [("expert", "naive"), ("expert", "intermediate"), ("intermediate", "naive")]


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "epoch_dims.csv"
    df = pd.read_csv(config.RESULTS_DIR / name)
    sig = df[df["sig"] == 1].copy()
    sig["peak_cc"] = pd.to_numeric(sig["peak_cc"], errors="coerce")
    print(f"{name}: {len(df)} dim-rows, {len(sig)} significant; "
          f"{df['animal'].nunique()} animals\n")
    print("Per-dimension held-out CC compared across epochs — rank-sum (Mann–Whitney U).")
    print("n_dim = pooled significant dims (pseudoreplicated); n_an = animals (honest unit).\n")
    for pair in PAIRS:
        g = sig[sig["pair"] == pair]
        if g.empty:
            continue
        print(f"  {pair}:")
        for a, b in CONTRASTS:
            va = g[g["epoch"] == a]["peak_cc"].dropna().to_numpy()
            vb = g[g["epoch"] == b]["peak_cc"].dropna().to_numpy()
            if va.size < 3 or vb.size < 3:
                continue
            U, p = stats.mannwhitneyu(va, vb, alternative="two-sided")
            n_an_a = g[g["epoch"] == a]["animal"].nunique()
            n_an_b = g[g["epoch"] == b]["animal"].nunique()
            star = " *" if p < 0.05 else ""
            print(f"    {a[:3]}({va.size}d/{n_an_a}an, cc={np.median(va):.3f}) vs "
                  f"{b[:3]}({vb.size}d/{n_an_b}an, cc={np.median(vb):.3f})  "
                  f"U-p={p:.3g}{star}")
    print("\nReminder: with dims-as-n the effective N is inflated ~5-15x over animals; "
          "p-values here are anti-conservative and not our inferential test.")


if __name__ == "__main__":
    main()
