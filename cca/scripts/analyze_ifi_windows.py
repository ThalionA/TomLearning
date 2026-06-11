"""Which lag-integration window gives the cleanest directionality?

Reads results/ifi_windows_bin{25,50}{,_fsincl}.csv. For every pair and every window
|lag| <= w, tests the across-animal IFI against 0 (paired t + signed-rank) and
reports the window that maximises the *reproducibility* of the directionality —
|mean / SEM| (the across-animal t-statistic). A narrow window misses an off-zero
peak; a wide one dilutes it with far-lag noise; the cleanest window is in between.

Run with no args for FS-excluded 25 ms; pass e.g. ``bin50`` / ``fsincl`` to switch.
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


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else ""
    binms = "50" if "bin50" in arg or "50" in arg else "25"
    suffix = "_fsincl" if "fsincl" in arg else ""
    path = config.RESULTS_DIR / f"ifi_windows_bin{binms}{suffix}.csv"
    if not path.is_file():
        print(f"missing {path.name} — run run_ifi_windows.py --bin-ms {binms} first")
        return
    df = pd.read_csv(path)
    bin_ms = int(df["bin_ms"].iloc[0])
    print(f"{path.name}: {df['animal'].nunique()} animals, "
          f"windows w=1..{int(df['window_bins'].max())} ({bin_ms} ms/bin)\n")
    print("=" * 96)
    print("Per pair: across-animal IFI by lag-integration window — | mean/SEM | (t-stat) "
          "picks the cleanest")
    print("=" * 96)

    best_per_pair = {}
    for pair in PAIRS:
        g = df[df["pair"] == pair]
        if g.empty:
            continue
        cells, scored = [], []
        for w in sorted(g["window_bins"].unique()):
            v = _num(g[g["window_bins"] == w]["ifi"]).dropna().to_numpy()
            if v.size < 3:
                continue
            sem = v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else np.nan
            tstat = v.mean() / sem if sem and np.isfinite(sem) and sem > 0 else 0.0
            tp = stats.ttest_1samp(v, 0).pvalue if v.size > 1 else np.nan
            scored.append((w, abs(tstat), v.mean(), tp, v.size))
            cells.append(f"w{w}({w*bin_ms}ms):{v.mean():+.3f}|t={tstat:+.1f}")
        if not scored:
            continue
        best = max(scored, key=lambda r: r[1])           # max |mean/SEM|
        best_per_pair[pair] = best
        print(f"\n  {pair:8s} (n={best[4]}) — cleanest = w{best[0]} "
              f"({best[0]*bin_ms} ms): IFI={best[2]:+.3f}, t={best[1]:.2f}, "
              f"t-p={best[3]:.2g}")
        print("     " + "  ".join(cells))

    print("\n" + "=" * 96)
    print("SUMMARY — cleanest integration window per pair (max across-animal t-stat)")
    print("=" * 96)
    for pair, (w, t, m, tp, n) in best_per_pair.items():
        star = "*" if (np.isfinite(tp) and tp < 0.05) else " "
        print(f"  {pair:8s} w={w:>2} ({w*bin_ms:>3} ms)  IFI={m:+.3f}  t={t:.2f}  p={tp:.2g}{star}")
    # global recommendation: the window most often best / best on the strong directional pairs
    if best_per_pair:
        ws = [b[0] for b in best_per_pair.values()]
        print(f"\n  median cleanest window across pairs: w={int(np.median(ws))} "
              f"({int(np.median(ws))*bin_ms} ms)")


if __name__ == "__main__":
    main()
