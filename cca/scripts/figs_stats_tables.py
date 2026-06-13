"""Stats-summary table FIGURES for the deck (the 'summary-table slides' view).

Renders clean colour-coded tables (green = significant) as PNGs into the vault
attachments, so the deck globs them like any figure. Two tables:
  * learning-slopes — per pair, d(Gini)/d(trial_frac) and d(IFI)/d(trial_frac),
    LEARNERS vs POOLED (all animals) side by side, signed-rank on per-animal slopes;
  * directionality  — per pair, the held-out IFI window-sweep result (mean IFI at the
    cleanest window, t, p) — the §3.2 directionality summary.

Run:  PYTHONPATH=src python scripts/figs_stats_tables.py [bin25|bin50]
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402
from figs_report import PAIRS, _slope_quant, _stars, _cohort  # noqa: E402
from tom_cca import config  # noqa: E402

ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"
SIG, NS = "#1a7f4b", "#7a7972"


def _cell(slope, p, n):
    star, _ = _stars(p)
    col = SIG if (np.isfinite(p) and p < 0.05) else NS
    if not np.isfinite(slope):
        return "—", NS
    arrow = " ↓" if slope < 0 else " ↑"
    return f"{slope:+.3f}{arrow}\n{star} · n={n}", col


def _table_fig(title, col_heads, rows, fname, note):
    """rows: list of (rowlabel, [(text, colour), ...])."""
    nr, nc = len(rows), len(col_heads)
    fig, ax = plt.subplots(figsize=(2.2 + 2.5 * nc, 1.1 + 0.62 * nr))
    ax.axis("off")
    ax.set_xlim(0, nc + 1); ax.set_ylim(0, nr + 1.6)
    ax.text((nc + 1) / 2, nr + 1.25, title, ha="center", va="center",
            fontsize=13, fontweight="bold", color="#23231f")
    for j, h in enumerate(col_heads):                       # header row
        ax.text(j + 1.5, nr + 0.5, h, ha="center", va="center",
                fontsize=10.5, fontweight="bold", color="#23231f")
    ax.plot([0.2, nc + 1], [nr + 0.1, nr + 0.1], color="#cfcec6", lw=1)
    for i, (rlab, cells) in enumerate(rows):
        y = nr - i - 0.5
        ax.text(0.7, y, rlab, ha="left", va="center", fontsize=10.5,
                fontweight="bold", color="#23231f")
        for j, (txt, col) in enumerate(cells):
            ax.text(j + 1.5, y, txt, ha="center", va="center", fontsize=9.5,
                    color=col, fontweight="bold", linespacing=1.3)
        if i % 2:                                           # zebra
            ax.axhspan(y - 0.5, y + 0.5, color="#f5f5f0", zorder=0)
    ax.text(0.7, -0.55, note, ha="left", va="center", fontsize=8.5, color="#5f5e5a")
    figstyle.save(fig, ATT / fname)
    print(f"wrote {fname}")


def learning_slopes(stem, tag):
    df = pd.read_csv(config.RESULTS_DIR / f"{stem}.csv")
    heads = ["Gini slope\n(learners)", "Gini slope\n(pooled, all 16)",
             "IFI slope\n(learners)", "IFI slope\n(pooled, all 16)"]
    rows = []
    for p in PAIRS:
        cells = []
        for metric in ("gini_x", "ifi"):
            for pool in (False, True):
                s, pv, n = _slope_quant(_cohort(df, pool), metric, "trial_frac", p)
                cells.append(_cell(s, pv, n))
        rows.append((p, cells))
    binms = "10" if "10" in stem else ("50" if "50" in stem else "25")
    _table_fig(f"Learning-axis slopes vs trial-fraction — {tag} ({binms} ms)",
               heads, rows, f"HCV1_stats_slopes_{tag}_bin{binms}.png",
               "signed-rank on per-animal slopes vs 0 · green = p<0.05 "
               "(*** p<.001  ** p<.01  * p<.05) · ↓/↑ = sign · pooled adds the 4 non-learners")


def directionality(stem, tag):
    """Held-out IFI window-sweep summary (cleanest window per pair)."""
    binms = "10" if "10" in stem else ("50" if "50" in stem else "25")
    p = config.RESULTS_DIR / f"ifi_windows_bin{binms}{'_fsincl' if tag == 'fsincl' else ''}.csv"
    if not p.is_file():
        print(f"  (missing {p.name} — skip directionality table)"); return
    d = pd.read_csv(p)
    from scipy import stats as st
    heads = ["cleanest\nwindow", "mean IFI\n(>0: 1st leads)", "t", "p\n(t vs 0, across animals)"]
    rows = []
    for pair in PAIRS:
        g = d[d["pair"] == pair]
        if g.empty:
            continue
        bt, bp, bw, bm = 0, np.nan, None, np.nan
        for w, gw in g.groupby("window_ms"):
            v = pd.to_numeric(gw["ifi"], errors="coerce").dropna().to_numpy()
            if v.size < 3 or np.std(v) == 0:
                continue
            t, pv = st.ttest_1samp(v, 0)
            if abs(t) > abs(bt):
                bt, bp, bw, bm = t, pv, int(w), float(np.mean(v))
        if bw is None:
            continue
        col = SIG if (np.isfinite(bp) and bp < 0.05) else NS
        star, _ = _stars(bp)
        rows.append((pair, [(f"±{bw} ms", "#23231f"),
                            (f"{bm:+.3f}", col), (f"{bt:.2f}", col), (star, col)]))
    _table_fig(f"Directionality — held-out IFI window sweep — {tag} ({binms} ms)",
               heads, rows, f"HCV1_stats_directionality_{tag}_bin{binms}.png",
               "per pair, the lag-integration window maximising |t| (held-out, segment-aware) · "
               "green = p<0.05 · CA1→RSC is the a-priori, robust flow")


def main():
    stem_excl = sys.argv[1] if len(sys.argv) > 1 else "trajectory_w15_bin25"
    learning_slopes(stem_excl, "fsexcl")
    fsincl = stem_excl + "_fsincl"
    if (config.RESULTS_DIR / f"{fsincl}.csv").is_file():
        learning_slopes(fsincl, "fsincl")
    directionality(stem_excl, "fsexcl")


if __name__ == "__main__":
    main()
