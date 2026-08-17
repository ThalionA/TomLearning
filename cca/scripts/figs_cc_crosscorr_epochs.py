"""2026-08-07 ask 2 figures — cross-correlograms of naive and expert, whole and split
by positive/negative-grouped (FF/FB) CCs.

Reads results/cc_crosscorr_epochs_bin10{,_fsincl}.csv (analyze_cc_crosscorr_epochs.py)
and its *_stats_* file for the panel titles.

Figure 1  HCV1_cc_crosscorr_epochs_<fs>_bin10
          Per pair, mean r(lag) over ALL of an animal's significant CCs, then across
          animals ± SEM: naive (blue) vs expert (red) bold with SEM bands,
          intermediate thin grey (mean only, no band). Title carries the expert −
          naive paired t on (i) the all-CC IFI at ±50 ms — the naive-vs-expert
          statistic — and (ii) the two halves of the curve-height difference: peak
          MINUS the |lag| ≥ 200 ms baseline (coupling-specific) and the baseline
          itself (a lag-independent offset = slow co-modulation, speed the obvious
          candidate). SEM bands are drawn only when n ≥ 3 animals.

Figure 2  HCV1_cc_crosscorr_epochs_fffb_<fs>_bin10
          The same, split by the whole-session FF/FB label (identical to slide 10 of
          CCA_update_20260807.pptx): FF solid, FB dashed, naive vs expert only.
          Signed curves of opposite-labelled CCs are never averaged together.

Frozen axes, so the SAME components are compared across epochs. r is IN-SAMPLE for
the whole session: a contrast statistic, not a coupling strength.
Positive lag = the FIRST-named area leads.

Usage: PYTHONPATH=src python scripts/figs_cc_crosscorr_epochs.py [fsincl]
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402
from analyze_cc_crosscorr_epochs import EPOCHS, PAIRS, cross_animal  # noqa: E402

figstyle.apply()
ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"
RES = Path(__file__).resolve().parents[1] / "results"
EPOCH_COLOUR = {"naive": "#2c6fbb", "intermediate": "#95a5a6", "expert": "#c0392b"}
HEADLINE_MS = 50


def _p(stats: pd.DataFrame, pair: str, metric: str, contrast: str, unit: str = "animals"):
    r = stats[(stats["pair"] == pair) & (stats["metric"] == metric) &
              (stats["contrast"] == contrast)]
    if "unit" in r.columns:
        r = r[r["unit"] == unit]
    return None if r.empty else r.iloc[0]


def _fmt(row, name, row_cc=None):
    """'name Δ p(animals) | p(CCs, n)' — the animals-as-n Δ and p, then the CCs-as-n
    p as the power check (n = significant CCs pooled over animals)."""
    if row is None:
        return f"{name} —"
    star = "*" if row["p"] < 0.05 else ""
    s = f"{name} Δ{row['delta']:+.3f} p={row['p']:.2g}{star}"
    if row_cc is not None:
        s += f" | CCs p={row_cc['p']:.2g}{'*' if row_cc['p'] < 0.05 else ''} (n={int(row_cc['n'])})"
    return s


def _draw(ax, curves, pair, epoch, group, ls, lw, alpha, label, band=True):
    lags, mean, sem, n = cross_animal(curves, pair, epoch, group)
    if n == 0:
        return 0
    c = EPOCH_COLOUR[epoch]
    ax.plot(lags, mean, ls, color=c, lw=lw, alpha=alpha, zorder=3,
            label=f"{label} (n={n})" + ("" if band else ", mean only"))
    if band and n >= 3:                    # an SEM from 2 animals is not a band
        ax.fill_between(lags, mean - sem, mean + sem, color=c, alpha=0.14, zorder=2)
    return n


def fig_all(curves: pd.DataFrame, stats: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        x_lead = pair.split("-")[0]
        drawn = 0
        for epoch in EPOCHS:
            mid = epoch == "intermediate"
            drawn += _draw(ax, curves, pair, epoch, "all", "-",
                           1.1 if mid else 2.0, 0.9 if mid else 1.0, epoch, band=not mid)
        if drawn == 0:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        ax.axvline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axhline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axvspan(-HEADLINE_MS, HEADLINE_MS, color="#000", alpha=0.03, zorder=0)
        c = "expert-naive"
        t1 = _fmt(_p(stats, pair, "ifi", c), "IFI@±50 exp−naive",
                  _p(stats, pair, "ifi", c, "ccs"))
        t2 = _fmt(_p(stats, pair, "peak_minus_far", c), "peak−baseline",
                  _p(stats, pair, "peak_minus_far", c, "ccs"))
        t3 = _fmt(_p(stats, pair, "far_r", c), "baseline offset",
                  _p(stats, pair, "far_r", c, "ccs"))
        ax.set_title(f"{pair}\n{t1}\nexp−naive: {t2}\n{t3}", fontsize=7.5)
        ax.set_xlabel(f"lag (ms)   +: {x_lead} leads", fontsize=8)
        ax.set_ylabel("mean r over sig. CCs (frozen, in-sample)", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper right", framealpha=0.6)
    ifi_an = stats[(stats["metric"] == "ifi") & (stats["contrast"] == "expert-naive")]
    if "unit" in ifi_an.columns:
        ifi_an = ifi_an[ifi_an["unit"] == "animals"]
    n_ifi_hits = int((ifi_an["p"] < 0.05).sum())
    fig.suptitle(f"Cross-correlograms by learning epoch, all significant CCs — {fs} | "
                 f"frozen axes (same components every epoch), per-animal-first mean ± SEM"
                 f" (bands only for n ≥ 3)\n"
                 f"IFI@±50 ms expert−naive: {n_ifi_hits}/{len(ifi_an)} pairs p<0.05 "
                 f"(animals-as-n) — the naive-vs-expert statistic. Curve HEIGHT is not: "
                 f"naive sits low at EVERY lag on this all-trials fit\n(a baseline offset "
                 f"= slow co-modulation, speed the obvious candidate; LP-independent; "
                 f"held-out refits show no strength change) — titles: Δ and p with "
                 f"animals as n, then p with significant CCs as n (power check; CCs are "
                 f"nested in animals)", fontsize=9.5)
    figstyle.save(fig, ATT / f"HCV1_cc_crosscorr_epochs_{fs}_bin10.png")
    print(f"wrote HCV1_cc_crosscorr_epochs_{fs}_bin10.png")


def fig_fffb(curves: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        x_lead = pair.split("-")[0]
        drawn = 0
        for group, ls in [("FF", "-"), ("FB", "--")]:
            for epoch in ("naive", "expert"):
                drawn += _draw(ax, curves, pair, epoch, group, ls, 2.0, 1.0,
                               f"{group} {epoch}")
        if drawn == 0:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        ax.axvline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axhline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axvspan(-HEADLINE_MS, HEADLINE_MS, color="#000", alpha=0.03, zorder=0)
        ax.set_title(pair, fontsize=10)
        ax.set_xlabel(f"lag (ms)   +: {x_lead} leads", fontsize=8)
        ax.set_ylabel("mean r over sig. CCs (frozen, in-sample)", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper right", framealpha=0.6, ncol=2)
    fig.suptitle(f"Cross-correlograms, naive vs expert, split by FF/FB label — {fs} | "
                 f"solid = +IFI-labelled (FF, first area leads), dashed = −IFI-labelled "
                 f"(FB); label from the whole-session fit at ±{HEADLINE_MS} ms (as on "
                 f"slide 10), so the asymmetry inside the shaded band is circular by "
                 f"construction\nper-animal-first mean ± SEM (bands only for n ≥ 3); "
                 f"frozen axes; the naive-vs-expert height caveat of the all-CC figure "
                 f"applies here too", fontsize=10)
    figstyle.save(fig, ATT / f"HCV1_cc_crosscorr_epochs_fffb_{fs}_bin10.png")
    print(f"wrote HCV1_cc_crosscorr_epochs_fffb_{fs}_bin10.png")


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    src = RES / f"cc_crosscorr_epochs_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/analyze_cc_crosscorr_epochs.py first")
    curves = pd.read_csv(src)
    stats = pd.read_csv(RES / f"cc_crosscorr_epochs_stats_bin10{suf}.csv")
    fig_all(curves, stats, fs)
    fig_fffb(curves, fs)


if __name__ == "__main__":
    main()
