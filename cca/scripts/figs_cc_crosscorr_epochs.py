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
          Tom's layout (2026-08-18): per pair TWO panels — naive vs expert for the
          FF-labelled CCs (first-named area leads at ±50 ms) and, next to it, naive vs
          expert for the FB-labelled CCs (second-named area leads). Labels from the
          whole-session fit (identical to slide 10 of CCA_update_20260807.pptx).
          Signed curves of opposite-labelled CCs are never averaged together. Titles
          carry item 2's per-label expert − naive IFI test (cc_label_track_stats_*).

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
import numpy as np
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


def _fmt(row, name, row_cc=None, with_n=False):
    """'name Δ p(animals) | CCs p' — the animals-as-n Δ and p, then the CCs-as-n p as
    the power check (n = significant CCs pooled over animals; printed once per panel)."""
    if row is None:
        return f"{name} —"
    star = "*" if row["p"] < 0.05 else ""
    s = f"{name} Δ{row['delta']:+.3f} p={row['p']:.2g}{star}"
    if row_cc is not None:
        s += f" | CCs p={row_cc['p']:.2g}{'*' if row_cc['p'] < 0.05 else ''}"
        if with_n:
            s += f" ({int(row_cc['n'])} CCs)"
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
        t1 = _fmt(_p(stats, pair, "ifi", c), "IFI@±50",
                  _p(stats, pair, "ifi", c, "ccs"), with_n=True)
        t2 = _fmt(_p(stats, pair, "peak_minus_far", c), "peak−base",
                  _p(stats, pair, "peak_minus_far", c, "ccs"))
        t3 = _fmt(_p(stats, pair, "far_r", c), "baseline",
                  _p(stats, pair, "far_r", c, "ccs"))
        ax.set_title(f"{pair} — expert − naive\n{t1}\n{t2}\n{t3}", fontsize=7.5)
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


def fig_fffb(curves: pd.DataFrame, fs: str, label_stats: pd.DataFrame | None = None):
    """Tom's 2026-08-18 layout: per pair, TWO panels side by side — the naive vs expert
    cross-correlograms of the FF-labelled CCs (first-named area leads at ±50 ms) in one,
    of the FB-labelled CCs (second-named area leads) in the other. 8 pairs x 2 = 16
    panels, two pairs per row. Titles carry item 2's per-label expert − naive paired t
    on the IFI at ±50 ms (`cc_label_track_stats_*`, animals-as-n).
    """
    fig, axes = figstyle.grid(2 * len(PAIRS), ncols=4)
    ls = label_stats.set_index("pair") if label_stats is not None else None
    for k, ax in enumerate(axes):
        pair, group = PAIRS[k // 2], ("FF", "FB")[k % 2]
        a1, a2 = pair.split("-")
        leader = a1 if group == "FF" else a2
        drawn = 0
        for epoch in ("naive", "expert"):
            drawn += _draw(ax, curves, pair, epoch, group, "-", 2.0, 1.0, epoch)
        if drawn == 0:
            ax.set_title(f"{pair} — {group} CCs\n(no data)", fontsize=9)
            ax.axis("off"); continue
        ax.axvline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axhline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axvspan(-HEADLINE_MS, HEADLINE_MS, color="#000", alpha=0.03, zorder=0)
        title = f"{pair} — {group}-labelled CCs ({leader} leads at ±{HEADLINE_MS} ms)"
        if ls is not None and pair in ls.index:
            r = ls.loc[pair]
            n, d, pv = r.get(f"n_{group}_ifi"), r.get(f"d_{group}_ifi"), r.get(f"p_{group}_ifi")
            if np.isfinite(pv):
                title += (f"\nexpert − naive IFI@±50: Δ{d:+.3f} p={pv:.2g}"
                          f"{'*' if pv < 0.05 else ''} (n={int(n)} animals)")
        ax.set_title(title, fontsize=8)
        ax.set_xlabel(f"lag (ms)   +: {a1} leads", fontsize=8)
        ax.set_ylabel(f"mean r over {group} CCs (frozen, in-sample)", fontsize=8)
        ax.legend(fontsize=7, loc="upper right", framealpha=0.6)
    fig.suptitle(f"Cross-correlograms, naive vs expert, ONE DIRECTION PER PANEL — {fs}\n"
                 f"per pair: left = FF-labelled CCs (+IFI, first-named area leads), right = "
                 f"FB-labelled CCs (−IFI, second-named area leads); label from the "
                 f"whole-session fit at ±{HEADLINE_MS} ms (as on slide 10) — the asymmetry "
                 f"inside the shaded band is circular by construction\n"
                 f"per-animal-first mean ± SEM (bands only for n ≥ 3); frozen axes; titles: "
                 f"item-2 paired t of the per-label IFI, animals-as-n; the naive-vs-expert "
                 f"HEIGHT caveat of the all-CC figure applies (baseline offset, not learning)",
                 fontsize=10)
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
    ls_csv = RES / f"cc_label_track_stats_bin10{suf}.csv"
    label_stats = pd.read_csv(ls_csv) if ls_csv.exists() else None
    fig_fffb(curves, fs, label_stats)


if __name__ == "__main__":
    main()
