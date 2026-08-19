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

Figure 3/4  HCV1_cc_crosscorr_epochs_FF_<fs>_bin10, HCV1_cc_crosscorr_epochs_FB_<fs>_bin10
          The same panels, one direction per FIGURE (8 panels each) — Theo, 2026-08-18.

Frozen axes, so the SAME components are compared across epochs. r is IN-SAMPLE for
the whole session: a contrast statistic, not a coupling strength.
Positive lag = the FIRST-named area leads.

Usage: PYTHONPATH=src python scripts/figs_cc_crosscorr_epochs.py [fsincl] [--label-col label|label_int|label_xv|label_loo]
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
from analyze_cc_crosscorr_epochs import (  # noqa: E402
    EPOCHS, LABEL_DESC, LABEL_TAGS, PAIRS, cross_animal, _label_col_from_argv)

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


def _label_panel(ax, curves, pair, group, ls, label_col="label"):
    """One panel: naive vs expert cross-correlograms of the ``group``-labelled CCs of
    ``pair``; title carries item 2's per-label expert − naive IFI test. Returns the
    number of curves drawn (0 = no data, panel switched off)."""
    a1, a2 = pair.split("-")
    leader = a1 if group == "FF" else a2
    drawn = 0
    for epoch in ("naive", "expert"):
        drawn += _draw(ax, curves, pair, epoch, group, "-", 2.0, 1.0, epoch)
    if drawn == 0:
        ax.set_title(f"{pair} — {group} CCs\n(no data)", fontsize=9)
        ax.axis("off"); return 0
    ax.axvline(0, color="k", lw=0.7, ls=":", zorder=0)
    ax.axhline(0, color="k", lw=0.7, ls=":", zorder=0)
    ax.axvspan(-HEADLINE_MS, HEADLINE_MS, color="#000", alpha=0.03, zorder=0)
    lab_desc = {"label": "whole-session label", "label_loo": "leave-epoch-out label",
                "label_int": "label from INTERMEDIATE trials only",
                "label_xv": "label from all trials OUTSIDE naive/expert"}[label_col]
    title = (f"{pair} — {group} CCs ({leader} leads at ±{HEADLINE_MS} ms)\n"
             f"{lab_desc}")
    if ls is not None:
        r = ls[(ls["pair"] == pair) & (ls["group"] == group) & (ls["metric"] == "ifi") &
               (ls["contrast"] == "expert-naive")]
        if not r.empty and np.isfinite(r.iloc[0]["p"]):
            r = r.iloc[0]
            title += (f"\nexpert − naive IFI@±50: Δ{r['delta']:+.3f} p={r['p']:.2g}"
                      f"{'*' if r['p'] < 0.05 else ''} (n={int(r['n'])} animals)")
    ax.set_title(title, fontsize=8)
    ax.set_xlabel(f"lag (ms)   +: {a1} leads", fontsize=8)
    ax.set_ylabel(f"mean r over {group} CCs (frozen, in-sample)", fontsize=8)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.6)
    return drawn


def fig_fffb(curves: pd.DataFrame, fs: str, label_stats: pd.DataFrame | None = None,
             label_col: str = "label"):
    """Tom's 2026-08-18 layout: per pair, TWO panels side by side — the naive vs expert
    cross-correlograms of the FF-labelled CCs (first-named area leads at ±50 ms) in one,
    of the FB-labelled CCs (second-named area leads) in the other. 8 pairs x 2 = 16
    panels, two pairs per row. :func:`fig_label` gives each direction its own figure.
    """
    fig, axes = figstyle.grid(2 * len(PAIRS), ncols=4)
    ls = _per_label_rows(label_stats)
    for k, ax in enumerate(axes):
        _label_panel(ax, curves, PAIRS[k // 2], ("FF", "FB")[k % 2], ls, label_col)
    tag = LABEL_TAGS[label_col]
    fig.suptitle(f"Cross-correlograms, naive vs expert, ONE DIRECTION PER PANEL — {fs}\n"
                 f"per pair: left = FF-labelled CCs (+IFI, first-named area leads), right = "
                 f"FB-labelled CCs (−IFI, second-named area leads)\n"
                 + _label_caveat(label_col), fontsize=10)
    figstyle.save(fig, ATT / f"HCV1_cc_crosscorr_epochs_fffb{tag}_{fs}_bin10.png")
    print(f"wrote HCV1_cc_crosscorr_epochs_fffb{tag}_{fs}_bin10.png")


def fig_label(curves: pd.DataFrame, fs: str, group: str,
              label_stats: pd.DataFrame | None = None, label_col: str = "label"):
    """One direction per FIGURE: 8 panels (one per pair), naive vs expert of the
    ``group``-labelled CCs only. HCV1_cc_crosscorr_epochs_{FF,FB}<tag>_<fs>_bin10,
    tag = "" (whole-session label) / _labint / _labxv / _loo."""
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    ls = _per_label_rows(label_stats)
    for ax, pair in zip(axes, PAIRS):
        _label_panel(ax, curves, pair, group, ls, label_col)
    what = ("FF-labelled CCs (+IFI: the FIRST-named area leads)" if group == "FF"
            else "FB-labelled CCs (−IFI: the SECOND-named area leads)")
    tag = LABEL_TAGS[label_col]
    fig.suptitle(f"Cross-correlograms, naive vs expert — {what} — {fs}\n"
                 + _label_caveat(label_col), fontsize=10)
    figstyle.save(fig, ATT / f"HCV1_cc_crosscorr_epochs_{group}{tag}_{fs}_bin10.png")
    print(f"wrote HCV1_cc_crosscorr_epochs_{group}{tag}_{fs}_bin10.png")


def _per_label_rows(stats):
    """The FF/FB rows of the tagged stats file (group != 'all'), or None."""
    if stats is None or "group" not in stats.columns:
        return None
    sub = stats[stats["group"].isin(["FF", "FB"])]
    return sub if not sub.empty else None


def _label_caveat(label_col: str) -> str:
    xv = label_col in ("label_int", "label_xv")
    circ = ("the label never saw the plotted trials, so the asymmetry inside the shaded "
            "band is NOT circular" if xv else
            "the asymmetry inside the shaded band is circular by construction")
    return (f"label = {LABEL_DESC[label_col]} — {circ}\n"
            f"per-animal-first mean ± SEM (bands only for n ≥ 3); frozen axes (weights "
            f"from ALL trials — only the FF/FB assignment is held out); titles: paired t "
            f"of the per-label IFI, animals-as-n; the naive-vs-expert HEIGHT caveat of the "
            f"all-CC figure applies (baseline offset, not learning)")


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    label_col = _label_col_from_argv(sys.argv[1:])
    tag = LABEL_TAGS[label_col]
    src = RES / f"cc_crosscorr_epochs{tag}_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/analyze_cc_crosscorr_epochs.py "
                 f"--label-col {label_col} first")
    curves = pd.read_csv(src)
    stats = pd.read_csv(RES / f"cc_crosscorr_epochs_stats{tag}_bin10{suf}.csv")
    if label_col == "label":              # the all-CC figure is label-independent
        fig_all(curves, stats, fs)
    fig_fffb(curves, fs, stats, label_col)
    fig_label(curves, fs, "FF", stats, label_col)
    fig_label(curves, fs, "FB", stats, label_col)


if __name__ == "__main__":
    main()
