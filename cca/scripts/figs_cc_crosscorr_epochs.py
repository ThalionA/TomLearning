"""2026-08-07 ask 2 figures — cross-correlograms of naive and expert, whole and split
by positive/negative-grouped (FF/FB) CCs.

Reads results/cc_crosscorr_epochs_bin10{,_fsincl}.csv (analyze_cc_crosscorr_epochs.py)
and its *_stats_* file for the panel titles.

Figure 1  HCV1_cc_crosscorr_epochs_<fs>_bin10
          Per pair, mean r(lag) over ALL of an animal's significant CCs, then across
          animals ± SEM: naive (blue) vs expert (red) bold, intermediate thin grey.
          Title carries the expert − naive paired t on the all-CC IFI at ±50 ms and on
          peak r, PLUS intermediate − naive on peak r — the built-in diagnostic that
          the curve-height difference is a first-trials/fit-set effect (intermediate,
          pre-LP, already sits with expert), not learning.

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


def _p(stats: pd.DataFrame, pair: str, metric: str, contrast: str):
    r = stats[(stats["pair"] == pair) & (stats["metric"] == metric) &
              (stats["contrast"] == contrast)]
    return None if r.empty else r.iloc[0]


def _fmt(row, name):
    if row is None:
        return f"{name} —"
    star = "*" if row["p"] < 0.05 else ""
    return f"{name} Δ{row['delta']:+.3f} p={row['p']:.2g}{star}"


def _draw(ax, curves, pair, epoch, group, ls, lw, alpha, label):
    lags, mean, sem, n = cross_animal(curves, pair, epoch, group)
    if n == 0:
        return 0
    c = EPOCH_COLOUR[epoch]
    ax.plot(lags, mean, ls, color=c, lw=lw, alpha=alpha, zorder=3,
            label=f"{label} (n={n})")
    if lw >= 1.6:
        ax.fill_between(lags, mean - sem, mean + sem, color=c, alpha=0.14, zorder=2)
    return n


def fig_all(curves: pd.DataFrame, stats: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        x_lead = pair.split("-")[0]
        drawn = 0
        for epoch in EPOCHS:
            lw, alpha = (1.1, 0.9) if epoch == "intermediate" else (2.0, 1.0)
            drawn += _draw(ax, curves, pair, epoch, "all", "-", lw, alpha, epoch)
        if drawn == 0:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        ax.axvline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axhline(0, color="k", lw=0.7, ls=":", zorder=0)
        ax.axvspan(-HEADLINE_MS, HEADLINE_MS, color="#000", alpha=0.03, zorder=0)
        t1 = _fmt(_p(stats, pair, "ifi", "expert-naive"), "IFI@±50 exp−naive")
        t2 = _fmt(_p(stats, pair, "peak_r", "expert-naive"), "peak r exp−naive")
        t3 = _fmt(_p(stats, pair, "peak_r", "intermediate-naive"), "peak r int−naive")
        ax.set_title(f"{pair}\n{t1}\n{t2} · {t3}", fontsize=8)
        ax.set_xlabel(f"lag (ms)   +: {x_lead} leads", fontsize=8)
        ax.set_ylabel("mean r over sig. CCs (frozen, in-sample)", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper right", framealpha=0.6)
    fig.suptitle(f"Cross-correlograms by learning epoch, all significant CCs — {fs} | "
                 f"frozen axes (same components every epoch), per-animal-first mean ± SEM\n"
                 f"⚠ curve HEIGHT: naive is uniquely low on this all-trials fit — "
                 f"intermediate (pre-LP) already equals expert, and the balanced-trial "
                 f"fit shows no rise — read height as 'first trials differ', not learning; "
                 f"IFI (shape) is the naive-vs-expert statistic and is null",
                 fontsize=10)
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
                 f"solid = +IFI-labelled (FF, {'first'} area leads), dashed = −IFI-labelled "
                 f"(FB); label from the whole-session fit at ±{HEADLINE_MS} ms (as on "
                 f"slide 10), so the asymmetry inside the shaded band is circular by "
                 f"construction\nper-animal-first mean ± SEM; frozen axes; the naive-"
                 f"vs-expert height caveat of the all-CC figure applies here too",
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
    fig_fffb(curves, fs)


if __name__ == "__main__":
    main()
