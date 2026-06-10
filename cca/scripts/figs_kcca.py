"""Figures for the nonlinear (kernel) CCA section of the report.

From results/kcca_metrics{,_fsincl}.csv and kcca_transition{,_fsincl}.csv, renders
into the ResearchVault attachments. Every bar carries the **per-animal points**
overlaid as dots and a significance annotation, matching the linear figures
(scripts/figs_report.py house style):
  * kcca_vs_linear  -- per-pair held-out KCCA vs LINEAR CC (session, all animals),
                       grouped bars + animal dots, paired-Wilcoxon star (superiority test)
  * kcca_levels     -- per-pair KCCA CC1, IFI (+ sign test vs 0), Gini; bars + dots
  * kcca_epochs     -- naive/intermediate/expert KCCA CC1, IFI, Gini (learners),
                       grouped bars + dots, expert−naive Wilcoxon star
  * kcca_transition -- cued vs uncued deltas (CC, IFI, Gini), bars + dots,
                       signed-rank star vs 0
Session = unit of analysis; dots = animals; bars = mean ± SEM over animals.
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
from tom_cca import config, mixed_effects, paired_stats  # noqa: E402

ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
EPOCHS = ["naive", "intermediate", "expert"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _sem(v):
    """(mean, sem, n) over the finite values."""
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if v.size == 0:
        return np.nan, np.nan, 0
    return (float(v.mean()),
            float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else 0.0,
            v.size)


def _per_animal(df, pair, col):
    """One value per animal (mean over that animal's rows in the given frame)."""
    return [float(_num(s[col]).mean()) for _, s in df[df["pair"] == pair].groupby("animal")]


def _star_vs0(vals):
    """Two-sided signed-rank p that the per-animal values differ from 0."""
    v = [x for x in vals if np.isfinite(x)]
    return paired_stats.wilcoxon_signed(v)[3] if len(v) >= 3 else np.nan


def _bar_points_sem(ax, pairs, series, ylabel, title, star_p=None,
                    color="#8e44ad", hline0=True):
    """Single-series bars (mean±SEM) + jittered per-animal dots + optional star.

    series: dict pair -> list of per-animal values. star_p: dict pair -> p (a star
    is drawn above the bar when p<0.05)."""
    xs = np.arange(len(pairs))
    means = [_sem(series.get(p, []))[0] for p in pairs]
    sems = [_sem(series.get(p, []))[1] for p in pairs]
    ax.bar(xs, means, yerr=sems, capsize=3, color=color, alpha=0.55, zorder=1)
    for i, p in enumerate(pairs):
        v = [x for x in series.get(p, []) if np.isfinite(x)]
        jit = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.3
        ax.scatter(np.full(len(v), i) + jit, v, s=12, color="#2c3e50",
                   zorder=2, alpha=0.8)
        if star_p is not None and np.isfinite(star_p.get(p, np.nan)) and star_p[p] < 0.05:
            top = max((means[i] or 0) + (sems[i] or 0), max(v) if v else 0)
            ax.text(i, top, "*", ha="center", va="bottom", fontsize=13)
    if hline0:
        ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(pairs, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel); ax.set_title(title, fontsize=9)


def _grouped_bars_points(ax, pairs, groups, ylabel, title, contrast_star=None,
                         hline0=False):
    """Grouped bars (one cluster per pair) + jittered per-animal dots per group.

    groups: list of (label, colour, series_dict[pair -> per-animal list]).
    contrast_star: dict pair -> p; a star is drawn above the cluster when p<0.05."""
    xs = np.arange(len(pairs))
    ng = len(groups)
    w = 0.8 / ng
    tops = np.full(len(pairs), -np.inf)
    for j, (lab, color, series) in enumerate(groups):
        off = (j - (ng - 1) / 2) * w
        means = [_sem(series.get(p, []))[0] for p in pairs]
        sems = [_sem(series.get(p, []))[1] for p in pairs]
        ax.bar(xs + off, means, w, yerr=sems, capsize=2, label=lab, color=color,
               alpha=0.75, zorder=1)
        for i, p in enumerate(pairs):
            v = [x for x in series.get(p, []) if np.isfinite(x)]
            jit = (np.random.default_rng(i * 13 + j).random(len(v)) - 0.5) * (w * 0.7)
            ax.scatter(np.full(len(v), i) + off + jit, v, s=9, color="#2c3e50",
                       zorder=2, alpha=0.7)
            if v:
                tops[i] = max(tops[i], (means[i] or 0) + (sems[i] or 0), max(v))
    if contrast_star:
        for i, p in enumerate(pairs):
            pv = contrast_star.get(p, np.nan)
            if np.isfinite(pv) and pv < 0.05 and np.isfinite(tops[i]):
                ax.text(i, tops[i], "*", ha="center", va="bottom", fontsize=12)
    if hline0:
        ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(pairs, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel(ylabel); ax.set_title(title, fontsize=9); ax.legend(fontsize=7)


def fig_kcca_vs_linear(df, tag):
    sess = df[df["scope"] == "session"]
    kser = {p: _per_animal(sess, p, "cc_kcca") for p in PAIRS}
    lser = {p: _per_animal(sess, p, "cc_linear") for p in PAIRS}
    star = {}                                           # paired Wilcoxon KCCA−linear
    for p in PAIRS:
        d = [a - b for a, b in zip(kser[p], lser[p]) if np.isfinite(a) and np.isfinite(b)]
        star[p] = paired_stats.wilcoxon_signed(d)[3] if len(d) >= 3 else np.nan
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    _grouped_bars_points(
        ax, PAIRS,
        [("kernel CCA", "#8e44ad", kser), ("linear CCA", "#2980b9", lser)],
        "held-out CC$_1$ (session)",
        f"Nonlinear vs linear coupling — {tag}\n(bars = mean±SEM over animals; "
        "dots = animals; * = paired Wilcoxon KCCA−linear p<0.05)",
        contrast_star=star, hline0=True)
    fig.tight_layout(); fig.savefig(ATT / f"HCV1_KCCA_{tag}_vs_linear.png", dpi=150)
    plt.close(fig)


def fig_kcca_levels(df, tag):
    sess = df[df["scope"] == "session"]
    panels = [("cc_kcca", "held-out KCCA CC$_1$", False),
              ("ifi", "IFI (>0: X leads Y)", True),
              ("gini_x", "Gini (structure coeff, X)", False)]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, (col, lab, do_star) in zip(axes, panels):
        series = {p: _per_animal(sess, p, col) for p in PAIRS}
        star = {p: _star_vs0(series[p]) for p in PAIRS} if do_star else None
        _bar_points_sem(ax, PAIRS, series, lab,
                        lab + (" (* = p<0.05 vs 0)" if do_star else ""),
                        star_p=star, hline0=(col == "ifi"))
    fig.suptitle(f"Kernel-CCA levels by pair — {tag} (session; dots = animals, "
                 "bar = mean ± SEM over animals)", fontsize=11)
    fig.tight_layout(); fig.savefig(ATT / f"HCV1_KCCA_{tag}_levels.png", dpi=150)
    plt.close(fig)


def fig_kcca_epochs(df, tag):
    ep = df[df["scope"].isin(EPOCHS)]
    if ep.empty:
        return
    panels = [("cc_kcca", "held-out KCCA CC$_1$"), ("ifi", "IFI"),
              ("gini_x", "Gini (structure coeff)")]
    cols = {"naive": "#3498db", "intermediate": "#9b59b6", "expert": "#e74c3c"}
    fig, axes = plt.subplots(3, 1, figsize=(11, 12))   # stacked → legible at note width
    for ax, (col, lab) in zip(axes, panels):
        groups = [(e, cols[e], {p: _per_animal(ep[ep["scope"] == e], p, col)
                                for p in PAIRS}) for e in EPOCHS]
        star = {}                                       # expert−naive per-animal Wilcoxon
        for p in PAIRS:
            g = ep[ep["pair"] == p]
            if g["animal"].nunique() >= 3:
                recs = [{"animal_id": r["animal"], "epoch": r["scope"],
                         "value": float(_num(pd.Series([r[col]]))[0])}
                        for _, r in g.iterrows()]
                star[p] = mixed_effects.collapsed_epoch_contrasts(recs) \
                    .get("expert-naive", {}).get("p", np.nan)
            else:
                star[p] = np.nan
        _grouped_bars_points(ax, PAIRS, groups, lab,
                             lab + "  (* = expert−naive Wilcoxon p<0.05)",
                             contrast_star=star, hline0=(col == "ifi"))
    fig.suptitle(f"Kernel-CCA across learning epochs — {tag} (learners; dots = animals, "
                 "bar = mean ± SEM over animals)", fontsize=11)
    fig.tight_layout(); fig.savefig(ATT / f"HCV1_KCCA_{tag}_epochs.png", dpi=150)
    plt.close(fig)


def fig_kcca_transition(tag, suf):
    path = config.RESULTS_DIR / f"kcca_transition{suf}.csv"
    if not path.is_file():
        return
    t = pd.read_csv(path)
    learn = t[t["learner"] == 1]
    panels = [("d_cc_kcca", "Δ KCCA CC$_1$"), ("d_ifi", "Δ IFI"), ("d_gini_x", "Δ Gini")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, (col, lab) in zip(axes, panels):
        series, star = {}, {}                           # one delta per animal per pair
        for p in PAIRS:
            vals = _num(learn[learn["pair"] == p][col]).dropna().tolist()
            series[p] = vals; star[p] = _star_vs0(vals)
        _bar_points_sem(ax, PAIRS, series, f"{lab} (cued − uncued)",
                        f"{lab}: uncued→cued (* = signed-rank p<0.05 vs 0)",
                        star_p=star, color="#16a085", hline0=True)
    fig.suptitle(f"Kernel-CCA uncued→cued — {tag} (learners; dots = animals, "
                 "bar = mean ± SEM over animals)", fontsize=11)
    fig.tight_layout(); fig.savefig(ATT / f"HCV1_KCCA_{tag}_transition.png", dpi=150)
    plt.close(fig)


def main():
    ATT.mkdir(parents=True, exist_ok=True)
    for tag in ("fsexcl", "fsincl"):
        suf = "" if tag == "fsexcl" else "_fsincl"
        path = config.RESULTS_DIR / f"kcca_metrics{suf}.csv"
        if not path.is_file():
            print(f"  (missing {path.name} — skip {tag})"); continue
        df = pd.read_csv(path)
        fig_kcca_vs_linear(df, tag)
        fig_kcca_levels(df, tag)
        fig_kcca_epochs(df, tag)
        fig_kcca_transition(tag, suf)
        print(f"  {tag}: KCCA figures written ({len(df)} rows)")
    print("figures ->", ATT)


if __name__ == "__main__":
    main()
