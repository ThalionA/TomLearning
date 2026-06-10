"""Figures for the nonlinear (kernel) CCA section of the report.

From results/kcca_metrics{,_fsincl}.csv and kcca_transition{,_fsincl}.csv, renders
into the ResearchVault attachments:
  * kcca_vs_linear  -- per-pair held-out KCCA vs LINEAR CC (session, all animals),
                       points + mean±SEM, paired-Wilcoxon star (the superiority test)
  * kcca_levels     -- per-pair KCCA CC1, significance rate, IFI (direction), Gini
  * kcca_epochs     -- naive/intermediate/expert KCCA CC1, IFI, Gini (learners)
  * kcca_transition -- cued vs uncued deltas (CC, IFI, Gini)
Every panel: labelled axes, title, legend; session = unit of analysis.
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
from tom_cca import config, paired_stats  # noqa: E402

ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
EPOCHS = ["naive", "intermediate", "expert"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _sem(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if v.size == 0:
        return np.nan, np.nan
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else 0.0


def _per_animal(df, pair, col):
    return [float(_num(s[col]).mean()) for _, s in df[df["pair"] == pair].groupby("animal")]


def fig_kcca_vs_linear(df, tag):
    sess = df[df["scope"] == "session"]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    xs = np.arange(len(PAIRS)); w = 0.38
    for off, col, lab, c in [(-w / 2, "cc_kcca", "kernel CCA", "#8e44ad"),
                             (+w / 2, "cc_linear", "linear CCA", "#2980b9")]:
        means, sems = [], []
        for p in PAIRS:
            m, s = _sem(_per_animal(sess, p, col)); means.append(m); sems.append(s)
        ax.bar(xs + off, means, w, yerr=sems, capsize=3, label=lab, color=c, alpha=0.85)
    for i, p in enumerate(PAIRS):                       # paired Wilcoxon star
        k = _per_animal(sess, p, "cc_kcca"); l = _per_animal(sess, p, "cc_linear")
        d = [a - b for a, b in zip(k, l) if np.isfinite(a) and np.isfinite(b)]
        if len(d) >= 3:
            pv = paired_stats.wilcoxon_signed(d)[3]
            if np.isfinite(pv) and pv < 0.05:
                ax.text(xs[i], max(np.nanmean(k), np.nanmean(l)) + 0.02, "*",
                        ha="center", fontsize=12)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("held-out CC$_1$ (session)")
    ax.set_title(f"Nonlinear vs linear coupling — kernel CCA vs linear CCA — {tag}\n"
                 "(bars = mean±SEM over animals; * = paired Wilcoxon p<0.05; "
                 "KCCA≈linear ⇒ subspace is linear)", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(ATT / f"HCV1_KCCA_{tag}_vs_linear.png", dpi=150)
    plt.close(fig)


def fig_kcca_levels(df, tag):
    sess = df[df["scope"] == "session"]
    panels = [("cc_kcca", "held-out KCCA CC$_1$"), ("ifi", "IFI (>0: X leads Y)"),
              ("gini_x", "Gini (structure coeff, X)")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, (col, lab) in zip(axes, panels):
        means, sems = [], []
        for p in PAIRS:
            m, s = _sem(_per_animal(sess, p, col)); means.append(m); sems.append(s)
        ax.bar(np.arange(len(PAIRS)), means, yerr=sems, capsize=3, color="#8e44ad",
               alpha=0.8)
        if col == "ifi":
            ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(range(len(PAIRS))); ax.set_xticklabels(PAIRS, rotation=45,
                                                             ha="right", fontsize=7)
        ax.set_ylabel(lab)
    fig.suptitle(f"Kernel-CCA levels by pair — {tag} (session, mean±SEM over animals)",
                 fontsize=11)
    fig.tight_layout(); fig.savefig(ATT / f"HCV1_KCCA_{tag}_levels.png", dpi=150)
    plt.close(fig)


def fig_kcca_epochs(df, tag):
    ep = df[df["scope"].isin(EPOCHS)]
    if ep.empty:
        return
    panels = [("cc_kcca", "held-out KCCA CC$_1$"), ("ifi", "IFI"),
              ("gini_x", "Gini (structure coeff)")]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    xs = np.arange(len(PAIRS)); w = 0.26
    cols = {"naive": "#3498db", "intermediate": "#9b59b6", "expert": "#e74c3c"}
    for ax, (col, lab) in zip(axes, panels):
        for j, e in enumerate(EPOCHS):
            sub = ep[ep["scope"] == e]
            means = [_sem(_per_animal(sub, p, col))[0] for p in PAIRS]
            sems = [_sem(_per_animal(sub, p, col))[1] for p in PAIRS]
            ax.bar(xs + (j - 1) * w, means, w, yerr=sems, capsize=2, label=e,
                   color=cols[e], alpha=0.85)
        if col == "ifi":
            ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(xs); ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel(lab); ax.legend(fontsize=7)
    fig.suptitle(f"Kernel-CCA across learning epochs — {tag} (learners, mean±SEM)",
                 fontsize=11)
    fig.tight_layout(); fig.savefig(ATT / f"HCV1_KCCA_{tag}_epochs.png", dpi=150)
    plt.close(fig)


def fig_kcca_transition(tag, suf):
    path = config.RESULTS_DIR / f"kcca_transition{suf}.csv"
    if not path.is_file():
        return
    t = pd.read_csv(path)
    learn = t[t["learner"] == 1]
    panels = [("d_cc_kcca", "Δ KCCA CC$_1$"), ("d_ifi", "Δ IFI"),
              ("d_gini_x", "Δ Gini")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, (col, lab) in zip(axes, panels):
        means, sems, stars = [], [], []
        for p in PAIRS:
            vals = _num(learn[learn["pair"] == p][col]).dropna().tolist()
            m, s = _sem(vals); means.append(m); sems.append(s)
            stars.append(paired_stats.wilcoxon_signed(vals)[3] if len(vals) >= 3 else np.nan)
        ax.bar(np.arange(len(PAIRS)), means, yerr=sems, capsize=3, color="#16a085",
               alpha=0.8)
        for i, (m, pv) in enumerate(zip(means, stars)):
            if np.isfinite(pv) and pv < 0.05 and np.isfinite(m):
                ax.text(i, m, "*", ha="center", va="bottom" if m >= 0 else "top",
                        fontsize=11)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(range(len(PAIRS))); ax.set_xticklabels(PAIRS, rotation=45,
                                                             ha="right", fontsize=7)
        ax.set_ylabel(f"{lab} (cued − uncued)")
    fig.suptitle(f"Kernel-CCA uncued→cued — {tag} (learners, * = signed-rank p<0.05)",
                 fontsize=11)
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
