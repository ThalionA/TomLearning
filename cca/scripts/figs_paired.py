"""Paired-comparison panels for Tom (smoothed 10 ms, FS-included by default).

Two modes:
  * ``epochs``     — naive vs trained(expert), from epoch_metrics_bin10*/epoch_dims_bin10*
  * ``transition`` — uncued vs cued, from transition_uncued_cued_bin10*/transition_dims_bin10*

Per metric, an 8-pair grid: each panel shows every animal's two values connected by a
line, group mean ± SEM, and the paired test (t + Wilcoxon). Metrics:
  cc1     — held-out CC of the dominant dim
  n_sig   — number of significant canonical dims
  ifi     — information-flow index (±50 ms), >0 ⇒ first area leads
  gini_x  — Gini of the area-X canonical weights (participation sparsity)
  mincc   — mean held-out CC over the MIN # of significant dims shared by the two
            conditions (a strength measure at common dimensionality)

Usage:  PYTHONPATH=src python scripts/figs_paired.py epochs|transition [fsexcl] [learners|all]
        (defaults: FS-included; epochs are learners-only by construction)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402
from tom_cca import paired_stats  # noqa: E402

figstyle.apply()
ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"
RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB", "RSC-SUB", "V1-RSC"]
SIG, NS = "#1a7f4b", "#7a7972"
METRICS = [("cc1", "held-out CC$_1$"), ("n_sig", "# sig dims"),
           ("ifi", "IFI (>0: 1st leads)"), ("gini_x", "Gini (area X weights)"),
           ("mincc", "mean CC over min # sig dims")]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _mincc(dims, animal, pair, ea, eb=None):
    """Mean peak-CC over the top-(min n_sig) dims. For epochs eb is the other epoch;
    for transition the two 'epochs' are the uncued/cued dim rows already split."""
    def topk_cc(sub):
        cc = _num(sub["peak_cc" if "peak_cc" in sub else "cc"]).to_numpy()
        sig = int((_num(sub["sig"]) == 1).sum()) if "sig" in sub else int(np.isfinite(cc).sum())
        cc = np.sort(cc[np.isfinite(cc)])[::-1]
        return cc, sig
    a = dims[(dims["animal"] == animal) & (dims["pair"] == pair) & (dims["epoch"] == ea)]
    b = dims[(dims["animal"] == animal) & (dims["pair"] == pair) & (dims["epoch"] == eb)]
    if a.empty or b.empty:
        return np.nan, np.nan
    cca, sa = topk_cc(a); ccb, sb = topk_cc(b)
    k = max(1, min(sa, sb))
    if cca.size < k or ccb.size < k:
        return np.nan, np.nan
    return float(cca[:k].mean()), float(ccb[:k].mean())


def _paired_panel(ax, pair, va, vb, labels):
    """va, vb: per-animal arrays (already matched/paired). Draw connected pairs."""
    m = np.isfinite(va) & np.isfinite(vb)
    va, vb = va[m], vb[m]
    n = va.size
    if n < 1:
        ax.set_title(f"{pair}\n(no data)", fontsize=10, color=NS); ax.axis("off"); return
    for a, b in zip(va, vb):
        ax.plot([0, 1], [a, b], "-", color="#b4b2a9", lw=0.8, alpha=0.7, zorder=1)
    ax.scatter(np.zeros(n), va, s=22, color="#2c3e50", zorder=3)
    ax.scatter(np.ones(n), vb, s=22, color="#2c3e50", zorder=3)
    for x, v in [(0, va), (1, vb)]:
        ax.plot([x - 0.18, x + 0.18], [v.mean(), v.mean()], color="#c0392b", lw=2.4, zorder=4)
        ax.errorbar(x, v.mean(), yerr=v.std(ddof=1) / np.sqrt(n) if n > 1 else 0,
                    color="#c0392b", capsize=3, zorder=4)
    p_t = st.ttest_rel(va, vb)[1] if n >= 2 else np.nan
    _, _, _, p_w = paired_stats.wilcoxon_signed((vb - va).tolist())
    col = SIG if (np.isfinite(p_t) and p_t < 0.05) else "#3a3a36"
    ax.set_title(f"{pair}\nt p={p_t:.3f} · W p={p_w:.2g} · n={n}", fontsize=10,
                 color=col, fontweight=("bold" if col == SIG else "normal"))
    ax.set_xlim(-0.5, 1.5); ax.set_xticks([0, 1]); ax.set_xticklabels(labels, fontsize=9)


def main():
    args = sys.argv[1:]
    mode = "transition" if "transition" in args else "epochs"
    suf = "_fsincl" if "fsexcl" not in args else ""          # FS-incl default (Tom)
    fstag = "fsincl" if suf else "fsexcl"
    cohort_all = "all" in args
    labels = ("uncued", "cued") if mode == "transition" else ("naïve", "trained")

    if mode == "epochs":
        em = pd.read_csv(RES / f"epoch_metrics_bin10{suf}.csv")
        ed = pd.read_csv(RES / f"epoch_dims_bin10{suf}.csv")

        def series(metric, pair):
            g = em[em["pair"] == pair]
            va, vb = [], []
            for an, gg in g.groupby("animal"):
                if metric == "mincc":
                    a, b = _mincc(ed, an, pair, "naive", "expert"); va.append(a); vb.append(b)
                else:
                    va.append(_num(gg[gg.epoch == "naive"][metric]).mean())
                    vb.append(_num(gg[gg.epoch == "expert"][metric]).mean())
            return np.array(va, float), np.array(vb, float)
        cohort_lbl = "learners"            # epochs are learner-only by construction
    else:
        t = pd.read_csv(RES / f"transition_uncued_cued_bin10{suf}.csv")
        td = pd.read_csv(RES / f"transition_dims_bin10{suf}.csv")
        if not cohort_all:
            t = t[t["learner"] == 1]
        cohort_lbl = "all animals" if cohort_all else "learners"

        def series(metric, pair):
            g = t[t["pair"] == pair]
            if metric == "mincc":
                va, vb = [], []
                for an in g["animal"].unique():
                    sub = td[(td.animal == an) & (td.pair == pair)]
                    if sub.empty:
                        continue
                    cu = np.sort(_num(sub["unc_cc"]).dropna().to_numpy())[::-1]
                    cc = np.sort(_num(sub["cue_cc"]).dropna().to_numpy())[::-1]
                    su = int((_num(sub["unc_sig"]) == 1).sum()); sc = int((_num(sub["cue_sig"]) == 1).sum())
                    k = max(1, min(su, sc))
                    if cu.size >= k and cc.size >= k:
                        va.append(cu[:k].mean()); vb.append(cc[:k].mean())
                return np.array(va, float), np.array(vb, float)
            return _num(g[f"unc_{metric}"]).to_numpy(), _num(g[f"cue_{metric}"]).to_numpy()

    metrics = METRICS if mode == "epochs" or not cohort_all else \
        [m for m in METRICS if m[0] in ("cc1", "mincc")]   # Tom's E = CC1 + mincc only
    for metric, mlab in metrics:
        fig, axes = figstyle.grid(len(PAIRS), ncols=4)
        for ax, p in zip(axes, PAIRS):
            va, vb = series(metric, p)
            _paired_panel(ax, p, va, vb, labels)
        fig.suptitle(f"{labels[0]} vs {labels[1]} — {mlab} — {fstag}, {cohort_lbl}, 10 ms smoothed "
                     "(red = mean±SEM; paired t & Wilcoxon)", fontsize=12)
        coh = "all" if cohort_all else "learn"
        out = ATT / f"HCV1_paired_{mode}_{metric}_{fstag}_{coh}_bin10.png"
        figstyle.save(fig, out)
        print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
