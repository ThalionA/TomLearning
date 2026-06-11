"""Report figures from the committed full-suite CSVs (no re-fitting).

Per FS condition (FS-excluded = trajectory_windows.csv; FS-included =
trajectory_windows_fsincl.csv) renders, into the ResearchVault attachments:
  * levels         -- per-pair metric levels, individual animal points + SEM
  * slopes_<axis>  -- ALL pairs x ALL metrics: median per-animal slope heatmap
                      with sign-test significance (shows every relationship)
  * gini_traj      -- Gini vs learning: across-animal mean +- SEM band + faint
                      per-animal lines (CA1-RSC / CA1-DG / CA1-V1)
  * ifi_traj       -- CA1->DG IFI vs learning, same band style
  * direction      -- per-pair mean IFI, points + SEM, sign test vs 0
  * rotation_floor -- cross-window rotation vs the split-half noise floor
Plus cross-cut comparisons: Gini-slope FS-excl vs FS-incl, and learners vs
non-learners. Every plot: axis labels with units, title, legend, quantified.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, paired_stats, trajectory  # noqa: E402
import figstyle  # noqa: E402

figstyle.apply()

ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")
RES = config.RESULTS_DIR
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
AXES = ["trial_frac", "performance", "lp_rel"]
SLOPE_METRICS = ["cc1", "n_sig", "mi_sig", "ifi", "gini_x", "gini_pearson_x",
                 "rot_x", "jac_x"]


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def _sem(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    return (np.nan, np.nan, 0) if v.size == 0 else (
        float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(v.size)) if v.size > 1 else 0.0,
        v.size)


def per_animal_level(df, pair, metric):
    g = df[df["pair"] == pair]
    return [float(_num(s[metric]).mean()) for _, s in g.groupby("animal")]


def per_animal_slopes(df, metric, axis):
    out = defaultdict(list)
    for (an, pr), g in df.groupby(["animal", "pair"]):
        x = _num(g[axis]).to_numpy(); y = _num(g[metric]).to_numpy()
        if np.sum(np.isfinite(x) & np.isfinite(y)) >= 4:
            s, _ = trajectory.linear_slope(x, y)
            if np.isfinite(s):
                out[pr].append(s)
    return out


def _bar_points_sem(ax, pairs, series, ylabel, title, star_p=None):
    """series: dict pair -> list of per-animal values."""
    xs = np.arange(len(pairs))
    means = [_sem(series.get(p, []))[0] for p in pairs]
    sems = [_sem(series.get(p, []))[1] for p in pairs]
    ax.bar(xs, means, yerr=sems, capsize=3, color="#7f8c8d", alpha=0.7, zorder=1)
    for i, p in enumerate(pairs):
        v = [x for x in series.get(p, []) if np.isfinite(x)]
        jit = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.3
        ax.scatter(np.full(len(v), i) + jit, v, s=12, color="#2c3e50",
                   zorder=2, alpha=0.8)
        if star_p is not None and p in star_p and np.isfinite(star_p[p]) and star_p[p] < 0.05:
            top = (means[i] or 0) + (sems[i] or 0)
            figstyle.star(ax, i, top, always=True)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(pairs, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10)


def _band(ax, df, pair, metric, axis, nbins=6, lo=0.0, hi=1.0):
    centres = (np.linspace(lo, hi, nbins + 1)[:-1] + np.linspace(lo, hi, nbins + 1)[1:]) / 2
    edges = np.linspace(lo, hi, nbins + 1)
    curves = []
    for an, g in df[df["pair"] == pair].groupby("animal"):
        x = _num(g[axis]).to_numpy(); y = _num(g[metric]).to_numpy()
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 3:
            continue
        binned = np.full(nbins, np.nan)
        for b in range(nbins):
            m = ok & (x >= edges[b]) & (x < edges[b + 1] if b < nbins - 1 else x <= edges[b + 1])
            if m.any():
                binned[b] = np.mean(y[m])
        curves.append(binned)
        ax.plot(centres, binned, color="#bdc3c7", lw=0.7, alpha=0.6, zorder=1)
    if not curves:
        return
    arr = np.array(curves)
    m = np.nanmean(arr, 0)
    n = np.sum(np.isfinite(arr), 0)
    sem = np.nanstd(arr, 0, ddof=1) / np.sqrt(np.maximum(n, 1))
    ax.plot(centres, m, color="#c0392b", lw=2, zorder=3, label="mean ± SEM")
    ax.fill_between(centres, m - sem, m + sem, color="#c0392b", alpha=0.25, zorder=2)


def _slope_quant(df, metric, axis, pair):
    sl = per_animal_slopes(df, metric, axis).get(pair, [])
    if len(sl) < 3:
        return np.nan, np.nan, len(sl)
    _, med, _, p = paired_stats.wilcoxon_signed(sl)
    return med, p, len(sl)


# --------------------------------------------------------------------------

def fig_levels(df, tag):
    learn = df[df["learner"] == 1]
    metrics = [("cc1", "held-out CC$_1$"), ("n_sig", "# sig dims"),
               ("mi_sig", "MI$_{sig}$ (nats)"), ("ifi", "IFI"),
               ("gini_x", "Gini (area X)")]
    fig, axes = plt.subplots(1, len(metrics), figsize=(3 * len(metrics), 3.8))
    for ax, (m, lab) in zip(axes, metrics):
        series = {p: per_animal_level(learn, p, m) for p in PAIRS}
        _bar_points_sem(ax, PAIRS, series, lab, lab)
    fig.suptitle(f"Subspace levels by pair — {tag} (points = animals, bar = mean ± SEM, learners)")
    figstyle.save(fig, ATT / f"HCV1_CCA_{tag}_levels.png")


def fig_slopes_heatmap(df, tag):
    learn = df[df["learner"] == 1]
    metrics = [m for m in SLOPE_METRICS if m in df.columns]
    for axis in AXES:
        M = np.full((len(PAIRS), len(metrics)), np.nan)
        P = np.full_like(M, np.nan)
        for i, p in enumerate(PAIRS):
            for j, m in enumerate(metrics):
                med, pp, n = _slope_quant(learn, m, axis, p)
                M[i, j] = med; P[i, j] = pp
        # normalise each metric column to its max abs for visual comparability
        norm = M / (np.nanmax(np.abs(M), 0, keepdims=True) + 1e-12)
        fig, ax = plt.subplots(figsize=(1.1 * len(metrics) + 2, 0.6 * len(PAIRS) + 2))
        im = ax.imshow(norm, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(len(metrics))); ax.set_xticklabels(metrics, rotation=45, ha="right")
        ax.set_yticks(range(len(PAIRS))); ax.set_yticklabels(PAIRS)
        for i in range(len(PAIRS)):
            for j in range(len(metrics)):
                if np.isfinite(M[i, j]):
                    star = "*" if (np.isfinite(P[i, j]) and P[i, j] < 0.05) else ""
                    ax.text(j, i, f"{M[i, j]:+.3f}{star}", ha="center", va="center",
                            fontsize=7, color="black")
        fig.colorbar(im, ax=ax, label="median per-animal slope (col-normalised)")
        ax.set_title(f"All relationships — d(metric)/d({axis}) — {tag}, learners\n"
                     "* = signed-rank p<0.05 (uncorrected)", fontsize=9)
        figstyle.save(fig, ATT / f"HCV1_CCA_{tag}_slopes_{axis}.png")


def _grid_traj(df, metric, axis, tag, fname, ylabel, suptitle, hline0=False):
    """Spacious 4×2 grid: one panel per area pair, mean±SEM band + faint per-animal."""
    learn = df[df["learner"] == 1]
    fig, axflat = figstyle.grid(len(PAIRS), ncols=2)
    for i, (ax, p) in enumerate(zip(axflat, PAIRS)):
        _band(ax, learn, p, metric, axis)
        if hline0:
            ax.axhline(0, color="k", lw=0.6)
        med, pv, n = _slope_quant(learn, metric, axis, p)
        star = " ✶" if (np.isfinite(pv) and pv < 0.05) else ""
        ax.set_title(f"{p}   slope={med:+.3f}  p={pv:.3g}  n={n}{star}")
        ax.margins(y=0.15)
        if i % 2 == 0:
            ax.set_ylabel(ylabel)
        if i >= len(PAIRS) - 2:
            ax.set_xlabel(axis)
    fig.suptitle(suptitle, fontsize=13)
    figstyle.save(fig, ATT / fname)


def fig_gini_traj(df, tag):
    _grid_traj(df, "gini_x", "trial_frac", tag, f"HCV1_CCA_{tag}_gini_traj.png",
               "Gini (area X)",
               f"Subspace sparsity vs learning — ALL pairs — {tag} "
               "(red = mean±SEM, faint = animals)")


def fig_ifi_traj(df, tag):
    _grid_traj(df, "ifi", "trial_frac", tag, f"HCV1_CCA_{tag}_ifi_traj.png",
               "IFI (>0: X leads Y)",
               f"Directionality vs learning — ALL pairs — {tag} "
               "(red = mean±SEM, faint = animals)", hline0=True)


def fig_transition():
    path = RES / "transition_uncued_cued.csv"
    if not path.is_file():
        return
    t = pd.read_csv(path)
    learn = t[t["learner"] == 1]
    metrics = [("d_cc1", "Δ CC$_1$"), ("d_n_sig", "Δ n_sig"),
               ("d_mi_sig", "Δ MI$_{sig}$"), ("d_ifi", "Δ IFI")]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for ax, (m, lab) in zip(axes.ravel(), metrics):
        ser = {p: _num(learn[learn["pair"] == p][m]).dropna().tolist() for p in PAIRS}
        star = {p: paired_stats.wilcoxon_signed(ser[p])[3] if len(ser[p]) >= 3 else np.nan
                for p in PAIRS}
        _bar_points_sem(ax, PAIRS, ser, f"{lab} (cued − uncued)",
                        f"{lab}: uncued→cued", star_p=star)
    fig.suptitle("Uncued→cued transition — ALL pairs (points = animals, bar = mean ± SEM, "
                 "* = signed-rank p<0.05)", fontsize=12)
    figstyle.save(fig, ATT / "HCV1_CCA_transition.png")


def fig_direction(df, tag):
    learn = df[df["learner"] == 1]
    series, star = {}, {}
    for p in PAIRS:
        per_an = [float(_num(s["ifi"]).mean()) for _, s in learn[learn["pair"] == p].groupby("animal")]
        series[p] = per_an
        star[p] = paired_stats.wilcoxon_signed([v for v in per_an if np.isfinite(v)])[3] \
            if sum(np.isfinite(per_an)) >= 3 else np.nan
    fig, ax = plt.subplots(figsize=(7, 4))
    _bar_points_sem(ax, PAIRS, series, "mean IFI (>0: first area leads)",
                    f"Directionality by pair — {tag} (* = p<0.05 vs 0)", star_p=star)
    figstyle.save(fig, ATT / f"HCV1_CCA_{tag}_direction.png")


def fig_rotation_floor(df, tag):
    learn = df[df["learner"] == 1]
    fig, ax = plt.subplots(figsize=(8, 4))
    xs = np.arange(len(PAIRS)); w = 0.38
    rot = [_sem(per_animal_level(learn, p, "rot_x")) for p in PAIRS]
    flo = [_sem(per_animal_level(learn, p, "sh_x")) for p in PAIRS]
    ax.bar(xs - w / 2, [r[0] for r in rot], w, yerr=[r[1] for r in rot], capsize=3,
           label="cross-window rotation", color="#2980b9", alpha=0.8)
    ax.bar(xs + w / 2, [f[0] for f in flo], w, yerr=[f[1] for f in flo], capsize=3,
           label="split-half noise floor", color="#95a5a6", alpha=0.8)
    ax.set_xticks(xs); ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("max principal angle (°)")
    ax.set_title(f"Rotation vs noise floor — {tag} (rotation>floor ⇒ genuine reorientation)",
                 fontsize=10)
    ax.legend(fontsize=8)
    figstyle.save(fig, ATT / f"HCV1_CCA_{tag}_rotation_floor.png")


def fig_gini_control(df, tag, axis="trial_frac"):
    """CCA-independent control: per-pair median slope of the canonical-weight Gini
    (gini_x) vs the raw cross-correlation coupling Gini (gini_pearson_x). If the
    de-sparsification is real, both fall together — Pearson is computed with no CCA."""
    if "gini_pearson_x" not in df.columns:
        return
    learn = df[df["learner"] == 1]
    fig, ax = plt.subplots(figsize=(8.5, 4))
    xs = np.arange(len(PAIRS)); w = 0.38
    for off, metric, lab, col in [(-w / 2, "gini_x", "Gini (CCA weights)", "#c0392b"),
                                  (+w / 2, "gini_pearson_x",
                                   "Gini (Pearson coupling — CCA-free)", "#8e44ad")]:
        slopes = per_animal_slopes(learn, metric, axis)     # once per metric
        meds, sems, stars = [], [], []
        for p in PAIRS:
            sl = slopes.get(p, [])
            if len(sl) >= 3:
                _, med, _, pv = paired_stats.wilcoxon_signed(sl)
            else:
                med, pv = np.nan, np.nan
            meds.append(med); sems.append(_sem(sl)[1]); stars.append(pv)
        ax.bar(xs + off, meds, w, yerr=sems, capsize=3, label=lab, color=col, alpha=0.85)
        for i, (m, pv) in enumerate(zip(meds, stars)):
            if np.isfinite(pv) and pv < 0.05 and np.isfinite(m):
                ax.text(xs[i] + off, m, "*", ha="center",
                        va="bottom" if m >= 0 else "top", fontsize=11)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"median per-animal slope d(Gini)/d({axis})")
    ax.set_title(f"De-sparsification control — CCA-weight Gini vs CCA-free Pearson Gini — "
                 f"{tag}, learners\n(both negative ⇒ not a CCA-weight artefact; shares the "
                 "residualisation, so not orthogonal to it; * = signed-rank p<0.05)", fontsize=8)
    ax.legend(fontsize=8)
    figstyle.save(fig, ATT / f"HCV1_CCA_{tag}_gini_control.png")


def fig_learner_vs_non(df, tag, metric="gini_x", axis="trial_frac"):
    fig, ax = plt.subplots(figsize=(8, 4))
    xs = np.arange(len(PAIRS)); w = 0.38
    for off, lv, lab, col in [(-w / 2, 1, "learners", "#c0392b"),
                              (+w / 2, 0, "non-learners", "#27ae60")]:
        sub = df[df["learner"] == lv]
        sl = per_animal_slopes(sub, metric, axis)
        means = [_sem(sl.get(p, []))[0] for p in PAIRS]
        sems = [_sem(sl.get(p, []))[1] for p in PAIRS]
        ax.bar(xs + off, means, w, yerr=sems, capsize=3, label=lab, color=col, alpha=0.8)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"median slope d({metric})/d({axis})")
    ax.set_title(f"Learners vs non-learners — {metric} slope — {tag}", fontsize=10)
    ax.legend(fontsize=8)
    figstyle.save(fig, ATT / f"HCV1_CCA_{tag}_learnervsnon_{metric}.png")


def make_all(csv, tag):
    df = pd.read_csv(csv)
    fig_levels(df, tag)
    fig_slopes_heatmap(df, tag)
    fig_gini_traj(df, tag)
    fig_ifi_traj(df, tag)
    fig_direction(df, tag)
    if "sh_x" in df.columns:
        fig_rotation_floor(df, tag)
    fig_gini_control(df, tag)
    fig_learner_vs_non(df, tag)
    print(f"  {tag}: figures written ({len(df)} rows)")


def main():
    ATT.mkdir(parents=True, exist_ok=True)
    # optional positional arg: trajectory CSV stem (default = committed; pass e.g.
    # "trajectory_w15_bin25" for the window=15 re-run)
    stem = sys.argv[1] if len(sys.argv) > 1 else "trajectory_windows"
    for suffix, tag in [("", "fsexcl"), ("_fsincl", "fsincl")]:
        csv = RES / f"{stem}{suffix}.csv"
        if csv.is_file():
            make_all(csv, tag)
        else:
            print(f"  (missing {csv.name} — skip {tag})")
    fig_transition()                       # all pairs; transition CSV is FS-excluded
    print(f"figures ({stem}) ->", ATT)


if __name__ == "__main__":
    main()
