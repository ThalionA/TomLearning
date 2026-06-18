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
        if star_p is not None and p in star_p:
            tier = _tier(star_p[p])
            if tier:
                top = (means[i] or 0) + (sems[i] or 0)
                ax.annotate(tier, (i, top), textcoords="offset points", xytext=(0, 5),
                            ha="center", va="bottom", color=SIG, fontweight="bold",
                            fontsize=14, zorder=20, clip_on=False)
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
    _, med, _, p = paired_stats.paired_t(sl)
    return med, p, len(sl)


# ---- cohort (pooled = all animals, vs learners-only) + clear-stats helpers ----

def _cohort(df, pool):
    """All animals when pool=True, else learners only."""
    return df if pool else df[df["learner"] == 1]


BINTAG = ""   # set in main() from the stem (e.g. "_bin10") so bins don't collide


def _ptag(tag, pool):
    return f"{tag}{'_pooled' if pool else ''}{BINTAG}"


def _clabel(df, pool):
    n = df["animal"].nunique() if pool else df[df["learner"] == 1]["animal"].nunique()
    return f"all {n} animals (learners + non)" if pool else f"{n} learners"


SIG = "#1a7f4b"      # significant (green)
NS = "#9a9990"       # n.s. (grey)


def _stars(p):
    """(asterisk string, colour) for a p-value; '' / grey when n.s. or undefined."""
    if not np.isfinite(p):
        return "n/a", NS
    if p < 0.001:
        return "*** p<0.001", SIG
    if p < 0.01:
        return f"** p={p:.3f}", SIG
    if p < 0.05:
        return f"* p={p:.3f}", SIG
    return f"n.s. p={p:.2f}", NS


def _tier(p):
    """Asterisk-tier only (for bar tops): '' when n.s./undefined."""
    if not np.isfinite(p) or p >= 0.05:
        return ""
    return "***" if p < 0.001 else ("**" if p < 0.01 else "*")


def _panel_stat(ax, pair, med, p, n, unit="slope"):
    """Non-occluding 2-line title: pair, then p (asterisk tier) · signed effect with an
    up/down arrow · n animals — all tinted green when significant. The test itself is
    named once in the figure suptitle, not repeated per panel."""
    star, col = _stars(p)
    arrow = " ↓" if med < 0 else (" ↑" if med > 0 else "")
    eff = f"{unit} {med:+.3f}{arrow}" if np.isfinite(med) else f"{unit} n/a"
    tcol = "#14532d" if col == SIG else "#3a3a36"
    ax.set_title(f"{pair}\n{star}  ·  {eff}  ·  n={n}",
                 fontsize=11, fontweight="bold", color=tcol, linespacing=1.5)


# --------------------------------------------------------------------------

def fig_levels(df, tag, pool=False):
    data = _cohort(df, pool)
    metrics = [("cc1", "held-out CC$_1$"), ("n_sig", "# sig dims"),
               ("mi_sig", "MI$_{sig}$ (nats)"), ("ifi", "IFI"),
               ("gini_x", "Gini (area X)")]
    fig, axes = plt.subplots(1, len(metrics), figsize=(3 * len(metrics), 3.8))
    for ax, (m, lab) in zip(axes, metrics):
        series = {p: per_animal_level(data, p, m) for p in PAIRS}
        _bar_points_sem(ax, PAIRS, series, lab, lab)
    fig.suptitle(f"Subspace levels by pair — {tag}, {_clabel(df, pool)} "
                 "(points = animals, bar = mean ± SEM)")
    figstyle.save(fig, ATT / f"HCV1_CCA_{_ptag(tag, pool)}_levels.png")


def fig_slopes_heatmap(df, tag, pool=False):
    data = _cohort(df, pool)
    metrics = [m for m in SLOPE_METRICS if m in df.columns]
    for axis in AXES:
        M = np.full((len(PAIRS), len(metrics)), np.nan)
        P = np.full_like(M, np.nan)
        for i, p in enumerate(PAIRS):
            for j, m in enumerate(metrics):
                med, pp, n = _slope_quant(data, m, axis, p)
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
                    tier = _tier(P[i, j])
                    sig = bool(tier)
                    col = "white" if abs(norm[i, j]) > 0.45 else "black"   # adaptive contrast
                    ax.text(j, i, f"{M[i, j]:+.3f}\n{tier}" if tier else f"{M[i, j]:+.3f}",
                            ha="center", va="center", fontsize=7.5, color=col,
                            fontweight=("bold" if sig else "normal"))
        fig.colorbar(im, ax=ax, label="median per-animal slope (col-normalised)")
        ax.set_title(f"All relationships — d(metric)/d({axis}) — {tag}, {_clabel(df, pool)}\n"
                     "paired t on per-animal slopes  (* p<.05  ** p<.01  *** p<.001, "
                     "uncorrected)", fontsize=9)
        figstyle.save(fig, ATT / f"HCV1_CCA_{_ptag(tag, pool)}_slopes_{axis}.png")


def _grid_traj(df, metric, axis, fname, ylabel, suptitle, hline0=False, pool=False):
    """Spacious 4×2 grid: one panel per area pair, mean±SEM band + faint per-animal.
    pool=True uses all animals (learners + non); else learners only. Per panel a clear
    stats caption: paired t test, p with asterisk tier, signed slope, n animals."""
    data = _cohort(df, pool)
    fig, axflat = figstyle.grid(len(PAIRS), ncols=2)
    for i, (ax, p) in enumerate(zip(axflat, PAIRS)):
        _band(ax, data, p, metric, axis)
        if hline0:
            ax.axhline(0, color="k", lw=0.6)
        med, pv, n = _slope_quant(data, metric, axis, p)
        _panel_stat(ax, p, med, pv, n)
        ax.margins(y=0.15)
        if i % 2 == 0:
            ax.set_ylabel(ylabel)
        if i >= len(PAIRS) - 2:
            ax.set_xlabel(axis)
    fig.suptitle(suptitle, fontsize=13)
    figstyle.save(fig, ATT / fname)


def fig_gini_traj(df, tag, pool=False):
    _grid_traj(df, "gini_x", "trial_frac", f"HCV1_CCA_{_ptag(tag, pool)}_gini_traj.png",
               "Gini (area X)",
               f"Subspace sparsity vs learning — {tag}, {_clabel(df, pool)}\n"
               "band = mean±SEM, faint = per-animal · per panel: paired t on per-animal "
               "slopes vs 0  (* p<.05  ** p<.01  *** p<.001)", pool=pool)


def fig_ifi_traj(df, tag, pool=False):
    _grid_traj(df, "ifi", "trial_frac", f"HCV1_CCA_{_ptag(tag, pool)}_ifi_traj.png",
               "IFI (>0: X leads Y)",
               f"Directionality vs learning — {tag}, {_clabel(df, pool)}\n"
               "band = mean±SEM, faint = per-animal · per panel: paired t on per-animal "
               "slopes vs 0  (* p<.05  ** p<.01  *** p<.001)", hline0=True, pool=pool)


def fig_transition(pool=False):
    path = RES / f"transition_uncued_cued{BINTAG}.csv"
    if not path.is_file() or path.stat().st_size < 50:
        return
    t = pd.read_csv(path)
    data = _cohort(t, pool)
    metrics = [("d_cc1", "Δ CC$_1$"), ("d_n_sig", "Δ n_sig"),
               ("d_mi_sig", "Δ MI$_{sig}$"), ("d_ifi", "Δ IFI")]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for ax, (m, lab) in zip(axes.ravel(), metrics):
        ser = {p: _num(data[data["pair"] == p][m]).dropna().tolist() for p in PAIRS}
        star = {p: paired_stats.paired_t(ser[p])[3] if len(ser[p]) >= 3 else np.nan
                for p in PAIRS}
        _bar_points_sem(ax, PAIRS, ser, f"{lab} (cued − uncued)",
                        f"{lab}: uncued→cued", star_p=star)
    fig.suptitle(f"Uncued→cued transition — {_clabel(t, pool)}\n(points = animals, bar = mean ± "
                 "SEM · paired t vs 0: * p<.05  ** p<.01  *** p<.001)", fontsize=12)
    figstyle.save(fig, ATT / f"HCV1_CCA_transition{BINTAG}{'_pooled' if pool else ''}.png")


def fig_direction(df, tag, pool=False):
    data = _cohort(df, pool)
    series, star = {}, {}
    for p in PAIRS:
        per_an = [float(_num(s["ifi"]).mean()) for _, s in data[data["pair"] == p].groupby("animal")]
        series[p] = per_an
        star[p] = paired_stats.paired_t([v for v in per_an if np.isfinite(v)])[3] \
            if sum(np.isfinite(per_an)) >= 3 else np.nan
    fig, ax = plt.subplots(figsize=(7, 4))
    _bar_points_sem(ax, PAIRS, series, "mean IFI (>0: first area leads)",
                    f"Directionality by pair — {tag}, {_clabel(df, pool)}\n"
                    "paired t vs 0: * p<.05  ** p<.01  *** p<.001", star_p=star)
    figstyle.save(fig, ATT / f"HCV1_CCA_{_ptag(tag, pool)}_direction.png")


def fig_rotation_floor(df, tag, pool=False):
    data = _cohort(df, pool)
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    xs = np.arange(len(PAIRS)); w = 0.38
    rot = [_sem(per_animal_level(data, p, "rot_x")) for p in PAIRS]
    flo = [_sem(per_animal_level(data, p, "sh_x")) for p in PAIRS]
    ax.bar(xs - w / 2, [r[0] for r in rot], w, yerr=[r[1] for r in rot], capsize=3,
           label="cross-window rotation", color="#2980b9", alpha=0.8)
    ax.bar(xs + w / 2, [f[0] for f in flo], w, yerr=[f[1] for f in flo], capsize=3,
           label="split-half noise floor", color="#95a5a6", alpha=0.8)
    for i, p in enumerate(PAIRS):     # paired rotation-vs-floor paired t per pair
        g = data[data["pair"] == p]
        diffs = [r - h for r, h in
                 ((_num(s["rot_x"]).mean(), _num(s["sh_x"]).mean()) for _, s in g.groupby("animal"))
                 if np.isfinite(r) and np.isfinite(h)]
        tier = _tier(paired_stats.paired_t(diffs)[3]) if len(diffs) >= 3 else ""
        if tier:
            top = max(rot[i][0] + rot[i][1], flo[i][0] + flo[i][1])
            ax.annotate(tier, (i, top), textcoords="offset points", xytext=(0, 5),
                        ha="center", color=SIG, fontweight="bold", fontsize=13,
                        zorder=20, clip_on=False)
    ax.set_xticks(xs); ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("max principal angle (°)")
    ax.set_title(f"Rotation vs noise floor — {tag}, {_clabel(df, pool)}\nrotation>floor ⇒ "
                 "reorientation (here rotation ≤ floor) · paired t rot vs floor: "
                 "* p<.05  ** p<.01  *** p<.001", fontsize=9)
    ax.legend(fontsize=8)
    figstyle.save(fig, ATT / f"HCV1_CCA_{_ptag(tag, pool)}_rotation_floor.png")


def fig_gini_control(df, tag, axis="trial_frac", pool=False):
    """CCA-independent control: per-pair median slope of the canonical-weight Gini
    (gini_x) vs the raw cross-correlation coupling Gini (gini_pearson_x). If the
    de-sparsification is real, both fall together — Pearson is computed with no CCA."""
    if "gini_pearson_x" not in df.columns:
        return
    data = _cohort(df, pool)
    fig, ax = plt.subplots(figsize=(8.5, 4.3))
    xs = np.arange(len(PAIRS)); w = 0.38
    for off, metric, lab, col in [(-w / 2, "gini_x", "Gini (CCA weights)", "#c0392b"),
                                  (+w / 2, "gini_pearson_x",
                                   "Gini (Pearson coupling — CCA-free)", "#8e44ad")]:
        slopes = per_animal_slopes(data, metric, axis)     # once per metric
        meds, sems, stars = [], [], []
        for p in PAIRS:
            sl = slopes.get(p, [])
            if len(sl) >= 3:
                _, med, _, pv = paired_stats.paired_t(sl)
            else:
                med, pv = np.nan, np.nan
            meds.append(med); sems.append(_sem(sl)[1]); stars.append(pv)
        ax.bar(xs + off, meds, w, yerr=sems, capsize=3, label=lab, color=col, alpha=0.85)
        for i, (m, pv) in enumerate(zip(meds, stars)):
            tier = _tier(pv)
            if tier and np.isfinite(m):
                ax.annotate(tier, (xs[i] + off, m), textcoords="offset points",
                            xytext=(0, 4 if m >= 0 else -12), ha="center", color=SIG,
                            fontweight="bold", fontsize=12, zorder=20, clip_on=False)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(PAIRS, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"median per-animal slope d(Gini)/d({axis})")
    ax.set_title(f"De-sparsification control — CCA-weight Gini vs CCA-free Pearson Gini — "
                 f"{tag}, {_clabel(df, pool)}\n(both negative ⇒ not a CCA-weight artefact; shares "
                 "the residualisation · paired t: * p<.05  ** p<.01  *** p<.001)", fontsize=8)
    ax.legend(fontsize=8)
    figstyle.save(fig, ATT / f"HCV1_CCA_{_ptag(tag, pool)}_gini_control.png")


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
    figstyle.save(fig, ATT / f"HCV1_CCA_{_ptag(tag, False)}_learnervsnon_{metric}.png")


def make_all(csv, tag):
    df = pd.read_csv(csv)
    for pool in (False, True):             # learners-only AND pooled (all animals)
        fig_levels(df, tag, pool)
        fig_slopes_heatmap(df, tag, pool)
        fig_gini_traj(df, tag, pool)
        fig_ifi_traj(df, tag, pool)
        fig_direction(df, tag, pool)
        if "sh_x" in df.columns:
            fig_rotation_floor(df, tag, pool)
        fig_gini_control(df, tag, pool=pool)
    fig_learner_vs_non(df, tag)            # the learner/non contrast — pooling N/A
    print(f"  {tag}: figures written ({len(df)} rows, learners + pooled)")


def main():
    ATT.mkdir(parents=True, exist_ok=True)
    # optional positional arg: trajectory CSV stem (default = committed; pass e.g.
    # "trajectory_w15_bin25" for the window=15 re-run)
    stem = sys.argv[1] if len(sys.argv) > 1 else "trajectory_windows"
    global BINTAG
    BINTAG = "_bin10" if "bin10" in stem else ("_bin50" if "bin50" in stem else "")
    for suffix, tag in [("", "fsexcl"), ("_fsincl", "fsincl")]:
        csv = RES / f"{stem}{suffix}.csv"
        if csv.is_file():
            make_all(csv, tag)
        else:
            print(f"  (missing {csv.name} — skip {tag})")
    for pool in (False, True):
        fig_transition(pool)               # all pairs; transition CSV is FS-excluded
    print(f"figures ({stem}) ->", ATT)


if __name__ == "__main__":
    main()
