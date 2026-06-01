"""Config-centric exploration of the spatial CCA sweep.

Complements ``summarise_sweep.py`` (which is organised *per area-pair*). This
script is organised *per config* -- it answers "which of the 66 hyperparameter
configurations are interesting, and how do the swept axes drive the results".

Reads ``figures/sweep_summary_spatial.csv`` (the long per-(config x pair) table
written by ``summarise_sweep.py``) and writes, into ``figures/sweep_explore/``:

  * sweep_redundancy.png      66x66 identical-result matrix -- how many of the
                              66 configs are genuinely distinct.
  * sweep_main_effects.png    mean CC / dCC / #significant pairs as a function
                              of the PC-count rule, split by CCA type.
  * sweep_config_ranking.png  every distinct config ranked by total significant
                              pairs (strength + IFI); committed config marked.
  * sweep_config_scatter.png  effect size vs significance, one point per config.
  * sweep_committed_context.png  where the committed config sits in the sweep
                              distribution of each headline metric.
  * sweep_pair_structure.png  per-pair: dCC across configs + robustness counts.

The IFI lag-scan window (5/10/20) leaves every headline metric here unchanged
(it only extends the scan range), so figures 2-6 use the lag=10 representative
of each (CCA type x PC-rule) cell -- 22 configs. Figure 1 shows all 66.

Run:  python scripts/plot_sweep_explore.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tom_cca import config  # noqa: E402

CSV_PATH = config.FIGURES_DIR / "sweep_summary_spatial.csv"
OUTDIR = config.FIGURES_DIR / "sweep_explore"
ALPHA = 0.05
COMMITTED = "res_samp15_lag10"
PAIRS = [f"{x}-{y}" for x, y in config.PAIRS]

# PC-count rules, ordered low -> high effective dimensionality within family.
KRULE_ORDER = ["fixed3", "fixed5", "fixed10", "fixed20", "fixed30",
               "samples15", "samples25", "samples40", "var75", "var85", "var95"]
RES_COLOUR = "#3c6fb0"
SIG_COLOUR = "#d4753a"


# --------------------------------------------------------------------------
# data access
# --------------------------------------------------------------------------
def _f(x):
    """CSV cell -> float, '' -> NaN."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def load_rows():
    if not CSV_PATH.exists():
        sys.exit(f"missing {CSV_PATH} -- run summarise_sweep.py first")
    with open(CSV_PATH) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"{CSV_PATH} is empty")
    return rows


def by_tag(rows):
    """{tag: {'cca','k_rule','max_lag', pair-> metric-dict}}."""
    out: dict[str, dict] = {}
    for r in rows:
        t = r["tag"]
        d = out.setdefault(t, {"cca": r["cca"], "k_rule": r["k_rule"],
                               "max_lag": int(_f(r["max_lag"])), "pairs": {}})
        d["pairs"][r["pair"]] = r
    return out


def config_aggregate(entry):
    """Per-config summary over the eight area-pairs."""
    cc, dcc, amf, gini = [], [], [], []
    n_sig_strength = n_sig_ifi = 0
    dcc_pos = 0
    for pair in PAIRS:
        m = entry["pairs"].get(pair, {})
        for key in ("cc_naive", "cc_expert"):
            v = _f(m.get(key))
            if np.isfinite(v):
                cc.append(v)
        d = _f(m.get("d_cc"))
        if np.isfinite(d):
            dcc.append(d)
            dcc_pos += d > 0
        a = _f(m.get("angle_minus_floor"))
        if np.isfinite(a):
            amf.append(a)
        g = _f(m.get("gini_expert"))
        if np.isfinite(g):
            gini.append(g)
        if _f(m.get("p_naive_vs_expert")) < ALPHA:
            n_sig_strength += 1
        if _f(m.get("p_ifi_w5")) < ALPHA:
            n_sig_ifi += 1
    mean = lambda v: float(np.mean(v)) if v else np.nan  # noqa: E731
    return {
        "mean_cc": mean(cc), "mean_dcc": mean(dcc),
        "mean_amf": mean(amf), "mean_gini": mean(gini),
        "n_sig_strength": n_sig_strength, "n_sig_ifi": n_sig_ifi,
        "n_sig_total": n_sig_strength + n_sig_ifi, "dcc_pos": dcc_pos,
    }


def short_label(cca, k_rule):
    """Compact config label, e.g. 'res samp15' / 'sig var85'."""
    c = "res" if cca == "residual" else "sig"
    k = (k_rule.replace("samples", "samp").replace("fixed", "fix"))
    return f"{c} {k}"


def _save(fig, name):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    path = OUTDIR / f"{name}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"saved {path}")


# --------------------------------------------------------------------------
# figure 1 -- redundancy: how many of the 66 configs are distinct
# --------------------------------------------------------------------------
def fig_redundancy(tags):
    """66x66 matrix: cell dark where two configs share an identical d_cc
    vector across the eight pairs (-> identical headline result)."""
    order = sorted(tags, key=lambda t: (
        0 if tags[t]["cca"] == "residual" else 1,
        KRULE_ORDER.index(tags[t]["k_rule"]), tags[t]["max_lag"]))
    vecs = {}
    for t in order:
        vecs[t] = tuple(round(_f(tags[t]["pairs"][p].get("d_cc")), 6)
                        for p in PAIRS)
    n = len(order)
    same = np.zeros((n, n), bool)
    for i, ti in enumerate(order):
        for j, tj in enumerate(order):
            same[i, j] = vecs[ti] == vecs[tj]
    n_distinct = len({vecs[t] for t in order})

    fig, ax = plt.subplots(figsize=(9.2, 8.4))
    ax.imshow(same, cmap="Greys", interpolation="nearest", vmin=0, vmax=1)
    split = sum(1 for t in order if tags[t]["cca"] == "residual")
    for s in (split - 0.5,):
        ax.axhline(s, color="#c0392b", lw=1.4)
        ax.axvline(s, color="#c0392b", lw=1.4)
    ax.set_xticks([split / 2, split + (n - split) / 2])
    ax.set_xticklabels(["residual configs", "signal configs"])
    ax.set_yticks([split / 2, split + (n - split) / 2])
    ax.set_yticklabels(["residual", "signal"], rotation=90, va="center")
    ax.set_title(
        f"Sweep redundancy -- {n} configs collapse to {n_distinct} distinct "
        f"headline results\n(dark = identical dCC across all 8 pairs;  lag "
        f"5/10/20 identical;  samp15/25/40 + fix30 share one k)", fontsize=10)
    _save(fig, "sweep_redundancy")
    return n_distinct


# --------------------------------------------------------------------------
# figure 2 -- main effects of the PC-count rule and CCA type
# --------------------------------------------------------------------------
def fig_main_effects(reps):
    """mean CC / dCC / #sig-strength / #sig-IFI vs PC-rule, split by CCA type."""
    metrics = [("mean_cc", "mean held-out CC", None),
               ("mean_dcc", "mean dCC (expert - naive)", 0.0),
               ("n_sig_strength", "# pairs, strength change p<.05  (of 8)", None),
               ("n_sig_ifi", "# pairs, IFI |lag|<=5 p<.05  (of 8)", None)]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.6))
    x = np.arange(len(KRULE_ORDER))
    for ax, (key, label, ref) in zip(axes.ravel(), metrics):
        for cca, colour in (("residual", RES_COLOUR), ("signal", SIG_COLOUR)):
            y = [config_aggregate(reps[(cca, k)])[key] if (cca, k) in reps
                 else np.nan for k in KRULE_ORDER]
            ax.plot(x, y, "o-", color=colour, lw=1.8, ms=5, label=cca)
        if ref is not None:
            ax.axhline(ref, color="k", ls="--", lw=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(KRULE_ORDER, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(label, fontsize=9)
        ax.grid(axis="y", alpha=0.25)
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle("Sweep main effects -- how the PC-count rule and CCA type "
                 "drive the result\n(lag-scan window omitted: 5/10/20 give "
                 "identical values for every metric shown)", fontsize=11)
    _save(fig, "sweep_main_effects")


# --------------------------------------------------------------------------
# figure 3 -- every distinct config ranked by total significant pairs
# --------------------------------------------------------------------------
def fig_config_ranking(reps):
    items = []
    for (cca, k), entry in reps.items():
        a = config_aggregate(entry)
        items.append((short_label(cca, k), cca, a["n_sig_strength"],
                      a["n_sig_ifi"], a["n_sig_total"],
                      entry["tag"] == COMMITTED))
    items.sort(key=lambda t: t[4])
    labels = [it[0] for it in items]
    y = np.arange(len(items))
    strength = np.array([it[2] for it in items])
    ifi = np.array([it[3] for it in items])

    fig, ax = plt.subplots(figsize=(9.5, 8.2))
    ax.barh(y, strength, color="#7b9fd0", label="strength change (of 8 pairs)")
    ax.barh(y, ifi, left=strength, color="#e0a05a",
            label="IFI |lag|<=5 (of 8 pairs)")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    for i, it in enumerate(items):
        if it[5]:
            ax.get_yticklabels()[i].set_color("#c0392b")
            ax.get_yticklabels()[i].set_fontweight("bold")
            ax.text(it[4] + 0.15, i, "<- committed", color="#c0392b",
                    va="center", fontsize=8)
    ax.set_xlabel("number of area-pairs with a significant effect")
    ax.set_ylabel("config  (lag=10 representative; residual vs signal)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("Sweep config ranking -- where the significant effects "
                 "concentrate\n(committed config in red)", fontsize=10)
    _save(fig, "sweep_config_ranking")


# --------------------------------------------------------------------------
# figure 4 -- effect size vs significance, one point per config
# --------------------------------------------------------------------------
def fig_config_scatter(reps):
    fig, ax = plt.subplots(figsize=(9.5, 7))
    pts = []
    for (cca, k), entry in reps.items():
        a = config_aggregate(entry)
        pts.append((a["mean_dcc"], a["n_sig_total"], cca,
                    short_label(cca, k), entry["tag"] == COMMITTED))
    for cca, colour in (("residual", RES_COLOUR), ("signal", SIG_COLOUR)):
        xs = [p[0] for p in pts if p[2] == cca]
        ys = [p[1] for p in pts if p[2] == cca]
        ax.scatter(xs, ys, s=70, color=colour, alpha=0.75, label=cca,
                   edgecolor="white", linewidth=0.6, zorder=3)
    ax.axvline(0, color="k", ls="--", lw=0.7)
    # Many configs collapse to the same (dCC, #sig) point -- group by position
    # so labels are not drawn on top of each other.
    groups: dict[tuple, list] = {}
    for x, y, _cca, lab, committed in pts:
        if not np.isfinite(x):
            continue
        groups.setdefault((round(x, 4), round(y)), []).append((lab, committed))
    ranked = sorted(groups, key=lambda k: -k[1])          # by #significant
    top = set(ranked[:4])
    for (x, y), members in groups.items():
        committed_here = any(c for _, c in members)
        if not (committed_here or (x, y) in top):
            continue
        names = [lab for lab, _ in members]
        txt = names[0] + (f" +{len(names) - 1}" if len(names) > 1 else "")
        if committed_here:
            txt = f"committed ({txt})"
        ax.annotate(txt, (x, y), textcoords="offset points", xytext=(8, 5),
                    fontsize=8, color="#c0392b" if committed_here else "0.3",
                    fontweight="bold" if committed_here else "normal")
    ax.set_xlabel("mean dCC across the 8 pairs  (expert - naive)")
    ax.set_ylabel("# pairs with a significant effect  (strength + IFI, of 16)")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    ax.set_title("Sweep configs -- effect size vs significance "
                 "(top-right = strongest, most significant)", fontsize=10)
    _save(fig, "sweep_config_scatter")


# --------------------------------------------------------------------------
# figure 5 -- where the committed config sits in the sweep distribution
# --------------------------------------------------------------------------
def fig_committed_context(reps):
    metrics = [("mean_cc", "mean held-out CC"),
               ("mean_dcc", "mean dCC (expert - naive)"),
               ("n_sig_strength", "# pairs strength p<.05"),
               ("n_sig_ifi", "# pairs IFI p<.05"),
               ("mean_amf", "mean angle - split-half floor"),
               ("mean_gini", "mean Gini (expert)")]
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.6))
    rng = np.random.default_rng(0)
    for ax, (key, label) in zip(axes.ravel(), metrics):
        vals, comm = [], np.nan
        for (cca, k), entry in reps.items():
            v = config_aggregate(entry)[key]
            vals.append(v)
            if entry["tag"] == COMMITTED:
                comm = v
        vals = np.array(vals, float)
        finite = vals[np.isfinite(vals)]
        jit = (rng.random(vals.size) - 0.5) * 0.3
        ax.scatter(vals, jit, s=32, color="0.6", alpha=0.7, zorder=2)
        if finite.size:
            ax.axvline(np.median(finite), color="k", ls="--", lw=0.9,
                       label="sweep median")
        if np.isfinite(comm):
            ax.scatter([comm], [0], s=130, color="#c0392b", marker="D",
                       zorder=4, label="committed")
            pct = 100.0 * np.mean(finite <= comm) if finite.size else np.nan
            ax.set_title(f"{label}\ncommitted at {pct:.0f}th percentile",
                         fontsize=8.5)
        else:
            ax.set_title(label, fontsize=8.5)
        ax.set_yticks([])
        ax.set_ylim(-0.4, 0.4)
    axes[0, 0].legend(frameon=False, fontsize=7, loc="upper right")
    fig.suptitle("Is the committed config representative? -- committed value "
                 "(red) vs the spread of all 22 distinct configs", fontsize=11)
    _save(fig, "sweep_committed_context")


# --------------------------------------------------------------------------
# figure 6 -- per-pair structure across the sweep
# --------------------------------------------------------------------------
def fig_pair_structure(reps):
    ordered = [reps[(c, k)] for c in ("residual", "signal")
               for k in KRULE_ORDER if (c, k) in reps]
    col_labels = [short_label(e["cca"], e["k_rule"]) for e in ordered]
    dcc = np.array([[_f(e["pairs"][p].get("d_cc")) for e in ordered]
                    for p in PAIRS], float)

    fig, axes = plt.subplots(2, 1, figsize=(11, 9),
                             gridspec_kw={"height_ratios": [1, 1.05]})

    # panel A -- dCC heatmap, pairs x configs
    vmax = np.nanmax(np.abs(dcc))
    im = axes[0].imshow(dcc, aspect="auto", cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax)
    axes[0].set_xticks(np.arange(len(ordered)))
    axes[0].set_xticklabels(col_labels, rotation=90, fontsize=7)
    axes[0].set_yticks(np.arange(len(PAIRS)))
    axes[0].set_yticklabels(PAIRS, fontsize=8)
    axes[0].axvline(len(ordered) / 2 - 0.5, color="k", lw=1.2)
    fig.colorbar(im, ax=axes[0], label="dCC (expert - naive)", pad=0.01)
    axes[0].set_title("dCC by area-pair across every distinct config "
                      "(red = strengthens with learning)", fontsize=10)

    # panel B -- per-pair robustness: fraction of configs passing each test
    tests = [("d_cc", "dCC > 0", lambda m: _f(m.get("d_cc")) > 0),
             ("p_str", "strength p<.05",
              lambda m: _f(m.get("p_naive_vs_expert")) < ALPHA),
             ("angle", "angle > floor",
              lambda m: _f(m.get("angle_minus_floor")) > 0),
             ("ifi", "IFI p<.05",
              lambda m: _f(m.get("p_ifi_w5")) < ALPHA)]
    n_cfg = len(ordered)
    x = np.arange(len(PAIRS))
    width = 0.2
    for i, (_k, label, fn) in enumerate(tests):
        frac = [np.mean([fn(e["pairs"][p]) for e in ordered]) for p in PAIRS]
        axes[1].bar(x + (i - 1.5) * width, frac, width, label=label)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(PAIRS, fontsize=8)
    axes[1].set_ylabel(f"fraction of the {n_cfg} distinct configs")
    axes[1].set_ylim(0, 1)
    axes[1].axhline(0.5, color="k", ls=":", lw=0.7)
    axes[1].legend(frameon=False, fontsize=8, ncol=4, loc="upper center")
    axes[1].set_title("Per-pair robustness -- how often each effect holds "
                      "across the sweep", fontsize=10)
    _save(fig, "sweep_pair_structure")


# --------------------------------------------------------------------------
def main():
    rows = load_rows()
    tags = by_tag(rows)
    for t, e in tags.items():
        e["tag"] = t
    # lag=10 representative of each (CCA type x PC-rule) cell -- 22 configs.
    reps = {(e["cca"], e["k_rule"]): e
            for e in tags.values() if e["max_lag"] == 10}

    n_distinct = fig_redundancy(tags)
    fig_main_effects(reps)
    fig_config_ranking(reps)
    fig_config_scatter(reps)
    fig_committed_context(reps)
    fig_pair_structure(reps)
    print(f"\nSweep exploration done -- {len(tags)} configs, "
          f"{len(reps)} distinct (CCA x PC-rule) cells, "
          f"{n_distinct} distinct headline results.")


if __name__ == "__main__":
    main()
