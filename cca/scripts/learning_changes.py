"""Paired naive -> expert tests on the landmark-arm cohort result.

For each ``(animal, pair, landmark)`` triple we compute the change from naive
to expert in two cell-level subspace-aggregated stats:

* ``delta_mi_sig``  -- Gaussian-MI surrogate summed over significant dims,
  expert minus naive. Communication-subspace strength change.
* ``delta_ifi_weighted`` -- CC-weighted mean IFI (high-CC dims dominate),
  expert minus naive. Directionality change.

Two test families:

1. **Per (pair, landmark)** -- Wilcoxon signed-rank on the per-animal delta.
   n = number of animals with both naive and expert cells. 8 pairs * 6
   landmarks = 48 tests, FDR(BH)-corrected within each contrast and stat.
2. **Per pair, pooling landmarks** -- same test on the pooled per-animal-per-
   landmark deltas. 8 tests per stat, FDR-corrected.

Writes:
* ``results/learning_changes_<tag>.csv`` -- long-form per-test table.
* ``figures/<tag>/learning_changes_perlandmark.png`` -- per-(pair, landmark)
  heatmap with FDR significance markers.
* ``figures/<tag>/learning_changes_pooled.png`` -- pooled per-pair bars.

Operates entirely on an existing landmark pkl -- no re-run needed.
"""

from __future__ import annotations

import argparse
import csv
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, subspace_stats  # noqa: E402


LANDMARKS = list(range(1, 7))
CONTRASTS = (("expert", "naive"), ("expert", "intermediate"),
             ("intermediate", "naive"))
STATS = ("mi_sig", "ifi_weighted")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tag", default="landmark50_res_samp15")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--fdr-q", type=float, default=0.05)
    return p.parse_args()


def _load(tag: str):
    pkl = config.RESULTS_DIR / f"{tag}.pkl"
    if not pkl.is_file():
        raise SystemExit(f"missing {pkl}")
    with open(pkl, "rb") as f:
        raw = pickle.load(f)
    fits = [r for r in raw["results"] if hasattr(r, "cells")]
    return fits, raw["cfg"]


def _ordered_pairs(fits) -> list[tuple[str, str]]:
    have = {(r.area_x, r.area_y) for r in fits}
    return [p for p in config.PAIRS if p in have]


def _per_cell_table(fits, alpha) -> dict:
    """(animal, pair, epoch, landmark) -> {stat: value}."""
    out: dict = {}
    for r in fits:
        for (epoch, lm), cell in r.cells.items():
            s = subspace_stats.cell_subspace_stats(cell, alpha=alpha)
            if s is None:
                continue
            out[(r.animal_id, (r.area_x, r.area_y), epoch, int(lm))] = {
                "mi_sig": s.mi_sig,
                "mi_all": s.mi_all,
                "n_sig": s.n_sig,
                "ifi_weighted": s.ifi_weighted,
                "ifi_mean_sig": s.ifi_mean_sig,
            }
    return out


def _paired_deltas(table, pair, lm, ep_a, ep_b, stat) -> list[float]:
    """Return per-animal (a - b) deltas for the given pair/landmark/contrast."""
    animals = sorted({k[0] for k in table if k[1] == pair and k[3] == lm})
    deltas = []
    for a in animals:
        ka = (a, pair, ep_a, lm)
        kb = (a, pair, ep_b, lm)
        if ka in table and kb in table:
            va, vb = table[ka][stat], table[kb][stat]
            if np.isfinite(va) and np.isfinite(vb):
                deltas.append(va - vb)
    return deltas


def _wilcoxon(deltas: list[float]):
    """Two-sided Wilcoxon signed-rank; returns (n, median, statistic, p)."""
    arr = np.asarray(deltas, dtype=float)
    arr = arr[np.isfinite(arr)]
    n = arr.size
    if n < 3 or np.all(arr == 0):
        return (n, float(np.nanmedian(arr) if n else np.nan),
                float("nan"), float("nan"))
    try:
        res = stats.wilcoxon(arr, zero_method="wilcox",
                             alternative="two-sided",
                             nan_policy="omit", method="auto")
        return (n, float(np.median(arr)), float(res.statistic), float(res.pvalue))
    except ValueError:                                  # all-zero / degenerate
        return (n, float(np.median(arr)), float("nan"), float("nan"))


def _fdr_bh(pvals: np.ndarray, q: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg FDR mask. Operates on finite p-values only."""
    p = np.asarray(pvals, dtype=float)
    finite = np.isfinite(p)
    if not np.any(finite):
        return np.zeros_like(p, dtype=bool)
    p_f = p[finite]
    n = p_f.size
    order = np.argsort(p_f)
    ranked = p_f[order]
    threshold = q * np.arange(1, n + 1) / n
    passed = ranked <= threshold
    cutoff = ranked[np.max(np.where(passed))] if np.any(passed) else -1
    out_f = p_f <= cutoff if cutoff >= 0 else np.zeros(n, dtype=bool)
    out = np.zeros_like(p, dtype=bool)
    out[finite] = out_f
    return out


def _run_tests(table, pairs, alpha, fdr_q) -> list[dict]:
    """Build per (pair, landmark, contrast, stat) test rows."""
    rows: list[dict] = []
    # First pass: per (pair, landmark, contrast, stat).
    for stat in STATS:
        for ep_a, ep_b in CONTRASTS:
            family_idx = []                             # row indices for FDR
            for pair in pairs:
                for lm in LANDMARKS:
                    deltas = _paired_deltas(table, pair, lm, ep_a, ep_b, stat)
                    n, med, w, p = _wilcoxon(deltas)
                    rows.append({
                        "scope": "per_landmark",
                        "pair": f"{pair[0]}-{pair[1]}",
                        "landmark": lm,
                        "stat": stat,
                        "contrast": f"{ep_a}-{ep_b}",
                        "n": n,
                        "median_delta": med,
                        "wilcoxon_W": w,
                        "p": p,
                    })
                    family_idx.append(len(rows) - 1)
            p_arr = np.array([rows[i]["p"] for i in family_idx])
            passed = _fdr_bh(p_arr, q=fdr_q)
            for i, pass_ in zip(family_idx, passed):
                rows[i]["p_fdr_pass"] = bool(pass_)
    # Second pass: per (pair, contrast, stat) pooling landmarks.
    for stat in STATS:
        for ep_a, ep_b in CONTRASTS:
            family_idx = []
            for pair in pairs:
                deltas = []
                for lm in LANDMARKS:
                    deltas.extend(_paired_deltas(table, pair, lm,
                                                 ep_a, ep_b, stat))
                n, med, w, p = _wilcoxon(deltas)
                rows.append({
                    "scope": "per_pair_pooled",
                    "pair": f"{pair[0]}-{pair[1]}",
                    "landmark": "",
                    "stat": stat,
                    "contrast": f"{ep_a}-{ep_b}",
                    "n": n,
                    "median_delta": med,
                    "wilcoxon_W": w,
                    "p": p,
                })
                family_idx.append(len(rows) - 1)
            p_arr = np.array([rows[i]["p"] for i in family_idx])
            passed = _fdr_bh(p_arr, q=fdr_q)
            for i, pass_ in zip(family_idx, passed):
                rows[i]["p_fdr_pass"] = bool(pass_)
    return rows


def _write_csv(rows: list[dict], path: Path):
    fieldnames = ["scope", "pair", "landmark", "stat", "contrast",
                  "n", "median_delta", "wilcoxon_W", "p", "p_fdr_pass"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def _plot_per_landmark(rows, pairs, out_path, stat_label, stat):
    contrast = "expert-naive"
    pair_labels = [f"{p[0]}-{p[1]}" for p in pairs]
    grid_delta = np.full((len(pairs), len(LANDMARKS)), np.nan)
    grid_passed = np.zeros((len(pairs), len(LANDMARKS)), dtype=bool)
    grid_n = np.zeros((len(pairs), len(LANDMARKS)), dtype=int)
    by_key = {(r["pair"], r["landmark"], r["stat"], r["contrast"]): r
              for r in rows if r["scope"] == "per_landmark"}
    for i, pair in enumerate(pairs):
        for j, lm in enumerate(LANDMARKS):
            r = by_key.get((f"{pair[0]}-{pair[1]}", lm, stat, contrast))
            if r is None:
                continue
            grid_delta[i, j] = r["median_delta"]
            grid_passed[i, j] = r.get("p_fdr_pass", False)
            grid_n[i, j] = r["n"]
    vmax = max(0.01, np.nanmax(np.abs(grid_delta)) if np.isfinite(grid_delta).any() else 0.01)
    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(grid_delta, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(len(pairs)):
        for j in range(len(LANDMARKS)):
            if grid_n[i, j] == 0:
                ax.text(j, i, "·", ha="center", va="center", color="grey")
            else:
                ann = (f"{grid_delta[i, j]:+.2f}"
                       + ("*" if grid_passed[i, j] else ""))
                ax.text(j, i, ann, ha="center", va="center", fontsize=8,
                        color="white" if abs(grid_delta[i, j]) > 0.6 * vmax
                        else "black")
    ax.set_xticks(range(len(LANDMARKS)), [f"L{l}" for l in LANDMARKS])
    ax.set_yticks(range(len(pairs)), pair_labels)
    ax.set_xlabel("Landmark id")
    ax.set_title(f"Per-(pair, landmark)  expert - naive  median Δ {stat_label}  "
                 f"(* = FDR < 0.05; n in cell, dot = no animals paired)")
    fig.colorbar(im, ax=ax, shrink=0.7, label=f"Δ {stat_label}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_pooled(rows, pairs, out_path):
    contrast = "expert-naive"
    pair_labels = [f"{p[0]}-{p[1]}" for p in pairs]
    by_key = {(r["pair"], r["stat"], r["contrast"]): r
              for r in rows if r["scope"] == "per_pair_pooled"}
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    for ax, stat, lab in zip(axes, STATS,
                             ("Δ MI_sig (nats)", "Δ IFI_weighted")):
        meds = []
        passed_mask = []
        ns = []
        for pair in pairs:
            r = by_key.get((f"{pair[0]}-{pair[1]}", stat, contrast))
            meds.append(r["median_delta"] if r else np.nan)
            passed_mask.append(r.get("p_fdr_pass", False) if r else False)
            ns.append(r["n"] if r else 0)
        xs = np.arange(len(pairs))
        bars = ax.bar(xs, meds,
                      color=["firebrick" if p else "steelblue"
                             for p in passed_mask])
        ax.axhline(0, color="grey", lw=0.8, ls=":")
        ax.set_xticks(xs, pair_labels, rotation=30)
        for x, n in zip(xs, ns):
            ax.text(x, ax.get_ylim()[1] * 0.95, f"n={n}",
                    ha="center", va="top", fontsize=8, color="grey")
        ax.set_ylabel(lab)
        ax.set_title(f"Pooled across landmarks  ({lab.split(' ')[1]} contrast: "
                     f"expert - naive)")
    fig.suptitle("Pooled per-pair learning change  "
                 "(red = FDR < 0.05 across the 8-pair family)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    args = parse_args()
    fits, cfg = _load(args.tag)
    table = _per_cell_table(fits, args.alpha)
    pairs = _ordered_pairs(fits)
    rows = _run_tests(table, pairs, alpha=args.alpha, fdr_q=args.fdr_q)

    fig_dir = config.FIGURES_DIR / args.tag
    fig_dir.mkdir(parents=True, exist_ok=True)
    csv_path = config.RESULTS_DIR / f"learning_changes_{args.tag}.csv"
    _write_csv(rows, csv_path)
    print(f"wrote {csv_path}")

    out1 = _plot_per_landmark(rows, pairs,
                              fig_dir / "learning_changes_perlandmark_mi.png",
                              "MI_sig (nats)", "mi_sig")
    out2 = _plot_per_landmark(rows, pairs,
                              fig_dir / "learning_changes_perlandmark_ifi.png",
                              "IFI_weighted", "ifi_weighted")
    out3 = _plot_pooled(rows, pairs,
                        fig_dir / "learning_changes_pooled.png")
    for p in (out1, out2, out3):
        print(f"wrote {p}")

    # Short text summary -- pooled survivors (FDR<0.05) per stat.
    print("\nPooled per-pair survivors (FDR < q):")
    for stat in STATS:
        for ep_a, ep_b in CONTRASTS:
            survs = [r for r in rows
                     if r["scope"] == "per_pair_pooled"
                     and r["stat"] == stat
                     and r["contrast"] == f"{ep_a}-{ep_b}"
                     and r.get("p_fdr_pass")]
            if survs:
                print(f"  {stat}  {ep_a}-{ep_b}:")
                for s in survs:
                    print(f"    {s['pair']}: n={s['n']}, "
                          f"med Δ={s['median_delta']:+.3f}, p={s['p']:.4f}")


if __name__ == "__main__":
    main()
