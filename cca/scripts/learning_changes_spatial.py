"""Paired naive -> expert tests on the SPATIAL-arm cohort result.

The spatial analogue of ``learning_changes.py``. The spatial arm fits one cell
per ``(animal, pair, epoch)`` over the full 200-bin corridor -- there is no
landmark dimension -- so there is a single test family: **per (pair, contrast)**
Wilcoxon signed-rank on the per-animal delta, FDR(BH)-corrected within each
(stat, contrast) family across the 8 pairs. This mirrors the landmark arm's
``per_pair_pooled`` scope.

Two cell-level stats (see ``subspace_stats.epoch_subspace_stats``):
* ``mi_sig``        -- Gaussian-MI surrogate over significant dims (strength).
* ``ifi_weighted``  -- CC-weighted mean IFI (directionality).

Both significance and strength are taken at lag 0 (the spatial arm's native D7
test), unlike the landmark arm's peak-across-lag readout -- so the spatial
numbers are comparable in form but not identical in definition.

Writes:
* ``results/learning_changes_spatial_<tag>.csv`` -- long-form per-test table.
* ``figures/spatial_learning/<tag>_pooled.png``  -- per-pair bars, expert-naive.

Operates entirely on an existing spatial stage-2 pkl -- no re-run needed.
"""

from __future__ import annotations

import argparse
import csv
import pickle
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, paired_stats, subspace_stats  # noqa: E402


CONTRASTS = (("expert", "naive"), ("expert", "intermediate"),
             ("intermediate", "naive"))
STATS = ("mi_sig", "ifi_weighted")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tag", default="res_samp15_lag10",
                   help="Spatial stage-2 config tag (no 'stage2_' prefix).")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--fdr-q", type=float, default=0.05)
    return p.parse_args()


def _load(tag: str):
    pkl = config.RESULTS_DIR / f"stage2_{tag}.pkl"
    if not pkl.is_file():
        raise SystemExit(f"missing {pkl}")
    with open(pkl, "rb") as f:
        raw = pickle.load(f)
    fits = [r for r in raw["results"] if hasattr(r, "epochs")]
    return fits, raw["cfg"]


def _ordered_pairs(fits) -> list[tuple[str, str]]:
    have = {(r.area_x, r.area_y) for r in fits}
    return [p for p in config.PAIRS if p in have]


def _per_cell_table(fits, alpha) -> dict:
    """(animal, pair, epoch) -> {stat: value}."""
    out: dict = {}
    for r in fits:
        for epoch, ep in r.epochs.items():
            s = subspace_stats.epoch_subspace_stats(ep, alpha=alpha)
            if s is None:
                continue
            out[(r.animal_id, (r.area_x, r.area_y), epoch)] = {
                "mi_sig": s.mi_sig,
                "ifi_weighted": s.ifi_weighted,
                "n_sig": s.n_sig,
            }
    return out


def _paired_deltas(table, pair, ep_a, ep_b, stat) -> list[float]:
    """Per-animal (a - b) deltas for the given pair/contrast."""
    animals = sorted({k[0] for k in table if k[1] == pair})
    deltas = []
    for a in animals:
        ka, kb = (a, pair, ep_a), (a, pair, ep_b)
        if ka in table and kb in table:
            va, vb = table[ka][stat], table[kb][stat]
            if np.isfinite(va) and np.isfinite(vb):
                deltas.append(va - vb)
    return deltas


def _run_tests(table, pairs, fdr_q) -> list[dict]:
    """Per (pair, contrast, stat); FDR within each (stat, contrast) family."""
    rows: list[dict] = []
    for stat in STATS:
        for ep_a, ep_b in CONTRASTS:
            family: list[dict] = []
            for pair in pairs:
                deltas = _paired_deltas(table, pair, ep_a, ep_b, stat)
                n, med, w, p = paired_stats.wilcoxon_signed(deltas)
                family.append({
                    "scope": "per_pair",
                    "pair": f"{pair[0]}-{pair[1]}",
                    "stat": stat,
                    "contrast": f"{ep_a}-{ep_b}",
                    "n": n,
                    "median_delta": med,
                    "wilcoxon_W": w,
                    "p": p,
                    "p_fdr_pass": False,
                })
            passed = paired_stats.fdr_bh(
                np.array([r["p"] for r in family]), q=fdr_q)
            for r, ok in zip(family, passed):
                r["p_fdr_pass"] = bool(ok)
            rows.extend(family)
    return rows


def _write_csv(rows: list[dict], path: Path):
    fieldnames = ["scope", "pair", "stat", "contrast", "n",
                  "median_delta", "wilcoxon_W", "p", "p_fdr_pass"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def _plot_pooled(rows, pairs, out_path, tag):
    contrast = "expert-naive"
    pair_labels = [f"{p[0]}-{p[1]}" for p in pairs]
    fig, axes = plt.subplots(1, len(STATS), figsize=(5 * len(STATS), 4.2),
                             squeeze=False)
    for ax, stat in zip(axes[0], STATS):
        by_pair = {r["pair"]: r for r in rows
                   if r["stat"] == stat and r["contrast"] == contrast}
        vals, colours = [], []
        for lab in pair_labels:
            r = by_pair.get(lab)
            v = r["median_delta"] if r and np.isfinite(r["median_delta"]) else 0.0
            vals.append(v)
            colours.append("#c0392b" if (r and r["p_fdr_pass"]) else "#95a5a6")
        ax.bar(range(len(pair_labels)), vals, color=colours)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(range(len(pair_labels)))
        ax.set_xticklabels(pair_labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(f"median Δ {stat} (expert - naive)")
        ax.set_title(stat)
    fig.suptitle(f"Spatial arm {tag}: naive→expert change "
                 "(red = FDR-significant, per-pair paired Wilcoxon)")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    fits, _cfg = _load(args.tag)
    pairs = _ordered_pairs(fits)
    table = _per_cell_table(fits, args.alpha)
    rows = _run_tests(table, pairs, fdr_q=args.fdr_q)

    csv_path = config.RESULTS_DIR / f"learning_changes_spatial_{args.tag}.csv"
    _write_csv(rows, csv_path)
    print(f"wrote {csv_path}  ({len(rows)} rows)")

    fig_path = config.FIGURES_DIR / "spatial_learning" / f"{args.tag}_pooled.png"
    _plot_pooled(rows, pairs, fig_path, args.tag)
    print(f"wrote {fig_path}")

    # Console summary: FDR survivors.
    survivors = [r for r in rows if r["p_fdr_pass"]]
    if survivors:
        print("FDR-significant (pair / stat / contrast / median_delta / p):")
        for r in sorted(survivors, key=lambda x: x["p"]):
            print(f"  {r['pair']:9s} {r['stat']:12s} {r['contrast']:20s} "
                  f"d={r['median_delta']:+.4f}  p={r['p']:.4g}")
    else:
        print("No FDR-significant naive->expert changes in the spatial arm "
              f"for {args.tag}.")


if __name__ == "__main__":
    main()
