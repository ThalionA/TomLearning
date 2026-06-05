"""Sweep summary for the landmark arm -- works on partial completion.

Scans ``results/`` for any ``landmark{bin}_{cca}_{krule}.pkl`` files that have
a matching ``.done`` marker. Computes per-(config, pair) summary stats in the
striatum-style convention (statistical unit = significant subspace dimension,
pooled across animals and landmarks within each pair) and writes:

  results/sweep_landmark_summary.csv    long-form per-(config, pair) table.
  figures/landmark_sweep/
    sweep_strength_heatmap.png  median held-out peak CC across sig dims,
                                configs as rows, pairs as columns.
    sweep_n_sig_heatmap.png     mean n_sig dims per cell.
    sweep_ifi_heatmap.png       median IFI across sig dims (full lag window).
    sweep_frac_sig_cells.png    fraction of cells with >=1 significant dim.
    sweep_lr_change.png         per-pair median (expert - naive) of CC and IFI
                                across configs (the "robustness of the
                                learning effect" view).

The `.done` marker filter excludes mid-run pkls that may still be growing.

Run:  python scripts/summarise_landmark_sweep.py
"""

from __future__ import annotations

import argparse
import csv
import pickle
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, subspace_stats  # noqa: E402


def _wilcoxon_vs0(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 6 or not np.any(v != 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(v).pvalue)
    except ValueError:
        return float("nan")


def _mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return float("nan")
    try:
        return float(stats.mannwhitneyu(a, b).pvalue)
    except ValueError:
        return float("nan")


def _stars(p, thresholds=(0.001, 0.01, 0.05)):
    if not np.isfinite(p):
        return ""
    for n, t in zip(("***", "**", "*"), thresholds):
        if p < t:
            return n
    return ""


EPOCHS = config.EPOCH_NAMES
ALPHA = 0.05
TAG_RE = re.compile(r"^landmark(\d+)_(res|sig)_(\w+)$")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--alpha", type=float, default=0.05)
    return p.parse_args()


def _find_done_configs() -> list[tuple[str, Path]]:
    """Tags whose .pkl has a corresponding .done marker, sorted naturally."""
    out = []
    for done in sorted(config.RESULTS_DIR.glob("landmark*.done")):
        tag = done.stem
        pkl = config.RESULTS_DIR / f"{tag}.pkl"
        if pkl.is_file():
            out.append((tag, pkl))
    return out


def _parse_tag(tag: str) -> tuple[int, str, str] | None:
    m = TAG_RE.match(tag)
    if not m:
        return None
    return int(m.group(1)), m.group(2), m.group(3)


def _per_pair_stats(fits, alpha: float) -> dict:
    """Per-pair, per-epoch summaries for one config.

    Stores pooled sig-dim values per epoch:
      cc:       held-out peak CC across sig dims (one value per sig dim)
      ifi:      full-window IFI (one value per sig dim)
      ifi_win:  dict {window -> list of IFI values} pooled across sig dims
      n_sig:    per-cell count of sig dims (one value per cell)
    """
    def _new_bucket():
        return {e: {"cc": [], "ifi": [], "n_sig_cells": [],
                    "ifi_win": defaultdict(list),
                    "n_cells_total": 0, "n_cells_with_sig": 0}
                for e in EPOCHS}

    by_pair: dict = defaultdict(_new_bucket)
    for r in fits:
        if not hasattr(r, "cells"):
            continue
        pair = (r.area_x, r.area_y)
        for (epoch, lm), cell in r.cells.items():
            if epoch not in EPOCHS:
                continue
            bucket = by_pair[pair][epoch]
            if cell.skip_reason is not None:
                continue
            s = subspace_stats.cell_subspace_stats(cell, alpha=alpha)
            if s is None:
                continue
            bucket["n_cells_total"] += 1
            bucket["n_sig_cells"].append(int(s.n_sig))
            if s.n_sig > 0:
                bucket["n_cells_with_sig"] += 1
                bucket["cc"].extend(s.peak_cc_per_dim[s.sig_mask].tolist())
                ifi_sig = cell.ifi_per_dim[s.sig_mask]
                bucket["ifi"].extend([float(x) for x in ifi_sig
                                      if np.isfinite(x)])
                # Per-window IFI: cell.ifi_windows has shape (n_dims, max_w).
                if cell.ifi_windows is not None and cell.ifi_windows.ndim == 2:
                    max_w = cell.ifi_windows.shape[1]
                    for w in range(1, max_w + 1):
                        vals = cell.ifi_windows[s.sig_mask, w - 1]
                        bucket["ifi_win"][w].extend(
                            [float(x) for x in vals if np.isfinite(x)])
    return by_pair


def _aggregate_rows(tag, pair, by_epoch) -> list[dict]:
    bin_ms, cca, krule = _parse_tag(tag)
    rows = []
    for epoch in EPOCHS:
        b = by_epoch[epoch]
        cc = np.array(b["cc"])
        ifi = np.array(b["ifi"])
        n_sig = np.array(b["n_sig_cells"])
        rows.append({
            "tag": tag,
            "bin_ms": bin_ms,
            "cca": cca,
            "krule": krule,
            "pair": f"{pair[0]}-{pair[1]}",
            "epoch": epoch,
            "n_cells_total": b["n_cells_total"],
            "n_cells_with_sig": b["n_cells_with_sig"],
            "frac_cells_sig": (b["n_cells_with_sig"] / b["n_cells_total"]
                               if b["n_cells_total"] > 0 else float("nan")),
            "mean_n_sig": float(np.mean(n_sig)) if n_sig.size else float("nan"),
            "median_cc_sig": float(np.median(cc)) if cc.size else float("nan"),
            "median_ifi_sig": float(np.median(ifi)) if ifi.size else float("nan"),
            "n_sig_dims_total": int(cc.size),
        })
    return rows


def _heatmap(values: np.ndarray, tags: list[str], pairs: list[str],
             title: str, label: str, cmap: str, vmin: float, vmax: float,
             centre_zero: bool = False, out_path: Path = None,
             cell_fmt: str = "{:.2f}"):
    fig, ax = plt.subplots(figsize=(max(8, 1.0 * len(pairs)),
                                    max(5, 0.35 * len(tags))))
    if centre_zero:
        m = max(abs(vmin), abs(vmax))
        im = ax.imshow(values, cmap=cmap, vmin=-m, vmax=m, aspect="auto")
    else:
        im = ax.imshow(values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            if np.isfinite(v):
                ax.text(j, i, cell_fmt.format(v), ha="center", va="center",
                        fontsize=7,
                        color="white" if (centre_zero and abs(v) > 0.55 * m)
                        or (not centre_zero and v > 0.6 * vmax)
                        else "black")
    ax.set_xticks(range(len(pairs)), pairs, rotation=30, ha="right")
    ax.set_yticks(range(len(tags)), tags, fontsize=8)
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.7, label=label)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    args = parse_args()
    configs = _find_done_configs()
    if not configs:
        raise SystemExit("no landmark*.done files in results/")
    print(f"found {len(configs)} done configs")

    rows: list[dict] = []
    for tag, pkl in configs:
        if _parse_tag(tag) is None:
            print(f"  skip (unparseable tag): {tag}")
            continue
        with open(pkl, "rb") as f:
            raw = pickle.load(f)
        fits = [r for r in raw["results"] if hasattr(r, "cells")]
        by_pair = _per_pair_stats(fits, alpha=args.alpha)
        for pair, by_epoch in sorted(by_pair.items(),
                                     key=lambda kv: config.PAIRS.index(kv[0])
                                     if kv[0] in config.PAIRS else 99):
            rows.extend(_aggregate_rows(tag, pair, by_epoch))
        print(f"  {tag}: {len(fits)} fits, {len(by_pair)} pairs")

    out_csv = config.RESULTS_DIR / "sweep_landmark_summary.csv"
    fieldnames = ["tag", "bin_ms", "cca", "krule", "pair", "epoch",
                  "n_cells_total", "n_cells_with_sig", "frac_cells_sig",
                  "mean_n_sig", "median_cc_sig", "median_ifi_sig",
                  "n_sig_dims_total"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_csv}  ({len(rows)} rows)")

    fig_dir = config.FIGURES_DIR / "landmark_sweep"
    fig_dir.mkdir(parents=True, exist_ok=True)

    tags = sorted({r["tag"] for r in rows})
    pairs = [f"{p[0]}-{p[1]}" for p in config.PAIRS
             if any(r["pair"] == f"{p[0]}-{p[1]}" for r in rows)]

    def _grid(metric: str, epoch: str) -> np.ndarray:
        g = np.full((len(tags), len(pairs)), np.nan)
        for r in rows:
            if r["epoch"] != epoch:
                continue
            i = tags.index(r["tag"])
            j = pairs.index(r["pair"])
            g[i, j] = r[metric]
        return g

    # Strength: median CC across sig dims, naive epoch (one heatmap per epoch).
    for epoch in EPOCHS:
        out = fig_dir / f"sweep_strength_{epoch}.png"
        _heatmap(_grid("median_cc_sig", epoch), tags, pairs,
                 title=f"Median held-out peak CC across significant dims  "
                       f"({epoch} epoch; pooled animals × landmarks)",
                 label="median CC", cmap="viridis", vmin=0, vmax=1,
                 out_path=out)
        print(f"wrote {out}")

    # n_sig per cell.
    for epoch in EPOCHS:
        out = fig_dir / f"sweep_n_sig_{epoch}.png"
        _heatmap(_grid("mean_n_sig", epoch), tags, pairs,
                 title=f"Mean # significant subspace dims per cell  "
                       f"({epoch} epoch)",
                 label="mean n_sig", cmap="magma", vmin=0,
                 vmax=max(1, float(np.nanmax(_grid("mean_n_sig", epoch)) or 1)),
                 out_path=out, cell_fmt="{:.1f}")
        print(f"wrote {out}")

    # IFI (expert epoch -- the one most likely to show learning effect).
    out = fig_dir / "sweep_ifi_expert.png"
    g = _grid("median_ifi_sig", "expert")
    vmax = float(np.nanmax(np.abs(g)) or 0.5)
    _heatmap(g, tags, pairs,
             title="Median IFI across significant dims  (expert epoch)",
             label="median IFI", cmap="RdBu_r",
             vmin=-vmax, vmax=vmax, centre_zero=True, out_path=out,
             cell_fmt="{:+.2f}")
    print(f"wrote {out}")

    # Learning change: median CC (expert - naive) per (config, pair).
    delta_cc = (_grid("median_cc_sig", "expert")
                - _grid("median_cc_sig", "naive"))
    out = fig_dir / "sweep_lr_change_cc.png"
    _heatmap(delta_cc, tags, pairs,
             title="Median CC change: expert - naive  "
                   "(robustness across configs)",
             label="Δ median CC", cmap="RdBu_r",
             vmin=-0.3, vmax=0.3, centre_zero=True, out_path=out,
             cell_fmt="{:+.2f}")
    print(f"wrote {out}")

    delta_dim = (_grid("mean_n_sig", "expert")
                 - _grid("mean_n_sig", "naive"))
    out = fig_dir / "sweep_lr_change_n_sig.png"
    _heatmap(delta_dim, tags, pairs,
             title="Mean n_sig change: expert - naive",
             label="Δ mean n_sig", cmap="RdBu_r",
             vmin=-2, vmax=2, centre_zero=True, out_path=out,
             cell_fmt="{:+.1f}")
    print(f"wrote {out}")

    # Frac cells with at least one sig dim (naive epoch baseline).
    out = fig_dir / "sweep_frac_sig_cells_naive.png"
    _heatmap(_grid("frac_cells_sig", "naive"), tags, pairs,
             title="Fraction of cells with ≥1 significant dim  (naive)",
             label="frac sig cells", cmap="viridis",
             vmin=0, vmax=1, out_path=out)
    print(f"wrote {out}")

    # --- New: per-config IFI-window heatmap + p-value tables ----------------
    # Need to re-load to get the per-epoch pooled values for stats.
    print("computing p-values + IFI-window grids...")
    per_config_pair_epoch: dict = {}        # (tag, pair, epoch) -> {cc, ifi, ifi_win}
    for tag, pkl in configs:
        if _parse_tag(tag) is None:
            continue
        with open(pkl, "rb") as f:
            raw = pickle.load(f)
        fits = [r for r in raw["results"] if hasattr(r, "cells")]
        bp = _per_pair_stats(fits, alpha=args.alpha)
        for pair, by_epoch in bp.items():
            for epoch in EPOCHS:
                per_config_pair_epoch[(tag, pair, epoch)] = by_epoch[epoch]

    _plot_ifi_window_per_pair(per_config_pair_epoch, tags, fig_dir)
    _plot_ifi_window_lr_change(per_config_pair_epoch, tags, fig_dir)
    _plot_pvalues_summary(per_config_pair_epoch, tags, fig_dir)


# ---------------------------------------------------------------------------
# New: IFI-window heatmap (per pair) + p-value summary
# ---------------------------------------------------------------------------
def _bin_ms_for_tag(tag):
    pt = _parse_tag(tag)
    return pt[0] if pt else 50


def _ifi_window_grid(per_cpe, tags, pair, epoch, max_ms=500):
    """Build a (configs x window_ms) grid of median IFI and p_vs_0."""
    # x-axis: physical lag window in ms. We sample at every distinct
    # window_ms appearing across the included configs.
    all_ms = sorted({_bin_ms_for_tag(t) * w
                     for t in tags
                     for w in per_cpe.get((t, pair, epoch),
                                          {"ifi_win": {}})["ifi_win"]})
    if not all_ms:
        return np.array([]), np.full((len(tags), 0), np.nan), \
               np.full((len(tags), 0), np.nan)
    med = np.full((len(tags), len(all_ms)), np.nan)
    pv0 = np.full((len(tags), len(all_ms)), np.nan)
    for i, tag in enumerate(tags):
        bin_ms = _bin_ms_for_tag(tag)
        cpe = per_cpe.get((tag, pair, epoch))
        if cpe is None:
            continue
        for w, vals in cpe["ifi_win"].items():
            win_ms = bin_ms * w
            if win_ms not in all_ms:
                continue
            j = all_ms.index(win_ms)
            vals = [v for v in vals if np.isfinite(v)]
            if vals:
                med[i, j] = float(np.median(vals))
                pv0[i, j] = _wilcoxon_vs0(np.array(vals))
    return np.array(all_ms), med, pv0


def _plot_ifi_window_per_pair(per_cpe, tags, fig_dir):
    """2x4 grid of pairs; each panel rows=configs, cols=lag window (ms),
    cell = median IFI in expert epoch, * = p_vs_0 < 0.05."""
    pairs_in_data = [p for p in config.PAIRS
                     if any((t, p, "naive") in per_cpe for t in tags)]
    fig, axes = plt.subplots(2, 4, figsize=(22, 8))
    axes = axes.ravel()
    vmax = 0.0
    grids = []
    for pair in pairs_in_data:
        ms, med, pv0 = _ifi_window_grid(per_cpe, tags, pair, "expert")
        grids.append((pair, ms, med, pv0))
        if med.size:
            vmax = max(vmax, float(np.nanmax(np.abs(med))))
    vmax = max(vmax, 0.1)
    for ax, (pair, ms, med, pv0) in zip(axes, grids):
        if med.size == 0:
            ax.set_visible(False)
            continue
        im = ax.imshow(med, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        for i in range(med.shape[0]):
            for j in range(med.shape[1]):
                v = med[i, j]
                p = pv0[i, j]
                if np.isfinite(v):
                    txt = f"{v:+.2f}{_stars(p)}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=6,
                            color="white" if abs(v) > 0.55 * vmax else "black")
        ax.set_xticks(range(len(ms)), [f"±{m}" for m in ms], fontsize=7)
        ax.set_yticks(range(len(tags)),
                      [t.replace("landmark", "") for t in tags], fontsize=6)
        ax.set_title(f"{pair[0]}-{pair[1]}", fontsize=10)
        ax.set_xlabel("|lag| window (ms)")
    for ax in axes[len(pairs_in_data):]:
        ax.set_visible(False)
    fig.suptitle("Sweep × IFI window (expert epoch) — median IFI across sig dims; "
                 "* p<0.05, ** p<0.01, *** p<0.001 (Wilcoxon vs 0)",
                 fontsize=11)
    fig.colorbar(im, ax=axes.tolist(), shrink=0.7, label="median IFI")
    out = fig_dir / "sweep_ifi_window_per_pair_expert.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def _plot_ifi_window_lr_change(per_cpe, tags, fig_dir):
    """Same 2x4 layout; cell = Δ median IFI (expert - naive),
    * = p_n_vs_e (Mann-Whitney on pooled sig-dim values)."""
    pairs_in_data = [p for p in config.PAIRS
                     if any((t, p, "naive") in per_cpe for t in tags)]
    fig, axes = plt.subplots(2, 4, figsize=(22, 8))
    axes = axes.ravel()
    vmax = 0.0
    grids = []
    for pair in pairs_in_data:
        ms_n, med_n, _ = _ifi_window_grid(per_cpe, tags, pair, "naive")
        ms_e, med_e, _ = _ifi_window_grid(per_cpe, tags, pair, "expert")
        if ms_n.size == 0 or ms_e.size == 0 or ms_n.tolist() != ms_e.tolist():
            grids.append((pair, np.array([]), np.array([]), np.array([])))
            continue
        delta = med_e - med_n
        # p value at each (config, window) via Mann-Whitney on the pooled vals.
        pmat = np.full_like(delta, np.nan)
        for i, tag in enumerate(tags):
            bin_ms = _bin_ms_for_tag(tag)
            for j, win_ms in enumerate(ms_n):
                w = int(win_ms // bin_ms)
                cpe_n = per_cpe.get((tag, pair, "naive"))
                cpe_e = per_cpe.get((tag, pair, "expert"))
                if cpe_n is None or cpe_e is None:
                    continue
                a = cpe_n["ifi_win"].get(w, [])
                b = cpe_e["ifi_win"].get(w, [])
                pmat[i, j] = _mwu(a, b)
        grids.append((pair, ms_n, delta, pmat))
        if delta.size:
            vmax = max(vmax, float(np.nanmax(np.abs(delta))))
    vmax = max(vmax, 0.1)
    for ax, (pair, ms, delta, pmat) in zip(axes, grids):
        if delta.size == 0:
            ax.set_visible(False)
            continue
        im = ax.imshow(delta, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        for i in range(delta.shape[0]):
            for j in range(delta.shape[1]):
                v = delta[i, j]
                p = pmat[i, j]
                if np.isfinite(v):
                    txt = f"{v:+.2f}{_stars(p)}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=6,
                            color="white" if abs(v) > 0.55 * vmax else "black")
        ax.set_xticks(range(len(ms)), [f"±{m}" for m in ms], fontsize=7)
        ax.set_yticks(range(len(tags)),
                      [t.replace("landmark", "") for t in tags], fontsize=6)
        ax.set_title(f"{pair[0]}-{pair[1]}", fontsize=10)
        ax.set_xlabel("|lag| window (ms)")
    for ax in axes[len(pairs_in_data):]:
        ax.set_visible(False)
    fig.suptitle("Sweep × IFI window — Δ median IFI (expert − naive); "
                 "* p<0.05, ** p<0.01, *** p<0.001 (Mann-Whitney on pooled sig dims)",
                 fontsize=11)
    fig.colorbar(im, ax=axes.tolist(), shrink=0.7, label="Δ median IFI")
    out = fig_dir / "sweep_ifi_window_lr_change.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def _plot_pvalues_summary(per_cpe, tags, fig_dir):
    """For each pair, a heatmap of p-values across the headline tests.

    Columns (6): CC: naive≠0, expert≠0, n vs e | IFI(maxwin): naive≠0, expert≠0, n vs e
    Rows: configs.
    Cell color = -log10(p), text = stars.
    """
    pairs_in_data = [p for p in config.PAIRS
                     if any((t, p, "naive") in per_cpe for t in tags)]
    col_labels = ["CC: n≠0", "CC: e≠0", "CC: n vs e",
                  "IFI: n≠0", "IFI: e≠0", "IFI: n vs e"]
    fig, axes = plt.subplots(2, 4, figsize=(22, 9))
    axes = axes.ravel()
    for ax, pair in zip(axes, pairs_in_data):
        pmat = np.full((len(tags), 6), np.nan)
        for i, tag in enumerate(tags):
            cpe_n = per_cpe.get((tag, pair, "naive"))
            cpe_e = per_cpe.get((tag, pair, "expert"))
            if cpe_n is None or cpe_e is None:
                continue
            pmat[i, 0] = _wilcoxon_vs0(cpe_n["cc"])
            pmat[i, 1] = _wilcoxon_vs0(cpe_e["cc"])
            pmat[i, 2] = _mwu(cpe_n["cc"], cpe_e["cc"])
            # IFI at the maximum available window for this config.
            max_w_n = max(cpe_n["ifi_win"]) if cpe_n["ifi_win"] else None
            max_w_e = max(cpe_e["ifi_win"]) if cpe_e["ifi_win"] else None
            if max_w_n:
                pmat[i, 3] = _wilcoxon_vs0(cpe_n["ifi_win"][max_w_n])
            if max_w_e:
                pmat[i, 4] = _wilcoxon_vs0(cpe_e["ifi_win"][max_w_e])
            if max_w_n and max_w_e and max_w_n == max_w_e:
                pmat[i, 5] = _mwu(cpe_n["ifi_win"][max_w_n],
                                  cpe_e["ifi_win"][max_w_e])
        # Color: -log10(p), capped at 4.
        with np.errstate(divide="ignore", invalid="ignore"):
            colour = -np.log10(np.where(np.isfinite(pmat) & (pmat > 0),
                                        pmat, 1.0))
        colour = np.clip(colour, 0, 4)
        im = ax.imshow(colour, cmap="OrRd", vmin=0, vmax=4, aspect="auto")
        for i in range(pmat.shape[0]):
            for j in range(pmat.shape[1]):
                p = pmat[i, j]
                if np.isfinite(p):
                    txt = f"{p:.3f}{_stars(p)}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=6,
                            color="white" if colour[i, j] > 2.5 else "black")
        ax.set_xticks(range(6), col_labels, fontsize=7, rotation=20, ha="right")
        ax.set_yticks(range(len(tags)),
                      [t.replace("landmark", "") for t in tags], fontsize=6)
        ax.set_title(f"{pair[0]}-{pair[1]}", fontsize=10)
    for ax in axes[len(pairs_in_data):]:
        ax.set_visible(False)
    fig.suptitle("Sweep × p-values  "
                 "(CC: Wilcoxon vs 0 / Mann-Whitney; IFI at max |lag| window same; "
                 "stars: * p<.05, ** p<.01, *** p<.001)", fontsize=11)
    fig.colorbar(im, ax=axes.tolist(), shrink=0.7, label="-log10(p)")
    out = fig_dir / "sweep_pvalues_summary.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
