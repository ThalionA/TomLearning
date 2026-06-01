"""Striatum-style Stage-2 figures for the Tom landmark arm.

Adapted from ``StriatumACC/Striatum project/cca/scripts/plot_stage2.py`` --
same convention throughout: **the statistical unit is the significant
canonical subspace dimension**, not the animal and not the (epoch, landmark)
cell. Each fit cell contributes its ``n_sig`` significant dims; the pooled
distribution per (pair, epoch) is the data.

Five figures, all written to ``figures/<tag>/``:

  coverage_<tag>.png            -- diagnostic (Tom-specific): per-(epoch,
                                   landmark) mean n_crossings per pair.
  stage2_comm_strength_<tag>.png  -- held-out peak CC across significant dims,
                                     box+jitter per epoch, per pair.
                                     * = vs-0 Wilcoxon p<alpha; ANOVA + Tukey
                                     across epochs in the panel title.
  stage2_subspace_dim_<tag>.png   -- # significant subspace dims per cell,
                                     box+jitter per epoch, per pair.
                                     Paired and unpaired naive vs expert.
  stage2_lag_curves_<tag>.png     -- mean held-out CC vs lag (ms), one line
                                     per epoch, per pair, pooled across all
                                     significant dims.
  stage2_ifi_winNN_<tag>.png      -- IFI distribution across significant dims
                                     per epoch, per pair. * = vs-0 Wilcoxon;
                                     ANOVA + Tukey. One figure per lag-window
                                     1..max.

Plus one Tom-specific extra:

  stage2_per_landmark_<tag>.png   -- per-pair x per-landmark heatmap of
                                     median held-out CC across sig dims per
                                     epoch, so landmark-specific structure is
                                     visible without picking from 6 panels.

Operates on the existing pkl -- no re-run needed.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import matplotlib  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, subspace_stats  # noqa: E402


EPOCHS = config.EPOCH_NAMES
EPOCH_LABEL = ["naive", "inter", "expert"]
EPOCH_COLOUR = config.EPOCH_COLOURS
ALPHA = 0.05
LANDMARK_IDS = list(range(1, 7))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tag", default="landmark50_res_samp15")
    return p.parse_args()


def _load(tag: str):
    pkl = config.RESULTS_DIR / f"{tag}.pkl"
    if not pkl.is_file():
        raise SystemExit(f"missing {pkl}; run run_temporal.py first")
    with open(pkl, "rb") as f:
        raw = pickle.load(f)
    fits = [r for r in raw["results"] if hasattr(r, "cells")]
    return fits, raw["cfg"]


def _ordered_pairs(fits) -> list[tuple[str, str]]:
    have = {(r.area_x, r.area_y) for r in fits}
    return [p for p in config.PAIRS if p in have]


# --- data access (per-pair, pooled over animals AND landmarks) --------------
def pair_cells(fits, area_x, area_y, epoch):
    """All (animal, landmark) cells for this (pair, epoch) that fit."""
    out = []
    for r in fits:
        if (r.area_x, r.area_y) != (area_x, area_y):
            continue
        for (ep, lm), cell in r.cells.items():
            if ep == epoch and cell.skip_reason is None:
                out.append((r.animal_id, int(lm), cell))
    return out


def sig_dim_values(cells, kind: str, window: int | None = None) -> np.ndarray:
    """Pool ``kind`` over significant dims across every (animal, landmark) cell.

    ``kind`` is "cc" (held-out peak CC), "ifi_win" (IFI at the given window), or
    "n_sig" (per-cell count, returned as one row per cell, not per dim).
    """
    rows = []
    for _, _, cell in cells:
        s = subspace_stats.cell_subspace_stats(cell, alpha=ALPHA)
        if s is None:
            continue
        if kind == "n_sig":
            rows.append(float(s.n_sig))
            continue
        if s.n_sig == 0:
            continue
        if kind == "cc":
            rows.extend(s.peak_cc_per_dim[s.sig_mask].tolist())
        elif kind == "ifi":
            ifi = cell.ifi_per_dim[s.sig_mask]
            rows.extend([float(x) for x in ifi if np.isfinite(x)])
        elif kind == "ifi_win":
            assert window is not None
            arr = cell.ifi_windows[s.sig_mask, window - 1]
            rows.extend([float(x) for x in arr if np.isfinite(x)])
        else:
            raise ValueError(kind)
    return np.array(rows, dtype=float)


def per_animal_means(cells, kind: str, window: int | None = None) -> dict[int, float]:
    """Per-animal mean over (landmarks × sig dims) -- one number per animal."""
    by_animal: dict[int, list[float]] = {}
    for animal_id, _, cell in cells:
        s = subspace_stats.cell_subspace_stats(cell, alpha=ALPHA)
        if s is None or (kind != "n_sig" and s.n_sig == 0):
            continue
        if kind == "n_sig":
            vals = [float(s.n_sig)]
        elif kind == "cc":
            vals = s.peak_cc_per_dim[s.sig_mask].tolist()
        elif kind == "ifi":
            vals = [float(x) for x in cell.ifi_per_dim[s.sig_mask]
                    if np.isfinite(x)]
        elif kind == "ifi_win":
            vals = [float(x) for x in cell.ifi_windows[s.sig_mask, window - 1]
                    if np.isfinite(x)]
        by_animal.setdefault(animal_id, []).extend(vals)
    return {a: float(np.nanmean(v)) for a, v in by_animal.items() if v}


# --- stats helpers (mirror plot_stage2 in striatum) -------------------------
def _wilcoxon_vs0(v: np.ndarray) -> float:
    v = v[np.isfinite(v)]
    if v.size < 6 or not np.any(v != 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(v).pvalue)
    except ValueError:
        return float("nan")


def _paired(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if a.size < 6 or not np.any(a != b):
        return float("nan")
    try:
        return float(stats.wilcoxon(a, b).pvalue)
    except ValueError:
        return float("nan")


def _unpaired(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return float("nan")
    try:
        return float(stats.mannwhitneyu(a, b).pvalue)
    except ValueError:
        return float("nan")


def _anova(per_epoch: list[np.ndarray]):
    groups = [v[np.isfinite(v)] for v in per_epoch]
    if min(g.size for g in groups) < 2:
        return float("nan"), (float("nan"),) * 3
    omnibus = float(stats.f_oneway(*groups).pvalue)
    th = stats.tukey_hsd(*groups).pvalue
    return omnibus, (float(th[0, 1]), float(th[1, 2]), float(th[0, 2]))


def _fmt_p(p):
    return f"{p:.3f}" if np.isfinite(p) else "n/a"


# --- plotting primitives ----------------------------------------------------
def _grid(n_pairs: int):
    rows = 2
    cols = (n_pairs + 1) // 2
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 6.8))
    return fig, np.atleast_1d(axes).ravel()


def _save(fig, name, tag, fig_dir):
    fig_dir.mkdir(parents=True, exist_ok=True)
    path = fig_dir / f"{name}_{tag}.png"
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _box_by_epoch(ax, per_epoch):
    for i, (epoch, v) in enumerate(zip(EPOCHS, per_epoch)):
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        colour = EPOCH_COLOUR[epoch]
        bp = ax.boxplot([v], positions=[i], widths=0.55, showfliers=False,
                        patch_artist=True, medianprops=dict(color="k"))
        bp["boxes"][0].set(facecolor=colour, alpha=0.35)
        jit = (np.random.default_rng(i).random(v.size) - 0.5) * 0.3
        ax.scatter(i + jit, v, s=10, color=colour, alpha=0.45, zorder=3)
    ax.set_xticks(range(len(EPOCHS)))
    ax.set_xticklabels(EPOCH_LABEL)
    ax.set_xlim(-0.6, len(EPOCHS) - 0.4)


def _star_vs0(ax, per_epoch):
    for i, v in enumerate(per_epoch):
        p = _wilcoxon_vs0(v)
        if np.isfinite(p) and p < ALPHA:
            ax.text(i, 0.93, "*", transform=ax.get_xaxis_transform(),
                    fontsize=14, ha="center")


def _stat_title(pair, n_animals, n_dims, anova_p, tukey):
    tk = "  ".join(f"{lab} {_fmt_p(p)}"
                   for lab, p in zip(("n-i", "i-e", "n-e"), tukey))
    return (f"{pair}  n={n_animals}  {n_dims}d\n"
            f"ANOVA p={_fmt_p(anova_p)}   Tukey {tk}")


# --- main figures -----------------------------------------------------------
def plot_comm_strength(fits, pairs, fig_dir, tag):
    fig, axes = _grid(len(pairs))
    for ax, (ax_x, ax_y) in zip(axes, pairs):
        per_epoch = [
            sig_dim_values(pair_cells(fits, ax_x, ax_y, e), "cc")
            for e in EPOCHS
        ]
        _box_by_epoch(ax, per_epoch)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_ylim(0, 1)
        _star_vs0(ax, per_epoch)
        anova_p, tukey = _anova(per_epoch)
        cells_naive = pair_cells(fits, ax_x, ax_y, "naive")
        n_animals = len({a for a, _, _ in cells_naive})
        n_dims = sum(v.size for v in per_epoch)
        ax.set_title(_stat_title(f"{ax_x}-{ax_y}", n_animals, n_dims,
                                 anova_p, tukey), fontsize=7.5)
    for ax in axes[: len(pairs):len(pairs) // 2 or 1]:
        ax.set_ylabel("held-out peak CC  (significant dims)")
    for ax in axes[len(pairs):]:
        ax.set_visible(False)
    fig.suptitle(f"Stage 2 -- communication strength (landmark arm, {tag}; "
                 f"learners)  '*' CC≠0 (Wilcoxon); ANOVA+Tukey over sig dims",
                 fontsize=10)
    return _save(fig, "stage2_comm_strength", tag, fig_dir)


def plot_subspace_dim(fits, pairs, fig_dir, tag):
    fig, axes = _grid(len(pairs))
    for ax, (ax_x, ax_y) in zip(axes, pairs):
        per_epoch = [
            sig_dim_values(pair_cells(fits, ax_x, ax_y, e), "n_sig")
            for e in EPOCHS
        ]
        _box_by_epoch(ax, per_epoch)
        # naive vs expert: paired across animals (per-animal mean), unpaired
        # over the cell-level n_sig values.
        naive_anim = per_animal_means(pair_cells(fits, ax_x, ax_y, "naive"), "n_sig")
        expert_anim = per_animal_means(pair_cells(fits, ax_x, ax_y, "expert"), "n_sig")
        common = sorted(set(naive_anim) & set(expert_anim))
        a_n = np.array([naive_anim[a] for a in common])
        a_e = np.array([expert_anim[a] for a in common])
        pp = _paired(a_n, a_e)
        up = _unpaired(per_epoch[0], per_epoch[2])
        ax.set_ylim(bottom=0)
        cells_naive = pair_cells(fits, ax_x, ax_y, "naive")
        n_animals = len({a for a, _, _ in cells_naive})
        ax.set_title(f"{ax_x}-{ax_y}  n={n_animals}\n"
                     f"n→e pair={_fmt_p(pp)} unp={_fmt_p(up)}", fontsize=8)
    for i, ax in enumerate(axes):
        if i < len(pairs) and i % (len(pairs) // 2 or 1) == 0:
            ax.set_ylabel("# significant subspace dims  (per cell)")
    for ax in axes[len(pairs):]:
        ax.set_visible(False)
    fig.suptitle(f"Stage 2 -- communication-subspace dimensionality "
                 f"(landmark arm, {tag}; learners)", fontsize=10)
    return _save(fig, "stage2_subspace_dim", tag, fig_dir)


def plot_lag_curves(fits, pairs, fig_dir, tag, cfg):
    fig, axes = _grid(len(pairs))
    bin_ms = int(cfg.temporal_bin_ms)
    for ax, (ax_x, ax_y) in zip(axes, pairs):
        n_dims_total = 0
        for epoch in EPOCHS:
            curves = []
            ref_lags = None
            for _, _, cell in pair_cells(fits, ax_x, ax_y, epoch):
                s = subspace_stats.cell_subspace_stats(cell, alpha=ALPHA)
                if s is None or s.n_sig == 0:
                    continue
                if ref_lags is None:
                    ref_lags = cell.lags
                if cell.lags.size != ref_lags.size:
                    continue
                for j in np.where(s.sig_mask)[0]:
                    curves.append(cell.lag_cc_per_dim[:, j])
            if not curves:
                continue
            arr = np.stack(curves, axis=0)
            n_dims_total += arr.shape[0]
            mean = np.nanmean(arr, axis=0)
            sem = np.nanstd(arr, axis=0) / np.sqrt(arr.shape[0])
            lag_ms = ref_lags * bin_ms
            ax.plot(lag_ms, mean, "-", color=EPOCH_COLOUR[epoch], lw=1.6,
                    label=epoch)
            ax.fill_between(lag_ms, mean - sem, mean + sem,
                            color=EPOCH_COLOUR[epoch], alpha=0.2)
        ax.axvline(0, color="k", lw=0.5, ls=":")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(f"{ax_x}-{ax_y}  ({n_dims_total} dims)", fontsize=10)
    for i, ax in enumerate(axes):
        if i < len(pairs) and i % (len(pairs) // 2 or 1) == 0:
            ax.set_ylabel("held-out CC  (significant dims)")
        if i >= len(pairs) // 2:
            ax.set_xlabel("lag (ms)  +ve: X leads Y")
    axes[0].legend(frameon=False, fontsize=8)
    for ax in axes[len(pairs):]:
        ax.set_visible(False)
    fig.suptitle(f"Stage 2 -- lagged-refit CC over significant subspace dims "
                 f"(landmark arm, {tag}; learner mean ± SEM)", fontsize=10)
    return _save(fig, "stage2_lag_curves", tag, fig_dir)


def plot_ifi_window(fits, pairs, fig_dir, tag, window: int, cfg):
    bin_ms = int(cfg.temporal_bin_ms)
    fig, axes = _grid(len(pairs))
    for ax, (ax_x, ax_y) in zip(axes, pairs):
        per_epoch = [
            sig_dim_values(pair_cells(fits, ax_x, ax_y, e), "ifi_win",
                           window=window)
            for e in EPOCHS
        ]
        _box_by_epoch(ax, per_epoch)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_ylim(-1.05, 1.05)
        _star_vs0(ax, per_epoch)
        anova_p, tukey = _anova(per_epoch)
        cells_naive = pair_cells(fits, ax_x, ax_y, "naive")
        n_animals = len({a for a, _, _ in cells_naive})
        n_dims = sum(v.size for v in per_epoch)
        ax.set_title(_stat_title(f"{ax_x}-{ax_y}", n_animals, n_dims,
                                 anova_p, tukey), fontsize=7.5)
    for i, ax in enumerate(axes):
        if i < len(pairs) and i % (len(pairs) // 2 or 1) == 0:
            ax.set_ylabel("IFI  (+ve: X leads Y)")
    for ax in axes[len(pairs):]:
        ax.set_visible(False)
    fig.suptitle(f"Stage 2 -- IFI, |lag| ≤ {window} bins (= {window*bin_ms} ms) "
                 f"(landmark arm, {tag}; learners)  '*' IFI≠0; ANOVA+Tukey",
                 fontsize=10)
    return _save(fig, f"stage2_ifi_win{window:02d}", tag, fig_dir)


# --- Tom-specific diagnostics ----------------------------------------------
def plot_coverage(fits, pairs, fig_dir, tag):
    fig, axes = _grid(len(pairs))
    vmax = 1.0
    grids = []
    for pair in pairs:
        grid = np.full((len(EPOCHS), len(LANDMARK_IDS)), np.nan)
        for ei, e in enumerate(EPOCHS):
            for li, lm in enumerate(LANDMARK_IDS):
                counts = [r.coverage.get((e, lm))
                          for r in fits if (r.area_x, r.area_y) == pair
                          and r.coverage.get((e, lm)) is not None]
                if counts:
                    grid[ei, li] = float(np.mean(counts))
        vmax = max(vmax, np.nanmax(grid) if np.isfinite(grid).any() else 1)
        grids.append(grid)
    for ax, pair, grid in zip(axes, pairs, grids):
        im = ax.imshow(grid, cmap="viridis", vmin=0, vmax=vmax, aspect="auto")
        for ei in range(len(EPOCHS)):
            for li in range(len(LANDMARK_IDS)):
                v = grid[ei, li]
                if np.isfinite(v):
                    ax.text(li, ei, f"{v:.0f}", ha="center", va="center",
                            color="white" if v < vmax * 0.5 else "black",
                            fontsize=8)
        ax.set_title(f"{pair[0]}-{pair[1]}", fontsize=9)
        ax.set_xticks(range(len(LANDMARK_IDS)),
                      [f"L{l}" for l in LANDMARK_IDS])
        ax.set_yticks(range(len(EPOCHS)), EPOCH_LABEL)
    for ax in axes[len(pairs):]:
        ax.set_visible(False)
    fig.suptitle(f"Coverage: mean n_crossings per (epoch, landmark)  "
                 f"({tag})", fontsize=10)
    fig.colorbar(im, ax=axes.tolist(), shrink=0.6, label="n_crossings")
    return _save(fig, "coverage", tag, fig_dir)


def plot_per_landmark_detail(fits, pairs, fig_dir, tag):
    """Median held-out CC per (pair x landmark x epoch) -- shows where the
    landmark axis matters, complementing the pooled main figures."""
    rows = []
    for pair in pairs:
        for li, lm in enumerate(LANDMARK_IDS):
            for ei, e in enumerate(EPOCHS):
                cells = [(a, l, c) for r in fits if (r.area_x, r.area_y) == pair
                         for (ep, l), c in r.cells.items()
                         if ep == e and l == lm and c.skip_reason is None
                         for a in [r.animal_id]]
                ccs = sig_dim_values(cells, "cc")
                med = float(np.median(ccs)) if ccs.size else np.nan
                rows.append((pair, lm, e, med, int(ccs.size)))
    fig, axes = plt.subplots(len(EPOCHS), 1, figsize=(11, 6), sharex=True)
    pair_labels = [f"{p[0]}-{p[1]}" for p in pairs]
    vmax = 0.7
    for ei, e in enumerate(EPOCHS):
        grid = np.full((len(pairs), len(LANDMARK_IDS)), np.nan)
        ns = np.zeros_like(grid, dtype=int)
        for (pair, lm, ep, med, n) in rows:
            if ep != e:
                continue
            i = pairs.index(pair)
            j = LANDMARK_IDS.index(lm)
            grid[i, j] = med
            ns[i, j] = n
        ax = axes[ei]
        im = ax.imshow(grid, cmap="viridis", vmin=0, vmax=vmax, aspect="auto")
        for i in range(len(pairs)):
            for j in range(len(LANDMARK_IDS)):
                v = grid[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}\n({ns[i, j]})", ha="center",
                            va="center", fontsize=7,
                            color="white" if v > 0.4 else "black")
        ax.set_yticks(range(len(pairs)), pair_labels, fontsize=8)
        ax.set_title(f"{e}", fontsize=10)
        if ei == len(EPOCHS) - 1:
            ax.set_xticks(range(len(LANDMARK_IDS)),
                          [f"L{l}" for l in LANDMARK_IDS])
        else:
            ax.set_xticks([])
    fig.suptitle(f"Per-landmark detail: median held-out peak CC across sig "
                 f"dims  (cell text: median (n_sig_dims))  ({tag})",
                 fontsize=10)
    fig.colorbar(im, ax=axes.tolist(), shrink=0.7, label="median CC over sig dims")
    return _save(fig, "stage2_per_landmark", tag, fig_dir)


def main():
    args = parse_args()
    fits, cfg = _load(args.tag)
    pairs = _ordered_pairs(fits)
    fig_dir = config.FIGURES_DIR / args.tag
    print(f"[plot_landmark] tag={args.tag}: {len(fits)} fits, "
          f"{len(pairs)} pairs -> {fig_dir}")

    outputs = [
        plot_coverage(fits, pairs, fig_dir, args.tag),
        plot_comm_strength(fits, pairs, fig_dir, args.tag),
        plot_subspace_dim(fits, pairs, fig_dir, args.tag),
        plot_lag_curves(fits, pairs, fig_dir, args.tag, cfg),
        plot_per_landmark_detail(fits, pairs, fig_dir, args.tag),
    ]
    # IFI windows: one figure per integration window.
    sample = next(iter(fits[0].cells.values()))
    while sample.skip_reason is not None:
        sample = next(iter(fits[0].cells.values()))
    max_window = sample.ifi_windows.shape[1]
    for window in range(1, max_window + 1):
        outputs.append(
            plot_ifi_window(fits, pairs, fig_dir, args.tag, window, cfg))
    for p in outputs:
        print(f"  wrote {p}")


if __name__ == "__main__":
    main()
