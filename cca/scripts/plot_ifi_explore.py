"""IFI deep-dive -- directionality across learning, and lag-window structure.

The sweep summary CSV pools naive + expert and reports IFI only at windows
1..10. This script goes back to the per-epoch data in the Stage 2 pkls and
focuses on the Information Flow Index (IFI):

  1. directionality across learning -- IFI at naive / intermediate / expert,
     per area-pair;
  2. lag-window dependence and the peak lead/lag, in bins and cm.

IFI is read from the *dominant canonical dimension* (dim 0) -- the pipeline's
own headline ``EpochAnalysis.ifi`` -- so a value is defined for every
(config, pair, epoch). Sign convention: +ve = X leads Y (X = first-listed
area of the pair). Uses the lag20 configs (22 of them: 2 CCA types x 11
PC-rules), which carry the full +/-20-bin window; for windows 1..10 the lag10
configs are identical. Learner cohort only.

Writes into figures/ifi_explore/:
  * ifi_across_learning.png   IFI vs epoch, per pair, one line per config
  * ifi_change_naive_expert.png  per-pair IFI(expert) - IFI(naive)
  * ifi_vs_window.png         IFI vs lag-integration window, naive vs expert
  * ifi_peak_lag.png          peak lead/lag (cm) per pair, naive vs expert

Run:  python scripts/plot_ifi_explore.py
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tom_cca import config  # noqa: E402

OUTDIR = config.FIGURES_DIR / "ifi_explore"
EPOCHS = list(config.EPOCH_NAMES)                 # naive, intermediate, expert
PAIRS = [f"{x}-{y}" for x, y in config.PAIRS]
BIN_CM = config.CORRIDOR_CM / config.N_BINS       # 2.5 cm per bin
STD_W = 10                                        # headline window: |lag|<=10
RES_COLOUR = "#3c6fb0"
SIG_COLOUR = "#d4753a"
EPOCH_COLOUR = {"naive": "#4c72b0", "intermediate": "#dd8452",
                "expert": "#55a868"}


# --------------------------------------------------------------------------
# load per-epoch IFI from the lag20 Stage 2 pkls
# --------------------------------------------------------------------------
def load():
    """data[tag] = {'cca','pairs': {pair: {epoch: {'ifi_w','peak'}}}}.

    ``ifi_w``  -- (n_windows,) IFI of dim 0, averaged over learner animals.
    ``peak``   -- list of dim-0 peak lags (cm), one per learner animal.
    """
    paths = sorted(config.RESULTS_DIR.glob("stage2_*_lag20.pkl"))
    if not paths:
        sys.exit("no stage2_*_lag20.pkl found -- run the sweep first")
    data: dict[str, dict] = {}
    for path in paths:
        tag = path.stem[len("stage2_"):]
        with open(path, "rb") as fh:
            results = pickle.load(fh)["results"]
        cca = "residual" if tag.startswith("res_") else "signal"
        pairs: dict[str, dict] = {}
        for (ax, ay), name in zip(config.PAIRS, PAIRS):
            rs = [r for r in results
                  if (r.area_x, r.area_y) == (ax, ay) and r.role == "learner"]
            ep: dict[str, dict] = {}
            for epoch in EPOCHS:
                ifi_rows, peaks = [], []
                for r in rs:
                    ea = r.epochs.get(epoch)
                    if ea is None:
                        continue
                    ifi_rows.append(np.asarray(ea.ifi_windows[0], float))
                    peaks.append(float(ea.peak_lag_per_dim[0]) * BIN_CM)
                if ifi_rows:
                    ep[epoch] = {
                        "ifi_w": np.nanmean(np.vstack(ifi_rows), axis=0),
                        "peak": peaks,
                    }
            if ep:
                pairs[name] = ep
        data[tag] = {"cca": cca, "pairs": pairs}
        del results
    return data


def _ifi_at(entry, pair, epoch, w=STD_W):
    """Dim-0 IFI of one config at lag-window ``w`` for one pair x epoch."""
    ep = entry["pairs"].get(pair, {}).get(epoch)
    if ep is None or w - 1 >= ep["ifi_w"].size:
        return np.nan
    return float(ep["ifi_w"][w - 1])


def _save(fig, name):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    path = OUTDIR / f"{name}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"saved {path}")


# --------------------------------------------------------------------------
# figure 1 -- IFI across learning, per pair
# --------------------------------------------------------------------------
def fig_across_learning(data):
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.8))
    x = np.arange(len(EPOCHS))
    seen: list[float] = []
    for ax, pair in zip(axes.ravel(), PAIRS):
        means = {"residual": [], "signal": []}
        for tag, entry in data.items():
            y = [_ifi_at(entry, pair, e) for e in EPOCHS]
            if all(np.isnan(v) for v in y):
                continue
            colour = RES_COLOUR if entry["cca"] == "residual" else SIG_COLOUR
            ax.plot(x, y, "-", color=colour, lw=0.7, alpha=0.3, zorder=1)
            means[entry["cca"]].append(y)
            seen += [v for v in y if np.isfinite(v)]
        for cca, colour in (("residual", RES_COLOUR), ("signal", SIG_COLOUR)):
            if means[cca]:
                m = np.nanmean(np.array(means[cca], float), axis=0)
                ax.plot(x, m, "-o", color=colour, lw=2.4, ms=5, zorder=3,
                        label=f"{cca} mean")
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(["naive", "inter", "expert"])
        ax.set_xlim(-0.3, len(EPOCHS) - 0.7)
        ax.set_title(pair, fontsize=10)
    # Shared symmetric y-limit, scaled to the data (IFI itself lives in
    # [-1, 1]; these effects are far smaller, so a fixed +/-1 axis would hide
    # the structure). Floor at +/-0.5 so small values are not exaggerated.
    lim = max(0.5, 1.15 * (max(abs(v) for v in seen) if seen else 0.5))
    for ax in axes.ravel():
        ax.set_ylim(-lim, lim)
    for ax in axes[:, 0]:
        ax.set_ylabel("IFI  (+ve: X leads Y)", fontsize=8)
    axes[0, 0].legend(frameon=False, fontsize=7, loc="upper left")
    fig.suptitle(f"IFI across learning, per area-pair  (|lag| <= {STD_W} bins; "
                 f"dim 0; one faint line per config, bold = CCA-type mean)",
                 fontsize=11)
    _save(fig, "ifi_across_learning")


# --------------------------------------------------------------------------
# figure 2 -- naive -> expert change in IFI, per pair
# --------------------------------------------------------------------------
def fig_change(data):
    fig, ax = plt.subplots(figsize=(11, 5.4))
    x = np.arange(len(PAIRS))
    for i, (cca, colour, off) in enumerate(
            (("residual", RES_COLOUR, -0.19), ("signal", SIG_COLOUR, 0.19))):
        for j, pair in enumerate(PAIRS):
            deltas = []
            for entry in data.values():
                if entry["cca"] != cca:
                    continue
                d = (_ifi_at(entry, pair, "expert")
                     - _ifi_at(entry, pair, "naive"))
                if np.isfinite(d):
                    deltas.append(d)
            if not deltas:
                continue
            bp = ax.boxplot([deltas], positions=[j + off], widths=0.32,
                            patch_artist=True, showfliers=False,
                            medianprops=dict(color="k"))
            bp["boxes"][0].set(facecolor=colour, alpha=0.55)
            jit = (np.random.default_rng(j + i).random(len(deltas)) - 0.5) * 0.2
            ax.scatter(j + off + jit, deltas, s=9, color=colour, alpha=0.6,
                       zorder=3)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(PAIRS)
    ax.set_ylabel("IFI(expert) - IFI(naive)   (+ve: flow shifts toward X)")
    handles = [plt.Line2D([], [], marker="s", ls="", color=RES_COLOUR,
                          label="residual"),
               plt.Line2D([], [], marker="s", ls="", color=SIG_COLOUR,
                          label="signal")]
    ax.legend(handles=handles, frameon=False, fontsize=9)
    ax.set_title("Change in information-flow direction across learning "
                 "(box = spread over the 22 configs; |lag| <= 10 bins, dim 0)",
                 fontsize=10)
    _save(fig, "ifi_change_naive_expert")


# --------------------------------------------------------------------------
# figure 3 -- IFI vs lag-integration window, naive vs expert
# --------------------------------------------------------------------------
def fig_vs_window(data):
    n_win = max(ep["ifi_w"].size
                for e in data.values() for p in e["pairs"].values()
                for ep in p.values())
    w = np.arange(1, n_win + 1)
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.8))
    seen: list[float] = []
    for ax, pair in zip(axes.ravel(), PAIRS):
        for cca, ls in (("residual", "--"), ("signal", "-")):
            for epoch in ("naive", "expert"):
                rows = []
                for entry in data.values():
                    if entry["cca"] != cca:
                        continue
                    ep = entry["pairs"].get(pair, {}).get(epoch)
                    if ep is not None and ep["ifi_w"].size == n_win:
                        rows.append(ep["ifi_w"])
                if not rows:
                    continue
                m = np.nanmean(np.array(rows, float), axis=0)
                ax.plot(w, m, ls, color=EPOCH_COLOUR[epoch], lw=1.8,
                        label=f"{epoch}, {cca}")
                seen += [v for v in m if np.isfinite(v)]
        ax.axhline(0, color="k", lw=0.7)
        ax.axvline(STD_W, color="0.6", ls=":", lw=0.8)
        ax.set_title(pair, fontsize=10)
    # Shared symmetric y-limit scaled to the data (see fig_across_learning).
    lim = max(0.5, 1.15 * (max(abs(v) for v in seen) if seen else 0.5))
    for ax in axes.ravel():
        ax.set_ylim(-lim, lim)
    for ax in axes[:, 0]:
        ax.set_ylabel("IFI  (+ve: X leads Y)", fontsize=8)
    for ax in axes[1, :]:
        ax.set_xlabel("lag window (bins; 1 bin = 2.5 cm)", fontsize=8)
    axes[0, 0].legend(frameon=False, fontsize=6.5, loc="upper right")
    fig.suptitle("IFI vs lag-integration window, naive vs expert  "
                 "(mean over configs; solid = signal, dashed = residual; "
                 "dotted = committed |lag|<=10)", fontsize=11)
    _save(fig, "ifi_vs_window")


# --------------------------------------------------------------------------
# figure 4 -- peak lead/lag, per pair, naive vs expert
# --------------------------------------------------------------------------
def fig_peak_lag(data):
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.8))
    span = STD_W * BIN_CM
    for ax, pair in zip(axes.ravel(), PAIRS):
        cols = []
        for epoch in ("naive", "expert"):
            vals = []
            for entry in data.values():
                ep = entry["pairs"].get(pair, {}).get(epoch)
                if ep is not None:
                    vals.extend(v for v in ep["peak"] if np.isfinite(v))
            cols.append(np.array(vals, float))
        for i, (epoch, vals) in enumerate(zip(("naive", "expert"), cols)):
            if vals.size == 0:
                continue
            bp = ax.boxplot([vals], positions=[i], widths=0.55,
                            patch_artist=True, showfliers=False, vert=True,
                            medianprops=dict(color="k"))
            bp["boxes"][0].set(facecolor=EPOCH_COLOUR[epoch], alpha=0.55)
            jit = (np.random.default_rng(i).random(vals.size) - 0.5) * 0.3
            ax.scatter(i + jit, vals, s=7, color=EPOCH_COLOUR[epoch],
                       alpha=0.4, zorder=3)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["naive", "expert"])
        ax.set_ylim(-span - 4, span + 4)
        ax.set_title(pair, fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel("peak lag (cm; +ve: X leads Y)", fontsize=8)
    fig.suptitle("Peak lead/lag of the dominant canonical dimension, per pair "
                 "(pooled over configs x animals; box = quartiles)",
                 fontsize=11)
    _save(fig, "ifi_peak_lag")


def main():
    data = load()
    n_res = sum(e["cca"] == "residual" for e in data.values())
    print(f"loaded {len(data)} lag20 configs "
          f"({n_res} residual, {len(data) - n_res} signal)")
    fig_across_learning(data)
    fig_change(data)
    fig_vs_window(data)
    fig_peak_lag(data)
    print("IFI deep-dive done.")


if __name__ == "__main__":
    main()
