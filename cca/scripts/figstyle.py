"""Shared figure style + saver for the CCA report figures.

Enforces the repo standard (a `.svg` + `.png` pair, PNG longest side <= 1600 px) and a
clean, readable house style: constrained layout so axis labels/titles never collide,
de-spined axes and no gridlines to cut clutter, restrained fonts, and a prominent
significance-star helper that always draws above the data. Usage:

    import figstyle
    figstyle.apply()                     # once, before creating figures
    ...
    figstyle.save(fig, ATT / "HCV1_X")   # writes HCV1_X.png (<=1600 px) + HCV1_X.svg
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAX_PX = 1600          # repo rule: PNG longest side <= 1600 px
STAR_COLOR = "#c0392b"
# Every figure is ALSO written here, so the project folder holds its own figures
# rather than them existing only inside the Obsidian vault. The temporal-arm scripts
# each hardcoded the vault attachments as their only destination, which meant a
# checkout of this repo had no temporal figures at all. Mirroring at the single
# shared saver fixes it for every script at once. (Gitignored, like the landmark and
# spatial figures — they regenerate from the scripts.)
REPO_FIGURES = Path(__file__).resolve().parents[1] / "figures"


def apply():
    """Install the shared rcParams (call once before building figures). Generous
    spacing + fonts so panels breathe even when Obsidian scales them to note width."""
    plt.rcParams.update({
        "figure.constrained_layout.use": True,      # auto-spacing, no label collisions
        "figure.constrained_layout.h_pad": 0.14,    # real whitespace between panels
        "figure.constrained_layout.w_pad": 0.14,
        "figure.constrained_layout.hspace": 0.12,
        "figure.constrained_layout.wspace": 0.10,
        "axes.spines.top": False,                   # de-clutter
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "axes.titlepad": 8,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "legend.frameon": False,
        "font.size": 12,
        "lines.linewidth": 1.8,
        "axes.linewidth": 0.9,
        "savefig.facecolor": "white",
    })


def grid(n, ncols=2, panel=(5.4, 3.5)):
    """Spacious per-pair panel grid: few columns + generous panel size so each
    panel has real room (default ~5.4×3.5 in). Returns ``(fig, axes_flat)``. For
    8 pairs this is a tall 4×2 — breathing room over compactness, by design."""
    import numpy as np
    nrows = int(np.ceil(n / max(1, ncols)))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(panel[0] * ncols, panel[1] * nrows),
                             squeeze=False)
    return fig, axes.ravel()


def save(fig, stem, max_px: int = MAX_PX, mirror: bool = True):
    """Write ``<stem>.png`` (longest side <= ``max_px``) AND ``<stem>.svg``, then close.

    ``stem`` may include or omit a ``.png`` suffix; it is stripped. The PNG dpi is
    chosen so the longest side is capped at ``max_px`` (never upscaled past 150).

    The pair is ALSO mirrored into :data:`REPO_FIGURES` (unless ``mirror=False``) so the
    project folder holds every figure, not just the vault. Pass ``mirror=False`` for a
    figure that is genuinely vault-only.
    """
    stem = str(stem)
    if stem.endswith(".png") or stem.endswith(".svg"):
        stem = stem[:-4]
    w_in, h_in = fig.get_size_inches()
    dpi = min(150.0, max_px / max(w_in, h_in))
    fig.savefig(f"{stem}.png", dpi=dpi)
    fig.savefig(f"{stem}.svg")
    if mirror:
        try:
            REPO_FIGURES.mkdir(parents=True, exist_ok=True)
            name = Path(stem).name
            if Path(stem).resolve().parent != REPO_FIGURES.resolve():
                # COPY the rendered files rather than calling savefig again.
                # constrained_layout re-solves on every draw, so a second render can
                # settle sub-pixel differently and the two copies would not be
                # byte-identical (4 of 24 were not). Copying also halves the cost.
                for ext in ("png", "svg"):
                    shutil.copyfile(f"{stem}.{ext}", REPO_FIGURES / f"{name}.{ext}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  [figstyle] repo mirror failed for {stem}: {exc}")
    plt.close(fig)


def star(ax, x, y, p=None, *, always: bool = False, size: int = 15):
    """Draw a bold, offset significance star ABOVE (x, y) when ``p < 0.05``.

    Offset upward in display points + high zorder + ``clip_on=False`` so it is never
    hidden behind bars or data points. Pass ``always=True`` to draw regardless of p."""
    import numpy as np
    ok = always or (p is not None and np.isfinite(p) and p < 0.05)
    if not ok or not np.isfinite(y):
        return
    ax.annotate("*", (x, y), xytext=(0, 5), textcoords="offset points",
                ha="center", va="bottom", fontsize=size, fontweight="bold",
                color=STAR_COLOR, clip_on=False, zorder=20)
