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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAX_PX = 1600          # repo rule: PNG longest side <= 1600 px
STAR_COLOR = "#c0392b"


def apply():
    """Install the shared rcParams (call once before building figures)."""
    plt.rcParams.update({
        "figure.constrained_layout.use": True,      # auto-spacing, no label collisions
        "figure.constrained_layout.h_pad": 0.06,
        "figure.constrained_layout.w_pad": 0.06,
        "figure.constrained_layout.hspace": 0.05,
        "figure.constrained_layout.wspace": 0.05,
        "axes.spines.top": False,                   # de-clutter
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "axes.titlepad": 4,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "font.size": 10,
        "savefig.facecolor": "white",
    })


def save(fig, stem, max_px: int = MAX_PX):
    """Write ``<stem>.png`` (longest side <= ``max_px``) AND ``<stem>.svg``, then close.

    ``stem`` may include or omit a ``.png`` suffix; it is stripped. The PNG dpi is
    chosen so the longest side is capped at ``max_px`` (never upscaled past 150)."""
    stem = str(stem)
    if stem.endswith(".png") or stem.endswith(".svg"):
        stem = stem[:-4]
    w_in, h_in = fig.get_size_inches()
    dpi = min(150.0, max_px / max(w_in, h_in))
    fig.savefig(f"{stem}.png", dpi=dpi)
    fig.savefig(f"{stem}.svg")
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
