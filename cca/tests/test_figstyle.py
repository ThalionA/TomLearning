"""figstyle.check_text_overflow — the render-time fence for clipped titles.

Promoted from ~/.claude/MISTAKES.md (`no-verification` x3, 2026-08-18): figures were
sent twice with panel titles cut off at the figure edge because only crops were
inspected. `figstyle.save` now measures every text artist against the figure box
and warns; these tests pin that it fires on an overflowing title and stays quiet on
a figure that fits."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import figstyle  # noqa: E402


def test_overflowing_suptitle_is_reported():
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.suptitle("x" * 400, fontsize=12)               # far wider than 4 inches
    issues = figstyle.check_text_overflow(fig)
    plt.close(fig)
    assert issues, "a 400-character suptitle on a 4-inch figure must be flagged"
    assert any("suptitle" in i or "xxx" in i for i in issues)


def test_overflowing_axes_title_is_reported_with_the_side():
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    axes[1].set_title("right-panel title that is much too long to fit " * 4, fontsize=9)
    issues = figstyle.check_text_overflow(fig)
    plt.close(fig)
    assert issues
    assert any("right" in i for i in issues)


def test_fitting_figure_is_quiet():
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot([0, 1], [0, 1])
    ax.set_title("short title"); ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.suptitle("fits fine")
    assert figstyle.check_text_overflow(fig) == []
    plt.close(fig)


def test_overlapping_panel_titles_are_reported():
    fig, axes = plt.subplots(1, 2, figsize=(5, 3))
    fig.subplots_adjust(wspace=0.02)
    for ax in axes:
        ax.set_title("a title long enough to run into its neighbour", fontsize=10)
    issues = figstyle.check_text_overflow(fig)
    plt.close(fig)
    assert any("overlap" in i for i in issues)
