"""Paired per-animal hypothesis-test helpers for learning-change analyses.

Two-sided Wilcoxon signed-rank on per-animal deltas, plus a Benjamini-Hochberg
FDR mask applied within a test family. Shared by the landmark- and spatial-arm
learning-change drivers so the statistics are defined in exactly one place.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def wilcoxon_signed(deltas) -> tuple[int, float, float, float]:
    """Two-sided Wilcoxon signed-rank on per-animal deltas.

    Returns ``(n, median, statistic, p)``. ``p`` is NaN when there are fewer
    than 3 finite, non-zero deltas (the signed-rank test is not meaningful) or
    when every delta is zero.
    """
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
        return (n, float(np.median(arr)), float(res.statistic),
                float(res.pvalue))
    except ValueError:                                  # degenerate / all-zero
        return (n, float(np.median(arr)), float("nan"), float("nan"))


def fdr_bh(pvals, q: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg FDR mask. Operates on finite p-values only.

    Non-finite p-values are excluded from the family (they never pass and do not
    count toward ``n``); the correction is over the finite tests only.
    """
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
