"""Per-unit spatial reliability from the trial × spatial-bin firing tensor.

Reliability of unit u at trial t = the mean Pearson correlation between the
unit's spatially-binned map on trial t and its maps on the neighbouring trials
within ``±half_window`` (session edges clipped: fewer neighbours, never padded,
never wrapped). Computed from ``Animal.spatial_fr`` directly — NOT from the
export's precomputed ``analysis_spatial/reliability`` fields, whose provenance
(smoothing, window, correlation type) we cannot verify.

Conventions:
- A trial pair correlates over the spatial bins finite in BOTH maps; pairs with
  fewer than ``min_bins`` such bins, or with zero variance on either side, are
  NaN and drop out of the mean (a silent trial does not drag reliability to 0).
- A trial whose every neighbour pair is NaN gets NaN reliability.
- ``epoch_mean_reliability`` takes 1-BASED trial ids (the temporal stream's
  ``trial_idx_50ms`` convention, measured on TF073: ids 1..214 for 214 spatial
  rows) and nan-means the per-trial reliability over them.
"""

from __future__ import annotations

import numpy as np


def _pair_corr(a: np.ndarray, b: np.ndarray, min_bins: int) -> np.ndarray:
    """Pearson r per unit between two (n_bins, n_units) maps, pairwise-complete.

    NaN where fewer than ``min_bins`` bins are finite in both maps or either
    side has zero variance over the shared bins.
    """
    finite = np.isfinite(a) & np.isfinite(b)
    n = finite.sum(axis=0)
    a0 = np.where(finite, a, 0.0)
    b0 = np.where(finite, b, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ma = a0.sum(axis=0) / n
        mb = b0.sum(axis=0) / n
        da = np.where(finite, a0 - ma, 0.0)
        db = np.where(finite, b0 - mb, 0.0)
        cov = (da * db).sum(axis=0)
        va = (da * da).sum(axis=0)
        vb = (db * db).sum(axis=0)
        r = cov / np.sqrt(va * vb)
    bad = (n < min_bins) | (va <= 0) | (vb <= 0)
    return np.where(bad, np.nan, r)


def trial_map_reliability(spatial_fr: np.ndarray, half_window: int = 2,
                          min_bins: int = 5) -> np.ndarray:
    """Sliding ±``half_window``-trial mean map correlation, per (trial, unit).

    ``spatial_fr`` is ``(n_trials, n_bins, n_units)`` (``Animal.spatial_fr``
    order). Returns ``(n_trials, n_units)``.
    """
    fr = np.asarray(spatial_fr, dtype=float)
    if fr.ndim != 3:
        raise ValueError(f"expected (n_trials, n_bins, n_units), got {fr.shape}")
    n_trials, _, n_units = fr.shape
    sums = np.zeros((n_trials, n_units))
    counts = np.zeros((n_trials, n_units), dtype=int)
    for d in range(1, half_window + 1):
        if d >= n_trials:
            break
        for t in range(n_trials - d):
            r = _pair_corr(fr[t], fr[t + d], min_bins)   # (n_units,)
            ok = np.isfinite(r)
            sums[t, ok] += r[ok]
            sums[t + d, ok] += r[ok]
            counts[t, ok] += 1
            counts[t + d, ok] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = sums / counts
    return np.where(counts > 0, rel, np.nan)


def epoch_mean_reliability(reliability: np.ndarray, trial_ids) -> np.ndarray:
    """Nan-mean per-unit reliability over 1-BASED trial ids.

    ``reliability`` is ``(n_trials, n_units)`` from
    :func:`trial_map_reliability`; ``trial_ids`` are the temporal stream's
    1-based trial ids (rows ``id - 1``). Returns ``(n_units,)``.
    """
    rows = np.asarray(list(trial_ids), dtype=int) - 1
    if rows.size and (rows.min() < 0 or rows.max() >= reliability.shape[0]):
        raise IndexError(
            f"trial ids {rows.min() + 1}..{rows.max() + 1} outside "
            f"1..{reliability.shape[0]}")
    with np.errstate(invalid="ignore"):
        return np.nanmean(reliability[rows], axis=0)
