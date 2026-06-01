"""Segment construction for Arm A (running-state temporal CCA).

Builds the segment list that gates the velocity-filtered lag analysis. Inputs
are the 50 ms-binned velocity mask and the per-bin trial index (NaN outside
the cued task). Outputs are the maximal contiguous valid runs **within a
single cued trial**, after morphological cleanup, that are long enough to
yield >= 5 paired samples at the most extreme lag.

A segment never crosses a trial boundary even if the velocity stays high
across the inter-trial gap -- the gap is a teleport/reset, not contiguous
neural data. See ``cca/UNDERSTANDING_temporal.md`` Sec. 5.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Segment:
    """One maximal contiguous valid run in a single trial.

    ``start`` is inclusive, ``stop`` is exclusive (slice-friendly). Both index
    into the 50 ms binned timeline.
    """

    trial: int
    start: int
    stop: int

    def __len__(self) -> int:
        return self.stop - self.start

    def slice(self) -> slice:
        return slice(self.start, self.stop)


def close_short_gaps(mask: np.ndarray, min_gap_close_bins: int) -> np.ndarray:
    """Fill any interior invalid run of length <= ``min_gap_close_bins``.

    "Interior" means the run is flanked by valid bins on both sides. Edge
    gaps (at the start or end of the trace) are left alone -- there is no
    flanking signal to interpolate over.
    """
    mask = np.asarray(mask, dtype=bool).copy()
    n = mask.size
    if n == 0 or min_gap_close_bins <= 0:
        return mask
    i = 0
    while i < n:
        if mask[i]:
            i += 1
            continue
        j = i
        while j < n and not mask[j]:
            j += 1
        gap_len = j - i
        # Only fill interior gaps: must have a valid bin on both sides.
        if i > 0 and j < n and gap_len <= min_gap_close_bins:
            mask[i:j] = True
        i = j
    return mask


def find_segments(mask: np.ndarray,
                  trial_idx: np.ndarray,
                  min_gap_close_bins: int,
                  min_run_bins: int) -> list[Segment]:
    """Find Segment runs after gap-closing, respecting trial boundaries.

    Parameters
    ----------
    mask : (n_bins,) bool
        Velocity-passing mask in the 50 ms grid.
    trial_idx : (n_bins,) float
        Trial index per 50 ms bin; NaN where the bin is outside any cued trial
        or spans a trial boundary.
    min_gap_close_bins : int
        Invalid runs <= this many bins are bridged (morphological close)
        before segmentation. Pass 0 to disable.
    min_run_bins : int
        Final runs strictly shorter than this are dropped.

    Returns
    -------
    list of :class:`Segment`
        Each segment carries the trial id and the inclusive/exclusive bin
        range. Trial boundaries always split, even if the closed mask would
        otherwise span them.
    """
    mask = np.asarray(mask, dtype=bool)
    trial_idx = np.asarray(trial_idx, dtype=float)
    if mask.shape != trial_idx.shape:
        raise ValueError(
            f"mask {mask.shape} and trial_idx {trial_idx.shape} must match"
        )
    if mask.ndim != 1:
        raise ValueError(f"expected 1-D arrays, got {mask.ndim}-D")

    closed = close_short_gaps(mask, min_gap_close_bins)
    n = closed.size
    segments: list[Segment] = []
    i = 0
    while i < n:
        # Skip invalid or no-trial bins.
        if not closed[i] or np.isnan(trial_idx[i]):
            i += 1
            continue
        # Walk forward while valid AND same trial.
        j = i
        t = trial_idx[i]
        while (j < n and closed[j]
               and not np.isnan(trial_idx[j])
               and trial_idx[j] == t):
            j += 1
        if (j - i) >= min_run_bins:
            segments.append(Segment(trial=int(t), start=i, stop=j))
        i = j
    return segments


def segment_paired_indices(segment: Segment, lag: int) -> tuple[np.ndarray, np.ndarray]:
    """Indices into ``segment`` for the lag-L pair list.

    At lag L > 0: pair bin t (in X) with bin t+L (in Y) -> X leads Y.
    At lag L < 0: pair bin t+|L| (in X) with bin t (in Y) -> Y leads X.

    Returns the two index arrays (X-bin indices, Y-bin indices) relative to
    the start of the segment, with both arrays the same length.
    """
    seg_len = len(segment)
    L = int(lag)
    if abs(L) >= seg_len:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    if L >= 0:
        x = np.arange(0, seg_len - L)
        y = np.arange(L, seg_len)
    else:
        x = np.arange(-L, seg_len)
        y = np.arange(0, seg_len + L)
    return x, y


def total_paired_samples(segments: list[Segment], lag: int) -> int:
    """Count the lag-L paired samples summed across all segments."""
    L = int(abs(lag))
    return sum(max(0, len(s) - L) for s in segments)
