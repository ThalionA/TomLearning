"""Landmark-entry detection + window alignment for Arm B.

The landmark arm takes 500 ms windows centred on each landmark entry event
(bins -5..+4, entry at index 5; spans -250 ms to +200 ms at 50 ms bins) and
fits CCA per ``(animal, pair, epoch, landmark_id)``. See
``cca/UNDERSTANDING_temporal.md`` Sec. 6.

Two stages:

1. :func:`find_entries` walks the 1 ms ``visual_landmark_binned`` stream and
   emits ``(entry_1ms_index, landmark_id)`` at every 0 -> non-zero transition.
   X -> Y transitions (different landmarks with no zero gap in between) are
   not counted as entries to Y -- only 0 -> Y transitions count.

2. :func:`align_windows` maps each entry to its 50 ms entry-bin, slices a
   window of ``(pre + 1 + post)`` bins around it, applies the per-window
   engagement gate (drop if mean window velocity < threshold), enforces the
   single-trial constraint (a window may not span a trial boundary), and
   groups the surviving windows by landmark_id.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LandmarkWindows:
    """All retained windows for one landmark id within an animal."""

    landmark_id: int
    spikes: np.ndarray           # (n_kept, n_window_bins, n_units) float
    trials: np.ndarray           # (n_kept,) int -- the trial each crossing came from
    entry_bins: np.ndarray       # (n_kept,) int -- 50 ms-bin index of entry
    mean_velocities: np.ndarray  # (n_kept,) float -- mean window velocity, cm/s

    @property
    def n_crossings(self) -> int:
        return self.spikes.shape[0]


def find_entries(landmark_id_1ms: np.ndarray) -> np.ndarray:
    """Detect 0 -> non-zero transitions in the 1 ms landmark stream.

    Parameters
    ----------
    landmark_id_1ms : (n_1ms,) array
        Landmark identity per 1 ms bin: 0 outside any landmark zone, otherwise
        the landmark id.

    Returns
    -------
    (n_entries, 2) int array
        Columns: (entry_1ms_index, landmark_id). Length 0 if there are no
        entries.
    """
    lid = np.asarray(landmark_id_1ms).ravel()
    if lid.size == 0:
        return np.zeros((0, 2), dtype=int)
    # An entry bin is non-zero AND the previous bin was zero. The first bin
    # counts as an entry if it is non-zero (treat the implicit pre-trace as
    # zero, matching the "before the recording started" intuition).
    prev_zero = np.concatenate(([True], lid[:-1] == 0))
    is_entry = (lid != 0) & prev_zero
    idxs = np.where(is_entry)[0]
    out = np.column_stack([idxs.astype(int), lid[idxs].astype(int)])
    return out


def _bin_index(idx_1ms: int, bin_ms: int) -> int:
    """1 ms index -> 50 ms bin index (integer floor)."""
    return int(idx_1ms) // int(bin_ms)


def align_windows(rebinned_spikes: np.ndarray,
                  vel_50ms: np.ndarray,
                  trial_idx_50ms: np.ndarray,
                  entries_1ms: np.ndarray,
                  bin_ms: int,
                  window_bins_pre: int,
                  window_bins_post: int,
                  vel_threshold_cm_s: float) -> dict[int, LandmarkWindows]:
    """Slice + gate landmark-centred windows; return per-landmark stacks.

    Parameters
    ----------
    rebinned_spikes : (n_bins_50ms, n_units) float
        50 ms-binned spike-count or rate matrix for the animal.
    vel_50ms : (n_bins_50ms,) float
        50 ms-binned mean velocity, cm/s.
    trial_idx_50ms : (n_bins_50ms,) float
        Per-bin trial index; NaN outside the cued task or at a trial-spanning
        bin (the all-or-nothing binning rule -- a 50 ms bin is in trial t only
        if every 1 ms bin in its span is in trial t).
    entries_1ms : (n_entries, 2) int
        Output of :func:`find_entries`.
    bin_ms : int
    window_bins_pre, window_bins_post : int
        Window spans ``[entry - pre, entry + 1 + post)`` in 50 ms bins.
    vel_threshold_cm_s : float
        Per-window mean velocity gate -- drop windows whose mean velocity is
        strictly below this.

    Returns
    -------
    dict landmark_id -> :class:`LandmarkWindows`
        Only landmarks with at least one kept crossing appear in the dict.
    """
    rebinned_spikes = np.asarray(rebinned_spikes, dtype=float)
    vel_50ms = np.asarray(vel_50ms, dtype=float)
    trial_idx_50ms = np.asarray(trial_idx_50ms, dtype=float)
    n_bins, n_units = rebinned_spikes.shape
    if vel_50ms.shape != (n_bins,) or trial_idx_50ms.shape != (n_bins,):
        raise ValueError(
            f"shape mismatch: spikes ({n_bins},{n_units}) vs vel {vel_50ms.shape} "
            f"vs trial {trial_idx_50ms.shape}"
        )
    pre, post = int(window_bins_pre), int(window_bins_post)
    win_len = pre + 1 + post

    buckets: dict[int, list[tuple[np.ndarray, int, int, float]]] = defaultdict(list)
    for row in np.asarray(entries_1ms, dtype=int):
        idx_1ms, lid = int(row[0]), int(row[1])
        entry_bin = _bin_index(idx_1ms, bin_ms)
        start, stop = entry_bin - pre, entry_bin + 1 + post
        if start < 0 or stop > n_bins:
            continue
        win_trials = trial_idx_50ms[start:stop]
        if np.any(np.isnan(win_trials)):
            continue
        # All bins in the window must come from the same trial.
        t0 = win_trials[0]
        if not np.all(win_trials == t0):
            continue
        win_vel = vel_50ms[start:stop]
        if np.any(np.isnan(win_vel)):
            continue
        mean_v = float(np.mean(win_vel))
        if mean_v < vel_threshold_cm_s:
            continue
        win_spikes = rebinned_spikes[start:stop, :]
        buckets[lid].append((win_spikes, int(t0), entry_bin, mean_v))

    out: dict[int, LandmarkWindows] = {}
    for lid, items in buckets.items():
        spikes = np.stack([it[0] for it in items], axis=0)  # (n, win_len, n_units)
        assert spikes.shape == (len(items), win_len, n_units)
        out[lid] = LandmarkWindows(
            landmark_id=lid,
            spikes=spikes,
            trials=np.array([it[1] for it in items], dtype=int),
            entry_bins=np.array([it[2] for it in items], dtype=int),
            mean_velocities=np.array([it[3] for it in items], dtype=float),
        )
    return out
