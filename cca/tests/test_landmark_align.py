"""Tests for tom_cca.landmark_align -- entry detection + window alignment.

The contract these tests pin down is the spec in
``cca/UNDERSTANDING_temporal.md`` Sec. 6.
"""

from __future__ import annotations

import numpy as np
import pytest

from tom_cca import landmark_align as la


# ---------------------------------------------------------------------------
# find_entries
# ---------------------------------------------------------------------------
def test_find_entries_detects_zero_to_nonzero_transitions():
    # Three landmark visits, ids 1, 2, 1
    lid = np.array([0, 0, 1, 1, 1, 0, 0, 2, 2, 0, 0, 1, 1])
    entries = la.find_entries(lid)
    # Entries at indices 2 (id 1), 7 (id 2), 11 (id 1).
    assert entries.shape == (3, 2)
    assert list(entries[:, 0]) == [2, 7, 11]
    assert list(entries[:, 1]) == [1, 2, 1]


def test_find_entries_does_not_count_nonzero_to_other_nonzero():
    # Landmark 1 directly transitions to landmark 2 with no zero gap.
    # That is NOT a new entry to landmark 2.
    lid = np.array([0, 1, 1, 2, 2, 0])
    entries = la.find_entries(lid)
    assert entries.shape == (1, 2)
    assert entries[0, 0] == 1 and entries[0, 1] == 1


def test_find_entries_first_bin_nonzero_counts_as_entry():
    # The implicit pre-trace is treated as zero, so a non-zero first bin is
    # an entry.
    lid = np.array([3, 3, 0, 0, 4])
    entries = la.find_entries(lid)
    assert list(entries[:, 0]) == [0, 4]
    assert list(entries[:, 1]) == [3, 4]


def test_find_entries_handles_all_zero():
    lid = np.zeros(20, dtype=int)
    entries = la.find_entries(lid)
    assert entries.shape == (0, 2)


def test_find_entries_empty_input():
    entries = la.find_entries(np.array([], dtype=int))
    assert entries.shape == (0, 2)


# ---------------------------------------------------------------------------
# align_windows
# ---------------------------------------------------------------------------
def _setup(n_bins=80, n_units=4, bin_ms=50):
    """Build a synthetic animal: 80 50 ms bins, 4 units. Trial 0 covers bins
    0..40, trial 1 covers bins 40..80; uniformly running at 5 cm/s."""
    spikes = np.arange(n_bins * n_units, dtype=float).reshape(n_bins, n_units)
    vel = np.full(n_bins, 5.0)
    trial = np.concatenate([np.zeros(40), np.ones(40)]).astype(float)
    return spikes, vel, trial, bin_ms


def test_align_windows_basic_slicing_and_grouping():
    spikes, vel, trial, bin_ms = _setup()
    # Two entries: landmark 1 at 1 ms index 500 (= bin 10), landmark 2 at
    # 1 ms index 1500 (= bin 30).
    entries = np.array([[500, 1], [1500, 2]])
    out = la.align_windows(spikes, vel, trial, entries,
                           bin_ms=bin_ms,
                           window_bins_pre=5, window_bins_post=4,
                           vel_threshold_cm_s=2.0)
    assert set(out.keys()) == {1, 2}
    # window length = 5 + 1 + 4 = 10 bins.
    assert out[1].spikes.shape == (1, 10, 4)
    # Slice for entry at bin 10: bins [10-5, 10+1+4) = [5, 15).
    np.testing.assert_array_equal(out[1].spikes[0], spikes[5:15, :])
    np.testing.assert_array_equal(out[2].spikes[0], spikes[25:35, :])
    assert out[1].trials[0] == 0
    assert out[2].trials[0] == 0


def test_align_windows_drops_windows_that_overrun_array():
    spikes, vel, trial, bin_ms = _setup(n_bins=20)
    # Entry at 1 ms index 100 (= bin 2). pre=5 -> start at -3, off the start.
    entries = np.array([[100, 1]])
    out = la.align_windows(spikes, vel, trial[:20], entries,
                           bin_ms=bin_ms,
                           window_bins_pre=5, window_bins_post=4,
                           vel_threshold_cm_s=2.0)
    assert out == {}


def test_align_windows_drops_windows_spanning_trial_boundary():
    spikes, vel, trial, bin_ms = _setup()
    # Trial flips at bin 40; entry at bin 38 means window [33, 43) spans both
    # trials -> dropped.
    entries = np.array([[38 * 50, 1]])
    out = la.align_windows(spikes, vel, trial, entries,
                           bin_ms=bin_ms,
                           window_bins_pre=5, window_bins_post=4,
                           vel_threshold_cm_s=2.0)
    assert out == {}


def test_align_windows_drops_window_below_velocity_threshold():
    spikes, vel, trial, bin_ms = _setup()
    # Make bins 5..15 slow (< threshold) so the window mean dips below 2 cm/s.
    vel[5:15] = 0.5
    entries = np.array([[500, 1]])   # entry at bin 10 -> window [5, 15)
    out = la.align_windows(spikes, vel, trial, entries,
                           bin_ms=bin_ms,
                           window_bins_pre=5, window_bins_post=4,
                           vel_threshold_cm_s=2.0)
    assert out == {}


def test_align_windows_per_window_gate_is_mean_not_min():
    spikes, vel, trial, bin_ms = _setup()
    # Most bins slow, a few fast. Mean should be above threshold -> kept.
    vel[5:15] = 0.0
    vel[7] = 30.0  # one fast bin lifts the mean to 3.0
    entries = np.array([[500, 1]])
    out = la.align_windows(spikes, vel, trial, entries,
                           bin_ms=bin_ms,
                           window_bins_pre=5, window_bins_post=4,
                           vel_threshold_cm_s=2.0)
    assert 1 in out and out[1].n_crossings == 1


def test_align_windows_stacks_multiple_crossings_per_landmark():
    spikes, vel, trial, bin_ms = _setup(n_bins=120)
    trial = np.concatenate([np.zeros(60), np.ones(60)]).astype(float)
    # Two visits to landmark 1: bins 10 and 50 (both in trial 0).
    entries = np.array([[10 * 50, 1], [50 * 50, 1]])
    out = la.align_windows(spikes, vel, trial, entries,
                           bin_ms=bin_ms,
                           window_bins_pre=5, window_bins_post=4,
                           vel_threshold_cm_s=2.0)
    assert out[1].n_crossings == 2
    assert out[1].spikes.shape == (2, 10, 4)
    assert list(out[1].trials) == [0, 0]


def test_align_windows_drops_window_with_any_nan_bin():
    spikes, vel, trial, bin_ms = _setup()
    # NaN inside trial_idx for one bin in the window.
    trial[8] = np.nan
    entries = np.array([[500, 1]])
    out = la.align_windows(spikes, vel, trial, entries,
                           bin_ms=bin_ms,
                           window_bins_pre=5, window_bins_post=4,
                           vel_threshold_cm_s=2.0)
    assert out == {}


def test_align_windows_raises_on_shape_mismatch():
    spikes = np.zeros((10, 3))
    vel = np.zeros(11)
    trial = np.zeros(10)
    with pytest.raises(ValueError):
        la.align_windows(spikes, vel, trial, np.zeros((0, 2), dtype=int),
                         bin_ms=50, window_bins_pre=2, window_bins_post=2,
                         vel_threshold_cm_s=1.0)
