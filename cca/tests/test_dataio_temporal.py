"""Tests for the temporal layer in tom_cca.dataio.

Covers:
* The 1 ms -> 50 ms binning helpers (rebin_spikes, bin_velocity,
  bin_trial_index) -- pure numpy, no I/O.
* An end-to-end integration check on the real TF028 file when present in
  ``HC_V1_data/`` (skipped otherwise). Confirms the lazy 1 ms stream loader
  composes correctly with the segment + landmark helpers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tom_cca import config, dataio, segments

CFG = config.DEFAULT


# ---------------------------------------------------------------------------
# rebin_spikes
# ---------------------------------------------------------------------------
def test_rebin_spikes_sums_within_bins():
    # 100 1 ms bins, 2 units. Unit 0 always 1, unit 1 always 2. At 50 ms,
    # each bin should hold sums of 50 and 100 respectively.
    spikes = np.zeros((100, 2), dtype=np.uint8)
    spikes[:, 0] = 1
    spikes[:, 1] = 2
    out = dataio.rebin_spikes(spikes, bin_ms=50)
    assert out.shape == (2, 2)
    np.testing.assert_array_equal(out[:, 0], [50, 50])
    np.testing.assert_array_equal(out[:, 1], [100, 100])
    assert out.dtype == np.float32


def test_rebin_spikes_drops_incomplete_final_bin():
    spikes = np.ones((125, 1), dtype=np.uint8)
    out = dataio.rebin_spikes(spikes, bin_ms=50)
    assert out.shape == (2, 1)        # 125 // 50 = 2; last 25 dropped
    np.testing.assert_array_equal(out[:, 0], [50, 50])


def test_rebin_spikes_empty():
    out = dataio.rebin_spikes(np.zeros((0, 3), dtype=np.uint8), bin_ms=50)
    assert out.shape == (0, 3)


# ---------------------------------------------------------------------------
# bin_velocity
# ---------------------------------------------------------------------------
def test_bin_velocity_takes_mean_per_bin():
    vel = np.concatenate([np.full(50, 1.0), np.full(50, 9.0)])
    out = dataio.bin_velocity(vel, bin_ms=50)
    np.testing.assert_array_equal(out, [1.0, 9.0])


def test_bin_velocity_handles_partial_bins_by_dropping():
    vel = np.full(75, 4.0)
    out = dataio.bin_velocity(vel, bin_ms=50)
    assert out.shape == (1,)
    assert out[0] == 4.0


# ---------------------------------------------------------------------------
# bin_trial_index
# ---------------------------------------------------------------------------
def test_bin_trial_index_all_same_returns_value():
    trial = np.full(100, 7.0)
    out = dataio.bin_trial_index(trial, bin_ms=50)
    np.testing.assert_array_equal(out, [7.0, 7.0])


def test_bin_trial_index_mixed_within_bin_returns_nan():
    # 50 ms bin spans 25 of trial 1 and 25 of trial 2 -> NaN.
    trial = np.concatenate([np.full(25, 1.0), np.full(25, 2.0),
                            np.full(50, 2.0)])
    out = dataio.bin_trial_index(trial, bin_ms=50)
    assert np.isnan(out[0])
    assert out[1] == 2.0


def test_bin_trial_index_any_nan_in_bin_returns_nan():
    trial = np.full(100, 3.0)
    trial[10] = np.nan
    out = dataio.bin_trial_index(trial, bin_ms=50)
    assert np.isnan(out[0])           # bin 0 has a NaN inside
    assert out[1] == 3.0              # bin 1 untouched


# ---------------------------------------------------------------------------
# Integration: lazy stream loader composed with segments + landmark_align
# ---------------------------------------------------------------------------
_REAL_TF = config.DATA_DIR / "TF028_export.mat"


@pytest.mark.skipif(not _REAL_TF.is_file(),
                    reason="TF028_export.mat not mounted in this environment")
def test_temporal_streams_load_and_segment_on_real_animal():
    animal = dataio.load_animal(_REAL_TF, animal_id=28,
                                region_labels=None) if False else None
    # Region labels live in cca_labels.json (Tom's exports use MATLAB `string`
    # which h5py cannot decode). Use the cohort loader so the companion is
    # consulted -- but only if the companion is present.
    companion = (config.DATA_DIR / config.COMPANION_FILE).is_file()
    if not companion:
        pytest.skip("cca_labels.json companion missing")
    animals = dataio.load_animals(config.DATA_DIR)
    animal = next((a for a in animals if a.animal_id == 28), None)
    if animal is None:
        pytest.skip("TF028 not in cohort")

    # Segments
    segs = dataio.temporal_segments(animal, CFG)
    assert len(segs) > 0, "expected at least one running-state segment"
    # Each segment respects the min-length constraint.
    mr = config.min_run_bins(CFG)
    assert all(len(s) >= mr for s in segs)

    # Area activity matches segment range.
    activity, idx = dataio.area_activity_50ms(animal, "CA1", CFG)
    assert activity.shape[1] == idx.size
    assert activity.shape[1] >= CFG.min_units
    # The flat 50 ms timeline length matches the binning -- spot-check that
    # the last segment fits inside the activity array.
    assert segs[-1].stop <= activity.shape[0]

    # Landmark windows: at least one landmark id should survive the gate.
    win_dict, idx2 = dataio.area_landmark_windows(animal, "CA1", CFG)
    assert len(win_dict) > 0
    win_len = config.landmark_window_bins(CFG)
    for lid, lw in win_dict.items():
        assert lw.spikes.shape == (lw.n_crossings, win_len, idx2.size)
        assert lw.n_crossings >= 1
        assert np.all(lw.mean_velocities >= CFG.velocity_thresh_cm_s)
