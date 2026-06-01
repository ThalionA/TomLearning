"""Tests for tom_cca.segments -- velocity-mask cleanup + segment construction.

The contract these tests pin down is the spec in
``cca/UNDERSTANDING_temporal.md`` Sec. 5.
"""

from __future__ import annotations

import numpy as np
import pytest

from tom_cca import segments


# ---------------------------------------------------------------------------
# close_short_gaps
# ---------------------------------------------------------------------------
def _mask(s: str) -> np.ndarray:
    """Compact mask helper: '.' = False, '#' = True."""
    return np.array([c == "#" for c in s], dtype=bool)


def test_close_gaps_fills_short_interior_runs():
    # gap of 1 between valid runs -> filled
    m = _mask("###.###")
    out = segments.close_short_gaps(m, min_gap_close_bins=2)
    assert np.array_equal(out, _mask("#######"))


def test_close_gaps_leaves_long_interior_runs_alone():
    # gap of 3 with min_gap_close=2 -> untouched
    m = _mask("###...###")
    out = segments.close_short_gaps(m, min_gap_close_bins=2)
    assert np.array_equal(out, m)


def test_close_gaps_leaves_edge_gaps_alone():
    # leading and trailing invalid runs have no flanking valid -> untouched
    m = _mask("...###...")
    out = segments.close_short_gaps(m, min_gap_close_bins=10)
    assert np.array_equal(out, m)


def test_close_gaps_min_zero_is_noop():
    m = _mask("###.###")
    out = segments.close_short_gaps(m, min_gap_close_bins=0)
    assert np.array_equal(out, m)


def test_close_gaps_idempotent():
    m = _mask("##.##.##")
    out1 = segments.close_short_gaps(m, min_gap_close_bins=1)
    out2 = segments.close_short_gaps(out1, min_gap_close_bins=1)
    assert np.array_equal(out1, out2)


# ---------------------------------------------------------------------------
# find_segments
# ---------------------------------------------------------------------------
def test_find_segments_single_trial_one_segment():
    mask = _mask("#" * 20)
    trial = np.zeros(20)
    segs = segments.find_segments(mask, trial,
                                  min_gap_close_bins=0, min_run_bins=5)
    assert len(segs) == 1
    assert segs[0] == segments.Segment(trial=0, start=0, stop=20)


def test_find_segments_drops_short_runs():
    # Three runs of lengths 5, 3, 8; min_run=5 -> two kept.
    mask = _mask("#####...###...########")
    trial = np.zeros(mask.size)
    segs = segments.find_segments(mask, trial,
                                  min_gap_close_bins=0, min_run_bins=5)
    assert [len(s) for s in segs] == [5, 8]
    assert segs[0].start == 0 and segs[0].stop == 5
    assert segs[1].start == 14 and segs[1].stop == 22


def test_find_segments_splits_at_trial_boundary_even_when_mask_continuous():
    # Mask is uniformly valid but trial flips midway.
    mask = _mask("#" * 20)
    trial = np.concatenate([np.zeros(10), np.ones(10)])
    segs = segments.find_segments(mask, trial,
                                  min_gap_close_bins=0, min_run_bins=5)
    assert len(segs) == 2
    assert segs[0].trial == 0 and (segs[0].start, segs[0].stop) == (0, 10)
    assert segs[1].trial == 1 and (segs[1].start, segs[1].stop) == (10, 20)


def test_find_segments_skips_nan_trial_bins():
    # Inter-trial NaN gap breaks the segment cleanly.
    mask = _mask("#" * 20)
    trial = np.concatenate([np.zeros(8), [np.nan, np.nan], np.ones(10)])
    segs = segments.find_segments(mask, trial,
                                  min_gap_close_bins=0, min_run_bins=5)
    assert len(segs) == 2
    assert segs[0].trial == 0 and (segs[0].start, segs[0].stop) == (0, 8)
    assert segs[1].trial == 1 and (segs[1].start, segs[1].stop) == (10, 20)


def test_find_segments_gap_close_does_not_bridge_trial_boundary():
    # Short invalid gap that, if naively closed, would span a trial boundary.
    # Gap is 1 bin wide (length 1), small enough for min_gap_close=2 to
    # consider closing -- but the trial flips inside the gap, so the segment
    # walk must still split.
    mask = _mask("#######.#######")          # gap of 1 at index 7
    trial = np.concatenate([np.zeros(8), np.ones(7)])  # flips at index 8
    segs = segments.find_segments(mask, trial,
                                  min_gap_close_bins=2, min_run_bins=4)
    # Closing fills index 7 with True, but trial walk still splits at the
    # trial change between bins 7 and 8.
    assert len(segs) == 2
    assert segs[0].trial == 0
    assert segs[1].trial == 1


def test_find_segments_no_valid_bins_returns_empty():
    mask = _mask("." * 20)
    trial = np.zeros(20)
    segs = segments.find_segments(mask, trial,
                                  min_gap_close_bins=0, min_run_bins=5)
    assert segs == []


def test_find_segments_min_run_enforced_after_close():
    # Bridged run length = 5 + 1 + 5 = 11; passes min_run=10.
    mask = _mask("#####.#####")
    trial = np.zeros(mask.size)
    segs = segments.find_segments(mask, trial,
                                  min_gap_close_bins=2, min_run_bins=10)
    assert len(segs) == 1
    assert len(segs[0]) == 11


def test_find_segments_raises_on_shape_mismatch():
    with pytest.raises(ValueError):
        segments.find_segments(_mask("###"), np.zeros(5),
                               min_gap_close_bins=0, min_run_bins=1)


# ---------------------------------------------------------------------------
# segment_paired_indices + total_paired_samples
# ---------------------------------------------------------------------------
def test_paired_indices_positive_lag_pairs_xt_with_yt_plus_lag():
    seg = segments.Segment(trial=0, start=100, stop=110)  # length 10
    x, y = segments.segment_paired_indices(seg, lag=3)
    # X-indices 0..6, Y-indices 3..9 -- both length 7.
    assert np.array_equal(x, np.arange(0, 7))
    assert np.array_equal(y, np.arange(3, 10))


def test_paired_indices_negative_lag_is_symmetric():
    seg = segments.Segment(trial=0, start=0, stop=10)
    x_pos, y_pos = segments.segment_paired_indices(seg, lag=3)
    x_neg, y_neg = segments.segment_paired_indices(seg, lag=-3)
    # Negative lag swaps X/Y roles: X starts later, Y starts earlier.
    assert np.array_equal(x_neg, y_pos)
    assert np.array_equal(y_neg, x_pos)


def test_paired_indices_lag_geq_segment_returns_empty():
    seg = segments.Segment(trial=0, start=0, stop=5)
    x, y = segments.segment_paired_indices(seg, lag=5)
    assert x.size == 0 and y.size == 0


def test_total_paired_samples_matches_sum_over_segments():
    segs = [segments.Segment(0, 0, 10), segments.Segment(0, 20, 28),
            segments.Segment(1, 50, 65)]
    # lengths: 10, 8, 15; at |L|=3: contributions 7, 5, 12 -> 24
    assert segments.total_paired_samples(segs, lag=3) == 24
    assert segments.total_paired_samples(segs, lag=-3) == 24
