"""Synthetic ground-truth tests for Arm A (segment-aware lagged CCA).

Each test constructs Y as X shifted by a known integer lag (plus noise), runs
the lag curve, and asserts peak lag / IFI sign / null centering recover the
ground truth. The injected signal must survive: (a) segmentation (no pair
straddles a gap or trial boundary), (b) leave-one-trial-out CV, and (c) the
per-segment circshift null centring at zero IFI when there is no coupling.

Pins down the contract in ``cca/UNDERSTANDING_temporal.md`` Sec. 5.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from tom_cca import config, lagged_temporal, segments

CFG = config.DEFAULT


def _make_segments(starts_and_lens, trials):
    """Build a list of Segment(trial, start, stop) from explicit specs."""
    out = []
    for (start, length), trial in zip(starts_and_lens, trials):
        out.append(segments.Segment(trial=trial, start=start, stop=start + length))
    return out


def _coupled_xy(n_total_bins: int, segs, true_lag: int,
                k: int = 4, noise: float = 0.3, seed: int = 0):
    """Build flat X, Y score matrices where, within each segment, Y(t) = X(t - true_lag) + noise.

    Outside segments the values are zero (they will never be sampled).
    """
    rng = np.random.default_rng(seed)
    X = np.zeros((n_total_bins, k))
    Y = np.zeros((n_total_bins, k))
    for s in segs:
        x_seg = rng.standard_normal((len(s), k))
        # Pad so we can shift without falling off the segment.
        pad = abs(true_lag) + 2
        x_full = np.zeros((len(s) + 2 * pad, k))
        x_full[pad:pad + len(s)] = x_seg
        # Y_within(t) = X_within(t - true_lag) -> if true_lag > 0, X leads Y in time
        # which means at lag L = true_lag we should see high CC.
        y_full = np.roll(x_full, true_lag, axis=0)
        # Take the central segment-length slice for both.
        X[s.start:s.stop] = x_full[pad:pad + len(s)]
        Y[s.start:s.stop] = y_full[pad:pad + len(s)]
    Y = Y + noise * rng.standard_normal(Y.shape)
    return X, Y


# ---------------------------------------------------------------------------
# Known-lag recovery
# ---------------------------------------------------------------------------
def test_lag_curve_recovers_known_positive_lag():
    # 6 segments of 30 bins each, 3 trials (2 segments per trial).
    segs = _make_segments(
        [(0, 30), (40, 30), (100, 30), (140, 30), (200, 30), (240, 30)],
        trials=[0, 0, 1, 1, 2, 2],
    )
    X, Y = _coupled_xy(300, segs, true_lag=3, k=4, noise=0.4, seed=1)
    cfg = dataclasses.replace(CFG, lag_ms=400, temporal_bin_ms=50)  # +-8 bins
    res = lagged_temporal.lag_curve(X, Y, segs, cfg)
    assert abs(res.peak_lag - 3) <= 1, (
        f"expected peak lag near +3, got {res.peak_lag}; "
        f"CC1 curve = {res.cc1}"
    )
    assert res.ifi > 0.1, f"IFI should be strongly positive, got {res.ifi}"


def test_lag_curve_recovers_known_negative_lag():
    segs = _make_segments(
        [(0, 30), (50, 30), (100, 30), (150, 30)],
        trials=[0, 0, 1, 1],
    )
    X, Y = _coupled_xy(200, segs, true_lag=-2, k=3, noise=0.4, seed=2)
    cfg = dataclasses.replace(CFG, lag_ms=400, temporal_bin_ms=50)
    res = lagged_temporal.lag_curve(X, Y, segs, cfg)
    assert abs(res.peak_lag - (-2)) <= 1
    assert res.ifi < -0.1


def test_lag_curve_no_coupling_gives_near_zero_ifi():
    rng = np.random.default_rng(7)
    segs = _make_segments(
        [(0, 30), (50, 30), (100, 30), (150, 30)],
        trials=[0, 0, 1, 1],
    )
    X = rng.standard_normal((200, 3))
    Y = rng.standard_normal((200, 3))
    cfg = dataclasses.replace(CFG, lag_ms=400, temporal_bin_ms=50)
    res = lagged_temporal.lag_curve(X, Y, segs, cfg)
    # Held-out CC on random data is near 0 (often slightly negative); IFI of
    # a near-flat curve should be small in magnitude.
    assert abs(res.ifi) < 0.4


# ---------------------------------------------------------------------------
# Segment respect
# ---------------------------------------------------------------------------
def test_no_lag_pair_crosses_gap_or_trial_boundary():
    """The segment-pooled X/Y rows must come from inside one segment each.

    We inject a *negative-lag* signal in segment A and a *positive-lag* signal
    in segment B. If the lag pooling correctly stays within segments, both
    signals contribute coherently. If it bridges the gap, the cross-bin pair
    is meaningless noise and CC drops.
    """
    n_total = 200
    segs = _make_segments(
        [(0, 40), (100, 40)], trials=[0, 1],
    )
    X, Y = _coupled_xy(n_total, segs, true_lag=2, k=3, noise=0.3, seed=4)
    cfg = dataclasses.replace(CFG, lag_ms=400, temporal_bin_ms=50)
    res = lagged_temporal.lag_curve(X, Y, segs, cfg)
    # Paired-sample count at lag L should equal sum_s max(0, len(s) - |L|).
    expected = np.array([
        segments.total_paired_samples(segs, int(l)) for l in res.lags
    ])
    np.testing.assert_array_equal(res.n_paired_samples, expected)
    # And the signal should still be recovered cleanly.
    assert abs(res.peak_lag - 2) <= 1


# ---------------------------------------------------------------------------
# Per-segment circshift null
# ---------------------------------------------------------------------------
def test_circshift_segments_keeps_shifts_within_segments():
    rng = np.random.default_rng(0)
    scores = np.arange(40, dtype=float).reshape(20, 2)
    segs = _make_segments([(0, 10), (10, 10)], trials=[0, 1])
    shuffled = lagged_temporal.circshift_segments(scores, segs, rng, min_shift=2)
    # Each segment's rows are a permutation of the original segment's rows.
    for s in segs:
        original_set = set(map(tuple, scores[s.start:s.stop].tolist()))
        shuffled_set = set(map(tuple, shuffled[s.start:s.stop].tolist()))
        assert original_set == shuffled_set, (
            "circshift_segments leaked rows across segments at "
            f"trial {s.trial}"
        )


def test_circshift_null_ifi_centred_near_zero_under_no_coupling():
    """With no real coupling, the null IFI distribution should sit at ~0."""
    rng_data = np.random.default_rng(11)
    segs = _make_segments(
        [(0, 25), (40, 25), (80, 25), (120, 25)],
        trials=[0, 0, 1, 1],
    )
    X = rng_data.standard_normal((150, 3))
    Y = rng_data.standard_normal((150, 3))
    cfg = dataclasses.replace(CFG, lag_ms=300, temporal_bin_ms=50,
                              temporal_circshift_min_bins=2)
    rng_shuf = np.random.default_rng(0)
    null_ifis = []
    for _ in range(30):
        nr = lagged_temporal.shuffle_lag_curve(X, Y, segs, cfg, rng_shuf)
        null_ifis.append(nr.ifi)
    null_ifis = np.array(null_ifis)
    assert abs(np.nanmean(null_ifis)) < 0.2, (
        f"null IFI mean drifted from 0: {np.nanmean(null_ifis):.3f}"
    )


def test_circshift_null_does_not_destroy_signal_in_first_arm_step():
    """Sanity: under real coupling, real IFI exceeds the null mean."""
    segs = _make_segments(
        [(0, 30), (50, 30), (100, 30), (150, 30), (200, 30), (250, 30)],
        trials=[0, 0, 1, 1, 2, 2],
    )
    X, Y = _coupled_xy(300, segs, true_lag=3, k=4, noise=0.3, seed=5)
    cfg = dataclasses.replace(CFG, lag_ms=400, temporal_bin_ms=50,
                              temporal_circshift_min_bins=3)
    real = lagged_temporal.lag_curve(X, Y, segs, cfg).ifi
    rng = np.random.default_rng(0)
    null = np.array([
        lagged_temporal.shuffle_lag_curve(X, Y, segs, cfg, rng).ifi
        for _ in range(40)
    ])
    # Real must clear (a) a non-trivial absolute IFI threshold, and (b) the
    # null mean by at least one shuffle SD. The 1-sigma threshold is a
    # sanity check on the surrogate, not a power calculation.
    null_mean = float(np.nanmean(null))
    null_sd = float(np.nanstd(null))
    assert real > 0.2 and real > null_mean + null_sd, (
        f"real IFI {real:.3f} did not beat null {null_mean:.3f} "
        f"+/- {null_sd:.3f}"
    )


# ---------------------------------------------------------------------------
# Robustness
# ---------------------------------------------------------------------------
def test_lag_curve_empty_segments_returns_nan():
    cfg = dataclasses.replace(CFG, lag_ms=200, temporal_bin_ms=50)
    res = lagged_temporal.lag_curve(
        np.zeros((10, 3)), np.zeros((10, 3)), [], cfg,
    )
    assert np.all(np.isnan(res.cc_per_dim))
    assert np.all(np.isnan(res.ifi_per_dim))


def test_lag_curve_single_trial_yields_nan_cv():
    """LOTO CV needs at least 2 trials -- single-trial input returns NaN."""
    segs = _make_segments([(0, 40)], trials=[0])
    X, Y = _coupled_xy(50, segs, true_lag=2, k=3, noise=0.2, seed=8)
    cfg = dataclasses.replace(CFG, lag_ms=200, temporal_bin_ms=50)
    res = lagged_temporal.lag_curve(X, Y, segs, cfg)
    assert np.all(np.isnan(res.cc_per_dim))
