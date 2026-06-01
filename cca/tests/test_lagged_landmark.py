"""Synthetic ground-truth tests for Arm B (window-bounded lagged CCA).

Each test builds a per-landmark window stack ``(n_crossings, win_len, k)``
where Y(bin) = X(bin - true_lag) + noise inside each window, runs the lag
curve, and checks peak lag / IFI sign. Also checks the per-window circshift
null centres at zero IFI under no coupling, and that the lag scan is clipped
to ``window - 5`` when the requested lag would exceed the window.

Pins down the contract in ``cca/UNDERSTANDING_temporal.md`` Sec. 6.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from tom_cca import config, lagged_landmark

CFG = config.DEFAULT


def _coupled_windows(n_crossings, win_len, true_lag, k=3, noise=0.3,
                     trials_per_crossing=None, seed=0):
    """Build per-window X, Y tensors with Y a lagged copy of X."""
    rng = np.random.default_rng(seed)
    pad = abs(true_lag) + 2
    X = np.zeros((n_crossings, win_len, k))
    Y = np.zeros((n_crossings, win_len, k))
    for c in range(n_crossings):
        x_full = rng.standard_normal((win_len + 2 * pad, k))
        y_full = np.roll(x_full, true_lag, axis=0)
        X[c] = x_full[pad:pad + win_len]
        Y[c] = y_full[pad:pad + win_len]
    Y = Y + noise * rng.standard_normal(Y.shape)
    if trials_per_crossing is None:
        trials_per_crossing = np.arange(n_crossings) % max(3, n_crossings // 4)
    return X, Y, np.asarray(trials_per_crossing, dtype=int)


# ---------------------------------------------------------------------------
# Known-lag recovery
# ---------------------------------------------------------------------------
def test_lag_curve_recovers_known_positive_lag_landmark_arm():
    # 30 crossings of 15-bin windows; lag = +2 bins.
    X, Y, trials = _coupled_windows(30, 15, true_lag=2, k=3, noise=0.4, seed=1)
    cfg = dataclasses.replace(CFG, landmark_lag_ms=400, temporal_bin_ms=50)
    res = lagged_landmark.lag_curve(X, Y, trials, landmark_id=1, cfg=cfg)
    assert abs(res.peak_lag - 2) <= 1
    assert res.ifi > 0.1


def test_lag_curve_recovers_known_negative_lag_landmark_arm():
    X, Y, trials = _coupled_windows(30, 15, true_lag=-3, k=3, noise=0.4, seed=2)
    cfg = dataclasses.replace(CFG, landmark_lag_ms=400, temporal_bin_ms=50)
    res = lagged_landmark.lag_curve(X, Y, trials, landmark_id=2, cfg=cfg)
    assert abs(res.peak_lag - (-3)) <= 1
    assert res.ifi < -0.1


def test_lag_curve_no_coupling_gives_near_zero_ifi():
    rng = np.random.default_rng(9)
    X = rng.standard_normal((25, 12, 3))
    Y = rng.standard_normal((25, 12, 3))
    trials = np.arange(25) % 5
    cfg = dataclasses.replace(CFG, landmark_lag_ms=300, temporal_bin_ms=50)
    res = lagged_landmark.lag_curve(X, Y, trials, landmark_id=3, cfg=cfg)
    assert abs(res.ifi) < 0.4


# ---------------------------------------------------------------------------
# Window-length bounds the lag scan
# ---------------------------------------------------------------------------
def test_lag_scan_clipped_to_window_minus_floor():
    """A 10-bin window with a >= 5 paired-sample floor caps |L| at 5."""
    X = np.zeros((10, 10, 3))
    Y = np.zeros((10, 10, 3))
    trials = np.arange(10) % 3
    cfg = dataclasses.replace(CFG, landmark_lag_ms=1000, temporal_bin_ms=50)
    res = lagged_landmark.lag_curve(X, Y, trials, landmark_id=1, cfg=cfg)
    assert res.lags.min() == -5 and res.lags.max() == 5


def test_short_window_returns_only_zero_lag():
    X = np.zeros((10, 5, 2))
    Y = np.zeros((10, 5, 2))
    trials = np.arange(10) % 3
    cfg = dataclasses.replace(CFG, landmark_lag_ms=200, temporal_bin_ms=50)
    res = lagged_landmark.lag_curve(X, Y, trials, landmark_id=1, cfg=cfg)
    # win=5 -> max_allowed = 0 -> lags = [0] only.
    assert list(res.lags) == [0]


# ---------------------------------------------------------------------------
# Per-window circshift null
# ---------------------------------------------------------------------------
def test_circshift_windows_keeps_shift_within_window():
    rng = np.random.default_rng(0)
    scores = np.tile(np.arange(10).reshape(10, 1), (1, 1))[None]  # (1, 10, 1)
    scores = np.broadcast_to(scores, (3, 10, 1)).copy().astype(float)
    shuffled = lagged_landmark.circshift_windows(scores, rng, min_shift=1)
    for c in range(3):
        # Each window's values are a permutation of [0..9].
        assert sorted(shuffled[c, :, 0].tolist()) == list(range(10))


def test_circshift_null_centred_near_zero_under_no_coupling():
    rng_data = np.random.default_rng(13)
    X = rng_data.standard_normal((25, 12, 3))
    Y = rng_data.standard_normal((25, 12, 3))
    trials = np.arange(25) % 5
    cfg = dataclasses.replace(CFG, landmark_lag_ms=300, temporal_bin_ms=50,
                              landmark_circshift_min_bins=1)
    rng_shuf = np.random.default_rng(0)
    null_ifis = []
    for _ in range(30):
        nr = lagged_landmark.shuffle_lag_curve(X, Y, trials, 1, cfg, rng_shuf)
        null_ifis.append(nr.ifi)
    null_ifis = np.array(null_ifis)
    assert abs(np.nanmean(null_ifis)) < 0.2


def test_real_ifi_beats_null_under_real_coupling():
    X, Y, trials = _coupled_windows(30, 15, true_lag=2, k=3, noise=0.3, seed=6)
    cfg = dataclasses.replace(CFG, landmark_lag_ms=300, temporal_bin_ms=50,
                              landmark_circshift_min_bins=4)
    real = lagged_landmark.lag_curve(X, Y, trials, 1, cfg).ifi
    rng = np.random.default_rng(0)
    null = np.array([
        lagged_landmark.shuffle_lag_curve(X, Y, trials, 1, cfg, rng).ifi
        for _ in range(40)
    ])
    null_mean = float(np.nanmean(null))
    null_sd = float(np.nanstd(null))
    assert real > 0.2 and real > null_mean + null_sd, (
        f"real IFI {real:.3f} did not beat null {null_mean:.3f} "
        f"+/- {null_sd:.3f}"
    )
