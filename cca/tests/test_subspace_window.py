"""Tests for tom_cca.subspace_window (full per-window subspace readout)."""

from __future__ import annotations

import numpy as np

from tom_cca import subspace_window as sw


def _grouped(n, n_groups, rng):
    return rng.integers(0, n_groups, n)


def test_recovers_shared_dims_and_significance():
    rng = np.random.default_rng(0)
    n = 6000
    lat = rng.normal(0, 1, (n, 2))               # two shared latents
    X = rng.normal(0, 1, (n, 8))
    Y = rng.normal(0, 1, (n, 8))
    X[:, :2] += 1.5 * lat
    Y[:, :2] += 1.5 * lat
    g = _grouped(n, 10, rng)
    r = sw.window_subspace(X, Y, g, k=6, max_lag=5, n_shuffles=30)
    assert r.cc[0] > 0.5                          # dominant CC recovered
    assert r.n_sig >= 2                           # two significant dims
    assert r.mi_sig > 0
    assert r.weights_x.shape[0] == 8 and r.member_x.dtype == bool


def test_independent_has_no_significant_dims():
    rng = np.random.default_rng(1)
    n = 6000
    X = rng.normal(0, 1, (n, 8))
    Y = rng.normal(0, 1, (n, 8))
    g = _grouped(n, 10, rng)
    r = sw.window_subspace(X, Y, g, k=6, max_lag=5, n_shuffles=40)
    assert r.n_sig <= 1                           # ~0 false positives
    assert r.mi_sig < 0.05 or r.n_sig == 0


def test_optimal_lag_recovers_known_shift():
    # Y is X's shared latent delayed by +3 bins -> optimal lag should be near +3
    rng = np.random.default_rng(2)
    n = 8000
    lat = rng.normal(0, 1, n)
    X = rng.normal(0, 0.3, (n, 5))
    Y = rng.normal(0, 0.3, (n, 5))
    X[:, 0] += lat
    Y[3:, 0] += lat[:-3]                          # Y lags X by 3 bins
    g = (np.arange(n) // 800)                     # contiguous trial blocks
    r = sw.window_subspace(X, Y, g, k=4, max_lag=8, n_shuffles=10)
    assert abs(r.optimal_lag - 3) <= 1
    assert np.isfinite(r.ifi)


def test_split_half_floor_small_for_stable_subspace():
    # A strong, stable shared subspace -> the two halves agree -> small angle.
    rng = np.random.default_rng(7)
    n = 9000
    lat = rng.normal(0, 1, (n, 3))            # 3 shared latents == the 3 compared dims
    X = rng.normal(0, 0.3, (n, 8)); X[:, :3] += 1.8 * lat
    Y = rng.normal(0, 0.3, (n, 8)); Y[:, :3] += 1.8 * lat
    g = np.arange(n) // 300
    r = sw.window_subspace(X, Y, g, k=6, max_lag=3, n_shuffles=10)
    assert np.isfinite(r.split_half_x)
    assert r.split_half_x < 45     # halves agree -> well below orthogonal


def test_split_half_floor_large_for_noise():
    rng = np.random.default_rng(8)
    n = 8000
    X = rng.normal(0, 1, (n, 8)); Y = rng.normal(0, 1, (n, 8))
    g = np.arange(n) // 400
    r = sw.window_subspace(X, Y, g, k=6, max_lag=3, n_shuffles=10)
    # no real shared structure -> dominant direction unstable -> large angle
    assert r.split_half_x > 45


def test_gini_and_members_present():
    rng = np.random.default_rng(3)
    n = 4000
    lat = rng.normal(0, 1, n)
    X = rng.normal(0, 1, (n, 10)); X[:, 0] += 2 * lat   # one unit dominates
    Y = rng.normal(0, 1, (n, 10)); Y[:, 0] += 2 * lat
    g = _grouped(n, 10, rng)
    r = sw.window_subspace(X, Y, g, k=5, max_lag=3, n_shuffles=10)
    assert 0.0 <= r.gini_x <= 1.0
    assert r.member_x.sum() >= 1                  # at least one member
