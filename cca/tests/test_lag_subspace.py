"""Tests for lagged communication subspaces (meeting 2026-07-28, items 2/3/4).

Ground truth: a latent that reaches Y a known number of bins AFTER it reaches X must
give (a) a held-out CC that peaks at exactly that lag, and (b) a lagged subspace that is
closer to the true one than the unlagged fit is.
"""

from __future__ import annotations

import numpy as np
import pytest

from tom_cca import lag_subspace


def _trials(n_trials=12, n_bins=120, rng_seed=0):
    """Trial ids for a flat (n_trials*n_bins, ...) sample stack."""
    return np.repeat(np.arange(n_trials), n_bins)


def _leading_pair(lag_true=4, n_trials=12, n_bins=120, n_units=6, noise=0.35,
                  seed=0):
    """X carries latent s(t); Y carries s(t - lag_true), i.e. X LEADS Y by lag_true.

    Pairing X[t] with Y[t + lag_true] therefore re-aligns them, so the held-out CC
    must peak at lag = +lag_true under the positive-lag-means-X-leads convention.
    """
    rng = np.random.default_rng(seed)
    groups = _trials(n_trials, n_bins)
    s = rng.standard_normal(groups.size)
    s_shift = np.empty_like(s)
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        v = s[idx]
        if lag_true == 0:
            s_shift[idx] = v                  # simultaneous, no lead
            continue
        sh = np.empty_like(v)
        sh[:lag_true] = rng.standard_normal(lag_true)
        sh[lag_true:] = v[:-lag_true]         # Y lags X within the trial
        s_shift[idx] = sh
    wx = rng.standard_normal(n_units)
    wy = rng.standard_normal(n_units)
    X = np.outer(s, wx) + noise * rng.standard_normal((groups.size, n_units))
    Y = np.outer(s_shift, wy) + noise * rng.standard_normal((groups.size, n_units))
    return X, Y, groups


# ---------------------------------------------------------------------------
# segment_lag — the alignment primitive
# ---------------------------------------------------------------------------
def test_zero_lag_is_identity():
    X, Y, groups = _leading_pair()
    Xp, Yp, _, _, gp = lag_subspace.segment_lag(X, Y, None, groups, 0)
    assert np.array_equal(Xp, X) and np.array_equal(Yp, Y)
    assert np.array_equal(gp, groups)


def test_positive_lag_pairs_x_with_later_y():
    groups = np.repeat([0, 1], 10)
    X = np.arange(20, dtype=float)[:, None]
    Y = np.arange(100, 120, dtype=float)[:, None]
    Xp, Yp, _, _, _ = lag_subspace.segment_lag(X, Y, None, groups, 3)
    # within trial 0: x bins 0..6 pair with y bins 3..9
    assert Xp[0, 0] == 0 and Yp[0, 0] == 103


def test_lag_never_crosses_a_trial_boundary():
    groups = np.repeat([0, 1], 10)
    X = np.arange(20, dtype=float)[:, None]
    Y = np.arange(20, dtype=float)[:, None]
    Xp, Yp, _, _, gp = lag_subspace.segment_lag(X, Y, None, groups, 3)
    # 7 usable pairs per trial, and no pair may mix trial 0's x with trial 1's y
    assert Xp.shape[0] == 14
    assert np.all((Xp[:, 0] < 10) == (gp == 0))
    assert np.all(np.abs(Yp[:, 0] - Xp[:, 0]) == 3)


def test_confound_is_aligned_to_each_area_separately():
    """Z must follow X's timepoints for X and Y's timepoints for Y — partialling a
    third area out of Y at the wrong time would inject a spurious lag."""
    groups = np.repeat([0], 10)
    X = np.arange(10, dtype=float)[:, None]
    Y = np.arange(10, dtype=float)[:, None]
    Z = np.arange(10, dtype=float)[:, None]
    Xp, Yp, Zx, Zy, _ = lag_subspace.segment_lag(X, Y, Z, groups, 2)
    assert np.array_equal(Zx, Xp)      # Z at X's times
    assert np.array_equal(Zy, Yp)      # Z at Y's times


def test_lag_too_long_for_every_trial_returns_none():
    groups = np.repeat([0, 1], 5)
    X = np.zeros((10, 2)); Y = np.zeros((10, 2))
    assert lag_subspace.segment_lag(X, Y, None, groups, 50) is None


# ---------------------------------------------------------------------------
# lagged_fit — leak-free CC + neuron-space weights at a fixed lag
# ---------------------------------------------------------------------------
def test_heldout_cc_peaks_at_the_true_lag():
    X, Y, groups = _leading_pair(lag_true=4)
    ccs = {}
    for lag in range(0, 8):
        fit = lag_subspace.lagged_fit(X, Y, groups, lag=lag, k=4, seed=0)
        ccs[lag] = fit.cc[0]
    assert max(ccs, key=ccs.get) == 4


def test_negative_lag_is_weaker_than_the_true_positive_lag():
    """Directionality sanity: X leads, so -4 must not beat +4."""
    X, Y, groups = _leading_pair(lag_true=4)
    pos = lag_subspace.lagged_fit(X, Y, groups, lag=4, k=4, seed=0).cc[0]
    neg = lag_subspace.lagged_fit(X, Y, groups, lag=-4, k=4, seed=0).cc[0]
    assert pos > neg


def test_weights_are_in_neuron_space():
    X, Y, groups = _leading_pair()
    fit = lag_subspace.lagged_fit(X, Y, groups, lag=2, k=4, seed=0)
    assert fit.wx.shape[0] == X.shape[1]
    assert fit.wy.shape[0] == Y.shape[1]
    assert fit.wx.shape[1] == fit.wy.shape[1]


def test_heldout_cc_is_not_inflated_on_unrelated_areas():
    rng = np.random.default_rng(3)
    groups = _trials()
    X = rng.standard_normal((groups.size, 6))
    Y = rng.standard_normal((groups.size, 6))
    fit = lag_subspace.lagged_fit(X, Y, groups, lag=0, k=4, seed=0)
    assert fit.cc[0] < 0.25          # honest CV on noise must stay near zero


def test_fit_returns_none_when_the_lag_is_unusable():
    groups = np.repeat([0, 1], 5)
    X = np.zeros((10, 3)); Y = np.zeros((10, 3))
    assert lag_subspace.lagged_fit(X, Y, groups, lag=50, k=2) is None


# ---------------------------------------------------------------------------
# subspace_angle — the item-3 similarity readout
# ---------------------------------------------------------------------------
def test_identical_subspaces_are_zero_degrees():
    w = np.random.default_rng(0).standard_normal((10, 3))
    assert lag_subspace.subspace_angle(w, w, 3) == pytest.approx(0.0, abs=1e-6)


def test_orthogonal_subspaces_are_ninety_degrees():
    wa = np.zeros((4, 2)); wa[0, 0] = wa[1, 1] = 1.0
    wb = np.zeros((4, 2)); wb[2, 0] = wb[3, 1] = 1.0
    assert lag_subspace.subspace_angle(wa, wb, 2) == pytest.approx(90.0, abs=1e-6)


def test_angle_is_symmetric():
    rng = np.random.default_rng(1)
    wa, wb = rng.standard_normal((10, 3)), rng.standard_normal((10, 3))
    assert lag_subspace.subspace_angle(wa, wb, 3) == pytest.approx(
        lag_subspace.subspace_angle(wb, wa, 3))


def test_angle_ignores_column_scaling():
    """Canonical weights have arbitrary scale; the SUBSPACE must not depend on it."""
    rng = np.random.default_rng(2)
    wa = rng.standard_normal((10, 3))
    assert lag_subspace.subspace_angle(wa, wa * np.array([1.0, -5.0, 0.2]), 3) == \
        pytest.approx(0.0, abs=1e-6)


def test_angle_is_nan_when_a_basis_is_empty():
    wa = np.random.default_rng(0).standard_normal((10, 3))
    assert np.isnan(lag_subspace.subspace_angle(wa, np.zeros((10, 0)), 3))


# ---------------------------------------------------------------------------
# lag_sweep — the efficient many-lag form used by the driver
# ---------------------------------------------------------------------------
def test_sweep_recovers_the_true_lag():
    X, Y, groups = _leading_pair(lag_true=4)
    fits = lag_subspace.lag_sweep(X, Y, groups, lags=range(-6, 7), k=4, seed=0)
    peak = max(fits, key=lambda g: fits[g].cc[0])
    assert peak == 4


def test_sweep_agrees_with_single_lag_fit_on_the_peak():
    """The shared-residualisation shortcut must not move the answer."""
    X, Y, groups = _leading_pair(lag_true=3)
    sweep = lag_subspace.lag_sweep(X, Y, groups, lags=[0, 3], k=4, seed=0)
    single = lag_subspace.lagged_fit(X, Y, groups, lag=3, k=4, seed=0)
    assert sweep[3].cc[0] == pytest.approx(single.cc[0], abs=0.05)
    assert sweep[3].cc[0] > sweep[0].cc[0]


def test_sweep_with_confound_still_recovers_the_lag():
    X, Y, groups = _leading_pair(lag_true=4, seed=5)
    rng = np.random.default_rng(9)
    Z = rng.standard_normal((groups.size, 4))
    X = X + Z @ rng.standard_normal((4, X.shape[1]))     # shared third-area drive
    Y = Y + Z @ rng.standard_normal((4, Y.shape[1]))
    fits = lag_subspace.lag_sweep(X, Y, groups, Z=Z, lags=range(-6, 7), k=4, seed=0)
    assert max(fits, key=lambda g: fits[g].cc[0]) == 4


def test_sweep_skips_unusable_lags_rather_than_failing():
    X, Y, groups = _leading_pair(n_bins=30)
    fits = lag_subspace.lag_sweep(X, Y, groups, lags=[0, 5, 500], k=4, seed=0)
    assert 0 in fits and 500 not in fits


def test_sweep_weights_have_consistent_neuron_dimension_across_lags():
    """Item 3 compares subspaces across lags — the bases must live in one space."""
    X, Y, groups = _leading_pair()
    fits = lag_subspace.lag_sweep(X, Y, groups, lags=[-4, 0, 4], k=4, seed=0)
    assert {f.wx.shape[0] for f in fits.values()} == {X.shape[1]}
    assert {f.wy.shape[0] for f in fits.values()} == {Y.shape[1]}


# ---------------------------------------------------------------------------
# split_half_floor — the noise floor the item-3 verdict is read against
# ---------------------------------------------------------------------------
def test_floor_returns_a_separate_value_per_area():
    """X and Y differ in unit count, so they must not share one floor."""
    X, Y, groups = _leading_pair(n_units=6)
    Y = Y[:, :3]
    fx, fy = lag_subspace.split_half_floor(X, Y, groups, k=3, d_use=2, seed=0)
    assert np.isfinite(fx) and np.isfinite(fy)
    assert fx != fy


def test_floor_is_an_angle_in_degrees():
    X, Y, groups = _leading_pair()
    fx, fy = lag_subspace.split_half_floor(X, Y, groups, k=4, d_use=2, seed=0)
    assert 0.0 <= fx <= 90.0 and 0.0 <= fy <= 90.0


def test_floor_is_lower_for_a_strongly_coupled_pair_than_for_noise():
    """A real shared latent makes the two halves agree; pure noise does not."""
    X, Y, groups = _leading_pair(lag_true=0, noise=0.05)
    rng = np.random.default_rng(7)
    Xn = rng.standard_normal((groups.size, 6))
    Yn = rng.standard_normal((groups.size, 6))
    sig_x, _ = lag_subspace.split_half_floor(X, Y, groups, k=4, d_use=2, seed=0)
    noise_x, _ = lag_subspace.split_half_floor(Xn, Yn, groups, k=4, d_use=2, seed=0)
    assert sig_x < noise_x


def test_floor_is_nan_with_too_few_trials():
    groups = np.repeat([0, 1], 50)
    X = np.zeros((100, 3)); Y = np.zeros((100, 3))
    fx, fy = lag_subspace.split_half_floor(X, Y, groups, k=2)
    assert np.isnan(fx) and np.isnan(fy)
