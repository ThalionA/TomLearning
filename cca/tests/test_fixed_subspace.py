"""Tests for the FIXED-subspace arm (meeting 2026-07-28, items 5 and 7).

The question these serve: every epoch comparison in the project so far REFITS the
subspace per epoch, so an epoch difference can be a difference in the fit rather than in
the activity. Here the subspace is identified once across all trials, frozen, and then
each epoch is projected through it — so one canonical component can be lagged across
time and compared naive vs expert on equal terms.
"""

from __future__ import annotations

import numpy as np
import pytest

from tom_cca import fixed_subspace


def _shifted(n_trials=10, n_bins=100, lag_true=3, noise=0.2, seed=0):
    """u carries s(t); v carries s(t - lag_true) => u LEADS v by lag_true bins."""
    rng = np.random.default_rng(seed)
    groups = np.repeat(np.arange(n_trials), n_bins)
    u = rng.standard_normal(groups.size)
    v = np.empty_like(u)
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        seg = u[idx]
        if lag_true == 0:
            v[idx] = seg
            continue
        out = np.empty_like(seg)
        out[:lag_true] = rng.standard_normal(lag_true)
        out[lag_true:] = seg[:-lag_true]
        v[idx] = out
    v = v + noise * rng.standard_normal(v.size)
    return u, v, groups


# ---------------------------------------------------------------------------
# variate_lag_curve — one component, lagged across time
# ---------------------------------------------------------------------------
def test_curve_peaks_at_the_true_lag():
    u, v, groups = _shifted(lag_true=3)
    lags = np.arange(-8, 9)
    r = fixed_subspace.variate_lag_curve(u, v, groups, lags)
    assert lags[int(np.nanargmax(r))] == 3


def test_curve_is_symmetric_for_simultaneous_signals():
    u, v, groups = _shifted(lag_true=0, noise=0.1)
    lags = np.arange(-5, 6)
    r = fixed_subspace.variate_lag_curve(u, v, groups, lags)
    assert lags[int(np.nanargmax(r))] == 0


def test_curve_never_pairs_across_a_trial_boundary():
    """Bin 0 of trial 2 must never be paired with the last bins of trial 1."""
    groups = np.repeat([0, 1], 6)
    u = np.array([0, 1, 2, 3, 4, 5, 100, 101, 102, 103, 104, 105], dtype=float)
    v = u.copy()
    r = fixed_subspace.variate_lag_curve(u, v, groups, [0])
    assert r[0] == pytest.approx(1.0)          # identical within trials
    # a lag that consumes a whole trial leaves nothing pairable
    assert np.isnan(fixed_subspace.variate_lag_curve(u, v, groups, [20])[0])


def test_curve_is_nan_for_a_constant_variate():
    groups = np.repeat([0, 1], 20)
    u = np.ones(40)
    v = np.random.default_rng(0).standard_normal(40)
    assert np.all(np.isnan(fixed_subspace.variate_lag_curve(u, v, groups, [0, 1])))


def test_uncorrelated_variates_stay_near_zero():
    rng = np.random.default_rng(4)
    groups = np.repeat(np.arange(10), 100)
    u = rng.standard_normal(groups.size)
    v = rng.standard_normal(groups.size)
    r = fixed_subspace.variate_lag_curve(u, v, groups, np.arange(-5, 6))
    assert np.nanmax(np.abs(r)) < 0.15


# ---------------------------------------------------------------------------
# balanced_trials — the fixed subspace must not be dominated by one epoch
# ---------------------------------------------------------------------------
def test_balanced_selection_takes_equal_trials_per_epoch():
    epoch_of_trial = {t: ("naive" if t < 10 else
                          "intermediate" if t < 25 else "expert")
                      for t in range(40)}
    sel = fixed_subspace.balanced_trials(epoch_of_trial,
                                         ["naive", "intermediate", "expert"], seed=0)
    counts = {e: sum(1 for t in sel if epoch_of_trial[t] == e)
              for e in ("naive", "intermediate", "expert")}
    assert len(set(counts.values())) == 1       # equal by construction
    assert counts["naive"] == 10                # capped by the smallest epoch


def test_balanced_selection_is_deterministic_for_a_seed():
    eot = {t: ("naive" if t < 12 else "expert") for t in range(30)}
    a = fixed_subspace.balanced_trials(eot, ["naive", "expert"], seed=3)
    b = fixed_subspace.balanced_trials(eot, ["naive", "expert"], seed=3)
    assert a == b


def test_balanced_selection_returns_empty_when_an_epoch_is_missing():
    eot = {t: "naive" for t in range(10)}
    assert fixed_subspace.balanced_trials(eot, ["naive", "expert"], seed=0) == []


# ---------------------------------------------------------------------------
# fit_fixed — identify once, project many times
# ---------------------------------------------------------------------------
def _paired_pops(n_trials=12, n_bins=100, n_units=6, noise=0.3, seed=0):
    rng = np.random.default_rng(seed)
    groups = np.repeat(np.arange(n_trials), n_bins)
    s = rng.standard_normal(groups.size)
    X = np.outer(s, rng.standard_normal(n_units)) + \
        noise * rng.standard_normal((groups.size, n_units))
    Y = np.outer(s, rng.standard_normal(n_units)) + \
        noise * rng.standard_normal((groups.size, n_units))
    return X, Y, groups


def test_fixed_fit_recovers_a_shared_latent():
    X, Y, groups = _paired_pops()
    fit = fixed_subspace.fit_fixed(X, Y, groups, k=4)
    u, v = fixed_subspace.project(X, Y, fit, dim=0)
    r = fixed_subspace.variate_lag_curve(u, v, groups, [0])[0]
    assert r > 0.7


def test_projection_is_the_same_map_for_every_subset():
    """The point of a FIXED subspace: two epochs must go through identical weights."""
    X, Y, groups = _paired_pops()
    fit = fixed_subspace.fit_fixed(X, Y, groups, k=4)
    half = groups < 6
    u_all, _ = fixed_subspace.project(X, Y, fit, dim=0)
    u_half, _ = fixed_subspace.project(X[half], Y[half], fit, dim=0)
    assert u_half == pytest.approx(u_all[half])


def test_fixed_fit_returns_none_on_degenerate_input():
    groups = np.repeat([0, 1], 5)
    X = np.zeros((10, 3)); Y = np.zeros((10, 3))
    assert fixed_subspace.fit_fixed(X, Y, groups, k=2) is None


def test_fixed_fit_restricted_to_a_trial_subset_uses_only_those_trials():
    X, Y, groups = _paired_pops(n_trials=12)
    sub = [0, 1, 2, 3]
    fit_a = fixed_subspace.fit_fixed(X, Y, groups, k=4, trials=sub)
    fit_b = fixed_subspace.fit_fixed(X[np.isin(groups, sub)],
                                     Y[np.isin(groups, sub)],
                                     groups[np.isin(groups, sub)], k=4)
    assert fit_a.wx == pytest.approx(fit_b.wx)


# ---------------------------------------------------------------------------
# curve_half_width — the integration window (meeting item 6)
# ---------------------------------------------------------------------------
def test_half_width_of_a_known_plateau():
    lags = np.arange(-5, 6, dtype=float)
    r = np.where(np.abs(lags) <= 2, 1.0, 0.1)     # >= half-max over lags -2..+2
    assert fixed_subspace.curve_half_width(lags, r) == pytest.approx(4.0)


def test_half_width_ignores_a_disconnected_far_lag_bump():
    """A noise excursion above half-max must not widen the window."""
    lags = np.arange(-6, 7, dtype=float)
    r = np.full(lags.size, 0.1)
    r[np.abs(lags) <= 1] = 1.0
    r[lags == 6] = 0.9                            # disconnected from the peak
    assert fixed_subspace.curve_half_width(lags, r) == pytest.approx(2.0)


def test_half_width_is_nan_for_a_non_positive_peak():
    lags = np.arange(-3, 4, dtype=float)
    assert np.isnan(fixed_subspace.curve_half_width(lags, np.full(lags.size, -0.2)))


def test_half_width_is_censored_at_the_swept_range():
    lags = np.arange(-4, 5, dtype=float)
    assert fixed_subspace.curve_half_width(lags, np.ones(lags.size)) == \
        pytest.approx(8.0)


def test_half_width_handles_unsorted_lags():
    lags = np.array([2.0, -2.0, 0.0, 1.0, -1.0])
    r = np.array([0.1, 0.1, 1.0, 1.0, 1.0])
    assert fixed_subspace.curve_half_width(lags, r) == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# side_peak — is the lag curve ringing rather than decaying?
# ---------------------------------------------------------------------------
def test_side_peak_finds_an_oscillation_period():
    lags = np.arange(-250, 251, 10, dtype=float)
    r = np.cos(2 * np.pi * lags / 130.0)          # ~7.7 Hz, theta
    off, ratio = fixed_subspace.side_peak(lags, r)
    assert off == pytest.approx(130.0, abs=10.0)
    assert ratio > 0.9


def test_side_peak_is_nan_for_a_decaying_curve():
    lags = np.arange(-100, 101, 10, dtype=float)
    r = np.exp(-(lags ** 2) / (2 * 30.0 ** 2))    # clean Gaussian, no ringing
    off, ratio = fixed_subspace.side_peak(lags, r)
    assert np.isnan(off) and np.isnan(ratio)


def test_side_peak_ignores_the_central_peak_itself():
    """A broad flat central lobe must not be reported as its own side peak."""
    lags = np.arange(-50, 51, 10, dtype=float)
    r = np.where(np.abs(lags) <= 20, 1.0, 0.05)
    off, _ = fixed_subspace.side_peak(lags, r)
    assert np.isnan(off) or off > 20


def test_side_peak_is_nan_for_a_non_positive_curve():
    lags = np.arange(-30, 31, 10, dtype=float)
    assert np.isnan(fixed_subspace.side_peak(lags, np.full(lags.size, -0.3))[0])


# ---------------------------------------------------------------------------
# side_peak_null — the ringing detector needs a calibrated null
# ---------------------------------------------------------------------------
def test_side_peak_fires_on_almost_all_noise_so_the_fraction_is_uninformative():
    """The reason the RATE of ringing must never be quoted as evidence."""
    lags = np.arange(-250, 251, 10, dtype=float)
    null = fixed_subspace.side_peak_null(lags, n_draws=300, seed=1)
    assert null.size > 0.9 * 300


def test_noise_side_peaks_are_broadly_spread_not_band_limited():
    lags = np.arange(-250, 251, 10, dtype=float)
    null = fixed_subspace.side_peak_null(lags, n_draws=800, seed=2)
    assert fixed_subspace.band_occupancy(null) < 0.6


def test_band_occupancy_of_a_true_oscillation_is_high():
    offsets = np.full(50, 140.0)              # 7.1 Hz, squarely in theta
    assert fixed_subspace.band_occupancy(offsets) == pytest.approx(1.0)


def test_band_occupancy_excludes_out_of_band_offsets():
    assert fixed_subspace.band_occupancy(np.array([500.0, 500.0])) == pytest.approx(0.0)


def test_band_occupancy_is_nan_without_usable_offsets():
    assert np.isnan(fixed_subspace.band_occupancy(np.array([np.nan, -5.0])))


# ---------------------------------------------------------------------------
# FixedFit exposes the PCA bases the CCA was fitted through
# ---------------------------------------------------------------------------
def test_fit_exposes_pca_bases_matching_the_weights():
    X, Y, groups = _paired_pops()
    fit = fixed_subspace.fit_fixed(X, Y, groups, k=4)
    assert fit.cx is not None and fit.cy is not None
    assert fit.cx.shape[0] == X.shape[1] and fit.cy.shape[0] == Y.shape[1]
    # scores built from the exposed basis must reproduce the projection exactly
    u_direct, _ = fixed_subspace.project(X, Y, fit, dim=0)
    Sx = (X - fit.x_mean) @ fit.cx
    coef = np.linalg.lstsq(fit.cx, fit.wx[:, 0], rcond=None)[0]
    assert (Sx @ coef) == pytest.approx(u_direct, abs=1e-8)


def test_pca_bases_are_orthonormal():
    X, Y, groups = _paired_pops()
    fit = fixed_subspace.fit_fixed(X, Y, groups, k=4)
    assert fit.cx.T @ fit.cx == pytest.approx(np.eye(fit.cx.shape[1]), abs=1e-8)


# ---------------------------------------------------------------------------
# frozen_perm_null — significance of the FROZEN fit's own canonical correlations
# ---------------------------------------------------------------------------
def test_perm_null_matches_cca_fit_on_the_observed_data():
    """The observed statistic must be the frozen fit's own canonical correlations.

    Canonical correlations are the singular values of Qx^T Qy for orthonormal bases of
    the centred scores, so this must agree with core.cca_fit to numerical precision.
    """
    from tom_cca import core
    X, Y, groups = _paired_pops(n_units=5)
    r_ref = core.cca_fit(X - X.mean(0), Y - Y.mean(0)).r
    obs = fixed_subspace.observed_canonical_r(X, Y)
    d = min(len(r_ref), len(obs))
    assert obs[:d] == pytest.approx(np.asarray(r_ref)[:d], abs=1e-8)


def test_perm_null_flags_a_real_shared_latent():
    X, Y, groups = _paired_pops(noise=0.2, n_units=5)
    res = fixed_subspace.frozen_perm_null(X, Y, n_shuffles=100, seed=0, correct=None)
    assert res.mask[0]
    assert res.p[0] < 0.05


def test_perm_null_is_sparse_on_pure_noise():
    rng = np.random.default_rng(11)
    X = rng.standard_normal((1200, 5))
    Y = rng.standard_normal((1200, 5))
    res = fixed_subspace.frozen_perm_null(X, Y, n_shuffles=100, seed=0, correct=None)
    assert res.mask.sum() <= 1              # alpha = 0.05 over 5 dims


def test_perm_null_compares_rank_k_to_rank_k():
    """No cross-fit alignment exists to get wrong: the null for dim k is the shuffled
    distribution of dim k, so the thresholds fall with rank like the observed spectrum.

    This is the structural fix for the bug where per-fold refit statistics were attached
    by bare index to a frozen fit's dimensions.
    """
    X, Y, groups = _paired_pops(n_units=6, noise=0.4)
    res = fixed_subspace.frozen_perm_null(X, Y, n_shuffles=60, seed=0)
    thr = res.threshold[np.isfinite(res.threshold)]
    assert thr[0] >= thr[-1]


def test_perm_null_p_has_the_permutation_floor():
    X, Y, groups = _paired_pops(noise=0.05, n_units=4)
    res = fixed_subspace.frozen_perm_null(X, Y, n_shuffles=50, seed=0, correct=None)
    assert res.p.min() >= 1 / 51 - 1e-12


def test_perm_null_is_deterministic_for_a_seed():
    X, Y, groups = _paired_pops()
    a = fixed_subspace.frozen_perm_null(X, Y, n_shuffles=40, seed=3)
    b = fixed_subspace.frozen_perm_null(X, Y, n_shuffles=40, seed=3)
    assert a.p == pytest.approx(b.p)


# ---------------------------------------------------------------------------
# trial_lag_moments / curve_from_moments — same answer, computed once
# ---------------------------------------------------------------------------
def test_moments_reproduce_variate_lag_curve_exactly():
    u, v, groups = _shifted(lag_true=3, noise=0.3)
    lags = np.arange(-8, 9)
    direct = fixed_subspace.variate_lag_curve(u, v, groups, lags)
    _, M = fixed_subspace.trial_lag_moments(u, v, groups, lags)
    pooled = fixed_subspace.curve_from_moments(M)
    ok = np.isfinite(direct) & np.isfinite(pooled)
    assert ok.sum() > 10
    assert pooled[ok] == pytest.approx(direct[ok], abs=1e-9)


def test_moments_subset_matches_direct_curve_on_that_subset():
    """An epoch curve from a trial mask must equal the direct computation on those
    trials — this is what licenses deriving every epoch and leave-one-out curve from
    one pass."""
    u, v, groups = _shifted(n_trials=12, lag_true=2, noise=0.3)
    lags = np.arange(-6, 7)
    trials, M = fixed_subspace.trial_lag_moments(u, v, groups, lags)
    sel = np.isin(trials, [0, 1, 2, 3])
    pooled = fixed_subspace.curve_from_moments(M, sel)
    m = np.isin(groups, [0, 1, 2, 3])
    direct = fixed_subspace.variate_lag_curve(u[m], v[m], groups[m], lags)
    ok = np.isfinite(direct) & np.isfinite(pooled)
    assert ok.sum() > 5
    assert pooled[ok] == pytest.approx(direct[ok], abs=1e-9)


def test_moments_complement_is_the_leave_epoch_out_curve():
    u, v, groups = _shifted(n_trials=12, noise=0.3)
    lags = np.arange(-5, 6)
    trials, M = fixed_subspace.trial_lag_moments(u, v, groups, lags)
    sel = np.isin(trials, [0, 1, 2])
    pooled = fixed_subspace.curve_from_moments(M, ~sel)
    m = ~np.isin(groups, [0, 1, 2])
    direct = fixed_subspace.variate_lag_curve(u[m], v[m], groups[m], lags)
    ok = np.isfinite(direct) & np.isfinite(pooled)
    assert pooled[ok] == pytest.approx(direct[ok], abs=1e-9)


def test_moments_are_additive_over_trials():
    u, v, groups = _shifted(n_trials=8)
    lags = np.arange(-4, 5)
    _, M = fixed_subspace.trial_lag_moments(u, v, groups, lags)
    assert M.sum(axis=0) == pytest.approx(M[:4].sum(0) + M[4:].sum(0))
