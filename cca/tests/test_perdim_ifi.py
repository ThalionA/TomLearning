"""Tests for the per-canonical-dimension Information Flow Index (meeting 2026-07-28,
item 1: "do different CCs have different IFI?").

Ground truth throughout: the sign convention is POSITIVE lag => X leads Y (see
``lagged.lag_slice``), so a curve whose mass sits at positive lags must give IFI = +1.
"""

from __future__ import annotations

import numpy as np
import pytest

from tom_cca import perdim_ifi


LAGS = np.arange(-5, 6)


# ---------------------------------------------------------------------------
# curve_ifi — the single-curve reduction
# ---------------------------------------------------------------------------
def test_symmetric_curve_gives_zero_ifi():
    cc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.9, 0.5, 0.4, 0.3, 0.2, 0.1])
    assert perdim_ifi.curve_ifi(LAGS, cc) == pytest.approx(0.0)


def test_all_mass_at_positive_lag_gives_plus_one():
    cc = np.where(LAGS > 0, 0.5, 0.0)
    assert perdim_ifi.curve_ifi(LAGS, cc) == pytest.approx(1.0)


def test_all_mass_at_negative_lag_gives_minus_one():
    cc = np.where(LAGS < 0, 0.5, 0.0)
    assert perdim_ifi.curve_ifi(LAGS, cc) == pytest.approx(-1.0)


def test_lag_zero_is_excluded_from_the_index():
    """Lag 0 is neither X-leads nor Y-leads; changing it must not move the IFI."""
    cc = np.where(LAGS > 0, 0.5, 0.0)
    bumped = cc.copy()
    bumped[LAGS == 0] = 99.0
    assert perdim_ifi.curve_ifi(LAGS, cc) == perdim_ifi.curve_ifi(LAGS, bumped)


def test_negative_cc_is_clipped_not_signed():
    """Held-out CC can go negative; it means 'no coupling', not 'reverse coupling'."""
    cc = np.where(LAGS > 0, 0.5, -0.5)
    assert perdim_ifi.curve_ifi(LAGS, cc) == pytest.approx(1.0)


def test_nan_lags_are_ignored():
    cc = np.where(LAGS > 0, 0.5, 0.0).astype(float)
    cc[0] = np.nan                      # a lag with too few segment-aware pairs
    assert perdim_ifi.curve_ifi(LAGS, cc) == pytest.approx(1.0)


def test_all_nan_curve_is_nan():
    assert np.isnan(perdim_ifi.curve_ifi(LAGS, np.full(LAGS.size, np.nan)))


def test_unsorted_input_gives_same_answer():
    """CSV row order is not guaranteed — the reduction must sort by lag itself."""
    cc = np.where(LAGS > 0, 0.5, 0.1)
    order = np.random.default_rng(0).permutation(LAGS.size)
    assert perdim_ifi.curve_ifi(LAGS[order], cc[order]) == pytest.approx(
        perdim_ifi.curve_ifi(LAGS, cc))


# ---------------------------------------------------------------------------
# curve_peak_lag
# ---------------------------------------------------------------------------
def test_peak_lag_finds_the_argmax():
    cc = np.zeros(LAGS.size)
    cc[LAGS == 3] = 1.0
    assert perdim_ifi.curve_peak_lag(LAGS, cc) == 3


def test_peak_lag_is_nan_when_all_nan():
    assert np.isnan(perdim_ifi.curve_peak_lag(LAGS, np.full(LAGS.size, np.nan)))


# ---------------------------------------------------------------------------
# ifi_by_dim — the long-format group-by
# ---------------------------------------------------------------------------
def _two_dim_frame():
    """dim 1 leads at POSITIVE lag, dim 2 at NEGATIVE lag — opposite directions."""
    lag = np.concatenate([LAGS, LAGS])
    dim = np.concatenate([np.ones(LAGS.size), 2 * np.ones(LAGS.size)])
    cc = np.concatenate([np.where(LAGS > 0, 0.5, 0.0),
                         np.where(LAGS < 0, 0.5, 0.0)])
    return lag, cc, dim


def test_ifi_by_dim_recovers_opposite_directions():
    lag, cc, dim = _two_dim_frame()
    dims, ifi, peaks = perdim_ifi.ifi_by_dim(lag, cc, dim)
    assert np.array_equal(dims, [1, 2])
    assert ifi[0] == pytest.approx(1.0)
    assert ifi[1] == pytest.approx(-1.0)
    assert peaks[0] > 0 and peaks[1] < 0


def test_ifi_by_dim_is_row_order_invariant():
    lag, cc, dim = _two_dim_frame()
    order = np.random.default_rng(1).permutation(lag.size)
    dims_a, ifi_a, _ = perdim_ifi.ifi_by_dim(lag, cc, dim)
    dims_b, ifi_b, _ = perdim_ifi.ifi_by_dim(lag[order], cc[order], dim[order])
    assert np.array_equal(dims_a, dims_b)
    assert ifi_a == pytest.approx(ifi_b)


def test_ifi_by_dim_returns_dims_in_ascending_order():
    lag, cc, dim = _two_dim_frame()
    dims, _, _ = perdim_ifi.ifi_by_dim(lag[::-1], cc[::-1], dim[::-1])
    assert np.array_equal(dims, np.sort(dims))


# ---------------------------------------------------------------------------
# sign_agreement — is the direction consistent across canonical rank?
# ---------------------------------------------------------------------------
def test_sign_agreement_all_same_sign_is_one():
    assert perdim_ifi.sign_agreement(np.array([0.3, 0.2, 0.5])) == pytest.approx(1.0)


def test_sign_agreement_all_opposite_is_zero():
    assert perdim_ifi.sign_agreement(np.array([0.3, -0.2, -0.5])) == pytest.approx(0.0)


def test_sign_agreement_half_is_half():
    assert perdim_ifi.sign_agreement(np.array([0.3, 0.2, -0.5])) == pytest.approx(0.5)


def test_sign_agreement_needs_a_reference_and_a_partner():
    assert np.isnan(perdim_ifi.sign_agreement(np.array([0.3])))
    assert np.isnan(perdim_ifi.sign_agreement(np.array([])))


def test_sign_agreement_ignores_exact_zero_ifi():
    """IFI == 0 is the 'no asymmetry' return value, not evidence either way."""
    assert perdim_ifi.sign_agreement(np.array([0.3, 0.0, 0.4])) == pytest.approx(1.0)


def test_sign_agreement_nan_reference_is_nan():
    assert np.isnan(perdim_ifi.sign_agreement(np.array([np.nan, 0.2, 0.5])))


# ---------------------------------------------------------------------------
# rank_split — CC1 vs the rest of the (significant) spectrum
# ---------------------------------------------------------------------------
def test_rank_split_averages_the_tail():
    ifi = np.array([0.6, 0.2, 0.4])
    sig = np.array([1, 1, 1])
    cc1, rest = perdim_ifi.rank_split(ifi, sig)
    assert cc1 == pytest.approx(0.6)
    assert rest == pytest.approx(0.3)


def test_rank_split_uses_significant_dims_only_for_the_tail():
    ifi = np.array([0.6, 0.2, 0.9])
    sig = np.array([1, 1, 0])            # dim 3 is not significant
    _, rest = perdim_ifi.rank_split(ifi, sig)
    assert rest == pytest.approx(0.2)


def test_rank_split_tail_is_nan_when_no_significant_tail():
    ifi = np.array([0.6, 0.2])
    sig = np.array([1, 0])
    cc1, rest = perdim_ifi.rank_split(ifi, sig)
    assert cc1 == pytest.approx(0.6)
    assert np.isnan(rest)


def test_rank_split_requires_cc1_regardless_of_its_own_flag():
    """CC1 is the reference by construction — it is reported even if flagged n.s.,
    so the contrast is never silently redefined onto a different dim."""
    ifi = np.array([0.6, 0.2])
    sig = np.array([0, 1])
    cc1, rest = perdim_ifi.rank_split(ifi, sig)
    assert cc1 == pytest.approx(0.6)
    assert rest == pytest.approx(0.2)


# ---------------------------------------------------------------------------
# ifi_windows_by_dim — IFI of each CC as a function of integration window
# (meeting item 1: "calculate IFI at different lags for each CC ... do CCs
#  have different signs?")
# ---------------------------------------------------------------------------
LAGB = np.arange(-5, 6)


def _two_dim_opposite():
    """dim 1 leads at POSITIVE lag, dim 2 at NEGATIVE lag — opposite signs."""
    lag = np.concatenate([LAGB, LAGB])
    dim = np.concatenate([np.ones(LAGB.size), 2 * np.ones(LAGB.size)])
    cc = np.concatenate([np.where(LAGB > 0, 0.5, 0.0),
                         np.where(LAGB < 0, 0.5, 0.0)])
    return lag, cc, dim


def test_windows_recover_opposite_signs_per_dim():
    lag, cc, dim = _two_dim_opposite()
    dims, wins, ifi, _ = perdim_ifi.ifi_windows_by_dim(lag, cc, dim)
    assert np.array_equal(dims, [1, 2])
    assert np.array_equal(wins, np.arange(1, 6))
    assert ifi[0, -1] == pytest.approx(1.0)      # dim 1 positive at widest window
    assert ifi[1, -1] == pytest.approx(-1.0)     # dim 2 negative


def test_windows_are_cumulative_in_lag():
    """Entry w-1 uses only |lag| <= w, so widening can change the value."""
    lag = np.arange(-3, 4)
    cc = np.zeros(lag.size)
    cc[lag == 3] = 1.0                            # mass only at the widest lag
    dim = np.ones(lag.size)
    _, wins, ifi, _ = perdim_ifi.ifi_windows_by_dim(lag, cc, dim)
    assert wins.tolist() == [1, 2, 3]
    assert ifi[0, 0] == pytest.approx(0.0)        # window 1 sees no mass
    assert ifi[0, 2] == pytest.approx(1.0)        # window 3 sees it


def test_degenerate_flag_marks_zero_from_no_signal_not_balance():
    """information_flow_index returns 0.0 when BOTH sides clip to zero. That is
    'no coupling', not 'balanced' — the flag keeps the two distinguishable."""
    lag = np.arange(-3, 4)
    cc = np.full(lag.size, -0.2)                  # all negative -> clips to 0
    dim = np.ones(lag.size)
    _, _, ifi, degen = perdim_ifi.ifi_windows_by_dim(lag, cc, dim)
    assert np.all(ifi[0] == 0.0)
    assert np.all(degen[0])


def test_genuinely_balanced_curve_is_zero_but_not_degenerate():
    lag = np.arange(-3, 4)
    cc = np.full(lag.size, 0.4)                   # symmetric, real signal
    dim = np.ones(lag.size)
    _, _, ifi, degen = perdim_ifi.ifi_windows_by_dim(lag, cc, dim)
    assert ifi[0, -1] == pytest.approx(0.0)
    assert not np.any(degen[0])


def test_windows_row_order_invariant():
    lag, cc, dim = _two_dim_opposite()
    order = np.random.default_rng(0).permutation(lag.size)
    a = perdim_ifi.ifi_windows_by_dim(lag, cc, dim)[2]
    b = perdim_ifi.ifi_windows_by_dim(lag[order], cc[order], dim[order])[2]
    assert a == pytest.approx(b, nan_ok=True)
