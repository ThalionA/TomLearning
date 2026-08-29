"""TDD for spatial_reliability.trial_map_reliability.

Reliability of unit u at trial t = the mean Pearson correlation between the
unit's spatially-binned map on trial t and its maps on the neighbouring trials
within ±half_window (session edges clipped: fewer neighbours, never padded).
"""
import numpy as np
import pytest

from tom_cca import spatial_reliability


def _fr(maps_by_trial):
    """Stack per-trial (n_bins,) maps for ONE unit -> (n_trials, n_bins, 1)."""
    return np.stack(maps_by_trial)[:, :, None]


def test_identical_maps_give_reliability_one():
    rng = np.random.default_rng(0)
    m = rng.normal(size=20)
    fr = _fr([m] * 5)
    rel = spatial_reliability.trial_map_reliability(fr, half_window=2)
    assert rel.shape == (5, 1)
    assert np.allclose(rel, 1.0)


def test_matches_bruteforce_corrcoef_on_random_data():
    rng = np.random.default_rng(1)
    n_trials, n_bins, n_units = 7, 30, 3
    fr = rng.normal(size=(n_trials, n_bins, n_units))
    rel = spatial_reliability.trial_map_reliability(fr, half_window=2)
    for t in range(n_trials):
        for u in range(n_units):
            neigh = [s for s in range(max(0, t - 2), min(n_trials, t + 3)) if s != t]
            expect = np.mean([np.corrcoef(fr[t, :, u], fr[s, :, u])[0, 1]
                              for s in neigh])
            assert rel[t, u] == pytest.approx(expect, abs=1e-12)


def test_edges_clip_to_existing_trials_only():
    # Trial 0 must average over trials {1, 2} only, NOT wrap or pad.
    rng = np.random.default_rng(2)
    fr = rng.normal(size=(4, 25, 1))
    rel = spatial_reliability.trial_map_reliability(fr, half_window=2)
    c = lambda a, b: np.corrcoef(fr[a, :, 0], fr[b, :, 0])[0, 1]
    assert rel[0, 0] == pytest.approx(np.mean([c(0, 1), c(0, 2)]), abs=1e-12)
    assert rel[3, 0] == pytest.approx(np.mean([c(3, 1), c(3, 2)]), abs=1e-12)


def test_half_window_one_uses_only_adjacent_trials():
    rng = np.random.default_rng(3)
    fr = rng.normal(size=(5, 25, 1))
    rel = spatial_reliability.trial_map_reliability(fr, half_window=1)
    c = lambda a, b: np.corrcoef(fr[a, :, 0], fr[b, :, 0])[0, 1]
    assert rel[2, 0] == pytest.approx(np.mean([c(2, 1), c(2, 3)]), abs=1e-12)
    assert rel[0, 0] == pytest.approx(c(0, 1), abs=1e-12)


def test_anticorrelated_neighbours_average_to_zero():
    rng = np.random.default_rng(4)
    m = rng.normal(size=20)
    # even trials: +m, odd trials: -m -> interior windows hold two +1 and two -1.
    fr = _fr([m if t % 2 == 0 else -m for t in range(6)])
    rel = spatial_reliability.trial_map_reliability(fr, half_window=2)
    assert rel[2, 0] == pytest.approx(0.0, abs=1e-12)
    assert rel[3, 0] == pytest.approx(0.0, abs=1e-12)


def test_flat_map_trial_is_nan_and_excluded_from_neighbours():
    rng = np.random.default_rng(5)
    m = rng.normal(size=20)
    fr = _fr([m, np.zeros(20), m, m, m])          # trial 1 silent (zero variance)
    rel = spatial_reliability.trial_map_reliability(fr, half_window=2)
    assert np.isnan(rel[1, 0])                    # its own reliability undefined
    # trial 0's neighbours are {1, 2}; the flat trial-1 pair is dropped, not 0.
    assert rel[0, 0] == pytest.approx(1.0, abs=1e-12)


def test_nan_bins_use_pairwise_complete_bins():
    rng = np.random.default_rng(6)
    a = rng.normal(size=20)
    b = rng.normal(size=20)
    b_nan = b.copy()
    b_nan[:5] = np.nan
    fr = _fr([a, b_nan])
    rel = spatial_reliability.trial_map_reliability(fr, half_window=2)
    expect = np.corrcoef(a[5:], b[5:])[0, 1]
    assert rel[0, 0] == pytest.approx(expect, abs=1e-12)
    assert rel[1, 0] == pytest.approx(expect, abs=1e-12)


def test_too_few_finite_bins_gives_nan():
    a = np.full(20, np.nan)
    a[:3] = [1.0, 2.0, 3.0]
    b = np.full(20, np.nan)
    b[:3] = [1.0, 2.5, 2.9]
    fr = _fr([a, b])
    rel = spatial_reliability.trial_map_reliability(fr, half_window=2, min_bins=10)
    assert np.isnan(rel).all()


def test_single_trial_returns_nan():
    fr = np.random.default_rng(7).normal(size=(1, 20, 2))
    rel = spatial_reliability.trial_map_reliability(fr, half_window=2)
    assert rel.shape == (1, 2)
    assert np.isnan(rel).all()


def test_epoch_mean_maps_one_based_trial_ids_onto_rows():
    # Temporal-stream trial ids are 1-BASED (measured on TF073: ids 1..214 for
    # 214 spatial rows). epoch_mean_reliability must subtract 1, not truncate.
    rel = np.arange(10, dtype=float)[:, None]     # rel[t, 0] == t
    out = spatial_reliability.epoch_mean_reliability(rel, trial_ids=[1, 2, 3])
    assert out[0] == pytest.approx(np.mean([0.0, 1.0, 2.0]))
    out_last = spatial_reliability.epoch_mean_reliability(rel, trial_ids=[10])
    assert out_last[0] == pytest.approx(9.0)


def test_epoch_mean_ignores_nan_trials():
    rel = np.array([[1.0], [np.nan], [3.0]])
    out = spatial_reliability.epoch_mean_reliability(rel, trial_ids=[1, 2, 3])
    assert out[0] == pytest.approx(2.0)
