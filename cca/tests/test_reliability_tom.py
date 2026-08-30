"""TDD for dataio.load_reliability_tom: Tom's precomputed moving-window
reliability with its shuffle null.

Real layout (measured on TF073_export.mat, 2026-08-30):
  analysis_spatial/reliability/units/reliability_moving_window              (n_trials, n_units)
  analysis_spatial/reliability/units/reliability_moving_window_shuffled_mean (same)
  analysis_spatial/reliability/units/reliability_moving_window_shuffled_sd   (same)
(trials, units) — same orientation as spatial_reliability outputs, transposed
relative to tuning_score.
"""
import h5py
import numpy as np
import pytest

from tom_cca import dataio


def _write(tmp_path, rel, mu, sd):
    p = tmp_path / "TF001_export.mat"
    with h5py.File(p, "w") as f:
        g = (f.create_group("analysis_spatial").create_group("reliability")
             .create_group("units"))
        g.create_dataset("reliability_moving_window", data=rel)
        g.create_dataset("reliability_moving_window_shuffled_mean", data=mu)
        g.create_dataset("reliability_moving_window_shuffled_sd", data=sd)
    return p


def test_z_matches_hand_computation(tmp_path):
    rng = np.random.default_rng(0)
    rel = rng.normal(size=(4, 3))
    mu = rng.normal(size=(4, 3))
    sd = rng.uniform(0.5, 2.0, size=(4, 3))
    out = dataio.load_reliability_tom(_write(tmp_path, rel, mu, sd))
    assert out["rel"].shape == (4, 3)
    np.testing.assert_allclose(out["z"], (rel - mu) / sd, atol=1e-12)


def test_zero_shuffle_sd_gives_nan_not_inf(tmp_path):
    rel = np.ones((2, 2)); mu = np.zeros((2, 2)); sd = np.zeros((2, 2))
    out = dataio.load_reliability_tom(_write(tmp_path, rel, mu, sd))
    assert np.isnan(out["z"]).all()


def test_nan_entries_propagate(tmp_path):
    rel = np.ones((2, 2)); rel[0, 1] = np.nan
    mu = np.zeros((2, 2)); sd = np.ones((2, 2))
    out = dataio.load_reliability_tom(_write(tmp_path, rel, mu, sd))
    assert np.isnan(out["z"][0, 1]) and np.isfinite(out["z"][1, 0])


def test_missing_group_raises_keyerror(tmp_path):
    p = tmp_path / "TF002_export.mat"
    with h5py.File(p, "w") as f:
        f.create_group("analysis_spatial")
    with pytest.raises(KeyError, match="reliability"):
        dataio.load_reliability_tom(p)
