"""TDD for dataio.load_tuning_scores: per-unit per-trial spatial tuning score
z-scored against the export's own 1000-shuffle null.

Real layout (measured on TF073_export.mat, 2026-08-28):
  analysis_spatial/tuning_score/score     (n_units, n_trials)
  analysis_spatial/tuning_score/shuffled  (n_units, n_trials, n_shuffles)
NB: transposed relative to reliability/units (trials, units).
"""
import h5py
import numpy as np
import pytest

from tom_cca import dataio


def _write(tmp_path, score, shuffled):
    p = tmp_path / "TF001_export.mat"
    with h5py.File(p, "w") as f:
        g = f.create_group("analysis_spatial").create_group("tuning_score")
        g.create_dataset("score", data=score)
        g.create_dataset("shuffled", data=shuffled)
    return p


def test_z_matches_hand_computation(tmp_path):
    rng = np.random.default_rng(0)
    score = rng.normal(size=(3, 4))
    shuffled = rng.normal(size=(3, 4, 50))
    out = dataio.load_tuning_scores(_write(tmp_path, score, shuffled))
    assert out["score"].shape == (3, 4)
    expect = (score - shuffled.mean(axis=2)) / shuffled.std(axis=2, ddof=0)
    np.testing.assert_allclose(out["z"], expect, atol=1e-12)


def test_zero_shuffle_sd_gives_nan_not_inf(tmp_path):
    score = np.ones((2, 3))
    shuffled = np.zeros((2, 3, 10))          # sd = 0 everywhere
    out = dataio.load_tuning_scores(_write(tmp_path, score, shuffled))
    assert np.isnan(out["z"]).all()


def test_nan_scores_propagate(tmp_path):
    rng = np.random.default_rng(1)
    score = rng.normal(size=(2, 3)); score[0, 1] = np.nan
    shuffled = rng.normal(size=(2, 3, 20))
    out = dataio.load_tuning_scores(_write(tmp_path, score, shuffled))
    assert np.isnan(out["z"][0, 1])
    assert np.isfinite(out["z"][1, 2])


def test_missing_group_raises_keyerror(tmp_path):
    p = tmp_path / "TF002_export.mat"
    with h5py.File(p, "w") as f:
        f.create_group("analysis_spatial")
    with pytest.raises(KeyError, match="tuning_score"):
        dataio.load_tuning_scores(p)
