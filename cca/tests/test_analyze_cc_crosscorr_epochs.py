"""Tests for scripts/analyze_cc_crosscorr_epochs.py — cross-correlograms (r vs lag)
per learning epoch, averaged over an animal's significant CCs, whole and split by
FF/FB label (2026-08-07 meeting ask 2).

Pinned here: per-animal-first aggregation, the significance gate, that the FF and FB
groups never mix, that the lag axis and its sign are passed through untouched, and
that the expert−naive contrast recovers a planted shift."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_cc_crosscorr_epochs as A  # noqa: E402

LAGS = np.arange(-50, 51, 10)


def _curve(animal, pair, dim, epoch, label, r_fn, sig=1):
    return [dict(animal=animal, learner=1, pair=pair, dim=dim, epoch=epoch, bin_ms=10,
                 lag_bins=int(l // 10), lag_ms=int(l), r=float(r_fn(l)), label=label,
                 label_loo=label, sig=sig)
            for l in LAGS]


def _table():
    rows = []
    for epoch in ("naive", "expert"):
        # animal A: two FF CCs at constant 0.4 / 0.2, one FB CC at 0.1
        rows += _curve("A", "CA1-RSC", 1, epoch, "FF", lambda l: 0.4)
        rows += _curve("A", "CA1-RSC", 2, epoch, "FF", lambda l: 0.2)
        rows += _curve("A", "CA1-RSC", 3, epoch, "FB", lambda l: 0.1)
        # a non-significant CC with a wild value: must never enter
        rows += _curve("A", "CA1-RSC", 4, epoch, "FF", lambda l: 9.0, sig=0)
        # animal B: one FF CC at 0.0
        rows += _curve("B", "CA1-RSC", 1, epoch, "FF", lambda l: 0.0)
    return pd.DataFrame(rows)


def test_build_curves_is_per_animal_first_and_gated():
    cur = A.build_curves(_table())
    a_all = cur[(cur["animal"] == "A") & (cur["group"] == "all") &
                (cur["epoch"] == "naive")].sort_values("lag_ms")
    # (0.4 + 0.2 + 0.1) / 3, the sig==0 CC excluded
    assert np.allclose(a_all["r_mean"], (0.4 + 0.2 + 0.1) / 3)
    assert (a_all["n_ccs"] == 3).all()
    a_ff = cur[(cur["animal"] == "A") & (cur["group"] == "FF") &
               (cur["epoch"] == "naive")]
    a_fb = cur[(cur["animal"] == "A") & (cur["group"] == "FB") &
               (cur["epoch"] == "naive")]
    assert np.allclose(a_ff["r_mean"], 0.3) and (a_ff["n_ccs"] == 2).all()
    assert np.allclose(a_fb["r_mean"], 0.1) and (a_fb["n_ccs"] == 1).all()
    # animal B has no FB CC -> no FB rows for B
    assert cur[(cur["animal"] == "B") & (cur["group"] == "FB")].empty


def test_build_curves_keeps_the_lag_axis_intact():
    cur = A.build_curves(_table())
    got = np.sort(cur["lag_ms"].unique())
    assert np.array_equal(got, LAGS)


def test_build_curves_uses_the_requested_label_column():
    tab = _table()
    # flip the LOO label of animal A's dim 1 to FB; with label_col="label_loo" the FF
    # group in A must then be dim 2 alone (0.2)
    tab.loc[(tab["animal"] == "A") & (tab["dim"] == 1), "label_loo"] = "FB"
    cur = A.build_curves(tab, label_col="label_loo")
    a_ff = cur[(cur["animal"] == "A") & (cur["group"] == "FF") &
               (cur["epoch"] == "naive")]
    assert np.allclose(a_ff["r_mean"], 0.2)


def test_cross_animal_mean_and_sem_are_animals_as_n():
    cur = A.build_curves(_table())
    lags, mean, sem, n = A.cross_animal(cur, pair="CA1-RSC", epoch="naive", group="all")
    # A: 0.2333, B: 0.0 -> mean 0.11667, n = 2, not 4 CCs
    assert n == 2
    assert np.allclose(mean, (0.7 / 3 + 0.0) / 2)
    assert np.allclose(sem, np.std([0.7 / 3, 0.0], ddof=1) / np.sqrt(2))
    assert np.array_equal(lags, LAGS)


def test_epoch_contrast_recovers_a_planted_shift():
    """Per (animal, pair, dim, epoch) reductions -> per-animal mean over sig CCs ->
    paired t expert−naive across animals. Plant +0.1 in every animal."""
    rows = []
    jitter = [0.0, 0.01, -0.01, 0.005, -0.005, 0.002]      # avoid a zero-variance t
    for i, an in enumerate("abcdef"):
        for dim in (1, 2):
            base = 0.05 * i
            rows.append(dict(animal=an, pair="CA1-RSC", dim=dim, epoch="naive",
                             ifi=base, peak_r=0.2, sig=1))
            rows.append(dict(animal=an, pair="CA1-RSC", dim=dim, epoch="expert",
                             ifi=base + 0.1 + jitter[i], peak_r=0.2, sig=1))
    red = pd.DataFrame(rows)
    out = A.epoch_contrast(red, metric="ifi", pairs=["CA1-RSC"])
    r = out.iloc[0]
    assert r["n"] == 6
    assert r["delta"] == pytest.approx(0.1 + np.mean(jitter))
    assert r["p"] < 1e-4
    # peak_r is identical in both epochs -> zero variance -> p is NaN, not 0
    out2 = A.epoch_contrast(red, metric="peak_r", pairs=["CA1-RSC"])
    assert np.isnan(out2.iloc[0]["p"])


def test_curve_metrics_separate_an_offset_from_a_peak():
    """A curve that is a flat 0.05 offset above another has the same peak-minus-
    baseline; a curve with a taller peak on the same baseline does not."""
    lags = np.arange(-250, 251, 10)
    base = lambda l: 0.10 * np.exp(-(l / 40.0) ** 2)                 # peak 0.10, far 0
    rows = []
    for l in lags:
        rows.append(dict(animal="A", pair="CA1-RSC", group="all", epoch="naive",
                         lag_ms=int(l), r_mean=base(l), n_ccs=3))
        rows.append(dict(animal="A", pair="CA1-RSC", group="all", epoch="expert",
                         lag_ms=int(l), r_mean=base(l) + 0.05, n_ccs=3))       # offset
        rows.append(dict(animal="A", pair="CA1-RSC", group="all", epoch="intermediate",
                         lag_ms=int(l), r_mean=1.5 * base(l), n_ccs=3))         # taller
    cur = pd.DataFrame(rows)
    m = A.curve_metrics(cur, far_ms=200).set_index("epoch")
    assert m.loc["naive", "peak_r"] == pytest.approx(0.10)
    assert m.loc["expert", "peak_r"] == pytest.approx(0.15)
    assert m.loc["expert", "far_r"] == pytest.approx(0.05, abs=1e-6)
    assert m.loc["expert", "peak_minus_far"] == pytest.approx(
        m.loc["naive", "peak_minus_far"], abs=1e-6)                       # offset cancels
    assert m.loc["intermediate", "peak_minus_far"] == pytest.approx(0.15, abs=1e-6)


def test_epoch_contrast_needs_both_epochs_per_animal():
    red = pd.DataFrame([
        dict(animal="a", pair="CA1-RSC", dim=1, epoch="naive", ifi=0.0, peak_r=0.1, sig=1),
        dict(animal="a", pair="CA1-RSC", dim=1, epoch="expert", ifi=0.2, peak_r=0.1, sig=1),
        dict(animal="b", pair="CA1-RSC", dim=1, epoch="naive", ifi=0.0, peak_r=0.1, sig=1),
        # animal b has no expert row -> dropped
        dict(animal="c", pair="CA1-RSC", dim=1, epoch="naive", ifi=0.1, peak_r=0.1, sig=1),
        dict(animal="c", pair="CA1-RSC", dim=1, epoch="expert", ifi=0.2, peak_r=0.1, sig=1),
        dict(animal="d", pair="CA1-RSC", dim=1, epoch="naive", ifi=0.1, peak_r=0.1, sig=1),
        dict(animal="d", pair="CA1-RSC", dim=1, epoch="expert", ifi=0.5, peak_r=0.1, sig=1),
    ])
    out = A.epoch_contrast(red, metric="ifi", pairs=["CA1-RSC"])
    assert out.iloc[0]["n"] == 3
