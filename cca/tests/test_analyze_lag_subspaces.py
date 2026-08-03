"""Tests for the ANALYSIS layer of the lagged-subspace arm.

This layer had no test coverage at all until 2026-08-03, and both bugs found by
adversarial verification lived here rather than in the primitives — one in which floor
column was subtracted from which angle, one in how the two areas were aggregated. Unit
tests on `lag_subspace` passed throughout both.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_lag_subspaces as A  # noqa: E402


def _frame(rows):
    base = dict(pair="CA1-DG", bin_ms=10, epoch="session")
    return pd.DataFrame([{**base, **r} for r in rows])


def _cell(animal, lag, ax, ay, fx, fy):
    return dict(animal=animal, lag_ms=lag, angle_x_cc1=ax, angle_y_cc1=ay,
                floor_x_cc1=fx, floor_y_cc1=fy, angle_x=ax, angle_y=ay,
                floor_x=fx, floor_y=fy)


def test_each_area_is_compared_against_its_own_floor():
    """X's angle must meet X's floor. Subtracting one shared floor from both was BUG 1."""
    rows = []
    for a in range(6):
        rows.append(_cell(a, 0, 10.0, 10.0, 10.0, 10.0))
        # X sits 20 deg above ITS floor; Y sits exactly AT its own, higher, floor.
        # Both floors stay under the estimability gate so both areas contribute.
        rows.append(_cell(a, 50, 40.0, 60.0, 20.0, 60.0))
    s = A.stability(_frame(rows), dims=1)
    got = s[s.lag_ms == 50]["angle_minus_floor"].iloc[0]
    assert got == pytest.approx(10.0)   # mean of (+20, 0); a shared X floor gives +30


def test_an_unestimable_area_is_excluded_not_averaged_in():
    """BUG 2 of 2026-08-03: gating on the POOLED floor let an area with no usable
    estimate contribute a large negative delta and drag the pair toward 'at floor'."""
    rows = []
    for a in range(6):
        rows.append(_cell(a, 0, 10.0, 10.0, 10.0, 10.0))
        # X estimable and rotating (+25); Y unestimable (floor 88) with a big negative term
        rows.append(_cell(a, 50, 45.0, 60.0, 20.0, 88.0))
    s = A.stability(_frame(rows), dims=1)
    r = s[s.lag_ms == 50].iloc[0]
    assert r["angle_minus_floor"] == pytest.approx(25.0)  # pooled would give (25-28)/2
    assert r["n_areas_used"] == 6                          # one area per animal, not two


def test_animal_drops_out_when_neither_area_is_estimable():
    rows = []
    for a in range(6):
        rows.append(_cell(a, 0, 10.0, 10.0, 10.0, 10.0))
        rows.append(_cell(a, 50, 80.0, 80.0, 85.0, 88.0))
    s = A.stability(_frame(rows), dims=1)
    r = s[s.lag_ms == 50].iloc[0]
    assert r["n_areas_used"] == 0 and not r["estimable"]


def test_gate_sensitivity_spans_the_threshold_and_restores_the_module_default():
    rows = []
    for a in range(6):
        rows.append(_cell(a, 0, 10.0, 10.0, 10.0, 10.0))
        rows.append(_cell(a, 50, 45.0, 60.0, 20.0, 75.0))
    before = A.NOT_ESTIMABLE_DEG
    sens = A.gate_sensitivity(_frame(rows), dims=1)
    assert len(sens) == len(A.GATE_SWEEP)
    assert A.NOT_ESTIMABLE_DEG == before      # must not leak the mutated global
    assert np.isinf(sens["gate_deg"].iloc[-1])


def test_bonferroni_is_across_distinct_absolute_lags():
    rng = np.random.default_rng(0)
    rows = []
    for a in range(6):
        for lag in (0, 50, -50, 100, -100):
            # jitter, else the within-animal deltas are identical, the paired-t variance
            # is zero and p comes back NaN rather than testing the correction
            rows.append(_cell(a, lag, 40.0 + rng.normal(0, 3), 40.0 + rng.normal(0, 3),
                              20.0, 20.0))
    s = A.stability(_frame(rows), dims=1)
    r = s[s.lag_ms == 50].iloc[0]
    assert r["p_bonf"] == pytest.approx(min(1.0, r["p"] * 2))   # |50| and |100|
