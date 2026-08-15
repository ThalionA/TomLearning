"""Tests for the OVERALL (ungrouped) IFI additions to scripts/analyze_cc_ifi_signs.py
— the 2026-08-07 meeting asks 1 and 3.

The heavy lifting is in tom_cca.cc_aggregate (tested separately); these pin the
wiring: which rows enter, which weight column is used, which window is tested,
and that the CC-strength-weighted mean is carried alongside the plain one."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_cc_ifi_signs as A  # noqa: E402


def _cell(animal, pair, window, ifis, peaks=None, sig=1, degenerate=0):
    peaks = peaks or [0.2] * len(ifis)
    return [dict(animal=animal, learner=1, pair=pair, dim=i + 1,
                 window_bins=window, window_ms=window * 10, ifi=v,
                 degenerate=degenerate, cc_peak=w, sig=sig)
            for i, (v, w) in enumerate(zip(ifis, peaks))]


def _windows_table():
    rows = []
    for w in (1, 5, 10):
        rows += _cell("a1", "CA1-RSC", w, [0.3, 0.1, -0.1], peaks=[0.3, 0.1, 0.05])
        rows += _cell("a2", "CA1-RSC", w, [0.2, 0.2])
        rows += _cell("a3", "CA1-RSC", w, [0.1])
        rows += _cell("a4", "CA1-RSC", w, [9.9], sig=0)          # must be ignored
        rows += _cell("a1", "CA1-DG", w, [0.05, -0.05])
        rows += _cell("a2", "CA1-DG", w, [-0.1, 0.1])
        rows += _cell("a3", "CA1-DG", w, [0.02])
        rows += _cell("a4", "CA1-DG", w, [0.0], degenerate=1)   # must be ignored
    return pd.DataFrame(rows)


def test_overall_by_animal_is_one_row_per_animal_pair_window():
    out = A.overall_by_animal(_windows_table())
    # a4 (non-sig in CA1-RSC, degenerate in CA1-DG) contributes nothing
    assert set(out["animal"]) == {"a1", "a2", "a3"}
    assert len(out) == 3 * 3 + 3 * 3          # 3 animals per pair x 3 windows
    a1 = out[(out["animal"] == "a1") & (out["pair"] == "CA1-RSC") &
             (out["window_bins"] == 5)].iloc[0]
    assert a1["mean_ifi"] == pytest.approx(0.1)
    assert a1["wmean_ifi"] == pytest.approx((0.3 * 0.3 + 0.1 * 0.1 - 0.1 * 0.05) / 0.45)
    assert a1["n_ccs"] == 3
    assert {"animal", "pair", "window_bins", "window_ms", "mean_ifi", "wmean_ifi",
            "n_ccs"} <= set(out.columns)


def test_overall_direction_tests_only_the_headline_window():
    pa = A.overall_by_animal(_windows_table())
    out = A.overall_direction(pa, window_bins=5)
    assert set(out["pair"]) == {"CA1-RSC", "CA1-DG"}
    rsc = out[out["pair"] == "CA1-RSC"].iloc[0]
    # per-animal means: a1 0.1, a2 0.2, a3 0.1 -> mean 0.1333, all positive
    assert rsc["n"] == 3
    assert rsc["mean"] == pytest.approx(0.13333, abs=1e-4)
    assert rsc["t"] > 0
    assert {"n", "mean", "sem", "t", "p", "bh_pass", "wmean", "wt", "wp"} <= set(out.columns)
    assert rsc["window_ms"] == 50


def test_overall_direction_skips_pairs_below_min_n():
    pa = A.overall_by_animal(_windows_table())
    out = A.overall_direction(pa, window_bins=5, min_n=4)
    assert out.empty                                  # 3 animals per pair < 4


def test_overall_recombines_from_sign_groups():
    """The all-CC mean equals the count-weighted mean of the +IFI and −IFI groups
    labelled at the SAME window — the identity the figure prints."""
    tab = _windows_table()
    pa = A.overall_by_animal(tab)
    at = tab[(tab["window_bins"] == 5) & (tab["sig"] == 1) & (tab["degenerate"] == 0)]
    for (an, pr), g in at.groupby(["animal", "pair"]):
        v = g["ifi"].to_numpy(float)
        v = v[np.isfinite(v)]
        expect = v.mean()
        got = pa[(pa["animal"] == an) & (pa["pair"] == pr) &
                 (pa["window_bins"] == 5)]["mean_ifi"].iloc[0]
        assert got == pytest.approx(expect)
