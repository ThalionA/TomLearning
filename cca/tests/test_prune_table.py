"""Synthetic ground-truth tests for the config-pruning aggregation."""

from __future__ import annotations

import numpy as np

from tom_cca import prune_table as pt


def _bucket(cc, ifi=None, ifi_win=None, n_sig_cells=None,
            n_cells_total=None, n_cells_with_sig=None):
    """Minimal bucket; defaults derived from the cc pool for convenience."""
    cc = list(cc)
    if ifi is None:
        ifi = [0.0] * len(cc)
    if ifi_win is None:
        ifi_win = {1: list(ifi)}
    if n_sig_cells is None:
        n_sig_cells = [1] * len(cc)
    if n_cells_total is None:
        n_cells_total = len(cc)
    if n_cells_with_sig is None:
        n_cells_with_sig = sum(1 for v in n_sig_cells if v > 0)
    return {"cc": cc, "ifi": ifi, "ifi_win": ifi_win,
            "n_sig_cells": n_sig_cells, "n_cells_total": n_cells_total,
            "n_cells_with_sig": n_cells_with_sig}


# --- parse_tag ------------------------------------------------------------

def test_parse_tag():
    assert pt.parse_tag("landmark50_res_samp15") == (50, "res", "samp15")
    assert pt.parse_tag("landmark25_sig_var95") == (25, "sig", "var95")
    assert pt.parse_tag("not_a_tag") is None


# --- cc_diag: overfitting diagnostics -------------------------------------

def test_cc_diag_known_values():
    d = pt.cc_diag([0.2, 0.4, 0.6, 0.8])
    assert d["n"] == 4
    assert d["max"] == 0.8
    assert np.isclose(d["median"], 0.5)
    assert d["frac_ge_099"] == 0.0


def test_cc_diag_flags_saturation():
    # Half the dims are essentially perfect held-out CC -> overfit tell.
    d = pt.cc_diag([1.0, 1.0, 0.3, 0.4])
    assert d["max"] == 1.0
    assert np.isclose(d["frac_ge_099"], 0.5)


def test_cc_diag_ignores_nonfinite_and_handles_empty():
    d = pt.cc_diag([np.nan, 0.5, np.inf])
    assert d["n"] == 1 and d["max"] == 0.5
    empty = pt.cc_diag([])
    assert empty["n"] == 0 and np.isnan(empty["median"])


# --- p-value wiring -------------------------------------------------------

def test_headline_pvalues_detects_offset_signal():
    rng = np.random.default_rng(0)
    naive = _bucket(list(0.05 + 0.02 * rng.standard_normal(40)))   # ~0
    expert = _bucket(list(0.6 + 0.05 * rng.standard_normal(40)))   # clearly >0
    p = pt.headline_pvalues(naive, expert)
    assert p["cc_p_expert_vs0"] < 1e-3          # expert CC clearly non-zero
    assert p["cc_p_naive_vs_expert"] < 1e-3     # epochs clearly differ


def test_headline_pvalues_nan_on_missing_epoch():
    p = pt.headline_pvalues(None, _bucket([0.5] * 10))
    assert all(np.isnan(v) for v in p.values())


def test_wilcoxon_and_mwu_underpowered_return_nan():
    assert np.isnan(pt.wilcoxon_vs0([0.1, 0.2]))        # < 6 samples
    assert np.isnan(pt.wilcoxon_vs0([0.0] * 10))        # all zero
    assert np.isnan(pt.mwu([0.1, 0.2], [0.3]))          # group < 3


# --- build_prune_rows -----------------------------------------------------

def _two_config_per_cpe():
    """A clean config and an overfit config, one pair, three epochs."""
    rng = np.random.default_rng(1)
    pair = ("CA1", "V1")
    per = {}
    # clean config: moderate CC, grows naive->expert
    for e, lo in (("naive", 0.3), ("intermediate", 0.4), ("expert", 0.55)):
        per[("landmark50_res_samp15", pair, e)] = _bucket(
            list(lo + 0.05 * rng.standard_normal(30)))
    # overfit config: held-out CC pinned at ~1.0
    for e in pt.EPOCHS:
        per[("landmark50_res_var95", pair, e)] = _bucket(
            list(np.clip(1.0 - 0.001 * rng.standard_normal(30), 0, 1)))
    return per


def test_build_prune_rows_shape_and_flags():
    rows = pt.build_prune_rows(_two_config_per_cpe())
    assert len(rows) == 2
    assert set(pt.prune_fieldnames()) == set(rows[0])  # exact column coverage
    by_tag = {r["tag"]: r for r in rows}
    assert by_tag["landmark50_res_var95"]["overfit_flag"] is True
    assert by_tag["landmark50_res_samp15"]["overfit_flag"] is False
    # clean config: expert median CC exceeds naive
    clean = by_tag["landmark50_res_samp15"]
    assert clean["delta_cc_expert_naive"] > 0
    assert np.isclose(clean["median_cc_naive"],
                      np.median(_two_config_per_cpe()[
                          ("landmark50_res_samp15", ("CA1", "V1"),
                           "naive")]["cc"]))


def test_rollup_counts():
    rows = pt.build_prune_rows(_two_config_per_cpe())
    roll = {r["tag"]: r for r in pt.rollup_rows(rows)}
    assert roll["landmark50_res_var95"]["n_overfit_pairs"] == 1
    assert roll["landmark50_res_var95"]["frac_overfit_pairs"] == 1.0
    assert roll["landmark50_res_samp15"]["n_overfit_pairs"] == 0
    assert np.isclose(roll["landmark50_res_var95"]["max_cc_any"], 1.0, atol=0.05)
    assert set(pt.rollup_fieldnames()) == set(
        next(iter(pt.rollup_rows(rows))))
