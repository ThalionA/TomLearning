"""Tests for tom_cca.cc_aggregate — averaging across an animal's canonical
components BEFORE averaging across animals (animals-as-n).

The whole point of the module is that an animal with five CCs contributes ONE
value to the cross-animal statistic, not five. Every test pins that or one of the
gates (significance, degeneracy, weights, NaNs)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tom_cca import cc_aggregate as agg


def _rows(animal, pair, vals, *, window=5, sig=1, degenerate=0, weight=None):
    """One row per canonical dim for one (animal, pair, window)."""
    out = []
    for i, v in enumerate(vals, start=1):
        r = dict(animal=animal, pair=pair, dim=i, window_bins=window, ifi=v,
                 sig=sig, degenerate=degenerate)
        if weight is not None:
            r["cc_peak"] = weight[i - 1]
        out.append(r)
    return out


# ---------------------------------------------------------------- per_animal_mean

def test_an_animal_with_many_ccs_counts_once():
    """Five CCs in animal A, one in animal B: the per-animal table has one row each,
    and A's row is the mean of its five, not five rows."""
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [0.5, 0.5, 0.5, 0.5, 0.5])
                       + _rows("B", "CA1-RSC", [-0.1]))
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"])
    assert len(out) == 2
    a = out[out["animal"] == "A"].iloc[0]
    b = out[out["animal"] == "B"].iloc[0]
    assert a["mean"] == pytest.approx(0.5) and a["n_ccs"] == 5
    assert b["mean"] == pytest.approx(-0.1) and b["n_ccs"] == 1


def test_non_significant_and_degenerate_rows_are_excluded():
    tab = pd.DataFrame(
        _rows("A", "CA1-RSC", [0.2, 0.4])                      # kept
        + _rows("A", "CA1-RSC", [9.0], sig=0)                  # not significant
        + _rows("A", "CA1-RSC", [-9.0], degenerate=1))         # degenerate
    # the three groups above share dims 1..; make dims unique so the row identity is
    # unambiguous (the function must not care about dim numbering anyway)
    tab["dim"] = range(1, len(tab) + 1)
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"])
    assert len(out) == 1
    assert out.iloc[0]["mean"] == pytest.approx(0.3)
    assert out.iloc[0]["n_ccs"] == 2


def test_gates_can_be_switched_off():
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [0.2, 0.4], sig=0))
    assert agg.per_animal_mean(tab, value="ifi", by=["window_bins"]).empty
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"], sig_only=False)
    assert out.iloc[0]["mean"] == pytest.approx(0.3)


def test_missing_sig_column_raises_rather_than_passing_every_rank():
    """A table gated upstream must SAY so (sig_only=False); silently including all
    ranks would defeat the 'significant CCs only' decision without a trace."""
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [0.2, 0.4])).drop(columns="sig")
    with pytest.raises(KeyError):
        agg.per_animal_mean(tab, value="ifi", by=["window_bins"])
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"], sig_only=False)
    assert out.iloc[0]["mean"] == pytest.approx(0.3)


def test_string_typed_sig_column_is_coerced_not_silently_emptied():
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [0.2, 0.4]))
    tab["sig"] = tab["sig"].astype(str)
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"])
    assert out.iloc[0]["mean"] == pytest.approx(0.3)


def test_weighted_mean_uses_clipped_weights():
    """Weights are held-out CC peaks; a negative one means 'no coupling' and must
    contribute nothing (clip at 0), not pull the mean the wrong way."""
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [1.0, 0.0, -1.0],
                             weight=[0.3, 0.1, -0.5]))
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"], weight="cc_peak")
    row = out.iloc[0]
    assert row["mean"] == pytest.approx(0.0)                    # plain mean
    assert row["wmean"] == pytest.approx((1.0 * 0.3 + 0.0 * 0.1) / 0.4)


def test_weighted_mean_is_nan_when_all_weights_clip_to_zero():
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [0.5, 0.5], weight=[-0.1, 0.0]))
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"], weight="cc_peak")
    assert np.isnan(out.iloc[0]["wmean"])
    assert out.iloc[0]["mean"] == pytest.approx(0.5)


def test_grouping_respects_every_by_column():
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [0.1, 0.3], window=1)
                       + _rows("A", "CA1-RSC", [0.5, 0.7], window=2)
                       + _rows("A", "CA1-DG", [-0.2], window=1))
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"])
    got = out.set_index(["pair", "window_bins"])["mean"]
    assert got[("CA1-RSC", 1)] == pytest.approx(0.2)
    assert got[("CA1-RSC", 2)] == pytest.approx(0.6)
    assert got[("CA1-DG", 1)] == pytest.approx(-0.2)


def test_nan_values_are_dropped_not_propagated():
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [0.2, np.nan, 0.4]))
    out = agg.per_animal_mean(tab, value="ifi", by=["window_bins"])
    assert out.iloc[0]["mean"] == pytest.approx(0.3)
    assert out.iloc[0]["n_ccs"] == 2


def test_all_nan_cell_is_dropped():
    tab = pd.DataFrame(_rows("A", "CA1-RSC", [np.nan, np.nan]))
    assert agg.per_animal_mean(tab, value="ifi", by=["window_bins"]).empty


def test_sign_groups_recombine_to_the_ungrouped_mean():
    """Identity that the figure relies on: the count-weighted combination of the
    +IFI group mean and the −IFI group mean IS the all-CC mean, per animal."""
    rng = np.random.default_rng(0)
    vals = rng.normal(size=9)
    tab = pd.DataFrame(_rows("A", "CA1-RSC", vals.tolist()))
    all_ = agg.per_animal_mean(tab, value="ifi", by=["window_bins"]).iloc[0]
    pos = agg.per_animal_mean(tab[tab["ifi"] > 0], value="ifi", by=["window_bins"]).iloc[0]
    neg = agg.per_animal_mean(tab[tab["ifi"] < 0], value="ifi", by=["window_bins"]).iloc[0]
    recombined = (pos["mean"] * pos["n_ccs"] + neg["mean"] * neg["n_ccs"]) / (
        pos["n_ccs"] + neg["n_ccs"])
    assert recombined == pytest.approx(all_["mean"])


# ------------------------------------------------------------- one_sample_by_pair

def _per_animal(pair, values):
    return pd.DataFrame([dict(animal=f"a{i}", pair=pair, mean=v)
                         for i, v in enumerate(values)])


def test_one_sample_recovers_a_planted_direction_and_a_null():
    pa = pd.concat([_per_animal("CA1-RSC", [0.10, 0.12, 0.09, 0.11, 0.13, 0.10]),
                    _per_animal("CA1-DG", [0.10, -0.10, 0.05, -0.05, 0.02, -0.02])],
                   ignore_index=True)
    out = agg.one_sample_by_pair(pa, value="mean", pairs=["CA1-RSC", "CA1-DG"])
    out = out.set_index("pair")
    assert out.loc["CA1-RSC", "n"] == 6
    assert out.loc["CA1-RSC", "mean"] == pytest.approx(0.108333, abs=1e-6)
    assert out.loc["CA1-RSC", "p"] < 1e-4
    assert out.loc["CA1-RSC", "t"] > 0
    assert out.loc["CA1-DG", "p"] > 0.5
    # BH across the pairs is a sensitivity column: the planted pair survives
    assert bool(out.loc["CA1-RSC", "bh_pass"]) is True
    assert bool(out.loc["CA1-DG", "bh_pass"]) is False


def test_one_sample_skips_pairs_with_too_few_animals():
    pa = pd.concat([_per_animal("CA1-RSC", [0.1, 0.2]),          # n=2 -> skipped
                    _per_animal("CA1-DG", [0.1, 0.2, 0.3])],
                   ignore_index=True)
    out = agg.one_sample_by_pair(pa, value="mean", pairs=["CA1-RSC", "CA1-DG"],
                                 min_n=3)
    assert out["pair"].tolist() == ["CA1-DG"]


def test_one_sample_refuses_more_than_one_row_per_animal():
    """Two windows' worth of per-animal rows passed straight in would silently be
    windows-as-n — exactly the pseudoreplication the helper exists to prevent."""
    pa = pd.concat([_per_animal("CA1-RSC", [0.1, 0.2, 0.3]),
                    _per_animal("CA1-RSC", [0.1, 0.2, 0.3])], ignore_index=True)
    with pytest.raises(ValueError):
        agg.one_sample_by_pair(pa, value="mean", pairs=["CA1-RSC"])


def test_one_sample_keeps_requested_pair_order():
    pa = pd.concat([_per_animal("V1-RSC", [0.1, 0.2, 0.3]),
                    _per_animal("CA1-RSC", [0.1, 0.2, 0.3])], ignore_index=True)
    out = agg.one_sample_by_pair(pa, value="mean", pairs=["CA1-RSC", "V1-RSC"])
    assert out["pair"].tolist() == ["CA1-RSC", "V1-RSC"]


def test_one_sample_wilcoxon_columns_are_additive():
    """wilcoxon=True appends w_stat/w_p; the default columns are unchanged, and the
    default output is byte-identical to the pre-flag behaviour."""
    pa = _per_animal("CA1-RSC", [0.20, 0.30, 0.25, 0.10, 0.15, 0.22])
    out0 = agg.one_sample_by_pair(pa, value="mean", pairs=["CA1-RSC"])
    out1 = agg.one_sample_by_pair(pa, value="mean", pairs=["CA1-RSC"], wilcoxon=True)
    assert list(out0.columns) == ["pair", "n", "mean", "sem", "t", "p", "bh_pass"]
    assert list(out1.columns) == ["pair", "n", "mean", "sem", "t", "p",
                                  "w_stat", "w_p", "bh_pass"]
    pd.testing.assert_frame_equal(out1[out0.columns], out0)
    assert np.isfinite(out1["w_p"].iloc[0])


def test_one_sample_wilcoxon_nan_below_three_nonzero():
    pa = _per_animal("CA1-RSC", [0.0, 0.0, 0.5])   # only one non-zero delta
    out = agg.one_sample_by_pair(pa, value="mean", pairs=["CA1-RSC"], wilcoxon=True)
    assert np.isnan(out["w_p"].iloc[0])
