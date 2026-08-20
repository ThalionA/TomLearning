"""Tests for scripts/analyze_trial12.py — the trial-1-vs-2 analysis wiring.

The statistics live in tom_cca (cc_aggregate, perdim_ifi, paired_stats) and are
tested there. These pin the wiring end-to-end on synthetic curve/trial tables with
a PLANTED trial-1-vs-2 difference: which rows enter, that dims collapse per animal
before any across-animal test, that the delta is trial1 - trial2 with the right
sign, that the adjacent-trial control measures scale rather than sign, and that a
gated (blank) trial drops the animal's delta instead of poisoning it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_trial12 as A  # noqa: E402

LAGS = list(range(-8, 9))


def _curve(ifi_strength, peak=0.4):
    """A lag curve whose positive side exceeds its negative side by
    ``ifi_strength`` (in units of the peak) — IFI sign follows its sign."""
    out = {}
    for lag in LAGS:
        base = peak * np.exp(-abs(lag) / 4.0)
        out[lag] = base * (1.0 + ifi_strength) if lag > 0 else base
    return out


def _rows(animal, pair, ordinal, ifi_strength, *, matched=0, dims=(1, 2),
          sig=(1, 1), r0=0.4, n_bins=600, blank=False):
    """Curve rows + trial rows for one (animal, pair, ordinal, matched) cell."""
    crows, trows = [], []
    for i, dim in enumerate(dims):
        c = _curve(ifi_strength if dim == 1 else ifi_strength * 0.5)
        meta = dict(animal=animal, learner=1, pair=pair, dim=dim, ordinal=ordinal,
                    trial_id=100 + ordinal, matched=matched, sig=sig[i],
                    p_perdim=0.01, r_frozen=0.5, n_bins=n_bins, n_fit_trials=20)
        trows.append({**meta, "r0": "" if blank else r0, "n_gaps": 2,
                      "vel_mean": 12.0, "vel_median": 11.0})
        if blank:
            continue
        for lag in LAGS:
            crows.append({**meta, "bin_ms": 10, "lag_bins": lag,
                          "lag_ms": lag * 10, "r": round(c[lag], 5)})
    return crows, trows


def _tables(spec, pair="CA1-RSC", matched_arms=(0, 1, 2)):
    """``spec`` = {animal: {ordinal: (ifi_strength, r0)}} -> (curves, trials)."""
    crows, trows = [], []
    for animal, ords in spec.items():
        for ordinal, val in ords.items():
            strength, r0 = val if isinstance(val, tuple) else (val, 0.4)
            for m in matched_arms:
                if m == 1 and ordinal not in (1, 2):
                    continue
                c, t = _rows(animal, pair, ordinal, strength, matched=m, r0=r0)
                crows += c; trows += t
    return pd.DataFrame(crows), pd.DataFrame(trows)


# ------------------------------------------------------------------ windows

def test_windows_table_is_per_ordinal_and_arm():
    curves, _ = _tables({1: {1: 0.5, 2: 0.0}})
    win = A.per_trial_windows(curves)
    assert {"ordinal", "matched", "window_bins", "ifi", "degenerate"} <= set(win.columns)
    # one row per (animal, pair, ordinal, matched, dim, window)
    n_arms = curves["matched"].nunique()
    assert win.groupby(["ordinal", "matched", "dim"]).ngroups == (
        curves.groupby(["ordinal", "matched", "dim"]).ngroups)
    assert win["window_bins"].max() == max(LAGS)
    assert n_arms == 3


def test_planted_ifi_sign_is_recovered_per_ordinal():
    """Trial 1 built X-leading, trial 2 balanced: IFI(1) > 0 ≈ IFI(2)."""
    curves, _ = _tables({1: {1: 0.8, 2: 0.0}})
    win = A.per_trial_windows(curves)
    hw = win[(win["window_bins"] == A.HEADLINE_W) & (win["dim"] == 1)
             & (win["matched"] == 1)]
    ifi = {int(o): float(v) for o, v in zip(hw["ordinal"], hw["ifi"])}
    assert ifi[1] > 0.15
    assert abs(ifi[2]) < 1e-9


def test_wide_band_is_disjoint_from_the_headline_window():
    curves, _ = _tables({1: {1: 0.5, 2: 0.0}})
    wide = A.wide_band(curves)
    assert not wide.empty
    # the wide band must ignore |lag| < WIDE_MIN — plant a curve that is flat
    # inside the headline window and asymmetric only outside it
    rows = []
    for lag in LAGS:
        r = 0.3 if abs(lag) < A.WIDE_MIN else (0.5 if lag > 0 else 0.1)
        rows.append(dict(animal=1, learner=1, pair="CA1-RSC", dim=1, ordinal=1,
                         matched=0, lag_bins=lag, r=r, sig=1))
    w = A.wide_band(pd.DataFrame(rows))
    assert float(w["ifi"].iloc[0]) > 0.5          # driven purely by the outer lags


# -------------------------------------------------------------------- deltas

def test_delta_is_trial1_minus_trial2_and_needs_both():
    pa = pd.DataFrame([
        dict(animal="a", pair="CA1-RSC", matched=1, window_bins=5, ifi=0.30),
        dict(animal="a", pair="CA1-RSC", matched=1, window_bins=5, ordinal=2, ifi=0.10),
        dict(animal="b", pair="CA1-RSC", matched=1, window_bins=5, ordinal=1, ifi=0.50),
    ])
    pa.loc[0, "ordinal"] = 1
    d = A.delta12(pa, "ifi")
    assert len(d) == 1                                   # animal b has no trial 2
    assert float(d["delta"].iloc[0]) == pytest_approx(0.20)      # 0.30 - 0.10
    assert float(d["delta_abs"].iloc[0]) == pytest_approx(0.20)


def test_a_gated_blank_trial_drops_that_animals_delta():
    """A trial below the bin gate writes a blank r0; that animal must not appear
    in the delta (and must be visible in the completeness table)."""
    crows, trows = [], []
    for animal in ("a", "b"):
        for ordinal in (1, 2):
            blank = (animal == "b" and ordinal == 1)
            c, t = _rows(animal, "CA1-RSC", ordinal, 0.4, matched=0, blank=blank)
            crows += c; trows += t
    trials = pd.DataFrame(trows)
    d = A.delta12(trials[trials["dim"] == 1], "r0")
    assert set(d["animal"]) == {"a"}
    comp = A.completeness(trials)
    row = comp[comp["matched"] == 0].iloc[0]
    assert row["finite_t1"] == 1 and row["finite_t2"] == 2 and row["n_deltas"] == 1


# ------------------------------------------------------------------- control

def test_control_z_measures_scale_not_sign():
    """z = d12 / (sqrt(2)*SD over ordinals 3..10): a 1-vs-2 delta the size of the
    adjacent-trial spread gives |z| ~ 1/sqrt(2)*1, and a huge one gives a big z —
    while a control band with the SAME signed mean does not cancel it."""
    rng = np.random.default_rng(0)
    ctrl = list(0.10 + 0.05 * rng.standard_normal(8))     # sd ~ 0.05
    rows = []
    for animal, d12 in (("a", 0.05), ("b", 0.50)):
        vals = {1: 0.30 + d12, 2: 0.30}
        vals.update({o: v for o, v in zip(A.CTRL_ORDS, ctrl)})
        for o, v in vals.items():
            rows.append(dict(animal=animal, pair="CA1-RSC", matched=2, ordinal=o,
                             ifi=v))
    cs = A.control_stats(pd.DataFrame(rows), "ifi")
    z = {r["animal"]: r["z"] for _, r in cs.iterrows()}
    sd = float(np.std(ctrl, ddof=1))
    assert z["a"] == pytest_approx(0.05 / (np.sqrt(2) * sd))
    assert abs(z["b"]) > abs(z["a"]) * 5                  # scale, not sign
    assert (cs["n_ctrl"] == 8).all()


def test_control_requires_enough_adjacent_ordinals():
    rows = [dict(animal="a", pair="CA1-RSC", matched=2, ordinal=o, ifi=0.2)
            for o in (1, 2, 3, 4)]                         # only 2 control ordinals
    assert A.control_stats(pd.DataFrame(rows), "ifi").empty


# ---------------------------------------------------------- animals-as-n unit

def test_allsig_collapses_dims_before_the_across_animal_test():
    """Two dims per animal must become ONE value per animal before any test —
    otherwise the n is dims, the most pseudoreplicated unit available."""
    curves, _ = _tables({a: {1: 0.6, 2: 0.0} for a in range(1, 7)})
    win = A.per_trial_windows(curves)
    from tom_cca import cc_aggregate
    allsig = cc_aggregate.per_animal_mean(
        win, value="ifi", by=["ordinal", "matched", "window_bins"], weight="cc_peak")
    cell = allsig[(allsig["ordinal"] == 1) & (allsig["matched"] == 1)
                  & (allsig["window_bins"] == A.HEADLINE_W)]
    assert len(cell) == 6                       # six animals, not twelve dim-rows
    assert (cell["n_ccs"] == 2).all()           # each collapsed two CCs
    d = A.delta12(allsig[allsig["window_bins"] == A.HEADLINE_W], "mean")
    assert len(d[d["matched"] == 1]) == 6
    assert (d[d["matched"] == 1]["delta"] > 0).all()      # planted direction


# -------------------------------------------- 2026-08-20 audit regression tests

def test_wide_band_cc_peak_is_the_all_lag_max():
    """`cc_peak` must mean the same thing in wide_band as in windows_table — both
    feed per_animal_mean(weight='cc_peak'). A band-local max is a different quantity
    under the same name and silently changed both wide-band results."""
    rows = []
    for lag in LAGS:
        # the true peak sits INSIDE the headline window, i.e. outside the wide band
        r = 0.9 if lag == 0 else (0.2 if lag > 0 else 0.05)
        rows.append(dict(animal=1, learner=1, pair="CA1-RSC", dim=1, ordinal=1,
                         matched=0, lag_bins=lag, r=r, sig=1))
    df = pd.DataFrame(rows)
    w = A.wide_band(df)
    assert float(w["cc_peak"].iloc[0]) == pytest_approx(0.9)   # not the band max 0.2
    wt = A.per_trial_windows(df)
    assert float(wt["cc_peak"].iloc[0]) == pytest_approx(0.9)  # same as windows_table


def test_dim_matched_keeps_only_dims_present_in_both_ordinals():
    """A delta must average the SAME dims in both trials. dim_matched drops a dim
    that is degenerate (or absent) in either ordinal, so the two means share support."""
    rows = []
    for ordinal in (1, 2):
        for dim in (1, 2, 3):
            # dim 3 is degenerate in trial 2 only -> must be dropped from BOTH
            degen = int(dim == 3 and ordinal == 2)
            rows.append(dict(animal="a", pair="CA1-RSC", matched=1, window_bins=5,
                             ordinal=ordinal, dim=dim, ifi=0.1 * dim,
                             degenerate=degen, sig=1, cc_peak=0.3))
    out = A.dim_matched(pd.DataFrame(rows), "ifi")
    assert sorted(out["dim"].unique()) == [1, 2]
    assert len(out) == 4                       # 2 dims x 2 ordinals
    # and a dim failing the sig gate goes too
    rows2 = [r | {"sig": 0} if r["dim"] == 2 else r for r in rows]
    out2 = A.dim_matched(pd.DataFrame(rows2), "ifi")
    assert sorted(out2["dim"].unique()) == [1]


def test_sweep_test_recovers_a_planted_window_effect_and_is_null_on_noise():
    """The nested-window test: a consistent delta across windows must be detected,
    and pure noise must not be."""
    rng = np.random.default_rng(7)
    wins = list(range(1, 11))
    rows = []
    for animal in range(12):
        for w in wins:
            # pair A: every animal shifts +0.30 at every window (real, consistent)
            rows.append(dict(animal=f"a{animal}", pair="CA1-RSC", matched=1,
                             window_bins=w, ordinal=1, dim=1, degenerate=0, sig=1,
                             ifi=0.30 + 0.05 * rng.standard_normal()))
            rows.append(dict(animal=f"a{animal}", pair="CA1-RSC", matched=1,
                             window_bins=w, ordinal=2, dim=1, degenerate=0, sig=1,
                             ifi=0.05 * rng.standard_normal()))
            # pair B: pure noise in both ordinals
            for o in (1, 2):
                rows.append(dict(animal=f"a{animal}", pair="V1-RSC", matched=1,
                                 window_bins=w, ordinal=o, dim=1, degenerate=0,
                                 sig=1, ifi=0.20 * rng.standard_normal()))
    win = pd.DataFrame(rows)
    out = A.sweep_test(win, pd.DataFrame(), n_perm=500, seed=1)
    got = {r["pair"]: r for _, r in out.iterrows()}
    assert got["CA1-RSC"]["p_perm"] < 0.01 and bool(got["CA1-RSC"]["bh_pass"])
    assert got["V1-RSC"]["p_perm"] > 0.05 and not bool(got["V1-RSC"]["bh_pass"])
    assert got["CA1-RSC"]["mean_at_best"] > 0.2


def pytest_approx(v, rel=1e-6):
    import pytest
    return pytest.approx(v, rel=rel)
