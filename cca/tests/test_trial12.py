"""Tests for the trial-1-vs-2 frozen-projection analysis (scripts/run_trial12.py).

Ground-truth synthetic data throughout: a session whose coupling is CONSTRUCTED to
differ between the first two running trials, projected through a reference subspace
fit on the later trials. The tests pin (1) the equivalence of the two frozen-fit
routes in src (fit_fixed(trials=) vs early_trials.reference_fit), (2) that
frozen_perm_null on the train rows tests the fit it claims to, (3) recovery of a
constructed trial-1-vs-2 difference with the right sign and lag, (4) the matched-
bins arms, the ordinal->trial_id mapping and the min-bins gate, and (5) the
single-trial moments/curve identity the analysis relies on.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from tom_cca import early_trials, fixed_subspace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_trial12  # noqa: E402


def _session(rng, n_trials=14, n_bins=400, nx=6, ny=6, lead_lag=2,
             uncoupled=(2,), trial_ids_start=0, n_bins_by_ordinal=None):
    """Synthetic (X, Y, trial_ids): a shared latent drives X in every trial; Y
    carries it at lag 0 in every trial EXCEPT ordinal 1 (where Y lags X by
    ``lead_lag`` bins => X leads) and the ``uncoupled`` ordinals (independent
    noise). Trial ids start at ``trial_ids_start`` and may be non-consecutive."""
    Wx = rng.standard_normal((1, nx))
    Wy = rng.standard_normal((1, ny))
    ids = trial_ids_start + np.arange(0, 3 * n_trials, 3)     # non-consecutive
    X_parts, Y_parts, tid_parts = [], [], []
    for o, t in enumerate(ids, start=1):
        nb = (n_bins_by_ordinal or {}).get(o, n_bins)
        s = rng.standard_normal((nb, 1))
        x = s @ Wx + 0.3 * rng.standard_normal((nb, nx))
        if o == 1:                        # X leads Y by lead_lag
            sl = np.vstack([np.zeros((lead_lag, 1)), s[:nb - lead_lag]])
            y = sl @ Wy + 0.3 * rng.standard_normal((nb, ny))
        elif o in uncoupled:              # no coupling at all
            y = rng.standard_normal((nb, ny))
        else:                             # zero-lag coupling (the reference regime)
            y = s @ Wy + 0.3 * rng.standard_normal((nb, ny))
        X_parts.append(x); Y_parts.append(y)
        tid_parts.append(np.full(nb, t))
    return (np.vstack(X_parts), np.vstack(Y_parts), np.concatenate(tid_parts))


def _fit(X, Y, tids, k=4, n_held=10):
    uniq = np.unique(tids)
    train_mask = np.isin(tids, uniq[n_held:])
    return early_trials.reference_fit(X, Y, train_mask, k=k), uniq, train_mask


# --------------------------------------------------------------- src equivalences

def test_fit_fixed_trials_equals_reference_fit_variates():
    """The two frozen-fit routes in src are the same map: fit_fixed(trials=train)
    projected per dim == reference_fit variates, to 1e-10."""
    rng = np.random.default_rng(0)
    X, Y, tids = _session(rng, n_trials=8, n_bins=250)
    uniq = np.unique(tids)
    train_ids = uniq[4:]
    fitA = early_trials.reference_fit(X, Y, np.isin(tids, train_ids), k=4)
    fitB = fixed_subspace.fit_fixed(X, Y, tids, Z=None, k=4, trials=train_ids)
    uA, vA = early_trials.variates(fitA, slice(None))
    for dim in range(min(uA.shape[1], fitB.wx.shape[1])):
        uB, vB = fixed_subspace.project(X, Y, fitB, dim=dim)
        np.testing.assert_allclose(uA[:, dim], uB, atol=1e-10)
        np.testing.assert_allclose(vA[:, dim], vB, atol=1e-10)


def test_frozen_perm_null_on_train_rows_tests_the_fit_itself():
    """r_obs of frozen_perm_null(Sx[train], Sy[train]) must equal model.r — the
    observed statistic and the null are the same estimator on the same fit."""
    rng = np.random.default_rng(1)
    X, Y, tids = _session(rng, n_trials=8, n_bins=250)
    fit, uniq, train_mask = _fit(X, Y, tids, k=4, n_held=4)
    res = fixed_subspace.frozen_perm_null(fit.Sx[train_mask], fit.Sy[train_mask],
                                          n_shuffles=20, seed=0)
    d = fit.model.r.size
    np.testing.assert_allclose(res.r_obs[:d], fit.model.r[:d], atol=1e-10)


def test_single_trial_moments_mask_equals_variate_lag_curve():
    """curve_from_moments with a one-trial mask == variate_lag_curve on that
    trial's bins — the identity that lets any analysis swap between the routes."""
    rng = np.random.default_rng(2)
    n = 900
    u = rng.standard_normal(n)
    v = np.roll(u, 2) + 0.5 * rng.standard_normal(n)
    groups = np.repeat([10, 20, 30], n // 3)
    lags = np.arange(-5, 6)
    trials, M = fixed_subspace.trial_lag_moments(u, v, groups, lags)
    for i, t in enumerate(trials):
        mask = np.zeros(trials.size, dtype=bool); mask[i] = True
        bins = groups == t
        np.testing.assert_allclose(
            fixed_subspace.curve_from_moments(M, mask),
            fixed_subspace.variate_lag_curve(u[bins], v[bins], groups[bins], lags),
            atol=1e-12, equal_nan=True)


# --------------------------------------------------------- ground-truth recovery

def test_constructed_trial12_difference_recovered_with_sign_and_lag():
    """Ordinal 1 is coupled with X leading by 2 bins, ordinal 2 uncoupled: the
    per-trial curves must peak at the true lag in trial 1 and stay flat in
    trial 2, held out of the fit."""
    rng = np.random.default_rng(3)
    X, Y, tids = _session(rng, lead_lag=2, uncoupled=(2,))
    fit, uniq, train_mask = _fit(X, Y, tids, k=4)
    u, v = early_trials.variates(fit, slice(None))
    lags = np.arange(-6, 7)
    crows, trows = run_trial12.pair_rows(
        u, v, tids, uniq, lags, d=1, sig_mask=np.array([True]),
        pvals=np.array([0.005]), r_frozen=fit.model.r[:1],
        velocity=None, stream_index=None, min_bins=50)
    c1 = {r["lag_bins"]: r["r"] for r in crows
          if r["ordinal"] == 1 and r["matched"] == 0}
    c2 = {r["lag_bins"]: r["r"] for r in crows
          if r["ordinal"] == 2 and r["matched"] == 0}
    assert max(c1, key=lambda k: c1[k]) == 2          # peak at the constructed lag
    assert c1[2] > 0.6                                # coupling recovered held-out
    assert max(c2.values()) < 0.3                     # none where none was built
    # positive lags dominate trial 1's curve => IFI would read X-leads
    pos = np.mean([v_ for k, v_ in c1.items() if k > 0])
    neg = np.mean([v_ for k, v_ in c1.items() if k < 0])
    assert pos > neg


def test_matched_arms_ordinal_mapping_and_gate():
    rng = np.random.default_rng(4)
    # ordinal 1 twice the bins of ordinal 2; ordinal 5 below any sensible gate
    X, Y, tids = _session(rng, n_bins_by_ordinal={1: 800, 2: 400, 5: 30},
                          trial_ids_start=7)
    fit, uniq, train_mask = _fit(X, Y, tids, k=4)
    u, v = early_trials.variates(fit, slice(None))
    lags = np.arange(-4, 5)
    crows, trows = run_trial12.pair_rows(
        u, v, tids, uniq, lags, d=2, sig_mask=np.array([True, False]),
        pvals=np.array([0.005, 0.5]), r_frozen=fit.model.r[:2],
        velocity=None, stream_index=None, min_bins=50)
    by = {(r["ordinal"], r["matched"], r["dim"]): r for r in trows}
    # ordinal -> trial_id maps through the NON-consecutive running-trial ids
    assert by[(1, 0, 1)]["trial_id"] == int(uniq[0]) == 7
    assert by[(3, 0, 1)]["trial_id"] == int(uniq[2])
    # raw arm keeps the full counts; pairwise arm cuts BOTH ordinals to min(n1, n2)
    assert by[(1, 0, 1)]["n_bins"] == 800 and by[(2, 0, 1)]["n_bins"] == 400
    assert by[(1, 1, 1)]["n_bins"] == by[(2, 1, 1)]["n_bins"] == 400
    assert (3, 1, 1) not in by                        # matched=1 exists only for 1 and 2
    # common arm: min over gated ordinals (30-bin trial 5 is EXCLUDED from the min)
    assert by[(1, 2, 1)]["n_bins"] == 400
    # the too-short trial is gated: blank r0, present in trials, absent from curves
    assert by[(5, 0, 1)]["r0"] == ""
    assert not any(r["ordinal"] == 5 for r in crows)
    # sig/p flow through per dim
    assert by[(1, 0, 1)]["sig"] == 1 and by[(1, 0, 2)]["sig"] == 0
