"""Equivalence guards for the 2026-08-19 de-duplication (audit/REPORT_2026-08-19.md §2).

Each ``_ref_*`` below is a VERBATIM copy of an implementation that existed in the tree
before the consolidation (file:line in the docstring). The tests assert that the single
shared implementation reproduces every old copy to floating precision on shared fixtures,
so the merge cannot change a number. Keep these references frozen — they are the record
of what the old code did, not code to maintain.
"""

from __future__ import annotations

import numpy as np
import pytest

from tom_cca import core, fixed_subspace, lag_subspace, lagged, lagpairs, paired_stats

RTOL = 1e-12
ATOL = 1e-12


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
def _trials(rng, n_trials=12, lo=4, hi=60):
    """Per-bin trial ids with uneven trial lengths (some shorter than |lag|+2)."""
    lens = rng.integers(lo, hi, size=n_trials)
    return np.repeat(np.arange(n_trials), lens)


@pytest.fixture
def flat_pair():
    rng = np.random.default_rng(3)
    groups = _trials(rng)
    n = groups.size
    X = rng.standard_normal((n, 9))
    Y = rng.standard_normal((n, 7))
    Z = rng.standard_normal((n, 4))
    return X, Y, Z, groups


# ---------------------------------------------------------------------------
# 2.2  BH-FDR: lagged.fdr_bh / fixed_subspace._fdr_bh / scripts -> paired_stats.fdr_bh
# ---------------------------------------------------------------------------
def _ref_fdr_bh_lagged(pvals, q: float = 0.05) -> np.ndarray:
    """Verbatim src/tom_cca/lagged.py:323 (pre-2026-08-19)."""
    p = np.asarray(pvals, dtype=float)
    finite = np.isfinite(p)
    if not np.any(finite):
        return np.zeros_like(p, dtype=bool)
    pf = p[finite]
    n = pf.size
    order = np.argsort(pf)
    ranked = pf[order]
    thresh = q * np.arange(1, n + 1) / n
    passed = ranked <= thresh
    k = np.max(np.flatnonzero(passed)) + 1 if np.any(passed) else 0
    cut = ranked[k - 1] if k else -np.inf
    out = np.zeros_like(p, dtype=bool)
    out[finite] = pf <= cut
    return out


def _ref_fdr_bh_fixed(pvals, q: float = 0.05) -> np.ndarray:
    """Verbatim src/tom_cca/fixed_subspace.py:340 (pre-2026-08-19)."""
    p = np.asarray(pvals, dtype=float)
    finite = np.isfinite(p)
    if not np.any(finite):
        return np.zeros_like(p, dtype=bool)
    pf = p[finite]
    n = pf.size
    ranked = np.sort(pf)
    passed = ranked <= q * np.arange(1, n + 1) / n
    k = np.max(np.flatnonzero(passed)) + 1 if np.any(passed) else 0
    cut = ranked[k - 1] if k else -np.inf
    out = np.zeros_like(p, dtype=bool)
    out[finite] = pf <= cut
    return out


@pytest.mark.parametrize("seed", range(20))
def test_fdr_bh_single_implementation(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(1, 40))
    p = rng.uniform(0, 1, n) ** rng.uniform(0.5, 4)      # mix of tiny and large p
    if n > 3:
        p[rng.integers(0, n, 2)] = np.nan                  # non-finite entries
    for q in (0.05, 0.1):
        new = paired_stats.fdr_bh(p, q)
        assert np.array_equal(new, _ref_fdr_bh_lagged(p, q))
        assert np.array_equal(new, _ref_fdr_bh_fixed(p, q))
        # the old names must still resolve to the same function
        assert np.array_equal(lagged.fdr_bh(p, q), new)
        assert np.array_equal(fixed_subspace._fdr_bh(p, q), new)


def test_fdr_bh_all_nan_and_empty():
    assert not paired_stats.fdr_bh(np.array([np.nan, np.nan])).any()
    assert paired_stats.fdr_bh(np.array([])).size == 0


# ---------------------------------------------------------------------------
# 2.3  PCA -> scores: driver _pca_scores x3, lag_subspace._scores, subspace_window,
#      early_trials, fixed_subspace.fit_fixed inline  ->  core.pca_fit_flat / pca_scores
# ---------------------------------------------------------------------------
def _ref_pca_scores_driver(M, k):
    """Verbatim scripts/run_lag_curves.py:83 == run_lag_cosine.py:82 (pre-2026-08-19)."""
    mu = M.mean(0)
    _, _, vt = np.linalg.svd(M - mu, full_matrices=False)
    return (M - mu) @ vt[:k].T


def _ref_scores_lag_subspace(M, k, train):
    """Verbatim src/tom_cca/lag_subspace.py:78 (pre-2026-08-19)."""
    k = int(min(k, M.shape[1], max(1, np.count_nonzero(train) - 1)))
    mu = M[train].mean(0)
    _, _, vt = np.linalg.svd(M[train] - mu, full_matrices=False)
    comp = vt[:k].T
    return (M - mu) @ comp, comp


def _ref_pca_fit_window(M, k):
    """Verbatim src/tom_cca/subspace_window.py:60 == early_trials.py:40 (pre-2026-08-19)."""
    mu = M.mean(0)
    _, _, vt = np.linalg.svd(M - mu, full_matrices=False)
    return mu, vt[:k].T


@pytest.mark.parametrize("seed", range(5))
def test_pca_scores_matches_every_old_copy(seed):
    rng = np.random.default_rng(seed)
    n, p, k = 400, 23, 30 if seed % 2 else 10          # k > p and k < p
    M = rng.standard_normal((n, p)) @ rng.standard_normal((p, p))
    # driver copy (all rows)
    S_new, comp_new = core.pca_scores(M, min(k, p))
    S_ref = _ref_pca_scores_driver(M, min(k, p))
    np.testing.assert_allclose(S_new, S_ref, rtol=RTOL, atol=ATOL)
    # window / early_trials copy
    mu_ref, comp_ref = _ref_pca_fit_window(M, min(k, p))
    st = core.pca_fit_flat(M, min(k, p))
    np.testing.assert_allclose(st.mean, mu_ref, rtol=RTOL, atol=ATOL)
    np.testing.assert_allclose(st.components, comp_ref, rtol=RTOL, atol=ATOL)
    # train-rows-only copy (lag_subspace)
    train = rng.uniform(size=n) < 0.7
    S_new, comp_new = core.pca_scores(M, k, train=train)
    S_ref, comp_ref = _ref_scores_lag_subspace(M, k, train)
    np.testing.assert_allclose(S_new, S_ref, rtol=RTOL, atol=ATOL)
    np.testing.assert_allclose(comp_new, comp_ref, rtol=RTOL, atol=ATOL)
    # the old module-level names still resolve
    S_ls, _ = lag_subspace._scores(M, k, train)
    np.testing.assert_allclose(S_ls, S_ref, rtol=RTOL, atol=ATOL)


def test_pca_scores_k_clip_rule():
    """k is clipped to min(k, n_units, n_train - 1); never raises for small train sets."""
    rng = np.random.default_rng(0)
    M = rng.standard_normal((6, 10))
    _, comp = core.pca_scores(M, 30)
    assert comp.shape == (10, 5)           # n_train - 1
    _, comp = core.pca_scores(M, 3)
    assert comp.shape == (10, 3)


# ---------------------------------------------------------------------------
# 2.1  Within-trial lag pairing: lagged._segment_lagged_pairs, lag_subspace.segment_lag,
#      fixed_subspace.variate_lag_curve / trial_lag_moments  ->  lagpairs.*
# ---------------------------------------------------------------------------
def _ref_segment_lagged_pairs(Sx, Sy, groups, lag):
    """Verbatim src/tom_cca/lagged.py:145 (pre-2026-08-19)."""
    Xs, Ys, gs = [], [], []
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        xt, yt = Sx[idx], Sy[idx]
        n = xt.shape[0]
        if n <= abs(lag) + 2:
            continue
        if lag >= 0:
            xp, yp = xt[: n - lag], yt[lag:]
        else:
            xp, yp = xt[-lag:], yt[: n + lag]
        Xs.append(xp); Ys.append(yp); gs.append(np.full(xp.shape[0], g))
    if not Xs:
        return None
    return np.vstack(Xs), np.vstack(Ys), np.concatenate(gs)


def _ref_segment_lag(X, Y, Z, groups, lag: int):
    """Verbatim src/tom_cca/lag_subspace.py:48 (pre-2026-08-19)."""
    lag = int(lag)
    Xs, Ys, Zxs, Zys, gs = [], [], [], [], []
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        n = idx.size
        if n <= abs(lag) + 2:
            continue
        if lag >= 0:
            xi, yi = idx[: n - lag], idx[lag:]
        else:
            xi, yi = idx[-lag:], idx[: n + lag]
        Xs.append(X[xi]); Ys.append(Y[yi]); gs.append(np.full(xi.size, g))
        if Z is not None:
            Zxs.append(Z[xi]); Zys.append(Z[yi])
    if not Xs:
        return None
    return (np.vstack(Xs), np.vstack(Ys),
            np.vstack(Zxs) if Zxs else None,
            np.vstack(Zys) if Zys else None,
            np.concatenate(gs))


def _ref_variate_lag_curve(u, v, groups, lags) -> np.ndarray:
    """Verbatim src/tom_cca/fixed_subspace.py:123 (pre-2026-08-19)."""
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    groups = np.asarray(groups)
    out = np.full(len(list(lags)), np.nan)
    for i, lag in enumerate(lags):
        lag = int(lag)
        us, vs = [], []
        for g in np.unique(groups):
            idx = np.where(groups == g)[0]
            n = idx.size
            if n <= abs(lag) + 2:
                continue
            xi, yi = (idx[: n - lag], idx[lag:]) if lag >= 0 else \
                     (idx[-lag:], idx[: n + lag])
            us.append(u[xi]); vs.append(v[yi])
        if not us:
            continue
        a, b = np.concatenate(us), np.concatenate(vs)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 3:
            continue
        a, b = a[ok], b[ok]
        if np.std(a) == 0 or np.std(b) == 0:
            continue
        out[i] = float(np.corrcoef(a, b)[0, 1])
    return out


def _ref_trial_lag_moments(u, v, groups, lags):
    """Verbatim src/tom_cca/fixed_subspace.py:356 (pre-2026-08-19)."""
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    groups = np.asarray(groups)
    lags = np.asarray(list(lags), dtype=int)
    trials = np.unique(groups)
    M = np.zeros((trials.size, lags.size, 6))
    idx_by_trial = [np.where(groups == t)[0] for t in trials]
    for ti, idx in enumerate(idx_by_trial):
        n = idx.size
        ut, vt = u[idx], v[idx]
        for li, lag in enumerate(lags):
            lag = int(lag)
            if n <= abs(lag) + 2:
                continue
            a = ut[: n - lag] if lag >= 0 else ut[-lag:]
            b = vt[lag:] if lag >= 0 else vt[: n + lag]
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < 1:
                continue
            a, b = a[ok], b[ok]
            M[ti, li] = (a.size, a.sum(), b.sum(), float(a @ b),
                         float(a @ a), float(b @ b))
    return trials, M


LAGS = list(range(-8, 9))


@pytest.mark.parametrize("lag", LAGS)
def test_lag_pair_indices_reproduce_segment_lagged_pairs(flat_pair, lag):
    X, Y, _, groups = flat_pair
    ref = _ref_segment_lagged_pairs(X, Y, groups, lag)
    new = lagged._segment_lagged_pairs(X, Y, groups, lag)
    if ref is None:
        assert new is None
        return
    for a, b in zip(ref, new):
        np.testing.assert_array_equal(a, b)
    # and the primitive itself yields the same rows in the same order
    ix, iy = lagpairs.lag_pair_indices(groups, lag)
    np.testing.assert_array_equal(X[ix], ref[0])
    np.testing.assert_array_equal(Y[iy], ref[1])
    np.testing.assert_array_equal(groups[ix], ref[2])
    np.testing.assert_array_equal(groups[iy], ref[2])   # never crosses a trial


@pytest.mark.parametrize("lag", LAGS)
@pytest.mark.parametrize("with_z", [True, False])
def test_segment_lag_unchanged(flat_pair, lag, with_z):
    X, Y, Z, groups = flat_pair
    Zc = Z if with_z else None
    ref = _ref_segment_lag(X, Y, Zc, groups, lag)
    new = lag_subspace.segment_lag(X, Y, Zc, groups, lag)
    if ref is None:
        assert new is None
        return
    for a, b in zip(ref, new):
        if a is None:
            assert b is None
        else:
            np.testing.assert_array_equal(a, b)


def test_variate_lag_curve_and_moments_unchanged(flat_pair):
    X, Y, _, groups = flat_pair
    rng = np.random.default_rng(11)
    u = X[:, 0] + 0.3 * np.roll(Y[:, 0], -2)
    v = Y[:, 0].copy()
    u[rng.integers(0, u.size, 5)] = np.nan                  # some non-finite bins
    ref = _ref_variate_lag_curve(u, v, groups, LAGS)
    new = fixed_subspace.variate_lag_curve(u, v, groups, LAGS)
    np.testing.assert_allclose(new, ref, rtol=RTOL, atol=ATOL, equal_nan=True)
    t_ref, M_ref = _ref_trial_lag_moments(u, v, groups, LAGS)
    t_new, M_new = fixed_subspace.trial_lag_moments(u, v, groups, LAGS)
    np.testing.assert_array_equal(t_new, t_ref)
    np.testing.assert_allclose(M_new, M_ref, rtol=RTOL, atol=ATOL)
    # the docstring claim at fixed_subspace.curve_from_moments: moments reproduce the curve
    curve = fixed_subspace.curve_from_moments(M_new)
    np.testing.assert_allclose(curve, new, rtol=1e-9, atol=1e-12, equal_nan=True)


def test_lag_pair_indices_guard_and_sign():
    """n <= |lag| + min_extra is skipped; positive lag pairs X[t] with Y[t+lag]."""
    groups = np.array([0] * 5 + [1] * 3)
    ix, iy = lagpairs.lag_pair_indices(groups, 2)
    np.testing.assert_array_equal(ix, [0, 1, 2])            # trial 1 (n=3) skipped
    np.testing.assert_array_equal(iy, [2, 3, 4])
    ix, iy = lagpairs.lag_pair_indices(groups, -2)
    np.testing.assert_array_equal(ix, [2, 3, 4])
    np.testing.assert_array_equal(iy, [0, 1, 2])
    ix, iy = lagpairs.lag_pair_indices(groups, 2, min_extra=0)
    np.testing.assert_array_equal(ix, [0, 1, 2, 5])
    np.testing.assert_array_equal(iy, [2, 3, 4, 7])
    ix, iy = lagpairs.lag_pair_indices(groups, 9)
    assert ix.size == 0 and iy.size == 0


# ---------------------------------------------------------------------------
# 2.7  IFI sides: the clipped-mean re-derivations in perdim_ifi / run_ifi_windows
#      -> lagged.ifi_sides
# ---------------------------------------------------------------------------
def test_ifi_sides_matches_inline_derivation_and_ifi():
    rng = np.random.default_rng(5)
    lags = np.arange(-10, 11)
    for _ in range(20):
        cc = rng.normal(0.05, 0.1, lags.size)
        pos_ref = np.clip(cc[lags > 0], 0.0, None)
        neg_ref = np.clip(cc[lags < 0], 0.0, None)
        pm_ref = np.nanmean(pos_ref) if np.any(np.isfinite(pos_ref)) else 0.0
        nm_ref = np.nanmean(neg_ref) if np.any(np.isfinite(neg_ref)) else 0.0
        pm, nm = lagged.ifi_sides(lags, cc)
        assert pm == pytest.approx(pm_ref, abs=0) and nm == pytest.approx(nm_ref, abs=0)
        ifi = lagged.information_flow_index(lags, cc)
        tot = pm + nm
        assert ifi == pytest.approx(0.0 if tot <= 0 else (pm - nm) / tot, abs=0)
