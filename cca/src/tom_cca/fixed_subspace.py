"""FIXED communication subspace — meeting 2026-07-28, items 5 and 7.

    "Go back to also doing one CC and lag it across time."
    "Check with fixed subspaces, identified across all trials, and compare naive/exp."

Every epoch contrast in the project so far REFITS CCA within each epoch, so an
epoch-to-epoch difference can be a difference in the *fit* (different weights, different
effective dimensionality, different overfitting) rather than in the *activity*. This
module identifies the subspace once across all trials, freezes it, and then projects
each epoch through the identical weights. One canonical component then has a single
time series per epoch that can be lagged across time and compared directly.

**On fitting across all trials.** The frozen fit sees every epoch, so the absolute
correlation it yields is optimistic — it is NOT a held-out number and must never be
quoted as one. What it is good for is the *contrast*: because both epochs pass through
the same weights, neither is favoured. That symmetry only holds if the fit is not
dominated by one epoch, which is why :func:`balanced_trials` equalises the trial count
per epoch before fitting. Use `lag_subspace` when a leak-free CC is what is wanted.

Sign convention, as everywhere in the temporal arm: positive lag => X leads Y.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import core, membership, partial


@dataclass
class FixedFit:
    """A subspace identified once and reused across epochs."""

    wx: np.ndarray          # (n_units_x, d) neuron-space canonical weights
    wy: np.ndarray          # (n_units_y, d)
    x_mean: np.ndarray      # (n_units_x,) centring used at fit time
    y_mean: np.ndarray      # (n_units_y,)
    r_fit: np.ndarray       # (d,) IN-SAMPLE canonical correlations — not held out
    n_trials: int
    # PCA bases the CCA was fitted through. Exposed so a caller can build the SAME
    # PC scores the fit used — needed to compute per-dimension significance on the
    # frozen fit without re-deriving (and possibly diverging from) the fit's own maths.
    cx: np.ndarray | None = None
    cy: np.ndarray | None = None


def balanced_trials(epoch_of_trial: dict, epochs, seed: int = 0) -> list:
    """Equal numbers of trials per epoch, so no epoch dominates the frozen fit.

    Takes the smallest epoch's trial count and draws that many from each, deterministic
    per ``seed``. Returns ``[]`` if any epoch is absent — an unbalanced fit would make
    the naive-vs-expert contrast asymmetric, which is exactly what this arm exists to
    avoid.
    """
    by_epoch = {e: sorted(t for t, lab in epoch_of_trial.items() if lab == e)
                for e in epochs}
    if any(len(v) == 0 for v in by_epoch.values()):
        return []
    n = min(len(v) for v in by_epoch.values())
    rng = np.random.default_rng(seed)
    out = []
    for e in epochs:
        pool = np.asarray(by_epoch[e])
        out.extend(sorted(rng.permutation(pool)[:n].tolist()))
    return sorted(out)


def fit_fixed(X, Y, groups, Z=None, k: int = 30, trials=None,
              seed: int = 0) -> FixedFit | None:
    """Identify the subspace once, on ``trials`` (default: all), and return the weights.

    ``X``/``Y`` are ``(n_samples, n_units)`` raw activity; ``Z`` the optional third-area
    confound, partialled out with in-sample coefficients (this is a descriptive fit by
    construction — see the module docstring).
    """
    del seed
    sel = np.ones(X.shape[0], dtype=bool) if trials is None else np.isin(groups, trials)
    if sel.sum() < max(10, k + 1):
        return None
    Xs, Ys = X[sel], Y[sel]
    Zs = Z[sel] if Z is not None else None
    xr = partial.partial_out(Xs, Zs) if Zs is not None else Xs
    yr = partial.partial_out(Ys, Zs) if Zs is not None else Ys
    kx = int(min(k, xr.shape[1], sel.sum() - 1))
    ky = int(min(k, yr.shape[1], sel.sum() - 1))
    if kx < 1 or ky < 1:
        return None
    x_mean, y_mean = xr.mean(0), yr.mean(0)
    _, sx, vtx = np.linalg.svd(xr - x_mean, full_matrices=False)
    _, sy, vty = np.linalg.svd(yr - y_mean, full_matrices=False)
    if not np.any(sx > 0) or not np.any(sy > 0):
        return None
    cx, cy = vtx[:kx].T, vty[:ky].T
    try:
        model = core.cca_fit((xr - x_mean) @ cx, (yr - y_mean) @ cy)
    except Exception:                                       # noqa: BLE001
        return None
    d = int(min(model.A.shape[1], model.B.shape[1]))
    if d < 1:
        return None
    return FixedFit(
        wx=membership.canonical_weight_scores(cx, model.A, d),
        wy=membership.canonical_weight_scores(cy, model.B, d),
        x_mean=x_mean, y_mean=y_mean,
        r_fit=np.asarray(model.r, dtype=float)[:d],
        n_trials=int(np.unique(groups[sel]).size),
        cx=cx, cy=cy)


def project(X, Y, fit: FixedFit, dim: int = 0):
    """``(u, v)`` — the canonical variate time series of one component ``dim``.

    Applies the frozen weights and the fit-time centring, so any subset of samples maps
    through exactly the same transform.
    """
    u = (np.asarray(X, dtype=float) - fit.x_mean) @ fit.wx[:, dim]
    v = (np.asarray(Y, dtype=float) - fit.y_mean) @ fit.wy[:, dim]
    return u, v


def variate_lag_curve(u, v, groups, lags) -> np.ndarray:
    """Pearson correlation of ``u[t]`` with ``v[t + lag]``, pooled WITHIN trials.

    The pairing never crosses a trial boundary, so the running-bin concatenation cannot
    manufacture a correlation at long lag. Returns NaN at a lag where no trial is long
    enough or either variate is constant.
    """
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


def curve_half_width(lags, r) -> float:
    """Width (in lag units) of the CONTIGUOUS region around the peak where the curve
    stays at or above half its peak value — the temporal integration window of the
    coupling.

    Contiguity matters: a far-lag noise excursion above half-max must not widen the
    window, so the region is grown outward from the peak and stops at the first bin
    that drops below the threshold. Returns NaN for a non-positive or all-NaN peak (a
    half-max width is meaningless when the peak correlation is <= 0), and a censored
    width equal to the swept range when the curve never falls below half-max.
    """
    lags = np.asarray(lags, dtype=float)
    r = np.asarray(r, dtype=float)
    ok = np.isfinite(lags) & np.isfinite(r)
    if not np.any(ok):
        return float("nan")
    lags, r = lags[ok], r[ok]
    order = np.argsort(lags, kind="stable")
    lags, r = lags[order], r[order]
    pk = int(np.argmax(r))
    if not np.isfinite(r[pk]) or r[pk] <= 0:
        return float("nan")
    half = r[pk] / 2.0
    lo = hi = pk
    while lo - 1 >= 0 and r[lo - 1] >= half:
        lo -= 1
    while hi + 1 < r.size and r[hi + 1] >= half:
        hi += 1
    return float(lags[hi] - lags[lo])


def side_peak(lags, r):
    """``(lag_offset, height_ratio)`` of the largest SECONDARY peak outside the central
    half-max region — the tell for an oscillatory lag curve.

    Intra-hippocampal canonical components are theta-modulated, so their lag curve rings
    at ~7-8 Hz rather than decaying. That matters for interpretation: on a ringing curve
    :func:`curve_half_width` measures the half-period of the oscillation, NOT a temporal
    integration window, and the IFI compares two lobes of a rhythm rather than a lead.

    Returns ``(nan, nan)`` when the curve has no secondary peak (a genuinely decaying
    curve, where the width does mean an integration window). ``height_ratio`` is the
    secondary peak's height as a fraction of the central peak.
    """
    lags = np.asarray(lags, dtype=float)
    r = np.asarray(r, dtype=float)
    ok = np.isfinite(lags) & np.isfinite(r)
    if ok.sum() < 5:
        return float("nan"), float("nan")
    lags, r = lags[ok], r[ok]
    order = np.argsort(lags, kind="stable")
    lags, r = lags[order], r[order]
    pk = int(np.argmax(r))
    if r[pk] <= 0:
        return float("nan"), float("nan")
    half = r[pk] / 2.0
    lo = hi = pk
    while lo - 1 >= 0 and r[lo - 1] >= half:
        lo -= 1
    while hi + 1 < r.size and r[hi + 1] >= half:
        hi += 1
    best_i, best_v = -1, -np.inf
    for i in range(1, r.size - 1):
        if lo <= i <= hi:
            continue
        # STRICT on both sides: a flat shoulder is not a ring, and with `>=` every
        # point of a plateau qualifies as a local maximum.
        if r[i] > r[i - 1] and r[i] > r[i + 1] and r[i] > best_v:
            best_i, best_v = i, r[i]
    if best_i < 0 or best_v <= 0:
        return float("nan"), float("nan")
    return float(abs(lags[best_i] - lags[pk])), float(best_v / r[pk])


def side_peak_null(lags, n_draws: int = 3000, ratio_gate: float = 0.5,
                   seed: int = 0) -> np.ndarray:
    """Side-peak offsets produced by NOISE alone, on the same lag grid.

    :func:`side_peak` fires on essentially every noisy curve — on a curve whose central
    peak sits near the noise floor almost any bin qualifies as a secondary maximum — so
    the FRACTION of curves with a side peak carries no information and must never be
    quoted as evidence of rhythmicity. What does carry information is WHERE the side
    peaks land: a real oscillation concentrates them at its period, whereas noise spreads
    them broadly across the grid. Compare an observed offset distribution against this
    null with a KS test and/or a band-occupancy contrast.
    """
    rng = np.random.default_rng(seed)
    lags = np.asarray(lags, dtype=float)
    out = []
    for _ in range(int(n_draws)):
        r = rng.standard_normal(lags.size) * 0.02 + 0.05
        ms, ratio = side_peak(lags, r)
        if np.isfinite(ms) and ratio > ratio_gate:
            out.append(ms)
    return np.asarray(out, dtype=float)


def band_occupancy(offsets_ms, lo_hz: float = 5.0, hi_hz: float = 12.0) -> float:
    """Fraction of side-peak offsets whose implied frequency falls in ``[lo_hz, hi_hz]``."""
    o = np.asarray(offsets_ms, dtype=float)
    o = o[np.isfinite(o) & (o > 0)]
    if o.size == 0:
        return float("nan")
    return float(np.mean((o >= 1000.0 / hi_hz) & (o <= 1000.0 / lo_hz)))


@dataclass
class FrozenSignificance:
    """Per-dimension significance of a FROZEN fit's canonical correlations."""

    mask: np.ndarray        # (d,) bool
    p: np.ndarray           # (d,) empirical p-values
    threshold: np.ndarray   # (d,) null quantile at alpha, per rank
    r_obs: np.ndarray       # (d,) observed canonical correlations


def _orth_basis(M):
    """Orthonormal basis for the column space of centred ``M`` (thin QR)."""
    M = np.asarray(M, dtype=float)
    q, _ = np.linalg.qr(M - M.mean(0))
    return q


def observed_canonical_r(X, Y) -> np.ndarray:
    """Canonical correlations of ``X`` vs ``Y`` = singular values of ``Qx.T @ Qy``."""
    qx, qy = _orth_basis(X), _orth_basis(Y)
    return np.clip(np.linalg.svd(qx.T @ qy, compute_uv=False), 0.0, 1.0)


def frozen_perm_null(X, Y, n_shuffles: int = 200, alpha: float = 0.05, seed: int = 0,
                     correct: str | None = "fdr",
                     fdr_dims: int | None = None) -> FrozenSignificance:
    """Circular-shift permutation test on a FROZEN fit's own canonical correlations.

    Dimension *k*'s observed canonical correlation is compared to the shuffled
    distribution of dimension *k*. Observed and null are the SAME estimator — the rank-k
    canonical correlation of a CCA fit to (real / shifted) data — so there is no
    cross-fit correspondence to get wrong.

    That matters here specifically. Gating a frozen fit with statistics taken from
    per-fold CCA refits (which average across fits BY RANK) mis-attributes them: on this
    data the fold ordering matched the frozen ordering in only 68-92 % of instances at
    ranks 2-5, producing both false admits and false excludes. This test cannot have that
    failure because it never mixes two fits.

    Fast because a circular shift is a PERMUTATION: for orthonormal bases, canonical
    correlations are the singular values of ``Qx.T @ Qy``, and ``QR(P @ Y) = P @ Qy`` for
    a permutation ``P`` — so the bases are computed ONCE and each surrogate costs one
    ``(d x n) @ (n x d)`` matmul instead of a full refit. ~20x cheaper, and exact.

    ``fdr_dims`` restricts the BH family to the leading dims; the permutation p floor is
    ``1/(n_shuffles+1)`` and BH needs ``alpha/d``, so the family must be small enough for
    the floor to clear it.
    """
    qx, qy = _orth_basis(X), _orth_basis(Y)
    r_obs = np.clip(np.linalg.svd(qx.T @ qy, compute_uv=False), 0.0, 1.0)
    d = r_obs.size
    n = qy.shape[0]
    if n_shuffles < 1 or d == 0:
        return FrozenSignificance(np.zeros(d, dtype=bool), np.ones(d),
                                  np.full(d, np.nan), r_obs)
    rng = np.random.default_rng(seed + 7)
    null = np.empty((n_shuffles, d))
    for s in range(n_shuffles):
        shift = int(rng.integers(1, max(2, n)))
        m = qx.T @ np.roll(qy, shift, axis=0)
        sv = np.linalg.svd(m, compute_uv=False)
        null[s, :sv.size] = np.clip(sv, 0.0, 1.0)
        if sv.size < d:
            null[s, sv.size:] = 0.0
    p = (1.0 + np.sum(null >= r_obs[None, :], axis=0)) / (1.0 + n_shuffles)
    thr = np.quantile(null, 1 - alpha, axis=0)
    if correct == "fdr":
        fam = np.zeros(d, dtype=bool)
        fam[: (d if fdr_dims is None else min(d, fdr_dims))] = True
        mask = np.zeros(d, dtype=bool)
        mask[fam] = _fdr_bh(p[fam], alpha)
    else:
        mask = p < alpha
    return FrozenSignificance(mask, p, thr, r_obs)


def _fdr_bh(pvals, q: float = 0.05) -> np.ndarray:
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


def trial_lag_moments(u, v, groups, lags):
    """Per-(trial, lag) sufficient statistics for the pooled within-trial correlation.

    Returns ``(trials, M)`` with ``M`` of shape ``(n_trials, n_lags, 6)`` holding
    ``[n, sum_u, sum_v, sum_uv, sum_uu, sum_vv]`` over the segment-aware pairs of that
    trial at that lag.

    These are ADDITIVE over trials, so any subset's pooled Pearson correlation follows
    from summing rows — see :func:`curve_from_moments`. Computing them once turns the
    whole-session curve, each epoch's curve and each leave-one-epoch-out curve into cheap
    aggregations instead of a fresh O(n_lags x n_trials) pass apiece. That is what makes
    the uncapped (full-trial) analysis tractable: the pairing work is done once.
    """
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


def curve_from_moments(M, trial_mask=None) -> np.ndarray:
    """Pooled within-trial Pearson correlation per lag from :func:`trial_lag_moments`.

    ``trial_mask`` selects which trials contribute (default all), so an epoch curve and a
    leave-that-epoch-out curve are two different masks over the same moments. Exactly
    reproduces :func:`variate_lag_curve` on the same trials.
    """
    M = np.asarray(M, dtype=float)
    S = M.sum(axis=0) if trial_mask is None else M[np.asarray(trial_mask, bool)].sum(0)
    n, su, sv, suv, suu, svv = (S[:, i] for i in range(6))
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = suv - su * sv / n
        vu = suu - su * su / n
        vv = svv - sv * sv / n
        r = cov / np.sqrt(vu * vv)
    r = np.where((n >= 3) & (vu > 0) & (vv > 0), r, np.nan)
    return np.clip(r, -1.0, 1.0)
