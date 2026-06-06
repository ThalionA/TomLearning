"""Full per-window subspace readout for the continuous-regime (p)CCA.

Given one window's residualised activity (the caller does the pCCA control), this
computes the readouts both reference papers report, not just the dominant CC:

* ``cc``          -- held-out canonical correlation per dimension (whole-group CV)
* ``n_sig``       -- number of significant dimensions (circular-shift surrogate)
* ``mi_sig``      -- Gaussian-MI surrogate over the significant dims (strength)
* ``ifi``         -- information-flow index from the lagged curve (direction)
* ``optimal_lag`` -- bin lag maximising the dominant-dim correlation (direction)
* ``gini_x/y``    -- subspace-weight sparsity per area (membership concentration)
* ``weights_x/y`` -- neuron-space canonical weights (for cross-window angles)
* ``member_x/y``  -- top-quartile member masks (for cross-window Jaccard overlap)

Reuses the tested primitives in ``core``, ``membership``, ``lagged``. The driver
computes principal-angle rotation and member overlap *between* windows from the
returned weights/masks. Lags act on the concatenated running bins (small lags vs
thousands of bins -> trial-boundary leakage is negligible); noted as a caveat.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import core, lagged, membership, subspace

CC_CLIP = 0.999


@dataclass
class WindowSubspace:
    cc: np.ndarray            # (d,) held-out per-dim canonical correlation
    n_sig: int
    mi_sig: float
    ifi: float
    optimal_lag: int
    lags: np.ndarray
    lag_cc1: np.ndarray
    gini_x: float
    gini_y: float
    weights_x: np.ndarray     # (n_units_x, d)
    weights_y: np.ndarray     # (n_units_y, d)
    member_x: np.ndarray      # (n_units_x,) bool
    member_y: np.ndarray      # (n_units_y,) bool
    n_units_x: int
    n_units_y: int
    split_half_x: float       # within-window split-half angle (deg) — rotation noise floor
    split_half_y: float


def _scores(M: np.ndarray, k: int):
    mu = M.mean(0)
    _, _, vt = np.linalg.svd(M - mu, full_matrices=False)
    comp = vt[:k].T                                  # (n_units, k)
    return (M - mu) @ comp, comp


def _heldout_perdim(Sx, Sy, groups, n_folds, seed):
    uniq = np.unique(groups)
    if uniq.size < n_folds:
        return None
    rng = np.random.default_rng(seed)
    folds = np.array_split(rng.permutation(uniq), n_folds)
    per = []
    for test_g in folds:
        te = np.isin(groups, test_g)
        tr = ~te
        if np.sum(tr) < 3:
            continue
        model = core.cca_fit(Sx[tr], Sy[tr])
        per.append(core.cca_score(Sx[te], Sy[te], model))
    if not per:
        return None
    d = min(p.size for p in per)
    return np.nanmean(np.array([p[:d] for p in per]), axis=0)


def _insample_perdim(Sx, Sy):
    model = core.cca_fit(Sx, Sy)
    return core.cca_score(Sx, Sy, model)


def _significance(Sx, Sy, cc_heldout, n_shuffles, alpha, seed):
    """Number of significant dims = held-out CC per dim exceeding a circular-shift
    threshold. The null is the distribution of the *dominant* canonical
    correlation under circular shift (the hardest bar); a dimension counts as
    significant only if its cross-validated CC clears that scalar threshold.
    Comparing the held-out CC (unbiased) to the dominant-dim shuffle avoids the
    overcounting of an in-sample, per-dim test.
    """
    d = cc_heldout.size
    if n_shuffles < 1 or d == 0:
        return np.zeros(d, dtype=bool)
    rng = np.random.default_rng(seed + 7)
    n = Sx.shape[0]
    null_top = np.empty(n_shuffles)
    for s in range(n_shuffles):
        shift = int(rng.integers(1, n))
        rs = _insample_perdim(Sx, np.roll(Sy, shift, axis=0))
        null_top[s] = rs[0] if rs.size else 0.0
    thr = float(np.quantile(null_top, 1 - alpha))
    return np.asarray(cc_heldout) > thr


def _lag_curve(Sx, Sy, lags):
    """Dominant-dim in-sample correlation as a function of integer bin lag.

    Positive lag = Y delayed relative to X (X leads). In-sample is sufficient
    for the *shape* used by IFI / optimal lag; both signs are equally biased.
    """
    n = Sx.shape[0]
    out = np.full(lags.size, np.nan)
    for i, L in enumerate(lags):
        if L >= 0:
            a, b = Sx[:n - L] if L else Sx, Sy[L:] if L else Sy
        else:
            a, b = Sx[-L:], Sy[:n + L]
        if a.shape[0] < 10:
            continue
        out[i] = _insample_perdim(a, b)[0]
    return out


def _max_angle(wa, wb, d_use):
    d = int(min(d_use, wa.shape[1], wb.shape[1]))
    if d < 1:
        return float("nan")
    qa, _ = np.linalg.qr(wa[:, :d])
    qb, _ = np.linalg.qr(wb[:, :d])
    return float(np.degrees(np.nanmax(subspace.principal_angles(qa, qb))))


def _half_weights(X, Y, mask, k, d_use):
    if np.sum(mask) < 3 * (d_use + 1):
        return None, None
    Sx, cx = _scores(X[mask], int(min(k, X.shape[1])))
    Sy, cy = _scores(Y[mask], int(min(k, Y.shape[1])))
    model = core.cca_fit(Sx, Sy)
    d = int(min(d_use, cx.shape[1], cy.shape[1], model.A.shape[1]))
    return (membership.canonical_weight_scores(cx, model.A, d),
            membership.canonical_weight_scores(cy, model.B, d))


def _split_half_angles(X, Y, groups, k, d_use, seed):
    """Rotation NOISE FLOOR: max principal angle between subspaces fit on two
    random halves of the window's trials. A genuine cross-window rotation only
    means reorientation if it exceeds this within-window floor."""
    uniq = np.unique(groups)
    if uniq.size < 4:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed + 11)
    perm = rng.permutation(uniq)
    h1 = perm[: uniq.size // 2]
    m1 = np.isin(groups, h1)
    wx1, wy1 = _half_weights(X, Y, m1, k, d_use)
    wx2, wy2 = _half_weights(X, Y, ~m1, k, d_use)
    if wx1 is None or wx2 is None:
        return float("nan"), float("nan")
    return _max_angle(wx1, wx2, d_use), _max_angle(wy1, wy2, d_use)


def window_subspace(X, Y, groups, k: int = 30, max_lag: int = 10,
                    n_shuffles: int = 30, alpha: float = 0.05,
                    n_folds: int = 5, member_q: float = 0.75,
                    seed: int = 0) -> WindowSubspace:
    """Full subspace readout for one window. ``X``/``Y`` are (n_samples, n_units),
    already residualised on the control area(s) by the caller."""
    kx = int(min(k, X.shape[1]))
    ky = int(min(k, Y.shape[1]))
    Sx, cx = _scores(X, kx)
    Sy, cy = _scores(Y, ky)
    d = min(kx, ky)

    cc = _heldout_perdim(Sx, Sy, groups, n_folds, seed)
    if cc is None:
        cc = np.full(d, np.nan)
    sig = _significance(Sx, Sy, np.nan_to_num(cc, nan=-1.0), n_shuffles, alpha, seed)
    sig = sig[:cc.size]
    n_sig = int(np.sum(sig))
    cc_clip = np.clip(cc, 0.0, CC_CLIP)
    if n_sig > 0:
        s = cc_clip[sig]
        s = s[np.isfinite(s) & (s > 0)]
        mi_sig = float(-0.5 * np.sum(np.log1p(-s ** 2)))
    else:
        mi_sig = 0.0

    lags = np.arange(-max_lag, max_lag + 1)
    lag_cc1 = _lag_curve(Sx, Sy, lags)
    ifi = lagged.information_flow_index(lags, lag_cc1)
    optimal_lag = int(lags[int(np.nanargmax(lag_cc1))]) if np.any(np.isfinite(lag_cc1)) else 0

    model = core.cca_fit(Sx, Sy)
    wx = membership.canonical_weight_scores(cx, model.A, d)   # (n_units_x, d)
    wy = membership.canonical_weight_scores(cy, model.B, d)
    contrib_x = membership.subspace_contribution(wx)
    contrib_y = membership.subspace_contribution(wy)
    sh_x, sh_y = _split_half_angles(X, Y, groups, k, d_use=3, seed=seed)
    return WindowSubspace(
        cc=cc, n_sig=n_sig, mi_sig=mi_sig, ifi=ifi, optimal_lag=optimal_lag,
        lags=lags, lag_cc1=lag_cc1,
        gini_x=membership.gini(contrib_x), gini_y=membership.gini(contrib_y),
        weights_x=wx, weights_y=wy,
        member_x=membership.member_mask(contrib_x, member_q),
        member_y=membership.member_mask(contrib_y, member_q),
        n_units_x=X.shape[1], n_units_y=Y.shape[1],
        split_half_x=sh_x, split_half_y=sh_y,
    )
