"""Lagged communication subspaces — the shared machinery for meeting 2026-07-28
items 2, 3 and 4.

Everything upstream reads canonical *correlations* at a lag (`lagged`,
`run_lag_curves`); nothing exports the *weights*, so the subspace itself has never been
compared across lags. This module fits CCA at a fixed segment-aware lag and returns
neuron-space weights alongside the held-out per-dim CC, which buys three things:

* **item 3** — principal angles between the lag-0 subspace and each lagged subspace,
  against a split-half noise floor: *how stable are the CCs across time lags?*
* **items 2/4** — a feedforward subspace (fit at +tau, X leads) and a feedback subspace
  (fit at -tau, Y leads) that can be tracked separately across learning. Defining FF/FB
  by the *fit lag* rather than by the sign of a per-dim IFI avoids the noise-splitting
  hazard found in `analyze_perdim_ifi` (the significant tail is dominated by dims sitting
  at the held-out CC floor, so their IFI sign is close to a coin flip).

Per the temporal-arm convention (`UNDERSTANDING_temporal.md`) this is a NEW module —
`subspace_window.window_subspace` is not modified.

Sign convention, inherited from `lagged.lag_slice`: **positive lag => X leads Y**.

Leakage policy matches the rest of the temporal arm: the confound regression, the PCA and
the CCA are all fit on the training trials of each fold only, and the CC is read on
held-out trials. The exported weights are the full-window (descriptive) fit — they
describe subspace *geometry*, and are never used to produce a performance number.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import core, membership, partial, subspace


@dataclass
class LagFit:
    """One CCA fit at a fixed lag."""

    lag: int
    cc: np.ndarray          # (d,) held-out per-dim canonical correlation
    wx: np.ndarray          # (n_units_x, d) neuron-space canonical weights
    wy: np.ndarray          # (n_units_y, d)
    n_samples: int          # segment-aware pairs surviving the lag


def segment_lag(X, Y, Z, groups, lag: int):
    """Pair ``X[t]`` with ``Y[t + lag]`` WITHIN each trial only.

    Returns ``(Xp, Yp, Zx, Zy, groups)`` or ``None`` if no trial is long enough. The
    confound is returned twice, aligned to each area's own timepoints: partialling a
    third area out of ``Y`` using ``Z`` at ``X``'s times would inject exactly the lag
    the analysis is trying to measure. ``Zx``/``Zy`` are ``None`` when ``Z`` is ``None``.
    """
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


def _scores(M, k, train):
    """PCA fit on training rows only, applied to all rows -> (scores, components)."""
    k = int(min(k, M.shape[1], max(1, np.count_nonzero(train) - 1)))
    mu = M[train].mean(0)
    _, _, vt = np.linalg.svd(M[train] - mu, full_matrices=False)
    comp = vt[:k].T                                  # (n_units, k)
    return (M - mu) @ comp, comp


def lagged_fit(X, Y, groups, Z=None, lag: int = 0, k: int = 30, n_folds: int = 5,
               seed: int = 0) -> LagFit | None:
    """Leak-free CCA at a fixed segment-aware ``lag``; weights in neuron space.

    ``X``/``Y`` are ``(n_samples, n_units)`` raw activity, ``groups`` the per-bin trial
    id, ``Z`` the optional third-area confound. Returns ``None`` when the lag leaves no
    usable trial or no fold is fittable.
    """
    packed = segment_lag(X, Y, Z, groups, lag)
    if packed is None:
        return None
    Xp, Yp, Zx, Zy, gp = packed
    uniq = np.unique(gp)
    if uniq.size < 2:
        return None
    rng = np.random.default_rng(seed)
    folds = [f for f in np.array_split(rng.permutation(uniq), n_folds) if f.size]

    per_fold = []
    for te in folds:
        tem = np.isin(gp, te)
        trm = ~tem
        if trm.sum() < max(5, k + 1) or tem.sum() < 3:
            continue
        xr = partial.partial_out_cv(Xp, Zx, trm) if Zx is not None else Xp
        yr = partial.partial_out_cv(Yp, Zy, trm) if Zy is not None else Yp
        sx, _ = _scores(xr, k, trm)
        sy, _ = _scores(yr, k, trm)
        model = core.cca_fit(sx[trm], sy[trm])
        r = core.cca_score(sx[tem], sy[tem], model)
        if r.size:
            per_fold.append(np.asarray(r, dtype=float))
    if not per_fold:
        return None
    d = min(r.size for r in per_fold)
    cc = np.nanmean([r[:d] for r in per_fold], axis=0)

    # descriptive full-window fit -> the geometry we compare across lags
    allrows = np.ones(Xp.shape[0], dtype=bool)
    xr = partial.partial_out(Xp, Zx) if Zx is not None else Xp
    yr = partial.partial_out(Yp, Zy) if Zy is not None else Yp
    sx, cx = _scores(xr, k, allrows)
    sy, cy = _scores(yr, k, allrows)
    model = core.cca_fit(sx, sy)
    d_w = int(min(d, model.A.shape[1], model.B.shape[1]))
    wx = membership.canonical_weight_scores(cx, model.A, d_w)
    wy = membership.canonical_weight_scores(cy, model.B, d_w)
    return LagFit(lag=int(lag), cc=cc, wx=wx, wy=wy, n_samples=int(Xp.shape[0]))


def lag_sweep(X, Y, groups, Z=None, lags=(0,), k: int = 30, n_folds: int = 5,
              seed: int = 0) -> dict[int, LagFit]:
    """:func:`lagged_fit` over many lags, sharing one residualisation AND one PCA basis
    per fold.

    Both the confound regression and the PCA are fit on the training trials of the
    UNLAGGED data and applied to all rows; lagging then only selects a subset of those
    rows, so neither set of coefficients changes and the fold is still leak-free. This
    is the same hoist `run_lag_curves` already makes (PCA once, CCA held-out per lag),
    and it matters three ways:

    * cost — the per-lag work drops to a CCA on ``k``-dimensional scores. Refitting the
      PCA inside the lag loop meant ``n_lags x n_folds x 2`` SVDs of a
      ``n_samples x n_units`` matrix, which dominated everything else.
    * comparability — every lag is expressed in ONE neuron->PC basis, so a cross-lag
      principal angle reflects the CCA rotation rather than a re-estimated PCA. With a
      per-lag PCA the item-3 readout would partly measure PCA jitter.
    * the confound is removed identically at every lag, so a cross-lag difference cannot
      be an artefact of re-fitting the third-area regression.

    Returns ``{lag: LagFit}``, skipping lags with no usable trial.
    """
    uniq = np.unique(groups)
    if uniq.size < 2:
        return {}
    rng = np.random.default_rng(seed)
    folds = [f for f in np.array_split(rng.permutation(uniq), n_folds) if f.size]

    prepared = []
    for te in folds:
        tem = np.isin(groups, te)
        trm = ~tem
        if trm.sum() < max(5, k + 1) or tem.sum() < 3:
            continue
        xr = partial.partial_out_cv(X, Z, trm) if Z is not None else X
        yr = partial.partial_out_cv(Y, Z, trm) if Z is not None else Y
        sx, _ = _scores(xr, k, trm)
        sy, _ = _scores(yr, k, trm)
        prepared.append((sx, sy, tem))
    if not prepared:
        return {}
    allrows = np.ones(X.shape[0], dtype=bool)
    xr_all = partial.partial_out(X, Z) if Z is not None else X
    yr_all = partial.partial_out(Y, Z) if Z is not None else Y
    sx_all, cx = _scores(xr_all, k, allrows)
    sy_all, cy = _scores(yr_all, k, allrows)

    out: dict[int, LagFit] = {}
    for lag in lags:
        per_fold = []
        for sx, sy, tem in prepared:
            packed = segment_lag(sx, sy, None, groups, lag)
            if packed is None:
                break
            xp, yp, _, _, gp = packed
            # the lag drops edge bins, so rebuild the fold mask on surviving rows
            tep = np.isin(gp, np.unique(groups[tem]))
            trp = ~tep
            if trp.sum() < max(5, k + 1) or tep.sum() < 3:
                continue
            model = core.cca_fit(xp[trp], yp[trp])
            r = core.cca_score(xp[tep], yp[tep], model)
            if r.size:
                per_fold.append(np.asarray(r, dtype=float))
        if not per_fold:
            continue
        d = min(r.size for r in per_fold)
        cc = np.nanmean([r[:d] for r in per_fold], axis=0)

        packed = segment_lag(sx_all, sy_all, None, groups, lag)
        if packed is None:
            continue
        xp, yp, _, _, _ = packed
        model = core.cca_fit(xp, yp)
        d_w = int(min(model.A.shape[1], model.B.shape[1]))
        out[int(lag)] = LagFit(
            lag=int(lag), cc=cc,
            wx=membership.canonical_weight_scores(cx, model.A, d_w),
            wy=membership.canonical_weight_scores(cy, model.B, d_w),
            n_samples=int(xp.shape[0]))
    return out


def subspace_angle(wa, wb, d_use: int) -> float:
    """Largest principal angle (degrees) between the spans of ``wa`` and ``wb``.

    0 = the same subspace, 90 = orthogonal. Uses the LARGEST angle (the worst-aligned
    direction) so a subspace that matches on its dominant axis but diverges on the rest
    is not scored as stable. Scale- and sign-invariant, since only the span matters.
    """
    wa, wb = np.asarray(wa, dtype=float), np.asarray(wb, dtype=float)
    d = int(min(d_use, wa.shape[1] if wa.ndim > 1 else 0,
                wb.shape[1] if wb.ndim > 1 else 0))
    if d < 1:
        return float("nan")
    qa, _ = np.linalg.qr(wa[:, :d])
    qb, _ = np.linalg.qr(wb[:, :d])
    ang = subspace.principal_angles(qa, qb)
    if ang.size == 0 or not np.any(np.isfinite(ang)):
        return float("nan")
    return float(np.degrees(np.nanmax(ang)))


def split_half_floor(X, Y, groups, Z=None, lag: int = 0, k: int = 30, d_use: int = 3,
                     seed: int = 0) -> tuple[float, float]:
    """Noise floor for :func:`subspace_angle`: the angle between two halves of the
    trials fit at the SAME lag. Returns ``(floor_x, floor_y)``.

    A cross-lag angle only means "the subspace moved" if it exceeds this. Splitting by
    whole trials keeps the comparison honest for autocorrelated within-trial bins.

    The two areas get SEPARATE floors on purpose: the floor rises with fewer units and
    fewer samples, and the pairs here are lopsided (CA1 carries many more units than
    SUB), so a single shared floor would set the bar too high for one area and too low
    for the other — and the floor is the threshold the item-3 verdict turns on.
    """
    nan2 = (float("nan"), float("nan"))
    packed = segment_lag(X, Y, Z, groups, lag)
    if packed is None:
        return nan2
    Xp, Yp, Zx, Zy, gp = packed
    uniq = np.unique(gp)
    if uniq.size < 4:
        return nan2
    rng = np.random.default_rng(seed + 11)
    perm = rng.permutation(uniq)
    halves = [np.isin(gp, perm[: perm.size // 2]), np.isin(gp, perm[perm.size // 2:])]
    ws = []
    for m in halves:
        if m.sum() < max(10, k + 1):
            return nan2
        xr = partial.partial_out(Xp[m], Zx[m]) if Zx is not None else Xp[m]
        yr = partial.partial_out(Yp[m], Zy[m]) if Zy is not None else Yp[m]
        allrows = np.ones(xr.shape[0], dtype=bool)
        sx, cx = _scores(xr, k, allrows)
        sy, cy = _scores(yr, k, allrows)
        model = core.cca_fit(sx, sy)
        d_w = int(min(model.A.shape[1], model.B.shape[1]))
        ws.append((membership.canonical_weight_scores(cx, model.A, d_w),
                   membership.canonical_weight_scores(cy, model.B, d_w)))
    return (subspace_angle(ws[0][0], ws[1][0], d_use),
            subspace_angle(ws[0][1], ws[1][1], d_use))
