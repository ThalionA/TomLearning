"""Segment-aware lagged CCA for Arm A (running-state temporal).

Refits CCA at every integer-bin lag in ``[-lag_bins, +lag_bins]``, where the
per-lag design matrix is built by pooling within-segment (X-bin t, Y-bin t+L)
pairs across all segments returned by :func:`tom_cca.segments.find_segments`.
No pair ever straddles a gap or a trial boundary.

Sign convention: positive lag => X leads Y (the X bin is earlier in time;
the matched Y bin is L bins later).

See ``cca/UNDERSTANDING_temporal.md`` Sec. 5.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from . import core, lagged, segments


@dataclass
class TemporalLagResult:
    """Lagged-CCA directionality for one (animal, pair, epoch), Arm A."""

    lags: np.ndarray              # (n_lags,) integer bin lags
    cc_per_dim: np.ndarray        # (n_lags, n_dims) held-out CC1..n at each lag
    ifi_per_dim: np.ndarray       # (n_dims,) IFI over the full lag range
    ifi_windows: np.ndarray       # (n_dims, max_window) IFI by lag window
    peak_lag_per_dim: np.ndarray  # (n_dims,) integer bin lag of per-dim CC peak
    n_paired_samples: np.ndarray  # (n_lags,) pooled sample count at each lag

    @property
    def cc1(self) -> np.ndarray:
        return self.cc_per_dim[:, 0]

    @property
    def ifi(self) -> float:
        return float(self.ifi_per_dim[0])

    @property
    def peak_lag(self) -> int:
        return int(self.peak_lag_per_dim[0])


def _segment_pairs(
    scores: np.ndarray,
    segments_: list[segments.Segment],
    lag: int,
    which: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Pool within-segment lag-L rows from a flat ``(n_bins, k)`` score matrix.

    Returns ``(pooled_scores, trial_labels)``. ``which`` is ``"x"`` (use the
    earlier-time bin of each pair) or ``"y"`` (the later-time bin), matching
    the sign convention.
    """
    rows: list[np.ndarray] = []
    trials: list[np.ndarray] = []
    L = int(lag)
    for s in segments_:
        n = len(s)
        if abs(L) >= n:
            continue
        if L >= 0:
            x_local = np.arange(0, n - L)
            y_local = np.arange(L, n)
        else:
            x_local = np.arange(-L, n)
            y_local = np.arange(0, n + L)
        local = x_local if which == "x" else y_local
        rows.append(scores[s.start + local])
        trials.append(np.full(local.size, s.trial, dtype=int))
    if not rows:
        return (np.zeros((0, scores.shape[1])), np.zeros(0, dtype=int))
    return np.concatenate(rows, axis=0), np.concatenate(trials, axis=0)


def _held_out_cc_lotot(
    x: np.ndarray, y: np.ndarray, trials: np.ndarray, n_dims_target: int,
) -> np.ndarray:
    """Leave-one-trial-out held-out CC per dimension.

    For each unique trial label, fit CCA on the other trials and score on the
    held-out trial; report the mean held-out per-dim correlation. Folds with
    fewer than 3 held-out paired samples are dropped.
    """
    unique = np.unique(trials)
    if unique.size < 2 or x.shape[0] < 6:
        return np.full(n_dims_target, np.nan)
    fold_cc: list[np.ndarray] = []
    for t in unique:
        train = trials != t
        test = ~train
        if not (np.any(train) and np.sum(test) >= 3):
            continue
        try:
            model = core.cca_fit(x[train], y[train])
        except ValueError:
            continue
        cc = core.cca_score(x[test], y[test], model)
        if cc.size == 0:
            continue
        fold_cc.append(cc)
    if not fold_cc:
        return np.full(n_dims_target, np.nan)
    # Pad each fold's CC vector to n_dims_target so we can average per-dim.
    out = np.full((len(fold_cc), n_dims_target), np.nan)
    for i, cc in enumerate(fold_cc):
        m = min(cc.size, n_dims_target)
        out[i, :m] = cc[:m]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN columns
        return np.nanmean(out, axis=0)


def lag_curve(
    scores_x: np.ndarray,
    scores_y: np.ndarray,
    segments_: list[segments.Segment],
    cfg,
    max_lag_bins: int | None = None,
) -> TemporalLagResult:
    """Refit CCA at every lag and summarise directionality, Arm A.

    Parameters
    ----------
    scores_x, scores_y : (n_bins, k) arrays
        Flat PCA-reduced score matrices for the two areas, indexed by 50 ms bin.
    segments_ : list of :class:`Segment`
        Segments restricting which (t, t+L) pairs are valid.
    cfg : Config
    max_lag_bins : int, optional
        Override ``cfg.lag_ms // cfg.temporal_bin_ms``.
    """
    from . import config as _cfg
    L = (max_lag_bins if max_lag_bins is not None
         else _cfg.lag_bins(cfg))
    lags = np.arange(-L, L + 1)

    # Number of dimensions at lag 0 anchors the result shape.
    x0, _ = _segment_pairs(scores_x, segments_, 0, "x")
    y0, _ = _segment_pairs(scores_y, segments_, 0, "y")
    if x0.shape[0] < 3:
        # Empty / too-fragmented epoch: return all-NaN.
        d = max(scores_x.shape[1], 1)
        return TemporalLagResult(
            lags=lags,
            cc_per_dim=np.full((lags.size, d), np.nan),
            ifi_per_dim=np.full(d, np.nan),
            ifi_windows=np.full((d, max(L, 1)), np.nan),
            peak_lag_per_dim=np.zeros(d, dtype=int),
            n_paired_samples=np.array(
                [segments.total_paired_samples(segments_, int(l)) for l in lags],
                dtype=int),
        )
    d = min(scores_x.shape[1], scores_y.shape[1])
    cc = np.full((lags.size, d), np.nan)
    n_paired = np.zeros(lags.size, dtype=int)
    for i, lag in enumerate(lags):
        x_pairs, trials_x = _segment_pairs(scores_x, segments_, int(lag), "x")
        y_pairs, trials_y = _segment_pairs(scores_y, segments_, int(lag), "y")
        # X / Y pairs share the same trial labels by construction (the same
        # segment contributes the same number of rows on both sides).
        assert np.array_equal(trials_x, trials_y)
        n_paired[i] = x_pairs.shape[0]
        if x_pairs.shape[0] < 6:
            continue
        cc_lag = _held_out_cc_lotot(x_pairs, y_pairs, trials_x, d)
        cc[i, :] = cc_lag

    ifi_per_dim = np.array(
        [lagged.information_flow_index(lags, cc[:, j]) for j in range(d)]
    )
    ifi_windows = np.array(
        [lagged.ifi_by_window(lags, cc[:, j]) for j in range(d)]
    )
    peak = np.array([
        int(lags[np.nanargmax(cc[:, j])]) if np.any(np.isfinite(cc[:, j])) else 0
        for j in range(d)
    ])
    return TemporalLagResult(
        lags=lags,
        cc_per_dim=cc,
        ifi_per_dim=ifi_per_dim,
        ifi_windows=ifi_windows,
        peak_lag_per_dim=peak,
        n_paired_samples=n_paired,
    )


# ---------------------------------------------------------------------------
# Segment-aware circshift surrogate
# ---------------------------------------------------------------------------
def circshift_segments(
    scores: np.ndarray,
    segments_: list[segments.Segment],
    rng: np.random.Generator,
    min_shift: int,
) -> np.ndarray:
    """Circularly shift each segment's rows independently, in place on a copy.

    For each segment of length n, draw a shift in ``[min_shift, n - 1]``
    (clamped to ``max(min_shift, 1)`` when the segment is too short for the
    full range) and ``np.roll`` that segment's slice of ``scores`` along
    axis 0.

    Within-segment shifts are essential: a global circshift across the
    concatenated scores would import the gap structure into the null and
    bias the significance test (the v4-era silent bug; see
    ``legacy/NOTES_temporal_v4.md``).
    """
    out = scores.copy()
    for s in segments_:
        n = len(s)
        if n <= 1:
            continue
        lo = max(int(min_shift), 1)
        hi = max(lo, n - 1)
        shift = int(rng.integers(lo, hi + 1)) if hi > lo else lo
        out[s.start:s.stop] = np.roll(out[s.start:s.stop], shift, axis=0)
    return out


def shuffle_lag_curve(
    scores_x: np.ndarray,
    scores_y: np.ndarray,
    segments_: list[segments.Segment],
    cfg,
    rng: np.random.Generator,
    max_lag_bins: int | None = None,
) -> TemporalLagResult:
    """One surrogate lag scan: per-segment circshift of Y, then refit."""
    shuffled_y = circshift_segments(
        scores_y, segments_, rng,
        min_shift=cfg.temporal_circshift_min_bins,
    )
    return lag_curve(scores_x, shuffled_y, segments_, cfg,
                     max_lag_bins=max_lag_bins)
