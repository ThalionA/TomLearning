"""Window-bounded lagged CCA for Arm B (landmark-centred temporal).

Refits CCA at every integer-bin lag in ``[-landmark_lag_bins, +landmark_lag_bins]``,
where the per-lag design matrix is built by pooling within-window (X-bin t,
Y-bin t+L) pairs across all kept crossings of one landmark. Window length
caps |L| at ``window_bins - 5`` (the minimum-paired-sample floor); the lag
range is bounded by the window, not the other way round.

Each call operates on a single per-landmark fit cell: PC scores per area
shaped ``(n_kept_crossings, n_window_bins, k)``. The per-(bin, unit) PSTH
residualisation (if requested) happens upstream in the pipeline before PCA.

Sign convention: positive lag => X leads Y. See
``cca/UNDERSTANDING_temporal.md`` Sec. 6.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import core, lagged, lagged_temporal


@dataclass
class LandmarkLagResult:
    """Lagged-CCA directionality for one (animal, pair, epoch, landmark)."""

    landmark_id: int
    lags: np.ndarray              # (n_lags,) integer bin lags
    cc_per_dim: np.ndarray        # (n_lags, n_dims)
    ifi_per_dim: np.ndarray       # (n_dims,)
    ifi_windows: np.ndarray       # (n_dims, max_window)
    peak_lag_per_dim: np.ndarray  # (n_dims,)
    n_paired_samples: np.ndarray  # (n_lags,) pooled (crossings * (win - |L|))
    n_crossings: int

    @property
    def cc1(self) -> np.ndarray:
        return self.cc_per_dim[:, 0]

    @property
    def ifi(self) -> float:
        return float(self.ifi_per_dim[0])

    @property
    def peak_lag(self) -> int:
        return int(self.peak_lag_per_dim[0])


def _window_pairs(
    scores: np.ndarray,            # (n_crossings, win_len, k)
    trials: np.ndarray,            # (n_crossings,)
    lag: int,
    which: str,                    # "x" or "y"
) -> tuple[np.ndarray, np.ndarray]:
    """Slice within each window for lag L, pool across crossings.

    Returns ``(pooled, trial_labels)`` where the trial label is broadcast to
    every row contributed by that crossing.
    """
    n_cr, win_len, k = scores.shape
    L = int(lag)
    if abs(L) >= win_len:
        return np.zeros((0, k)), np.zeros(0, dtype=int)
    if L >= 0:
        if which == "x":
            sliced = scores[:, : win_len - L, :]
        else:
            sliced = scores[:, L:, :]
    else:
        if which == "x":
            sliced = scores[:, -L:, :]
        else:
            sliced = scores[:, : win_len + L, :]
    pooled = sliced.reshape(n_cr * sliced.shape[1], k)
    trial_labels = np.repeat(trials, sliced.shape[1])
    return pooled, trial_labels


def lag_curve(
    scores_x: np.ndarray,
    scores_y: np.ndarray,
    trials: np.ndarray,
    landmark_id: int,
    cfg,
    max_lag_bins: int | None = None,
) -> LandmarkLagResult:
    """Refit CCA at every lag, pooling within-window pairs across crossings.

    Parameters
    ----------
    scores_x, scores_y : (n_crossings, win_len, k) arrays
    trials : (n_crossings,) int -- the original trial label per crossing
    landmark_id : int
    cfg : Config
    max_lag_bins : int, optional
        Override ``cfg.landmark_lag_ms // cfg.temporal_bin_ms``.
    """
    from . import config as _cfg
    n_cr, win_len, _ = scores_x.shape
    if scores_y.shape != scores_x.shape:
        raise ValueError("scores_x and scores_y must share shape")
    if trials.shape != (n_cr,):
        raise ValueError("trials shape mismatch")
    requested = (max_lag_bins if max_lag_bins is not None
                 else _cfg.landmark_lag_bins(cfg))
    # Bound the lag by what the window allows (need >= 5 paired samples).
    max_allowed = max(0, win_len - 5)
    L = int(min(requested, max_allowed))
    lags = np.arange(-L, L + 1)
    d = min(scores_x.shape[2], scores_y.shape[2])
    cc = np.full((lags.size, d), np.nan)
    n_paired = np.zeros(lags.size, dtype=int)

    for i, lag in enumerate(lags):
        x_pairs, trials_x = _window_pairs(scores_x, trials, int(lag), "x")
        y_pairs, trials_y = _window_pairs(scores_y, trials, int(lag), "y")
        assert np.array_equal(trials_x, trials_y)
        n_paired[i] = x_pairs.shape[0]
        if x_pairs.shape[0] < 6:
            continue
        cc[i, :] = lagged_temporal._held_out_cc_lotot(
            x_pairs, y_pairs, trials_x, d
        )

    ifi_per_dim = np.array(
        [lagged.information_flow_index(lags, cc[:, j]) for j in range(d)]
    )
    ifi_windows = np.array(
        [lagged.ifi_by_window(lags, cc[:, j]) for j in range(d)]
    ) if L > 0 else np.full((d, 1), np.nan)
    peak = np.array([
        int(lags[np.nanargmax(cc[:, j])]) if np.any(np.isfinite(cc[:, j])) else 0
        for j in range(d)
    ])
    return LandmarkLagResult(
        landmark_id=int(landmark_id),
        lags=lags,
        cc_per_dim=cc,
        ifi_per_dim=ifi_per_dim,
        ifi_windows=ifi_windows,
        peak_lag_per_dim=peak,
        n_paired_samples=n_paired,
        n_crossings=int(n_cr),
    )


# ---------------------------------------------------------------------------
# Per-window circshift surrogate
# ---------------------------------------------------------------------------
def circshift_windows(
    scores: np.ndarray,         # (n_crossings, win_len, k)
    rng: np.random.Generator,
    min_shift: int,
) -> np.ndarray:
    """Circularly shift each window's bin axis independently.

    For each crossing's window, draw a shift in ``[min_shift, win_len - 1]``
    and roll that window along the bin axis. Windows are short (10 bins) but
    pooling many windows per surrogate is what gives the null its statistical
    power -- the small per-window range is fine.
    """
    n_cr, win_len, _ = scores.shape
    out = scores.copy()
    lo = max(int(min_shift), 1)
    hi = max(lo, win_len - 1)
    if hi <= lo:
        # Window too short to shuffle meaningfully -- return unchanged.
        return out
    for i in range(n_cr):
        s = int(rng.integers(lo, hi + 1))
        out[i] = np.roll(out[i], s, axis=0)
    return out


def shuffle_lag_curve(
    scores_x: np.ndarray,
    scores_y: np.ndarray,
    trials: np.ndarray,
    landmark_id: int,
    cfg,
    rng: np.random.Generator,
    max_lag_bins: int | None = None,
) -> LandmarkLagResult:
    """One surrogate lag scan: per-window circshift of Y, then refit."""
    shuffled_y = circshift_windows(
        scores_y, rng,
        min_shift=cfg.landmark_circshift_min_bins,
    )
    return lag_curve(scores_x, shuffled_y, trials, landmark_id, cfg,
                     max_lag_bins=max_lag_bins)
