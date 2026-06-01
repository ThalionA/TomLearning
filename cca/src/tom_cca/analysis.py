"""Stage 2 analysis: lagged CCA + surrogate significance per (animal, pair).

:func:`analyse_pair` is pure compute over the small per-epoch PC-score tensors
in a :class:`~tom_cca.pipeline.PreparedPair`, so it parallelises cleanly
across processes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import (config, core, lagged, lagged_landmark, lagged_temporal,
               pipeline, surrogate)


@dataclass
class EpochAnalysis:
    """Full Stage-2 result for one (animal, pair, epoch). All CC values are
    held-out (5-fold cross-validated)."""

    epoch: str
    k: int
    n_samples: int
    # Lag-0 CCA.
    held_out_cc: np.ndarray        # (d,) held-out canonical correlations
    in_sample_cc: np.ndarray       # (d,) in-sample CC (diagnostic only)
    cca_full: core.CCAResult       # lag-0 full-data fit (A, B for Stage 3)
    # Lagged CCA / directionality (D6) — held-out, all canonical dimensions.
    lags: np.ndarray               # (2*max_lag+1,) integer bin lags
    lag_cc_per_dim: np.ndarray     # (n_lags, d) held-out CC per lag per dim
    ifi_per_dim: np.ndarray        # (d,) IFI per dimension
    ifi_windows: np.ndarray        # (d, max_window) IFI by lag window (point 4)
    peak_lag_per_dim: np.ndarray   # (d,)
    # Subspace significance (D7) — held-out-CC permutation test, per dimension.
    null_held_out: np.ndarray      # (n_shuffles, d) trial-permutation null
    p_per_dim: np.ndarray          # (d,) per-dimension held-out p-values
    n_significant: int             # significant communication-subspace dims

    # Convenience accessors for the dominant canonical dimension.
    @property
    def lag_cc1(self) -> np.ndarray:
        return self.lag_cc_per_dim[:, 0]

    @property
    def ifi(self) -> float:
        return float(self.ifi_per_dim[0])

    @property
    def peak_lag(self) -> int:
        return int(self.peak_lag_per_dim[0])

    @property
    def p_cc1(self) -> float:
        return float(self.p_per_dim[0])


@dataclass
class PairAnalysis:
    """Stage-2 result for one (animal, area-pair) across all epochs."""

    animal_id: int
    area_x: str
    area_y: str
    role: str
    lp: int
    k: int
    n_units_x: int
    n_units_y: int
    unit_index_x: np.ndarray
    unit_index_y: np.ndarray
    pca_x: dict[str, core.PCAState]      # epoch -> PCA basis (for Stage 3)
    pca_y: dict[str, core.PCAState]
    epochs: dict[str, EpochAnalysis]


def analyse_pair(prepared: pipeline.PreparedPair, cfg=config.DEFAULT) -> PairAnalysis:
    """Run the lagged CCA and surrogate nulls for one prepared area pair."""
    epochs: dict[str, EpochAnalysis] = {}
    for epoch in config.EPOCH_NAMES:
        scores_x = prepared.scores_x[epoch]
        scores_y = prepared.scores_y[epoch]

        # Lag-0 held-out CCA -- unbiased effect size + the full A, B.
        cv0 = core.cca_cv(scores_x, scores_y, cfg)
        # Held-out lagged curve, all canonical dimensions -- directionality.
        lag = lagged.lag_curve(scores_x, scores_y, cfg, held_out=True)
        # Significance: held-out-CC permutation test, per dimension (the
        # in-sample test over-called -- spectrum shift; see surrogate.py).
        null = surrogate.build_null(scores_x, scores_y, cv0.held_out_r, cfg)

        epochs[epoch] = EpochAnalysis(
            epoch=epoch,
            k=cv0.k,
            n_samples=cv0.n_samples,
            held_out_cc=cv0.held_out_r,
            in_sample_cc=cv0.in_sample_r,
            cca_full=cv0.full,
            lags=lag.lags,
            lag_cc_per_dim=lag.cc_per_dim,
            ifi_per_dim=lag.ifi_per_dim,
            ifi_windows=lag.ifi_windows,
            peak_lag_per_dim=lag.peak_lag_per_dim,
            null_held_out=null.null_held_out,
            p_per_dim=null.p_per_dim,
            n_significant=null.n_significant,
        )

    return PairAnalysis(
        animal_id=prepared.animal_id,
        area_x=prepared.area_x,
        area_y=prepared.area_y,
        role=prepared.role,
        lp=prepared.lp,
        k=prepared.k,
        n_units_x=prepared.n_units_x,
        n_units_y=prepared.n_units_y,
        unit_index_x=prepared.unit_index_x,
        unit_index_y=prepared.unit_index_y,
        pca_x=prepared.pca_x,
        pca_y=prepared.pca_y,
        epochs=epochs,
    )


# ===========================================================================
# Arm A -- running-state temporal analysis
# ===========================================================================
@dataclass
class TemporalEpochAnalysis:
    """Arm A Stage-2 result for one (animal, pair, epoch)."""

    epoch: str
    k: int
    n_segments: int
    n_bins_in_segments: int
    n_paired_samples_at_lag0: int
    lags: np.ndarray
    lag_cc_per_dim: np.ndarray
    ifi_per_dim: np.ndarray
    ifi_windows: np.ndarray
    peak_lag_per_dim: np.ndarray
    n_paired_samples: np.ndarray
    # Surrogate (per-segment circshift on Y, full lag scan per shuffle).
    null_ifi: np.ndarray              # (n_shuffles,) per-dim-1 IFI under null
    null_peak_cc_per_dim: np.ndarray  # (n_shuffles, d) peak CC per dim per shuffle
    p_ifi: float                      # P(null IFI >= real IFI), +1/+1 corrected
    p_peak_cc_per_dim: np.ndarray     # (d,) P(null peak CC >= real peak CC)

    @property
    def ifi(self) -> float:
        return float(self.ifi_per_dim[0])

    @property
    def peak_lag(self) -> int:
        return int(self.peak_lag_per_dim[0])


@dataclass
class TemporalPairAnalysis:
    """Arm A Stage-2 result across all epochs."""

    animal_id: int
    area_x: str
    area_y: str
    role: str
    lp: int
    k: int
    n_units_x: int
    n_units_y: int
    unit_index_x: np.ndarray
    unit_index_y: np.ndarray
    epochs: dict[str, TemporalEpochAnalysis]


def _shuffle_lag_temporal(
    scores_x, scores_y, segs, cfg, n_shuffles: int, seed: int,
):
    """Run ``n_shuffles`` per-segment-circshift lag scans on Y."""
    null_ifi = np.full(n_shuffles, np.nan)
    null_peak = None
    for s in range(n_shuffles):
        rng = np.random.default_rng(seed + s + 1)
        nr = lagged_temporal.shuffle_lag_curve(scores_x, scores_y, segs, cfg, rng)
        null_ifi[s] = nr.ifi
        peak = np.nanmax(nr.cc_per_dim, axis=0)
        if null_peak is None:
            null_peak = np.full((n_shuffles, peak.size), np.nan)
        m = min(peak.size, null_peak.shape[1])
        null_peak[s, :m] = peak[:m]
    if null_peak is None:
        null_peak = np.zeros((n_shuffles, 1))
    return null_ifi, null_peak


def _p_value(real: float, null: np.ndarray) -> float:
    null = null[np.isfinite(null)]
    if null.size == 0:
        return float("nan")
    return float((1 + np.sum(null >= real)) / (1 + null.size))


def analyse_pair_temporal(
    prepared: pipeline.PreparedTemporalPair, cfg=config.DEFAULT,
) -> TemporalPairAnalysis:
    """Arm A Stage 2: lag scan + per-segment circshift null per epoch."""
    epochs: dict[str, TemporalEpochAnalysis] = {}
    for epoch in config.EPOCH_NAMES:
        scores_x = prepared.scores_x[epoch]
        scores_y = prepared.scores_y[epoch]
        segs = prepared.segments[epoch]

        real = lagged_temporal.lag_curve(scores_x, scores_y, segs, cfg)
        null_ifi, null_peak = _shuffle_lag_temporal(
            scores_x, scores_y, segs, cfg, cfg.n_shuffles, cfg.surrogate_seed,
        )
        real_peak = np.nanmax(real.cc_per_dim, axis=0)
        d = real.cc_per_dim.shape[1]
        p_peak = np.array([
            _p_value(real_peak[j] if np.isfinite(real_peak[j]) else float("-inf"),
                     null_peak[:, j] if j < null_peak.shape[1]
                     else np.full(cfg.n_shuffles, np.nan))
            for j in range(d)
        ])
        epochs[epoch] = TemporalEpochAnalysis(
            epoch=epoch,
            k=prepared.k,
            n_segments=len(segs),
            n_bins_in_segments=sum(len(s) for s in segs),
            n_paired_samples_at_lag0=int(real.n_paired_samples[
                len(real.lags) // 2]),
            lags=real.lags,
            lag_cc_per_dim=real.cc_per_dim,
            ifi_per_dim=real.ifi_per_dim,
            ifi_windows=real.ifi_windows,
            peak_lag_per_dim=real.peak_lag_per_dim,
            n_paired_samples=real.n_paired_samples,
            null_ifi=null_ifi,
            null_peak_cc_per_dim=null_peak,
            p_ifi=_p_value(real.ifi, null_ifi),
            p_peak_cc_per_dim=p_peak,
        )
    return TemporalPairAnalysis(
        animal_id=prepared.animal_id,
        area_x=prepared.area_x, area_y=prepared.area_y,
        role=prepared.role, lp=prepared.lp, k=prepared.k,
        n_units_x=prepared.n_units_x, n_units_y=prepared.n_units_y,
        unit_index_x=prepared.unit_index_x,
        unit_index_y=prepared.unit_index_y,
        epochs=epochs,
    )


# ===========================================================================
# Arm B -- landmark-centred analysis
# ===========================================================================
@dataclass
class LandmarkCellAnalysis:
    """Arm B Stage-2 result for one (animal, pair, epoch, landmark)."""

    epoch: str
    landmark_id: int
    n_crossings: int
    skip_reason: str | None
    # Populated only when skip_reason is None:
    k: int = 0
    lags: np.ndarray | None = None
    lag_cc_per_dim: np.ndarray | None = None
    ifi_per_dim: np.ndarray | None = None
    ifi_windows: np.ndarray | None = None
    peak_lag_per_dim: np.ndarray | None = None
    n_paired_samples: np.ndarray | None = None
    null_ifi: np.ndarray | None = None
    null_peak_cc_per_dim: np.ndarray | None = None
    p_ifi: float = float("nan")
    p_peak_cc_per_dim: np.ndarray | None = None

    @property
    def ifi(self) -> float:
        return (float("nan") if self.ifi_per_dim is None
                else float(self.ifi_per_dim[0]))

    @property
    def peak_lag(self) -> int:
        return (0 if self.peak_lag_per_dim is None
                else int(self.peak_lag_per_dim[0]))


@dataclass
class LandmarkPairAnalysis:
    """Arm B Stage-2 result across all (epoch, landmark) cells."""

    animal_id: int
    area_x: str
    area_y: str
    role: str
    lp: int
    n_units_x: int
    n_units_y: int
    unit_index_x: np.ndarray
    unit_index_y: np.ndarray
    cells: dict[tuple[str, int], LandmarkCellAnalysis]
    coverage: dict[tuple[str, int], int]


def _shuffle_lag_landmark(
    scores_x, scores_y, trials, landmark_id, cfg, n_shuffles: int, seed: int,
):
    null_ifi = np.full(n_shuffles, np.nan)
    null_peak = None
    for s in range(n_shuffles):
        rng = np.random.default_rng(seed + s + 1)
        nr = lagged_landmark.shuffle_lag_curve(
            scores_x, scores_y, trials, landmark_id, cfg, rng,
        )
        null_ifi[s] = nr.ifi
        peak = np.nanmax(nr.cc_per_dim, axis=0)
        if null_peak is None:
            null_peak = np.full((n_shuffles, peak.size), np.nan)
        m = min(peak.size, null_peak.shape[1])
        null_peak[s, :m] = peak[:m]
    if null_peak is None:
        null_peak = np.zeros((n_shuffles, 1))
    return null_ifi, null_peak


def analyse_pair_landmark(
    prepared: pipeline.PreparedLandmarkPair, cfg=config.DEFAULT,
) -> LandmarkPairAnalysis:
    """Arm B Stage 2: per-(epoch, landmark) lag scan + per-window circshift null."""
    cells: dict[tuple[str, int], LandmarkCellAnalysis] = {}
    for key, prep_cell in prepared.cells.items():
        if prep_cell.skip_reason is not None:
            cells[key] = LandmarkCellAnalysis(
                epoch=prep_cell.epoch,
                landmark_id=prep_cell.landmark_id,
                n_crossings=prep_cell.n_crossings,
                skip_reason=prep_cell.skip_reason,
            )
            continue
        real = lagged_landmark.lag_curve(
            prep_cell.scores_x, prep_cell.scores_y,
            prep_cell.trials, prep_cell.landmark_id, cfg,
        )
        null_ifi, null_peak = _shuffle_lag_landmark(
            prep_cell.scores_x, prep_cell.scores_y, prep_cell.trials,
            prep_cell.landmark_id, cfg, cfg.n_shuffles, cfg.surrogate_seed,
        )
        real_peak = np.nanmax(real.cc_per_dim, axis=0)
        d = real.cc_per_dim.shape[1]
        p_peak = np.array([
            _p_value(real_peak[j] if np.isfinite(real_peak[j]) else float("-inf"),
                     null_peak[:, j] if j < null_peak.shape[1]
                     else np.full(cfg.n_shuffles, np.nan))
            for j in range(d)
        ])
        cells[key] = LandmarkCellAnalysis(
            epoch=prep_cell.epoch,
            landmark_id=prep_cell.landmark_id,
            n_crossings=prep_cell.n_crossings,
            skip_reason=None,
            k=prep_cell.k,
            lags=real.lags,
            lag_cc_per_dim=real.cc_per_dim,
            ifi_per_dim=real.ifi_per_dim,
            ifi_windows=real.ifi_windows,
            peak_lag_per_dim=real.peak_lag_per_dim,
            n_paired_samples=real.n_paired_samples,
            null_ifi=null_ifi,
            null_peak_cc_per_dim=null_peak,
            p_ifi=_p_value(real.ifi, null_ifi),
            p_peak_cc_per_dim=p_peak,
        )
    return LandmarkPairAnalysis(
        animal_id=prepared.animal_id,
        area_x=prepared.area_x, area_y=prepared.area_y,
        role=prepared.role, lp=prepared.lp,
        n_units_x=len(prepared.unit_index_x),
        n_units_y=len(prepared.unit_index_y),
        unit_index_x=prepared.unit_index_x,
        unit_index_y=prepared.unit_index_y,
        cells=cells,
        coverage=prepared.coverage,
    )
