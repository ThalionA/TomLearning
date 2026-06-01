"""Per-(animal, pair) orchestration of the CCA pipeline -- Tom edition.

Two-step structure so the data-heavy preparation can be separated from the
parallelisable compute:

* :func:`prepare_pair` -- load, FS-select, residualise, per-epoch PCA. Produces
  small per-epoch PC-score tensors. Runs in the main process.
* :func:`fit_pair` (Stage 1) and :func:`tom_cca.analysis.analyse_pair`
  (Stage 2) -- pure compute on those small score tensors.

Spatial-only (no temporal binning), per the user's instruction.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from . import config, core, dataio, partial, segments as segments_mod


@dataclass
class PreparedPair:
    """Per-epoch PCA-reduced residual scores for one (animal, area-pair)."""

    animal_id: int
    area_x: str
    area_y: str
    role: str                              # "learner" | "nonlearner"
    lp: int
    k: int                                 # PCs per area, fixed over epochs
    n_units_x: int
    n_units_y: int
    unit_index_x: np.ndarray
    unit_index_y: np.ndarray
    scores_x: dict[str, np.ndarray]        # epoch -> (n_trials, n_bins, k)
    scores_y: dict[str, np.ndarray]
    pca_x: dict[str, core.PCAState]        # epoch -> PCA basis for X
    pca_y: dict[str, core.PCAState]


@dataclass
class SkippedPair:
    """A (animal, pair) that could not be prepared, with the reason."""

    animal_id: int
    area_x: str
    area_y: str
    reason: str


@dataclass
class PairFit:
    """Stage-1 cross-validated CCA for one (animal, area-pair) over all epochs."""

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
    pca_x_by_epoch: dict[str, core.PCAState]
    pca_y_by_epoch: dict[str, core.PCAState]
    epochs: dict[str, core.CVResult]


def config_label(cfg) -> str:
    """One-line human description of an analysis config."""
    cca = "residual" if cfg.subtract_trial_mean else "signal"
    fs = "FS-excl" if cfg.exclude_fast_spiking else "FS-incl"
    z = "z-on" if cfg.zscore_units else "z-off"
    return f"spatial, {cca}, {fs}, {z}, k={cfg.k_mode}"


def prepare_pair(
    animal: dataio.Animal,
    area_x: str,
    area_y: str,
    entry: dataio.CohortEntry,
    cfg=config.DEFAULT,
) -> PreparedPair | SkippedPair:
    """Build per-epoch PCA-reduced residual scores for one area pair.

    Steps: area tensors (FS units excluded) -> whole-engaged-period z-scoring
    -> per-epoch residualisation -> k from the sample budget, capped at the
    smallest per-epoch numerical rank -> per-epoch PCA.
    """
    tensor_x, idx_x = dataio.area_tensor(animal, area_x, cfg)
    tensor_y, idx_y = dataio.area_tensor(animal, area_y, cfg)
    if cfg.zscore_units:
        tensor_x = core.zscore_units(tensor_x)
        tensor_y = core.zscore_units(tensor_y)

    n_units_x = len(idx_x)
    n_units_y = len(idx_y)
    if n_units_x < cfg.min_units or n_units_y < cfg.min_units:
        return SkippedPair(
            animal.animal_id, area_x, area_y,
            f"too few units ({area_x}={n_units_x}, {area_y}={n_units_y}; "
            f"min {cfg.min_units})",
        )

    n_use = tensor_x.shape[0]
    windows = dataio.epoch_windows(entry.lp, n_use, cfg)
    if windows is None:
        return SkippedPair(
            animal.animal_id, area_x, area_y,
            f"no valid epochs (lp={entry.lp}, usable trials={n_use})",
        )

    res_x = {e: _residual(tensor_x[idx], cfg) for e, idx in windows.items()}
    res_y = {e: _residual(tensor_y[idx], cfg) for e, idx in windows.items()}

    n_valid = min(core.n_valid_samples(res_x[e]) for e in config.EPOCH_NAMES)
    min_rank = min(
        min(core.numerical_rank(res_x[e]) for e in config.EPOCH_NAMES),
        min(core.numerical_rank(res_y[e]) for e in config.EPOCH_NAMES),
    )
    k = core.choose_k(n_units_x, n_units_y, n_valid, cfg, max_rank=min_rank,
                      variance_k=_variance_k([res_x, res_y], cfg))

    scores_x, scores_y = {}, {}
    pca_x, pca_y = {}, {}
    for epoch in config.EPOCH_NAMES:
        px = core.pca_fit(res_x[epoch], k)
        py = core.pca_fit(res_y[epoch], k)
        pca_x[epoch] = px
        pca_y[epoch] = py
        scores_x[epoch] = core.pca_transform(res_x[epoch], px)
        scores_y[epoch] = core.pca_transform(res_y[epoch], py)

    return PreparedPair(
        animal_id=animal.animal_id,
        area_x=area_x, area_y=area_y,
        role=entry.role, lp=entry.lp, k=k,
        n_units_x=n_units_x, n_units_y=n_units_y,
        unit_index_x=idx_x, unit_index_y=idx_y,
        scores_x=scores_x, scores_y=scores_y,
        pca_x=pca_x, pca_y=pca_y,
    )


def _residual_tensors(
    animal: dataio.Animal, area: str, entry: dataio.CohortEntry, cfg
) -> tuple[dict[str, np.ndarray], np.ndarray] | None:
    """Per-epoch residualised (z-scored) neuron tensors for one area, or None."""
    tensor, idx = dataio.area_tensor(animal, area, cfg)
    if len(idx) < cfg.min_units:
        return None
    if cfg.zscore_units:
        tensor = core.zscore_units(tensor)
    windows = dataio.epoch_windows(entry.lp, tensor.shape[0], cfg)
    if windows is None:
        return None
    res = {e: _residual(tensor[e_idx], cfg) for e, e_idx in windows.items()}
    return res, idx


def prepare_pair_partial(
    animal: dataio.Animal,
    area_x: str,
    area_y: str,
    entry: dataio.CohortEntry,
    cfg=config.DEFAULT,
) -> PreparedPair | SkippedPair:
    """Like :func:`prepare_pair`, but with every other recorded area removed.

    Each other area's PC scores are regressed out of X's and Y's residualised
    neuron tensors *before* the per-epoch PCA -- neuron-level partialling so
    the PCA basis and Stage-3 back-projection stay in X's / Y's own neuron
    space. Returns a SkippedPair if X or Y is unusable or the animal has no
    other area to condition on.
    """
    rx = _residual_tensors(animal, area_x, entry, cfg)
    ry = _residual_tensors(animal, area_y, entry, cfg)
    if rx is None or ry is None:
        return SkippedPair(animal.animal_id, area_x, area_y,
                           "X or Y unusable for partial preparation")
    res_x, idx_x = rx
    res_y, idx_y = ry

    confounds = []
    for area_z in config.AREAS:
        if area_z in (area_x, area_y):
            continue
        sz = prepare_area(animal, area_z, entry, cfg)
        if sz is not None:
            confounds.append(sz)
    if not confounds:
        return SkippedPair(animal.animal_id, area_x, area_y,
                           "no other area to partial out")

    res_xp, res_yp = {}, {}
    for epoch in config.EPOCH_NAMES:
        z = np.concatenate([c[epoch] for c in confounds], axis=-1)
        res_xp[epoch] = partial.partial_out_tensor(res_x[epoch], z)
        res_yp[epoch] = partial.partial_out_tensor(res_y[epoch], z)

    n_units_x, n_units_y = len(idx_x), len(idx_y)
    n_valid = min(core.n_valid_samples(res_xp[e]) for e in config.EPOCH_NAMES)
    min_rank = min(
        min(core.numerical_rank(res_xp[e]) for e in config.EPOCH_NAMES),
        min(core.numerical_rank(res_yp[e]) for e in config.EPOCH_NAMES),
    )
    k = core.choose_k(n_units_x, n_units_y, n_valid, cfg, max_rank=min_rank,
                      variance_k=_variance_k([res_xp, res_yp], cfg))

    scores_x, scores_y = {}, {}
    pca_x, pca_y = {}, {}
    for epoch in config.EPOCH_NAMES:
        px = core.pca_fit(res_xp[epoch], k)
        py = core.pca_fit(res_yp[epoch], k)
        pca_x[epoch] = px
        pca_y[epoch] = py
        scores_x[epoch] = core.pca_transform(res_xp[epoch], px)
        scores_y[epoch] = core.pca_transform(res_yp[epoch], py)

    return PreparedPair(
        animal_id=animal.animal_id,
        area_x=area_x, area_y=area_y,
        role=entry.role, lp=entry.lp, k=k,
        n_units_x=n_units_x, n_units_y=n_units_y,
        unit_index_x=idx_x, unit_index_y=idx_y,
        scores_x=scores_x, scores_y=scores_y,
        pca_x=pca_x, pca_y=pca_y,
    )


def fit_pair(
    animal: dataio.Animal,
    area_x: str,
    area_y: str,
    entry: dataio.CohortEntry,
    cfg=config.DEFAULT,
) -> PairFit | SkippedPair:
    """Stage-1 cross-validated CCA for one area pair (lag 0, no surrogates)."""
    prepared = prepare_pair(animal, area_x, area_y, entry, cfg)
    if isinstance(prepared, SkippedPair):
        return prepared
    epochs = {
        e: core.cca_cv(prepared.scores_x[e], prepared.scores_y[e], cfg)
        for e in config.EPOCH_NAMES
    }
    return PairFit(
        animal_id=prepared.animal_id,
        area_x=area_x, area_y=area_y,
        role=prepared.role, lp=prepared.lp, k=prepared.k,
        n_units_x=prepared.n_units_x, n_units_y=prepared.n_units_y,
        unit_index_x=prepared.unit_index_x,
        unit_index_y=prepared.unit_index_y,
        pca_x_by_epoch=prepared.pca_x,
        pca_y_by_epoch=prepared.pca_y,
        epochs=epochs,
    )


def prepare_area(
    animal: dataio.Animal, area: str, entry: dataio.CohortEntry, cfg=config.DEFAULT
) -> dict[str, np.ndarray] | None:
    """Per-epoch residual PC scores for a single area (used by partial CCA)."""
    tensor, idx = dataio.area_tensor(animal, area, cfg)
    if len(idx) < cfg.min_units:
        return None
    if cfg.zscore_units:
        tensor = core.zscore_units(tensor)
    windows = dataio.epoch_windows(entry.lp, tensor.shape[0], cfg)
    if windows is None:
        return None
    res = {e: _residual(tensor[e_idx], cfg) for e, e_idx in windows.items()}
    n_units = len(idx)
    n_valid = min(core.n_valid_samples(res[e]) for e in config.EPOCH_NAMES)
    min_rank = min(core.numerical_rank(res[e]) for e in config.EPOCH_NAMES)
    k = core.choose_k(n_units, n_units, n_valid, cfg, max_rank=min_rank,
                      variance_k=_variance_k([res], cfg))
    return {
        e: core.pca_transform(res[e], core.pca_fit(res[e], k))
        for e in config.EPOCH_NAMES
    }


def _residual(tensor: np.ndarray, cfg) -> np.ndarray:
    """Per-epoch centring: residualise or, for the signal variant, remove the
    per-unit grand mean. Units are z-scored upstream so no per-unit scaling
    happens here. Missing samples stay NaN and are dropped later, at fit time.
    """
    if cfg.subtract_trial_mean:
        return core.residualise(tensor)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        grand = np.nanmean(tensor, axis=(0, 1), keepdims=True)
    return tensor - grand


# ===========================================================================
# Arm A -- running-state temporal (per (animal, pair, epoch))
# ===========================================================================
@dataclass
class PreparedTemporalPair:
    """Arm A: per-epoch PCA-reduced score matrices + segment lists.

    Scores are flat ``(n_bins_50ms, k)`` matrices indexed by the global 50 ms
    timeline; the ``segments`` field restricts which bins enter the lag scan.
    """

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
    scores_x: dict[str, np.ndarray]                  # epoch -> (n_bins_50ms, k)
    scores_y: dict[str, np.ndarray]
    pca_x: dict[str, core.PCAState]
    pca_y: dict[str, core.PCAState]
    segments: dict[str, list[segments_mod.Segment]]  # epoch -> segments in epoch


def _epoch_rows(activity: np.ndarray, segs: list[segments_mod.Segment]) -> np.ndarray:
    """Vertical stack of the 50 ms rows inside the given segments."""
    if not segs:
        return np.zeros((0, activity.shape[1]), dtype=activity.dtype)
    return np.concatenate([activity[s.start:s.stop] for s in segs], axis=0)


def prepare_pair_temporal(
    animal: dataio.Animal,
    area_x: str,
    area_y: str,
    entry: dataio.CohortEntry,
    cfg=config.DEFAULT,
) -> PreparedTemporalPair | SkippedPair:
    """Arm A: build per-epoch score matrices for one area pair.

    Steps: per-area 50 ms activity (z-scored on running+cued reference) ->
    segment construction (velocity-passing, trial-bounded) -> per-epoch
    segment filter -> k from min epoch sample budget -> PCA per area per
    epoch -> project full 50 ms timeline to scores.

    Signal-only (ragged segment lengths preclude per-bin PSTH).
    """
    act_x, idx_x = dataio.area_activity_50ms(animal, area_x, cfg)
    act_y, idx_y = dataio.area_activity_50ms(animal, area_y, cfg)
    n_x, n_y = len(idx_x), len(idx_y)
    if n_x < cfg.min_units or n_y < cfg.min_units:
        return SkippedPair(
            animal.animal_id, area_x, area_y,
            f"too few units ({area_x}={n_x}, {area_y}={n_y}; min {cfg.min_units})",
        )

    all_segs = dataio.temporal_segments(animal, cfg)
    if not all_segs:
        return SkippedPair(animal.animal_id, area_x, area_y,
                           "no segments after velocity filter + morph cleanup")

    windows = dataio.epoch_windows(entry.lp, animal.n_trials, cfg)
    if windows is None:
        return SkippedPair(animal.animal_id, area_x, area_y,
                           f"no valid epochs (lp={entry.lp}, "
                           f"usable trials={animal.n_trials})")

    per_epoch_segs: dict[str, list[segments_mod.Segment]] = {}
    per_epoch_samples: dict[str, int] = {}
    for epoch, trial_arr in windows.items():
        trial_set = set(int(t) for t in trial_arr.tolist())
        ep = [s for s in all_segs if s.trial in trial_set]
        per_epoch_segs[epoch] = ep
        per_epoch_samples[epoch] = sum(len(s) for s in ep)

    # Need enough segment bins in every epoch to support PCA + a lag scan.
    min_samples_floor = max(20, 3 * config.lag_bins(cfg))
    if any(per_epoch_samples[e] < min_samples_floor for e in config.EPOCH_NAMES):
        return SkippedPair(
            animal.animal_id, area_x, area_y,
            f"too few segment bins per epoch: {per_epoch_samples} "
            f"(floor {min_samples_floor})",
        )

    n_valid = min(per_epoch_samples.values())
    min_rank = min(
        min(core.numerical_rank(_epoch_rows(act_x, per_epoch_segs[e]))
            for e in config.EPOCH_NAMES),
        min(core.numerical_rank(_epoch_rows(act_y, per_epoch_segs[e]))
            for e in config.EPOCH_NAMES),
    )
    k = core.choose_k(n_x, n_y, n_valid, cfg, max_rank=min_rank, variance_k=None)

    scores_x, scores_y, pca_x, pca_y = {}, {}, {}, {}
    for epoch in config.EPOCH_NAMES:
        ep_x = _epoch_rows(act_x, per_epoch_segs[epoch])
        ep_y = _epoch_rows(act_y, per_epoch_segs[epoch])
        # core.pca_fit wraps the flatten/finite-row handling; give it a
        # (1, n_rows, n_units) "trial of one" tensor.
        px = core.pca_fit(ep_x[None, :, :], k)
        py = core.pca_fit(ep_y[None, :, :], k)
        pca_x[epoch] = px
        pca_y[epoch] = py
        # Project the full 50 ms timeline so the lag scan can index with the
        # global segment starts/stops without per-epoch coordinate gymnastics.
        scores_x[epoch] = (act_x - px.mean) @ px.components
        scores_y[epoch] = (act_y - py.mean) @ py.components

    return PreparedTemporalPair(
        animal_id=animal.animal_id, area_x=area_x, area_y=area_y,
        role=entry.role, lp=entry.lp, k=k,
        n_units_x=n_x, n_units_y=n_y,
        unit_index_x=idx_x, unit_index_y=idx_y,
        scores_x=scores_x, scores_y=scores_y,
        pca_x=pca_x, pca_y=pca_y,
        segments=per_epoch_segs,
    )


# ===========================================================================
# Arm B -- landmark-centred (per (animal, pair, epoch, landmark))
# ===========================================================================
@dataclass
class PreparedLandmarkCell:
    """One (epoch, landmark) fit cell for Arm B."""

    epoch: str
    landmark_id: int
    n_crossings: int
    skip_reason: str | None
    # Only populated when skip_reason is None:
    k: int = 0
    n_units_x: int = 0
    n_units_y: int = 0
    scores_x: np.ndarray | None = None   # (n_crossings, win_len, k)
    scores_y: np.ndarray | None = None
    trials: np.ndarray | None = None     # (n_crossings,) int
    pca_x: core.PCAState | None = None
    pca_y: core.PCAState | None = None


@dataclass
class PreparedLandmarkPair:
    """Arm B: per-(epoch, landmark) prepared cells."""

    animal_id: int
    area_x: str
    area_y: str
    role: str
    lp: int
    unit_index_x: np.ndarray
    unit_index_y: np.ndarray
    cells: dict[tuple[str, int], PreparedLandmarkCell]
    # Coverage summary -- n_kept_crossings per (epoch, landmark) and skip-reason.
    coverage: dict[tuple[str, int], int]


def _residualise_windows(tensor: np.ndarray, subtract_trial_mean: bool) -> np.ndarray:
    """Centre a (n_crossings, win_len, n_units) tensor.

    With ``subtract_trial_mean=True`` (the residual-CCA default) the per-
    ``(bin, unit)`` mean across crossings is subtracted -- removes the
    landmark-evoked response. Otherwise just the per-unit grand mean is
    removed (signal CCA).
    """
    if subtract_trial_mean:
        return tensor - tensor.mean(axis=0, keepdims=True)
    return tensor - tensor.mean(axis=(0, 1), keepdims=True)


def prepare_pair_landmark(
    animal: dataio.Animal,
    area_x: str,
    area_y: str,
    entry: dataio.CohortEntry,
    cfg=config.DEFAULT,
) -> PreparedLandmarkPair | SkippedPair:
    """Arm B: build per-(epoch, landmark) score tensors for one area pair.

    Steps: per-area landmark windows (already gated by per-window velocity)
    -> epoch filter on each window's trial -> per-landmark per-epoch
    residualisation + PCA -> score tensor (n_crossings, win_len, k).

    Cells with fewer than ``cfg.min_crossings_per_landmark`` kept crossings
    are emitted with ``skip_reason="low_crossings"`` instead of a fit.
    """
    win_x, idx_x = dataio.area_landmark_windows(animal, area_x, cfg)
    win_y, idx_y = dataio.area_landmark_windows(animal, area_y, cfg)
    n_x, n_y = len(idx_x), len(idx_y)
    if n_x < cfg.min_units or n_y < cfg.min_units:
        return SkippedPair(
            animal.animal_id, area_x, area_y,
            f"too few units ({area_x}={n_x}, {area_y}={n_y}; min {cfg.min_units})",
        )

    windows = dataio.epoch_windows(entry.lp, animal.n_trials, cfg)
    if windows is None:
        return SkippedPair(animal.animal_id, area_x, area_y,
                           f"no valid epochs (lp={entry.lp})")

    landmark_ids = sorted(set(win_x.keys()) & set(win_y.keys()))
    if not landmark_ids:
        return SkippedPair(animal.animal_id, area_x, area_y,
                           "no landmark passes the velocity gate")

    cells: dict[tuple[str, int], PreparedLandmarkCell] = {}
    coverage: dict[tuple[str, int], int] = {}
    win_len = config.landmark_window_bins(cfg)

    for lm in landmark_ids:
        x_lm = win_x[lm]
        y_lm = win_y[lm]
        if not np.array_equal(x_lm.trials, y_lm.trials):
            # Different gating would produce different trial labels -- a bug
            # somewhere upstream. Surface it loudly.
            raise RuntimeError(
                f"landmark {lm}: X and Y disagree on trial labels for kept crossings"
            )
        for epoch, trial_arr in windows.items():
            trial_set = set(int(t) for t in trial_arr.tolist())
            mask = np.array([int(t) in trial_set for t in x_lm.trials])
            n_kept = int(np.sum(mask))
            coverage[(epoch, lm)] = n_kept
            if n_kept < cfg.min_crossings_per_landmark:
                cells[(epoch, lm)] = PreparedLandmarkCell(
                    epoch=epoch, landmark_id=lm, n_crossings=n_kept,
                    skip_reason="low_crossings",
                )
                continue
            X_raw = x_lm.spikes[mask]                  # (n_kept, win_len, n_x)
            Y_raw = y_lm.spikes[mask]
            X_res = _residualise_windows(X_raw, cfg.subtract_trial_mean)
            Y_res = _residualise_windows(Y_raw, cfg.subtract_trial_mean)
            n_samples = n_kept * win_len
            min_rank = min(core.numerical_rank(X_res), core.numerical_rank(Y_res))
            k = core.choose_k(n_x, n_y, n_samples, cfg, max_rank=min_rank)
            if k < 1:
                cells[(epoch, lm)] = PreparedLandmarkCell(
                    epoch=epoch, landmark_id=lm, n_crossings=n_kept,
                    skip_reason="rank_deficient",
                )
                continue
            try:
                px = core.pca_fit(X_res, k)
                py = core.pca_fit(Y_res, k)
            except Exception as e:                     # pragma: no cover
                cells[(epoch, lm)] = PreparedLandmarkCell(
                    epoch=epoch, landmark_id=lm, n_crossings=n_kept,
                    skip_reason=f"pca_failed: {type(e).__name__}",
                )
                continue
            sx = core.pca_transform(X_res, px)        # (n_kept, win_len, k)
            sy = core.pca_transform(Y_res, py)
            cells[(epoch, lm)] = PreparedLandmarkCell(
                epoch=epoch, landmark_id=lm, n_crossings=n_kept,
                skip_reason=None,
                k=k, n_units_x=n_x, n_units_y=n_y,
                scores_x=sx, scores_y=sy,
                trials=x_lm.trials[mask].astype(int),
                pca_x=px, pca_y=py,
            )

    return PreparedLandmarkPair(
        animal_id=animal.animal_id, area_x=area_x, area_y=area_y,
        role=entry.role, lp=entry.lp,
        unit_index_x=idx_x, unit_index_y=idx_y,
        cells=cells, coverage=coverage,
    )


def _variance_k(area_residuals, cfg) -> int | None:
    """PCs reaching ``cfg.k_variance`` cumulative variance, pooled over epochs
    and taken as the symmetric (min) value across the supplied areas. Returns
    None unless ``cfg.k_mode == "variance"``.
    """
    if cfg.k_mode != "variance":
        return None
    ks = []
    for res in area_residuals:
        pooled = np.concatenate(
            [r.reshape(-1, r.shape[-1]) for r in res.values()])
        ks.append(core.k_for_variance(pooled, cfg.k_variance))
    return min(ks)
