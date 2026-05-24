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

from . import config, core, dataio, partial


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
