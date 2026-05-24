"""Data IO and cohort logic for the Tom-learning CCA pipeline.

Loads Tom's per-animal exports (``TF*_export.mat``, MATLAB v7.3 / HDF5) and the
cohort learning-point file (``animal_behaviour.mat``), classifies animals into
learners and yoked non-learners, and builds per-area activity tensors with
fast-spiking units excluded (Tom convention: FS is only flagged in
V1/RSC/CA1/CA3 -- see :mod:`tom_cca.config`).

Conventions
-----------
The MATLAB pipeline stores the spatial firing rate as
``analysis_spatial.firing.cued.freq`` of shape ``(n_units, n_trials, n_bins)``
(see HC_V1_Code/legacy/CCA_HC_V1_spatial_v2.m). h5py reads MATLAB v7.3 datasets
with dimensions reversed, so the on-disk read yields ``(n_bins, n_trials,
n_units)``; the loader transposes that to ``(n_trials, n_bins, n_units)`` --
the layout the rest of the pipeline expects.

Animal ids are 1-indexed (matching the MATLAB convention). Each TF*_export.mat
holds one animal; the loader infers the animal id from the filename digits
when available, otherwise from the file ordering (after a stable sort).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from . import config


# ---------------------------------------------------------------------------
# Raw per-animal record
# ---------------------------------------------------------------------------
@dataclass
class Animal:
    """One animal's preprocessed data, as loaded from a TF*_export.mat."""

    animal_id: int                       # 1-indexed
    spatial_fr: np.ndarray               # (n_trials, n_bins, n_units), float
    area_masks: dict[str, np.ndarray]    # area -> (n_units,) bool
    fs_mask: np.ndarray                  # (n_units,) bool -- raw FS flag
    n_trials: int
    filename: str = ""

    @property
    def n_units(self) -> int:
        return self.spatial_fr.shape[2]

    @property
    def n_bins(self) -> int:
        return self.spatial_fr.shape[1]


# ---------------------------------------------------------------------------
# .mat helpers
# ---------------------------------------------------------------------------
def _h5_string(node) -> str:
    """Decode a MATLAB v7.3 cellstr / char array into a Python string."""
    arr = np.asarray(node).ravel()
    if arr.dtype.kind in ("U", "S"):
        return "".join(arr.astype(str).tolist())
    return "".join(chr(int(v)) for v in arr)


def _read_region_labels(file_handle, units_group) -> list[str]:
    """Read ``units.regions_label`` as a list of strings (one per area row).

    In MATLAB v7.3 a cellstr is stored as an array of object references; each
    reference points to a uint16 (or uint8) char array.
    """
    raw = units_group["regions_label"]
    labels: list[str] = []
    # Either a flat array of refs (cellstr) or a single string.
    arr = np.asarray(raw)
    if arr.dtype == object or arr.dtype.kind == "O" or h5py.check_dtype(ref=raw.dtype):
        for ref in arr.ravel():
            labels.append(_h5_string(file_handle[ref]))
    else:
        labels.append(_h5_string(arr))
    return labels


def _read_units(file_handle, units_group) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Read ``(area_idx, region_labels, fs_flag)`` from ``units``.

    Returns
    -------
    area_idx : ndarray (n_areas, n_units) bool
    region_labels : list of str, length n_areas
    fs_flag : ndarray (n_units,) bool -- True for fast-spiking
    """
    idx = np.asarray(units_group["idx"]).astype(bool)
    # h5py reads MATLAB matrices transposed; idx is (n_areas, n_units) on disk
    # but stored as (n_units, n_areas) by h5py. Detect and align by matching
    # the regions_label length on the area axis.
    labels = _read_region_labels(file_handle, units_group)
    if idx.shape[0] != len(labels):
        idx = idx.T
    if idx.shape[0] != len(labels):
        raise ValueError(
            f"units.idx shape {idx.shape} does not match {len(labels)} areas"
        )

    if "idx_fs" in units_group:
        fs = np.asarray(units_group["idx_fs"]).astype(bool).ravel()
    else:
        fs = np.zeros(idx.shape[1], dtype=bool)
    if fs.size != idx.shape[1]:
        raise ValueError(
            f"units.idx_fs size {fs.size} does not match {idx.shape[1]} units"
        )
    return idx, labels, fs


def _read_spatial_fr(file_handle, analysis_spatial, field: str) -> np.ndarray:
    """Read ``analysis_spatial.firing.cued.<field>`` as (n_trials, n_bins, n_units).

    MATLAB stores ``freq`` as (n_units, n_trials, n_bins); h5py reads it with
    the axes reversed, so the on-disk read is (n_bins, n_trials, n_units). The
    return shape is the (trials, bins, units) layout the pipeline expects.
    """
    firing = analysis_spatial["firing"]["cued"]
    if field not in firing:
        raise KeyError(f"analysis_spatial.firing.cued.{field} not found")
    raw = np.asarray(firing[field], dtype=float)
    if raw.ndim != 3:
        raise ValueError(f"expected 3-D spatial firing array, got {raw.shape}")
    # h5py read order is reversed vs. MATLAB; freq -> (n_bins, n_trials, n_units)
    # in Python, which we transpose to (n_trials, n_bins, n_units).
    return np.transpose(raw, (1, 0, 2))


# ---------------------------------------------------------------------------
# Per-animal loader
# ---------------------------------------------------------------------------
def _infer_animal_id(filename: str, fallback: int) -> int:
    """Extract the digits from a TF<id>_export.mat filename; fall back to index."""
    m = re.search(r"\d+", filename)
    if m:
        return int(m.group())
    return fallback


def load_animal(path: str | Path, animal_id: int | None = None,
                spatial_field: str = "freq") -> Animal:
    """Load one ``TF*_export.mat`` file."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"animal file not found: {path}")
    aid = animal_id if animal_id is not None else _infer_animal_id(path.name, 1)
    with h5py.File(path, "r") as f:
        if "units" not in f or "analysis_spatial" not in f:
            raise ValueError(
                f"{path.name} missing required fields (units, analysis_spatial)"
            )
        area_idx, region_labels, fs_flag = _read_units(f, f["units"])
        spatial_fr = _read_spatial_fr(f, f["analysis_spatial"], spatial_field)

    n_units = spatial_fr.shape[2]
    if area_idx.shape[1] != n_units:
        raise ValueError(
            f"{path.name}: units.idx has {area_idx.shape[1]} units but "
            f"spatial firing has {n_units}"
        )

    masks = {area: np.zeros(n_units, dtype=bool) for area in config.AREAS}
    for row, label in enumerate(region_labels):
        if label in masks:
            masks[label] |= area_idx[row]
    return Animal(
        animal_id=aid,
        spatial_fr=spatial_fr,
        area_masks=masks,
        fs_mask=fs_flag,
        n_trials=spatial_fr.shape[0],
        filename=path.name,
    )


# ---------------------------------------------------------------------------
# Cohort loader (all TF*_export.mat in a directory)
# ---------------------------------------------------------------------------
def _read_behaviour_file(path: Path) -> dict:
    """Read ``animal_behaviour.mat`` -- LP per animal.

    Returns a flat lookup dict ``('period_experienced', animal_id) -> int LP``.
    """
    lookup: dict = {}
    if not path.is_file():
        return lookup
    with h5py.File(path, "r") as f:
        if "animal_id" not in f or "period_experienced" not in f:
            return lookup
        animal_ids = np.asarray(f["animal_id"]).astype(int).ravel()
        # MATLAB saves ``period_experienced`` as (n_animals, 2) -- column 1 is
        # the LP, column 2 is the corresponding ratio (HC_V1_Code/HC_V1_temporal.m
        # line 37). h5py reads it with the axes reversed -> (2, n_animals).
        pe = np.asarray(f["period_experienced"], dtype=float)
        if pe.ndim == 1:
            pe = pe[:, None]
        # Prefer the h5py-transposed shape (n_cols, n_animals); fall back to
        # the in-MATLAB orientation (n_animals, n_cols). When the matrix is
        # square the orientation is ambiguous -- trust the MATLAB-via-h5py
        # convention (column-major on disk -> transposed read).
        if pe.shape[1] == animal_ids.size:
            lps = pe[0, :]
        elif pe.shape[0] == animal_ids.size:
            lps = pe[:, 0]
        else:
            raise ValueError(
                f"period_experienced shape {pe.shape} does not match "
                f"{animal_ids.size} animals"
            )
        for aid, lp in zip(animal_ids.tolist(), lps.tolist()):
            lookup[("period_experienced", int(aid))] = (
                int(round(lp)) if np.isfinite(lp) else None
            )
        if "zscored_lick_errors" in f:
            # Optional per-trial z, indexed by animal_id (object-ref array).
            zarr = f["zscored_lick_errors"]
            for i, aid in enumerate(animal_ids.tolist()):
                ref = zarr[i, 0] if zarr.ndim == 2 else zarr[i]
                lookup[("zscored_lick_errors", int(aid))] = (
                    np.asarray(f[ref], dtype=float).ravel()
                )
    return lookup


def load_animals(data_dir: str | Path | None = None,
                 spatial_field: str = "freq") -> list[Animal]:
    """Load every ``TF*_export.mat`` in ``data_dir`` (default config.DATA_DIR).

    Each file becomes one :class:`Animal`. The accompanying
    ``animal_behaviour.mat`` (if present) supplies the cohort learning points
    and per-trial behaviour.
    """
    data_dir = Path(data_dir) if data_dir is not None else config.DATA_DIR
    if not data_dir.is_dir():
        raise FileNotFoundError(f"data directory not found: {data_dir}")
    files = sorted(data_dir.glob(config.DATA_FILE_PATTERN))
    if not files:
        raise FileNotFoundError(
            f"no {config.DATA_FILE_PATTERN} files in {data_dir}"
        )
    animals: list[Animal] = []
    for i, path in enumerate(files, start=1):
        aid = _infer_animal_id(path.name, i)
        animals.append(load_animal(path, animal_id=aid,
                                   spatial_field=spatial_field))
    return animals


# ---------------------------------------------------------------------------
# Learning point + cohort classification
# ---------------------------------------------------------------------------
# Learning points come entirely from Tom's cohort behaviour file
# (``animal_behaviour.mat``; ``period_experienced(:, 1)`` indexed by
# ``animal_id``). The user committed to "whatever LP Tom is giving" -- there is
# no Python-side LP detection and no LP-criterion sweep. Animals without a
# recorded LP are dropped from the cohort (committed: learners only).


@dataclass
class CohortEntry:
    """Cohort role and learning point for one animal."""

    animal_id: int
    role: str             # "learner" | "nonlearner"
    lp: int               # learning point used for epochs (real or yoked)
    raw_lp: int | None    # recorded LP from Tom's file, or None


def _per_animal_lp(animal: Animal, behaviour_lookup: dict | None) -> int | None:
    """The LP Tom's cohort file records for this animal, or None."""
    if behaviour_lookup is None:
        return None
    return behaviour_lookup.get(("period_experienced", animal.animal_id))


def classify_cohort(animals: list[Animal], cfg,
                    behaviour_lookup: dict | None = None
                    ) -> dict[int, CohortEntry]:
    """Return the LEARNER subset, keyed by ``animal_id``.

    A learner has a recorded LP in Tom's cohort behaviour file and is not in
    ``cfg.manual_nonlearners``. Animals without a recorded LP are dropped
    (committed: "only work with the learners" -- no yoked non-learner branch).
    """
    entries: dict[int, CohortEntry] = {}
    for a in animals:
        lp = _per_animal_lp(a, behaviour_lookup)
        if lp is None or a.animal_id in cfg.manual_nonlearners:
            continue
        entries[a.animal_id] = CohortEntry(a.animal_id, "learner", lp, lp)
    return entries


# ---------------------------------------------------------------------------
# Epochs
# ---------------------------------------------------------------------------
def n_usable_trials(animal: Animal) -> int:
    """Number of trials usable for epoch analysis.

    Tom's pipeline does not provide a disengagement / change-point trial in
    the per-animal export, so the full trial count is taken as usable. Adapt
    here if a per-animal truncation index is added later.
    """
    return int(animal.n_trials)


def epoch_windows(lp: int, n_usable: int, cfg) -> dict[str, np.ndarray] | None:
    """0-indexed trial indices for the three learning epochs.

    naive        = trials 1..trials_per_epoch       (0-indexed 0..e-1)
    intermediate = trials_per_epoch trials ending at lp  (0-indexed lp-e .. lp-1)
    expert       = trials_per_epoch trials after lp  (0-indexed lp .. lp+e-1)

    Returns None if naive and intermediate would overlap (lp < 2 * e) or the
    expert window would exceed the usable trials.
    """
    e = cfg.trials_per_epoch
    if lp < 2 * e:
        return None
    if lp + e > n_usable:
        return None
    return {
        "naive": np.arange(0, e),
        "intermediate": np.arange(lp - e, lp),
        "expert": np.arange(lp, lp + e),
    }


# ---------------------------------------------------------------------------
# Per-area activity tensors
# ---------------------------------------------------------------------------
def select_units(animal: Animal, area: str, cfg) -> np.ndarray:
    """Indices of the units in ``area`` that survive selection.

    Keeps units flagged for the area. When ``cfg.exclude_fast_spiking`` is set,
    drops fast-spiking units in the four areas Tom flags FS in (V1/RSC/CA1/CA3,
    :data:`config.FS_AREAS`); other areas are unaffected.
    """
    keep = animal.area_masks[area].copy()
    if cfg.exclude_fast_spiking and area in config.FS_AREAS:
        keep &= ~animal.fs_mask
    return np.where(keep)[0]


def area_tensor(animal: Animal, area: str, cfg) -> tuple[np.ndarray, np.ndarray]:
    """Build one area's activity tensor.

    Returns ``(activity, unit_indices)`` where ``activity`` is
    ``(n_usable_trials, n_bins, n_units_kept)``.
    """
    idx = select_units(animal, area, cfg)
    n_use = n_usable_trials(animal)
    return animal.spatial_fr[:n_use][:, :, idx], idx
