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

import json
import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from . import config, landmark_align, segments


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
    # Path back to the .mat file, used by the temporal arms to lazy-load the
    # 1 ms ``binned_spikes`` and ``data_behaviour`` streams on demand. Set by
    # load_animal; None when an Animal is constructed in-memory (test fixtures).
    streams_path: Path | None = None

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


def _read_units(file_handle, units_group,
                region_labels: list[str] | None = None
                ) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Read ``(area_idx, region_labels, fs_flag)`` from ``units``.

    Parameters
    ----------
    region_labels : list of str, optional
        Area names, one per row of ``units.idx``. When supplied (from the
        JSON companion file -- see :data:`config.COMPANION_FILE`), the labels
        are used verbatim and ``units.regions_label`` is not read from the
        .mat file. This is the path taken for Tom's real exports, where
        ``regions_label`` is a MATLAB ``string`` array h5py cannot decode.
        When ``None`` the labels are read from the file (cellstr layout --
        used by the synthetic test fixtures and any cellstr-typed export).

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
    if region_labels is not None:
        labels = list(region_labels)
    else:
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
                spatial_field: str = "freq",
                region_labels: list[str] | None = None) -> Animal:
    """Load one ``TF*_export.mat`` file.

    ``region_labels`` (one per ``units.idx`` row), when supplied, overrides the
    in-file ``units.regions_label`` -- the path used for Tom's real exports,
    whose ``regions_label`` is a MATLAB ``string`` array (see
    :func:`_read_units`). :func:`load_animals` supplies it from the JSON
    companion file automatically.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"animal file not found: {path}")
    aid = animal_id if animal_id is not None else _infer_animal_id(path.name, 1)
    with h5py.File(path, "r") as f:
        if "units" not in f or "analysis_spatial" not in f:
            raise ValueError(
                f"{path.name} missing required fields (units, analysis_spatial)"
            )
        area_idx, region_labels, fs_flag = _read_units(
            f, f["units"], region_labels=region_labels)
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
        streams_path=path,
    )


# ---------------------------------------------------------------------------
# Cohort loader (all TF*_export.mat in a directory)
# ---------------------------------------------------------------------------
def _read_companion(data_dir: str | Path) -> dict | None:
    """Parse the JSON companion file in ``data_dir``, or None if it is absent.

    The companion file (:data:`config.COMPANION_FILE`, written by
    ``scripts/export_cca_labels.m``) carries the region labels and learning
    points that Tom's exports store as MATLAB ``string`` arrays -- a type
    h5py cannot decode. When present it is the authoritative source for both.

    Returns
    -------
    dict or None
        ``{"regions": {animal_id: [area, ...]},
           "lp":      {animal_id: int | None}}``  or ``None`` when no
        companion file exists in ``data_dir``.
    """
    path = Path(data_dir) / config.COMPANION_FILE
    if not path.is_file():
        return None
    data = json.loads(path.read_text())
    regions: dict[int, list[str]] = {}
    for entry in data.get("animals", []):
        regions[int(entry["id"])] = [str(r) for r in entry["regions"]]
    lp: dict[int, int | None] = {}
    for entry in data.get("behaviour", []):
        val = entry.get("lp")
        lp[int(entry["id"])] = (
            int(round(val)) if val is not None and np.isfinite(val) else None
        )
    return {"regions": regions, "lp": lp}


def _read_behaviour_file(path: Path) -> dict:
    """Read learning points -- one per animal.

    Returns a flat lookup dict ``('period_experienced', animal_id) -> int LP``.
    When a JSON companion file sits beside ``path`` it is used (Tom's real
    cohort); otherwise the learning points are read from ``animal_behaviour.mat``
    directly (the numeric layout the synthetic test fixtures write).
    """
    path = Path(path)
    companion = _read_companion(path.parent)
    if companion is not None:
        return {("period_experienced", aid): lp
                for aid, lp in companion["lp"].items()}

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
    # When the JSON companion is present it supplies region labels (Tom's
    # exports store regions_label as an undecodable MATLAB `string`); when
    # absent each file's labels are read directly (cellstr test fixtures).
    companion = _read_companion(data_dir)
    animals: list[Animal] = []
    for i, path in enumerate(files, start=1):
        aid = _infer_animal_id(path.name, i)
        region_labels = None
        if companion is not None:
            region_labels = companion["regions"].get(aid)
            if region_labels is None:
                raise ValueError(
                    f"{path.name}: animal {aid} is missing from "
                    f"{config.COMPANION_FILE} -- re-run "
                    f"scripts/export_cca_labels.m to refresh it"
                )
        animals.append(load_animal(path, animal_id=aid,
                                   spatial_field=spatial_field,
                                   region_labels=region_labels))
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


# ===========================================================================
# Temporal arms (Arm A: running-state full-traversal; Arm B: landmark-centred)
# ===========================================================================
# 1 ms streams are read on demand from the TF*_export.mat and immediately
# binned to ``cfg.temporal_bin_ms`` (default 50 ms). Only the 50 ms result is
# cached -- the raw 1 ms ``binned_spikes`` is large (~1.6 GB per animal).
# Cache holds one animal; cleared when a different animal is requested. See
# ``cca/UNDERSTANDING_temporal.md`` Sec. 3 for the data contract.


def rebin_spikes(spikes_1ms, bin_ms: int,
                 chunk_out_bins: int = 2000,
                 max_1ms: int | None = None,
                 gaussian_sd_ms: float = 0.0) -> np.ndarray:
    """Sum 1 ms spike counts into ``bin_ms``-wide bins.

    ``(n_1ms, n_units) -> (n_bins, n_units)`` float32. An incomplete final
    bin is dropped.

    If ``gaussian_sd_ms > 0`` the 1 ms spike train of each unit is first
    convolved with a Gaussian of that s.d. (in ms = 1 ms samples), yielding a
    smoothed spike-density, *then* integrated into the bins — the Gonzalez &
    Buzsáki 2026 preprocessing ("convolved the spike times of each unit with a
    Gaussian kernel with a 2.5 ms s.d. ... yielding smoothed spiking activity";
    z-scoring is applied downstream). ``gaussian_filter1d`` conserves total mass,
    so smoothed bins remain spike-count-scaled. A ±4σ halo per chunk removes
    chunk-boundary edge effects (exact to numerical precision for σ≪chunk).

    Accepts either a numpy array or an h5py Dataset for ``spikes_1ms`` -- both
    support ``[start:stop]`` slicing and ``.shape``. Rebinning is done in
    chunks of ``chunk_out_bins`` bins so the intermediate float32 array stays
    bounded. ``max_1ms`` caps the input length (cross-stream alignment).
    """
    bin_ms = int(bin_ms)
    n_1ms = int(spikes_1ms.shape[0])
    if max_1ms is not None:
        n_1ms = min(n_1ms, int(max_1ms))
    n_units = int(spikes_1ms.shape[1])
    n_bins = n_1ms // bin_ms
    if n_bins == 0:
        return np.zeros((0, n_units), dtype=np.float32)
    smooth = float(gaussian_sd_ms) > 0
    if smooth:
        from scipy.ndimage import gaussian_filter1d
        halo = int(np.ceil(4 * float(gaussian_sd_ms)))    # 1 ms samples
    out = np.empty((n_bins, n_units), dtype=np.float32)
    for start in range(0, n_bins, int(chunk_out_bins)):
        stop = min(start + int(chunk_out_bins), n_bins)
        lo, hi = start * bin_ms, stop * bin_ms
        if smooth:
            a, b = max(0, lo - halo), min(n_1ms, hi + halo)
            sm = gaussian_filter1d(np.asarray(spikes_1ms[a:b]).astype(np.float32),
                                   sigma=float(gaussian_sd_ms), axis=0, mode="nearest")
            chunk = sm[lo - a: lo - a + (hi - lo)]
        else:
            chunk = np.asarray(spikes_1ms[lo:hi]).astype(np.float32)
        out[start:stop] = chunk.reshape(stop - start, bin_ms, n_units).sum(axis=1)
    return out


def bin_velocity(vel_1ms: np.ndarray, bin_ms: int) -> np.ndarray:
    """Mean velocity per 50 ms bin. ``(n_1ms,) -> (n_bins,)``."""
    vel = np.asarray(vel_1ms, dtype=float).ravel()
    n_bins = vel.size // int(bin_ms)
    if n_bins == 0:
        return np.zeros(0, dtype=float)
    return vel[: n_bins * int(bin_ms)].reshape(n_bins, int(bin_ms)).mean(axis=1)


def bin_trial_index(trial_1ms: np.ndarray, bin_ms: int) -> np.ndarray:
    """All-or-nothing 50 ms trial index.

    A 50 ms bin is in trial t iff *every* 1 ms bin in its span is in trial t.
    Returns NaN for bins that are outside any trial or that span a trial
    boundary (the conservative rule that keeps segments from straddling
    trial boundaries -- see ``segments.find_segments``).
    """
    trial = np.asarray(trial_1ms, dtype=float).ravel()
    n_bins = trial.size // int(bin_ms)
    if n_bins == 0:
        return np.zeros(0, dtype=float)
    blk = trial[: n_bins * int(bin_ms)].reshape(n_bins, int(bin_ms))
    first = blk[:, 0]
    # consistent <=> every cell in the row equals the first column AND is
    # finite. NaN == NaN is False so a NaN-in-first-column row collapses to
    # False here, which is exactly the desired "drop" behaviour.
    consistent = np.all(blk == first[:, None], axis=1)
    consistent &= np.all(np.isfinite(blk), axis=1)
    return np.where(consistent, first, np.nan)


@dataclass
class _TemporalStreams:
    """Binned 50 ms streams for one animal."""

    spikes_50ms: np.ndarray         # (n_bins, n_units) float32
    vel_50ms: np.ndarray            # (n_bins,) float
    trial_idx_50ms: np.ndarray      # (n_bins,) float, NaN where out-of-trial
    entries_1ms: np.ndarray         # (n_entries, 2) int -- (1ms_idx, landmark_id)


_TEMPORAL_CACHE: dict = {}


def _load_temporal_streams(animal: Animal, cfg) -> _TemporalStreams:
    """Read + bin the 1 ms streams for ``animal``. Cached per animal.

    Cache holds one animal at a time. Re-binning into ``cfg.temporal_bin_ms``
    happens on every reload (cache key includes the bin size).
    """
    key = (animal.animal_id, int(cfg.temporal_bin_ms), float(getattr(cfg, "gaussian_sd_ms", 0.0)))
    cached = _TEMPORAL_CACHE.get(key)
    if cached is not None:
        return cached
    _TEMPORAL_CACHE.clear()                                # 1-slot cache
    if animal.streams_path is None:
        raise ValueError(
            f"animal {animal.animal_id} has no streams_path -- cannot load "
            f"1 ms streams (in-memory Animal fixtures cannot drive the "
            f"temporal arms; use load_animal to attach the .mat path)"
        )
    bin_ms = int(cfg.temporal_bin_ms)
    with h5py.File(animal.streams_path, "r") as f:
        ds = f["binned_spikes"]                         # h5py Dataset, lazy
        beh = f["data_behaviour"]
        vel_1ms = np.asarray(beh["velocity_binned_gf"]).ravel()
        trial_1ms = np.asarray(beh["trial_binned_cued"]).ravel()
        landmark_1ms = np.asarray(
            beh["visual_landmark_binned"]).ravel().astype(int)
        # Align lengths -- the four streams should be 1:1 by construction;
        # trim to the common length defensively.
        n = min(int(ds.shape[0]), vel_1ms.size, trial_1ms.size, landmark_1ms.size)
        vel_1ms = vel_1ms[:n]
        trial_1ms = trial_1ms[:n]
        landmark_1ms = landmark_1ms[:n]
        # Chunked read directly from the h5py dataset -- the full 1.6 GB uint8
        # binned_spikes is never materialised in memory.
        spikes_50ms = rebin_spikes(ds, bin_ms, max_1ms=n,
                                   gaussian_sd_ms=getattr(cfg, "gaussian_sd_ms", 0.0))

    streams = _TemporalStreams(
        spikes_50ms=spikes_50ms,
        vel_50ms=bin_velocity(vel_1ms, bin_ms),
        trial_idx_50ms=bin_trial_index(trial_1ms, bin_ms),
        entries_1ms=landmark_align.find_entries(landmark_1ms),
    )
    _TEMPORAL_CACHE[key] = streams
    return streams


def _zscore_engaged(spikes_area: np.ndarray, vel_50ms: np.ndarray,
                    trial_idx_50ms: np.ndarray, vel_threshold: float,
                    require_running: bool) -> np.ndarray:
    """Per-unit z-score using the engaged-period reference.

    Reference = bins in any cued trial (require_running=False) or in any cued
    trial AND running >= threshold (require_running=True). The choice differs
    by arm: Arm A residualises against running-state distribution; Arm B uses
    the wider in-cued reference so the landmark windows include the full
    behavioural envelope around landmarks.
    """
    in_trial = ~np.isnan(trial_idx_50ms)
    if require_running:
        ref_mask = in_trial & (vel_50ms >= vel_threshold)
    else:
        ref_mask = in_trial
    if not np.any(ref_mask):
        return spikes_area
    ref = spikes_area[ref_mask]
    mu = ref.mean(axis=0)
    sd = ref.std(axis=0)
    sd[sd == 0] = 1.0
    return (spikes_area - mu) / sd


def temporal_segments(animal: Animal, cfg) -> list[segments.Segment]:
    """Velocity-passing segments for Arm A. Area-agnostic (same for X and Y)."""
    streams = _load_temporal_streams(animal, cfg)
    vel_mask = streams.vel_50ms >= cfg.velocity_thresh_cm_s
    return segments.find_segments(
        vel_mask,
        streams.trial_idx_50ms,
        min_gap_close_bins=config.min_gap_close_bins(cfg),
        min_run_bins=config.min_run_bins(cfg),
    )


def area_activity_50ms(animal: Animal, area: str, cfg
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Arm A: per-area 50 ms activity matrix, z-scored over running+cued bins.

    Returns ``(spikes_50ms[:, kept_units], kept_unit_indices)``.
    """
    streams = _load_temporal_streams(animal, cfg)
    idx = select_units(animal, area, cfg)
    spikes_area = streams.spikes_50ms[:, idx]
    if cfg.zscore_units:
        spikes_area = _zscore_engaged(
            spikes_area, streams.vel_50ms, streams.trial_idx_50ms,
            cfg.velocity_thresh_cm_s, require_running=True,
        )
    return spikes_area, idx


def area_landmark_windows(animal: Animal, area: str, cfg
                          ) -> tuple[dict[int, landmark_align.LandmarkWindows],
                                     np.ndarray]:
    """Arm B: per-landmark window stacks for one area's units.

    Returns ``(dict landmark_id -> LandmarkWindows, kept_unit_indices)``.
    """
    streams = _load_temporal_streams(animal, cfg)
    idx = select_units(animal, area, cfg)
    spikes_area = streams.spikes_50ms[:, idx]
    if cfg.zscore_units:
        spikes_area = _zscore_engaged(
            spikes_area, streams.vel_50ms, streams.trial_idx_50ms,
            cfg.velocity_thresh_cm_s, require_running=False,
        )
    windows = landmark_align.align_windows(
        spikes_area,
        streams.vel_50ms,
        streams.trial_idx_50ms,
        streams.entries_1ms,
        bin_ms=int(cfg.temporal_bin_ms),
        window_bins_pre=config.landmark_window_bins_pre(cfg),
        window_bins_post=config.landmark_window_bins_post(cfg),
        vel_threshold_cm_s=cfg.velocity_thresh_cm_s,
    )
    return windows, idx
