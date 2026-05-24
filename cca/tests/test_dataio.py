"""Tests for the Tom dataio module.

Two layers:
1. Unit tests for the cohort / epoch logic that don't touch .mat files.
2. Round-trip tests that write a synthetic ``TF*_export.mat`` (and a matching
   ``animal_behaviour.mat``) to a tmp directory, load it with
   :func:`tom_cca.dataio.load_animal`, and check the in-memory layout.
"""

from __future__ import annotations

import dataclasses

import h5py
import numpy as np
import pytest

from tom_cca import config, dataio

CFG = config.DEFAULT


# ---------------------------------------------------------------------------
# Synthetic animal (in-memory)
# ---------------------------------------------------------------------------
def make_animal(animal_id, n_trials, units_per_area, fs_units=()):
    """Build an Animal directly (no .mat round-trip)."""
    n_units = sum(units_per_area.values())
    rng = np.random.default_rng(animal_id)
    spatial_fr = rng.standard_normal((n_trials, config.N_BINS, n_units))
    masks = {a: np.zeros(n_units, dtype=bool) for a in config.AREAS}
    start = 0
    for area in config.AREAS:
        count = units_per_area.get(area, 0)
        masks[area][start : start + count] = True
        start += count
    fs_mask = np.zeros(n_units, dtype=bool)
    for u in fs_units:
        fs_mask[u] = True
    return dataio.Animal(
        animal_id=animal_id,
        spatial_fr=spatial_fr,
        area_masks=masks,
        fs_mask=fs_mask,
        n_trials=n_trials,
    )


# ---------------------------------------------------------------------------
# Cohort classification -- learners only, LP from Tom's behaviour file
# ---------------------------------------------------------------------------
def test_classify_cohort_returns_only_learners_with_recorded_lp():
    a1 = make_animal(1, 130, {"CA1": 10})
    a2 = make_animal(2, 130, {"CA1": 10})
    a3 = make_animal(3, 130, {"CA1": 10})
    # Animal 3 has no recorded LP -- dropped (no yoking).
    behaviour = {
        ("period_experienced", 1): 40,
        ("period_experienced", 2): 60,
        ("period_experienced", 3): None,
    }
    entries = dataio.classify_cohort([a1, a2, a3], CFG, behaviour)
    assert set(entries) == {1, 2}
    assert entries[1].role == "learner" and entries[1].lp == 40
    assert entries[2].role == "learner" and entries[2].lp == 60


def test_classify_cohort_drops_animals_without_any_recorded_lp():
    a = make_animal(1, 100, {"CA1": 10})
    # No behaviour entry at all -> not a learner.
    entries = dataio.classify_cohort([a], CFG, behaviour_lookup={})
    assert entries == {}


def test_classify_cohort_drops_when_behaviour_file_missing():
    a = make_animal(1, 100, {"CA1": 10})
    entries = dataio.classify_cohort([a], CFG, behaviour_lookup=None)
    assert entries == {}


def test_classify_cohort_honours_manual_nonlearners():
    a = make_animal(7, 130, {"CA1": 10})
    other = make_animal(1, 130, {"CA1": 10})
    cfg = dataclasses.replace(CFG, manual_nonlearners=frozenset({7}))
    behaviour = {("period_experienced", 7): 30, ("period_experienced", 1): 50}
    entries = dataio.classify_cohort([a, other], cfg, behaviour)
    # Animal 7 is forced out even though it has a recorded LP.
    assert set(entries) == {1}
    assert entries[1].lp == 50


# ---------------------------------------------------------------------------
# Usable trials / epoch windows
# ---------------------------------------------------------------------------
def test_n_usable_trials_uses_full_trial_count():
    a = make_animal(1, 90, {"CA1": 5})
    assert dataio.n_usable_trials(a) == 90


def test_epoch_windows_are_disjoint_and_correct():
    win = dataio.epoch_windows(lp=42, n_usable=100, cfg=CFG)
    assert win is not None
    assert set(win) == {"naive", "intermediate", "expert"}
    assert np.array_equal(win["naive"], np.arange(0, 10))
    assert np.array_equal(win["intermediate"], np.arange(32, 42))
    assert np.array_equal(win["expert"], np.arange(42, 52))


def test_epoch_windows_boundary_lp_equals_twice_trials_per_epoch():
    win = dataio.epoch_windows(lp=20, n_usable=100, cfg=CFG)
    assert win is not None
    assert np.array_equal(win["intermediate"], np.arange(10, 20))
    assert np.array_equal(win["expert"], np.arange(20, 30))


def test_epoch_windows_reject_overlap_when_lp_too_early():
    assert dataio.epoch_windows(lp=19, n_usable=100, cfg=CFG) is None


def test_epoch_windows_reject_when_expert_exceeds_trials():
    assert dataio.epoch_windows(lp=42, n_usable=48, cfg=CFG) is None


# ---------------------------------------------------------------------------
# Unit selection / area tensor
# ---------------------------------------------------------------------------
def test_select_units_excludes_fast_spiking_in_fs_areas_only():
    # FS flag on unit 2 (V1 -> FS area, dropped) and unit 12 (DG -> kept).
    a = make_animal(1, 60, {"V1": 10, "DG": 10}, fs_units=(2, 12))
    assert 2 not in dataio.select_units(a, "V1", CFG).tolist()
    assert 12 in dataio.select_units(a, "DG", CFG).tolist()


def test_select_units_keeps_fast_spiking_when_disabled():
    a = make_animal(1, 60, {"V1": 10}, fs_units=(3,))
    cfg = dataclasses.replace(CFG, exclude_fast_spiking=False)
    assert dataio.select_units(a, "V1", cfg).size == 10


def test_area_tensor_shape_and_truncation():
    a = make_animal(1, 80, {"V1": 12, "DG": 8}, fs_units=(0,))
    tensor, idx = dataio.area_tensor(a, "V1", CFG)
    assert tensor.shape == (80, config.N_BINS, 11)
    assert idx.size == 11


# ---------------------------------------------------------------------------
# Committed defaults the user pinned
# ---------------------------------------------------------------------------
def test_committed_defaults():
    assert CFG.zscore_units is True
    assert CFG.min_units == 5
    assert CFG.exclude_fast_spiking is True
    assert not hasattr(CFG, "lp_min_consecutive")     # LP detection removed


# ---------------------------------------------------------------------------
# .mat round-trip
# ---------------------------------------------------------------------------
def _write_tom_mat(path, *, n_units, n_trials, n_bins, region_for_unit,
                   fs_units=()):
    """Write a v7.3 .mat file mimicking Tom's TF*_export.mat layout.

    ``region_for_unit`` is a list of length n_units giving each unit's area.
    Stores:
      - ``units.idx``    : (n_areas, n_units) uint8
      - ``units.regions_label`` : cellstr of area names (one per area row)
      - ``units.idx_fs`` : (n_units,) uint8
      - ``analysis_spatial.firing.cued.freq`` : (n_units, n_trials, n_bins)
    """
    regions = list(dict.fromkeys(region_for_unit))   # preserve first-seen order
    area_idx = np.zeros((len(regions), n_units), dtype=np.uint8)
    for i, area in enumerate(region_for_unit):
        area_idx[regions.index(area), i] = 1
    fs_flag = np.zeros(n_units, dtype=np.uint8)
    for u in fs_units:
        fs_flag[u] = 1

    rng = np.random.default_rng(0)
    freq = rng.standard_normal((n_units, n_trials, n_bins)).astype(np.float64)

    with h5py.File(path, "w") as f:
        units = f.create_group("units")
        # h5py writes arrays as-is; MATLAB stores idx as (n_areas, n_units), so
        # writing the *transposed* array makes the on-disk shape match Tom's.
        units.create_dataset("idx", data=area_idx.T)
        units.create_dataset("idx_fs", data=fs_flag)

        ref_dtype = h5py.special_dtype(ref=h5py.Reference)
        label_refs = units.create_dataset("regions_label", (len(regions), 1),
                                          dtype=ref_dtype)
        for i, area in enumerate(regions):
            ds = f.create_dataset(f"_label_{i}",
                                  data=np.frombuffer(area.encode("utf-16-le"),
                                                     dtype=np.uint16))
            label_refs[i, 0] = ds.ref

        # analysis_spatial.firing.cued.freq stored with MATLAB axis order
        # (n_units, n_trials, n_bins) -- which means h5py writes
        # (n_bins, n_trials, n_units) on disk.
        ana = f.create_group("analysis_spatial")
        firing = ana.create_group("firing")
        cued = firing.create_group("cued")
        cued.create_dataset("freq", data=np.transpose(freq, (2, 1, 0)))


def _write_behaviour_mat(path, animal_ids, learning_points):
    """Write a minimal ``animal_behaviour.mat`` with period_experienced + animal_id."""
    pe = np.column_stack([np.asarray(learning_points, dtype=float),
                          np.zeros(len(learning_points))])
    with h5py.File(path, "w") as f:
        # MATLAB axis-reverses on write -- a (n, 2) on disk becomes (2, n) here.
        f.create_dataset("period_experienced", data=pe.T)
        f.create_dataset("animal_id", data=np.asarray(animal_ids, dtype=np.int64))


def test_load_animal_round_trip(tmp_path):
    n_units, n_trials, n_bins = 20, 50, config.N_BINS
    regions = (["CA1"] * 10) + (["V1"] * 10)
    path = tmp_path / "TF03_export.mat"
    _write_tom_mat(path, n_units=n_units, n_trials=n_trials, n_bins=n_bins,
                   region_for_unit=regions, fs_units=(2, 11))

    a = dataio.load_animal(path)
    assert a.animal_id == 3
    assert a.n_trials == n_trials
    assert a.n_units == n_units
    assert a.n_bins == n_bins
    assert a.spatial_fr.shape == (n_trials, n_bins, n_units)
    assert int(a.area_masks["CA1"].sum()) == 10
    assert int(a.area_masks["V1"].sum()) == 10
    # FS in CA1 (an FS area) gets dropped at select time; FS in V1 too.
    assert dataio.select_units(a, "CA1", CFG).size == 9      # unit 2 dropped
    assert dataio.select_units(a, "V1", CFG).size == 9       # unit 11 dropped


def test_load_animals_with_behaviour_uses_recorded_lp(tmp_path):
    for i in (1, 2, 3):
        _write_tom_mat(tmp_path / f"TF{i:02d}_export.mat", n_units=10,
                       n_trials=80, n_bins=config.N_BINS,
                       region_for_unit=(["CA1"] * 6) + (["V1"] * 4))
    _write_behaviour_mat(tmp_path / config.LEARNING_FILE,
                         animal_ids=[1, 2, 3], learning_points=[35, 45, 55])

    animals = dataio.load_animals(tmp_path)
    behaviour = dataio._read_behaviour_file(tmp_path / config.LEARNING_FILE)
    entries = dataio.classify_cohort(animals, CFG, behaviour)
    assert set(entries) == {1, 2, 3}
    assert entries[1].lp == 35
    assert entries[2].lp == 45
    assert entries[3].lp == 55


def test_load_animal_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        dataio.load_animal(tmp_path / "does_not_exist.mat")


def test_load_animals_with_no_matches_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        dataio.load_animals(tmp_path)
