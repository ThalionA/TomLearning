"""Running-state preprocessing shared by every temporal driver.

One animal → one :class:`RunningSession`: the engaged running bins (in-trial and
velocity ≥ threshold) of every area with enough units, z-scored, plus the per-bin trial
id. Drivers then ask for a pair (``X, Y, Z``) and run their own method on it.

Before 2026-08-19 this block — load the 1 ms streams, build the running mask, check the
trial count, select areas with ≥ ``min_units``, apply the sample cap, concatenate the
"other" areas into the confound ``Z`` — was copy-pasted into 13 drivers, with the sample
cap implemented four times (two different selection rules) and a silent
``except Exception: continue`` that made an animal vanish from the cohort without a
message (audit/REPORT_2026-08-19.md §2.4, §2.5, §4.2). Now every skip has a reason.

Preprocessing contract (HANDOFF.md §2, unchanged):
  1. 1 ms spikes → Gaussian ``cfg.gaussian_sd_ms`` → ``cfg.temporal_bin_ms`` bins
     (``dataio.area_activity_50ms``, z-scored on the engaged reference).
  2. Running bins = in-trial AND velocity ≥ ``cfg.velocity_thresh_cm_s``.
  3. An area needs ≥ ``cfg.min_units`` units; FS exclusion per ``cfg.exclude_fast_spiking``.
  4. Optional cap: consecutive leading block per trial, earliest trials first
     (:func:`cap_running_bins`; default uncapped since 2026-08-17).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import config, core, dataio, partial


@dataclass
class RunningSession:
    """Engaged running bins of one animal, every qualifying area, trial-labelled."""

    animal_id: int
    learner: bool
    lp: int | None                       # learning point (trial index) if a learner
    areas: dict[str, np.ndarray]         # area -> (n_bins, n_units) z-scored activity
    trial_ids: np.ndarray                # (n_bins,) trial id per kept running bin
    n_running_total: int                 # running bins before the cap
    cap_index: np.ndarray                # indices (into the running-bin axis) kept
    unit_counts: dict[str, int] = field(default_factory=dict)
    trials_full: np.ndarray | None = None   # unique running-trial ids BEFORE the cap
    velocity: np.ndarray | None = None      # (n_bins,) cm/s per kept running bin

    @property
    def n_bins(self) -> int:
        return int(self.trial_ids.size)

    @property
    def trials(self) -> np.ndarray:
        return np.unique(self.trial_ids)

    def subset(self, sel) -> RunningSession:
        """The same session restricted to running-bin indices / mask ``sel``."""
        sel = np.asarray(sel)
        if sel.dtype == bool:
            sel = np.flatnonzero(sel)
        return RunningSession(
            animal_id=self.animal_id, learner=self.learner, lp=self.lp,
            areas={a: m[sel] for a, m in self.areas.items()},
            trial_ids=self.trial_ids[sel], n_running_total=self.n_running_total,
            cap_index=self.cap_index[sel], unit_counts=dict(self.unit_counts),
            trials_full=self.trials_full,
            velocity=self.velocity[sel] if self.velocity is not None else None)

    def pair(self, ax: str, ay: str):
        """``(X, Y, Z)`` for one pair — ``Z`` = the other present areas, concatenated in
        ``config.AREAS`` order (``None`` when no third area). ``None`` if either area is
        missing. This is the confound every driver partials out."""
        if ax not in self.areas or ay not in self.areas:
            return None
        others = [m for a, m in self.areas.items() if a not in (ax, ay)]
        Z = np.concatenate(others, axis=1) if others else None
        return self.areas[ax], self.areas[ay], Z


def cap_running_bins(trial_ids, max_samples: int = 0, max_per_trial: int = 0) -> np.ndarray:
    """Indices keeping a CONSECUTIVE leading block (≤ ``max_per_trial`` bins) of each
    trial, earliest trials first, stopping once ``max_samples`` is reached.

    ``0`` for either limit means "no limit"; with both 0 every bin is kept (the default
    since 2026-08-17). Consecutive blocks preserve the bin adjacency the segment-aware
    lag pairing needs — an even stride would destroy it (GOTCHAS.md). ⚠ Any cap takes the
    EARLIEST trials: at 12 000 bins / 600 per trial that is the first ~20 trials.

    Reproduces all three pre-2026-08-19 driver versions: ``run_lag_curves._capped_index``
    (this signature), ``run_lag_subspaces`` / ``run_lag_cosine._capped_index(12000, 600)``
    and ``run_fixed_subspace._cap_bins(600)`` ≡ ``(0, 600)`` — see
    ``tests/test_dedup_equivalence.py``.
    """
    trial_ids = np.asarray(trial_ids)
    if not max_samples and not max_per_trial:
        return np.arange(trial_ids.size)
    keep, total = [], 0
    for t in np.unique(trial_ids):
        pos = np.flatnonzero(trial_ids == t)               # temporal order within trial
        if max_per_trial:
            pos = pos[:max_per_trial]
        keep.append(pos)
        total += pos.size
        if max_samples and total >= max_samples:
            break
    if not keep:
        return np.zeros(0, dtype=int)
    return np.sort(np.concatenate(keep))


def load_running_session(animal, cfg, entries: dict, *, max_samples: int = 0,
                         max_per_trial: int = 0, min_trials: int = config.TEMPORAL.min_trials,
                         ) -> RunningSession | str:
    """Build the :class:`RunningSession` for ``animal`` or return a one-line skip reason.

    ``entries`` is the learner cohort from :func:`dataio.classify_cohort` (keyed by
    animal id; ``entries[id].lp`` is the learning point). Areas with fewer than
    ``cfg.min_units`` units are omitted; the caller decides which pairs are therefore
    unavailable via :meth:`RunningSession.pair`.
    """
    try:
        streams = dataio._load_temporal_streams(animal, cfg)
    except Exception as e:                                  # noqa: BLE001
        return f"stream load failed: {type(e).__name__}: {e}"
    run = (~np.isnan(streams.trial_idx_50ms)) & (streams.vel_50ms >= cfg.velocity_thresh_cm_s)
    trial_ids_full = streams.trial_idx_50ms[run]
    n_trials = np.unique(trial_ids_full).size
    if n_trials < min_trials:
        return f"only {n_trials} running trials (< {min_trials})"
    cap = cap_running_bins(trial_ids_full, max_samples, max_per_trial)
    areas, counts = {}, {}
    for area in config.AREAS:
        m, idx = dataio.area_activity_50ms(animal, area, cfg)
        counts[area] = int(len(idx))
        if len(idx) >= cfg.min_units:
            areas[area] = m[run][cap]
    entry = entries.get(animal.animal_id)
    lp = int(entry.lp) if entry is not None and getattr(entry, "lp", None) is not None else None
    return RunningSession(
        animal_id=int(animal.animal_id), learner=entry is not None, lp=lp,
        areas=areas, trial_ids=np.asarray(trial_ids_full[cap]),
        n_running_total=int(trial_ids_full.size), cap_index=cap, unit_counts=counts,
        trials_full=np.unique(trial_ids_full),
        velocity=np.asarray(streams.vel_50ms[run][cap], dtype=float))


def residual_pca_scores(X, Y, Z, k: int):
    """Refit-per-lag family prep: partial the third areas out of both (in-sample, all
    rows), then PCA to ``k`` per area → ``(Sx, Sy)``. Used by ``run_lag_curves``,
    ``run_lag_cosine``, ``run_ifi_windows``; the frozen-axes drivers residualise on
    their own row subset and call :func:`core.pca_fit_flat` via ``fixed_subspace``.
    NOTE (HANDOFF §2.4 correction, 2026-08-19): the confound coefficients here see every
    row, including the folds later held out — unsupervised w.r.t. the X–Y relation but
    not train-only. ``lag_subspace.lag_sweep`` does it per fold.
    """
    Xr = partial.partial_out(X, Z) if Z is not None else X
    Yr = partial.partial_out(Y, Z) if Z is not None else Y
    Sx, _ = core.pca_scores(Xr, min(k, Xr.shape[1]))
    Sy, _ = core.pca_scores(Yr, min(k, Yr.shape[1]))
    return Sx, Sy


def epoch_of_trial(lp: int, trials, cfg, epochs=config.EPOCH_NAMES) -> dict[int, str] | None:
    """Map trial id → epoch name for the three-epoch design (``dataio.epoch_windows`` on
    the running trials present). ``None`` if the learning point leaves no room for the
    epochs. Shared by the frozen-axes drivers (was two inline copies)."""
    trials = np.unique(np.asarray(trials))
    ew = dataio.epoch_windows(int(lp), trials.size, cfg)
    if ew is None:
        return None
    out: dict[int, str] = {}
    for epoch in epochs:
        pos = np.asarray(ew[epoch])
        pos = pos[(pos >= 0) & (pos < trials.size)]
        for t in trials[pos]:
            out[int(t)] = epoch
    return out or None
