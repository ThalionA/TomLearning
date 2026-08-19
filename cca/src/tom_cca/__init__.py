"""Inter-areal communication subspaces across learning — Tom (hippocampal/cortical) cohort.

Two layers live here:

**Temporal arm (PRIMARY — carries the write-up; `HANDOFF.md`, `STATE.md` §3.0).**
Running-state CCA between area pairs at 10 ms resolution, with lag curves, IFI, frozen-axes
epoch contrasts and lagged subspaces. The method layer, in pipeline order:

    config        dataset contract (AREAS, PAIRS_TEMPORAL, FS_AREAS, paths) + TEMPORAL defaults
    dataio        loaders for THIS dataset's .mat exports (spikes, velocity, trials, learning points)
    preprocess    RunningSession: running mask, unit selection, sample cap, pair + confound (X, Y, Z)
    partial       third-area partial-out
    core          pca_fit_flat / pca_scores, cca_fit / cca_score
    lagpairs      THE within-trial lag pairer (X[t] ↔ Y[t+lag]; positive lag ⇒ X leads)
    lagged        held-out per-dim lag curves, per-dim circular-shift null, IFI (+ ifi_sides)
    lag_subspace  refit-per-lag subspaces (FF/FB), principal angles, split-half floor
    fixed_subspace frozen-axes fits, variate lag curves / moments, frozen permutation null
    perdim_ifi    per-canonical-dim IFI re-reductions of on-disk curves
    cc_aggregate  animals-as-n aggregation across CCs; paired_stats  paired t / Wilcoxon / BH-FDR
    membership, subspace, mixed_effects  shared readouts (Gini, principal angles, LMM)

Drivers: ``scripts/run_{lag_curves,lag_subspaces,lag_cosine,fixed_subspace,cc_label_track}.py``
(all on ``scripts/_common.py``); see ``audit/REPORT_2026-08-19.md`` for the consolidation.

**Spatial + landmark arms (superseded, retained for the record — `STATE.md` §2).**
``analysis``, ``pipeline``, ``stage3``, ``sweep``, ``surrogate``, ``crosspair``,
``lagged_temporal``, ``lagged_landmark``, ``landmark_align``, ``segments``, ``prune_table``,
``subspace_stats``, ``subspace_window``, ``kernel_cca``, ``kcca_window``, ``trajectory``,
``early_trials``. Originally a port of the striatum CCA pipeline
(``StriatumACC/Striatum project/cca/``).
"""

from . import (
    analysis,
    cc_aggregate,
    config,
    core,
    crosspair,
    dataio,
    fixed_subspace,
    lag_subspace,
    lagged,
    lagpairs,
    membership,
    paired_stats,
    partial,
    perdim_ifi,
    pipeline,
    preprocess,
    stage3,
    subspace,
    surrogate,
    sweep,
)

__all__ = [
    "analysis",
    "cc_aggregate",
    "config",
    "core",
    "crosspair",
    "dataio",
    "fixed_subspace",
    "lag_subspace",
    "lagged",
    "lagpairs",
    "membership",
    "paired_stats",
    "partial",
    "perdim_ifi",
    "pipeline",
    "preprocess",
    "stage3",
    "subspace",
    "surrogate",
    "sweep",
]
