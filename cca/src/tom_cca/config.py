"""Configuration for the Tom-learning CCA pipeline.

Single source of truth for paths, area definitions and the analysis knobs --
the Tom-side adaptation of the striatum_cca package. The major data-level
differences from the striatum pipeline are:

* one .mat file per animal (TF*_export.mat) instead of one .mat for the cohort;
* spatial firing rate lives in ``analysis_spatial.firing.cued.freq`` --
  ``(n_units, n_trials, n_bins)``, 200 bins x 2.5 cm = 500 cm corridor;
* learning points come from a separate ``animal_behaviour.mat``
  (``period_experienced(:, 1)`` indexed by ``animal_id``);
* areas are CA1/V1/DG/CA3/RSC/SUB; pair list is the eight Tom hard-codes
  in HC_V1_Code/HC_V1_temporal.m and HC_V1_Code/legacy/CCA_HC_V1_spatial_v2.m;
* FS exclusion uses ``units.idx_fs`` and is applied to V1/RSC/CA1/CA3 only --
  the same convention Tom's MATLAB pipeline uses (no FS flags in DG/SUB).

The pipeline always excludes FS units (the user's instruction); this is fixed,
not swept.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# --- Paths -------------------------------------------------------------------
# This file lives at  TomLearning/cca/src/tom_cca/config.py
_REPO_ROOT = Path(__file__).resolve().parents[3]          # ".../TomLearning"
# Per-animal exports + cohort learning-point file. Tom's MATLAB pipeline points
# at HC_V1_data/ inside the repo; the loader accepts a custom path too.
DATA_DIR = _REPO_ROOT / "HC_V1_data"
DATA_FILE_PATTERN = "TF*_export.mat"
LEARNING_FILE = "animal_behaviour.mat"
CCA_DIR = _REPO_ROOT / "cca"
RESULTS_DIR = CCA_DIR / "results"
FIGURES_DIR = CCA_DIR / "figures"

# --- Areas and pairs ---------------------------------------------------------
# Six areas appear in the Tom cohort. The first element of each pair tuple is
# "X" -- it fixes the sign convention for the Information Flow Index
# (positive IFI => X leads Y), matching how the striatum pipeline reports IFI.
# The pair list is byte-identical to the hard-coded list in
# HC_V1_Code/HC_V1_temporal.m and the spatial CCA scripts.
AREAS: tuple[str, ...] = ("CA1", "V1", "DG", "CA3", "RSC", "SUB")
PAIRS: tuple[tuple[str, str], ...] = (
    ("CA1", "V1"),
    ("CA1", "DG"),
    ("CA1", "CA3"),
    ("CA1", "RSC"),
    ("CA1", "SUB"),
    ("V1", "RSC"),
    ("RSC", "SUB"),
    ("CA3", "DG"),
)

# --- FS exclusion convention (Tom MATLAB) ------------------------------------
# Tom's pipeline only flags fast-spiking units in V1, RSC, CA1, CA3 -- the
# four areas with the histology / spike-width discriminability for it. DG and
# SUB have no FS flags. The Tom MATLAB scripts apply the FS mask only to those
# four areas; we copy that convention verbatim.
FS_AREAS: tuple[str, ...] = ("V1", "RSC", "CA1", "CA3")

# --- Spatial geometry --------------------------------------------------------
# Tom's track is 500 cm long, binned into 200 spatial bins at 2.5 cm each
# (see HC_V1_Code/legacy/CCA_HC_V1_spatial_v2.m).
N_BINS = 200
CORRIDOR_CM = 500.0


def bin_size_cm(n_bins: int) -> float:
    """Spatial bin width in cm for a given bin count."""
    return CORRIDOR_CM / n_bins


# --- Epoch names (ordered naive -> intermediate -> expert) -------------------
# Three-epoch design, matching the striatum pipeline post-round-10
# (naive = trials 1-10, intermediate = the 10 trials ending at LP,
# expert = the 10 trials after LP).
EPOCH_NAMES: tuple[str, str, str] = ("naive", "intermediate", "expert")

# Epoch plotting colours, kept byte-identical to the striatum pipeline so the
# Tom and striatum figures are visually consistent.
EPOCH_COLOURS: dict[str, tuple[float, float, float]] = {
    "naive": (0.298, 0.447, 0.690),          # blue
    "intermediate": (0.867, 0.518, 0.322),   # orange
    "expert": (0.333, 0.776, 0.333),         # green
}


@dataclass(frozen=True)
class Config:
    """Tunable analysis parameters. Frozen so a run's settings are immutable."""

    # Learning point detection (project rule).
    lp_z_threshold: float = -2.0
    lp_window: int = 10
    lp_min_consecutive: int = 7
    # Animals forced to non-learner regardless of a detected LP. Empty for Tom
    # by default -- there is no analogue of the striatum animal-8 artefact.
    manual_nonlearners: frozenset[int] = frozenset()

    # Epochs.
    trials_per_epoch: int = 10

    # Unit inclusion. FS-excluded is FIXED here (the user's request) -- it is
    # not on the sweep grid. The flag stays so the FS-included variant can be
    # run by overriding the config explicitly (matching striatum's
    # `--include-fs` switch).
    min_units: int = 6
    exclude_fast_spiking: bool = True

    # Residualisation (D2 in the striatum spec): subtract the per-(bin, unit)
    # trial mean. Default True; signal-CCA variant flips this off.
    subtract_trial_mean: bool = True
    # Z-score each unit by its std over the entire engaged period, applied to
    # the raw activity *before* residualisation and epoch slicing (striatum
    # round 7 convention).
    zscore_units: bool = True

    # Source field in analysis_spatial.firing.cued. "freq" = raw firing rate
    # (Hz) -- the natural input that lets the pipeline apply its own z-score.
    # "freq_z" = the file's pre-z-scored variant; selecting it makes the
    # pipeline's `zscore_units` redundant but is provided for symmetry with
    # Tom's MATLAB scripts.
    spatial_field: str = "freq"

    # PCA / k rule (D4). k_mode selects how many PCs to keep per area:
    #   "samples"  -- k = floor(n_samples / samples_per_pc)   [default]
    #   "fixed"    -- k = k_fixed
    #   "variance" -- k = #PCs reaching k_variance cumulative variance
    # All are capped by the smaller area's unit count, k_cap and the per-epoch
    # numerical rank.
    k_mode: str = "samples"
    k_fixed: int = 10
    k_variance: float = 0.90
    samples_per_pc: int = 15
    k_cap: int = 30

    # Cross-validation: 5-fold over whole trials.
    n_folds: int = 5
    cv_seed: int = 0

    # Lagged CCA / directionality. 10 bins at 2.5 cm = +/-25 cm -- matches the
    # striatum pipeline's spatial window.
    max_lag_bins: int = 10

    # Surrogate null: per-dimension held-out-CC permutation test.
    # null_type "circshift" -- per-trial circular shift of the bin axis by
    # >= circshift_min_bins (Gonzalez et al.; committed default in the
    # striatum pipeline). "trials" -- permute trial correspondence (H&H).
    n_shuffles: int = 200
    surrogate_seed: int = 0
    null_type: str = "circshift"
    circshift_min_bins: int = 15

    # Parallelism: processes for the cohort run.
    n_jobs: int = 4


DEFAULT = Config()
