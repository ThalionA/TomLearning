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
import os
from pathlib import Path

# --- Paths -------------------------------------------------------------------
# This file lives at  TomLearning/cca/src/tom_cca/config.py
_REPO_ROOT = Path(__file__).resolve().parents[3]          # ".../TomLearning"
# Per-animal exports + cohort learning-point file. Tom's MATLAB pipeline points
# at HC_V1_data/ inside the repo; the loader accepts a custom path too.
DATA_DIR = _REPO_ROOT / "HC_V1_data"
DATA_FILE_PATTERN = "TF*_export.mat"
LEARNING_FILE = "animal_behaviour.mat"
# Plain-JSON companion file holding the region labels and learning points that
# Tom's exports store as MATLAB `string` arrays (a type h5py cannot decode).
# Generated once by scripts/export_cca_labels.m; consumed by dataio. When this
# file is present alongside the .mat exports it is the authoritative source for
# region labels and LPs; when absent the loader falls back to reading those
# fields directly from the .mat files (works for cellstr-typed test fixtures).
COMPANION_FILE = "cca_labels.json"
CCA_DIR = _REPO_ROOT / "cca"
RESULTS_DIR = CCA_DIR / "results"
FIGURES_DIR = CCA_DIR / "figures"

# Temporal-arm figures are written to FIGURES_DIR and MIRRORED (byte-identical copy)
# into the research vault so Obsidian notes can embed them. Set the environment
# variable CCA_FIG_MIRROR to another directory, or to "" to disable the mirror —
# e.g. on a machine without the vault. (Before 2026-08-19 the vault path was
# hard-coded in 24 figure scripts, 7 of them as /Users/theoamvr/...)
_mirror_env = os.environ.get("CCA_FIG_MIRROR")
FIGURE_MIRROR_DIR: Path | None = (
    (Path(_mirror_env).expanduser() if _mirror_env else None) if _mirror_env is not None
    else Path.home() / "Documents" / "ResearchVault" / "attachments"
)
# Figure filename prefix for this dataset (``HCV1_<stem>_<fs>_bin10.svg``).
FIG_PREFIX = "HCV1"

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

# The SAME eight pairs in the panel/table order every temporal-arm script has used
# since 2026-06 (figures, md tables, CSV iteration). ``PAIRS`` above keeps the MATLAB
# order for the spatial arm. Before 2026-08-19 this list was re-declared in 51 scripts.
PAIRS_TEMPORAL: tuple[tuple[str, str], ...] = (
    ("CA1", "RSC"), ("CA1", "CA3"), ("CA1", "DG"), ("CA1", "V1"),
    ("CA3", "DG"), ("CA1", "SUB"), ("RSC", "SUB"), ("V1", "RSC"),
)
PAIR_NAMES: tuple[str, ...] = tuple(f"{x}-{y}" for x, y in PAIRS_TEMPORAL)

# --- Surrogate / significance ------------------------------------------------
# Shuffle count for the circular-shift surrogate used to test canonical-dimension
# significance (95th-percentile threshold). Centralised so all continuous-regime
# drivers (run_trajectory / run_transition / run_epochs) stay consistent; 100
# matches the reference papers (20 is too few for a stable 95th percentile).
SURROGATE_SHUFFLES: int = 100

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
    """Tunable analysis parameters. Frozen so a run's settings are immutable.

    Committed-by-user-instruction (NOT on the sweep grid):
      * ``zscore_units = True`` -- whole-engaged-period per-unit z-scoring
      * ``min_units = 5``       -- minimum units per area per pair
      * ``exclude_fast_spiking = True`` -- FS-excluded (Tom convention)
      * spatial bin width: 200 bins x 2.5 cm (no rebin)
      * learning point: whatever ``animal_behaviour.mat`` provides
        (no Python-side LP detection; no LP-criterion sweep)
    """

    # Animals forced to non-learner regardless of the recorded LP. Empty for
    # Tom -- if a known LP artefact surfaces, add the animal_id here.
    manual_nonlearners: frozenset[int] = frozenset()

    # Epochs.
    trials_per_epoch: int = 10

    # Unit inclusion. FS-excluded is FIXED (user instruction). Override the
    # `exclude_fast_spiking` flag explicitly to run the FS-included variant
    # (matching striatum's `--include-fs` switch).
    min_units: int = 5
    exclude_fast_spiking: bool = True

    # Residualisation: subtract the per-(bin, unit) trial mean. Default True;
    # signal-CCA variant flips this off (the one remaining sweep axis on the
    # "what is communication" question).
    subtract_trial_mean: bool = True
    # Z-score each unit by its std over the entire engaged period, applied to
    # the raw activity *before* residualisation and epoch slicing. FIXED ON
    # (user instruction); not on the sweep grid.
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

    # ---- Temporal arms (added 2026-05-26) ----------------------------------
    # bin_mode selects which arm to run. "spatial" preserves the current spatial
    # CCA path byte-identically. "temporal_runstate" is Arm A (running-state
    # full-traversal, signal-only, segment-fragmented). "landmark" is Arm B
    # (landmark-centred per-landmark, residual + signal). See
    # cca/UNDERSTANDING_temporal.md for the full design.
    bin_mode: str = "spatial"
    # Locked at 50 ms (Tom MATLAB-v6 alignment). Not on the temporal sweep grid.
    temporal_bin_ms: int = 50
    # Gaussian s.d. (ms) for convolving the 1 ms spike train before binning, to a
    # smoothed spike-density. 0 = raw summed counts (legacy). Gonzalez & Buzsáki
    # 2026 use 2.5 ms; set via the drivers' --smooth-ms flag.
    gaussian_sd_ms: float = 0.0
    # Lag ranges expressed in ms so the physical directionality window stays
    # constant if bin width ever changes.
    lag_ms: int = 500              # Arm A: +-500 ms = +-10 bins at 50 ms
    landmark_lag_ms: int = 250     # Arm B: +-250 ms forced by 500 ms window
    # Velocity gate: applied to the 50 ms-binned mean velocity. Threshold and
    # morphological cleanup parameters; the min surviving run length is
    # lag_ms + min_run_extra_ms so every segment yields >= 5 paired samples at
    # the most extreme lag.
    velocity_thresh_cm_s: float = 2.0
    min_gap_close_ms: int = 100    # close invalid runs <= this many ms
    min_run_extra_ms: int = 250    # min run = lag_ms + min_run_extra_ms
    # Arm B landmark window: 5 bins pre + entry bin + 4 bins post = 10 bins
    # spanning -250 ms to +200 ms relative to entry.
    # Window specified in ms so bin-width changes (e.g. 25 ms) don't silently
    # change the physical extent of the analysis window. Defaults span
    # -250 ms to +200 ms relative to entry -- 10 bins at 50 ms, 20 bins at 25 ms.
    landmark_window_pre_ms: int = 250
    landmark_window_post_ms: int = 200
    min_crossings_per_landmark: int = 5
    # Per-segment circshift floor (Arm A): minimum shift in bins so the null
    # beats neural autocorrelation. 5 bins = 250 ms, larger than the typical
    # cortex/hippocampus autocorrelation timescale.
    temporal_circshift_min_bins: int = 5
    # Per-window circshift floor (Arm B): a 1-bin shift on a 10-bin window
    # leaves the rolled version highly correlated with the original (only one
    # bin moved), so the floor should be a meaningful fraction of the window.
    # 3 bins = 150 ms, large enough to break within-window auto-structure
    # while leaving 7 valid shift values (3..9) for the random draw.
    landmark_circshift_min_bins: int = 3


DEFAULT = Config()


@dataclass(frozen=True)
class TemporalDefaults:
    """The canonical running-state (temporal-arm) configuration — ONE object.

    These are the values every live temporal driver ran with (HANDOFF.md §2); until
    2026-08-19 each script re-declared them as module constants with three different
    ``--bin-ms`` defaults and two shuffle counts. ``Config`` above still carries the
    older per-field knobs the spatial/landmark arms read (``temporal_bin_ms = 50``
    is that legacy default, NOT the temporal arm's).
    """

    bin_ms: int = 10            # Gaussian-smoothed spike density, 10 ms bins
    smooth_ms: float = 2.5      # sigma of the 1 ms smoothing kernel (Gonzalez & Buzsáki)
    k: int = 30                 # PCs per area before CCA
    n_folds: int = 5            # CV folds, split by WHOLE trials
    n_shuffles: int = 200       # circular-shift surrogates (p floor 1/201 < 0.05/10)
    fdr_dims: int = 10          # BH family = the leading dims (p-floor arithmetic)
    max_lag_bins: int = 25      # +/-250 ms at 10 ms
    label_w_bins: int = 5       # FF/FB label window: IFI over |lag| <= 5 bins (+/-50 ms)
    min_trials: int = 6         # an animal needs > n_folds running trials
    velocity_thresh_cm_s: float = 2.0
    min_units: int = 5


TEMPORAL = TemporalDefaults()


# ---- Derived temporal quantities -------------------------------------------
# These live as module-level helpers (not @property) to keep Config a pure
# frozen data record.
def lag_bins(cfg: Config) -> int:
    """Arm A lag range in 50 ms bins (one-sided; full scan is [-lag, +lag])."""
    return cfg.lag_ms // cfg.temporal_bin_ms


def landmark_lag_bins(cfg: Config) -> int:
    """Arm B lag range in 50 ms bins (one-sided)."""
    return cfg.landmark_lag_ms // cfg.temporal_bin_ms


def min_run_bins(cfg: Config) -> int:
    """Arm A: minimum length of a retained velocity-passing segment, in bins.

    Sized so that even at +-lag_bins, at least 5 paired samples survive.
    """
    return (cfg.lag_ms + cfg.min_run_extra_ms) // cfg.temporal_bin_ms


def min_gap_close_bins(cfg: Config) -> int:
    """Arm A: invalid runs of length <= this many bins are bridged."""
    return cfg.min_gap_close_ms // cfg.temporal_bin_ms


def landmark_window_bins_pre(cfg: Config) -> int:
    """Arm B: number of bins before the entry bin."""
    return cfg.landmark_window_pre_ms // cfg.temporal_bin_ms


def landmark_window_bins_post(cfg: Config) -> int:
    """Arm B: number of bins after the entry bin."""
    return cfg.landmark_window_post_ms // cfg.temporal_bin_ms


def landmark_window_bins(cfg: Config) -> int:
    """Arm B: total window length in bins (pre + entry-bin + post)."""
    return landmark_window_bins_pre(cfg) + 1 + landmark_window_bins_post(cfg)


def landmark_entry_index(cfg: Config) -> int:
    """Arm B: index within the window where the entry bin lives."""
    return landmark_window_bins_pre(cfg)
