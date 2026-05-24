"""Parameter sweep for the Tom-learning CCA pipeline.

The user's commitments (NOT on the sweep grid):

* spatial analysis only -- no temporal axis;
* Tom's eight hard-coded pairs;
* FS-excluded fixed;
* bin width fixed at 200 bins (Tom's native preprocessing);
* z-scoring fixed ON;
* min_units fixed at 5;
* learning point taken from Tom's ``animal_behaviour.mat`` (no LP-criterion
  sweep, no Python-side LP detection);
* learners only.

Swept axes:

    CCA type        x  AXIS_CCA   (residual / signal)
    PC-count rule   x  AXIS_KRULE (samples 15/25/40, fixed 3/5/10/20/30,
                                   variance 75/85/95 %)
    IFI lag window  x  AXIS_LAG   (5 / 10 / 20 bins  ==  +/-12.5 / 25 / 50 cm)

Total = 2 * 11 * 3 = 66 configs.

``build_sweep(name)`` returns a list of ``(tag, cfg)`` pairs.
"""

from __future__ import annotations

import dataclasses

from . import config

# --- sweep axes (edit to trim the grid) -------------------------------------
AXIS_CCA = ((True, "res"), (False, "sig"))
AXIS_KRULE = (
    ("samp15", {"k_mode": "samples", "samples_per_pc": 15}),
    ("samp25", {"k_mode": "samples", "samples_per_pc": 25}),
    ("samp40", {"k_mode": "samples", "samples_per_pc": 40}),
    ("fix03", {"k_mode": "fixed", "k_fixed": 3}),
    ("fix05", {"k_mode": "fixed", "k_fixed": 5}),
    ("fix10", {"k_mode": "fixed", "k_fixed": 10}),
    ("fix20", {"k_mode": "fixed", "k_fixed": 20}),
    ("fix30", {"k_mode": "fixed", "k_fixed": 30}),
    ("var75", {"k_mode": "variance", "k_variance": 0.75}),
    ("var85", {"k_mode": "variance", "k_variance": 0.85}),
    ("var95", {"k_mode": "variance", "k_variance": 0.95}),
)
# IFI lag-scan window in bins. 2.5 cm/bin: +/-5 = +/-12.5 cm, +/-10 = +/-25 cm,
# +/-20 = +/-50 cm. The IFI is reported at every sub-window of the scan via
# ``analysis.EpochAnalysis.ifi_windows`` (shape (n_dims, max_lag_bins)), so
# every sweep config implicitly carries all sub-windows of its scan.
AXIS_LAG = (5, 10, 20)


def _spatial() -> list[tuple[str, config.Config]]:
    """The Cartesian product of every sweep axis."""
    base = config.DEFAULT
    out: list[tuple[str, config.Config]] = []
    for resid, ctag in AXIS_CCA:
        for krtag, kover in AXIS_KRULE:
            for lag in AXIS_LAG:
                cfg = dataclasses.replace(
                    base,
                    subtract_trial_mean=resid,
                    max_lag_bins=lag,
                    **kover,
                )
                tag = f"{ctag}_{krtag}_lag{lag:02d}"
                out.append((tag, cfg))
    return out


def build_sweep(name: str) -> list[tuple[str, config.Config]]:
    """``(tag, cfg)`` pairs for sweep ``name``. ``name`` must be ``"spatial"``."""
    if name == "spatial":
        return _spatial()
    raise ValueError(f"unknown sweep: {name!r}")
