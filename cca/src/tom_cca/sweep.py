"""Parameter sweep for the Tom-learning CCA pipeline.

The user's brief:

* spatial analysis only -- no temporal axis;
* Tom's eight hard-coded pairs;
* FS-excluded is fixed (not on the sweep grid);
* bin width is fixed at 200 bins (Tom's native preprocessing -- no rebin);
* sweep every other knob -- residual/signal CCA, z-scoring on/off, min-units,
  LP criterion, PC-count rule.

That leaves the following axes (compared to the striatum sweep we have
dropped ``AXIS_BINS`` and ``AXIS_FS``):

    CCA type        x  AXIS_CCA       (residual / signal)
    z-scoring       x  AXIS_Z         (on / off)
    min units       x  AXIS_MIN_UNITS (4 / 6 / 10)
    LP criterion    x  AXIS_LP        (7 / 8 consecutive)
    PC-count rule   x  AXIS_KRULE     (samples 15/25/40, fixed 3/5/10/20/30,
                                       variance 75/85/95 %)

Total = 2 * 2 * 3 * 2 * 11 = 264 configs.

``build_sweep(name)`` returns a list of ``(tag, cfg)`` pairs. Defined in one
place so run_stage2.py, run_stage3.py and summarise_sweep.py all see the same
grid.
"""

from __future__ import annotations

import dataclasses

from . import config


# --- sweep axes (edit to trim the grid) -------------------------------------
AXIS_CCA = ((True, "res"), (False, "sig"))
AXIS_Z = ((True, "z1"), (False, "z0"))
AXIS_MIN_UNITS = (4, 6, 10)
AXIS_LP = (7, 8)
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


def _spatial() -> list[tuple[str, config.Config]]:
    """The Cartesian product of every sweep axis. FS is fixed-excluded."""
    base = config.DEFAULT
    out: list[tuple[str, config.Config]] = []
    for resid, ctag in AXIS_CCA:
        for z, ztag in AXIS_Z:
            for mu in AXIS_MIN_UNITS:
                for lpc in AXIS_LP:
                    for krtag, kover in AXIS_KRULE:
                        cfg = dataclasses.replace(
                            base, subtract_trial_mean=resid,
                            exclude_fast_spiking=True,    # fixed
                            zscore_units=z, min_units=mu,
                            lp_min_consecutive=lpc, **kover)
                        tag = (f"{ctag}_{ztag}_mu{mu:02d}"
                               f"_lp{lpc}_{krtag}")
                        out.append((tag, cfg))
    return out


def build_sweep(name: str) -> list[tuple[str, config.Config]]:
    """``(tag, cfg)`` pairs for sweep ``name``.

    ``name`` must be ``"spatial"`` -- there is no temporal arm for Tom.
    """
    if name == "spatial":
        return _spatial()
    raise ValueError(f"unknown sweep: {name!r}")
