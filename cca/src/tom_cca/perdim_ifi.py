"""Per-canonical-dimension Information Flow Index.

Meeting 2026-07-28 (Nathalie/Tom), item 1: *do different CCs have different IFI?*

The shipped directionality readout (`run_ifi_windows`, STATE.md §3.0 finding 3) summarises
only the dominant canonical dimension. But `run_lag_curves` already exports the held-out,
segment-aware lag curve for **every** canonical dimension, so the question is answerable by
re-reducing data that is already on disk — no refits.

**Sign convention.** Positive lag => the FIRST-named area leads (`lagged.lag_slice`), so
positive IFI => X leads Y. Inherited unchanged from `lagged.information_flow_index`.

**Interpretive caveat — read before using `dim` as an identity.**
`lagged.heldout_lag_curve_flat_perdim` refits CCA independently at every lag, so `dim = 2`
at lag -100 ms and `dim = 2` at lag +100 ms are the second-ranked canonical dimension *of
two different fits*, not a tracked subspace. Everything here is therefore a statement about
**canonical rank**, not about a persistent component. Whether a stronger reading is available
is exactly what the cross-lag subspace-stability analysis (meeting item 3) tests.
"""

from __future__ import annotations

import numpy as np

from . import lagged


def curve_ifi(lag: np.ndarray, cc: np.ndarray) -> float:
    """IFI of one (lag, cc) curve, tolerant of unsorted rows and NaN lags.

    Sorts by lag (CSV row order is not guaranteed), drops lags whose held-out CC is
    non-finite (too few segment-aware pairs for any fold), then delegates to
    :func:`lagged.information_flow_index`. Returns NaN if nothing finite survives on
    either side of lag 0 — distinct from the 0.0 that function returns for a curve that
    is finite but has no asymmetry.
    """
    lag = np.asarray(lag, dtype=float)
    cc = np.asarray(cc, dtype=float)
    keep = np.isfinite(lag) & np.isfinite(cc)
    if not np.any(keep):
        return float("nan")
    lag, cc = lag[keep], cc[keep]
    order = np.argsort(lag, kind="stable")
    lag, cc = lag[order], cc[order]
    if not (np.any(lag > 0) and np.any(lag < 0)):
        return float("nan")
    return float(lagged.information_flow_index(lag, cc))


def curve_peak_lag(lag: np.ndarray, cc: np.ndarray) -> float:
    """Lag at which the held-out CC peaks. NaN if the curve is entirely non-finite."""
    lag = np.asarray(lag, dtype=float)
    cc = np.asarray(cc, dtype=float)
    keep = np.isfinite(lag) & np.isfinite(cc)
    if not np.any(keep):
        return float("nan")
    return float(lag[keep][int(np.nanargmax(cc[keep]))])


def ifi_by_dim(lag: np.ndarray, cc: np.ndarray, dim: np.ndarray):
    """Reduce a long-format ``(lag, cc, dim)`` block to one IFI per canonical rank.

    Returns ``(dims, ifi, peak_lag)`` with ``dims`` ascending. Intended to be applied
    within one (animal, pair) group; grouping is the caller's job.
    """
    lag = np.asarray(lag, dtype=float)
    cc = np.asarray(cc, dtype=float)
    dim = np.asarray(dim)
    dims = np.unique(dim[np.isfinite(np.asarray(dim, dtype=float))])
    ifi = np.full(dims.size, np.nan)
    peak = np.full(dims.size, np.nan)
    for i, d in enumerate(dims):
        m = dim == d
        ifi[i] = curve_ifi(lag[m], cc[m])
        peak[i] = curve_peak_lag(lag[m], cc[m])
    return dims, ifi, peak


def sign_agreement(ifi: np.ndarray) -> float:
    """Fraction of the non-dominant dims whose IFI shares CC1's sign.

    The dims-differ question in its crudest form: 1.0 = every canonical dimension points
    the same way as CC1, 0.5 = chance, 0.0 = the spectrum systematically reverses.
    ``ifi[0]`` is the CC1 reference. Dims with an exactly-zero IFI are excluded — that is
    :func:`lagged.information_flow_index`'s "no asymmetry" return, not evidence of a
    direction. NaN if the reference is non-finite or no comparable dim survives.
    """
    ifi = np.asarray(ifi, dtype=float)
    if ifi.size < 2 or not np.isfinite(ifi[0]) or ifi[0] == 0:
        return float("nan")
    tail = ifi[1:]
    tail = tail[np.isfinite(tail) & (tail != 0)]
    if tail.size == 0:
        return float("nan")
    return float(np.mean(np.sign(tail) == np.sign(ifi[0])))


def rank_split(ifi: np.ndarray, sig: np.ndarray) -> tuple[float, float]:
    """``(ifi_cc1, ifi_rest)`` — the dominant dim vs the mean of the significant tail.

    CC1 is the reference *by construction* and is returned whatever its own significance
    flag says, so the contrast can never be silently redefined onto a different dimension.
    The tail averages only dims flagged significant against the circular-shift null; NaN
    when no significant dim beyond the first survives.
    """
    ifi = np.asarray(ifi, dtype=float)
    sig = np.asarray(sig)
    if ifi.size == 0:
        return float("nan"), float("nan")
    cc1 = float(ifi[0])
    tail = ifi[1:][(sig[1:] == 1) & np.isfinite(ifi[1:])]
    rest = float(np.mean(tail)) if tail.size else float("nan")
    return cc1, rest


def ifi_windows_by_dim(lag, cc, dim, max_window: int | None = None):
    """IFI of EACH canonical dimension as a function of integration window.

    Meeting item 1 as asked: *"calculate IFI at different lags for each CC, for each
    area, and then look at them — are there CCs that have different signs?"*

    ``lag`` must be in BIN units (not ms): the window grid is one entry per integer
    ``w = 1 .. max|lag|``, and entry ``w`` integrates only lags with ``|lag| <= w``
    (:func:`lagged.ifi_by_window`). Passing ms would create one window per millisecond,
    almost all of them duplicates.

    Returns ``(dims, windows, ifi, degenerate)``:
      * ``dims``       ascending canonical ranks present
      * ``windows``    ``1 .. max_window`` in BINS
      * ``ifi``        ``(n_dims, n_windows)``
      * ``degenerate`` ``(n_dims, n_windows)`` bool — True where the IFI is 0.0 only
        because BOTH sides clipped to zero. :func:`lagged.information_flow_index`
        returns 0.0 for "no coupling either way" and for "genuinely balanced", and those
        must not be pooled: a held-out CC that is negative at every lag is noise, not a
        symmetric flow. High canonical ranks sit at the CC floor, so this is common.
    """
    lag = np.asarray(lag, dtype=float)
    cc = np.asarray(cc, dtype=float)
    dim = np.asarray(dim)
    finite_dim = np.asarray(dim, dtype=float)
    dims = np.unique(dim[np.isfinite(finite_dim)])
    max_w = int(max_window if max_window is not None
                else (np.nanmax(np.abs(lag)) if lag.size else 0))
    windows = np.arange(1, max_w + 1)
    ifi = np.full((dims.size, windows.size), np.nan)
    degen = np.zeros((dims.size, windows.size), dtype=bool)
    for i, d in enumerate(dims):
        m = (dim == d) & np.isfinite(lag) & np.isfinite(cc)
        if not np.any(m):
            continue
        lg, c = lag[m], cc[m]
        order = np.argsort(lg, kind="stable")
        lg, c = lg[order], c[order]
        for j, w in enumerate(windows):
            sel = np.abs(lg) <= w
            if not (np.any(lg[sel] > 0) and np.any(lg[sel] < 0)):
                continue
            ifi[i, j] = lagged.information_flow_index(lg[sel], c[sel])
            pm, nm = lagged.ifi_sides(lg[sel], c[sel])
            degen[i, j] = (pm + nm) <= 0
    return dims, windows, ifi, degen
