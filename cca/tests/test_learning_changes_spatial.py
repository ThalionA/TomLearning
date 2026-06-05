"""Tests for the spatial-arm paired learning-change machinery.

Covers two new pieces:
* ``subspace_stats.epoch_subspace_stats`` -- spatial analogue of
  ``cell_subspace_stats`` that reads the lag-0 ``held_out_cc`` + ``p_per_dim``
  (the spatial arm's native significance) instead of the peak-across-lag
  quantities used by the landmark/temporal cells.
* ``paired_stats`` -- the per-animal Wilcoxon signed-rank + Benjamini-Hochberg
  FDR helpers used to test naive -> expert changes.
"""

from __future__ import annotations

import numpy as np

from tom_cca import paired_stats, subspace_stats


class _MockEpoch:
    """Duck-typed spatial EpochAnalysis (only the read fields)."""

    def __init__(self, held_out_cc, p_per_dim, ifi_per_dim):
        self.held_out_cc = np.asarray(held_out_cc, dtype=float)
        self.p_per_dim = np.asarray(p_per_dim, dtype=float)
        self.ifi_per_dim = np.asarray(ifi_per_dim, dtype=float)


class _MockCell:
    """Duck-typed landmark cell for the cross-check."""

    def __init__(self, lag_cc_per_dim, p_peak_cc_per_dim, ifi_per_dim):
        self.lag_cc_per_dim = np.asarray(lag_cc_per_dim, dtype=float)
        self.p_peak_cc_per_dim = np.asarray(p_peak_cc_per_dim, dtype=float)
        self.ifi_per_dim = np.asarray(ifi_per_dim, dtype=float)
        self.skip_reason = None


# --------------------------------------------------------------------------
# epoch_subspace_stats
# --------------------------------------------------------------------------

def test_epoch_stats_n_sig_and_mi_sig():
    ep = _MockEpoch(held_out_cc=[0.6, 0.3, 0.05],
                    p_per_dim=[0.001, 0.2, 0.9],
                    ifi_per_dim=[0.5, -0.2, 0.1])
    s = subspace_stats.epoch_subspace_stats(ep, alpha=0.05)
    assert s.n_sig == 1
    assert list(s.sig_mask) == [True, False, False]
    # Only dim 0 is significant: MI = -0.5 * log(1 - 0.6**2).
    expected_mi_sig = -0.5 * np.log1p(-0.6 ** 2)
    assert abs(s.mi_sig - expected_mi_sig) < 1e-9


def test_epoch_stats_ifi_weighted_uses_all_finite_dims():
    ep = _MockEpoch(held_out_cc=[0.6, 0.3, 0.05],
                    p_per_dim=[0.001, 0.2, 0.9],
                    ifi_per_dim=[0.5, -0.2, 0.1])
    s = subspace_stats.epoch_subspace_stats(ep, alpha=0.05)
    # CC-weighted over ALL finite dims, not just significant ones.
    num = 0.5 * 0.6 + (-0.2) * 0.3 + 0.1 * 0.05
    den = 0.6 + 0.3 + 0.05
    assert abs(s.ifi_weighted - num / den) < 1e-9


def test_epoch_stats_mi_sig_zero_when_none_significant():
    ep = _MockEpoch(held_out_cc=[0.5, 0.5],
                    p_per_dim=[0.5, 0.9],
                    ifi_per_dim=[0.1, 0.2])
    s = subspace_stats.epoch_subspace_stats(ep, alpha=0.05)
    assert s.n_sig == 0
    assert s.mi_sig == 0.0
    assert s.mi_all > 0


def test_epoch_stats_none_when_no_cc():
    class _Empty:
        held_out_cc = None
    assert subspace_stats.epoch_subspace_stats(_Empty()) is None


def test_epoch_stats_matches_cell_stats_on_equivalent_inputs():
    # When a landmark cell's peak-across-lag CC equals an epoch's lag-0 CC and
    # the p-values match, both helpers must agree (shared aggregation core).
    cc = [0.7, 0.4, 0.2, 0.05]
    p = [0.001, 0.03, 0.2, 0.8]
    ifi = [0.3, -0.1, 0.2, 0.0]
    # lag_cc with the peak (over lags) equal to `cc` on the centre row.
    lag_cc = np.tile(np.array(cc) * 0.5, (5, 1))
    lag_cc[2] = cc                                   # peak row == cc
    ep = _MockEpoch(cc, p, ifi)
    cell = _MockCell(lag_cc, p, ifi)
    se = subspace_stats.epoch_subspace_stats(ep, alpha=0.05)
    sc = subspace_stats.cell_subspace_stats(cell, alpha=0.05)
    assert se.n_sig == sc.n_sig
    assert abs(se.mi_sig - sc.mi_sig) < 1e-12
    assert abs(se.ifi_weighted - sc.ifi_weighted) < 1e-12


# --------------------------------------------------------------------------
# paired_stats
# --------------------------------------------------------------------------

def test_fdr_bh_known_pvalues():
    p = np.array([0.001, 0.04, 0.2, 0.5])
    mask = paired_stats.fdr_bh(p, q=0.05)
    assert list(mask) == [True, False, False, False]


def test_fdr_bh_all_pass_when_tiny():
    p = np.array([0.0001, 0.0002, 0.0003])
    mask = paired_stats.fdr_bh(p, q=0.05)
    assert list(mask) == [True, True, True]


def test_fdr_bh_handles_nan():
    p = np.array([0.001, np.nan, 0.9])
    mask = paired_stats.fdr_bh(p, q=0.05)
    assert mask[0] and not mask[2]
    assert not mask[1]


def test_wilcoxon_signed_detects_positive_shift():
    n, med, w, p = paired_stats.wilcoxon_signed([0.4, 0.5, 0.6, 0.7, 0.8])
    assert n == 5
    assert abs(med - 0.6) < 1e-9
    assert np.isfinite(p) and p < 0.1


def test_wilcoxon_signed_too_few():
    n, med, w, p = paired_stats.wilcoxon_signed([0.5, 0.6])
    assert n == 2
    assert np.isnan(p)


def test_wilcoxon_signed_all_zero():
    n, med, w, p = paired_stats.wilcoxon_signed([0.0, 0.0, 0.0, 0.0])
    assert np.isnan(p)
