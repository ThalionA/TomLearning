"""
Python prototype + synthetic-data tests for the v4 temporal CCA algorithm.

This is a *shadow implementation* used to validate the algorithm before porting
to MATLAB. It is not part of the MATLAB pipeline; the MATLAB port lives in
HC_V1_Code/ (CCA_HC_V1_temporal_v4.m + v4_*.m helpers + tests/test_temporal_v4.m).

Synthetic ground truths exercised:
  T1 - Lagged AR(1) coupling: v lags u by 1 bin -> per-trial IFI is
       significantly negative (u leads v -> r_pos > r_neg -> IFI < 0).
  T2 - No coupling: per-trial CC is not significantly above shuffle CC.
  T3 - RZ-locked coupling: bursty cross-region coupling only in the last
       250 ms before RZ entry -> RZ-aligned CC peaks in the pre-RZ bins,
       and RZ-aligned per-trial IFI sign matches the lag of the burst.
  T4 - Block segmentation: a trial with two valid blocks separated by a
       velocity drop yields a weighted IFI that respects block lengths.

Run:  python3 tests/test_temporal_v4_proto.py
"""

from __future__ import annotations

import math
import sys
import numpy as np
from numpy.linalg import LinAlgError
from scipy.linalg import qr

RNG_GLOBAL = np.random.default_rng(0)


# -----------------------------------------------------------------------------
# Core helpers (Python mirror of the MATLAB helpers)
# -----------------------------------------------------------------------------

def canoncorr_first(X: np.ndarray, Y: np.ndarray):
    """First canonical pair via QR + SVD (matches MATLAB's canoncorr semantics).

    Returns A (k1, 1), B (k2, 1), r (scalar). Raises if rank-deficient.
    """
    n, k1 = X.shape
    _, k2 = Y.shape
    if n <= k1 + k2:
        raise ValueError("n must exceed k1 + k2 for stable canoncorr")
    Xc = X - X.mean(axis=0, keepdims=True)
    Yc = Y - Y.mean(axis=0, keepdims=True)
    Qx, Rx = qr(Xc, mode="economic")
    Qy, Ry = qr(Yc, mode="economic")
    if np.linalg.matrix_rank(Rx) < k1 or np.linalg.matrix_rank(Ry) < k2:
        raise LinAlgError("rank deficient")
    M = Qx.T @ Qy
    U, s, Vt = np.linalg.svd(M, full_matrices=False)
    r = float(s[0])
    a = np.linalg.solve(Rx, U[:, 0]) * math.sqrt(n - 1)
    b = np.linalg.solve(Ry, Vt.T[:, 0]) * math.sqrt(n - 1)
    return a.reshape(-1, 1), b.reshape(-1, 1), r


def lag_corr_simple(u: np.ndarray, v: np.ndarray, lags) -> np.ndarray:
    """Lag correlation; lag L>=0 -> corr(u[:n-L], v[L:n]) (u leads v at +L)."""
    u = np.asarray(u).ravel()
    v = np.asarray(v).ravel()
    n = len(u)
    out = np.full(len(lags), np.nan)
    for il, L in enumerate(lags):
        if L >= 0:
            a, b = u[: n - L], v[L:]
        else:
            a, b = u[-L:], v[: n + L]
        if len(a) > 5 and np.std(a) > 0 and np.std(b) > 0:
            out[il] = np.corrcoef(a, b)[0, 1]
    return out


def ifi_from_lags(r_lags: np.ndarray, lags) -> float:
    """IFI = (r_neg - r_pos) / (r_neg + r_pos). u leads v -> r_pos large -> IFI<0."""
    lags = np.asarray(lags)
    r_neg = np.nanmean(r_lags[lags < 0])
    r_pos = np.nanmean(r_lags[lags > 0])
    if (r_neg + r_pos) > 0.001:
        return float((r_neg - r_pos) / (r_neg + r_pos))
    return float("nan")


def contiguous_blocks(is_valid: np.ndarray, min_block: int):
    """Return list of (start, stop_exclusive) for runs of True with length >= min_block."""
    is_valid = np.asarray(is_valid, dtype=bool)
    blocks = []
    i = 0
    n = len(is_valid)
    while i < n:
        if not is_valid[i]:
            i += 1
            continue
        j = i
        while j < n and is_valid[j]:
            j += 1
        if (j - i) >= min_block:
            blocks.append((i, j))
        i = j
    return blocks


def per_trial_cca(X_valid: np.ndarray, Y_valid: np.ndarray):
    """Single canoncorr fit on all valid bins of one trial. Returns (r, A, B) or None."""
    if X_valid.shape[0] < X_valid.shape[1] + Y_valid.shape[1] + 5:
        return None
    try:
        A, B, r = canoncorr_first(X_valid, Y_valid)
    except (LinAlgError, ValueError):
        return None
    return r, A, B


def per_trial_ifi_weighted(
    X_trial: np.ndarray,
    Y_trial: np.ndarray,
    is_valid: np.ndarray,
    A: np.ndarray,
    B: np.ndarray,
    lags,
    min_block_bins: int,
):
    """Per-block IFI on canonical projections, weighted average by block length.

    Returns (ifi_trial, info). ifi_trial is NaN if no eligible block.
    """
    Xc = X_trial - X_trial.mean(axis=0, keepdims=True)
    Yc = Y_trial - Y_trial.mean(axis=0, keepdims=True)
    u_full = (Xc @ A).ravel()
    v_full = (Yc @ B).ravel()

    blocks = contiguous_blocks(is_valid, min_block_bins)
    if not blocks:
        return float("nan"), {"n_blocks": 0, "block_lens": []}

    weights, ifis = [], []
    for s, e in blocks:
        u_b = u_full[s:e]
        v_b = v_full[s:e]
        if len(u_b) < (max(abs(L) for L in lags) + 5):
            continue
        rl = lag_corr_simple(u_b, v_b, lags)
        ifi_b = ifi_from_lags(rl, lags)
        if not np.isnan(ifi_b):
            weights.append(e - s)
            ifis.append(ifi_b)
    if not ifis:
        return float("nan"), {"n_blocks": len(blocks), "block_lens": [e - s for s, e in blocks]}
    w = np.asarray(weights, dtype=float)
    ifi_trial = float(np.sum(w * np.asarray(ifis)) / np.sum(w))
    return ifi_trial, {"n_blocks": len(blocks), "block_lens": [e - s for s, e in blocks]}


def rz_align_pre(
    X_trial: np.ndarray,
    Y_trial: np.ndarray,
    pos_trial: np.ndarray,
    is_valid_trial: np.ndarray,
    rz_entry_cm: float,
    n_rz_bins: int,
):
    """Take the n_rz_bins immediately *before* the first bin where pos >= rz_entry_cm.

    Returns (X_seg, Y_seg) of shape (n_rz_bins, k) each, or (None, None) if the
    segment is shorter than n_rz_bins or contains any invalid bins.
    """
    above = np.where(pos_trial >= rz_entry_cm)[0]
    if len(above) == 0:
        return None, None
    i_rz = int(above[0])
    s = i_rz - n_rz_bins
    if s < 0:
        return None, None
    seg = slice(s, i_rz)
    if not np.all(is_valid_trial[seg]):
        return None, None
    return X_trial[seg, :], Y_trial[seg, :]


# -----------------------------------------------------------------------------
# Synthetic generators
# -----------------------------------------------------------------------------

def make_lagged_ar1(n_bins, k1, k2, lag, coupling=0.7, ar=0.6, rng=None):
    """Generate two regions where Y(t) is driven by X(t - lag) plus AR(1) noise.

    lag > 0: u leads v (X drives Y in the future) -> per-trial IFI should be < 0.
    lag < 0: v leads u                              -> IFI > 0.
    """
    rng = rng or np.random.default_rng()
    nx = rng.standard_normal((n_bins, k1))
    X = np.zeros_like(nx)
    X[0] = nx[0]
    for t in range(1, n_bins):
        X[t] = ar * X[t - 1] + nx[t]
    ny = rng.standard_normal((n_bins, k2))
    Y = np.zeros_like(ny)
    # Coupling matrix from X -> Y
    W = rng.standard_normal((k1, k2)) * (coupling / math.sqrt(k1))
    for t in range(n_bins):
        src_t = t - lag
        if 0 <= src_t < n_bins:
            Y[t] = X[src_t] @ W + 0.5 * ny[t]
        else:
            Y[t] = 0.5 * ny[t]
        if t > 0:
            Y[t] += ar * Y[t - 1] - ar * (Y[t - 1])  # keep AR mild on Y
            Y[t] = 0.5 * Y[t] + 0.5 * (ar * Y[t - 1])
    return X, Y


def make_uncoupled(n_bins, k1, k2, ar=0.6, rng=None):
    rng = rng or np.random.default_rng()
    def ar1(k):
        n = rng.standard_normal((n_bins, k))
        x = np.zeros_like(n)
        x[0] = n[0]
        for t in range(1, n_bins):
            x[t] = ar * x[t - 1] + n[t]
        return x
    return ar1(k1), ar1(k2)


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

def assert_true(cond, msg):
    if not cond:
        print(f"FAIL: {msg}")
        sys.exit(1)
    print(f"  ok: {msg}")


def test_lagged_coupling():
    print("\n[T1] Lagged AR(1) coupling, lag=+1 (u leads v)")
    rng = np.random.default_rng(1)
    k1, k2, n_bins = 3, 3, 200
    n_trials = 12
    ifis = []
    rs = []
    for _ in range(n_trials):
        X, Y = make_lagged_ar1(n_bins, k1, k2, lag=1, coupling=0.9, rng=rng)
        out = per_trial_cca(X, Y)
        assert_true(out is not None, "per_trial_cca returns a fit")
        r, A, B = out
        rs.append(r)
        is_valid = np.ones(n_bins, dtype=bool)
        ifi, _ = per_trial_ifi_weighted(X, Y, is_valid, A, B,
                                         lags=range(-3, 4), min_block_bins=12)
        ifis.append(ifi)
    mean_ifi = float(np.nanmean(ifis))
    mean_r = float(np.nanmean(rs))
    print(f"    mean per-trial CC = {mean_r:.3f}")
    print(f"    mean per-trial IFI = {mean_ifi:.3f}  (expect < 0)")
    assert_true(mean_r > 0.5, "per-trial CC > 0.5 with strong coupling")
    assert_true(mean_ifi < -0.05, "per-trial IFI < -0.05 (u leads v -> IFI negative)")


def test_uncoupled_null():
    print("\n[T2] Uncoupled regions, real CC vs shuffle CC")
    rng = np.random.default_rng(2)
    k1, k2, n_bins = 3, 3, 200
    n_trials = 20
    rs_real, rs_shuff = [], []
    for _ in range(n_trials):
        X, Y = make_uncoupled(n_bins, k1, k2, rng=rng)
        out = per_trial_cca(X, Y)
        if out is None: continue
        r, _, _ = out
        rs_real.append(r)
        # shuffle: permute Y rows
        perm = rng.permutation(n_bins)
        out_s = per_trial_cca(X, Y[perm])
        if out_s is None: continue
        rs_shuff.append(out_s[0])
    mean_r = np.mean(rs_real); mean_s = np.mean(rs_shuff)
    print(f"    real CC = {mean_r:.3f}, shuffle CC = {mean_s:.3f}")
    # Both will be non-zero due to small-sample bias of canoncorr;
    # with 200 samples and k=3+3 the bias should be similar and small.
    assert_true(abs(mean_r - mean_s) < 0.10,
                "real CC ~ shuffle CC for uncoupled regions (within 0.10)")


def test_block_weighted_ifi():
    print("\n[T3] Block weighted IFI: short noisy block + long signal block")
    rng = np.random.default_rng(3)
    k1, k2 = 3, 3
    # Long block: lag=+1, strong coupling -> IFI < 0
    X_long, Y_long = make_lagged_ar1(180, k1, k2, lag=1, coupling=0.9, rng=rng)
    # Short noisy block: uncoupled
    X_short, Y_short = make_uncoupled(15, k1, k2, rng=rng)
    # Stitch with a 10-bin invalid gap in between
    gap = 10
    X = np.vstack([X_short, np.zeros((gap, k1)), X_long])
    Y = np.vstack([Y_short, np.zeros((gap, k2)), Y_long])
    is_valid = np.concatenate([np.ones(15, bool), np.zeros(gap, bool), np.ones(180, bool)])
    # Fit on all valid bins (concat; per-trial CCA uses only valid ones)
    X_valid = X[is_valid]; Y_valid = Y[is_valid]
    out = per_trial_cca(X_valid, Y_valid)
    assert_true(out is not None, "per_trial_cca handles concatenated valid bins")
    r, A, B = out
    ifi, info = per_trial_ifi_weighted(X, Y, is_valid, A, B,
                                        lags=range(-3, 4), min_block_bins=12)
    print(f"    blocks: {info['block_lens']}, weighted IFI = {ifi:.3f}")
    assert_true(info["n_blocks"] == 2, "two contiguous blocks identified")
    # Short block (15 bins) above min_block (12), but its IFI may dominate noise.
    # With weights [15, 180], the 180-bin block dominates -> IFI < 0.
    assert_true(ifi < 0, "weighted IFI dominated by long lagged block")

    # Sanity: drop short block below min_block -> only long block contributes.
    is_valid2 = is_valid.copy(); is_valid2[:15] = False
    ifi2, info2 = per_trial_ifi_weighted(X, Y, is_valid2, A, B,
                                          lags=range(-3, 4), min_block_bins=12)
    print(f"    after dropping short: blocks={info2['block_lens']}, IFI={ifi2:.3f}")
    assert_true(info2["n_blocks"] == 1 and info2["block_lens"][0] == 180,
                "min_block_bins filter drops short block")


def test_rz_alignment():
    print("\n[T4] RZ-pre alignment: coupling burst only in last 250 ms before RZ")
    print("     (compare real CC to shuffle CC at each aligned bin to remove"
          " small-sample bias)")
    rng = np.random.default_rng(4)
    k1, k2 = 3, 3
    n_bins = 400
    rz_entry_cm = 200.0
    pos = np.linspace(0, 250, n_bins)
    i_rz = int(np.where(pos >= rz_entry_cm)[0][0])
    print(f"    rz_entry bin = {i_rz}")

    n_rz_bins = 20
    n_trials = 60                     # more trials so n > k1+k2 with margin
    X_aligned, Y_aligned, pre_ifis = [], [], []

    for itr in range(n_trials):
        X, Y = make_uncoupled(n_bins, k1, k2, rng=rng)
        burst = slice(i_rz - 10, i_rz)
        Xb = X[burst]
        W = rng.standard_normal((k1, k2)) * 0.95 / math.sqrt(k1)
        # u leads v by 1 bin within the burst
        for j, t in enumerate(range(burst.start + 1, burst.stop)):
            Y[t] = Xb[j] @ W + 0.3 * rng.standard_normal(k2)

        out = per_trial_cca(X, Y)
        if out is None: continue
        r, A, B = out
        is_valid = np.ones(n_bins, bool)
        Xs, Ys = rz_align_pre(X, Y, pos, is_valid, rz_entry_cm, n_rz_bins)
        if Xs is None: continue
        X_aligned.append(Xs); Y_aligned.append(Ys)

        Xc = Xs - Xs.mean(0, keepdims=True)
        Yc = Ys - Ys.mean(0, keepdims=True)
        u = (Xc @ A).ravel(); v = (Yc @ B).ravel()
        rl = lag_corr_simple(u, v, range(-3, 4))
        pre_ifi = ifi_from_lags(rl, range(-3, 4))
        if not np.isnan(pre_ifi):
            pre_ifis.append(pre_ifi)

    X_aligned = np.stack(X_aligned, axis=0)
    Y_aligned = np.stack(Y_aligned, axis=0)
    n_valid = X_aligned.shape[0]
    print(f"    valid trials = {n_valid}")
    assert_true(n_valid >= 40, "most trials yield a pre-RZ segment")

    # Per-aligned-bin CC across trials (real and shuffle)
    n_shuf = 20
    cc_real = np.full(n_rz_bins, np.nan)
    cc_shuf = np.full(n_rz_bins, np.nan)
    for b in range(n_rz_bins):
        Xb = X_aligned[:, b, :]
        Yb = Y_aligned[:, b, :]
        out = per_trial_cca(Xb, Yb)
        if out is not None: cc_real[b] = out[0]
        rs = []
        for _ in range(n_shuf):
            perm = rng.permutation(n_valid)
            o = per_trial_cca(Xb, Yb[perm])
            if o is not None: rs.append(o[0])
        if rs: cc_shuf[b] = np.mean(rs)
    cc_excess = cc_real - cc_shuf
    print(f"    real - shuffle CC, last 5 bins: {cc_excess[-5:]}")
    print(f"    real - shuffle CC, first 5 bins: {cc_excess[:5]}")
    assert_true(np.nanmean(cc_excess[-5:]) > np.nanmean(cc_excess[:5]) + 0.05,
                "real-minus-shuffle CC is higher at the end of pre-RZ window")
    mean_pre_ifi = float(np.nanmean(pre_ifis))
    print(f"    mean pre-RZ IFI = {mean_pre_ifi:.3f} (expect < 0; lag=+1)")
    assert_true(mean_pre_ifi < 0.0, "pre-RZ per-trial IFI is negative on average")


if __name__ == "__main__":
    test_lagged_coupling()
    test_uncoupled_null()
    test_block_weighted_ifi()
    test_rz_alignment()
    print("\nAll prototype tests passed.")
