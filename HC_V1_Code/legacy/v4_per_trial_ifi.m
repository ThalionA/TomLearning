function [ifi_trial, info] = v4_per_trial_ifi(X_trial, Y_trial, is_valid, A, B, lags, min_block_bins)
%V4_PER_TRIAL_IFI  Per-trial IFI as length-weighted mean of per-block IFIs.
%   [ifi_trial, info] = v4_per_trial_ifi(X_trial, Y_trial, is_valid, A, B, lags, min_block_bins)
%
%   Computes IFI on the canonical projections u = Xc * A, v = Yc * B over
%   each contiguous run of valid bins (length >= min_block_bins) within the
%   trial, then averages the per-block IFIs weighted by block length.
%
%   X_trial, Y_trial : [n_tr_bins x k1] / [n_tr_bins x k2]
%   is_valid         : [n_tr_bins x 1] logical, true where bin is in
%                      base_valid (cued & vel >= min_speed_cms).
%   A, B             : canonical loadings from v4_per_trial_cca for this
%                      trial (k1 x 1, k2 x 1).
%   lags             : vector of lag offsets (e.g. -3:3).
%   min_block_bins   : scalar, minimum block length to contribute.
%
%   Returns:
%     ifi_trial : weighted-mean IFI across eligible blocks. NaN if no
%                 eligible block produced a finite IFI.
%     info      : struct with fields:
%                   .blocks      - Nx2 [start stop] of all blocks meeting
%                                  the length filter
%                   .block_lens  - 1xN block lengths
%                   .block_ifis  - 1xN per-block IFIs (NaN if a block
%                                  could not yield a finite IFI)
%                   .block_used  - 1xN logical, true if the block's IFI
%                                  was included in the weighted mean
%                   .total_used_bins - sum of lengths over used blocks
%                   .reason      - 'ok' / 'no_eligible_block' /
%                                  'no_finite_ifi'
%
%   The trial-level (A, B) basis is intentional: we want IFI to measure
%   lag asymmetry along the direction of maximum trial-level coupling,
%   not a per-block redirection that would inflate variance.

    info = struct('blocks', zeros(0, 2), 'block_lens', [], ...
                  'block_ifis', [], 'block_used', false(1,0), ...
                  'total_used_bins', 0, 'reason', '');

    n = size(X_trial, 1);
    if n == 0 || any(isnan(A)) || any(isnan(B))
        info.reason = 'no_eligible_block'; ifi_trial = NaN; return;
    end

    % Canonical projections over the full trial (zero-meaned over the
    % trial). Per-block centring would bias the lag-corr; trial-mean
    % centring is the consistent choice with v3.
    Xc = X_trial - repmat(mean(X_trial, 1), n, 1);
    Yc = Y_trial - repmat(mean(Y_trial, 1), n, 1);
    u_full = Xc * A;
    v_full = Yc * B;

    blocks = v4_contiguous_blocks(is_valid, min_block_bins);
    if isempty(blocks)
        info.reason = 'no_eligible_block'; ifi_trial = NaN; return;
    end
    info.blocks     = blocks;
    n_blk = size(blocks, 1);
    info.block_lens = blocks(:, 2)' - blocks(:, 1)' + 1;
    info.block_ifis = nan(1, n_blk);
    info.block_used = false(1, n_blk);

    % A block needs at least max(|lag|) + 6 bins for the most extreme lag
    % overlap to have > 5 samples (matches v4_lag_corr's > 5 check).
    min_for_lag = max(abs(lags)) + 6;

    for ib = 1:n_blk
        s = blocks(ib, 1); e = blocks(ib, 2);
        if (e - s + 1) < min_for_lag, continue; end
        u_b = u_full(s:e);
        v_b = v_full(s:e);
        rl  = v4_lag_corr(u_b, v_b, lags);
        ifi_b = v4_ifi_from_lags(rl, lags);
        info.block_ifis(ib) = ifi_b;
        if isfinite(ifi_b)
            info.block_used(ib) = true;
        end
    end

    used = info.block_used;
    if ~any(used)
        info.reason = 'no_finite_ifi'; ifi_trial = NaN; return;
    end
    w = info.block_lens(used);
    info.total_used_bins = sum(w);
    ifi_trial = sum(w .* info.block_ifis(used)) / sum(w);
    info.reason = 'ok';
end
