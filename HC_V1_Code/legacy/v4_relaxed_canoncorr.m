function [r, ok, info] = v4_relaxed_canoncorr(X, Y, min_extra)
%V4_RELAXED_CANONCORR  canoncorr with a configurable sample-margin guard.
%   [r, ok, info] = v4_relaxed_canoncorr(X, Y, min_extra)
%
%   Same rank check as v4_per_trial_cca, but the sample-size guard is
%   n > k1 + k2 + min_extra (default 4). Used by the RZ-per-bin path,
%   where n = trials_in_epoch is small and we accept some upward bias in
%   absolute CC, relying on real-minus-shuffle as the interpretable
%   quantity.
%
%   Inputs:
%     X         : [n x k1]
%     Y         : [n x k2]
%     min_extra : (default 4) minimum n - (k1 + k2)
%
%   Returns:
%     r    : first canonical correlation (NaN if skipped).
%     ok   : true if canoncorr ran successfully.
%     info : struct with fields .reason ('ok' / 'insufficient_samples' /
%            'rank_deficient' / 'canoncorr_error'), .n_samples, .k1, .k2.

    if nargin < 3 || isempty(min_extra), min_extra = 4; end

    info = struct('reason', '', ...
                  'n_samples', size(X, 1), ...
                  'k1', size(X, 2), ...
                  'k2', size(Y, 2));
    n  = info.n_samples;
    k1 = info.k1;
    k2 = info.k2;
    r  = NaN; ok = false;

    if n <= (k1 + k2 + min_extra)
        info.reason = 'insufficient_samples'; return;
    end
    Xc = X - repmat(mean(X, 1), n, 1);
    Yc = Y - repmat(mean(Y, 1), n, 1);
    if rank(Xc) < k1 || rank(Yc) < k2
        info.reason = 'rank_deficient'; return;
    end
    try
        [~, ~, rall] = canoncorr(X, Y);
    catch
        info.reason = 'canoncorr_error'; return;
    end
    r = rall(1); ok = true; info.reason = 'ok';
end
