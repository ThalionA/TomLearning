function [r, A, B, info] = v4_per_trial_cca(X_valid, Y_valid)
%V4_PER_TRIAL_CCA  Single canoncorr fit on the valid bins of one trial.
%   [r, A, B, info] = v4_per_trial_cca(X_valid, Y_valid)
%
%   X_valid : [n_valid x k1] PCA-reduced timecourse for region 1.
%   Y_valid : [n_valid x k2] PCA-reduced timecourse for region 2.
%
%   Returns:
%     r    : first canonical correlation (scalar). NaN if the fit was skipped.
%     A,B  : canonical loadings (k1 x 1 and k2 x 1). NaN-filled if skipped.
%     info : struct with fields:
%              .ok        - logical, true if canoncorr ran successfully
%              .reason    - char, reason for skip (or 'ok')
%              .n_samples - number of samples used
%              .k1, .k2   - dimensionality of inputs
%
%   Skip rules (match the spirit of v3's per-bin checks but applied at the
%   trial level):
%     - n_samples <= k1 + k2 + 4   ('insufficient_samples')
%     - rank(Xc) < k1 or rank(Yc) < k2 ('rank_deficient')
%     - canoncorr throws ('canoncorr_error')

    info = struct('ok', false, 'reason', '', ...
                  'n_samples', size(X_valid, 1), ...
                  'k1', size(X_valid, 2), ...
                  'k2', size(Y_valid, 2));
    k1 = info.k1; k2 = info.k2; n = info.n_samples;
    A = nan(k1, 1); B = nan(k2, 1); r = NaN;

    if n <= (k1 + k2 + 4)
        info.reason = 'insufficient_samples'; return;
    end
    Xc = X_valid - repmat(mean(X_valid, 1), n, 1);
    Yc = Y_valid - repmat(mean(Y_valid, 1), n, 1);
    if rank(Xc) < k1 || rank(Yc) < k2
        info.reason = 'rank_deficient'; return;
    end
    try
        [Aall, Ball, rall] = canoncorr(X_valid, Y_valid);
    catch
        info.reason = 'canoncorr_error'; return;
    end
    A = Aall(:, 1); B = Ball(:, 1); r = rall(1);
    info.ok = true; info.reason = 'ok';
end
