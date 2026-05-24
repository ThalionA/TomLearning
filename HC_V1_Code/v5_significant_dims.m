function sig_mask = v5_significant_dims(cc_observed, cc_null, alpha)
%V5_SIGNIFICANT_DIMS  Per-dimension significance vs a shuffle null.
%
%   SIG_MASK = V5_SIGNIFICANT_DIMS(CC_OBSERVED, CC_NULL) returns a
%   logical mask (r x 1) marking each canonical dimension as
%   significant if the observed correlation exceeds the
%   (1 - alpha) quantile of its per-dimension null distribution.
%
%   SIG_MASK = V5_SIGNIFICANT_DIMS(CC_OBSERVED, CC_NULL, ALPHA)
%   overrides the default alpha = 0.05 (i.e. 95th percentile threshold).
%
%   This follows the Gonzalez et al. paper's procedure for counting
%   significant subspace dimensions: shuffle the trial axis of one
%   population, refit, record the full canonical-correlation vector
%   per shuffle, and threshold each dimension at the 95th percentile
%   of its null. The v5 shuffle helper (v5_shuffle_null) already
%   returns the full (n_shuffles x r) cc_cv matrix, so this is a
%   thin wrapper that turns it into a per-dim significance mask.
%
%   Inputs
%       cc_observed  (r x 1) double. Observed canonical correlations
%                    for one fit (typically cc_cv from v5_cca_fit_cv).
%       cc_null      (n_shuffles x r) double. Per-shuffle, per-dim
%                    canonical correlations from v5_shuffle_null.cc_cv.
%       alpha        scalar in (0, 1). Default 0.05.
%
%   Output
%       sig_mask     (r x 1) logical. True where cc_observed(i) exceeds
%                    the (1 - alpha) quantile of cc_null(:, i).
%
%   Notes
%       - Strict inequality (`>`), matching the paper.
%       - NaN entries in cc_null are excluded from the quantile.
%       - If cc_null has fewer columns than cc_observed (e.g. caller
%         passed an outdated null), an error is raised.

    if nargin < 3 || isempty(alpha), alpha = 0.05; end
    if ~isscalar(alpha) || alpha <= 0 || alpha >= 1
        error('v5_significant_dims:badAlpha', ...
            'alpha must be a scalar in (0, 1); got %g', alpha);
    end
    if ~isvector(cc_observed)
        error('v5_significant_dims:badShape', 'cc_observed must be a vector');
    end
    if ~ismatrix(cc_null)
        error('v5_significant_dims:badShape', 'cc_null must be 2D');
    end
    cc_observed = double(cc_observed(:));
    r = numel(cc_observed);
    if size(cc_null, 2) ~= r
        error('v5_significant_dims:badShape', ...
            'cc_null has %d columns; cc_observed has %d entries', ...
            size(cc_null, 2), r);
    end

    sig_mask = false(r, 1);
    for i = 1:r
        col = cc_null(:, i);
        col = col(isfinite(col));
        if isempty(col)
            % No null samples for this dim — cannot decide. Mark as
            % not significant to be conservative.
            continue;
        end
        thresh = quantile(col, 1 - alpha);
        sig_mask(i) = cc_observed(i) > thresh;
    end
end
