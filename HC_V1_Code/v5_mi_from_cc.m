function mi = v5_mi_from_cc(cc, sig_mask)
%V5_MI_FROM_CC  Gaussian-MI approximation from canonical correlations.
%
%   MI = V5_MI_FROM_CC(CC) computes the mutual-information surrogate
%   used in the Gonzalez et al. pCCA paper (Adapted for our spec —
%   plain CCA + per-period PSTH residualisation, but the MI formula is
%   the same):
%
%       MI = -0.5 * sum_i log(1 - cc_i^2)    in nats
%
%   This is the exact MI for two jointly Gaussian populations whose
%   canonical correlations are cc_i. Used as a single-number summary
%   of shared information per fit (animal x area-pair x epoch).
%
%   MI = V5_MI_FROM_CC(CC, SIG_MASK) sums only over the dimensions
%   where SIG_MASK is true. Useful for restricting MI to significant
%   canonical pairs identified via v5_significant_dims. Default
%   SIG_MASK = true(size(CC)).
%
%   Inputs
%       cc        (r x 1) double. Canonical correlations in [0, 1].
%       sig_mask  (r x 1) logical, optional. Default all true.
%
%   Output
%       mi        scalar double, in nats. >= 0.
%
%   Notes
%       - cc values are clamped to [0, 1 - eps] to avoid log(0).
%       - Negative cc (which can arise from per-fold CV under sign
%         confusion) are treated as their absolute value, since the
%         formula uses cc^2.
%       - Empty sig_mask (all false) yields MI = 0.

    if ~isvector(cc)
        error('v5_mi_from_cc:badShape', 'cc must be a vector');
    end
    cc = double(cc(:));
    r = numel(cc);
    if nargin < 2 || isempty(sig_mask)
        sig_mask = true(r, 1);
    else
        if ~isvector(sig_mask) || numel(sig_mask) ~= r
            error('v5_mi_from_cc:badShape', ...
                'sig_mask length %d must equal cc length %d', ...
                numel(sig_mask), r);
        end
        sig_mask = logical(sig_mask(:));
    end

    if ~any(sig_mask)
        mi = 0;
        return;
    end

    % Use cc^2 so sign of cc is irrelevant. Clamp to avoid log(0) when
    % cc = 1 (perfect correlation, e.g. degenerate identical input).
    cc2 = cc(sig_mask) .^ 2;
    cc2 = min(cc2, 1 - eps);
    mi = -0.5 * sum(log(1 - cc2));
end
