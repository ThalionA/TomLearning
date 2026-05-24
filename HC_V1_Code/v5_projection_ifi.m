function [ifi, lag_corr] = v5_projection_ifi(u, v, max_lag, min_paired_samples)
%V5_PROJECTION_IFI  Per-trial projection IFI from lag-correlation asymmetry.
%
%   [IFI, LAG_CORR] = V5_PROJECTION_IFI(U, V) defaults: max_lag=8,
%   min_paired_samples=5.
%
%   For each trial: compute lag correlations r(tau) =
%   corr(u_trial(t), v_trial(t+tau)) for tau in -max_lag..+max_lag,
%   then form the signed asymmetry index:
%
%       IFI = (sum_{tau>0} r(tau) - sum_{tau<0} r(tau)) /
%             (sum_{tau>0} |r(tau)| + sum_{tau<0} |r(tau)|)
%
%   Bounded in [-1, 1]. Positive => u leads v (consistent with
%   lagged_refit_ifi sign convention). Lags with fewer than
%   min_paired_samples after trim are dropped (NaN) rather than zero.
%
%   Stage 8 of the H&H-adapted CCA pipeline.
%
%   Inputs
%       u, v                (n_trials x n_bins) projected canonical variates
%       max_lag             max |tau| in bins. Default 8.
%       min_paired_samples  drop a lag if fewer paired samples remain
%                           after trim. Default 5.
%
%   Outputs
%       ifi       (n_trials x 1). Per-trial IFI; NaN where insufficient.
%       lag_corr  (n_trials x (2*max_lag + 1)). Per-trial lag-corr curve
%                 with index (tau + max_lag + 1) for tau in -max_lag..max_lag.
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_projection_ifi.m
%   Python reference at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::projection_ifi

    if nargin < 3 || isempty(max_lag),            max_lag = 8; end
    if nargin < 4 || isempty(min_paired_samples), min_paired_samples = 5; end

    if ~ismatrix(u) || ~ismatrix(v)
        error('v5_projection_ifi:badShape', 'u, v must be 2D');
    end
    if ~isequal(size(u), size(v))
        error('v5_projection_ifi:badShape', ...
            'u, v shape mismatch: %s vs %s', mat2str(size(u)), mat2str(size(v)));
    end
    [n_trials, n_bins] = size(u);
    if max_lag < 1
        error('v5_projection_ifi:badShape', 'max_lag must be >= 1');
    end
    if max_lag >= n_bins
        error('v5_projection_ifi:badShape', ...
            'max_lag=%d must be < n_bins=%d', max_lag, n_bins);
    end

    n_lag_axis = 2 * max_lag + 1;
    lag_corr = NaN(n_trials, n_lag_axis);
    ifi = NaN(n_trials, 1);

    for t = 1:n_trials
        ut = u(t, :);
        vt = v(t, :);
        for tau = -max_lag:max_lag
            if tau >= 0
                a = ut(1:n_bins - tau);
                b = vt(tau + 1:end);
            else
                a = ut(-tau + 1:end);
                b = vt(1:n_bins + tau);
            end
            if numel(a) < min_paired_samples, continue; end
            std_a = std(a);
            std_b = std(b);
            if std_a < 1e-15 || std_b < 1e-15, continue; end
            c = corrcoef(a, b);
            lag_corr(t, tau + max_lag + 1) = c(1, 2);
        end
        r_neg = lag_corr(t, 1:max_lag);
        r_pos = lag_corr(t, max_lag + 2:end);
        if any(isnan(r_neg)) || any(isnan(r_pos)), continue; end
        sum_pos = sum(r_pos);
        sum_neg = sum(r_neg);
        denom = sum(abs(r_pos)) + sum(abs(r_neg));
        if denom > 1e-15
            ifi(t) = (sum_pos - sum_neg) / denom;
        else
            ifi(t) = 0;
        end
    end
end
