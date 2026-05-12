function S_res = v5_residualise(S, condition_labels, period_mask)
%V5_RESIDUALISE  Per-(unit, bin) conditional-PSTH subtraction within a period mask.
%
%   S_RES = V5_RESIDUALISE(S, CONDITION_LABELS) subtracts, for each unique
%   value in CONDITION_LABELS, the across-trial mean of S for that group
%   from every trial in the group, leaving S_RES whose per-group, per-
%   (bin, unit) mean is exactly zero. This is Stage 2 of the H&H-adapted
%   CCA pipeline (see ResearchVault/Methods/CCA_HH_Adapted.md §2).
%
%   S_RES = V5_RESIDUALISE(S, CONDITION_LABELS, PERIOD_MASK) restricts the
%   subtraction to bins where PERIOD_MASK is true. Bins where the mask is
%   false are returned unchanged.
%
%   Inputs
%       S                 (n_trials × n_bins × n_units) double. Trial-
%                         aligned neural data (already binned).
%       condition_labels  (n_trials × 1) numeric. Grouping variable
%                         defining the residualisation target. May
%                         contain NaN; NaNs are bucketed as their own
%                         group.
%       period_mask       (n_bins × 1) logical, optional. True for bins
%                         to residualise; false leaves bins untouched.
%                         Default: residualise every bin.
%
%   Output
%       S_res             (n_trials × n_bins × n_units) double.
%
%   Notes
%       - The input is not mutated.
%       - Singleton groups (one trial) collapse that trial to zeros
%         (its mean = itself).
%       - Per-group bin/unit residual mean should be < 1e-12 for double-
%         precision inputs.
%
%   Cross-project mirror: this file must stay byte-identical (modulo
%   path-only differences) with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_residualise.m
%   Update both in lockstep. The Python reference lives in
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::residualise

    if nargin < 3 || isempty(period_mask)
        period_mask = [];
    else
        period_mask = logical(period_mask(:));
    end

    sz = size(S);
    if numel(sz) ~= 3
        error('v5_residualise:badShape', ...
            'S must be 3D (n_trials × n_bins × n_units); got size [%s]', ...
            num2str(sz));
    end
    n_trials = sz(1);
    n_bins   = sz(2);
    n_units  = sz(3);
    if numel(condition_labels) ~= n_trials
        error('v5_residualise:badShape', ...
            'condition_labels length %d != n_trials %d', ...
            numel(condition_labels), n_trials);
    end
    if ~isempty(period_mask) && numel(period_mask) ~= n_bins
        error('v5_residualise:badMask', ...
            'period_mask length %d != n_bins %d', ...
            numel(period_mask), n_bins);
    end

    condition_labels = condition_labels(:);  % force column
    S_res = double(S);                       % working copy, double precision

    % Determine which bin indices to operate on.
    if isempty(period_mask)
        bin_sel = 1:n_bins;
    else
        bin_sel = find(period_mask);
    end
    if isempty(bin_sel)
        return;  % nothing to do
    end

    % Identify NaN-labeled trials (own bucket).
    if isfloat(condition_labels)
        nan_mask = isnan(condition_labels);
    else
        nan_mask = false(size(condition_labels));
    end
    finite_labels = condition_labels(~nan_mask);
    if isempty(finite_labels)
        unique_groups = [];
    else
        unique_groups = unique(finite_labels);
    end

    % Subtract group mean within bin_sel for each group of finite labels.
    for g_idx = 1:numel(unique_groups)
        g = unique_groups(g_idx);
        trial_mask = condition_labels == g;
        if ~any(trial_mask), continue; end
        group_mean = mean(S(trial_mask, bin_sel, :), 1);   % 1 × |bin_sel| × n_units
        S_res(trial_mask, bin_sel, :) = ...
            S(trial_mask, bin_sel, :) - group_mean;
    end

    % NaN group (if any).
    if any(nan_mask)
        group_mean = mean(S(nan_mask, bin_sel, :), 1);
        S_res(nan_mask, bin_sel, :) = ...
            S(nan_mask, bin_sel, :) - group_mean;
    end
end
