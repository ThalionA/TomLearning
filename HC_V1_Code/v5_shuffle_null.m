function null = v5_shuffle_null(X, Y, n_shuffles, samples_per_trial, n_splits, seed, shuffle_folds, fold_seed)
%V5_SHUFFLE_NULL  H&H trial-shuffle null distribution for cc_cv.
%
%   NULL = V5_SHUFFLE_NULL(X, Y) builds the null with defaults:
%       n_shuffles=50, samples_per_trial=1, n_splits=10,
%       seed=0, shuffle_folds=true, fold_seed=0.
%
%   For each of N_SHUFFLES permutations, permute the trial axis of Y
%   (preserving within-trial sample order), refit CCA on (X, Y_shuffled)
%   with the same CV protocol as the real fit, and record the resulting
%   cc_cv vector. Threshold = mean(cc1_null) + 3 * sd(cc1_null) per H&H.
%
%   Stage 5 of the H&H-adapted CCA pipeline.
%
%   Inputs
%       X                   (n_samples × k_x)
%       Y                   (n_samples × k_y). Same n_samples, same trial layout.
%       n_shuffles          number of trial-axis permutations. Default 50.
%       samples_per_trial   bins per trial. Must divide n_samples. Default 1.
%       n_splits            forwarded to v5_cca_fit_cv. Default 10.
%       seed                master seed for the shuffle permutations.
%                           NOTE: MATLAB's Mersenne Twister produces
%                           DIFFERENT permutation sequences than
%                           numpy.default_rng for the same seed. The
%                           cross-language test uses a degenerate input
%                           where this doesn't matter.
%       shuffle_folds       forwarded to v5_cca_fit_cv. Default true.
%       fold_seed           forwarded to v5_cca_fit_cv. Default 0.
%
%   Output
%       null struct with fields:
%           cc_cv      (n_shuffles × r) per-shuffle out-of-fold CCs
%           cc1_null   (n_shuffles × 1) top-canonical CC per shuffle
%           threshold  scalar — mean(cc1_null) + 3*sd(cc1_null)
%           mean       scalar — mean(cc1_null), for shuffle-correction
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_shuffle_null.m
%   Python reference at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::shuffle_null

    if nargin < 3 || isempty(n_shuffles),         n_shuffles = 50; end
    if nargin < 4 || isempty(samples_per_trial),  samples_per_trial = 1; end
    if nargin < 5 || isempty(n_splits),           n_splits = 10; end
    if nargin < 6 || isempty(seed),               seed = 0; end
    if nargin < 7 || isempty(shuffle_folds),      shuffle_folds = true; end
    if nargin < 8 || isempty(fold_seed),          fold_seed = 0; end

    if ~ismatrix(X) || ~ismatrix(Y)
        error('v5_shuffle_null:badShape', 'X, Y must be 2D');
    end
    n = size(X, 1);
    if size(Y, 1) ~= n
        error('v5_shuffle_null:badShape', ...
            'X, Y row count mismatch: %d vs %d', n, size(Y, 1));
    end
    if samples_per_trial < 1
        error('v5_shuffle_null:badShape', ...
            'samples_per_trial must be >= 1; got %d', samples_per_trial);
    end
    if mod(n, samples_per_trial) ~= 0
        error('v5_shuffle_null:badShape', ...
            'n_samples=%d not divisible by samples_per_trial=%d', ...
            n, samples_per_trial);
    end
    n_trials = n / samples_per_trial;
    if n_trials < 2
        error('v5_shuffle_null:badShape', ...
            'need >= 2 trials to shuffle; got %d', n_trials);
    end
    if n_shuffles < 1
        error('v5_shuffle_null:badShape', ...
            'n_shuffles must be >= 1; got %d', n_shuffles);
    end

    k_x = size(X, 2);
    k_y = size(Y, 2);
    r = min(k_x, k_y);

    % Seed the master RNG for permutation generation.
    rng_state = rng(seed, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>

    cc_cv_all = NaN(n_shuffles, r);
    for s = 1:n_shuffles
        perm = randperm(n_trials);
        Y_shuf = trial_permute(Y, perm, samples_per_trial);
        res = v5_cca_fit_cv(X, Y_shuf, n_splits, shuffle_folds, fold_seed);
        cc_cv_all(s, :) = res.cc_cv.';
    end

    cc1_null = cc_cv_all(:, 1);
    mean_null = mean(cc1_null);
    if n_shuffles > 1
        sd_null = std(cc1_null);   % ddof=1 by default in MATLAB std
    else
        sd_null = 0;
    end
    threshold = mean_null + 3 * sd_null;

    null = struct( ...
        'cc_cv', cc_cv_all, ...
        'cc1_null', cc1_null, ...
        'threshold', threshold, ...
        'mean', mean_null);
end

function Y_shuf = trial_permute(Y, perm, samples_per_trial)
% TRIAL_PERMUTE  Reorder Y's rows so trial-block i lands at position perm(i).
%   Within-trial sample order is preserved.
    spt = samples_per_trial;
    block_idx = bsxfun(@plus, (perm(:) - 1) * spt, 1:spt);
    block_idx = reshape(block_idx.', [], 1);
    Y_shuf = Y(block_idx, :);
end
