function res = v5_cca_fit_cv(X, Y, n_splits, shuffle, seed)
%V5_CCA_FIT_CV  K-fold cross-validated CCA.
%
%   RES = V5_CCA_FIT_CV(X, Y) runs CCA on the full sample for the
%   canonical loadings (A, B, cc), then runs 10-fold CV to estimate
%   out-of-fold canonical correlations cc_cv. Default n_splits=10,
%   shuffle=true, seed=0.
%
%   RES = V5_CCA_FIT_CV(X, Y, N_SPLITS, SHUFFLE, SEED) overrides defaults.
%   SHUFFLE=false uses contiguous-block folds in sample order (the
%   cross-language test relies on this).
%
%   Per fold:
%       1. Fit (A_fold, B_fold) on the training rows.
%       2. Project the test rows through (A_fold, B_fold).
%       3. Compute Pearson corr of paired test variates per canonical pair.
%   Returns the mean over folds (cc_cv) plus the per-fold breakdown.
%
%   Inputs
%       X         (n_samples × k_x)
%       Y         (n_samples × k_y)
%       n_splits  positive integer ≥ 2, ≤ n_samples. Default 10.
%       shuffle   logical, default true.
%       seed      integer, default 0. Ignored if shuffle=false.
%
%   Output
%       res struct, same field set as v5_cca_fit, plus cc_cv and
%       cc_cv_per_fold filled in.
%
%   Notes
%       - Per-fold test correlations are NOT clamped to [0, 1] —
%         negative values can arise from overfit + sign confusion and
%         are reported honestly.
%       - When a fold's test set has zero variance in a canonical
%         direction (typical for very small folds), the per-fold cc is
%         set to 0 rather than NaN.
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_cca_fit_cv.m
%   Python reference at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::cca_fit_cv

    if nargin < 3 || isempty(n_splits), n_splits = 10; end
    if nargin < 4 || isempty(shuffle),  shuffle  = true; end
    if nargin < 5 || isempty(seed),     seed     = 0; end

    [n, k_x] = size(X);
    if size(Y, 1) ~= n
        error('v5_cca_fit_cv:badShape', ...
            'X, Y row count mismatch: %d vs %d', n, size(Y, 1));
    end
    k_y = size(Y, 2);
    r = min(k_x, k_y);

    if n_splits < 2
        error('v5_cca_fit_cv:badSplits', 'n_splits must be >= 2');
    end
    if n_splits > n
        error('v5_cca_fit_cv:badSplits', ...
            'n_splits=%d > n_samples=%d', n_splits, n);
    end

    % Full-sample fit.
    full = v5_cca_fit(X, Y);

    % Fold indices. Match Python's _make_fold_indices exactly when
    % shuffle=false (contiguous blocks in 1..n order). When shuffle=true
    % we just use a MATLAB-side permutation; the cross-language test
    % only runs the shuffle=false path so this is fine.
    if shuffle
        rng_state = rng(seed, 'twister');
        cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
        idx = randperm(n);
    else
        idx = 1:n;
    end
    fold_sizes = repmat(floor(n / n_splits), 1, n_splits);
    fold_sizes(1:mod(n, n_splits)) = fold_sizes(1:mod(n, n_splits)) + 1;

    cc_cv_per_fold = NaN(n_splits, r);
    start = 1;
    for f = 1:n_splits
        sz = fold_sizes(f);
        test_idx  = idx(start : start + sz - 1);
        train_idx = [idx(1:start-1), idx(start+sz:end)];
        start = start + sz;

        fit_f = v5_cca_fit(X(train_idx, :), Y(train_idx, :));
        U_te = (X(test_idx, :) - fit_f.x_mean) * fit_f.A;   % (n_test × r)
        V_te = (Y(test_idx, :) - fit_f.y_mean) * fit_f.B;
        for i = 1:r
            u = U_te(:, i);  v = V_te(:, i);
            if numel(u) < 2, continue; end
            std_u = std(u);  std_v = std(v);
            if std_u < 1e-15 || std_v < 1e-15
                cc_cv_per_fold(f, i) = 0;
            else
                c = corrcoef(u, v);
                cc_cv_per_fold(f, i) = c(1, 2);
            end
        end
    end
    cc_cv = mean(cc_cv_per_fold, 1, 'omitnan').';   % (r × 1)

    res = full;
    res.cc_cv = cc_cv;
    res.cc_cv_per_fold = cc_cv_per_fold;
end
