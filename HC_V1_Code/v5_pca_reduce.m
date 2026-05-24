function [P, state] = v5_pca_reduce(S, k)
%V5_PCA_REDUCE  SVD-based PCA on flattened (trials × bins, units) matrix.
%
%   [P, STATE] = V5_PCA_REDUCE(S, K) fits PCA on the (n_trials × n_bins,
%   n_units) flattened view of S and returns the top-K PC scores reshaped
%   back to (n_trials × n_bins × K), plus a STATE struct sufficient to
%   apply the same projection to held-out data via v5_pca_apply.
%
%   Stage 3 of the spec (ResearchVault/Methods/CCA_HH_Adapted.md §2).
%
%   Inputs
%       S   (n_trials × n_bins × n_units) double. Residualised data
%           from Stage 2, or raw if the caller chooses to skip Stage 2.
%       K   positive integer ≤ min(n_trials*n_bins, n_units).
%
%   Outputs
%       P       (n_trials × n_bins × K) double. PC scores.
%       state   struct with fields:
%           mean                      (1 × n_units)
%           components                (K × n_units)   rows are PCs
%           explained_variance_ratio  (K × 1)
%
%   Sign convention: sklearn's `svd_flip` on V — for each PC vector,
%   find the index of the largest absolute element and force its sign
%   positive. This is what makes Python and MATLAB agree element-for-
%   element on the cross-language test.
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_pca_reduce.m
%   The Python reference is at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::pca_reduce

    sz = size(S);
    if numel(sz) ~= 3
        error('v5_pca_reduce:badShape', ...
            'S must be 3D (n_trials × n_bins × n_units); got size [%s]', ...
            num2str(sz));
    end
    n_trials = sz(1);
    n_bins   = sz(2);
    n_units  = sz(3);
    n_samples = n_trials * n_bins;

    if ~isscalar(k) || k ~= round(k) || k < 1
        error('v5_pca_reduce:badK', 'k must be a positive integer; got %s', ...
            mat2str(k));
    end
    k_max = min(n_samples, n_units);
    if k > k_max
        error('v5_pca_reduce:badK', ...
            'k=%d exceeds min(n_samples=%d, n_units=%d) = %d', ...
            k, n_samples, n_units, k_max);
    end

    X = reshape(double(S), n_samples, n_units);
    mu = mean(X, 1);                            % (1 × n_units)
    Xc = X - mu;

    [U, Sigma, V] = svd(Xc, 'econ');
    sigma = diag(Sigma);                        % column vector of singular values

    % svd_flip: make max-abs entry of each column of V positive. V has
    % shape (n_units × r); each column is one PC vector. Flip both V's
    % column and U's column to keep U*Sigma*V' unchanged.
    r = size(V, 2);
    signs = ones(1, r);
    for i = 1:r
        [~, j] = max(abs(V(:, i)));
        s = sign(V(j, i));
        if s == 0, s = 1; end
        signs(i) = s;
    end
    V = V .* signs;
    U = U .* signs;

    % Variance explained ratio over the full singular spectrum.
    total_var = sum(sigma .^ 2);
    if total_var > 0
        evr_all = (sigma .^ 2) / total_var;
    else
        evr_all = zeros(size(sigma));
    end
    evr_k = evr_all(1:k);

    components = V(:, 1:k).';                   % (k × n_units), rows are PCs
    P_flat = U(:, 1:k) .* sigma(1:k).';         % (n_samples × k)
    P = reshape(P_flat, n_trials, n_bins, k);

    state = struct( ...
        'mean', mu, ...
        'components', components, ...
        'explained_variance_ratio', evr_k);
end
