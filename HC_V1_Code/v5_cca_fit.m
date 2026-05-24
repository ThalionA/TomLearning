function res = v5_cca_fit(X, Y)
%V5_CCA_FIT  Closed-form CCA via QR + SVD.
%
%   RES = V5_CCA_FIT(X, Y) fits CCA on (X, Y) and returns a struct
%   with fields A, B, cc, cc_cv (NaN here — see v5_cca_fit_cv for CV),
%   cc_cv_per_fold (NaN), x_mean, y_mean.
%
%   Stage 4 inner of the H&H-adapted CCA pipeline. Mean-centres both
%   arrays internally, performs QR on each, then SVDs the cross-product
%   of the orthonormal bases. The returned loadings (A, B) satisfy:
%       corr((X − x_mean) * A(:, i), (Y − y_mean) * B(:, i)) = cc(i)
%   to numerical precision, for each i = 1..r where r = min(k_x, k_y).
%
%   Sign convention: max-abs entry of each column of A is positive; B
%   carries the same sign. Matches the Python reference.
%
%   Inputs
%       X   (n_samples × k_x) double. n_samples ≥ k_x required.
%       Y   (n_samples × k_y) double. Same n_samples as X; n ≥ k_y.
%
%   Output
%       res struct with fields:
%           A               (k_x × r) canonical loadings
%           B               (k_y × r) canonical loadings
%           cc              (r × 1)  canonical correlations
%           cc_cv           (r × 1)  NaN — set by v5_cca_fit_cv
%           cc_cv_per_fold  (n_splits × r) NaN — set by v5_cca_fit_cv
%           x_mean          (1 × k_x) per-feature mean of X
%           y_mean          (1 × k_y) per-feature mean of Y
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_cca_fit.m
%   Python reference at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::cca_fit

    if ~ismatrix(X) || ~ismatrix(Y)
        error('v5_cca_fit:badShape', 'X, Y must be 2D');
    end
    [n, k_x] = size(X);
    [n_y, k_y] = size(Y);
    if n ~= n_y
        error('v5_cca_fit:badShape', ...
            'X and Y must have same n_samples; got %d vs %d', n, n_y);
    end
    if n < max(k_x, k_y)
        error('v5_cca_fit:badShape', ...
            'n_samples=%d must be >= max(k_x=%d, k_y=%d)', n, k_x, k_y);
    end

    X = double(X);
    Y = double(Y);
    x_mean = mean(X, 1);                        % (1 × k_x)
    y_mean = mean(Y, 1);                        % (1 × k_y)
    Xc = X - x_mean;
    Yc = Y - y_mean;

    % Economy QR (MATLAB returns reduced when input is tall).
    [Qx, Rx] = qr(Xc, 0);                       % Qx: n×k_x, Rx: k_x×k_x
    [Qy, Ry] = qr(Yc, 0);                       % Qy: n×k_y, Ry: k_y×k_y
    % Canonicalise QR so diag(R) >= 0. Required for cross-language
    % determinism since neither MATLAB nor numpy enforce this natively;
    % without it, sign differences in R propagate to A, B downstream.
    sx = sign(diag(Rx));  sx(sx == 0) = 1;
    sy = sign(diag(Ry));  sy(sy == 0) = 1;
    Qx = Qx .* sx.';
    Rx = Rx .* sx;
    Qy = Qy .* sy.';
    Ry = Ry .* sy;

    % Condition-number diagnostic: warn when Rx or Ry is ill-conditioned.
    % This is the root cause of non-finite canonical loadings (A, B) and
    % the downstream NaN principal-angles bug.
    cond_rx = cond(Rx);
    cond_ry = cond(Ry);
    if cond_rx > 1e10
        warning('v5_cca_fit:illConditioned', ...
            'cond(Rx) = %.2e — X-side loadings may be unreliable.', cond_rx);
    end
    if cond_ry > 1e10
        warning('v5_cca_fit:illConditioned', ...
            'cond(Ry) = %.2e — Y-side loadings may be unreliable.', cond_ry);
    end

    M = Qx.' * Qy;                              % (k_x × k_y)
    [U, S, V] = svd(M, 'econ');
    sigma = diag(S);                            % column of singular values
    r = min(k_x, k_y);
    cc = max(0, min(1, sigma(1:r)));            % numerical clamp to [0,1]

    % Canonical loadings: solve Rx * A = U(:, 1:r), Ry * B = V(:, 1:r).
    % Use pinv (Moore-Penrose pseudo-inverse) instead of mldivide so that
    % near-singular Rx / Ry give a minimum-norm solution instead of NaN/Inf
    % entries. On well-conditioned inputs pinv and `\` agree at machine
    % precision (verified by the cross-language reference test). Degenerate
    % PCs (those without enough variance in the training fold to be well-
    % defined) produce small loadings, not garbage propagating downstream.
    A = pinv(Rx) * U(:, 1:r);                   % (k_x × r)
    B = pinv(Ry) * V(:, 1:r);                   % (k_y × r)

    % Sign convention: make max-abs of each column of A positive.
    % Guard: if a column has any non-finite entry (NaN/Inf from a rank-
    % deficient upstream that pinv couldn't tame), skip the sign fix.
    % Multiplying by NaN/Inf would poison the entire column; leaving it
    % alone lets v5_principal_angles' sanitise_basis drop it downstream
    % rather than letting non-finite values propagate silently.
    for i = 1:r
        col_abs = abs(A(:, i));
        if any(~isfinite(col_abs))
            continue;
        end
        [~, j] = max(col_abs);
        s = sign(A(j, i));
        if s == 0, s = 1; end
        A(:, i) = A(:, i) * s;
        B(:, i) = B(:, i) * s;
    end

    res = struct( ...
        'A', A, ...
        'B', B, ...
        'cc', cc, ...
        'cc_cv', NaN(r, 1), ...
        'cc_cv_per_fold', [], ...
        'x_mean', x_mean, ...
        'y_mean', y_mean);
end
