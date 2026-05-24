function res = v5_lagged_refit_ifi(X, Y, max_lag, central_window, smooth_size, smooth_sigma)
%V5_LAGGED_REFIT_IFI  H&H lagged-refit direction index.
%
%   RES = V5_LAGGED_REFIT_IFI(X, Y) defaults: max_lag=8, central_window=12,
%   smooth_size=3, smooth_sigma=1.
%
%   For each positive lag tau in 1..max_lag:
%       - X-leads-Y branch: pair X(t, b, :) with Y(t, b+tau, :);
%         refit CCA on the flattened (n_trials*(n_bins-tau)) matrices.
%         Record cc1.
%       - Y-leads-X branch: symmetric, pair Y(t,b,:) with X(t,b+tau,:).
%   The zero-lag fit is shared. Both cc curves are smoothed with a
%   Gaussian kernel. IFI is the (xy - yx) / (xy + yx) ratio averaged
%   over the central positive-lag window. Positive IFI => X leads Y.
%
%   Stage 7 of the H&H-adapted CCA pipeline.
%
%   Inputs
%       X, Y            (n_trials x n_bins x k) tensors. Same trial and
%                       bin axes; feature dims may differ.
%       max_lag         max lag in bins. Default 8 = 200 ms at 25 ms.
%       central_window  max tau used for IFI averaging. Default 12.
%       smooth_size     Gaussian kernel size (odd integer). Default 3.
%       smooth_sigma    Gaussian sigma. Default 1.
%
%   Output
%       res struct with fields:
%           cc_xy           (max_lag + 1) raw cc curve, X-leads-Y branch
%           cc_yx           (max_lag + 1) raw cc curve, Y-leads-X branch
%           cc_xy_smooth    smoothed cc_xy
%           cc_yx_smooth    smoothed cc_yx
%           ifi             scalar in [-1, 1]
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_lagged_refit_ifi.m
%   Python reference at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::lagged_refit_ifi

    if nargin < 3 || isempty(max_lag),        max_lag = 8;        end
    if nargin < 4 || isempty(central_window), central_window = 12; end
    if nargin < 5 || isempty(smooth_size),    smooth_size = 3;    end
    if nargin < 6 || isempty(smooth_sigma),   smooth_sigma = 1.0; end

    if ndims(X) ~= 3 || ndims(Y) ~= 3
        error('v5_lagged_refit_ifi:badShape', ...
            'X, Y must be 3D; got %s, %s', mat2str(size(X)), mat2str(size(Y)));
    end
    if ~isequal(size(X, 1), size(Y, 1)) || ~isequal(size(X, 2), size(Y, 2))
        error('v5_lagged_refit_ifi:badShape', ...
            'X, Y must share (n_trials, n_bins); got %s vs %s', ...
            mat2str(size(X)), mat2str(size(Y)));
    end
    [n_trials, n_bins, k_x] = size(X);
    k_y = size(Y, 3);
    if max_lag < 1
        error('v5_lagged_refit_ifi:badShape', 'max_lag must be >= 1');
    end
    if max_lag >= n_bins
        error('v5_lagged_refit_ifi:badShape', ...
            'max_lag=%d must be < n_bins=%d', max_lag, n_bins);
    end
    n_samples_lagged = n_trials * (n_bins - max_lag);
    if n_samples_lagged < max(k_x, k_y)
        error('v5_lagged_refit_ifi:badShape', ...
            'lagged sample count %d < max(k_x=%d, k_y=%d)', ...
            n_samples_lagged, k_x, k_y);
    end

    cc_xy = zeros(max_lag + 1, 1);
    cc_yx = zeros(max_lag + 1, 1);

    % Zero-lag fit, shared between directions.
    X_flat = reshape(permute(X, [2 1 3]), n_trials * n_bins, k_x);
    Y_flat = reshape(permute(Y, [2 1 3]), n_trials * n_bins, k_y);
    res0 = v5_cca_fit(X_flat, Y_flat);
    cc_xy(1) = res0.cc(1);
    cc_yx(1) = res0.cc(1);

    for tau = 1:max_lag
        % X-leads-Y: X(:, 1:end-tau, :) paired with Y(:, tau+1:end, :)
        X_part = X(:, 1:end-tau, :);
        Y_part = Y(:, tau+1:end, :);
        X_lag = reshape(permute(X_part, [2 1 3]), n_trials * (n_bins - tau), k_x);
        Y_lag = reshape(permute(Y_part, [2 1 3]), n_trials * (n_bins - tau), k_y);
        r_xy = v5_cca_fit(X_lag, Y_lag);
        cc_xy(tau + 1) = r_xy.cc(1);
        % Y-leads-X: Y(:, 1:end-tau, :) paired with X(:, tau+1:end, :)
        Y_part2 = Y(:, 1:end-tau, :);
        X_part2 = X(:, tau+1:end, :);
        Y_lag2 = reshape(permute(Y_part2, [2 1 3]), n_trials * (n_bins - tau), k_y);
        X_lag2 = reshape(permute(X_part2, [2 1 3]), n_trials * (n_bins - tau), k_x);
        r_yx = v5_cca_fit(X_lag2, Y_lag2);
        cc_yx(tau + 1) = r_yx.cc(1);
    end

    kernel = gaussian_kernel(smooth_size, smooth_sigma);
    cc_xy_smooth = gaussian_smooth(cc_xy, kernel);
    cc_yx_smooth = gaussian_smooth(cc_yx, kernel);

    cw = min(central_window, max_lag);
    % Indices: cc_xy(2:cw+1) are taus 1..cw (since index 1 is tau=0).
    xy_mean = mean(cc_xy_smooth(2:cw + 1));
    yx_mean = mean(cc_yx_smooth(2:cw + 1));
    denom = xy_mean + yx_mean;
    if denom > 1e-15
        ifi = (xy_mean - yx_mean) / denom;
    else
        ifi = 0;
    end

    res = struct( ...
        'cc_xy', cc_xy, ...
        'cc_yx', cc_yx, ...
        'cc_xy_smooth', cc_xy_smooth, ...
        'cc_yx_smooth', cc_yx_smooth, ...
        'ifi', ifi);
end

function k = gaussian_kernel(size_in, sigma)
% 1D Gaussian kernel normalised to sum to 1.
    if size_in < 1 || mod(size_in, 2) ~= 1
        error('gaussian_kernel:badSize', 'size must be a positive odd integer');
    end
    half = (size_in - 1) / 2;
    x = (-half:half).';
    k = exp(-(x .^ 2) / (2 * sigma ^ 2));
    k = k / sum(k);
end

function out = gaussian_smooth(y, kernel)
% Convolve y with kernel, preserving length via reflect padding (matches
% numpy's behavior in the Python reference).
    y = y(:);
    n = numel(y);
    half = (numel(kernel) - 1) / 2;
    % Reflect padding: y(half+1:-1:2) at start, y(end-1:-1:end-half) at end
    if n >= half + 1
        pad_l = y(half+1:-1:2);
        pad_r = y(end-1:-1:end-half);
    else
        pad_l = zeros(half, 1);
        pad_r = zeros(half, 1);
    end
    padded = [pad_l; y; pad_r];
    full = conv(padded, kernel, 'valid');
    out = full(1:n);
end
