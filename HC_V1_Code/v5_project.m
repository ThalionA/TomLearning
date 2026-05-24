function Q = v5_project(P, mu, loadings)
%V5_PROJECT  Apply (P - mu) * loadings, preserving leading dims of P.
%
%   Q = V5_PROJECT(P, MU, LOADINGS) takes P of any shape ending in
%   k_in features, subtracts the per-feature mean MU (length k_in),
%   and applies the LOADINGS matrix (k_in x k_out). The output Q has
%   the same leading dimensions as P, with the last dim replaced by k_out.
%
%   Stage 6 of the H&H-adapted CCA pipeline. The canonical use is
%   projecting per-trial PC scores through CCA loadings to get per-
%   trial canonical variates: u(t) = (P_X - x_mean) * A.
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_project.m
%   Python reference at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::project

    sz = size(P);
    if numel(sz) < 1
        error('v5_project:badShape', 'P must have at least 1 dim');
    end
    if ~isvector(mu)
        error('v5_project:badShape', 'mu must be 1D');
    end
    if ~ismatrix(loadings)
        error('v5_project:badShape', 'loadings must be 2D');
    end
    k_in = sz(end);
    if numel(mu) ~= k_in
        error('v5_project:badShape', ...
            'mu has %d features; P has %d', numel(mu), k_in);
    end
    if size(loadings, 1) ~= k_in
        error('v5_project:badShape', ...
            'loadings has %d rows; P has %d features', size(loadings, 1), k_in);
    end
    k_out = size(loadings, 2);

    n_samples = prod(sz(1:end-1));
    flat = reshape(double(P), n_samples, k_in);
    flat_q = (flat - mu(:).') * loadings;
    out_shape = [sz(1:end-1), k_out];
    if numel(out_shape) == 1
        out_shape = [out_shape, 1];   % MATLAB needs 2D minimum for reshape
    end
    Q = reshape(flat_q, out_shape);
end
