function P = v5_pca_apply(S, state)
%V5_PCA_APPLY  Project new data through an existing PCA fit.
%
%   P = V5_PCA_APPLY(S, STATE) applies the mean+components in STATE
%   (returned by V5_PCA_REDUCE) to a new (held-out) data tensor S.
%
%   Inputs
%       S       (n_trials × n_bins × n_units) double. Must have the
%               same n_units as the data used to fit STATE.
%       state   struct with fields `mean`, `components`,
%               `explained_variance_ratio`. Returned by v5_pca_reduce.
%
%   Output
%       P       (n_trials × n_bins × K) double, where K = size(state.components, 1).
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_pca_apply.m
%   Python reference at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::pca_apply

    sz = size(S);
    if numel(sz) ~= 3
        error('v5_pca_apply:badShape', ...
            'S must be 3D (n_trials × n_bins × n_units); got size [%s]', ...
            num2str(sz));
    end
    n_trials = sz(1);
    n_bins   = sz(2);
    n_units  = sz(3);

    if n_units ~= size(state.mean, 2)
        error('v5_pca_apply:unitMismatch', ...
            'S has %d units; state was fit on %d', ...
            n_units, size(state.mean, 2));
    end

    k = size(state.components, 1);
    X = reshape(double(S), n_trials * n_bins, n_units);
    Xc = X - state.mean;
    P_flat = Xc * state.components.';            % (n_samples × k)
    P = reshape(P_flat, n_trials, n_bins, k);
end
