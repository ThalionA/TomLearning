function W_cells = v5_cell_weights(A_pc, pca_state)
%V5_CELL_WEIGHTS  Back-project PC-space CCA loadings to cell space.
%
%   W_CELLS = V5_CELL_WEIGHTS(A_PC, PCA_STATE) composes the PCA basis
%   stored in PCA_STATE with the PC-space CCA loadings A_PC to produce
%   per-cell weight vectors in the original neuron-by-neuron coordinate
%   system. This is what makes principal-angle comparisons between
%   different fits meaningful: PCA bases differ across fits (each fit
%   has its own mean and components), but cell space is shared.
%
%   Construction
%       The CCA canonical variate is
%           u(t) = (P_x(t) - x_mean_pc) * A_pc
%       where P_x(t) is the per-sample PC-score vector. Substituting
%       P_x(t) = (X(t) - mu_pca) * components.' gives
%           u(t) = (X(t) - mu_combined) * (components.' * A_pc)
%       so the cell-space weight vector is `components.' * A_pc`. The
%       PCA mean shows up as a per-canonical-direction offset, which
%       does not affect column-space comparisons (principal angles,
%       directionality, weight-magnitude distributions).
%
%   Inputs
%       A_pc       (k_pc x r) double. CCA loadings in PC space; r is
%                  the number of canonical pairs.
%       pca_state  struct returned by v5_pca_reduce with field
%                  `components` of shape (k_pc x n_units) — rows are PCs.
%
%   Output
%       W_cells    (n_units x r) double. Per-cell weight for each
%                  canonical dimension.
%
%   Notes
%       - Used in the v5 spatial / temporal drivers to compare epoch-
%         to-epoch subspace rotations via v5_principal_angles. Comparing
%         A_pc directly is meaningless when k_pc equals the fitted rank,
%         because both bases trivially span the entire PC ambient.
%       - Cross-project mirror: this helper is identical across the
%         Hippocampus-V1 (TomLearning) and Striatum drivers. The Python
%         equivalent is `components.T @ A_pc` written inline in IBL's
%         src/cca/pipelines.py.

    if ~ismatrix(A_pc)
        error('v5_cell_weights:badShape', 'A_pc must be 2D');
    end
    if ~isstruct(pca_state) || ~isfield(pca_state, 'components')
        error('v5_cell_weights:badShape', ...
            'pca_state must be a struct with a `components` field');
    end
    components = pca_state.components;
    if ~ismatrix(components)
        error('v5_cell_weights:badShape', ...
            'pca_state.components must be 2D (k_pc x n_units)');
    end
    if size(components, 1) ~= size(A_pc, 1)
        error('v5_cell_weights:badShape', ...
            'components has %d PCs; A_pc has %d rows', ...
            size(components, 1), size(A_pc, 1));
    end

    W_cells = components.' * A_pc;
end
