function angles = v5_principal_angles(A1, A2)
%V5_PRINCIPAL_ANGLES  Principal angles between two subspaces.
%
%   ANGLES = V5_PRINCIPAL_ANGLES(A1, A2) returns the principal angles
%   (radians, in [0, pi/2]) between the subspaces spanned by the
%   columns of A1 and A2.
%
%   - All zero => identical subspaces.
%   - All pi/2 => fully orthogonal subspaces.
%
%   Implementation: QR-orthonormalise both bases, then SVD the cross-
%   product Q1' * Q2; the singular values are cosines of the principal
%   angles. Clipped to [-1, 1] before arccos to absorb numerical
%   overshoots.
%
%   Degenerate columns (non-finite or near-zero norm) are silently
%   dropped before the QR step.  This handles the common case where
%   upstream CCA loadings from an ill-conditioned Rx\U solve contain
%   Inf/NaN entries.  If *all* columns are degenerate on either side,
%   a scalar NaN is returned.
%
%   Diagnostic per spec §5 — used to quantify how much the canonical
%   subspace itself rotates across paired fits (e.g. naive vs post-LP).
%
%   Cross-project mirror: must stay byte-identical with
%       /Users/theoamvr/Desktop/Experiments/StriatumACC/Striatum project/v5_principal_angles.m
%   Python reference at
%       /Users/theoamvr/Desktop/Experiments/IBL/src/cca/core.py::principal_angles

    if ~ismatrix(A1) || ~ismatrix(A2)
        error('v5_principal_angles:badShape', 'A1, A2 must be 2D');
    end
    if size(A1, 1) ~= size(A2, 1)
        error('v5_principal_angles:badShape', ...
            'A1, A2 must have same row count (ambient dim); got %d vs %d', ...
            size(A1, 1), size(A2, 1));
    end

    % Guard: drop columns that are non-finite or effectively zero-norm.
    A1 = sanitise_basis(A1);
    A2 = sanitise_basis(A2);
    if isempty(A1) || isempty(A2)
        angles = NaN;
        return;
    end

    [Q1, ~] = qr(A1, 0);
    [Q2, ~] = qr(A2, 0);
    s = svd(Q1.' * Q2);
    s = max(-1, min(1, s));
    angles = acos(s);
end


function B = sanitise_basis(A)
%SANITISE_BASIS  Remove columns with non-finite or near-zero norm.
    norms = vecnorm(A);
    keep = isfinite(norms) & (norms > 1e-12);
    B = A(:, keep);
end
