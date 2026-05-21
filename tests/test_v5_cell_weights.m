function test_v5_cell_weights()
%TEST_V5_CELL_WEIGHTS  Synthetic-data tests for v5_cell_weights.
%
%   Covers basic correctness (composition matches the algebra) plus the
%   regression case that motivated the helper: principal angles between
%   per-epoch CCA loadings must NOT trivially collapse to zero. Before
%   this fix, the driver passed (k x k) PC-space loadings into
%   v5_principal_angles; both bases trivially spanned the entire k-D
%   ambient and the diagnostic always returned ~0. The cell-space
%   composition recovers the actual rotation between two fits.

    fprintf('=== v5_cell_weights tests ===\n');
    results = {};
    results(end+1, :) = run_case('identity_pca_passthrough',     @t_identity);
    results(end+1, :) = run_case('shape_correctness',            @t_shape);
    results(end+1, :) = run_case('composition_matches_algebra',  @t_composition);
    results(end+1, :) = run_case('cell_space_principal_angles',  @t_principal_angles);
    results(end+1, :) = run_case('shape_mismatch_raises',        @t_shape_mismatch);

    n_pass = sum(cellfun(@(x) x, results(:, 2)));
    n_total = size(results, 1);
    fprintf('\n%d/%d passed\n', n_pass, n_total);
    if n_pass < n_total
        error('test_v5_cell_weights:fail', 'Some tests failed; see output above.');
    end
end

% -------------------------------------------------------------------- %
function row = run_case(name, fn)
    try
        fn();
        fprintf('  PASS: %s\n', name);
        row = {name, true};
    catch err
        fprintf('  FAIL: %s: %s\n', name, err.message);
        row = {name, false};
    end
end

% -------------------------------------------------------------------- %
function t_identity()
    % When components is the identity, cell weights equal A_pc.
    A_pc = randn(5, 3);
    state = struct('components', eye(5));
    W = v5_cell_weights(A_pc, state);
    assert(max(abs(W(:) - A_pc(:))) < 1e-12);
end

function t_shape()
    % (k_pc x r) loadings + (k_pc x n_units) components -> (n_units x r).
    A_pc = randn(4, 2);
    state = struct('components', randn(4, 7));
    W = v5_cell_weights(A_pc, state);
    assert(isequal(size(W), [7 2]));
end

function t_composition()
    % Direct check: W_cells == components.' * A_pc.
    rng_state = rng(0, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    A_pc = randn(6, 4);
    comps = randn(6, 10);
    state = struct('components', comps);
    W = v5_cell_weights(A_pc, state);
    assert(max(abs(W(:) - reshape(comps.' * A_pc, [], 1))) < 1e-12);
end

function t_principal_angles()
    % REGRESSION: principal angles between two CCA fits must be non-trivial
    % when computed in cell space. Build two epochs that share a known
    % shared latent in cell space but with different PCA bases per epoch.
    %
    % Setup. Two epochs, each with n_units = 12 cells. We construct a
    % rotation in cell space that swaps one shared latent component with a
    % new orthogonal direction. The principal angle between the two CCA
    % subspaces in CELL space should equal the planted rotation angle.
    rng_state = rng(7, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>

    n_units = 12;
    theta = pi / 6;                  % planted rotation

    % Cell-space loading for epoch 1: standard basis e1.
    W_cells_1 = zeros(n_units, 1);
    W_cells_1(1) = 1.0;

    % Cell-space loading for epoch 2: rotated within (e1, e2) plane by theta.
    W_cells_2 = zeros(n_units, 1);
    W_cells_2(1) = cos(theta);
    W_cells_2(2) = sin(theta);

    % Each epoch has its own PCA basis (arbitrary orthonormal). The A_pc
    % loading in each epoch is whatever maps the PC basis to its cell-space
    % weight: A_pc = components * W_cells (since components.' * A_pc = W_cells
    % when components has orthonormal rows).
    [Q1, ~] = qr(randn(n_units), 0);
    [Q2, ~] = qr(randn(n_units), 0);
    % Use first k_pc = 6 rows as a PCA basis (k_pc x n_units).
    k_pc = 6;
    comps_1 = Q1(:, 1:k_pc).';
    comps_2 = Q2(:, 1:k_pc).';
    state_1 = struct('components', comps_1);
    state_2 = struct('components', comps_2);
    A_pc_1 = comps_1 * W_cells_1;
    A_pc_2 = comps_2 * W_cells_2;

    % Round-trip: v5_cell_weights should recover W_cells_1 and W_cells_2.
    W_back_1 = v5_cell_weights(A_pc_1, state_1);
    W_back_2 = v5_cell_weights(A_pc_2, state_2);
    assert(max(abs(W_back_1 - W_cells_1)) < 1e-12, 'round-trip 1');
    assert(max(abs(W_back_2 - W_cells_2)) < 1e-12, 'round-trip 2');

    % Now the regression check: principal angles between A_pc_1 and A_pc_2
    % directly are MEANINGLESS, but in cell space they recover theta.
    angles_cell = v5_principal_angles(W_back_1, W_back_2);
    assert(abs(angles_cell(1) - theta) < 1e-10, ...
        'cell-space angle should be %g; got %g', theta, angles_cell(1));
    assert(angles_cell(1) > 0.1, ...
        'cell-space angle must not collapse to ~0');
end

function t_shape_mismatch()
    raised = false;
    try
        v5_cell_weights(zeros(3, 2), struct('components', zeros(4, 5)));
    catch err
        raised = strcmp(err.identifier, 'v5_cell_weights:badShape');
    end
    assert(raised, 'expected v5_cell_weights:badShape');
end
