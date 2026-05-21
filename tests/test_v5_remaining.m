function test_v5_remaining()
%TEST_V5_REMAINING  Tests for v5_project, v5_lagged_refit_ifi,
%                   v5_projection_ifi, v5_principal_angles.
%
%   Mirrors the Python tests in
%       /Users/theoamvr/Desktop/Experiments/IBL/tests/test_cca_core.py
%   for the same primitives. Cross-language reference values are
%   intentionally degenerate (perfect coupling or known geometry)
%   since the per-primitive numerics are already underwritten by the
%   stage-4 cca_fit cross-lang test.

    fprintf('=== v5_project / lagged_refit_ifi / projection_ifi / principal_angles ===\n');
    results = {};

    % project
    results(end+1, :) = run_case('project_basic',                @t_project_basic);
    results(end+1, :) = run_case('project_shape_mismatch',       @t_project_shape_mismatch);

    % lagged_refit_ifi
    results(end+1, :) = run_case('lagged_zero_simultaneous',     @t_lagged_zero);
    results(end+1, :) = run_case('lagged_positive_when_x_leads', @t_lagged_x_leads);
    results(end+1, :) = run_case('lagged_negative_when_y_leads', @t_lagged_y_leads);
    results(end+1, :) = run_case('lagged_validation',            @t_lagged_validation);

    % projection_ifi
    results(end+1, :) = run_case('projection_zero_when_symmetric',  @t_proj_zero);
    results(end+1, :) = run_case('projection_positive_when_u_leads',@t_proj_u_leads);
    results(end+1, :) = run_case('projection_negative_when_v_leads',@t_proj_v_leads);
    results(end+1, :) = run_case('projection_shapes',               @t_proj_shapes);

    % principal_angles
    results(end+1, :) = run_case('pa_identical_subspaces',  @t_pa_identical);
    results(end+1, :) = run_case('pa_orthogonal_subspaces', @t_pa_orthogonal);
    results(end+1, :) = run_case('pa_known_rotation',       @t_pa_rotation);
    results(end+1, :) = run_case('pa_validation',           @t_pa_validation);
    results(end+1, :) = run_case('pa_cross_lang_reference', @t_pa_cross_lang);
    results(end+1, :) = run_case('pa_degenerate_zeros',     @t_pa_degenerate_zeros);
    results(end+1, :) = run_case('pa_degenerate_nan_inf',   @t_pa_degenerate_nan_inf);

    n_pass = sum(cellfun(@(x) x, results(:, 2)));
    n_total = size(results, 1);
    fprintf('\n%d/%d passed\n', n_pass, n_total);
    if n_pass < n_total
        error('test_v5_remaining:fail', 'Some tests failed; see output above.');
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
function t_project_basic()
    % P = reshape(0..23, [2 3 4]), mean = [1 2 3 4], loadings = eye(4)(:, 1:2)
    % Q should equal (P - mean)(..., 1:2)
    P = reshape(0:23, [4 3 2]);  P = permute(P, [3 2 1]);  % numpy C-order layout
    mu = [1.0, 2.0, 3.0, 4.0];
    loadings = eye(4);
    loadings = loadings(:, 1:2);
    Q = v5_project(P, mu, loadings);
    assert(isequal(size(Q), [2 3 2]));
    expected = P - reshape(mu, [1 1 4]);
    expected = expected(:, :, 1:2);
    assert(max(abs(Q(:) - expected(:))) < 1e-12);
end

function t_project_shape_mismatch()
    raised = false;
    try
        v5_project(zeros(2, 4), zeros(3, 1), eye(4));
    catch err
        raised = strcmp(err.identifier, 'v5_project:badShape');
    end
    assert(raised);
end

function t_lagged_zero()
    rng_state = rng(0, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n_trials = 30; n_bins = 30;
    L = randn(n_trials, n_bins, 2);
    X = L + 0.05 * randn(n_trials, n_bins, 2);
    Y = L + 0.05 * randn(n_trials, n_bins, 2);
    res = v5_lagged_refit_ifi(X, Y, 5);
    fprintf('    ifi (simultaneous): %.4f\n', res.ifi);
    assert(abs(res.ifi) < 0.05, 'ifi: %g', res.ifi);
end

function t_lagged_x_leads()
    rng_state = rng(1, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n_trials = 40; n_bins = 30; lead = 2;
    X = randn(n_trials, n_bins, 2);
    Y = zeros(size(X));
    Y(:, lead+1:end, :) = X(:, 1:end-lead, :);
    Y = Y + 0.1 * randn(size(Y));
    res = v5_lagged_refit_ifi(X, Y, 5);
    fprintf('    ifi (X leads): %.4f\n', res.ifi);
    assert(res.ifi > 0.3);
end

function t_lagged_y_leads()
    rng_state = rng(2, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n_trials = 40; n_bins = 30; lead = 2;
    Y = randn(n_trials, n_bins, 2);
    X = zeros(size(Y));
    X(:, lead+1:end, :) = Y(:, 1:end-lead, :);
    X = X + 0.1 * randn(size(X));
    res = v5_lagged_refit_ifi(X, Y, 5);
    fprintf('    ifi (Y leads): %.4f\n', res.ifi);
    assert(res.ifi < -0.3);
end

function t_lagged_validation()
    raised = false;
    try
        v5_lagged_refit_ifi(zeros(2, 5, 3), zeros(2, 4, 3), 2);
    catch err
        raised = strcmp(err.identifier, 'v5_lagged_refit_ifi:badShape');
    end
    assert(raised);
end

function t_proj_zero()
    rng_state = rng(3, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    u = randn(10, 30);
    v = u;
    [ifi, ~] = v5_projection_ifi(u, v, 5);
    finite_ifi = ifi(~isnan(ifi));
    fprintf('    ifi (u=v): max abs = %.4f\n', max(abs(finite_ifi)));
    assert(max(abs(finite_ifi)) < 0.05);
end

function t_proj_u_leads()
    rng_state = rng(4, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    u = randn(20, 30);
    v = zeros(size(u));
    v(:, 3:end) = u(:, 1:end-2);
    v = v + 0.05 * randn(size(v));
    [ifi, ~] = v5_projection_ifi(u, v, 5);
    finite_ifi = ifi(~isnan(ifi));
    fprintf('    ifi (u leads): mean = %.4f\n', mean(finite_ifi));
    assert(mean(finite_ifi) > 0.3);
end

function t_proj_v_leads()
    rng_state = rng(5, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    v = randn(20, 30);
    u = zeros(size(v));
    u(:, 3:end) = v(:, 1:end-2);
    u = u + 0.05 * randn(size(u));
    [ifi, ~] = v5_projection_ifi(u, v, 5);
    finite_ifi = ifi(~isnan(ifi));
    fprintf('    ifi (v leads): mean = %.4f\n', mean(finite_ifi));
    assert(mean(finite_ifi) < -0.3);
end

function t_proj_shapes()
    rng_state = rng(6, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    u = randn(5, 20);
    v = randn(5, 20);
    [ifi, lag_corr] = v5_projection_ifi(u, v, 4);
    assert(numel(ifi) == 5);
    assert(isequal(size(lag_corr), [5, 9]));
end

function t_pa_identical()
    rng_state = rng(8, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    A = randn(10, 3);
    angles = v5_principal_angles(A, A);
    assert(all(angles < 1e-7), 'identical-subspace angles should be ~0; max=%g', max(angles));
end

function t_pa_orthogonal()
    A1 = zeros(6, 2);  A1(1, 1) = 1;  A1(2, 2) = 1;
    A2 = zeros(6, 2);  A2(3, 1) = 1;  A2(4, 2) = 1;
    angles = v5_principal_angles(A1, A2);
    assert(max(abs(angles - pi/2)) < 1e-10);
end

function t_pa_rotation()
    A1 = [1.0; 0.0; 0.0];
    theta = pi / 6;
    A2 = [cos(theta); sin(theta); 0.0];
    angles = v5_principal_angles(A1, A2);
    assert(abs(angles(1) - theta) < 1e-10);
end

function t_pa_validation()
    raised = false;
    try
        v5_principal_angles(zeros(5, 2), zeros(4, 2));
    catch err
        raised = strcmp(err.identifier, 'v5_principal_angles:badShape');
    end
    assert(raised);
end

function t_pa_cross_lang()
    % CROSS-LANG: deterministic geometry. Two 2D subspaces in R^3 with a
    % known 45° rotation between their second basis vectors. Expected
    % angles: [0, pi/4] (first basis vectors identical, second 45° apart).
    A1 = [1.0, 0.0; 0.0, 1.0; 0.0, 0.0];
    A2 = [1.0, 0.0; 0.0, cos(pi/4); 0.0, sin(pi/4)];
    angles = v5_principal_angles(A1, A2);
    expected = [0.0; pi/4];
    assert(max(abs(angles - expected)) < 1e-10, ...
        'angles: %s, expected: %s', mat2str(angles), mat2str(expected));
end

function t_pa_degenerate_zeros()
    % All-zero columns should be dropped; with nothing left, return NaN.
    A1 = zeros(5, 3);
    A2 = randn(5, 3);
    angles = v5_principal_angles(A1, A2);
    assert(isscalar(angles) && isnan(angles), ...
        'all-zero A1 should yield scalar NaN');
    % Same for both sides zero.
    angles2 = v5_principal_angles(zeros(4, 2), zeros(4, 2));
    assert(isscalar(angles2) && isnan(angles2));
end

function t_pa_degenerate_nan_inf()
    % Columns with NaN/Inf entries (from ill-conditioned Rx\U) should be
    % dropped. If some columns survive, angles are computed on the clean subset.
    A1 = [1 NaN; 0 Inf; 0 0; 0 0];     % col-2 is degenerate
    A2 = [1 0;   0 1;   0 0; 0 0];     % both cols clean
    angles = v5_principal_angles(A1, A2);
    % Only col-1 of A1 survives → compare 1D subspaces.
    % A1_clean = [1;0;0;0], A2_clean = [1 0; 0 1; 0 0; 0 0]
    % Principal angle between span([1;0;0;0]) and span([1,0;0,1;0,0;0,0])
    % The first col of A2 is identical → angle should be 0.
    assert(~any(isnan(angles)), 'clean columns should yield finite angles');
    assert(angles(1) < 1e-10, 'first angle should be ~0; got %g', angles(1));
end
