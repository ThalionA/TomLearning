function test_v5_cca_fit()
%TEST_V5_CCA_FIT  Synthetic-data tests for v5_cca_fit and v5_cca_fit_cv.
%
%   Mirrors the Python tests in
%       /Users/theoamvr/Desktop/Experiments/IBL/tests/test_cca_core.py::TestCCAFit
%
%   Cross-language reference values are hardcoded from a Python run of
%   core.cca_fit(X, Y) and core.cca_fit_cv(X, Y, n_splits=4, shuffle=false).

    fprintf('=== v5_cca_fit / v5_cca_fit_cv tests ===\n');
    results = {};
    results(end+1, :) = run_case('identity_yields_unit_correlations', @t_identity);
    results(end+1, :) = run_case('coupled_latent',                    @t_coupled_latent);
    results(end+1, :) = run_case('independent_populations_near_zero', @t_independent);
    results(end+1, :) = run_case('sign_convention_a_max_positive',    @t_sign_convention);
    results(end+1, :) = run_case('cc_matches_manual_corr',            @t_cc_matches_manual);
    results(end+1, :) = run_case('cv_shapes_and_coupled',             @t_cv_shapes);
    results(end+1, :) = run_case('cv_independent_near_zero',          @t_cv_independent);
    results(end+1, :) = run_case('cv_no_shuffle_deterministic',       @t_cv_no_shuffle);
    results(end+1, :) = run_case('shape_mismatch_raises',             @t_shape_mismatch);
    results(end+1, :) = run_case('cross_language_reference_case',     @t_cross_lang);
    results(end+1, :) = run_case('finite_inputs_unchanged_by_guard',  @t_finite_unchanged);

    n_pass = sum(cellfun(@(x) x, results(:, 2)));
    n_total = size(results, 1);
    fprintf('\n%d/%d passed\n', n_pass, n_total);
    if n_pass < n_total
        error('test_v5_cca_fit:fail', 'Some tests failed; see output above.');
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
    rng_state = rng(0, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(50, 3);
    res = v5_cca_fit(X, X);
    assert(numel(res.cc) == 3);
    assert(max(abs(res.cc - 1)) < 1e-10, 'identity CCs not ~1: %s', mat2str(res.cc));
end

function t_coupled_latent()
    rng_state = rng(1, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n = 200;
    L = randn(n, 2);
    noise = 0.2;
    X = [L(:, 1) + noise * randn(n, 1), ...
         L(:, 2) + noise * randn(n, 1), ...
         randn(n, 1), randn(n, 1)];
    Y = [L(:, 1) + noise * randn(n, 1), ...
         L(:, 2) + noise * randn(n, 1), ...
         randn(n, 1), randn(n, 1)];
    res = v5_cca_fit(X, Y);
    assert(res.cc(1) > 0.85, 'cc1: %g', res.cc(1));   % MATLAB rng ≠ numpy rng so margins
    assert(res.cc(2) > 0.7,  'cc2: %g', res.cc(2));   % are looser than Python.
    assert(res.cc(3) < 0.5,  'cc3: %g', res.cc(3));
end

function t_independent()
    rng_state = rng(2, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(200, 4);
    Y = randn(200, 4);
    res = v5_cca_fit(X, Y);
    assert(res.cc(1) < 0.3, 'independent cc1: %g', res.cc(1));
end

function t_sign_convention()
    rng_state = rng(3, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(100, 4);
    Y = X * randn(4, 4) + 0.1 * randn(100, 4);
    res = v5_cca_fit(X, Y);
    for i = 1:size(res.A, 2)
        [~, j] = max(abs(res.A(:, i)));
        assert(res.A(j, i) >= 0, 'col %d max-abs negative', i);
    end
end

function t_cc_matches_manual()
    rng_state = rng(4, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(100, 4);
    Y = randn(100, 4);
    res = v5_cca_fit(X, Y);
    U = (X - res.x_mean) * res.A;
    V = (Y - res.y_mean) * res.B;
    for i = 1:numel(res.cc)
        c = corrcoef(U(:, i), V(:, i));
        assert(abs(c(1,2) - res.cc(i)) < 1e-10, 'col %d', i);
    end
end

function t_cv_shapes()
    rng_state = rng(5, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(100, 3);
    Y = X + 0.5 * randn(100, 3);
    res = v5_cca_fit_cv(X, Y, 10, true, 0);
    assert(isequal(size(res.A), [3 3]));
    assert(numel(res.cc_cv) == 3);
    assert(isequal(size(res.cc_cv_per_fold), [10 3]));
    assert(res.cc_cv(1) > 0.7, 'cc_cv(1): %g', res.cc_cv(1));
end

function t_cv_independent()
    rng_state = rng(6, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(200, 3);
    Y = randn(200, 3);
    res = v5_cca_fit_cv(X, Y, 10, true, 0);
    assert(abs(res.cc_cv(1)) < 0.2, 'cv cc1: %g', res.cc_cv(1));
end

function t_cv_no_shuffle()
    rng_state = rng(7, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(100, 3);
    Y = randn(100, 3);
    res1 = v5_cca_fit_cv(X, Y, 5, false);
    res2 = v5_cca_fit_cv(X, Y, 5, false);
    assert(isequal(res1.cc_cv_per_fold, res2.cc_cv_per_fold));
end

function t_shape_mismatch()
    raised = false;
    try
        v5_cca_fit(zeros(10, 2), zeros(9, 2));
    catch err
        raised = strcmp(err.identifier, 'v5_cca_fit:badShape');
    end
    assert(raised);
end

function t_finite_unchanged()
    % REGRESSION: the NaN guard in _cca_sign_fix must be a no-op for
    % finite inputs. The cross-language reference test above already
    % covers byte-identical numerics; this redundantly asserts the guard
    % path doesn't perturb a well-conditioned fit so the cross-language
    % invariant is preserved for any new test author who modifies the
    % guard.
    rng_state = rng(11, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(80, 4);
    Y = X * randn(4, 4) + 0.2 * randn(80, 4);
    res = v5_cca_fit(X, Y);
    for i = 1:size(res.A, 2)
        col = res.A(:, i);
        [~, j] = max(abs(col));
        assert(col(j) >= 0, 'sign convention violated on col %d', i);
        assert(all(isfinite(col)), 'finite input produced non-finite A col %d', i);
    end
end

function t_cross_lang()
    % CROSS-LANG: deterministic 8×2 X, Y matrices. Hardcoded outputs
    % match Python's core.cca_fit(X, Y) and core.cca_fit_cv(X, Y, n_splits=4, shuffle=false).
    X = [ 1.0,  0.5; -1.0,  0.5; ...
          2.0,  1.0; -2.0,  1.0; ...
          0.5, -0.5; -0.5, -0.5; ...
          1.5, -1.0; -1.5, -1.0];
    Y = [ 2.0,  0.3; -2.0,  0.4; ...
          4.0,  0.2; -4.0,  0.5; ...
          1.0, -0.1; -1.0, -0.2; ...
          3.0, -0.3; -3.0, -0.4];

    res = v5_cca_fit(X, Y);

    expected_cc = [1.0; 0.9480269968318369];
    expected_A = [ 0.2581988897471611,  -1.7753745435194553e-17; ...
                   8.545928625897863e-18, 0.447213595499958];
    expected_B = [ 0.12909944487358055,   0.018595200082640036; ...
                   0.0,                   1.1157120049584026];
    expected_x_mean = [0.0, 0.0];
    expected_y_mean = [0.0, 0.05];

    assert(max(abs(res.cc - expected_cc)) < 1e-10, 'cc mismatch: %s', mat2str(res.cc));
    assert(max(abs(res.A(:) - expected_A(:))) < 1e-10, 'A mismatch');
    assert(max(abs(res.B(:) - expected_B(:))) < 1e-10, 'B mismatch');
    assert(max(abs(res.x_mean - expected_x_mean)) < 1e-12, 'x_mean mismatch');
    assert(max(abs(res.y_mean - expected_y_mean)) < 1e-12, 'y_mean mismatch');

    % CV without shuffle: 4 folds of size 2 each. CC2 is 0 in every fold
    % because the 2-sample test set has zero variance in the second
    % canonical direction → guarded to 0.
    res_cv = v5_cca_fit_cv(X, Y, 4, false);
    expected_cv_per_fold = [1.0, 0.0; 1.0, 0.0; 1.0, 0.0; 1.0, 0.0];
    expected_cc_cv = [1.0; 0.0];
    assert(max(abs(res_cv.cc_cv_per_fold(:) - expected_cv_per_fold(:))) < 1e-10, ...
        'cc_cv_per_fold mismatch: %s', mat2str(res_cv.cc_cv_per_fold));
    assert(max(abs(res_cv.cc_cv - expected_cc_cv)) < 1e-10, 'cc_cv mismatch');
end
