function test_v5_pca_reduce()
%TEST_V5_PCA_REDUCE  Synthetic-data tests for v5_pca_reduce and v5_pca_apply.
%
%   Mirrors the Python tests in
%       /Users/theoamvr/Desktop/Experiments/IBL/tests/test_cca_core.py::TestPCAReduce
%
%   Cross-language reference values are hardcoded from a Python run of
%   src.cca.core.pca_reduce on the same deterministic 6×1×4 input. Both
%   languages should agree to within 1e-10.

    fprintf('=== v5_pca_reduce / v5_pca_apply tests ===\n');
    results = {};
    results(end+1, :) = run_case('shape_and_evr',                       @t_shape_and_evr);
    results(end+1, :) = run_case('against_manual_svd',                  @t_against_manual_svd);
    results(end+1, :) = run_case('k_above_rank_trails_to_zero',         @t_k_above_rank);
    results(end+1, :) = run_case('pca_apply_round_trip',                @t_apply_round_trip);
    results(end+1, :) = run_case('pca_apply_to_new_data',               @t_apply_new_data);
    results(end+1, :) = run_case('pca_apply_unit_count_mismatch_raises', @t_apply_unit_mismatch);
    results(end+1, :) = run_case('pca_reduce_bad_k_raises',             @t_bad_k);
    results(end+1, :) = run_case('input_unchanged',                     @t_input_unchanged);
    results(end+1, :) = run_case('cross_language_reference_case',       @t_cross_lang_reference);

    n_pass = sum(cellfun(@(x) x, results(:, 2)));
    n_total = size(results, 1);
    fprintf('\n%d/%d passed\n', n_pass, n_total);
    if n_pass < n_total
        error('test_v5_pca_reduce:fail', 'Some tests failed; see output above.');
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
function t_shape_and_evr()
    rng_state = rng(0, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    S = randn(20, 8, 5);
    [P, state] = v5_pca_reduce(S, 3);
    assert(isequal(size(P), [20 8 3]));
    assert(isequal(size(state.components), [3 5]));
    assert(isequal(size(state.mean), [1 5]));
    assert(numel(state.explained_variance_ratio) == 3);
    s = sum(state.explained_variance_ratio);
    assert(s > 0 && s <= 1 + 1e-12);
    assert(all(state.explained_variance_ratio >= 0));
    assert(all(diff(state.explained_variance_ratio) <= 1e-12));
end

function t_against_manual_svd()
    rng_state = rng(1, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    S = randn(15, 10, 4);
    X = reshape(S, [], 4);
    Xc = X - mean(X, 1);
    [U, Sigma, V] = svd(Xc, 'econ');
    sigma = diag(Sigma);
    % Same svd_flip as the primitive.
    r = size(V, 2);
    signs = ones(1, r);
    for i = 1:r
        [~, j] = max(abs(V(:, i)));
        s = sign(V(j, i));
        if s == 0, s = 1; end
        signs(i) = s;
    end
    V = V .* signs; U = U .* signs;
    P_expected = reshape(U(:, 1:3) .* sigma(1:3).', 15, 10, 3);
    evr_expected = (sigma(1:3) .^ 2) / sum(sigma .^ 2);
    [P, state] = v5_pca_reduce(S, 3);
    assert(max(abs(P(:) - P_expected(:))) < 1e-12);
    assert(max(abs(state.explained_variance_ratio - evr_expected)) < 1e-12);
end

function t_k_above_rank()
    rng_state = rng(2, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    latent = randn(30, 2);
    loadings = randn(2, 5);
    X = latent * loadings;
    S = reshape(X, 6, 5, 5);
    [P, state] = v5_pca_reduce(S, 3);
    assert(sum(state.explained_variance_ratio(1:2)) > 1 - 1e-10);
    assert(state.explained_variance_ratio(3) < 1e-10);
    assert(max(abs(reshape(P(:, :, 3), [], 1))) < 1e-8);
end

function t_apply_round_trip()
    rng_state = rng(3, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    S = randn(25, 6, 4);
    [P_fit, state] = v5_pca_reduce(S, 2);
    P_apply = v5_pca_apply(S, state);
    assert(max(abs(P_apply(:) - P_fit(:))) < 1e-12);
end

function t_apply_new_data()
    rng_state = rng(4, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    S_train = randn(20, 5, 4);
    S_test = randn(10, 5, 4);
    [~, state] = v5_pca_reduce(S_train, 3);
    P_test = v5_pca_apply(S_test, state);
    X_test = reshape(S_test, [], 4);
    expected = reshape((X_test - state.mean) * state.components.', 10, 5, 3);
    assert(max(abs(P_test(:) - expected(:))) < 1e-12);
end

function t_apply_unit_mismatch()
    rng_state = rng(5, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    S_train = randn(10, 4, 5);
    S_bad = randn(10, 4, 7);
    [~, state] = v5_pca_reduce(S_train, 2);
    raised = false;
    try
        v5_pca_apply(S_bad, state);
    catch err
        raised = strcmp(err.identifier, 'v5_pca_apply:unitMismatch');
    end
    assert(raised, 'expected v5_pca_apply:unitMismatch error');
end

function t_bad_k()
    raised = false;
    try
        v5_pca_reduce(zeros(10, 5, 4), 0);
    catch err
        raised = strcmp(err.identifier, 'v5_pca_reduce:badK');
    end
    assert(raised, 'expected badK for k=0');
    raised = false;
    try
        v5_pca_reduce(zeros(10, 5, 4), 5);
    catch err
        raised = strcmp(err.identifier, 'v5_pca_reduce:badK');
    end
    assert(raised, 'expected badK for k > n_units');
end

function t_input_unchanged()
    rng_state = rng(6, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    S = randn(15, 4, 6);
    S_orig = S;
    [~, ~] = v5_pca_reduce(S, 3);
    assert(isequal(S, S_orig));
end

function t_cross_lang_reference()
    % CROSS-LANG: deterministic 6×1×4 input, expected outputs hardcoded
    % from a Python run of src.cca.core.pca_reduce(S, k=2).
    % Tolerance 1e-10 absorbs LAPACK implementation differences.
    S = zeros(6, 1, 4);
    S(1, 1, :) = [ 1.0,  2.0,  0.5, -1.0];
    S(2, 1, :) = [ 2.0,  4.0,  1.5,  0.0];
    S(3, 1, :) = [ 3.0,  6.0,  0.5,  2.0];
    S(4, 1, :) = [ 4.0,  8.0,  2.5, -3.0];
    S(5, 1, :) = [ 5.0, 10.0,  1.5,  4.0];
    S(6, 1, :) = [ 6.0, 12.0,  3.5, -2.0];
    [P, state] = v5_pca_reduce(S, 2);

    expected_mean = [3.5, 7.0, 1.6666666666666667, 0.0];
    expected_evr = [0.7149918629210316; 0.2809271997111383];
    expected_components_row1 = [0.43647285879065656, 0.8729457175813128, ...
                                0.21771313333807601, 0.0076294998312245715];
    expected_components_row2 = [0.020188404092723603, 0.0403768081854472, ...
                                -0.2363831965828122, 0.9706106974859882];
    % P flat (C-order, varying last axis fastest): trial-by-trial, then PCs.
    % Shape is (6, 1, 2) so flat-C-order = [t0_pc0, t0_pc1, t1_pc0, t1_pc1, ...]
    expected_P_flat_C = [...
        -5.717538890275518, -0.9471853526317526, ...
        -3.309831963152937, -0.11201583126495801, ...
        -1.3299218028752808, 2.1665307807534475, ...
         1.249721258598031, -3.0583470793784984, ...
         3.2677789180318095, 4.073253020069847, ...
         5.839792479673897, -2.1222355375480864];
    % Reshape C-order to MATLAB layout: numpy (6,1,2) → MATLAB (2,1,6) reshape → permute(... ,[3 2 1])
    expected_P = permute(reshape(expected_P_flat_C, [2 1 6]), [3 2 1]);

    assert(max(abs(state.mean - expected_mean)) < 1e-12, ...
        'mean mismatch; max diff = %g', max(abs(state.mean - expected_mean)));
    assert(max(abs(state.explained_variance_ratio - expected_evr)) < 1e-12, ...
        'evr mismatch; max diff = %g', max(abs(state.explained_variance_ratio - expected_evr)));
    assert(max(abs(state.components(1, :) - expected_components_row1)) < 1e-10, ...
        'components row 1 mismatch');
    assert(max(abs(state.components(2, :) - expected_components_row2)) < 1e-10, ...
        'components row 2 mismatch');
    assert(max(abs(P(:) - expected_P(:))) < 1e-10, ...
        'P mismatch; max diff = %g', max(abs(P(:) - expected_P(:))));
end
