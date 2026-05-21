function test_v5_shuffle_null()
%TEST_V5_SHUFFLE_NULL  Synthetic-data tests for v5_shuffle_null.
%
%   Mirrors the Python tests in
%       /Users/theoamvr/Desktop/Experiments/IBL/tests/test_cca_core.py
%   (the shuffle_null section to be added).
%
%   The seeded permutation path uses MATLAB's Mersenne Twister, which
%   produces different sequences than numpy.default_rng. So the
%   "deterministic given seed" property is tested WITHIN MATLAB (same
%   seed twice gives same answer) but the cross-language deterministic
%   test uses degenerate inputs where the permutation choice doesn't
%   matter.

    fprintf('=== v5_shuffle_null tests ===\n');
    results = {};
    results(end+1, :) = run_case('coupled_clears_threshold',         @t_coupled);
    results(end+1, :) = run_case('independent_does_not_clear',       @t_independent);
    results(end+1, :) = run_case('trial_perm_preserves_marginals',   @t_trial_perm_marginals);
    results(end+1, :) = run_case('same_seed_gives_same_null',        @t_deterministic);
    results(end+1, :) = run_case('different_seed_gives_different',   @t_different_seed);
    results(end+1, :) = run_case('threshold_formula',                @t_threshold);
    results(end+1, :) = run_case('shapes',                           @t_shapes);
    results(end+1, :) = run_case('validation_errors',                @t_validation);
    results(end+1, :) = run_case('cross_language_degenerate_case',   @t_cross_lang);

    n_pass = sum(cellfun(@(x) x, results(:, 2)));
    n_total = size(results, 1);
    fprintf('\n%d/%d passed\n', n_pass, n_total);
    if n_pass < n_total
        error('test_v5_shuffle_null:fail', 'Some tests failed; see output above.');
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
function t_coupled()
    rng_state = rng(0, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n_trials = 40; n_bins = 10;
    L = randn(n_trials, n_bins, 2);
    X = cat(3, L, 0.2 * randn(n_trials, n_bins, 2));
    Y = cat(3, L + 0.2 * randn(n_trials, n_bins, 2), 0.2 * randn(n_trials, n_bins, 2));
    X_flat = reshape(permute(X, [2 1 3]), n_trials*n_bins, 4);
    Y_flat = reshape(permute(Y, [2 1 3]), n_trials*n_bins, 4);
    % Use the spec's flatten convention: (trial, bin, feat) → (trial*bin, feat) with bin fast.
    % Above permute puts bins first → reshape walks bins within trial, then next trial — matches numpy C-order.
    real = v5_cca_fit_cv(X_flat, Y_flat, 10, true, 0);
    null = v5_shuffle_null(X_flat, Y_flat, 20, n_bins, 10, 0, true, 0);
    fprintf('    real cc_cv(1)=%.3f, null mean=%.3f, threshold=%.3f\n', ...
        real.cc_cv(1), null.mean, null.threshold);
    assert(real.cc_cv(1) > null.threshold, ...
        'real %g should clear threshold %g', real.cc_cv(1), null.threshold);
end

function t_independent()
    rng_state = rng(1, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n_trials = 40; n_bins = 10;
    X = randn(n_trials, n_bins, 3);
    Y = randn(n_trials, n_bins, 3);
    X_flat = reshape(permute(X, [2 1 3]), n_trials*n_bins, 3);
    Y_flat = reshape(permute(Y, [2 1 3]), n_trials*n_bins, 3);
    real = v5_cca_fit_cv(X_flat, Y_flat, 10, true, 0);
    null = v5_shuffle_null(X_flat, Y_flat, 20, n_bins, 10, 0, true, 0);
    fprintf('    real cc_cv(1)=%.3f, null mean=%.3f, threshold=%.3f\n', ...
        real.cc_cv(1), null.mean, null.threshold);
    assert(real.cc_cv(1) < null.threshold);
end

function t_trial_perm_marginals()
    rng_state = rng(2, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    % 6 trials × 5 bins flatten — call shuffle_null with n_shuffles=1 and
    % inspect via the public output indirectly: marginals are unchanged.
    Y = randn(30, 4);
    % Manually trial-permute via the same logic and check.
    perm = randperm(6);
    spt = 5;
    block_idx = bsxfun(@plus, (perm(:) - 1) * spt, 1:spt);
    block_idx = reshape(block_idx.', [], 1);
    Y_shuf = Y(block_idx, :);
    assert(size(Y_shuf, 1) == size(Y, 1));
    % Same set of rows: sorted rows are equal.
    assert(isequal(sortrows(Y_shuf), sortrows(Y)));
    % Each trial-block is preserved as a contiguous unit.
    for new_trial = 1:6
        new_block = Y_shuf((new_trial-1)*spt + 1 : new_trial*spt, :);
        original  = Y((perm(new_trial)-1)*spt + 1 : perm(new_trial)*spt, :);
        assert(isequal(new_block, original));
    end
end

function t_deterministic()
    rng_state = rng(3, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(50, 3);  Y = randn(50, 3);
    n1 = v5_shuffle_null(X, Y, 10, 5, 10, 42, true, 0);
    n2 = v5_shuffle_null(X, Y, 10, 5, 10, 42, true, 0);
    assert(isequal(n1.cc_cv, n2.cc_cv));
    assert(n1.threshold == n2.threshold);
end

function t_different_seed()
    rng_state = rng(4, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(50, 3);  Y = randn(50, 3);
    n1 = v5_shuffle_null(X, Y, 10, 5, 10, 1, true, 0);
    n2 = v5_shuffle_null(X, Y, 10, 5, 10, 2, true, 0);
    assert(~isequal(n1.cc_cv, n2.cc_cv));
end

function t_threshold()
    rng_state = rng(5, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(30, 2);  Y = randn(30, 2);
    n = v5_shuffle_null(X, Y, 20, 3, 10, 0, true, 0);
    expected = mean(n.cc1_null) + 3 * std(n.cc1_null);
    assert(abs(n.threshold - expected) < 1e-12);
end

function t_shapes()
    rng_state = rng(6, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    X = randn(40, 3);  Y = randn(40, 4);
    n = v5_shuffle_null(X, Y, 8, 4, 10, 0, true, 0);
    assert(isequal(size(n.cc_cv), [8, min(3, 4)]));
    assert(numel(n.cc1_null) == 8);
end

function t_validation()
    raised = false;
    try
        v5_shuffle_null(zeros(10, 2), zeros(10, 2), 5, 3);
    catch err
        raised = strcmp(err.identifier, 'v5_shuffle_null:badShape');
    end
    assert(raised, 'expected v5_shuffle_null:badShape');
end

function t_cross_lang()
    % CROSS-LANG (degenerate): with X = Y = column vector of 1..16,
    % 2 trials × 8 bins, ANY trial permutation gives cc=1 (since X and Y
    % are perfectly correlated). Tests the threshold formula and output
    % shape agree, not the RNG path (which differs between languages).
    X = (1:16).';
    Y = X;
    null = v5_shuffle_null(X, Y, 2, 8, 2, 0, false, 0);
    assert(isequal(size(null.cc_cv), [2, 1]));
    assert(max(abs(null.cc_cv - 1)) < 1e-10);
    assert(abs(null.mean - 1) < 1e-10);
    % With std=0 across identical 1s, threshold = mean + 3*0 = 1
    assert(abs(null.threshold - 1) < 1e-10);
end
