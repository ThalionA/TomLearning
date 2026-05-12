function test_v5_residualise()
%TEST_V5_RESIDUALISE  Synthetic-data tests for v5_residualise.
%
%   Mirrors the Python tests in
%       /Users/theoamvr/Desktop/Experiments/IBL/tests/test_cca_core.py::TestResidualise
%
%   Includes the cross-language reference case (test 8 / `cross_language_reference_case`).
%   The reference output was computed once in Python with seed-free,
%   np.arange-based input so MATLAB and Python produce the same number.
%
%   Run from the TomLearning repo root:
%       cd HC_V1_Code; addpath(pwd); cd ..
%       cd tests; test_v5_residualise

    fprintf('=== v5_residualise tests ===\n');
    results = {};
    results(end+1, :) = run_case('planted_psth_recovered_exactly',     @t_planted_psth_recovered_exactly);
    results(end+1, :) = run_case('zero_psth_is_noop',                  @t_zero_psth_is_noop);
    results(end+1, :) = run_case('period_mask_restricts_residualisation', @t_period_mask);
    results(end+1, :) = run_case('multiple_groups_each_zeroed',        @t_multiple_groups);
    results(end+1, :) = run_case('singleton_group_is_zeroed',          @t_singleton_group);
    results(end+1, :) = run_case('shape_mismatch_raises',              @t_shape_mismatch);
    results(end+1, :) = run_case('input_unchanged',                    @t_input_unchanged);
    results(end+1, :) = run_case('cross_language_reference_case',      @t_cross_lang_reference);
    results(end+1, :) = run_case('nan_labels_handled',                 @t_nan_labels);

    n_pass = sum(cellfun(@(x) x, results(:, 2)));
    n_total = size(results, 1);
    fprintf('\n%d/%d passed\n', n_pass, n_total);
    if n_pass < n_total
        error('test_v5_residualise:fail', 'Some tests failed; see output above.');
    end
end

% -------------------------------------------------------------------- %
% Test harness                                                         %
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
% Individual test cases                                                %
% -------------------------------------------------------------------- %

function t_planted_psth_recovered_exactly()
    % Bake in a known per-condition PSTH; residualisation removes it exactly.
    rng_state = rng(0, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n_trials = 60; n_bins = 20; n_units = 8;
    labels = [zeros(30, 1); ones(30, 1)];
    psth = zeros(2, n_bins, n_units);
    psth(1, 6, 1) = 1.0;        % group 0, bin 6 (Python bin 5), unit 1
    psth(2, 11, 4) = -2.0;      % group 1, bin 11 (Python bin 10), unit 4
    noise = randn(n_trials, n_bins, n_units);
    S = noise + psth(labels + 1, :, :);
    S_res = v5_residualise(S, labels);
    for g = 0:1
        m = labels == g;
        gm = mean(S_res(m, :, :), 1);
        assert(max(abs(gm(:))) < 1e-12, ...
            'group %d residual mean: %g', g, max(abs(gm(:))));
    end
    % Manual check: residual = S − group_mean(S) per group.
    expected = S;
    for g = 0:1
        m = labels == g;
        expected(m, :, :) = S(m, :, :) - mean(S(m, :, :), 1);
    end
    assert(max(abs(S_res(:) - expected(:))) < 1e-12);
end

function t_zero_psth_is_noop()
    % If the input is already mean-zero per condition, residualise is a no-op.
    rng_state = rng(1, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    raw = randn(40, 15, 5);
    S = raw - mean(raw, 1);
    S_res = v5_residualise(S, zeros(40, 1));
    assert(max(abs(S_res(:) - S(:))) < 1e-12);
    assert(max(abs(reshape(mean(S_res, 1), [], 1))) < 1e-12);
end

function t_period_mask()
    % Inside the period mask: per-group bin means are zero. Outside: untouched.
    rng_state = rng(2, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    labels = [zeros(15, 1); ones(15, 1)];
    psth = zeros(2, 10, 4);
    psth(1, :, :) = 3.0;
    psth(2, :, :) = -1.0;
    S = randn(30, 10, 4) + psth(labels + 1, :, :);
    mask = false(10, 1);  mask(4:7) = true;     % MATLAB 1-indexed: bins 4..7 = Python 3..6 (inclusive Python 3,4,5,6)
    S_res = v5_residualise(S, labels, mask);
    for g = 0:1
        m = labels == g;
        inside = mean(S_res(m, mask, :), 1);
        assert(max(abs(inside(:))) < 1e-12, 'inside mask, group %d', g);
    end
    out_S = S(:, ~mask, :);
    out_R = S_res(:, ~mask, :);
    assert(isequal(out_S, out_R), 'outside-mask bins changed');
end

function t_multiple_groups()
    rng_state = rng(3, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    labels = randi([0, 4], 100, 1);
    S = randn(100, 8, 6);
    S_res = v5_residualise(S, labels);
    groups = unique(labels);
    for k = 1:numel(groups)
        m = labels == groups(k);
        gm = mean(S_res(m, :, :), 1);
        assert(max(abs(gm(:))) < 1e-12, 'group %d', groups(k));
    end
end

function t_singleton_group()
    rng_state = rng(4, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    S = randn(5, 4, 3);
    labels = [0; 0; 0; 0; 1];   % trial 5 alone in its group
    S_res = v5_residualise(S, labels);
    assert(max(abs(reshape(S_res(5, :, :), [], 1))) < 1e-12);
end

function t_shape_mismatch()
    raised = false;
    try
        v5_residualise(zeros(10, 5, 3), zeros(9, 1));
    catch err
        raised = strcmp(err.identifier, 'v5_residualise:badShape');
    end
    assert(raised, 'expected v5_residualise:badShape error');
end

function t_input_unchanged()
    % Workspace semantics in MATLAB are copy-on-write, but explicitly check.
    rng_state = rng(5, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    S = randn(20, 6, 4);
    S_orig = S;
    [~] = v5_residualise(S, zeros(20, 1));
    assert(isequal(S, S_orig));
end

function t_cross_lang_reference()
    % CROSS-LANG: deterministic input, output must match Python's
    % core.residualise(np.arange(24).reshape(4,3,2), [0,0,1,1])
    % to within 1e-12.
    flat_in_C = 0:23;
    % numpy C-order (last axis varies fastest, shape 4×3×2) becomes
    % MATLAB (reshape into [2 3 4], then permute to [4 3 2]).
    S = permute(reshape(flat_in_C, [2 3 4]), [3 2 1]);
    labels = [0; 0; 1; 1];
    S_res = v5_residualise(S, labels);
    expected_flat_C = [...
        -3, -3, -3, -3, -3, -3, ...
         3,  3,  3,  3,  3,  3, ...
        -3, -3, -3, -3, -3, -3, ...
         3,  3,  3,  3,  3,  3];
    expected = permute(reshape(expected_flat_C, [2 3 4]), [3 2 1]);
    assert(max(abs(S_res(:) - expected(:))) < 1e-12, ...
        'cross-lang mismatch; max diff = %g', max(abs(S_res(:) - expected(:))));
end

function t_nan_labels()
    rng_state = rng(6, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    labels = [0; 0; NaN; NaN; 1; 1];
    S = randn(6, 4, 3);
    S_res = v5_residualise(S, labels);
    nan_mask = isnan(labels);
    gm = mean(S_res(nan_mask, :, :), 1);
    assert(max(abs(gm(:))) < 1e-12);
end
