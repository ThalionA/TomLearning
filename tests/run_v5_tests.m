function results = run_v5_tests()
%RUN_V5_TESTS  Run every v5 test file and report a single pass/fail summary.
%
%   From the TomLearning repo root:
%       cd HC_V1_Code; addpath(pwd); cd ../tests; run_v5_tests
%
%   Or as a one-liner:
%       cd /Users/theoamvr/Desktop/Experiments/TomLearning && \
%         matlab -batch "addpath('HC_V1_Code'); cd tests; run_v5_tests"
%
%   This runner exists so the user can run every test in one call rather
%   than enumerating them. It reports failures concisely and exits with a
%   non-zero error if any test failed (so it composes with CI / shell
%   pipelines if the user later wires one up).

    test_files = {
        'test_v5_residualise',     ...
        'test_v5_pca_reduce',      ...
        'test_v5_cca_fit',         ...
        'test_v5_shuffle_null',    ...
        'test_v5_remaining',       ...
        'test_v5_cell_weights',    ...
        'test_v5_aux_analyses',    ...
    };

    results = struct('file', {}, 'ok', {}, 'err', {});
    fprintf('\n');
    for i = 1:numel(test_files)
        name = test_files{i};
        fprintf('### %s ###\n', name);
        try
            feval(name);
            results(end+1) = struct('file', name, 'ok', true, 'err', '');  %#ok<AGROW>
        catch err
            fprintf('  *** %s threw: %s\n', name, err.message);
            results(end+1) = struct('file', name, 'ok', false, 'err', err.message);  %#ok<AGROW>
        end
        fprintf('\n');
    end

    n_ok = sum([results.ok]);
    n_total = numel(results);
    fprintf('======================================================\n');
    fprintf('OVERALL: %d/%d test files passed.\n', n_ok, n_total);
    fprintf('======================================================\n');
    if n_ok < n_total
        for i = 1:n_total
            if ~results(i).ok
                fprintf('  FAIL %s: %s\n', results(i).file, results(i).err);
            end
        end
        error('run_v5_tests:fail', '%d/%d test files failed', ...
            n_total - n_ok, n_total);
    end
end
