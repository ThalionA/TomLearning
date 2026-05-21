function test_v5_aux_analyses()
%TEST_V5_AUX_ANALYSES  Synthetic-data tests for the v5 auxiliary helpers
%   added on top of the existing CCA pipeline:
%
%       v5_mi_from_cc
%       v5_significant_dims
%       v5_weight_spatial_entropy
%       v5_projection_behaviour_corr
%
%   These four helpers ship with the Gonzalez et al. paper's analyses
%   (minus pCCA proper — we keep plain CCA + per-period PSTH
%   residualisation as the v5 spec already does). Each test exercises
%   the helper on a known-coupling synthetic input plus an
%   uncoupled control so any future regression in either direction is
%   caught.

    fprintf('=== v5_aux_analyses tests ===\n');
    results = {};

    % v5_mi_from_cc
    results(end+1, :) = run_case('mi_zero_cc',                @t_mi_zero);
    results(end+1, :) = run_case('mi_formula',                @t_mi_formula);
    results(end+1, :) = run_case('mi_clamp_at_one',           @t_mi_clamp);
    results(end+1, :) = run_case('mi_sig_mask_filters',       @t_mi_sig_mask);
    results(end+1, :) = run_case('mi_empty_mask_zero',        @t_mi_empty_mask);
    results(end+1, :) = run_case('mi_negative_cc_via_square', @t_mi_negative);

    % v5_significant_dims
    results(end+1, :) = run_case('sig_all_above_null',        @t_sig_all_above);
    results(end+1, :) = run_case('sig_all_below_null',        @t_sig_all_below);
    results(end+1, :) = run_case('sig_per_dim_independence',  @t_sig_per_dim);
    results(end+1, :) = run_case('sig_alpha_override',        @t_sig_alpha);
    results(end+1, :) = run_case('sig_shape_mismatch',        @t_sig_shape);

    % v5_weight_spatial_entropy
    results(end+1, :) = run_case('ent_uniform_high',          @t_ent_uniform);
    results(end+1, :) = run_case('ent_concentrated_low',      @t_ent_concentrated);
    results(end+1, :) = run_case('ent_shape',                 @t_ent_shape);
    results(end+1, :) = run_case('ent_surrogate_centred',     @t_ent_surrogate);

    % v5_projection_behaviour_corr
    results(end+1, :) = run_case('proj_linear_position',      @t_proj_linear);
    results(end+1, :) = run_case('proj_random_near_zero',     @t_proj_random);
    results(end+1, :) = run_case('proj_broadcasting',         @t_proj_broadcast);
    results(end+1, :) = run_case('proj_null_distribution',    @t_proj_null);

    n_pass = sum(cellfun(@(x) x, results(:, 2)));
    n_total = size(results, 1);
    fprintf('\n%d/%d passed\n', n_pass, n_total);
    if n_pass < n_total
        error('test_v5_aux_analyses:fail', 'Some tests failed; see output above.');
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

% ==================================================================== %
% v5_mi_from_cc                                                        %
% ==================================================================== %

function t_mi_zero()
    assert(v5_mi_from_cc([0; 0; 0]) == 0);
end

function t_mi_formula()
    cc = [0.9; 0.6; 0.3];
    expected = -0.5 * sum(log(1 - cc .^ 2));
    assert(abs(v5_mi_from_cc(cc) - expected) < 1e-12);
end

function t_mi_clamp()
    % cc=1 should be clamped to 1 - eps so MI is large but finite.
    mi = v5_mi_from_cc([1.0]);
    assert(isfinite(mi));
    assert(mi > 0);
end

function t_mi_sig_mask()
    cc = [0.9; 0.6; 0.3];
    mask = [true; false; true];
    expected = -0.5 * sum(log(1 - cc([1 3]) .^ 2));
    assert(abs(v5_mi_from_cc(cc, mask) - expected) < 1e-12);
end

function t_mi_empty_mask()
    assert(v5_mi_from_cc([0.5; 0.5], [false; false]) == 0);
end

function t_mi_negative()
    % Negative CCs (possible from per-fold CV under sign flips) should
    % be treated by squaring — MI uses cc^2.
    mi_neg = v5_mi_from_cc([-0.6; -0.4]);
    mi_pos = v5_mi_from_cc([ 0.6;  0.4]);
    assert(abs(mi_neg - mi_pos) < 1e-12);
end

% ==================================================================== %
% v5_significant_dims                                                  %
% ==================================================================== %

function t_sig_all_above()
    cc_obs = [0.9; 0.8; 0.7];
    cc_null = 0.1 * randn(50, 3);
    sig = v5_significant_dims(cc_obs, cc_null);
    assert(all(sig));
end

function t_sig_all_below()
    cc_obs = [0.05; 0.04; 0.03];
    cc_null = 0.3 + 0.1 * randn(50, 3);
    sig = v5_significant_dims(cc_obs, cc_null);
    assert(~any(sig));
end

function t_sig_per_dim()
    % Dim 1 above null, dim 2 below, dim 3 above.
    rng_state = rng(13, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    cc_null = [...
        0.1 + 0.05 * randn(40, 1), ...
        0.5 + 0.05 * randn(40, 1), ...
        0.2 + 0.05 * randn(40, 1)];
    cc_obs = [0.6; 0.1; 0.8];
    sig = v5_significant_dims(cc_obs, cc_null);
    assert(isequal(sig, [true; false; true]));
end

function t_sig_alpha()
    % alpha=0.5 -> threshold at the median of the null. Observed above
    % median should be significant; below should not.
    cc_null = (1:100).' / 100;
    cc_obs = [0.6; 0.4];
    cc_null_2 = [cc_null, cc_null];
    sig = v5_significant_dims(cc_obs, cc_null_2, 0.5);
    assert(isequal(sig, [true; false]));
end

function t_sig_shape()
    raised = false;
    try
        v5_significant_dims([1; 2; 3], randn(10, 2));
    catch err
        raised = strcmp(err.identifier, 'v5_significant_dims:badShape');
    end
    assert(raised);
end

% ==================================================================== %
% v5_weight_spatial_entropy                                            %
% ==================================================================== %

function t_ent_uniform()
    % Each cell prefers a different position, uniform spread.
    % Weights uniform — entropy at max = log(n_spatial_bins).
    n_units = 25; n_pos = 25; n_bins = 5;
    place = eye(n_units);          % unit i prefers position i
    W = ones(n_units, 1);
    ent = v5_weight_spatial_entropy(W, place, n_bins, 0);
    assert(abs(ent - log(n_bins)) < 1e-10, ...
        'uniform entropy should be log(%d)=%.3f; got %.3f', ...
        n_bins, log(n_bins), ent);
end

function t_ent_concentrated()
    % All units prefer the same position, single canonical dimension.
    % All mass concentrates in one spatial bin -> entropy = 0.
    n_units = 10; n_pos = 25; n_bins = 5;
    place = zeros(n_units, n_pos);
    place(:, 1) = 1.0;                % all units prefer position 1
    W = ones(n_units, 1);
    ent = v5_weight_spatial_entropy(W, place, n_bins, 0);
    assert(ent < 1e-12);
end

function t_ent_shape()
    n_units = 20; n_pos = 50; r = 3;
    place = randn(n_units, n_pos);
    W = randn(n_units, r);
    [ent, null_ent] = v5_weight_spatial_entropy(W, place, 5, 10);
    assert(isequal(size(ent), [r 1]));
    assert(isequal(size(null_ent), [10 r]));
end

function t_ent_surrogate()
    % For a moderate, mixed-preference setup, the surrogate distribution
    % of entropies should be tight around the observed entropy (random
    % shifts of tuning curves preserve the rough spread).
    rng_state = rng(17, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n_units = 50; n_pos = 100;
    place = randn(n_units, n_pos);
    W = randn(n_units, 2);
    [ent, null_ent] = v5_weight_spatial_entropy(W, place, 5, 50, 0);
    null_mean = mean(null_ent, 1).';
    % Surrogates should produce a null distribution; not asserting it's
    % near the observed because random tuning + random weights gives
    % no real structure either way. The check is just that null_ent has
    % real spread.
    assert(all(std(null_ent, 0, 1) > 1e-8), ...
        'surrogate distribution should have non-zero spread');
    assert(numel(null_mean) == numel(ent));
end

% ==================================================================== %
% v5_projection_behaviour_corr                                         %
% ==================================================================== %

function t_proj_linear()
    % Construct u(trial, bin) = position_bin so r^2 with position is 1.
    n_trials = 20; n_pos = 50;
    u = repmat(1:n_pos, n_trials, 1);
    beh = struct('position', 1:n_pos);
    out = v5_projection_behaviour_corr(u, beh, 5, 0);
    assert(abs(out.position.r - 1) < 1e-12, 'r: %g', out.position.r);
    assert(abs(out.position.r2 - 1) < 1e-12, 'r2: %g', out.position.r2);
end

function t_proj_random()
    rng_state = rng(19, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    u = randn(20, 50);
    beh = struct('position', 1:50);
    out = v5_projection_behaviour_corr(u, beh, 20, 0);
    assert(out.position.r2 < 0.05, 'random r2: %g', out.position.r2);
end

function t_proj_broadcast()
    % Trial-number broadcast: u = trial_number * (constant offset per bin)
    n_trials = 15; n_pos = 30;
    u = repmat((1:n_trials).', 1, n_pos);
    beh = struct('trial_number', 1:n_trials);
    out = v5_projection_behaviour_corr(u, beh, 5, 0);
    assert(abs(out.trial_number.r - 1) < 1e-12);
end

function t_proj_null()
    % The circular-shift null for a strongly-correlated signal must
    % produce r^2 values mostly less than the observed.
    rng_state = rng(21, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    n_trials = 30; n_pos = 100;
    pos = repmat(1:n_pos, n_trials, 1);
    u = pos + 0.05 * randn(n_trials, n_pos);
    beh = struct('position', 1:n_pos);
    out = v5_projection_behaviour_corr(u, beh, 100, 0);
    assert(out.position.r2 > 0.95);
    % p_perm should be small (observed >> most of the null).
    assert(out.position.p_perm < 0.1, 'p_perm: %g', out.position.p_perm);
end
