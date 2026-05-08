%% test_temporal_v4.m
% Synthetic-ground-truth tests for the v4 temporal CCA helpers.
%
% Run from MATLAB:
%   addpath('../HC_V1_Code'); test_temporal_v4
% Or:
%   run('/Users/theoamvr/Desktop/Experiments/TomLearning/tests/test_temporal_v4.m')
%
% The Python prototype tests/test_temporal_v4_proto.py exercises the same
% scenarios. This MATLAB script is the authoritative test for the actual
% code that ships in HC_V1_Code/.
%
% Tests:
%   T1 - Lagged AR(1) coupling (lag=+1) -> per-trial CC > 0.5,
%        per-trial IFI < 0 (u leads v -> negative).
%   T2 - Uncoupled regions -> real CC ~ shuffle CC.
%   T3 - Two valid blocks of unequal length -> weighted IFI dominated by
%        the long lagged block; min_block_bins drops the short one.
%   T4 - RZ-pre alignment with a coupling burst near RZ entry -> real-
%        minus-shuffle CC is higher in late aligned bins than early.
%        Pre-RZ per-trial IFI is negative on average.

clear; clc;

% Path to the helpers (assumes this file lives in tests/).
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'HC_V1_Code'));

rng(0);

n_pass = 0; n_fail = 0;

%% T1 - lagged AR(1) coupling
fprintf('\n[T1] Lagged AR(1) coupling, lag=+1 (u leads v)\n');
k1 = 3; k2 = 3; n_bins = 200;
n_trials = 12;
ifis = nan(n_trials, 1);
rs   = nan(n_trials, 1);
for it = 1:n_trials
    [X, Y] = make_lagged_ar1(n_bins, k1, k2, 1, 0.9, 0.6);
    [r, A, B, info] = v4_per_trial_cca(X, Y);
    if ~info.ok, continue; end
    rs(it) = r;
    is_valid = true(n_bins, 1);
    [ifi, ~] = v4_per_trial_ifi(X, Y, is_valid, A, B, -3:3, 12);
    ifis(it) = ifi;
end
mean_r = mean(rs, 'omitnan');
mean_ifi = mean(ifis, 'omitnan');
fprintf('    mean per-trial CC  = %.3f\n', mean_r);
fprintf('    mean per-trial IFI = %.3f (expect < 0)\n', mean_ifi);
[n_pass, n_fail] = check(mean_r > 0.5, ...
    'per-trial CC > 0.5 with strong coupling', n_pass, n_fail);
[n_pass, n_fail] = check(mean_ifi < -0.05, ...
    'per-trial IFI < -0.05 (u leads v)', n_pass, n_fail);

%% T2 - uncoupled
fprintf('\n[T2] Uncoupled regions: real CC vs shuffle CC\n');
n_trials = 20;
rs_real = nan(n_trials, 1);
rs_sh   = nan(n_trials, 1);
for it = 1:n_trials
    [X, Y] = make_uncoupled(n_bins, k1, k2, 0.6);
    [r, ~, ~, info] = v4_per_trial_cca(X, Y);
    if info.ok, rs_real(it) = r; end
    perm = randperm(n_bins);
    [rs, ~, ~, info_s] = v4_per_trial_cca(X, Y(perm, :));
    if info_s.ok, rs_sh(it) = rs; end
end
mean_real = mean(rs_real, 'omitnan');
mean_sh   = mean(rs_sh, 'omitnan');
fprintf('    real CC = %.3f, shuffle CC = %.3f\n', mean_real, mean_sh);
[n_pass, n_fail] = check(abs(mean_real - mean_sh) < 0.10, ...
    'real CC ~ shuffle CC for uncoupled regions (within 0.10)', n_pass, n_fail);

%% T3 - block weighted IFI
fprintf('\n[T3] Block weighted IFI: short noisy + long lagged block\n');
[Xl, Yl] = make_lagged_ar1(180, k1, k2, 1, 0.9, 0.6);
[Xs, Ys] = make_uncoupled(15, k1, k2, 0.6);
gap = 10;
X_full = [Xs; zeros(gap, k1); Xl];
Y_full = [Ys; zeros(gap, k2); Yl];
is_valid = [true(15, 1); false(gap, 1); true(180, 1)];

X_valid = X_full(is_valid, :);
Y_valid = Y_full(is_valid, :);
[r, A, B, info_cc] = v4_per_trial_cca(X_valid, Y_valid);
[n_pass, n_fail] = check(info_cc.ok, ...
    'per_trial_cca handles concatenated valid bins', n_pass, n_fail);

[ifi, info] = v4_per_trial_ifi(X_full, Y_full, is_valid, A, B, -3:3, 12);
fprintf('    blocks lengths: %s, weighted IFI = %.3f\n', ...
    mat2str(info.block_lens), ifi);
[n_pass, n_fail] = check(size(info.blocks, 1) == 2, ...
    'two contiguous blocks identified', n_pass, n_fail);
[n_pass, n_fail] = check(ifi < 0, ...
    'weighted IFI dominated by long lagged block (IFI < 0)', n_pass, n_fail);

% Drop short block via min_block_bins (set to 16 -> only 180-bin block survives)
[ifi2, info2] = v4_per_trial_ifi(X_full, Y_full, is_valid, A, B, -3:3, 16);
fprintf('    min_block=16 -> blocks=%s, IFI=%.3f\n', ...
    mat2str(info2.block_lens), ifi2);
[n_pass, n_fail] = check(size(info2.blocks, 1) == 1 && info2.block_lens(1) == 180, ...
    'min_block_bins filter drops short block', n_pass, n_fail);

%% T4 - RZ alignment
fprintf('\n[T4] RZ-pre alignment with coupling burst near RZ entry\n');
n_bins_rz = 400;
rz_entry_cm = 200;
pos = linspace(0, 250, n_bins_rz)';
i_rz = find(pos >= rz_entry_cm, 1, 'first');
fprintf('    rz_entry bin = %d\n', i_rz);
n_rz = 20;
n_trials_rz = 60;

X_aligned = nan(n_rz, k1, n_trials_rz);
Y_aligned = nan(n_rz, k2, n_trials_rz);
pre_ifis  = nan(n_trials_rz, 1);
n_valid = 0;
for it = 1:n_trials_rz
    [X, Y] = make_uncoupled(n_bins_rz, k1, k2, 0.6);
    burst_start = i_rz - 10;
    burst_end   = i_rz - 1;
    Xb = X(burst_start:burst_end, :);
    W = randn(k1, k2) * (0.95 / sqrt(k1));
    for j = 1:(burst_end - burst_start)
        t = burst_start + j;
        Y(t, :) = Xb(j, :) * W + 0.3 * randn(1, k2);
    end
    [r, A, B, info] = v4_per_trial_cca(X, Y);
    if ~info.ok, continue; end
    is_valid = true(n_bins_rz, 1);
    [Xs, Ys, ok_rz] = v4_rz_align_pre(X, Y, pos, is_valid, rz_entry_cm, n_rz);
    if ~ok_rz, continue; end
    n_valid = n_valid + 1;
    X_aligned(:, :, n_valid) = Xs;
    Y_aligned(:, :, n_valid) = Ys;
    Xc = Xs - repmat(mean(Xs, 1), n_rz, 1);
    Yc = Ys - repmat(mean(Ys, 1), n_rz, 1);
    u = Xc * A; v = Yc * B;
    pre_ifis(n_valid) = v4_ifi_from_lags(v4_lag_corr(u, v, -3:3), -3:3);
end
X_aligned = X_aligned(:, :, 1:n_valid);
Y_aligned = Y_aligned(:, :, 1:n_valid);
fprintf('    valid trials = %d\n', n_valid);
[n_pass, n_fail] = check(n_valid >= 40, ...
    'most trials yield a pre-RZ segment', n_pass, n_fail);

% Per-aligned-bin CC across trials, real vs shuffle
n_shuf = 20;
cc_real = nan(1, n_rz);
cc_sh   = nan(1, n_rz);
for b = 1:n_rz
    Xb = squeeze(X_aligned(b, :, :))';
    Yb = squeeze(Y_aligned(b, :, :))';
    [r, ~, ~, info] = v4_per_trial_cca(Xb, Yb);
    if info.ok, cc_real(b) = r; end
    rs_iter = nan(n_shuf, 1);
    for ish = 1:n_shuf
        perm = randperm(n_valid);
        [rs, ~, ~, info_s] = v4_per_trial_cca(Xb, Yb(perm, :));
        if info_s.ok, rs_iter(ish) = rs; end
    end
    cc_sh(b) = mean(rs_iter, 'omitnan');
end
cc_excess = cc_real - cc_sh;
fprintf('    real - shuffle CC, last 5 bins: %s\n', mat2str(cc_excess(end-4:end), 3));
fprintf('    real - shuffle CC, first 5 bins: %s\n', mat2str(cc_excess(1:5), 3));
[n_pass, n_fail] = check(mean(cc_excess(end-4:end), 'omitnan') > mean(cc_excess(1:5), 'omitnan') + 0.05, ...
    'real-minus-shuffle CC is higher near RZ entry', n_pass, n_fail);
mean_pre_ifi = mean(pre_ifis, 'omitnan');
fprintf('    mean pre-RZ IFI = %.3f (expect < 0)\n', mean_pre_ifi);
[n_pass, n_fail] = check(mean_pre_ifi < 0.0, ...
    'pre-RZ per-trial IFI is negative on average', n_pass, n_fail);

%% Summary
fprintf('\n=== %d passed, %d failed ===\n', n_pass, n_fail);
if n_fail > 0
    error('test_temporal_v4: %d test(s) failed', n_fail);
end


%% --- Helpers used only inside this test ---

function [n_pass, n_fail] = check(cond, msg, n_pass, n_fail)
    if cond
        fprintf('  ok: %s\n', msg);
        n_pass = n_pass + 1;
    else
        fprintf('  FAIL: %s\n', msg);
        n_fail = n_fail + 1;
    end
end

function [X, Y] = make_lagged_ar1(n_bins, k1, k2, lag, coupling, ar)
% Two regions where Y(t) is driven by X(t-lag) plus AR(1) noise.
    nx = randn(n_bins, k1);
    X = zeros(n_bins, k1);
    X(1, :) = nx(1, :);
    for t = 2:n_bins
        X(t, :) = ar * X(t-1, :) + nx(t, :);
    end
    ny = randn(n_bins, k2);
    Y = zeros(n_bins, k2);
    W = randn(k1, k2) * (coupling / sqrt(k1));
    for t = 1:n_bins
        src_t = t - lag;
        if src_t >= 1 && src_t <= n_bins
            Y(t, :) = X(src_t, :) * W + 0.5 * ny(t, :);
        else
            Y(t, :) = 0.5 * ny(t, :);
        end
        if t > 1
            Y(t, :) = 0.5 * Y(t, :) + 0.5 * (ar * Y(t-1, :));
        end
    end
end

function [X, Y] = make_uncoupled(n_bins, k1, k2, ar)
    X = ar1(n_bins, k1, ar);
    Y = ar1(n_bins, k2, ar);
end

function x = ar1(n_bins, k, ar)
    n = randn(n_bins, k);
    x = zeros(n_bins, k);
    x(1, :) = n(1, :);
    for t = 2:n_bins
        x(t, :) = ar * x(t-1, :) + n(t, :);
    end
end
