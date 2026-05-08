%% cca_simulation.m
% Sanity check and parameter sweep for CCA and Kernel CCA using synthetic data.
% This script simulates two views (X and Y) with shared latent dynamics,
% and tests how temporal CCA degrades under various conditions:
%   1. Noise
%   2. Nonlinearities (Linear CCA vs Kernel CCA)
%   3. Bin width (Timescale smoothing)
%   4. Lag (Temporal offsets between views)
%   5. Dimensionality (True vs Over/Under-fitted CCA)
%   6. Sample Size (Number of bins per window)

clear; clc; close all;

% Reproducibility (CLAUDE.md scientific correctness rule).
rng(42);

%% Setup Plotting
% Resolve project root from this script's location so the script runs
% from any cwd (CLAUDE.md no-hardcoded-paths rule).
this_file = mfilename('fullpath');
[script_dir, ~] = fileparts(this_file);
figures_dir = fullfile(script_dir, '..', 'HC_V1_figures', 'CCA_Simulations');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

% Set default plotting parameters
set(groot, 'DefaultAxesFontSize', 12);
set(groot, 'DefaultLineLineWidth', 2);

% A helper `top_cc(X, Y, kappa)` is defined at the bottom of this file
% and returns [r_lin, r_kcca] -- the top-1 canonical correlation under
% linear CCA and RBF kernel CCA respectively. Sigma uses the median
% heuristic on the combined (X, Y) pairwise distances; kappa is the
% ridge (1e-3 is a safe default).

%% 0. KCCA / linear-CCA self-validation (S8)
% Sanity check the helper against three constructed cases:
%   (1) X and Y identical            -> CC1 = 1
%   (2) X random, Y = sin(X) + noise -> kernel >> linear
%   (3) X random, Y independent      -> CC1 close to bias floor
% Errors out if any of these violate expectations.
fprintf('--- Self-check: top_cc helper ---\n');
[ok_msg, ok] = check_top_cc();
fprintf('  %s\n', ok_msg);
if ~ok
    error('Self-check on top_cc failed; aborting before sweeps.');
end

%% 1. Effect of Noise (linear CCA + kernel CCA)
fprintf('--- Running Noise Sweep ---\n');
noise_levels = 0:0.5:5;
N = 500; % samples
D_true = 3;
cc_noise   = zeros(length(noise_levels), 1);   % linear CC1
cc_noise_k = zeros(length(noise_levels), 1);   % kernel CC1
cc_noise_lin_sh = zeros(length(noise_levels), 1);   % linear shuffle
cc_noise_kcca_sh = zeros(length(noise_levels), 1);  % kernel shuffle

% True latent
Z = randn(N, D_true);
W_x = randn(D_true, 10);
W_y = randn(D_true, 15);

for i = 1:length(noise_levels)
    noise = noise_levels(i);
    X = Z * W_x + noise * randn(N, 10);
    Y = Z * W_y + noise * randn(N, 15);

    [r_lin, r_kcca] = top_cc(X, Y, 1e-3);
    cc_noise(i)   = r_lin;
    cc_noise_k(i) = r_kcca;

    % Shuffle null: independently permute Y rows.
    perm = randperm(N);
    [r_lin_sh, r_kcca_sh] = top_cc(X, Y(perm,:), 1e-3);
    cc_noise_lin_sh(i)  = r_lin_sh;
    cc_noise_kcca_sh(i) = r_kcca_sh;
end

figure('Name', 'Effect of Noise', 'Position', [100 100 600 400]);
plot(noise_levels, cc_noise,           '-o', 'Color', [0.2 0.6 0.8], 'DisplayName','Linear CCA real'); hold on;
plot(noise_levels, cc_noise_lin_sh,    '--', 'Color', [0.2 0.6 0.8], 'DisplayName','Linear CCA shuffle');
plot(noise_levels, cc_noise_k,         '-s', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel CCA real');
plot(noise_levels, cc_noise_kcca_sh,   '--', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel CCA shuffle');
xlabel('Noise Standard Deviation');
ylabel('CC1');
title(sprintf('Sweep 1: Noise (N=%d, D_{true}=%d)', N, D_true));
legend('Location','best');
ylim([0 1]); grid on;
saveas(gcf, fullfile(figures_dir, 'sim_1_noise.png'));

%% 2. Nonlinearities — linear vs kernel CCA across families (S7)
% Identity (control), tanh (mild), squared (severe), sine+cosine (oscillatory).
% View 1 is always a linear projection of the 1D latent. View 2 applies
% the chosen nonlinearity to z before mixing. We compare linear and
% kernel CCA, with a shuffle null overlaid as a hatched bar.
fprintf('--- Running Nonlinearity sweep ---\n');
N_nonlin = 500;
Z_nonlin = linspace(-3, 3, N_nonlin)' + 0.5*randn(N_nonlin, 1);
X_nonlin = Z_nonlin * randn(1, 5) + 0.2*randn(N_nonlin, 5);

nonlin_funcs = {
    'identity',   @(z) z;
    'tanh',       @(z) tanh(z);
    'squared',    @(z) z.^2;
    'sin+cos',    @(z) [sin(z), cos(z)] * randn(2,1);  % 1D combined
};

n_nl = size(nonlin_funcs, 1);
res_lin    = nan(n_nl, 1);
res_lin_sh = nan(n_nl, 1);
res_k      = nan(n_nl, 1);
res_k_sh   = nan(n_nl, 1);

for inl = 1:n_nl
    name = nonlin_funcs{inl, 1};
    f    = nonlin_funcs{inl, 2};
    Y_nl = f(Z_nonlin) * randn(1, 5) + 0.2*randn(N_nonlin, 5);

    [r_lin, r_kcca] = top_cc(X_nonlin, Y_nl, 1e-3);
    res_lin(inl) = r_lin;
    res_k(inl)   = r_kcca;

    perm = randperm(N_nonlin);
    [r_lin_sh, r_kcca_sh] = top_cc(X_nonlin, Y_nl(perm,:), 1e-3);
    res_lin_sh(inl) = r_lin_sh;
    res_k_sh(inl)   = r_kcca_sh;
    fprintf('  %-10s: lin=%.3f (sh %.3f)  kernel=%.3f (sh %.3f)\n', ...
        name, res_lin(inl), res_lin_sh(inl), res_k(inl), res_k_sh(inl));
end

figure('Name', 'Nonlinearity sweep', 'Position', [150 150 750 420]);
xs = 1:n_nl;
group = [res_lin, res_lin_sh, res_k, res_k_sh];
b = bar(xs, group);
b(1).FaceColor = [0.2 0.6 0.8]; b(1).DisplayName = 'Linear real';
b(2).FaceColor = [0.4 0.7 0.9]; b(2).DisplayName = 'Linear shuffle';
b(3).FaceColor = [0.8 0.4 0.2]; b(3).DisplayName = 'Kernel real';
b(4).FaceColor = [0.9 0.6 0.4]; b(4).DisplayName = 'Kernel shuffle';
set(gca, 'XTickLabel', nonlin_funcs(:,1));
ylabel('CC1'); ylim([0 1]);
title('Sweep 2: nonlinearity family (linear vs kernel CCA, real vs shuffle)');
legend('Location','best'); grid on;
saveas(gcf, fullfile(figures_dir, 'sim_2_nonlinearity.png'));

%% 3. Timescale (Bin Width Smoothing) — linear + kernel
fprintf('--- Running Timescale Sweep ---\n');
% Simulate fine timescale (e.g. 1ms bins) and then smooth by different window sizes
N_fine = 5000;
Z_fine = filter(ones(50,1)/50, 1, randn(N_fine, D_true)); % auto-correlated latent

X_fine = Z_fine * W_x + 2*randn(N_fine, 10);
Y_fine = Z_fine * W_y + 2*randn(N_fine, 15);

smooth_windows = [1, 5, 10, 25, 50, 100];
% S5: hold sample count fixed across smoothing levels by always taking
% n_fixed evenly-spaced samples *after* smoothing. Previously, larger w
% both smoothed and reduced the sample count, conflating SNR gain with
% loss of statistical power.
n_fixed = 500;
cc_smooth     = nan(length(smooth_windows), 1);
cc_smooth_k   = nan(length(smooth_windows), 1);
cc_smooth_sh  = nan(length(smooth_windows), 1);
cc_smooth_ksh = nan(length(smooth_windows), 1);

for i = 1:length(smooth_windows)
    w = smooth_windows(i);
    X_sm = movmean(X_fine, w, 1);
    Y_sm = movmean(Y_fine, w, 1);
    % Sample n_fixed evenly-spaced bins, regardless of w.
    idx = round(linspace(1, N_fine, n_fixed));
    Xs = X_sm(idx,:); Ys = Y_sm(idx,:);
    if size(Xs,1) <= size(Xs,2) + size(Ys,2) + 5, continue; end

    [r_lin, r_kcca] = top_cc(Xs, Ys, 1e-3);
    cc_smooth(i)   = r_lin;
    cc_smooth_k(i) = r_kcca;

    perm = randperm(size(Ys,1));
    [r_lin_sh, r_kcca_sh] = top_cc(Xs, Ys(perm,:), 1e-3);
    cc_smooth_sh(i)  = r_lin_sh;
    cc_smooth_ksh(i) = r_kcca_sh;
end

figure('Name', 'Effect of Timescale', 'Position', [200 200 600 400]);
plot(smooth_windows, cc_smooth,     '-o', 'Color', [0.3 0.7 0.3], 'DisplayName','Linear real'); hold on;
plot(smooth_windows, cc_smooth_sh,  '--', 'Color', [0.3 0.7 0.3], 'DisplayName','Linear shuffle');
plot(smooth_windows, cc_smooth_k,   '-s', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel real');
plot(smooth_windows, cc_smooth_ksh, '--', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel shuffle');
xlabel('Smoothing window / bin size (samples)');
ylabel('CC1');
title(sprintf('Sweep 3: timescale smoothing (sample count fixed at n=%d)', n_fixed));
legend('Location','best'); ylim([0 1]); grid on;
saveas(gcf, fullfile(figures_dir, 'sim_3_timescale.png'));

%% 4. Lag effect (linear + kernel)
fprintf('--- Running Lag Sweep ---\n');
% Shift Y relative to X by a known true lag, in *native* fine-time samples.
true_lag_native = 25;
Y_lagged = circshift(Y_fine, true_lag_native);

% Subsample by factor `subs`. Sweep the lag axis in native samples (then
% convert to subsampled units for the actual fit).
subs = 10;
X_sub = X_fine(1:subs:end, :);
Y_sub = Y_lagged(1:subs:end, :);

lags_native = -50:5:50;             % imposed lag in native samples
cc_lag      = nan(length(lags_native), 1);
cc_lag_k    = nan(length(lags_native), 1);

for i = 1:length(lags_native)
    L_native = lags_native(i);
    L = round(L_native / subs);     % shift in subsampled units
    if L >= 0
        X_w = X_sub(1:end-L, :);
        Y_w = Y_sub(1+L:end, :);
    else
        X_w = X_sub(1-L:end, :);
        Y_w = Y_sub(1:end+L, :);
    end
    if size(X_w, 1) <= size(X_sub, 2) + size(Y_sub, 2) + 5, continue; end
    [r_lin, r_kcca] = top_cc(X_w, Y_w, 1e-3);
    cc_lag(i)   = r_lin;
    cc_lag_k(i) = r_kcca;
end

figure('Name', 'Effect of Lag', 'Position', [250 250 600 400]);
plot(lags_native, cc_lag,   '-o', 'Color', [0.2 0.6 0.8], 'DisplayName','Linear CCA'); hold on;
plot(lags_native, cc_lag_k, '-s', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel CCA');
xline(true_lag_native, 'k--', 'True lag', 'HandleVisibility','off');
xlabel('Imposed lag (native samples; subsampling=10)');
ylabel('CC1');
title('Sweep 4: recovering a known lag (axis in native units)');
legend('Location','best'); ylim([0 1]); grid on;
saveas(gcf, fullfile(figures_dir, 'sim_4_lag.png'));

%% 5. Dimensionality (linear + kernel; with shuffle null)
fprintf('--- Running Dimensionality Sweep ---\n');
dims_to_test = 1:10;
cc_dims      = nan(length(dims_to_test), 1);   % linear real CC1
cc_dims_k    = nan(length(dims_to_test), 1);   % kernel real CC1
cc_dims_sh   = nan(length(dims_to_test), 1);   % linear shuffle
cc_dims_ksh  = nan(length(dims_to_test), 1);   % kernel shuffle

[~, scoreX] = pca(X_fine);
[~, scoreY] = pca(Y_fine);

for i = 1:length(dims_to_test)
    d = dims_to_test(i);
    k_x = min(d, size(scoreX,2));
    k_y = min(d, size(scoreY,2));
    Xs = scoreX(:, 1:k_x);
    Ys = scoreY(:, 1:k_y);
    [r_lin, r_kcca] = top_cc(Xs, Ys, 1e-3);
    cc_dims(i)   = r_lin;
    cc_dims_k(i) = r_kcca;
    perm = randperm(size(Ys,1));
    [r_lin_sh, r_kcca_sh] = top_cc(Xs, Ys(perm,:), 1e-3);
    cc_dims_sh(i)  = r_lin_sh;
    cc_dims_ksh(i) = r_kcca_sh;
end

figure('Name', 'Effect of Dimensionality', 'Position', [300 300 600 400]);
plot(dims_to_test, cc_dims,    '-o', 'Color', [0.2 0.6 0.8], 'DisplayName','Linear real'); hold on;
plot(dims_to_test, cc_dims_sh, '--', 'Color', [0.2 0.6 0.8], 'DisplayName','Linear shuffle');
plot(dims_to_test, cc_dims_k,  '-s', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel real');
plot(dims_to_test, cc_dims_ksh,'--', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel shuffle');
xline(D_true, 'k--', 'True latent dim', 'HandleVisibility','off');
xlabel('Number of PCA components retained before CCA');
ylabel('CC1');
title('Sweep 5: dimensionality (real vs shuffle)');
legend('Location','best'); ylim([0 1]); grid on;
saveas(gcf, fullfile(figures_dir, 'sim_5_dimensionality.png'));

%% 6. Sample Size — bias floor on pure noise (linear + kernel)
fprintf('--- Running Sample Size Sweep (pure noise) ---\n');
sample_sizes = [15, 20, 30, 50, 100, 200, 500];
num_iterations = 20;

k_X = 10; k_Y = 15;

cc_samples_mean   = nan(length(sample_sizes), 1);
cc_samples_std    = nan(length(sample_sizes), 1);
cc_samples_k_mean = nan(length(sample_sizes), 1);
cc_samples_k_std  = nan(length(sample_sizes), 1);

for i = 1:length(sample_sizes)
    n_samp = sample_sizes(i);
    if n_samp <= k_X + k_Y, continue; end
    r_iters_lin = nan(num_iterations, 1);
    r_iters_k   = nan(num_iterations, 1);
    for it = 1:num_iterations
        X_noise = randn(n_samp, k_X);
        Y_noise = randn(n_samp, k_Y);
        [r_lin, r_kcca] = top_cc(X_noise, Y_noise, 1e-3);
        r_iters_lin(it) = r_lin;
        r_iters_k(it)   = r_kcca;
    end
    cc_samples_mean(i)   = mean(r_iters_lin, 'omitnan');
    cc_samples_std(i)    = std(r_iters_lin, 'omitnan');
    cc_samples_k_mean(i) = mean(r_iters_k,   'omitnan');
    cc_samples_k_std(i)  = std(r_iters_k,    'omitnan');
end

figure('Name', 'Effect of Sample Size — pure noise', 'Position', [350 350 600 400]);
errorbar(sample_sizes, cc_samples_mean,   cc_samples_std,   '-o', 'Color', [0.2 0.6 0.8], 'DisplayName','Linear (noise)'); hold on;
errorbar(sample_sizes, cc_samples_k_mean, cc_samples_k_std, '-s', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel (noise)');
xline(k_X + k_Y, 'r--', 'k1+k2', 'HandleVisibility','off');
xlabel('Sample size (bins per window)');
ylabel('Spurious CC1');
title('Sweep 6: bias floor on pure noise');
legend('Location','best'); ylim([0 1]); grid on;
saveas(gcf, fullfile(figures_dir, 'sim_6_samplesize.png'));

%% 6b. Sample Size — actual generative model (P1 S2 + S4 reference line)
% Same x-axis as sweep 6, but X and Y are generated from a shared latent
% Z (D_true=3) with additive noise. Plots real and shuffle CC1 for both
% linear and kernel CCA, plus the asymptotic linear-CCA reference line
% estimated empirically with N=1e5.
fprintf('--- Running Sample Size Sweep (generative model) ---\n');
sample_sizes_g = [30, 50, 100, 200, 500, 1000, 2000];
gen_noise_amp = 1.5;

% Asymptotic linear CC1 reference (large-N empirical estimate).
N_ref = 1e5;
Z_ref = randn(N_ref, D_true);
X_ref = Z_ref * W_x + gen_noise_amp * randn(N_ref, k_X);
Y_ref = Z_ref * W_y + gen_noise_amp * randn(N_ref, k_Y);
[~, ~, r_ref] = canoncorr(X_ref, Y_ref);
asymptotic_cc1 = r_ref(1);
fprintf('  asymptotic linear CC1 (noise=%.1f, N=%g) = %.3f\n', gen_noise_amp, N_ref, asymptotic_cc1);
clear X_ref Y_ref Z_ref;

cc_g_lin     = nan(length(sample_sizes_g), 1);
cc_g_lin_sh  = nan(length(sample_sizes_g), 1);
cc_g_k       = nan(length(sample_sizes_g), 1);
cc_g_k_sh    = nan(length(sample_sizes_g), 1);

for i = 1:length(sample_sizes_g)
    n_samp = sample_sizes_g(i);
    r_iters_lin    = nan(num_iterations, 1);
    r_iters_lin_sh = nan(num_iterations, 1);
    r_iters_k      = nan(num_iterations, 1);
    r_iters_k_sh   = nan(num_iterations, 1);
    for it = 1:num_iterations
        Z_g  = randn(n_samp, D_true);
        Xg = Z_g * W_x + gen_noise_amp * randn(n_samp, k_X);
        Yg = Z_g * W_y + gen_noise_amp * randn(n_samp, k_Y);
        [r_lin, r_kcca] = top_cc(Xg, Yg, 1e-3);
        r_iters_lin(it) = r_lin;
        r_iters_k(it)   = r_kcca;
        perm = randperm(n_samp);
        [r_lin_sh, r_kcca_sh] = top_cc(Xg, Yg(perm,:), 1e-3);
        r_iters_lin_sh(it) = r_lin_sh;
        r_iters_k_sh(it)   = r_kcca_sh;
    end
    cc_g_lin(i)    = mean(r_iters_lin,    'omitnan');
    cc_g_lin_sh(i) = mean(r_iters_lin_sh, 'omitnan');
    cc_g_k(i)      = mean(r_iters_k,      'omitnan');
    cc_g_k_sh(i)   = mean(r_iters_k_sh,   'omitnan');
end

figure('Name', 'Effect of Sample Size — generative model', 'Position', [400 400 600 400]);
plot(sample_sizes_g, cc_g_lin,    '-o', 'Color', [0.2 0.6 0.8], 'DisplayName','Linear real'); hold on;
plot(sample_sizes_g, cc_g_lin_sh, '--', 'Color', [0.2 0.6 0.8], 'DisplayName','Linear shuffle');
plot(sample_sizes_g, cc_g_k,      '-s', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel real');
plot(sample_sizes_g, cc_g_k_sh,   '--', 'Color', [0.8 0.4 0.2], 'DisplayName','Kernel shuffle');
xline(k_X + k_Y, 'r--', 'k1+k2', 'HandleVisibility','off');
yline(asymptotic_cc1, 'k:', sprintf('Asymptotic CC1 = %.2f', asymptotic_cc1), ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
xlabel('Sample size');
ylabel('CC1');
title(sprintf('Sweep 6b: generative model, noise=%.1f, D_{true}=%d', gen_noise_amp, D_true));
legend('Location','best'); ylim([0 1]); grid on;
saveas(gcf, fullfile(figures_dir, 'sim_6b_samplesize_generative.png'));

% Persist underlying numerical data alongside figures (CLAUDE.md viz rule).
save(fullfile(figures_dir, 'sim_data.mat'), ...
     'noise_levels','cc_noise','cc_noise_k','cc_noise_lin_sh','cc_noise_kcca_sh', ...
     'smooth_windows','cc_smooth','cc_smooth_k','cc_smooth_sh','cc_smooth_ksh', ...
     'lags_native','cc_lag','cc_lag_k','true_lag_native','subs', ...
     'dims_to_test','cc_dims','cc_dims_k','cc_dims_sh','cc_dims_ksh','D_true', ...
     'sample_sizes','cc_samples_mean','cc_samples_std','cc_samples_k_mean','cc_samples_k_std', ...
     'sample_sizes_g','cc_g_lin','cc_g_lin_sh','cc_g_k','cc_g_k_sh','gen_noise_amp','asymptotic_cc1', ...
     'res_lin','res_lin_sh','res_k','res_k_sh', ...
     '-v7.3');

fprintf('Simulations complete. Plots saved to %s\n', figures_dir);

%% Helper Functions

function [r] = kernel_cca(X, Y, kernel_type, sigma, kappa)
    % kernel_cca: Computes regularized Kernel Canonical Correlation Analysis
    % X, Y: N_samples x N_features matrices
    
    N = size(X, 1);
    
    % 1. Compute Kernel Matrices (Gram Matrices)
    if strcmp(kernel_type, 'rbf')
        distX = pdist2(X, X, 'squaredeuclidean');
        distY = pdist2(Y, Y, 'squaredeuclidean');
        Kx = exp(-distX / (2 * sigma^2));
        Ky = exp(-distY / (2 * sigma^2));
    elseif strcmp(kernel_type, 'linear')
        Kx = X * X';
        Ky = Y * Y';
    else
        error('Unsupported kernel type. Use ''rbf'' or ''linear''.');
    end
    
    % 2. Center the Kernel Matrices in Feature Space
    H = eye(N) - (1/N) * ones(N, N);
    Kx_c = H * Kx * H;
    Ky_c = H * Ky * H;
    
    % 3. Regularization (Crucial for small sample sizes in sliding windows)
    I = eye(N);
    Rx = Kx_c * Kx_c + kappa * Kx_c + 1e-8 * I; % Added slight ridge for stability
    Ry = Ky_c * Ky_c + kappa * Ky_c + 1e-8 * I;
    
    % 4. Solve the Generalized Eigenvalue Problem
    % Equation: [0, Kx_c*Ky_c; Ky_c*Kx_c, 0] * w = lambda * [Rx, 0; 0, Ry] * w
    Z = zeros(N, N);
    A = [Z, Kx_c * Ky_c; Ky_c * Kx_c, Z];
    B = [Rx, Z; Z, Ry];
    
    % Ensure B is symmetric positive definite
    B = (B + B') / 2; 
    
    try
        [~, D] = eig(A, B);
        eigenvalues = diag(D);

        % Extract real eigenvalues, sort them descending
        real_idx = imag(eigenvalues) == 0;
        evals_real = real(eigenvalues(real_idx));
        [sorted_evals, ~] = sort(evals_real, 'descend');

        % The canonical correlations are the positive eigenvalues
        r = sorted_evals(sorted_evals > 1e-4);
        r(r > 1) = 1; % Cap at 1 to prevent floating point overshoot
    catch
        r = nan; % Fail gracefully if matrix is poorly conditioned
    end
end

function [msg, pass] = check_top_cc()
% Self-validation for the `top_cc` helper.
% Returns msg = one-line summary, pass = bool.
    rng(7);
    N = 400;

    % (1) Identical views -> CC1 should be ~1.0 (linear and kernel).
    X1 = randn(N, 5);
    [r_lin1, r_kcca1] = top_cc(X1, X1, 1e-3);

    % (2) Strong nonlinear coupling: kernel should beat linear by >=0.2.
    z = linspace(-3, 3, N)' + 0.3*randn(N,1);
    X2 = z * randn(1,5) + 0.1*randn(N,5);
    Y2 = sin(z) * randn(1,5) + cos(z) * randn(1,5) + 0.1*randn(N,5);
    [r_lin2, r_kcca2] = top_cc(X2, Y2, 1e-3);

    % (3) Independent views: CC1 should be small (bias floor only).
    X3 = randn(N, 5);
    Y3 = randn(N, 5);
    [r_lin3, r_kcca3] = top_cc(X3, Y3, 1e-3);

    pass = (r_lin1 > 0.95) && (r_kcca1 > 0.85) ...
        && (r_kcca2 - r_lin2 > 0.15) ...
        && (r_lin3 < 0.5);

    msg = sprintf(['(1) identical: lin=%.3f kcca=%.3f | ' ...
                   '(2) sin/cos: lin=%.3f kcca=%.3f (kcca-lin=%+.2f) | ' ...
                   '(3) indep: lin=%.3f kcca=%.3f | pass=%d'], ...
                  r_lin1, r_kcca1, r_lin2, r_kcca2, r_kcca2 - r_lin2, ...
                  r_lin3, r_kcca3, pass);
end

function [r_lin, r_kcca] = top_cc(X, Y, kappa)
% top_cc: returns the top-1 canonical correlation under linear CCA and
% RBF kernel CCA. Sigma uses the median heuristic on stacked (X, Y).
% NaN if either fit fails or the inputs are degenerate.
%
% Used by every sweep in this script.
    r_lin  = NaN;
    r_kcca = NaN;
    if size(X,1) ~= size(Y,1) || size(X,1) <= size(X,2) + size(Y,2)
        return;
    end

    % --- Linear ---
    try
        [~,~,r] = canoncorr(X, Y);
        if ~isempty(r), r_lin = r(1); end
    catch
    end

    % --- Kernel (RBF) ---
    try
        % Median heuristic on combined data; fall back to per-view mean.
        d_combined = pdist([X, Y]);
        sigma = median(d_combined);
        if ~isfinite(sigma) || sigma == 0
            sigma = mean([median(pdist(X)), median(pdist(Y))]);
        end
        rk = kernel_cca(X, Y, 'rbf', sigma, kappa);
        if ~isempty(rk) && all(isfinite(rk(1)))
            r_kcca = rk(1);
        end
    catch
    end
end
