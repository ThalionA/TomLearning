%% cca_animal_sanity_check.m
% Detailed sanity checks on a single animal (Animal 5: TF052_export.mat)
% Performs linear CCA on a single-trial basis and plots diagnostics.

clear; clc; close all;

% Reproducibility: every random shuffle / split below uses this seed.
rng(42);

%% Setup Parameters
% Resolve project root from this script's location so the script runs
% from any cwd. Repo layout is .../TomLearning/HC_V1_Code/this_script.m
% so the project root is the parent of the script's folder.
this_file = mfilename('fullpath');
[script_dir, ~] = fileparts(this_file);
base_dir       = fullfile(script_dir, '..');
data_dir       = fullfile(base_dir, 'HC_V1_data');
figures_dir    = fullfile(base_dir, 'HC_V1_figures', 'CCA_Animal5');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

filename       = 'TF052_export.mat'; % Animal 5
target_bin_ms  = 25;
min_speed_cms  = 2;
half_window_bins = 10; % +-10 bins = 21 bin window = 525 ms (at 25ms bins)
window_bins    = -half_window_bins:half_window_bins;
n_window       = numel(window_bins);

% Areas to test
area1 = 'CA1';
area2 = 'V1';

%% 0. Synthetic-control TDD test (P1 C3)
% Run the same sliding-window CCA path on synthetic two-region data with
% a known shared latent. Confirms (a) the pipeline recovers a real CC1
% above the shuffle null when one exists, and (b) it does NOT find one
% when X and Y are independent. Errors out before touching real data
% if either fails -- catches regressions early.
run_synthetic_control = false;
if run_synthetic_control
    [synth_pass, synth_msg] = synthetic_cca_check(half_window_bins, 5, 5);
    fprintf('[Synthetic control] %s\n', synth_msg);
    if ~synth_pass
        error('Synthetic-control test failed; aborting before real data.');
    end
end

%% 1. Load Data
fullpath = fullfile(data_dir, filename);
fprintf('Loading data for %s...\n', filename);
D = load(fullpath);

if ~isfield(D, 'units') || ~isfield(D, 'binned_spikes')
    error('Missing required fields in data file.');
end

units = D.units;
bs_native = D.binned_spikes; 
native_bin_s = double(D.params_main.bin_size);
rebin_factor = round((target_bin_ms / 1000) / native_bin_s);

fprintf('Rebinning to %d ms (factor = %d)...\n', target_bin_ms, rebin_factor);

%% 2. Data Preparation (Binning)
bs_layout_units_first = (size(bs_native, 1) == length(units.unit_id));
if bs_layout_units_first
    nUnits = size(bs_native, 1);
    T      = size(bs_native, 2);
else
    nUnits = size(bs_native, 2);
    T      = size(bs_native, 1);
end

T_use  = floor(T / rebin_factor) * rebin_factor;
T_re   = T_use / rebin_factor;

spikes_re = zeros(T_re, nUnits);
for u = 1:nUnits
    if bs_layout_units_first
        row = double(bs_native(u, 1:T_use));
    else
        row = double(bs_native(1:T_use, u))';
    end
    spikes_re(:, u) = sum(reshape(row, rebin_factor, T_re), 1)';
end

vel_gf = double(D.data_behaviour.velocity_binned_gf(1:T_use));
tr_cued = double(D.data_behaviour.trial_binned_cued(1:T_use));
pos_cm = double(D.data_behaviour.pos_binned_gf(1:T_use));
mask_cued = logical(D.analysis_behaviour.masks.tunnel_cued(1:T_use));

% Lick rate handling
if isfield(D.data_behaviour, 'lick_rate_binned_gf')
    licks_gf = double(D.data_behaviour.lick_rate_binned_gf);
elseif isfield(D.data_behaviour, 'licks_binned_gf')
    licks_gf = double(D.data_behaviour.licks_binned_gf);
    licks_gf = licks_gf / native_bin_s; % convert 0-1 counts to Hz
elseif isfield(D.data_behaviour, 'licks_binned')
    licks_gf = double(D.data_behaviour.licks_binned);
    licks_gf = licks_gf / native_bin_s; % convert 0-1 counts to Hz
else
    licks_gf = nan(1, T);
end

% Smooth with 50ms Gaussian kernel
if ~all(isnan(licks_gf))
    sigma_s = 50 / 1000;
    sigma_bins = sigma_s / native_bin_s;
    x_win = -ceil(4*sigma_bins):ceil(4*sigma_bins);
    gauss_kernel = exp(-(x_win.^2) / (2 * sigma_bins^2));
    gauss_kernel = gauss_kernel / sum(gauss_kernel);
    licks_gf = conv(licks_gf, gauss_kernel, 'same');
end
licks_gf = licks_gf(1:T_use);

% Lick ratio per trial handling
if isfield(D.analysis_behaviour, 'lick_ratio')
    lick_ratio_all = D.analysis_behaviour.lick_ratio(:);
else
    lick_ratio_all = nan(max(tr_cued), 1);
end

vel_re   = mean(reshape(vel_gf,   rebin_factor, T_re), 1)';
tc_re    = mode(reshape(tr_cued,  rebin_factor, T_re), 1)';
pos_re   = mean(reshape(pos_cm,   rebin_factor, T_re), 1)';
mc_re    = all(reshape(mask_cued, rebin_factor, T_re), 1)';
licks_re = mean(reshape(licks_gf, rebin_factor, T_re), 1)';

% Validity mask (cued period AND min speed)
base_valid = mc_re & (vel_re >= min_speed_cms);

%% 2.1 Contiguous-valid-block-length distribution (P1 A4)
% How long are the unbroken stretches of (cued tunnel) AND (>=min_speed)?
% This is what the strict-contiguity sliding window can actually use.
% Compute per trial so cross-trial gaps don't inflate the lengths.
trial_ids = unique(tc_re(tc_re > 0));
block_lengths_bins = [];
for tt = 1:length(trial_ids)
    bv = base_valid(tc_re == trial_ids(tt));
    bv = bv(:);
    % Find runs of true.
    d = diff([0; bv; 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    block_lengths_bins = [block_lengths_bins; ends - starts + 1]; %#ok<AGROW>
end
block_lengths_ms = block_lengths_bins * target_bin_ms;

figure('Name','Contiguous valid-block lengths','Position',[120 120 700 380]);
subplot(1,2,1);
histogram(block_lengths_bins, 'BinWidth', 1, 'FaceColor', [0.2 0.5 0.8]);
xline(2*half_window_bins+1, 'r--', 'LineWidth', 1.5);
xlabel('Block length (50 ms bins)'); ylabel('Count');
title(sprintf('Per-trial contiguous valid blocks  (n=%d)', numel(block_lengths_bins)));
legend({'Blocks', sprintf('Window size = %d', 2*half_window_bins+1)}, 'Location','best');
grid on;
subplot(1,2,2);
histogram(block_lengths_ms, 'BinWidth', 100, 'FaceColor', [0.2 0.5 0.8]);
xline((2*half_window_bins+1)*target_bin_ms, 'r--', 'LineWidth', 1.5);
xlabel('Block length (ms)'); ylabel('Count');
title('Same in ms');
grid on;
saveas(gcf, fullfile(figures_dir, '00_block_lengths.png'));
fprintf('Block-length stats: median=%.0f bins (%.0f ms), 90th=%.0f bins, %.1f%% blocks >= window size.\n', ...
    median(block_lengths_bins), median(block_lengths_ms), ...
    quantile(block_lengths_bins, 0.9), ...
    100 * mean(block_lengths_bins >= (2*half_window_bins+1)));

%% 3. Load Pure Spatial Data
if isfield(D.analysis_spatial, 'firing') && isfield(D.analysis_spatial.firing, 'cued')
    raw_sp = D.analysis_spatial.firing.cued.freq_z;
    sp_data = permute(raw_sp, [1, 3, 2]); % [units x bins x trials]
else
    error('Spatial firing data not found.');
end

%% 4. Extract Regions and Plot PCA Variance
fprintf('Extracting %s and %s...\n', area1, area2);
if isfield(units, 'region')
    unit_regions = units.region(:);
elseif isfield(units, 'regions_label')
    unit_regions = units.regions_label(:);
end

%% 4.0 Per-region unit stats (P2 A5)
% Show how dense the spike data actually is per region, before PCA. Sparse
% regions can fool PCA into picking near-empty PCs.
all_region_names = unique(unit_regions);
fprintf('\n[Per-region stats]\n');
fprintf('  %-8s %-8s %-12s %-12s\n', 'region', 'n_units', 'rate (Hz)', 'frac_zero');
unit_stats = struct();
for ir = 1:length(all_region_names)
    nm = char(all_region_names{ir});
    if isempty(nm), continue; end
    msk = strcmp(unit_regions, nm);
    n_u = sum(msk);
    if n_u == 0, continue; end
    region_spikes = spikes_re(:, msk);              % [T_re x n_u]
    rate_hz   = mean(region_spikes(:)) / (target_bin_ms/1000);
    frac_zero = mean(region_spikes(:) == 0);
    unit_stats.(matlab.lang.makeValidName(nm)) = struct('n', n_u, 'rate_hz', rate_hz, 'frac_zero', frac_zero);
    fprintf('  %-8s %-8d %-12.3f %-12.3f\n', nm, n_u, rate_hz, frac_zero);
end

% Bar plot of rate and zero-fraction per region.
figure('Name','Per-region unit stats','Position',[140 140 800 360]);
reg_names_kept = fieldnames(unit_stats);

% Above is awkward; use explicit loop.
nR = numel(reg_names_kept);
n_units_vec = zeros(nR,1); rate_vec = zeros(nR,1); zero_frac_vec = zeros(nR,1);
for ir = 1:nR
    n_units_vec(ir)  = unit_stats.(reg_names_kept{ir}).n;
    rate_vec(ir)     = unit_stats.(reg_names_kept{ir}).rate_hz;
    zero_frac_vec(ir)= unit_stats.(reg_names_kept{ir}).frac_zero;
end
subplot(1,3,1); bar(n_units_vec, 'FaceColor', [0.2 0.4 0.7]);
xticks(1:nR); xticklabels(reg_names_kept); ylabel('# units');
title('Units per region'); grid on;
subplot(1,3,2); bar(rate_vec, 'FaceColor', [0.2 0.6 0.4]);
xticks(1:nR); xticklabels(reg_names_kept); ylabel('Mean rate (Hz)');
title(sprintf('Firing rate (post %d ms rebin)', target_bin_ms)); grid on;
subplot(1,3,3); bar(zero_frac_vec, 'FaceColor', [0.7 0.4 0.2]);
xticks(1:nR); xticklabels(reg_names_kept); ylabel('Fraction zero bins');
title('Sparsity (fraction of zero-spike bins)'); grid on;
saveas(gcf, fullfile(figures_dir, '00b_region_stats.png'));

% Extract full PCA space (Temporal)
[X_all, coeff_X, exp_X, u_mask_X] = get_pca_data(area1, unit_regions, spikes_re, base_valid);
[Y_all, coeff_Y, exp_Y, u_mask_Y] = get_pca_data(area2, unit_regions, spikes_re, base_valid);

% Extract full PCA space (Spatial)
[X_sp_all, ~, ~, ~] = get_spatial_pca_data(area1, unit_regions, sp_data, 5);
[Y_sp_all, ~, ~, ~] = get_spatial_pca_data(area2, unit_regions, sp_data, 5);

% Define dimensionalities to compare
% CRITICAL: CCA requires N_samples > k_x + k_y. With a window of n_window, we must cap K.
max_k = floor(n_window / 3);
fprintf('NOTE: Window size is %d bins. Capping max PCA components to %d (1/3 of bins) to prevent overfitting.\n', n_window, max_k);

k5_X = min(5, max_k); k5_Y = min(5, max_k);
fprintf('Final CCA K threshold: 5 (capped to %d if necessary)\n', max_k);

figure('Name', 'PCA Explained Variance', 'Position', [150 150 800 400]);
subplot(1,2,1);
plot(exp_X, '-o', 'LineWidth', 2);
xlabel('Principal Component'); ylabel('Cumulative Variance (%)');
title(sprintf('%s (N=%d units)', area1, sum(u_mask_X)));
grid on;

subplot(1,2,2);
plot(exp_Y, '-o', 'LineWidth', 2);
xlabel('Principal Component'); ylabel('Cumulative Variance (%)');
title(sprintf('%s (N=%d units)', area2, sum(u_mask_Y)));
grid on;
saveas(gcf, fullfile(figures_dir, '01_pca_variance.png'));

%% 5. Sliding-Window CCA on Trials (Temporal vs Spatial)
fprintf('Running CCA on trials...\n');
trials_all = unique(tc_re(tc_re > 0));
trials_to_plot = [10, 20, 30, 40, 50]; % Plot traces for these trials

% Drop the head of each trial: animals typically don't fully engage with
% the cued tunnel for the first few seconds. Trim these as a generic
% "warm-up" guard. Tune `trial_warmup_ms` if your protocol differs.
trial_warmup_ms = 5000;
bins_to_remove  = ceil(trial_warmup_ms / target_bin_ms);

n_sp_bins = 200;
bin_size_cm = 2.5;
x_axis_cm = (1:n_sp_bins) * bin_size_cm;

% Data structures for scatter plots
scatter_vel  = [];
scatter_lick = [];
scatter_cc1  = [];
scatter_cc1_sh = [];   % paired shuffle CC1 per window (for distribution plot)

% Data structures for Method Comparison Scatters
comp_t_cc = [];
comp_p_cc = [];
comp_t_ifi = [];
comp_p_ifi = [];

% Per-trial held-out CC1 (P0 A3): one scalar per trial.
% We pool all valid windows in a trial, randomly split 70/30, fit canoncorr
% on the train half, project the test half, report the held-out canonical
% correlation. This is the honest small-sample-bias-corrected number.
heldout_cc1 = nan(length(trials_all), 1);
heldout_cc1_sh = nan(length(trials_all), 1);
trainset_cc1 = nan(length(trials_all), 1);  % reference: training-set CC1

% Data structures for full spatial profiles
all_cc_temp_sp = nan(length(trials_all), n_sp_bins);
all_cc_temp_sp_sh = nan(length(trials_all), n_sp_bins);  % shuffle null
all_ifi_temp_sp = nan(length(trials_all), n_sp_bins);
all_cc_pure_sp = nan(length(trials_all), n_sp_bins);
all_ifi_pure_sp = nan(length(trials_all), n_sp_bins);

fprintf('\nProcessing %d trials...\n', length(trials_all));
for i = 1:length(trials_all)
    tr = trials_all(i);
    if mod(i, 10) == 0
        fprintf('  trial %d / %d\n', i, length(trials_all));
    end
    trial_idx = find(tc_re == tr);
    
    % Truncate first 5 seconds completely
    if length(trial_idx) > bins_to_remove
        trial_idx = trial_idx((bins_to_remove + 1):end);
    else
        continue; % Trial is too short
    end
    
    n_tr_bins = length(trial_idx);
    if n_tr_bins < n_window + 5
        continue;
    end
    
    X_tr = X_all(trial_idx, :);
    Y_tr = Y_all(trial_idx, :);
    vel_tr = vel_re(trial_idx);
    pos_tr = pos_re(trial_idx);
    lick_tr = licks_re(trial_idx);
    valid_tr = base_valid(trial_idx);
    
    % ----------------------------------------------------
    % A. TEMPORAL METHOD
    % ----------------------------------------------------
    cc_5  = nan(n_tr_bins, 1);
    cc_sh_5_mean = nan(n_tr_bins, 1);
    cc_sh_5_sem  = nan(n_tr_bins, 1);
    ifi_tr= nan(n_tr_bins, 1);
    
    for ic = (half_window_bins+1):(n_tr_bins - half_window_bins)
        idx_w = (ic - half_window_bins):(ic + half_window_bins);
        if ~all(valid_tr(idx_w)), continue; end
        
        X_w = X_tr(idx_w, :);
        Y_w = Y_tr(idx_w, :);
        
        % Check rank using max dimensionalities
        if rank(X_w(:, 1:k5_X)) < k5_X || rank(Y_w(:, 1:k5_Y)) < k5_Y, continue; end
        
        % CC at K=5
        [A, B, r] = canoncorr(X_w(:, 1:k5_X), Y_w(:, 1:k5_Y));
        if ~isempty(r)
            cc_5(ic) = r(1);
            
            % Information Flow Index (IFI)
            u = (X_w(:, 1:k5_X) - mean(X_w(:, 1:k5_X))) * A(:,1);
            v = (Y_w(:, 1:k5_Y) - mean(Y_w(:, 1:k5_Y))) * B(:,1);
            lags = -3:3;
            r_lags = nan(1, length(lags));
            for il = 1:length(lags)
                L = lags(il);
                if L >= 0
                    a = u(1:(n_window-L)); b = v((1+L):n_window);
                else
                    a = u((1-L):n_window); b = v(1:(n_window+L));
                end
                if numel(a) > 5, r_lags(il) = corr(a, b); end
            end
            % Signed-lag IFI with symmetric denominator guard.
            % Using signed averages (not |r_lags|) preserves directionality:
            % strongly anti-correlated negative-lag and positive-lag pairs
            % no longer alias to the same sign. Guard on |neg+pos| so the
            % ratio doesn't blow up when neg ~= -pos.
            r_neg = abs(mean(r_lags(lags < 0), 'omitnan'));
            r_pos = abs(mean(r_lags(lags > 0), 'omitnan'));
            denom = r_neg + r_pos;
            if abs(denom) > 1e-3
                ifi_tr(ic) = (r_neg - r_pos) / denom;
            end

            % Shuffle Control (20 iterations of random permutation)
            n_sh = 20;
            r_sh_5_iter = nan(n_sh, 1);
            for sh = 1:n_sh
                perm = randperm(n_window);
                [~,~,rs5] = canoncorr(X_w(:, 1:k5_X), Y_w(perm, 1:k5_Y));
                if ~isempty(rs5), r_sh_5_iter(sh) = rs5(1); end
            end
            cc_sh_5_mean(ic) = mean(r_sh_5_iter, 'omitnan');
            cc_sh_5_sem(ic)  = std(r_sh_5_iter, 'omitnan') / sqrt(sum(~isnan(r_sh_5_iter)));
        end
    end
    
    % Accumulate valid scatter data (using K=5 as the main CC measure)
    valid_cc = ~isnan(cc_5);
    if any(valid_cc)
        scatter_cc1     = [scatter_cc1;    cc_5(valid_cc)];
        scatter_cc1_sh  = [scatter_cc1_sh; cc_sh_5_mean(valid_cc)];
        scatter_vel     = [scatter_vel;    vel_tr(valid_cc)];
        scatter_lick    = [scatter_lick;   lick_tr(valid_cc)];
    end

    % Posthoc Spatial Binning of Temporal CC1/IFI (real and shuffle)
    cc_temp_sp    = nan(1, n_sp_bins);
    cc_temp_sp_sh = nan(1, n_sp_bins);
    ifi_temp_sp   = nan(1, n_sp_bins);
    pos_bin_idx = ceil(pos_tr / bin_size_cm);
    pos_bin_idx(pos_bin_idx < 1) = 1;
    pos_bin_idx(pos_bin_idx > n_sp_bins) = n_sp_bins;

    for b = 1:n_sp_bins
        idx_b = (pos_bin_idx == b);
        if any(idx_b)
            cc_temp_sp(b)    = mean(cc_5(idx_b),         'omitnan');
            cc_temp_sp_sh(b) = mean(cc_sh_5_mean(idx_b), 'omitnan');
            ifi_temp_sp(b)   = mean(ifi_tr(idx_b),       'omitnan');
        end
    end

    % --- Held-out CC1 (P0 A3) ---
    % Pool all valid window samples for this trial, fit canoncorr on a
    % random 70%, project the 30% held out, report Pearson correlation
    % of the projected canonical variates. Compare against shuffle.
    [heldout_cc1(i), trainset_cc1(i), heldout_cc1_sh(i)] = ...
        compute_heldout_cc1(X_tr, Y_tr, valid_tr, half_window_bins, ...
                            k5_X, k5_Y, 0.7, 20);
    
    % ----------------------------------------------------
    % B. PURE SPATIAL METHOD (Multi-Trial Pooled)
    % ----------------------------------------------------
    cc_pure_sp = nan(1, n_sp_bins);
    ifi_pure_sp = nan(1, n_sp_bins);
    
    tr_win = (tr - 3) : (tr + 3);
    tr_win = tr_win(tr_win >= 1 & tr_win <= size(X_sp_all, 3));
    
    for b = 4:(n_sp_bins - 3)
        b_win = (b - 3) : (b + 3);
        
        X_w_sp = X_sp_all(:, b_win, tr_win); % [5 x bins x trials]
        Y_w_sp = Y_sp_all(:, b_win, tr_win);
        
        x = reshape(X_w_sp, 5, [])'; % [samples x 5]
        y = reshape(Y_w_sp, 5, [])'; 
        
        valid_s = ~any(isnan(x), 2) & ~any(isnan(y), 2);
        
        if sum(valid_s) > 12 && rank(x(valid_s,:)) >= 5 && rank(y(valid_s,:)) >= 5
            [A, B, r] = canoncorr(x(valid_s,:), y(valid_s,:));
            cc_pure_sp(b) = r(1);
            
            u = (x(valid_s,:) - mean(x(valid_s,:))) * A(:,1);
            v = (y(valid_s,:) - mean(y(valid_s,:))) * B(:,1);
            lags = -3:3;
            r_lags = nan(1, length(lags));
            n_samp = sum(valid_s);
            for il = 1:length(lags)
                L = lags(il);
                if L >= 0
                    aa = u(1:(n_samp-L)); bb = v((1+L):n_samp);
                else
                    aa = u((1-L):n_samp); bb = v(1:(n_samp+L));
                end
                if numel(aa) > 5, r_lags(il) = corr(aa, bb); end
            end
            % Signed IFI with symmetric guard (matches temporal method).
            r_neg = abs(mean(r_lags(lags < 0), 'omitnan'));
            r_pos = abs(mean(r_lags(lags > 0), 'omitnan'));
            denom = r_neg + r_pos;
            if abs(denom) > 1e-3
                ifi_pure_sp(b) = (r_neg - r_pos) / denom;
            end
        end
    end
    
    % Method Comparison Trial Averages
    if sum(~isnan(cc_temp_sp)) > 10 && sum(~isnan(cc_pure_sp)) > 10
        comp_t_cc = [comp_t_cc; mean(cc_temp_sp, 'omitnan')];
        comp_p_cc = [comp_p_cc; mean(cc_pure_sp, 'omitnan')];
        comp_t_ifi = [comp_t_ifi; mean(ifi_temp_sp, 'omitnan')];
        comp_p_ifi = [comp_p_ifi; mean(ifi_pure_sp, 'omitnan')];
    end
    
    % Store full spatial curves
    all_cc_temp_sp(i, :)    = cc_temp_sp;
    all_cc_temp_sp_sh(i, :) = cc_temp_sp_sh;
    all_ifi_temp_sp(i, :)   = ifi_temp_sp;
    all_cc_pure_sp(i, :)    = cc_pure_sp;
    all_ifi_pure_sp(i, :)   = ifi_pure_sp;
    
    % Plot Single Trial Method Comparison if requested
    if ismember(tr, trials_to_plot)
        time_ax = (1:n_tr_bins) * target_bin_ms / 1000;
        
        figure('Name', sprintf('Trial %d Method Comparison', tr), 'Position', [200 200 800 800]);
        
        subplot(3,1,1);
        plot(time_ax, pos_tr, 'k', 'LineWidth', 2);
        ylabel('Position (cm)'); xlabel('Time (s)');
        title(sprintf('Trial %d Temporal Context (First 5s cropped)', tr));
        grid on;
        
        subplot(3,1,2);
        plot(x_axis_cm, cc_temp_sp,    'Color', [0 0.6 0], 'LineWidth', 2, ...
            'DisplayName', 'Temporal real');
        hold on;
        plot(x_axis_cm, cc_temp_sp_sh, 'Color', [0 0.6 0], 'LineStyle', ':', 'LineWidth', 1.5, ...
            'DisplayName', 'Temporal shuffle');
        plot(x_axis_cm, cc_pure_sp,    'Color', [0 0 0.8], 'LineStyle', '--', 'LineWidth', 2, ...
            'DisplayName', 'Pure spatial real');
        ylim([0 1.05]);
        ylabel('CC1'); xlabel('Position (cm)');
        legend('Location', 'best');
        title(sprintf('Trial %d: CC1 by position (real vs shuffle)', tr));
        grid on;
        
        subplot(3,1,3);
        plot(x_axis_cm, ifi_temp_sp, 'm', 'LineWidth', 2, 'DisplayName', 'Temporal IFI');
        hold on;
        plot(x_axis_cm, ifi_pure_sp, 'Color', [0.5 0 0.5], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Pure Spatial IFI');
        ylim([-1 1]);
        yline(0, 'k--');
        ylabel('IFI'); xlabel('Position (cm)');
        legend('Location', 'best');
        title('IFI Method Comparison');
        grid on;
        
        saveas(gcf, fullfile(figures_dir, sprintf('02_trial_%d_comparison.png', tr)));
        fprintf('  Processed and saved method comparison for Trial %d\n', tr);
    end
end

%% 6. Summary Scatters
fprintf('Generating summary scatter plots...\n');

% Scatter: Temporal Metrics vs Behavior
figure('Name', 'Temporal CC1 vs Behavior', 'Position', [250 250 800 400]);

subplot(1,2,1);
scatter(scatter_vel, scatter_cc1, 10, 'filled', 'MarkerFaceAlpha', 0.1); hold on;
robust_fit_line(scatter_vel, scatter_cc1);
xlabel('Velocity (cm/s)'); ylabel('Temporal CC1 (K=5)');
title('Per-Bin CC1 vs Velocity');
grid on;

subplot(1,2,2);
scatter(scatter_lick, scatter_cc1, 10, 'filled', 'MarkerFaceAlpha', 0.1); hold on;
robust_fit_line(scatter_lick, scatter_cc1);
xlabel('Lick Rate (Hz)'); ylabel('Temporal CC1 (K=5)');
title('Per-Bin CC1 vs Lick Rate');
grid on;

saveas(gcf, fullfile(figures_dir, '03_behavior_scatters.png'));

% Scatter: Method Comparison
figure('Name', 'Method Comparison Scatter', 'Position', [300 300 1000 450]);

subplot(1,2,1);
scatter(comp_t_cc, comp_p_cc, 40, 'filled', 'MarkerFaceAlpha', 0.8); hold on;
robust_fit_line(comp_t_cc, comp_p_cc);
min_val = min([comp_t_cc; comp_p_cc]); max_val = max([comp_t_cc; comp_p_cc]);
plot([min_val, max_val], [min_val, max_val], 'k--', 'LineWidth', 1);
xlabel('Temporal Method (Trial Avg CC1)'); ylabel('Spatial Method (Trial Avg CC1)');
title('Per-Trial CC1 Comparison');
grid on;

subplot(1,2,2);
scatter(comp_t_ifi, comp_p_ifi, 40, 'filled', 'MarkerFaceAlpha', 0.8); hold on;
robust_fit_line(comp_t_ifi, comp_p_ifi);
min_val = min([comp_t_ifi; comp_p_ifi]); max_val = max([comp_t_ifi; comp_p_ifi]);
plot([min_val, max_val], [min_val, max_val], 'k--', 'LineWidth', 1);
xlabel('Temporal Method (Trial Avg IFI)'); ylabel('Spatial Method (Trial Avg IFI)');
title('Per-Trial IFI Comparison');
grid on;

saveas(gcf, fullfile(figures_dir, '04_method_comparison.png'));

% -----------------------------------------------------------
% Plot: Trial-Averaged Spatial Profiles
% -----------------------------------------------------------
figure('Name', 'Trial-Averaged Spatial Profiles', 'Position', [350 350 800 600]);

% CC1 Profile (real and shuffle)
subplot(2,1,1);
mu_cc_t  = mean(all_cc_temp_sp, 1, 'omitnan');
sem_cc_t = std(all_cc_temp_sp, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_cc_temp_sp), 1));
mu_cc_t_sh  = mean(all_cc_temp_sp_sh, 1, 'omitnan');
sem_cc_t_sh = std(all_cc_temp_sp_sh, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_cc_temp_sp_sh), 1));
mu_cc_p  = mean(all_cc_pure_sp, 1, 'omitnan');
sem_cc_p = std(all_cc_pure_sp, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_cc_pure_sp), 1));

valid_t = ~isnan(mu_cc_t);
if any(valid_t)
    x_t = x_axis_cm(valid_t); m_t = mu_cc_t(valid_t); s_t = sem_cc_t(valid_t);
    patch([x_t, fliplr(x_t)], [m_t + s_t, fliplr(m_t - s_t)], [0 0.6 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
    plot(x_t, m_t, 'Color', [0 0.6 0], 'LineWidth', 2, 'DisplayName', 'Temporal real');
end

valid_t_sh = ~isnan(mu_cc_t_sh);
if any(valid_t_sh)
    x_ts = x_axis_cm(valid_t_sh); m_ts = mu_cc_t_sh(valid_t_sh); s_ts = sem_cc_t_sh(valid_t_sh);
    patch([x_ts, fliplr(x_ts)], [m_ts + s_ts, fliplr(m_ts - s_ts)], [0 0.6 0], 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
    plot(x_ts, m_ts, 'Color', [0 0.6 0], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Temporal shuffle');
end

valid_p = ~isnan(mu_cc_p);
if any(valid_p)
    x_p = x_axis_cm(valid_p); m_p = mu_cc_p(valid_p); s_p = sem_cc_p(valid_p);
    patch([x_p, fliplr(x_p)], [m_p + s_p, fliplr(m_p - s_p)], [0 0 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
    plot(x_p, m_p, 'Color', [0 0 0.8], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Pure spatial real');
end
ylabel('CC1'); xlabel('Position (cm)');
title('Trial-averaged CC1 profile (real vs shuffle)');
legend('Location', 'best');
grid on;

% IFI Profile
subplot(2,1,2);
mu_ifi_t = mean(all_ifi_temp_sp, 1, 'omitnan');
sem_ifi_t = std(all_ifi_temp_sp, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_ifi_temp_sp), 1));
mu_ifi_p = mean(all_ifi_pure_sp, 1, 'omitnan');
sem_ifi_p = std(all_ifi_pure_sp, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_ifi_pure_sp), 1));

valid_it = ~isnan(mu_ifi_t);
if any(valid_it)
    x_it = x_axis_cm(valid_it); m_it = mu_ifi_t(valid_it); s_it = sem_ifi_t(valid_it);
    patch([x_it, fliplr(x_it)], [m_it + s_it, fliplr(m_it - s_it)], 'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
    plot(x_it, m_it, 'm', 'LineWidth', 2, 'DisplayName', 'Temporal IFI');
end

valid_ip = ~isnan(mu_ifi_p);
if any(valid_ip)
    x_ip = x_axis_cm(valid_ip); m_ip = mu_ifi_p(valid_ip); s_ip = sem_ifi_p(valid_ip);
    patch([x_ip, fliplr(x_ip)], [m_ip + s_ip, fliplr(m_ip - s_ip)], [0.5 0 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
    plot(x_ip, m_ip, 'Color', [0.5 0 0.5], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Pure Spatial IFI');
end
ylabel('IFI'); xlabel('Position (cm)');
title('Trial-Averaged IFI Profile');
legend('Location', 'best');
grid on;

saveas(gcf, fullfile(figures_dir, '05_spatial_profiles.png'));

% -----------------------------------------------------------
% Plot: Bin-by-Bin Pooled Scatter
% -----------------------------------------------------------
figure('Name', 'Bin-by-Bin Pooled Scatter', 'Position', [400 400 1000 450]);

flat_cc_t = all_cc_temp_sp(:);
flat_cc_p = all_cc_pure_sp(:);
valid_flat_cc = ~isnan(flat_cc_t) & ~isnan(flat_cc_p);

subplot(1,2,1);
scatter(flat_cc_t(valid_flat_cc), flat_cc_p(valid_flat_cc), 5, 'filled', 'MarkerFaceAlpha', 0.1); hold on;
robust_fit_line(flat_cc_t(valid_flat_cc), flat_cc_p(valid_flat_cc));
min_val = min([flat_cc_t(valid_flat_cc); flat_cc_p(valid_flat_cc)]);
max_val = max([flat_cc_t(valid_flat_cc); flat_cc_p(valid_flat_cc)]);
plot([min_val, max_val], [min_val, max_val], 'k--', 'LineWidth', 1);
xlabel('Temporal Method (Bin CC1)'); ylabel('Spatial Method (Bin CC1)');
title('Bin-by-Bin Pooled CC1');
grid on;

flat_ifi_t = all_ifi_temp_sp(:);
flat_ifi_p = all_ifi_pure_sp(:);
valid_flat_ifi = ~isnan(flat_ifi_t) & ~isnan(flat_ifi_p);

subplot(1,2,2);
scatter(flat_ifi_t(valid_flat_ifi), flat_ifi_p(valid_flat_ifi), 5, 'filled', 'MarkerFaceAlpha', 0.1); hold on;
robust_fit_line(flat_ifi_t(valid_flat_ifi), flat_ifi_p(valid_flat_ifi));
min_val = min([flat_ifi_t(valid_flat_ifi); flat_ifi_p(valid_flat_ifi)]);
max_val = max([flat_ifi_t(valid_flat_ifi); flat_ifi_p(valid_flat_ifi)]);
plot([min_val, max_val], [min_val, max_val], 'k--', 'LineWidth', 1);
xlabel('Temporal Method (Bin IFI)'); ylabel('Spatial Method (Bin IFI)');
title('Bin-by-Bin Pooled IFI');
grid on;

saveas(gcf, fullfile(figures_dir, '06_bin_by_bin_scatter.png'));

% -----------------------------------------------------------
% Plot 07: Pooled CC1 distribution — real vs shuffle (P0 A1, P1 A6)
% Histograms over every valid sliding window in the session for the
% selected pair. The shuffle distribution shows the small-sample
% bias floor; if real overlaps it, there's no genuine coupling
% above the bias.
% -----------------------------------------------------------
figure('Name', 'CC1 distribution: real vs shuffle', 'Position', [450 450 700 400]);
hold on;
edges = linspace(0, 1, 41);
if ~isempty(scatter_cc1)
    histogram(scatter_cc1,    edges, 'Normalization','pdf', ...
              'FaceColor',[0 0.6 0],     'FaceAlpha', 0.55, 'EdgeColor','none', ...
              'DisplayName','Real');
end
if ~isempty(scatter_cc1_sh)
    histogram(scatter_cc1_sh, edges, 'Normalization','pdf', ...
              'FaceColor',[0.5 0.5 0.5], 'FaceAlpha', 0.55, 'EdgeColor','none', ...
              'DisplayName','Shuffle');
end
xline(median(scatter_cc1,    'omitnan'), '-',  'Color', [0 0.4 0],   'LineWidth', 1.5, 'HandleVisibility','off');
xline(median(scatter_cc1_sh, 'omitnan'), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel(sprintf('CC1 per window (%s vs %s)', area1, area2));
ylabel('PDF (windows pooled across trials)');
title(sprintf('Per-window CC1 distribution: real (n=%d) vs shuffle  [bias floor]', sum(~isnan(scatter_cc1))));
legend('Location','best');
xlim([0 1]); grid on;
saveas(gcf, fullfile(figures_dir, '07_cc1_distribution.png'));

% -----------------------------------------------------------
% Plot 08: Per-trial held-out CC1 (P0 A3)
% For each trial we pool all valid window samples, fit canoncorr on a
% random 70%, and evaluate canonical correlation on the held-out 30%.
% This is the bias-corrected number — the per-window mean overstates it.
% -----------------------------------------------------------
figure('Name', 'Held-out CC1 per trial', 'Position', [500 500 800 450]);
hold on;
n_trials_local = length(trials_all);
plot(trials_all, trainset_cc1,    '-o', 'Color', [0 0.6 0],     'LineWidth', 1.2, 'MarkerSize', 4, ...
     'DisplayName','Training-set CC1 (biased high)');
plot(trials_all, heldout_cc1,     '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'MarkerSize', 5, ...
     'DisplayName','Held-out CC1 (honest)');
plot(trials_all, heldout_cc1_sh,  ':',  'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, ...
     'DisplayName','Held-out CC1 (shuffle null)');
xlabel('Trial number');
ylabel('CC1 (canonical correlation)');
title(sprintf('Held-out CC1 per trial (%s vs %s, k=%d, 70/30 split)', area1, area2, k5_X));
legend('Location','best');
ylim([-0.2 1]); grid on;
saveas(gcf, fullfile(figures_dir, '08_heldout_cc1.png'));

% -----------------------------------------------------------
% Plot 09: Window-size sensitivity sweep (P1 A7)
% Repeat the sliding-window CCA at multiple window sizes, plotting
% mean CC1 (real) and mean CC1 (shuffle) vs window size. If both
% increase together as windows shrink, that's the small-sample
% bias floor; if real grows faster, you have signal.
% -----------------------------------------------------------
sweep_half_w   = [5, 10, 15, 20];           % bins
sweep_n_trials = min(30, length(trials_all)); % subset for speed
sweep_trials   = trials_all(1:sweep_n_trials);
sweep_real_mu  = nan(length(sweep_half_w), 1);
sweep_real_sem = nan(length(sweep_half_w), 1);
sweep_shuf_mu  = nan(length(sweep_half_w), 1);
sweep_shuf_sem = nan(length(sweep_half_w), 1);
fprintf('\nWindow-size sweep over %d trials...\n', sweep_n_trials);
for iw = 1:length(sweep_half_w)
    h = sweep_half_w(iw);
    n_w_local = 2*h + 1;
    if n_w_local <= k5_X + k5_Y + 1
        fprintf('  half_w=%d skipped (too small for k1+k2=%d).\n', h, k5_X+k5_Y);
        continue;
    end
    cc_pool    = [];
    cc_sh_pool = [];
    for tr_i = 1:sweep_n_trials
        tr = sweep_trials(tr_i);
        ti = find(tc_re == tr);
        if length(ti) > bins_to_remove, ti = ti((bins_to_remove+1):end); end
        if length(ti) < n_w_local + 5, continue; end
        Xt = X_all(ti, 1:k5_X);
        Yt = Y_all(ti, 1:k5_Y);
        vt = base_valid(ti);
        for ic = (h+1):(length(ti) - h)
            iw_w = (ic-h):(ic+h);
            if ~all(vt(iw_w)), continue; end
            Xw = Xt(iw_w, :); Yw = Yt(iw_w, :);
            if rank(Xw - mean(Xw,1)) < k5_X || rank(Yw - mean(Yw,1)) < k5_Y, continue; end
            try
                [~,~,r_real] = canoncorr(Xw, Yw);
            catch
                continue;
            end
            cc_pool(end+1,1) = r_real(1); %#ok<AGROW>
            % Single shuffle (random permutation) per window for this sweep.
            try
                [~,~,r_sh] = canoncorr(Xw, Yw(randperm(n_w_local), :));
                cc_sh_pool(end+1,1) = r_sh(1); %#ok<AGROW>
            catch
            end
        end
    end
    if ~isempty(cc_pool)
        sweep_real_mu(iw)  = mean(cc_pool, 'omitnan');
        sweep_real_sem(iw) = std(cc_pool, 'omitnan') / sqrt(numel(cc_pool));
    end
    if ~isempty(cc_sh_pool)
        sweep_shuf_mu(iw)  = mean(cc_sh_pool, 'omitnan');
        sweep_shuf_sem(iw) = std(cc_sh_pool, 'omitnan') / sqrt(numel(cc_sh_pool));
    end
    fprintf('  half_w=%2d (%d bins, %d ms): real=%.3f sh=%.3f n=%d\n', ...
        h, n_w_local, n_w_local*target_bin_ms, sweep_real_mu(iw), sweep_shuf_mu(iw), numel(cc_pool));
end

figure('Name','Window-size sweep','Position',[550 550 700 400]);
hold on;
errorbar(2*sweep_half_w+1, sweep_real_mu, sweep_real_sem, '-o', ...
    'Color', [0 0.6 0], 'LineWidth', 2, 'DisplayName','Real');
errorbar(2*sweep_half_w+1, sweep_shuf_mu, sweep_shuf_sem, '--s', ...
    'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName','Shuffle');
xlabel(sprintf('Window size (bins; bin = %d ms)', target_bin_ms));
ylabel('Mean CC1 across windows');
title(sprintf('Window-size sensitivity (%s vs %s, k=%d)', area1, area2, k5_X));
legend('Location','best'); grid on;
ylim([0 1]);
saveas(gcf, fullfile(figures_dir, '09_window_sweep.png'));

% Save the underlying numerical data alongside the figures (CLAUDE.md viz rule).
save(fullfile(figures_dir, 'sanity_data.mat'), ...
     'trials_all', 'scatter_cc1', 'scatter_cc1_sh', 'scatter_vel', 'scatter_lick', ...
     'all_cc_temp_sp', 'all_cc_temp_sp_sh', 'all_ifi_temp_sp', ...
     'all_cc_pure_sp', 'all_ifi_pure_sp', ...
     'heldout_cc1', 'heldout_cc1_sh', 'trainset_cc1', ...
     'comp_t_cc', 'comp_p_cc', 'comp_t_ifi', 'comp_p_ifi', ...
     'block_lengths_bins', 'block_lengths_ms', ...
     'sweep_half_w', 'sweep_real_mu', 'sweep_real_sem', 'sweep_shuf_mu', 'sweep_shuf_sem', ...
     'area1', 'area2', 'target_bin_ms', 'half_window_bins', ...
     '-v7.3');

fprintf('All sanity checks complete. Plots saved to %s\n', figures_dir);

%% Local Helper Functions
function [X_pca, coeff, cum_var, u_mask] = get_pca_data(area, regions, spikes, valid_mask)
    u_mask = strcmp(regions, area);
    if sum(u_mask) < 3
        error('Not enough units in %s', area);
    end

    X = spikes(:, u_mask);
    X_valid = X(valid_mask, :);

    % Z-score against valid-window stats. After subtracting `mu` and
    % dividing by `sd`, X_z is already mean-zero in the valid window, so
    % we don't subtract the mean a second time before projecting.
    mu = mean(X_valid, 1, 'omitnan');
    sd = std(X_valid, 0, 1, 'omitnan'); sd(sd==0) = 1;
    T = size(X, 1);
    X_z = (X - repmat(mu, T, 1)) ./ repmat(sd, T, 1);

    [coeff, ~, ~, ~, explained] = pca(X_valid);
    cum_var = cumsum(explained);

    X_pca = X_z * coeff;
end

function [X_pca_sp, coeff, cum_var, u_mask] = get_spatial_pca_data(area, regions, sp_data, k)
    u_mask = strcmp(regions, area);
    if sum(u_mask) < 3
        error('Not enough units in %s', area);
    end
    X = sp_data(u_mask, :, :); % [units x bins x trials]
    [nU, nB, nT] = size(X);
    X_flat = reshape(X, nU, [])'; % [bins*trials x units]
    valid_idx = ~any(isnan(X_flat), 2);
    
    mu = mean(X_flat(valid_idx, :), 1, 'omitnan');
    sd = std(X_flat(valid_idx, :), 0, 1, 'omitnan'); sd(sd==0) = 1;
    
    % Handle potential div by 0
    X_z = nan(size(X_flat));
    X_z(valid_idx, :) = (X_flat(valid_idx, :) - mu) ./ sd;
    
    [coeff, ~, ~, ~, explained] = pca(X_z(valid_idx, :));
    cum_var = cumsum(explained);
    
    k_use = min(k, size(coeff, 2));
    X_pca_flat = nan(size(X_z, 1), k_use);
    X_pca_flat(valid_idx, :) = X_z(valid_idx, :) * coeff(:, 1:k_use);
    
    X_pca_sp = reshape(X_pca_flat', k_use, nB, nT); % [k x bins x trials]
end

function robust_fit_line(x, y)
% Drop-in replacement for `lsline` that filters NaN/Inf before fitting.
% Plots the OLS regression line over the visible x-range.
    x = x(:); y = y(:);
    keep = isfinite(x) & isfinite(y);
    if sum(keep) < 5, return; end
    p = polyfit(x(keep), y(keep), 1);
    xs = linspace(min(x(keep)), max(x(keep)), 100);
    plot(xs, polyval(p, xs), 'r-', 'LineWidth', 1.2, 'HandleVisibility','off');
end

function [pass, msg] = synthetic_cca_check(half_w, k_x, k_y)
% synthetic_cca_check: TDD-style smoke test for the sliding-window CCA path.
%
% Generates two regions (X, Y) over T=2000 samples:
%   coupled case: shared latent z(t) (3 dims, smoothed), X = z*Wx + noise,
%                 Y = z*Wy + noise -> expected CC1 ~ 0.6-0.8.
%   null case:    independent Gaussians in X and Y -> expected CC1 ~ shuffle.
%
% Then runs a sliding-window CCA on each and checks:
%   (1) coupled mean CC1 > coupled shuffle CC1 + 0.05
%   (2) null    mean CC1 < null shuffle CC1 + 0.05  (no spurious signal)
%
% Returns pass=true if both hold, false otherwise. msg is a one-line summary.
    n_w = 2*half_w + 1;
    T = 2000;
    nX = 8; nY = 10;
    rng(123);  % deterministic for the test

    % Smoothed shared latent (3 dims).
    Z = filter(ones(20,1)/20, 1, randn(T, 3));
    Wx = randn(3, nX);
    Wy = randn(3, nY);
    noise_amp = 1.5;

    Xc = Z * Wx + noise_amp * randn(T, nX);
    Yc = Z * Wy + noise_amp * randn(T, nY);

    Xn = randn(T, nX);
    Yn = randn(T, nY);

    function vals = run_sliding(X, Y)
        % PCA-reduce each side to k_x / k_y components and slide.
        [~, sX] = pca(X); [~, sY] = pca(Y);
        sX = sX(:, 1:min(k_x, size(sX,2)));
        sY = sY(:, 1:min(k_y, size(sY,2)));
        cc = []; cc_sh = [];
        for ic = (half_w+1):(size(sX,1) - half_w)
            iw = (ic-half_w):(ic+half_w);
            Xw = sX(iw,:); Yw = sY(iw,:);
            if rank(Xw - mean(Xw,1)) < size(Xw,2) || rank(Yw - mean(Yw,1)) < size(Yw,2), continue; end
            try, [~,~,r] = canoncorr(Xw, Yw); cc(end+1,1) = r(1); catch, end %#ok<AGROW>
            try, [~,~,r] = canoncorr(Xw, Yw(randperm(n_w),:)); cc_sh(end+1,1) = r(1); catch, end %#ok<AGROW>
        end
        vals = struct('mu_real', mean(cc,'omitnan'), 'mu_sh', mean(cc_sh,'omitnan'), ...
                      'n', numel(cc));
    end

    rc = run_sliding(Xc, Yc);
    rn = run_sliding(Xn, Yn);

    sep_c = rc.mu_real - rc.mu_sh;
    sep_n = rn.mu_real - rn.mu_sh;
    pass = (sep_c > 0.05) && (abs(sep_n) < 0.05);

    msg = sprintf(['coupled real=%.3f shuf=%.3f sep=%+.3f  |  ' ...
                   'null real=%.3f shuf=%.3f sep=%+.3f  | n=%d windows | pass=%d'], ...
                  rc.mu_real, rc.mu_sh, sep_c, ...
                  rn.mu_real, rn.mu_sh, sep_n, rc.n, pass);
end

function [cc_test, cc_train, cc_test_sh] = compute_heldout_cc1(X_tr, Y_tr, valid_tr, half_w, k_x, k_y, train_frac, n_sh)
% compute_heldout_cc1: honest CC1 estimate via held-out canonical-variate correlation.
%
% For one trial:
%   1. Find all valid sliding-window centers and pool their bins into a
%      single (n_samples x k) matrix per region.
%   2. Randomly split rows 70/30 (or train_frac).
%   3. Fit canoncorr(X_train, Y_train) -> A, B.
%   4. Project X_test, Y_test with A(:,1), B(:,1) and report
%      corr(X_test*A(:,1), Y_test*B(:,1)) -- the held-out canonical
%      correlation. This is the bias-corrected version of CC1.
%   5. Repeat n_sh times with Y_train rows randomly permuted, average,
%      report as the shuffle null.
%
% Returns NaN if too few valid samples or if rank deficient.
    cc_test    = nan;
    cc_train   = nan;
    cc_test_sh = nan;

    n = size(X_tr, 1);
    if n < (2*half_w + 5), return; end

    % Pool every valid bin (fully-contiguous valid window centered there).
    keep = false(n, 1);
    for ic = (half_w+1):(n - half_w)
        if all(valid_tr((ic-half_w):(ic+half_w))), keep(ic) = true; end
    end

    Xv = X_tr(keep, 1:k_x);
    Yv = Y_tr(keep, 1:k_y);
    n_valid = size(Xv, 1);
    if n_valid < (k_x + k_y + 10), return; end

    n_train = max(k_x + k_y + 5, round(train_frac * n_valid));
    if n_train >= n_valid, n_train = n_valid - max(5, ceil(0.1*n_valid)); end

    perm = randperm(n_valid);
    idx_train = perm(1:n_train);
    idx_test  = perm(n_train+1:end);

    Xt = Xv(idx_train, :); Yt = Yv(idx_train, :);
    Xs = Xv(idx_test,  :); Ys = Yv(idx_test,  :);

    if rank(Xt - mean(Xt)) < k_x || rank(Yt - mean(Yt)) < k_y, return; end

    try
        [A, B, r] = canoncorr(Xt, Yt);
    catch
        return;
    end
    cc_train = r(1);

    % Use the SAME centering as canoncorr (column means of training set)
    u_test = (Xs - mean(Xt, 1)) * A(:,1);
    v_test = (Ys - mean(Yt, 1)) * B(:,1);
    if numel(u_test) > 5
        cc_test = corr(u_test, v_test);
    end

    % Shuffle null: permute Y rows in the training set, refit, project test.
    r_sh = nan(n_sh, 1);
    for is = 1:n_sh
        Yt_p = Yt(randperm(n_train), :);
        if rank(Yt_p - mean(Yt_p)) < k_y, continue; end
        try
            [As, Bs, ~] = canoncorr(Xt, Yt_p);
        catch
            continue;
        end
        us = (Xs - mean(Xt, 1))   * As(:,1);
        vs = (Ys - mean(Yt_p, 1)) * Bs(:,1);
        if numel(us) > 5, r_sh(is) = corr(us, vs); end
    end
    cc_test_sh = mean(r_sh, 'omitnan');
end
