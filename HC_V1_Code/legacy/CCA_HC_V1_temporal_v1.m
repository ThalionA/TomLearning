%% CCA_HC_V1_temporal_v1.m
% Temporal CCA Analysis with Masked Cross-Correlation
%
% DESCRIPTION:
%   Mirrors the spatial v2 pipeline (CCA_HC_V1_spatial_v2.m) but operates on
%   binned_spikes (time domain) rather than analysis_spatial.firing.cued.freq_z
%   (space domain). For every region pair and every epoch (early / pre-learn /
%   post-learn) we:
%     1) Re-bin binned_spikes from 1 ms (native) to 50 ms (target).
%     2) Build idx_valid = tunnel_cued AND velocity_gf >= 2 cm/s
%        AND trial_binned_cued in epoch trials.
%     3) PCA-reduce each region (>=90% variance, capped to >= num_ccs_analyze).
%     4) Fit canoncorr() on the validity-restricted samples.
%     5) Project the full timecourse onto the canonical components.
%     6) Compute Tom's masked cross-correlation across lags = -5:5 (i.e. +-250 ms),
%        using the centered-mask trick so no lagged sample crosses a discontinuity.
%     7) Shuffle control: permute valid-sample order of region B before xcorr.
%
% Lag interpretation (matches Tom's note):
%   negative lag  ==>  region B (=area2) leads region A (=area1)
%   positive lag  ==>  region A (=area1) leads region B (=area2)

%% 1. SETUP & PARAMETERS
clear; clc; close all;

% --- Paths (same layout as spatial v2) ---
base_dir       = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir       = fullfile(base_dir, data_subfolder);
file_pattern   = 'TF*_export.mat';
learning_file  = 'animal_behaviour.mat';

current_date = datestr(now, 'yyyy_mm_dd');
save_path    = fullfile(data_dir, sprintf('Temporal_CCA_Results_%s.mat', current_date));

% --- Temporal binning ---
target_bin_ms          = 50;        % rebin binned_spikes to 50 ms
lags                   = -5:5;      % +/- 5 bins  =>  +/- 250 ms
n_lags                 = numel(lags);

% --- PCA / CCA params (match spatial v2 conventions) ---
pca_variance_threshold = 90;
num_ccs_analyze        = 3;
n_shuffles             = 50;
min_units_per_region   = 5;
min_speed_cms          = 2;         % matches params_main.behaviour.minspeed
n_trials_epoch         = 10;        % early/pre/post epoch length (trials)

% --- Region pairs (same as spatial v2) ---
area_pairs_to_analyze = {'CA1','V1';  'CA1','DG'; 'CA1','CA3'; 'CA1','RSC'; ...
                         'CA1','SUB'; 'V1','RSC'; 'RSC','SUB'; 'CA3','DG'};
n_pairs = size(area_pairs_to_analyze, 1);

% Set to true to ignore saved files and re-run the analysis.
force_reprocess = true;

% Open a parpool if one is not already running. (The xcorr shuffle loop is parfor.)
if isempty(gcp('nocreate'))
    try, parpool; catch, warning('Could not start parpool; continuing serial.'); end
end

%% 2. LOAD & YOKE LEARNING POINTS  (same logic as spatial v2)
file_list = dir(fullfile(data_dir, file_pattern));
n_animals = length(file_list);

lp_path = fullfile(data_dir, learning_file);
if exist(lp_path, 'file')
    dat_lp = load(lp_path);
    learning_points = dat_lp.period_experienced(:, 1);
    if isfield(dat_lp, 'animal_id')
        [~, sorting_idx] = sort(dat_lp.animal_id);
        learning_points  = learning_points(sorting_idx);
    end
else
    warning('Learning file not found. Alignment will be skipped.');
    learning_points = nan(n_animals, 1);
end
if length(learning_points) < n_animals
    learning_points(end+1:n_animals) = nan;
end

is_learner   = ~isnan(learning_points);
mean_lp      = round(mean(learning_points(is_learner), 'omitnan'));
analysis_lp  = learning_points;
analysis_lp(~is_learner) = mean_lp;

fprintf('Identified %d Learners and %d Non-Learners. Yoked LP = %d\n', ...
    sum(is_learner), sum(~is_learner), mean_lp);

%% 3. INITIALISE RESULTS STRUCTURE
group_results = struct( ...
    'pair_name',       cell(n_pairs, 1), ...
    'xc_early',        cell(n_pairs, 1), ...   % {animal} -> [num_ccs x n_lags]
    'xc_pre',          cell(n_pairs, 1), ...
    'xc_post',         cell(n_pairs, 1), ...
    'xc_early_shuff',  cell(n_pairs, 1), ...
    'xc_pre_shuff',    cell(n_pairs, 1), ...
    'xc_post_shuff',   cell(n_pairs, 1), ...
    'cca_r_early',     cell(n_pairs, 1), ...   % {animal} -> [num_ccs x 1]
    'cca_r_pre',       cell(n_pairs, 1), ...
    'cca_r_post',      cell(n_pairs, 1), ...
    'n_samples_early', cell(n_pairs, 1), ...
    'n_samples_pre',   cell(n_pairs, 1), ...
    'n_samples_post',  cell(n_pairs, 1));
for ipair = 1:n_pairs
    group_results(ipair).pair_name = sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});
end

%% 4. MAIN ANALYSIS LOOP
current_config = struct( ...
    'target_bin_ms',          target_bin_ms, ...
    'lags',                   lags, ...
    'pca_variance_threshold', pca_variance_threshold, ...
    'num_ccs_analyze',        num_ccs_analyze, ...
    'n_shuffles',             n_shuffles, ...
    'min_units_per_region',   min_units_per_region, ...
    'min_speed_cms',          min_speed_cms, ...
    'n_trials_epoch',         n_trials_epoch);

existing_files = dir(fullfile(data_dir, 'Temporal_CCA_Results_*.mat'));

if ~isempty(existing_files) && ~force_reprocess
    [~, latest_idx] = max([existing_files.datenum]);
    load_target = fullfile(data_dir, existing_files(latest_idx).name);
    load(load_target, 'group_results', 'is_learner', 'analysis_lp', 'saved_config');
    fprintf('Loaded existing temporal results: %s\n', existing_files(latest_idx).name);
else
    if force_reprocess
        fprintf('Force reprocess enabled. Starting temporal CCA from scratch...\n');
    end

    for ianimal = 1:n_animals
        filename = file_list(ianimal).name;
        fullpath = fullfile(data_dir, filename);
        fprintf('\nProcessing Animal %d/%d: %s\n', ianimal, n_animals, filename);

        try
            D = load(fullpath);
            req_fields = {'units','binned_spikes','data_behaviour','analysis_behaviour','params_main'};
            if ~all(isfield(D, req_fields))
                fprintf('  Missing required fields. Skipping.\n'); continue;
            end

            units = D.units;

            % --- FS unit exclusion (same idiom as spatial v2) ---
            keep_mask = true(length(units.unit_id), 1);
            if isfield(units, 'idx_fs')
                target_fs_areas = {'V1','RSC','CA1','CA3'};
                is_fs = logical(units.idx_fs);
                for r = 1:length(target_fs_areas)
                    fs_in_area = strcmp(units.region, target_fs_areas{r}) & is_fs;
                    keep_mask(fs_in_area) = false;
                end
            end
            % Zero-out the rows of units.idx that correspond to excluded units
            % (matches spatial v2's behaviour exactly).
            units.idx(~keep_mask, :) = 0;

            % --- Determine native bin size (s) and rebin factor ---
            native_bin_s = double(D.params_main.bin_size);
            target_bin_s = target_bin_ms / 1000;
            rebin_factor = round(target_bin_s / native_bin_s);
            if rebin_factor < 1, rebin_factor = 1; end
            fprintf('  native bin = %.4f s, rebin factor = %d -> target bin = %.3f s\n', ...
                native_bin_s, rebin_factor, rebin_factor*native_bin_s);

            % --- Re-bin binned_spikes from [units x T] to [T_re x units] ---
            spikes = double(D.binned_spikes);          % MATLAB layout: [nUnits x T]
            if size(spikes, 1) == length(units.unit_id) && size(spikes, 2) ~= length(units.unit_id)
                % already [units x time]
            else
                spikes = spikes';                      % fall-back if layout is [T x units]
            end
            nUnits = size(spikes, 1);
            T      = size(spikes, 2);
            T_use  = floor(T / rebin_factor) * rebin_factor;
            spikes = spikes(:, 1:T_use);
            % Sum spikes across each rebin block: reshape to [nUnits x rebin_factor x T_re]
            T_re = T_use / rebin_factor;
            spikes_re = squeeze(sum(reshape(spikes, nUnits, rebin_factor, T_re), 2));   % [nUnits x T_re]
            spikes_re = spikes_re';                                                     % [T_re x nUnits]

            % --- Re-bin behavioural masks/vectors (truncate then reshape) ---
            mask_cued  = logical(D.analysis_behaviour.masks.tunnel_cued(:));
            vel_gf     = double(D.data_behaviour.velocity_binned_gf(:));
            tr_cued    = double(D.data_behaviour.trial_binned_cued(:));

            mask_cued = mask_cued(1:T_use);
            vel_gf    = vel_gf(1:T_use);
            tr_cued   = tr_cued(1:T_use);

            mc_re = all(reshape(mask_cued, rebin_factor, T_re), 1)';                     % [T_re x 1]
            vel_re = mean(reshape(vel_gf,  rebin_factor, T_re), 1)';                     % mean speed in 50 ms bin
            ms_re  = vel_re >= min_speed_cms;
            % For the trial label of a 50 ms bin, take the mode across the block.
            % (In practice every native ms inside the same bin shares the same trial,
            % except at the trial boundary. mode picks the dominant trial.)
            tc_re  = mode(reshape(tr_cued, rebin_factor, T_re), 1)';                      % [T_re x 1]

            base_valid = mc_re & ms_re;            % "in cued tunnel + moving"

            if sum(base_valid) < 100
                fprintf('  Too few valid samples (%d). Skipping.\n', sum(base_valid));
                continue;
            end

            % --- Z-score each unit using valid samples only ---
            mu  = mean(spikes_re(base_valid, :), 1, 'omitnan');
            sd  = std( spikes_re(base_valid, :), 0, 1, 'omitnan');
            sd(sd == 0) = 1;
            spikes_z = (spikes_re - mu) ./ sd;

            % --- Build per-region PCA-reduced timecourses ---
            % Mirror spatial v2: pick units by strcmp on units.region (per-unit area name)
            % rather than relying on the row order of units.idx.
            animal_areas = unique(units.regions_label);
            AreaActivity = struct();
            for ia = 1:length(animal_areas)
                area = char(animal_areas{ia});
                if isempty(area), continue; end

                % Per-unit logical index for this area, AND with the FS keep_mask.
                u_logical = strcmp(units.region, area) & keep_mask(:);
                u_logical = u_logical(:);
                if sum(u_logical) < min_units_per_region, continue; end

                X = spikes_z(:, u_logical);                  % [T_re x nUnitsRegion]
                if size(X, 2) < 2, continue; end
                X_valid = X(base_valid, :);
                if size(X_valid, 1) < 50, continue; end
                if any(~isfinite(X_valid(:))), continue; end

                [coeff, ~, ~, ~, explained] = pca(X_valid);
                cum_var = cumsum(explained);
                k = find(cum_var >= pca_variance_threshold, 1);
                if isempty(k), k = size(coeff, 2); end
                k = max(num_ccs_analyze, k);
                k = min(k, size(coeff, 2));

                % Project ALL timepoints (centered with the SAME mean used for PCA).
                X_centered = X - mean(X_valid, 1);
                scores_all = X_centered * coeff(:, 1:k);     % [T_re x k]

                fld = matlab.lang.makeValidName(area);
                AreaActivity.(fld).data = scores_all;
                AreaActivity.(fld).k    = k;
                AreaActivity.(fld).orig_name = area;
            end

            % --- Define epochs by trial range (same as spatial v2) ---
            lp_animal  = analysis_lp(ianimal);
            num_trials = max(tc_re);
            if isnan(lp_animal) || lp_animal <= n_trials_epoch || (lp_animal + n_trials_epoch - 1) > num_trials
                fprintf('  LP=%g unusable for animal (num_trials=%d). Skipping.\n', lp_animal, num_trials);
                continue;
            end
            tr_early = 1:n_trials_epoch;
            tr_pre   = (lp_animal - n_trials_epoch):(lp_animal - 1);
            tr_post  = lp_animal:(lp_animal + n_trials_epoch - 1);

            % --- Pair-wise temporal CCA + masked xcorr ---
            for ipair = 1:n_pairs
                a1 = matlab.lang.makeValidName(area_pairs_to_analyze{ipair, 1});
                a2 = matlab.lang.makeValidName(area_pairs_to_analyze{ipair, 2});
                if ~isfield(AreaActivity, a1) || ~isfield(AreaActivity, a2), continue; end

                X = AreaActivity.(a1).data;     % [T_re x kA]
                Y = AreaActivity.(a2).data;     % [T_re x kB]

                xc = struct();
                xc.early       = nan(num_ccs_analyze, n_lags);
                xc.pre         = nan(num_ccs_analyze, n_lags);
                xc.post        = nan(num_ccs_analyze, n_lags);
                xc.early_shuff = nan(num_ccs_analyze, n_lags);
                xc.pre_shuff   = nan(num_ccs_analyze, n_lags);
                xc.post_shuff  = nan(num_ccs_analyze, n_lags);
                rs = struct('early',nan(num_ccs_analyze,1),'pre',nan(num_ccs_analyze,1),'post',nan(num_ccs_analyze,1));
                ns = struct('early',0,'pre',0,'post',0);

                for iep = 1:3
                    switch iep
                        case 1, tr_set = tr_early; ep_name = 'early';
                        case 2, tr_set = tr_pre;   ep_name = 'pre';
                        case 3, tr_set = tr_post;  ep_name = 'post';
                    end

                    valid = base_valid & ismember(tc_re, tr_set);
                    ns.(ep_name) = sum(valid);
                    if ns.(ep_name) < (max(lags) - min(lags) + 50), continue; end

                    Xv = X(valid, :);
                    Yv = Y(valid, :);
                    if any(~isfinite(Xv(:))) || any(~isfinite(Yv(:))), continue; end

                    % --- Fit CCA on the valid samples ---
                    try
                        [A, B, r] = canoncorr(Xv, Yv);
                    catch
                        continue;
                    end
                    n_use = min([num_ccs_analyze, size(A, 2), size(B, 2)]);
                    rs.(ep_name)(1:n_use) = r(1:n_use);

                    % Project full timecourse onto canonical components.
                    Xc = X - mean(Xv, 1);
                    Yc = Y - mean(Yv, 1);
                    Ax = Xc * A(:, 1:n_use);    % [T_re x n_use]   (region A canonical)
                    By = Yc * B(:, 1:n_use);    % [T_re x n_use]   (region B canonical)

                    % --- Build centered validity mask (Tom's recipe) ---
                    idx_centered = valid;
                    for L = lags
                        idx_centered = idx_centered & circshift(valid, L);
                    end
                    n_centered = sum(idx_centered);
                    if n_centered < 50, continue; end

                    % --- Real masked cross-correlation across CCs ---
                    xc_real = nan(n_use, n_lags);
                    for ic = 1:n_use
                        a_ref = Ax(idx_centered, ic);
                        for il = 1:n_lags
                            L = lags(il);
                            idx_lag = circshift(idx_centered, L);
                            b_lag   = By(idx_lag, ic);
                            xc_real(ic, il) = corr(a_ref, b_lag, 'rows', 'complete');
                        end
                    end

                    % --- Shuffle controls: scramble valid-sample order of region B ---
                    valid_idx_list = find(valid);
                    n_valid = numel(valid_idx_list);
                    xc_sh = nan(n_use, n_lags, n_shuffles);
                    parfor ish = 1:n_shuffles
                        offset  = randi(max(n_valid-1, 1));
                        permuted = mod((1:n_valid)-1 + offset, n_valid) + 1;
                        By_s = By;
                        By_s(valid_idx_list, :) = By(valid_idx_list(permuted), :);
                        xc_sh_iter = nan(n_use, n_lags);
                        for ic = 1:n_use
                            a_ref = Ax(idx_centered, ic);
                            for il = 1:n_lags
                                L = lags(il);
                                idx_lag = circshift(idx_centered, L);
                                xc_sh_iter(ic, il) = corr(a_ref, By_s(idx_lag, ic), 'rows', 'complete');
                            end
                        end
                        xc_sh(:, :, ish) = xc_sh_iter;
                    end
                    xc_sh_mean = mean(xc_sh, 3, 'omitnan');

                    % Pad to num_ccs_analyze
                    pad_real = nan(num_ccs_analyze, n_lags); pad_real(1:n_use, :) = xc_real;
                    pad_sh   = nan(num_ccs_analyze, n_lags); pad_sh(1:n_use, :)   = xc_sh_mean;

                    xc.(ep_name)              = pad_real;
                    xc.([ep_name '_shuff'])   = pad_sh;
                end

                % --- Store ---
                group_results(ipair).xc_early{ianimal}        = xc.early;
                group_results(ipair).xc_pre{ianimal}          = xc.pre;
                group_results(ipair).xc_post{ianimal}         = xc.post;
                group_results(ipair).xc_early_shuff{ianimal}  = xc.early_shuff;
                group_results(ipair).xc_pre_shuff{ianimal}    = xc.pre_shuff;
                group_results(ipair).xc_post_shuff{ianimal}   = xc.post_shuff;
                group_results(ipair).cca_r_early{ianimal}     = rs.early;
                group_results(ipair).cca_r_pre{ianimal}       = rs.pre;
                group_results(ipair).cca_r_post{ianimal}      = rs.post;
                group_results(ipair).n_samples_early{ianimal} = ns.early;
                group_results(ipair).n_samples_pre{ianimal}   = ns.pre;
                group_results(ipair).n_samples_post{ianimal}  = ns.post;
            end
        catch ME
            fprintf('  Error processing %s: %s\n', filename, ME.message);
        end
    end

    saved_config = current_config;
    save(save_path, 'group_results', 'is_learner', 'analysis_lp', 'saved_config', '-v7.3');
    fprintf('Saved temporal results to %s\n', save_path);
end

%% 5. PADDING (so cell arrays line up with is_learner length)
n_total_animals = length(is_learner);
for ipair = 1:n_pairs
    fields = fieldnames(group_results);
    for f = 1:length(fields)
        if iscell(group_results(ipair).(fields{f})) && length(group_results(ipair).(fields{f})) < n_total_animals
            group_results(ipair).(fields{f}){n_total_animals} = [];
        end
    end
end

lag_ms = lags * target_bin_ms;

%% 6. PLOTS
fprintf('\nGenerating temporal CCA plots...\n');

% --- A. xcorr curves per pair, Learners overlay (early/pre/post) ---
figure('Name','Temporal Cross-Correlation (CC1) - Learners','Color','w','Position',[100 100 1400 800]);
tiledlayout('flow','TileSpacing','compact');
col_e = [0.6 0.7 0.95];      % early - light
col_p = [0.35 0.5 0.85];     % pre   - mid
col_x = [0.10 0.25 0.6 ];    % post  - dark
for ipair = 1:n_pairs
    nexttile; hold on;
    [mu_e, se_e] = collect_xc(group_results(ipair).xc_early, is_learner, 1);
    [mu_p, se_p] = collect_xc(group_results(ipair).xc_pre,   is_learner, 1);
    [mu_x, se_x] = collect_xc(group_results(ipair).xc_post,  is_learner, 1);
    [mu_es,~]   = collect_xc(group_results(ipair).xc_early_shuff, is_learner, 1);
    [mu_xs,~]   = collect_xc(group_results(ipair).xc_post_shuff,  is_learner, 1);

    if ~isempty(mu_e), errorbar_shaded(lag_ms, mu_e, se_e, col_e); end
    if ~isempty(mu_p), errorbar_shaded(lag_ms, mu_p, se_p, col_p); end
    if ~isempty(mu_x), errorbar_shaded(lag_ms, mu_x, se_x, col_x); end
    if ~isempty(mu_es), plot(lag_ms, mu_es, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); end
    if ~isempty(mu_xs), plot(lag_ms, mu_xs, ':',  'Color', [0.3 0.3 0.3], 'LineWidth', 1); end

    xline(0, 'k:'); yline(0, 'k:');
    xlabel(sprintf('Lag (ms)   %s leads <- 0 -> %s leads', area_pairs_to_analyze{ipair, 2}, area_pairs_to_analyze{ipair, 1}));
    ylabel('xcorr (CC1)');
    title(group_results(ipair).pair_name);
end
legend({'Early','Pre','Post','Early shuff','Post shuff'}, 'Location','best');
sgtitle('Learners: Temporal xcorr (centered-mask)');
try, saveas(gcf, fullfile(data_dir, 'Temporal_xcorr_Learners.png')); end

% --- A2. xcorr curves per pair, Non-Learners overlay (early/pre/post) ---
if any(~is_learner)
    figure('Name','Temporal Cross-Correlation (CC1) - Non-Learners','Color','w','Position',[100 100 1400 800]);
    tiledlayout('flow','TileSpacing','compact');
    col_eN = [0.95 0.7 0.55];
    col_pN = [0.85 0.5 0.3];
    col_xN = [0.6  0.2 0.1];
    for ipair = 1:n_pairs
        nexttile; hold on;
        [mu_e, se_e] = collect_xc(group_results(ipair).xc_early, ~is_learner, 1);
        [mu_p, se_p] = collect_xc(group_results(ipair).xc_pre,   ~is_learner, 1);
        [mu_x, se_x] = collect_xc(group_results(ipair).xc_post,  ~is_learner, 1);

        if ~isempty(mu_e), errorbar_shaded(lag_ms, mu_e, se_e, col_eN); end
        if ~isempty(mu_p), errorbar_shaded(lag_ms, mu_p, se_p, col_pN); end
        if ~isempty(mu_x), errorbar_shaded(lag_ms, mu_x, se_x, col_xN); end

        xline(0, 'k:'); yline(0, 'k:');
        xlabel('Lag (ms)'); ylabel('xcorr (CC1)');
        title(group_results(ipair).pair_name);
    end
    sgtitle('Non-Learners: Temporal xcorr');
    try, saveas(gcf, fullfile(data_dir, 'Temporal_xcorr_NonLearners.png')); end
end

% --- B. Peak |xcorr| per epoch (Learners vs Non-Learners) ---
figure('Name','Peak |xcorr| per Epoch (CC1)','Color','w','Position',[100 100 1400 800]);
tiledlayout('flow','TileSpacing','compact');
for ipair = 1:n_pairs
    nexttile;
    pk_e = peak_per_animal(group_results(ipair).xc_early, 1);
    pk_p = peak_per_animal(group_results(ipair).xc_pre,   1);
    pk_x = peak_per_animal(group_results(ipair).xc_post,  1);
    plot_peaks_bars(pk_e, pk_p, pk_x, is_learner, group_results(ipair).pair_name, 'Peak |xcorr|');
end
sgtitle('Peak |xcorr| (CC1) per epoch');
try, saveas(gcf, fullfile(data_dir, 'Temporal_peak_xcorr_per_epoch.png')); end

% --- C. Peak lag (ms) per epoch ---
figure('Name','Peak lag per epoch (CC1)','Color','w','Position',[100 100 1400 800]);
tiledlayout('flow','TileSpacing','compact');
for ipair = 1:n_pairs
    nexttile;
    plg_e = peaklag_per_animal(group_results(ipair).xc_early, 1, lag_ms);
    plg_p = peaklag_per_animal(group_results(ipair).xc_pre,   1, lag_ms);
    plg_x = peaklag_per_animal(group_results(ipair).xc_post,  1, lag_ms);
    plot_peaks_bars(plg_e, plg_p, plg_x, is_learner, group_results(ipair).pair_name, 'Peak lag (ms)');
    yline(0, 'k--');
end
sgtitle('Peak lag (CC1) per epoch  (neg = area2 leads, pos = area1 leads)');
try, saveas(gcf, fullfile(data_dir, 'Temporal_peak_lag_per_epoch.png')); end

% --- D. CCA r1 (zero-lag canonical correlation) per epoch ---
figure('Name','CCA r_1 per Epoch','Color','w','Position',[100 100 1400 800]);
tiledlayout('flow','TileSpacing','compact');
for ipair = 1:n_pairs
    nexttile;
    r_e = scalar_per_animal(group_results(ipair).cca_r_early, 1);
    r_p = scalar_per_animal(group_results(ipair).cca_r_pre,   1);
    r_x = scalar_per_animal(group_results(ipair).cca_r_post,  1);
    plot_peaks_bars(r_e, r_p, r_x, is_learner, group_results(ipair).pair_name, 'CCA r_1');
end
sgtitle('Canonical correlation r_1 per epoch');
try, saveas(gcf, fullfile(data_dir, 'Temporal_CCA_r1_per_epoch.png')); end

fprintf('Done.\n');

%% =========================================================================
%% LOCAL HELPERS
%% =========================================================================
function [mu, se] = collect_xc(C, mask, cc)
    mat = [];
    for i = 1:length(C)
        if i > length(mask) || ~mask(i), continue; end
        if isempty(C{i}), continue; end
        if size(C{i}, 1) < cc, continue; end
        row = C{i}(cc, :);
        if all(isnan(row)), continue; end
        mat(end+1, :) = row; %#ok<AGROW>
    end
    if isempty(mat), mu = []; se = []; return; end
    n  = sum(~isnan(mat), 1);
    mu = mean(mat, 1, 'omitnan');
    se = std(mat, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
end

function v = peak_per_animal(C, cc)
    v = nan(length(C), 1);
    for i = 1:length(C)
        if isempty(C{i}) || size(C{i}, 1) < cc, continue; end
        row = C{i}(cc, :);
        if all(isnan(row)), continue; end
        v(i) = max(abs(row), [], 'omitnan');
    end
end

function v = peaklag_per_animal(C, cc, lag_ms)
    v = nan(length(C), 1);
    for i = 1:length(C)
        if isempty(C{i}) || size(C{i}, 1) < cc, continue; end
        row = C{i}(cc, :);
        if all(isnan(row)), continue; end
        [~, k] = max(abs(row), [], 'omitnan');
        v(i) = lag_ms(k);
    end
end

function v = scalar_per_animal(C, cc)
    v = nan(length(C), 1);
    for i = 1:length(C)
        if isempty(C{i}) || numel(C{i}) < cc, continue; end
        v(i) = C{i}(cc);
    end
end

function plot_peaks_bars(e_vals, p_vals, x_vals, is_learner, ttl, ylab)
    hold on;
    means_L = [mean(e_vals(is_learner), 'omitnan'), ...
               mean(p_vals(is_learner), 'omitnan'), ...
               mean(x_vals(is_learner), 'omitnan')];
    sem_L   = [sem_(e_vals(is_learner)), sem_(p_vals(is_learner)), sem_(x_vals(is_learner))];
    errorbar(1:3, means_L, sem_L, 'o-', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10], 'DisplayName', 'Learners');

    if any(~is_learner)
        means_NL = [mean(e_vals(~is_learner), 'omitnan'), ...
                    mean(p_vals(~is_learner), 'omitnan'), ...
                    mean(x_vals(~is_learner), 'omitnan')];
        sem_NL   = [sem_(e_vals(~is_learner)), sem_(p_vals(~is_learner)), sem_(x_vals(~is_learner))];
        errorbar(1:3, means_NL, sem_NL, 's-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Non-Learners');
    end
    xticks(1:3); xticklabels({'Early','Pre','Post'}); xlim([0.5 3.5]);
    ylabel(ylab); title(ttl);
    legend('Location','best');
end

function s = sem_(v)
    v = v(~isnan(v));
    if isempty(v), s = NaN; else, s = std(v) / sqrt(numel(v)); end
end

function errorbar_shaded(x, y, e, c)
    % Lightweight shaded errorbar (avoids dependency on shadedErrorBar.m).
    x = x(:)'; y = y(:)'; e = e(:)';
    fill([x, fliplr(x)], [y+e, fliplr(y-e)], c, 'EdgeColor','none','FaceAlpha',0.2);
    plot(x, y, 'Color', c, 'LineWidth', 1.6);
end
