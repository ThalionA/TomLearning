%% CCA_HC_V1_Temporal_Continuous_Windowed.m
% Continuous Temporal CCA with Sliding Trial Windows and Masked Lags
%
% DESCRIPTION:
% Operates on continuous binned time, applying a sliding window defined by 
% trial labels to dynamically calculate CC1 and the Information Flow Index (IFI).

%% 1. SETUP & PARAMETERS
clear; clc; close all;

base_dir       = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir       = fullfile(base_dir, data_subfolder);
file_pattern   = 'TF*_export.mat';
learning_file  = 'animal_behaviour.mat';

current_date = datestr(now, 'yyyy_mm_dd');
save_path    = fullfile(data_dir, sprintf('Temporal_Continuous_CCA_%s.mat', current_date));

% --- Temporal & Analysis Parameters ---
target_bin_ms          = 50;        
lags                   = -5:5;      % Lags to analyze (e.g., +/- 250ms)
pca_variance_threshold = 90;
n_trials_window        = -3:3;      % Sliding window for CCA fit
n_shuffles             = 50;
min_units_per_region   = 5;
min_speed_cms          = 2;         
min_samples_window     = 30;        % Minimum valid bins required in a window to fit CCA

area_pairs_to_analyze = {'CA1','V1';  'CA1','DG'; 'CA1','CA3'; 'CA1','RSC'; ...
                         'CA1','SUB'; 'V1','RSC'; 'RSC','SUB'; 'CA3','DG'};
n_pairs = size(area_pairs_to_analyze, 1);
force_reprocess = true;

if isempty(gcp('nocreate'))
    parpool;
end

%% 2. LOAD & YOKE LEARNING POINTS
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

%% 3. INITIALIZE RESULTS
group_results = struct( ...
    'pair_name', cell(n_pairs, 1), ...
    'trial_corr_early', cell(n_pairs, 1), 'trial_corr_pre', cell(n_pairs, 1), 'trial_corr_post', cell(n_pairs, 1), ...
    'trial_corr_early_shuff', cell(n_pairs, 1), 'trial_corr_pre_shuff', cell(n_pairs, 1), 'trial_corr_post_shuff', cell(n_pairs, 1), ...
    'trial_precession_early_idx', cell(n_pairs, 1), 'trial_precession_pre_idx', cell(n_pairs, 1), 'trial_precession_post_idx', cell(n_pairs, 1), ...
    'trial_precession_early_idx_shuff', cell(n_pairs, 1), 'trial_precession_pre_idx_shuff', cell(n_pairs, 1), 'trial_precession_post_idx_shuff', cell(n_pairs, 1));

for ipair = 1:n_pairs
    group_results(ipair).pair_name = sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});
end

%% 4. MAIN PIPELINE
for ianimal = 1:n_animals
    filename = file_list(ianimal).name;
    fullpath = fullfile(data_dir, filename);
    fprintf('\nProcessing Animal %d/%d: %s\n', ianimal, n_animals, filename);

    try
        D = load(fullpath);
        if ~isfield(D, 'units') || ~isfield(D, 'binned_spikes'), continue; end
        units = D.units;

        % 4.1 Filter FS units
        keep_mask = true(length(units.unit_id), 1);
        if isfield(units, 'idx_fs')
            target_fs_areas = {'V1','RSC','CA1','CA3'};
            is_fs = logical(units.idx_fs);
            for r = 1:length(target_fs_areas)
                fs_in_area = strcmp(units.region, target_fs_areas{r}) & is_fs;
                keep_mask(fs_in_area) = false;
            end
        end

        % 4.2 Re-binning to target_bin_ms
        native_bin_s = double(D.params_main.bin_size);
        rebin_factor = round((target_bin_ms / 1000) / native_bin_s);
        if rebin_factor < 1, rebin_factor = 1; end
        
        spikes = double(D.binned_spikes);
        if size(spikes, 1) ~= length(units.unit_id), spikes = spikes'; end
        
        nUnits = size(spikes, 1);
        T      = size(spikes, 2);
        T_use  = floor(T / rebin_factor) * rebin_factor;
        T_re   = T_use / rebin_factor;
        
        spikes_re = squeeze(sum(reshape(spikes(:, 1:T_use), nUnits, rebin_factor, T_re), 2))';
        
        mask_cued = logical(D.analysis_behaviour.masks.tunnel_cued(1:T_use));
        vel_gf    = double(D.data_behaviour.velocity_binned_gf(1:T_use));
        tr_cued   = double(D.data_behaviour.trial_binned_cued(1:T_use));

        mc_re  = all(reshape(mask_cued, rebin_factor, T_re), 1)';
        vel_re = mean(reshape(vel_gf,  rebin_factor, T_re), 1)';
        tc_re  = mode(reshape(tr_cued, rebin_factor, T_re), 1)';
        
        base_valid = mc_re & (vel_re >= min_speed_cms);

        % Z-score on valid points
        mu = mean(spikes_re(base_valid, :), 1, 'omitnan');
        sd = std(spikes_re(base_valid, :), 0, 1, 'omitnan'); sd(sd==0) = 1;
        spikes_z = (spikes_re - mu) ./ sd;

        % 4.3 PCA Generation (Continuous)
        animal_areas = unique(units.regions_label);
        AreaActivity = struct();
        for ia = 1:length(animal_areas)
            area = char(animal_areas{ia});
            if isempty(area), continue; end
            
            u_logical = strcmp(units.region, area) & keep_mask(:);
            if sum(u_logical) < min_units_per_region, continue; end
            
            X = spikes_z(:, u_logical);
            X_valid = X(base_valid, :);
            if size(X_valid, 1) < min_samples_window, continue; end
            
            [coeff, ~, ~, ~, explained] = pca(X_valid);
            cum_var = cumsum(explained);
            k = max(1, find(cum_var >= pca_variance_threshold, 1));
            
            % Center and project continuous timeline
            AreaActivity.(area).data = (X - mean(X_valid, 1)) * coeff(:, 1:k);
            AreaActivity.(area).k = k;
        end

        % 4.4 Sliding Window Trial Loop
        num_trials = max(tc_re);
        
        for ipair = 1:n_pairs
            a1 = area_pairs_to_analyze{ipair, 1}; a2 = area_pairs_to_analyze{ipair, 2};
            if ~isfield(AreaActivity, a1) || ~isfield(AreaActivity, a2), continue; end
            
            X_full = AreaActivity.(a1).data; k1 = AreaActivity.(a1).k;
            Y_full = AreaActivity.(a2).data; k2 = AreaActivity.(a2).k;
            
            cca_tr        = nan(1, num_trials);
            ifi_tr        = nan(1, num_trials);
            cca_tr_shuff  = nan(1, num_trials);
            ifi_tr_shuff  = nan(1, num_trials);
            
            for t = 1:num_trials
                % Define window dynamically
                win_bounds = (t + n_trials_window);
                valid_window = base_valid & ismember(tc_re, win_bounds(win_bounds >= 1 & win_bounds <= num_trials));
                
                if sum(valid_window) < max(k1, k2) + max(lags) + 5, continue; end
                
                X_win = X_full(valid_window, :);
                Y_win = Y_full(valid_window, :);
                
                % Fit CCA on active window
                try
                    [A, B, r] = canoncorr(X_win, Y_win);
                    cca_tr(t) = r(1);
                catch, continue; end
                
                % Project canonical variates for full continuous timeline
                U_full = (X_full - mean(X_win)) * A(:,1);
                V_full = (Y_full - mean(Y_win)) * B(:,1);
                
                % Real Masked Lag Correlation
                r_lags = calc_masked_lags(U_full, V_full, valid_window, lags);
                ifi_tr(t) = calc_ifi(r_lags, lags);
                
                % Shuffle Control (Continuous Circular Shift preserves autocorrelation)
                r_sh = nan(n_shuffles, 1); p_sh = nan(n_shuffles, 1);
                parfor ish = 1:n_shuffles
                    % Shift Y continuous array randomly to break anatomical coupling
                    shift_offset = randi([100, length(valid_window)-100]);
                    Y_full_shuff = circshift(Y_full, shift_offset);
                    Y_win_shuff  = Y_full_shuff(valid_window, :);
                    
                    try
                        [~, Bs, rs] = canoncorr(X_win, Y_win_shuff);
                        r_sh(ish) = rs(1);
                        V_full_shuff = (Y_full_shuff - mean(Y_win_shuff)) * Bs(:,1);
                        r_lags_sh = calc_masked_lags(U_full, V_full_shuff, valid_window, lags);
                        p_sh(ish) = calc_ifi(r_lags_sh, lags);
                    catch
                        continue;
                    end
                end
                cca_tr_shuff(t) = mean(r_sh, 'omitnan');
                ifi_tr_shuff(t) = mean(p_sh, 'omitnan');
            end
            
            % 4.5 Extract Epochs via Yoked LP
            lp = analysis_lp(ianimal);
            if ~isnan(lp) && lp > 10 && (lp + 9) <= num_trials
                idx_early = 1:10; idx_pre = (lp - 10) : (lp - 1); idx_post = lp : (lp + 9);
                
                group_results(ipair).trial_corr_early{ianimal} = cca_tr(idx_early);
                group_results(ipair).trial_corr_pre{ianimal}   = cca_tr(idx_pre);
                group_results(ipair).trial_corr_post{ianimal}  = cca_tr(idx_post);
                
                group_results(ipair).trial_precession_early_idx{ianimal} = ifi_tr(idx_early);
                group_results(ipair).trial_precession_pre_idx{ianimal}   = ifi_tr(idx_pre);
                group_results(ipair).trial_precession_post_idx{ianimal}  = ifi_tr(idx_post);
                
                group_results(ipair).trial_corr_early_shuff{ianimal} = cca_tr_shuff(idx_early);
                group_results(ipair).trial_corr_pre_shuff{ianimal}   = cca_tr_shuff(idx_pre);
                group_results(ipair).trial_corr_post_shuff{ianimal}  = cca_tr_shuff(idx_post);
                
                group_results(ipair).trial_precession_early_idx_shuff{ianimal} = ifi_tr_shuff(idx_early);
                group_results(ipair).trial_precession_pre_idx_shuff{ianimal}   = ifi_tr_shuff(idx_pre);
                group_results(ipair).trial_precession_post_idx_shuff{ianimal}  = ifi_tr_shuff(idx_post);
            end
        end
    catch ME
        fprintf('  Error processing %s: %s\n', filename, ME.message);
    end
end

save(save_path, 'group_results', 'is_learner', 'analysis_lp', '-v7.3');
fprintf('\nFinished! Results formatted for drop-in use with v2 plotting scripts.\n');

%% LOCAL HELPERS

function r_lags = calc_masked_lags(U, V, valid_window, lags)
    % Calculates cross-correlation iteratively only where shifting does not cross a validity boundary.
    r_lags = nan(1, length(lags));
    for il = 1:length(lags)
        L = lags(il);
        
        % Shifting V by -L aligns V(t+L) with U(t). 
        % If L > 0, we check if U precedes V (Region 1 leads).
        V_shifted = circshift(V, -L);
        mask_shifted = circshift(valid_window, -L);
        
        valid_pairs = valid_window & mask_shifted;
        if sum(valid_pairs) > 10
            r_lags(il) = corr(U(valid_pairs), V_shifted(valid_pairs), 'rows','complete');
        end
    end
end

function ifi = calc_ifi(r_lags, lags)
    % Extracts IFI based strictly on defined Eq 1.
    neg_lags = r_lags(lags < 0);
    pos_lags = r_lags(lags > 0);
    
    r_neg = mean(neg_lags, 'omitnan');
    r_pos = mean(pos_lags, 'omitnan');
    
    if (r_neg + r_pos) > 0.001
        ifi = (r_neg - r_pos) / (r_neg + r_pos);
    else
        ifi = nan;
    end
end