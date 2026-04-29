%% CCA_HC_V1_Spatial_Final_Optimized.m
% Final Spatial CCA Analysis with Dual Precession Metrics & Yoked Controls
%
% DESCRIPTION:
% Performs CCA to quantify shared variance and extracts Precession metrics.
% Splits animals into Learners and Non-Learners.
% Uses a 2-way Mixed Repeated Measures ANOVA for epoch comparisons.
%% 1. SETUP & PARAMETERS
clear; clc; close all;
% --- Paths ---
base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir = fullfile(base_dir, data_subfolder);
file_pattern = 'TF*_export.mat';
learning_file = 'animal_behaviour.mat';
% --- Save Path ---
current_date = datestr(now, 'yyyy_mm_dd');
save_path = fullfile(data_dir, sprintf('Spatial_CCA_Results_%s.mat', current_date));
% --- PCA Parameters ---
pca_selection_method = 'variance'; 
pca_variance_threshold = 90;       
n_components_reduced = 4;          
% --- Analysis Parameters ---
num_ccs_analyze = 3;         
n_trials_window = -3:3;      
n_bins_window = -3:3;        
n_shuffles = 50;             
max_shift_bins = 2;          
min_units_per_region = 5;
% Spatial Parameters
n_bins = 200;                
bin_size_cm = 2.5;           
track_length_cm = n_bins * bin_size_cm; 
landmarks_cm = [0, 50; 100, 125; 175, 200; 250, 275; 325, 350; 400, 425]; 
% --- Pairs to Analyze ---
area_pairs_to_analyze = {'CA1', 'V1'; 'CA1', 'DG'; 'CA1', 'CA3'; 'CA1', 'RSC'; ...
                         'CA1', 'SUB'; 'V1', 'RSC'; 'RSC', 'SUB'; 'CA3', 'DG'};
n_pairs = size(area_pairs_to_analyze, 1);
if isempty(gcp('nocreate'))
    parpool; 
end

% Set to true to ignore saved files and re-run CCA
force_reprocess = true; 
%% 2. LOAD & YOKE LEARNING POINTS 
file_list = dir(fullfile(data_dir, file_pattern));
n_animals = length(file_list);

lp_path = fullfile(data_dir, learning_file);
if exist(lp_path, 'file')
    dat_lp = load(lp_path);
    learning_points = dat_lp.period_experienced(:, 1);
    if isfield(dat_lp, 'animal_id')
        [~, sorting_idx] = sort(dat_lp.animal_id);
        learning_points = learning_points(sorting_idx);
    end
else
    warning('Learning file not found. Alignment will be skipped.');
    learning_points = nan(n_animals, 1);
end

if length(learning_points) < n_animals
    learning_points(end+1:n_animals) = nan;
end

is_learner = ~isnan(learning_points);
if any(is_learner)
    mean_lp = mean(learning_points(is_learner));
    sem_lp = sem(learning_points(is_learner));
else
    mean_lp = NaN; 
end

analysis_lp = learning_points;
analysis_lp(~is_learner) = round(mean_lp); 

fprintf('Identified %d Learners and %d Non-Learners. Mean LP used for yoking: %d\n', ...
    sum(is_learner), sum(~is_learner), mean_lp);

%% 3. INITIALIZE RESULTS STRUCTURE
group_results = struct('pair_name', cell(n_pairs, 1), ...
                       'all_bins_corr', cell(n_pairs, 1), ...
                       'all_bins_corr_shuff', cell(n_pairs, 1), ...
                       'all_bins_precession_idx', cell(n_pairs, 1), ...
                       'all_bins_precession_idx_shuff', cell(n_pairs, 1), ...
                       'all_bins_precession_curve', cell(n_pairs, 1), ...
                       'all_bins_precession_curve_shuff', cell(n_pairs, 1), ...
                       'trial_corr_early', cell(n_pairs, 1), ...
                       'trial_corr_pre', cell(n_pairs, 1), ...
                       'trial_corr_post', cell(n_pairs, 1), ...
                       'trial_corr_early_shuff', cell(n_pairs, 1), ...
                       'trial_corr_pre_shuff', cell(n_pairs, 1), ...
                       'trial_corr_post_shuff', cell(n_pairs, 1), ...
                       'trial_precession_early_idx', cell(n_pairs, 1), ... 
                       'trial_precession_pre_idx', cell(n_pairs, 1), ...   
                       'trial_precession_post_idx', cell(n_pairs, 1), ...
                       'trial_precession_early_idx_shuff', cell(n_pairs, 1), ...
                       'trial_precession_pre_idx_shuff', cell(n_pairs, 1), ...
                       'trial_precession_post_idx_shuff', cell(n_pairs, 1), ...
                       'trial_precession_early_curve', cell(n_pairs, 1), ... 
                       'trial_precession_pre_curve', cell(n_pairs, 1), ...   
                       'trial_precession_post_curve', cell(n_pairs, 1), ...
                       'trial_precession_early_curve_shuff', cell(n_pairs, 1), ...
                       'trial_precession_pre_curve_shuff', cell(n_pairs, 1), ...
                       'trial_precession_post_curve_shuff', cell(n_pairs, 1));
                   
for ipair = 1:n_pairs
    group_results(ipair).pair_name = sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});
end
%% 4. MAIN ANALYSIS LOOP & FILE LOADING

% Bundle current parameters for tracking and comparison
current_config = struct(...
    'pca_selection_method', pca_selection_method, ...
    'pca_variance_threshold', pca_variance_threshold, ...
    'n_components_reduced', n_components_reduced, ...
    'num_ccs_analyze', num_ccs_analyze, ...
    'max_shift_bins', max_shift_bins, ...
    'min_units_per_region', min_units_per_region, ...
    'n_bins', n_bins, ...
    'n_shuffles', n_shuffles);

existing_files = dir(fullfile(data_dir, 'Spatial_CCA_Results_*.mat'));

% --> MODIFIED: Bypass loading if force_reprocess is true
if ~isempty(existing_files) && ~force_reprocess
    
    % Load the most recently modified file
    [~, latest_idx] = max([existing_files.datenum]);
    load_target = fullfile(data_dir, existing_files(latest_idx).name);
    
    % Check if previous run saved its configuration
    vars_in_file = whos('-file', load_target);
    has_config = ismember('saved_config', {vars_in_file.name});
    
    if has_config
        load(load_target, 'group_results', 'is_learner', 'analysis_lp', 'saved_config');
        fprintf('Found existing results. Loading: %s\n', existing_files(latest_idx).name);
        
        fprintf('\n--- Loaded Configuration ---\n');
        fields = fieldnames(saved_config);
        mismatch_found = false;
        for i = 1:numel(fields)
            f = fields{i};
            val_saved = saved_config.(f);
            
            % Print formatted value
            if ischar(val_saved)
                fprintf('  %-25s : %s\n', f, val_saved);
            else
                fprintf('  %-25s : %g\n', f, val_saved);
            end
            
            % Check for mismatches against active script parameters
            if isfield(current_config, f) && ~isequal(saved_config.(f), current_config.(f))
                val_curr = current_config.(f);
                if ischar(val_curr)
                    fprintf('      *** WARNING: Active script is set to ''%s'' ***\n', val_curr);
                else
                    fprintf('      *** WARNING: Active script is set to %g ***\n', val_curr);
                end
                mismatch_found = true;
            end
        end
        fprintf('----------------------------\n');
        
        if mismatch_found
            warning('Loaded results parameters DO NOT MATCH active script parameters. Plotting will proceed using LOADED data.');
        end
    else
        load(load_target, 'group_results', 'is_learner', 'analysis_lp');
        fprintf('Found existing results (Legacy file, no config saved). Loading: %s\n', existing_files(latest_idx).name);
    end
else
    
    if force_reprocess
        fprintf('Force reprocess enabled. Starting CCA analysis from scratch...\n');
    else
        fprintf('No existing results found. Starting CCA analysis...\n');
    end
    
    n_shifts = 2 * max_shift_bins + 1; 
    for ianimal = 1:n_animals
        filename = file_list(ianimal).name;
        fullpath = fullfile(data_dir, filename);
        fprintf('\nProcessing Animal %d/%d: %s\n', ianimal, n_animals, filename);
        
        try
            D = load(fullpath);
            if ~isfield(D, 'units') || ~isfield(D, 'analysis_spatial'), continue; end
            units = D.units;
            
            keep_mask = true(length(units.unit_id), 1);
            if isfield(units, 'idx_fs')
                target_fs_areas = {'V1', 'RSC', 'CA1', 'CA3'};
                is_fs = logical(units.idx_fs);
                for r = 1:length(target_fs_areas)
                    fs_in_area = strcmp(units.region, target_fs_areas{r}) & is_fs;
                    keep_mask(fs_in_area) = false;
                end
            end
            units.idx(~keep_mask, :) = 0;
            
            if isfield(D.analysis_spatial, 'firing') && isfield(D.analysis_spatial.firing, 'cued')
                raw_spatial = D.analysis_spatial.firing.cued.freq_z;
                num_raw_nans = sum(isnan(raw_spatial(:)));
                if num_raw_nans > 0
                    fprintf('  [NaN Tracker] Animal %d (%s): Raw firing data contains %d NaNs.\n', ...
                        ianimal, filename, num_raw_nans);
                end
                animal_data = permute(raw_spatial, [1, 3, 2]);
                [~, curr_n_bins, num_trials] = size(animal_data);
                if curr_n_bins ~= n_bins
                    warning('Bin mismatch for %s. Skipping.', filename); continue;
                end
            else, continue; end
            
            animal_areas = unique(units.regions_label);
            AreaActivity = struct();
            
            for ia = 1:length(animal_areas)
                area_name = animal_areas{ia};
                u_idx = units.idx(strcmp(units.regions_label, area_name), :);
                if sum(u_idx) < min_units_per_region, continue; end
                
                area_dat = animal_data(u_idx, :, :);
                reshaped_dat = reshape(area_dat, sum(u_idx), [])'; 
                
                if size(reshaped_dat, 2) >= 2 
                    [~, scores, ~, ~, explained] = pca(reshaped_dat);
                    if strcmp(pca_selection_method, 'variance')
                        cum_var = cumsum(explained);
                        n_comps = find(cum_var >= pca_variance_threshold, 1);
                        if isempty(n_comps), n_comps = size(scores, 2); end
                    else, n_comps = n_components_reduced; end
                    
                    n_comps = max(n_comps, num_ccs_analyze); 
                    n_comps = min(n_comps, size(scores, 2));
                    
                    AreaActivity.(area_name).data = reshape(scores(:, 1:n_comps)', n_comps, n_bins, num_trials);
                    AreaActivity.(area_name).n_comps = n_comps;
                end
            end
            
            for ipair = 1:n_pairs
                a1 = area_pairs_to_analyze{ipair, 1};
                a2 = area_pairs_to_analyze{ipair, 2};
                
                if ~isfield(AreaActivity, a1) || ~isfield(AreaActivity, a2), continue; end
                
                d1 = AreaActivity.(a1).data; d2 = AreaActivity.(a2).data;
                nc1 = AreaActivity.(a1).n_comps; nc2 = AreaActivity.(a2).n_comps;
                
                cca_tr = nan(num_ccs_analyze, num_trials);
                cca_tr_shuff = nan(num_ccs_analyze, num_trials); 
                prec_tr_idx = nan(num_ccs_analyze, num_trials);
                prec_tr_idx_shuff = nan(num_ccs_analyze, num_trials); 
                prec_tr_curve = nan(num_ccs_analyze, n_shifts, num_trials);
                prec_tr_curve_shuff = nan(num_ccs_analyze, n_shifts, num_trials);
                
                cca_bin = nan(num_ccs_analyze, n_bins);
                cca_bin_shuff = nan(num_ccs_analyze, n_bins); 
                prec_bin_idx = nan(num_ccs_analyze, n_bins);
                prec_bin_idx_shuff = nan(num_ccs_analyze, n_bins);
                prec_bin_curve = nan(num_ccs_analyze, n_shifts, n_bins);
                prec_bin_curve_shuff = nan(num_ccs_analyze, n_shifts, n_bins);
                
                % =========================================================
                % ANALYSIS 1: TRIAL-WISE
                % =========================================================
                for t = 1:num_trials
                    win = t + n_trials_window; win = win(win>=1 & win<=num_trials);
                    if isempty(win), continue; end
                    
                    D1_local = d1(:, :, win); D2_local = d2(:, :, win);
                    x = reshape(D1_local, nc1, []); y = reshape(D2_local, nc2, []);
                    
                    if check_dims(x, y, max(nc1, nc2) + 5)
                        [~,~,r] = canoncorr(x', y');
                        cca_tr(1:min(length(r), num_ccs_analyze), t) = r(1:min(length(r), num_ccs_analyze));
                        [c_curve, c_idx] = calc_precession(D1_local, D2_local, max_shift_bins, num_ccs_analyze, nc1, nc2);
                        prec_tr_curve(:, :, t) = c_curve;
                        prec_tr_idx(:, t) = c_idx;
                        r_sh_iter = nan(num_ccs_analyze, n_shuffles);
                        p_sh_idx_iter = nan(num_ccs_analyze, n_shuffles);
                        p_sh_curve_iter = nan(num_ccs_analyze, n_shifts, n_shuffles);
                        n_bins_local = size(D1_local, 2);
                        parfor ish = 1:n_shuffles
                            perm_idx = randperm(n_bins_local);
                            D1_s = D1_local(:, perm_idx, :); x_s = reshape(D1_s, nc1, []);
                            [~,~,rs] = canoncorr(x_s', y');
                            r_sh_iter(:, ish) = rs(1:min(length(rs), num_ccs_analyze));
                            [c_curve_sh, c_idx_sh] = calc_precession(D1_s, D2_local, max_shift_bins, num_ccs_analyze, nc1, nc2);
                            p_sh_curve_iter(:, :, ish) = c_curve_sh;
                            p_sh_idx_iter(:, ish) = c_idx_sh;
                        end
                        cca_tr_shuff(:, t) = mean(r_sh_iter, 2, 'omitnan');
                        prec_tr_curve_shuff(:, :, t) = mean(p_sh_curve_iter, 3, 'omitnan');
                        prec_tr_idx_shuff(:, t) = mean(p_sh_idx_iter, 2, 'omitnan');
                    else
                        if any(isnan(x(:))) || any(isnan(y(:)))
                            fprintf('  [NaN Tracker] Trial %d skipped (Pair %d): NaNs present in PCA scores.\n', t, ipair);
                        else
                            fprintf('  [NaN Tracker] Trial %d skipped (Pair %d): Insufficient samples (Size: %d).\n', t, ipair, size(x, 2));
                        end
                    end
                end
                
                % =========================================================
                % ANALYSIS 2: BIN-WISE
                % =========================================================
                for b = 1:n_bins
                    win = b + n_bins_window; win = win(win>=1 & win<=n_bins);
                    if isempty(win), continue; end
                    
                    D1_local = d1(:, win, :); D2_local = d2(:, win, :);
                    x = reshape(D1_local, nc1, []); y = reshape(D2_local, nc2, []);
                    
                    if check_dims(x, y, max(nc1, nc2) + 5)
                        [~,~,r] = canoncorr(x', y');
                        cca_bin(1:min(length(r), num_ccs_analyze), b) = r(1:min(length(r), num_ccs_analyze));
                        [c_curve, c_idx] = calc_precession(D1_local, D2_local, max_shift_bins, num_ccs_analyze, nc1, nc2);
                        prec_bin_curve(:, :, b) = c_curve;
                        prec_bin_idx(:, b) = c_idx;
                        
                        r_sh_iter = nan(num_ccs_analyze, n_shuffles);
                        p_sh_idx_iter = nan(num_ccs_analyze, n_shuffles);
                        p_sh_curve_iter = nan(num_ccs_analyze, n_shifts, n_shuffles);
                        n_tr_local = size(D1_local, 3);
                        
                        parfor ish = 1:n_shuffles
                             perm_idx = randperm(n_tr_local);
                             D1_s = D1_local(:, :, perm_idx); x_s = reshape(D1_s, nc1, []);
                             [~,~,rs] = canoncorr(x_s', y');
                             r_sh_iter(:, ish) = rs(1:min(length(rs), num_ccs_analyze));
                             [c_curve_sh, c_idx_sh] = calc_precession(D1_s, D2_local, max_shift_bins, num_ccs_analyze, nc1, nc2);
                             p_sh_curve_iter(:, :, ish) = c_curve_sh;
                             p_sh_idx_iter(:, ish) = c_idx_sh;
                        end
                        cca_bin_shuff(:, b) = mean(r_sh_iter, 2, 'omitnan');
                        prec_bin_curve_shuff(:, :, b) = mean(p_sh_curve_iter, 3, 'omitnan');
                        prec_bin_idx_shuff(:, b) = mean(p_sh_idx_iter, 2, 'omitnan');
                    end
                end
                
                % --- F. Store Results ---
                group_results(ipair).all_bins_corr{ianimal} = cca_bin;
                group_results(ipair).all_bins_corr_shuff{ianimal} = cca_bin_shuff;
                group_results(ipair).all_bins_precession_idx{ianimal} = prec_bin_idx;
                group_results(ipair).all_bins_precession_idx_shuff{ianimal} = prec_bin_idx_shuff; 
                group_results(ipair).all_bins_precession_curve{ianimal} = prec_bin_curve;
                group_results(ipair).all_bins_precession_curve_shuff{ianimal} = prec_bin_curve_shuff; 
                
                % Process Epochs using Yoked Analysis LP
                lp = analysis_lp(ianimal);
                if ~isnan(lp) && lp > 10 && (lp + 9) <= num_trials
                    idx_early = 1:10; idx_pre = (lp - 10) : (lp - 1); idx_post = lp : (lp + 9);
                    
                    get_cols_2d = @(data, cols) data(:, cols);
                    get_cols_3d = @(data, cols) data(:, :, cols);
                    
                    group_results(ipair).trial_corr_early{ianimal} = get_cols_2d(cca_tr, idx_early);
                    group_results(ipair).trial_corr_pre{ianimal}   = get_cols_2d(cca_tr, idx_pre);
                    group_results(ipair).trial_corr_post{ianimal}  = get_cols_2d(cca_tr, idx_post);
                    group_results(ipair).trial_corr_early_shuff{ianimal} = get_cols_2d(cca_tr_shuff, idx_early);
                    group_results(ipair).trial_corr_pre_shuff{ianimal}   = get_cols_2d(cca_tr_shuff, idx_pre);
                    group_results(ipair).trial_corr_post_shuff{ianimal}  = get_cols_2d(cca_tr_shuff, idx_post);
                    
                    group_results(ipair).trial_precession_early_idx{ianimal} = get_cols_2d(prec_tr_idx, idx_early);
                    group_results(ipair).trial_precession_pre_idx{ianimal}   = get_cols_2d(prec_tr_idx, idx_pre);
                    group_results(ipair).trial_precession_post_idx{ianimal}  = get_cols_2d(prec_tr_idx, idx_post);
                    group_results(ipair).trial_precession_early_idx_shuff{ianimal} = get_cols_2d(prec_tr_idx_shuff, idx_early);
                    group_results(ipair).trial_precession_pre_idx_shuff{ianimal}   = get_cols_2d(prec_tr_idx_shuff, idx_pre);
                    group_results(ipair).trial_precession_post_idx_shuff{ianimal}  = get_cols_2d(prec_tr_idx_shuff, idx_post);
                    
                    group_results(ipair).trial_precession_early_curve{ianimal} = get_cols_3d(prec_tr_curve, idx_early);
                    group_results(ipair).trial_precession_pre_curve{ianimal}   = get_cols_3d(prec_tr_curve, idx_pre);
                    group_results(ipair).trial_precession_post_curve{ianimal}  = get_cols_3d(prec_tr_curve, idx_post);
                    group_results(ipair).trial_precession_early_curve_shuff{ianimal} = get_cols_3d(prec_tr_curve_shuff, idx_early);
                    group_results(ipair).trial_precession_pre_curve_shuff{ianimal}   = get_cols_3d(prec_tr_curve_shuff, idx_pre);
                    group_results(ipair).trial_precession_post_curve_shuff{ianimal}  = get_cols_3d(prec_tr_curve_shuff, idx_post);
                end
            end 
        catch ME
            fprintf('Error processing %s: %s\n', filename, ME.message);
        end
    end
    saved_config = current_config;
    save(save_path, 'group_results', 'is_learner', 'analysis_lp', 'saved_config', '-v7.3');
    fprintf('Finished analysis and saving! \n')
end

%% 5. PLOTTING (Summary Figures)
fprintf('\nGenerating Plots...\n');
lags = -max_shift_bins:max_shift_bins;

% --- NEW FIX: Pad cell arrays to prevent out-of-bounds indexing ---
n_total_animals = length(is_learner); 
for ipair = 1:n_pairs
    fields = fieldnames(group_results);
    for f = 1:length(fields)
        if iscell(group_results(ipair).(fields{f}))
            if length(group_results(ipair).(fields{f})) < n_total_animals
                % Assigning an empty array to the last index forces MATLAB 
                % to pad the missing intermediate indices with []
                group_results(ipair).(fields{f}){n_total_animals} = [];
            end
        end
    end
end

% % --- A1. Spatial Correlation (Bin-wise) ---
% figure('Name', 'Group Spatial Correlation', 'Color', 'w', 'Position', [100 100 1200 800]);
% tiledlayout('flow');
% for ipair = 1:n_pairs
%     nexttile
%     [mu, se] = aggregate_cells(group_results(ipair).all_bins_corr, n_bins, 1); 
%     [mu_s, ~] = aggregate_cells(group_results(ipair).all_bins_corr_shuff, n_bins, 1);
%     if isempty(mu), continue; end
%     plot(1:n_bins, mu_s, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
%     shadedErrorBar(1:n_bins, mu, se, 'lineProps', {'Color', 'b'});
%     xline((landmarks_cm(:) / track_length_cm) * n_bins); 
%     title(group_results(ipair).pair_name); xlabel('Spatial Bin'); ylabel('Correlation (CC1)'); xlim([1 n_bins]);
% end
% save_to_svg(fullfile(data_dir, 'Spatial_Corr_Group'));
% 
% % --- A2. Spatial Precession (Average Index) ---
% figure('Name', 'Group Spatial Precession Index', 'Color', 'w', 'Position', [100 100 1200 800]);
% tiledlayout('flow');
% for ipair = 1:n_pairs
%     nexttile
%     [mu, se] = aggregate_cells(group_results(ipair).all_bins_precession_idx, n_bins, 1);
%     [mu_s, ~] = aggregate_cells(group_results(ipair).all_bins_precession_idx_shuff, n_bins, 1); 
%     if isempty(mu), continue; end
%     plot(1:n_bins, mu_s, 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'LineStyle', '--'); hold on;
%     shadedErrorBar(1:n_bins, mu, se, 'lineProps', {'Color', 'r'});
%     xline((landmarks_cm(:) / track_length_cm) * n_bins); yline(0, '--k');
%     title(group_results(ipair).pair_name); xlim([1 n_bins]); ylim([-0.2 0.2]); ylabel('Precession Idx (CC1)');
% end
% save_to_svg(fullfile(data_dir, 'Spatial_Precession_Idx_Group'));

% % --- A3. Spatial Precession (Lag vs Space Heatmap) ---
% figure('Name', 'Group Spatial Lag-Tuning', 'Color', 'w', 'Position', [100 100 1200 800]);
% tiledlayout('flow');
% for ipair = 1:n_pairs
%     nexttile
%     mu_real = aggregate_shift_heatmaps(group_results(ipair).all_bins_precession_curve, n_bins, 1);
%     mu_shuff = aggregate_shift_heatmaps(group_results(ipair).all_bins_precession_curve_shuff, n_bins, 1);
%     if isempty(mu_real), continue; end
%     imagesc(1:n_bins, lags, mu_real - mu_shuff); axis xy; colormap(gca, 'parula');
%     c = colorbar; c.Label.String = '\Delta Corr (Real-Shuff)';
%     xline((landmarks_cm(:) / track_length_cm) * n_bins, 'w--'); yline(0, 'w:');
%     title(group_results(ipair).pair_name); xlabel('Spatial Bin');
%     ylabel(sprintf('%s Leads <-Lag-> %s Leads', area_pairs_to_analyze{ipair,1}, area_pairs_to_analyze{ipair,2}));
% end
% save_to_svg(fullfile(data_dir, 'Spatial_Precession_Heatmaps'));

% --- B1. Epoch Correlation (Bar Graph with Mixed rmANOVA) ---
fprintf('\n--- Epoch Correlation Mixed rmANOVA Results ---\n');
figure('Name', 'Epoch Correlation (Learners vs Non)', 'Color', 'w', 'Position', [100 100 1400 800]);
tiledlayout('flow');
for ipair = 1:n_pairs
    e_vals = extract_animal_means(group_results(ipair).trial_corr_early, 1);
    p_vals = extract_animal_means(group_results(ipair).trial_corr_pre, 1);
    x_vals = extract_animal_means(group_results(ipair).trial_corr_post, 1);
    
    if all(isnan(e_vals)), continue; end
    
    nexttile;
    fprintf('\nPair: %s\n', group_results(ipair).pair_name);
    plot_grouped_bars_with_rmanova(e_vals, p_vals, x_vals, is_learner, group_results(ipair).pair_name, 'Correlation (CC1)', false);
end
save_to_svg(fullfile(data_dir, 'Epoch_Corr_Bars_Split'));

% --- B2. Epoch Precession Index (Bar Graph with Mixed rmANOVA) ---
fprintf('\n--- Epoch Precession Index Mixed rmANOVA Results ---\n');
figure('Name', 'Epoch Precession (Learners vs Non)', 'Color', 'w', 'Position', [100 100 1400 800]);
tiledlayout('flow');
for ipair = 1:n_pairs
    e_vals = extract_animal_means(group_results(ipair).trial_precession_early_idx, 1);
    p_vals = extract_animal_means(group_results(ipair).trial_precession_pre_idx, 1);
    x_vals = extract_animal_means(group_results(ipair).trial_precession_post_idx, 1);
    
    if all(isnan(e_vals)), continue; end
    
    nexttile;
    fprintf('\nPair: %s\n', group_results(ipair).pair_name);
    plot_grouped_bars_with_rmanova(e_vals, p_vals, x_vals, is_learner, group_results(ipair).pair_name, 'Precession Index', true);
end
save_to_svg(fullfile(data_dir, 'Epoch_Precession_Bars_Split'));

% --- C1. Continuous Trial-by-Trial Precession (3 Epochs Linked) ---
fprintf('\n--- Continuous Trial Precession Stats (Mixed rmANOVA, Uncorrected Local Tests) ---\n');
figure('Name', 'Continuous Epoch Precession', 'Color', 'w', 'Position', [150 150 1500 800]);
t = tiledlayout('flow');
for ipair = 1:n_pairs
    early_mat = extract_epoch_trials(group_results(ipair).trial_precession_early_idx, n_animals);
    pre_mat   = extract_epoch_trials(group_results(ipair).trial_precession_pre_idx, n_animals);
    post_mat  = extract_epoch_trials(group_results(ipair).trial_precession_post_idx, n_animals);
    
    if isempty(early_mat), continue; end
    
    nexttile; hold on;
    
    plot_epoch_line(1:10, early_mat(is_learner, :), [0.85 0.33 0.1], '-');
    plot_epoch_line(11:20, pre_mat(is_learner, :), [0.85 0.33 0.1], '-');
    plot_epoch_line(21:30, post_mat(is_learner, :), [0.85 0.33 0.1], '-');
    
    if any(~is_learner)
        plot_epoch_line(1:10, early_mat(~is_learner, :), [0.5 0.5 0.5], '-');
        plot_epoch_line(11:20, pre_mat(~is_learner, :), [0.5 0.5 0.5], '-');
        plot_epoch_line(21:30, post_mat(~is_learner, :), [0.5 0.5 0.5], '-');
    end
    
    xline(10.5, ':k'); xline(20.5, ':k'); yline(0, '-k');
    xticks([5, 15, 25]); xticklabels({'Early', 'Pre-Learn', 'Post-Learn'});
    xlim([1 30]);
    title(sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2}));
    ylabel('Trial-by-Trial Precession');
    
    % --- Stats Setup ---
    full_mat = [early_mat, pre_mat, post_mat]; % [n_animals x 30]
    
    yl = ylim; offset = (yl(2) - yl(1)) * 0.05;
    y_L_zero = yl(2) + offset;      
    y_NL_zero = yl(2) + 2*offset;   
    y_diff = yl(2) + 3*offset;    
    
    % --- 1. Uncorrected Local T-Tests ---
    for tr = 1:30
        dat_L = full_mat(is_learner, tr);
        if sum(~isnan(dat_L)) > 2
            [~, p_L] = ttest(dat_L, 0);
            if p_L < 0.05, plot(tr, y_L_zero, '.', 'Color', [0.85 0.33 0.1], 'MarkerSize', 10); end
        end
        
        if any(~is_learner)
            dat_NL = full_mat(~is_learner, tr);
            if sum(~isnan(dat_NL)) > 2
                [~, p_NL] = ttest(dat_NL, 0);
                if p_NL < 0.05, plot(tr, y_NL_zero, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 10); end
            end
            
            if sum(~isnan(dat_L)) > 2 && sum(~isnan(dat_NL)) > 2
                [~, p_diff] = ttest2(dat_L, dat_NL);
                if p_diff < 0.05, plot(tr, y_diff, '*', 'Color', 'k', 'MarkerSize', 5); end
            end
        end
    end
    
    % Add legend proxies
    text(30.5, y_L_zero, 'L \neq 0', 'Color', [0.85 0.33 0.1], 'FontSize', 8, 'VerticalAlignment', 'middle');
    if any(~is_learner)
        text(30.5, y_NL_zero, 'NL \neq 0', 'Color', [0.5 0.5 0.5], 'FontSize', 8, 'VerticalAlignment', 'middle');
        text(30.5, y_diff, 'L \neq NL', 'Color', 'k', 'FontSize', 8, 'VerticalAlignment', 'middle');
    end
    ylim([yl(1), yl(2) + 4.5*offset]); 
    
    % --- 2. 30-Trial Mixed Repeated Measures ANOVA ---
    fprintf('\nPair: %s\n', area_pairs_to_analyze{ipair, 1});
    
    % Find complete cases (animals with zero NaNs across all 30 trials)
    complete_cases = all(~isnan(full_mat), 2);
    n_complete = sum(complete_cases);
    
    if n_complete > max(2, length(unique(is_learner(complete_cases))))
        % Build Table
        varNames = arrayfun(@(x) sprintf('Trial%d', x), 1:30, 'UniformOutput', false);
        varNames{end+1} = 'LearnerGroup';
        
        T_data = [full_mat(complete_cases, :), is_learner(complete_cases)];
        T = array2table(T_data, 'VariableNames', varNames);
        T.LearnerGroup = categorical(T.LearnerGroup);
        
        Meas = table((1:30)', 'VariableNames', {'Trial'});
        
        try
            rm = fitrm(T, 'Trial1-Trial30 ~ LearnerGroup', 'WithinDesign', Meas);
            
            % Within-subjects (Trial, Interaction)
            ranovatbl = ranova(rm);
            p_trial = ranovatbl.pValue(1);
            p_int   = ranovatbl.pValue(2);
            
            % Between-subjects (Group)
            betweentbl = anova(rm);
            p_group = betweentbl.pValue(2); 
            
            disp('--- Within-Subjects (Trial & Interaction) ---');
            disp(ranovatbl);
            disp('--- Between-Subjects (Group) ---');
            disp(betweentbl);
            
            sig_text = {sprintf('rmANOVA (n=%d complete):', n_complete)};
            has_sig = false;
            
            if p_trial < 0.05, sig_text{end+1} = sprintf('Trial: p=%.3f', p_trial); has_sig = true; end
            if p_group < 0.05, sig_text{end+1} = sprintf('Group: p=%.3f', p_group); has_sig = true; end
            if p_int < 0.05,   sig_text{end+1} = sprintf('Tr x Grp: p=%.3f', p_int); has_sig = true; end
            if ~has_sig, sig_text{end+1} = 'n.s.'; end
            
            curr_yl = ylim; curr_xl = xlim;
            text(curr_xl(1) + 1, curr_yl(2) - 0.5*offset, ...
                 strjoin(sig_text, '\n'), 'VerticalAlignment', 'top', ...
                 'FontSize', 8, 'EdgeColor', 'k', 'BackgroundColor', 'w', 'Margin', 3);
                 
        catch ME
            fprintf('rmANOVA failed: %s\n', ME.message);
        end
    else
        fprintf('Not enough complete cases (n=%d) across all 30 trials for strict rmANOVA.\n', n_complete);
    end
end
save_to_svg(fullfile(data_dir, 'Continuous_Trial_Precession_Split'));

% --- C2. Continuous Trial-by-Trial CC1 (Correlation) ---
fprintf('\n--- Continuous Trial Correlation Stats (Mixed rmANOVA, Uncorrected Local Tests) ---\n');
figure('Name', 'Continuous Epoch Correlation', 'Color', 'w', 'Position', [150 150 1500 800]);
t = tiledlayout('flow');
for ipair = 1:n_pairs
    early_mat = extract_epoch_trials(group_results(ipair).trial_corr_early, n_animals);
    pre_mat   = extract_epoch_trials(group_results(ipair).trial_corr_pre, n_animals);
    post_mat  = extract_epoch_trials(group_results(ipair).trial_corr_post, n_animals);
    
    if isempty(early_mat), continue; end
    
    nexttile; hold on;
    
    plot_epoch_line(1:10, early_mat(is_learner, :), [0.85 0.33 0.1], '-');
    plot_epoch_line(11:20, pre_mat(is_learner, :), [0.85 0.33 0.1], '-');
    plot_epoch_line(21:30, post_mat(is_learner, :), [0.85 0.33 0.1], '-');
    
    if any(~is_learner)
        plot_epoch_line(1:10, early_mat(~is_learner, :), [0.5 0.5 0.5], '-');
        plot_epoch_line(11:20, pre_mat(~is_learner, :), [0.5 0.5 0.5], '-');
        plot_epoch_line(21:30, post_mat(~is_learner, :), [0.5 0.5 0.5], '-');
    end
    
    xline(10.5, ':k'); xline(20.5, ':k');
    xticks([5, 15, 25]); xticklabels({'Early', 'Pre-Learn', 'Post-Learn'});
    xlim([1 30]);
    title(sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2}));
    ylabel('Trial-by-Trial Correlation (CC1)');
    
    % --- Stats Setup ---
    full_mat = [early_mat, pre_mat, post_mat]; % [n_animals x 30]
    
    yl = ylim; offset = (yl(2) - yl(1)) * 0.05;
    y_diff = yl(2) + offset;    
    
    % --- 1. Uncorrected Local T-Tests ---
    if any(~is_learner)
        for tr = 1:30
            dat_L = full_mat(is_learner, tr);
            dat_NL = full_mat(~is_learner, tr);
            if sum(~isnan(dat_L)) > 2 && sum(~isnan(dat_NL)) > 2
                [~, p_diff] = ttest2(dat_L, dat_NL);
                if p_diff < 0.05, plot(tr, y_diff, '*', 'Color', 'k', 'MarkerSize', 5); end
            end
        end
        text(30.5, y_diff, 'L \neq NL', 'Color', 'k', 'FontSize', 8, 'VerticalAlignment', 'middle');
    end
    
    ylim([yl(1), yl(2) + 2.5*offset]); 
    
    % --- 2. 30-Trial Mixed Repeated Measures ANOVA ---
    fprintf('\nPair: %s\n', area_pairs_to_analyze{ipair, 1});
    complete_cases = all(~isnan(full_mat), 2);
    n_complete = sum(complete_cases);
    
    if n_complete > max(2, length(unique(is_learner(complete_cases))))
        varNames = arrayfun(@(x) sprintf('Trial%d', x), 1:30, 'UniformOutput', false);
        varNames{end+1} = 'LearnerGroup';
        
        T_data = [full_mat(complete_cases, :), is_learner(complete_cases)];
        T = array2table(T_data, 'VariableNames', varNames);
        T.LearnerGroup = categorical(T.LearnerGroup);
        
        Meas = table((1:30)', 'VariableNames', {'Trial'});
        
        try
            rm = fitrm(T, 'Trial1-Trial30 ~ LearnerGroup', 'WithinDesign', Meas);
            
            ranovatbl = ranova(rm);
            p_trial = ranovatbl.pValue(1);
            p_int   = ranovatbl.pValue(2);
            
            betweentbl = anova(rm);
            p_group = betweentbl.pValue(2); 
            
            disp('--- Within-Subjects (Trial & Interaction) ---');
            disp(ranovatbl);
            disp('--- Between-Subjects (Group) ---');
            disp(betweentbl);
            
            sig_text = {sprintf('rmANOVA (n=%d complete):', n_complete)};
            has_sig = false;
            
            if p_trial < 0.05, sig_text{end+1} = sprintf('Trial: p=%.3f', p_trial); has_sig = true; end
            if p_group < 0.05, sig_text{end+1} = sprintf('Group: p=%.3f', p_group); has_sig = true; end
            if p_int < 0.05,   sig_text{end+1} = sprintf('Tr x Grp: p=%.3f', p_int); has_sig = true; end
            if ~has_sig, sig_text{end+1} = 'n.s.'; end
            
            curr_yl = ylim; curr_xl = xlim;
            text(curr_xl(1) + 1, curr_yl(2) - 0.5*offset, ...
                 strjoin(sig_text, '\n'), 'VerticalAlignment', 'top', ...
                 'FontSize', 8, 'EdgeColor', 'k', 'BackgroundColor', 'w', 'Margin', 3);
        catch ME
            fprintf('rmANOVA failed: %s\n', ME.message);
        end
    else
        fprintf('Not enough complete cases (n=%d) across all 30 trials for strict rmANOVA.\n', n_complete);
    end
end
save_to_svg(fullfile(data_dir, 'Continuous_Trial_Corr_Split'));
%

% % --- D. Dual-Encoded Network Visualization (Learners vs Non-Learners) ---
% fprintf('\nGenerating Network Graphs (Real vs Shuffled Filtered)...\n');
% layout_def.names = {'CA1', 'V1', 'DG', 'CA3', 'RSC', 'SUB'};
% layout_def.x     = [7.0,  8.0,  3.0,  5.0,  2.0,  2.0];
% layout_def.y     = [5.0,  9.0,  4.0,  2.0,  9.0,  6.5];
% 
% % Preallocate tracking structures
% net = struct();
% metrics = {'L_early_cc', 'L_post_cc', 'L_early_ifi', 'L_post_ifi', ...
%            'NL_early_cc', 'NL_post_cc', 'NL_early_ifi', 'NL_post_ifi', ...
%            'L_early_pval', 'L_post_pval', 'NL_early_pval', 'NL_post_pval'};
% for m = metrics, net.(m{1}) = nan(n_pairs, 1); end
% 
% % Extract mean state values and run paired t-tests (Real vs Shuffled)
% for ipair = 1:n_pairs
%     e_cc = extract_animal_means(group_results(ipair).trial_corr_early, 1);
%     x_cc = extract_animal_means(group_results(ipair).trial_corr_post, 1);
% 
%     e_ifi   = extract_animal_means(group_results(ipair).trial_precession_early_idx, 1);
%     e_ifi_s = extract_animal_means(group_results(ipair).trial_precession_early_idx_shuff, 1);
%     x_ifi   = extract_animal_means(group_results(ipair).trial_precession_post_idx, 1);
%     x_ifi_s = extract_animal_means(group_results(ipair).trial_precession_post_idx_shuff, 1);
% 
%     % Means
%     net.L_early_cc(ipair)  = mean(e_cc(is_learner), 'omitnan');
%     net.L_post_cc(ipair)   = mean(x_cc(is_learner), 'omitnan');
%     net.L_early_ifi(ipair) = mean(e_ifi(is_learner), 'omitnan');
%     net.L_post_ifi(ipair)  = mean(x_ifi(is_learner), 'omitnan');
% 
%     % Significance testing: Paired Real vs Shuffled
%     valid_L_e = ~isnan(e_ifi(is_learner)) & ~isnan(e_ifi_s(is_learner));
%     if sum(valid_L_e) > 2
%         [~, net.L_early_pval(ipair)] = ttest(e_ifi(is_learner(valid_L_e)), e_ifi_s(is_learner(valid_L_e))); 
%     end
% 
%     valid_L_x = ~isnan(x_ifi(is_learner)) & ~isnan(x_ifi_s(is_learner));
%     if sum(valid_L_x) > 2
%         [~, net.L_post_pval(ipair)] = ttest(x_ifi(is_learner(valid_L_x)), x_ifi_s(is_learner(valid_L_x))); 
%     end
% 
%     if any(~is_learner)
%         net.NL_early_cc(ipair)  = mean(e_cc(~is_learner), 'omitnan');
%         net.NL_post_cc(ipair)   = mean(x_cc(~is_learner), 'omitnan');
%         net.NL_early_ifi(ipair) = mean(e_ifi(~is_learner), 'omitnan');
%         net.NL_post_ifi(ipair)  = mean(x_ifi(~is_learner), 'omitnan');
% 
%         valid_NL_e = ~isnan(e_ifi(~is_learner)) & ~isnan(e_ifi_s(~is_learner));
%         if sum(valid_NL_e) > 2
%             [~, net.NL_early_pval(ipair)] = ttest(e_ifi(~is_learner(valid_NL_e)), e_ifi_s(~is_learner(valid_NL_e))); 
%         end
% 
%         valid_NL_x = ~isnan(x_ifi(~is_learner)) & ~isnan(x_ifi_s(~is_learner));
%         if sum(valid_NL_x) > 2
%             [~, net.NL_post_pval(ipair)] = ttest(x_ifi(~is_learner(valid_NL_x)), x_ifi_s(~is_learner(valid_NL_x))); 
%         end
%     end
% end
% 
% % Global Maxima for Scaling
% global_max_cc = max(max([net.L_early_cc, net.L_post_cc, net.NL_early_cc, net.NL_post_cc]));
% global_max_ifi = max(max(abs([net.L_early_ifi, net.L_post_ifi, net.NL_early_ifi, net.NL_post_ifi])));
% if global_max_cc == 0 || isnan(global_max_cc), global_max_cc = 1; end
% if global_max_ifi == 0 || isnan(global_max_ifi), global_max_ifi = 1; end
% 
% figure('Name', 'Network Evolution', 'Color', 'w', 'Position', [100 100 1000 800]);
% if any(~is_learner)
%     t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% else
%     t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% end
% 
% % Plot Learners
% nexttile(t); 
% plot_dual_network(area_pairs_to_analyze, net.L_early_cc, net.L_early_ifi, net.L_early_pval, 'Learners: Naive (Early)', layout_def, global_max_cc, global_max_ifi);
% nexttile(t); 
% plot_dual_network(area_pairs_to_analyze, net.L_post_cc, net.L_post_ifi, net.L_post_pval, 'Learners: Expert (Post)', layout_def, global_max_cc, global_max_ifi);
% 
% % Plot Non-Learners
% if any(~is_learner)
%     nexttile(t); 
%     plot_dual_network(area_pairs_to_analyze, net.NL_early_cc, net.NL_early_ifi, net.NL_early_pval, 'Non-Learners: Naive (Early)', layout_def, global_max_cc, global_max_ifi);
%     nexttile(t); 
%     plot_dual_network(area_pairs_to_analyze, net.NL_post_cc, net.NL_post_ifi, net.NL_post_pval, 'Non-Learners: Expert (Post)', layout_def, global_max_cc, global_max_ifi);
% end
% save_to_svg(fullfile(data_dir, 'Network_Evolution_Dual'));
% fprintf('Done.\n');


%% 5B. EXTENDED PLOTTING: MULTIPLE CCs (REAL VS SHUFFLED)
fprintf('\n--- Extended Quantifications: Real vs Shuffled (%d CCs) ---\n', num_ccs_analyze);
for g_idx = 1:2
    if g_idx == 1
        mask = is_learner; group_label = 'Learners';
    else
        mask = ~is_learner; group_label = 'Non-Learners';
    end
    
    if sum(mask) == 0, continue; end
    
    % Initialize Figures
    fig_cc = figure('Name', sprintf('CC Errorbars - %s', group_label), 'Color', 'w', 'Position', [100 100 1600 900]);
    t_cc = tiledlayout('flow', 'TileSpacing', 'compact');
    
    fig_idx = figure('Name', sprintf('Precession Idx Errorbars - %s', group_label), 'Color', 'w', 'Position', [150 150 1600 900]);
    t_idx = tiledlayout('flow', 'TileSpacing', 'compact');
    
    fig_cc_crv = figure('Name', sprintf('CC Trial Curves - %s', group_label), 'Color', 'w', 'Position', [200 200 1600 900]);
    t_cc_crv = tiledlayout('flow', 'TileSpacing', 'compact');
    
    fig_pi_crv = figure('Name', sprintf('Precession Trial Curves - %s', group_label), 'Color', 'w', 'Position', [250 250 1600 900]);
    t_pi_crv = tiledlayout('flow', 'TileSpacing', 'compact');
    
    for cc = 1:num_ccs_analyze
        for ipair = 1:n_pairs
            pair_name = group_results(ipair).pair_name;
            title_str = sprintf('%s (CC%d)', pair_name, cc);
            
            % --- Extract CC Mean Data ---
            e_cc_r = extract_animal_means(group_results(ipair).trial_corr_early, cc);
            e_cc_s = extract_animal_means(group_results(ipair).trial_corr_early_shuff, cc);
            x_cc_r = extract_animal_means(group_results(ipair).trial_corr_post, cc);
            x_cc_s = extract_animal_means(group_results(ipair).trial_corr_post_shuff, cc);
            
            % --- Extract Precession Idx Mean Data ---
            e_pi_r = extract_animal_means(group_results(ipair).trial_precession_early_idx, cc);
            e_pi_s = extract_animal_means(group_results(ipair).trial_precession_early_idx_shuff, cc);
            x_pi_r = extract_animal_means(group_results(ipair).trial_precession_post_idx, cc);
            x_pi_s = extract_animal_means(group_results(ipair).trial_precession_post_idx_shuff, cc);
            
            % --- Extract Trial-by-Trial Data (Using new CC-aware helper) ---
            e_cc_tr_r = extract_epoch_trials_cc(group_results(ipair).trial_corr_early, n_animals, cc);
            e_cc_tr_s = extract_epoch_trials_cc(group_results(ipair).trial_corr_early_shuff, n_animals, cc);
            x_cc_tr_r = extract_epoch_trials_cc(group_results(ipair).trial_corr_post, n_animals, cc);
            x_cc_tr_s = extract_epoch_trials_cc(group_results(ipair).trial_corr_post_shuff, n_animals, cc);
            
            e_pi_tr_r = extract_epoch_trials_cc(group_results(ipair).trial_precession_early_idx, n_animals, cc);
            e_pi_tr_s = extract_epoch_trials_cc(group_results(ipair).trial_precession_early_idx_shuff, n_animals, cc);
            x_pi_tr_r = extract_epoch_trials_cc(group_results(ipair).trial_precession_post_idx, n_animals, cc);
            x_pi_tr_s = extract_epoch_trials_cc(group_results(ipair).trial_precession_post_idx_shuff, n_animals, cc);
            
            % --- Plot CC Bars ---
            figure(fig_cc); nexttile(t_cc);
            plot_real_shuff_bars_with_stats(e_cc_r(mask), e_cc_s(mask), x_cc_r(mask), x_cc_s(mask), title_str, 'Correlation (CC)');
                
            % --- Plot Precession Bars ---
            figure(fig_idx); nexttile(t_idx);
            plot_real_shuff_bars_with_stats(e_pi_r(mask), e_pi_s(mask), x_pi_r(mask), x_pi_s(mask), title_str, 'Precession Index');
                
            % --- Plot CC Trial Curves ---
            figure(fig_cc_crv); nexttile(t_cc_crv);
            plot_continuous_trials_with_shuff(e_cc_tr_r(mask,:), e_cc_tr_s(mask,:), x_cc_tr_r(mask,:), x_cc_tr_s(mask,:), title_str, 'Correlation (CC)');
            
            % --- Plot Precession Trial Curves ---
            figure(fig_pi_crv); nexttile(t_pi_crv);
            plot_continuous_trials_with_shuff(e_pi_tr_r(mask,:), e_pi_tr_s(mask,:), x_pi_tr_r(mask,:), x_pi_tr_s(mask,:), title_str, 'Precession Index');
        end
    end
    
    save_to_svg(fullfile(data_dir, sprintf('Extended_CC_Bars_%s', group_label)));
    save_to_svg(fullfile(data_dir, sprintf('Extended_Prec_Bars_%s', group_label)));
    save_to_svg(fullfile(data_dir, sprintf('Extended_CC_TrialCurves_%s', group_label)));
    save_to_svg(fullfile(data_dir, sprintf('Extended_Prec_TrialCurves_%s', group_label)));
end

%% 5C. COMBINED CONTINUOUS CURVES, ERROR BARS & NETWORKS (REAL VS SHUFFLED)
fprintf('\n--- Generating Combined Curves, Error Bars & Networks (%d CCs) ---\n', num_ccs_analyze);

% Define layout for the network plots
layout_def.names = {'CA1', 'V1', 'DG', 'CA3', 'RSC', 'SUB'};
layout_def.x     = [7.0,  8.0,  3.0,  5.0,  2.0,  2.0];
layout_def.y     = [5.0,  9.0,  4.0,  2.0,  9.0,  6.5];

for g_idx = 1:2
    if g_idx == 1
        mask = is_learner; group_label = 'Learners';
        c_real_e = [0.4 0.6 0.8]; % Light Blue (Naive)
        c_real_x = [0.1 0.3 0.6]; % Dark Blue (Expert)
    else
        mask = ~is_learner; group_label = 'Non-Learners';
        c_real_e = [0.9 0.5 0.3]; % Light Rust (Naive)
        c_real_x = [0.7 0.2 0.1]; % Dark Rust (Expert)
    end
    c_shuff_e = [0.7 0.7 0.7]; % Light Gray
    c_shuff_x = [0.4 0.4 0.4]; % Dark Gray
    
    if sum(mask) == 0, continue; end
    
    for cc = 1:num_ccs_analyze
        % Initialize Correlation Combined Figure
        fig_cc_comb = figure('Name', sprintf('Combined_CC_%d_%s', cc, group_label), 'Color', 'w', 'Position', [100 100 1000 200*n_pairs]);
        t_cc_comb = tiledlayout(n_pairs, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Initialize Precession Combined Figure
        fig_pi_comb = figure('Name', sprintf('Combined_Prec_%d_%s', cc, group_label), 'Color', 'w', 'Position', [150 150 1000 200*n_pairs]);
        t_pi_comb = tiledlayout(n_pairs, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Preallocate arrays for the Network Plot
        net_e_cc = nan(n_pairs, 1); net_x_cc = nan(n_pairs, 1);
        net_e_ifi = nan(n_pairs, 1); net_x_ifi = nan(n_pairs, 1);
        net_e_cc_pval = nan(n_pairs, 1); net_x_cc_pval = nan(n_pairs, 1);
        net_e_ifi_pval = nan(n_pairs, 1); net_x_ifi_pval = nan(n_pairs, 1);
        
        for ipair = 1:n_pairs
            pair_name = group_results(ipair).pair_name;
            title_str = sprintf('%s (CC%d)', pair_name, cc);
            
            % Extract CC Data
            e_cc_r = extract_animal_means(group_results(ipair).trial_corr_early, cc);
            e_cc_s = extract_animal_means(group_results(ipair).trial_corr_early_shuff, cc);
            x_cc_r = extract_animal_means(group_results(ipair).trial_corr_post, cc);
            x_cc_s = extract_animal_means(group_results(ipair).trial_corr_post_shuff, cc);
            
            e_cc_tr_r = extract_epoch_trials_cc(group_results(ipair).trial_corr_early, n_animals, cc);
            e_cc_tr_s = extract_epoch_trials_cc(group_results(ipair).trial_corr_early_shuff, n_animals, cc);
            x_cc_tr_r = extract_epoch_trials_cc(group_results(ipair).trial_corr_post, n_animals, cc);
            x_cc_tr_s = extract_epoch_trials_cc(group_results(ipair).trial_corr_post_shuff, n_animals, cc);
            
            % Extract Precession Data
            e_pi_r = extract_animal_means(group_results(ipair).trial_precession_early_idx, cc);
            e_pi_s = extract_animal_means(group_results(ipair).trial_precession_early_idx_shuff, cc);
            x_pi_r = extract_animal_means(group_results(ipair).trial_precession_post_idx, cc);
            x_pi_s = extract_animal_means(group_results(ipair).trial_precession_post_idx_shuff, cc);
            
            e_pi_tr_r = extract_epoch_trials_cc(group_results(ipair).trial_precession_early_idx, n_animals, cc);
            e_pi_tr_s = extract_epoch_trials_cc(group_results(ipair).trial_precession_early_idx_shuff, n_animals, cc);
            x_pi_tr_r = extract_epoch_trials_cc(group_results(ipair).trial_precession_post_idx, n_animals, cc);
            x_pi_tr_s = extract_epoch_trials_cc(group_results(ipair).trial_precession_post_idx_shuff, n_animals, cc);
            
            % Calculate Valid N's for the current pair and group
            n_e = sum(~isnan(e_cc_r(mask)) & ~isnan(e_cc_s(mask)));
            n_x = sum(~isnan(x_cc_r(mask)) & ~isnan(x_cc_s(mask)));
            annotated_title = sprintf('%s (n_{naive}=%d, n_{expert}=%d)', title_str, n_e, n_x);
            
            % Compute CC Paired T-Tests (Real vs Shuff) for edge existence
            if n_e > 2, [~, p_cc_e] = ttest(e_cc_r(mask), e_cc_s(mask)); else, p_cc_e = nan; end
            if n_x > 2, [~, p_cc_x] = ttest(x_cc_r(mask), x_cc_s(mask)); else, p_cc_x = nan; end
            
            % Compute Precession Paired T-Tests (Real vs Shuff) for arrows/asterisks
            n_e_pi = sum(~isnan(e_pi_r(mask)) & ~isnan(e_pi_s(mask)));
            n_x_pi = sum(~isnan(x_pi_r(mask)) & ~isnan(x_pi_s(mask)));
            if n_e_pi > 2, [~, p_ifi_e] = ttest(e_pi_r(mask), e_pi_s(mask)); else, p_ifi_e = nan; end
            if n_x_pi > 2, [~, p_ifi_x] = ttest(x_pi_r(mask), x_pi_s(mask)); else, p_ifi_x = nan; end
            
            % Store Network Metrics
            net_e_cc(ipair) = mean(e_cc_r(mask), 'omitnan');
            net_x_cc(ipair) = mean(x_cc_r(mask), 'omitnan');
            net_e_ifi(ipair) = mean(e_pi_r(mask), 'omitnan');
            net_x_ifi(ipair) = mean(x_pi_r(mask), 'omitnan');
            net_e_cc_pval(ipair) = p_cc_e;
            net_x_cc_pval(ipair) = p_cc_x;
            net_e_ifi_pval(ipair) = p_ifi_e;
            net_x_ifi_pval(ipair) = p_ifi_x;
            
            % --- Render CC Combined ---
            figure(fig_cc_comb); 
            nexttile(t_cc_comb);
            plot_continuous_trials_with_shuff(e_cc_tr_r(mask,:), e_cc_tr_s(mask,:), x_cc_tr_r(mask,:), x_cc_tr_s(mask,:), sprintf('%s Trial Curve', annotated_title), 'Correlation (CC)', c_real_e, c_real_x, c_shuff_e, c_shuff_x);
            if ipair ~= n_pairs, xlabel(''); end
            
            nexttile(t_cc_comb);
            plot_real_shuff_bars_with_stats(e_cc_r(mask), e_cc_s(mask), x_cc_r(mask), x_cc_s(mask), sprintf('%s Averages', annotated_title), 'Correlation (CC)', c_real_e, c_real_x, c_shuff_e, c_shuff_x);
            
            % --- Render Precession Combined ---
            figure(fig_pi_comb); 
            nexttile(t_pi_comb);
            plot_continuous_trials_with_shuff(e_pi_tr_r(mask,:), e_pi_tr_s(mask,:), x_pi_tr_r(mask,:), x_pi_tr_s(mask,:), sprintf('%s Trial Curve', annotated_title), 'Precession Index', c_real_e, c_real_x, c_shuff_e, c_shuff_x);
            if ipair ~= n_pairs, xlabel(''); end
            
            nexttile(t_pi_comb);
            plot_real_shuff_bars_with_stats(e_pi_r(mask), e_pi_s(mask), x_pi_r(mask), x_pi_s(mask), sprintf('%s Averages', annotated_title), 'Precession Index', c_real_e, c_real_x, c_shuff_e, c_shuff_x);
        end
        
        % Add Legends to Combined Figures
        L = gobjects(4,1);
        L(1) = plot(nan, nan, 'Color', c_real_e, 'LineWidth', 2); hold on;
        L(2) = plot(nan, nan, 'Color', c_real_x, 'LineWidth', 2);
        L(3) = plot(nan, nan, 'Color', c_shuff_e, 'LineWidth', 2, 'LineStyle', '--');
        L(4) = plot(nan, nan, 'Color', c_shuff_x, 'LineWidth', 2, 'LineStyle', '--');
        
        figure(fig_cc_comb);
        lgd1 = legend(L, {'Real Naive', 'Real Expert', 'Shuff Naive', 'Shuff Expert'}, 'Orientation', 'horizontal', 'NumColumns', 4);
        lgd1.Layout.Tile = 'north';
        save_to_svg(fullfile(data_dir, sprintf('Combined_CC_%d_%s', cc, group_label)));
        
        figure(fig_pi_comb);
        lgd2 = legend(L, {'Real Naive', 'Real Expert', 'Shuff Naive', 'Shuff Expert'}, 'Orientation', 'horizontal', 'NumColumns', 4);
        lgd2.Layout.Tile = 'north';
        save_to_svg(fullfile(data_dir, sprintf('Combined_Prec_%d_%s', cc, group_label)));
        
        
        % --- Render Filtered Network for this CC and Group ---
        fig_net = figure('Name', sprintf('Network_CC%d_%s', cc, group_label), 'Color', 'w', 'Position', [200 200 1000 500]);
        t_net = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
        
        global_max_cc = max([net_e_cc; net_x_cc], [], 'omitnan');
        global_max_ifi = max(abs([net_e_ifi; net_x_ifi]), [], 'omitnan');
        if isempty(global_max_cc) || global_max_cc == 0 || isnan(global_max_cc), global_max_cc = 1; end
        if isempty(global_max_ifi) || global_max_ifi == 0 || isnan(global_max_ifi), global_max_ifi = 1; end
        
        nexttile(t_net);
        plot_dual_network(area_pairs_to_analyze, net_e_cc, net_e_ifi, net_e_cc_pval, net_e_ifi_pval, sprintf('%s: Naive (CC%d)', group_label, cc), layout_def, global_max_cc, global_max_ifi);
        
        nexttile(t_net);
        plot_dual_network(area_pairs_to_analyze, net_x_cc, net_x_ifi, net_x_cc_pval, net_x_ifi_pval, sprintf('%s: Expert (CC%d)', group_label, cc), layout_def, global_max_cc, global_max_ifi);
        
        save_to_svg(fullfile(data_dir, sprintf('Network_Filtered_CC%d_%s', cc, group_label)));
    end
end

%% 6. LOCAL HELPERS

function mat = extract_epoch_trials(cell_data, n_animals)
    mat = nan(n_animals, 10);
    for i = 1:n_animals
        d = cell_data{i};
        if ~isempty(d) && size(d, 2) == 10
            mat(i, :) = d(1, :); 
        end
    end
end

function plot_epoch_line(x_range, data_mat, color, linestyle)
    valid_mask = sum(~isnan(data_mat), 2) > 0;
    valid_data = data_mat(valid_mask, :);
    if isempty(valid_data), return; end
    
    mu = mean(valid_data, 1, 'omitnan');
    se = std(valid_data, 0, 1, 'omitnan') ./ sqrt(size(valid_data, 1));
    
    if strcmp(linestyle, '--')
        plot(x_range, mu, 'Color', color, 'LineStyle', linestyle, 'LineWidth', 2);
        plot(x_range, mu+se, 'Color', color, 'LineStyle', ':', 'LineWidth', 0.5);
        plot(x_range, mu-se, 'Color', color, 'LineStyle', ':', 'LineWidth', 0.5);
    else
        shadedErrorBar(x_range, mu, se, 'lineProps', {'Color', color, 'LineStyle', linestyle, 'LineWidth', 2});
    end
end

function [r_shifts, p_idx] = calc_precession(D1, D2, max_shift, num_ccs, nc1, nc2)
    shifts = -max_shift : max_shift;
    r_shifts = zeros(num_ccs, length(shifts));
    for is = 1:length(shifts)
        shift = shifts(is);
        if shift < 0      
             s1 = D1(:, 1:end+shift, :); s2 = D2(:, 1-shift:end, :);
        elseif shift > 0  
             s1 = D1(:, 1+shift:end, :); s2 = D2(:, 1:end-shift, :);
        else
             s1 = D1; s2 = D2;
        end
        sx = reshape(s1, nc1, []); sy = reshape(s2, nc2, []);
        if size(sx,2) > max(nc1, nc2) + 2
             [~,~,rs] = canoncorr(sx', sy');
             r_shifts(1:min(length(rs), num_ccs), is) = rs(1:min(length(rs), num_ccs));
        end
    end
    idx_neg = 1:max_shift; idx_pos = max_shift+2 : length(shifts);
    p_idx = nan(num_ccs, 1);
    for cc = 1:num_ccs
        neg_val = mean(r_shifts(cc, idx_neg), 'omitnan');
        pos_val = mean(r_shifts(cc, idx_pos), 'omitnan');
        if (neg_val + pos_val) > 0.001
            p_idx(cc) = (neg_val - pos_val) / (neg_val + pos_val);
        else
            p_idx(cc) = nan;
            % Inject tracking here
            fprintf('  [NaN Tracker] Precession NaN assigned: Correlation sum too low (%.5f)\n', neg_val + pos_val);
        end
    end
end

function raw_means = extract_animal_means(cell_data, cc_idx)
    n_animals = length(cell_data); raw_means = nan(n_animals, 1);
    for i = 1:n_animals, if ~isempty(cell_data{i}), raw_means(i) = mean(cell_data{i}(cc_idx, :), 2, 'omitnan'); end; end
end

function save_to_svg(fig_name)
    fig = gcf; 
    fprintf('Saving %s.svg...\n', fig_name);
    try 
        % Modern export: completely preserves screen aspect ratios and layout
        exportgraphics(fig, [fig_name '.svg'], 'ContentType', 'vector', 'BackgroundColor', 'w'); 
    catch 
        % Fallback for older MATLAB versions: locks the paper size to screen size
        set(fig, 'PaperPositionMode', 'auto'); 
        print(fig, '-dpng', [fig_name '.png']); 
    end
end

function ok = check_dims(X, Y, min_samples)
    if isempty(X) || isempty(Y), ok = false; return; end
    n_samples = size(X, 2); if any(isnan(X(:))) || any(isnan(Y(:))), ok = false; return; end
    ok = n_samples > min_samples; 
end

function plot_dual_network(pairs, cc_vals, ifi_vals, pval_cc, pval_ifi, titleStr, layout, max_cc, max_ifi)
    regions = unique(pairs(:)); 
    sources = []; targets = []; 
    weights_cc = []; weights_ifi = []; weights_pval_ifi = [];
    
    n_edges = size(pairs, 1);
    for i = 1:n_edges
        u = find(strcmp(regions, pairs{i, 1})); 
        v = find(strcmp(regions, pairs{i, 2}));
        
        cc = cc_vals(i);
        ifi = ifi_vals(i);
        p_cc = pval_cc(i);
        p_ifi = pval_ifi(i);
        
        % STRUCTURAL FILTER: Keep edge only if CC is significant vs Shuffled
        if isnan(cc) || isnan(p_cc) || p_cc >= 0.05
            continue; 
        end
        
        % Direction defined strictly by IFI
        if ifi >= 0
            s = u; t = v; 
        else
            s = v; t = u; 
        end
    
        sources = [sources; s]; 
        targets = [targets; t]; 
        weights_cc = [weights_cc; cc];
        weights_ifi = [weights_ifi; abs(ifi)];
        weights_pval_ifi = [weights_pval_ifi; p_ifi];
    end
    
    NodeTable = table(regions, 'VariableNames', {'Name'});
    x_coords = zeros(1, numel(regions)); y_coords = zeros(1, numel(regions));
    
    for i = 1:numel(regions)
        idx = find(strcmp(layout.names, regions{i}));
        if ~isempty(idx), x_coords(i) = layout.x(idx); y_coords(i) = layout.y(idx); end
    end
    
    if isempty(sources)
        G_empty = graph([], [], [], NodeTable);
        p_empty = plot(G_empty, 'XData', x_coords, 'YData', y_coords);
        p_empty.NodeColor = [0.85 0.85 0.85]; p_empty.MarkerSize = 30; p_empty.NodeFontSize = 12; p_empty.NodeFontWeight = 'bold';
        p_empty.EdgeColor = 'none';
        axis off; xlim([0 10]); ylim([0 10]); 
        title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
        return;
    end
    
    % --- 1. Plot UNDIRECTED BASE GRAPH (For all significant CCs) ---
    G_un = graph(sources, targets, weights_cc, NodeTable);
    
    p_un = plot(G_un, 'XData', x_coords, 'YData', y_coords);
    p_un.NodeColor = [0.85 0.85 0.85]; p_un.MarkerSize = 30; p_un.NodeFontSize = 12; p_un.NodeFontWeight = 'bold';
    p_un.EdgeColor = [0.85 0.85 0.85]; 
    
    norm_cc_un = G_un.Edges.Weight / max_cc;
    p_un.LineWidth = 1 + (norm_cc_un * 8); 
    
    hold on;
    
    % --- 2. Plot DIRECTED OVERLAY GRAPH (Only if Precession is significant) ---
    sig_idx = weights_pval_ifi < 0.05;
    
    if any(sig_idx)
        sig_sources = sources(sig_idx);
        sig_targets = targets(sig_idx);
        sig_cc      = weights_cc(sig_idx);
        sig_ifi     = weights_ifi(sig_idx);
        
        EdgeTable_dir = table([sig_sources, sig_targets], sig_cc, sig_ifi, ...
            'VariableNames', {'EndNodes', 'Weight', 'IFI'});
            
        G_dir = digraph(EdgeTable_dir, NodeTable);
        
        p_dir = plot(G_dir, 'XData', x_coords, 'YData', y_coords);
        p_dir.NodeColor = 'none'; 
        p_dir.EdgeLabel = {};     
        p_dir.NodeLabel = {};     
        
        norm_cc_dir = G_dir.Edges.Weight / max_cc;
        p_dir.LineWidth = 1 + (norm_cc_dir * 8); 
        
        norm_ifi_dir = min(G_dir.Edges.IFI / max_ifi, 1); 
        base_color = [0.7 0.7 0.7];   
        target_color = [0.1 0.4 0.6]; 
        
        edge_colors_dir = zeros(length(norm_ifi_dir), 3);
        for i = 1:length(norm_ifi_dir)
            edge_colors_dir(i, :) = base_color + (target_color - base_color) * norm_ifi_dir(i);
        end
        
        p_dir.EdgeColor = edge_colors_dir; 
        p_dir.ArrowSize = 15;
        
        % Draw asterisks aligned to directed edges
        for i = 1:numedges(G_dir)
            s_name = G_dir.Edges.EndNodes(i, 1);
            t_name = G_dir.Edges.EndNodes(i, 2);
            
            s_idx = findnode(G_dir, s_name);
            t_idx = findnode(G_dir, t_name);
            
            x1 = x_coords(s_idx); y1 = y_coords(s_idx);
            x2 = x_coords(t_idx); y2 = y_coords(t_idx);
            
            x_mid = (x1 + x2) / 2;
            y_mid = (y1 + y2) / 2;
            
            dx = x2 - x1; dy = y2 - y1;
            len = sqrt(dx^2 + dy^2);
            nx = -dy / len * 0.15; 
            ny = dx / len * 0.15;  
            
            text(x_mid + nx, y_mid + ny, '*', 'FontSize', 22, 'Color', edge_colors_dir(i, :), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
        end
    end
    
    axis off; 
    xlim([0 10]); ylim([0 10]);
    title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
end

function plot_grouped_bars_with_rmanova(e_vals, p_vals, x_vals, is_learner, pair_name, y_label, test_zero)
    L_e = e_vals(is_learner); L_p = p_vals(is_learner); L_x = x_vals(is_learner);
    NL_e = e_vals(~is_learner); NL_p = p_vals(~is_learner); NL_x = x_vals(~is_learner);
    
    means_L = [mean(L_e, 'omitnan'), mean(L_p, 'omitnan'), mean(L_x, 'omitnan')];
    sems_L  = [std(L_e, 'omitnan')/sqrt(sum(~isnan(L_e))), std(L_p, 'omitnan')/sqrt(sum(~isnan(L_p))), std(L_x, 'omitnan')/sqrt(sum(~isnan(L_x)))];
    
    means_NL = [mean(NL_e, 'omitnan'), mean(NL_p, 'omitnan'), mean(NL_x, 'omitnan')];
    sems_NL  = [std(NL_e, 'omitnan')/sqrt(sum(~isnan(NL_e))), std(NL_p, 'omitnan')/sqrt(sum(~isnan(NL_p))), std(NL_x, 'omitnan')/sqrt(sum(~isnan(NL_x)))];
    
    hold on;
    b = bar(1:3, [means_L; means_NL]', 'grouped');
    b(1).FaceColor = [0.85 0.33 0.1]; 
    b(2).FaceColor = [0.6 0.6 0.6];   
    
    x_offset = [b(1).XEndPoints; b(2).XEndPoints]';
    errorbar(x_offset(:,1), means_L, sems_L, 'k.', 'LineWidth', 1.5);
    errorbar(x_offset(:,2), means_NL, sems_NL, 'k.', 'LineWidth', 1.5);
    
    yline(0, '-k', 'LineWidth', 0.5);
    xticks(1:3); xticklabels({'Early', 'Pre', 'Post'});
    title(pair_name); ylabel(y_label);
    
    yl = ylim; offset = (yl(2) - yl(1)) * 0.05; 
    min_y_needed = yl(1); max_y_needed = yl(2);
    
    % --- 1-Sample T-tests against 0 ---
    if test_zero
        p_L = nan(1,3); p_NL = nan(1,3);
        if sum(~isnan(L_e)) > 2, [~, p_L(1)] = ttest(L_e, 0); end
        if sum(~isnan(L_p)) > 2, [~, p_L(2)] = ttest(L_p, 0); end
        if sum(~isnan(L_x)) > 2, [~, p_L(3)] = ttest(L_x, 0); end
        
        if sum(~isnan(NL_e)) > 2, [~, p_NL(1)] = ttest(NL_e, 0); end
        if sum(~isnan(NL_p)) > 2, [~, p_NL(2)] = ttest(NL_p, 0); end
        if sum(~isnan(NL_x)) > 2, [~, p_NL(3)] = ttest(NL_x, 0); end
        
        for idx = 1:3
            % Learners
            if p_L(idx) < 0.05
                if means_L(idx) >= 0
                    y_pos = means_L(idx) + sems_L(idx) + offset; vert_align = 'bottom';
                else
                    y_pos = means_L(idx) - sems_L(idx) - offset; vert_align = 'top'; 
                end
                if y_pos > max_y_needed, max_y_needed = y_pos + offset; end
                if y_pos < min_y_needed, min_y_needed = y_pos - offset; end
                text(x_offset(idx, 1), y_pos, '*', 'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', vert_align, 'Color', [0.85 0.33 0.1]);
            end
            % Non-Learners
            if p_NL(idx) < 0.05
                if means_NL(idx) >= 0
                    y_pos = means_NL(idx) + sems_NL(idx) + offset; vert_align = 'bottom';
                else
                    y_pos = means_NL(idx) - sems_NL(idx) - offset; vert_align = 'top'; 
                end
                if y_pos > max_y_needed, max_y_needed = y_pos + offset; end
                if y_pos < min_y_needed, min_y_needed = y_pos - offset; end
                text(x_offset(idx, 2), y_pos, '*', 'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', vert_align, 'Color', [0.4 0.4 0.4]);
            end
        end
    end
    
    ylim([min(yl(1), min_y_needed), max(yl(2), max_y_needed)]); 
    
    % --- Mixed Repeated Measures ANOVA ---
    valid_mask = ~isnan(e_vals) & ~isnan(p_vals) & ~isnan(x_vals);
    n_complete = sum(valid_mask);
    
    if n_complete > max(2, length(unique(is_learner(valid_mask))))
        T = table(e_vals(valid_mask), p_vals(valid_mask), x_vals(valid_mask), ...
            categorical(is_learner(valid_mask), [1 0], {'Learner', 'NonLearner'}), ...
            'VariableNames', {'Early', 'Pre', 'Post', 'LearnerGroup'});
        
        Meas = table([1 2 3]', 'VariableNames', {'Epoch'});
        
        try
            rm = fitrm(T, 'Early-Post ~ LearnerGroup', 'WithinDesign', Meas);
            
            % Within-Subjects
            ranovatbl = ranova(rm);
            p_epoch = ranovatbl.pValue(1);
            p_int   = ranovatbl.pValue(2);
            
            % Between-Subjects
            betweentbl = anova(rm);
            p_group = betweentbl.pValue(2);
            
            % Print to Console
            disp('--- Within-Subjects (Epoch & Interaction) ---');
            disp(ranovatbl);
            disp('--- Between-Subjects (LearnerGroup) ---');
            disp(betweentbl);
            
            % Build Significance Text Box
            sig_text = {sprintf('rmANOVA (n=%d):', n_complete)};
            has_sig = false;
            if p_epoch < 0.05, sig_text{end+1} = sprintf('Epoch: p=%.3f', p_epoch); has_sig = true; end
            if p_group < 0.05, sig_text{end+1} = sprintf('Group: p=%.3f', p_group); has_sig = true; end
            if p_int < 0.05,   sig_text{end+1} = sprintf('Ep x Grp: p=%.3f', p_int); has_sig = true; end
            if ~has_sig, sig_text{end+1} = 'n.s.'; end
            
            curr_yl = ylim; curr_xl = xlim;
            y_range = curr_yl(2) - curr_yl(1);
            ylim([curr_yl(1), curr_yl(2) + y_range*0.25]); 
            new_yl = ylim;
            
            text(curr_xl(1) + 0.03*(curr_xl(2)-curr_xl(1)), new_yl(2) - 0.03*(new_yl(2)-new_yl(1)), ...
                 strjoin(sig_text, '\n'), 'VerticalAlignment', 'top', ...
                 'FontSize', 9, 'EdgeColor', 'k', 'BackgroundColor', 'w', 'Margin', 3);
                 
            % Run Post-hocs
            if p_int < 0.05
                fprintf('  -> Interaction Significant. Post-hoc Epoch across Groups:\n');
                mc = multcompare(rm, 'LearnerGroup', 'By', 'Epoch');
                disp(mc);
            else
                if p_epoch < 0.05
                    fprintf('  -> Main Effect of Epoch Significant. Post-hoc:\n');
                    mc = multcompare(rm, 'Epoch');
                    disp(mc);
                end
                if p_group < 0.05
                    fprintf('  -> Main Effect of Group Significant. Post-hoc:\n');
                    mc = multcompare(rm, 'LearnerGroup');
                    disp(mc);
                end
            end
            
        catch ME
            fprintf('rmANOVA failed: %s\n', ME.message);
        end
    else
        fprintf('Not enough complete cases (n=%d) for strict rmANOVA.\n', n_complete);
    end
end

function plot_real_shuff_bars_with_stats(e_r, e_s, x_r, x_s, title_str, y_label, c_re, c_rx, c_se, c_sx)
    % Averages and SEMs
    means = [mean(e_r, 'omitnan'), mean(e_s, 'omitnan'); mean(x_r, 'omitnan'), mean(x_s, 'omitnan')];
    sems  = [std(e_r, 'omitnan')/sqrt(sum(~isnan(e_r))), std(e_s, 'omitnan')/sqrt(sum(~isnan(e_s))); ...
             std(x_r, 'omitnan')/sqrt(sum(~isnan(x_r))), std(x_s, 'omitnan')/sqrt(sum(~isnan(x_s)))];
             
    hold on;
    b = bar(1:2, means, 'grouped');
    
    % Apply Flat Shading to synchronize bar colors with curve colors
    b(1).FaceColor = 'flat';
    b(1).CData = [c_re; c_rx]; 
    
    b(2).FaceColor = 'flat';
    b(2).CData = [c_se; c_sx];
    
    % Add error bars manually
    x_offset = [b(1).XEndPoints; b(2).XEndPoints]';
    errorbar(x_offset(:,1), means(:,1), sems(:,1), 'k.', 'LineWidth', 1.5);
    errorbar(x_offset(:,2), means(:,2), sems(:,2), 'k.', 'LineWidth', 1.5);
    
    yline(0, '-k', 'LineWidth', 0.5);
    xticks(1:2); xticklabels({'Naive (Early)', 'Expert (Post)'});
    ylabel(y_label); title(title_str);
    
    % --- Statistics (Paired T-Tests) ---
    valid_e = ~isnan(e_r) & ~isnan(e_s);
    valid_x = ~isnan(x_r) & ~isnan(x_s);
    valid_ep = ~isnan(e_r) & ~isnan(x_r); % Real Early vs Real Post
    
    p_e_shuff = nan; p_x_shuff = nan; p_epoch = nan;
    if sum(valid_e) > 2, [~, p_e_shuff] = ttest(e_r(valid_e), e_s(valid_e)); end
    if sum(valid_x) > 2, [~, p_x_shuff] = ttest(x_r(valid_x), x_s(valid_x)); end
    if sum(valid_ep) > 2, [~, p_epoch] = ttest(x_r(valid_ep), e_r(valid_ep)); end
    
    yl = ylim; 
    
    % Only calculate offset if bars have actual height to prevent scaling bugs
    if yl(2) == 0 && yl(1) == 0, yl = [-1 1]; end 
    offset = (yl(2) - yl(1)) * 0.08;
    max_y = max([means(1,:) + sems(1,:), means(2,:) + sems(2,:)]);
    
    % Add connecting lines/stars for Real vs Shuff within epochs
    if p_e_shuff < 0.05
        plot(x_offset(1,:), [1 1]*(max_y+offset), '-k');
        text(mean(x_offset(1,:)), max_y+offset*1.4, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
    end
    if p_x_shuff < 0.05
        plot(x_offset(2,:), [1 1]*(max_y+offset), '-k');
        text(mean(x_offset(2,:)), max_y+offset*1.4, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
    end
    
    % Add text overlay for Real Epoch shift (Early vs Post)
    text_y = max_y + offset * 3;
    if p_epoch < 0.05
        text(1.5, text_y, sprintf('Early vs Post: p = %.3f', p_epoch), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    else
        text(1.5, text_y, 'Early vs Post: n.s.', ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    ylim([min(yl(1), 0 - offset), text_y + offset*2]);
    % Standardize tick format to ensure identical text bounding boxes across figures
    ytickformat('%.2f');
    
    % Force the axes to update their margin properties immediately
    drawnow;
end

function plot_continuous_trials_with_shuff(e_r, e_s, x_r, x_s, title_str, y_label, c_re, c_rx, c_se, c_sx)
    hold on;
    % Plot Shuffled Data (Grays, Dashed)
    plot_epoch_line(1:10, e_s, c_se, '--'); 
    plot_epoch_line(11:20, x_s, c_sx, '--'); 
    
    % Plot Real Data (Colors, Solid)
    plot_epoch_line(1:10, e_r, c_re, '-');  
    plot_epoch_line(11:20, x_r, c_rx, '-');  
    
    xline(10.5, 'k:'); 
    yline(0, '-k', 'LineWidth', 0.5);
    
    xticks([5, 15]); 
    xticklabels({'Naive', 'Expert'});
    xlim([1 20]);
    title(title_str); ylabel(y_label);

    % Standardize tick format to ensure identical text bounding boxes across figures
    ytickformat('%.2f');
    
    % Force the axes to update their margin properties immediately
    drawnow;
end

function mat = extract_epoch_trials_cc(cell_data, n_animals, cc_idx)
    % Extracts trial-by-trial data dynamically based on the requested CC
    mat = nan(n_animals, 10);
    for i = 1:n_animals
        d = cell_data{i};
        if ~isempty(d) && size(d, 1) >= cc_idx && size(d, 2) == 10
            mat(i, :) = d(cc_idx, :); 
        end
    end
end

