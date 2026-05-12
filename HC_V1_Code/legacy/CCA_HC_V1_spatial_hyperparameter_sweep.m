%% CCA_HC_V1_Hyperparameter_Sweep.m
% Spatial CCA Analysis - Hyperparameter Sweep
%
% DESCRIPTION:
% Sweeps across PCA dimensionality reduction strategies and Precession Lag
% limits to assess the robustness of the CCA spatial/trial-wise pipeline.
%
% SWEEP PARAMETERS:
% - PCA Method: Variance (50%, 60%, 75%, 90%, 95%) AND Fixed PCs (3, 4, 5)
% - Max Shift Lags: 2, 3, 4, 5 bins
%
% OUTPUTS:
% - Generates a separate .mat file and set of SVGs for each configuration.
% - Saves all outputs to a 'Sweep_Results' subfolder.

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
save_dir = fullfile(data_dir, sprintf('Sweep_Results_%s', current_date));
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

% --- Generate Sweep Configurations ---
var_thresholds = [50, 60, 75, 90, 95];
fixed_pcs = [3, 4, 5];
shift_bins = [2, 3, 4, 5];

configs = [];
idx = 1;
% 1. Variance configurations
for v = var_thresholds
    for s = shift_bins
        configs(idx).method = 'variance';
        configs(idx).val = v;
        configs(idx).shift = s;
        configs(idx).name = sprintf('Var%d_Shift%d', v, s);
        idx = idx + 1;
    end
end
% 2. Fixed configurations
for f = fixed_pcs
    for s = shift_bins
        configs(idx).method = 'fixed';
        configs(idx).val = f;
        configs(idx).shift = s;
        configs(idx).name = sprintf('Fix%d_Shift%d', f, s);
        idx = idx + 1;
    end
end
num_configs = length(configs);
fprintf('Initialized Sweep with %d configurations.\n', num_configs);

% --- Static Analysis Parameters ---
num_ccs_analyze = 2;         
n_trials_window = -3:3;      
n_bins_window = -3:3;        
n_shuffles = 50;             
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

%% 2. LOAD LEARNING POINTS 
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
    learning_points = [];
end

% Pre-load file list
file_list = dir(fullfile(data_dir, file_pattern));
n_animals = length(file_list);

%% 3. MASTER SWEEP LOOP
for icfg = 1:num_configs
    current_cfg = configs(icfg);
    fprintf('\n=================================================================\n');
    fprintf('--- Running Configuration %d/%d: %s ---\n', icfg, num_configs, current_cfg.name);
    fprintf('=================================================================\n');
    
    % Initialize Results Structure for this config
    group_results = struct('pair_name', cell(n_pairs, 1), ...
                           'all_bins_corr', cell(n_pairs, 1), ...
                           'all_bins_corr_shuff', cell(n_pairs, 1), ...
                           'all_bins_precession', cell(n_pairs, 1), ...
                           'all_bins_precession_shuff', cell(n_pairs, 1), ...
                           ...
                           'trial_corr_early', cell(n_pairs, 1), ...
                           'trial_corr_pre', cell(n_pairs, 1), ...
                           'trial_corr_post', cell(n_pairs, 1), ...
                           'trial_corr_early_shuff', cell(n_pairs, 1), ...
                           'trial_corr_pre_shuff', cell(n_pairs, 1), ...
                           'trial_corr_post_shuff', cell(n_pairs, 1), ...
                           ...
                           'trial_precession_early', cell(n_pairs, 1), ... 
                           'trial_precession_pre', cell(n_pairs, 1), ...   
                           'trial_precession_post', cell(n_pairs, 1), ...
                           'trial_precession_early_shuff', cell(n_pairs, 1), ...
                           'trial_precession_pre_shuff', cell(n_pairs, 1), ...
                           'trial_precession_post_shuff', cell(n_pairs, 1));
                       
    for ipair = 1:n_pairs
        group_results(ipair).pair_name = sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});
    end

    % --- ANIMAL LOOP ---
    for ianimal = 1:n_animals
        filename = file_list(ianimal).name;
        fullpath = fullfile(data_dir, filename);
        
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
                 raw_spatial = D.analysis_spatial.firing.cued.freq;
                 animal_data = permute(raw_spatial, [1, 3, 2]); 
                 [~, curr_n_bins, num_trials] = size(animal_data);
                 if curr_n_bins ~= n_bins, continue; end
            else, continue; end
            
            animal_areas = unique(units.regions_label);
            AreaActivity = struct();
            
            % --- PCA Reduction using Config ---
            for ia = 1:length(animal_areas)
                area_name = animal_areas{ia};
                u_idx = units.idx(strcmp(units.regions_label, area_name), :);
                if sum(u_idx) < min_units_per_region, continue; end
                
                area_dat = animal_data(u_idx, :, :);
                reshaped_dat = reshape(area_dat, sum(u_idx), [])'; 
                
                if size(reshaped_dat, 2) >= 2 
                    [~, scores, ~, ~, explained] = pca(reshaped_dat);
                    if strcmp(current_cfg.method, 'variance')
                        cum_var = cumsum(explained);
                        n_comps = find(cum_var >= current_cfg.val, 1);
                        if isempty(n_comps), n_comps = size(scores, 2); end
                    else
                        n_comps = current_cfg.val; 
                    end
                    
                    n_comps = max(n_comps, num_ccs_analyze); 
                    n_comps = min(n_comps, size(scores, 2));
                    
                    AreaActivity.(area_name).data = reshape(scores(:, 1:n_comps)', n_comps, n_bins, num_trials);
                    AreaActivity.(area_name).n_comps = n_comps;
                end
            end
            
            % --- CCA Loop ---
            for ipair = 1:n_pairs
                a1 = area_pairs_to_analyze{ipair, 1};
                a2 = area_pairs_to_analyze{ipair, 2};
                
                if ~isfield(AreaActivity, a1) || ~isfield(AreaActivity, a2), continue; end
                
                d1 = AreaActivity.(a1).data; d2 = AreaActivity.(a2).data;
                nc1 = AreaActivity.(a1).n_comps; nc2 = AreaActivity.(a2).n_comps;
                
                % Preallocate arrays
                cca_tr = nan(num_ccs_analyze, num_trials);
                cca_tr_shuff = nan(num_ccs_analyze, num_trials); 
                prec_tr = nan(num_ccs_analyze, num_trials);
                prec_tr_shuff = nan(num_ccs_analyze, num_trials); 
                
                cca_bin = nan(num_ccs_analyze, n_bins);
                cca_bin_shuff = nan(num_ccs_analyze, n_bins); 
                prec_bin = nan(num_ccs_analyze, n_bins);
                prec_bin_shuff = nan(num_ccs_analyze, n_bins);
                
                % 1. TRIAL-WISE
                for t = 1:num_trials
                    win = t + n_trials_window; win = win(win>=1 & win<=num_trials);
                    if isempty(win), continue; end
                    
                    D1_local = d1(:, :, win); D2_local = d2(:, :, win);
                    x = reshape(D1_local, nc1, []); y = reshape(D2_local, nc2, []);
                    
                    if check_dims(x, y, max(nc1, nc2) + 5)
                        [~,~,r] = canoncorr(x', y');
                        cca_tr(1:min(length(r), num_ccs_analyze), t) = r(1:min(length(r), num_ccs_analyze));
                        prec_tr(:, t) = calc_precession(D1_local, D2_local, current_cfg.shift, num_ccs_analyze, nc1, nc2);
                        
                        r_sh_iter = nan(num_ccs_analyze, n_shuffles);
                        p_sh_iter = nan(num_ccs_analyze, n_shuffles);
                        n_bins_local = size(D1_local, 2);
                        
                        parfor ish = 1:n_shuffles
                             perm_idx = randperm(n_bins_local);
                             D1_s = D1_local(:, perm_idx, :); x_s = reshape(D1_s, nc1, []);
                             [~,~,rs] = canoncorr(x_s', y');
                             r_sh_iter(:, ish) = rs(1:min(length(rs), num_ccs_analyze));
                             p_sh_iter(:, ish) = calc_precession(D1_s, D2_local, current_cfg.shift, num_ccs_analyze, nc1, nc2);
                        end
                        cca_tr_shuff(:, t) = mean(r_sh_iter, 2, 'omitnan');
                        prec_tr_shuff(:, t) = mean(p_sh_iter, 2, 'omitnan');
                    end
                end
                
                % 2. BIN-WISE
                for b = 1:n_bins
                    win = b + n_bins_window; win = win(win>=1 & win<=n_bins);
                    if isempty(win), continue; end
                    
                    D1_local = d1(:, win, :); D2_local = d2(:, win, :);
                    x = reshape(D1_local, nc1, []); y = reshape(D2_local, nc2, []);
                    
                    if check_dims(x, y, max(nc1, nc2) + 5)
                        [~,~,r] = canoncorr(x', y');
                        cca_bin(1:min(length(r), num_ccs_analyze), b) = r(1:min(length(r), num_ccs_analyze));
                        prec_bin(:, b) = calc_precession(D1_local, D2_local, current_cfg.shift, num_ccs_analyze, nc1, nc2);
                        
                        r_sh_iter = nan(num_ccs_analyze, n_shuffles);
                        p_sh_iter = nan(num_ccs_analyze, n_shuffles);
                        n_tr_local = size(D1_local, 3);
                        
                        parfor ish = 1:n_shuffles
                             perm_idx = randperm(n_tr_local);
                             D1_s = D1_local(:, :, perm_idx); x_s = reshape(D1_s, nc1, []);
                             [~,~,rs] = canoncorr(x_s', y');
                             r_sh_iter(:, ish) = rs(1:min(length(rs), num_ccs_analyze));
                             p_sh_iter(:, ish) = calc_precession(D1_s, D2_local, current_cfg.shift, num_ccs_analyze, nc1, nc2);
                        end
                        cca_bin_shuff(:, b) = mean(r_sh_iter, 2, 'omitnan');
                        prec_bin_shuff(:, b) = mean(p_sh_iter, 2, 'omitnan');
                    end
                end
                
                % STORE
                group_results(ipair).all_bins_corr{ianimal} = cca_bin;
                group_results(ipair).all_bins_corr_shuff{ianimal} = cca_bin_shuff;
                group_results(ipair).all_bins_precession{ianimal} = prec_bin;
                group_results(ipair).all_bins_precession_shuff{ianimal} = prec_bin_shuff; 
                
                if ianimal <= length(learning_points) && ~isnan(learning_points(ianimal))
                    lp = learning_points(ianimal);
                    if lp > 10 && (lp + 9) <= num_trials
                        idx_early = 1:10; idx_pre = (lp - 10) : (lp - 1); idx_post = lp : (lp + 9);
                        get_cols = @(data, cols) data(:, cols);
                        
                        group_results(ipair).trial_corr_early{ianimal} = get_cols(cca_tr, idx_early);
                        group_results(ipair).trial_corr_pre{ianimal}   = get_cols(cca_tr, idx_pre);
                        group_results(ipair).trial_corr_post{ianimal}  = get_cols(cca_tr, idx_post);
                        group_results(ipair).trial_corr_early_shuff{ianimal} = get_cols(cca_tr_shuff, idx_early);
                        group_results(ipair).trial_corr_pre_shuff{ianimal}   = get_cols(cca_tr_shuff, idx_pre);
                        group_results(ipair).trial_corr_post_shuff{ianimal}  = get_cols(cca_tr_shuff, idx_post);
                        
                        group_results(ipair).trial_precession_early{ianimal} = get_cols(prec_tr, idx_early);
                        group_results(ipair).trial_precession_pre{ianimal}   = get_cols(prec_tr, idx_pre);
                        group_results(ipair).trial_precession_post{ianimal}  = get_cols(prec_tr, idx_post);
                        group_results(ipair).trial_precession_early_shuff{ianimal} = get_cols(prec_tr_shuff, idx_early);
                        group_results(ipair).trial_precession_pre_shuff{ianimal}   = get_cols(prec_tr_shuff, idx_pre);
                        group_results(ipair).trial_precession_post_shuff{ianimal}  = get_cols(prec_tr_shuff, idx_post);
                    end
                end
            end 
        catch ME
            fprintf('  Failed animal %d: %s\n', ianimal, ME.message);
        end
    end
    
    % Save .mat for this configuration
    mat_filename = fullfile(save_dir, sprintf('CCA_%s.mat', current_cfg.name));
    save(mat_filename, 'group_results', '-v7.3');
    
    % =====================================================================
    % PLOTTING FOR CURRENT CONFIG
    % =====================================================================
    fprintf('  Generating Plots for %s...\n', current_cfg.name);
    
    % A. Spatial Correlation
    figure('Name', sprintf('Spatial Corr (%s)', current_cfg.name), 'Color', 'w', 'Position', [100 100 1200 800], 'Visible', 'off');
    tiledlayout('flow');
    for ipair = 1:n_pairs
        nexttile;
        [mu, se] = aggregate_cells(group_results(ipair).all_bins_corr, n_bins, 1); 
        [mu_s, ~] = aggregate_cells(group_results(ipair).all_bins_corr_shuff, n_bins, 1);
        if isempty(mu), continue; end
        plot(1:n_bins, mu_s, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
        shadedErrorBar_local(1:n_bins, mu, se, 'b');
        xline((landmarks_cm(:) / track_length_cm) * n_bins); 
        title(group_results(ipair).pair_name); xlabel('Spatial Bin'); ylabel('Correlation (CC1)'); xlim([1 n_bins]);
    end
    save_to_svg(fullfile(save_dir, sprintf('Spatial_Corr_%s', current_cfg.name)));

    % B. Spatial Precession
    figure('Name', sprintf('Spatial Precession (%s)', current_cfg.name), 'Color', 'w', 'Position', [100 100 1200 800], 'Visible', 'off');
    tiledlayout('flow');
    for ipair = 1:n_pairs
        nexttile;
        [mu, se] = aggregate_cells(group_results(ipair).all_bins_precession, n_bins, 1);
        [mu_s, ~] = aggregate_cells(group_results(ipair).all_bins_precession_shuff, n_bins, 1); 
        if isempty(mu), continue; end
        plot(1:n_bins, mu_s, 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'LineStyle', '--'); hold on;
        shadedErrorBar_local(1:n_bins, mu, se, 'r');
        xline((landmarks_cm(:) / track_length_cm) * n_bins); yline(0, '--k');
        title(group_results(ipair).pair_name); xlim([1 n_bins]); ylim([-0.2 0.2]); ylabel('Precession Idx');
    end
    save_to_svg(fullfile(save_dir, sprintf('Spatial_Precession_%s', current_cfg.name)));

    % C. Epoch Correlation
    figure('Name', sprintf('Epoch Corr (%s)', current_cfg.name), 'Color', 'w', 'Visible', 'off');
    tiledlayout('flow');
    for ipair = 1:n_pairs
        e_vals = extract_animal_means(group_results(ipair).trial_corr_early, 1);
        p_vals = extract_animal_means(group_results(ipair).trial_corr_pre, 1);
        x_vals = extract_animal_means(group_results(ipair).trial_corr_post, 1);
        es_vals = extract_animal_means(group_results(ipair).trial_corr_early_shuff, 1);
        ps_vals = extract_animal_means(group_results(ipair).trial_corr_pre_shuff, 1);
        xs_vals = extract_animal_means(group_results(ipair).trial_corr_post_shuff, 1);
        
        if all(isnan(e_vals)), continue; end
        nexttile; hold on;
        means = [mean(e_vals, 'omitnan'), mean(p_vals, 'omitnan'), mean(x_vals, 'omitnan')];
        sems  = [std(e_vals, 'omitnan')/sqrt(sum(~isnan(e_vals))), std(p_vals, 'omitnan')/sqrt(sum(~isnan(p_vals))), std(x_vals, 'omitnan')/sqrt(sum(~isnan(x_vals)))];
        shuffs = [mean(es_vals, 'omitnan'), mean(ps_vals, 'omitnan'), mean(xs_vals, 'omitnan')];
        
        bar(1:3, means, 'FaceColor', [0.8 0.8 0.8]);
        errorbar(1:3, means, sems, 'k.', 'LineWidth', 1.5);
        plot(1:3, shuffs, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', '--');
        xticklabels({'Early', 'Pre', 'Post'}); title(group_results(ipair).pair_name); ylabel('Correlation');
        
        y_stats = [e_vals; p_vals; x_vals];
        g_stats = [ones(length(e_vals),1); 2*ones(length(p_vals),1); 3*ones(length(x_vals),1)];
        valid_mask = ~isnan(y_stats);
        if sum(valid_mask) > 3
            [pval_anova, ~, stats] = anova1(y_stats(valid_mask), g_stats(valid_mask), 'off');
            if pval_anova < 0.05
                c = multcompare(stats, 'Display', 'off'); pairs = {}; pvals_mc = [];
                for k = 1:size(c,1), if c(k,6) < 0.05, pairs{end+1} = [c(k,1), c(k,2)]; pvals_mc(end+1) = c(k,6); end, end
                if ~isempty(pairs), sigstar(pairs, pvals_mc); end
            end
        end
    end
    save_to_svg(fullfile(save_dir, sprintf('Epoch_Corr_%s', current_cfg.name)));

    % D. Epoch Precession
    figure('Name', sprintf('Epoch Precession (%s)', current_cfg.name), 'Color', 'w', 'Visible', 'off');
    tiledlayout('flow');
    for ipair = 1:n_pairs
        e_vals = extract_animal_means(group_results(ipair).trial_precession_early, 1);
        p_vals = extract_animal_means(group_results(ipair).trial_precession_pre, 1);
        x_vals = extract_animal_means(group_results(ipair).trial_precession_post, 1);
        es_vals = extract_animal_means(group_results(ipair).trial_precession_early_shuff, 1);
        ps_vals = extract_animal_means(group_results(ipair).trial_precession_pre_shuff, 1);
        xs_vals = extract_animal_means(group_results(ipair).trial_precession_post_shuff, 1);
        
        if all(isnan(e_vals)), continue; end
        nexttile; hold on;
        means = [mean(e_vals, 'omitnan'), mean(p_vals, 'omitnan'), mean(x_vals, 'omitnan')];
        sems  = [std(e_vals, 'omitnan')/sqrt(sum(~isnan(e_vals))), std(p_vals, 'omitnan')/sqrt(sum(~isnan(p_vals))), std(x_vals, 'omitnan')/sqrt(sum(~isnan(x_vals)))];
        shuffs = [mean(es_vals, 'omitnan'), mean(ps_vals, 'omitnan'), mean(xs_vals, 'omitnan')];
        
        bar(1:3, means, 'FaceColor', [0.8 0.5 0.5]);
        errorbar(1:3, means, sems, 'k.', 'LineWidth', 1.5);
        plot(1:3, shuffs, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', '--'); yline(0, '-k', 'LineWidth', 0.5);
        xticklabels({'Early', 'Pre', 'Post'}); title(group_results(ipair).pair_name); ylabel('Precession Index');
        
        y_stats = [e_vals; p_vals; x_vals];
        g_stats = [ones(length(e_vals),1); 2*ones(length(p_vals),1); 3*ones(length(x_vals),1)];
        valid_mask = ~isnan(y_stats);
        if sum(valid_mask) > 3
            [pval_anova, ~, stats] = anova1(y_stats(valid_mask), g_stats(valid_mask), 'off');
            if pval_anova < 0.05
                c = multcompare(stats, 'Display', 'off'); pairs = {}; pvals_mc = [];
                for k = 1:size(c,1), if c(k,6) < 0.05, pairs{end+1} = [c(k,1), c(k,2)]; pvals_mc(end+1) = c(k,6); end, end
                if ~isempty(pairs), sigstar(pairs, pvals_mc); end
            end
        end
        p_zero = nan(1, 3);
        if sum(~isnan(e_vals)) > 2, [~, p_zero(1)] = ttest(e_vals, 0); end
        if sum(~isnan(p_vals)) > 2, [~, p_zero(2)] = ttest(p_vals, 0); end
        if sum(~isnan(x_vals)) > 2, [~, p_zero(3)] = ttest(x_vals, 0); end
        
        yl = ylim; offset = (yl(2) - yl(1)) * 0.05; min_y_needed = yl(1); 
        for idx = 1:3
            if p_zero(idx) < 0.05
                if means(idx) >= 0, y_pos = means(idx) + sems(idx) + offset; vert_align = 'bottom';
                else, y_pos = means(idx) - sems(idx) - offset; vert_align = 'top'; if y_pos < min_y_needed, min_y_needed = y_pos - offset; end
                end
                text(idx, y_pos, '*', 'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', vert_align);
            end
        end
        if min_y_needed < yl(1), ylim([min_y_needed, yl(2)]); end
    end
    save_to_svg(fullfile(save_dir, sprintf('Epoch_Precession_%s', current_cfg.name)));

    % E. Network Visualization
    layout_def.names = {'CA1', 'V1', 'DG', 'CA3', 'RSC', 'SUB'};
    layout_def.x     = [6.8,  8.5,  3.5,  5.0,  1.5,  3.0];
    layout_def.y     = [6.0,  9.5,  5.5,  2.5,  9.5,  7.5];

    net_data.early_cc = zeros(n_pairs, 1); net_data.early_ifi = zeros(n_pairs, 1);
    net_data.post_cc  = zeros(n_pairs, 1); net_data.post_ifi  = zeros(n_pairs, 1);

    for ipair = 1:n_pairs
        net_data.early_cc(ipair) = mean(extract_animal_means(group_results(ipair).trial_corr_early, 1), 'omitnan');
        net_data.post_cc(ipair) = mean(extract_animal_means(group_results(ipair).trial_corr_post, 1), 'omitnan');
        net_data.early_ifi(ipair) = mean(extract_animal_means(group_results(ipair).trial_precession_early, 1), 'omitnan');
        net_data.post_ifi(ipair) = mean(extract_animal_means(group_results(ipair).trial_precession_post, 1), 'omitnan');
    end

    figure('Name', sprintf('Network (%s)', current_cfg.name), 'Color', 'w', 'Position', [100 100 1600 500], 'Visible', 'off');
    t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'none');
    nexttile(t); plot_directed_network(area_pairs_to_analyze, net_data.early_cc, net_data.early_ifi, 'Early Trials', false, layout_def);
    nexttile(t); plot_directed_network(area_pairs_to_analyze, net_data.post_cc, net_data.post_ifi, 'Expert Trials (Post)', false, layout_def);
    nexttile(t); plot_directed_network(area_pairs_to_analyze, net_data.post_cc - net_data.early_cc, net_data.post_ifi - net_data.early_ifi, 'Difference (Expert - Early)', true, layout_def);
    save_to_svg(fullfile(save_dir, sprintf('Network_%s', current_cfg.name)));

    % Prevent OOM by closing hidden figures generated in this loop
    close all; 
end
fprintf('\n>>> Hyperparameter Sweep Complete. All results saved to: %s <<<\n', save_dir);

%% 4. LOCAL HELPERS

function p_idx = calc_precession(D1, D2, max_shift, num_ccs, nc1, nc2)
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
        else, p_idx(cc) = nan; end
    end
end

function plot_directed_network(pairs, cc_vals, ifi_vals, titleStr, isDiff, layout)
    regions = unique(pairs(:)); sources = []; targets = []; weights = []; edge_colors = [];
    n_edges = size(pairs, 1);
    for i = 1:n_edges
        u = find(strcmp(regions, pairs{i, 1})); v = find(strcmp(regions, pairs{i, 2}));
        cc = cc_vals(i); ifi = ifi_vals(i);
        if isnan(cc) || isnan(ifi), continue; end
        if ifi >= 0, s = u; t = v; else, s = v; t = u; end
        sources = [sources; s]; targets = [targets; t]; weights = [weights; abs(cc)];
        if isDiff
            if cc >= 0, edge_colors = [edge_colors; 0.85, 0.33, 0.10]; 
            else, edge_colors = [edge_colors; 0.00, 0.45, 0.74]; end
        else, edge_colors = [edge_colors; 0.2, 0.2, 0.2]; end
    end
    G = digraph(sources, targets, weights, numel(regions)); p = plot(G);
    x_coords = zeros(1, numel(regions)); y_coords = zeros(1, numel(regions));
    for i = 1:numel(regions)
        idx = find(strcmp(layout.names, regions{i}));
        if ~isempty(idx), x_coords(i) = layout.x(idx); y_coords(i) = layout.y(idx);
        else, x_coords(i) = 0; y_coords(i) = 0; end
    end
    p.XData = x_coords; p.YData = y_coords; axis equal; axis manual; xlim([0 10]); ylim([0 11]); 
    labelnode(p, 1:numel(regions), regions);
    p.NodeColor = [0.9 0.9 0.9]; p.MarkerSize = 18; p.NodeFontSize = 12; p.NodeFontWeight = 'bold';
    max_w = max(G.Edges.Weight);
    if isempty(max_w) || max_w == 0, norm_w = zeros(size(G.Edges.Weight)); else, norm_w = G.Edges.Weight / max_w; end
    p.LineWidth = (norm_w * 7) + 0.5; p.EdgeColor = edge_colors; p.ArrowSize = 12;
    if ~isDiff, p.EdgeAlpha = 0.8; end
    axis off; title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
end

function [mu, se] = aggregate_cells(cell_data, n_bins, cc_idx)
    valid_data = cell_data(~cellfun(@isempty, cell_data)); n_animals = length(valid_data);
    if n_animals == 0, mu = []; se = []; return; end
    stack = nan(n_animals, n_bins);
    for i = 1:n_animals, d = valid_data{i}; if size(d, 2) == n_bins, stack(i, :) = d(cc_idx, :); end; end
    mu = mean(stack, 1, 'omitnan'); se = std(stack, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(stack), 1));
end

function raw_means = extract_animal_means(cell_data, cc_idx)
    n_animals = length(cell_data); raw_means = nan(n_animals, 1);
    for i = 1:n_animals, if ~isempty(cell_data{i}), raw_means(i) = mean(cell_data{i}(cc_idx, :), 2, 'omitnan'); end; end
end

function save_to_svg(fig_name)
    fig = gcf; set(fig, 'Renderer', 'painters');
    try print(fig, '-dsvg', [fig_name '.svg']); catch, print(fig, '-dpng', [fig_name '.png']); end
end

function shadedErrorBar_local(x, y, err, color)
    fill([x, fliplr(x)], [y+err, fliplr(y-err)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on; plot(x, y, 'Color', color, 'LineWidth', 1.5);
end

function ok = check_dims(X, Y, min_samples)
    if isempty(X) || isempty(Y), ok = false; return; end
    n_samples = size(X, 2); if any(isnan(X(:))) || any(isnan(Y(:))), ok = false; return; end
    ok = n_samples > min_samples; 
end