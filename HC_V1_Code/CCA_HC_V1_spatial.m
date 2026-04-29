%% CCA_HC_V1_Spatial_Final_Optimized.m
% Final Spatial CCA Analysis with Shuffled Precession Controls
%
% DESCRIPTION:
% This script performs Canonical Correlation Analysis (CCA) to quantify shared 
% variance and lead/lag relationships (Precession).
%
% 1. TRIAL-WISE ANALYSIS (Spatial Map Similarity):
%    - Control: Shuffle Spatial Bins (destroys spatial topology).
%
% 2. BIN-WISE ANALYSIS (Co-fluctuation):
%    - Control: Shuffle Trial Identities (destroys co-fluctuation).
%
% 3. PRECESSION ANALYSIS (Lead/Lag):
%    - Logic: Shift Area A relative to B.
%    - Control: Same as above (Shuffle Bins for Trial-wise, Trials for Bin-wise).
%
% PARAMETERS:
% - Spatial Resolution: 200 Bins (2.5 cm/bin)
% - Shuffles: 50 (Parallelized)
%
% OUTPUTS:
% - Incremental saving of .mat files
% - SVG Figures

%% 1. SETUP & PARAMETERS
clear; clc;

% --- Paths ---
base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir = fullfile(base_dir, data_subfolder);
file_pattern = 'TF*_export.mat';
learning_file = 'animal_behaviour.mat';

% --- Save Path ---
current_date = datestr(now, 'yyyy_mm_dd');
save_path = fullfile(data_dir, sprintf('Spatial_CCA_Results_%s.mat', current_date));

% --- Parameters ---
n_components_reduced = 4;    % PCA components to keep
num_ccs_analyze = 2;         % Number of CCA components to store/plot
n_trials_window = -3:3;      % Sliding window for Trial-wise analysis
n_bins_window = -3:3;        % Sliding window for Bin-wise analysis
n_shuffles = 50;             % Number of shuffles
max_shift_bins = 3;          % Max shift for precession
min_units_per_region = 5;

% Spatial Parameters
n_bins = 200;                % Fixed bin count
bin_size_cm = 2.5;           % cm per bin
track_length_cm = n_bins * bin_size_cm; 
landmarks_cm = [0, 50; 100, 125; 175, 200; 250, 275; 325, 350; 400, 425]; 

% --- Pairs to Analyze ---
area_pairs_to_analyze = {'CA1', 'V1'; 'CA1', 'DG'; 'CA1', 'CA3'; 'CA1', 'RSC'; ...
                         'CA1', 'SUB'; 'V1', 'RSC'; 'RSC', 'SUB'; 'CA3', 'DG'};
n_pairs = size(area_pairs_to_analyze, 1);

% --- Parallel Pool ---
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

%% 3. INITIALIZE RESULTS STRUCTURE
% Structure: group_results(pair_idx).field{animal_idx} = [Matrix]
group_results = struct('pair_name', cell(n_pairs, 1), ...
                       'all_bins_corr', cell(n_pairs, 1), ...        % [nCCs x nBins]
                       'all_bins_corr_shuff', cell(n_pairs, 1), ...  % [nCCs x nBins]
                       'all_bins_precession', cell(n_pairs, 1), ...  % [nCCs x nBins]
                       'all_bins_precession_shuff', cell(n_pairs, 1), ... % <--- ADDED (Bin Shuffle)
                       ...
                       'trial_corr_early', cell(n_pairs, 1), ...     % [nCCs x 10]
                       'trial_corr_pre', cell(n_pairs, 1), ...
                       'trial_corr_post', cell(n_pairs, 1), ...
                       ...
                       'trial_corr_early_shuff', cell(n_pairs, 1), ...
                       'trial_corr_pre_shuff', cell(n_pairs, 1), ...
                       'trial_corr_post_shuff', cell(n_pairs, 1), ...
                       ...
                       'trial_precession_early', cell(n_pairs, 1), ... 
                       'trial_precession_pre', cell(n_pairs, 1), ...   
                       'trial_precession_post', cell(n_pairs, 1), ...
                       ...
                       'trial_precession_early_shuff', cell(n_pairs, 1), ... % <--- ADDED (Epoch Shuffle)
                       'trial_precession_pre_shuff', cell(n_pairs, 1), ...   % <--- ADDED
                       'trial_precession_post_shuff', cell(n_pairs, 1));     % <--- ADDED

for ipair = 1:n_pairs
    group_results(ipair).pair_name = sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});
end
all_lick_errors = {};

%% 4. MAIN ANALYSIS LOOP
file_list = dir(fullfile(data_dir, file_pattern));
n_animals = length(file_list);
fprintf('Found %d animals. Starting analysis...\n', n_animals);

for ianimal = 1:n_animals
    filename = file_list(ianimal).name;
    fullpath = fullfile(data_dir, filename);
    fprintf('\nProcessing Animal %d/%d: %s\n', ianimal, n_animals, filename);
    
    try
        D = load(fullpath);
        
        % --- A. Load & Filter Units ---
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
        
        % --- B. Get Spatial Data ---
        if isfield(D.analysis_spatial, 'firing') && isfield(D.analysis_spatial.firing, 'cued')
             raw_spatial = D.analysis_spatial.firing.cued.freq;
             % Permute to [Units x Bins x Trials]
             animal_data = permute(raw_spatial, [1, 3, 2]);
             [~, curr_n_bins, num_trials] = size(animal_data);
             
             if curr_n_bins ~= n_bins
                 warning('Bin mismatch for %s. Found %d, expected %d. Skipping.', filename, curr_n_bins, n_bins);
                 continue;
             end
        else
            continue;
        end
        
        % --- C. Licking Data ---
        if isfield(D.analysis_behaviour, 'spatial_licking_gf')
            licks = D.analysis_behaviour.spatial_licking_gf; 
            bin_indices = repmat(1:n_bins, num_trials, 1);
            target_bin = n_bins; 
            current_lick_error = mean(licks .* (bin_indices - target_bin).^2, 2);
            z_lick_error = (current_lick_error - mean(current_lick_error)) ./ std(current_lick_error);
            all_lick_errors{ianimal} = z_lick_error;
        else
            all_lick_errors{ianimal} = nan(num_trials, 1);
        end
        
        % --- D. PCA Reduction ---
        animal_areas = unique(units.regions_label);
        AreaActivity = struct();
        
        for ia = 1:length(animal_areas)
            area_name = animal_areas{ia};
            u_idx = units.idx(strcmp(units.regions_label, area_name), :);
            if sum(u_idx) < min_units_per_region, continue; end
            
            area_dat = animal_data(u_idx, :, :);
            reshaped_dat = reshape(area_dat, sum(u_idx), [])'; 
            
            if size(reshaped_dat, 2) >= n_components_reduced
                [~, scores, ~] = pca(reshaped_dat, 'NumComponents', n_components_reduced);
                reduced = reshape(scores', n_components_reduced, n_bins, num_trials);
                AreaActivity.(area_name) = reduced;
            end
        end
        
        % --- E. CCA Analysis Loop ---
        for ipair = 1:n_pairs
            a1 = area_pairs_to_analyze{ipair, 1};
            a2 = area_pairs_to_analyze{ipair, 2};
            
            if ~isfield(AreaActivity, a1) || ~isfield(AreaActivity, a2)
                continue; 
            end
            
            d1 = AreaActivity.(a1);
            d2 = AreaActivity.(a2);
            
            % Preallocate Local Results
            cca_tr = nan(num_ccs_analyze, num_trials);
            cca_tr_shuff = nan(num_ccs_analyze, num_trials); 
            
            prec_tr = nan(num_ccs_analyze, num_trials);
            prec_tr_shuff = nan(num_ccs_analyze, num_trials); % <--- Added
            
            cca_bin = nan(num_ccs_analyze, n_bins);
            cca_bin_shuff = nan(num_ccs_analyze, n_bins); 
            
            prec_bin = nan(num_ccs_analyze, n_bins);
            prec_bin_shuff = nan(num_ccs_analyze, n_bins); % <--- Added
            
            % =========================================================
            % ANALYSIS 1: TRIAL-WISE (Are maps similar?)
            % Control: Shuffle BINS (Space)
            % =========================================================
            for t = 1:num_trials
                win = t + n_trials_window;
                win = win(win>=1 & win<=num_trials);
                if isempty(win), continue; end
                
                % Collapse time window: [Comps x Bins x nWin]
                D1_local = d1(:, :, win); 
                D2_local = d2(:, :, win);
                
                x = reshape(D1_local, n_components_reduced, []);
                y = reshape(D2_local, n_components_reduced, []);
                
                if check_dims(x, y, n_components_reduced + 5)
                    % 1. Real CCA
                    [~,~,r] = canoncorr(x', y');
                    cca_tr(1:min(length(r), num_ccs_analyze), t) = r(1:min(length(r), num_ccs_analyze));
                    
                    % 3. Real Precession
                    prec_tr(:, t) = calc_precession(D1_local, D2_local, max_shift_bins, num_ccs_analyze, n_components_reduced);
                    
                    % --- SHUFFLE BLOCK ---
                    r_sh_iter = nan(num_ccs_analyze, n_shuffles);
                    p_sh_iter = nan(num_ccs_analyze, n_shuffles);
                    
                    n_bins_local = size(D1_local, 2);
                    
                    parfor ish = 1:n_shuffles
                         % Shuffle SPATIAL BINS (Dimension 2)
                         perm_idx = randperm(n_bins_local);
                         
                         % Use permuted data for Precession (keep structure)
                         D1_s = D1_local(:, perm_idx, :); 
                         
                         % Use permuted data for CCA (flattened)
                         x_s = reshape(D1_s, n_components_reduced, []);
                         
                         % A. Shuffled CCA
                         [~,~,rs] = canoncorr(x_s', y');
                         r_sh_iter(:, ish) = rs(1:min(length(rs), num_ccs_analyze));
                         
                         % B. Shuffled Precession
                         p_sh_iter(:, ish) = calc_precession(D1_s, D2_local, max_shift_bins, num_ccs_analyze, n_components_reduced);
                    end
                    cca_tr_shuff(:, t) = mean(r_sh_iter, 2, 'omitnan');
                    prec_tr_shuff(:, t) = mean(p_sh_iter, 2, 'omitnan');
                end
            end
            
            % =========================================================
            % ANALYSIS 2: BIN-WISE (Do they co-fluctuate?)
            % Control: Shuffle TRIALS
            % =========================================================
            for b = 1:n_bins
                win = b + n_bins_window;
                win = win(win>=1 & win<=n_bins);
                if isempty(win), continue; end
                
                D1_local = d1(:, win, :);
                D2_local = d2(:, win, :);
                
                x = reshape(D1_local, n_components_reduced, []);
                y = reshape(D2_local, n_components_reduced, []);
                
                if check_dims(x, y, n_components_reduced + 5)
                    % 1. Real CCA
                    [~,~,r] = canoncorr(x', y');
                    cca_bin(1:min(length(r), num_ccs_analyze), b) = r(1:min(length(r), num_ccs_analyze));
                    
                    % 3. Real Precession
                    prec_bin(:, b) = calc_precession(D1_local, D2_local, max_shift_bins, num_ccs_analyze, n_components_reduced);
                    
                    % --- SHUFFLE BLOCK ---
                    r_sh_iter = nan(num_ccs_analyze, n_shuffles);
                    p_sh_iter = nan(num_ccs_analyze, n_shuffles);
                    
                    n_tr_local = size(D1_local, 3);
                    
                    parfor ish = 1:n_shuffles
                         % Shuffle TRIALS (Dimension 3)
                         perm_idx = randperm(n_tr_local);
                         D1_s = D1_local(:, :, perm_idx);
                         
                         x_s = reshape(D1_s, n_components_reduced, []);
                         
                         % A. Shuffled CCA
                         [~,~,rs] = canoncorr(x_s', y');
                         r_sh_iter(:, ish) = rs(1:min(length(rs), num_ccs_analyze));
                         
                         % B. Shuffled Precession
                         p_sh_iter(:, ish) = calc_precession(D1_s, D2_local, max_shift_bins, num_ccs_analyze, n_components_reduced);
                    end
                    cca_bin_shuff(:, b) = mean(r_sh_iter, 2, 'omitnan');
                    prec_bin_shuff(:, b) = mean(p_sh_iter, 2, 'omitnan');
                end
            end
            
            % --- F. Store Results ---
            
            % 1. Spatial Bin Results
            group_results(ipair).all_bins_corr{ianimal} = cca_bin;
            group_results(ipair).all_bins_corr_shuff{ianimal} = cca_bin_shuff;
            group_results(ipair).all_bins_precession{ianimal} = prec_bin;
            group_results(ipair).all_bins_precession_shuff{ianimal} = prec_bin_shuff; % <--- Store
            
            % 2. Trial Results (Align by Learning Point)
            if ianimal <= length(learning_points) && ~isnan(learning_points(ianimal))
                lp = learning_points(ianimal);
                if lp > 10 && (lp + 9) <= num_trials
                    idx_early = 1:10;
                    idx_pre   = (lp - 10) : (lp - 1);
                    idx_post  = lp : (lp + 9);
                    get_cols = @(data, cols) data(:, cols);
                    
                    % Correlation
                    group_results(ipair).trial_corr_early{ianimal} = get_cols(cca_tr, idx_early);
                    group_results(ipair).trial_corr_pre{ianimal}   = get_cols(cca_tr, idx_pre);
                    group_results(ipair).trial_corr_post{ianimal}  = get_cols(cca_tr, idx_post);
                    
                    group_results(ipair).trial_corr_early_shuff{ianimal} = get_cols(cca_tr_shuff, idx_early);
                    group_results(ipair).trial_corr_pre_shuff{ianimal}   = get_cols(cca_tr_shuff, idx_pre);
                    group_results(ipair).trial_corr_post_shuff{ianimal}  = get_cols(cca_tr_shuff, idx_post);
                    
                    % Precession (Real)
                    group_results(ipair).trial_precession_early{ianimal} = get_cols(prec_tr, idx_early);
                    group_results(ipair).trial_precession_pre{ianimal}   = get_cols(prec_tr, idx_pre);
                    group_results(ipair).trial_precession_post{ianimal}  = get_cols(prec_tr, idx_post);
                    
                    % Precession (Shuffle) <--- STORED HERE
                    group_results(ipair).trial_precession_early_shuff{ianimal} = get_cols(prec_tr_shuff, idx_early);
                    group_results(ipair).trial_precession_pre_shuff{ianimal}   = get_cols(prec_tr_shuff, idx_pre);
                    group_results(ipair).trial_precession_post_shuff{ianimal}  = get_cols(prec_tr_shuff, idx_post);
                end
            end
        end % End Pairs
        
        save(save_path, 'group_results', 'all_lick_errors', '-v7.3');
    catch ME
        fprintf('Error processing %s: %s\n', filename, ME.message);
    end
end
fprintf('\nAnalysis Complete. Final Save...\n');
save(save_path, 'group_results', 'all_lick_errors', '-v7.3');

%% 5. PLOTTING (Summary Figures)
fprintf('Generating Plots...\n');

% --- A. Spatial Correlation (Bin-wise) ---
figure('Name', 'Group Spatial Correlation', 'Color', 'w', 'Position', [100 100 1200 800]);
tiledlayout('flow');
for ipair = 1:n_pairs
    nexttile
    [mu, se] = aggregate_cells(group_results(ipair).all_bins_corr, n_bins, 1); 
    [mu_s, ~] = aggregate_cells(group_results(ipair).all_bins_corr_shuff, n_bins, 1);
    
    if isempty(mu), continue; end
    
    x = 1:n_bins;
    plot(x, mu_s, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
    shadedErrorBar_local(x, mu, se, 'b');
    
    lm_bins = (landmarks_cm(:) / track_length_cm) * n_bins;
    xline(lm_bins); 
    title(group_results(ipair).pair_name);
    xlabel('Spatial Bin'); ylabel('Correlation (CC1)');
    xlim([1 n_bins]);
end
linkaxes
save_to_svg(fullfile(data_dir, 'Spatial_Corr_Group'));

% --- B. Spatial Precession ---
figure('Name', 'Group Spatial Precession', 'Color', 'w', 'Position', [100 100 1200 800]);
tiledlayout('flow');
for ipair = 1:n_pairs
    nexttile
    [mu, se] = aggregate_cells(group_results(ipair).all_bins_precession, n_bins, 1);
    [mu_s, ~] = aggregate_cells(group_results(ipair).all_bins_precession_shuff, n_bins, 1); % Shuffle
    
    if isempty(mu), continue; end
    
    x = 1:n_bins;
    % Plot Shuffle as grey line
    plot(x, mu_s, 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'LineStyle', '--'); hold on;
    shadedErrorBar_local(x, mu, se, 'r');
    
    lm_bins = (landmarks_cm(:) / track_length_cm) * n_bins;
    xline(lm_bins); 
    yline(0, '--k');
    title(group_results(ipair).pair_name);
    xlim([1 n_bins]); ylim([-0.2 0.2]);
    ylabel('Precession Idx (CC1)');
end
linkaxes
save_to_svg(fullfile(data_dir, 'Spatial_Precession_Group'));

% --- C. Trial Correlation (Epochs) ---
figure('Name', 'Group Epoch Correlation', 'Color', 'w');
tiledlayout('flow');
for ipair = 1:n_pairs
    [m_early, s_early] = aggregate_epochs(group_results(ipair).trial_corr_early);
    [m_pre, s_pre]     = aggregate_epochs(group_results(ipair).trial_corr_pre);
    [m_post, s_post]   = aggregate_epochs(group_results(ipair).trial_corr_post);
    
    [ms_early, ~] = aggregate_epochs(group_results(ipair).trial_corr_early_shuff);
    [ms_pre, ~]   = aggregate_epochs(group_results(ipair).trial_corr_pre_shuff);
    [ms_post, ~]  = aggregate_epochs(group_results(ipair).trial_corr_post_shuff);
    
    if isnan(m_early), continue; end
    
    nexttile
    means = [m_early, m_pre, m_post];
    sems  = [s_early, s_pre, s_post];
    shuffs = [ms_early, ms_pre, ms_post];
    
    bar(1:3, means, 'FaceColor', [0.8 0.8 0.8]); hold on;
    errorbar(1:3, means, sems, 'k.', 'LineWidth', 1.5);
    plot(1:3, shuffs, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', '--');
    
    xticklabels({'Early', 'Pre', 'Post'});
    title(group_results(ipair).pair_name);
    ylabel('Correlation (CC1)');
end
linkaxes
save_to_svg(fullfile(data_dir, 'Trial_Corr_Epochs'));

% --- D. Trial Precession (Epochs) ---
figure('Name', 'Group Epoch Precession', 'Color', 'w');
tiledlayout('flow');
for ipair = 1:n_pairs
    % Real Data
    [m_early, s_early] = aggregate_epochs(group_results(ipair).trial_precession_early);
    [m_pre, s_pre]     = aggregate_epochs(group_results(ipair).trial_precession_pre);
    [m_post, s_post]   = aggregate_epochs(group_results(ipair).trial_precession_post);
    
    % Shuffle Data
    [ms_early, ~] = aggregate_epochs(group_results(ipair).trial_precession_early_shuff);
    [ms_pre, ~]   = aggregate_epochs(group_results(ipair).trial_precession_pre_shuff);
    [ms_post, ~]  = aggregate_epochs(group_results(ipair).trial_precession_post_shuff);
    
    if isnan(m_early), continue; end
    
    nexttile
    means = [m_early, m_pre, m_post];
    sems  = [s_early, s_pre, s_post];
    shuffs = [ms_early, ms_pre, ms_post];
    
    bar(1:3, means, 'FaceColor', [0.8 0.5 0.5]); hold on;
    errorbar(1:3, means, sems, 'k.', 'LineWidth', 1.5);
    
    % Add Shuffle Level (Expectation ~0)
    plot(1:3, shuffs, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', '--');
    yline(0, '-k', 'LineWidth', 0.5);
    
    xticklabels({'Early', 'Pre', 'Post'});
    title(group_results(ipair).pair_name);
    ylabel('Precession Index');
end
linkaxes
save_to_svg(fullfile(data_dir, 'Trial_Precession_Epochs'));
fprintf('Done.\n');

%% E. NETWORK VISUALIZATION (Directed Graph)
fprintf('Generating Network Graphs...\n');

% 1. Setup Layout (Anatomically Calibrated "C-Shape")
% We use a 0-10 scale to allow spacing from the edges.
% Order: {'CA1', 'V1', 'DG', 'CA3', 'RSC', 'SUB'}
layout_def.names = {'CA1', 'V1', 'DG', 'CA3', 'RSC', 'SUB'};
layout_def.x     = [6.8,  8.5,  3.5,  5.0,  1.5,  3.0];
layout_def.y     = [6.0,  9.5,  5.5,  2.5,  9.5,  7.5];

% 2. Extract Data using aggregate_epochs
net_data.early_cc  = zeros(n_pairs, 1);
net_data.early_ifi = zeros(n_pairs, 1);
net_data.post_cc   = zeros(n_pairs, 1);
net_data.post_ifi  = zeros(n_pairs, 1);

for ipair = 1:n_pairs
    % -- Early --
    % We use [mu, ~] to discard the SEM, as we only need the mean for edge width
    [net_data.early_cc(ipair), ~]  = aggregate_epochs(group_results(ipair).trial_corr_early);
    [net_data.early_ifi(ipair), ~] = aggregate_epochs(group_results(ipair).trial_precession_early);
    
    % -- Post (Expert) --
    [net_data.post_cc(ipair), ~]   = aggregate_epochs(group_results(ipair).trial_corr_post);
    [net_data.post_ifi(ipair), ~]  = aggregate_epochs(group_results(ipair).trial_precession_post);
end

% 3. Generate Plot
figure('Name', 'Network Connectivity', 'Color', 'w', 'Position', [100 100 1600 500]);
t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'none');

% --- Subplot 1: Early ---
nexttile(t);
plot_directed_network(area_pairs_to_analyze, ...
    net_data.early_cc, net_data.early_ifi, ...
    'Early Trials', false, layout_def);

% --- Subplot 2: Expert (Post) ---
nexttile(t);
plot_directed_network(area_pairs_to_analyze, ...
    net_data.post_cc, net_data.post_ifi, ...
    'Expert Trials (Post)', false, layout_def);

% --- Subplot 3: Difference (Post - Early) ---
diff_cc  = net_data.post_cc - net_data.early_cc;
diff_ifi = net_data.post_ifi - net_data.early_ifi;

nexttile(t);
plot_directed_network(area_pairs_to_analyze, ...
    diff_cc, diff_ifi, ...
    'Difference (Expert - Early)', true, layout_def);

% Save
save_to_svg(fullfile(data_dir, 'Network_Connectivity_Group'));
fprintf('Done.\n');


% --- LOCAL PLOTTING FUNCTION ---
function plot_directed_network(pairs, cc_vals, ifi_vals, titleStr, isDiff, layout)
    % Get unique regions from the pairs list
    regions = unique(pairs(:));
    
    sources = [];
    targets = [];
    weights = [];
    edge_colors = [];
    
    n_edges = size(pairs, 1);
    
    for i = 1:n_edges
        u_name = pairs{i, 1};
        v_name = pairs{i, 2};
        
        % Map string names to numeric indices 1..N
        u = find(strcmp(regions, u_name));
        v = find(strcmp(regions, v_name));
        
        cc = cc_vals(i);
        ifi = ifi_vals(i);
        
        % Check for NaNs (in case aggregate_epochs returned NaN for missing data)
        if isnan(cc) || isnan(ifi)
            continue; 
        end

        % Direction Logic: Positive IFI = u->v, Negative IFI = v->u
        if ifi >= 0
            s = u; t = v;
        else
            s = v; t = u;
        end
        
        sources = [sources; s];
        targets = [targets; t];
        weights = [weights; abs(cc)];
        
        % Color Logic
        if isDiff
            if cc >= 0
                edge_colors = [edge_colors; 0.85, 0.33, 0.10]; % Orange (Increase)
            else
                edge_colors = [edge_colors; 0.00, 0.45, 0.74]; % Blue (Decrease)
            end
        else
            edge_colors = [edge_colors; 0.2, 0.2, 0.2]; % Dark Grey
        end
    end

    % Create Graph
    G = digraph(sources, targets, weights, numel(regions));
    p = plot(G);
    
    % --- Apply Anatomical Layout ---
    x_coords = zeros(1, numel(regions));
    y_coords = zeros(1, numel(regions));
    
    for i = 1:numel(regions)
        reg = regions{i};
        idx = find(strcmp(layout.names, reg));
        if ~isempty(idx)
            x_coords(i) = layout.x(idx);
            y_coords(i) = layout.y(idx);
        else
            x_coords(i) = 0; y_coords(i) = 0; % Fallback
        end
    end
    p.XData = x_coords;
    p.YData = y_coords;
    
    % --- CRITICAL FIX: Lock Aspect Ratio & Set Bounds ---
    axis equal;    % Ensures physical distances are consistent (no stretching)
    axis manual;   % Prevents auto-scaling to data limits
    xlim([0 10]);  % Set wider frame than data range (1.5 to 8.5)
    ylim([0 11]);  % Set taller frame
    
    % --- Styling ---
    labelnode(p, 1:numel(regions), regions);
    p.NodeColor = [0.9 0.9 0.9];
    p.MarkerSize = 18;
    p.NodeFontSize = 12;
    p.NodeFontWeight = 'bold';
    
    % Scale Line Width (Prevent 0-width error)
    max_w = max(G.Edges.Weight);
    if isempty(max_w) || max_w == 0
        norm_w = zeros(size(G.Edges.Weight));
    else
        norm_w = G.Edges.Weight / max_w;
    end
    p.LineWidth = (norm_w * 7) + 0.5; 
    
    p.EdgeColor = edge_colors;
    p.ArrowSize = 12;
    
    if ~isDiff
        p.EdgeAlpha = 0.8;
    end
    
    axis off;
    title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
end

%% 6. LOCAL HELPERS
function p_idx = calc_precession(D1, D2, max_shift, num_ccs, n_comps)
    shifts = -max_shift : max_shift;
    r_shifts = zeros(num_ccs, length(shifts));
    
    for is = 1:length(shifts)
        shift = shifts(is);
        if shift < 0      % A Left -> A(1:end-k) vs B(1+k:end)
             s1 = D1(:, 1:end+shift, :); 
             s2 = D2(:, 1-shift:end, :);
        elseif shift > 0  % A Right -> A(1+k:end) vs B(1:end-k)
             s1 = D1(:, 1+shift:end, :);
             s2 = D2(:, 1:end-shift, :);
        else
             s1 = D1; s2 = D2;
        end
        
        sx = reshape(s1, n_comps, []);
        sy = reshape(s2, n_comps, []);
        
        if size(sx,2) > n_comps + 2
             [~,~,rs] = canoncorr(sx', sy');
             r_shifts(1:min(length(rs), num_ccs), is) = rs(1:min(length(rs), num_ccs));
        end
    end
    
    idx_neg = 1:max_shift; 
    idx_pos = max_shift+2 : length(shifts);
    p_idx = nan(num_ccs, 1);
    
    for cc = 1:num_ccs
        neg_val = mean(r_shifts(cc, idx_neg), 'omitnan');
        pos_val = mean(r_shifts(cc, idx_pos), 'omitnan');
        % Simple difference for index
        if (neg_val + pos_val) > 0.001 
            p_idx(cc) = (neg_val - pos_val) / (neg_val + pos_val);
        else
            p_idx(cc) = nan;
        end
    end
end

function [mu, se] = aggregate_cells(cell_data, n_bins, cc_idx)
    valid_data = cell_data(~cellfun(@isempty, cell_data));
    n_animals = length(valid_data);
    if n_animals == 0
        mu = []; se = []; return;
    end
    stack = nan(n_animals, n_bins);
    for i = 1:n_animals
        d = valid_data{i};
        if size(d, 2) == n_bins
            stack(i, :) = d(cc_idx, :);
        end
    end
    mu = mean(stack, 1, 'omitnan');
    se = std(stack, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(stack), 1));
end

function [mu, se] = aggregate_epochs(cell_data)
    valid_data = cell_data(~cellfun(@isempty, cell_data));
    if isempty(valid_data)
        mu = nan; se = nan; return;
    end
    animal_means = [];
    for i = 1:length(valid_data)
        d = valid_data{i}; 
        animal_means(end+1) = mean(d(1, :), 2, 'omitnan'); % CC1 only
    end
    mu = mean(animal_means, 'omitnan');
    se = std(animal_means, 'omitnan') / sqrt(length(animal_means));
end

function save_to_svg(fig_name)
    fig = gcf;
    set(fig, 'Renderer', 'painters'); 
    fprintf('Saving %s.svg...\n', fig_name);
    try
        print(fig, '-dsvg', [fig_name '.svg']);
    catch
        print(fig, '-dpng', [fig_name '.png']); 
    end
end

function shadedErrorBar_local(x, y, err, color)
    fill([x, fliplr(x)], [y+err, fliplr(y-err)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on;
    plot(x, y, 'Color', color, 'LineWidth', 1.5);
end

function ok = check_dims(X, Y, min_samples)
    if isempty(X) || isempty(Y), ok = false; return; end
    if size(X, 2) ~= size(Y, 2), ok = false; return; end
    n_samples = size(X, 2);
    if any(isnan(X(:))) || any(isnan(Y(:))), ok = false; return; end
    ok = n_samples > min_samples; 
end