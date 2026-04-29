%% HC_V1_temporal_v6.m
% Final Optimized Temporal CCA Analysis (Dual Mode)
%
% Modes:
% 1. NO FILTER: All trials in cued tunnel.
% 2. TRIAL FILTER: Only trials where animal is active (>2cm/s) for >50% of bins.
%
% Features:
% - Detailed Progress Tracking (File, Animal, Region Pair)
% - Parallelized Shuffles
% - Caching of results
% - SVG Export

%% 1. SETUP & PARAMETERS
clear; clc;

% --- Paths ---
base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir = fullfile(base_dir, data_subfolder);
file_pattern = 'TF*_export.mat';
learning_file = 'animal_behaviour.mat';

% --- Results Caching ---
current_date = datetime('today', 'Format', 'uuuu-MM-dd');
results_filename = sprintf('Analysis_Results_Dual_%s.mat', current_date);
results_path = fullfile(data_dir, results_filename);

% --- Analysis Parameters ---
bin_ms = 50;                  % Time bin size
n_components_reduced = 4;     % PCA components
n_window = -3:3;              % CCA Sliding window (+/- trials)
max_shift_bins = 5;           % Precession/Lag max shift
n_min_animals = 3;            % Min animals to analyze a pair
min_units_per_region = 5;     % Keep only areas with >= 5 units
n_shuffles = 100;             % Number of shuffles

% --- Velocity Parameters ---
min_speed_threshold = 2;      % cm/s definition of "active"
trial_min_active_frac = 0.5;  % Trial must be active for >50% bins (For Filtered Mode)

% --- Parallelization Setup ---
if isempty(gcp('nocreate'))
    parpool; 
end

% --- Plotting Parameters ---
n_pad = 10;                   % Trials for Naive/Pre/Expert alignment

% Load Learning Points
lp_path = fullfile(data_dir, learning_file);
if exist(lp_path, 'file')
    dat_lp = load(lp_path);
    learning_points = dat_lp.period_experienced(:, 1);
    
    % Sort learning points to match file list (alphabetical TF number)
    if isfield(dat_lp, 'animal_id')
        [~, sorting_idx] = sort(dat_lp.animal_id);
        learning_points = learning_points(sorting_idx);
    end
else
    warning('Learning file "%s" not found. Plots will not be aligned.', lp_path);
    learning_points = [];
end

%% 2. CHECK FOR CACHED RESULTS
if exist(results_path, 'file')
    fprintf('>>> Found cached results: %s\n', results_filename);
    fprintf('>>> Loading...\n');
    load(results_path);
    fprintf('>>> Data loaded. Skipping analysis steps.\n');
else
    fprintf('>>> No cached results found. Starting fresh analysis.\n');
    
    %% 3. PRE-PROCESSING (Load, Bin, PCA, Calculate Activity)
    file_list = dir(fullfile(data_dir, file_pattern));
    n_files = length(file_list);
    
    Cohort = struct('Filename', {}, 'Units', {}, 'Trials', {}, 'Reduced', {});
    
    fprintf('\n=== STEP 1: PRE-PROCESSING ===\n');
    
    for i_animal = 1:n_files
        filename = file_list(i_animal).name;
        fullpath = fullfile(data_dir, filename);
        
        fprintf('[Pre-Proc] File %d/%d: %-20s ... ', i_animal, n_files, filename);
        
        try 
            dat = load(fullpath);
            if ~isfield(dat, 'data_behaviour') || ~isfield(dat, 'binned_spikes')
                fprintf('Skipped (Missing Data)\n');
                continue;
            end
            
            beh = dat.data_behaviour;
            spikes = dat.binned_spikes;
            units = dat.units;
            
            % --- A. Filter Fast Spiking Units ---
            if isfield(units, 'idx_fs')
                target_fs_areas = {'V1', 'RSC', 'CA1', 'CA3'};
                fs_vec = logical(units.idx_fs);
                for t_area = 1:length(target_fs_areas)
                    row_mask = strcmp(units.regions_label, target_fs_areas{t_area});
                    if any(row_mask)
                        units.idx(row_mask, fs_vec) = 0; 
                    end
                end
            end
            
            Cohort(i_animal).Filename = filename;
            Cohort(i_animal).Units = units;
            
            % --- B. Context Filtering Only ---
            valid_mask = ~isnan(beh.trial_binned_cued);
            
            % Slice Data to Task Bounds
            task_start = find(valid_mask, 1, 'first');
            task_end = find(valid_mask, 1, 'last');
            
            if isempty(task_start) || isempty(task_end)
                fprintf('Skipped (No Valid Task Data)\n');
                continue;
            end
            
            raw_spikes_t = spikes(:, task_start:task_end);
            trial_vec = beh.trial_binned_cued(task_start:task_end);
            
            % Get Velocity
            n_total_bins = size(spikes, 2);
            if isfield(beh, 'velocity_binned_gf')
                full_vel = beh.velocity_binned_gf;
            else
                full_vel = beh.velocity_gf;
            end
            v_gf = full_vel(1:min(length(full_vel), n_total_bins));
            v_gf_sliced = v_gf(task_start:task_end);
            
            [n_units, ~] = size(raw_spikes_t);
            num_trials = max(trial_vec);
            
            % Pre-allocate Trials
            Cohort(i_animal).Trials = repmat(struct('binned_spikes', [], ...
                'frac_active', 0, 'Reduced', struct()), num_trials, 1);
            
            % --- C. Binning Loop ---
            for itrial = 1:num_trials
                t_idx = (trial_vec == itrial);
                if ~any(t_idx), continue; end
                
                % Calculate Activity Metric
                trial_vels = v_gf_sliced(t_idx);
                frac_active = mean(trial_vels >= min_speed_threshold);
                
                t_spikes = raw_spikes_t(:, t_idx);
                n_samps = size(t_spikes, 2);
                time_bins = ceil((1:n_samps) / bin_ms)'; 
                n_new = time_bins(end);
                
                % Vectorized Binning
                subs_u = repmat((1:n_units)', 1, n_samps);
                subs_t = repmat(time_bins', n_units, 1);
                
                binned_spikes = accumarray([subs_u(:), subs_t(:)], t_spikes(:), [n_units, n_new]) ./ bin_ms * 1000;
                
                Cohort(i_animal).Trials(itrial).binned_spikes = binned_spikes;
                Cohort(i_animal).Trials(itrial).frac_active = frac_active;
            end
            
            % --- D. PCA Reduction ---
            regions = unique(units.regions_label);
            for r = 1:length(regions)
                reg = regions{r};
                reg_idx = units.idx(strcmp(units.regions_label, reg), :);
                if sum(reg_idx) < min_units_per_region, continue; end
                
                data_cells = {Cohort(i_animal).Trials.binned_spikes};
                valid_trials = find(~cellfun(@isempty, data_cells));
                if isempty(valid_trials), continue; end
                
                temp_data = cell(1, length(valid_trials));
                lens = zeros(1, num_trials);
                
                for k = 1:length(valid_trials)
                    t = valid_trials(k);
                    d = Cohort(i_animal).Trials(t).binned_spikes(reg_idx, :);
                    temp_data{k} = d;
                    lens(t) = size(d, 2);
                end
                all_data = [temp_data{:}];
                
                [n_units_reg, n_time_total] = size(all_data);
                actual_components = min([n_components_reduced, n_units_reg, n_time_total]);
                
                if actual_components > 0 && n_time_total > actual_components + 5
                    [~, score, ~] = pca(all_data', 'NumComponents', actual_components);
                    score = score'; 
                    curr = 1;
                    for k = 1:length(valid_trials)
                        t = valid_trials(k);
                        Cohort(i_animal).Trials(t).Reduced.(reg) = score(:, curr:curr+lens(t)-1);
                        curr = curr + lens(t);
                    end
                end
            end
            
            fprintf('Done.\n');
            
        catch ME
            fprintf('Error: %s\n', ME.message);
        end
    end
    fprintf('Pre-processing complete.\n');

    %% 4. RUN ANALYSES (Dual Mode)
    modes = [false, true];
    mode_names = {'All_Trials', 'Trial_Filtered'};
    
    GroupCCA_All = [];
    GroupCCA_Filtered = [];
    
    for m = 1:2
        use_filter = modes(m);
        curr_name = mode_names{m};
        
        fprintf('\n==============================================\n');
        fprintf('>>> STARTING ANALYSIS MODE: %s \n', curr_name);
        fprintf('==============================================\n');
        
        % Run the CCA Engine
        ResultStruct = run_cca_engine(Cohort, n_components_reduced, max_shift_bins, ...
                                      n_window, n_shuffles, use_filter, trial_min_active_frac, curr_name);
                                  
        if use_filter
            GroupCCA_Filtered = ResultStruct;
        else
            GroupCCA_All = ResultStruct;
        end
    end

    %% 5. SAVE RESULTS
    fprintf('\n>>> Saving results to %s...\n', results_path);
    save(results_path, 'GroupCCA_All', 'GroupCCA_Filtered', 'Cohort', 'learning_points', '-v7.3');
end

%% 7. PLOTTING (Dual Mode)
fprintf('\n=== STEP 4: PLOTTING RESULTS ===\n');
if isempty(learning_points)
    warning('No learning points loaded. Skipping plots.');
    return;
end

plot_cca_results(GroupCCA_All, learning_points, n_pad, 'All_Trials', data_dir);
plot_cca_results(GroupCCA_Filtered, learning_points, n_pad, 'Trial_Filtered', data_dir);

fprintf('>>> Analysis & Plotting Complete.\n');


%% 8. LOCAL HELPER FUNCTIONS

% --- MAIN CCA ENGINE ---
function GroupCCA = run_cca_engine(Cohort, n_comps, max_shift, n_win, n_shuff, use_filter, min_act, mode_name)
    valid_pairs = {'CA1', 'V1'; 'CA1', 'DG'; 'CA1', 'CA3'; 'CA1', 'RSC'; 'CA1', 'SUB'; 'V1', 'RSC'; 'RSC', 'SUB'; 'CA3', 'DG'};
    n_pairs = size(valid_pairs, 1);
    
    % Initialize Struct
    GroupCCA = struct('Name', cell(1, n_pairs), 'AnimalIndex', [], 'Reduced', struct(), 'Full', struct());
    for ipair = 1:n_pairs
        GroupCCA(ipair).Name = sprintf('%s-%s', valid_pairs{ipair,1}, valid_pairs{ipair,2});
        GroupCCA(ipair).Reduced.TrialCorr = {}; GroupCCA(ipair).Reduced.Precession = {};
        GroupCCA(ipair).Full.TrialCorr = {};    GroupCCA(ipair).Full.Precession = {};
        GroupCCA(ipair).Reduced.ShuffleCorr = {}; GroupCCA(ipair).Reduced.ShufflePrecession = {};
        GroupCCA(ipair).Full.ShuffleCorr = {};    GroupCCA(ipair).Full.ShufflePrecession = {};
    end

    all_regs = {};
    for i=1:length(Cohort), if ~isempty(Cohort(i).Units), all_regs=[all_regs, unique(Cohort(i).Units.regions_label)]; end; end
    unique_regions = unique(all_regs);
    n_regs = length(unique_regions);
    n_animals = length(Cohort);

    for i_animal = 1:n_animals
        animal = Cohort(i_animal);
        if isempty(animal.Trials), continue; end
        
        fprintf('\n[%s] Processing Animal %d/%d (%s):\n', mode_name, i_animal, n_animals, animal.Filename);
        
        reg_map = struct();
        for r = 1:n_regs
            reg = unique_regions{r};
            reg_map.(reg) = animal.Units.idx(strcmp(animal.Units.regions_label, reg), :);
        end
        
        for ipair = 1:n_pairs
            a1 = valid_pairs{ipair, 1}; 
            a2 = valid_pairs{ipair, 2};
            
            fprintf('   > Pair %d/%d (%s-%s): ', ipair, n_pairs, a1, a2);
            
            if ~isfield(animal.Trials(1).Reduced, a1) || ~isfield(animal.Trials(1).Reduced, a2)
                fprintf('Skipped (Missing Region)\n');
                continue;
            end
            
            n_trials = length(animal.Trials);
            rho_red = nan(1, n_trials); rho_full = nan(1, n_trials);
            prec_red = nan(1, n_trials); prec_full = nan(1, n_trials);
            shuff_rho_red = nan(1, n_trials); shuff_rho_full = nan(1, n_trials);
            shuff_prec_red = nan(1, n_trials); shuff_prec_full = nan(1, n_trials);
            
            has_valid_data = false;
            
            % Process Trials
            for t = 1:n_trials
                % Filter check
                if use_filter && (animal.Trials(t).frac_active < min_act)
                    continue; 
                end
                
                win = t + n_win; 
                win = win(win>=1 & win<=n_trials);
                
                X_red = []; Y_red = []; X_full = []; Y_full = [];
                
                for w = win
                   if isempty(animal.Trials(w).binned_spikes), continue; end
                   if use_filter && (animal.Trials(w).frac_active < min_act), continue; end
                   
                   if isfield(animal.Trials(w).Reduced, a1) && isfield(animal.Trials(w).Reduced, a2)
                       X_red = [X_red, animal.Trials(w).Reduced.(a1)];
                       Y_red = [Y_red, animal.Trials(w).Reduced.(a2)];
                       s = animal.Trials(w).binned_spikes;
                       X_full = [X_full, s(reg_map.(a1), :)];
                       Y_full = [Y_full, s(reg_map.(a2), :)];
                   end
                end
                
                if isempty(X_red) || isempty(Y_red), continue; end
                
                % Process Red/Full (One dot per trial)
                ok_r = check_dims(X_red, Y_red, n_comps + max_shift + 2);
                ok_f = check_dims(X_full, Y_full, size(X_full,1)+size(Y_full,1) + 2);
                
                if ok_r
                    [r_real, p_real] = calc_cca_precession(X_red, Y_red, n_comps, max_shift);
                    rho_red(t) = r_real; prec_red(t) = p_real;
                    [shuff_rho_red(t), shuff_prec_red(t)] = run_shuffles_par(X_red, Y_red, n_comps, max_shift, n_shuff);
                end
                
                if ok_f
                    n_vars = size(X_full,1) + size(Y_full,1);
                    [r_real, p_real] = calc_cca_precession(X_full, Y_full, n_vars, max_shift);
                    rho_full(t) = r_real; prec_full(t) = p_real;
                    [shuff_rho_full(t), shuff_prec_full(t)] = run_shuffles_par(X_full, Y_full, n_vars, max_shift, n_shuff);
                end
                
                if ok_r || ok_f, has_valid_data = true; fprintf('.'); end
            end
            
            if has_valid_data
                GroupCCA(ipair).AnimalIndex(end+1) = i_animal;
                GroupCCA(ipair).Reduced.TrialCorr{end+1} = rho_red;
                GroupCCA(ipair).Reduced.Precession{end+1} = prec_red;
                GroupCCA(ipair).Reduced.ShuffleCorr{end+1} = shuff_rho_red;
                GroupCCA(ipair).Reduced.ShufflePrecession{end+1} = shuff_prec_red;
                GroupCCA(ipair).Full.TrialCorr{end+1} = rho_full;
                GroupCCA(ipair).Full.Precession{end+1} = prec_full;
                GroupCCA(ipair).Full.ShuffleCorr{end+1} = shuff_rho_full;
                GroupCCA(ipair).Full.ShufflePrecession{end+1} = shuff_prec_full;
                fprintf(' Done.\n');
            else
                fprintf(' No valid trials found.\n');
            end
        end
    end
end

% --- Parallel Shuffle Wrapper ---
function [m_r, m_p] = run_shuffles_par(X, Y, n_comps, max_shift, n_shuff)
    r_shuffs = nan(n_shuff, 1);
    p_shuffs = nan(n_shuff, 1);
    
    parfor ish = 1:n_shuff
        shift_amount = randi([10, size(Y, 2)-10]);
        Y_shuff = circshift(Y, shift_amount, 2);
        [r_s, p_s] = calc_cca_precession_local(X, Y_shuff, n_comps, max_shift);
        r_shuffs(ish) = r_s;
        p_shuffs(ish) = p_s;
    end
    m_r = mean(r_shuffs, 'omitnan');
    m_p = mean(p_shuffs, 'omitnan');
end

% --- PLOTTING FUNCTION ---
function plot_cca_results(G, lps, npad, tag, save_dir)
    % CCA Plot
    f = figure('Name', ['CCA Aligned - ' tag], 'Color', 'w');
    tiledlayout('flow', 'TileSpacing', 'compact');
    for i = 1:length(G)
        if isempty(G(i).AnimalIndex), continue; end
        idxs = G(i).AnimalIndex;
        if max(idxs) <= length(lps)
            [r_n, r_p, r_e] = align_data_cells(G(i).Reduced.TrialCorr, lps(idxs), npad);
            [f_n, f_p, f_e] = align_data_cells(G(i).Full.TrialCorr, lps(idxs), npad);
            [rs_n, rs_p, rs_e] = align_data_cells(G(i).Reduced.ShuffleCorr, lps(idxs), npad);
            
            n_v = sum(~isnan(mean(r_n, 2)));
            if n_v > 0
                nexttile; hold on;
                plot_shaded_error(rs_n, rs_p, rs_e, [0.7 0.7 0.7], npad);
                plot_shaded_error(r_n, r_p, r_e, 'b', npad); 
                plot_shaded_error(f_n, f_p, f_e, 'r', npad); 
                title(sprintf('%s (n=%d)', G(i).Name, n_v));
                ylabel('Correlation'); ylim([0 0.8]);
            end
        end
    end
    save_to_svg(fullfile(save_dir, ['CCA_Aligned_' tag]));

    % Precession Plot
    f2 = figure('Name', ['Precession Aligned - ' tag], 'Color', 'w');
    tiledlayout('flow', 'TileSpacing', 'compact');
    for i = 1:length(G)
        if isempty(G(i).AnimalIndex), continue; end
        idxs = G(i).AnimalIndex;
        if max(idxs) <= length(lps)
            [r_n, r_p, r_e] = align_data_cells(G(i).Reduced.Precession, lps(idxs), npad);
            [rs_n, rs_p, rs_e] = align_data_cells(G(i).Reduced.ShufflePrecession, lps(idxs), npad);
            
            n_v = sum(~isnan(mean(r_n, 2)));
            if n_v > 0
                nexttile; hold on;
                plot_shaded_error(rs_n, rs_p, rs_e, [0.7 0.7 0.7], npad);
                plot_shaded_error(r_n, r_p, r_e, 'b', npad); 
                title(sprintf('%s (n=%d)', G(i).Name, n_v));
                ylabel('Precession Idx'); yline(0, '--k'); ylim([-0.2 0.2]); 
            end
        end
    end
    save_to_svg(fullfile(save_dir, ['Precession_Aligned_' tag]));
end

% --- Internal shadedErrorBar ---
function plot_shaded_error(na, pr, ex, col, npad)
    gap = 2;
    x_na = 1:npad; 
    x_pr = x_na(end)+gap : x_na(end)+gap+npad-1;
    x_ex = x_pr(end)+gap : x_pr(end)+gap+npad-1;
    
    do_shaded_plot(x_na, na, col);
    do_shaded_plot(x_pr, pr, col);
    do_shaded_plot(x_ex, ex, col);
    
    xticks([mean(x_na) mean(x_pr) mean(x_ex)]);
    xticklabels({'Naive','Pre-Expert','Expert'});
end

function do_shaded_plot(x, data, col)
    if any(~isnan(mean(data, 1)))
        y_mean = mean(data, 1, 'omitnan');
        y_std = std(data, 0, 1, 'omitnan');
        n_obs = sum(~isnan(data), 1);
        y_sem = y_std ./ sqrt(n_obs);
        valid = ~isnan(y_mean);
        if any(valid)
            upper = y_mean + y_sem;
            lower = y_mean - y_sem;
            fill([x(valid) fliplr(x(valid))], [upper(valid) fliplr(lower(valid))], ...
                 col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            plot(x(valid), y_mean(valid), 'Color', col, 'LineWidth', 2);
        end
    end
end

function save_to_svg(fig_name)
    fig = gcf;
    set(fig, 'Renderer', 'painters'); 
    [fpath, fname, ~] = fileparts(fig_name);
    full_path = fullfile(fpath, [fname '.svg']);
    fprintf('Saving figure to %s...\n', full_path);
    print(fig, '-dsvg', full_path);
end

function ok = check_dims(X, Y, min_samples)
    if isempty(X) || isempty(Y), ok = false; return; end
    if size(X, 2) ~= size(Y, 2), ok = false; return; end
    n_samples = size(X, 2);
    if any(isnan(X(:))) || any(isnan(Y(:))), ok = false; return; end
    ok = n_samples > min_samples; 
end

% Function copy for internal usage (avoids scoping issues in parfor)
function [rho, prec_str] = calc_cca_precession_local(X, Y, n_comps, max_shift)
    [~,~,r] = canoncorr(X', Y');
    rho = r(1);
    r_shifts = nan(1, 2*max_shift+1);
    for ishift = -max_shift:max_shift
        if ishift < 0
            Xc = X(:, 1-ishift:end); Yc = Y(:, 1:end+ishift);
        elseif ishift > 0
            Xc = X(:, 1:end-ishift); Yc = Y(:, 1+ishift:end);
        else
            Xc = X; Yc = Y;
        end
        if size(Xc, 2) > n_comps + 2
             [~,~,rs] = canoncorr(Xc', Yc');
             r_shifts(ishift + max_shift + 1) = rs(1);
        end
    end
    y_leads = mean(r_shifts(1:max_shift), 'omitnan');
    x_leads = mean(r_shifts(max_shift+2:end), 'omitnan');
    prec_str = (x_leads - y_leads) ./ (x_leads + y_leads);
end

% Duplicated for main scope access if needed
function [rho, prec_str] = calc_cca_precession(X, Y, n_comps, max_shift)
    [rho, prec_str] = calc_cca_precession_local(X, Y, n_comps, max_shift);
end

function [naive, pre, expert, count] = align_data_cells(cells, lps, npad)
    n_anim = length(cells);
    naive = nan(n_anim, npad); 
    pre = nan(n_anim, npad); 
    expert = nan(n_anim, npad);
    count = 0;
    for k=1:n_anim
        trace = cells{k};
        if isempty(trace), continue; end
        if length(trace) >= npad
            naive(k, :) = trace(1:npad); 
        end
        lp = lps(k);
        if ~isnan(lp) && lp > npad && (lp + npad - 1) <= length(trace)
            pre(k, :) = trace(lp-npad : lp-1);
            expert(k, :) = trace(lp : lp+npad-1);
            count = count + 1;
        end
    end
end