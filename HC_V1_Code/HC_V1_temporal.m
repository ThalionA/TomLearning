%% HC_V1_temporal.m 
% Refactored and Optimized Temporal CCA Analysis
%
% Steps:
% 1. Setup & Parameters
% 2. Pre-processing (Load, Filter FS, Bin, PCA)
% 3. Region Pair Discovery (Hardcoded)
% 4. CCA Analysis (Region-Region w/ Shuffles & Region-Behavior)
% 5. Plotting (Aligned & Consolidated) & Export

%% 1. SETUP & PARAMETERS
clear; clc;

% --- Paths & Files ---
data_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
file_pattern = 'TF*_export.mat';
processed_file = fullfile(data_dir, 'Cohort_Processed.mat');
learning_file = 'animal_behaviour.mat';

% --- Analysis Parameters ---
bin_ms = 50;                  % Time bin size
n_components_reduced = 4;     % PCA components
n_window = -3:3;              % CCA Sliding window (+/- trials)
max_shift_bins = 5;           % Precession/Lag max shift
n_min_animals = 3;            % Min animals to analyze a pair
min_units_per_region = 5;     % Keep only areas with >= 5 units
n_shuffles = 50;              % Number of shuffled controls

% --- Plotting Parameters ---
n_pad = 10;                   % Trials for Early/Pre/Post alignment
target_behavior = 'velocity'; % 'velocity' or 'licks'

% Load Learning Points
lp_path = fullfile(data_dir, learning_file);
if exist(lp_path, 'file')
    dat_lp = load(lp_path);
    learning_points = dat_lp.period_experienced(:, 1);
    animal_idx = dat_lp.animal_id;

    [sorted_animal_idx, sorting_idx] = sort(animal_idx);
    sorted_learning_points = learning_points(sorting_idx);
else
    warning('File "%s" not found. Cannot proceed with alignment plots.', learning_file);
    learning_points = [];
end

% Initialize Master Structure
Cohort = struct('Filename', {}, 'Units', {}, 'Trials', {}, 'Correlations', {}, 'Reduced', {});

%% 2. PRE-PROCESSING (Load, Filter, Bin, PCA)
file_list = dir(fullfile(data_dir, file_pattern));
n_files = length(file_list);

fprintf('>>> S T E P  1: Pre-processing %d files...\n', n_files);

for i_animal = 1:n_files
    filename = file_list(i_animal).name;
    fullpath = fullfile(data_dir, filename);
    fprintf('   Processing %d/%d: %s\n', i_animal, n_files, filename);
    
    try 
        dat = load(fullpath);
        if ~isfield(dat, 'data_behaviour') || ~isfield(dat, 'binned_spikes')
            warning('Skipping %s (Missing vars)', filename);
            continue;
        end
        
        beh = dat.data_behaviour;
        spikes = dat.binned_spikes;
        units = dat.units;
        
        % --- 1. Filter Fast Spiking Units (V1, RSC, CA1, CA3) ---
        if isfield(units, 'idx_fs')
            target_fs_areas = {'V1', 'RSC', 'CA1', 'CA3'};
            fs_vec = logical(units.idx_fs); % Logical vector of FS units
            
            for t_area = 1:length(target_fs_areas)
                area_lbl = target_fs_areas{t_area};
                row_mask = strcmp(units.regions_label, area_lbl);
                
                % Zero out the indices in units.idx for FS units in these areas.
                % This effectively removes them from the region definition 
                % without resizing the global spike matrix.
                if any(row_mask)
                    units.idx(row_mask, fs_vec) = 0; 
                end
            end
        end

        % Store basic info
        Cohort(i_animal).Filename = filename;
        Cohort(i_animal).Units = units;
        
        % --- Define Task Period ---
        task_start = find(~isnan(beh.trial_binned_cued), 1);
        task_end = find(~isnan(beh.trial_binned_cued), 1, 'last');
        
        if isempty(task_start) || isempty(task_end)
            warning('No valid trial data for %s', filename);
            continue;
        end
        
        % Slice Data
        raw_spikes_t = spikes(:, task_start:task_end);
        trial_vec = beh.trial_binned_cued(task_start:task_end);
        pos = beh.pos_binned(task_start:task_end);
        licks = beh.licks_binned(task_start:task_end);
        vel = beh.velocity_binned(task_start:task_end);
        
        [n_units, ~] = size(raw_spikes_t);
        num_trials = trial_vec(end);
        num_areas = length(units.regions_label);
        
        % Pre-allocate Simple Correlations
        Cohort(i_animal).Correlations.Velocity = nan(num_areas, num_trials);
        Cohort(i_animal).Correlations.Licks = nan(num_areas, num_trials);
        
        % Pre-allocate Trials Struct Array
        empty_trial = struct('binned_spikes', [], 'binned_pos', [], ...
                             'binned_licks', [], 'binned_vel', [], ...
                             'Reduced', struct());
        Cohort(i_animal).Trials = repmat(empty_trial, num_trials, 1);
        
        % --- Trial Loop (Binning & Simple Corr) ---
        for itrial = 1:num_trials
            t_idx = trial_vec == itrial;
            if ~any(t_idx), continue; end
            
            % Optimized Binning
            t_spikes = raw_spikes_t(:, t_idx);
            t_pos = pos(t_idx);
            t_licks = licks(t_idx);
            t_vel = vel(t_idx);
            
            n_samps = length(t_pos);
            time_bins = ceil((1:n_samps) / bin_ms)'; 
            n_new = time_bins(end);
            
            % Vectorized Spike Binning
            subs_u = repmat((1:n_units)', 1, n_samps);
            subs_t = repmat(time_bins', n_units, 1);
            
            if isempty(subs_u), continue; end
            
            binned_spikes = accumarray([subs_u(:), subs_t(:)], t_spikes(:), [n_units, n_new]) ./ bin_ms * 1000;
            
            % Behavior Binning
            b_pos = accumarray(time_bins, t_pos, [n_new, 1], @mean);
            b_licks = accumarray(time_bins, t_licks, [n_new, 1], @sum);
            b_vel = accumarray(time_bins, t_vel, [n_new, 1], @mean);
            
            % Store
            Cohort(i_animal).Trials(itrial).binned_spikes = binned_spikes;
            Cohort(i_animal).Trials(itrial).binned_pos = b_pos;
            Cohort(i_animal).Trials(itrial).binned_licks = b_licks;
            Cohort(i_animal).Trials(itrial).binned_vel = b_vel;
            
            % Simple Correlation
            for iarea = 1:num_areas
                a_idx = units.idx(iarea, :);
                % 2. Check for Min Units (for simple corr)
                if sum(a_idx) < min_units_per_region
                    continue; 
                end
                
                act = mean(binned_spikes(a_idx, :), 1)';
                if ~isempty(act) && ~isempty(b_vel) && length(act) == length(b_vel)
                     Cohort(i_animal).Correlations.Velocity(iarea, itrial) = corr(act, b_vel, 'Rows', 'complete');
                     Cohort(i_animal).Correlations.Licks(iarea, itrial) = corr(act, b_licks, 'Rows', 'complete');
                end
            end
        end
        
        % --- PCA Reduction ---
        regions = unique(units.regions_label);
        Cohort(i_animal).Reduced = struct();
        
        for r = 1:length(regions)
            reg = regions{r};
            reg_idx = units.idx(strcmp(units.regions_label, reg), :);
            
            % 2. Keep only areas/animals with at least 5 units
            if sum(reg_idx) < min_units_per_region
                continue;
            end
            
            % Collect data
            data_cells = {Cohort(i_animal).Trials.binned_spikes};
            valid_trials = cellfun(@(x) ~isempty(x), data_cells);
            
            if ~any(valid_trials), continue; end
            
            lens = zeros(1, num_trials);
            temp_data_cell = cell(1, num_trials);
            
            for t = 1:num_trials
                 if valid_trials(t)
                     d = Cohort(i_animal).Trials(t).binned_spikes(reg_idx, :);
                     temp_data_cell{t} = d;
                     lens(t) = size(d, 2);
                 end
            end
            all_data = [temp_data_cell{:}]; 
            
            % PCA Robustness Check
            [n_units_reg, n_time_total] = size(all_data);
            actual_components = min([n_components_reduced, n_units_reg, n_time_total]);
            
            if actual_components > 0
                [~, score, ~] = pca(all_data', 'NumComponents', actual_components);
                score = score'; 
                
                curr = 1;
                for t = 1:num_trials
                    if lens(t) > 0
                        Cohort(i_animal).Trials(t).Reduced.(reg) = score(:, curr:curr+lens(t)-1);
                        curr = curr + lens(t);
                    end
                end
            end
        end
        
    catch ME
        warning('Error processing %s: %s', filename, ME.message);
    end
end
fprintf('Pre-processing complete.\n');

%% 3. REGION & PAIR DISCOVERY
% Hardcoded Pairs
valid_pairs = {
  'CA1', 'V1';
  'CA1', 'DG';
  'CA1', 'CA3';
  'CA1', 'RSC';
  'CA1', 'SUB';
  'V1', 'RSC';
  'RSC', 'SUB';
  'CA3', 'DG'};
n_pairs = size(valid_pairs, 1);
fprintf('Using %d hardcoded region pairs.\n', n_pairs);

%% 4. ANALYSIS: REGION-REGION CCA (With Shuffles)
fprintf('>>> S T E P  3: Running Region-Region CCA...\n');
warning('off', 'stats:canoncorr:NotFullRank'); 

GroupCCA = struct('Name', cell(1, n_pairs), 'AnimalIndex', [], 'Reduced', struct(), 'Full', struct());

for i = 1:n_pairs
    GroupCCA(i).Name = sprintf('%s-%s', valid_pairs{i,1}, valid_pairs{i,2});
    % Data storage
    GroupCCA(i).Reduced.TrialCorr = {}; 
    GroupCCA(i).Reduced.Precession = {};
    GroupCCA(i).Full.TrialCorr = {};    
    GroupCCA(i).Full.Precession = {};
    
    % Shuffle storage (Storing the MEAN of the 50 shuffles per trial)
    GroupCCA(i).Reduced.ShuffleCorr = {}; 
    GroupCCA(i).Reduced.ShufflePrecession = {};
    GroupCCA(i).Full.ShuffleCorr = {};    
    GroupCCA(i).Full.ShufflePrecession = {};
end

% Helper list for mapping
all_regs = {};
for i=1:length(Cohort), if ~isempty(Cohort(i).Units), all_regs=[all_regs, unique(Cohort(i).Units.regions_label)]; end; end
unique_regions = unique(all_regs);
n_regs = length(unique_regions);

for i_animal = 1:length(Cohort)
    animal = Cohort(i_animal);
    if isempty(animal.Units), continue; end
    
    n_trials = length(animal.Trials);
    if n_trials == 0, continue; end
    
    % Cache indices for Full Space
    reg_map = struct();
    for r = 1:n_regs
        reg = unique_regions{r};
        reg_map.(reg) = animal.Units.idx(strcmp(animal.Units.regions_label, reg), :);
    end
    
    for ipair = 1:n_pairs
        a1 = valid_pairs{ipair, 1}; 
        a2 = valid_pairs{ipair, 2};
        
        if ~isfield(animal.Trials(1).Reduced, a1) || ~isfield(animal.Trials(1).Reduced, a2)
            continue;
        end
        
        cca_rho = nan(2, n_trials); 
        cca_shift = nan(2, n_trials, 2*max_shift_bins+1);
        
        % Shuffles: [2 spaces x n_trials x n_shuffles]
        cca_rho_shuff = nan(2, n_trials, n_shuffles);
        cca_shift_shuff = nan(2, n_trials, 2*max_shift_bins+1, n_shuffles);
        
        has_valid_data = false;
        
        for t = 1:n_trials
            win = t + n_window; 
            win = win(win>=1 & win<=n_trials);
            
            X_red_c = cell(1, length(win)); Y_red_c = cell(1, length(win));
            X_full_c = cell(1, length(win)); Y_full_c = cell(1, length(win));
            valid_win_idx = false(1, length(win));
            
            for k = 1:length(win)
                w = win(k);
                if isempty(animal.Trials(w).binned_spikes), continue; end
                
                if isfield(animal.Trials(w).Reduced, a1) && isfield(animal.Trials(w).Reduced, a2)
                    X_red_c{k} = animal.Trials(w).Reduced.(a1);
                    Y_red_c{k} = animal.Trials(w).Reduced.(a2);
                    
                    s = animal.Trials(w).binned_spikes;
                    X_full_c{k} = s(reg_map.(a1), :);
                    Y_full_c{k} = s(reg_map.(a2), :);
                    valid_win_idx(k) = true;
                end
            end
            
            if ~any(valid_win_idx), continue; end
            
            X_red = [X_red_c{:}]; Y_red = [Y_red_c{:}];
            X_full = [X_full_c{:}]; Y_full = [Y_full_c{:}];
            
            % --- REDUCED CCA ---
            if check_dims(X_red, Y_red, n_components_reduced + max_shift_bins + 2)
                % Real
                [~,~,r] = canoncorr(X_red', Y_red');
                cca_rho(1, t) = r(1);
                
                for s_i = 1:(2*max_shift_bins+1)
                    s_val = s_i - (max_shift_bins + 1);
                    [Xc, Yc] = shift_and_cut(X_red, Y_red, s_val);
                    if check_dims(Xc, Yc, n_components_reduced + 2)
                        [~,~,r_l] = canoncorr(Xc', Yc');
                        cca_shift(1, t, s_i) = r_l(1);
                    end
                end
                
                % Shuffles (Circular shift Y to break temporal structure)
                for ish = 1:n_shuffles
                    shift_amount = randi([10, size(Y_red, 2)-10]); % Random circular shift
                    Y_red_s = circshift(Y_red, shift_amount, 2);
                    
                    if check_dims(X_red, Y_red_s, n_components_reduced + max_shift_bins + 2)
                        [~,~,rs] = canoncorr(X_red', Y_red_s');
                        cca_rho_shuff(1, t, ish) = rs(1);
                        
                        % Shuffle Precession
                         for s_i = 1:(2*max_shift_bins+1)
                            s_val = s_i - (max_shift_bins + 1);
                            [Xc, Yc] = shift_and_cut(X_red, Y_red_s, s_val);
                            if check_dims(Xc, Yc, n_components_reduced + 2)
                                [~,~,rls] = canoncorr(Xc', Yc');
                                cca_shift_shuff(1, t, s_i, ish) = rls(1);
                            end
                         end
                    end
                end
                has_valid_data = true;
            end
            
            % --- FULL CCA ---
            n_vars = size(X_full,1) + size(Y_full,1);
            if check_dims(X_full, Y_full, n_vars + 2)
                 % Real
                 [~,~,r] = canoncorr(X_full', Y_full');
                 cca_rho(2, t) = r(1);
                 
                 for s_i = 1:(2*max_shift_bins+1)
                    s_val = s_i - (max_shift_bins + 1);
                    [Xc, Yc] = shift_and_cut(X_full, Y_full, s_val);
                    if check_dims(Xc, Yc, n_vars + 2)
                        [~,~,r_l] = canoncorr(Xc', Yc');
                        cca_shift(2, t, s_i) = r_l(1);
                    end
                 end
                 
                 % Shuffles
                 for ish = 1:n_shuffles
                     shift_amount = randi([10, size(Y_full, 2)-10]);
                     Y_full_s = circshift(Y_full, shift_amount, 2);
                     
                     if check_dims(X_full, Y_full_s, n_vars + 2)
                         [~,~,rs] = canoncorr(X_full', Y_full_s');
                         cca_rho_shuff(2, t, ish) = rs(1);
                         
                         % Shuffle Precession
                         for s_i = 1:(2*max_shift_bins+1)
                            s_val = s_i - (max_shift_bins + 1);
                            [Xc, Yc] = shift_and_cut(X_full, Y_full_s, s_val);
                            if check_dims(Xc, Yc, n_vars + 2)
                                [~,~,rls] = canoncorr(Xc', Yc');
                                cca_shift_shuff(2, t, s_i, ish) = rls(1);
                            end
                         end
                     end
                 end
                 has_valid_data = true;
            end
        end
        
        if has_valid_data
            % 1. Process Real Data
            a2_lead_r = squeeze(mean(cca_shift(1, :, 1:max_shift_bins), 3, 'omitnan'));
            a1_lead_r = squeeze(mean(cca_shift(1, :, max_shift_bins+2:end), 3, 'omitnan'));
            prec_str_r = (a1_lead_r - a2_lead_r) ./ (a1_lead_r + a2_lead_r);
            
            a2_lead_f = squeeze(mean(cca_shift(2, :, 1:max_shift_bins), 3, 'omitnan'));
            a1_lead_f = squeeze(mean(cca_shift(2, :, max_shift_bins+2:end), 3, 'omitnan'));
            prec_str_f = (a1_lead_f - a2_lead_f) ./ (a1_lead_f + a2_lead_f);
            
            % 2. Process Shuffle Data
            % Calculate mean correlation across shuffles
            shuff_mean_rho_r = mean(squeeze(cca_rho_shuff(1,:,:)), 2, 'omitnan')';
            shuff_mean_rho_f = mean(squeeze(cca_rho_shuff(2,:,:)), 2, 'omitnan')';
            
            % Process Shuffle Precession
            % cca_shift_shuff is [2 x Trials x ShiftBins x Shuffles]
            a2_lead_rs = squeeze(mean(cca_shift_shuff(1, :, 1:max_shift_bins, :), 3, 'omitnan'));
            a1_lead_rs = squeeze(mean(cca_shift_shuff(1, :, max_shift_bins+2:end, :), 3, 'omitnan'));
            % Calculate precession per shuffle, then average
            prec_str_rs_all = (a1_lead_rs - a2_lead_rs) ./ (a1_lead_rs + a2_lead_rs);
            shuff_mean_prec_r = mean(prec_str_rs_all, 2, 'omitnan')';
            
            a2_lead_fs = squeeze(mean(cca_shift_shuff(2, :, 1:max_shift_bins, :), 3, 'omitnan'));
            a1_lead_fs = squeeze(mean(cca_shift_shuff(2, :, max_shift_bins+2:end, :), 3, 'omitnan'));
            prec_str_fs_all = (a1_lead_fs - a2_lead_fs) ./ (a1_lead_fs + a2_lead_fs);
            shuff_mean_prec_f = mean(prec_str_fs_all, 2, 'omitnan')';

            % Store
            GroupCCA(ipair).AnimalIndex(end+1) = i_animal;
            
            GroupCCA(ipair).Reduced.TrialCorr{end+1} = cca_rho(1, :);
            GroupCCA(ipair).Reduced.Precession{end+1} = prec_str_r;
            GroupCCA(ipair).Reduced.ShuffleCorr{end+1} = shuff_mean_rho_r;
            GroupCCA(ipair).Reduced.ShufflePrecession{end+1} = shuff_mean_prec_r;
            
            GroupCCA(ipair).Full.TrialCorr{end+1} = cca_rho(2, :);
            GroupCCA(ipair).Full.Precession{end+1} = prec_str_f;
            GroupCCA(ipair).Full.ShuffleCorr{end+1} = shuff_mean_rho_f;
            GroupCCA(ipair).Full.ShufflePrecession{end+1} = shuff_mean_prec_f;
        end
        fprintf('pair %d/%d done \n', ipair, n_pairs)
    end
    fprintf('animal %d/%d done \n', i_animal, length(Cohort))
end
fprintf('Region-Region CCA Done.\n');

%% 5. ANALYSIS: REGION-BEHAVIOR CCA
fprintf('>>> S T E P  3b: Running Region-Behavior CCA...\n');
BehavCCA = struct('Name', {}, 'AnimalIndex', {}, 'TrialCorr', {});
for i=1:n_regs 
    BehavCCA(i).Name = unique_regions{i}; 
    BehavCCA(i).AnimalIndex=[]; 
    BehavCCA(i).TrialCorr={}; 
end

for i_animal = 1:length(Cohort)
    animal = Cohort(i_animal);
    if isempty(animal.Units), continue; end
    n_trials = length(animal.Trials);
    
    reg_map = struct();
    for r=1:n_regs, reg = unique_regions{r}; reg_map.(reg) = animal.Units.idx(strcmp(animal.Units.regions_label, reg), :); end
    
    for i_reg = 1:n_regs
        reg = unique_regions{i_reg};
        if n_trials == 0 || ~isfield(animal.Trials(1).Reduced, reg), continue; end
        
        idx_units = reg_map.(reg);
        n_u = sum(idx_units);
        if n_u < min_units_per_region, continue; end % Apply limit here too
        
        cca_rho = nan(1, n_trials);
        has_data = false;
        
        for t = 1:n_trials
             win = t + n_window; win = win(win>=1 & win<=n_trials);
             
             X_c = cell(1, length(win)); Y_c = cell(1, length(win));
             valid_win = false(1, length(win));
             for k = 1:length(win)
                 w = win(k);
                 s = animal.Trials(w).binned_spikes;
                 if isempty(s), continue; end
                 
                 if strcmp(target_behavior, 'velocity')
                     b = animal.Trials(w).binned_vel;
                 else
                     b = animal.Trials(w).binned_licks;
                 end
                 
                 x = s(idx_units, :);
                 min_l = min(size(x, 2), length(b));
                 if min_l > 5
                     X_c{k} = x(:, 1:min_l);
                     Y_c{k} = b(1:min_l)';
                     valid_win(k) = true;
                 end
             end
             
             if ~any(valid_win), continue; end
             X_pool = [X_c{:}]; Y_pool = [Y_c{:}];
             
             if check_dims(X_pool, Y_pool, n_u + 5)
                 [~,~,r] = canoncorr(X_pool', Y_pool');
                 cca_rho(t) = r(1);
                 has_data = true;
             end
        end
        
        if has_data
            BehavCCA(i_reg).AnimalIndex(end+1) = i_animal;
            BehavCCA(i_reg).TrialCorr{end+1} = cca_rho;
        end
    end
end
fprintf('Behavior CCA Done.\n');

%% 6. PLOTTING
fprintf('>>> S T E P  4: Plotting Results...\n');
if isempty(learning_points)
    warning('No learning points loaded. Skipping plots.');
    return;
end

% --- Plot 1: Simple Velocity Correlation (Aligned) ---
f1 = figure('Name', 'Simple Corr Aligned', 'Color', 'w');
tiledlayout('flow', 'TileSpacing', 'compact');
for r = 1:n_regs
    reg = unique_regions{r};
    data_mat = nan(length(Cohort), 100); 
    valid_id = [];
    
    for i=1:length(Cohort)
        if ~isempty(Cohort(i).Correlations) && ~isempty(Cohort(i).Filename)
            idx = find(strcmp(Cohort(i).Units.regions_label, reg));
            if ~isempty(idx)
                tr = mean(Cohort(i).Correlations.Velocity(idx, :), 1, 'omitnan');
                if length(tr) > size(data_mat, 2), data_mat(i, length(tr)) = 0; end
                data_mat(i, 1:length(tr)) = tr;
                valid_id(end+1) = i; 
            end
        end
    end
    
    if ~isempty(valid_id)
        [d_early, d_pre, d_post, n_valid] = align_data(data_mat(valid_id, :), learning_points(valid_id), n_pad);
        if n_valid > 0
            nexttile; hold on;
            plot_aligned_set(d_early, d_pre, d_post, [0.2 0.2 0.2], n_pad);
            title(sprintf('%s (n=%d)', reg, n_valid));
            ylim([-0.2 0.6]); yline(0, '--k');
        end
    end
end
sgtitle('Simple Activity-Velocity Correlation');
print(f1, 'Simple_Corr_Aligned.eps', '-depsc', '-painters'); % 3. Export EPS

% --- Plot 2: Region-Region CCA (Aligned with Shuffles) ---
f2 = figure('Name', 'All Paired CCA Aligned', 'Color', 'w');
t = tiledlayout('flow', 'TileSpacing', 'compact');
h_leg = []; 
for i = 1:n_pairs
    if isempty(GroupCCA(i).AnimalIndex), continue; end
    
    idxs = GroupCCA(i).AnimalIndex;
    if max(idxs) <= length(learning_points)
        % Align Real Data
        [r_early, r_pre, r_post, ~] = align_data_cells(GroupCCA(i).Reduced.TrialCorr, learning_points(idxs), n_pad);
        [f_early, f_pre, f_post, n_v] = align_data_cells(GroupCCA(i).Full.TrialCorr, learning_points(idxs), n_pad);
        
        % Align Shuffle Data
        [rs_early, rs_pre, rs_post, ~] = align_data_cells(GroupCCA(i).Reduced.ShuffleCorr, learning_points(idxs), n_pad);
        
        if n_v > 0
            nexttile; hold on;
            
            % Plot Shuffles (Gray)
            plot_aligned_set(rs_early, rs_pre, rs_post, [0.7 0.7 0.7], n_pad);
            
            % Plot Real
            plot_aligned_set(r_early, r_pre, r_post, 'b', n_pad); 
            plot_aligned_set(f_early, f_pre, f_post, 'r', n_pad); 
            
            title(sprintf('%s (n=%d)', GroupCCA(i).Name, n_v));
            ylabel('Correlation');
            
            if isempty(h_leg)
                h_leg(1) = plot(nan, nan, 'Color', 'b', 'LineWidth', 2);
                h_leg(2) = plot(nan, nan, 'Color', 'r', 'LineWidth', 2);
                h_leg(3) = fill(nan, nan, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            end
        end
    end
end
linkaxes(findall(gcf,'type','axes'), 'xy'); 
if ~isempty(h_leg)
    legend(h_leg, 'Reduced', 'Full', 'Shuffle'); 
end
print(f2, 'Paired_CCA_Aligned.eps', '-depsc', '-painters'); % 3. Export EPS

% --- Plot 3: Region-Region Precession (Aligned with Shuffles) ---
f3 = figure('Name', 'All Paired Precession Aligned', 'Color', 'w');
t = tiledlayout('flow', 'TileSpacing', 'compact');
h_leg = []; 
for i = 1:n_pairs
    if isempty(GroupCCA(i).AnimalIndex), continue; end
    
    idxs = GroupCCA(i).AnimalIndex;
    if max(idxs) <= length(learning_points)
        % Align Real
        [r_early, r_pre, r_post, ~] = align_data_cells(GroupCCA(i).Reduced.Precession, learning_points(idxs), n_pad);
        [f_early, f_pre, f_post, n_v] = align_data_cells(GroupCCA(i).Full.Precession, learning_points(idxs), n_pad);
        
        % Align Shuffle
        [rs_early, rs_pre, rs_post, ~] = align_data_cells(GroupCCA(i).Reduced.ShufflePrecession, learning_points(idxs), n_pad);
        
        if n_v > 0
            nexttile; hold on;
            % Plot Shuffle
            plot_aligned_set(rs_early, rs_pre, rs_post, [0.7 0.7 0.7], n_pad);
            
            % Plot Real
            plot_aligned_set(r_early, r_pre, r_post, 'b', n_pad); 
            plot_aligned_set(f_early, f_pre, f_post, 'r', n_pad); 
            
            title(sprintf('%s (n=%d)', GroupCCA(i).Name, n_v));
            ylabel('Precession Idx');
            yline(0, '--k');
            ylim([-0.2 0.2]); 
            
            if isempty(h_leg)
                h_leg(1) = plot(nan, nan, 'Color', 'b', 'LineWidth', 2);
                h_leg(2) = plot(nan, nan, 'Color', 'r', 'LineWidth', 2);
                h_leg(3) = fill(nan, nan, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            end
        end
    end
end
linkaxes(findall(gcf,'type','axes'), 'xy'); 
if ~isempty(h_leg)
    legend(h_leg, 'Reduced', 'Full', 'Shuffle'); 
end
print(f3, 'Paired_Precession_Aligned.eps', '-depsc', '-painters'); % 3. Export EPS

% --- Plot 4: Region-Behavior CCA (Aligned) ---
f4 = figure('Name', 'Behavior CCA Aligned', 'Color', 'w');
tiledlayout('flow', 'TileSpacing', 'compact');
for i = 1:n_regs
    idxs = BehavCCA(i).AnimalIndex;
    if isempty(idxs), continue; end
    
    if max(idxs) <= length(learning_points)
        [d_early, d_pre, d_post, n_v] = align_data_cells(BehavCCA(i).TrialCorr, learning_points(idxs), n_pad);
        
        if n_v > 0
            nexttile; hold on;
            plot_aligned_set(d_early, d_pre, d_post, [0.85 0.325 0.098], n_pad);
            title(sprintf('%s (n=%d)', BehavCCA(i).Name, n_v));
            ylim([0 0.8]);
        end
    end
end
sgtitle(['Full Space CCA vs ' target_behavior]);
print(f4, 'Behavior_CCA_Aligned.eps', '-depsc', '-painters'); % 3. Export EPS

fprintf('>>> Analysis Complete.\n');

%% 7. HELPER FUNCTIONS
function ok = check_dims(X, Y, min_samples)
    if isempty(X) || isempty(Y), ok = false; return; end
    if size(X, 2) ~= size(Y, 2), ok = false; return; end
    n_samples = size(X, 2);
    n_feats = size(X, 1) + size(Y, 1);
    if any(isnan(X(:))) || any(isnan(Y(:))), ok = false; return; end
    ok = n_samples > n_feats + 2; 
    if nargin >= 3, ok = ok && (n_samples >= min_samples); end
end
function [Xc, Yc] = shift_and_cut(X, Y, s)
    if s < 0, Xc = X(:, 1-s:end); Yc = Y(:, 1:end+s);
    elseif s > 0, Xc = X(:, 1:end-s); Yc = Y(:, 1+s:end);
    else, Xc = X; Yc = Y; end
end
function [mat, max_len] = pad_trials(cell_arr)
    if isempty(cell_arr), mat=[]; max_len=0; return; end
    max_len = max(cellfun(@length, cell_arr));
    mat = nan(length(cell_arr), max_len);
    for k=1:length(cell_arr), d=cell_arr{k}; mat(k, 1:length(d))=d; end
end
function plot_shaded(x, data, col)
    if isempty(data), return; end
    mu = mean(data,1,'omitnan');
    se = std(data,0,1,'omitnan')./sqrt(sum(~isnan(data),1));
    valid = ~isnan(mu);
    if ~any(valid), return; end
    x_poly = [x(valid) fliplr(x(valid))];
    y_poly = [mu(valid)+se(valid) fliplr(mu(valid)-se(valid))];
    fill(x_poly, y_poly, col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(x(valid), mu(valid), 'Color', col, 'LineWidth', 2);
end
function [early, pre, post, count] = align_data(mat, lps, npad)
    [n_anim, n_t] = size(mat);
    early = mat(:, 1:npad);
    pre = nan(n_anim, npad); post = nan(n_anim, npad);
    count = 0;
    for k=1:n_anim
        lp = lps(k);
        if ~isnan(lp) && lp>npad && (lp+npad-1)<=n_t
            pre(k,:) = mat(k, lp-npad:lp-1);
            post(k,:) = mat(k, lp:lp+npad-1);
            count = count+1;
        end
    end
end
function [early, pre, post, count] = align_data_cells(cells, lps, npad)
    n_anim = length(cells);
    early = nan(n_anim, npad); pre = nan(n_anim, npad); post = nan(n_anim, npad);
    count = 0;
    for k=1:n_anim
        trace = cells{k};
        if length(trace)>=npad, early(k,:)=trace(1:npad); end
        lp = lps(k);
        if ~isnan(lp) && lp>npad && (lp+npad-1)<=length(trace)
            pre(k,:) = trace(lp-npad:lp-1);
            post(k,:) = trace(lp:lp+npad-1);
            count = count+1;
        end
    end
end
function plot_aligned_set(e, pr, po, col, npad)
    gap=2;
    xe = 1:npad; 
    xp = xe(end)+gap : xe(end)+gap+npad-1;
    xpo = xp(end)+gap : xp(end)+gap+npad-1;
    plot_shaded(xe, e, col);
    plot_shaded(xp, pr, col);
    plot_shaded(xpo, po, col);
    xticks([mean(xe) mean(xp) mean(xpo)]);
    xticklabels({'Early','Pre','Post'});
end