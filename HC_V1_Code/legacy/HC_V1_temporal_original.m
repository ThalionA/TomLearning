%% 1. SETUP & PARAMETERS
clear; clc;

% --- Paths & Files ---
data_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
file_list = dir(fullfile(data_dir, 'TF*_Theo.mat'));
save_path = fullfile(data_dir, 'Cohort_Processed.mat');

% --- Analysis Parameters ---
bin_ms = 50;                % Time bin size
n_components_reduced = 4;   % PCA components
n_window = -3:3;            % CCA Sliding window (+/- trials)
max_shift_bins = 5;         % Precession/Lag max shift
n_min_animals = 3;          % Min animals to analyze a pair

% --- Plotting Parameters ---
n_pad = 10;                 % Trials for Early/Pre/Post alignment
smooth_win = 3;             % Smoothing for simple correlation plots

% Load Learning Points (Prerequisite for Plotting)
if exist('period_experienced.mat', 'file')
    load('period_experienced.mat'); 
    learning_points = period_experienced(:, 1);
else
    error('File "period_experienced.mat" not found. Cannot proceed with alignment.');
end

% Initialize Master Structure
Cohort = struct();

%% 2. MAIN PRE-PROCESSING LOOP (Load, Bin, Simple Corr, PCA)
fprintf('Starting Pre-processing on %d files...\n', length(file_list));

for i_animal = 1:length(file_list)
    filename = file_list(i_animal).name;
    fullpath = fullfile(data_dir, filename);
    fprintf('Processing %d/%d: %s\n', i_animal, length(file_list), filename);
    
    % --- A. Load Data ---
    loaded_data = load(fullpath);
    if ~isfield(loaded_data, 'data_behaviour') || ~isfield(loaded_data, 'binned_spikes')
        warning('Skipping %s (Missing vars)', filename); continue;
    end
    
    beh = loaded_data.data_behaviour;
    spikes = loaded_data.binned_spikes;
    units = loaded_data.units;
    
    Cohort(i_animal).Filename = filename;
    Cohort(i_animal).Units = units;
    
    % --- B. Define Task Period & Binning ---
    task_start = find(~isnan(beh.trial_binned_cued), 1);
    task_end = find(~isnan(beh.trial_binned_cued), 1, 'last');
    
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
    
    % Initialize PCA Storage
    Cohort(i_animal).Reduced = struct();
    
    Trials = struct(); 
    
    % --- C. Trial Loop (Binning & Simple Corr) ---
    for itrial = 1:num_trials
        t_idx = trial_vec == itrial;
        if ~any(t_idx), continue; end
        
        % 1. Optimized Binning
        t_raw_spikes = raw_spikes_t(:, t_idx);
        n_samps = length(pos(t_idx));
        time_bins = ceil((1:n_samps) / bin_ms)';
        n_new = time_bins(end);
        
        % Vectorized Spike Binning
        unit_subs = repmat((1:n_units)', 1, n_samps);
        time_subs = repmat(time_bins', n_units, 1);
        binned_spikes = accumarray([unit_subs(:), time_subs(:)], t_raw_spikes(:), [n_units, n_new]);
        
        % Behavior Binning
        b_pos = accumarray(time_bins, pos(t_idx), [], @mean);
        b_licks = accumarray(time_bins, licks(t_idx), [], @sum);
        b_vel = accumarray(time_bins, vel(t_idx), [], @mean);
        
        Trials(itrial).binned_spikes = binned_spikes;
        Trials(itrial).binned_pos = b_pos;
        Trials(itrial).binned_licks = b_licks;
        Trials(itrial).binned_vel = b_vel;
        
        % 2. Simple Correlation (Area Mean vs Behavior)
        for iarea = 1:num_areas
            a_idx = units.idx(iarea, :);
            if any(a_idx)
                act = mean(binned_spikes(a_idx, :), 1)';
                Cohort(i_animal).Correlations.Velocity(iarea, itrial) = corr(act, b_vel, 'Rows', 'complete');
                Cohort(i_animal).Correlations.Licks(iarea, itrial) = corr(act, b_licks, 'Rows', 'complete');
            end
        end
    end
    Cohort(i_animal).Trials = Trials;
    
    % --- D. PCA Reduction (Per Area, Concatenated Trials) ---
    regions = unique(units.regions_label);
    for r = 1:length(regions)
        reg = regions{r};
        reg_idx = units.idx(strcmp(units.regions_label, reg), :);
        if ~any(reg_idx), continue; end
        
        % Concatenate trials
        all_data = [];
        lens = [];
        for t=1:num_trials
            if ~isempty(Trials(t).binned_spikes)
                d = Trials(t).binned_spikes(reg_idx, :);
                all_data = [all_data, d]; %#ok<AGROW>
                lens(t) = size(d, 2);
            else
                lens(t) = 0;
            end
        end
        
        % PCA
        if size(all_data, 1) >= n_components_reduced
            [~, score, ~] = pca(all_data', 'NumComponents', n_components_reduced);
            score = score'; % [PC x Time]
            
            % Distribute back
            curr = 1;
            for t=1:num_trials
                if lens(t)>0
                    Cohort(i_animal).Trials(t).Reduced.(reg) = score(:, curr:curr+lens(t)-1);
                    curr = curr + lens(t);
                end
            end
        end
    end
end
fprintf('Pre-processing complete.\n');

%% 3. REGION & PAIR DISCOVERY
% Collect all unique regions
all_regs = {};
animal_regions = cell(length(Cohort), 1);
for i = 1:length(Cohort)
    if isfield(Cohort(i).Units, 'regions_label')
        u_l = unique(Cohort(i).Units.regions_label);
        animal_regions{i} = u_l;
        all_regs = [all_regs, u_l]; %#ok<AGROW>
    end
end
unique_regions = unique(all_regs);
n_regs = length(unique_regions);

% Find Pairs shared by N animals
valid_pairs = cell(0, 2);
for i = 1:n_regs
    for j = i+1:n_regs
        cnt = 0;
        for k = 1:length(Cohort)
            if ismember(unique_regions{i}, animal_regions{k}) && ismember(unique_regions{j}, animal_regions{k})
                cnt = cnt + 1;
            end
        end
        if cnt >= n_min_animals
            valid_pairs(end+1, :) = {unique_regions{i}, unique_regions{j}}; %#ok<AGROW>
        end
    end
end
n_pairs = size(valid_pairs, 1);
fprintf('Found %d valid region pairs.\n', n_pairs);

%% 4. ANALYSIS: REGION-REGION CCA (Reduced & Full)
warning off; % Suppress rank deficiency warnings
GroupCCA = struct();

for ipair = 1:n_pairs
    pair_name = sprintf('%s-%s', valid_pairs{ipair, 1}, valid_pairs{ipair, 2});
    GroupCCA(ipair).Name = pair_name;
    GroupCCA(ipair).Reduced.TrialCorr = {}; GroupCCA(ipair).Reduced.Precession = {};
    GroupCCA(ipair).Full.TrialCorr = {};    GroupCCA(ipair).Full.Precession = {};
    GroupCCA(ipair).AnimalIndex = []; % Store ID immediately
end

fprintf('Running Region-Region CCA...\n');

for i_animal = 1:length(Cohort)
    n_trials = length(Cohort(i_animal).Trials);
    units = Cohort(i_animal).Units;
    if isempty(units), continue; end
    
    % Cache indices for Full Space
    reg_map = struct();
    for r=1:n_regs, reg = unique_regions{r}; reg_map.(reg) = units.idx(strcmp(units.regions_label, reg), :); end
    
    for ipair = 1:n_pairs
        a1 = valid_pairs{ipair, 1}; a2 = valid_pairs{ipair, 2};
        
        % Check data existence
        if ~isfield(Cohort(i_animal).Trials(1).Reduced, a1) || ~isfield(Cohort(i_animal).Trials(1).Reduced, a2)
            continue;
        end
        
        cca_rho = nan(2, n_trials); % 1=Red, 2=Full
        cca_shift = nan(2, n_trials, 2*max_shift_bins+1);
        
        for t = 1:n_trials
            % Sliding Window
            win = t + n_window; win = win(win>=1 & win<=n_trials);
            
            X_red = []; Y_red = []; X_full = []; Y_full = [];
            for w = win
                % Reduced
                X_red = [X_red, Cohort(i_animal).Trials(w).Reduced.(a1)]; %#ok<AGROW>
                Y_red = [Y_red, Cohort(i_animal).Trials(w).Reduced.(a2)]; %#ok<AGROW>
                % Full
                s = Cohort(i_animal).Trials(w).binned_spikes;
                if ~isempty(s)
                    X_full = [X_full, s(reg_map.(a1), :)]; %#ok<AGROW>
                    Y_full = [Y_full, s(reg_map.(a2), :)]; %#ok<AGROW>
                end
            end
            
            % --- REDUCED CCA ---
            if size(X_red, 2) > (n_components_reduced + max_shift_bins + 5)
                try
                    [~,~,r] = canoncorr(X_red', Y_red');
                    cca_rho(1, t) = r(1);
                    % Shifts
                    for s_i = 1:size(cca_shift,3)
                        s_val = s_i - (max_shift_bins + 1);
                        [Xc, Yc] = shift_and_cut(X_red, Y_red, s_val);
                        if size(Xc,2) > n_components_reduced
                            [~,~,r_l] = canoncorr(Xc', Yc');
                            cca_shift(1, t, s_i) = r_l(1);
                        end
                    end
                catch; end
            end
            
            % --- FULL CCA ---
            n_vars = size(X_full,1) + size(Y_full,1);
            if size(X_full, 2) > n_vars + 5
                try
                    [~,~,r] = canoncorr(X_full', Y_full');
                    cca_rho(2, t) = r(1);
                    % Shifts
                    for s_i = 1:size(cca_shift,3)
                        s_val = s_i - (max_shift_bins + 1);
                        [Xc, Yc] = shift_and_cut(X_full, Y_full, s_val);
                        if size(Xc,2) > n_vars
                            [~,~,r_l] = canoncorr(Xc', Yc');
                            cca_shift(2, t, s_i) = r_l(1);
                        end
                    end
                catch; end
            end
        end
        
        % Store Results
        if ~all(isnan(cca_rho(:)))
            GroupCCA(ipair).AnimalIndex(end+1) = i_animal;
            GroupCCA(ipair).Reduced.TrialCorr{end+1} = cca_rho(1, :);
            GroupCCA(ipair).Reduced.Precession{end+1} = squeeze(mean(cca_shift(1,:,:), 2, 'omitnan'));
            GroupCCA(ipair).Full.TrialCorr{end+1} = cca_rho(2, :);
            GroupCCA(ipair).Full.Precession{end+1} = squeeze(mean(cca_shift(2,:,:), 2, 'omitnan'));
        end
    end
end
fprintf('Region-Region CCA Done.\n');

%% 5. ANALYSIS: REGION-BEHAVIOR CCA (Full Space)
target_behavior = 'velocity'; 
BehavCCA = struct();
for i=1:n_regs, BehavCCA(i).Name = unique_regions{i}; BehavCCA(i).AnimalIndex=[]; BehavCCA(i).TrialCorr={}; end

fprintf('Running Region-Behavior CCA (%s)...\n', target_behavior);

for i_animal = 1:length(Cohort)
    n_trials = length(Cohort(i_animal).Trials);
    units = Cohort(i_animal).Units;
    if isempty(units), continue; end
    
    % Cache indices
    reg_map = struct();
    for r=1:n_regs, reg = unique_regions{r}; reg_map.(reg) = units.idx(strcmp(units.regions_label, reg), :); end
    
    for i_reg = 1:n_regs
        reg = unique_regions{i_reg};
        if ~isfield(Cohort(i_animal).Trials(1).Reduced, reg), continue; end
        
        idx_units = reg_map.(reg);
        n_u = sum(idx_units);
        cca_rho = nan(1, n_trials);
        
        for t = 1:n_trials
            win = t + n_window; win = win(win>=1 & win<=n_trials);
            
            X_pool = []; Y_pool = [];
            for w = win
                s = Cohort(i_animal).Trials(w).binned_spikes;
                if isempty(s), continue; end
                
                if strcmp(target_behavior, 'velocity')
                    b = Cohort(i_animal).Trials(w).binned_vel;
                else
                    b = Cohort(i_animal).Trials(w).binned_licks;
                end
                
                x = s(idx_units, :);
                min_l = min(size(x, 2), length(b));
                if min_l > 0
                    X_pool = [X_pool, x(:, 1:min_l)]; %#ok<AGROW>
                    Y_pool = [Y_pool, b(1:min_l)'];   %#ok<AGROW>
                end
            end
            
            % Full CCA Constraints
            if size(X_pool, 2) > (n_u + 5)
                try
                    [~,~,r] = canoncorr(X_pool', Y_pool');
                    cca_rho(t) = r(1);
                catch; end
            end
        end
        
        if ~all(isnan(cca_rho))
            BehavCCA(i_reg).AnimalIndex(end+1) = i_animal;
            BehavCCA(i_reg).TrialCorr{end+1} = cca_rho;
        end
    end
end
fprintf('Behavior CCA Done.\n');

%% 6. PLOTTING

% --- Plot 1: Simple Velocity Correlation (Aligned) ---
figure('Name', 'Simple Corr Aligned', 'Color', 'w');
t = tiledlayout('flow', 'TileSpacing', 'compact');
for r = 1:n_regs
    reg = unique_regions{r};
    
    % Collect Data from Cohort
    data_mat = nan(length(Cohort), 100); % Arbitrary max
    valid_id = [];
    for i=1:length(Cohort)
        u_l = Cohort(i).Units.regions_label;
        idx = find(strcmp(u_l, reg));
        if ~isempty(idx)
            tr = Cohort(i).Correlations.Velocity(idx, :);
            data_mat(i, 1:length(tr)) = tr;
            valid_id(end+1) = i; %#ok<AGROW>
        end
    end
    
    if isempty(valid_id), continue; end
    
    % Prepare Aligned Data
    [d_early, d_pre, d_post, n_valid] = align_data(data_mat(valid_id, :), learning_points(valid_id), n_pad);
    
    if n_valid > 0
        nexttile; hold on;
        plot_aligned_set(d_early, d_pre, d_post, [0.2 0.2 0.2], n_pad);
        title(sprintf('%s (n=%d)', reg, n_valid));
        ylim([-0.2 0.6]); yline(0, '--k');
    end
end
sgtitle('Simple Activity-Velocity Correlation');

% --- Plot 2: Region-Region CCA (Comparison Reduced vs Full) ---
for i = 1:n_pairs
    [r_mat, max_t] = pad_trials(GroupCCA(i).Reduced.TrialCorr);
    [f_mat, ~]     = pad_trials(GroupCCA(i).Full.TrialCorr);
    
    if isempty(r_mat), continue; end
    
    figure('Name', ['Comparison: ' GroupCCA(i).Name], 'Color', 'w', 'Position', [100 100 800 300]);
    subplot(1,2,1); hold on;
    plot_shaded(1:max_t, r_mat, 'b');
    plot_shaded(1:max_t, f_mat, 'r');
    title('Trial Evolution'); xlabel('Trial'); ylabel('CC1');
    legend('Reduced', 'Full');
    
    subplot(1,2,2); hold on;
    % Precession Plot
    r_prec = pad_shifts(GroupCCA(i).Reduced.Precession);
    f_prec = pad_shifts(GroupCCA(i).Full.Precession);
    lags = -max_shift_bins:max_shift_bins;
    plot_shaded(lags, r_prec, 'b');
    plot_shaded(lags, f_prec, 'r');
    title('Precession'); xlabel('Lag'); xline(0,'--k');
end

% --- Plot 3: Region-Region CCA (Aligned) ---
for i = 1:n_pairs
    idxs = GroupCCA(i).AnimalIndex;
    if isempty(idxs), continue; end
    
    [r_early, r_pre, r_post, ~] = align_data_cells(GroupCCA(i).Reduced.TrialCorr, learning_points(idxs), n_pad);
    [f_early, f_pre, f_post, n_v] = align_data_cells(GroupCCA(i).Full.TrialCorr, learning_points(idxs), n_pad);
    
    figure('Name', ['Aligned: ' GroupCCA(i).Name], 'Color', 'w', 'Position', [100 100 600 300]);
    hold on;
    plot_aligned_set(r_early, r_pre, r_post, [0 0.447 0.741], n_pad); % Blue
    plot_aligned_set(f_early, f_pre, f_post, [0.85 0.325 0.098], n_pad); % Red
    title(sprintf('%s (n=%d)', GroupCCA(i).Name, n_v));
    legend('Reduced', 'Full');
end

% --- Plot 4: Region-Behavior CCA (Full Space, Aligned) ---
figure('Name', 'Behavior CCA Aligned', 'Color', 'w');
t = tiledlayout('flow', 'TileSpacing', 'compact');
for i = 1:n_regs
    idxs = BehavCCA(i).AnimalIndex;
    if isempty(idxs), continue; end
    
    [d_early, d_pre, d_post, n_v] = align_data_cells(BehavCCA(i).TrialCorr, learning_points(idxs), n_pad);
    
    if n_v > 0
        nexttile; hold on;
        plot_aligned_set(d_early, d_pre, d_post, [0.85 0.325 0.098], n_pad);
        title(sprintf('%s (n=%d)', BehavCCA(i).Name, n_v));
        ylim([0 0.8]);
    end
end
sgtitle(['Full Space CCA vs ' target_behavior]);


%% 7. HELPER FUNCTIONS
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

function mat = pad_shifts(cell_arr)
    if isempty(cell_arr), mat=[]; return; end
    try mat = cell2mat(cell_arr'); catch, mat=nan; end
end

function plot_shaded(x, data, col)
    if isempty(data), return; end
    mu = mean(data,1,'omitnan');
    se = std(data,0,1,'omitnan')./sqrt(sum(~isnan(data),1));
    valid = ~isnan(mu);
    if ~any(valid), return; end
    fill([x(valid) fliplr(x(valid))], [mu(valid)+se(valid) fliplr(mu(valid)-se(valid))], ...
         col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(x(valid), mu(valid), 'Color', col, 'LineWidth', 2);
end

% --- Alignment Helpers ---
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