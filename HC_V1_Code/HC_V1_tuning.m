%% Spatial_Tuning_Evolution_Epochs.m
% 1. Aligns data to Learning Points (Early, Pre, Post).
% 2. Quantifies % significant bins (BOTH DIRECTIONS: < 2.5th or > 97.5th).
% 3. Classifies units as "Tuned" if > 5% of their bins are significant.
% 4. Post-hoc selection of the best examples, plotting Mean +/- SEM of shuffles.

%% 1. SETUP & PARAMETERS
% clear; clc; close all;

% --- Paths ---
base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir = fullfile(base_dir, data_subfolder);
file_pattern = 'TF*_export.mat';
learning_file = 'animal_behaviour.mat';

% --- Analysis Parameters ---
n_shuffles = 100;
shuffle_method = 'permute'; % Options: 'permute' OR 'circshift'
prc_lower = 2.5;  % Two-tailed lower bound
prc_upper = 97.5; % Two-tailed upper bound
tuning_threshold_pct = 5.0; % >5% of bins must be significant

n_bins = 200;
min_units_per_region = 5;
target_areas = {'CA1', 'CA3', 'DG', 'SUB', 'V1', 'RSC'};
n_areas = length(target_areas);

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
    error('Learning file not found. Epoch alignment requires LP data.');
end

if length(learning_points) < n_animals
    learning_points(end+1:n_animals) = nan;
end

is_learner = ~isnan(learning_points);
mean_lp = round(mean(learning_points(is_learner)));
analysis_lp = learning_points;
analysis_lp(~is_learner) = mean_lp; 

%% 3. INITIALIZE RESULTS STRUCTURE
results_pct_bins = struct();
results_pct_tuned = struct();
for a = 1:n_areas
    results_pct_bins.(target_areas{a}) = nan(n_animals, 30);
    results_pct_tuned.(target_areas{a}) = nan(n_animals, 30);
end

% Master pool to store all tuned units for post-hoc "best of" selection
% (Note: se_shuff removed, as we are now plotting the full 95% CI bounds directly)
candidate_pool = struct('v', {}, 'mu_shuff', {}, 'p_lo', {}, 'p_hi', {}, 'area', {}, 'trial', {}, 'pct_sig', {}, 'animal', {});

%% 4. MAIN LOOP
fprintf('\nStarting single-unit spatial significance analysis (Absolute Z-Threshold Contiguous Field Method)...\n');

for ianimal = 1:n_animals
    filename = file_list(ianimal).name;
    fullpath = fullfile(data_dir, filename);
    fprintf('Processing Animal %d/%d: %s\n', ianimal, n_animals, filename);
    
    lp = analysis_lp(ianimal);
    
    try
        D = load(fullpath);
        if ~isfield(D, 'units') || ~isfield(D, 'analysis_spatial'), continue; end
        
        if isfield(D.analysis_spatial, 'firing') && isfield(D.analysis_spatial.firing, 'cued')
            raw_spatial = D.analysis_spatial.firing.cued.freq_z;
            animal_data = permute(raw_spatial, [1, 3, 2]); % [cells, bins, trials]
            [~, curr_n_bins, num_trials] = size(animal_data);
            
            if curr_n_bins ~= n_bins || isnan(lp) || lp < 11 || (lp + 9) > num_trials
                continue; 
            end
        else
            continue;
        end
        
        % Target Epoch Extraction
        idx_early = 1:10;
        idx_pre   = (lp - 10) : (lp - 1);
        idx_post  = lp : (lp + 9);
        target_trials = [idx_early, idx_pre, idx_post]; % 30 trials total
        
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
        
        animal_areas = unique(units.regions_label);
        
        for ia = 1:n_areas
            area_name = target_areas{ia};
            if ~ismember(area_name, animal_areas), continue; end
            
            u_idx = units.idx(strcmp(units.regions_label, area_name), :);
            if sum(u_idx) < min_units_per_region, continue; end
            
            area_dat = animal_data(u_idx, :, target_trials); % [units x bins x 30]
            n_units = size(area_dat, 1);
            
            unit_field_width = nan(n_units, 30);
            unit_is_tuned = nan(n_units, 30);
            
            for u = 1:n_units
                for t_idx = 1:30
                    v = squeeze(area_dat(u, :, t_idx)); % 1 x n_bins
                    valid_bins = sum(~isnan(v));
                    
                    if valid_bins == 0 || var(v, 0, 'omitnan') == 0
                        continue;
                    end
                    
                    % 1. Define active bins using an ABSOLUTE Z-score threshold
                    z_active_threshold = 1.0; 
                    real_active_mask = v > z_active_threshold;
                    
                    % 2. Find longest contiguous sequence of active bins
                    cc_real = bwconncomp(real_active_mask);
                    if cc_real.NumObjects == 0
                        real_max_width = 0;
                    else
                        real_max_width = max(cellfun(@length, cc_real.PixelIdxList));
                    end
                    
                    % 3. Build Null Distribution of Contiguous Widths via Shuffling
                    shuff_max_widths = zeros(n_shuffles, 1);
                    for s = 1:n_shuffles
                        v_shuff = v(randperm(n_bins));
                        cc_shuff = bwconncomp(v_shuff > z_active_threshold);
                        if cc_shuff.NumObjects > 0
                            shuff_max_widths(s) = max(cellfun(@length, cc_shuff.PixelIdxList));
                        end
                    end
                    
                    % 4. Tuning Classification
                    p95_shuff_width = prctile(shuff_max_widths, 95);
                    
                    % Must beat the shuffle AND have a minimum meaningful width
                    is_tuned_flag = (real_max_width > p95_shuff_width) && (real_max_width > 2);
                    
                    unit_field_width(u, t_idx) = real_max_width;
                    unit_is_tuned(u, t_idx) = is_tuned_flag;
                    
                    % Save tuned units for post-hoc selection
                    if is_tuned_flag
                        candidate_pool(end+1).v = v;
                        candidate_pool(end).active_thresh = z_active_threshold;
                        candidate_pool(end).real_width = real_max_width;
                        candidate_pool(end).shuff_thresh = p95_shuff_width;
                        candidate_pool(end).area = area_name;
                        candidate_pool(end).trial = target_trials(t_idx);
                        candidate_pool(end).animal = ianimal;
                    end
                end
            end
            
            % Save animal means
            results_pct_bins.(area_name)(ianimal, :) = mean(unit_field_width, 1, 'omitnan'); 
            results_pct_tuned.(area_name)(ianimal, :) = mean(unit_is_tuned, 1, 'omitnan') * 100;
        end
    catch ME
        fprintf('Error processing %s: %s\n', filename, ME.message);
    end
end

fprintf('\nAnalysis complete. Selecting best examples and plotting...\n');

%% 5. POST-HOC EXAMPLE SELECTION & PLOTTING
if ~isempty(candidate_pool)
    % Sort by largest contiguous field width to find the best examples
    all_widths = [candidate_pool.real_width];
    [~, sort_idx] = sort(all_widths, 'descend');
    sorted_pool = candidate_pool(sort_idx);
    
    [~, unique_idx] = unique({sorted_pool.area}, 'stable');
    best_examples = sorted_pool(unique_idx);
    
    n_examples = length(best_examples);
    
    figure('Name', 'Best Example Units (Contiguous Field Method)', 'Color', 'w', 'Position', [100, 100, 1200, 800]);
    t = tiledlayout('flow', 'TileSpacing', 'compact');
    
    for i = 1:n_examples
        nexttile; hold on;
        ex = best_examples(i); 
        
        v = ex.v;
        x = 1:n_bins;
        
        % Plot Real Data
        plot(x, v, 'k', 'LineWidth', 1.5);
        
        % Highlight active threshold
        yline(ex.active_thresh, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        
        % Highlight the active bins in red
        sig_bins = v > ex.active_thresh;
        plot(x(sig_bins), v(sig_bins), 'r.', 'MarkerSize', 12);
        
        title(sprintf('Area: %s | Animal: %d | Trial: %d\nMax Field Width: %d bins (Chance limit: %.1f)', ...
            ex.area, ex.animal, ex.trial, ex.real_width, ex.shuff_thresh));
        xlabel('Spatial Bin'); ylabel('Z-Scored Firing Rate');
        xlim([1 n_bins]);
        
        y_min = min(v); y_max = max(v);
        y_range = y_max - y_min;
        if y_range == 0, y_range = 1; end
        ylim([y_min - 0.1*y_range, y_max + 0.1*y_range]);
    end
    title(t, 'Top Spatially Tuned Units by Contiguous Field Width');
end
%% 6. PLOTTING HELPER
plot_3_epochs = @(data_struct, y_label, fig_title) build_epoch_plot(data_struct, target_areas, n_areas, is_learner, y_label, fig_title);

plot_3_epochs(results_pct_bins, '% Significant Bins', 'Evolution of Significant Bins');
plot_3_epochs(results_pct_tuned, '% Tuned Units', 'Evolution of Spatially Tuned Units');

%% LOCAL FUNCTIONS
function build_epoch_plot(data_struct, areas, n_areas, is_learner, y_label, fig_title)
    figure('Name', fig_title, 'Color', 'w', 'Position', [150, 150, 1600, 900]);
    tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    for a = 1:n_areas
        area_name = areas{a};
        data_mat = data_struct.(area_name); 
        
        nexttile; hold on;
        
        valid_mask = sum(~isnan(data_mat), 2) > 0;
        if sum(valid_mask) == 0
            title(sprintf('%s (No Data)', area_name)); continue;
        end
        
        idx_L = valid_mask & is_learner;
        idx_NL = valid_mask & ~is_learner;
        
        if sum(idx_NL) > 0
            mu_NL = mean(data_mat(idx_NL, :), 1, 'omitnan');
            se_NL = std(data_mat(idx_NL, :), 0, 1, 'omitnan') ./ sqrt(sum(idx_NL));
            plot_epoch_line(1:10, mu_NL(1:10), se_NL(1:10), [0.5 0.5 0.5]);
            plot_epoch_line(11:20, mu_NL(11:20), se_NL(11:20), [0.5 0.5 0.5]);
            plot_epoch_line(21:30, mu_NL(21:30), se_NL(21:30), [0.5 0.5 0.5]);
        end
        
        if sum(idx_L) > 0
            mu_L = mean(data_mat(idx_L, :), 1, 'omitnan');
            se_L = std(data_mat(idx_L, :), 0, 1, 'omitnan') ./ sqrt(sum(idx_L));
            plot_epoch_line(1:10, mu_L(1:10), se_L(1:10), [0.85 0.33 0.1]);
            plot_epoch_line(11:20, mu_L(11:20), se_L(11:20), [0.85 0.33 0.1]);
            plot_epoch_line(21:30, mu_L(21:30), se_L(21:30), [0.85 0.33 0.1]);
        end
        
        xline(10.5, 'k:'); xline(20.5, 'k:');
        xticks([5, 15, 25]); xticklabels({'Early', 'Pre', 'Post'});
        xlim([1 30]);
        ylabel(y_label);
        title(sprintf('%s (L=%d, NL=%d)', area_name, sum(idx_L), sum(idx_NL)));
        
        y_max = max(mean(data_mat(valid_mask,:), 1, 'omitnan')) * 1.5;
        if isnan(y_max) || y_max == 0, y_max = 10; end
        ylim([0, y_max]);
    end
end

function plot_epoch_line(x_range, mu, se, color)
    fill([x_range, fliplr(x_range)], [mu+se, fliplr(mu-se)], color, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(x_range, mu, '-', 'Color', color, 'LineWidth', 2);
end