%% Standalone_Dimensionality_Analysis.m
clear; clc; close all;

% --- Paths & Parameters ---
base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir = fullfile(base_dir, data_subfolder);
file_pattern = 'TF*_export.mat';

min_units_per_region = 5;
n_bins = 200;
target_variance = 90; % Evaluate PCs required to hit this variance threshold

file_list = dir(fullfile(data_dir, file_pattern));
n_animals = length(file_list);

dim_data = struct();
top_units_activity = struct(); % New struct for top unit activity

fprintf('Running standalone PCA dimensionality and top unit extraction...\n');

% --- Data Extraction Loop ---
for ianimal = 1:n_animals
    filename = file_list(ianimal).name;
    fullpath = fullfile(data_dir, filename);
    
    try
        D = load(fullpath);
        if ~isfield(D, 'units') || ~isfield(D, 'analysis_spatial')
            fprintf('  Skipping animal %d (%s): Missing required fields\n', ianimal, filename);
            continue; 
        end
        units = D.units;
        
        % Filter units (matching sweep logic)
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
             [~, curr_n_bins, ~] = size(animal_data);
             if curr_n_bins ~= n_bins
                 fprintf('  Skipping animal %d (%s): Incorrect bin size\n', ianimal, filename);
                 continue; 
             end
        else
            fprintf('  Skipping animal %d (%s): Missing spatial firing data\n', ianimal, filename);
            continue; 
        end
        
        animal_areas = unique(units.regions_label);
        
        for ia = 1:length(animal_areas)
            area_name = animal_areas{ia};
            u_idx = units.idx(strcmp(units.regions_label, area_name), :);
            num_valid_units = sum(u_idx);
            
            if num_valid_units < min_units_per_region, continue; end
            
            area_dat = animal_data(u_idx, :, :); % [units, bins, trials]
            reshaped_dat = reshape(area_dat, num_valid_units, [])'; 
            
            if size(reshaped_dat, 2) >= 2 
                % Compute PCA (extract coeffs this time)
                [coeffs, ~, ~, ~, explained] = pca(reshaped_dat);
                cum_var = cumsum(explained);
                n_comps = find(cum_var >= target_variance, 1);
                if isempty(n_comps), n_comps = size(reshaped_dat, 2); end
                
                % --- 1. Track Dimensionality ---
                if ~isfield(dim_data, area_name)
                    dim_data.(area_name).n_comps = [];
                    dim_data.(area_name).n_units = [];
                    dim_data.(area_name).norm_dim = [];
                end
                
                dim_data.(area_name).n_comps(end+1) = n_comps;
                dim_data.(area_name).n_units(end+1) = num_valid_units;
                dim_data.(area_name).norm_dim(end+1) = n_comps / num_valid_units;
                
                % --- 2. Track Top 3 Units Activity ---
                % Find top 3 units based on absolute PC1 loadings
                [~, sorted_idx] = sort(abs(coeffs(:, 1)), 'descend');
                n_top = min(3, length(sorted_idx));
                top3_idx = sorted_idx(1:n_top);
                
                % Get spatial activity: trial-average, then average across the 3 units
                top3_dat = area_dat(top3_idx, :, :); % [3, 200, trials]
                top3_trial_avg = mean(top3_dat, 3, 'omitnan'); % [3, 200]
                top3_overall_avg = mean(top3_trial_avg, 1, 'omitnan'); % [1, 200]
                
                if ~isfield(top_units_activity, area_name)
                    top_units_activity.(area_name) = [];
                end
                % Store the 1x200 spatial curve for this animal/area
                top_units_activity.(area_name)(end+1, :) = top3_overall_avg;
            end
        end
        
    catch ME
        fprintf('  Failed animal %d (%s): %s\n', ianimal, filename, ME.message);
    end
    
    fprintf('Finished processing animal %d/%d (%s)\n', ianimal, n_animals, filename);
end

% =========================================================================
% --- Plotting 1: Normalized Dimensionality Summary ---
% =========================================================================
figure('Name', 'Normalized Dimensionality Summary', 'Color', 'w', 'Position', [200 200 800 500]);
hold on;
areas_tracked = fieldnames(dim_data);
n_areas = length(areas_tracked);

means = zeros(n_areas, 1);
sems  = zeros(n_areas, 1);

for ia = 1:n_areas
    vals = dim_data.(areas_tracked{ia}).norm_dim;
    if ~isempty(vals)
        means(ia) = mean(vals);
        sems(ia)  = std(vals) / sqrt(length(vals));
    end
end

bar(1:n_areas, means, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1);

for ia = 1:n_areas
    vals = dim_data.(areas_tracked{ia}).norm_dim;
    if isempty(vals), continue; end
    x_jitter = ia + (rand(length(vals), 1) - 0.5) * 0.25;
    scatter(x_jitter, vals, 35, 'k', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'w');
end

errorbar(1:n_areas, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);

xticks(1:n_areas); xticklabels(areas_tracked);
ylabel(sprintf('Normalized Dimensionality\n(PCs for %d%% Var / Total Units)', target_variance));
title('Intrinsic Dimensionality Scaled by Unit Yield');
xlim([0.4, n_areas + 0.6]);
ylim([0 max(means + sems)*1.2]); 
box off;

% =========================================================================
% --- Plotting 2: Top Units Average Spatial Activity ---
% =========================================================================
figure('Name', 'Top Units Spatial Activity', 'Color', 'w', 'Position', [250 250 1200 800]);
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for ia = 1:n_areas
    area = areas_tracked{ia};
    data = top_units_activity.(area); % [n_valid_animals x n_bins]
    
    if isempty(data) || size(data, 1) < 2, continue; end
    
    % Calculate cross-animal mean and SEM
    mu = mean(data, 1, 'omitnan');
    se = std(data, 0, 1, 'omitnan') ./ sqrt(size(data, 1));
    
    nexttile;
    hold on;
    x = 1:n_bins;
    
    % Plot shaded SEM
    fill([x, fliplr(x)], [mu+se, fliplr(mu-se)], [0.85 0.33 0.1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Plot mean line
    plot(x, mu, 'Color', [0.85 0.33 0.1], 'LineWidth', 2);
    
    title(sprintf('%s (n=%d animals)', area, size(data, 1)));
    xlabel('Spatial Bin');
    ylabel('Norm. Firing Rate');
    xlim([1 n_bins]);
    box off;
end

fprintf('\nExtraction and plotting complete.\n');