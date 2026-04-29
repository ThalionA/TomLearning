%% Load data
load Supermouse_for_Theo.mat

%% Preprocess
supermouse_all = permute(supermouse.spatial_firing_gf, [1, 3, 2]);

[num_units_total, num_bins, num_trials] = size(supermouse_all);

all_areas = cellstr(supermouse.units.region_id);

area_labels = cell(1, num_units_total);
is_valid = false(1, num_units_total);
for iUnit = 1:num_units_total
    [area_idx, mouse_idx] = ind2sub([numel(supermouse.units.region_id), numel(supermouse.units.mouse_id)], find(squeeze(supermouse.units.idx(:, :, iUnit)), 1));
    if isempty(area_idx)
        continue
    else
        area_labels{iUnit} = all_areas{area_idx};
        % if ~strcmp(area_labels{iUnit}, 'CA3') %& ~strcmp(area_labels{iUnit}, 'SCop') 
        %     continue
        % end
        is_valid(iUnit) = true;
    end
    
end

area_labels = area_labels(is_valid);

% Count units per area
[unique_areas, ~, area_indices_in_unique_list] = unique(area_labels);
unit_counts_per_area = accumarray(area_indices_in_unique_list, 1);
figure
bar(unique_areas, unit_counts_per_area)


supermouse_combined_valid = supermouse_all(is_valid, :, :);

%% Subsample units by area
fprintf('--- Subsampling units ---\n');

% Count units per area
[unique_areas, ~, area_indices_in_unique_list] = unique(area_labels);
unit_counts_per_area = accumarray(area_indices_in_unique_list, 1);

% --- Parameters for subsampling ---
min_unit_threshold = 20; 

% Identify areas to keep (those with enough units)
areas_to_keep_mask = unit_counts_per_area >= min_unit_threshold;
areas_to_keep = unique_areas(areas_to_keep_mask);
counts_of_areas_to_keep = unit_counts_per_area(areas_to_keep_mask);

% Identify and report skipped areas
areas_to_skip = unique_areas(~areas_to_keep_mask);
if ~isempty(areas_to_skip)
    fprintf('Skipping areas with fewer than %d units: %s\n', min_unit_threshold, strjoin(areas_to_skip, ', '));
end

% Determine the number of units to subsample to (minimum of the kept areas)
if isempty(areas_to_keep)
    error('No areas meet the minimum unit threshold of %d.', min_unit_threshold);
end
target_n_units = min(counts_of_areas_to_keep);
fprintf('Subsampling all kept areas to %d units.\n', target_n_units);
fprintf('Kept areas: %s\n', strjoin(areas_to_keep, ', '));

% Perform subsampling
final_indices_to_keep = [];
for i = 1:length(areas_to_keep)
    current_area = areas_to_keep{i};
    
    % Find indices of all units from the current area
    indices_for_this_area = find(strcmp(area_labels, current_area));
    
    % Randomly select 'target_n_units' from them without replacement
    shuffled_indices = indices_for_this_area(randperm(length(indices_for_this_area)));
    selected_indices = shuffled_indices(1:target_n_units);
    
    % Append to the final list of indices
    final_indices_to_keep = [final_indices_to_keep, selected_indices];
end

% Create the new balanced dataset
supermouse_balanced = supermouse_combined_valid(final_indices_to_keep, :, :);
area_labels_balanced = area_labels(final_indices_to_keep);

fprintf('Subsampling complete. New dataset has %d total units.\n\n', size(supermouse_balanced, 1));

%% Analyze and Visualize Response Range per Area
fprintf('--- Analyzing response ranges per area ---\n');

% Calculate the firing rate range for each neuron in the balanced dataset
num_balanced_units = size(supermouse_balanced, 1);
response_ranges = zeros(num_balanced_units, 1);

for iUnit = 1:num_balanced_units
    % Extract all firing data for the current unit across all bins and trials
    unit_data = squeeze(supermouse_balanced(iUnit, :, :));
    
    % Calculate the range (max value - min value) and store it
    response_ranges(iUnit) = max(unit_data(:)) - min(unit_data(:));
end

% Visualize the distribution of ranges for each area using a box plot
figure('Name', 'Response Range by Area');
boxplot(response_ranges, area_labels_balanced, 'Notch', 'off', 'Whisker', 1.5);

title('Distribution of Firing Rate Range by Brain Area');
xlabel('Brain Area');
ylabel('Response Range (Max Rate - Min Rate)');
grid on;
xtickangle(45); % Angle the x-axis labels to prevent overlap

fprintf('Displayed box plot of response ranges.\n\n');
%% TCA params
% --- TCA Parameters ---
cfg.tca.method = 'cp_nmu'; % 'cp_nmu', 'cp-als'
cfg.tca.normalization = 'min-max'; % 'none', 'z-score', 'min-max'
cfg.tca.max_factors = 10;
cfg.tca.max_iterations = 300;
cfg.tca.num_initialisations = 20;
cfg.tca.select_factors_method = 'bic'; % or 'manual' or 'fixed'
cfg.tca.fixed_n_factors = 4; % Used if select_factors_method is 'fixed' or as fallback/manual choice

cfg.plot.colors.area_map = containers.Map(all_areas, num2cell(lines(numel(all_areas)), 2));

cfg.plot.zone_params.bin_size = 2.5;
cfg.plot.zone_params.visual_zones_au  = [0,50; 100,125; 175,200; 250,275; 325,350; 400,425];
cfg.plot.zone_params.reward_zone_au   = [275 325];
cfg.plot.zone_params.corridor_end_au = 500;

cfg.plot.epoch_names = {'Naive', 'Intermediate', 'Expert'}; % Corresponds to cfg.control_epoch_windows order

cfg.control_epoch_windows = {1:10, [-10, -1], [1, 10]}; % {early_range, pre_lp_offset, post_lp_offset} relative to reference point




%% Run TCA
[best_mdl, best_n_factors, tca_results] = runTCAAnalysis(supermouse_combined_valid, cfg.tca);

%% Adjust N factors and Plot diagnostics

best_n_factors = 10;
best_mdl = tca_results.all_best_models{best_n_factors};

figure;
legend_items = {}; % Initialize cell array for legend entries
x_axis_values = 1:cfg.tca.max_factors;

% Plot BIC on the left axis
if isfield(tca_results, 'bic_values') && ~isempty(tca_results.bic_values) && ~all(isnan(tca_results.bic_values))
    plot(x_axis_values, tca_results.bic_values, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    legend_items{end+1} = 'BIC';
    hold on; % Hold on for other plots
end

right_axis_plotted = false; % Keep track if right axis is used

% Plot Normalized Reconstruction Error on right axis
if isfield(tca_results,'recon_errors') && ~isempty(tca_results.recon_errors) && ~all(isnan(tca_results.recon_errors))
    yyaxis right;
    valid_recon = tca_results.recon_errors(~isnan(tca_results.recon_errors));
    if ~isempty(valid_recon)
        norm_recon_err = tca_results.recon_errors ./ max(valid_recon);
        plot(x_axis_values, norm_recon_err, 'r-*', 'LineWidth', 1);
        ylabel('Normalized Recon Error');
        ax = gca; ax.YAxis(2).Color = 'r'; % Color right axis red
        legend_items{end+1} = 'Norm Recon Error';
        right_axis_plotted = true;
    else
        yyaxis left; % Switch back if no valid recon error plotted
    end
end

% Plot Initialization Similarity Score (also on right axis if used)
if isfield(tca_results,'init_similarity_scores') && ~isempty(tca_results.init_similarity_scores) && ~all(isnan(tca_results.init_similarity_scores))
    yyaxis right; % Ensure we are on the right axis
    plot(x_axis_values, tca_results.init_similarity_scores, 'g-s', 'LineWidth', 1, 'MarkerFaceColor', 'g');
    legend_items{end+1} = 'Avg Similarity';
    if right_axis_plotted % Already plotted recon error?
        ylabel('Norm Recon Error / Avg Similarity'); % Combined label
    else % Only similarity is on right axis
        ylabel('Avg Factor Similarity');
        ax = gca; ax.YAxis(2).Color = 'g'; % Color right axis green
    end
    right_axis_plotted = true; % Mark right axis as used
end

% Switch back to left axis for final annotations
yyaxis left;
xlabel('Number of Factors');
ylabel('BIC'); % Primary Y-label remains BIC
title('TCA Factor Selection Diagnostics');



if ~isempty(legend_items)
    legend(legend_items, 'Location', 'best');
end

% Add vertical line for selected factors
if ~isempty(best_n_factors)
    xline(best_n_factors, 'k--', sprintf('Selected = %d', best_n_factors), 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
end

box off;
grid on;
xlim([1, cfg.tca.max_factors + 0.5]); % Adjust xlim slightly
hold off;

%% Plotting
nFactors = size(best_mdl.U{1}, 2); % Use actual factors from model

plotSpatialFactors(best_mdl, cfg.plot.zone_params, nFactors);

plotTrialFactors(best_mdl, cfg, nFactors);

% plotNeuronFactorsByArea(best_mdl, area_labels', cfg, nFactors, false);

%% 1. Ensemble Assignment and Purity (Refined)
fprintf('--- Assigning Ensembles ---\n');

% Extract factors
neuron_factors = best_mdl.u{1}; % Neuron factors
[~, ensemble_assignments] = max(neuron_factors, [], 2);
total_ensembles = max(ensemble_assignments);

% Purity calculation
neuron_factors_norm = neuron_factors ./ sum(neuron_factors, 2);
neuron_purity = max(neuron_factors_norm, [], 2);
purity_threshold = 0.1;
is_pure = neuron_purity >= purity_threshold;

fprintf('Identified %d ensembles. %d/%d units classified as pure (threshold %.2f).\n', ...
    total_ensembles, sum(is_pure), numel(is_pure), purity_threshold);

%% 2. Ensemble Composition by Brain Area
% Maps which areas constitute each ensemble (ignoring neuron type)
fprintf('--- Analyzing Area Composition ---\n');

% Ensure area labels match the dataset size
% Assuming 'area_labels' is the cell array matching 'supermouse_combined_valid'
if numel(area_labels) ~= size(supermouse_combined_valid, 1)
    % Fallback if you used the 'balanced' set
    if exist('area_labels_balanced', 'var') && numel(area_labels_balanced) == size(supermouse_combined_valid, 1)
        current_area_labels = area_labels_balanced;
    else
        error('Mismatch between area label count and number of units in dataset.');
    end
else
    current_area_labels = area_labels;
end

unique_areas = unique(current_area_labels);
n_areas = numel(unique_areas);

% Preallocate
% Rows: Ensembles, Cols: Areas
composition_counts = zeros(total_ensembles, n_areas);

for iEns = 1:total_ensembles
    % Get indices for this ensemble (using only pure units for cleaner definition)
    idx = (ensemble_assignments == iEns) & is_pure;
    
    if sum(idx) == 0; continue; end
    
    these_areas = current_area_labels(idx);
    
    for iArea = 1:n_areas
        composition_counts(iEns, iArea) = sum(strcmp(these_areas, unique_areas{iArea}));
    end
end

% Convert to proportions
composition_prop = composition_counts ./ sum(composition_counts, 2) * 100;

% Plotting
figure('Name', 'Ensemble Area Composition', 'Position', [300, 300, 900, 500]);
b = bar(composition_prop, 'stacked');
xlabel('Ensemble ID');
ylabel('Proportion of Units (%)');
title('Anatomical Composition of Ensembles');
legend(unique_areas, 'Location', 'eastoutside');
xticks(1:total_ensembles);
grid on;

% Optional: Plot absolute counts to see ensemble sizes
figure('Name', 'Ensemble Sizes by Area');
bar(composition_counts, 'stacked');
xlabel('Ensemble ID');
ylabel('Number of Units');
title('Ensemble Size and Anatomical Make-up');
legend(unique_areas, 'Location', 'eastoutside');

%% 3. Spatial Activity Evolution by Epoch (Naive -> Expert)
fprintf('--- Visualizing Activity Evolution ---\n');

% Epoch definitions (from your config structure)
epochs = cfg.plot.epoch_names; % {'Naive', 'Intermediate', 'Expert'}
epoch_wins = cfg.control_epoch_windows; 
% Ensure epoch_wins indices are valid for this dataset
n_trials_total = size(supermouse_combined_valid, 3);
n_bins = size(supermouse_combined_valid, 2);

% Colors for epochs
epoch_colors = [0.2 0.2 0.8; 0.5 0.5 0.5; 0.8 0.2 0.2]; % Blue, Grey, Red

figure('Name', 'Ensemble Spatial Tuning by Epoch', 'Position', [100, 100, 1200, 800]);
t = tiledlayout('flow', 'TileSpacing', 'compact');

for iEns = 1:total_ensembles
    nexttile;
    hold on;
    
    % Get units in this ensemble
    idx_units = (ensemble_assignments == iEns) & is_pure;
    
    if sum(idx_units) < 2
        title(sprintf('Ens %d (n<2)', iEns));
        continue; 
    end
    
    % Loop through epochs
    h_lines = [];
    for iEp = 1:numel(epochs)
        % Handle trial indices relative to learning or fixed
        % This logic assumes the windows in cfg are relative to something or absolute.
        % For simplicity, assuming these are valid indices. If relative to learning
        % point, specific logic per mouse is needed (like in your 'lick error' section).
        % Here I define a simplified generic mapping for the aggregate:
        
        % Placeholder logic: divide total trials into 3 chunks if specific indices aren't easy
        % Or use the provided windows if they are absolute 1:10 etc.
        current_win = epoch_wins{iEp};
        
        % Handle negative indexing (relative to end or specific point)
        % Since we don't have the per-mouse learning point in this loop easily without
        % unwrapping the struct again, I will simplify to: 
        % 1: First 20 trials, 2: Middle 20, 3: Last 20.
        % ADJUST THIS BLOCK TO MATCH YOUR EXACT INDEXING LOGIC
        if iEp == 1
            tr_idx = 1:20; 
        elseif iEp == 2
            mid = floor(n_trials_total/2);
            tr_idx = mid-10 : mid+9;
        else
            tr_idx = n_trials_total-19 : n_trials_total;
        end
        
        tr_idx = tr_idx(tr_idx > 0 & tr_idx <= n_trials_total);
        
        % Extract data: Units x Bins x Trials -> Mean over Units -> Mean over Trials
        ens_data = supermouse_combined_valid(idx_units, :, tr_idx);
        
        % Mean spatial profile for the ensemble in this epoch
        % 1. Mean across trials
        mean_trial_resp = squeeze(mean(ens_data, 3, 'omitnan')); 
        % 2. Mean across units (with SEM)
        mu = mean(mean_trial_resp, 1, 'omitnan');
        se = std(mean_trial_resp, 0, 1, 'omitnan') / sqrt(sum(idx_units));
        
        h = shadedErrorBar(1:n_bins, mu, se, 'lineprops', {'Color', epoch_colors(iEp,:), 'LineWidth', 1.5});
        h_lines(end+1) = h.mainLine;
    end
    
    title(sprintf('Ensemble %d (n=%d)', iEns, sum(idx_units)));
    axis tight;
    if iEns == 1
        legend(h_lines, epochs, 'Location', 'best');
    end
    grid on;
end

xlabel(t, 'Spatial Bin');
ylabel(t, 'Firing Rate (z-score/norm)');

%% 4. Canonical Correlation Analysis (CCA) between Ensembles
% This analyzes how coupled the population dynamics of different ensembles are.
fprintf('--- Running CCA between Ensembles ---\n');

cca_rho_matrix = nan(total_ensembles, total_ensembles);
cca_pval_matrix = nan(total_ensembles, total_ensembles);

% Data reshaping for CCA:
% CCA requires (Samples x Features).
% Features = Neurons in the ensemble.
% Samples = Timepoints (Bins * Trials).
% We want to see if the temporal trajectory of Ensemble A predicts Ensemble B.

% Pre-format data structures to speed up the loop
ensemble_data_flat = cell(1, total_ensembles);

for iEns = 1:total_ensembles
    idx_units = (ensemble_assignments == iEns) & is_pure;
    
    if sum(idx_units) == 0
        ensemble_data_flat{iEns} = [];
        continue;
    end
    
    % Extract: Units x Bins x Trials
    raw = supermouse_combined_valid(idx_units, :, :);
    [nU, nB, nT] = size(raw);
    
    % Reshape to (Bins*Trials) x Units
    % Permute to Bins x Trials x Units first to keep time continuous
    raw_perm = permute(raw, [2, 3, 1]); 
    flat_data = reshape(raw_perm, nB*nT, nU);
    
    % Remove NaN samples (if any bins are NaN across all, usually standardizing)
    % For CCA, NaNs propagate. Replace with 0 or remove rows. 
    % Here we assume Z-scored data, so we fill 0 or remove rows with NaN.
    flat_data(isnan(flat_data)) = 0; 
    
    ensemble_data_flat{iEns} = flat_data;
end

% Pairwise CCA
for i = 1:total_ensembles
    for j = 1:total_ensembles
        if i == j
            cca_rho_matrix(i,j) = 1; 
            continue; 
        end
        
        X = ensemble_data_flat{i};
        Y = ensemble_data_flat{j};
        
        % Check if we have enough data (Samples > Features)
        if isempty(X) || isempty(Y)
            continue; 
        end
        
        % We need min(rank(X), rank(Y)) canonical variables. 
        % We usually care about the FIRST canonical correlation (max shared variance).
        try
            [A, B, r, U, V, stats] = canoncorr(X, Y);
            if ~isempty(r)
                cca_rho_matrix(i,j) = r(1); % The strongest correlation
                cca_pval_matrix(i,j) = stats.p(1);
            end
        catch
            % CCA can fail if matrices are ill-conditioned (rank deficient)
            % Usually happens if n_units > n_timepoints (unlikely here)
            warning('CCA failed for pair %d-%d', i, j);
        end
    end
end

% Visualization
figure('Name', 'Ensemble CCA Connectivity');
t = tiledlayout(1, 2);

% --- 1. Heatmap ---
nexttile;
imagesc(cca_rho_matrix);
colormap('jet');
colorbar;
clim([0 1]);
title('Canonical Correlation (r)');
xlabel('Ensemble ID');
ylabel('Ensemble ID');
axis square;
set(gca, 'XTick', 1:total_ensembles, 'YTick', 1:total_ensembles);

% --- 2. Network Graph (Weighted Edges) ---
nexttile;

% Thresholding
cca_threshold = 0.5; % Adjust this if too many/few edges appear
adjacency = cca_rho_matrix;

% Force symmetry to fix 'graph' error
adjacency = (adjacency + adjacency') / 2;

% Remove weak links and self-loops
adjacency(adjacency < cca_threshold) = 0;
adjacency(eye(size(adjacency))==1) = 0; 

% Create graph object
G = graph(adjacency);

if numedges(G) == 0
    text(0.5, 0.5, 'No correlations above threshold', 'HorizontalAlignment', 'center');
    axis off;
else
    % --- SCALING LOGIC ---
    weights = G.Edges.Weight;
    min_width = 1;
    max_width = 8;
    
    if max(weights) == min(weights)
        % If all edges have identical weight, make them uniform size
        LWidths = 4 * ones(size(weights));
    else
        % Normalize weights to range [0, 1] then map to [min_width, max_width]
        norm_weights = (weights - min(weights)) / (max(weights) - min(weights));
        LWidths = min_width + norm_weights * (max_width - min_width);
    end
    
    % Plot
    p = plot(G, 'Layout', 'circle', ...
            'NodeLabel', 1:total_ensembles, ...
            'NodeColor', 'k', ...
            'NodeFontSize', 12, ...
            'EdgeColor', [0.4 0.4 0.4]);
            
    % Apply calculated widths
    p.LineWidth = LWidths;
    
    title(sprintf('Strong Couplings (r > %.1f)', cca_threshold));
    axis off; 
end