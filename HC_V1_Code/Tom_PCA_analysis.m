animal_id = 'TF077';

load([animal_id '_export_colab.mat'])

col_naive = [0.466, 0.674, 0.188];
col_experienced = [0.494, 0.184, 0.556];

%% Setup and PCA per region

% Assume spatialFiring_gf is already defined
[numUnits, numTrials, numBins] = size(spatialFiring_gf);
numComponents = 8;
numRegions = size(idx_units, 1);  % idx_units is logical: [numRegions x numUnits]
components_actually_used = nan(1, numRegions);

% Define a color map for the regions (one color per region)
cmap_regions = lines(numRegions);

% Preallocate cell arrays to store results per region
regionCumsum = cell(numRegions, 1);
effectiveDim = zeros(numRegions, 1);
regionMeanTraj = cell(numRegions, 1);

for r = 1:numRegions
    % Get indices for the current region's units
    region_units = find(idx_units(r, :));
    if isempty(region_units)
        warning('No units found for region %s', regions_label(r));
        continue;
    end
    
    % Subset the data for this region
    regionFiring = spatialFiring_gf_z(region_units, :, :);
    [nUnitsRegion, numTrials, numBins] = size(regionFiring);
    
    % Rearrange dimensions: units x bins x trials, then reshape to units x (bins*numTrials)
    activity = permute(regionFiring, [1, 3, 2]);
    activity_reshaped = reshape(activity, nUnitsRegion, numBins*numTrials);
    
    % Run PCA (observations in rows)
    components_actually_used(r) = min([numComponents, numel(region_units)-1]);
    [~, region_activity_reduced, ~, ~, explained, ~] = pca(activity_reshaped', 'NumComponents', components_actually_used(r));
    
    % Calculate cumulative explained variance
    cumsum_explained = cumsum(explained);
    regionCumsum{r} = cumsum_explained;
    
    % Effective dimensionality: (# components to reach 90% variance)/(# units)
    idx_90 = find(cumsum_explained >= 90, 1);
    effectiveDim(r) = idx_90 / nUnitsRegion;
    
    % Reshape PCA scores to: numBins x numTrials x numComponents
    score_reshaped = reshape(region_activity_reduced, [numBins, numTrials, components_actually_used(r)]);
    % Compute mean trajectory (average across trials)
    mean_traj = squeeze(mean(score_reshaped, 2, 'omitnan'));
    regionMeanTraj{r} = mean_traj;
end

% Figure 1: Cumulative Explained Variance (one subplot per region)
figure;
for r = 1:numRegions
    subplot(numRegions, 1, r);
    cumsum_explained = regionCumsum{r};
    if isempty(cumsum_explained)
        continue;
    end
    plot(cumsum_explained, 'LineWidth', 1, 'Color', cmap_regions(r, :));
    hold on;
    % Mark the 90% threshold
    idx_90 = find(cumsum_explained >= 90, 1);
    plot([idx_90, idx_90], [0, cumsum_explained(idx_90)], '--', 'Color', [0.6 0.6 0.6]);
    plot([0, idx_90], [90, 90], '--', 'Color', [0.6 0.6 0.6]);
    axis tight;
    ylabel('Explained Variance (%)');
    title(sprintf('%s: EffDim = %.3f', regions_label(r), effectiveDim(r)));
end
xlabel('Component #');

% Figure 2: Average Trajectories (one subplot per region)
figure;
t = tiledlayout('flow', 'TileSpacing', 'compact');
for r = 1:numRegions
    nexttile
    % Adjust the subplot arrangement as desired
    mean_traj = regionMeanTraj{r};
    if isempty(mean_traj)
        continue;
    end
    if numComponents >= 3
        plot3(mean_traj(:,1), mean_traj(:,2), mean_traj(:,3), 'Color', cmap_regions(r, :), 'LineWidth', 2);
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
        grid on;
        title(sprintf('Avg Traj: %s', regions_label(r)));
        view(3);
    else
        plot(mean_traj(:,1), mean_traj(:,2), 'Color', cmap_regions(r, :), 'LineWidth', 2);
        xlabel('PC1'); ylabel('PC2');
        title(sprintf('Avg Traj: %s', regions_label(r)));
    end
end

% Figure 3: Projections on each PCA Dimension (tiles: numRegions x numComponents)
figure;
for r = 1:numRegions
    mean_traj = regionMeanTraj{r};  % size: [numBins x numComponents]
    for comp = 1:components_actually_used(r)
        subplot(numRegions, numComponents, (r-1)*numComponents + comp);
        % Plot the mean projection for this region and component
        plot(mean_traj(:, comp), 'Color', cmap_regions(r, :), 'LineWidth', 2);
        xlabel('Bin'); ylabel(sprintf('PC%d', comp));
        title(sprintf('%s - PC%d', regions_label(r), comp));
        hold on;
        % Lightly shade the reward zone:
        % Assume idx_reward_zone is a vector of bin indices that are in the reward zone
        yl = ylim;
        x_patch = [idx_reward_zone(1), idx_reward_zone(end), idx_reward_zone(end), idx_reward_zone(1)];
        y_patch = [yl(1), yl(1), yl(2), yl(2)];
        patch(x_patch, y_patch, [0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        uistack(gca, 'bottom');  % Ensure the patch is behind the line plot
    end
end


%% PCA and Data Processing per Region for Two Trial Epochs
% Assumes the following variables are defined:
%   spatialFiring_gf, idx_units, regions_label, period_naive, period_experienced, idx_reward_zone

[numUnits, numTrials, numBins] = size(spatialFiring_gf);
numComponents = 8;
numRegions = size(idx_units, 1);  % idx_units is logical: [numRegions x numUnits]
components_actually_used = nan(1, numRegions);

% Define a color map for the regions (one color per region)
cmap_regions = lines(numRegions);

% Preallocate cell arrays to store results per region
regionCumsum = cell(numRegions, 1);
effectiveDim = zeros(numRegions, 1);
regionMeanTraj = cell(numRegions, 1);

for r = 1:numRegions
    % Get indices for the current region's units
    region_units = find(idx_units(r, :));
    if isempty(region_units)
        warning('No units found for region %s', regions_label(r));
        continue;
    end
    
    % Subset the data for this region
    regionFiring = spatialFiring_gf_z(region_units, :, :);
    [nUnitsRegion, numTrials, numBins] = size(regionFiring);
    
    % Rearrange dimensions: units x bins x trials
    activity = permute(regionFiring, [1, 3, 2]);
    % Reshape to: (nUnitsRegion) x (numBins*numTrials)
    activity_reshaped = reshape(activity, nUnitsRegion, numBins*numTrials);
    
    % Run PCA (observations in rows)
    components_actually_used(r) = min([numComponents, numel(region_units)-1]);

    [~, region_activity_reduced, ~, ~, explained, ~] = pca(activity_reshaped', 'NumComponents', components_actually_used(r));
    
    % Calculate cumulative explained variance
    cumsum_explained = cumsum(explained);
    regionCumsum{r} = cumsum_explained;
    
    % Compute effective dimensionality:
    % (# components required to reach 90% variance) / (# units in the region)
    idx_90 = find(cumsum_explained >= 90, 1);
    effectiveDim(r) = idx_90 / nUnitsRegion;
    
    % Reshape PCA scores to: numBins x numTrials x numComponents
    score_reshaped = reshape(region_activity_reduced, [numBins, numTrials, components_actually_used(r)]);
    
    % Compute mean trajectories for each epoch:
    mean_traj_naive = squeeze(mean(score_reshaped(:, period_naive, :), 2, 'omitnan'));
    mean_traj_experienced = squeeze(mean(score_reshaped(:, period_experienced, :), 2, 'omitnan'));
    
    % Store in a structure for this region
    regionMeanTraj{r}.naive = mean_traj_naive;
    regionMeanTraj{r}.experienced = mean_traj_experienced;
end


% Figure 2: Average Trajectories per Region (Naive vs. Experienced)
% The naive trajectories are plotted in a paler version of the region’s color.
figure;
t = tiledlayout('flow', 'TileSpacing', 'compact');
for r = 1:numRegions
    nexttile;
    % Check that regionMeanTraj for this region exists.
    if isempty(regionMeanTraj{r})
        continue;
    end
    
    % Retrieve the mean trajectories for the two epochs.
    mean_traj_naive = regionMeanTraj{r}.naive;
    mean_traj_experienced = regionMeanTraj{r}.experienced;
    
    % Get the region's color and compute a paler version.
    col = cmap_regions(r, :);
    pale_col = col + (1 - col) * 0.5;  % 50% blend with white
    
    % Plot the trajectories
    if numComponents >= 3
        % Plot naive (pale) and experienced (bold) trajectories in 3D
        plot3(mean_traj_naive(:,1), mean_traj_naive(:,2), mean_traj_naive(:,3), ...
              'Color', pale_col, 'LineWidth', 2);
        hold on;
        plot3(mean_traj_experienced(:,1), mean_traj_experienced(:,2), mean_traj_experienced(:,3), ...
              'Color', col, 'LineWidth', 2);
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
        grid on;
        view(3);
    else
        % 2D plot for 2 components
        plot(mean_traj_naive(:,1), mean_traj_naive(:,2), 'Color', pale_col, 'LineWidth', 2);
        hold on;
        plot(mean_traj_experienced(:,1), mean_traj_experienced(:,2), 'Color', col, 'LineWidth', 2);
        xlabel('PC1'); ylabel('PC2');
    end
    title(sprintf('Avg Traj: %s', regions_label(r)));
end
t.Padding = 'compact';
t.TileSpacing = 'compact';

% Figure 3: Projections on Each PCA Dimension (using tiledLayout)
% Each tile shows the mean projection (over bins) of one PCA component for both epochs.
% The reward zone (using idx_reward_zone) is lightly shaded.
figure;
% t = tiledlayout(numRegions, numComponents, 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:numRegions
    % Retrieve the mean trajectories for naive and experienced epochs.
    mean_traj_naive = regionMeanTraj{r}.naive;  % size: [numBins x numComponents]
    mean_traj_experienced = regionMeanTraj{r}.experienced;
    
    % Compute common y-axis limits for this region (across all components)
    all_data = [mean_traj_naive(:); mean_traj_experienced(:)];
    y_min = min(all_data);
    y_max = max(all_data);
    
    % Get the region's color and compute a paler version.
    col = cmap_regions(r, :);
    pale_col = col + (1 - col) * 0.5;  % 50% blend with white for a lighter tone
    
    for comp = 1:components_actually_used(r)
        % ax = nexttile;
        ax = subplot(numRegions, numComponents, (r-1)*numComponents + comp);
        % Plot naive projection (pale)
        plot(ax, mean_traj_naive(:, comp), 'Color', pale_col, 'LineWidth', 2);
        hold(ax, 'on');
        % Plot experienced projection (bold)
        plot(ax, mean_traj_experienced(:, comp), 'Color', col, 'LineWidth', 2);
        xlabel(ax, 'Bin');
        ylabel(ax, sprintf('PC%d', comp));
        title(ax, sprintf('%s - PC%d', regions_label(r), comp));
        
        % Lightly shade the reward zone:
        % Use the common y-axis limits for this row
        x_patch = [idx_reward_zone(1), idx_reward_zone(end), idx_reward_zone(end), idx_reward_zone(1)];
        y_patch = [y_min, y_min, y_max, y_max];
        p = patch(ax, x_patch, y_patch, [0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        uistack(p, 'bottom');  % Ensure the patch is behind the line plots
        
        % Set the common y-axis limits for this tile
        ylim(ax, [y_min, y_max]);
    end
end

%% Coding Dimension Analysis: Difference between Naive and Experienced

% Assumptions:
% spatialFiring_gf: [numUnits x numTrials x numBins] firing rate data
% idx_units: [numRegions x numUnits] logical matrix for unit assignment to regions
% regions_label: 1 x numRegions string array with region labels
% period_naive: vector of trial indices for the naive epoch
% period_experienced: vector of trial indices for the experienced epoch
% idx_reward_zone: vector of bin indices that fall within the reward zone

numRegions = size(idx_units, 1);
regionCodingProjection = cell(numRegions, 1);

% Define a color map for the regions (one color per region)
cmap_regions = lines(numRegions);

for r = 1:numRegions
    % Get indices for the current region's units
    region_units = find(idx_units(r, :));
    if isempty(region_units)
        warning('No units found for region %s', regions_label(r));
        continue;
    end
    
    % Extract firing data for this region (nUnitsRegion x numTrials x numBins)
    regionFiring = spatialFiring_gf_z(region_units, :, :);
    [nUnitsRegion, numTrials, numBins] = size(regionFiring);
    
    % Compute the average activity across spatial bins for each trial.
    % This produces a matrix of size [nUnitsRegion x numTrials] used to compute the coding direction.
    activity_avg = squeeze(mean(regionFiring, 3));
    
    % Compute mean responses for naive and experienced epochs
    mu_naive = mean(activity_avg(:, period_naive), 2, 'omitnan');
    mu_experienced = mean(activity_avg(:, period_experienced), 2, 'omitnan');
    
    % Compute the coding direction (difference vector) and normalize it to be a unit vector.
    CD = mu_experienced - mu_naive;
    CD = CD / norm(CD);
    
    % Project the activity from each spatial bin onto the coding direction.
    % For each bin, extract the activity (neurons x trials) and compute the dot product.
    projections_bins = zeros(numBins, numTrials);
    for b = 1:numBins
        activity_bin = squeeze(regionFiring(:, :, b));  % [nUnitsRegion x numTrials]
        projections_bins(b, :) = CD' * activity_bin;
    end
    
    % Compute the mean projection trajectory (across trials) for naive and experienced epochs.
    mean_proj_naive = mean(projections_bins(:, period_naive), 2, 'omitnan');
    mean_proj_experienced = mean(projections_bins(:, period_experienced), 2, 'omitnan');
    
    % Store the results in a structure for this region.
    regionCodingProjection{r}.naive = mean_proj_naive;
    regionCodingProjection{r}.experienced = mean_proj_experienced;
end

% Plotting: Average Projection Trajectories for Each Area Using tiledLayout

% This figure will show, for each region, the average projection (across spatial bins)
% for both naive and experienced epochs in separate tiles.
figure;
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:numRegions
    nexttile;
    if isempty(regionCodingProjection{r})
        continue;
    end
    
    % Get the region's color and create a paler version for the naive epoch.
    col = cmap_regions(r, :);
    pale_col = col + (1 - col) * 0.5;  % 50% blend with white
    
    % Retrieve the average projections for this region.
    mean_proj_naive = regionCodingProjection{r}.naive;
    mean_proj_experienced = regionCodingProjection{r}.experienced;
    
    % Plot the naive (pale) and experienced (bold) trajectories.
    plot(mean_proj_naive, 'Color', pale_col, 'LineWidth', 2);
    hold on;
    plot(mean_proj_experienced, 'Color', col, 'LineWidth', 2);
    xlabel('Spatial Bin');
    ylabel('Projection (a.u.)');
    title(sprintf('Coding Dim: %s', regions_label(r)));

    axis tight
    
    % Optionally, lightly shade the reward zone.
    yl = ylim;
    x_patch = [idx_reward_zone(1), idx_reward_zone(end), idx_reward_zone(end), idx_reward_zone(1)];
    y_patch = [yl(1), yl(1), yl(2), yl(2)];
    patch(x_patch, y_patch, [0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    uistack(gca, 'bottom');  % Ensure the patch stays behind the line plots
end
linkaxes


%% Coding Dimension Analysis: Difference between Naive and Experienced
% Assumptions:
% spatialFiring_gf_z: [numUnits x numTrials x numBins] z-scored firing rate data
% idx_units: [numRegions x numUnits] logical matrix for unit assignment to regions
% regions_label: 1 x numRegions string array with region labels
% period_naive: vector of trial indices for the naive epoch
% period_experienced: vector of trial indices for the experienced epoch
% idx_reward_zone: vector of bin indices that fall within the reward zone

numRegions = size(idx_units, 1);
regionCodingProjection = cell(numRegions, 1);

% Define a color map for the regions (one color per region)
cmap_regions = lines(numRegions);

for r = 1:numRegions
    % Get indices for the current region's units
    region_units = find(idx_units(r, :));
    if isempty(region_units)
        warning('No units found for region %s', regions_label(r));
        continue;
    end
    
    % Extract firing data for this region (nUnitsRegion x numTrials x numBins)
    regionFiring = spatialFiring_gf_z(region_units, :, :);
    [nUnitsRegion, numTrials, numBins] = size(regionFiring);
    
    % Compute the average activity across spatial bins for each trial.
    % This yields a matrix of size [nUnitsRegion x numTrials]
    activity_avg = squeeze(mean(regionFiring, 3));
    
    % Compute mean responses for naive and experienced epochs
    mu_naive = mean(activity_avg(:, period_naive), 2, 'omitnan');
    mu_experienced = mean(activity_avg(:, period_experienced), 2, 'omitnan');
    
    % Compute the coding direction (difference vector) and normalize it.
    CD = mu_experienced - mu_naive;
    CD = CD / norm(CD);
    
    % Project the activity from each spatial bin onto the coding direction.
    % For each bin, extract activity (nUnitsRegion x numTrials) and compute the dot product.
    projections_bins = zeros(numBins, numTrials);
    for b = 1:numBins
        activity_bin = squeeze(regionFiring(:, :, b));  % [nUnitsRegion x numTrials]
        projections_bins(b, :) = CD' * activity_bin;
    end
    
    % Now, compute the average projection per trial (averaging across spatial bins)
    trial_projection = mean(projections_bins, 1);  % 1 x numTrials
    
    % (Optionally, you could still compute spatial profiles for each epoch)
    mean_proj_naive = mean(projections_bins(:, period_naive), 2, 'omitnan');
    mean_proj_experienced = mean(projections_bins(:, period_experienced), 2, 'omitnan');
    
    % Store results for this region
    regionCodingProjection{r}.trial = trial_projection;
    regionCodingProjection{r}.naive = mean_proj_naive;
    regionCodingProjection{r}.experienced = mean_proj_experienced;
end

% Compute Global Y-Axis Limits Across Regions (for trial projections only)
global_trial_proj = [];
for r = 1:numRegions
    if ~isempty(regionCodingProjection{r})
        global_trial_proj = [global_trial_proj; regionCodingProjection{r}.trial(:)];
    end
end
global_ymin = min(global_trial_proj);
global_ymax = max(global_trial_proj);

% Plotting: Progression of Average Projection Across Trials for Each Area
% This figure shows, for each region, a line plot of the average projection (averaged over all spatial bins)
% as a function of trial number. The naive and experienced trial periods are shaded in green and magenta.
figure;
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:numRegions
    ax = nexttile;
    if isempty(regionCodingProjection{r})
        continue;
    end
    
    % Retrieve the trial projection for this region (ensure it's a column vector)
    trial_proj = regionCodingProjection{r}.trial(:);  % [numTrials x 1]
    
    % Plot the progression of the average projection across trials.
    plot(ax, 1:numTrials, trial_proj, 'k', 'LineWidth', 2);
    hold(ax, 'on');
    xlabel(ax, 'Trial');
    ylabel(ax, 'Avg Projection (a.u.)');
    title(ax, sprintf('Coding Dim Progression: %s', regions_label(r)));
    axis tight
    % Set common y-axis limits for all region tiles
    ylim(ax, [global_ymin, global_ymax]);
    
    % Shade the naive trial period with the specified green color.
    if ~isempty(period_naive)
        x_naive = [min(period_naive) max(period_naive) max(period_naive) min(period_naive)];
        y_patch = [global_ymin, global_ymin, global_ymax, global_ymax];
        patch(ax, x_naive, y_patch, col_naive, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    
    % Shade the experienced trial period with the specified magenta color.
    if ~isempty(period_experienced)
        x_exp = [min(period_experienced) max(period_experienced) max(period_experienced) min(period_experienced)];
        y_patch = [global_ymin, global_ymin, global_ymax, global_ymax];
        patch(ax, x_exp, y_patch, col_experienced, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    
    % Replot the trial projection on top so it is visible.
    plot(ax, 1:numTrials, trial_proj, 'k', 'LineWidth', 2);
end

% Additional Tiles (e.g., lick ratio, licks, occupancy)
% These tiles will not have the global y-axis limits enforced.
nexttile
hold on
mov_window = 10;
shadedErrorBar(1:numTrials, movmean(lick_ratio, mov_window), movmean(lick_ratio, mov_window)/sqrt(mov_window))
title('lick ratio')
axis tight
if ~isempty(period_naive)
    x_naive = [min(period_naive) max(period_naive) max(period_naive) min(period_naive)];
    yl = ylim;
    patch(x_naive, [yl(1) yl(1) yl(2) yl(2)], col_naive, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
if ~isempty(period_experienced)
    x_exp = [min(period_experienced) max(period_experienced) max(period_experienced) min(period_experienced)];
    yl = ylim;
    patch(x_exp, [yl(1) yl(1) yl(2) yl(2)], col_experienced, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

nexttile
imagesc(spatial_licking_gf)
title('licks')

nexttile
imagesc(spatialOccupation_time_gf)
title('occupancy')

% Figure: Correlation of Each Area's CD Projection and movmean(lick_ratio)

mov_window = 5;
lick_ratio_smoothed = movmean(lick_ratio, mov_window);  % smooth the lick ratio

figure;
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:numRegions
    ax = nexttile;
    if isempty(regionCodingProjection{r})
        continue;
    end
    
    % Retrieve the CD projection for this region (as a column vector)
    trial_proj = regionCodingProjection{r}.trial(:);  % [numTrials x 1]
    
    % Scatter plot: x = CD projection, y = smoothed lick ratio
    scatter(ax, trial_proj, lick_ratio_smoothed, 36, 'filled');
    xlabel(ax, 'CD Projection');
    ylabel(ax, 'Lick Ratio (movmean)');
    title(ax, sprintf('Area: %s', regions_label(r)));
    hold(ax, 'on');
    
    % Add a least-squares line
    lsline
    
    % Compute and display Pearson correlation coefficient and p-value
    [R, P] = corr(trial_proj, lick_ratio_smoothed', 'rows', 'complete');
    x_text = min(trial_proj) + 0.05*(max(trial_proj)-min(trial_proj));
    y_text = max(lick_ratio_smoothed) - 0.1*(max(lick_ratio_smoothed)-min(lick_ratio_smoothed));
    text(ax, x_text, y_text, sprintf('r = %.2f, p = %.3f', R, P), 'FontSize', 8, 'BackgroundColor', 'w');
end

%% Cross-Correlation Analysis: CD Projection vs. Lick Ratio
% Assumptions:
% - regionCodingProjection{r}.trial contains the CD projection per trial (1 x numTrials)
% - lick_ratio_smoothed is defined (length = numTrials)
% - regions_label is a 1 x numRegions string array

% Set maximum lag (in trials) for cross-correlation
maxlag = 20;

figure;
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:numRegions
    ax = nexttile;
    if isempty(regionCodingProjection{r})
        continue;
    end
    
    % Get the CD projection for this region (as a column vector)
    trial_proj = regionCodingProjection{r}.trial(:);  % [numTrials x 1]
    
    % Compute the cross-correlation between the CD projection and the lick_ratio
    % 'coeff' normalizes the correlation so that it is scaled as a correlation coefficient.
    [xc, lags] = xcorr(trial_proj, lick_ratio_smoothed, maxlag, 'coeff');
    
    % Plot the cross-correlation as a stem plot
    stem(ax, lags, xc, 'filled', 'Color', cmap_regions(r, :));
    xlabel(ax, 'Lag (trials)');
    ylabel(ax, 'Cross-correlation');
    title(ax, sprintf('Cross-Correlation: %s', regions_label(r)));
    
end

linkaxes


%% TCA Analysis per Region

if ~exist([animal_id, '_TCA_results.mat'], "file")
    % TCA parameters (using z-scored data)
    tca_type = 'cp_nmu';
    normalisation_to_apply = 'min-max';
    cross_validation_fold = 5;
    maxNumFactors = 8;
    max_iterations = 50;
    
    % Preallocate a cell array to store TCA results per region
    regionTCAResults = cell(numRegions, 1);
    
    for r = 1:numRegions
        % Get indices for the current region's units
        region_units = find(idx_units(r, :));
        if isempty(region_units) || numel(region_units) < 3
            warning('No units found for region %s', regions_label(r));
            continue;
        end
        
        % Extract the z-scored firing data for the current region: [nUnits x numTrials x numBins]
        regionFiring = spatialFiring_gf(region_units, :, :);
        
        % Permute dimensions to [nUnits x numBins x numTrials] for TCA
        regionData = permute(regionFiring, [1, 3, 2]);
        
        % Run TCA with cross-validation
        [best_mdl, variance_explained, mean_cv_errors, sem_cv_errors] = ...
            tca_with_cv(regionData, tca_type, normalisation_to_apply, cross_validation_fold, maxNumFactors, max_iterations);
        
        best_numFactors = length(variance_explained);
        
        % Store the results in a structure for this region
        regionTCAResults{r}.best_mdl = best_mdl;
        regionTCAResults{r}.variance_explained = variance_explained;
        regionTCAResults{r}.mean_cv_errors = mean_cv_errors;
        regionTCAResults{r}.sem_cv_errors = sem_cv_errors;
        regionTCAResults{r}.best_numFactors = best_numFactors;

    end
    save([animal_id, '_TCA_results.mat'], "regionTCAResults")
else
    load([animal_id, '_TCA_results.mat'])
end




%% Plot TCA

for r = 1:numRegions
    best_mdl = regionTCAResults{r}.best_mdl;
    variance_explained = regionTCAResults{r}.variance_explained;
    mean_cv_errors = regionTCAResults{r}.mean_cv_errors;
    sem_cv_errors = regionTCAResults{r}.sem_cv_errors;
    best_numFactors = regionTCAResults{r}.best_numFactors;
    %---------------------------
    % Plot 1: CV Error vs. Number of Factors
    %---------------------------
    figure;
    shadedErrorBar(1:maxNumFactors, mean_cv_errors, sem_cv_errors);
    xlabel('Number of Factors');
    ylabel('Mean CV Reconstruction Error');
    title(sprintf('CV Error vs. Number of Factors: %s', regions_label(r)));
    
    %---------------------------
    % Plot 2: Spatial Factors (Factors in U{2})
    %---------------------------
    figure;
    t_spatial = tiledlayout(best_numFactors, 1, "TileSpacing", "compact");
    for iFactor = 1:best_numFactors
        nexttile
        plot(best_mdl.U{2}(:, iFactor), 'LineWidth', 1.5);
        axis tight
        if iFactor ~= best_numFactors
            xticks([]);  % Remove x-tick labels for all but the last subplot
        end
        ylabel(sprintf('Comp %d', iFactor));
        title(sprintf('Variance Explained: %.2f%%', variance_explained(iFactor)*100));
    end
    xlabel(t_spatial, 'Spatial Bin');
    sgtitle(sprintf('Spatial Factors: %s', regions_label(r)));
    
    %---------------------------
    % Plot 3: Trial Factors (Factors in U{3})
    %---------------------------
    figure;
    t_trial = tiledlayout(best_numFactors, 1, "TileSpacing", "compact");
    for iFactor = 1:best_numFactors
        nexttile
        shadedErrorBar(1:numTrials, ...
            movmean(best_mdl.U{3}(:, iFactor), 10), ...
            movstd(best_mdl.U{3}(:, iFactor), 10)/sqrt(10));
        axis tight
        if iFactor ~= best_numFactors
            xticks([]);
        end
        ylabel(sprintf('Comp %d', iFactor));
        title(sprintf('Variance Explained: %.2f%%', variance_explained(iFactor)*100));
    end
    xlabel(t_trial, 'Trial #');
    sgtitle(sprintf('Trial Factors: %s', regions_label(r)));
    
end

%% Combined TCA: Spatial Factors and Naive vs Experienced (Trial Factors)
numRegions = size(regions_label, 2);
for r = 1:numRegions
    % Skip if no TCA result for this region
    if isempty(regionTCAResults{r})
        continue;
    end
    
    % Retrieve TCA results for this region
    best_mdl         = regionTCAResults{r}.best_mdl;
    variance_explained = regionTCAResults{r}.variance_explained;
    best_numFactors  = regionTCAResults{r}.best_numFactors;
    
    % Scale figure height based on the number of factors
    baseHeight = 150;  % adjust this value as needed (in pixels per factor)
    figHeight = best_numFactors * baseHeight;
    
    % Create a new figure with a tiled layout:
    % Rows: each factor; Columns: 2 (left = spatial factor, right = naive vs. experienced)
    figure('Position', [100, 100, 640, figHeight]);
    t = tiledlayout(best_numFactors, 2, "TileSpacing", "compact", "Padding", "compact");
    
    for iFactor = 1:best_numFactors
        %-------------------------
        % Left Column: Spatial Factor Plot (U{2})
        %-------------------------
        ax1 = nexttile(t, 2*(iFactor-1)+1);
        plot(ax1, best_mdl.U{2}(:, iFactor), 'LineWidth', 1.5);
        axis(ax1, 'tight');
        ylabel(ax1, sprintf('Comp %d', iFactor));
        title(ax1, 'Spatial');
        
        %-------------------------
        % Right Column: Naive vs. Experienced Error Bar Plot (from U{3})
        %-------------------------
        ax2 = nexttile(t, 2*(iFactor-1)+2);
        % Plot the error bar plot using your custom function
        my_errorbar_plot( best_mdl.U{3}(period_naive, iFactor), best_mdl.U{3}(period_experienced, iFactor) );
        xticklabels(ax2, {'naive', 'experienced'});
        ylabel(ax2, 'Trial Factor');
        title(ax2, sprintf('Naive vs Exp: Comp %d', iFactor));
        
        % Perform a t-test between naive and experienced trial factors for this component
        [~, p_val] = ttest2(best_mdl.U{3}(period_naive, iFactor), best_mdl.U{3}(period_experienced, iFactor));
        
        % Add significance marker using sigstar (make sure sigstar is in your path)
        hold(ax2, 'on');
        if p_val < 0.05
            sigstar({[1, 2]}, p_val);
        end
    end
    title(t, regions_label{r});  % if regions_label is a cell array; use regions_label(r) if it's a string array
end

%% Plot raw and processed licks


% Exclude data before the trial index becomes 1
valid_idx = trial_binned >= 1;

pos_valid = pos_binned_gf(valid_idx);

max_allowed_position = floor(max(pos_valid)/100)*100;

pos_valid(pos_valid < 0) = 0;
pos_valid(pos_valid > max_allowed_position) = max_allowed_position;
licks_valid = licks_binned(valid_idx);
trial_valid = trial_binned(valid_idx);


% Set maximum allowed position and bin size
spatial_bin_size = 2.5;

% Define bin edges from 0 to 500 (inclusive) in steps of 2.5.
% This creates 201 edges and hence 200 bins.
edges = 0:spatial_bin_size:max_allowed_position;

% Use discretize to assign each position to a bin.
% With the default behavior, bins are [edge(i), edge(i+1)) except for the last bin which is inclusive.
spatial_bin_index = discretize(pos_valid, edges);

% Determine the number of bins and trials.
num_bins = max(spatial_bin_index); % exactly, since we set edges from 0 to 500 in steps of 2.5.
num_trials = max(trial_valid);

% Use accumarray to sum the lick counts for each spatial bin and trial.
licks_matrix = accumarray([spatial_bin_index(:), trial_valid(:)], licks_valid(:), [num_bins, num_trials], @sum, 0);

figure
subplot(1, 2, 1)
imagesc(licks_matrix')

subplot(1, 2, 2)
imagesc(spatial_licking_gf)


%% Plot lick ratio and lick error
% Assume the following variables are defined (all are 1 x T vectors):
%   pos_binned_gf : continuous position (cm) at each timepoint
%   licks_binned  : integer number of licks at each timepoint
%   trial_binned  : integer trial identifier at each timepoint

% Determine the total number of trials from trial_binned.
num_trials = max(trial_binned);

% Initialize a cell array to hold lick positions for each trial.
lick_positions = cell(num_trials, 1);

% Loop over each trial.
for trial = 1:num_trials
    % Find indices for the current trial.
    idx = trial_binned == trial;
    
    % Get the positions and lick counts for this trial.
    trial_positions = pos_binned_gf(idx);
    trial_lick_counts = licks_binned(idx);
    
    % Initialize an empty array for lick positions in this trial.
    trial_lick_positions = [];
    
    % Loop through each timepoint in the trial.
    for i = 1:length(trial_positions)
        % If there are any licks at this timepoint...
        if trial_lick_counts(i) > 0
            % Replicate the position according to the number of licks and append.
            trial_lick_positions = [trial_lick_positions, repmat(trial_positions(i), 1, trial_lick_counts(i))];
        end
    end
    
    % Store the lick positions for the current trial in the cell array.
    lick_positions{trial} = trial_lick_positions;
end


nTrials = numel(lick_positions);

% Preallocate output vectors
trial_lick_errors = nan(nTrials, 1);
shuffled_lick_error_means = nan(nTrials, 1);
shuffled_lick_error_stds = nan(nTrials, 1);
chance_quantile_all = nan(nTrials, 1);

rz_start = idx_reward_zone(1) + 10 + 1;  % reward zone start position in cm

% Loop over each trial
for itrial = 1:nTrials
    % Get lick positions for the current trial.
    % Note: ensure that lick_positions{itrial} is a row or column vector.
    l_trial = lick_positions{itrial};
    
    if ~isempty(l_trial)
    % Calculate lick precision metrics for this trial.
    [trial_lick_errors(itrial), shuffled_lick_error_means(itrial), shuffled_lick_error_stds(itrial), chance_quantile_all(itrial)] = ...
        calculate_lick_precision(l_trial, rz_start);
    end
end

% Compute z-scored lick errors per trial.
zscored_lick_errors = (trial_lick_errors - shuffled_lick_error_means) ./ shuffled_lick_error_stds;

mov_window = 20;

% Plot the z-scored lick errors and the lick ratio.
figure;
subplot(1, 2, 1)
shadedErrorBar(1:nTrials, movmean(zscored_lick_errors, mov_window, 'omitmissing'), movstd(zscored_lick_errors, mov_window, 'omitmissing')/sqrt(mov_window));
xlabel('Trial');
ylabel('Z-scored Lick Error');
title('Z-scored Lick Errors');
axis tight
if ~isempty(period_naive)
    x_naive = [min(period_naive) max(period_naive) max(period_naive) min(period_naive)];
    yl = ylim;
    patch(x_naive, [yl(1) yl(1) yl(2) yl(2)], col_naive, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
if ~isempty(period_experienced)
    x_exp = [min(period_experienced) max(period_experienced) max(period_experienced) min(period_experienced)];
    yl = ylim;
    patch(x_exp, [yl(1) yl(1) yl(2) yl(2)], col_experienced, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

subplot(1, 2, 2)
shadedErrorBar(1:nTrials, movmean(lick_ratio, mov_window, 'omitmissing'), movstd(lick_ratio, mov_window, 'omitmissing')/sqrt(mov_window));
xlabel('Trial');
ylabel('Lick Ratio');
title('Lick Ratio');
axis tight
if ~isempty(period_naive)
    x_naive = [min(period_naive) max(period_naive) max(period_naive) min(period_naive)];
    yl = ylim;
    patch(x_naive, [yl(1) yl(1) yl(2) yl(2)], col_naive, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
if ~isempty(period_experienced)
    x_exp = [min(period_experienced) max(period_experienced) max(period_experienced) min(period_experienced)];
    yl = ylim;
    patch(x_exp, [yl(1) yl(1) yl(2) yl(2)], col_experienced, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end