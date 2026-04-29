%% Load data
load("data_for_Theo.mat")
load('period_experienced.mat')
load('licking_behaviour.mat')
%% Compare behavioural metrics
num_shuffles_licks = 100;
all_zcored_lick_error = cell(1, length(spatial_licking_gf));
zscore_learning_points = nan(1, length(spatial_licking_gf));
learning_points = period_experienced(:, 1);
for ianimal = 1:length(spatial_licking_gf)
    current_licks = spatial_licking_gf{ianimal};
    current_lick_ratio = lick_ratio{ianimal};
    [num_trials, num_bins] = size(current_licks);
    current_lick_error = mean(current_licks .* (repmat(1:num_bins, num_trials, 1) - 100).^2, 2);
    lick_error_distribution = nan(size(current_licks, 1), num_shuffles_licks);
    for ishuffle = 1:num_shuffles_licks
        lick_error_distribution(:, ishuffle) = mean(current_licks(:, randperm(num_bins, num_bins)) .* (repmat(1:num_bins, num_trials, 1) - 100).^2, 2);
    end
    z_scored_lick_error = (current_lick_error - mean(lick_error_distribution, 2))./std(lick_error_distribution, [], 2);
    zscore_learning = find(movsum(z_scored_lick_error <= -2, [0 9], 1, "omitmissing") >= 8, 10) + 1;
    
    if ~isempty(zscore_learning)
        zscore_learning(zscore_learning < 11) = nan;
        zscore_learning_point = zscore_learning(find(~isnan(zscore_learning), 1));
    else
        zscore_learning_point = nan;
    end
    figure
    subplot(1, 2, 1)
    imagesc(current_licks)
    subplot(1, 2, 2)
    plot(movmean(current_lick_ratio, 10))
    xline(learning_points(ianimal), 'b');
    hold on
    yyaxis("right")
    plot(movmean(z_scored_lick_error, 10, 1, "omitmissing"))
    xline(zscore_learning_point, 'r')
    sgtitle(sprintf('animal %d', ianimal))
    all_zcored_lick_error{ianimal} = z_scored_lick_error;
    zscore_learning_points(ianimal) = zscore_learning_point;
end
%% Parameters
n_animals = length(analysis_spatial);
area_activity = cell(1, n_animals);
n_components_reduced = 4;
n_trials_window = -3:3;
n_bins_window = -3:3;
n_shuffles = 10;
max_shift_bins = 3;
landmarks = [0, 50;...
             100, 125;...
         	 175, 200;...
             250, 275;...
         	 325, 350;...
         	 400, 425]/2.5;
% --- Define all area pairs to analyze ---
area_pairs_to_analyze = {'CA1', 'V1'; ...
                         'CA1', 'DG'; ...
                         'CA1', 'CA3'; ...
                         'CA1', 'RSC'; ...
                         'CA1', 'SUB'; ...
                         'V1', 'RSC'; ...
                         'RSC', 'SUB'; ...
                         'CA3', 'DG'};
n_pairs = size(area_pairs_to_analyze, 1);
%% Create Figure Handles for each pair
% We create one set of figures (t, u, v) for each pair
t_handles = cell(n_pairs, 1); % Trial correlation
s_handles = cell(n_pairs, 1); % Trial correlation vs behaviour
u_handles = cell(n_pairs, 1); % Precession
v_handles = cell(n_pairs, 1); % Spatial bin correlation
w_handles = cell(n_pairs, 1); % Precession in space
for ipair = 1:n_pairs
    pair_name = sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});
    
    figure('Name', ['Trial Correlation: ' pair_name]);
    t_handles{ipair} = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    figure('Name', ['CCA vs Behaviour: ' pair_name]);
    s_handles{ipair} = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    
    figure('Name', ['Precession: ' pair_name]);
    u_handles{ipair} = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    
    figure('Name', ['Spatial Correlation: ' pair_name]);
    v_handles{ipair} = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    figure('Name', ['Spatial Precession: ' pair_name]);
    w_handles{ipair} = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
end
%% Create struct to hold all animal results
group_results = cell(n_pairs, 1);
for ipair = 1:n_pairs
    group_results{ipair}.pair_name = sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});
    
    % --- CC1 Data ---
    group_results{ipair}.all_bins_corr = [];
    group_results{ipair}.all_bins_precession = [];
    group_results{ipair}.trial_corr_early = [];
    group_results{ipair}.trial_corr_pre = [];
    group_results{ipair}.trial_corr_post = [];
    group_results{ipair}.trial_precession_early = [];
    group_results{ipair}.trial_precession_pre = [];
    group_results{ipair}.trial_precession_post = [];
    
    % --- ADDED: CC2 Data ---
    group_results{ipair}.all_bins_corr_cc2 = [];
    group_results{ipair}.all_bins_precession_cc2 = [];
    group_results{ipair}.trial_corr_early_cc2 = [];
    group_results{ipair}.trial_corr_pre_cc2 = [];
    group_results{ipair}.trial_corr_post_cc2 = [];
    group_results{ipair}.trial_precession_early_cc2 = [];
    group_results{ipair}.trial_precession_pre_cc2 = [];
    group_results{ipair}.trial_precession_post_cc2 = [];
    % -----------------------
end
%% Run CCA for all animals and all defined pairs
for ianimal = 1:n_animals
    animal_data = permute(analysis_spatial{ianimal}.firing.cued.freq, [1, 3, 2]);
    [~, num_bins, num_trials] = size(animal_data);
    animal_areas = units{ianimal}.regions_label;
    n_areas = numel(animal_areas);
    unit_areas = units{ianimal}.idx;
    
    % First, run PCA for all areas for this animal
    for iarea = 1:n_areas
        current_area_activity = animal_data(unit_areas(iarea, :), :, :);
        n_units_area = size(current_area_activity, 1);
        area_activity{ianimal}.(animal_areas(iarea)).spatial_binned_fr = current_area_activity;
        if n_units_area >= n_components_reduced
            [~,current_area_scores,~, ~, explained, ~] = pca(current_area_activity(:, :)', NumComponents=n_components_reduced);
            area_activity_reduced = reshape(current_area_scores', n_components_reduced, num_bins, num_trials);
            area_activity{ianimal}.(animal_areas(iarea)).area_activity_reduced = area_activity_reduced;
            fprintf('%s pca for animal %d with %d dimensions explains %.2f total variance\n', animal_areas(iarea), ianimal, n_components_reduced, sum(explained(1:n_components_reduced)))
        end
    end
    
    % --- Generalized CCA Loop ---
    % Now, loop through each requested pair and run the analysis
    for ipair = 1:n_pairs
        area1_name = area_pairs_to_analyze{ipair, 1};
        area2_name = area_pairs_to_analyze{ipair, 2};
        
        % Get the correct figure handle for this pair
        current_t = t_handles{ipair};
        current_s = s_handles{ipair};
        current_u = u_handles{ipair};
        current_v = v_handles{ipair};
        current_w = w_handles{ipair};
        
        % Check if this animal has reduced data for both areas in the pair
        if isfield(area_activity{ianimal}, area1_name) && isfield(area_activity{ianimal}, area2_name) && ...
           isfield(area_activity{ianimal}.(area1_name), 'area_activity_reduced') && isfield(area_activity{ianimal}.(area2_name), 'area_activity_reduced') 
            
            % Pre-allocate results matrices (now inside the pair loop)
            cca_all_trials = nan(n_components_reduced, num_trials);
            cca_shuffled_trials = nan(n_components_reduced, num_trials, n_shuffles);
            cca_shifted_trials = nan(n_components_reduced, num_trials, 2*max_shift_bins+1);
            cca_shited_shuffled_trials = nan(n_components_reduced, num_trials, 2*max_shift_bins+1, n_shuffles); % (Typo 'shited' kept from original)
            cca_all_bins = nan(n_components_reduced, num_bins);
            cca_shuffled_bins = nan(n_components_reduced, num_bins, n_shuffles);
    
            % --- 1. CCA by Trial ---
            for itrial = 1:num_trials
                trial_window = itrial + n_trials_window;
                trial_window = trial_window(trial_window >= 1 & trial_window <= num_trials);
                
                % Robust NaN check: find bins that are valid in *both* areas
                nan_bins_1 = squeeze(any(isnan(area_activity{ianimal}.(area1_name).area_activity_reduced(:, :, trial_window)), [1, 3]));
                nan_bins_2 = squeeze(any(isnan(area_activity{ianimal}.(area2_name).area_activity_reduced(:, :, trial_window)), [1, 3]));
                valid_bins = ~nan_bins_1 & ~nan_bins_2;
                
                area1_reduced_pooled = area_activity{ianimal}.(area1_name).area_activity_reduced(:, valid_bins, trial_window);
                area2_reduced_pooled = area_activity{ianimal}.(area2_name).area_activity_reduced(:, valid_bins, trial_window);
                
                try
                    [~,~,r,~,~] = canoncorr(area1_reduced_pooled(:, :)', area2_reduced_pooled(:, :)');
                    n_r = min(length(r), n_components_reduced);
                    cca_all_trials(1:n_r, itrial) = r(1:n_r);
                    
                    for ishuffle = 1:n_shuffles
                        % Shuffle bins within the trial window
                        area1_reduced_shuffled = area1_reduced_pooled(:, randperm(sum(valid_bins), sum(valid_bins)), :);
                        [~,~,r,~,~] = canoncorr(area1_reduced_shuffled(:, :)', area2_reduced_pooled(:, :)');
                        n_r = min(length(r), n_components_reduced);
                        cca_shuffled_trials(1:n_r, itrial, ishuffle) = r(1:n_r);
                    end
                    
                    % --- 2. Shifted CCA (for Precession) ---
                    for ishift = -max_shift_bins:max_shift_bins
                        % area1_reduced_shifted = circshift(area1_reduced_pooled, ishift, 2); % Shift bins
                        % remove tails from shifted data
                        if ishift < 0
                            area1_reduced_shifted_cut = area1_reduced_pooled(:, 1-ishift:end, :);
                            area2_reduced_cut = area2_reduced_pooled(:, 1:end+ishift, :);
                        elseif ishift > 0
                            area1_reduced_shifted_cut = area1_reduced_pooled(:, 1:end-ishift, :);
                            area2_reduced_cut = area2_reduced_pooled(:, 1+ishift:end, :);
                        else
                            area1_reduced_shifted_cut = area1_reduced_pooled;
                            area2_reduced_cut = area2_reduced_pooled;
                        end
                        [~,~,r,~,~] = canoncorr(area1_reduced_shifted_cut(:, :)', area2_reduced_cut(:, :)');
                        n_r = min(length(r), n_components_reduced);
                        cca_shifted_trials(1:n_r, itrial, ishift+1+max_shift_bins) = r(1:n_r);
                        
                        for ishuffle = 1:n_shuffles
                            area1_shuffled = area1_reduced_pooled(:, randperm(sum(valid_bins), sum(valid_bins)), :);
                            if ishift < 0
                                area1_shifted_shuffled_cut = area1_shuffled(:, 1-ishift:end, :);
                            elseif ishift > 0
                                area1_shifted_shuffled_cut = area1_shuffled(:, 1:end-ishift, :);
                            else
                                area1_shifted_shuffled_cut = area1_shuffled;
                            end
                            [~,~,r,~,~] = canoncorr(area1_shifted_shuffled_cut(:, :)', area2_reduced_cut(:, :)');
                            n_r = min(length(r), n_components_reduced);
                            cca_shited_shuffled_trials(1:n_r, itrial, ishift+1+max_shift_bins, ishuffle) = r(1:n_r);
                        end
                    end
                catch
                    warning('%s-%s CCA failed on trial %d for animal %d', area1_name, area2_name, itrial, ianimal);
                end
            end
            
            % --- Plot Trial-wise CCA (MODIFIED) ---
            nexttile(current_t)
            plot(cca_all_trials(1, :)', 'LineWidth', 1, 'Color', 'b') % CC1
            hold on
            shadedErrorBar(1:num_trials, squeeze(mean(cca_shuffled_trials(1, :, :), 3, "omitmissing")), squeeze(sem(cca_shuffled_trials(1, :, :), 3)), 'lineProps', {'Color', 'b', 'LineStyle', '--'})
            plot(cca_all_trials(2, :)', 'LineWidth', 1, 'Color', 'r') % CC2
            shadedErrorBar(1:num_trials, squeeze(mean(cca_shuffled_trials(2, :, :), 3, "omitmissing")), squeeze(sem(cca_shuffled_trials(2, :, :), 3)), 'lineProps', {'Color', 'r', 'LineStyle', '--'})
            legend({'CC1', 'CC1 Shuffle', 'CC2', 'CC2 Shuffle'}, 'Location', 'best')
            axis tight
            yyaxis("right")
            plot(movmean(-all_zcored_lick_error{ianimal}, 10, 'omitmissing'), 'LineWidth', 1, 'Color', 'k');
            title(sprintf('%s-%s - animal %d', area1_name, area2_name, ianimal))
            xline(zscore_learning_points(ianimal))
            
            % --- Plot CCA vs behaviour scatter (MODIFIED) ---
            nexttile(current_s)
            scatter(cca_all_trials(1, :)', movmean(-all_zcored_lick_error{ianimal}, 10, 'omitmissing')', 'b') % CC1
            hold on
            scatter(cca_all_trials(2, :)', movmean(-all_zcored_lick_error{ianimal}, 10, 'omitmissing')', 'r') % CC2
            lsline
            axis tight
            [rho1, pval1] = corr(cca_all_trials(1, :)', movmean(-all_zcored_lick_error{ianimal}, 10, 'omitmissing'), 'Rows','complete');
            [rho2, pval2] = corr(cca_all_trials(2, :)', movmean(-all_zcored_lick_error{ianimal}, 10, 'omitmissing'), 'Rows','complete');
            title(sprintf('A%d (CC1: r=%.2f, p=%.2f | CC2: r=%.2f, p=%.2f)', ianimal, rho1, pval1, rho2, pval2))
            legend({'CC1', 'CC2'})

            
            % --- Calculate and Plot Precession (MODIFIED) ---
            % CC1
            area2_leads_area1 = squeeze(mean(cca_shifted_trials(1, :, 1:max_shift_bins), 3)); % Negative shifts
            area1_leads_area2 = squeeze(mean(cca_shifted_trials(1, :, max_shift_bins+2:end), 3)); % Positive shifts
            precession_strength = (area1_leads_area2 - area2_leads_area1)./(area1_leads_area2 + area2_leads_area1);
            
            area2_leads_area1_shuffled = squeeze(mean(cca_shited_shuffled_trials(1, :, 1:max_shift_bins, :), 3));
            area1_leads_area2_shuffled = squeeze(mean(cca_shited_shuffled_trials(1, :, max_shift_bins+2:end, :), 3));
            precession_strength_shuffled = (area1_leads_area2_shuffled - area2_leads_area1_shuffled)./(area1_leads_area2_shuffled + area2_leads_area1_shuffled);
            
            % --- ADDED: CC2 Calculation ---
            area2_leads_area1_cc2 = squeeze(mean(cca_shifted_trials(2, :, 1:max_shift_bins), 3)); % Negative shifts
            area1_leads_area2_cc2 = squeeze(mean(cca_shifted_trials(2, :, max_shift_bins+2:end), 3)); % Positive shifts
            precession_strength_cc2 = (area1_leads_area2_cc2 - area2_leads_area1_cc2)./(area1_leads_area2_cc2 + area2_leads_area1_cc2);
            
            area2_leads_area1_shuffled_cc2 = squeeze(mean(cca_shited_shuffled_trials(2, :, 1:max_shift_bins, :), 3));
            area1_leads_area2_shuffled_cc2 = squeeze(mean(cca_shited_shuffled_trials(2, :, max_shift_bins+2:end, :), 3));
            precession_strength_shuffled_cc2 = (area1_leads_area2_shuffled_cc2 - area2_leads_area1_shuffled_cc2)./(area1_leads_area2_shuffled_cc2 + area2_leads_area1_shuffled_cc2);
            % -------------------------------

            nexttile(current_u)
            plot(precession_strength, 'LineWidth', 1, 'Color', 'b') % CC1
            hold on
            shadedErrorBar(1:num_trials, squeeze(mean(precession_strength_shuffled, 2, "omitmissing")), squeeze(sem(precession_strength_shuffled, 2)), 'lineProps', {'Color', 'b', 'LineStyle', '--'})
            plot(precession_strength_cc2, 'LineWidth', 1, 'Color', 'r') % CC2
            shadedErrorBar(1:num_trials, squeeze(mean(precession_strength_shuffled_cc2, 2, "omitmissing")), squeeze(sem(precession_strength_shuffled_cc2, 2)), 'lineProps', {'Color', 'r', 'LineStyle', '--'})
            legend({'CC1', 'CC1 Shuffle', 'CC2', 'CC2 Shuffle'}, 'Location', 'best')
            
            title(sprintf('%s-%s - animal %d', area1_name, area2_name, ianimal))
            yline(0, '--')
            xline(zscore_learning_points(ianimal))
            axis tight
            ylim([-0.4, 0.4])
            yyaxis('right')
            plot(movmean(-all_zcored_lick_error{ianimal}, 10, 'omitmissing'), 'LineWidth', 1, 'Color', 'k');
            
            % --- 3. CCA by Spatial Bin ---
            
            % --- Rename old matrix for clarity ---
            cca_shuffled_bins_by_trial = nan(n_components_reduced, num_bins, n_shuffles);
            % --- Add new matrices for spatial precession ---
            cca_shifted_bins = nan(n_components_reduced, num_bins, 2*max_shift_bins+1);
            cca_shited_shuffled_bins_by_bin = nan(n_components_reduced, num_bins, 2*max_shift_bins+1, n_shuffles); % (Typo 'shited' kept)
            % Robust NaN check: find trials that are valid in *both* areas
            nan_trials_1 = squeeze(any(isnan(area_activity{ianimal}.(area1_name).area_activity_reduced), [1, 2]));
            nan_trials_2 = squeeze(any(isnan(area_activity{ianimal}.(area2_name).area_activity_reduced), [1, 2]));
            valid_trials = ~nan_trials_1 & ~nan_trials_2;
            
            for ibin = 1:num_bins
                bin_window_indices = ibin + n_bins_window;
                bin_window_indices = bin_window_indices(bin_window_indices >= 1 & bin_window_indices <= num_bins);
                
                area1_reduced_pooled = area_activity{ianimal}.(area1_name).area_activity_reduced(:, bin_window_indices, valid_trials);
                area2_reduced_pooled = area_activity{ianimal}.(area2_name).area_activity_reduced(:, bin_window_indices, valid_trials);
                
                num_bins_in_window = size(area1_reduced_pooled, 2);
                num_valid_trials_in_window = size(area1_reduced_pooled, 3);
                
                % Check for sufficient data
                if num_bins_in_window < (n_components_reduced + max_shift_bins) || num_valid_trials_in_window < n_components_reduced
                    warning('Skipping bin %d for animal %d (%s-%s): not enough valid data', ibin, ianimal, area1_name, area2_name);
                    continue;
                end
                try
                    % --- Analysis 3a: Correlation vs. Bin (for 'v' plot) ---
                    [~,~,r,~,~] = canoncorr(area1_reduced_pooled(:, :)', area2_reduced_pooled(:, :)');
                    n_r = min(length(r), n_components_reduced);
                    cca_all_bins(1:n_r, ibin) = r(1:n_r);
                    
                    % Shuffle by TRIAL (dim 3)
                    for ishuffle = 1:n_shuffles
                        area1_reduced_shuffled_by_trial = area1_reduced_pooled(:, :, randperm(num_valid_trials_in_window));
                        [~,~,r,~,~] = canoncorr(area1_reduced_shuffled_by_trial(:, :)', area2_reduced_pooled(:, :)');
                        n_r = min(length(r), n_components_reduced);
                        cca_shuffled_bins_by_trial(1:n_r, ibin, ishuffle) = r(1:n_r);
                    end
                    % --- Analysis 3b: Precession vs. Bin (for 'w' plot) ---
                    for ishift = -max_shift_bins:max_shift_bins
                        % Slicing logic on DIM 2 (bins)
                        if ishift < 0
                            bins_1 = (1-ishift):num_bins_in_window;
                            bins_2 = 1:(num_bins_in_window+ishift);
                        elseif ishift > 0
                            bins_1 = 1:(num_bins_in_window-ishift);
                            bins_2 = (1+ishift):num_bins_in_window;
                        else % ishift == 0
                            bins_1 = 1:num_bins_in_window;
                            bins_2 = 1:num_bins_in_window;
                        end
                        area1_data_for_cca = area1_reduced_pooled(:, bins_1, :);
                        area2_data_for_cca = area2_reduced_pooled(:, bins_2, :);
                        % Run CCA on real shifted data
                        [~,~,r,~,~] = canoncorr(area1_data_for_cca(:, :)', area2_data_for_cca(:, :)');
                        n_r = min(length(r), n_components_reduced);
                        cca_shifted_bins(1:n_r, ibin, ishift+1+max_shift_bins) = r(1:n_r);
                        
                        % Shuffle by BIN (dim 2)
                        for ishuffle = 1:n_shuffles
                            % Shuffle bins (dim 2)
                            area1_shuffled_pooled_by_bin = area1_reduced_pooled(:, randperm(num_bins_in_window), :);
                            % Apply same slice to shuffled data
                            area1_shuffled_for_cca = area1_shuffled_pooled_by_bin(:, bins_1, :);
                            [~,~,r,~,~] = canoncorr(area1_shuffled_for_cca(:, :)', area2_data_for_cca(:, :)');
                            n_r = min(length(r), n_components_reduced);
                            cca_shited_shuffled_bins_by_bin(1:n_r, ibin, ishift+1+max_shift_bins, ishuffle) = r(1:n_r);
                        end
                    end
                catch ME
                    warning('%s-%s CCA failed on bin %d for animal %d: %s', area1_name, area2_name, ibin, ianimal, ME.message);
                end 
            end
            
            % --- Plot Spatial Bin CCA (v plot) (MODIFIED) ---
            nexttile(current_v)
            plot(cca_all_bins(1, :)', 'LineWidth', 1, 'Color', 'b') % CC1
            hold on
            shadedErrorBar(1:num_bins, squeeze(mean(cca_shuffled_bins_by_trial(1, :, :), 3, "omitmissing")), squeeze(sem(cca_shuffled_bins_by_trial(1, :, :), 3)), 'lineProps', {'Color', 'b', 'LineStyle', '--'})
            plot(cca_all_bins(2, :)', 'LineWidth', 1, 'Color', 'r') % CC2
            shadedErrorBar(1:num_bins, squeeze(mean(cca_shuffled_bins_by_trial(2, :, :), 3, "omitmissing")), squeeze(sem(cca_shuffled_bins_by_trial(2, :, :), 3)), 'lineProps', {'Color', 'r', 'LineStyle', '--'})
            legend({'CC1', 'CC1 Shuffle', 'CC2', 'CC2 Shuffle'}, 'Location', 'best')
            title(sprintf('%s-%s - animal %d', area1_name, area2_name, ianimal))
            axis tight
            ylim([0, 0.7])
            xline(landmarks(:))
            
            % --- Calculate and Plot Spatial Precession (w plot) (MODIFIED) ---
            % (CC1)
            area2_leads_area1 = squeeze(mean(cca_shifted_bins(1, :, 1:max_shift_bins), 3));
            area1_leads_area2 = squeeze(mean(cca_shifted_bins(1, :, max_shift_bins+2:end), 3));
            spatial_precession_strength = (area1_leads_area2 - area2_leads_area1)./(area1_leads_area2 + area1_leads_area2);
            
            area2_leads_area1_shuffled = squeeze(mean(cca_shited_shuffled_bins_by_bin(1, :, 1:max_shift_bins, :), 3));
            area1_leads_area2_shuffled = squeeze(mean(cca_shited_shuffled_bins_by_bin(1, :, max_shift_bins+2:end, :), 3));
            spatial_precession_shuffled = (area1_leads_area2_shuffled - area2_leads_area1_shuffled)./(area1_leads_area2_shuffled + area2_leads_area1_shuffled);
            
            % --- ADDED: CC2 Calculation ---
            area2_leads_area1_cc2 = squeeze(mean(cca_shifted_bins(2, :, 1:max_shift_bins), 3));
            area1_leads_area2_cc2 = squeeze(mean(cca_shifted_bins(2, :, max_shift_bins+2:end), 3));
            spatial_precession_strength_cc2 = (area1_leads_area2_cc2 - area2_leads_area1_cc2)./(area1_leads_area2_cc2 + area1_leads_area2_cc2);
            
            area2_leads_area1_shuffled_cc2 = squeeze(mean(cca_shited_shuffled_bins_by_bin(2, :, 1:max_shift_bins, :), 3));
            area1_leads_area2_shuffled_cc2 = squeeze(mean(cca_shited_shuffled_bins_by_bin(2, :, max_shift_bins+2:end, :), 3));
            spatial_precession_shuffled_cc2 = (area1_leads_area2_shuffled_cc2 - area2_leads_area1_shuffled_cc2)./(area1_leads_area2_shuffled_cc2 + area2_leads_area1_shuffled_cc2);
            % -------------------------------
            
            nexttile(current_w) % Plot on the new 'w' figure
            plot(spatial_precession_strength, 'LineWidth', 1, 'Color', 'b') % CC1
            hold on
            shadedErrorBar(1:num_bins, squeeze(mean(spatial_precession_shuffled, 2, "omitmissing")), squeeze(sem(spatial_precession_shuffled, 2)), 'lineProps', {'Color', 'b', 'LineStyle', '--'})
            plot(spatial_precession_strength_cc2, 'LineWidth', 1, 'Color', 'r') % CC2
            shadedErrorBar(1:num_bins, squeeze(mean(spatial_precession_shuffled_cc2, 2, "omitmissing")), squeeze(sem(spatial_precession_shuffled_cc2, 2)), 'lineProps', {'Color', 'r', 'LineStyle', '--'})
            legend({'CC1', 'CC1 Shuffle', 'CC2', 'CC2 Shuffle'}, 'Location', 'best')
            
            title(sprintf('%s-%s - animal %d', area1_name, area2_name, ianimal))
            yline(0, '--')
            xline(landmarks(:))
            axis tight
            ylim([-0.3, 0.3])
            
            % --- 1. Accumulate bin data (all animals) (MODIFIED) ---
            % (CC1)
            group_results{ipair}.all_bins_corr = [group_results{ipair}.all_bins_corr, cca_all_bins(1, :)'];
            group_results{ipair}.all_bins_precession = [group_results{ipair}.all_bins_precession, spatial_precession_strength'];
            % --- ADDED: (CC2) ---
            group_results{ipair}.all_bins_corr_cc2 = [group_results{ipair}.all_bins_corr_cc2, cca_all_bins(2, :)'];
            group_results{ipair}.all_bins_precession_cc2 = [group_results{ipair}.all_bins_precession_cc2, spatial_precession_strength_cc2'];
            % --------------------

            % --- 2. Accumulate trial data (only for animals with a valid learning point) (MODIFIED) ---
            current_zscore_learning_point = zscore_learning_points(ianimal);
            
            % Check if learning point is valid and windows are within trial bounds
            if ~isnan(current_zscore_learning_point) && ...
               current_zscore_learning_point > 10 && ... % Ensures 'pre' window is valid
               (current_zscore_learning_point + 9) <= num_trials % Ensures 'post' window is valid
                
                % Define indices
                idx_early = 1:10;
                idx_pre   = (current_zscore_learning_point - 10) : (current_zscore_learning_point - 1);
                idx_post  = current_zscore_learning_point : (current_zscore_learning_point + 9);
                
                % Extract correlation data (CC1)
                corr_early = cca_all_trials(1, idx_early)';
                corr_pre   = cca_all_trials(1, idx_pre)';
                corr_post  = cca_all_trials(1, idx_post)';
                
                % Extract precession data (CC1)
                prec_early = precession_strength(idx_early)';
                prec_pre   = precession_strength(idx_pre)';
                prec_post  = precession_strength(idx_post)';
                
                % --- ADDED: Extract (CC2) data ---
                corr_early_cc2 = cca_all_trials(2, idx_early)';
                corr_pre_cc2   = cca_all_trials(2, idx_pre)';
                corr_post_cc2  = cca_all_trials(2, idx_post)';
                
                prec_early_cc2 = precession_strength_cc2(idx_early)';
                prec_pre_cc2   = precession_strength_cc2(idx_pre)';
                prec_post_cc2  = precession_strength_cc2(idx_post)';
                % ---------------------------------
                
                % Store data [10 x n_learning_animals] (CC1)
                group_results{ipair}.trial_corr_early = [group_results{ipair}.trial_corr_early, corr_early];
                group_results{ipair}.trial_corr_pre   = [group_results{ipair}.trial_corr_pre,   corr_pre];
                group_results{ipair}.trial_corr_post  = [group_results{ipair}.trial_corr_post,  corr_post];
                
                group_results{ipair}.trial_precession_early = [group_results{ipair}.trial_precession_early, prec_early];
                group_results{ipair}.trial_precession_pre   = [group_results{ipair}.trial_precession_pre,   prec_pre];
                group_results{ipair}.trial_precession_post  = [group_results{ipair}.trial_precession_post,  prec_post];
            
                % --- ADDED: Store (CC2) data ---
                group_results{ipair}.trial_corr_early_cc2 = [group_results{ipair}.trial_corr_early_cc2, corr_early_cc2];
                group_results{ipair}.trial_corr_pre_cc2   = [group_results{ipair}.trial_corr_pre_cc2,   corr_pre_cc2];
                group_results{ipair}.trial_corr_post_cc2  = [group_results{ipair}.trial_corr_post_cc2,  corr_post_cc2];
                
                group_results{ipair}.trial_precession_early_cc2 = [group_results{ipair}.trial_precession_early_cc2, prec_early_cc2];
                group_results{ipair}.trial_precession_pre_cc2   = [group_results{ipair}.trial_precession_pre_cc2,   prec_pre_cc2];
                group_results{ipair}.trial_precession_post_cc2  = [group_results{ipair}.trial_precession_post_cc2,  prec_post_cc2];
                % ---------------------------------
            
            end
            
        end % End of "if animal has data" block
    end % End of pair loop
end % End of animal loop
% Add axis labels to all figures
for ipair = 1:n_pairs
    area1_name = area_pairs_to_analyze{ipair, 1};
    area2_name = area_pairs_to_analyze{ipair, 2};
    
    xlabel(t_handles{ipair}, 'trials')
    ylabel(t_handles{ipair}, 'correlation')
    xlabel(s_handles{ipair}, 'cca')
    ylabel(s_handles{ipair}, 'lick ratio')
    
    xlabel(u_handles{ipair}, 'trials')
    ylabel(u_handles{ipair}, sprintf('%s-leads vs %s-leads', area1_name, area2_name))
    
    xlabel(v_handles{ipair}, 'spatial bins')
    ylabel(v_handles{ipair}, 'correlation')
    xlabel(w_handles{ipair}, 'spatial bins')
    ylabel(w_handles{ipair}, sprintf('%s-leads vs %s-leads', area1_name, area2_name))
end
%% Group-level analysis 
figure('Name', 'Group Spatial Correlation')
tiledlayout('flow')
for ipair = 1:n_pairs
    nexttile
    shadedErrorBar(1:num_bins, mean(group_results{ipair}.all_bins_corr, 2, 'omitmissing'), sem(group_results{ipair}.all_bins_corr, 2), 'lineProps', {'Color', 'b'})
    hold on
    % shadedErrorBar(1:num_bins, mean(group_results{ipair}.all_bins_corr_cc2, 2, 'omitmissing'), sem(group_results{ipair}.all_bins_corr_cc2, 2), 'lineProps', {'Color', 'r'})
    xline(landmarks(:))
    title(group_results{ipair}.pair_name)
    % legend({'CC1', 'CC2'}, 'Location', 'best')
end

figure('Name', 'Group Spatial Precession')
tiledlayout('flow')
for ipair = 1:n_pairs
    nexttile
    shadedErrorBar(1:num_bins, mean(group_results{ipair}.all_bins_precession, 2, 'omitmissing'), sem(group_results{ipair}.all_bins_precession, 2), 'lineProps', {'Color', 'b'})
    hold on
    % shadedErrorBar(1:num_bins, mean(group_results{ipair}.all_bins_precession_cc2, 2, 'omitmissing'), sem(group_results{ipair}.all_bins_precession_cc2, 2), 'lineProps', {'Color', 'r'})
    xline(landmarks(:))
    title(group_results{ipair}.pair_name)
    yline(0, '--')
    ylim([-0.3 0.3])
    % legend({'CC1', 'CC2'}, 'Location', 'best')
end

figure('Name', 'Group Trial Correlation')
tiledlayout('flow')
for ipair = 1:n_pairs
    nexttile
    hold on
    shadedErrorBar(1:10, mean(group_results{ipair}.trial_corr_early, 2, 'omitmissing'), sem(group_results{ipair}.trial_corr_early, 2), 'lineProps', {'Color', 'b'})
    shadedErrorBar(11:20, mean(group_results{ipair}.trial_corr_pre, 2, 'omitmissing'), sem(group_results{ipair}.trial_corr_pre, 2), 'lineProps', {'Color', 'b'})
    shadedErrorBar(21:30, mean(group_results{ipair}.trial_corr_post, 2, 'omitmissing'), sem(group_results{ipair}.trial_corr_post, 2), 'lineProps', {'Color', 'b'})
    
    % shadedErrorBar(1:10, mean(group_results{ipair}.trial_corr_early_cc2, 2, 'omitmissing'), sem(group_results{ipair}.trial_corr_early_cc2, 2), 'lineProps', {'Color', 'r'})
    % shadedErrorBar(11:20, mean(group_results{ipair}.trial_corr_pre_cc2, 2, 'omitmissing'), sem(group_results{ipair}.trial_corr_pre_cc2, 2), 'lineProps', {'Color', 'r'})
    % shadedErrorBar(21:30, mean(group_results{ipair}.trial_corr_post_cc2, 2, 'omitmissing'), sem(group_results{ipair}.trial_corr_post_cc2, 2), 'lineProps', {'Color', 'r'})
    
    xline([0.5, 10.5, 20.5], 'k', {'First', 'Pre', 'Post'})
    title(group_results{ipair}.pair_name)
    ylim([0, 0.6])
    % legend({'CC1', 'CC2'}, 'Location', 'best')
end

figure('Name', 'Group Trial Precession')
tiledlayout('flow')
for ipair = 1:n_pairs
    nexttile
    hold on
    shadedErrorBar(1:10, mean(group_results{ipair}.trial_precession_early, 2, 'omitmissing'), sem(group_results{ipair}.trial_precession_early, 2), 'lineProps', {'Color', 'b'})
    shadedErrorBar(11:20, mean(group_results{ipair}.trial_precession_pre, 2, 'omitmissing'), sem(group_results{ipair}.trial_precession_pre, 2), 'lineProps', {'Color', 'b'})
    shadedErrorBar(21:30, mean(group_results{ipair}.trial_precession_post, 2, 'omitmissing'), sem(group_results{ipair}.trial_precession_post, 2), 'lineProps', {'Color', 'b'})
    
    % shadedErrorBar(1:10, mean(group_results{ipair}.trial_precession_early_cc2, 2, 'omitmissing'), sem(group_results{ipair}.trial_precession_early_cc2, 2), 'lineProps', {'Color', 'r'})
    % shadedErrorBar(11:20, mean(group_results{ipair}.trial_precession_pre_cc2, 2, 'omitmissing'), sem(group_results{ipair}.trial_precession_pre_cc2, 2), 'lineProps', {'Color', 'r'})
    % shadedErrorBar(21:30, mean(group_results{ipair}.trial_precession_post_cc2, 2, 'omitmissing'), sem(group_results{ipair}.trial_precession_post_cc2, 2), 'lineProps', {'Color', 'r'})
    % 
    xline([0.5, 10.5, 20.5], 'k', {'First', 'Pre', 'Post'})
    title(group_results{ipair}.pair_name)
    ylim([-0.3, 0.3])
    % legend({'CC1', 'CC2'}, 'Location', 'best')
    yline(0, '--')
end
%% Epoch average CCA
figure('Name', 'Group Epoch Correlation (bar)')
tiledlayout('flow')
for ipair = 1:n_pairs
    nexttile
    hold on
    % Combine CC1 and CC2 data
    data_cc1 = [mean(group_results{ipair}.trial_corr_early, 1, 'omitmissing')', ...
                mean(group_results{ipair}.trial_corr_pre, 1, 'omitmissing')', ...
                mean(group_results{ipair}.trial_corr_post, 1, 'omitmissing')'];
    data_cc2 = [mean(group_results{ipair}.trial_corr_early_cc2, 1, 'omitmissing')', ...
                mean(group_results{ipair}.trial_corr_pre_cc2, 1, 'omitmissing')', ...
                mean(group_results{ipair}.trial_corr_post_cc2, 1, 'omitmissing')'];
            
    my_errorbar_plot([data_cc1])
    % xticks(1:6)
    % xticklabels({'E-CC1', 'P-CC1', 'Po-CC1', 'E-CC2', 'P-CC2', 'Po-CC2'})

    xticks(1:3)
    xticklabels({'Early', 'Pre', 'Post'})
end

figure('Name', 'Group Epoch Precession (bar)')
tiledlayout('flow')
for ipair = 1:n_pairs
    nexttile
    hold on
    % Combine CC1 and CC2 data
    data_cc1 = [mean(group_results{ipair}.trial_precession_early, 1, 'omitmissing')', ...
                mean(group_results{ipair}.trial_precession_pre, 1, 'omitmissing')', ...
                mean(group_results{ipair}.trial_precession_post, 1, 'omitmissing')'];
    data_cc2 = [mean(group_results{ipair}.trial_precession_early_cc2, 1, 'omitmissing')', ...
                mean(group_results{ipair}.trial_precession_pre_cc2, 1, 'omitmissing')', ...
                mean(group_results{ipair}.trial_precession_post_cc2, 1, 'omitmissing')'];
            
    % my_errorbar_plot([data_cc1, data_cc2])
    my_errorbar_plot([data_cc1])
    % xticks(1:6)
    % xticklabels({'E-CC1', 'P-CC1', 'Po-CC1', 'E-CC2', 'P-CC2', 'Po-CC2'})

    xticks(1:3)
    xticklabels({'Early', 'Pre', 'Post'})

    title(group_results{ipair}.pair_name)
end
%% (MODIFIED)
% First plot: Correlation data
figure('Name', 'Group Epoch Corr Stats (plot)')
tiledlayout('flow')
for ipair = 1:n_pairs
    nexttile
    hold on
    
    % Prepare data (CC1)
    data_early_cc1 = group_results{ipair}.trial_corr_early;
    data_pre_cc1 = group_results{ipair}.trial_corr_pre;
    data_post_cc1 = group_results{ipair}.trial_corr_post;
    
    % Get means for plotting (CC1 ONLY)
    plot_data = [mean(data_early_cc1, 1, 'omitmissing')', ...
                 mean(data_pre_cc1, 1, 'omitmissing')', ...
                 mean(data_post_cc1, 1, 'omitmissing')'];
    
    my_errorbar_plot(plot_data)
    
    % --- Perform repeated measures ANOVA (ON CC1 DATA ONLY) ---
    data_early_means = mean(data_early_cc1, 1, 'omitmissing')';
    data_pre_means = mean(data_pre_cc1, 1, 'omitmissing')';
    data_post_means = mean(data_post_cc1, 1, 'omitmissing')';
    
    % Organize data: rows = subjects, columns = time points
    rm_data = [data_early_means, data_pre_means, data_post_means];
    
    % Remove rows with any NaN values
    valid_rows = ~any(isnan(rm_data), 2);
    rm_data_clean = rm_data(valid_rows, :);
    
    if size(rm_data_clean, 1) >= 3  % Need at least 3 subjects
        % Create table for repeated measures
        t = array2table(rm_data_clean, 'VariableNames', {'Early', 'Pre', 'Post'});
        
        % Fit repeated measures model
        Time = [1 2 3]';
        rm = fitrm(t, 'Early-Post~1', 'WithinDesign', Time);
        
        % Run repeated measures ANOVA
        ranova_results = ranova(rm);
        p_value = ranova_results.pValue(1);
        
        % --- IMPROVED: Post-hoc pairwise comparisons ---
        mct = multcompare(rm, 'Time', 'ComparisonType', 'bonferroni');
        
        % Extract p-values for sigstar
        % Find rows for 1vs2, 2vs3, 1vs3
        p_early_pre = mct.pValue((mct.Time_1 == 1 & mct.Time_2 == 2));
        p_pre_post = mct.pValue((mct.Time_1 == 2 & mct.Time_2 == 3));
        p_early_post = mct.pValue((mct.Time_1 == 1 & mct.Time_2 == 3));
        
        % Use sigstar to display significant comparisons (on CC1 positions)
        groups = {};
        pvals = [];
        
        if ~isempty(p_early_pre) & p_early_pre < 0.05
            groups{end+1} = [1, 2];
            pvals(end+1) = p_early_pre;
        end
        
        if ~isempty(p_pre_post) & p_pre_post < 0.05
            groups{end+1} = [2, 3];
            pvals(end+1) = p_pre_post;
        end
        
        if ~isempty(p_early_post) & p_early_post < 0.05
            groups{end+1} = [1, 3];
            pvals(end+1) = p_early_post;
        end
        
        if ~isempty(groups)
            sigstar(groups, pvals);
        end
        
        % Update title with ANOVA p-value
        title(sprintf('%s (RM-ANOVA: p=%.4f)', group_results{ipair}.pair_name, p_value))
        
    else
        title(sprintf('%s (insufficient data)', group_results{ipair}.pair_name))
    end
    
    % --- FIX: Correct xticks for 3 columns ---
    xticks(1:3)
    xticklabels({'Early', 'Pre', 'Post'})
end

% Second plot: Precession data
figure('Name', 'Group Epoch Precession Stats (plot)')
tiledlayout('flow')
for ipair = 1:n_pairs
    nexttile
    hold on
    
    % Prepare data (CC1)
    data_early_cc1 = group_results{ipair}.trial_precession_early;
    data_pre_cc1 = group_results{ipair}.trial_precession_pre;
    data_post_cc1 = group_results{ipair}.trial_precession_post;
    
    % Get means for plotting (CC1 ONLY)
    plot_data = [mean(data_early_cc1, 1, 'omitmissing')', ...
                 mean(data_pre_cc1, 1, 'omitmissing')', ...
                 mean(data_post_cc1, 1, 'omitmissing')'];
    
    my_errorbar_plot(plot_data)
    
    % --- Perform repeated measures ANOVA (ON CC1 DATA ONLY) ---
    data_early_means = mean(data_early_cc1, 1, 'omitmissing')';
    data_pre_means = mean(data_pre_cc1, 1, 'omitmissing')';
    data_post_means = mean(data_post_cc1, 1, 'omitmissing')';
    
    rm_data = [data_early_means, data_pre_means, data_post_means];
    valid_rows = ~any(isnan(rm_data), 2);
    rm_data_clean = rm_data(valid_rows, :);
    
    if size(rm_data_clean, 1) >= 3
        t = array2table(rm_data_clean, 'VariableNames', {'Early', 'Pre', 'Post'});
        Time = [1 2 3]';
        rm = fitrm(t, 'Early-Post~1', 'WithinDesign', Time);
        ranova_results = ranova(rm);
        p_value = ranova_results.pValue(1);
        
        % --- IMPROVED: Post-hoc pairwise comparisons ---
        mct = multcompare(rm, 'Time', 'ComparisonType', 'bonferroni');
        
        % Extract p-values for sigstar
        p_early_pre = mct.pValue((mct.Time_1 == 1 & mct.Time_2 == 2));
        p_pre_post = mct.pValue((mct.Time_1 == 2 & mct.Time_2 == 3));
        p_early_post = mct.pValue((mct.Time_1 == 1 & mct.Time_2 == 3));
        
        % Use sigstar to display significant comparisons
        groups = {};
        pvals = [];
        
        if ~isempty(p_early_pre) & p_early_pre < 0.05
            groups{end+1} = [1, 2];
            pvals(end+1) = p_early_pre;
        end
        
        if ~isempty(p_pre_post) & p_pre_post < 0.05
            groups{end+1} = [2, 3];
            pvals(end+1) = p_pre_post;
        end
        
        if ~isempty(p_early_post) & p_early_post < 0.05
            groups{end+1} = [1, 3];
            pvals(end+1) = p_early_post;
        end
        
        if ~isempty(groups)
            sigstar(groups, pvals);
        end
        
        title(sprintf('%s (RM-ANOVA: p=%.4f)', group_results{ipair}.pair_name, p_value))
    else
        title(sprintf('%s (insufficient data)', group_results{ipair}.pair_name))
    end
    
    % --- FIX: Correct xticks for 3 columns ---
    xticks(1:3)
    xticklabels({'Early', 'Pre', 'Post'})
end
%% Group-level Statistics (Trial-Aligned) (CC1)
fprintf('--- Group Statistics (CC1) (n=%d animals) ---\n', size(group_results{1}.trial_corr_early, 2));
for ipair = 1:n_pairs
    pair_name = group_results{ipair}.pair_name;
    
    % --- Correlation Stats ---
    % Take the mean correlation in each epoch for each animal
    early_corr_means = mean(group_results{ipair}.trial_corr_early, 1, 'omitmissing');
    pre_corr_means = mean(group_results{ipair}.trial_corr_pre, 1, 'omitmissing');
    post_corr_means = mean(group_results{ipair}.trial_corr_post, 1, 'omitmissing');
    
    % Paired test: Early vs Pre-learning 
    [p_corr, ~, stats_corr] = signrank(early_corr_means, pre_corr_means);
    fprintf('%s Correlation (Early vs Pre): p = %.4f (W = %d)\n', ...
        pair_name, p_corr, stats_corr.signedrank);
    % Paired test: Pre-learning vs. Post-learning
    [p_corr, ~, stats_corr] = signrank(pre_corr_means, post_corr_means);
    fprintf('%s Correlation (Pre vs Post): p = %.4f (W = %d)\n', ...
        pair_name, p_corr, stats_corr.signedrank);
    % Paired test: Early vs. Post-learning
    [p_corr, ~, stats_corr] = signrank(early_corr_means, post_corr_means);
    fprintf('%s Correlation (Early vs Post): p = %.4f (W = %d)\n', ...
        pair_name, p_corr, stats_corr.signedrank);
    % --- Precession Stats ---
    early_prec_means = mean(group_results{ipair}.trial_precession_early, 1, 'omitmissing');
    pre_prec_means = mean(group_results{ipair}.trial_precession_pre, 1, 'omitmissing');
    post_prec_means = mean(group_results{ipair}.trial_precession_post, 1, 'omitmissing');
    % Paired test: Early vs. Pre-learning
    [p_prec, ~, stats_prec] = signrank(early_prec_means, pre_prec_means);
    fprintf('%s Precession (Early vs Pre): p = %.4f (W = %d)\n', ...
        pair_name, p_prec, stats_prec.signedrank);
    
    % Paired test: Pre-learning vs. Post-learning
    [p_prec, ~, stats_prec] = signrank(pre_prec_means, post_prec_means);
    fprintf('%s Precession (Pre vs Post): p = %.4f (W = %d)\n', ...
        pair_name, p_prec, stats_prec.signedrank);
    % Paired test: Early vs. Post-learning
    [p_prec, ~, stats_prec] = signrank(early_prec_means, post_prec_means);
    fprintf('%s Precession (Early vs Post): p = %.4f (W = %d)\n', ...
        pair_name, p_prec, stats_prec.signedrank);
       
end

%% Group-level Statistics (Trial-Aligned) (CC2) ---
fprintf('--- Group Statistics (CC2) (n=%d animals) ---\n', size(group_results{1}.trial_corr_early_cc2, 2));
for ipair = 1:n_pairs
    pair_name = group_results{ipair}.pair_name;
    
    % --- Correlation Stats ---
    early_corr_means = mean(group_results{ipair}.trial_corr_early_cc2, 1, 'omitmissing');
    pre_corr_means = mean(group_results{ipair}.trial_corr_pre_cc2, 1, 'omitmissing');
    post_corr_means = mean(group_results{ipair}.trial_corr_post_cc2, 1, 'omitmissing');
    
    [p_corr, ~, stats_corr] = signrank(early_corr_means, pre_corr_means);
    fprintf('%s CC2 Correlation (Early vs Pre): p = %.4f (W = %d)\n', ...
        pair_name, p_corr, stats_corr.signedrank);
    [p_corr, ~, stats_corr] = signrank(pre_corr_means, post_corr_means);
    fprintf('%s CC2 Correlation (Pre vs Post): p = %.4f (W = %d)\n', ...
        pair_name, p_corr, stats_corr.signedrank);
    [p_corr, ~, stats_corr] = signrank(early_corr_means, post_corr_means);
    fprintf('%s CC2 Correlation (Early vs Post): p = %.4f (W = %d)\n', ...
        pair_name, p_corr, stats_corr.signedrank);
        
    % --- Precession Stats ---
    early_prec_means = mean(group_results{ipair}.trial_precession_early_cc2, 1, 'omitmissing');
    pre_prec_means = mean(group_results{ipair}.trial_precession_pre_cc2, 1, 'omitmissing');
    post_prec_means = mean(group_results{ipair}.trial_precession_post_cc2, 1, 'omitmissing');

    [p_prec, ~, stats_prec] = signrank(early_prec_means, pre_prec_means);
    fprintf('%s CC2 Precession (Early vs Pre): p = %.4f (W = %d)\n', ...
        pair_name, p_prec, stats_prec.signedrank);
    [p_prec, ~, stats_prec] = signrank(pre_prec_means, post_prec_means);
    fprintf('%s CC2 Precession (Pre vs Post): p = %.4f (W = %d)\n', ...
        pair_name, p_prec, stats_prec.signedrank);
    [p_prec, ~, stats_prec] = signrank(early_prec_means, post_prec_means);
    fprintf('%s CC2 Precession (Early vs Post): p = %.4f (W = %d)\n', ...
        pair_name, p_prec, stats_prec.signedrank);
       
end
% -----------------------------------------------------------

%% Save Results
disp('Saving results...')
my_save_path = '/Users/theoamvr/Desktop/Experiments/TomLearning';
% 1. Save the group data structure
save(fullfile(my_save_path, 'cca_group_results.mat'), 'group_results');
% 2. Save all the figures
all_fig_handles = [t_handles; s_handles; u_handles; v_handles; w_handles];
all_fig_names = cellfun(@(h) h.Parent.Name, all_fig_handles, 'UniformOutput', false);
for iFig = 1:length(all_fig_handles)
    fig_name = matlab.lang.makeValidName(all_fig_names{iFig}); % Clean up name for file
    savefig(all_fig_handles{iFig}.Parent, fullfile(my_save_path, [fig_name '.fig']));
    exportgraphics(all_fig_handles{iFig}.Parent, fullfile(my_save_path, [fig_name '.png']), 'Resolution', 300);
end
disp('Done.');