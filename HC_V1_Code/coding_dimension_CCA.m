%% Setup and Parameters
% Parameters for xcorr
max_lag_bins = 10; 
x_axis_vals = 1:30; 

% Dimensionality Reduction Parameters
n_pca_components = 5; % User requested 3-5. 3 is safest for 20 trials.
min_units = 5;        % Minimum neurons required to attempt analysis

% Storage
aligned_zero_lag = cell(1, n_pairs);
aligned_peak_lag = cell(1, n_pairs);

%% Main Analysis Loop
for ianimal = 1:n_animals
    
    % --- Safety Check: Learning Point ---
    if isnan(zscore_learning_points(ianimal)) || zscore_learning_points(ianimal) <= 11
        continue;
    end
    
    % --- Data Prep ---
    % Get data: (Units x Bins x Trials)
    animal_data = permute(analysis_spatial{ianimal}.firing.cued.freq_z, [1, 3, 2]);
    [~, num_bins, num_trials] = size(animal_data);
    
    animal_areas = units{ianimal}.regions_label;
    n_areas = numel(animal_areas);
    unit_areas = units{ianimal}.idx;
    
    % Define trial indices for Axis Definition
    idx_naive = 1:10;
    idx_expert = (zscore_learning_points(ianimal)+1):(zscore_learning_points(ianimal)+10);
    
    current_animal_proj = struct();
    
    % --- 1. PCA -> LDA Projection ---
    for iarea = 1:n_areas
        current_area_name = animal_areas(iarea);
        current_area_activity = animal_data(unit_areas(iarea, :), :, :);
        n_units_area = size(current_area_activity, 1);
        
        % Check: Enough neurons?
        if n_units_area < min_units
            fprintf('Skipping %s (Animal %d): Only %d units.\n', current_area_name, ianimal, n_units_area);
            continue;
        end
        
        % A. Prepare Data for Dimensionality Reduction
        % We use the 20 definition trials to find the subspace
        % Average over time bins to get 1 point per trial
        data_naive = permute(mean(current_area_activity(:, :, idx_naive), 2, 'omitmissing'), [3, 1, 2]); % 10 x Units
        data_expert = permute(mean(current_area_activity(:, :, idx_expert), 2, 'omitmissing'), [3, 1, 2]); % 10 x Units
        
        X_def = [data_naive; data_expert]; % (20 x Units)
        Y_def = [ones(10, 1); 2*ones(10, 1)]; % Labels
        
        % B. PCA Reduction
        % Calculate PCA on the 20 trials (Naive + Expert)
        % We want to find the dominant modes of variance in these states
        
        % Handle case where we requested 3 components but have e.g. 5 units
        % (PCA requires N_obs-1 max, but here N_units limits the rank usually)
        k_components = min([n_pca_components, n_units_area, size(X_def, 1)-1]);
        
        [coeff, score, ~] = pca(X_def, 'NumComponents', k_components);
        % coeff: (Units x k) -> The projection matrix to PC space
        % score: (20 x k)    -> The data in PC space
        
        % C. LDA on PC Scores
        % Now we find the separation axis inside the PC space
        try
            % We use 'linear' because we are now Low-D (k <= 5) and N=20.
            % Covariance is invertible.
            Mdl = fitcdiscr(score, Y_def, 'DiscrimType', 'linear');
            
            % Get the linear coefficients (normal vector to the boundary)
            % This is a (k x 1) vector
            lda_weights_pc = Mdl.Coeffs(1, 2).Linear; 
            
        catch
            warning('LDA failed for area %s. Skipping.', current_area_name);
            continue;
        end
        
        % D. Global Projection
        % We need to project ALL data (not just the 20 trials)
        % Path: Raw (Units) -> PCA (k) -> LDA (1)
        
        % 1. Flatten all data for projection: (Units x TimePoints)
        activity_to_project = reshape(current_area_activity, n_units_area, num_bins*num_trials);
        
        % 2. Project to PC space: (k x Units)' * (Units x TimePoints) -> (k x TimePoints)
        data_pc = coeff' * activity_to_project;
        
        % 3. Project to LDA axis: (1 x k) * (k x TimePoints) -> (1 x TimePoints)
        projected_flat = lda_weights_pc' * data_pc;
        
        % 4. Reshape back
        projected_activity = reshape(projected_flat, num_bins, num_trials);
        
        current_animal_proj.(string(current_area_name)) = projected_activity;
    end
    
    % --- 2. Pairwise Correlation (xcorr) ---
    % (Logic remains identical to previous version)
    
    for ipair = 1:n_pairs
        area1_name = area_pairs_to_analyze{ipair, 1};
        area2_name = area_pairs_to_analyze{ipair, 2};
        
        if ~isfield(current_animal_proj, area1_name) || ~isfield(current_animal_proj, area2_name)
            continue;
        end
        
        data1 = current_animal_proj.(area1_name);
        data2 = current_animal_proj.(area2_name);
        
        pair_zero_lag = nan(num_trials, 1);
        pair_peak_lag = nan(num_trials, 1);
        
        for t = 1:num_trials
            if any(isnan(data1(:, t))) || any(isnan(data2(:, t)))
                continue;
            end
            
            [xc, lags] = xcorr(data1(:, t), data2(:, t), max_lag_bins, 'coeff');
            
            % Zero Lag
            pair_zero_lag(t) = xc(lags == 0);
            
            % Peak Lag
            [~, idx_max] = max(abs(xc)); 
            pair_peak_lag(t) = lags(idx_max);
        end
        
        % --- 3. Align Data ---
        lp = zscore_learning_points(ianimal);
        idx_transition = (lp - 10):(lp - 1);
        idx_expert_plot = (lp + 1):(lp + 10);
        idx_naive_plot = 1:10;
        
        aligned_zero_lag{ipair}(ianimal, :) = [ ...
             pair_zero_lag(idx_naive_plot)', ...
             pair_zero_lag(idx_transition)', ...
             pair_zero_lag(idx_expert_plot)' ];
         
        aligned_peak_lag{ipair}(ianimal, :) = [ ...
             pair_peak_lag(idx_naive_plot)', ...
             pair_peak_lag(idx_transition)', ...
             pair_peak_lag(idx_expert_plot)' ];
    end
end

%% Visualization
% Figure 1: Zero Lag Correlation
figure('Name', 'PCA-LDA Projected Correlation', 'Color', 'w');
for ipair = 1:n_pairs
    area1_name = area_pairs_to_analyze{ipair, 1};
    area2_name = area_pairs_to_analyze{ipair, 2};
    
    subplot(2, 4, ipair)
    hold on
    
    data_matrix = aligned_zero_lag{ipair};
    data_matrix(all(isnan(data_matrix), 2), :) = [];
    
    if ~isempty(data_matrix)
        shadedErrorBar(x_axis_vals, mean(data_matrix, 1), sem(data_matrix, 1), 'lineProps', 'b');
    end
    
    xline(10.5, '--k'); xline(20.5, '--k');
    axis tight
    title(sprintf('%s vs %s', area1_name, area2_name));
    ylabel('Correlation');
end

% Figure 2: Peak Lag
figure('Name', 'PCA-LDA Lag Dynamics', 'Color', 'w');
for ipair = 1:n_pairs
    area1_name = area_pairs_to_analyze{ipair, 1};
    area2_name = area_pairs_to_analyze{ipair, 2};
    
    subplot(2, 4, ipair)
    hold on
    
    lag_matrix = aligned_peak_lag{ipair};
    lag_matrix(all(isnan(lag_matrix), 2), :) = [];
    
    if ~isempty(lag_matrix)
        shadedErrorBar(x_axis_vals, mean(lag_matrix, 1), sem(lag_matrix, 1), 'lineProps', 'r');
    end
    
    yline(0, '-', 'Color', [0.5 0.5 0.5]); 
    xline(10.5, '--k'); xline(20.5, '--k');
    axis tight
    ylim([-max_lag_bins max_lag_bins]);
    ylabel('Lag (Bins)');
    title(sprintf('Lag: %s vs %s', area1_name, area2_name));
end