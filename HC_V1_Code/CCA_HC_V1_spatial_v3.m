%% CCA_HC_V1_spatial_v3.m
% Spatial CCA driver, v5 spec.
%
% Composes the v5_* primitives into the H&H-adapted CCA pipeline for the
% Hippocampus-V1 project, spatial alignment mode (200 position bins × 2.5
% cm/bin = 5 m corridor). One fit per (animal × area-pair × epoch) with
% epochs = {naive, pre-LP, post-LP} as 10-trial windows.
%
% Spec: /Users/theoamvr/Documents/ResearchVault/Methods/CCA_HH_Adapted.md
%       See §6.1 (project-specific deviations) and §6.1.x (data IO contract).
%
% Per spec §3.2, k is fixed within an animal × area-pair across the three
% epochs. Per spec §6.1, region tagging derives unit→area from
% `units.idx` + `units.regions_label`; FS units are excluded from cortical
% areas.

clear; clc; close all;

%% 1. CONFIG

cfg = struct();

% Paths
cfg.base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
cfg.data_dir = fullfile(cfg.base_dir, 'HC_V1_data');
cfg.file_pattern = 'TF*_export.mat';
cfg.learning_file = 'animal_behaviour.mat';

% Output
cfg.current_date = datestr(now, 'yyyy_mm_dd');
cfg.save_path = fullfile(cfg.data_dir, ...
    sprintf('Spatial_CCA_v3_%s.mat', cfg.current_date));

% Spatial binning
cfg.n_position_bins = 200;
cfg.bin_size_cm = 2.5;

% Epoch definition (spec §6.1)
cfg.n_trials_per_epoch = 10;     % naive = first 10; pre-LP = 10 before LP; post-LP = 10 after LP

% Unit selection
cfg.min_units_per_region = 5;
cfg.exclude_fs_areas = {'V1', 'RSC', 'CA1', 'CA3'};   % exclude FS units in these areas

% PCA / CCA k rule (spec §3.2)
cfg.k_cap = 30;
cfg.k_samples_rule = 50;        % k <= n_samples / k_samples_rule per side
cfg.k_variance_target = 0.80;    % optional: explained-variance mode, capped at k_cap
cfg.k_mode = 'variance';          % 'samples' (default, per spec) | 'variance'

% Shuffle null (spec §3.1)
cfg.n_shuffles = 25;
cfg.cv_splits = 5;

% Lagged-refit IFI (spec §3.1, adapted for spatial)
% Note: spatial-mode "lag" is position-bin shift, not time. Legacy v2 used
% ±3 position bins (= ±7.5 cm). The spec's ±200 ms is for temporal mode;
% the spatial analogue lives here as max_lag in position bins.
cfg.max_lag_bins = 5;
cfg.central_window_bins = 3;
cfg.smooth_size = 3;
cfg.smooth_sigma = 1.0;

% Projection IFI
cfg.proj_max_lag_bins = 5;
cfg.proj_min_paired_samples = 5;

% Auxiliary analyses (Gonzalez et al. paper, adapted to v5 + plain CCA).
cfg.sig_alpha = 0.05;             % per-dim significance via 95th percentile of null
cfg.n_spatial_bins_entropy = 5;   % paper convention for weight spatial entropy
cfg.n_surrogates_entropy = 100;
cfg.n_surrogates_proj_corr = 100;
cfg.aux_seed = 0;                 % RNG seed for circular-shift surrogates

% Area pairs (per spec §6.1)
cfg.area_pairs = {
    'CA1', 'V1';
    'CA1', 'CA3';
    'CA1', 'DG';
    'CA1', 'RSC';
    'CA1', 'SUB';
    'V1',  'RSC';
    'RSC', 'SUB';
    'CA3', 'DG';
};
n_pairs = size(cfg.area_pairs, 1);

% Reproducibility
cfg.cca_fold_seed = 0;
cfg.shuffle_seed = 0;

% Make the path visible
addpath(fullfile(cfg.base_dir, 'HC_V1_Code'));

fprintf('CCA_HC_V1_spatial_v3 — %s\n', cfg.current_date);
fprintf('  Output: %s\n', cfg.save_path);

%% 2. LOAD ANIMAL FILE LIST + LEARNING POINTS

file_list = dir(fullfile(cfg.data_dir, cfg.file_pattern));
n_animals = length(file_list);
fprintf('  Found %d animal files.\n', n_animals);

learning_points = load_learning_points(...
    fullfile(cfg.data_dir, cfg.learning_file), n_animals);
is_learner = ~isnan(learning_points);
mean_lp = round(mean(learning_points(is_learner), 'omitnan'));
% Yoked: assign mean LP to non-learners.
analysis_lp = learning_points;
analysis_lp(~is_learner) = mean_lp;
fprintf('  Learners: %d/%d. Mean LP (yoke): %d\n', ...
    sum(is_learner), n_animals, mean_lp);

%% 3. INITIALISE RESULTS

group_results = init_results_struct(cfg, n_animals, n_pairs);

%% 4. PER-ANIMAL LOOP

for ianimal = 1:n_animals
    filename = file_list(ianimal).name;
    fullpath = fullfile(cfg.data_dir, filename);
    fprintf('\n[%d/%d] Processing %s (LP=%d, learner=%d)\n', ...
        ianimal, n_animals, filename, analysis_lp(ianimal), is_learner(ianimal));

    try
        animal = load_animal(fullpath, cfg);
    catch err
        warning('  Failed to load %s: %s', filename, err.message);
        continue;
    end
    if isempty(animal.area_tensors)
        fprintf('  No valid areas after unit selection. Skipping.\n');
        continue;
    end

    % Define epoch trial indices for this animal.
    n_trials = animal.n_trials;
    lp = analysis_lp(ianimal);
    epochs = epoch_indices_for_animal(n_trials, lp, cfg.n_trials_per_epoch);
    if isempty(epochs)
        fprintf('  No valid epochs (trial count too low or LP out of range). Skipping.\n');
        continue;
    end

    % Per area-pair loop.
    for ipair = 1:n_pairs
        a1 = cfg.area_pairs{ipair, 1};
        a2 = cfg.area_pairs{ipair, 2};
        if ~isfield(animal.area_tensors, a1) || ~isfield(animal.area_tensors, a2)
            continue;
        end
        S_x = animal.area_tensors.(a1);   % (n_trials × n_pos_bins × n_units_a1)
        S_y = animal.area_tensors.(a2);   % (n_trials × n_pos_bins × n_units_a2)
        n_units_x = size(S_x, 3);
        n_units_y = size(S_y, 3);

        % Decide k for this animal × pair, fixed across the three epochs.
        % Use the smallest epoch's sample count to set the cap.
        smallest_epoch_size = min(cellfun(@(e) numel(e), {epochs.trials_idx}));
        n_samples_smallest = smallest_epoch_size * cfg.n_position_bins;
        k_chosen = choose_k(cfg, n_units_x, n_units_y, n_samples_smallest, ...
                            S_x, S_y, epochs(1).trials_idx);

        % Per-epoch fit. Cache the CELL-SPACE canonical loadings from the
        % previous epoch for the principal-angles diagnostic. Comparing
        % res.A (PC-space, k×k) directly is meaningless — both bases
        % trivially span the entire k-D ambient and the diagnostic
        % collapses to zero. Composing PCA components with A puts the
        % comparison in cell space, which is shared across epochs.
        prev_W_cells_A = []; prev_W_cells_B = [];
        speed_per_trial_bin = animal.speed_cm_s;
        for iep = 1:numel(epochs)
            ep = epochs(iep);
            tr_idx = ep.trials_idx;
            S_x_ep = S_x(tr_idx, :, :);
            S_y_ep = S_y(tr_idx, :, :);
            speed_ep = speed_per_trial_bin(tr_idx, :);

            res = fit_one(S_x_ep, S_y_ep, k_chosen, cfg, speed_ep, tr_idx);

            % Principal angles in cell space.
            r_dims = size(res.W_cells_A, 2);
            if ~isempty(prev_W_cells_A)
                res.principal_angles_A = v5_principal_angles( ...
                    prev_W_cells_A, res.W_cells_A);
                res.principal_angles_B = v5_principal_angles( ...
                    prev_W_cells_B, res.W_cells_B);
            else
                res.principal_angles_A = nan(r_dims, 1);
                res.principal_angles_B = nan(r_dims, 1);
            end
            prev_W_cells_A = res.W_cells_A;
            prev_W_cells_B = res.W_cells_B;

            % Store into group_results.
            group_results(ipair).per_epoch(ianimal, iep) = res;
        end
        fprintf('  pair %s-%s: k=%d, samples/var=%.1f, ', ...
            a1, a2, k_chosen, n_samples_smallest / (2 * k_chosen));
        cc_vec = arrayfun(@(e) ...
            group_results(ipair).per_epoch(ianimal, e.iep).cc_cv(1), ...
            struct('iep', num2cell(1:numel(epochs))));
        fprintf('cc_cv across epochs: %s\n', mat2str(cc_vec, 3));
    end
end

%% 5. SAVE

save(cfg.save_path, 'group_results', 'cfg', 'learning_points', ...
    'analysis_lp', 'is_learner', '-v7.3');
fprintf('\nSaved to %s\n', cfg.save_path);


%% Helpers

function learning_points = load_learning_points(lp_path, n_animals)
    if exist(lp_path, 'file')
        dat_lp = load(lp_path);
        learning_points = dat_lp.period_experienced(:, 1);
        if isfield(dat_lp, 'animal_id')
            [~, sorting_idx] = sort(dat_lp.animal_id);
            learning_points = learning_points(sorting_idx);
        end
        if length(learning_points) < n_animals
            learning_points(end+1:n_animals) = NaN;
        end
    else
        warning('Learning file not found at %s. All NaN.', lp_path);
        learning_points = nan(n_animals, 1);
    end
end


function animal = load_animal(fullpath, cfg)
% Load one TF*_export.mat and build (n_trials × n_pos_bins × n_units) per area.
    D = load(fullpath);
    if ~isfield(D, 'units') || ~isfield(D, 'analysis_spatial') || ...
       ~isfield(D.analysis_spatial, 'firing') || ...
       ~isfield(D.analysis_spatial.firing, 'cued') || ...
       ~isfield(D.analysis_spatial.firing.cued, 'freq_z')
        error('Missing units / analysis_spatial.firing.cued.freq_z');
    end
    units = D.units;
    raw_spatial = D.analysis_spatial.firing.cued.freq_z;   % (n_units × n_trials × n_pos_bins)
    [n_units, n_trials, n_pos_bins] = size(raw_spatial);
    if n_pos_bins ~= cfg.n_position_bins
        error('Position-bin count %d != cfg.n_position_bins=%d', ...
              n_pos_bins, cfg.n_position_bins);
    end
    % Precondition: freq_z must be finite. v5_residualise / v5_pca_reduce /
    % v5_cca_fit do not guard against NaN; a single NaN entry would
    % propagate through the pipeline and silently produce non-finite CCA
    % loadings. The freq_z field is the gap-filled+z-scored version and
    % should be NaN-free in well-formed exports; assert it explicitly.
    if any(~isfinite(raw_spatial(:)))
        n_bad = sum(~isfinite(raw_spatial(:)));
        error('load_animal:nonFiniteData', ...
            ['analysis_spatial.firing.cued.freq_z contains %d non-finite ', ...
             'entries — pipeline assumes clean (gap-filled, z-scored) input. ', ...
             'Did the export step run with the wrong field?'], n_bad);
    end

    % Per-unit region via units.idx + units.regions_label (spec §6.1 deviation).
    if ~isfield(units, 'idx') || ~isfield(units, 'regions_label')
        error('units.idx / units.regions_label missing');
    end
    region_labels = cellstr(units.regions_label);
    idx_logical = logical(units.idx);
    if size(idx_logical, 2) ~= n_units
        error('units.idx column count %d != n_units %d', ...
              size(idx_logical, 2), n_units);
    end

    % Per-unit FS keep mask.
    keep_mask = true(n_units, 1);
    if isfield(units, 'idx_fs')
        is_fs = logical(units.idx_fs(:));
        for r = 1:numel(cfg.exclude_fs_areas)
            area_name = cfg.exclude_fs_areas{r};
            row_mask = strcmp(region_labels, area_name);
            if any(row_mask)
                in_area = idx_logical(row_mask, :).' & is_fs;
                keep_mask(in_area) = false;
            end
        end
    end

    % Build per-area tensor (n_trials × n_pos_bins × n_units_in_area).
    % raw_spatial is (n_units × n_trials × n_pos_bins) → permute to
    % (n_trials × n_pos_bins × n_units) then slice on the unit axis.
    activity_tnu = permute(raw_spatial, [2 3 1]);   % (n_trials × n_pos × n_units)

    area_tensors = struct();
    area_kept_counts = struct();
    for r = 1:numel(region_labels)
        area_name = region_labels{r};
        if isempty(area_name), continue; end
        u_logical = idx_logical(r, :).' & keep_mask;
        if sum(u_logical) < cfg.min_units_per_region, continue; end
        area_tensors.(area_name) = activity_tnu(:, :, u_logical);
        area_kept_counts.(area_name) = sum(u_logical);
    end

    % Per-trial-per-bin speed (cm / s), derived from occupancy time.
    % `analysis_behaviour.spatialOccupation.cued.time_gf` is the gap-
    % filled time per bin per trial in seconds — bin_width / time gives
    % the bin-average speed. We fall back to NaN if the field is absent
    % (e.g. legacy export) so the downstream behaviour-correlation step
    % just skips speed.
    if isfield(D, 'analysis_behaviour') && ...
       isfield(D.analysis_behaviour, 'spatialOccupation') && ...
       isfield(D.analysis_behaviour.spatialOccupation, 'cued') && ...
       isfield(D.analysis_behaviour.spatialOccupation.cued, 'time_gf')
        time_per_bin = D.analysis_behaviour.spatialOccupation.cued.time_gf;
        % Expected shape: (n_pos_bins x n_trials). Transpose to
        % (n_trials x n_pos_bins) for consistency with our tensors.
        if size(time_per_bin, 1) == n_pos_bins && size(time_per_bin, 2) == n_trials
            time_per_bin = time_per_bin.';
        elseif size(time_per_bin, 1) ~= n_trials || size(time_per_bin, 2) ~= n_pos_bins
            warning('load_animal:speedShape', ...
                'time_gf has shape %s; expected %d-by-%d or its transpose. Speed unavailable.', ...
                mat2str(size(time_per_bin)), n_trials, n_pos_bins);
            time_per_bin = nan(n_trials, n_pos_bins);
        end
        % Speed (cm/s) = bin_width / time. Guard against zero / negative.
        bin_width_cm = cfg.bin_size_cm;
        speed_cm_s = bin_width_cm ./ time_per_bin;
        speed_cm_s(time_per_bin <= 0) = NaN;
    else
        speed_cm_s = nan(n_trials, n_pos_bins);
    end

    animal = struct( ...
        'n_trials', n_trials, ...
        'n_pos_bins', n_pos_bins, ...
        'area_tensors', area_tensors, ...
        'area_kept_counts', area_kept_counts, ...
        'speed_cm_s', speed_cm_s);
end


function epochs = epoch_indices_for_animal(n_trials, lp, n_per_epoch)
% Return struct array with .name and .trials_idx for the three epochs.
% Naive = trials 1..n_per_epoch.
% Pre-LP = trials lp-n_per_epoch+1..lp.
% Post-LP = trials lp+1..lp+n_per_epoch.
    epochs = struct('name', {}, 'trials_idx', {}, 'iep', {});
    if n_trials < 3 * n_per_epoch || lp < n_per_epoch || lp + n_per_epoch > n_trials
        return;
    end
    epochs(1) = struct('name', 'naive',   'trials_idx', 1:n_per_epoch, 'iep', 1);
    epochs(2) = struct('name', 'pre_lp',  'trials_idx', (lp - n_per_epoch + 1):lp, 'iep', 2);
    epochs(3) = struct('name', 'post_lp', 'trials_idx', (lp + 1):(lp + n_per_epoch), 'iep', 3);
end


function k = choose_k(cfg, n_units_x, n_units_y, n_samples, S_x_full, S_y_full, ref_trials)
% Choose k per spec §3.2.
%   k = min(units_smaller_area, n_samples/k_samples_rule, k_cap)
% In 'variance' mode, k is set to capture cfg.k_variance_target of variance
% per side on the reference trials, capped by the same upper limit.
    k_units = min(n_units_x, n_units_y);
    k_samples = floor(n_samples / cfg.k_samples_rule);
    k_default = min([k_units, k_samples, cfg.k_cap]);

    switch lower(cfg.k_mode)
        case 'samples'
            k = max(1, k_default);
        case 'variance'
            % Fit PCA on the reference (typically smallest) epoch and pick the
            % smaller of the X- and Y-side k that hits the variance target.
            S_x_ref = S_x_full(ref_trials, :, :);
            S_y_ref = S_y_full(ref_trials, :, :);
            S_x_res = v5_residualise(S_x_ref, ones(size(S_x_ref, 1), 1));
            S_y_res = v5_residualise(S_y_ref, ones(size(S_y_ref, 1), 1));
            % PCA at the full available k, then pick the cumulative-variance threshold.
            [~, state_x] = v5_pca_reduce(S_x_res, min(k_default, size(S_x_ref, 3)));
            [~, state_y] = v5_pca_reduce(S_y_res, min(k_default, size(S_y_ref, 3)));
            cum_x = cumsum(state_x.explained_variance_ratio);
            cum_y = cumsum(state_y.explained_variance_ratio);
            k_x = find(cum_x >= cfg.k_variance_target, 1, 'first');
            k_y = find(cum_y >= cfg.k_variance_target, 1, 'first');
            if isempty(k_x), k_x = numel(cum_x); end
            if isempty(k_y), k_y = numel(cum_y); end
            k = min([k_x, k_y, k_default]);
            k = max(1, k);
        otherwise
            error('Unknown k_mode: %s', cfg.k_mode);
    end
end


function res = fit_one(S_x_3d, S_y_3d, k, cfg, speed_cm_s, trial_idx_in_session)
% One CCA fit for one (animal × area-pair × epoch). Composes the v5
% primitives into the full Stage 2-8 sequence on this slice, then adds
% the auxiliary analyses (MI, per-dim significance, cell-space weights,
% weight spatial entropy, projection-vs-behaviour correlations).
    [n_trials, n_pos_bins, n_units_x] = size(S_x_3d);
    n_units_y = size(S_y_3d, 3);

    % Stage 2 — residualise per (unit, position-bin) PSTH.
    % condition_labels = ones means subtract grand mean across all trials
    % per (bin, unit), i.e. the position-tuning curve.
    cl = ones(n_trials, 1);
    S_x_res = v5_residualise(S_x_3d, cl);
    S_y_res = v5_residualise(S_y_3d, cl);

    % Stage 3 — PCA per area.
    [P_x_3d, state_x] = v5_pca_reduce(S_x_res, k);
    [P_y_3d, state_y] = v5_pca_reduce(S_y_res, k);

    % Flatten to (n_samples, k) for Stages 4-5.
    P_x_flat = reshape(permute(P_x_3d, [2 1 3]), n_trials * n_pos_bins, k);
    P_y_flat = reshape(permute(P_y_3d, [2 1 3]), n_trials * n_pos_bins, k);

    % Stage 4 — CV CCA.
    fit_res = v5_cca_fit_cv(P_x_flat, P_y_flat, cfg.cv_splits, true, cfg.cca_fold_seed);

    % Stage 5 — shuffle null. samples_per_trial = n_pos_bins so the
    % permutation preserves within-trial position structure.
    null = v5_shuffle_null(P_x_flat, P_y_flat, cfg.n_shuffles, n_pos_bins, ...
                            cfg.cv_splits, cfg.shuffle_seed, true, cfg.cca_fold_seed);

    % Stage 7 — lagged-refit IFI on the 3D PC tensors.
    lagged = v5_lagged_refit_ifi(P_x_3d, P_y_3d, ...
        cfg.max_lag_bins, cfg.central_window_bins, ...
        cfg.smooth_size, cfg.smooth_sigma);

    % Stage 6+8 — per-trial projection IFI. Project each trial through
    % the full-sample (A, B) loadings, then run projection_ifi on the
    % (n_trials × n_pos_bins) scalar variates.
    u_3d = v5_project(P_x_3d, fit_res.x_mean, fit_res.A(:, 1));   % (n_trials × n_pos × 1)
    v_3d = v5_project(P_y_3d, fit_res.y_mean, fit_res.B(:, 1));
    u_first = squeeze_to_trial_bin(u_3d, n_trials, n_pos_bins);
    v_first = squeeze_to_trial_bin(v_3d, n_trials, n_pos_bins);
    [ifi_proj_pt, lag_corr_pt] = v5_projection_ifi(u_first, v_first, ...
        cfg.proj_max_lag_bins, cfg.proj_min_paired_samples);

    % --- Auxiliary analyses (Gonzalez et al., adapted) ----------------

    % Cell-space weights via PCA composition. These give the per-cell
    % weight for each canonical dimension and are the basis on which
    % cross-epoch principal angles are computed in the driver.
    W_cells_A = v5_cell_weights(fit_res.A, state_x);   % (n_units_x × r)
    W_cells_B = v5_cell_weights(fit_res.B, state_y);   % (n_units_y × r)

    % Per-dimension significance via the full per-dim null distribution
    % already returned by v5_shuffle_null. Use cc_cv_observed = fit_res.cc_cv.
    sig_mask = v5_significant_dims(fit_res.cc_cv, null.cc_cv, cfg.sig_alpha);

    % MI = -0.5 * sum log(1 - cc^2). Report both versions: over all dims
    % and over the significant subset.
    mi_all = v5_mi_from_cc(fit_res.cc_cv);
    mi_sig = v5_mi_from_cc(fit_res.cc_cv, sig_mask);

    % Weight spatial entropy. Use the trial-averaged residualised PSTH
    % as each cell's tuning curve (n_units × n_pos_bins). This is what
    % the entropy helper bins by preferred position.
    tuning_x = squeeze(mean(S_x_res, 1)).';   % (n_pos_bins × n_units) → transpose
    tuning_y = squeeze(mean(S_y_res, 1)).';
    if size(tuning_x, 2) ~= n_pos_bins, tuning_x = tuning_x.'; end
    if size(tuning_y, 2) ~= n_pos_bins, tuning_y = tuning_y.'; end
    [ent_x, ent_x_null] = v5_weight_spatial_entropy( ...
        W_cells_A, tuning_x, cfg.n_spatial_bins_entropy, ...
        cfg.n_surrogates_entropy, cfg.aux_seed);
    [ent_y, ent_y_null] = v5_weight_spatial_entropy( ...
        W_cells_B, tuning_y, cfg.n_spatial_bins_entropy, ...
        cfg.n_surrogates_entropy, cfg.aux_seed);

    % Projection-vs-behaviour correlations for each significant dim.
    % Project all r canonical pairs; correlate u_d and v_d separately
    % against position, speed, trial_number. (Trial number = within-
    % epoch index; absolute session indices in trial_idx_in_session.)
    r_dims = numel(fit_res.cc_cv);
    behaviour = struct( ...
        'position', 1:n_pos_bins, ...
        'trial_number', trial_idx_in_session(:));
    if any(isfinite(speed_cm_s(:)))
        behaviour.speed = speed_cm_s;
    end
    proj_corr_u = cell(r_dims, 1);
    proj_corr_v = cell(r_dims, 1);
    for d = 1:r_dims
        u_d = v5_project(P_x_3d, fit_res.x_mean, fit_res.A(:, d));
        v_d = v5_project(P_y_3d, fit_res.y_mean, fit_res.B(:, d));
        u_d = squeeze_to_trial_bin(u_d, n_trials, n_pos_bins);
        v_d = squeeze_to_trial_bin(v_d, n_trials, n_pos_bins);
        proj_corr_u{d} = v5_projection_behaviour_corr( ...
            u_d, behaviour, cfg.n_surrogates_proj_corr, cfg.aux_seed + d);
        proj_corr_v{d} = v5_projection_behaviour_corr( ...
            v_d, behaviour, cfg.n_surrogates_proj_corr, cfg.aux_seed + d + r_dims);
    end

    % --- Bundle ------------------------------------------------------
    cc_excess = fit_res.cc_cv(1) - null.mean;
    res = struct( ...
        'A', fit_res.A, ...
        'B', fit_res.B, ...
        'x_mean', fit_res.x_mean, ...
        'y_mean', fit_res.y_mean, ...
        'cc_full', fit_res.cc, ...
        'cc_cv', fit_res.cc_cv, ...
        'cc_cv_per_fold', fit_res.cc_cv_per_fold, ...
        'cc_excess', cc_excess, ...
        'null_threshold', null.threshold, ...
        'null_mean', null.mean, ...
        'null_cc1', null.cc1_null, ...
        'null_cc_cv', null.cc_cv, ...
        'ifi_lagged', lagged.ifi, ...
        'ifi_lagged_curve_xy', lagged.cc_xy, ...
        'ifi_lagged_curve_yx', lagged.cc_yx, ...
        'ifi_lagged_curve_xy_smooth', lagged.cc_xy_smooth, ...
        'ifi_lagged_curve_yx_smooth', lagged.cc_yx_smooth, ...
        'ifi_proj_per_trial', ifi_proj_pt, ...
        'ifi_proj_lag_corr', lag_corr_pt, ...
        'k_used', k, ...
        'var_x_at_k', sum(state_x.explained_variance_ratio), ...
        'var_y_at_k', sum(state_y.explained_variance_ratio), ...
        'n_samples', n_trials * n_pos_bins, ...
        'samples_per_var', (n_trials * n_pos_bins) / (2 * k), ...
        'W_cells_A', W_cells_A, ...
        'W_cells_B', W_cells_B, ...
        'sig_dims', sig_mask, ...
        'mi_all', mi_all, ...
        'mi_sig', mi_sig, ...
        'weight_entropy_A', ent_x, ...
        'weight_entropy_A_null', ent_x_null, ...
        'weight_entropy_B', ent_y, ...
        'weight_entropy_B_null', ent_y_null, ...
        'proj_corr_u', {proj_corr_u}, ...
        'proj_corr_v', {proj_corr_v}, ...
        'principal_angles_A', [], ...
        'principal_angles_B', []);
end


function y = squeeze_to_trial_bin(P_one_dim, n_trials, n_pos_bins)
% SQUEEZE_TO_TRIAL_BIN  Reshape a (n_trials × n_pos × 1) tensor to (n_trials × n_pos).
%   Handles the awkward edge case where MATLAB's `squeeze` collapses singleton
%   trial or bin dimensions when they happen to equal 1.
    y = squeeze(P_one_dim);
    if isvector(y)
        if n_trials > 1
            y = reshape(y, n_trials, n_pos_bins);
        else
            y = y(:).';
        end
    end
end


function group_results = init_results_struct(cfg, n_animals, n_pairs)
    empty_res = struct( ...
        'A', [], 'B', [], 'x_mean', [], 'y_mean', [], ...
        'cc_full', [], 'cc_cv', [], 'cc_cv_per_fold', [], 'cc_excess', NaN, ...
        'null_threshold', NaN, 'null_mean', NaN, 'null_cc1', [], ...
        'null_cc_cv', [], ...
        'ifi_lagged', NaN, 'ifi_lagged_curve_xy', [], 'ifi_lagged_curve_yx', [], ...
        'ifi_lagged_curve_xy_smooth', [], 'ifi_lagged_curve_yx_smooth', [], ...
        'ifi_proj_per_trial', [], 'ifi_proj_lag_corr', [], ...
        'k_used', NaN, 'var_x_at_k', NaN, 'var_y_at_k', NaN, ...
        'n_samples', NaN, 'samples_per_var', NaN, ...
        'W_cells_A', [], 'W_cells_B', [], ...
        'sig_dims', [], ...
        'mi_all', NaN, 'mi_sig', NaN, ...
        'weight_entropy_A', [], 'weight_entropy_A_null', [], ...
        'weight_entropy_B', [], 'weight_entropy_B_null', [], ...
        'proj_corr_u', {{}}, 'proj_corr_v', {{}}, ...
        'principal_angles_A', [], 'principal_angles_B', []);
    group_results = struct();
    for ipair = 1:n_pairs
        group_results(ipair).pair_name = sprintf('%s-%s', ...
            cfg.area_pairs{ipair, 1}, cfg.area_pairs{ipair, 2});
        group_results(ipair).area_x = cfg.area_pairs{ipair, 1};
        group_results(ipair).area_y = cfg.area_pairs{ipair, 2};
        group_results(ipair).per_epoch = repmat(empty_res, n_animals, 3);
    end
end
