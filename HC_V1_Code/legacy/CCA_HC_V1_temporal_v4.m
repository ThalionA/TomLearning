%% CCA_HC_V1_temporal_v4.m
% Per-trial Temporal CCA + RZ-pre per-bin Temporal CCA.
%
% DESIGN (vs. v3)
% ---------------
% v3 fits canoncorr on a sliding TIME-BIN window WITHIN each trial and
% gives a time-resolved CC1/IFI trace. This script (v4) instead fits
% canoncorr ONCE per trial on all valid bins of the trial, and adds an
% "RZ-pre per-bin" analysis where trials become n at each aligned bin
% in the last 500 ms BEFORE reward-zone entry.
%
% PIPELINE PER ANIMAL
%   1. Re-bin spikes / behaviour from native (~1 ms) to target_bin_ms.
%   2. Per region, z-score units across base_valid bins, then PCA-reduce
%      (>=90% var, floor 3 PCs, capped at max_k_per_region). [v3-identical]
%   3. For each trial in the epoch set:
%        (A) Per-trial CC: take the trial's base_valid bins, fit canoncorr
%            once -> r(1). Shuffle = random permutation of Y rows.
%        (B) Per-trial IFI: project full trial onto trial-level (A,B);
%            split into contiguous valid blocks (>= min_block_bins); per
%            block compute lag corr / IFI; weighted-mean by block length.
%            Shuffle = per-block row-shift on v.
%        (C) RZ-pre per-trial IFI: project the last 500 ms before RZ entry
%            onto trial-level (A,B); compute lag corr / IFI on that
%            segment. (Per-trial output; epoch averages handled below.)
%        (D) RZ-pre alignment: stash the (n_rz_bins x k) X/Y segments
%            for trials with a clean pre-RZ window. Per-bin CCA across
%            trials happens after the trial loop, once per epoch.
%   4. Aggregate per pair x epoch:
%        - trial_cc1 / trial_ifi: epoch-mean of per-trial scalars
%        - rz_cc1: at each aligned bin, fit canoncorr across the epoch's
%                  trials (n = trials_in_epoch) -> CC trace of length n_rz_bins
%        - rz_ifi: epoch-mean of per-trial pre-RZ IFI
%
% NO sliding window inside trials. NO spatial post-hoc binning. Those
% remain in v3 if you want them.

%% 1. SETUP & PARAMETERS
clear; clc; close all;

base_dir       = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir       = fullfile(base_dir, data_subfolder);
file_pattern   = 'TF*_export.mat';
learning_file  = 'animal_behaviour.mat';

current_date = datestr(now, 'yyyy_mm_dd');
save_path    = fullfile(data_dir, sprintf('Temporal_CCA_v4_%s.mat', current_date));

% --- Binning ---
target_bin_ms          = 25;        % matches v3
min_speed_cms          = 2;         % matches v3

% --- PCA / CCA ---
pca_variance_threshold = 90;        % matches spatial v2 / v3
min_units_per_region   = 5;
max_k_per_region       = 15;         % matches v3

% --- Per-trial IFI (lag analysis) ---
% lags = -3:3 at 25 ms binning -> +-75 ms, same as v3.
lags                   = -3:3;
min_block_bins         = 12;        % per design discussion

% --- Shuffling ---
n_shuffles             = 20;        % more than v3 since this script is faster
                                    % (no sliding window per trial)

% --- RZ pre-window ---
rz_entry_cm            = 200;       % matches v3
rz_window_ms           = 500;
n_rz_bins              = round(rz_window_ms / target_bin_ms);   % 20 at 25 ms

% --- RZ-per-bin small-sample tension ---
% At each aligned bin we run canoncorr with n = trials_in_epoch (<=10).
% Two knobs control whether the fit is even attempted:
%
%   max_k_rz_per_region : cap each region's PC dimensionality SPECIFICALLY
%       for the cross-trial RZ-per-bin path. The per-trial CC and per-trial
%       IFI still use the full max_k_per_region. With cap=3 and n_tr=10,
%       the outer guard n > k1+k2+min_extra becomes 10 > 6+2 = 8 -> passes.
%
%   rz_min_extra_samples : how many extra samples beyond k1+k2 to require.
%       Smaller = more permissive. We rely on real-minus-shuffle to absorb
%       the small-sample bias.
%
% NOTE: the previous default (cap = max_k_per_region = 15, min_extra = 2)
% caused the outer guard to fail for nearly every (pair, epoch); see
% HC_V1_Code/PLAN_v4_fixes.md.
max_k_rz_per_region    = 3;
rz_min_extra_samples   = 2;

% --- Epochs ---
n_trials_epoch         = 10;
epoch_names            = {'early', 'pre', 'post'};

% --- Region pairs ---
area_pairs_to_analyze = {'CA1','V1';  'CA1','DG'; 'CA1','CA3'; 'CA1','RSC'; ...
                         'CA1','SUB'; 'V1','RSC'; 'RSC','SUB'; 'CA3','DG'};
n_pairs = size(area_pairs_to_analyze, 1);

% Silence the canoncorr full-rank warning; v4_per_trial_cca does an
% explicit rank check.
warning('off', 'stats:canoncorr:NotFullRank');

% Helpers live in this folder; make sure they're on the path.
addpath(fileparts(mfilename('fullpath')));

%% 2. LOAD & YOKE LEARNING POINTS
file_list = dir(fullfile(data_dir, file_pattern));
n_animals = length(file_list);

lp_path = fullfile(data_dir, learning_file);
if exist(lp_path, 'file')
    dat_lp = load(lp_path);
    learning_points = dat_lp.period_experienced(:, 1);
    if isfield(dat_lp, 'animal_id')
        [~, sorting_idx] = sort(dat_lp.animal_id);
        learning_points  = learning_points(sorting_idx);
    end
else
    learning_points = nan(n_animals, 1);
end
if length(learning_points) < n_animals
    learning_points(end+1:n_animals) = nan;
end

is_learner   = ~isnan(learning_points);
mean_lp      = round(mean(learning_points(is_learner), 'omitnan'));
analysis_lp  = learning_points;
analysis_lp(~is_learner) = mean_lp;

fprintf('Identified %d Learners and %d Non-Learners. Yoked LP = %d\n', ...
    sum(is_learner), sum(~is_learner), mean_lp);

%% 3. INITIALIZE RESULTS
group_results = struct();
for ipair = 1:n_pairs
    group_results(ipair).pair_name = sprintf('%s-%s', ...
        area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});

    for e = 1:numel(epoch_names)
        ep = epoch_names{e};

        % Per-trial CC and IFI scalars (one per animal, mean over epoch trials).
        group_results(ipair).(['trial_cc1_'    ep])    = nan(n_animals, 1);
        group_results(ipair).(['trial_cc1_sh_' ep])    = nan(n_animals, 1);
        group_results(ipair).(['trial_ifi_'    ep])    = nan(n_animals, 1);
        group_results(ipair).(['trial_ifi_sh_' ep])    = nan(n_animals, 1);

        % Per-trial vectors (cell array indexed by animal). For pairing
        % against spatial v2's trial-level outputs.
        group_results(ipair).(['trial_cc1_pertrial_'    ep]) = cell(n_animals, 1);
        group_results(ipair).(['trial_cc1_sh_pertrial_' ep]) = cell(n_animals, 1);
        group_results(ipair).(['trial_ifi_pertrial_'    ep]) = cell(n_animals, 1);
        group_results(ipair).(['trial_ifi_sh_pertrial_' ep]) = cell(n_animals, 1);
        group_results(ipair).(['trial_id_pertrial_'     ep]) = cell(n_animals, 1);

        % RZ-pre per-bin CC trace (one per animal x aligned bin).
        group_results(ipair).(['rz_cc1_'    ep])       = nan(n_animals, n_rz_bins);
        group_results(ipair).(['rz_cc1_sh_' ep])       = nan(n_animals, n_rz_bins);

        % RZ-pre per-trial IFI scalar (epoch-mean of per-trial values).
        group_results(ipair).(['rz_ifi_'    ep])       = nan(n_animals, 1);
        group_results(ipair).(['rz_ifi_sh_' ep])       = nan(n_animals, 1);
        group_results(ipair).(['rz_ifi_pertrial_'    ep]) = cell(n_animals, 1);
        group_results(ipair).(['rz_ifi_sh_pertrial_' ep]) = cell(n_animals, 1);
    end
end

%% 4. MAIN PIPELINE (serial)
for ianimal = 1:n_animals
    filename = file_list(ianimal).name;
    fullpath = fullfile(data_dir, filename);
    fprintf('\n=== Animal %d/%d: %s ===\n', ianimal, n_animals, filename);
    t_animal = tic;

    D = load(fullpath);
    if ~isfield(D, 'units') || ~isfield(D, 'binned_spikes')
        fprintf('  -> Missing fields. Skipping.\n'); continue;
    end
    units = D.units;

    % --- Per-unit region array (use units.idx + units.regions_label) ---
    % Mirrors spatial_v2's tagging. units.region is unreliable: in 11 of
    % 12 RSC-containing animals it tags zero RSC units even though
    % units.idx marks tens-to-hundreds. See HC_V1_Code/PLAN_v4_fixes.md
    % and the probe_units_and_v4.m output. v4_unit_regions returns the
    % per-unit region label and the canonical region name list.
    try
        [unit_regions, animal_areas, ru_info] = v4_unit_regions(units);
    catch ME
        fprintf('  -> v4_unit_regions failed: %s. Skipping.\n', ME.message);
        continue;
    end
    if ru_info.n_multi > 0
        fprintf('  -> %d/%d units in >1 idx row; using last-match assignment.\n', ...
            ru_info.n_multi, ru_info.n_units);
    end
    if ru_info.n_unassigned > 0
        fprintf('  -> %d/%d units have no idx-row membership (will be skipped).\n', ...
            ru_info.n_unassigned, ru_info.n_units);
    end

    % --- FS unit exclusion (v3-identical) ---
    keep_mask = true(length(units.unit_id), 1);
    if isfield(units, 'idx_fs')
        target_fs_areas = {'V1','RSC','CA1','CA3'};
        is_fs = logical(units.idx_fs(:));
        if length(is_fs) == length(unit_regions)
            for r = 1:length(target_fs_areas)
                fs_in_area = strcmp(unit_regions, target_fs_areas{r}) & is_fs;
                keep_mask(fs_in_area) = false;
            end
        end
    end

    % --- Re-binning (v3-identical, per-unit to keep memory bounded) ---
    native_bin_s = double(D.params_main.bin_size);
    rebin_factor = round((target_bin_ms / 1000) / native_bin_s);
    if rebin_factor < 1, rebin_factor = 1; end
    fprintf('  native bin = %.4f s, rebin factor = %d (target %d ms)\n', ...
        native_bin_s, rebin_factor, target_bin_ms);

    bs_native = D.binned_spikes;
    bs_layout_units_first = (size(bs_native, 1) == length(units.unit_id));
    if bs_layout_units_first
        nUnits = size(bs_native, 1);
        T      = size(bs_native, 2);
    else
        nUnits = size(bs_native, 2);
        T      = size(bs_native, 1);
    end
    T_use = floor(T / rebin_factor) * rebin_factor;
    T_re  = T_use / rebin_factor;

    spikes_re = zeros(T_re, nUnits);
    for u = 1:nUnits
        if bs_layout_units_first
            row = double(bs_native(u, 1:T_use));
        else
            row = double(bs_native(1:T_use, u))';
        end
        spikes_re(:, u) = sum(reshape(row, rebin_factor, T_re), 1)';
    end
    clear bs_native;

    mask_cued = logical(D.analysis_behaviour.masks.tunnel_cued(:));
    if isfield(D.data_behaviour, 'velocity_binned_gf')
        vel_gf = double(D.data_behaviour.velocity_binned_gf(:));
    else
        vel_gf = double(D.data_behaviour.velocity_gf(:));
    end
    tr_cued = double(D.data_behaviour.trial_binned_cued(:));
    if isfield(D.data_behaviour, 'pos_binned_gf')
        pos_cm = double(D.data_behaviour.pos_binned_gf(:));
    elseif isfield(D.data_behaviour, 'pos_binned')
        pos_cm = double(D.data_behaviour.pos_binned(:));
    else
        pos_cm = nan(size(tr_cued));
    end
    clear D;

    mask_cued = mask_cued(1:T_use);
    vel_gf    = vel_gf(1:T_use);
    tr_cued   = tr_cued(1:T_use);
    pos_cm    = pos_cm(1:T_use);

    mc_re   = all(reshape(mask_cued, rebin_factor, T_re), 1)';
    vel_re  = mean(reshape(vel_gf,   rebin_factor, T_re), 1)';
    tc_re   = mode(reshape(tr_cued,  rebin_factor, T_re), 1)';
    pos_re  = mean(reshape(pos_cm,   rebin_factor, T_re), 1)';

    base_valid = mc_re & (vel_re >= min_speed_cms);

    % --- Z-score using valid samples ---
    mu = mean(spikes_re(base_valid, :), 1, 'omitnan');
    sd = std( spikes_re(base_valid, :), 0, 1, 'omitnan'); sd(sd==0) = 1;
    spikes_z = (spikes_re - repmat(mu, T_re, 1)) ./ repmat(sd, T_re, 1);

    % --- Per-region PCA (v3-identical) ---
    AreaActivity = struct();
    for ia = 1:length(animal_areas)
        area = char(animal_areas{ia});
        if isempty(area), continue; end

        u_logical = strcmp(unit_regions, area) & keep_mask(:);
        if sum(u_logical) < min_units_per_region, continue; end

        X = spikes_z(:, u_logical);
        X_valid = X(base_valid, :);
        if size(X_valid, 1) < 50, continue; end
        if any(~isfinite(X_valid(:))), continue; end

        [coeff, ~, ~, ~, explained] = pca(X_valid);
        cum_var = cumsum(explained);
        k_pca = find(cum_var >= pca_variance_threshold, 1);
        if isempty(k_pca), k_pca = size(coeff, 2); end
        k = min([max(3, k_pca), max_k_per_region, size(coeff, 2)]);

        X_centered = X - repmat(mean(X_valid, 1), T_re, 1);
        AreaActivity.(area).data = X_centered * coeff(:, 1:k);
        AreaActivity.(area).k    = k;
        fprintf('    [PCA] %s: %d units -> k=%d (k_pca@90%%=%d, var captured = %.1f%%)\n', ...
            area, sum(u_logical), k, k_pca, cum_var(k));
    end

    % --- Trial set / epochs ---
    num_trials = max(tc_re);
    lp = analysis_lp(ianimal);
    if isnan(lp) || lp <= n_trials_epoch || (lp + n_trials_epoch - 1) > num_trials
        fprintf('  LP=%g unusable for this animal (num_trials=%d). Skipping.\n', lp, num_trials);
        continue;
    end

    epoch_trials = struct( ...
        'early', 1:n_trials_epoch, ...
        'pre',   (lp - n_trials_epoch):(lp - 1), ...
        'post',  lp:(lp + n_trials_epoch - 1));
    epoch_trial_set = unique([epoch_trials.early, epoch_trials.pre, epoch_trials.post]);

    % Per-pair caches across trials (one row per trial actually run).
    pair_cache = repmat(struct('trial_id', [], ...
                               'cc1', [],  'cc1_sh', [], ...
                               'ifi', [],  'ifi_sh', [], ...
                               'rz_ifi', [], 'rz_ifi_sh', [], ...
                               'rz_X', {{}}, 'rz_Y', {{}}, ...
                               'rz_trial_id', []), n_pairs, 1);

    for it = 1:numel(epoch_trial_set)
        tr_id     = epoch_trial_set(it);
        trial_in  = (tc_re == tr_id);
        trial_idx = find(trial_in);
        n_tr_bins = numel(trial_idx);
        if n_tr_bins < min_block_bins
            continue;
        end
        is_valid_tr = base_valid(trial_idx);
        pos_tr      = pos_re(trial_idx);

        for ipair = 1:n_pairs
            a1 = area_pairs_to_analyze{ipair, 1};
            a2 = area_pairs_to_analyze{ipair, 2};
            if ~isfield(AreaActivity, a1) || ~isfield(AreaActivity, a2), continue; end

            X_full = AreaActivity.(a1).data;
            Y_full = AreaActivity.(a2).data;

            X_tr = X_full(trial_idx, :);
            Y_tr = Y_full(trial_idx, :);

            % --- (A) Per-trial CC ---
            X_valid = X_tr(is_valid_tr, :);
            Y_valid = Y_tr(is_valid_tr, :);
            [r_real, A_load, B_load, info_cc] = v4_per_trial_cca(X_valid, Y_valid);

            cc_real = NaN; cc_sh = NaN;
            ifi_real = NaN; ifi_sh = NaN;
            rz_ifi_real = NaN; rz_ifi_sh = NaN;
            rz_X_seg = []; rz_Y_seg = [];

            if info_cc.ok
                cc_real = r_real;

                % --- Shuffle CC: permute Y_valid rows ---
                r_sh_iter = nan(n_shuffles, 1);
                n_v = size(Y_valid, 1);
                for ish = 1:n_shuffles
                    perm = randperm(n_v);
                    [r_s, ~, ~, info_s] = v4_per_trial_cca(X_valid, Y_valid(perm, :));
                    if info_s.ok, r_sh_iter(ish) = r_s; end
                end
                cc_sh = mean(r_sh_iter, 'omitnan');

                % --- (B) Per-trial IFI (per-block weighted) ---
                [ifi_real, info_ifi] = v4_per_trial_ifi(X_tr, Y_tr, is_valid_tr, ...
                                                        A_load, B_load, lags, min_block_bins);

                % Shuffle IFI: per-block row-shift on v (preserves block
                % autocorrelation, breaks alignment). Build full v_shuf,
                % then compute weighted IFI the same way.
                if ~isempty(info_ifi.blocks)
                    ifi_sh_iter = nan(n_shuffles, 1);
                    for ish = 1:n_shuffles
                        ifi_sh_iter(ish) = local_block_shifted_ifi( ...
                            X_tr, Y_tr, is_valid_tr, A_load, B_load, ...
                            lags, info_ifi.blocks, info_ifi.block_used);
                    end
                    ifi_sh = mean(ifi_sh_iter, 'omitnan');
                end

                % --- (D) RZ pre-segment for cross-trial fit ---
                [rz_X_seg, rz_Y_seg, ok_rz] = v4_rz_align_pre( ...
                    X_tr, Y_tr, pos_tr, is_valid_tr, rz_entry_cm, n_rz_bins);

                % --- (C) RZ pre IFI on this trial's segment ---
                if ok_rz
                    Xc = rz_X_seg - repmat(mean(rz_X_seg, 1), n_rz_bins, 1);
                    Yc = rz_Y_seg - repmat(mean(rz_Y_seg, 1), n_rz_bins, 1);
                    u_rz = Xc * A_load;
                    v_rz = Yc * B_load;
                    rz_ifi_real = v4_ifi_from_lags(v4_lag_corr(u_rz, v_rz, lags), lags);

                    % Shuffle RZ IFI: shift v_rz relative to u_rz, then
                    % trim both so the surviving overlap is NaN-free.
                    sh_lo = max(abs(lags)) + 1;
                    sh_hi = n_rz_bins - max(abs(lags)) - 2;
                    if sh_hi >= sh_lo
                        ifi_rz_sh_iter = nan(n_shuffles, 1);
                        for ish = 1:n_shuffles
                            k_sh = randi([sh_lo, sh_hi]);
                            if rand() < 0.5, k_sh = -k_sh; end
                            [u_a, v_a] = local_shift_align(u_rz, v_rz, k_sh);
                            ifi_rz_sh_iter(ish) = v4_ifi_from_lags( ...
                                v4_lag_corr(u_a, v_a, lags), lags);
                        end
                        rz_ifi_sh = mean(ifi_rz_sh_iter, 'omitnan');
                    end
                end
            else
                if it == 1
                    fprintf('    [skip] %s-%s tr=%d: %s (n=%d, k=%d+%d)\n', ...
                        a1, a2, tr_id, info_cc.reason, info_cc.n_samples, ...
                        info_cc.k1, info_cc.k2);
                end
            end

            % --- Stash ---
            pair_cache(ipair).trial_id(end+1, 1)  = tr_id;
            pair_cache(ipair).cc1(end+1, 1)       = cc_real;
            pair_cache(ipair).cc1_sh(end+1, 1)    = cc_sh;
            pair_cache(ipair).ifi(end+1, 1)       = ifi_real;
            pair_cache(ipair).ifi_sh(end+1, 1)    = ifi_sh;
            pair_cache(ipair).rz_ifi(end+1, 1)    = rz_ifi_real;
            pair_cache(ipair).rz_ifi_sh(end+1, 1) = rz_ifi_sh;
            if ~isempty(rz_X_seg)
                % Cap PCs for the cross-trial RZ-per-bin fit only. The
                % per-trial pre-RZ IFI computation above used the full
                % rz_X_seg / rz_Y_seg already.
                k_cap_x = min(size(rz_X_seg, 2), max_k_rz_per_region);
                k_cap_y = min(size(rz_Y_seg, 2), max_k_rz_per_region);
                pair_cache(ipair).rz_X{end+1}        = rz_X_seg(:, 1:k_cap_x);
                pair_cache(ipair).rz_Y{end+1}        = rz_Y_seg(:, 1:k_cap_y);
                pair_cache(ipair).rz_trial_id(end+1, 1) = tr_id;
            end
        end

        if mod(it, 5) == 0 || it == numel(epoch_trial_set)
            fprintf('    trial %d / %d done (t=%.1f s)\n', ...
                it, numel(epoch_trial_set), toc(t_animal));
        end
    end

    % --- Reduce per-pair x epoch ---
    for ipair = 1:n_pairs
        if isempty(pair_cache(ipair).trial_id), continue; end
        tids = pair_cache(ipair).trial_id;

        for e = 1:numel(epoch_names)
            ep = epoch_names{e};
            sel = ismember(tids, epoch_trials.(ep));
            if ~any(sel), continue; end

            group_results(ipair).(['trial_cc1_'    ep])(ianimal) = mean(pair_cache(ipair).cc1(sel),     'omitnan');
            group_results(ipair).(['trial_cc1_sh_' ep])(ianimal) = mean(pair_cache(ipair).cc1_sh(sel),  'omitnan');
            group_results(ipair).(['trial_ifi_'    ep])(ianimal) = mean(pair_cache(ipair).ifi(sel),     'omitnan');
            group_results(ipair).(['trial_ifi_sh_' ep])(ianimal) = mean(pair_cache(ipair).ifi_sh(sel),  'omitnan');
            group_results(ipair).(['rz_ifi_'    ep])(ianimal)    = mean(pair_cache(ipair).rz_ifi(sel),    'omitnan');
            group_results(ipair).(['rz_ifi_sh_' ep])(ianimal)    = mean(pair_cache(ipair).rz_ifi_sh(sel), 'omitnan');

            % Per-trial vectors for paired comparisons against spatial v2.
            group_results(ipair).(['trial_cc1_pertrial_'    ep]){ianimal} = pair_cache(ipair).cc1(sel);
            group_results(ipair).(['trial_cc1_sh_pertrial_' ep]){ianimal} = pair_cache(ipair).cc1_sh(sel);
            group_results(ipair).(['trial_ifi_pertrial_'    ep]){ianimal} = pair_cache(ipair).ifi(sel);
            group_results(ipair).(['trial_ifi_sh_pertrial_' ep]){ianimal} = pair_cache(ipair).ifi_sh(sel);
            group_results(ipair).(['trial_id_pertrial_'     ep]){ianimal} = pair_cache(ipair).trial_id(sel);
            group_results(ipair).(['rz_ifi_pertrial_'    ep]){ianimal}    = pair_cache(ipair).rz_ifi(sel);
            group_results(ipair).(['rz_ifi_sh_pertrial_' ep]){ianimal}    = pair_cache(ipair).rz_ifi_sh(sel);

            % --- RZ per-bin CC across this epoch's trials ---
            sel_rz = ismember(pair_cache(ipair).rz_trial_id, epoch_trials.(ep));
            if any(sel_rz)
                X_stack = cat(3, pair_cache(ipair).rz_X{sel_rz});   % [n_rz_bins x k1 x n_tr]
                Y_stack = cat(3, pair_cache(ipair).rz_Y{sel_rz});   % [n_rz_bins x k2 x n_tr]
                n_tr_rz = size(X_stack, 3);
                k1 = size(X_stack, 2);
                k2 = size(Y_stack, 2);
                if n_tr_rz > (k1 + k2 + rz_min_extra_samples)
                    [cc_trace, cc_trace_sh] = local_rz_per_bin_cc( ...
                        X_stack, Y_stack, n_shuffles, rz_min_extra_samples);
                    group_results(ipair).(['rz_cc1_'    ep])(ianimal, :) = cc_trace;
                    group_results(ipair).(['rz_cc1_sh_' ep])(ianimal, :) = cc_trace_sh;
                else
                    fprintf('    [rz-skip] %s ep=%s: n_tr=%d, k1+k2=%d (need n>k1+k2+%d)\n', ...
                        group_results(ipair).pair_name, ep, n_tr_rz, k1+k2, rz_min_extra_samples);
                end
            else
                fprintf('    [rz-empty] %s ep=%s: 0 RZ-stash trials in this epoch\n', ...
                    group_results(ipair).pair_name, ep);
            end
        end
    end

    % --- Save partial progress ---
    save(save_path, 'group_results', 'is_learner', 'analysis_lp', ...
                    'target_bin_ms', 'lags', 'min_block_bins', ...
                    'rz_entry_cm', 'n_rz_bins', 'n_shuffles', ...
                    'max_k_per_region', '-v7.3');
    fprintf('  saved partial results -> %s  (animal took %.1f s)\n', save_path, toc(t_animal));
end

fprintf('\nFinished! Results saved to: %s\n', save_path);

%% LOCAL HELPERS

function ifi = local_block_shifted_ifi(X_tr, Y_tr, is_valid_tr, A, B, lags, blocks, block_used_real)
% Per-block row-shift on the canonical projections for the IFI null.
% Reuses the same blocks identified for the real IFI; if block was not
% used in the real IFI it isn't used here either.
    n = size(X_tr, 1);
    Xc = X_tr - repmat(mean(X_tr, 1), n, 1);
    Yc = Y_tr - repmat(mean(Y_tr, 1), n, 1);
    u_full = Xc * A;
    v_full = Yc * B;

    n_blk = size(blocks, 1);
    block_lens = blocks(:, 2)' - blocks(:, 1)' + 1;
    block_ifis_sh = nan(1, n_blk);

    min_for_lag = max(abs(lags)) + 6;

    for ib = 1:n_blk
        if ~block_used_real(ib), continue; end
        s = blocks(ib, 1); e = blocks(ib, 2); L = e - s + 1;
        if L < min_for_lag, continue; end

        u_b = u_full(s:e);
        v_b = v_full(s:e);

        sh_lo = max(abs(lags)) + 1;
        sh_hi = L - max(abs(lags)) - 2;
        if sh_hi < sh_lo, continue; end
        k_sh = randi([sh_lo, sh_hi]);
        if rand() < 0.5, k_sh = -k_sh; end
        [u_a, v_a] = local_shift_align(u_b, v_b, k_sh);
        rl_sh = v4_lag_corr(u_a, v_a, lags);
        block_ifis_sh(ib) = v4_ifi_from_lags(rl_sh, lags);
    end

    used = isfinite(block_ifis_sh);
    if ~any(used)
        ifi = NaN; return;
    end
    w = block_lens(used);
    ifi = sum(w .* block_ifis_sh(used)) / sum(w);
end

function [u_aligned, v_aligned] = local_shift_align(u, v, k)
% Non-circular shift of v relative to u by k rows, then trim BOTH so the
% surviving overlap is the same length and contains no padding.
%   k > 0: u "leads" by k -> u_aligned = u(k+1:end), v_aligned = v(1:end-k)
%   k < 0: v "leads" by |k| -> symmetric
% Both outputs have length n - |k|, both are NaN-free if inputs were.
    u = u(:); v = v(:);
    n = numel(u);
    if k >= 0
        u_aligned = u((k+1):n);
        v_aligned = v(1:(n-k));
    else
        u_aligned = u(1:(n+k));
        v_aligned = v((-k+1):n);
    end
end

function [cc_trace, cc_trace_sh] = local_rz_per_bin_cc(X_stack, Y_stack, n_shuffles, min_extra)
% At each aligned bin b in 1..n_rz_bins, fit canoncorr across the trial
% axis (n = n_tr_rz). Return CC trace and shuffle CC trace.
%
% Uses v4_relaxed_canoncorr (via local_relaxed_canoncorr) to allow a
% smaller n>k1+k2 margin than the default v4_per_trial_cca check, since
% we deliberately operate in a small-sample regime here and rely on
% shuffle correction.
    n_rz = size(X_stack, 1);
    n_tr = size(X_stack, 3);
    cc_trace    = nan(1, n_rz);
    cc_trace_sh = nan(1, n_rz);

    for b = 1:n_rz
        Xb = squeeze(X_stack(b, :, :))';   % [n_tr x k1]
        Yb = squeeze(Y_stack(b, :, :))';   % [n_tr x k2]

        [r, ok] = local_relaxed_canoncorr(Xb, Yb, min_extra);
        if ok, cc_trace(b) = r; end

        rs_iter = nan(n_shuffles, 1);
        for ish = 1:n_shuffles
            perm = randperm(n_tr);
            [r_s, ok_s] = local_relaxed_canoncorr(Xb, Yb(perm, :), min_extra);
            if ok_s, rs_iter(ish) = r_s; end
        end
        cc_trace_sh(b) = mean(rs_iter, 'omitnan');
    end
end

function [r, ok] = local_relaxed_canoncorr(X, Y, min_extra)
% Thin wrapper around v4_relaxed_canoncorr (standalone helper) so the
% RZ-per-bin path can be unit-tested directly.
    [r, ok, ~] = v4_relaxed_canoncorr(X, Y, min_extra);
end
