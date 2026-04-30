%% CCA_HC_V1_temporal_v3.m
% Per-bin Temporal CCA (within-trial sliding window)
%
% DESIGN
% ------
% v2 fit canoncorr on a sliding TRIAL window (concatenating bins from 7
% trials of unequal length). v3 fits canoncorr on a sliding TIME-BIN
% window WITHIN each trial. Per trial we get a time-resolved CC1 trace
% and a time-resolved IFI trace.
%
% PIPELINE PER ANIMAL
%   1. Re-bin spikes / behaviour from native (~1 ms) to 5 ms bins.
%   2. Per region, z-score units across valid bins, then PCA-reduce
%      (>=90% var, floor 3 PCs, capped at available PCs).
%   3. For each trial:
%        - extract that trial's PCA-reduced timecourse for each region
%        - slide a +-10 bin window (21 bins, 105 ms) along the trial
%        - at each center bin: fit canoncorr, compute lagged corr
%          (lags = -5..+5 bins = -25..+25 ms), compute IFI
%        - 5 lightweight non-circular shuffles per bin for a null
%   4. Aggregate per trial three ways:
%        (a) trial-mean: one CC1 / one IFI per trial
%        (b) spatial-binned: re-bin the trial's CC1/IFI trace into
%            n_pos_bins bins of position (cm)
%        (c) reward-zone-aligned: -500..+500 ms around RZ entry
%            (placeholder: RZ entry at first bin where pos >= rz_entry_cm)
%   5. Group results per epoch (early / pre-learning / post-learning),
%      yoked to LP. Each epoch averages within trial first, then across
%      the 10 trials of the epoch.
%
% Notes:
%   - Serial only (no parfor) so progress is observable.
%   - No circshift anywhere. Shuffle uses a non-circular row-shift on the
%     in-window Y, dropping the overhang.

%% 1. SETUP & PARAMETERS
clear; clc; close all;

base_dir       = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir       = fullfile(base_dir, data_subfolder);
file_pattern   = 'TF*_export.mat';
learning_file  = 'animal_behaviour.mat';

current_date = datestr(now, 'yyyy_mm_dd');
save_path    = fullfile(data_dir, sprintf('Temporal_CCA_v3_%s.mat', current_date));

% --- Binning ---
target_bin_ms          = 50;         % 50 ms bins (smooth velocity, dense PCA)
half_window_bins       = 5;          % +-5 bins -> 11-bin window (550 ms)
window_bins            = -half_window_bins:half_window_bins;
n_window               = numel(window_bins);

% --- Stride for sliding-window centers ---
% At 50 ms binning, every bin is a meaningful step. Set to 1 for full
% resolution; bump up if speed becomes an issue.
stride                 = 1;

% --- Lags for IFI ---
% With n_window = 11, lag L uses (n_window - |L|) paired samples. We
% need that to be safely > 5 for corr to be meaningful: keep |L| <= 3
% so the most extreme lag still has 8 pairs.
lags                   = -3:3;       % +-150 ms
max_lag_abs            = max(abs(lags));

% --- PCA / CCA ---
pca_variance_threshold = 90;
min_units_per_region   = 5;
min_speed_cms          = 2;
n_shuffles             = 5;          % light null per bin

% --- Shuffle type for the per-bin null ---
%   'rowshift' : non-circular row shift of Y_win by k_sh rows. Preserves
%                each region's autocorrelation; breaks only the
%                cross-region temporal alignment. Recommended for IFI
%                because the small-sample bias in canoncorr is matched
%                between real and null.
%   'random'   : independent random permutation of Y_win's rows.
%                Destroys Y's autocorrelation as well as the alignment.
%                Stronger null; tests "any relationship at all" rather
%                than "any temporal alignment." Faster, simpler.
shuffle_type           = 'rowshift';

% --- Hard cap on PCA dimensionality per region ---
% Window is 11 bins; canoncorr requires n_window > k1 + k2 and full rank
% Xc, Yc. With k1=k2=3 the constraint is 11 > 6 with rank 3 in 11 samples
% (generically fine even for sparse regions like CA3 with 8 units at 50 ms).
% Bump to 4 if you want more signal capacity at the cost of more
% rank-deficient windows.
max_k_per_region       = 3;          % k1 + k2 <= 6 << n_window=11

% --- Epochs ---
n_trials_epoch         = 10;

% --- Spatial alignment (post-hoc binning) ---
n_pos_bins             = 200;         % spatial bins along trial (2.5 cm each)
track_length_cm        = 500;        % from spatial code: 200 bins x 2.5 cm
pos_edges              = linspace(0, track_length_cm, n_pos_bins + 1);

% --- Reward-zone alignment (PLACEHOLDER) ---
% Spatial code defines landmarks_cm = [0,50; 100,125; 175,200; 250,275;
% 325,350; 400,425]. Reward zone is at 200 cm.
rz_entry_cm            = 200;        
rz_window_ms           = 500;        % +- 500 ms around RZ entry
rz_window_bins_half    = round(rz_window_ms / target_bin_ms);

% --- Region pairs ---
area_pairs_to_analyze = {'CA1','V1';  'CA1','DG'; 'CA1','CA3'; 'CA1','RSC'; ...
                         'CA1','SUB'; 'V1','RSC'; 'RSC','SUB'; 'CA3','DG'};
n_pairs = size(area_pairs_to_analyze, 1);

% Silence the canoncorr full-rank warning; we do an explicit rank check.
warning('off', 'stats:canoncorr:NotFullRank');

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
% For each pair, per epoch, store (animal x time/space) traces and scalars.
% We pre-allocate the "epoch matrices" lazily inside the animal loop.
epoch_names = {'early','pre','post'};

group_results = struct();
for ipair = 1:n_pairs
    group_results(ipair).pair_name = sprintf('%s-%s', ...
        area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});

    for e = 1:numel(epoch_names)
        ep = epoch_names{e};
        % Trial-mean: one scalar per animal (mean over the epoch's trials).
        group_results(ipair).(['trial_mean_cc1_' ep])      = nan(n_animals, 1);
        group_results(ipair).(['trial_mean_cc1_sh_' ep])   = nan(n_animals, 1);
        group_results(ipair).(['trial_mean_ifi_' ep])      = nan(n_animals, 1);
        group_results(ipair).(['trial_mean_ifi_sh_' ep])   = nan(n_animals, 1);

        % Spatial-binned: n_pos_bins per animal (mean over epoch trials of
        % the trial's spatial trace).
        group_results(ipair).(['spatial_cc1_' ep])         = nan(n_animals, n_pos_bins);
        group_results(ipair).(['spatial_cc1_sh_' ep])      = nan(n_animals, n_pos_bins);
        group_results(ipair).(['spatial_ifi_' ep])         = nan(n_animals, n_pos_bins);
        group_results(ipair).(['spatial_ifi_sh_' ep])      = nan(n_animals, n_pos_bins);

        % RZ-entry-aligned: 2*rz_window_bins_half + 1 per animal.
        n_rz = 2*rz_window_bins_half + 1;
        group_results(ipair).(['rz_cc1_' ep])              = nan(n_animals, n_rz);
        group_results(ipair).(['rz_cc1_sh_' ep])           = nan(n_animals, n_rz);
        group_results(ipair).(['rz_ifi_' ep])              = nan(n_animals, n_rz);
        group_results(ipair).(['rz_ifi_sh_' ep])           = nan(n_animals, n_rz);
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

    % --- Per-unit region array ---
    if isfield(units, 'region')
        unit_regions = units.region(:);
    elseif isfield(units, 'regions_label') && length(units.regions_label) == length(units.unit_id)
        unit_regions = units.regions_label(:);
    else
        fprintf('  -> Could not find a valid per-unit region array. Skipping.\n'); continue;
    end

    if isfield(units, 'regions_label') && length(units.regions_label) < length(units.unit_id)
        animal_areas = units.regions_label(:);
    else
        animal_areas = unique(unit_regions);
    end

    % --- FS unit exclusion ---
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

    % --- Re-binning to 5 ms ---
    native_bin_s = double(D.params_main.bin_size);
    rebin_factor = round((target_bin_ms / 1000) / native_bin_s);
    if rebin_factor < 1, rebin_factor = 1; end
    fprintf('  native bin = %.4f s, rebin factor = %d (target %d ms)\n', ...
        native_bin_s, rebin_factor, target_bin_ms);

    % Keep binned_spikes in its native (uint8) type and only convert
    % per-unit slices. Calling double(D.binned_spikes) up-front allocates
    % an enormous matrix (e.g. 9.9M x 455 -> 36 GB) which causes severe
    % swap on big animals. We do the rebin per-unit, converting just one
    % unit's row at a time.
    bs_native = D.binned_spikes;       % whatever uint type it is
    bs_layout_units_first = (size(bs_native, 1) == length(units.unit_id));

    if bs_layout_units_first
        nUnits = size(bs_native, 1);
        T      = size(bs_native, 2);
    else
        nUnits = size(bs_native, 2);
        T      = size(bs_native, 1);
    end

    T_use  = floor(T / rebin_factor) * rebin_factor;
    T_re   = T_use / rebin_factor;

    spikes_re = zeros(T_re, nUnits);
    for u = 1:nUnits
        if bs_layout_units_first
            row = double(bs_native(u, 1:T_use));   % [1 x T_use]
        else
            row = double(bs_native(1:T_use, u))';  % [1 x T_use]
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

    % Position (try both common field names; default to NaN if absent).
    if isfield(D.data_behaviour, 'pos_binned_gf')
        pos_cm = double(D.data_behaviour.pos_binned_gf(:));
    elseif isfield(D.data_behaviour, 'pos_binned')
        pos_cm = double(D.data_behaviour.pos_binned(:));
    else
        pos_cm = nan(size(tr_cued));
    end

    % Free the loaded raw data: we have everything we need in
    % spikes_re / mask_cued / vel_gf / tr_cued / pos_cm / units / unit_regions.
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

    % --- Per-region PCA ---
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
        % Cap k at max_k_per_region so the 21-bin window has enough samples
        % per dimension. Floor at 3.
        k = min([max(3, k_pca), max_k_per_region, size(coeff, 2)]);

        X_centered = X - repmat(mean(X_valid, 1), T_re, 1);
        AreaActivity.(area).data = X_centered * coeff(:, 1:k);   % [T_re x k]
        AreaActivity.(area).k    = k;
        fprintf('    [PCA] %s: %d units -> k=%d (k_pca@90%%=%d, var captured = %.1f%%)\n', ...
            area, sum(u_logical), k, k_pca, cum_var(k));
    end

    % --- Trial loop: extract per-trial bins, compute per-bin CCA traces ---
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

    % Per pair, accumulate per-trial outputs (one row per trial actually run).
    pair_cache = repmat(struct('trial_id',[], 'tmean_cc1',[], 'tmean_cc1_sh',[], ...
                               'tmean_ifi',[], 'tmean_ifi_sh',[], ...
                               'spat_cc1',[], 'spat_cc1_sh',[], ...
                               'spat_ifi',[], 'spat_ifi_sh',[], ...
                               'rz_cc1',[],  'rz_cc1_sh',[], ...
                               'rz_ifi',[],  'rz_ifi_sh',[]), n_pairs, 1);

    for it = 1:numel(epoch_trial_set)
        tr_id = epoch_trial_set(it);

        % STRICT CONTIGUITY: take ALL bins of the trial (not only valid
        % ones), then require that every bin inside the sliding window
        % satisfies base_valid. No stitching across velocity-drop gaps.
        trial_in_tr   = (tc_re == tr_id);
        trial_idx     = find(trial_in_tr);          % all bins of the trial
        n_tr_bins     = numel(trial_idx);
        is_valid_tr   = base_valid(trial_idx);      % per-bin validity

        if n_tr_bins < (n_window + 5)
            continue;
        end

        % Per-pair CCA over the trial's bins.
        for ipair = 1:n_pairs
            a1 = area_pairs_to_analyze{ipair, 1};
            a2 = area_pairs_to_analyze{ipair, 2};
            if ~isfield(AreaActivity, a1) || ~isfield(AreaActivity, a2), continue; end

            X_full = AreaActivity.(a1).data;     % [T_re x k1]
            Y_full = AreaActivity.(a2).data;     % [T_re x k2]
            k1 = AreaActivity.(a1).k;
            k2 = AreaActivity.(a2).k;

            % Trial timecourse over ALL trial bins (not just valid).
            X_tr   = X_full(trial_idx, :);        % [n_tr_bins x k1]
            Y_tr   = Y_full(trial_idx, :);        % [n_tr_bins x k2]
            pos_tr = pos_re(trial_idx);           % [n_tr_bins x 1]

            % Per-bin output traces (one entry per bin in the trial; NaN
            % for centers whose +-half_window contains any invalid bin).
            cc1_tr      = nan(n_tr_bins, 1);
            ifi_tr      = nan(n_tr_bins, 1);
            cc1_tr_sh   = nan(n_tr_bins, 1);
            ifi_tr_sh   = nan(n_tr_bins, 1);

            % Diagnostic counters (per trial x pair).
            n_attempted = 0; n_skip_contig = 0; n_skip_dim = 0;
            n_skip_rank = 0; n_skip_canoncorr = 0; n_ok = 0;
            n_sh_ok = 0;

            for ic = (half_window_bins+1):stride:(n_tr_bins - half_window_bins)
                n_attempted = n_attempted + 1;
                idx_w = (ic - half_window_bins):(ic + half_window_bins);

                % Strict contiguity: skip if any bin in the window is invalid.
                if ~all(is_valid_tr(idx_w)), n_skip_contig = n_skip_contig + 1; continue; end

                Xw = X_tr(idx_w, :);
                Yw = Y_tr(idx_w, :);

                % Rank check.
                Xc = Xw - repmat(mean(Xw,1), n_window, 1);
                Yc = Yw - repmat(mean(Yw,1), n_window, 1);
                if n_window <= size(Xw,2) + size(Yw,2)
                    n_skip_dim = n_skip_dim + 1; continue;
                end
                if rank(Xc) < size(Xw,2) || rank(Yc) < size(Yw,2)
                    n_skip_rank = n_skip_rank + 1; continue;
                end

                try
                    [A, B, r] = canoncorr(Xw, Yw);
                catch
                    n_skip_canoncorr = n_skip_canoncorr + 1; continue;
                end
                cc1_tr(ic) = r(1);
                n_ok = n_ok + 1;

                % Lag corr / IFI within this window's projections.
                u = Xc * A(:,1);                  % [n_window x 1]
                v = Yc * B(:,1);                  % [n_window x 1]
                r_lags = lag_corr_simple(u, v, lags);
                ifi_tr(ic) = ifi_from_lags(r_lags, lags);

                % --- Lightweight shuffle null ---
                % Two modes (set by shuffle_type at top of script):
                %   'rowshift' : non-circular row shift of Y_win, preserving
                %                Y's autocorrelation. Need shift in
                %                [max_lag_abs+1, n_window - max(k1,k2) - 3].
                %   'random'   : independent random permutation of Y_win's
                %                rows. No trimming, no shift constraint.
                shuffle_can_run = false;
                if strcmpi(shuffle_type, 'rowshift')
                    sh_lo = max_lag_abs + 1;
                    sh_hi = n_window - max(k1, k2) - 3;
                    shuffle_can_run = (sh_hi >= sh_lo);
                else  % 'random'
                    shuffle_can_run = (n_window >= max(k1, k2) + 3);
                end

                if shuffle_can_run
                    r_sh = nan(n_shuffles,1); p_sh = nan(n_shuffles,1);
                    for ish = 1:n_shuffles

                        % --- Build the (Xk, Yk) pair for this shuffle ---
                        if strcmpi(shuffle_type, 'rowshift')
                            k_sh = randi([sh_lo, sh_hi]);
                            if rand() < 0.5, k_sh = -k_sh; end
                            if k_sh >= 0
                                ix = (k_sh+1):n_window;
                                iy = 1:(n_window - k_sh);
                            else
                                ix = 1:(n_window + k_sh);
                                iy = (-k_sh + 1):n_window;
                            end
                            n_keep = numel(ix);
                            if n_keep < (max(k1,k2) + 3), continue; end
                            Xk = Xw(ix, :);
                            Yk = Yw(iy, :);
                        else  % 'random'
                            perm = randperm(n_window);
                            Xk = Xw;
                            Yk = Yw(perm, :);
                            n_keep = n_window;
                        end

                        % --- Rank check on the null pair ---
                        Xk_c = Xk - repmat(mean(Xk,1), n_keep, 1);
                        Yk_c = Yk - repmat(mean(Yk,1), n_keep, 1);
                        if rank(Xk_c) < size(Xk,2) || rank(Yk_c) < size(Yk,2), continue; end

                        try
                            [As, Bs, rs] = canoncorr(Xk, Yk);
                        catch
                            continue;
                        end
                        r_sh(ish) = rs(1);

                        % IFI null: project the SAME (real) window onto the
                        % null canonical directions, then compute lag corr.
                        % For 'random', mean(Yk) == mean(Yw) (permutation
                        % preserves the mean), so this reduces to the
                        % standard centering.
                        u_sh = (Xw - repmat(mean(Xk,1), n_window, 1)) * As(:,1);
                        v_sh = (Yw - repmat(mean(Yk,1), n_window, 1)) * Bs(:,1);
                        rl_sh = lag_corr_simple(u_sh, v_sh, lags);
                        p_sh(ish) = ifi_from_lags(rl_sh, lags);
                    end
                    cc1_tr_sh(ic) = mean(r_sh, 'omitnan');
                    ifi_tr_sh(ic) = mean(p_sh, 'omitnan');
                    if isfinite(cc1_tr_sh(ic)), n_sh_ok = n_sh_ok + 1; end
                end
            end

            % Diagnostic for the first trial of each animal x pair.
            if it == 1
                fprintf('    [diag] %s-%s tr=%d k1=%d k2=%d  attempted=%d  contig=%d  dim=%d  rank=%d  cc=%d  ok=%d  sh_ok=%d\n', ...
                    a1, a2, tr_id, k1, k2, n_attempted, n_skip_contig, n_skip_dim, n_skip_rank, n_skip_canoncorr, n_ok, n_sh_ok);
            end

            % --- Aggregate per trial ---
            % (a) trial mean
            tmean_cc1     = mean(cc1_tr,    'omitnan');
            tmean_cc1_sh  = mean(cc1_tr_sh, 'omitnan');
            tmean_ifi     = mean(ifi_tr,    'omitnan');
            tmean_ifi_sh  = mean(ifi_tr_sh, 'omitnan');

            % (b) spatial bin: bin trial bins by position, average within bin.
            spat_cc1    = bin_by_position(cc1_tr,    pos_tr, pos_edges);
            spat_cc1_sh = bin_by_position(cc1_tr_sh, pos_tr, pos_edges);
            spat_ifi    = bin_by_position(ifi_tr,    pos_tr, pos_edges);
            spat_ifi_sh = bin_by_position(ifi_tr_sh, pos_tr, pos_edges);

            % (c) RZ-entry-aligned: index the +- rz_window around the first
            %     bin where pos crosses rz_entry_cm. Pad with NaN if window
            %     extends beyond trial bounds.
            rz_cc1    = rz_align(cc1_tr,    pos_tr, rz_entry_cm, rz_window_bins_half);
            rz_cc1_sh = rz_align(cc1_tr_sh, pos_tr, rz_entry_cm, rz_window_bins_half);
            rz_ifi    = rz_align(ifi_tr,    pos_tr, rz_entry_cm, rz_window_bins_half);
            rz_ifi_sh = rz_align(ifi_tr_sh, pos_tr, rz_entry_cm, rz_window_bins_half);

            % --- Stash ---
            pair_cache(ipair).trial_id(end+1)        = tr_id;
            pair_cache(ipair).tmean_cc1(end+1,1)     = tmean_cc1;
            pair_cache(ipair).tmean_cc1_sh(end+1,1)  = tmean_cc1_sh;
            pair_cache(ipair).tmean_ifi(end+1,1)     = tmean_ifi;
            pair_cache(ipair).tmean_ifi_sh(end+1,1)  = tmean_ifi_sh;

            pair_cache(ipair).spat_cc1(end+1,:)      = spat_cc1;
            pair_cache(ipair).spat_cc1_sh(end+1,:)   = spat_cc1_sh;
            pair_cache(ipair).spat_ifi(end+1,:)      = spat_ifi;
            pair_cache(ipair).spat_ifi_sh(end+1,:)   = spat_ifi_sh;

            pair_cache(ipair).rz_cc1(end+1,:)        = rz_cc1;
            pair_cache(ipair).rz_cc1_sh(end+1,:)     = rz_cc1_sh;
            pair_cache(ipair).rz_ifi(end+1,:)        = rz_ifi;
            pair_cache(ipair).rz_ifi_sh(end+1,:)     = rz_ifi_sh;
        end

        if mod(it, 5) == 0 || it == numel(epoch_trial_set)
            fprintf('    trial %d / %d done (t=%.1f s)\n', it, numel(epoch_trial_set), toc(t_animal));
        end
    end

    % --- Reduce per-pair-per-epoch over trials ---
    for ipair = 1:n_pairs
        if isempty(pair_cache(ipair).trial_id), continue; end

        tids = pair_cache(ipair).trial_id;
        for e = 1:numel(epoch_names)
            ep = epoch_names{e};
            sel = ismember(tids, epoch_trials.(ep));
            if ~any(sel), continue; end

            % Trial-mean (mean across the epoch's trials of the per-trial mean).
            group_results(ipair).(['trial_mean_cc1_' ep])(ianimal)    = mean(pair_cache(ipair).tmean_cc1(sel),    'omitnan');
            group_results(ipair).(['trial_mean_cc1_sh_' ep])(ianimal) = mean(pair_cache(ipair).tmean_cc1_sh(sel), 'omitnan');
            group_results(ipair).(['trial_mean_ifi_' ep])(ianimal)    = mean(pair_cache(ipair).tmean_ifi(sel),    'omitnan');
            group_results(ipair).(['trial_mean_ifi_sh_' ep])(ianimal) = mean(pair_cache(ipair).tmean_ifi_sh(sel), 'omitnan');

            % Spatial.
            group_results(ipair).(['spatial_cc1_' ep])(ianimal,:)    = mean(pair_cache(ipair).spat_cc1(sel,:),    1, 'omitnan');
            group_results(ipair).(['spatial_cc1_sh_' ep])(ianimal,:) = mean(pair_cache(ipair).spat_cc1_sh(sel,:), 1, 'omitnan');
            group_results(ipair).(['spatial_ifi_' ep])(ianimal,:)    = mean(pair_cache(ipair).spat_ifi(sel,:),    1, 'omitnan');
            group_results(ipair).(['spatial_ifi_sh_' ep])(ianimal,:) = mean(pair_cache(ipair).spat_ifi_sh(sel,:), 1, 'omitnan');

            % RZ-aligned.
            group_results(ipair).(['rz_cc1_' ep])(ianimal,:)    = mean(pair_cache(ipair).rz_cc1(sel,:),    1, 'omitnan');
            group_results(ipair).(['rz_cc1_sh_' ep])(ianimal,:) = mean(pair_cache(ipair).rz_cc1_sh(sel,:), 1, 'omitnan');
            group_results(ipair).(['rz_ifi_' ep])(ianimal,:)    = mean(pair_cache(ipair).rz_ifi(sel,:),    1, 'omitnan');
            group_results(ipair).(['rz_ifi_sh_' ep])(ianimal,:) = mean(pair_cache(ipair).rz_ifi_sh(sel,:), 1, 'omitnan');
        end
    end

    % --- Save partial progress every animal so a hang doesn't lose work ---
    save(save_path, 'group_results', 'is_learner', 'analysis_lp', ...
                    'pos_edges', 'rz_window_bins_half', 'target_bin_ms', ...
                    'lags', 'half_window_bins', 'stride', 'rz_entry_cm', ...
                    'shuffle_type', 'n_shuffles', 'max_k_per_region', '-v7.3');
    fprintf('  saved partial results -> %s  (animal took %.1f s)\n', save_path, toc(t_animal));
end

fprintf('\nFinished! Results saved to: %s\n', save_path);

%% LOCAL HELPERS

function r_lags = lag_corr_simple(u, v, lags)
% In-window lag correlations (no validity mask: window is contiguous bins
% within a single trial, so all positions are valid by construction).
%   L >= 0:  corr(u(1:end-L), v(1+L:end))    -> u leads v
%   L <  0:  corr(u(1-L:end), v(1:end+L))    -> v leads u
    n = length(u); u = u(:); v = v(:);
    r_lags = nan(1, length(lags));
    for il = 1:length(lags)
        L = lags(il);
        if L >= 0
            a = u(1:(n-L));      b = v((1+L):n);
        else
            a = u((1-L):n);      b = v(1:(n+L));
        end
        if numel(a) > 5
            r_lags(il) = corr(a, b, 'rows', 'complete');
        end
    end
end

function ifi = ifi_from_lags(r_lags, lags)
    r_neg = mean(r_lags(lags < 0), 'omitnan');
    r_pos = mean(r_lags(lags > 0), 'omitnan');
    if (r_neg + r_pos) > 0.001
        ifi = (r_neg - r_pos) / (r_neg + r_pos);
    else
        ifi = NaN;
    end
end

function out = bin_by_position(x_per_bin, pos_per_bin, edges)
% Average x_per_bin into the spatial bins defined by edges.
% Returns a 1 x (numel(edges)-1) vector with NaN where no bins fell in
% that spatial bin.
    n_bins = numel(edges) - 1;
    out = nan(1, n_bins);
    if all(isnan(x_per_bin)) || all(isnan(pos_per_bin)), return; end
    bin_idx = discretize(pos_per_bin, edges);
    for b = 1:n_bins
        sel = (bin_idx == b) & ~isnan(x_per_bin(:));
        if any(sel)
            out(b) = mean(x_per_bin(sel), 'omitnan');
        end
    end
end

function out = rz_align(x_per_bin, pos_per_bin, rz_cm, half_bins)
% Align x_per_bin to reward-zone entry. RZ entry = first bin where
% pos_per_bin >= rz_cm. Returns a 1 x (2*half_bins+1) vector. NaN where
% the alignment window extends past trial bounds, or if RZ entry not
% found in this trial.
    n_out = 2*half_bins + 1;
    out = nan(1, n_out);
    if all(isnan(pos_per_bin)), return; end
    rz_idx = find(pos_per_bin(:) >= rz_cm, 1, 'first');
    if isempty(rz_idx), return; end

    src_lo = rz_idx - half_bins;
    src_hi = rz_idx + half_bins;
    n_tr = numel(x_per_bin);

    dst_lo = 1; dst_hi = n_out;
    if src_lo < 1
        dst_lo = 1 - src_lo + 1;
        src_lo = 1;
    end
    if src_hi > n_tr
        dst_hi = dst_hi - (src_hi - n_tr);
        src_hi = n_tr;
    end
    if src_lo > src_hi || dst_lo > dst_hi, return; end
    out(dst_lo:dst_hi) = x_per_bin(src_lo:src_hi);
end
