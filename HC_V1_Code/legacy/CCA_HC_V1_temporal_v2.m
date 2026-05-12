%% CCA_HC_V1_Temporal_Continuous_Windowed.m
% Continuous Temporal CCA with Sliding Trial Windows and Masked Lags
%
% DESCRIPTION:
% Operates on continuous binned time, applying a sliding window defined by 
% trial labels to dynamically calculate CC1 and the Information Flow Index (IFI).

%% 1. SETUP & PARAMETERS
clear; clc; close all;

base_dir       = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir       = fullfile(base_dir, data_subfolder);
file_pattern   = 'TF*_export.mat';
learning_file  = 'animal_behaviour.mat';

current_date = datestr(now, 'yyyy_mm_dd');
save_path    = fullfile(data_dir, sprintf('Temporal_Continuous_CCA_%s.mat', current_date));

% --- Temporal & Analysis Parameters ---
target_bin_ms          = 50;        
lags                   = -5:5;      
pca_variance_threshold = 90;
n_trials_window        = -3:3;      
n_shuffles             = 25;
min_units_per_region   = 5;
min_speed_cms          = 2;         
min_samples_window     = 20;        

area_pairs_to_analyze = {'CA1','V1';  'CA1','DG'; 'CA1','CA3'; 'CA1','RSC'; ...
                         'CA1','SUB'; 'V1','RSC'; 'RSC','SUB'; 'CA3','DG'};
n_pairs = size(area_pairs_to_analyze, 1);

if isempty(gcp('nocreate'))
    parpool;
end

% canoncorr will print "X is not full rank" whenever a trial window happens
% to be near rank-deficient. We do an explicit rank check below and NaN out
% those windows, so silence that specific warning to keep logs readable.
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
group_results = struct( ...
    'pair_name', cell(n_pairs, 1), ...
    'trial_corr_early', cell(n_pairs, 1), 'trial_corr_pre', cell(n_pairs, 1), 'trial_corr_post', cell(n_pairs, 1), ...
    'trial_corr_early_shuff', cell(n_pairs, 1), 'trial_corr_pre_shuff', cell(n_pairs, 1), 'trial_corr_post_shuff', cell(n_pairs, 1), ...
    'trial_precession_early_idx', cell(n_pairs, 1), 'trial_precession_pre_idx', cell(n_pairs, 1), 'trial_precession_post_idx', cell(n_pairs, 1), ...
    'trial_precession_early_idx_shuff', cell(n_pairs, 1), 'trial_precession_pre_idx_shuff', cell(n_pairs, 1), 'trial_precession_post_idx_shuff', cell(n_pairs, 1));

for ipair = 1:n_pairs
    group_results(ipair).pair_name = sprintf('%s-%s', area_pairs_to_analyze{ipair, 1}, area_pairs_to_analyze{ipair, 2});
end

%% 4. MAIN PIPELINE
for ianimal = 1:n_animals
    filename = file_list(ianimal).name;
    fullpath = fullfile(data_dir, filename);
    fprintf('\nProcessing Animal %d/%d: %s\n', ianimal, n_animals, filename);

    D = load(fullpath);
    if ~isfield(D, 'units') || ~isfield(D, 'binned_spikes')
        fprintf('  -> Missing fields. Skipping.\n'); continue; 
    end
    units = D.units;

    % Identify the true per-unit array vs the unique labels
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

    % 4.1 Filter FS units
    keep_mask = true(length(units.unit_id), 1);
    if isfield(units, 'idx_fs')
        target_fs_areas = {'V1','RSC','CA1','CA3'};
        is_fs = logical(units.idx_fs(:)); 
        
        if length(is_fs) == length(unit_regions)
            for r = 1:length(target_fs_areas)
                fs_in_area = strcmp(unit_regions, target_fs_areas{r}) & is_fs;
                keep_mask(fs_in_area) = false;
            end
        else
            fprintf('  -> [Warning] idx_fs length (%d) != units.region length (%d). Skipping FS filter.\n', length(is_fs), length(unit_regions));
        end
    end

    % 4.2 Re-binning
    native_bin_s = double(D.params_main.bin_size);
    rebin_factor = round((target_bin_ms / 1000) / native_bin_s);
    if rebin_factor < 1, rebin_factor = 1; end
    
    spikes = double(D.binned_spikes);
    if size(spikes, 1) ~= length(units.unit_id), spikes = spikes'; end
    
    nUnits = size(spikes, 1);
    T      = size(spikes, 2);
    T_use  = floor(T / rebin_factor) * rebin_factor;
    T_re   = T_use / rebin_factor;
    
    spikes_re = nan(T_re, nUnits);
    for u = 1:nUnits
        spikes_re(:, u) = sum(reshape(spikes(u, 1:T_use), rebin_factor, T_re), 1)';
    end
    
    mask_cued = logical(D.analysis_behaviour.masks.tunnel_cued(:));
    if isfield(D.data_behaviour, 'velocity_binned_gf')
        vel_gf = double(D.data_behaviour.velocity_binned_gf(:));
    else
        vel_gf = double(D.data_behaviour.velocity_gf(:));
    end
    tr_cued = double(D.data_behaviour.trial_binned_cued(:));

    mask_cued = mask_cued(1:T_use);
    vel_gf    = vel_gf(1:T_use);
    tr_cued   = tr_cued(1:T_use);

    mc_re  = all(reshape(mask_cued, rebin_factor, T_re), 1)';
    vel_re = mean(reshape(vel_gf,  rebin_factor, T_re), 1)';
    tc_re  = mode(reshape(tr_cued, rebin_factor, T_re), 1)';
    
    base_valid = mc_re & (vel_re >= min_speed_cms);

    mu = mean(spikes_re(base_valid, :), 1, 'omitnan');
    sd = std(spikes_re(base_valid, :), 0, 1, 'omitnan'); sd(sd==0) = 1;
    
    spikes_z = (spikes_re - repmat(mu, size(spikes_re, 1), 1)) ./ repmat(sd, size(spikes_re, 1), 1);

    % 4.3 PCA Generation
    AreaActivity = struct();
    for ia = 1:length(animal_areas)
        area = char(animal_areas{ia});
        if isempty(area), continue; end
        
        u_logical = strcmp(unit_regions, area) & keep_mask(:);
        if sum(u_logical) < min_units_per_region, continue; end
        
        X = spikes_z(:, u_logical);
        X_valid = X(base_valid, :);
        if size(X_valid, 1) < min_samples_window, continue; end
        
        [coeff, ~, ~, ~, explained] = pca(X_valid);
        cum_var = cumsum(explained);
        k_pca = find(cum_var >= pca_variance_threshold, 1);
        if isempty(k_pca), k_pca = size(coeff, 2); end
        k = min(max(3, k_pca), size(coeff, 2));  % floor at 3, cap at available PCs
        
        X_centered = X - repmat(mean(X_valid, 1), size(X, 1), 1);
        AreaActivity.(area).data = X_centered * coeff(:, 1:k);
        AreaActivity.(area).k = k;
    end

    % 4.4 Sliding Window Trial Loop
    num_trials = max(tc_re);
    
    for ipair = 1:n_pairs
        a1 = area_pairs_to_analyze{ipair, 1}; a2 = area_pairs_to_analyze{ipair, 2};
        if ~isfield(AreaActivity, a1) || ~isfield(AreaActivity, a2), continue; end
        
        X_full = AreaActivity.(a1).data; k1 = AreaActivity.(a1).k;
        Y_full = AreaActivity.(a2).data; k2 = AreaActivity.(a2).k;
        T_sess = size(X_full, 1);

        cca_tr        = nan(1, num_trials);
        ifi_tr        = nan(1, num_trials);
        cca_tr_shuff  = nan(1, num_trials);
        ifi_tr_shuff  = nan(1, num_trials);

        % parfor over trials — each trial window is independent.
        % The shuffle loop inside is a serial for (parfor cannot be nested).
        %
        % Memory note: do NOT circshift Y_full inside the inner shuffle loop.
        % Y_full is the full-session matrix [T_sess x k2] and copying it
        % n_shuffles * num_trials times per pair causes severe memory
        % pressure and effectively hangs the script. The correct null is to
        % shuffle the in-window block of Y, which is what v1 did.
        max_lag_abs = max(abs(lags));
        min_window_samples = max(k1, k2) + max_lag_abs + 5;

        parfor t = 1:num_trials
            win_bounds   = t + n_trials_window;
            win_bounds   = win_bounds(win_bounds >= 1 & win_bounds <= num_trials);
            valid_window = base_valid & ismember(tc_re, win_bounds);

            t_cca = nan; t_ifi = nan; t_cca_sh = nan; t_ifi_sh = nan;

            n_win = sum(valid_window);
            if n_win >= min_window_samples
                X_win = X_full(valid_window, :);
                Y_win = Y_full(valid_window, :);

                % --- Explicit rank check (silences canoncorr:NotFullRank
                %     for the dominant cause: short / near-degenerate windows).
                Xc_chk = X_win - repmat(mean(X_win,1), n_win, 1);
                Yc_chk = Y_win - repmat(mean(Y_win,1), n_win, 1);
                rank_ok = (n_win > size(X_win,2) + size(Y_win,2)) && ...
                          (rank(Xc_chk) == size(X_win,2)) && ...
                          (rank(Yc_chk) == size(Y_win,2));

                U_loc = []; V_loc = [];
                if rank_ok
                    try
                        [A, B, r] = canoncorr(X_win, Y_win);
                        t_cca = r(1);
                        U_loc = (X_full - repmat(mean(X_win, 1), T_sess, 1)) * A(:,1);
                        V_loc = (Y_full - repmat(mean(Y_win, 1), T_sess, 1)) * B(:,1);
                        t_ifi = calc_ifi(calc_masked_lags(U_loc, V_loc, valid_window, lags), lags);
                    catch
                    end
                end

                % --- Shuffle control: NON-CIRCULAR shift of the in-window Y.
                % We misalign X_win and Y_win by k_sh valid-samples (drop the
                % overhanging rows on both sides) and refit canoncorr on the
                % misaligned pair. The fitted Bs is then used to project the
                % full Y_full onto its first canonical direction, exactly as
                % in the real path. Lag correlations use the original
                % valid_window so the only thing different from the real path
                % is the misaligned fit.
                %
                % Why non-circular: a circshift would create a wrap junction
                % (last row glued to first row) which is not a real signal.
                % A plain shift drops the overhanging rows instead.
                if ~isempty(U_loc) && n_win > 2*max_lag_abs + 20
                    r_sh = nan(n_shuffles, 1); p_sh = nan(n_shuffles, 1);

                    % Shift range: at least max_lag_abs+1 so the lag pattern
                    % is broken, at most ~half the window so the kept block
                    % is still long enough to fit canoncorr.
                    sh_lo = max_lag_abs + 1;
                    sh_hi = max(sh_lo+1, floor(n_win/2));

                    for ish = 1:n_shuffles
                        % Random sign so the shift can go either direction.
                        k_sh = randi([sh_lo, sh_hi]);
                        if rand() < 0.5, k_sh = -k_sh; end

                        % Misalign by k_sh valid-samples and drop overhang.
                        if k_sh >= 0
                            ix = (k_sh+1):n_win;
                            iy = 1:(n_win - k_sh);
                        else
                            ix = 1:(n_win + k_sh);
                            iy = (-k_sh + 1):n_win;
                        end
                        n_keep = numel(ix);
                        if n_keep < (max(k1,k2) + max_lag_abs + 5), continue; end

                        X_keep = X_win(ix, :);
                        Y_keep = Y_win(iy, :);

                        % Rank check on the misaligned block.
                        Xk_c = X_keep - repmat(mean(X_keep,1), n_keep, 1);
                        Yk_c = Y_keep - repmat(mean(Y_keep,1), n_keep, 1);
                        if rank(Xk_c) < size(X_keep,2) || rank(Yk_c) < size(Y_keep,2)
                            continue;
                        end

                        try
                            [~, Bs, rs] = canoncorr(X_keep, Y_keep);
                            r_sh(ish)   = rs(1);

                            % Project full-session Y onto Bs(:,1), centered
                            % the same way as the real path.
                            V_loc_sh = (Y_full - repmat(mean(Y_keep,1), T_sess, 1)) * Bs(:,1);
                            p_sh(ish) = calc_ifi(calc_masked_lags(U_loc, V_loc_sh, valid_window, lags), lags);
                        catch
                        end
                    end
                    t_cca_sh = mean(r_sh, 'omitnan');
                    t_ifi_sh = mean(p_sh, 'omitnan');
                end
            end

            cca_tr(t)       = t_cca;
            ifi_tr(t)       = t_ifi;
            cca_tr_shuff(t) = t_cca_sh;
            ifi_tr_shuff(t) = t_ifi_sh;
        end
        
        % 4.5 Extract Epochs via Yoked LP
        lp = analysis_lp(ianimal);
        if ~isnan(lp) && lp > 10 && (lp + 9) <= num_trials
            idx_early = 1:10; idx_pre = (lp - 10) : (lp - 1); idx_post = lp : (lp + 9);
            
            group_results(ipair).trial_corr_early{ianimal} = cca_tr(idx_early);
            group_results(ipair).trial_corr_pre{ianimal}   = cca_tr(idx_pre);
            group_results(ipair).trial_corr_post{ianimal}  = cca_tr(idx_post);
            
            group_results(ipair).trial_precession_early_idx{ianimal} = ifi_tr(idx_early);
            group_results(ipair).trial_precession_pre_idx{ianimal}   = ifi_tr(idx_pre);
            group_results(ipair).trial_precession_post_idx{ianimal}  = ifi_tr(idx_post);
            
            group_results(ipair).trial_corr_early_shuff{ianimal} = cca_tr_shuff(idx_early);
            group_results(ipair).trial_corr_pre_shuff{ianimal}   = cca_tr_shuff(idx_pre);
            group_results(ipair).trial_corr_post_shuff{ianimal}  = cca_tr_shuff(idx_post);
            
            group_results(ipair).trial_precession_early_idx_shuff{ianimal} = ifi_tr_shuff(idx_early);
            group_results(ipair).trial_precession_pre_idx_shuff{ianimal}   = ifi_tr_shuff(idx_pre);
            group_results(ipair).trial_precession_post_idx_shuff{ianimal}  = ifi_tr_shuff(idx_post);
        end
    end
end

%% 5. PADDING
n_total_animals = length(is_learner);
for ipair = 1:n_pairs
    fields = fieldnames(group_results);
    for f = 1:length(fields)
        if iscell(group_results(ipair).(fields{f})) && length(group_results(ipair).(fields{f})) < n_total_animals
            group_results(ipair).(fields{f}){n_total_animals} = [];
        end
    end
end

save(save_path, 'group_results', 'is_learner', 'analysis_lp', '-v7.3');
fprintf('\nFinished! Results saved to: %s\n', save_path);

%% =========================================================================
%% LOCAL HELPERS
%% =========================================================================

function r_lags = calc_masked_lags(U, V, valid_window, lags)
% Masked lag correlations using direct index arithmetic (no circshift).
% For lag L: correlates U[t] with V[t+L], keeping only pairs where both
% indices lie within the session and both timepoints are valid.
%   L < 0  =>  V leads U  (used for r_neg in IFI)
%   L > 0  =>  U leads V  (used for r_pos in IFI)
    T = length(U);
    U = U(:); V = V(:);
    valid_window = logical(valid_window(:));
    r_lags = nan(1, length(lags));
    for il = 1:length(lags)
        L = lags(il);
        if L >= 0
            t_u = 1:(T - L);
            t_v = (1 + L):T;
        else
            t_u = (1 - L):T;   % (1+|L|) to T
            t_v = 1:(T + L);   % 1 to (T-|L|)
        end
        valid_pairs = valid_window(t_u) & valid_window(t_v);
        if sum(valid_pairs) > 10
            r_lags(il) = corr(U(t_u(valid_pairs)), V(t_v(valid_pairs)), 'rows', 'complete');
        end
    end
end

function ifi = calc_ifi(r_lags, lags)
    neg_lags = r_lags(lags < 0);
    pos_lags = r_lags(lags > 0);
    
    r_neg = mean(neg_lags, 'omitnan');
    r_pos = mean(pos_lags, 'omitnan');
    
    if (r_neg + r_pos) > 0.001
        ifi = (r_neg - r_pos) / (r_neg + r_pos);
    else
        ifi = nan;
    end
end