%% CCA_HC_V1_temporal_v3_plot.m
% Plotting / verification script for v3 temporal CCA results.
%
% Loads the latest Temporal_CCA_v3_<date>.mat and produces:
%   1. Trial-mean bars (CC1, IFI) per pair, early/pre/post,
%      learners vs non-learners side by side, with shuffle null overlay.
%   2. Spatial trajectories (CC1, IFI) per pair, one line per epoch,
%      learners vs non-learners as separate panels per pair.
%   3. RZ-entry-aligned trajectories (CC1, IFI) per pair, one line per epoch.
%   4. Diagnostic / sanity panel: fraction of finite animals per pair x
%      epoch, real vs shuffle scatter for CC1 trial-mean, n contributing.
%
% All figures are saved as .svg into the same data folder.

%% 1. Locate latest results file
clear; clc; close all;

base_dir       = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir       = fullfile(base_dir, data_subfolder);

files = dir(fullfile(data_dir, 'Temporal_CCA_v3_*.mat'));
if isempty(files)
    error('No Temporal_CCA_v3_*.mat files found in %s', data_dir);
end
[~, idx] = max([files.datenum]);
results_path = fullfile(data_dir, files(idx).name);
fprintf('Loading %s\n', results_path);
S = load(results_path);

group_results       = S.group_results;
is_learner          = logical(S.is_learner(:));
analysis_lp         = S.analysis_lp;
pos_edges           = S.pos_edges(:)';
target_bin_ms       = double(S.target_bin_ms);
half_window_bins    = double(S.half_window_bins);
rz_window_bins_half = double(S.rz_window_bins_half);
lags                = double(S.lags(:)');

n_pairs   = numel(group_results);
n_animals = numel(is_learner);

% Build x axes for spatial and RZ plots.
pos_centers = (pos_edges(1:end-1) + pos_edges(2:end)) / 2;            % cm
rz_t_ms     = (-rz_window_bins_half:rz_window_bins_half) * target_bin_ms;

% Landmark positions (from spatial code).
landmarks_cm = [0,50; 100,125; 175,200; 250,275; 325,350; 400,425];
rz_entry_cm  = 400;  % start of last landmark = reward zone

% --- Style ---
col_L  = [0.85 0.33 0.10];   % Learners
col_NL = [0.30 0.30 0.30];   % Non-Learners
col_real = [0.10 0.45 0.85];
col_shuf = [0.55 0.55 0.55];
ep_names = {'Early','Pre','Post'};
ep_keys  = {'early','pre','post'};
ep_marker = {'o','s','^'};

fprintf('Found %d pairs, %d animals, %d learners.\n', n_pairs, n_animals, sum(is_learner));

%% A) TRIAL-MEAN BARS (CC1, IFI) — learners vs non-learners
for which_metric = 1:2
    if which_metric == 1
        metric_real = 'trial_mean_cc1_';
        metric_shuf = 'trial_mean_cc1_sh_';
        metric_label = 'CC1';
        ylab = 'CC1 (real - shuffle overlay)';
    else
        metric_real = 'trial_mean_ifi_';
        metric_shuf = 'trial_mean_ifi_sh_';
        metric_label = 'IFI';
        ylab = 'IFI (real - shuffle overlay)';
    end

    figure('Name', sprintf('Trial-mean %s bars', metric_label), ...
           'Color','w','Position',[100 100 1500 800]);
    tiledlayout('flow','TileSpacing','compact','Padding','compact');

    for ipair = 1:n_pairs
        nexttile; hold on;
        gr = group_results(ipair);

        % Real values per epoch [animals]
        re = gr.([metric_real 'early']);
        rp = gr.([metric_real 'pre']);
        rx = gr.([metric_real 'post']);

        % Shuffle nulls
        se = gr.([metric_shuf 'early']);
        sp = gr.([metric_shuf 'pre']);
        sx = gr.([metric_shuf 'post']);

        plot_bars_pair(re, rp, rx, se, sp, sx, is_learner, ...
                       col_L, col_NL, col_shuf);

        title(gr.pair_name, 'Interpreter','none');
        ylabel(ylab);
        xticks(1:3); xticklabels(ep_names); xlim([0.5 3.5]);
        if ipair == 1
            legend({'Learners (real)','Learners (shuf)', ...
                    'NonL (real)','NonL (shuf)'}, ...
                    'Location','best','Interpreter','none');
        end
    end
    sgtitle(sprintf('%s: trial-mean per epoch  (real circles, shuffle dashed)', metric_label));
    save_fig(gcf, fullfile(data_dir, sprintf('TemporalCCAv3_TrialMean_%s.svg', metric_label)));
end

%% B) SPATIAL TRAJECTORIES (CC1, IFI) per pair, one panel per pair
for which_metric = 1:2
    if which_metric == 1
        metric_real = 'spatial_cc1_';
        metric_shuf = 'spatial_cc1_sh_';
        metric_label = 'CC1';
    else
        metric_real = 'spatial_ifi_';
        metric_shuf = 'spatial_ifi_sh_';
        metric_label = 'IFI';
    end

    figure('Name', sprintf('Spatial %s', metric_label), ...
           'Color','w','Position',[100 100 1600 900]);
    tiledlayout('flow','TileSpacing','compact','Padding','compact');

    for ipair = 1:n_pairs
        nexttile; hold on;
        gr = group_results(ipair);

        for iep = 1:3
            mat   = gr.([metric_real ep_keys{iep}]);   % [animals x n_pos_bins]
            mat_s = gr.([metric_shuf ep_keys{iep}]);
            plot_spatial_curve(pos_centers, mat,   is_learner,  col_for_epoch(iep,'L'),  ep_names{iep},  '-');
            plot_spatial_curve(pos_centers, mat_s, is_learner,  col_for_epoch(iep,'L'),  '',             ':');
        end
        % Landmarks as faint vertical bands.
        for il = 1:size(landmarks_cm,1)
            patch([landmarks_cm(il,1) landmarks_cm(il,2) landmarks_cm(il,2) landmarks_cm(il,1)], ...
                  [-1 -1 1 1], [0.9 0.9 0.95], 'EdgeColor','none', 'FaceAlpha',0.3, 'HandleVisibility','off');
        end
        xline(rz_entry_cm, 'k--', 'HandleVisibility','off');
        if which_metric == 1
            ylim([0 1]);
        else
            yline(0,'k:','HandleVisibility','off');
            ylim([-1 1]);
        end
        title(sprintf('%s (Learners)', gr.pair_name), 'Interpreter','none');
        xlabel('Position (cm)'); ylabel(metric_label);
    end
    sgtitle(sprintf('%s along position - Learners (solid: real, dotted: shuffle)', metric_label));
    save_fig(gcf, fullfile(data_dir, sprintf('TemporalCCAv3_Spatial_%s_Learners.svg', metric_label)));

    % Non-learners variant (separate figure for clarity).
    if any(~is_learner)
        figure('Name', sprintf('Spatial %s NL', metric_label), ...
               'Color','w','Position',[100 100 1600 900]);
        tiledlayout('flow','TileSpacing','compact','Padding','compact');
        for ipair = 1:n_pairs
            nexttile; hold on;
            gr = group_results(ipair);
            for iep = 1:3
                mat   = gr.([metric_real ep_keys{iep}]);
                mat_s = gr.([metric_shuf ep_keys{iep}]);
                plot_spatial_curve(pos_centers, mat,   ~is_learner, col_for_epoch(iep,'NL'), ep_names{iep}, '-');
                plot_spatial_curve(pos_centers, mat_s, ~is_learner, col_for_epoch(iep,'NL'), '',            ':');
            end
            for il = 1:size(landmarks_cm,1)
                patch([landmarks_cm(il,1) landmarks_cm(il,2) landmarks_cm(il,2) landmarks_cm(il,1)], ...
                      [-1 -1 1 1], [0.95 0.92 0.88], 'EdgeColor','none', 'FaceAlpha',0.3, 'HandleVisibility','off');
            end
            xline(rz_entry_cm, 'k--', 'HandleVisibility','off');
            if which_metric == 1
                ylim([0 1]);
            else
                yline(0,'k:','HandleVisibility','off');
                ylim([-1 1]);
            end
            title(sprintf('%s (Non-Learners)', gr.pair_name), 'Interpreter','none');
            xlabel('Position (cm)'); ylabel(metric_label);
        end
        sgtitle(sprintf('%s along position - Non-Learners', metric_label));
        save_fig(gcf, fullfile(data_dir, sprintf('TemporalCCAv3_Spatial_%s_NonLearners.svg', metric_label)));
    end
end

%% C) RZ-ENTRY-ALIGNED TRAJECTORIES (CC1, IFI)
for which_metric = 1:2
    if which_metric == 1
        metric_real = 'rz_cc1_';
        metric_shuf = 'rz_cc1_sh_';
        metric_label = 'CC1';
    else
        metric_real = 'rz_ifi_';
        metric_shuf = 'rz_ifi_sh_';
        metric_label = 'IFI';
    end

    figure('Name', sprintf('RZ-aligned %s', metric_label), ...
           'Color','w','Position',[100 100 1600 900]);
    tiledlayout('flow','TileSpacing','compact','Padding','compact');

    for ipair = 1:n_pairs
        nexttile; hold on;
        gr = group_results(ipair);
        for iep = 1:3
            mat   = gr.([metric_real ep_keys{iep}]);
            mat_s = gr.([metric_shuf ep_keys{iep}]);
            plot_spatial_curve(rz_t_ms, mat,   is_learner,  col_for_epoch(iep,'L'),  ep_names{iep}, '-');
            plot_spatial_curve(rz_t_ms, mat_s, is_learner,  col_for_epoch(iep,'L'),  '',            ':');
        end
        xline(0, 'k--', 'HandleVisibility','off');
        if which_metric == 1
            ylim([0 1]);
        else
            yline(0,'k:','HandleVisibility','off');
            ylim([-1 1]);
        end
        title(sprintf('%s (Learners)', gr.pair_name), 'Interpreter','none');
        xlabel('Time from RZ entry (ms)'); ylabel(metric_label);
    end
    sgtitle(sprintf('%s around RZ entry - Learners', metric_label));
    save_fig(gcf, fullfile(data_dir, sprintf('TemporalCCAv3_RZ_%s_Learners.svg', metric_label)));

    if any(~is_learner)
        figure('Name', sprintf('RZ-aligned %s NL', metric_label), ...
               'Color','w','Position',[100 100 1600 900]);
        tiledlayout('flow','TileSpacing','compact','Padding','compact');
        for ipair = 1:n_pairs
            nexttile; hold on;
            gr = group_results(ipair);
            for iep = 1:3
                mat   = gr.([metric_real ep_keys{iep}]);
                mat_s = gr.([metric_shuf ep_keys{iep}]);
                plot_spatial_curve(rz_t_ms, mat,   ~is_learner, col_for_epoch(iep,'NL'), ep_names{iep}, '-');
                plot_spatial_curve(rz_t_ms, mat_s, ~is_learner, col_for_epoch(iep,'NL'), '',            ':');
            end
            xline(0, 'k--', 'HandleVisibility','off');
            if which_metric == 1
                ylim([0 1]);
            else
                yline(0,'k:','HandleVisibility','off');
                ylim([-1 1]);
            end
            title(sprintf('%s (Non-Learners)', gr.pair_name), 'Interpreter','none');
            xlabel('Time from RZ entry (ms)'); ylabel(metric_label);
        end
        sgtitle(sprintf('%s around RZ entry - Non-Learners', metric_label));
        save_fig(gcf, fullfile(data_dir, sprintf('TemporalCCAv3_RZ_%s_NonLearners.svg', metric_label)));
    end
end

%% D) DIAGNOSTIC / SANITY PANEL
figure('Name','Diagnostic','Color','w','Position',[100 100 1600 900]);
tiledlayout(2, 3, 'TileSpacing','compact','Padding','compact');

% --- D1. Fraction of finite animals per pair x epoch (CC1 trial-mean) ---
nexttile;
finite_frac_cc = nan(n_pairs, 3);
for ipair = 1:n_pairs
    gr = group_results(ipair);
    for iep = 1:3
        v = gr.(['trial_mean_cc1_' ep_keys{iep}]);
        finite_frac_cc(ipair, iep) = mean(isfinite(v));
    end
end
imagesc(finite_frac_cc, [0 1]); colorbar;
yticks(1:n_pairs); yticklabels({group_results.pair_name});
xticks(1:3); xticklabels(ep_names);
title('Frac. animals with finite CC1 trial-mean');

% --- D2. Same for IFI ---
nexttile;
finite_frac_ifi = nan(n_pairs, 3);
for ipair = 1:n_pairs
    gr = group_results(ipair);
    for iep = 1:3
        v = gr.(['trial_mean_ifi_' ep_keys{iep}]);
        finite_frac_ifi(ipair, iep) = mean(isfinite(v));
    end
end
imagesc(finite_frac_ifi, [0 1]); colorbar;
yticks(1:n_pairs); yticklabels({group_results.pair_name});
xticks(1:3); xticklabels(ep_names);
title('Frac. animals with finite IFI trial-mean');

% --- D3. n animals contributing per pair (any epoch) ---
nexttile;
n_contrib_L  = zeros(n_pairs, 1);
n_contrib_NL = zeros(n_pairs, 1);
for ipair = 1:n_pairs
    gr = group_results(ipair);
    has_data = false(n_animals, 1);
    for iep = 1:3
        v = gr.(['trial_mean_cc1_' ep_keys{iep}]);
        has_data = has_data | isfinite(v(:));
    end
    n_contrib_L(ipair)  = sum(has_data &  is_learner);
    n_contrib_NL(ipair) = sum(has_data & ~is_learner);
end
bar([n_contrib_L, n_contrib_NL], 'stacked');
xticks(1:n_pairs); xticklabels({group_results.pair_name}); xtickangle(45);
ylabel('# animals'); legend({'Learners','Non-Learners'}, 'Location','best');
title('Animals contributing per pair (any epoch)');

% --- D4. Real vs shuffle scatter for CC1 trial-mean (post epoch) ---
nexttile;
r_all = []; s_all = []; pair_id = [];
for ipair = 1:n_pairs
    r = group_results(ipair).trial_mean_cc1_post(:);
    s = group_results(ipair).trial_mean_cc1_sh_post(:);
    fin = isfinite(r) & isfinite(s);
    r_all = [r_all; r(fin)]; %#ok<AGROW>
    s_all = [s_all; s(fin)]; %#ok<AGROW>
    pair_id = [pair_id; repmat(ipair, sum(fin), 1)]; %#ok<AGROW>
end
if ~isempty(r_all)
    scatter(s_all, r_all, 30, pair_id, 'filled');
    lim = [min([r_all; s_all])-0.05, max([r_all; s_all])+0.05];
    plot(lim, lim, 'k--', 'HandleVisibility','off');
    xlim(lim); ylim(lim);
    xlabel('CC1 shuffle (post)'); ylabel('CC1 real (post)');
    title('Real vs shuffle CC1 (post epoch, color=pair)');
    cb = colorbar; cb.Ticks = 1:n_pairs; cb.TickLabels = {group_results.pair_name};
else
    title('No finite real-shuffle pairs (post)');
end

% --- D5. Same for IFI ---
nexttile;
r_all = []; s_all = []; pair_id = [];
for ipair = 1:n_pairs
    r = group_results(ipair).trial_mean_ifi_post(:);
    s = group_results(ipair).trial_mean_ifi_sh_post(:);
    fin = isfinite(r) & isfinite(s);
    r_all = [r_all; r(fin)]; %#ok<AGROW>
    s_all = [s_all; s(fin)]; %#ok<AGROW>
    pair_id = [pair_id; repmat(ipair, sum(fin), 1)]; %#ok<AGROW>
end
if ~isempty(r_all)
    scatter(s_all, r_all, 30, pair_id, 'filled');
    lim = [min([r_all; s_all])-0.05, max([r_all; s_all])+0.05];
    plot(lim, lim, 'k--', 'HandleVisibility','off');
    xlim(lim); ylim(lim);
    xlabel('IFI shuffle (post)'); ylabel('IFI real (post)');
    title('Real vs shuffle IFI (post epoch, color=pair)');
    cb = colorbar; cb.Ticks = 1:n_pairs; cb.TickLabels = {group_results.pair_name};
else
    title('No finite real-shuffle pairs (post)');
end

% --- D6. Per-pair: LP per animal, with whether each animal contributed ---
nexttile;
contrib_per_animal = false(n_animals, n_pairs);
for ipair = 1:n_pairs
    v = group_results(ipair).trial_mean_cc1_post(:);
    contrib_per_animal(:, ipair) = isfinite(v);
end
imagesc(contrib_per_animal'); colormap(gca, [0.95 0.95 0.95; 0.2 0.6 0.3]);
yticks(1:n_pairs); yticklabels({group_results.pair_name});
xlabel('Animal #'); xticks(1:n_animals);
title('Contribution matrix (green = animal x pair has data)');

if isfield(S, 'shuffle_type')
    sgtitle(sprintf('Diagnostic panel  (shuffle: %s, n=%d, k cap=%d, bin=%d ms, window=+/-%d bins)', ...
        S.shuffle_type, S.n_shuffles, S.max_k_per_region, target_bin_ms, half_window_bins));
else
    sgtitle('Diagnostic panel');
end
save_fig(gcf, fullfile(data_dir, 'TemporalCCAv3_Diagnostic.svg'));

fprintf('\nAll figures saved to %s\n', data_dir);

%% =========================================================================
%% LOCAL HELPERS
%% =========================================================================
function plot_bars_pair(re, rp, rx, se, sp, sx, is_learner, col_L, col_NL, col_shuf)
% Three-epoch errorbar lines for learners and non-learners, with the
% shuffle null overlaid as dashed.
    re = re(:); rp = rp(:); rx = rx(:);
    se = se(:); sp = sp(:); sx = sx(:);

    L = logical(is_learner(:));
    means_L_real = [mean(re(L),'omitnan'), mean(rp(L),'omitnan'), mean(rx(L),'omitnan')];
    sem_L_real   = [sem_(re(L)), sem_(rp(L)), sem_(rx(L))];
    means_L_shuf = [mean(se(L),'omitnan'), mean(sp(L),'omitnan'), mean(sx(L),'omitnan')];
    sem_L_shuf   = [sem_(se(L)), sem_(sp(L)), sem_(sx(L))];

    errorbar(1:3, means_L_real, sem_L_real, 'o-',  'Color', col_L, 'LineWidth', 1.5, 'CapSize', 6);
    errorbar(1:3, means_L_shuf, sem_L_shuf, 'o--', 'Color', col_L, 'LineWidth', 1.0, 'CapSize', 4, 'MarkerSize', 4);

    if any(~L)
        means_NL_real = [mean(re(~L),'omitnan'), mean(rp(~L),'omitnan'), mean(rx(~L),'omitnan')];
        sem_NL_real   = [sem_(re(~L)), sem_(rp(~L)), sem_(rx(~L))];
        means_NL_shuf = [mean(se(~L),'omitnan'), mean(sp(~L),'omitnan'), mean(sx(~L),'omitnan')];
        sem_NL_shuf   = [sem_(se(~L)), sem_(sp(~L)), sem_(sx(~L))];
        errorbar(1:3, means_NL_real, sem_NL_real, 's-',  'Color', col_NL, 'LineWidth', 1.5, 'CapSize', 6);
        errorbar(1:3, means_NL_shuf, sem_NL_shuf, 's--', 'Color', col_NL, 'LineWidth', 1.0, 'CapSize', 4, 'MarkerSize', 4);
    end
    grid on;
end

function plot_spatial_curve(x, mat, mask, col, name, style)
% Mean +/- SEM line for the subset of animals defined by mask.
% mat: [n_animals x n_x]
    mask = logical(mask(:));
    if size(mat,1) ~= numel(mask)
        % MATLAB save / h5py read may transpose; try the other orientation.
        if size(mat,2) == numel(mask), mat = mat'; end
    end
    sub = mat(mask, :);
    if isempty(sub) || all(~isfinite(sub(:))), return; end
    n  = sum(isfinite(sub), 1);
    mu = mean(sub, 1, 'omitnan');
    se = std(sub, 0, 1, 'omitnan') ./ sqrt(max(n,1));
    keep = isfinite(mu);
    if ~any(keep), return; end
    if isempty(name)
        plot(x(keep), mu(keep), style, 'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    else
        fill([x(keep), fliplr(x(keep))], [mu(keep)+se(keep), fliplr(mu(keep)-se(keep))], ...
             col, 'EdgeColor','none','FaceAlpha',0.15, 'HandleVisibility','off');
        plot(x(keep), mu(keep), style, 'Color', col, 'LineWidth', 1.5, 'DisplayName', name);
    end
end

function s = sem_(v)
    v = v(isfinite(v));
    if isempty(v), s = NaN; else, s = std(v) / sqrt(numel(v)); end
end

function c = col_for_epoch(iep, group)
% Three shades for early/pre/post in either learner or non-learner palette.
    if strcmpi(group, 'L')
        cmap = [0.65 0.78 1.00; 0.30 0.55 0.95; 0.05 0.25 0.60];
    else
        cmap = [0.95 0.78 0.65; 0.85 0.50 0.25; 0.55 0.20 0.05];
    end
    c = cmap(iep, :);
end

function save_fig(h, path)
    try
        % Vector SVG export.
        set(h, 'Renderer', 'painters');
        saveas(h, path);
        fprintf('  saved %s\n', path);
    catch ME
        warning('Failed to save %s: %s', path, ME.message);
    end
end
