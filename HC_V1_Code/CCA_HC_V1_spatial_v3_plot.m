%% CCA_HC_V1_spatial_v3_plot.m
% Diagnostic and group-summary plots for the spatial v5 CCA output.
%
% Loads the most-recent `Spatial_CCA_v3_<date>.mat` from HC_V1_data/ and
% produces a suite of SVG figures (one per topic), all saved into
% HC_V1_data/ with the `SpatialCCAv3_` prefix to match legacy plot-script
% conventions.

clear; clc; close all;

%% 1. LOAD LATEST RESULTS

base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_dir = fullfile(base_dir, 'HC_V1_data');

result_files = dir(fullfile(data_dir, 'Spatial_CCA_v3_*.mat'));
if isempty(result_files)
    error('No Spatial_CCA_v3_*.mat found in %s', data_dir);
end
[~, latest_idx] = max([result_files.datenum]);
load_path = fullfile(data_dir, result_files(latest_idx).name);
fprintf('Loading: %s\n', load_path);
loaded = load(load_path);
group_results  = loaded.group_results;
cfg            = loaded.cfg;
is_learner     = logical(loaded.is_learner(:));
analysis_lp    = loaded.analysis_lp;       %#ok<NASGU>

current_date = datestr(now, 'yyyy_mm_dd');
epoch_names = {'naive', 'pre-LP', 'post-LP'};
n_pairs   = length(group_results);
n_animals = size(group_results(1).per_epoch, 1);
n_epochs  = size(group_results(1).per_epoch, 2);

pair_names = cell(1, n_pairs);
for ip = 1:n_pairs, pair_names{ip} = group_results(ip).pair_name; end

%% 2. EXTRACT TIDY ARRAYS  (n_animals × n_pairs × n_epochs)

[cc_full, cc_cv, cc_excess, ifi_lagged, ifi_proj_mean, ...
 null_mean, null_thr, k_used, samples_per_var, ...
 var_x, var_y, prin_angle_A_top, prin_angle_B_top, ...
 prin_angle_A_mean, prin_angle_B_mean] = deal(nan(n_animals, n_pairs, n_epochs));

% IFI per-trial pooled across animals × trials per (pair × epoch)
ifi_proj_all = cell(n_pairs, n_epochs);

for ip = 1:n_pairs
    for ia = 1:n_animals
        for ie = 1:n_epochs
            r = group_results(ip).per_epoch(ia, ie);
            if isempty(r.cc_cv), continue; end
            cc_full(ia, ip, ie)        = r.cc_full(1);
            cc_cv(ia, ip, ie)          = r.cc_cv(1);
            cc_excess(ia, ip, ie)      = r.cc_excess;
            ifi_lagged(ia, ip, ie)     = r.ifi_lagged;
            ifi_proj_mean(ia, ip, ie)  = mean(r.ifi_proj_per_trial, 'omitnan');
            null_mean(ia, ip, ie)      = r.null_mean;
            null_thr(ia, ip, ie)       = r.null_threshold;
            k_used(ia, ip, ie)         = r.k_used;
            samples_per_var(ia, ip, ie) = r.samples_per_var;
            var_x(ia, ip, ie)          = r.var_x_at_k;
            var_y(ia, ip, ie)          = r.var_y_at_k;
            if ~isempty(r.principal_angles_A) && ~all(isnan(r.principal_angles_A))
                prin_angle_A_top(ia, ip, ie)  = r.principal_angles_A(1);
                prin_angle_A_mean(ia, ip, ie) = mean(r.principal_angles_A, 'omitnan');
            end
            if ~isempty(r.principal_angles_B) && ~all(isnan(r.principal_angles_B))
                prin_angle_B_top(ia, ip, ie)  = r.principal_angles_B(1);
                prin_angle_B_mean(ia, ip, ie) = mean(r.principal_angles_B, 'omitnan');
            end
            % Pool per-trial IFIs
            ifi_proj_all{ip, ie} = [ifi_proj_all{ip, ie}; r.ifi_proj_per_trial(:)];
        end
    end
end

%% 3. PLOT — per-pair epoch summaries (learner split)

plot_metric_per_pair(cc_cv, pair_names, epoch_names, is_learner, ...
    'CC1 (cross-validated)', sprintf('Spatial CCA v3 — CC1 cv (%s)', current_date), ...
    fullfile(data_dir, sprintf('SpatialCCAv3_CC1_cv_%s.svg', current_date)));

plot_metric_per_pair(cc_excess, pair_names, epoch_names, is_learner, ...
    'CC1 excess (real − shuffle mean)', sprintf('Spatial CCA v3 — CC1 excess (%s)', current_date), ...
    fullfile(data_dir, sprintf('SpatialCCAv3_CC1_excess_%s.svg', current_date)));

plot_metric_per_pair(ifi_lagged, pair_names, epoch_names, is_learner, ...
    'IFI (lagged-refit, position bins)', ...
    sprintf('Spatial CCA v3 — IFI lagged: positive ⇒ X leads Y (%s)', current_date), ...
    fullfile(data_dir, sprintf('SpatialCCAv3_IFI_lagged_%s.svg', current_date)));

plot_metric_per_pair(ifi_proj_mean, pair_names, epoch_names, is_learner, ...
    'IFI (projection, per-trial mean)', ...
    sprintf('Spatial CCA v3 — IFI projection (per-trial mean) (%s)', current_date), ...
    fullfile(data_dir, sprintf('SpatialCCAv3_IFI_proj_%s.svg', current_date)));

%% 4. PLOT — principal-angle subspace alignment

% Principal angles are stored at epoch ie>=2 (relative to ie-1).
% Plot two transitions: naive→pre-LP (ie=2) and pre-LP→post-LP (ie=3).
plot_principal_angles( ...
    prin_angle_A_top(:, :, 2:end), prin_angle_B_top(:, :, 2:end), ...
    prin_angle_A_mean(:, :, 2:end), prin_angle_B_mean(:, :, 2:end), ...
    pair_names, {'naive → pre-LP', 'pre-LP → post-LP'}, is_learner, ...
    fullfile(data_dir, sprintf('SpatialCCAv3_PrincipalAngles_%s.svg', current_date)));

%% 5. PLOT — lagged-refit CC curves per pair × epoch

plot_lagged_cc_curves(group_results, cfg, ...
    fullfile(data_dir, sprintf('SpatialCCAv3_LaggedCCcurves_%s.svg', current_date)));

%% 6. PLOT — per-trial IFI distributions pooled across animals

plot_ifi_proj_distributions(ifi_proj_all, pair_names, epoch_names, ...
    fullfile(data_dir, sprintf('SpatialCCAv3_IFI_proj_distributions_%s.svg', current_date)));

%% 7. PLOT — diagnostics (k_used, samples/var, variance explained)

plot_diagnostics(k_used, samples_per_var, var_x, var_y, pair_names, ...
    fullfile(data_dir, sprintf('SpatialCCAv3_Diagnostics_%s.svg', current_date)));

%% 8. PLOT — significance summary (cc_cv vs null threshold)

plot_significance_summary(cc_cv, null_thr, null_mean, pair_names, epoch_names, ...
    is_learner, ...
    fullfile(data_dir, sprintf('SpatialCCAv3_Significance_%s.svg', current_date)));

fprintf('\nAll plots saved to %s\n', data_dir);


%% Plot functions

function plot_metric_per_pair(values, pair_names, epoch_names, is_learner, ...
                              ylabel_text, title_text, save_path)
% Bar plot per area-pair: epochs on x-axis, learner / non-learner split.
% values: (n_animals × n_pairs × n_epochs)
    n_pairs = numel(pair_names);
    n_epochs = numel(epoch_names);
    [n_rows, n_cols] = grid_for(n_pairs);

    fig = figure('Color', 'w', 'Position', [50 50 220 * n_cols 200 * n_rows]);
    set(fig, 'PaperPositionMode', 'auto');

    learner_color    = [0.20 0.40 0.80];
    nonlearner_color = [0.80 0.35 0.20];

    for ip = 1:n_pairs
        subplot(n_rows, n_cols, ip);
        hold on;
        x_centers = 1:n_epochs;
        bar_w = 0.35;
        for ie = 1:n_epochs
            v_learn    = values(is_learner,  ip, ie);
            v_nonlearn = values(~is_learner, ip, ie);
            m_l = mean(v_learn,    'omitnan');
            m_n = mean(v_nonlearn, 'omitnan');
            s_l = sem(v_learn);
            s_n = sem(v_nonlearn);
            bar(ie - bar_w/2, m_l, bar_w, 'FaceColor', learner_color, ...
                'EdgeColor', 'none');
            bar(ie + bar_w/2, m_n, bar_w, 'FaceColor', nonlearner_color, ...
                'EdgeColor', 'none');
            errorbar(ie - bar_w/2, m_l, s_l, 'Color', 'k', 'LineWidth', 0.5, ...
                'CapSize', 0);
            errorbar(ie + bar_w/2, m_n, s_n, 'Color', 'k', 'LineWidth', 0.5, ...
                'CapSize', 0);
            % per-animal points
            scatter(repmat(ie - bar_w/2, sum(~isnan(v_learn)), 1), ...
                v_learn(~isnan(v_learn)), 10, 'k', 'filled', ...
                'MarkerFaceAlpha', 0.35);
            scatter(repmat(ie + bar_w/2, sum(~isnan(v_nonlearn)), 1), ...
                v_nonlearn(~isnan(v_nonlearn)), 10, 'k', 'filled', ...
                'MarkerFaceAlpha', 0.35);
        end
        xlim([0.5, n_epochs + 0.5]);
        xticks(x_centers);
        xticklabels(epoch_names);
        ylabel(ylabel_text, 'Interpreter', 'tex');
        title(pair_names{ip}, 'FontWeight', 'normal');
        if ip == 1
            legend({'learner', 'non-learner'}, ...
                'Location', 'best', 'Box', 'off', 'FontSize', 7);
        end
        grid on; box off;
    end
    sgtitle(title_text, 'Interpreter', 'none');
    saveas(fig, save_path);
    fprintf('  saved %s\n', save_path);
end


function plot_principal_angles(angle_A_top, angle_B_top, angle_A_mean, angle_B_mean, ...
                                pair_names, transition_names, is_learner, save_path)
% Two-panel layout: angle_A and angle_B side by side; bars are per-pair
% mean across animals, with the top angle (most aligned direction) and the
% mean angle plotted.
    n_pairs = numel(pair_names);
    n_trans = numel(transition_names);
    fig = figure('Color', 'w', 'Position', [80 80 1100 360]);
    set(fig, 'PaperPositionMode', 'auto');

    learner_color    = [0.20 0.40 0.80];
    nonlearner_color = [0.80 0.35 0.20];

    for side = 1:2
        if side == 1
            top  = angle_A_top;
            avg  = angle_A_mean;     %#ok<NASGU>
            ttl  = 'Principal angles — X side (A)';
        else
            top  = angle_B_top;
            avg  = angle_B_mean;     %#ok<NASGU>
            ttl  = 'Principal angles — Y side (B)';
        end
        subplot(1, 2, side); hold on;
        for it = 1:n_trans
            for ip = 1:n_pairs
                v_l = top(is_learner,  ip, it);
                v_n = top(~is_learner, ip, it);
                x_l = (ip - 1) * (n_trans + 1) + it;
                bar(x_l, mean(v_l, 'omitnan'), 0.3, ...
                    'FaceColor', learner_color, 'EdgeColor', 'none');
                bar(x_l + 0.3, mean(v_n, 'omitnan'), 0.3, ...
                    'FaceColor', nonlearner_color, 'EdgeColor', 'none');
            end
        end
        yline(pi/2, 'k--', 'orthogonal (π/2)', 'FontSize', 7, 'LineWidth', 0.5);
        yline(0,    'k:',  'identical (0)',    'FontSize', 7, 'LineWidth', 0.5);
        xticks((1:n_pairs) * (n_trans + 1) - n_trans/2 - 0.5);
        xticklabels(pair_names);
        xtickangle(30);
        ylabel('principal angle (rad)');
        title(ttl, 'FontWeight', 'normal');
        ylim([0, pi/2 + 0.05]);
        if side == 1
            legend({'learner', 'non-learner'}, 'Location', 'best', 'Box', 'off');
        end
        grid on; box off;
    end
    sgtitle(sprintf('Subspace alignment across epoch transitions (top angle, %s)', ...
        strjoin(transition_names, ', ')), 'Interpreter', 'none');
    saveas(fig, save_path);
    fprintf('  saved %s\n', save_path);
end


function plot_lagged_cc_curves(group_results, cfg, save_path)
% Per pair: mean (across animals) cc_xy and cc_yx vs lag, one line per epoch.
    n_pairs = numel(group_results);
    n_epochs = size(group_results(1).per_epoch, 2);
    max_lag = cfg.max_lag_bins;
    lag_axis_bins = 0:max_lag;
    lag_axis_cm = lag_axis_bins * cfg.bin_size_cm;

    [n_rows, n_cols] = grid_for(n_pairs);
    fig = figure('Color', 'w', 'Position', [80 80 240 * n_cols 220 * n_rows]);
    set(fig, 'PaperPositionMode', 'auto');

    cmap = [0.20 0.50 0.80;  % naive
            0.85 0.55 0.20;  % pre-LP
            0.40 0.70 0.30]; % post-LP
    epoch_labels = {'naive', 'pre-LP', 'post-LP'};

    for ip = 1:n_pairs
        subplot(n_rows, n_cols, ip); hold on;
        for ie = 1:n_epochs
            n_animals_ep = size(group_results(ip).per_epoch, 1);
            xy = nan(n_animals_ep, max_lag + 1);
            yx = nan(n_animals_ep, max_lag + 1);
            for ia = 1:n_animals_ep
                r = group_results(ip).per_epoch(ia, ie);
                if isempty(r.ifi_lagged_curve_xy_smooth), continue; end
                xy(ia, :) = r.ifi_lagged_curve_xy_smooth(:).';
                yx(ia, :) = r.ifi_lagged_curve_yx_smooth(:).';
            end
            xy_mean = mean(xy, 1, 'omitnan');
            yx_mean = mean(yx, 1, 'omitnan');
            plot(lag_axis_cm, xy_mean, '-',  'Color', cmap(ie, :), 'LineWidth', 1.2);
            plot(lag_axis_cm, yx_mean, '--', 'Color', cmap(ie, :), 'LineWidth', 1.2);
        end
        xlabel('|spatial lag| (cm)');
        ylabel('CC1');
        title(group_results(ip).pair_name, 'FontWeight', 'normal');
        if ip == 1
            % solid = X-leads-Y, dashed = Y-leads-X
            h1 = plot(nan, nan, '-k', 'LineWidth', 1.2);
            h2 = plot(nan, nan, '--k', 'LineWidth', 1.2);
            l = legend([h1, h2], {'X leads Y', 'Y leads X'}, ...
                'Box', 'off', 'Location', 'best', 'FontSize', 7);
            l.AutoUpdate = 'off';
        end
        grid on; box off;
    end
    sgtitle('Lagged-refit CC curves (mean across animals, by epoch)', ...
        'Interpreter', 'none');
    saveas(fig, save_path);
    fprintf('  saved %s\n', save_path);
end


function plot_ifi_proj_distributions(ifi_pooled, pair_names, epoch_names, save_path)
% Per pair: violin / kernel-density of pooled per-trial IFI per epoch.
    n_pairs = numel(pair_names);
    n_epochs = numel(epoch_names);
    [n_rows, n_cols] = grid_for(n_pairs);
    fig = figure('Color', 'w', 'Position', [80 80 240 * n_cols 200 * n_rows]);
    set(fig, 'PaperPositionMode', 'auto');

    cmap = [0.20 0.50 0.80; 0.85 0.55 0.20; 0.40 0.70 0.30];

    for ip = 1:n_pairs
        subplot(n_rows, n_cols, ip); hold on;
        for ie = 1:n_epochs
            v = ifi_pooled{ip, ie};
            v = v(~isnan(v));
            if isempty(v), continue; end
            % Plot histogram outlines.
            edges = -1:0.1:1;
            [counts, ~] = histcounts(v, edges);
            centers = (edges(1:end-1) + edges(2:end)) / 2;
            plot(centers, counts / sum(counts), '-', ...
                'Color', cmap(ie, :), 'LineWidth', 1.2);
        end
        xline(0, 'k:');
        xlabel('IFI (projection, per trial)');
        ylabel('relative frequency');
        title(pair_names{ip}, 'FontWeight', 'normal');
        xlim([-1, 1]);
        if ip == 1
            legend(epoch_names, 'Box', 'off', 'Location', 'best', 'FontSize', 7);
        end
        grid on; box off;
    end
    sgtitle('Per-trial projection IFI distributions (pooled across animals)', ...
        'Interpreter', 'none');
    saveas(fig, save_path);
    fprintf('  saved %s\n', save_path);
end


function plot_diagnostics(k_used, samples_per_var, var_x, var_y, pair_names, save_path)
% Four-panel diagnostics: k_used histogram, samples/var histogram, var_x and
% var_y scatter vs k.
    fig = figure('Color', 'w', 'Position', [80 80 1100 700]);
    set(fig, 'PaperPositionMode', 'auto');

    k_flat = k_used(:);    k_flat = k_flat(~isnan(k_flat));
    spv_flat = samples_per_var(:);   spv_flat = spv_flat(~isnan(spv_flat));
    vx_flat = var_x(:);    vx_flat_valid = ~isnan(vx_flat);
    vy_flat = var_y(:);    vy_flat_valid = ~isnan(vy_flat);
    k_for_var = k_used(:);

    subplot(2, 2, 1);
    histogram(k_flat, 'BinMethod', 'integers', 'FaceColor', [0.5 0.5 0.5]);
    xlabel('k_{used} per side');
    ylabel('count of (animal × pair × epoch)');
    title('Distribution of k', 'FontWeight', 'normal');
    grid on; box off;

    subplot(2, 2, 2);
    histogram(spv_flat, 30, 'FaceColor', [0.5 0.5 0.5]);
    xline(50, 'r--', 'H&H stability threshold (50)', 'FontSize', 8);
    xline(25, 'r:',  'below 25 = sample-limited', 'FontSize', 8);
    xlabel('samples per variable');
    ylabel('count');
    title('Distribution of samples/var', 'FontWeight', 'normal');
    grid on; box off;

    subplot(2, 2, 3);
    scatter(k_for_var(vx_flat_valid), vx_flat(vx_flat_valid), 18, ...
        [0.20 0.40 0.80], 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('k_{used}');
    ylabel('total variance explained, X side');
    title('Variance captured per side (X)', 'FontWeight', 'normal');
    ylim([0, 1.02]);
    grid on; box off;

    subplot(2, 2, 4);
    scatter(k_for_var(vy_flat_valid), vy_flat(vy_flat_valid), 18, ...
        [0.80 0.35 0.20], 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('k_{used}');
    ylabel('total variance explained, Y side');
    title('Variance captured per side (Y)', 'FontWeight', 'normal');
    ylim([0, 1.02]);
    grid on; box off;

    sgtitle('Spatial CCA v3 diagnostics', 'Interpreter', 'none');
    saveas(fig, save_path);
    fprintf('  saved %s\n', save_path);
end


function plot_significance_summary(cc_cv, null_thr, null_mean, pair_names, ...
                                    epoch_names, is_learner, save_path)
% Per pair × epoch: dot plot of (cc_cv − null_thr). Positive = real fit
% clears H&H significance bar. Colored by learner status.
    n_pairs = numel(pair_names);
    n_epochs = numel(epoch_names);
    [n_rows, n_cols] = grid_for(n_pairs);
    fig = figure('Color', 'w', 'Position', [80 80 240 * n_cols 200 * n_rows]);
    set(fig, 'PaperPositionMode', 'auto');

    for ip = 1:n_pairs
        subplot(n_rows, n_cols, ip); hold on;
        for ie = 1:n_epochs
            v_l = cc_cv(is_learner,  ip, ie) - null_thr(is_learner,  ip, ie);
            v_n = cc_cv(~is_learner, ip, ie) - null_thr(~is_learner, ip, ie);
            jitter = 0.07;
            x_l = ie - 0.15 + jitter * (rand(numel(v_l), 1) - 0.5);
            x_n = ie + 0.15 + jitter * (rand(numel(v_n), 1) - 0.5);
            scatter(x_l, v_l, 18, [0.20 0.40 0.80], 'filled', 'MarkerFaceAlpha', 0.6);
            scatter(x_n, v_n, 18, [0.80 0.35 0.20], 'filled', 'MarkerFaceAlpha', 0.6);
        end
        yline(0, 'k--');
        xlim([0.5, n_epochs + 0.5]);
        xticks(1:n_epochs);
        xticklabels(epoch_names);
        ylabel('cc_{cv} − threshold');
        title(pair_names{ip}, 'FontWeight', 'normal');
        if ip == 1
            legend({'learner', 'non-learner'}, 'Box', 'off', ...
                'Location', 'best', 'FontSize', 7);
        end
        grid on; box off;
    end
    sgtitle('Significance: real CC1 minus shuffle threshold (3σ bar)', ...
        'Interpreter', 'none');
    saveas(fig, save_path);
    fprintf('  saved %s\n', save_path);
end


function s = sem(x)
    x = x(~isnan(x));
    if isempty(x)
        s = NaN;
    else
        s = std(x) / sqrt(numel(x));
    end
end


function [n_rows, n_cols] = grid_for(n_panels)
% Pick a reasonable subplot grid for n_panels.
    n_cols = ceil(sqrt(n_panels));
    n_rows = ceil(n_panels / n_cols);
end
