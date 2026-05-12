%% compare_excess_scatter.m
% Shuffle-corrected paired scatter (Figures 3 + 4).
%
% Loads compare_summary.mat (produced by compare_temporal_spatial.m) and
% plots excess CC1 = real - shuffle  for both temporal and spatial. The
% same-shape figures as compare_temporal_spatial.m's Figure 1 + 2 but
% with the small-sample inflation removed.
%
% Why this matters
%   The absolute CC1 from canoncorr is upward-biased when sample size is
%   small relative to k1 + k2, especially when the input matrix is
%   smoothed (spatial pipeline averages firing rates per position bin
%   across trials, dramatically reducing effective n). Real - shuffle
%   absorbs that bias: a method-driven inflation appears equally in real
%   and shuffle, and cancels. Biological coupling does not appear in the
%   shuffle, so the excess preserves it.
%
% IFI excess is also computed (real - shuffle) for symmetry, but is less
% commonly used: IFI is already a normalized asymmetry index, so its
% shuffle is centered near 0 by construction. Subtracting it doesn't add
% much. We plot it anyway in case the per-pair shuffle has structure.

clear; clc; close all;

base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning';
data_dir = fullfile(base_dir, 'HC_V1_data');
fig_dir  = fullfile(base_dir, 'HC_V1_figures', 'CCA_Compare');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

S = load(fullfile(data_dir, 'compare_summary.mat'));
ct = S.compare_table;

% Shuffle-corrected per-trial values
t_cc_x  = ct.temporal_cc1 - ct.temporal_cc1_sh;
s_cc_x  = ct.spatial_cc1  - ct.spatial_cc1_sh;
t_ifi_x = ct.temporal_ifi - ct.temporal_ifi_sh;
s_ifi_x = ct.spatial_ifi  - ct.spatial_ifi_sh;

epochs = {'early','pre','post'};
pairs  = unique(ct.pair, 'stable');
n_pairs = numel(pairs);
cmap   = lines(n_pairs);

%% Figure 3 — Combined excess scatter (3 epochs x {CC1, IFI})
fig3 = figure('Name', 'Excess CC1 / IFI (combined)', 'Color', 'w', ...
              'Position', [100 100 1500 800]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Temporal vs Spatial — shuffle-corrected (real - shuffle)');

metrics = {'cc1', 'ifi'};
metric_labels = {'\Delta CC1 (real - shuffle)', '\Delta IFI (real - shuffle)'};
for im = 1:numel(metrics)
    for ie = 1:numel(epochs)
        ax = nexttile;
        hold(ax, 'on');
        ep = epochs{ie};
        sel_ep = strcmp(ct.epoch, ep);

        for ipair = 1:n_pairs
            sel = sel_ep & strcmp(ct.pair, pairs{ipair});
            if ~any(sel), continue; end
            if im == 1
                xs = s_cc_x(sel);  ys = t_cc_x(sel);
            else
                xs = s_ifi_x(sel); ys = t_ifi_x(sel);
            end
            scatter(ax, xs, ys, 14, cmap(ipair, :), 'filled', ...
                'MarkerFaceAlpha', 0.5, 'DisplayName', pairs{ipair});
        end

        % Identity line through origin
        xl = xlim(ax); yl = ylim(ax);
        lo = min([xl yl]); hi = max([xl yl]);
        plot(ax, [lo hi], [lo hi], 'k--', 'HandleVisibility', 'off');
        plot(ax, [lo hi], [0 0],   'k:',  'HandleVisibility', 'off');
        plot(ax, [0 0],   [lo hi], 'k:',  'HandleVisibility', 'off');

        if im == 1
            xs_all = s_cc_x(sel_ep);  ys_all = t_cc_x(sel_ep);
        else
            xs_all = s_ifi_x(sel_ep); ys_all = t_ifi_x(sel_ep);
        end
        ok = isfinite(xs_all) & isfinite(ys_all);
        if sum(ok) >= 3
            rho = corr(xs_all(ok), ys_all(ok), 'type', 'Spearman');
        else
            rho = NaN;
        end
        title(ax, sprintf('%s — %s   \\rho_S = %.2f (n=%d)', ...
            metric_labels{im}, ep, rho, sum(ok)));
        xlabel(ax, sprintf('Spatial %s', metric_labels{im}));
        ylabel(ax, sprintf('Temporal %s', metric_labels{im}));
        grid(ax, 'on'); axis(ax, 'square');

        if im == 1 && ie == numel(epochs)
            legend(ax, 'Location', 'eastoutside', 'NumColumns', 1, ...
                'FontSize', 7);
        end
    end
end
saveas(fig3, fullfile(fig_dir, 'TemporalVsSpatial_excess_combined.svg'));

%% Figure 4 — Per-pair excess scatter (n_pairs x 6)
fig4 = figure('Name', 'Excess CC1 / IFI (per pair)', 'Color', 'w', ...
              'Position', [100 100 1700 200 + 200 * n_pairs]);
tl4 = tiledlayout(n_pairs, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl4, 'Temporal vs Spatial (excess) — per pair, per epoch');

col_L  = [0.1 0.3 0.7];
col_NL = [0.7 0.2 0.1];
for ipair = 1:n_pairs
    name = pairs{ipair};
    sel_pair = strcmp(ct.pair, name);
    for im = 1:numel(metrics)
        for ie = 1:numel(epochs)
            ax = nexttile;
            hold(ax, 'on');
            ep = epochs{ie};
            sel = sel_pair & strcmp(ct.epoch, ep);

            if im == 1
                xs = s_cc_x(sel);  ys = t_cc_x(sel);
            else
                xs = s_ifi_x(sel); ys = t_ifi_x(sel);
            end
            ll = logical(ct.is_learner(sel));

            scatter(ax, xs(ll),  ys(ll),  18, col_L,  'filled', ...
                'MarkerFaceAlpha', 0.6);
            scatter(ax, xs(~ll), ys(~ll), 18, col_NL, 'filled', ...
                'MarkerFaceAlpha', 0.6);

            xl = xlim(ax); yl = ylim(ax);
            lo = min([xl yl]); hi = max([xl yl]);
            plot(ax, [lo hi], [lo hi], 'k--', 'HandleVisibility', 'off');
            plot(ax, [lo hi], [0 0],   'k:',  'HandleVisibility', 'off');
            plot(ax, [0 0],   [lo hi], 'k:',  'HandleVisibility', 'off');

            ok = isfinite(xs) & isfinite(ys);
            if sum(ok) >= 3
                rho = corr(xs(ok), ys(ok), 'type', 'Spearman');
            else
                rho = NaN;
            end
            title(ax, sprintf('%s — %s — %s   \\rho=%.2f (n=%d)', ...
                name, metric_labels{im}, ep, rho, sum(ok)), ...
                'FontSize', 8);
            grid(ax, 'on'); axis(ax, 'square');
            if ipair == n_pairs
                xlabel(ax, sprintf('Spatial %s', metric_labels{im}));
            end
            if ie == 1 && im == 1
                ylabel(ax, sprintf('%s\nTemporal', name), 'FontWeight', 'bold');
            elseif ie == 1
                ylabel(ax, sprintf('Temporal %s', metric_labels{im}));
            end
        end
    end
end
saveas(fig4, fullfile(fig_dir, 'TemporalVsSpatial_excess_per_pair.svg'));

fprintf('\nFigures 3 + 4 saved to %s\n', fig_dir);
