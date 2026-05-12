%% compare_marginals.m
% Marginal distributions per learner / non-learner group on each metric
% axis (Figure 5).
%
% Why this matters
%   The paired scatter answers "do the two methods agree?" — that is
%   largely a property of the data + pipeline and shouldn't depend much
%   on learning state. Group differences live in the MARGINAL
%   distributions on each axis: does temporal CC1 shift between learners
%   and non-learners on this pair x epoch? Does spatial CC1 shift?
%   Those can be very different stories from the same scatter.
%
% Layout
%   One row per pair, four columns per epoch:
%     col 1: violin / strip of temporal CC1, learner vs non-learner
%     col 2: spatial CC1
%     col 3: temporal IFI
%     col 4: spatial IFI
%   3 epochs x 4 metrics = 12 columns. Or split into 3 figures (one per
%   epoch) for legibility.
%
%   Default: one figure per epoch (saner aspect ratio).
%
% Stats
%   Two-sample Kolmogorov-Smirnov test (kstest2) for distributional
%   shift. Mann-Whitney U (ranksum) for median shift. Both reported per
%   panel; whichever is smaller (more conservative) governs the * mark.

clear; clc; close all;

base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning';
data_dir = fullfile(base_dir, 'HC_V1_data');
fig_dir  = fullfile(base_dir, 'HC_V1_figures', 'CCA_Compare');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

S = load(fullfile(data_dir, 'compare_summary.mat'));
ct = S.compare_table;

epochs = {'early','pre','post'};
pairs  = unique(ct.pair, 'stable');
n_pairs = numel(pairs);

metrics = {'temporal_cc1', 'spatial_cc1', 'temporal_ifi', 'spatial_ifi'};
metric_labels = {'Temporal CC1', 'Spatial CC1', 'Temporal IFI', 'Spatial IFI'};
n_metrics = numel(metrics);

col_L  = [0.1 0.3 0.7];
col_NL = [0.7 0.2 0.1];

for ie = 1:numel(epochs)
    ep = epochs{ie};
    sel_ep = strcmp(ct.epoch, ep);

    fig = figure('Name', sprintf('Marginals — %s', ep), 'Color', 'w', ...
                 'Position', [100 100 200 + 220*n_metrics, 200 + 180*n_pairs]);
    tl = tiledlayout(n_pairs, n_metrics, ...
                     'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Marginal distributions — %s epoch', ep));

    for ipair = 1:n_pairs
        name = pairs{ipair};
        sel = sel_ep & strcmp(ct.pair, name);
        ll = logical(ct.is_learner(sel));

        for im = 1:n_metrics
            ax = nexttile;
            hold(ax, 'on');
            metric = metrics{im};
            v = ct.(metric)(sel);

            v_L  = v(ll  & isfinite(v));
            v_NL = v(~ll & isfinite(v));

            % Strip + box: jittered points + summary box.
            xj_L  = 1 + 0.10*randn(numel(v_L),  1);
            xj_NL = 2 + 0.10*randn(numel(v_NL), 1);
            scatter(ax, xj_L,  v_L,  10, col_L,  'filled', ...
                'MarkerFaceAlpha', 0.5);
            scatter(ax, xj_NL, v_NL, 10, col_NL, 'filled', ...
                'MarkerFaceAlpha', 0.5);

            % Median lines
            if ~isempty(v_L),  plot(ax, [0.7 1.3], [median(v_L)  median(v_L)],  '-', 'Color', col_L,  'LineWidth', 2); end
            if ~isempty(v_NL), plot(ax, [1.7 2.3], [median(v_NL) median(v_NL)], '-', 'Color', col_NL, 'LineWidth', 2); end

            % Stats: KS + ranksum
            p_ks = NaN; p_rs = NaN;
            if numel(v_L) >= 3 && numel(v_NL) >= 3
                try, [~, p_ks] = kstest2(v_L, v_NL); catch, end
                try, p_rs = ranksum(v_L, v_NL);      catch, end
            end
            p_min = min(p_ks, p_rs);
            sig = '';
            if ~isnan(p_min)
                if     p_min < 0.001, sig = '***';
                elseif p_min < 0.01,  sig = '**';
                elseif p_min < 0.05,  sig = '*';
                end
            end

            xticks(ax, [1 2]);
            xticklabels(ax, {sprintf('L (%d)', numel(v_L)), sprintf('NL (%d)', numel(v_NL))});
            xlim(ax, [0.4 2.6]);
            ylabel(ax, metric_labels{im}, 'FontSize', 8);

            ttl = sprintf('%s   p_{KS}=%.3f  p_{RS}=%.3f %s', ...
                name, p_ks, p_rs, sig);
            title(ax, ttl, 'FontSize', 7, 'FontWeight', 'normal');
            grid(ax, 'on');
        end
    end

    saveas(fig, fullfile(fig_dir, sprintf('TemporalVsSpatial_marginals_%s.svg', ep)));
end
fprintf('\nMarginal-distribution figures saved to %s\n', fig_dir);
