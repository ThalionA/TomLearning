%% CCA_HC_V1_temporal_v4_plot.m
% Plotting / verification script for v4 temporal CCA results.
%
% Loads the latest Temporal_CCA_v4_<date>.mat and produces:
%   A. Per-trial CC1 bars per pair (Learners vs Non-Learners), early/pre/post,
%      with rmANOVA + 1-sample t-test against the shuffle null on each bar.
%   B. Per-trial IFI bars (same layout, with t-test against zero).
%   C. RZ-pre per-bin CC trace, REAL - SHUFFLE as the primary line (raw real
%      / shuffle on a separate second figure since the absolute level is
%      biased high by the small n vs k1+k2 fit).
%   D. RZ-pre IFI scalar bars (per-trial IFI on the pre-RZ segment, epoch
%      mean per animal), rmANOVA + 1-sample t-test against zero.
%   E. EXCESS per-trial CC (real - shuffle) bars, rmANOVA. This is often
%      more interpretable than raw CC1 in regimes with non-trivial shuffle bias.
%   F. Diagnostic panel: finite-fraction heatmaps, real-vs-shuffle scatters,
%      animal x pair contribution matrix.
%
% All figures are saved as .svg into the same data folder.

%% 1. Locate latest results file
clear; clc; close all;

base_dir       = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_subfolder = 'HC_V1_data';
data_dir       = fullfile(base_dir, data_subfolder);

addpath(fileparts(mfilename('fullpath')));    % so v4_*.m helpers are reachable

files = dir(fullfile(data_dir, 'Temporal_CCA_v4_*.mat'));
if isempty(files)
    error('No Temporal_CCA_v4_*.mat files found in %s', data_dir);
end
[~, idx] = max([files.datenum]);
results_path = fullfile(data_dir, files(idx).name);
fprintf('Loading %s\n', results_path);
S = load(results_path);

group_results  = S.group_results;
is_learner     = logical(S.is_learner(:));
analysis_lp    = S.analysis_lp;
target_bin_ms  = double(S.target_bin_ms);
n_rz_bins      = double(S.n_rz_bins);
lags           = double(S.lags(:)');
n_shuffles     = double(S.n_shuffles);
min_block_bins = double(S.min_block_bins);
max_k          = double(S.max_k_per_region);

n_pairs   = numel(group_results);
n_animals = numel(is_learner);

% RZ x-axis: segment runs (i_rz - n_rz_bins) : (i_rz - 1). The k-th element
% of the segment (k=1..n_rz_bins) sits (k - n_rz_bins - 1) bins before
% RZ entry, so times are -n_rz_bins*bin_ms .. -bin_ms (exclusive of 0).
rz_t_ms = ((1:n_rz_bins) - n_rz_bins - 1) * target_bin_ms;

% --- Style ---
col_L    = [0.85 0.33 0.10];
col_NL   = [0.30 0.30 0.30];
ep_names = {'Early','Pre','Post'};
ep_keys  = {'early','pre','post'};

fprintf('Pairs=%d, animals=%d, learners=%d. RZ window: %g..%g ms (%d bins).\n', ...
    n_pairs, n_animals, sum(is_learner), rz_t_ms(1), rz_t_ms(end), n_rz_bins);

%% A) PER-TRIAL CC1 BARS WITH rmANOVA
fprintf('\n=== Section A: per-trial CC1 bars (rmANOVA) ===\n');
figure('Name','Per-trial CC1 (rmANOVA)','Color','w','Position',[100 100 1500 850]);
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for ipair = 1:n_pairs
    nexttile; hold on;
    gr = group_results(ipair);
    fprintf('\n--- %s : per-trial CC1 ---\n', gr.pair_name);
    plot_bars_with_rmanova( ...
        gr.trial_cc1_early,  gr.trial_cc1_pre,  gr.trial_cc1_post, ...
        gr.trial_cc1_sh_early, gr.trial_cc1_sh_pre, gr.trial_cc1_sh_post, ...
        is_learner, gr.pair_name, 'CC1', col_L, col_NL, 'shuffle');
end
sgtitle('Per-trial CC1 (mean over epoch trials), * = real > shuffle (paired t-test, p<0.05)');
save_fig(gcf, fullfile(data_dir, 'TemporalCCAv4_TrialCC1_Bars.svg'));

%% B) PER-TRIAL IFI BARS WITH rmANOVA
fprintf('\n=== Section B: per-trial IFI bars (rmANOVA) ===\n');
figure('Name','Per-trial IFI (rmANOVA)','Color','w','Position',[100 100 1500 850]);
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for ipair = 1:n_pairs
    nexttile; hold on;
    gr = group_results(ipair);
    fprintf('\n--- %s : per-trial IFI ---\n', gr.pair_name);
    plot_bars_with_rmanova( ...
        gr.trial_ifi_early,  gr.trial_ifi_pre,  gr.trial_ifi_post, ...
        gr.trial_ifi_sh_early, gr.trial_ifi_sh_pre, gr.trial_ifi_sh_post, ...
        is_learner, gr.pair_name, 'IFI', col_L, col_NL, 'zero');
end
sgtitle('Per-trial IFI (block-weighted mean), * = different from 0 (1-sample t-test, p<0.05)');
save_fig(gcf, fullfile(data_dir, 'TemporalCCAv4_TrialIFI_Bars.svg'));

%% C) RZ-PRE PER-BIN CC TRACE
% Primary: real - shuffle. Secondary: raw real and shuffle separately.
fprintf('\n=== Section C: RZ-pre per-bin CC trace ===\n');

% Primary figure: real - shuffle, learners and non-learners overlaid.
figure('Name','RZ-pre CC trace (real - shuffle)','Color','w','Position',[100 100 1600 900]);
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for ipair = 1:n_pairs
    nexttile; hold on;
    gr = group_results(ipair);
    for iep = 1:3
        rl  = gr.(['rz_cc1_'    ep_keys{iep}]);   % [animals x n_rz_bins]
        rls = gr.(['rz_cc1_sh_' ep_keys{iep}]);
        excess = rl - rls;
        plot_trace(rz_t_ms, excess, is_learner, col_for_epoch(iep,'L'), ...
                   sprintf('%s L', ep_names{iep}), '-');
        if any(~is_learner)
            plot_trace(rz_t_ms, excess, ~is_learner, col_for_epoch(iep,'NL'), ...
                       sprintf('%s NL', ep_names{iep}), '--');
        end
    end
    yline(0, 'k:', 'HandleVisibility','off');
    xline(0, 'k--', 'HandleVisibility','off');
    title(gr.pair_name, 'Interpreter','none');
    xlabel('Time before RZ entry (ms)');
    ylabel('CC1 real - shuffle');
    if ipair == 1
        legend('Location','best','Box','off','FontSize',7,'Interpreter','none');
    end
end
sgtitle(sprintf('RZ-pre per-bin CC1 (real - shuffle); n trials per epoch up to %d, k1+k2 up to %d', ...
    10, 2*max_k));
save_fig(gcf, fullfile(data_dir, 'TemporalCCAv4_RZ_CC_excess.svg'));

% Secondary figure: raw real (solid) and shuffle (dotted) for inspection.
figure('Name','RZ-pre CC trace (raw real & shuffle)','Color','w','Position',[100 100 1600 900]);
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for ipair = 1:n_pairs
    nexttile; hold on;
    gr = group_results(ipair);
    for iep = 1:3
        rl  = gr.(['rz_cc1_'    ep_keys{iep}]);
        rls = gr.(['rz_cc1_sh_' ep_keys{iep}]);
        plot_trace(rz_t_ms, rl,  is_learner, col_for_epoch(iep,'L'),  ...
                   sprintf('%s L real', ep_names{iep}), '-');
        plot_trace(rz_t_ms, rls, is_learner, col_for_epoch(iep,'L'),  '', ':');
    end
    xline(0, 'k--', 'HandleVisibility','off');
    title(sprintf('%s (Learners)', gr.pair_name), 'Interpreter','none');
    xlabel('Time before RZ entry (ms)'); ylabel('CC1');
    ylim([0 1]);
end
sgtitle('RZ-pre per-bin raw CC1 (solid: real, dotted: shuffle) - Learners');
save_fig(gcf, fullfile(data_dir, 'TemporalCCAv4_RZ_CC_raw_Learners.svg'));

if any(~is_learner)
    figure('Name','RZ-pre CC raw NL','Color','w','Position',[100 100 1600 900]);
    tiledlayout('flow','TileSpacing','compact','Padding','compact');
    for ipair = 1:n_pairs
        nexttile; hold on;
        gr = group_results(ipair);
        for iep = 1:3
            rl  = gr.(['rz_cc1_'    ep_keys{iep}]);
            rls = gr.(['rz_cc1_sh_' ep_keys{iep}]);
            plot_trace(rz_t_ms, rl,  ~is_learner, col_for_epoch(iep,'NL'), ...
                       sprintf('%s NL real', ep_names{iep}), '-');
            plot_trace(rz_t_ms, rls, ~is_learner, col_for_epoch(iep,'NL'), '', ':');
        end
        xline(0, 'k--', 'HandleVisibility','off');
        title(sprintf('%s (Non-Learners)', gr.pair_name), 'Interpreter','none');
        xlabel('Time before RZ entry (ms)'); ylabel('CC1');
        ylim([0 1]);
    end
    sgtitle('RZ-pre per-bin raw CC1 (solid: real, dotted: shuffle) - Non-Learners');
    save_fig(gcf, fullfile(data_dir, 'TemporalCCAv4_RZ_CC_raw_NonLearners.svg'));
end

%% D) RZ-PRE IFI SCALAR BARS (rmANOVA)
fprintf('\n=== Section D: RZ-pre IFI bars (rmANOVA) ===\n');
figure('Name','RZ-pre IFI (rmANOVA)','Color','w','Position',[100 100 1500 850]);
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for ipair = 1:n_pairs
    nexttile; hold on;
    gr = group_results(ipair);
    fprintf('\n--- %s : RZ-pre IFI ---\n', gr.pair_name);
    plot_bars_with_rmanova( ...
        gr.rz_ifi_early,  gr.rz_ifi_pre,  gr.rz_ifi_post, ...
        gr.rz_ifi_sh_early, gr.rz_ifi_sh_pre, gr.rz_ifi_sh_post, ...
        is_learner, gr.pair_name, 'RZ-pre IFI', col_L, col_NL, 'zero');
end
sgtitle('RZ-pre per-trial IFI (epoch mean), * = different from 0');
save_fig(gcf, fullfile(data_dir, 'TemporalCCAv4_RZIFI_Bars.svg'));

%% E) EXCESS (real - shuffle) PER-TRIAL CC1 BARS (rmANOVA)
fprintf('\n=== Section E: per-trial CC1 excess (real - shuffle) bars ===\n');
figure('Name','Per-trial CC1 excess','Color','w','Position',[100 100 1500 850]);
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for ipair = 1:n_pairs
    nexttile; hold on;
    gr = group_results(ipair);
    fprintf('\n--- %s : per-trial CC1 excess ---\n', gr.pair_name);
    e = gr.trial_cc1_early - gr.trial_cc1_sh_early;
    p = gr.trial_cc1_pre   - gr.trial_cc1_sh_pre;
    x = gr.trial_cc1_post  - gr.trial_cc1_sh_post;
    plot_bars_with_rmanova(e, p, x, [], [], [], is_learner, gr.pair_name, ...
        'CC1 (real-shuf)', col_L, col_NL, 'zero');
end
sgtitle('Per-trial CC1 EXCESS over shuffle, * = excess different from 0');
save_fig(gcf, fullfile(data_dir, 'TemporalCCAv4_TrialCC1_Excess_Bars.svg'));

%% F) DIAGNOSTIC PANEL
fprintf('\n=== Section F: diagnostic ===\n');
figure('Name','v4 diagnostic','Color','w','Position',[100 100 1700 950]);
tiledlayout(2, 3, 'TileSpacing','compact','Padding','compact');

% F1 - finite fraction: per-trial CC1
nexttile;
mat = nan(n_pairs, 3);
for ipair = 1:n_pairs
    for iep = 1:3
        v = group_results(ipair).(['trial_cc1_' ep_keys{iep}]);
        mat(ipair, iep) = mean(isfinite(v));
    end
end
imagesc(mat, [0 1]); colorbar;
yticks(1:n_pairs); yticklabels({group_results.pair_name});
xticks(1:3); xticklabels(ep_names);
title('Frac. animals w/ finite per-trial CC1');

% F2 - finite fraction: per-trial IFI
nexttile;
mat = nan(n_pairs, 3);
for ipair = 1:n_pairs
    for iep = 1:3
        v = group_results(ipair).(['trial_ifi_' ep_keys{iep}]);
        mat(ipair, iep) = mean(isfinite(v));
    end
end
imagesc(mat, [0 1]); colorbar;
yticks(1:n_pairs); yticklabels({group_results.pair_name});
xticks(1:3); xticklabels(ep_names);
title('Frac. animals w/ finite per-trial IFI');

% F3 - finite fraction of RZ CC bins, post epoch
nexttile;
mat = nan(n_pairs, n_animals);
for ipair = 1:n_pairs
    rl = group_results(ipair).rz_cc1_post;     % [animals x n_rz_bins]
    if ~isempty(rl)
        mat(ipair, :) = mean(isfinite(rl), 2)';
    end
end
imagesc(mat, [0 1]); colorbar;
yticks(1:n_pairs); yticklabels({group_results.pair_name});
xlabel('Animal #'); xticks(1:n_animals);
title('Frac. RZ-bins finite (post epoch)');

% F4 - real vs shuffle CC scatter (post)
nexttile; hold on;
r_all = []; s_all = []; pair_id = [];
for ipair = 1:n_pairs
    r = group_results(ipair).trial_cc1_post(:);
    s = group_results(ipair).trial_cc1_sh_post(:);
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
    title('Real vs shuffle CC1 (post)');
    cb = colorbar; cb.Ticks = 1:n_pairs; cb.TickLabels = {group_results.pair_name};
end

% F5 - real vs shuffle IFI scatter (post)
nexttile; hold on;
r_all = []; s_all = []; pair_id = [];
for ipair = 1:n_pairs
    r = group_results(ipair).trial_ifi_post(:);
    s = group_results(ipair).trial_ifi_sh_post(:);
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
    title('Real vs shuffle IFI (post)');
    cb = colorbar; cb.Ticks = 1:n_pairs; cb.TickLabels = {group_results.pair_name};
end

% F6 - animal x pair contribution matrix
nexttile;
contrib = false(n_animals, n_pairs);
for ipair = 1:n_pairs
    contrib(:, ipair) = isfinite(group_results(ipair).trial_cc1_post(:));
end
imagesc(contrib'); colormap(gca, [0.95 0.95 0.95; 0.2 0.6 0.3]);
yticks(1:n_pairs); yticklabels({group_results.pair_name});
xlabel('Animal #'); xticks(1:n_animals);
title('Contribution matrix (post epoch, green = finite)');

sgtitle(sprintf(['v4 diagnostic | bin=%d ms, n_shuf=%d, k cap=%d, ' ...
                 'min_block=%d, lags=[%d..%d], RZ pre %d bins'], ...
    target_bin_ms, n_shuffles, max_k, min_block_bins, lags(1), lags(end), n_rz_bins));
save_fig(gcf, fullfile(data_dir, 'TemporalCCAv4_Diagnostic.svg'));

fprintf('\nAll v4 figures saved to %s\n', data_dir);


%% LOCAL HELPERS

function plot_bars_with_rmanova(re, rp, rx, se, sp, sx, is_learner, ...
                                 pair_name, y_label, col_L, col_NL, sig_mode)
% Grouped bars (Learner vs Non-Learner) over Early/Pre/Post with:
%   - per-bar 1-sample t-test (* if p<0.05) against either the shuffle null
%     (sig_mode='shuffle': paired t-test real vs shuf) or zero
%     (sig_mode='zero')
%   - mixed rmANOVA (Epoch within, Group between) with results overlaid as a text box.
% se/sp/sx may be empty; if empty AND sig_mode=='shuffle' the per-bar test is skipped.
    re = re(:); rp = rp(:); rx = rx(:);
    L  = logical(is_learner(:));
    NL = ~L;
    grp = {L, NL};
    cols = {col_L, col_NL};

    means_real = nan(2, 3); sems_real = nan(2, 3);
    for g = 1:2
        v_e = re(grp{g}); v_p = rp(grp{g}); v_x = rx(grp{g});
        means_real(g, :) = [mean(v_e,'omitnan'), mean(v_p,'omitnan'), mean(v_x,'omitnan')];
        sems_real (g, :) = [sem_(v_e), sem_(v_p), sem_(v_x)];
    end

    hold on;
    b = bar(1:3, means_real', 'grouped');
    b(1).FaceColor = col_L;
    b(2).FaceColor = col_NL;
    x_off = [b(1).XEndPoints; b(2).XEndPoints]';

    errorbar(x_off(:,1), means_real(1,:), sems_real(1,:), 'k.', 'LineWidth', 1.2);
    errorbar(x_off(:,2), means_real(2,:), sems_real(2,:), 'k.', 'LineWidth', 1.2);
    yline(0, '-k', 'LineWidth', 0.5);
    xticks(1:3); xticklabels({'Early','Pre','Post'});
    title(pair_name, 'Interpreter','none'); ylabel(y_label);

    yl = ylim; offset = (yl(2) - yl(1)) * 0.05;
    min_y = yl(1); max_y = yl(2);

    use_shuf = strcmpi(sig_mode, 'shuffle') && ~isempty(se) && ~isempty(sp) && ~isempty(sx);
    if use_shuf
        se = se(:); sp = sp(:); sx = sx(:);
    end

    p_g = nan(2, 3);
    for g = 1:2
        v = {re(grp{g}), rp(grp{g}), rx(grp{g})};
        if use_shuf
            vs = {se(grp{g}), sp(grp{g}), sx(grp{g})};
        end
        for ie = 1:3
            if use_shuf
                pair_fin = isfinite(v{ie}) & isfinite(vs{ie});
                if sum(pair_fin) > 2
                    [~, p_g(g, ie)] = ttest(v{ie}(pair_fin), vs{ie}(pair_fin));
                end
            else
                fin = isfinite(v{ie});
                if sum(fin) > 2
                    [~, p_g(g, ie)] = ttest(v{ie}(fin), 0);
                end
            end
            if ~isnan(p_g(g, ie)) && p_g(g, ie) < 0.05
                if means_real(g, ie) >= 0
                    yp = means_real(g, ie) + sems_real(g, ie) + offset;
                    va = 'bottom';
                else
                    yp = means_real(g, ie) - sems_real(g, ie) - offset;
                    va = 'top';
                end
                if yp > max_y, max_y = yp + offset; end
                if yp < min_y, min_y = yp - offset; end
                text(x_off(ie, g), yp, '*', 'FontSize', 18, ...
                     'HorizontalAlignment','center','VerticalAlignment',va, 'Color', cols{g});
            end
        end
    end

    ylim([min(yl(1), min_y), max(yl(2), max_y)]);

    % --- Mixed rmANOVA (Epoch within, Group between) ---
    valid_mask = isfinite(re) & isfinite(rp) & isfinite(rx);
    n_complete = sum(valid_mask);
    if n_complete > max(2, length(unique(is_learner(valid_mask))))
        T = table(re(valid_mask), rp(valid_mask), rx(valid_mask), ...
                  categorical(is_learner(valid_mask), [1 0], {'Learner','NonLearner'}), ...
                  'VariableNames', {'Early','Pre','Post','LearnerGroup'});
        Meas = table([1 2 3]', 'VariableNames', {'Epoch'});
        try
            rm = fitrm(T, 'Early-Post ~ LearnerGroup', 'WithinDesign', Meas);
            ranovatbl  = ranova(rm);
            betweentbl = anova(rm);
            p_epoch = ranovatbl.pValue(1);
            p_int   = ranovatbl.pValue(2);
            p_group = betweentbl.pValue(2);

            disp('--- Within-Subjects (Epoch & Interaction) ---'); disp(ranovatbl);
            disp('--- Between-Subjects (LearnerGroup) ---');       disp(betweentbl);

            sig_text = {sprintf('rmANOVA (n=%d):', n_complete)};
            has_sig = false;
            if p_epoch < 0.05, sig_text{end+1} = sprintf('Epoch: p=%.3f', p_epoch);   has_sig = true; end
            if p_group < 0.05, sig_text{end+1} = sprintf('Group: p=%.3f', p_group);   has_sig = true; end
            if p_int   < 0.05, sig_text{end+1} = sprintf('Ep x Grp: p=%.3f', p_int);  has_sig = true; end
            if ~has_sig, sig_text{end+1} = 'n.s.'; end

            curr_yl = ylim; curr_xl = xlim;
            yr = curr_yl(2) - curr_yl(1);
            ylim([curr_yl(1), curr_yl(2) + yr*0.25]);
            new_yl = ylim;
            text(curr_xl(1) + 0.03*(curr_xl(2)-curr_xl(1)), ...
                 new_yl(2) - 0.03*(new_yl(2)-new_yl(1)), ...
                 strjoin(sig_text, '\n'), 'VerticalAlignment','top', ...
                 'FontSize', 9, 'EdgeColor','k', 'BackgroundColor','w', 'Margin', 3);

            % Post-hocs (printed only)
            if p_int < 0.05
                fprintf('  -> interaction sig. Post-hoc Group | Epoch:\n');
                disp(multcompare(rm, 'LearnerGroup', 'By', 'Epoch'));
            else
                if p_epoch < 0.05
                    fprintf('  -> Epoch main effect sig. Post-hoc Epoch:\n');
                    disp(multcompare(rm, 'Epoch'));
                end
                if p_group < 0.05
                    fprintf('  -> Group main effect sig. Post-hoc Group:\n');
                    disp(multcompare(rm, 'LearnerGroup'));
                end
            end
        catch ME
            fprintf('rmANOVA failed: %s\n', ME.message);
        end
    else
        fprintf('Not enough complete cases (n=%d) for rmANOVA.\n', n_complete);
    end
end

function plot_trace(x, mat, mask, col, name, style)
% Mean +/- SEM across animals (mat: [animals x n_x]).
    mask = logical(mask(:));
    if size(mat, 1) ~= numel(mask)
        if size(mat, 2) == numel(mask), mat = mat'; end
    end
    sub = mat(mask, :);
    if isempty(sub) || all(~isfinite(sub(:))), return; end
    n  = sum(isfinite(sub), 1);
    mu = mean(sub, 1, 'omitnan');
    se = std(sub, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
    keep = isfinite(mu);
    if ~any(keep), return; end
    if isempty(name)
        plot(x(keep), mu(keep), style, 'Color', col, 'LineWidth', 1.0, 'HandleVisibility','off');
    else
        fill([x(keep), fliplr(x(keep))], [mu(keep)+se(keep), fliplr(mu(keep)-se(keep))], ...
             col, 'EdgeColor','none','FaceAlpha',0.15, 'HandleVisibility','off');
        plot(x(keep), mu(keep), style, 'Color', col, 'LineWidth', 1.4, 'DisplayName', name);
    end
end

function s = sem_(v)
    v = v(isfinite(v));
    if isempty(v), s = NaN; else, s = std(v) / sqrt(numel(v)); end
end

function c = col_for_epoch(iep, group)
    if strcmpi(group, 'L')
        cmap = [0.65 0.78 1.00; 0.30 0.55 0.95; 0.05 0.25 0.60];
    else
        cmap = [0.95 0.78 0.65; 0.85 0.50 0.25; 0.55 0.20 0.05];
    end
    c = cmap(iep, :);
end

function save_fig(h, path)
    try
        set(h, 'Renderer', 'painters');
        saveas(h, path);
        fprintf('  saved %s\n', path);
    catch ME
        warning('Failed to save %s: %s', path, ME.message);
    end
end
