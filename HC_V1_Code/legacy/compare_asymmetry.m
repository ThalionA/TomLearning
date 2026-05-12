%% compare_asymmetry.m
% Per-pair summary of how the two methods diverge (Figure 6).
%
% Shows:
%   (a) Bar chart of mean(temporal_cc1) - mean(spatial_cc1) per pair x
%       epoch. Both real and shuffle-corrected (excess) variants. A bar
%       above zero means temporal CC1 is higher; below zero means
%       spatial CC1 is higher.
%   (b) IFI sign-agreement rate per pair x epoch. For each (animal,
%       trial), compute sign(temporal_ifi) and sign(spatial_ifi); report
%       fraction of trials where the signs match (excluding any zeros /
%       NaNs).
%   (c) IFI magnitude correlation per pair x epoch (Spearman of
%       |temporal_ifi| vs |spatial_ifi|). Probes whether the two methods
%       agree on which trials have stronger asymmetry (regardless of
%       direction).

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
n_epochs = numel(epochs);

% Excess
t_cc_x = ct.temporal_cc1 - ct.temporal_cc1_sh;
s_cc_x = ct.spatial_cc1  - ct.spatial_cc1_sh;

%% Aggregate per pair x epoch
asym_real    = nan(n_pairs, n_epochs);
asym_excess  = nan(n_pairs, n_epochs);
sign_agree   = nan(n_pairs, n_epochs);
mag_rho      = nan(n_pairs, n_epochs);
n_per_cell   = zeros(n_pairs, n_epochs);

for ipair = 1:n_pairs
    for ie = 1:n_epochs
        sel = strcmp(ct.pair, pairs{ipair}) & strcmp(ct.epoch, epochs{ie});
        if ~any(sel), continue; end

        tcc  = ct.temporal_cc1(sel);
        scc  = ct.spatial_cc1(sel);
        tcx  = t_cc_x(sel);
        scx  = s_cc_x(sel);
        tifi = ct.temporal_ifi(sel);
        sifi = ct.spatial_ifi(sel);

        ok_cc = isfinite(tcc) & isfinite(scc);
        if any(ok_cc)
            asym_real(ipair, ie)   = mean(tcc(ok_cc) - scc(ok_cc));
        end
        ok_x = isfinite(tcx) & isfinite(scx);
        if any(ok_x)
            asym_excess(ipair, ie) = mean(tcx(ok_x) - scx(ok_x));
        end

        ok_ifi = isfinite(tifi) & isfinite(sifi) & abs(tifi) > 0 & abs(sifi) > 0;
        if any(ok_ifi)
            agree = sign(tifi(ok_ifi)) == sign(sifi(ok_ifi));
            sign_agree(ipair, ie) = mean(agree);
            mag_rho(ipair, ie) = corr(abs(tifi(ok_ifi)), abs(sifi(ok_ifi)), 'type', 'Spearman');
        end
        n_per_cell(ipair, ie) = sum(ok_cc);
    end
end

%% Figure 6 — four-panel summary
fig6 = figure('Name', 'Asymmetry summary', 'Color', 'w', ...
              'Position', [100 100 1500 900]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Method-asymmetry summary across pairs');

% Panel A: real CC1 asymmetry
ax = nexttile; bar(ax, asym_real); hold(ax, 'on');
yline(ax, 0, 'k-');
xticks(ax, 1:n_pairs); xticklabels(ax, pairs); xtickangle(ax, 30);
ylabel(ax, 'mean(temporal CC1) - mean(spatial CC1)');
title(ax, 'A.  CC1 asymmetry (real values)');
legend(ax, epochs, 'Location', 'best');
grid(ax, 'on');

% Panel B: excess CC1 asymmetry
ax = nexttile; bar(ax, asym_excess); hold(ax, 'on');
yline(ax, 0, 'k-');
xticks(ax, 1:n_pairs); xticklabels(ax, pairs); xtickangle(ax, 30);
ylabel(ax, 'mean(temporal \DeltaCC1) - mean(spatial \DeltaCC1)');
title(ax, 'B.  CC1 asymmetry (shuffle-corrected)');
legend(ax, epochs, 'Location', 'best');
grid(ax, 'on');

% Panel C: IFI sign agreement
ax = nexttile; bar(ax, sign_agree); hold(ax, 'on');
yline(ax, 0.5, 'k--');
xticks(ax, 1:n_pairs); xticklabels(ax, pairs); xtickangle(ax, 30);
ylabel(ax, 'P(sign agreement)  [chance = 0.5]');
title(ax, 'C.  IFI sign agreement between methods');
ylim(ax, [0 1]);
legend(ax, epochs, 'Location', 'best');
grid(ax, 'on');

% Panel D: IFI magnitude correlation
ax = nexttile; bar(ax, mag_rho); hold(ax, 'on');
yline(ax, 0, 'k-');
xticks(ax, 1:n_pairs); xticklabels(ax, pairs); xtickangle(ax, 30);
ylabel(ax, 'Spearman \rho ( |IFI_T|, |IFI_S| )');
title(ax, 'D.  IFI magnitude correlation between methods');
legend(ax, epochs, 'Location', 'best');
grid(ax, 'on');

saveas(fig6, fullfile(fig_dir, 'TemporalVsSpatial_asymmetry_summary.svg'));

%% Console table
fprintf('\n=== Asymmetry summary table ===\n');
fprintf('%-15s %-6s %12s %14s %14s %12s %5s\n', ...
    'pair', 'epoch', 'asym_real', 'asym_excess', 'sign_agree', 'mag_rho', 'n');
for ipair = 1:n_pairs
    for ie = 1:n_epochs
        fprintf('%-15s %-6s %12.3f %14.3f %14.2f %12.2f %5d\n', ...
            pairs{ipair}, epochs{ie}, ...
            asym_real(ipair, ie), asym_excess(ipair, ie), ...
            sign_agree(ipair, ie), mag_rho(ipair, ie), n_per_cell(ipair, ie));
    end
end

% Save the numbers behind the figure
save(fullfile(data_dir, 'compare_asymmetry.mat'), ...
    'pairs', 'epochs', 'asym_real', 'asym_excess', 'sign_agree', ...
    'mag_rho', 'n_per_cell');
fprintf('\nFigure 6 saved to %s\n', fig_dir);
