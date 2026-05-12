%% compare_temporal_spatial.m
% Trial-by-trial paired comparison of Temporal CCA v4 vs Spatial CCA v2.
%
% For each (animal, pair, epoch, trial_id) the script extracts:
%   - temporal CC1   : group_results(iT).trial_cc1_pertrial_<ep>{ia}
%   - spatial  CC1   : group_results(iS).trial_corr_<ep>{ia}
%   - temporal IFI   : group_results(iT).trial_ifi_pertrial_<ep>{ia}
%   - spatial  IFI   : group_results(iS).trial_precession_<ep>_idx{ia}
%
% Trials are paired by trial_id (temporal stores it explicitly; spatial
% reconstructs from analysis_lp + epoch position). Pairs that appear in
% only one .mat are reported but skipped.
%
% Outputs:
%   1. Combined paired scatter (3 epochs x {CC1, IFI} = 6 panels). All
%      pairs co-plotted, color-coded by pair. Identity line. Spearman rho
%      reported per panel.
%   2. Per-pair paired scatter (n_pairs rows x {CC1 early/pre/post, IFI
%      early/pre/post} = 6 cols). Color-coded by learner / non-learner.
%      Per-panel Spearman rho.
%   3. compare_summary.mat with the paired arrays so the figures can be
%      regenerated without reloading the source results.
%
% Conventions:
%   - Temporal IFI is signed lag asymmetry along TIME (-3:3 bins of 25 ms).
%   - Spatial  IFI is signed lag asymmetry along POSITION (-3:3 bins of 2.5 cm).
%   They should agree in SIGN when within-trial speed is roughly constant
%   but their magnitudes are not directly comparable.

clear; clc; close all;

base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning';
data_dir = fullfile(base_dir, 'HC_V1_data');
fig_dir  = fullfile(base_dir, 'HC_V1_figures', 'CCA_Compare');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

addpath(fullfile(base_dir, 'HC_V1_Code'));

%% Discover latest result files
T_files = dir(fullfile(data_dir, 'Temporal_CCA_v4_*.mat'));
S_files = dir(fullfile(data_dir, 'Spatial_CCA_Results_*.mat'));
assert(~isempty(T_files), 'No Temporal_CCA_v4_*.mat found in %s', data_dir);
assert(~isempty(S_files), 'No Spatial_CCA_Results_*.mat found in %s', data_dir);
[~, iT] = max([T_files.datenum]); T_path = fullfile(T_files(iT).folder, T_files(iT).name);
[~, iS] = max([S_files.datenum]); S_path = fullfile(S_files(iS).folder, S_files(iS).name);
fprintf('Loading %s\n        %s\n', T_path, S_path);
T = load(T_path);
S = load(S_path);

%% Reconcile pair sets
pair_names_T = arrayfun(@(g) g.pair_name, T.group_results, 'uni', 0);
pair_names_S = arrayfun(@(g) g.pair_name, S.group_results, 'uni', 0);
common = intersect(pair_names_T, pair_names_S, 'stable');
only_T = setdiff(pair_names_T, common);
only_S = setdiff(pair_names_S, common);
if ~isempty(only_T)
    fprintf('Pairs only in temporal:  %s\n', strjoin(only_T, ', '));
end
if ~isempty(only_S)
    fprintf('Pairs only in spatial:   %s\n', strjoin(only_S, ', '));
end
fprintf('Pairs in both (will compare): %s\n\n', strjoin(common, ', '));
n_pairs = numel(common);

%% Sanity: confirm the temporal results have per-trial fields
required_T_fields = {};
for ie = 1:3
    ep = {'early','pre','post'}; ep = ep{ie};
    required_T_fields = [required_T_fields, {['trial_cc1_pertrial_' ep], ...
                                              ['trial_ifi_pertrial_' ep], ...
                                              ['trial_id_pertrial_'  ep]}]; %#ok<AGROW>
end
have_T = isfield(T.group_results, required_T_fields);
if ~all(have_T)
    missing = required_T_fields(~have_T);
    error(['Temporal results .mat is missing per-trial fields: %s\n' ...
           'Re-run CCA_HC_V1_temporal_v4.m with the current version of ' ...
           'the script (it must include the trial_*_pertrial_* and ' ...
           'trial_id_pertrial_* assignments).'], strjoin(missing, ', '));
end

required_S_fields = {};
for ie = 1:3
    ep = {'early','pre','post'}; ep = ep{ie};
    required_S_fields = [required_S_fields, ...
        {['trial_corr_' ep], ['trial_precession_' ep '_idx']}]; %#ok<AGROW>
end
have_S = isfield(S.group_results, required_S_fields);
if ~all(have_S)
    missing = required_S_fields(~have_S);
    error('Spatial results .mat is missing fields: %s', strjoin(missing, ', '));
end

% Skip-reason counters so a 0-row outcome is debuggable.
skip_counts = struct('lp_nan', 0, 'sp_empty', 0, 'tp_empty', 0, ...
                     'sp_len_mismatch', 0, 'no_trial_intersection', 0);

% Animal count and learner mask
assert(numel(T.analysis_lp) == numel(S.analysis_lp), ...
    'analysis_lp length differs between temporal and spatial.');
n_animals = numel(T.analysis_lp);
assert(all(T.analysis_lp == S.analysis_lp | (isnan(T.analysis_lp) & isnan(S.analysis_lp))), ...
    'analysis_lp values differ between temporal and spatial — check animal sorting.');
is_learner = T.is_learner;
analysis_lp = T.analysis_lp;

epochs = {'early','pre','post'};

%% Build paired arrays (long form)
% Each row is a single (animal, pair, epoch, trial_id) datum.
rows_pair      = {};
rows_animal    = [];
rows_is_learner= [];
rows_epoch     = {};
rows_trial_id  = [];
rows_t_cc      = [];
rows_t_ifi     = [];
rows_s_cc      = [];
rows_s_ifi     = [];
rows_t_cc_sh   = [];
rows_t_ifi_sh  = [];
rows_s_cc_sh   = [];
rows_s_ifi_sh  = [];

for ipair = 1:n_pairs
    name = common{ipair};
    idxT = find(strcmp(pair_names_T, name), 1);
    idxS = find(strcmp(pair_names_S, name), 1);

    for ie = 1:numel(epochs)
        ep = epochs{ie};
        for ia = 1:n_animals
            lp = analysis_lp(ia);
            % Skip animals with unusable LP (matches v4 logic).
            if isnan(lp)
                skip_counts.lp_nan = skip_counts.lp_nan + 1; continue;
            end
            switch ep
                case 'early', sp_ids = 1:10;
                case 'pre',   sp_ids = (lp-10):(lp-1);
                case 'post',  sp_ids = lp:(lp+9);
            end

            % Spatial side (positional). Real + shuffle.
            sp_cc       = S.group_results(idxS).(['trial_corr_' ep]){ia};
            sp_ifi      = S.group_results(idxS).(['trial_precession_' ep '_idx']){ia};
            sp_cc_sh    = S.group_results(idxS).(['trial_corr_' ep '_shuff']){ia};
            sp_ifi_sh   = S.group_results(idxS).(['trial_precession_' ep '_idx_shuff']){ia};
            if isempty(sp_cc) || isempty(sp_ifi)
                skip_counts.sp_empty = skip_counts.sp_empty + 1; continue;
            end
            % spatial stores [num_ccs x n_trials]; we want CC1 row.
            sp_cc       = sp_cc(1, :);
            sp_ifi      = sp_ifi(1, :);
            if ~isempty(sp_cc_sh),  sp_cc_sh  = sp_cc_sh(1, :);  else sp_cc_sh  = nan(size(sp_cc));  end
            if ~isempty(sp_ifi_sh), sp_ifi_sh = sp_ifi_sh(1, :); else sp_ifi_sh = nan(size(sp_ifi)); end
            if numel(sp_cc) ~= numel(sp_ids)
                skip_counts.sp_len_mismatch = skip_counts.sp_len_mismatch + 1;
                fprintf('  [warn] spatial len=%d, expected=%d (animal %d, ep %s, pair %s)\n', ...
                    numel(sp_cc), numel(sp_ids), ia, ep, name);
                continue;
            end

            % Temporal side (explicit trial_id). Real + shuffle.
            tp_cc      = T.group_results(idxT).(['trial_cc1_pertrial_'    ep]){ia};
            tp_ifi     = T.group_results(idxT).(['trial_ifi_pertrial_'    ep]){ia};
            tp_cc_sh   = T.group_results(idxT).(['trial_cc1_sh_pertrial_' ep]){ia};
            tp_ifi_sh  = T.group_results(idxT).(['trial_ifi_sh_pertrial_' ep]){ia};
            tp_ids     = T.group_results(idxT).(['trial_id_pertrial_'     ep]){ia};
            if isempty(tp_cc)
                skip_counts.tp_empty = skip_counts.tp_empty + 1; continue;
            end
            tp_cc     = tp_cc(:);     tp_ifi     = tp_ifi(:);
            tp_ids    = tp_ids(:);
            if ~isempty(tp_cc_sh),  tp_cc_sh  = tp_cc_sh(:);  else tp_cc_sh  = nan(size(tp_cc));  end
            if ~isempty(tp_ifi_sh), tp_ifi_sh = tp_ifi_sh(:); else tp_ifi_sh = nan(size(tp_ifi)); end

            % Pair by trial_id (intersection only).
            [common_ids, ~, ~] = intersect(sp_ids(:), tp_ids(:), 'stable');
            if isempty(common_ids)
                skip_counts.no_trial_intersection = skip_counts.no_trial_intersection + 1;
                if skip_counts.no_trial_intersection <= 3
                    fprintf('  [warn] no trial-id intersection (animal %d, ep %s, pair %s)\n', ia, ep, name);
                    fprintf('         spatial sp_ids   = %s\n', mat2str(sp_ids(:)'));
                    fprintf('         temporal tp_ids  = %s\n', mat2str(tp_ids(:)'));
                end
                continue;
            end
            for cid = common_ids(:)'
                sp_pos = find(sp_ids == cid, 1);
                tp_pos = find(tp_ids == cid, 1);
                rows_pair{end+1, 1}      = name;
                rows_animal(end+1, 1)    = ia;
                rows_is_learner(end+1,1) = logical(is_learner(ia));
                rows_epoch{end+1, 1}     = ep;
                rows_trial_id(end+1,1)   = cid;
                rows_t_cc(end+1, 1)      = tp_cc(tp_pos);
                rows_t_ifi(end+1, 1)     = tp_ifi(tp_pos);
                rows_s_cc(end+1, 1)      = sp_cc(sp_pos);
                rows_s_ifi(end+1, 1)     = sp_ifi(sp_pos);
                rows_t_cc_sh(end+1, 1)   = tp_cc_sh(tp_pos);
                rows_t_ifi_sh(end+1, 1)  = tp_ifi_sh(tp_pos);
                rows_s_cc_sh(end+1, 1)   = sp_cc_sh(sp_pos);
                rows_s_ifi_sh(end+1, 1)  = sp_ifi_sh(sp_pos);
            end
        end
    end
end
fprintf('Total paired rows: %d\n', numel(rows_pair));
fprintf('Skip counts: lp_nan=%d, sp_empty=%d, tp_empty=%d, sp_len_mismatch=%d, no_trial_intersection=%d\n', ...
    skip_counts.lp_nan, skip_counts.sp_empty, skip_counts.tp_empty, ...
    skip_counts.sp_len_mismatch, skip_counts.no_trial_intersection);
if numel(rows_pair) == 0
    error(['compare_temporal_spatial: 0 paired rows.\n' ...
           'See skip counts above. If tp_empty dominates, the temporal ' ...
           'results were generated by an older v4 script that did not ' ...
           'save per-trial vectors — re-run CCA_HC_V1_temporal_v4.m.\n' ...
           'If no_trial_intersection dominates, the trial_id conventions ' ...
           'between the two pipelines have diverged.']);
end

%% Save the long-form table for downstream regeneration
compare_table = struct( ...
    'pair',       {rows_pair}, ...
    'animal',     rows_animal, ...
    'is_learner', rows_is_learner, ...
    'epoch',      {rows_epoch}, ...
    'trial_id',   rows_trial_id, ...
    'temporal_cc1',     rows_t_cc, ...
    'temporal_ifi',     rows_t_ifi, ...
    'spatial_cc1',      rows_s_cc, ...
    'spatial_ifi',      rows_s_ifi, ...
    'temporal_cc1_sh',  rows_t_cc_sh, ...
    'temporal_ifi_sh',  rows_t_ifi_sh, ...
    'spatial_cc1_sh',   rows_s_cc_sh, ...
    'spatial_ifi_sh',   rows_s_ifi_sh);
save(fullfile(data_dir, 'compare_summary.mat'), 'compare_table', '-v7.3');

%% Figure 1 — Combined paired scatter (3 epochs x {CC1, IFI})
% All pairs co-plotted, color = pair name.
cmap = lines(n_pairs);
fig1 = figure('Name', 'Temporal vs Spatial (combined)', 'Color', 'w', ...
              'Position', [100 100 1500 800]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Temporal CCA v4 vs Spatial CCA v2 — paired per (animal, trial)');

metrics = {'cc1', 'ifi'};
metric_labels = {'CC1', 'IFI'};
for im = 1:numel(metrics)
    for ie = 1:numel(epochs)
        ax = nexttile;
        hold(ax, 'on');
        ep = epochs{ie};

        sel_ep = strcmp(rows_epoch, ep);
        for ipair = 1:n_pairs
            sel = sel_ep & strcmp(rows_pair, common{ipair});
            if ~any(sel), continue; end
            if im == 1
                xs = rows_s_cc(sel); ys = rows_t_cc(sel);
            else
                xs = rows_s_ifi(sel); ys = rows_t_ifi(sel);
            end
            scatter(ax, xs, ys, 14, cmap(ipair, :), 'filled', ...
                'MarkerFaceAlpha', 0.5, 'DisplayName', common{ipair});
        end

        % Identity line
        xl = xlim(ax); yl = ylim(ax);
        lo = min([xl yl]); hi = max([xl yl]);
        plot(ax, [lo hi], [lo hi], 'k--', 'HandleVisibility', 'off');

        % Spearman rho (all paired points in this panel)
        if im == 1
            xs_all = rows_s_cc(sel_ep);  ys_all = rows_t_cc(sel_ep);
        else
            xs_all = rows_s_ifi(sel_ep); ys_all = rows_t_ifi(sel_ep);
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
saveas(fig1, fullfile(fig_dir, 'TemporalVsSpatial_combined.svg'));

%% Figure 2 — Per-pair paired scatter (n_pairs rows x 6 cols)
fig2 = figure('Name', 'Temporal vs Spatial (per pair)', 'Color', 'w', ...
              'Position', [100 100 1700 200 + 200 * n_pairs]);
tl2 = tiledlayout(n_pairs, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Temporal vs Spatial — per pair, per epoch');

col_L  = [0.1 0.3 0.7];
col_NL = [0.7 0.2 0.1];
for ipair = 1:n_pairs
    name = common{ipair};
    sel_pair = strcmp(rows_pair, name);
    for im = 1:numel(metrics)
        for ie = 1:numel(epochs)
            ax = nexttile;
            hold(ax, 'on');
            ep = epochs{ie};
            sel = sel_pair & strcmp(rows_epoch, ep);

            if im == 1
                xs = rows_s_cc(sel);  ys = rows_t_cc(sel);
            else
                xs = rows_s_ifi(sel); ys = rows_t_ifi(sel);
            end
            ll = logical(rows_is_learner(sel));   % cast back to logical
                                                   % (growing from [] coerces
                                                   % to double).

            scatter(ax, xs(ll),  ys(ll),  18, col_L,  'filled', ...
                'MarkerFaceAlpha', 0.6, 'DisplayName', 'Learner');
            scatter(ax, xs(~ll), ys(~ll), 18, col_NL, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'DisplayName', 'Non-learner');

            xl = xlim(ax); yl = ylim(ax);
            lo = min([xl yl]); hi = max([xl yl]);
            plot(ax, [lo hi], [lo hi], 'k--', 'HandleVisibility', 'off');

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
saveas(fig2, fullfile(fig_dir, 'TemporalVsSpatial_per_pair.svg'));

%% Console summary table
fprintf('\n=== Spearman rho summary ===\n');
fprintf('%-15s %-6s %8s %8s %8s\n', 'pair', 'epoch', 'rho_CC1', 'rho_IFI', 'n');
for ipair = 1:n_pairs
    for ie = 1:numel(epochs)
        ep = epochs{ie};
        sel = strcmp(rows_pair, common{ipair}) & strcmp(rows_epoch, ep);
        ok_cc  = sel & isfinite(rows_t_cc) & isfinite(rows_s_cc);
        ok_ifi = sel & isfinite(rows_t_ifi) & isfinite(rows_s_ifi);
        rho_cc  = NaN; rho_ifi = NaN;
        if sum(ok_cc)  >= 3, rho_cc  = corr(rows_s_cc(ok_cc),   rows_t_cc(ok_cc),   'type','Spearman'); end
        if sum(ok_ifi) >= 3, rho_ifi = corr(rows_s_ifi(ok_ifi), rows_t_ifi(ok_ifi), 'type','Spearman'); end
        fprintf('%-15s %-6s %8.2f %8.2f %8d\n', common{ipair}, ep, rho_cc, rho_ifi, sum(ok_cc));
    end
end
fprintf('\nFigures saved to %s\n', fig_dir);
fprintf('Long-form table:  %s\n', fullfile(data_dir, 'compare_summary.mat'));
