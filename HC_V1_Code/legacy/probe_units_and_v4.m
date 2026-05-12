%% probe_units_and_v4.m
% Two probes for the temporal-vs-spatial RSC mystery.
%   Probe 1: inspect units.region / units.regions_label / units.idx in each
%            TF*_export file. Are the per-unit and idx-row region tags
%            consistent? RSC is the suspect.
%   Probe 2: in the v4 results .mat, tally which (animal x pair) cells
%            are all-NaN, broken out by whether the pair contains RSC.
%
% Writes a plain-text report to probe_units_and_v4_output.txt.
% Run from the repo root.

clear; clc;

base_dir   = '/Users/theoamvr/Desktop/Experiments/TomLearning';
data_dir   = fullfile(base_dir, 'HC_V1_data');
out_path   = fullfile(base_dir, 'probe_units_and_v4_output.txt');
v4_results = fullfile(data_dir, 'Temporal_CCA_v4_2026_05_05.mat');

fid = fopen(out_path, 'w');
both(fid, '=== probe_units_and_v4.m  %s ===\n\n', datestr(now));

%% Probe 1: per-animal units struct
both(fid, '--- Probe 1: units.region vs units.idx vs units.regions_label ---\n\n');

files = dir(fullfile(data_dir, 'TF*_export.mat'));
target_areas = {'CA1','V1','DG','CA3','RSC','SUB'};

for ii = 1:numel(files)
    fp = fullfile(files(ii).folder, files(ii).name);
    both(fid, 'Animal: %s\n', files(ii).name);
    try
        D = load(fp, 'units');
        u = D.units;
    catch ME
        both(fid, '  [load failed: %s]\n\n', ME.message);
        continue;
    end

    has_region = isfield(u, 'region');
    has_idx    = isfield(u, 'idx');
    has_rl     = isfield(u, 'regions_label');
    has_uid    = isfield(u, 'unit_id');

    nU = NaN;
    if has_uid, nU = numel(u.unit_id); end
    both(fid, '  n_units = %d   has region=%d  idx=%d  regions_label=%d\n', ...
         nU, has_region, has_idx, has_rl);

    rl = {};
    if has_rl
        rl = cellstr(u.regions_label);
        both(fid, '  regions_label (%d): { %s }\n', numel(rl), strjoin(rl, ', '));
    end

    reg = {};
    if has_region
        reg = cellstr(u.region(:));
        both(fid, '  size(units.region) = [%d %d]\n', size(u.region,1), size(u.region,2));
    end

    if has_idx
        both(fid, '  size(units.idx)    = [%d %d]\n', size(u.idx,1), size(u.idx,2));
    end

    % Per-area count comparison
    both(fid, '  per-area counts:\n');
    both(fid, '    %6s | region | idx  | match?\n', '');
    for k = 1:numel(target_areas)
        a = target_areas{k};
        c_reg = NaN; c_idx = NaN;
        if has_region
            c_reg = sum(strcmp(reg, a));
        end
        if has_idx && has_rl
            row = strcmp(rl, a);
            if any(row)
                c_idx = sum(logical(u.idx(row, :)));
            else
                c_idx = -1; % area not in regions_label at all
            end
        end
        if isnan(c_reg) || isnan(c_idx)
            mark = '?';
        elseif c_idx == -1
            mark = 'no-row';
        elseif c_reg == c_idx
            mark = 'OK';
        else
            mark = 'DIFFER';
        end
        both(fid, '    %6s | %5s  | %5s | %s\n', a, num2str(c_reg), num2str(c_idx), mark);
    end

    % Multi-membership: how many units appear in >1 idx row?
    if has_idx
        col_sum = sum(logical(u.idx), 1);
        both(fid, '  units in >1 idx row : %d / %d\n', sum(col_sum > 1), numel(col_sum));
    end

    both(fid, '\n');
end

%% Probe 2: v4 results, NaN audit by pair
both(fid, '--- Probe 2: Temporal_CCA_v4_2026_05_05.mat NaN audit ---\n\n');

if ~exist(v4_results, 'file')
    both(fid, '  results file not found: %s\n', v4_results);
else
    R = load(v4_results);
    if ~isfield(R, 'group_results')
        both(fid, '  group_results field missing\n');
    else
        gr = R.group_results;
        n_pairs = numel(gr);
        epoch_names = {'early','pre','post'};
        % find n_animals from any trial_cc1 field
        n_animals = NaN;
        for ee = 1:numel(epoch_names)
            f = ['trial_cc1_' epoch_names{ee}];
            if isfield(gr(1), f)
                n_animals = numel(gr(1).(f));
                break;
            end
        end
        both(fid, '  n_pairs   = %d\n', n_pairs);
        both(fid, '  n_animals = %d\n', n_animals);
        both(fid, '  epochs    = %s\n\n', strjoin(epoch_names, ', '));

        both(fid, '  pair                    | trial_cc1 NaN per epoch (early/pre/post) | rz_cc1 NaN per epoch | contains_RSC\n');
        for p = 1:n_pairs
            name = gr(p).pair_name;
            has_rsc = any(strcmpi(strsplit(name,'-'), 'RSC'));
            tc = '';
            rc = '';
            for ee = 1:numel(epoch_names)
                fcc  = ['trial_cc1_'  epoch_names{ee}];
                frcc = ['rz_cc1_'     epoch_names{ee}];
                v_tc  = NaN; v_rcc = NaN;
                if isfield(gr(p), fcc),  v_tc  = sum(isnan(gr(p).(fcc)(:))); end
                if isfield(gr(p), frcc)
                    rcc_arr = gr(p).(frcc);
                    if size(rcc_arr,1) == n_animals
                        v_rcc = sum(all(isnan(rcc_arr), 2));
                    else
                        v_rcc = sum(isnan(rcc_arr(:)));
                    end
                end
                tc = sprintf('%s%2d/%d ', tc, v_tc, n_animals);
                rc = sprintf('%s%2d/%d ', rc, v_rcc, n_animals);
            end
            both(fid, '  %-23s | %s| %s| %s\n', name, tc, rc, mat2str(has_rsc));
        end
    end
end

both(fid, '\n=== done ===\n');
fclose(fid);
fprintf('\nWrote %s\n', out_path);


%% Local helper: write a printf to BOTH stdout and a file handle.
function both(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end
