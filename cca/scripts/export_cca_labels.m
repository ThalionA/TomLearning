function export_cca_labels(data_dir, lp_column)
% EXPORT_CCA_LABELS  Decode MATLAB `string`-typed label fields into a plain
% JSON companion file that the tom_cca Python pipeline can read.
%
% Why this exists
% ---------------
% Tom's TF*_export.mat files store units.regions_label (and animal_behaviour.mat
% stores animal_id) as MATLAB `string` arrays. The Python pipeline reads .mat
% files with h5py, and h5py cannot decode the MATLAB `string` type -- so the
% pipeline cannot map the columns of units.idx to area names and cannot run.
% MATLAB handles `string` natively, so this script reads the labels here and
% writes HC_V1_data/cca_labels.json, which dataio.py then consumes.
%
% Usage (run once, from anywhere in MATLAB with this file on the path):
%   export_cca_labels                          % data_dir = ../../HC_V1_data
%   export_cca_labels('/path/to/HC_V1_data')   % explicit data directory
%   export_cca_labels('/path/to/HC_V1_data', 1)% LP = period_experienced(:,1)
%
% Re-run whenever the exports or animal_behaviour.mat change.
%
% lp_column selects which column of period_experienced is the learning point.
% It defaults to 1 (the column Tom's MATLAB pipeline uses, per HC_V1_temporal.m).
% period_experienced in these files has 10 columns -- the script prints all of
% them for the first animal so you can confirm column 1 is correct.

    if nargin < 1 || isempty(data_dir)
        here = fileparts(mfilename('fullpath'));      % .../cca/scripts
        data_dir = fullfile(here, '..', '..', 'HC_V1_data');
    end
    if nargin < 2 || isempty(lp_column)
        lp_column = 1;
    end
    data_dir = char(data_dir);
    assert(isfolder(data_dir), 'data_dir not found: %s', data_dir);

    % --- per-animal region labels ------------------------------------------
    files = dir(fullfile(data_dir, 'TF*_export.mat'));
    assert(~isempty(files), 'no TF*_export.mat files in %s', data_dir);
    fprintf('Region labels (%d animals):\n', numel(files));
    animals = struct('id', {}, 'file', {}, 'regions', {});
    for i = 1:numel(files)
        fn = files(i).name;
        S = load(fullfile(data_dir, fn), 'units');
        rl = cellstr(string(S.units.regions_label));   % `string` -> cellstr
        rl = reshape(rl, 1, []);                       % force a row
        tok = regexp(fn, '\d+', 'match', 'once');
        animals(i).id = str2double(tok);
        animals(i).file = fn;
        animals(i).regions = rl;
        fprintf('  %-22s id=%-4d  %s\n', fn, animals(i).id, strjoin(rl, ', '));
    end

    % --- cohort behaviour (learning points) --------------------------------
    behaviour = struct('id', {}, 'lp', {}, 'period_experienced', {});
    bpath = fullfile(data_dir, 'animal_behaviour.mat');
    if isfile(bpath)
        B = load(bpath);
        ids = local_ids_to_int(B.animal_id);
        pe = B.period_experienced;                     % animals x columns
        if size(pe, 1) ~= numel(ids) && size(pe, 2) == numel(ids)
            pe = pe.';                                 % orient: rows = animals
        end
        assert(size(pe, 1) == numel(ids), ...
            'period_experienced has %d rows but animal_id has %d entries', ...
            size(pe, 1), numel(ids));
        assert(lp_column >= 1 && lp_column <= size(pe, 2), ...
            'lp_column %d is outside 1..%d', lp_column, size(pe, 2));
        for i = 1:numel(ids)
            behaviour(i).id = ids(i);
            behaviour(i).lp = pe(i, lp_column);        % NaN -> JSON null
            behaviour(i).period_experienced = pe(i, :);
        end
        fprintf(['\nBehaviour: %d animals; period_experienced is %d x %d.\n' ...
                 'LP taken from column %d. First animal, all columns: %s\n'], ...
                numel(ids), size(pe, 1), size(pe, 2), lp_column, ...
                mat2str(pe(1, :)));
    else
        warning('animal_behaviour.mat not found -- behaviour block left empty');
    end

    % --- write the JSON companion file -------------------------------------
    out = struct();
    out.schema = 'tom_cca_labels_v1';
    out.generated = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    out.lp_column = lp_column;
    out.animals = animals;
    out.behaviour = behaviour;

    outpath = fullfile(data_dir, 'cca_labels.json');
    try
        txt = jsonencode(out, 'PrettyPrint', true);    % R2021a+
    catch
        txt = jsonencode(out);
    end
    fid = fopen(outpath, 'w');
    assert(fid ~= -1, 'could not open %s for writing', outpath);
    fwrite(fid, txt, 'char');
    fclose(fid);
    fprintf('\nWrote %s\n', outpath);
end


function ids = local_ids_to_int(aid)
% Normalise an animal_id field (string / char / cell / numeric) to integer ids
% by pulling the digits out of each entry.
    if isnumeric(aid)
        ids = double(aid(:)).';
        return;
    end
    s = reshape(string(aid), 1, []);
    ids = zeros(1, numel(s));
    for i = 1:numel(s)
        tok = regexp(char(s(i)), '\d+', 'match', 'once');
        if isempty(tok)
            ids(i) = NaN;
        else
            ids(i) = str2double(tok);
        end
    end
end
