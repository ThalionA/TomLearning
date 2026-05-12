function [unit_regions, region_names, info] = v4_unit_regions(units)
%V4_UNIT_REGIONS  Build a per-unit region label cellstr from a units struct.
%   [unit_regions, region_names, info] = v4_unit_regions(units)
%
%   Uses units.idx (a [num_areas x num_units] logical) addressed by
%   units.regions_label (a list of region-name strings) as the source of
%   truth. This mirrors what CCA_HC_V1_spatial_v2.m does and is more
%   reliable than reading units.region directly: in 11 of 12 RSC-containing
%   animals, units.region tags ZERO RSC units even though units.idx marks
%   tens-to-hundreds of them. See HC_V1_Code/PLAN_v4_fixes.md and the
%   probe_units_and_v4.m output for context.
%
%   Inputs:
%     units : struct with fields .unit_id, .regions_label, .idx
%
%   Outputs:
%     unit_regions : cellstr (n_units x 1). Empty string '' for any unit
%                    that does not appear in any idx row.
%     region_names : cellstr of region names from units.regions_label.
%     info         : struct with diagnostic fields:
%                      .n_units
%                      .n_regions
%                      .n_unassigned    : units with no idx-row membership
%                      .n_multi         : units with membership in >1 idx row
%                                         (we use last-match assignment; warn
%                                         if any).
%
%   Errors out cleanly if required fields are missing.

    if ~isfield(units, 'unit_id')
        error('v4_unit_regions:missing_unit_id', 'units.unit_id is required.');
    end
    if ~isfield(units, 'regions_label')
        error('v4_unit_regions:missing_regions_label', 'units.regions_label is required.');
    end
    if ~isfield(units, 'idx')
        error('v4_unit_regions:missing_idx', 'units.idx is required.');
    end

    n_units      = length(units.unit_id);
    region_names = cellstr(units.regions_label);
    region_names = region_names(:);                 % column
    n_regions    = numel(region_names);

    idx_logical = logical(units.idx);
    if size(idx_logical, 1) ~= n_regions
        error('v4_unit_regions:size_mismatch', ...
            'units.idx has %d rows but regions_label has %d entries.', ...
            size(idx_logical, 1), n_regions);
    end
    if size(idx_logical, 2) ~= n_units
        error('v4_unit_regions:size_mismatch', ...
            'units.idx has %d columns but unit_id has %d entries.', ...
            size(idx_logical, 2), n_units);
    end

    unit_regions = repmat({''}, n_units, 1);
    col_sum = sum(idx_logical, 1);

    % Last-match assignment: if a unit happens to appear in multiple idx
    % rows, the last loop iteration wins. Probe shows multi-membership is
    % 0/N for every animal we have, so this is normally a strict one-hot.
    for r = 1:n_regions
        unit_regions(idx_logical(r, :)) = region_names(r);
    end

    info = struct();
    info.n_units       = n_units;
    info.n_regions     = n_regions;
    info.n_unassigned  = sum(col_sum == 0);
    info.n_multi       = sum(col_sum > 1);
end
