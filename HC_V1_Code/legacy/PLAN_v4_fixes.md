# Plan: v4 fixes (RSC tagging + RZ-per-bin universal failure)

Two independent bugs uncovered by `probe_units_and_v4.m` on 2026-05-08.

## Bug 1 — RSC tagging in v4 reads the wrong field

### Symptom
`Temporal_CCA_v4_2026_05_05.mat`: every RSC-containing pair (CA1-RSC,
V1-RSC, RSC-SUB) is 16/16 NaN across `trial_cc1_*`, `trial_ifi_*`,
`rz_*`. Spatial v2's results are not affected.

### Cause
`CCA_HC_V1_temporal_v4.m` line 164:
```matlab
if isfield(units, 'region')
    unit_regions = units.region(:);
elseif isfield(units, 'regions_label') && length(units.regions_label) == length(units.unit_id)
    unit_regions = units.regions_label(:);
else
    fprintf('  -> Could not find a valid per-unit region array. Skipping.\n'); continue;
end
```
This reads the per-unit `units.region` string array. For 11 of the 12
RSC-containing animals, that field assigns 0 RSC units even though
`units.idx` (a `[num_areas x num_units]` logical, addressed by
`units.regions_label`) marks 10–225 units in the RSC row. TF070 is the
only animal where `units.region` correctly tags RSC. Multi-membership
across idx rows is 0/N for every animal probed, so the idx scheme is a
one-hot encoding that is a strict superset of the per-unit `region`
field.

Spatial v2 doesn't have this bug because it indexes through
`units.idx` + `units.regions_label` (line 235 of
`CCA_HC_V1_spatial_v2.m`).

### Fix
Replace the region-tag block in `CCA_HC_V1_temporal_v4.m` (≈ lines
163-170) with the spatial-style derivation:

```matlab
% --- Per-unit region array (mirrors spatial_v2: build from idx + regions_label) ---
if ~isfield(units, 'idx') || ~isfield(units, 'regions_label')
    fprintf('  -> Missing units.idx or units.regions_label. Skipping.\n'); continue;
end
rl = cellstr(units.regions_label);
unit_regions = repmat({''}, length(units.unit_id), 1);
idx_logical = logical(units.idx);
% Sanity: every unit should belong to at most one row (probe confirms this).
multi = sum(idx_logical, 1) > 1;
if any(multi)
    fprintf('  -> %d/%d units in >1 idx row; using last-match assignment.\n', ...
        sum(multi), length(unit_regions));
end
for r = 1:numel(rl)
    in_region = idx_logical(r, :);
    unit_regions(in_region) = rl(r);
end
```

Notes:
- The `units.region` fallback is dropped entirely. We trust `units.idx`
  + `units.regions_label` as the single source of truth, matching
  spatial v2.
- A diagnostic line fires if any animal ever has multi-membership; in
  the current dataset this is 0/N for every animal but we want the warning
  in case a future export changes that.
- Region-name list (`animal_areas` derivation later in the script) should
  become `unique(rl)` rather than `unique(unit_regions)`. Otherwise empty
  strings can sneak in if some idx columns are all-false.

### Test
Add a new check to `tests/test_temporal_v4.m`:

> Build a fake `units` struct with `regions_label = {'CA1','RSC'}`,
> `idx = [1 1 0 0 0; 0 0 1 1 1]`, and `region = {'CA1';'CA1';'';'';''}`
> (mimicking the RSC-blind data). Run the new tagging block on it and
> assert `unit_regions == {'CA1';'CA1';'RSC';'RSC';'RSC'}`.

Add a second positive case with multi-membership to confirm the warning
fires and last-match assignment is consistent.

## Bug 2 — RZ-per-bin guard too tight given typical PC counts

### Symptom
`rz_cc1_*` is 16/16 NaN for every pair × every animal, including
non-RSC pairs. No `[rz-skip]` lines were printed because the
short-circuit at line 441 (`if any(sel_rz)`) eats the case where no
trials in the epoch produced an RZ stash.

### Cause
`CCA_HC_V1_temporal_v4.m` lines 447-455:
```matlab
if n_tr_rz > (k1 + k2 + rz_min_extra_samples)
    [cc_trace, cc_trace_sh] = local_rz_per_bin_cc( ...
        X_stack, Y_stack, n_shuffles, rz_min_extra_samples);
    ...
else
    fprintf('    [rz-skip] %s ep=%s: n_tr=%d, k1+k2=%d (need n>k1+k2+%d)\n', ...
```
With `n_trials_epoch = 10` (so n_tr_rz <= 10) and `max_k_per_region = 15`,
k1 + k2 typically lands at 6-10. The guard `10 > 8 + 2` is false, so
nearly every (pair, epoch) gets skipped. The internal helper
`local_relaxed_canoncorr` is wired correctly and would run if reached;
it's the OUTER guard that fails.

### Fix
Two complementary changes (both needed):

1. **Cap PC dimensionality on the RZ path.** Add a parameter
   `max_k_rz_per_region` (default 3) and apply it when stashing
   `rz_X_seg / rz_Y_seg`. Concretely: at line 366-367, pass a k-cap into
   `v4_rz_align_pre` (or truncate the returned segments), so each stash
   has at most `max_k_rz_per_region` columns. With cap = 3 and n_tr = 10,
   the guard becomes `10 > 6 + 2 = 8`, comfortably true.

2. **Make the silent-skip case audible.** Add a print when
   `~any(sel_rz)` (i.e. no stash for this epoch), so the run log
   distinguishes "no RZ trials available" from "RZ trials existed but
   guard failed":
   ```matlab
   if any(sel_rz)
       ...
   else
       fprintf('    [rz-empty] %s ep=%s: 0 RZ-stash trials in this epoch\n', ...
           group_results(ipair).pair_name, ep);
   end
   ```

Optional, lower priority:

3. **Consider lowering `max_k_per_region` overall to 10** to match
   spatial v2. This reduces the global PC count and makes the RZ guard
   easier to satisfy without a separate RZ-only cap. Trade-off: the
   per-trial CC for high-PC regions (e.g. V1 with 100+ units) could
   change. We have v3 + v4 baselines to compare against if we make this
   change, so it's reversible.

### Test
Add a new check to `tests/test_temporal_v4.m`:

> Generate synthetic RZ stash with n_tr = 10, k1 = k2 = 3, with a
> known coupling at one aligned bin. Run the RZ-per-bin path and assert
> the cc_trace has finite values at every bin (no NaN due to the guard)
> and that the coupled bin has higher real-minus-shuffle CC than the
> uncoupled bins.

## Implementation order

Both fixes are independent. Suggest in this order:

1. RSC tagging fix in `CCA_HC_V1_temporal_v4.m`.
2. New tagging test in `tests/test_temporal_v4.m`.
3. RZ k-cap parameter + stash truncation + silent-skip print.
4. New RZ test in `tests/test_temporal_v4.m`.
5. Run `tests/test_temporal_v4.m` locally; confirm 9/9 (was 7).
6. Re-run `CCA_HC_V1_temporal_v4.m` on real data → new
   `Temporal_CCA_v4_<today>.mat`.
7. Re-run `CCA_HC_V1_temporal_v4_plot.m` and inspect.

Once 1-7 are clean, the paired temporal-vs-spatial scatter (per-animal
× per-pair × per-trial) becomes meaningful. Building it before the
fixes would be misleading because the temporal side has no RSC data
and no RZ-per-bin data.

## Open questions

- Is there an upstream reason `units.region` is RSC-blind in 11/12
  animals? Worth flagging to whoever produced the exports — they may
  want to fix it at source so future scripts don't have to choose.
- Should we add `unique(units.region)` printout to the v4 run log so
  this kind of upstream tagging drift surfaces immediately?
- Default value for `max_k_rz_per_region`: 3 is conservative.
  Alternative: `min(3, max_k_per_region)` so it tracks the global cap.
  3 is simpler and matches the v3 floor (k_pca floor 3, cap
  max_k_per_region).
