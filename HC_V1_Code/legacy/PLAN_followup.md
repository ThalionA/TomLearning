# Plan — three follow-up investigations after Fig 6

Triggered by `TemporalVsSpatial_asymmetry_summary.svg` (2026-05-08).
Order is from cheapest / most diagnostic to most invasive.

## Investigation 1 — visual inspection of Figure 3

**Question.** In Panel B of Fig 6, *almost every pair × epoch* is
negative (spatial excess CC1 > temporal excess CC1). Is this a uniform
trial-by-trial shift, or are a few outlier trials per pair pushing the
mean? Same magnitude either way, but very different interpretation.

**Action.** Open
`HC_V1_figures/CCA_Compare/TemporalVsSpatial_excess_combined.svg`. For
each of the six panels (CC1 × 3 epochs, IFI × 3 epochs), look at:

- The cloud shape: are points roughly aligned along *y = x* (methods
  agree on rank) but offset below the line (uniform spatial bias)? Or
  are most points on the diagonal with a tail of outliers below?
- The distribution along each axis: is spatial-excess CC1 uniformly
  higher across the full range, or only in a specific quadrant?
- Per-pair colour clusters: do the colours that were "spatial-biased"
  in Fig 6B form distinct clusters in the lower-right of the panels?

**Decision rule.**
- *Uniform shift* (cloud parallel to identity, offset below) → the
  spatial-bias is a real method asymmetry on every trial. Strong case
  for the apples-to-apples re-run.
- *Outlier-driven* (most points on identity, sparse cloud below) → the
  Fig 6B pattern is fragile. Investigate which animals/trials drive
  the outliers before drawing conclusions.

No code change. User opens the SVG, reports what they see, I interpret.

## Investigation 2 — apples-to-apples re-run

**Question.** Does the Panel B uniformity persist when both pipelines
use matched parameters? If yes → fundamental method asymmetry. If no →
the spatial PCA cap (k=3) vs temporal PCA cap (k=15) was driving it.

**Plan.**

1. Edit `CCA_HC_V1_temporal_v4.m` parameters (single block at top):
   - `pca_variance_threshold = 90;` (unchanged)
   - **`max_k_per_region = 3;`** *(was 15)*
   - **`min_units_per_region = 5;`** (unchanged, but verify)
   - **`n_shuffles = 25;`** *(was 20)*
   - **`save_path = fullfile(data_dir, sprintf('Temporal_CCA_v4_matched_%s.mat', current_date));`** *(was Temporal_CCA_v4_<date>)*
   - Also force the per-region PCA to use exactly k=3 floor (currently
     `max(3, k_pca)` which already gives k=3 when k_pca <3); cap at 3
     ensures it never exceeds.

2. Edit `CCA_HC_V1_spatial_v2.m` parameters:
   - **`n_shuffles = 25;`** (already 25)
   - `n_components_reduced = 3;` `pca_selection_method = 'fixed';` (already)
   - **`save_path = ...Spatial_CCA_Results_matched_...;`**

3. Run both pipelines.
4. Add a `MATCHED = true` flag to the four compare scripts so they
   auto-discover `*_matched_*.mat` instead of the default. Re-run them.
5. Compare the new Fig 6 to the original Fig 6.

**Cleanup.** After matched run, revert parameter changes (or keep them
on a branch / parallel set of files, since user prefers main).
Suggestion: copy the modified main scripts to
`CCA_HC_V1_temporal_v4_matched.m` and `_spatial_v2_matched.m`,
hard-coded with matched params. Doesn't require touching the originals.

I'll write the matched copies and a `compare_*_matched.m` wrapper set
when we get to this step.

## Investigation 3 — CA1-CA3 early case study

**Question.** Panel B of Fig 6 has exactly one positive bar:
CA1-CA3, early epoch, ≈ +0.08. Is this real fast-time-scale coupling
(theta / ripple) that temporal CCA detects but spatial doesn't, or a
fluke from a small `n`?

**Plan.** Write `case_study_CA1_CA3_early.m` that:

1. Loads `compare_summary.mat`.
2. Filters to CA1-CA3 early (one cell in the table).
3. Reports:
   - n animals contributing, n trials per animal.
   - Per-animal mean (real CC1, excess CC1, IFI) for both methods.
   - Pairwise temporal vs spatial scatter just for this cell, with
     animals labelled.
4. Optionally (more invasive, deferred unless useful): re-load the raw
   `TF*_export.mat` for a representative animal, recompute the lag-
   correlation curves (temporal and spatial) for one trial in the
   early epoch, and plot them side by side. This shows the actual
   shape of the cross-region timing structure rather than just the
   reduced IFI scalar.

**Decision rule.**
- If multiple animals contribute and trial-by-trial agreement is high
  for CC1 in this cell, the result is robust → CA1-CA3 early is a
  positive control for temporal CCA's added value.
- If 1-2 animals dominate and rest are NaN, the bar is noise → revisit
  whether n_animals × n_trials is sufficient anywhere.

## Suggested execution order

1. **Inv 1** (today): user looks at Fig 3, reports. I interpret.
2. **Inv 3** (next): write the CA1-CA3 case-study script. Runs in
   minutes. Adds clarity on whether temporal CCA has any unique
   advantage anywhere.
3. **Inv 2** (last): set up matched pipelines. Slowest, most
   invasive, but the definitive test of whether the panel B uniformity
   is methodological or fundamental.

Each step is independent — if inv 1 already gives a clear answer
(e.g., outlier-driven), inv 2 might become unnecessary.
