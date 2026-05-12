# Temporal CCA v4 vs Spatial CCA v2 — methodology and how to read the comparison

## What each method estimates

Both fit canonical correlation between region-1 and region-2 PCA-reduced
population activity, but the input matrix is different.

### Temporal CCA v4 (per trial)

For one trial, `X` is `[n_valid_bins x k1]` where each row is a 25 ms
time bin's region-1 PC scores. `Y` is the matching `[n_valid_bins x k2]`
for region 2. `canoncorr(X, Y)` finds linear combinations of region-1
PCs and region-2 PCs whose **moment-to-moment time courses** are
maximally correlated across the trial. The first canonical correlation
CC1 is one scalar per trial.

The per-trial IFI then projects the full trial onto that trial's
canonical loadings (A, B), computes lagged correlation between u(t)
and v(t) at lags = -3:3 bins (= ±75 ms), and reports
`(r_neg − r_pos) / (r_neg + r_pos)`. Sign convention: IFI < 0 when
region 1 (u) leads region 2 (v) in time.

### Spatial CCA v2 (per trial)

For one trial, `X` is `[200 x k1]` where each row is a 2.5 cm position
bin's region-1 PC scores (the population firing-rate vector at that
position on the linearised track). `Y` is `[200 x k2]` for region 2.
`canoncorr(X, Y)` finds linear combinations whose **position-tuning
profiles** along the track are maximally correlated. CC1 is one scalar
per trial.

The spatial precession index is the spatial analogue of IFI: lagged
correlation between u(position) and v(position) at lags = -3:3 position
bins (= ±7.5 cm). Sign convention: precession_idx < 0 when region 1
leads region 2 in space.

## When the methods agree, when they disagree

Within a corridor-running trial, position is roughly a monotonic
function of time, so the two analyses see closely related data
matrices. They will tend to agree, but they answer different questions
and can disentangle in informative ways.

| What's true about a pair                                    | Temporal CC1 | Spatial CC1 |
|-------------------------------------------------------------|:------------:|:-----------:|
| Shared **moment-to-moment** state (arousal, locomotion ramp), no co-tuning | high   | low |
| Shared **position-tuned ensembles** that drift trial-to-trial in timing  | low    | high |
| Both: place fields aligned AND co-fluctuating in time       | high   | high |
| Neither: independent populations                            | ≈ shuffle | ≈ shuffle |

For lead/lag (IFI vs precession_idx):

- **Same sign** is expected when within-trial speed is approximately
  stationary — being ahead in time is equivalent to being ahead in
  position, up to a velocity factor.
- **Opposite sign** is a red flag. It can happen if the same neural lead
  manifests as a position-shift but not a time-shift (or vice versa)
  due to non-stationary speed or asymmetric trial structure (slow at
  start, fast at end).
- **Magnitudes are not directly comparable.** Temporal IFI is unitless
  but defined at ±3 bins of 25 ms; spatial IFI is unitless at ±3 bins of
  2.5 cm. Comparing magnitudes across methods conflates the binning.

## How to read the paired scatter

For each (animal, pair, epoch, trial_id) we extract the two scalars
from each method and pair them by trial_id. The combined scatter
(`TemporalVsSpatial_combined.svg`) overlays all pairs in six panels
(2 metrics × 3 epochs). The per-pair scatter
(`TemporalVsSpatial_per_pair.svg`) breaks it down per pair × epoch with
a learner / non-learner split.

Read each panel as:

- **High Spearman ρ, near identity line:** the methods agree both
  rank-wise and absolutely. Most reassuring case — a trial that scores
  high on temporal CC1 also scores high on spatial CC1.
- **High Spearman ρ, off identity:** rank-wise agreement but a constant
  offset or scale difference. Often expected for IFI (binning).
- **Low Spearman ρ, no clear scatter pattern:** the metrics are
  measuring different things on these data. Worth investigating with the
  decision tree above.
- **Anti-correlated (ρ < 0):** something to dig into. Either a sign-
  convention mismatch, non-stationary speed, or a real
  spatial-vs-temporal dissociation.

## What the comparison **cannot** tell us

- It does not prove the pipeline is "correct" — only that the two
  views of the data are mutually consistent (or not). A bug in both
  pipelines (same direction) will show as agreement.
- It does not capture the **shuffle null**. The interpretable quantity
  for absolute CC1 in the small-sample regime is real - shuffle, not
  real alone. The paired scatter shows real CC1 vs real CC1; we can
  optionally add shuffle-corrected versions in a second pass.
- It does not handle epoch-mean comparisons, which would tell us if
  the two methods agree on **trial-averaged** structure for an animal.
  The current script is per-trial only.

## Known asymmetries between the two pipelines (as of 2026-05-08)

- **PCA cap.** Temporal v4: floor 3, cap `max_k_per_region = 15`,
  threshold 90% variance. Spatial v2 (current working tree): fixed k = 3
  (floor 3, cap 10, threshold 90% — but `pca_selection_method = 'fixed'`
  selects exactly 3). So spatial uses fewer PCs by construction in the
  current parameterization. Effect: spatial CC1 may be slightly lower
  than it would be with k = 15, but is more reproducible across animals.
- **n_shuffles.** Temporal: 20. Spatial: 25.
- **min_units_per_region.** Both 5 (after the latest spatial edit).
- **RSC tagging.** v4 reads `units.idx` + `units.regions_label` (after
  the 2026-05-08 fix in `v4_unit_regions.m`). Spatial v2 always did.
  Now equivalent.

If we want a strict apples-to-apples comparison rather than "what each
method's chosen parameterization says", we'd run both with matched
`max_k_per_region` (e.g. cap 3 for both) and matched n_shuffles. That's
a separate experiment worth queuing if the current comparison shows
disagreement.

## Files

- `HC_V1_Code/compare_temporal_spatial.m` — the comparison script.
- `HC_V1_data/compare_summary.mat` — long-form pairing table, regenerated
  on each run.
- `HC_V1_figures/CCA_Compare/TemporalVsSpatial_combined.svg`
- `HC_V1_figures/CCA_Compare/TemporalVsSpatial_per_pair.svg`
