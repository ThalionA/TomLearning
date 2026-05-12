# METHODS — HC↔V1 cross-region CCA pipelines

Detailed methodology for the spatial and temporal CCA analyses
implemented in `HC_V1_Code/`. Aimed at the level of a paper Methods
section: enough to replicate, with explicit parameter values, sample
gating, shuffle controls, and statistics. Last updated 2026-05-08.

## 1. Data

### Animals and recordings

Each animal contributes one `TF<id>_export.mat` file under
`HC_V1_data/`. Each export contains:

- `binned_spikes` : native-resolution (~1 ms) spike count matrix
  `[n_units × n_time_bins]`.
- `units.unit_id`, `units.regions_label`, `units.idx`, `units.idx_fs`:
  per-unit metadata. `units.idx` is a `[n_areas × n_units]` logical
  membership matrix indexed by `units.regions_label`. `units.idx_fs`
  flags fast-spiking units. **Per-unit `units.region` is unreliable
  for some areas (notably RSC) and is no longer read by either
  pipeline; see `v4_unit_regions.m`.**
- `data_behaviour.{velocity_gf, pos_binned, trial_binned_cued, ...}`:
  behavioural traces aligned to `binned_spikes`.
- `analysis_behaviour.masks.tunnel_cued` : logical mask, true on bins
  inside the cued portion of the corridor.
- `analysis_spatial.firing.cued.freq_z` : trial-averaged
  spatially-binned firing rate, z-scored. Used by spatial v2 only.

Animals' learning points are pulled from `animal_behaviour.mat`
(`period_experienced(:, 1)`). Non-learners get a yoked LP equal to the
mean LP of learners, so all animals contribute three epochs of equal
trial count.

### Region pairs

Both pipelines analyse the same eight pairs:
```
CA1-V1, CA1-DG, CA1-CA3, CA1-RSC, CA1-SUB, V1-RSC, RSC-SUB, CA3-DG
```

### Epochs

Per animal, three 10-trial epochs:
- `early` : trials 1..10 (naïve)
- `pre`   : trials (LP - 10)..(LP - 1) (just before learning)
- `post`  : trials LP..(LP + 9) (just after learning)

## 2. Pre-processing common to both pipelines

### Region tagging (`v4_unit_regions.m`)

Per-unit region label is built from `units.idx` × `units.regions_label`,
*not* `units.region`. For each row `r` of `units.idx`, all units with
`idx(r, :) == true` are assigned label `regions_label{r}`. Multi-row
membership has been 0/N in every animal probed; if it ever becomes
non-zero, a warning fires and the last-matching row wins.

Spatial v2 has always used this scheme. Temporal v4 was switched to it
on 2026-05-08; before that it read the per-unit `units.region` field,
which under-counted RSC severely (see `PLAN_v4_fixes.md`).

### Fast-spiking unit exclusion

For target areas `{CA1, V1, CA3, RSC}`, units with `idx_fs == true` are
excluded from PCA. Other areas (DG, SUB) keep all units. Identical in
both pipelines, matches v3 convention.

### Per-region PCA

Each region's z-scored unit-by-time matrix (temporal) or
unit-by-position matrix (spatial) is run through PCA. The number of
components retained, `k`, is selected by:

- **Temporal v4**: smallest `k` such that cumulative variance ≥ 90 %,
  floored at 3, capped at `max_k_per_region = 15`. Different `k` per
  region per animal.
- **Spatial v2**: `pca_selection_method = 'fixed'`, `k = 3` (current
  setting; was previously variance-driven). Same `k` for every region
  every animal. Floor `num_ccs_analyze = 1`, cap
  `max_k_per_region = 10`.

A region with fewer than `min_units_per_region = 5` units (after FS
exclusion) is dropped — no PCA fit, and any pair containing it is
skipped.

## 3. Temporal CCA v4 — `CCA_HC_V1_temporal_v4.m`

### Inputs

For each trial, region 1's PCA-reduced matrix `X` is `[n_valid_bins × k1]`
where each row is a 25 ms time bin's PC scores. Valid bins are those
where:

- `mask_cued == true` (animal in the cued corridor segment), and
- `velocity ≥ 2 cm/s` (animal is moving)

Region 2's `Y` is constructed identically.

### Per-trial CC1 (`v4_per_trial_cca`)

`canoncorr(X, Y)` is fit once per trial, returning the first canonical
correlation `r(1)` and the canonical loadings `A`, `B`. Skip rules:

- `n ≤ k1 + k2 + 4` → `info.reason = 'insufficient_samples'`
- `rank(X_centered) < k1` or `rank(Y_centered) < k2` → `'rank_deficient'`
- `canoncorr` throws → `'canoncorr_error'`

Output: one scalar per trial.

### Per-trial IFI (`v4_per_trial_ifi`, `v4_lag_corr`, `v4_ifi_from_lags`)

Project the trial's full timecourse onto the trial-level loadings:
`u = X_full * A`, `v = Y_full * B`. Split the trial into contiguous
runs of valid bins (`v4_contiguous_blocks`); discard any run shorter
than `min_block_bins = 12` (300 ms). For each surviving run, compute
lagged correlation between `u` and `v` at lags = -3:3 (= ±75 ms at
25 ms binning) using `v4_lag_corr`. Compute per-block IFI from the
lagged-correlation curve (see §6 below). Final per-trial IFI is the
length-weighted mean of per-block IFIs.

### RZ-pre per-bin CC1

For each trial, `v4_rz_align_pre` extracts the 20 bins (= 500 ms)
immediately before the first bin where `pos ≥ 200 cm` (reward zone
entry). The segment is rejected if position never crosses 200 cm, if
the segment would extend before bin 1, or if any bin in the segment is
invalid.

The retrieved segment (`[20 × k1]` and `[20 × k2]`) is *truncated to
`max_k_rz_per_region = 3`* columns before being stashed for the cross-
trial fit. The per-trial pre-RZ IFI computation (next subsection) uses
the full untruncated segment.

After the trial loop, for each (pair, epoch), the stashed segments are
stacked along trial axis: `X_stack` is `[20 × 3 × n_tr_rz]`,
`Y_stack` is `[20 × 3 × n_tr_rz]`. At each aligned bin `b`:

- `Xb = squeeze(X_stack(b, :, :))'` is `[n_tr_rz × 3]`.
- `Yb = squeeze(Y_stack(b, :, :))'` is `[n_tr_rz × 3]`.
- `v4_relaxed_canoncorr(Xb, Yb, rz_min_extra_samples = 2)` fits CCA
  with the relaxed margin `n > k1 + k2 + 2` instead of the default
  `n > k1 + k2 + 4`. Real CC1 is one scalar per aligned bin per
  (animal, pair, epoch); the per-bin shuffle null permutes trial labels
  on `Y` independently per bin (`n_shuffles = 20`).

The interpretable quantity here is **real - shuffle**, since at
`n_tr_rz ≤ 10` the absolute CC1 is biased substantially upward.

### RZ-pre per-trial IFI

For each trial with a valid pre-RZ segment, project onto the trial-
level `(A, B)`, compute lagged correlation on the 20-bin segment at
lags = -3:3, and apply the IFI formula. One scalar per trial.

### Shuffle nulls (per-trial)

- **CC1 null**: random permutation of `Y_valid` rows, refit
  `canoncorr`. `n_shuffles = 20`.
- **IFI null**: per-block row-shift on canonical projections —
  `v` is shifted relative to `u` by `±k` bins (`k` random in
  `[max(|lag|)+1, L - max(|lag|) - 2]`), then both are trimmed so the
  surviving overlap is the same length and NaN-free
  (`local_shift_align`). This preserves block autocorrelation while
  breaking alignment. Apply per-block, weight by block length, mean
  over `n_shuffles = 20`.

## 4. Spatial CCA v2 — `CCA_HC_V1_spatial_v2.m`

### Inputs

For each trial, region 1's matrix `D1` is built from
`analysis_spatial.firing.cued.freq_z` permuted to `[n_units × n_bins ×
n_trials]`, then PCA-reduced and reshaped to `[k1 × n_bins × n_trials]`
where `n_bins = 200` (positions of 2.5 cm each). For a single trial,
`D1(:, :, t)` is `[k1 × 200]`. Region 2 likewise.

### Per-trial CC1

For each trial, `canoncorr(D1(:, :, t)', D2(:, :, t)')` returns CC1.
The first canonical correlation is one scalar per trial.

### Per-trial precession_idx (`calc_precession`)

For shifts in `-max_shift_bins : max_shift_bins` (= -3:3 spatial bins =
±7.5 cm), shift `D1` relative to `D2` along the position axis, refit
`canoncorr` on the overlap, and record CC1 at each shift. The
precession_idx (the spatial analogue of IFI) is computed from the
seven-element lag-CC curve via the canonical bounded formula in §6.

### Per-spatial-bin sliding window

A separate analysis fits `canoncorr` in a sliding spatial window
(`n_bins_window = -3:3`) at each of the 200 position bins. This produces
a position-tuned CC1 trace per trial. Used in `cca_bin` outputs;
peripheral to the temporal-vs-spatial comparison.

### Shuffle nulls

For each trial, `D2` is randomly shifted along the trial axis by
`±max_shift_bins`, refit. `n_shuffles = 25`.

## 5. Region-tagging fix history

The temporal v4 pipeline read the per-unit `units.region` field until
2026-05-08. In 11 of 12 RSC-containing animals, that field reports
**zero** RSC units even when `units.idx` (cross-referenced against
`units.regions_label`) marks tens to hundreds of RSC units. Result: all
RSC-containing pairs were silently dropped from temporal v4 results.
Confirmed by `probe_units_and_v4.m` on 2026-05-08; fixed in
`v4_unit_regions.m` by switching to the `units.idx` × `units.regions_label`
scheme that spatial v2 has always used.

## 6. IFI / precession_idx definition

Both methods compute an asymmetry index from a 7-element lagged
correlation curve `r(lag)` for `lag ∈ {-3, -2, -1, 0, 1, 2, 3}`.

The canonical bounded form is
```
IFI = (|r_neg| - |r_pos|) / (|r_neg| + |r_pos|)
```
where `r_neg = mean(r(lag < 0))` and `r_pos = mean(r(lag > 0))`. This
is implemented in `v4_ifi_from_lags.m` (temporal) and `calc_precession`
in `CCA_HC_V1_spatial_v2.m` (spatial). Returns NaN if
`|r_neg| + |r_pos| < 1e-3`.

**Sign convention**: IFI > 0 ↔ region 2 leads region 1 (in time for
temporal, in space for spatial); IFI < 0 ↔ region 1 leads region 2.

**Earlier (buggy) form**, in use until 2026-05-08:
```
IFI_old = (r_neg - r_pos) / (r_neg + r_pos)
```
without absolute values, with an asymmetric guard
`(r_neg + r_pos) > 0.001` (temporal v4 / v3) or symmetric guard
`abs(r_neg + r_pos) > eps_tol` (spatial v2). Either form produces
`|IFI| > 1` when `r_neg` and `r_pos` have opposite signs (e.g.
`r_neg = 0.5, r_pos = -0.3` → IFI = 4.0). All v4 / spatial v2 outputs
generated before 2026-05-08 should be considered superseded.

The bounded form preserves sign in the all-positive regime (the common
case) but interprets opposite-sign cases as magnitude asymmetry rather
than signed correlation difference.

## 7. Temporal vs spatial comparison

### Pairing procedure (`compare_temporal_spatial.m`)

For each `(animal, pair, epoch)`:

- Spatial trial vector `trial_corr_<ep>{ia}` is positionally indexed:
  position `i` corresponds to trial id `idx_<ep>(i)` where
  `idx_early = 1:10`, `idx_pre = (lp-10):(lp-1)`,
  `idx_post = lp:(lp+9)`.
- Temporal trial vector `trial_cc1_pertrial_<ep>{ia}` is paired with
  `trial_id_pertrial_<ep>{ia}` (explicit trial ids).
- Trial ids are intersected; matched (real, real) pairs are appended
  to the long-form table along with shuffle counterparts.

The full long-form table is saved to
`HC_V1_data/compare_summary.mat`.

### Statistics

Per panel:
- **Spearman ρ**: rank-correlation between paired (spatial, temporal)
  values, robust to monotone transformations and outliers.
- **K-S two-sample test** (`kstest2`): in marginal-distribution panels,
  tests whether learner and non-learner distributions differ in any
  way (location or shape).
- **Wilcoxon rank-sum** (`ranksum`): tests location shift between
  learner and non-learner distributions.

No multiple-comparison correction is applied at the figure level —
significance marks are descriptive. For paper claims, apply
Bonferroni or FDR over the relevant family of tests.

## 8. Known parameter asymmetries between the two pipelines

As of 2026-05-08:

| Parameter            | Temporal v4 | Spatial v2 |
|----------------------|:-----------:|:----------:|
| `max_k_per_region`   | 15          | 10         |
| `pca_selection`      | variance ≥ 90% | fixed = 3  |
| `min_units_per_region` | 5         | 5          |
| `n_shuffles` (per-trial) | 20      | 25         |
| `max_k_rz_per_region` | 3 (RZ-only) | n/a       |

Effective per-region `k` ends up at 3–15 for temporal, exactly 3 for
spatial. This biases the comparison: for pairs where one region has
few useful PCs, the methods are roughly symmetric; where regions have
many useful PCs, temporal CCA has more flexibility and may report
higher CC1 (or, with small effective `n`, lower if dimensionality
overwhelms the sample-size guard).

## 9. Apples-to-apples control protocol

To rule out parameter-driven asymmetries, re-run both pipelines with:

- `max_k_per_region = 3`, `pca_selection_method = 'fixed'` in
  `CCA_HC_V1_temporal_v4.m` (a one-line patch in §1 SETUP and a
  one-line tweak in the per-region PCA block).
- `max_k_per_region = 3` in `CCA_HC_V1_spatial_v2.m` (already at fixed
  `k = 3` if `pca_selection_method = 'fixed'`).
- `n_shuffles = 25` in both.
- `min_units_per_region = 5` in both (already aligned).

Save under different output filenames
(`Temporal_CCA_v4_matched_<date>.mat`,
`Spatial_CCA_Results_matched_<date>.mat`) and run
`compare_temporal_spatial.m` against those, plus all four downstream
comparison scripts. If the qualitative pattern (RSC pairs spatial-
biased, CA3 pairs temporal-biased) survives, it's biological. If it
collapses, it was the PCA cap.

This control is queued, not yet run.

## 10. Files

### Source

- `HC_V1_Code/CCA_HC_V1_temporal_v4.m` — main temporal pipeline.
- `HC_V1_Code/CCA_HC_V1_spatial_v2.m` — main spatial pipeline.
- `HC_V1_Code/v4_*.m` — temporal v4 helpers (CCA, IFI, lag corr,
  contiguous blocks, RZ alignment, relaxed canoncorr, unit regions).
- `HC_V1_Code/compare_temporal_spatial.m` — pairing + Figures 1, 2.
- `HC_V1_Code/compare_excess_scatter.m` — Figures 3, 4
  (shuffle-corrected).
- `HC_V1_Code/compare_marginals.m` — Figure 5 (per-epoch marginals).
- `HC_V1_Code/compare_asymmetry.m` — Figure 6 (per-pair summary bars).

### Tests

- `tests/test_temporal_v4.m` — synthetic ground-truth tests for v4
  helpers. T1-T7 covering AR(1) coupling, uncoupled null, multi-block
  weighted IFI, RZ alignment, RSC tagging fix, relaxed canoncorr
  guard, IFI bound under all sign patterns.

### Outputs

- `HC_V1_data/Temporal_CCA_v4_<date>.mat`
- `HC_V1_data/Spatial_CCA_Results_<date>.mat`
- `HC_V1_data/compare_summary.mat` — long-form pairing table.
- `HC_V1_data/compare_asymmetry.mat` — Figure 6 numerical content.

### Figures

All under `HC_V1_figures/CCA_Compare/`:
- `TemporalVsSpatial_combined.svg` (Figure 1)
- `TemporalVsSpatial_per_pair.svg` (Figure 2)
- `TemporalVsSpatial_excess_combined.svg` (Figure 3)
- `TemporalVsSpatial_excess_per_pair.svg` (Figure 4)
- `TemporalVsSpatial_marginals_{early,pre,post}.svg` (Figure 5 ×3)
- `TemporalVsSpatial_asymmetry_summary.svg` (Figure 6)

### Documentation

- `HC_V1_Code/METHODS.md` — this file.
- `HC_V1_Code/FIGURE_LEGENDS.md` — paper-style legends for Figures 1-6.
- `HC_V1_Code/NOTES_temporal_v4.md` — v4 design log.
- `HC_V1_Code/NOTES_temporal_vs_spatial.md` — earlier comparison notes
  (will be retired once findings are folded into METHODS.md).
- `HC_V1_Code/PLAN_v4_fixes.md` — RSC + RZ fix plan, 2026-05-08.
- `HC_V1_Code/NOTES.md` — sanity-audit / simulation P0-P3 list.
