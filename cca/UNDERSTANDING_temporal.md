# Tom-learning CCA -- Temporal arms (Addendum)

Status: **design ratified 2026-05-26; scaffold + tests pending**
Companion to `cca/UNDERSTANDING.md` (the spatial arm).

---

## 1. Scope

Two new analysis arms, both Tom-only for now, both at 50 ms binning:

* **Arm A -- running-state full-traversal.** Full corridor traversals re-binned
  to 50 ms. Velocity-filtered (>= 2 cm/s) per bin; the resulting fragmented
  timeseries is segmented into contiguous valid runs and the lag analysis is
  segment-restricted. Signal-only CCA (ragged trial lengths preclude residualisation).
* **Arm B -- landmark-centred per-landmark.** 500 ms windows centred on each
  landmark entry event (bins -5..+4, entry at index 5). Per-window engagement
  gate (mean velocity >= 2 cm/s). Fit per `(animal, pair, epoch, landmark_id)`
  -- six landmark identities, each its own fit cell. Residual CCA is the
  committed default (PSTH per `(bin, unit)` within each landmark's slice);
  signal CCA stays on the sweep grid.

The two arms answer different scientific questions and coexist; neither
supersedes the other and neither supersedes the spatial arm.

---

## 2. Why two arms

| Question | Arm |
|---|---|
| Is there moment-to-moment HC<->cortex communication during running, distinct from co-modulation by shared evoked structure? | **A** -- full-traversal signal CCA on the residual fluctuations *cannot* be formed (no per-time-bin PSTH across ragged trials), so the signal-CCA reading carries both evoked and trial-to-trial structure. The question it answers is therefore necessarily about total moment-to-moment co-variation. |
| Is there event-locked HC<->cortex communication around landmarks, above and beyond what each area's shared landmark response can explain? | **B** -- aligned windows make a per-`(bin, unit)` PSTH well-defined, so residual CCA isolates trial-to-trial fluctuations around the landmark response. This is the temporal analogue of the spatial arm's "communication" definition. |

The spatial arm answers a third question (position-tuning-aligned communication
across the corridor) and is untouched.

---

## 3. Data contract

Per-animal `TF*_export.mat` (h5py-readable; verified on TF028 2026-05-26):

* `binned_spikes` `(n_units, n_1ms_bins)` uint8 -- 1 ms spike counts on a flat
  timeline. Re-binned to 50 ms by summation; an incomplete final bin is dropped.
* `data_behaviour.velocity_binned_gf` `(n_1ms_bins, 1)` float64 -- gap-filled
  velocity at 1 ms. Aggregated to 50 ms by mean.
* `data_behaviour.trial_binned_cued` `(n_1ms_bins, 1)` float64 -- per-1ms-bin
  cued trial index; NaN outside the cued task. Aggregated to 50 ms by mode
  (segment boundary if mode is NaN).
* `data_behaviour.visual_landmark_binned` `(n_1ms_bins, 1)` float64 -- landmark
  identity 1..6, 0 outside any landmark zone. Aggregated to 50 ms by the first
  non-zero value within the 1 ms span.
* `units.idx` + `cca_labels.json` -- area-membership per unit (unchanged from
  spatial path).

Probe of TF028: 7.92M 1ms bins, 137 cued trials, 6 landmark IDs, ~8.7
landmark entries per cued trial. Landmark identity is consistent across the
cohort (verification deferred to the smoke test).

---

## 4. Locked design parameters (both arms)

| Parameter | Value | Notes |
|---|---|---|
| `bin_ms` | 50 | Locked, not swept. Tom's MATLAB v6 default; departure from v5 cross-project lock at 25 ms is deliberate and noted in the methods writeup. |
| `min_units` | 5 per area | Existing spatial-arm value. |
| `velocity_thresh_cm_s` | 2.0 | Threshold applied to 50 ms-binned mean velocity. |
| `n_shuffles` (sweep) | 100 | |
| `n_shuffles` (committed) | 200 | |
| CV | 10-fold leave-one-trial-out | Trials are the natural exchangeable unit; within-trial segments / crossings are autocorrelated and must stay together in a fold. |
| `cca_type` (Arm A) | signal-only | Forced. No residualisation possible with ragged trial lengths. |
| `cca_type` (Arm B) | residual + signal | Both on sweep grid; residual is committed default. |

---

## 5. Arm A -- running-state full-traversal

### Data preparation

1. Bin `binned_spikes` to 50 ms per cued trial. Trials with fewer than `max_lag_bins + 5` valid bins are dropped at the diagnostic stage, not the loader stage, so the skip-reason is recorded.
2. Aggregate `velocity_binned_gf` to 50 ms (mean).
3. Per-unit z-score over all cued+running 1 ms bins **before** rebinning, so the scale matches what the analysis sees. Whole-engaged-period reference window (matches spatial-arm convention).

### Velocity mask + segments

1. Threshold the 50 ms binned velocity at 2 cm/s -> raw boolean mask.
2. Morphological close: fill any invalid run of length <= `min_gap_close_bins = 2` (= 100 ms). Removes threshold-flicker artefacts.
3. Morphological open: drop any valid run of length < `min_run_bins = max_lag_bins + 5 = 15` (= 750 ms). Ensures every retained segment yields >= 5 paired samples at the most extreme lag.
4. Segments are maximal contiguous valid runs **within a single cued trial**. A segment never crosses a trial boundary even if velocity stays high across the inter-trial gap (the gap is a teleport/reset, not contiguous neural data).

### Fit unit and CCA

* Fit unit: `(animal, pair, epoch)`. Three epochs, same definition as spatial: naive = trials 1..10, intermediate = (LP-9)..LP, expert = (LP+1)..(LP+10).
* Activity matrix per area: vertical stack of all retained 50 ms bins across all segments in the epoch -- shape `(n_valid_bins_in_epoch, n_units_in_area)`.
* CCA: signal-only (no residualisation), PCA-reduced per area, k chosen by `k_mode` (sweep axis), symmetric within pair, fixed across epochs.

### Lag analysis (IFI_lagged, segment-aware)

* Lag scan: integer bin lags in `[-10, +10]` (= +-500 ms at 50 ms).
* At each lag L, build the design matrix by pooling within-segment lag-L pairs across all segments: a segment of length n contributes n - |L| pairs. No pair ever straddles a gap or trial boundary.
* Refit CCA on the pooled per-lag design (IFI_lagged, H&H-style). Held-out CC via the same 10-fold leave-one-trial-out CV (the trial label is carried with each pair).
* IFI from the per-dim lag curve using the existing `information_flow_index` formula (clipped at 0, normalised). Reported per dimension; headline = dim 1.

### Surrogate null (segment-aware circshift)

* For each surrogate, shift each segment's Y rows independently by a random amount drawn from `[1, segment_length - 1]` and wrapped within the segment.
* Refit CCA at each lag per surrogate; per-dim null distribution.
* The within-segment shift is essential: a global circshift across the
  concatenated array imports the gap structure into the null and biases the
  significance test (a v4-era silent bug, see `legacy/NOTES_temporal_v4.md`).

### Sweep grid + committed config

* Sweep: 11 configs over the k-rule axis (samples 15/25/40, fixed 3/5/10/20/30, variance 75/85/95). Signal-only is forced; bin width is locked at 50; lag range is locked at +-500 ms. Tag: `temp50_{krule}`.
* Committed: `temp50_samp15`. Mirrors the spatial committed (`res_samp15_lag10`) on the PC-count axis.

---

## 6. Arm B -- landmark-centred per-landmark

### Landmark entry detection

* Walk `visual_landmark_binned` at 1 ms; an "entry" is any 1 ms bin where the value transitions from 0 to a non-zero landmark ID i. Repeat entries to the same landmark within one trial are kept separately (each is its own crossing).
* Entry time -> entry bin in the 50 ms grid (the 50 ms bin containing the 1 ms entry index).

### Window definition

* Window: 10 bins centred on entry -- 5 bins pre, entry bin at index 5, 4 bins post. Spans -250 ms to +200 ms relative to entry (asymmetric by one bin; the alternative -200/+250 is equivalent under sign convention).
* A window is dropped if it overruns the cued task boundary or the trial boundary at either end.
* Per-window engagement gate: window is dropped if mean velocity over its 10 bins is < 2 cm/s. Per-bin fragmentation inside windows is not handled; windows are all-or-nothing.

### Fit unit and CCA

* Fit unit: `(animal, pair, epoch, landmark_id)`. With 6 landmarks, that is 6x the per-epoch fit count of the spatial / Arm A path.
* Per-fit-cell activity per area: rectangular `(n_kept_crossings, 10, n_units)` tensor.
* Residualisation (committed default): subtract the per-`(bin, unit)` trial-mean across the kept crossings, within this landmark's own slice. Signal-CCA variant available on the sweep axis.
* `min_crossings_per_landmark = 5` per fit cell. Cells below the floor are skipped with reason `low_crossings`.
* k chosen by `k_mode` (sweep axis), capped by smaller area unit count and `k_cap = 30`; the per-cell sample count (~145) keeps `k_samples = floor(145/15) = 9` above the typical unit-count cap for most pairs.

### Lag analysis (IFI_lagged, window-bounded)

* Lag scan: integer bin lags in `[-5, +5]` (= +-250 ms at 50 ms). Bounded by window length: a 10-bin window with the >= 5 paired-sample requirement caps |L| at 5 bins.
* At each lag L, slice within each kept window to get the (n_kept_crossings, 10 - |L|, k) sub-tensor for X and Y; pool the (X-bin, Y-bin+L) pairs across all kept crossings. Refit CCA on the pooled per-lag design. Held-out CC via 10-fold leave-one-trial-out (a held-out trial removes all its crossings).
* IFI per dim from the per-dim lag curve, same formula as spatial / Arm A.

### Surrogate null (per-window circshift)

* For each surrogate, shift each window's Y rows independently by a random amount in `[1, 9]` (>= 1 bin shift; 10-bin window). Refit CCA at each lag.

### Sweep grid + committed config

* Sweep: 22 configs over (cca_type x k-rule) = 2 x 11. Tag: `landmark50_{cca}_{krule}`.
* Committed: `landmark50_res_samp15`.

### Skip-reason diagnostics

Every fit cell records one of:

* `ok` -- fit completed
* `low_crossings` -- < 5 kept crossings after the velocity gate
* `low_units` -- one or both areas have < `min_units` after FS exclusion
* `rank_deficient` -- per-epoch PCA produces fewer dims than needed for the requested k
* `samples_per_pc_not_met` -- when `k_mode = "samples"` and `n_samples < samples_per_pc`

Per-`(animal, pair, epoch, landmark_id)` coverage table written to the result pkl and rolled up into the sweep summary.

---

## 7. Diagnostics per fit cell (both arms)

Carried on every result regardless of arm:

* `n_samples_used`, `k1`, `k2`, `samples_per_dim`
* Mean spikes/bin per area, fraction of zero bins per area (sparsity check)
* Mean and SD of retained-bin velocity per epoch (Arm A) or mean window-velocity per epoch (Arm B), for the epoch-confound check
* Arm A: `n_segments`, segment-length histogram (`[<150, 150-300, 300-600, >600]` ms buckets), fraction of bins lost to the velocity filter, `n_paired_samples` per lag.
* Arm B: `n_kept_crossings` per `(epoch, landmark_id)`, per-fit skip-reason.

---

## 8. Architecture

Code lives in `cca/src/tom_cca/`. Spatial path is untouched.

* `config.py` -- new fields: `bin_mode` (`"spatial" | "temporal_runstate" | "landmark"`), `temporal_bin_ms = 50`, `lag_ms` (replaces `max_lag_bins` for temporal arms), `velocity_thresh_cm_s = 2.0`, `min_gap_close_ms = 100`, `min_run_extra_ms = 250`, `landmark_window_bins_pre = 5`, `landmark_window_bins_post = 5`, `min_crossings_per_landmark = 5`.
* `dataio.py` -- new per-animal cache for the 1 ms triple `(binned_spikes, velocity_binned_gf, trial_binned_cued, visual_landmark_binned)`; helpers `_rebin_spikes`, `_bin_velocity`, `_bin_trial_index`, `_bin_landmark_id`, `_velocity_mask`, `_landmark_entries`.
* `segments.py` (new) -- `find_segments(mask_50ms, trial_idx_50ms, min_gap_close_bins, min_run_bins) -> list[Segment]`. Pure, fully testable.
* `landmark_align.py` (new) -- `align_windows(rebinned_spikes, entries, window_bins_pre, window_bins_post) -> dict[landmark_id, (n_crossings, n_bins, n_units)]`, plus the per-window mean-velocity engagement gate.
* `lagged_temporal.py` (new) -- segment-aware `lag_pairs`, segment-pooled `lag_curve`, segment-aware circshift surrogate.
* `lagged_landmark.py` (new) -- window-bounded `lag_pairs`, window-pooled `lag_curve`, per-window circshift surrogate.
* `pipeline.py` -- branches on `cfg.bin_mode` only inside `prepare_pair` and `fit_pair`. Spatial branch byte-identical to the current code.
* `analysis.py` -- Arm B loops over `landmark_id` within each `(animal, pair, epoch)` cell.
* `sweep.py` -- `build_sweep("temporal_runstate")`, `build_sweep("temporal_landmark")` -- new builders; spatial builder unchanged.

The spatial `lagged.py` is **not modified**. Per the `feedback_v4_conventions` rule (never modify prior versions for incompatible extensions), the temporal lag machinery lives in its own modules.

---

## 9. Tests (TDD, before any implementation)

Per CLAUDE.md TDD policy. Tests on synthetic ground truth precede the implementations they validate.

* `tests/test_segments.py` -- morphological cleanup idempotence; no segment shorter than `min_run` survives; segment never crosses a trial boundary.
* `tests/test_landmark_align.py` -- correct entry detection (0->n_id transitions only, not n_id_a->n_id_b); window slicing bounds; per-window velocity-gate behaviour.
* `tests/test_lagged_temporal.py` -- known-lag recovery on a fragmented synthetic signal (inject a 75 ms lead, sprinkle gaps, confirm peak lag = +75 ms / IFI sign correct); no lag-pair ever crosses a gap or trial boundary; segment-aware circshift null is symmetric around zero IFI under no-coupling synthetic data.
* `tests/test_lagged_landmark.py` -- known-lag recovery on aligned synthetic windows; residualisation removes an injected PSTH (post-residual CC drops by the expected amount); `min_crossings_per_landmark` skip behaviour; per-window circshift null symmetric.
* `tests/test_pipeline_bin_mode.py` -- spatial path produces byte-identical results before and after the bin_mode addition (regression guard).

---

## 10. Won't-do this round

* **Bin-width sweep.** Locked at 50 ms; v5 lock-at-25 ms tension acknowledged but Tom's existing MATLAB-v6 alignment wins for now.
* **Striatum port.** Both arms are Tom-only for this iteration. Once results land, the landmark arm may or may not map to Striatum's task structure (the landmark concept is HC-corridor-specific).
* **Across-landmark pooling fallback for Arm B.** When a per-landmark cell is below the crossings floor, it is skipped with `low_crossings` -- not pooled with adjacent landmarks. The diagnostics table makes coverage visible.
* **Residual CCA for Arm A.** Impossible -- ragged trial lengths preclude a per-(bin, unit) PSTH.
* **Per-bin velocity filtering inside Arm B windows.** Per-window all-or-nothing engagement gate only.
* **GPU backend.** CPU is fine.

---

## 11. Open items

* Confirm landmark ID consistency across the cohort (TF028 has 1..6; verified
  during the smoke test).
* `min_gap_close_ms` and `min_run_extra_ms` may need tuning after the first
  real-data smoke -- review the per-animal segment-length histogram.
* Sweep execution is a multi-hour native job. The Cowork sandbox caps shell
  commands at 45 s and kills background jobs, so the sweep is handed off to
  Theo to run natively after the pipeline passes its smoke test.
* Methods doc (`ResearchVault/Methods/CCA_HH_Adapted.md`) needs a Tom-50ms
  departure note added by hand -- not mounted in the container.

---

## 12. Provenance

Decisions ratified 2026-05-26 in the design session that produced this doc.
Decisions log (chronological, this session):

* Q1 bin width: 10 / 25 / 50 ms sweep -> revised to 50 ms locked.
* Q2 velocity mask: morphological cleanup (close 100 ms, drop short runs).
* Q3 IFI flavour: refit CCA per lag (IFI_lagged, H&H-style).
* Q4 epoch confound: diagnostics only, no sample-count matching.
* Q5 lag range: +- 500 ms (Arm A); +- 250 ms forced by window length (Arm B).
* Q6 n_shuffles: 100 sweep / 200 committed.
* Q7 Arm coexistence: both arms run, neither supersedes the other.
* Q8 Landmark window: 500 ms (10 bins), per-landmark fits.
* Q9 Sample budget: ~145 samples per fit confirmed sufficient (k bound by
  unit count, not sample count).
