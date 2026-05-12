# NOTES — Temporal CCA v4 (session 2026-05-05)

Branch: `temporal-v4` (branched off `main`). Files added in this session, no
existing files modified.

## What v4 is

Per-trial Temporal CCA + RZ-pre per-bin Temporal CCA. Departs from v3's
sliding within-trial design:

- **Per-trial CC**: one canoncorr fit per trial, on all valid bins of that
  trial (>=2 cm/s, cued). One scalar per trial.
- **Per-trial IFI**: project the trial's full timecourse onto the
  trial-level (A, B); split into contiguous valid blocks; per-block IFI;
  weighted-mean by block length.
- **RZ-pre per-bin CC**: at each timebin in the last 500 ms before RZ
  entry, fit canoncorr across the epoch's trials (n = trials_in_epoch,
  k1 + k2 up to 8). Real - shuffle is the interpretable quantity here.
- **RZ-pre IFI**: per-trial scalar computed on the pre-RZ segment only,
  using trial-level (A, B). Epoch-mean per animal.

## Design decisions (locked in this session)

1. **CCA basis for IFI projections**: trial-level (A, B). Same fit for
   per-trial IFI and RZ-pre IFI. Per-block re-fits rejected as
   variance-inflating.
2. **Lag range**: `lags = -3:3` at 25 ms binning -> +-75 ms. Same as v3.
3. **Min block length for IFI**: 12 bins (300 ms). Drops short noisy
   blocks where lag corr is unstable.
4. **RZ window**: pre-only, `last 500 ms`. Segment ends one bin BEFORE
   RZ entry. Animals stop to consume reward at RZ entry, so the post-RZ
   side is contaminated.
5. **PCA basis**: global per-region PCA on all valid session bins, k from
   90% var floor 3 cap 4. Identical to v3.
6. **Shuffles**: 20 per trial. Real-CC null = random permutation of Y
   rows within trial. IFI null = per-block row-shift on canonical
   projections (preserves block autocorrelation, breaks alignment).

## Files added

```
HC_V1_Code/
  v4_per_trial_cca.m          single canoncorr per trial w/ rank+sample guards
  v4_per_trial_ifi.m          per-block IFI, length-weighted mean
  v4_rz_align_pre.m           extract last n_rz_bins before pos>=rz_entry_cm
  v4_contiguous_blocks.m      runs of true >= min_block, returns Nx2 [start stop]
  v4_lag_corr.m               lagged correlation (promoted from v3 local)
  v4_ifi_from_lags.m          IFI = (r_neg - r_pos) / (r_neg + r_pos)
  CCA_HC_V1_temporal_v4.m     main pipeline, mirrors v3 loading/PCA
  CCA_HC_V1_temporal_v4_plot.m  six figures with rmANOVA on bars

tests/
  test_temporal_v4.m          MATLAB synthetic-data tests (run locally)
  test_temporal_v4_proto.py   Python shadow used to validate algorithm
                              (passes 4/4 ground-truth scenarios)
```

## Outputs (group_results fields, per pair x epoch)

| Field                    | Shape                          | Notes |
|--------------------------|--------------------------------|-------|
| `trial_cc1_<ep>`         | n_animals x 1                  | Epoch mean of per-trial CC |
| `trial_cc1_sh_<ep>`      | n_animals x 1                  | Same, shuffle null         |
| `trial_ifi_<ep>`         | n_animals x 1                  | Epoch mean of weighted IFI |
| `trial_ifi_sh_<ep>`      | n_animals x 1                  | Same, block-row-shift null |
| `rz_cc1_<ep>`            | n_animals x n_rz_bins (=20)    | Per-aligned-bin CC trace   |
| `rz_cc1_sh_<ep>`         | n_animals x n_rz_bins          | Same, trial-permutation null |
| `rz_ifi_<ep>`            | n_animals x 1                  | Pre-RZ IFI scalar          |
| `rz_ifi_sh_<ep>`         | n_animals x 1                  | Same, shift-trim null      |

Save file: `Temporal_CCA_v4_<date>.mat`. Partial save after each animal.

## Gotchas

- **Small-sample tension at RZ-per-bin CC.** n = trials_in_epoch <= 10
  vs k1 + k2 up to 8. The default `n > k1+k2+4` guard would skip every
  bin. Mitigated by a relaxed guard (`rz_min_extra_samples = 2`, knob in
  the main script) and by reporting **real - shuffle** as the primary
  quantity (the absolute CC is biased high). If most pairs/epochs end up
  printing `[rz-skip]`, lower `max_k_per_region` for the RZ path or
  enlarge the epoch.
- **Canonical-projection shuffle and NaN padding.** First draft used a
  NaN-padded row-shift; v4_lag_corr's `std(b) > 0` check returned false
  on NaN-containing vectors and silently produced all-NaN IFIs. Fixed
  with `local_shift_align`: shifts v relative to u and trims BOTH so the
  surviving overlap is the same length and NaN-free.
- **Cell-field initialisation in `repmat(struct(...))`.** Need
  `'rz_X', {{}}` (double braces) so each replicated struct gets its own
  empty cell array, not 0 fields. Otherwise `pair_cache(ipair).rz_X{end+1}`
  errors.

## Validation

`tests/test_temporal_v4_proto.py` exercises four synthetic scenarios:

1. Lagged AR(1) coupling (lag=+1) -> per-trial CC = 0.579, IFI = -0.482
   (expect IFI < 0 since u leads v).
2. Uncoupled regions -> real CC ~ shuffle CC (within 0.10).
3. Two valid blocks (15 + 180 bins) -> weighted IFI = -0.410 (long block
   dominates). Setting min_block_bins=16 drops the short block and
   leaves IFI = -0.351 from the long block alone.
4. RZ-locked coupling (10-bin burst before RZ entry) -> real-minus-shuffle
   CC higher in last 5 aligned bins than first 5; pre-RZ IFI = -0.209.

All four pass. Same scenarios are encoded in
`tests/test_temporal_v4.m` for local MATLAB execution before running on
real data.

## Open / next steps

- [ ] Run `tests/test_temporal_v4.m` locally; assert all 7 checks pass.
- [ ] Run `CCA_HC_V1_temporal_v4.m` on real data.
- [ ] Run `CCA_HC_V1_temporal_v4_plot.m` and inspect figures.
- [ ] Decide whether to keep `rz_min_extra_samples = 2` or tighten further;
      depends on how many `[rz-skip]` lines appear in the run log.
- [ ] Commit on branch `temporal-v4` and merge once happy.
- [ ] (deferred) `tests/` directory now exists — could move other future
      tests there. Currently only v4 tests live there.

## Convention notes for future v* iterations

- Helper functions go in standalone `.m` files (`v4_*.m`) so tests can
  call them directly without hoisting locals out of a script.
- Synthetic ground-truth tests come BEFORE running on real data, per
  CLAUDE.md TDD policy. Python prototype is acceptable as an algorithmic
  sanity check; the MATLAB test is the authoritative pre-flight.
- Don't modify v3 (or any prior version). New design = new file. Old
  versions are reference points for reviewers who want to see what
  changed.
