# Appendix D — temporal-arm structural deep read (raw, 2026-08-19)

Read-only audit of the method layer + drivers by an Opus subagent; items marked
**[verified]** were re-read in source by the main session. Paths relative to `cca/`.
One correction to the raw output is applied in §A6 (see note there).

## A. Drift risk

### A1. Lag pairing — 9 implementations, identical arithmetic, different guards [verified 5]

| impl | shape | short-segment guard |
|---|---|---|
| `lagged.lag_slice:23` | 3-D (trials,bins,k) | raises if \|lag\|≥n_bins |
| `lagged._segment_lagged_pairs:145` | flat + `groups` | skip if `n <= |lag|+2` (`:159`) |
| `lagged_temporal._segment_pairs:48` | flat + `Segment` list | skip if `|L| >= n` (`:66`) |
| `lagged_landmark._window_pairs:52` | (crossings, win, k) | `|L| >= win_len` |
| `segments.segment_paired_indices:124` | index-only | `|L| >= seg_len` — **no production caller** |
| `lag_subspace.segment_lag:48` | flat + `groups` + Z twice | `n <= |lag|+2` (`:61`) |
| `subspace_window._lag_curve_perdim:138` | flat, **crosses trial boundaries** | `< 10 rows` (`:151`) |
| `fixed_subspace.variate_lag_curve:123` | 1-D variates + `groups` | `n <= |lag|+2` (`:140`) |
| `fixed_subspace.trial_lag_moments:356` | same, as moments | `n <= |lag|+2` (`:381`) |

Sign convention consistent everywhere (`x[:n-L] / y[L:]` for L≥0). `fixed_subspace.py:394-400`
asserts `curve_from_moments` ≡ `variate_lag_curve` — no test checks it.

### A2. CCA fitting
`core.cca_fit:156` (covariance route, `_COV_EIG_FLOOR=1e-10` at `:142`) is the only CCA on the
temporal path; `core._cca_fit_svd:228` is the pinned reference (`tests/test_core.py:295-326`);
`fixed_subspace.observed_canonical_r:281` / `frozen_perm_null:312-313` compute canonical r as singular
values of `Qx.T@Qy` — a different estimator of the same quantity, and both reach
`cc_label_track_*.csv` (`run_cc_label_track.py:166` `r_fit`; `:255` `cc_heldout`). `kernel_cca.kcca_fit:73` off-path.

### A3. PCA → scores — 7 implementations
`core.pca_fit/pca_transform:102,123` (NaN-safe) · `subspace_window._pca_fit/_scores:60,71` ·
`lag_subspace._scores:78` · `early_trials._pca_fit:40` · inline SVD `fixed_subspace.fit_fixed:91-95` ·
driver `_pca_scores` at `run_lag_curves.py:83` ≡ `run_lag_cosine.py:82` (md5 identical [verified]),
`run_ifi_windows.py:51` (variant).

### A4. Partial-out — 1 implementation, 4 leakage levels at call sites [verified run_lag_curves]
`partial.partial_out_cv:17` (real), `partial_out:34`, `partial_out_tensor:43`, `partial_cca_cv:57`
(legacy in-sample Z); `kcca_window._residualise:48` a copy. Call sites: `run_lag_curves.py:168` and
`run_lag_cosine.py:158` in-sample over all rows before a held-out curve; `lag_subspace.lag_sweep:171`
per-fold; `run_fixed_subspace.py:140` fit-trials; `run_cc_label_track.py:165` all rows (by design).

### A5. Fold construction — 4 rules
1. `core.trial_folds:383` interleaved `perm[i::n_folds]` (spatial, `surrogate.build_null`, `partial_cca_cv`).
2. `np.array_split(rng.permutation(uniq), n_folds)`: `lagged.py:188`, `lag_subspace.py:103,163`, `subspace_window.py:90`.
3. Leave-one-trial-out: `lagged_temporal._held_out_cc_lotot:81-116` (archived pipeline — **not** 5-fold).
4. Split-half by trials: `subspace.split_half_angles:44`, `subspace_window._split_half_angles:179`,
   `lag_subspace.split_half_floor:240`, inline `run_lag_cosine.py:172-175`.
Within `lag_subspace`: `lagged_fit` builds the fold mask once (`:107`), `lag_sweep` rebuilds it on
lag-surviving rows (`:192-193`).

### A6. Null distributions — 8 implementations, 3 live [verified 2]

| impl | shuffled | granularity | min shift | statistic |
|---|---|---|---|---|
| `surrogate.build_null:70` + `permute_trials:29`/`circshift_bins:34` | Y tensor | per-trial | `cfg.circshift_min_bins=15` | 5-fold held-out CC refit (spatial) |
| `lagged_temporal.circshift_segments:197` | Y scores | per segment | `temporal_circshift_min_bins=5` | LOTO lag scan (archived) |
| `lagged_landmark.circshift_windows:153` | Y windows | per window | 3 | LOTO (archived) |
| **`lagged.perdim_significance:238`** (`np.roll` `:287`) | Y scores | global | 1 | held-out per-dim @lag 0 (`run_lag_curves`) |
| **`subspace_window._significance:116`** (`:132`) | Y scores | global | 1 | in-sample dominant dim, no p (`run_lag_subspaces` via `window_subspace`) |
| **`fixed_subspace.frozen_perm_null:287`** (`:323`) | orthonormal `Qy` | global | 1 | in-sample rank-for-rank r (`run_cc_label_track`, `run_lag_cosine`) |
| `kcca_window._significance:85` | Y scores | global | 1 | held-out KCCA |
| `fixed_subspace.side_peak_null:232` | synthetic | — | — | different question |

The three live nulls ignore `cfg.temporal_circshift_min_bins` (`config.py:213` — dead config) and use a
single global roll. **Main-session note:** the raw output called this "the exact v4-era silent bug";
that overstates it. `circshift_segments`' docstring describes importing *gap structure* into a
per-segment lag-scan null; a single global roll over ~3e5 concatenated running bins at lag 0 is a valid
surrogate (P(shift < 50 bins) ≈ 2e-4 uncapped, ≈ 4e-3 at the 12 k cap). Treat as config/doc drift +
three constructions, not a result-changing bug.

### A7. IFI
`lagged.information_flow_index:42` canonical (`total<=0 → 0.0`); `lagged.ifi_by_window:54`;
`perdim_ifi.curve_ifi:28` returns **NaN** for empty sides (`:45-46`); `perdim_ifi.ifi_windows_by_dim:159-162`
re-computes the clipped means inline for `degenerate`; `run_ifi_windows.py:117-118` third inline copy;
`analyze_cc_label_track.py:73` → `lagged.information_flow_index`, `analyze_fixed_subspace.py:56` →
`perdim_ifi.curve_ifi` on the same class of curve; `run_cc_label_track.py:203,214,232` on the bin axis,
`analyze_cc_label_track.py:64,73` re-derives the label window on the ms axis.

### A8. FDR — 4 identical BH implementations [verified 3]
`paired_stats.fdr_bh:71`, `lagged.fdr_bh:323` (docstring: deliberate), `fixed_subspace._fdr_bh:340`,
`scripts/learning_changes.py:123`. p-floor handling explicit (`lagged.py:263-270`, `fixed_subspace.py:308-310`)
via `fdr_dims`; `FDR_DIMS=10` hard-coded at `run_lag_curves.py:54`, `run_lag_cosine.py:61`,
`run_cc_label_track.py:67`. `subspace_window._significance` applies no correction and emits no p.

### A9. Per-animal aggregation
`cc_aggregate.per_animal_mean:26` + `one_sample_by_pair:72` (sig gate `:48-51`, degenerate drop `:52-53`,
duplicate-(pair,animal) raise `:87-90`) used by `analyze_cc_ifi_signs.py:180,196`,
`analyze_cc_crosscorr_epochs.py:97,163`. `analyze_cc_label_track._per_animal:79` pivot_table with none of
the guards; `analyze_fixed_subspace.epoch_contrasts:74-87` by hand; `analyze_lag_subspaces._paired:46` thin.
Dims-as-n power checks (deliberate, documented): `analyze_cc_ifi_signs.overall_direction_ccs:205`,
`analyze_cc_crosscorr_epochs.epoch_contrast_ccs:187`.

### A10. Running-bin selection — 2 definitions
`dataio.temporal_segments:624` (velocity mask + `close_short_gaps` + `min_run_bins` + `Segment`s; sole
caller `pipeline.prepare_pair_temporal:360`, archived) vs inline `run = ~isnan(trial_idx) & (vel >= thr)`
in `run_lag_curves.py:144`, `run_lag_subspaces.py:123`, `run_lag_cosine.py:139`, `run_fixed_subspace.py:96`,
`run_cc_label_track.py:130`, `run_ifi_windows.py:85`, `run_epochs.py:72`, `run_trajectory.py:138`,
`run_early_trials.py:124`, `run_kcca.py:100`, `run_trajectory_bins.py:63`, both prototypes. **Every current
result uses the driver version** (no gap-closing, no min run). Velocity threshold and z-score reference
(`dataio._zscore_engaged:599`) single-sourced.

### A11. Sample cap — 4 implementations, 2 selections [verified]

| impl | totals | CLI | surfaced |
|---|---|---|---|
| `run_lag_curves.py:89` | 0/0 (uncapped since 2026-08-15) | yes `:72-77` | printed `:121-126`, not a column |
| `run_lag_subspaces.py:80` | 12000/600 (`:48-49`) | no | no |
| `run_lag_cosine.py:88` | 12000/600 (`:58-59`) | no | no |
| `run_fixed_subspace.py:58` `_cap_bins` | per-trial 600, no total | no | no |
| `run_cc_label_track.py:138` | none (`:68-75`) | — | — |

The three 12 k versions select the same bins (ascending trials, leading 600 of each → first ~20 trials);
`_cap_bins` selects a different sample (all trials, first 600 bins each) feeding the same naive-vs-expert
contrast as `run_cc_label_track` (uncapped). Unrelated: `run_trajectory.py:40` 6000,
`run_trajectory_bins.py:25,86-87` 8000 even stride, `run_transition.py:36` 150 000 ms.

### A12. FS inclusion — no drift
`dataio.select_units:425-435` only; drivers flip `exclude_fast_spiking=not args.include_fs`. Caveat:
FS-included changes 4 of 6 areas (DG/SUB carry no flags); only the `_fsincl` suffix records the condition.

## B. API inconsistency
- PC count: `k` (core/subspace_window/lag_subspace/fixed_subspace) vs `cfg.k_cap/k_fixed/k_mode` vs driver `K=30`; `n_components` only in `kernel_cca`.
- Lag units: `max_lag` bins vs `max_lag_bins` vs explicit `lags` vs `cfg.lag_ms` (ms) vs `cfg.max_lag_bins` (spatial); analyze scripts on `lag_ms`, run scripts on `lag_bins`; `perdim_ifi.ifi_windows_by_dim:123-127` warns lag must be bins.
- Bin width: `cfg.temporal_bin_ms` vs `bin_ms` param vs `--bin-ms` vs spatial `bin_size_cm`.
- Trial grouping: `groups` vs `trials` vs `segments_` vs `trial_ids`.
- Lag-curve return: `LagResult` / `TemporalLagResult` / `LandmarkLagResult` / `(lags, cc)` / `dict[int, LagFit]` / bare `ndarray`.
- Significance return: `PerDimSignificance(mask,p,threshold,null_mode)` / `FrozenSignificance(mask,p,threshold,r_obs)` / `NullResult(...)` / bare bool array.
- Shuffles: `cfg.n_shuffles=200` vs `SURROGATE_SHUFFLES=100` vs `window_subspace(n_shuffles=30)` vs `--n-shuffles 200`; `run_lag_subspaces.py:206` passes 100, `run_lag_curves.py:67` 200.
- Seeds: `cfg.cv_seed/surrogate_seed` vs `seed=0` + ad-hoc `seed+7` (`lagged.py:282`, `subspace_window.py:127`, `fixed_subspace.py:319`) and `seed+11` (`subspace_window.py:191`, `lag_subspace.py:261`).
- Angles: `subspace.principal_angles:30` (rad, all) vs `lag_subspace.subspace_angle:220` (deg, max) vs `subspace_window._max_angle:159` (duplicate).
- Pairs: `config.PAIRS` tuples (one order) vs driver tuples vs analyze strings.

## C. Dataset coupling inside `src/`
`config.py` — right place: `DATA_DIR`/pattern/learning file `:31-33`, `COMPANION_FILE :40`, `AREAS :51`,
`PAIRS :52-61`, `FS_AREAS :68`, `N_BINS/CORRIDOR_CM :80-81`, `EPOCH_NAMES :93`, `EPOCH_COLOURS :97-101`,
`spatial_field :144`, `temporal_bin_ms=50 :186`, landmark window `:202-208`.
`dataio.py` — .mat keys `:80,114,129,147,267-273,293,572-577`, `TF<id>` regex `:162`, area masks keyed by `AREAS :200`.
Learning-point / epoch design in the method layer: `dataio.classify_cohort:369-384`, `epoch_windows:400-419`;
`pipeline.py:381,398`, `analysis.py:209` iterate `EPOCH_NAMES`.
Dataset facts in method docstrings/defaults: `lag_subspace.py:249-251` (CA1 vs SUB unit counts),
`fixed_subspace.py:191-196` + `band_occupancy:255` (theta 5–12 Hz), `core.py:25-30` (NaN rates).
`bin10` not in `src/`.

## D. Minimal port set

| stage | functions | module |
|---|---|---|
| load | `load_animals`, `load_animal`, `_read_behaviour_file`, `_read_companion`, `_load_temporal_streams`, `rebin_spikes`, `bin_velocity`, `bin_trial_index`, `classify_cohort`, `epoch_windows`, `select_units`, `area_activity_50ms`, `_zscore_engaged` | `dataio` |
| preprocess | inline running mask + `_capped_index` (**drivers only**), `partial.partial_out(_cv)`, `_pca_scores` (**drivers only**) | `partial` + drivers |
| fit | `cca_fit`, `cca_score` | `core` |
| lag curves | refit: `_segment_lagged_pairs`, `heldout_lag_curve_flat_perdim`, `heldout_lag_curve_flat` (`lagged`) · subspace: `segment_lag`, `lag_sweep`, `subspace_angle`, `split_half_floor` (`lag_subspace`) · frozen: `balanced_trials`, `fit_fixed`, `project`, `trial_lag_moments`, `curve_from_moments`, `curve_half_width`, `side_peak` (`fixed_subspace`) | 3 |
| nulls | `perdim_significance` + `fdr_bh` (`lagged`) / `frozen_perm_null` (`fixed_subspace`) | 2 |
| IFI | `information_flow_index`, `ifi_by_window` (`lagged`); `curve_ifi`, `curve_peak_lag`, `ifi_windows_by_dim`, `sign_agreement`, `rank_split` (`perdim_ifi`) | 2 |
| aggregate | `per_animal_mean`, `one_sample_by_pair` (`cc_aggregate`); `paired_t`, `wilcoxon_signed`, `fdr_bh` (`paired_stats`) | 2 |

Not on the path: `lagged_temporal`, `lagged_landmark`, `segments`, `pipeline.prepare_pair_temporal/landmark`,
`analysis.analyse_pair_temporal/landmark`, `surrogate`, `core.{residualise,zscore_units,pca_fit,pca_transform,cca_cv,trial_folds,choose_k,numerical_rank,k_for_variance,cca_in_sample}`,
`lagged.{lag_slice,lag_curve,LagResult}`, `subspace.{canonical_weights,split_half_angles}`,
`subspace_window.window_subspace` (only `n_sig` consumed at `run_lag_subspaces.py:204-207`),
`partial.{partial_out_tensor,partial_cca_cv}`.

## E. Robustness smells
- **Silent defaults not in outputs**: caps (A11); `_COV_EIG_FLOOR=1e-10` (`core.py:142`) truncates d silently; `CC_CLIP=0.999` (`subspace_window.py:29`); `NOT_ESTIMABLE_DEG=70` (`analyze_lag_subspaces.py:54`); `LEAD_DIMS=5` (`analyze_cc_ifi_signs.py:48`); `far_ms=200` (`analyze_cc_crosscorr_epochs.py:124`); `TAU_BINS=5`/`ANGLE_DIMS=3` (`run_lag_subspaces.py:50-51`); `FDR_DIMS`, `K`, `N_FOLDS`; `min_units`, `velocity_thresh` never in a row.
- **Mislabelled output [verified]**: `run_cc_label_track.py:255` writes `res.r_obs` (in-sample) as `cc_heldout`; consumed `analyze_cc_label_track.py:71`.
- **Long functions**: `run_cc_label_track.main:95` (170 lines), `run_lag_cosine.main:111` (162), `run_lag_curves.main:110` (107), `run_fixed_subspace.main:66` (107), `analyze_cc_ifi_signs.main:479` (105), `pipeline.prepare_pair_landmark:474` (98), `subspace_window.window_subspace:203` (87), `pipeline.prepare_pair_temporal:335` (86), `lagged.perdim_significance:238` (83), `lag_subspace.lag_sweep:137` (81), `run_lag_subspaces._do_pairs:186` (76).
- **Hard-coded verdicts**: `analyze_fixed_subspace._md_ringing:197-202`; `analyze_lag_subspaces.stability:57-74` docstring; `analyze_cc_ifi_signs._mixing_verdict:413` now data-driven (fixed 08-17).
- **Broad `except`** (49 total): drivers `except Exception: continue` around stream load (`run_lag_curves.py:140-143`, `run_lag_subspaces.py:119-122`, `run_lag_cosine.py:135-138`, `run_fixed_subspace.py:92-95`, `run_cc_label_track.py:126-129`); `run_lag_subspaces.py:208-209` failed `window_subspace` → `n_sig=0`; `lagged.py:296` failed surrogate → empty row, denominator shifts (`:307-308`); `fixed_subspace.py:98`; `kcca_window.py:73,101,160`; bare `assert` for alignment invariants `lagged_temporal.py:167`, `lagged_landmark.py:120`.
- **Numerical**: IFI 0.0 ambiguity (`lagged.py:49-50`) — disambiguated only in `perdim_ifi.ifi_windows_by_dim:163`; `subspace_window._significance:131` `integers(1,n)` raises at n≤1; `core.cca_score:296` all-NaN for <3 finite → `nanmean` over folds hides "no fold fittable"; `analyze_cc_crosscorr_epochs.curve_metrics:145-147` `peak=max(r)` no positivity guard; `paired_stats.paired_t:54` allows n=2 vs `cc_aggregate` `min_n=3` vs `>=2` in `analyze_cc_label_track.py:111,130`, `analyze_fixed_subspace.py:86,148`, `analyze_lag_subspaces.py:48`. Guarded well: `curve_from_moments:404-409`, `membership.gini:156`, `partial.partial_out_cv:27-30`.
- **Tests**: `test_core.py:278-326` pins `cca_fit` to `_cca_fit_svd` (deliberate); `test_early_trials.py:54-55` pins `_gini_of_corr`. Gaps: no test for any `_capped_index`; none that the four `fdr_bh` agree; none that `curve_from_moments` ≡ `variate_lag_curve`; none that the lag pairers agree on a shared fixture.
