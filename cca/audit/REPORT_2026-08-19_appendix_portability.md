# Appendix C — dataset / machine coupling (raw sweep, 2026-08-19)

Scope `src/tom_cca/`, `scripts/`, `tests/`; excludes `__pycache__`, `results/`, `figures/`.

## 1. Absolute / machine paths

Literal `/Users/theoamvr/...`: `scripts/figs_report.py:35`, `figs_units.py:37`, `figs_early_trials.py:29`,
`figs_ifi_windows.py:25`, `figs_kcca.py:34`, `figs_rotation_cc1.py:30`, `analyze_epochs.py:32`
(all `ATT = Path("/Users/theoamvr/Documents/ResearchVault/attachments")`); `run_review_batch.sh:6`.

`Path.home()/"Documents"/"ResearchVault"/"attachments"`: `figs_lag_subspaces_epochs.py:37`,
`figs_cc_ifi_signs.py:40`, `figs_paired.py:37`, `figs_lag_curves.py:29`, `figs_cc_label_track.py:43`,
`figs_perdim_ifi.py:31`, `figs_lag_subspaces.py:33`, `figs_lag_cosine.py:37`, `figs_integration_windows.py:41`,
`figs_cc_crosscorr_epochs.py:51`, `figs_ifi_windows_epochs.py:37`, `figs_trajectory_dims.py:31`,
`figs_stats_tables.py:30`, `figs_area_gini.py:33`, `figs_trajectory_bins.py:29`, `figs_fixed_subspace.py:34`;
`make_figure_deck.py:25,27`, `make_smoothed_deck.py:22,23` (→ `.../Projects/Hippocampus-V1`).
26 scripts contain the `attachments` literal; only `figstyle.py:31` derives a repo-relative path
(`REPO_FIGURES = parents[1]/"figures"`, mirror at `:164-166`).

## 2. Dataset names / .mat keys

Centralised: `config.py:31` `DATA_DIR = _REPO_ROOT/"HC_V1_data"` (`_REPO_ROOT = parents[3]`),
`:32` `DATA_FILE_PATTERN = "TF*_export.mat"`, `:33` `LEARNING_FILE = "animal_behaviour.mat"`
(**never referenced**), `:40` `COMPANION_FILE = "cca_labels.json"`.

HDF5 keys, all in `dataio.py`: `:80` `regions_label`; `:114` `units.idx`; `:130` `idx_fs`;
`:147-149` `analysis_spatial.firing.cued.freq/freq_z`; `:185` requires `units` + `analysis_spatial`;
`:267-273` `animal_id`, `period_experienced`; `:293-295` `zscored_lick_errors`; `:572-577`
`binned_spikes`, `data_behaviour`, `velocity_binned_gf`, `trial_binned_cued`, `visual_landmark_binned`;
`:238-243` companion-JSON schema.

Keys leaking outside dataio: `run_transition.py:49`, `run_kcca_transition.py:47`
(`f["data_behaviour"]["world_binned"]`).

`"animal_behaviour.mat"` re-spelled (calling private `dataio._read_behaviour_file`) in:
`run_transition.py:89`, `run_ifi_windows.py:66`, `run_temporal.py:210`, `run_lag_subspaces.py:103`,
`run_trajectory.py:108`, `prototype_kcca.py:104`, `run_trajectory_bins.py:50`, `run_kcca.py:78`,
`run_fixed_subspace.py:74`, `run_early_trials.py:99`, `run_lag_curves.py:118`, `run_kcca_transition.py:68`,
`run_epochs.py:56`, `prototype_continuous_pcca.py:75`, `run_lag_cosine.py:119`, `run_cc_label_track.py:103`.

Test coupling: `tests/test_dataio_temporal.py:96` `config.DATA_DIR/"TF028_export.mat"` (skips if absent).

## 3. Area / pair hard-coding

`config.py:51` `AREAS`, `:52-61` `PAIRS` (8), `:68` `FS_AREAS`. Imported by only
`plot_sweep_explore.py:49`, `plot_ifi_explore.py:46`. Tuple-form redefinitions (different order): 13
drivers (`run_ifi_windows.py:33`, `run_transition.py:42`, `run_lag_subspaces.py:54`, `run_trajectory.py:43`,
`run_trajectory_bins.py:27`, `run_kcca.py:39`, `run_fixed_subspace.py:42`, `run_lag_curves.py:52`,
`run_kcca_transition.py:34`, `run_epochs.py:29`, `run_early_trials.py:47`, `run_cc_label_track.py:77`,
`run_lag_cosine.py:62`); 5-pair subsets `prototype_kcca.py:39`, `prototype_continuous_pcca.py:34`.
String-form `"CA1-RSC"` in ~35 analyze/figs scripts. `figs_area_gini.py:35` redefines `AREAS` in another
order. `analyze_bin10_full.py:30` takes `PAIRS` from `figs_report`.
**`src/` is clean: zero area-name literals outside `config.py`; only `dataio.select_units:434` consults `FS_AREAS`.**

## 4. Animal-specific

`config.py:119` `manual_nonlearners = frozenset()` (only per-animal hook; consumed `dataio.py:381`).
No literal animal-id lists anywhere. Cohort = any animal with a recorded LP (`dataio.classify_cohort:369-384`).
Cohort-size assumptions: `analyze_ifi_windows_epochs.py:53` and `analyze_cc_label_track.py:49`
`EXPECTED_ANIMALS = 12`; `mixed_effects.py:27` `MIN_ANIMALS = 3`; `figs_perdim_ifi.py:36` `MIN_FRAC_ANIMALS = 0.5`.
`learners` filter re-implemented in `run_temporal.py:212`, `prototype_kcca.py:106`,
`prototype_continuous_pcca.py:77`, `committed_ifi.py:90`, `compare_nulls.py:51`.
Animal ids 1-indexed (MATLAB), inferred from filename digits (`dataio._infer_animal_id:161`).

## 5. Magic numbers

- **12000 cap**: `run_lag_subspaces.py:48`, `run_lag_cosine.py:58`; siblings `run_trajectory.py:40`=6000,
  `run_trajectory_bins.py:25`=8000, `run_transition.py:36`=150 000 ms, `run_lag_curves.py:50`=0.
  `MAX_BINS_PER_TRIAL=600` in `run_lag_subspaces.py:49`, `run_lag_cosine.py:59`, `run_fixed_subspace.py:40`.
- **k=30**: src defaults `kcca_window.py:107`, `early_trials.py:46,123`, `trajectory.py:60`,
  `lag_subspace.py:87,137,240`, `fixed_subspace.py:70`, `subspace_window.py:203`; scripts
  `run_ifi_windows.py:31`, `run_trajectory.py:35`, `run_trajectory_bins.py:25`, `run_lag_subspaces.py:46`,
  `prototype_kcca.py:33`, `run_lag_curves.py:43`, `run_fixed_subspace.py:37`, `prototype_continuous_pcca.py:32`,
  `run_early_trials.py:37`, `run_cc_label_track.py:62`, `run_lag_cosine.py:57`. Divergent: `run_transition.py:31` K=15,
  `run_kcca.py:31` 20, `run_kcca_transition.py:25` 15. Config has `k_cap=30` (`:156`), `k_fixed=10` (`:153`).
- **σ=2.5 ms**: `default=2.5` only `run_lag_subspaces.py:69`, `run_fixed_subspace.py:52`; `default=0.0`
  (help text says 2.5) `run_transition.py:80`, `run_ifi_windows.py:47`, `run_trajectory.py:55`,
  `run_lag_curves.py:66`, `run_early_trials.py:65-66`. Config `gaussian_sd_ms=0.0` (`:190`).
  Collision: 2.5 also = cm per spatial bin (`config.py:9,78,112,162`, `core.py:26`, `sweep.py:20,48`).
- **bin_ms**: config `temporal_bin_ms=50` "Locked at 50 ms" (`:186`); drivers default 10
  (`run_lag_subspaces:67`, `run_trajectory_bins:35`, `run_fixed_subspace:50`, `run_lag_curves:61`,
  `run_cc_label_track:86`, `run_lag_cosine:73`), 25 (`run_ifi_windows:41`, `run_transition:75`,
  `run_trajectory:49`, `run_kcca:47`, `run_epochs:38`), 50 (`run_temporal:50`); only
  `run_early_trials.py:64` reads `config.DEFAULT.temporal_bin_ms`.
- **min_units=5**: centralised (`config.py:127`), respected everywhere.
- **velocity 2.0 cm/s**: centralised (`config.py:199`), respected everywhere.
- **n_shuffles**: `config.n_shuffles=200` (`:170`, spatial), `config.SURROGATE_SHUFFLES=100` (`:75`) used
  by `run_transition:33`, `run_trajectory:38`, `run_trajectory_bins:25`, `run_epochs:26`,
  `run_early_trials:40`, `run_lag_subspaces:206`; local 30/10/20 in `run_kcca:35`, `run_kcca_transition:29`,
  `prototype_kcca:38`; `subspace_window.py:204` default 30; `fixed_subspace.py:287` default 200.
- **n_folds=5**: `config.py:159`; redeclared in 10 drivers; `run_transition:32`, `run_kcca_transition:30` = 4.
- **alpha=0.05**: never in config; src defaults in `lagged:239,323`, `subspace_stats:111,129,147`,
  `prune_table:219`, `subspace_window:204`, `surrogate:75`, `fixed_subspace:287,340`, `paired_stats:71`;
  script `ALPHA` in 10 scripts.
- **MAX_LAG**: 8/10/8/4/10/25 across `run_transition:34`, `run_trajectory:37`, `run_kcca:34`,
  `run_kcca_transition:28`, `run_early_trials:39`, `run_fixed_subspace:38`.

## 6. Output-location hard-coding

`config.py:41-43` `CCA_DIR`, `RESULTS_DIR`, `FIGURES_DIR` — used by ~25 scripts; bypassed by all 26
figs/deck scripts (vault primary, `figstyle.save(mirror=True)` copies back).
`HCV1_` prefix in 26 files (`figs_cc_crosscorr_epochs` ×10, `figs_report` ×9, `figs_cc_label_track` ×9, …).
`bin10` baked: `figs_lag_subspaces_epochs.py:90-134`, `figs_cc_label_track.py:100-217`,
`figs_cc_ifi_signs.py:154-213`, `figs_perdim_ifi.py:105-149`, `analyze_bin10_full.py:14,32`;
`analyze_trajectory_dims.py:69` bakes `w15`; `plot_ifi_fs.py:131` `win10`.

## 7. `config.py` — covered vs bypassed

Covered and respected: `min_units`, `velocity_thresh_cm_s`, `FS_AREAS`, `RESULTS_DIR` (run/analyze),
ms→bins helpers, `exclude_fast_spiking`, `zscore_units`.
Covered but bypassed: `PAIRS` (~48 scripts), `AREAS` (1), `LEARNING_FILE` (never used), `n_folds` (12),
`n_shuffles`/`SURROGATE_SHUFFLES` (competing), `k_cap` (20 sites), `temporal_bin_ms` (17 drivers override
with 3 defaults), `FIGURES_DIR` (26 figs scripts), `gaussian_sd_ms`.
Not covered: `alpha`, sample caps, temporal `MAX_LAG`, `HCV1_` prefix, attachment destination,
`EXPECTED_ANIMALS`.

## 8. `dataio.py` — the data contract

Docstring (`:1-21`): `TF*_export.mat` MATLAB v7.3/HDF5; spatial FR at `analysis_spatial.firing.cued.freq`,
on-disk `(n_units, n_trials, n_bins)` → loader `(n_trials, n_bins, n_units)`; animal ids 1-indexed from
filename digits; one animal per file; FS flagged only in V1/RSC/CA1/CA3.

`Animal` (`:39-64`): `animal_id`, `spatial_fr`, `area_masks`, `fs_mask`, `n_trials`, `filename`,
`streams_path` (None for fixtures ⇒ temporal arms unusable).

| function | line | contract |
|---|---|---|
| `load_animal(path, animal_id=None, spatial_field="freq", region_labels=None)` | 169 | needs `units` + `analysis_spatial`; `region_labels` overrides undecodable MATLAB string array |
| `load_animals(data_dir=None, spatial_field="freq")` | 304 | globs `DATA_FILE_PATTERN`; reads `cca_labels.json` |
| `classify_cohort(animals, cfg, behaviour_lookup=None)` | 369 | learners only = has LP and not in `manual_nonlearners` |
| `n_usable_trials(animal)` | 390 | full `n_trials` (no disengagement index in Tom's export) |
| `epoch_windows(lp, n_usable, cfg)` | 400 | naive `0..e-1`, intermediate `lp-e..lp-1`, expert `lp..lp+e-1`; None if `lp<2e` or `lp+e>n_usable` |
| `select_units(animal, area, cfg)` | 425 | area mask minus FS (FS only if `area in FS_AREAS`) |
| `area_tensor(animal, area, cfg)` | 438 | `(activity, unit_indices)` |
| `rebin_spikes(spikes_1ms, bin_ms, chunk_out_bins=2000, max_1ms=None, gaussian_sd_ms=0.0)` | 459 | `(n_1ms, n_units) → (n_bins, n_units)` float32; Gaussian before binning |
| `bin_velocity(vel_1ms, bin_ms)` | 509 | mean per bin; cm/s assumed |
| `bin_trial_index(trial_1ms, bin_ms)` | 518 | all-or-nothing bin→trial; NaN at boundaries |
| `temporal_segments(animal, cfg)` | 624 | Arm-A segments (archived pipeline only) |
| `area_activity_50ms(animal, area, cfg)` | 636 | `(spikes[:, kept], kept_idx)` z-scored over running & cued bins |
| `area_landmark_windows(animal, area, cfg)` | 653 | Arm-B windows |
| `_read_behaviour_file(path)` (private, called by 16 scripts) | 249 | `{('period_experienced', id): LP}` |
| `_load_temporal_streams` (private) | 553-597 | 1-slot cache keyed `(animal_id, bin_ms, sd_ms)`; reads `binned_spikes` (~1.6 GB/animal) + behaviour streams |
