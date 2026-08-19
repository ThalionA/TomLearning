# Appendix B — `scripts/` map (raw sweep, 2026-08-19)

Kind: RUN / ANALYZE / FIGS / PLOT / EXPLORE / PROTO / UTIL / BATCH / MATLAB.
refs = line-count mentions of the basename in HANDOFF (H), README (R), STATE (S),
PROJECT_LOG (PL), GOTCHAS (G), NOTES (N), OPPORTUNITIES (O), results/*.md (res),
scripts/*.sh (sh), other scripts (py), tests (T). PL-date = nearest dated heading above the
newest PROJECT_LOG mention. 🔴 = zero references by name. 🟠 = newest PL mention < 2026-07-01
or none. **Caution:** "zero refs by name" ≠ dead — `analyze_cc_label_anova.py` is 🔴 yet its
CSV feeds `figs_cc_label_track.py` (verified). Check outputs before archiving.

## TEMPORAL arm

| file | LOC | purpose | kind | src imports | script imports | refs | PL-date | git |
|---|---|---|---|---|---|---|---|---|
| `run_lag_curves.py` | 220 | Held-out per-dim lagged CC curves, segment-aware, refit per lag/fold | RUN | config, dataio, lagged, partial | — | 23 | 2026-08-17 | 2026-08-17 |
| `run_lag_subspaces.py` | 265 | Lag-0 vs lagged principal angle; FF(+τ)/FB(−τ) fits | RUN | config, dataio, lag_subspace, membership, subspace_window | — | 13 | 2026-08-03 | 2026-08-03 |
| `run_cc_label_track.py` | 268 | Fit once on all running bins, freeze axes, label FF/FB, project epochs | RUN | config, dataio, fixed_subspace, lagged, partial | — | 11 | 2026-08-07 | 2026-08-19 |
| `run_fixed_subspace.py` | 176 | Epoch-balanced frozen subspace, one CC lagged per epoch | RUN | config, dataio, fixed_subspace, partial | — | 8 | 2026-07-29 | 2026-08-03 |
| `run_lag_cosine.py` | 276 | cos of each canonical vector to itself across lags; rank-swap rate | RUN | config, core, dataio, fixed_subspace, lag_subspace, partial | — | 8 | 2026-08-07 | 2026-08-06 |
| `run_ifi_windows.py` | 135 | Held-out lag curve → IFI per window (early) | RUN | config, dataio, lagged, partial | — | 7 | 2026-08-15 | 2026-06-16 |
| `run_transition.py` | 254 | Uncued→cued subspace readout | RUN | config, dataio, membership, paired_stats, subspace, subspace_window | — | 13 | 2026-07-29 | 2026-07-28 |
| `run_temporal.py` | 232 | Original Arm-A/B committed driver (`--arm runstate|landmark`) | RUN | analysis, config, dataio, pipeline, sweep | — | 2 🟠 | none | 2026-06-01 |
| `analyze_cc_ifi_signs.py` | 587 | Item 1 — per-CC IFI over window sweep; overall IFI test | ANALYZE | cc_aggregate, perdim_ifi, paired_stats | — | 6 | 2026-08-15 | 2026-08-17 |
| `analyze_cc_crosscorr_epochs.py` | 446 | Ask 2 — cross-correlograms naive/expert, all/FF/FB | ANALYZE | cc_aggregate, paired_stats | — | 5 | 2026-08-15 | 2026-08-19 |
| `analyze_lag_subspaces.py` | 436 | Items 2/3/4 — angle vs split-half floor; FF−FB | ANALYZE | paired_stats | — | 9 | 2026-08-03 | 2026-08-03 |
| `analyze_fixed_subspace.py` | 264 | Items 5/6/7 — frozen CC1 curve metrics naive vs expert | ANALYZE | fixed_subspace, paired_stats, perdim_ifi | — | 6 | 2026-07-29 | 2026-08-03 |
| `analyze_cc_label_track.py` | 258 | Item 2 — FF/FB labels across epochs | ANALYZE | lagged, paired_stats, perdim_ifi | — | 4 🟠 | none | 2026-08-06 |
| `analyze_cc_label_anova.py` | 227 | Item 2 factorial RM-ANOVA | ANALYZE | none | — | 0 🔴 (but output consumed by `figs_cc_label_track`) | none | 2026-08-07 |
| `analyze_perdim_ifi.py` | 222 | Item 1 — per-dim IFI from lag_curves | ANALYZE | mixed_effects, paired_stats, perdim_ifi | — | 4 | 2026-07-29 | 2026-07-28 |
| `analyze_ifi_windows_epochs.py` | 189 | Item 4 — window vs IFI by epoch | ANALYZE | paired_stats, perdim_ifi | — | 3 | 2026-08-07 | 2026-08-07 |
| `analyze_ifi.py` | 165 | IFI battery (animals- and dims-as-n) | ANALYZE | config | — | 4 | 2026-08-07 | 2026-06-11 |
| `analyze_bin10_full.py` | 301 | → `results/bin10_tables.md` | ANALYZE | config, paired_stats | `import figs_report as F` | 5 🟠 | 2026-06-17 | 2026-06-18 |
| `analyze_ifi_windows.py` | 87 | Which window maximises IFI reproducibility | ANALYZE | config | — | 3 | 2026-08-07 | 2026-06-11 |
| `figs_cc_crosscorr_epochs.py` | 246 | Ask-2 figures | FIGS | none | figstyle, `from analyze_cc_crosscorr_epochs import …` | 2 | 2026-08-15 | 2026-08-19 |
| `figs_cc_label_track.py` | 229 | Item-2 figures (+ ANOVA heatmap) | FIGS | none | figstyle | 1 🟠 | none | 2026-08-07 |
| `figs_cc_ifi_signs.py` | 222 | Item-1 figures | FIGS | cc_aggregate | figstyle | 1 🟠 | none | 2026-08-17 |
| `figs_integration_windows.py` | 182 | Half-max width anatomy | FIGS | fixed_subspace | figstyle | 1 🟠 | none | 2026-08-03 |
| `figs_lag_cosine.py` | 165 | cos_same vs cos_best | FIGS | none | figstyle | 1 🟠 | none | 2026-08-06 |
| `figs_perdim_ifi.py` | 158 | IFI vs canonical rank | FIGS | perdim_ifi | figstyle | 0 🔴 by name (figure `HCV1_perdim_ifi_*` regenerated 08-17 per PROJECT_LOG) | none | 2026-07-28 |
| `figs_ifi_windows_epochs.py` | 150 | Item-4 figures | FIGS | none | figstyle | 1 🟠 | none | 2026-08-07 |
| `figs_lag_subspaces_epochs.py` | 145 | FF/FB CC1 across epochs | FIGS | none | figstyle | 1 🟠 | none | 2026-08-03 |
| `figs_lag_subspaces.py` | 143 | Angle vs floor; FF/FB panels | FIGS | none | figstyle | 1 🟠 | none | 2026-08-04 |
| `figs_fixed_subspace.py` | 122 | Frozen-subspace lag curves per epoch | FIGS | none | figstyle | 1 🟠 | none | 2026-08-04 |
| `figs_ifi_windows.py` | 89 | IFI vs window per pair (early) | FIGS | config | figstyle | 1 🟠 | 2026-06-11 | 2026-06-13 |
| `figs_lag_curves.py` | 88 | Held-out lag CC curves per pair | FIGS | none | figstyle | 4 | 2026-08-17 | 2026-06-18 |
| `merge_lag_curve_parts.py` | 43 | Merge `lag_curves_*_part*.csv` | UTIL | none | — | 3 | 2026-08-17 | 2026-08-17 |
| `run_lag_curves_uncapped_batch.sh` | 35 | Uncapped batch, 3 procs × 2 FS | BATCH | — | run_lag_curves, merge | 3 | 2026-08-17 | 2026-08-17 |

## SPATIAL arm (all 🟠 unless noted)

`run_trajectory.py` (243, PL 07-28) · `run_early_trials.py` (223) · `run_epochs.py` (167, PL 07-28) ·
`run_committed.py` (155) · `run_stage2.py` (152) · `run_trajectory_bins.py` (134) · `run_partial.py` (116) ·
`run_stage3.py` (79) · `summarise_sweep.py` (359) · `analyze_trajectory.py` (189) ·
`learning_changes_spatial.py` (192) · `analyze_early_trials.py` (187) · `analyze_epochs.py` (133) ·
`analyze_trajectory_dims.py` (131) · `compare_nulls.py` (120) · `analyze_learning_vs_time.py` (110) ·
`analyze_dims_as_n.py` (68) · `committed_ifi.py` (163) · `figs_report.py` (423, PL 07-29, 20 refs) ·
`figs_paired.py` (154) · `figs_units.py` (143) · `figs_stats_tables.py` (133, imports `figs_report`) ·
`figs_area_gini.py` (131) · `figs_early_trials.py` (125) · `figs_trajectory_dims.py` (112, imports
`analyze_trajectory_dims`) · `figs_rotation_cc1.py` (103) · `figs_trajectory_bins.py` (95) ·
`plot_sweep_explore.py` (403) · `plot_stage2.py` (319) · `plot_ifi_explore.py` (293) · `plot_stage3.py` (245) ·
**🔴** `plot_parcoords.py` (221) · **🔴** `plot_common_units.py` (188) · **🔴** `plot_subspace_similarity.py` (151) ·
**🔴** `plot_ifi_fs.py` (144) · **🔴** `plot_partial.py` (126) · `prototype_continuous_pcca.py` (132) ·
`run_review_batch.sh` (30).

## LANDMARK arm (all 🟠)

`summarise_landmark_sweep.py` (547) · `plot_landmark.py` (489, "adapted from plot_stage2") ·
`learning_changes.py` (326) · `build_landmark_index.py` (258) · `learning_changes_mixed.py` (139) ·
`build_prune_table.py` (102) · `regen_all_landmark_figs.sh` (50).

## KCCA arm (all 🟠)

`figs_kcca.py` (227) · `prototype_kcca.py` (180) · `analyze_kcca.py` (172) · `run_kcca.py` (156) ·
`run_kcca_transition.py` (154).

## SHARED

`figstyle.py` (188, imported by 23 scripts, 29 refs) · **🔴** `make_figure_deck.py` (179) ·
`make_smoothed_deck.py` (159) · `export_cca_labels.m` (120, MATLAB; produces `cca_labels.json`).

## Duplication candidates (file:line in the sweep; summarised)

- **D1 `PAIRS` ×51**: `config.PAIRS` used by 2 scripts; 29 string-form, 15 tuple-form copies; `prototype_kcca` 5-pair subset.
- **D2 `K` ×14** with values 30/20/15. No shared constant.
- **D3 `MAX_SAMPLES=12000`/cap logic ×4** (`run_lag_cosine:58,88`, `run_lag_subspaces:48,80`, `run_lag_curves:89`, `run_fixed_subspace:58`) + inline in `run_trajectory:170`.
- **D4 divergent constants**: `MAX_LAG` 25/10/8/4; `TAU_BINS=5` ×2; `LABEL_W=5` ×2 with "must match" comment; `BIN_MS=10` ×3.
- **D5 results/ATT dirs**: `RES = parents[1]/"results"` ×20 (config.RESULTS_DIR used by 3); `ATT` absolute ×7 vs `Path.home()` ×17; neither in config.
- **D6 `sys.path.insert(0, …/src)`** in ~80 scripts.
- **D7 FS suffix ×36 in two idioms**; polarity inversion in `figs_area_gini.py:64`.
- **D8 inline stats helpers**: `_wilcoxon_vs0` ×7, `_paired` ×5, `_unpaired` ×2, local `_fdr_bh` (`learning_changes.py:123`), star formatters ×9, `_sem` ×8, 23 scripts calling scipy directly.
- **D9 bar/panel plotters**: `figs_kcca._bar_points_sem` ≈ `figs_report._bar_points_sem`.
- **D10 plot_* vs figs_* overlap**: four IFI-vs-window renderings (`plot_ifi_explore`, `figs_ifi_windows`, `figs_ifi_windows_epochs`, `committed_ifi`, `plot_ifi_fs`); three rotation-vs-floor; `plot_landmark` ≈ `plot_stage2` (same helper set); `make_figure_deck` ≈ `make_smoothed_deck`.
- **D11 argparse**: 29 parsers; `--include-fs` ×15, `--bin-ms` ×13, `--smooth-ms` ×11, `--max-lag` ×9, `--tag` ×9, `--out` ×8; no parent parser.
- **D12 per-cell prep loop** (`if len(idx) >= cfg.min_units:` + stream load + area prep) copied in 14 drivers.
- **D13 per-script CSV/pkl loaders** ×13 + ad-hoc `pd.read_csv(RES / f"…{suffix}.csv")` in every analyze/figs script.
