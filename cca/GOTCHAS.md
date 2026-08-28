# GOTCHAS

One-line entries for non-obvious bugs, so they are not reintroduced.

- **Trial 1 is not behaviourally comparable to trial 2: it is LONGER, SLOWER and more
  FRAGMENTED.** Measured 2026-08-20 across 16 animals, both FS: trial 1 has **+650 running
  bins (≈ +6.5 s, p = 5.5×10⁻⁶, 12/13 animals)**, runs **−3.3 cm/s slower**, and has **+6.8
  more running gaps** than trial 2; speed and gaps survive BH in **8/8** pairs, duration in
  7/8 (CA3-DG n = 5 fails). Any per-trial contrast in this project must bin-match before it
  can be read as communication. ⚠ And bin-matching fixes ONLY duration: taking the leading
  N bins equalises `n_bins` exactly (Δ = 0.000, 71/71 animal-pairs) but leaves speed
  **undiminished or worse** (−3.30 raw → −3.51 matched) and gaps largely intact, because it
  compares roughly the first 80 % of trial 1 against the whole of trial 2 — and the phase it
  over-weights is the trial-onset phase already documented (HANDOFF §6) as carrying a doubled
  speed confound. `analyze_trial12.py` therefore tests the three covariates on EVERY arm; it
  originally tested them only on the raw arm, with a comment asserting matching removed them,
  which hid this entirely.

- **A frozen CCA projection can collapse on a single trial and read |r| ≈ 1 at EVERY rank.**
  Animal 70, ordinals 1 and 9, pairs CA1-DG / CA1-SUB / CA1-V1 (2026-08-20): CC₁ = −0.941
  with **CC₈ = −0.996 on the same trial**, on that animal's LONGEST trials (5773/6448 bins, so
  not a short-trial clip), while its other three pairs read a normal 0.02–0.25. Genuine
  held-out coupling DECAYS with rank; a flat |r| ≈ 1 profile is a rank/collinearity degeneracy
  of the projection, not coupling. 131/15070 rows, and those cells alone inflated the delta SD
  ~10× and made the strength result formally uninformative for exactly those three pairs.
  Gate on the tail: `analyze_trial12.degenerate_trials` drops a cell when mean |r₀| over dims
  6–10 exceeds 0.5.

- **`frozen_perm_null` on a whole-session frozen fit is NOT a dimension selector — it passes
  essentially everything.** Measured 2026-08-20: **100.00 %** of (animal, pair, dim) cells pass
  FS-included and **99.56 %** FS-excluded, with ~99 % sitting exactly at the permutation p-floor
  1/201. So "average over all SIGNIFICANT CCs" is arithmetically "average over all 10 CCs", and
  any sentence contrasting the two is describing a distinction the data does not contain. The
  null is well calibrated under a total null; what it cannot do is separate ranks when hundreds
  of thousands of train bins make every rank beat a circular shift. Report the pass RATE
  whenever a `sig` gate is invoked — and note that a low per-animal-pair `n_sig` usually means
  FEWER DIMS WERE AVAILABLE, not that dims failed the test.

- **Cumulative integration windows cannot be tested with a per-window BH — it manufactures
  survivors.** Windows are nested (|lag| ≤ w), so 25 per-window tests are near-duplicates:
  applying BH across pairs at each window separately yields "5 BH-surviving" window results on
  the trial-1-vs-2 data, none of which survives the defensible statistic — **max |t| across
  windows against a sign-flip permutation null on each animal's whole delta-vs-window profile**,
  which preserves the nesting exactly (`analyze_trial12.sweep_test`, 0 BH survivors either FS).
  Related: seed the permutation PER GROUP (`zlib.crc32` of the group key, never `hash()`, which
  Python randomises per process) — a single stream consumed in group order silently changes
  every later group's p when any upstream row is dropped, and moved two results across 0.05.

- **On the ALL-TRIALS frozen fit (`run_cc_label_track`), curve HEIGHT is lower in the naive
  epoch at EVERY lag — a baseline offset, not learning.** Found 2026-08-15, adversarially
  verified: peak r "rises" expert − naive 8/8 pairs (both FS), but ~60 % of it is present at
  |lag| ≥ 200 ms (slow co-modulation on whole-session axes; speed is the candidate); peak
  MINUS baseline expert − naive is null in every pair with n > 4; the rise is LP-independent
  (ρ ≈ 0); and the held-out per-epoch REFIT arms show no naive deficit (`epoch_metrics_bin10`
  CC₁ 0.149/0.139/0.145). Consequence: `cc_label_track`'s in-sample `r`, `peak_r`, curve
  height — anything with units of correlation — must not be read as a naive→expert change.
  Report `peak_minus_far` next to `peak_r` (`analyze_cc_crosscorr_epochs.curve_metrics`).
  IFI is insensitive to scaling but NOT to an additive offset — fine only when null. Item 2's
  `cc_label_track_stats` per-label `p_*_peak_r` columns carry this caveat.

- **Cumulative integration windows make a sign-split at ±50 ms circular at EVERY window.**
  `HCV1_cc_ifi_windows_*` groups CCs by IFI sign at |lag| ≤ 50 ms and plots them out to
  ±250 ms; every wider window still contains those lags, which carry the largest CCs, so the
  groups "staying apart" out to 250 ms is the same circularity leaking outward — on the
  disjoint lags 60–250 ms the +/− gap collapses to ~0 in every pair (2026-08-15). Only a
  disjoint-lag comparison is non-circular.

- **A driver's sample cap can silently turn a session question into a first-20-trials
  question.** `run_lag_curves` defaulted to 12 000 bins = the first ≤ 600 bins of the first
  ~20 trials from 2026-07-28 to 2026-08-17; item 1 and the 08-07 asks 1/3 were answered on
  that without anyone saying so, and their CC₁ disagreed with the whole-session §B at
  r ≈ 0.4 per animal-pair (magnitudes 3–5× inflated). Since 2026-08-17 the default is
  uncapped and CC₁ reproduces §B to 3 dp. `run_lag_subspaces` / `run_lag_cosine` still cap
  — check `_capped_index` before trusting any "session-level" claim from a capped driver.

- **Never hard-code a verdict string next to a test.** "Sign-mixing is indistinguishable
  from a coin flip" and "nothing survives FDR" were fixed strings in `analyze_cc_ifi_signs`;
  the whole-session data flipped both and the md would have kept printing the old verdict.
  Same class as the interpretive suptitles the verifiers flagged. Derive the sentence from
  the numbers (`_mixing_verdict`) or leave interpretation to the log.

- **A thin SVD of the n×k data matrix is the wrong CCA at n ~ 3e5.** `core.cca_fit` took
  1.66 s per fit at 300 k rows (×1255 fits per pair = 35 min/pair, >12 h per run); the
  k×k covariance route gives identical results in 97 ms. Fixed 2026-08-17 (`d3ffec4`),
  pinned to `_cca_fit_svd`. If a run is unexpectedly slow, profile before parallelising —
  the first batch spent more time in `system` (BLAS thread spin) than in `user`.

- **A figure's SEM band and its title's p-value must be the same unit.** The +/− group
  lines on `HCV1_cc_ifi_windows_*` were mean ± SEM over (animal, dim) rows — dims-as-n —
  until 2026-08-15. Collapse each animal's CCs first (`cc_aggregate.per_animal_mean`),
  then average across animals; the means barely move but the bands are honest.

- **Canonical dimensions are NOT comparable across two CCA fits — 82 % of the time the
  best match at another lag is a different dimension.** Measured (`run_lag_cosine`,
  2026-08-07): split-half |cos| of a canonical vector at a FIXED lag is 0.59 (CC1), 0.39,
  0.26, 0.22 by CC4, against a 0.146 random-vector baseline in 30-D. Only CC1 has real
  identity. Consequence: never attach one fit's per-dim statistics to another fit's
  dimensions by bare index. This caused two separate bugs — `d38a833` (significance from
  `window_subspace` attached to `heldout_lag_curve_flat_perdim`'s dims; 19 % of "significant"
  dims had a NEGATIVE held-out CC) and `1bae90e` (per-fold refit statistics averaged BY RANK
  then attached to a frozen fit; `cc_heldout` rose with rank in 38/38 cells). Compute
  significance in the same driver, on the same fit, from the same scores — or use frozen
  axes, where no correspondence is needed.

- **A permutation p cannot go below 1/(n_shuffles+1), so BH can be arithmetically
  impossible.** BH needs the best of `d` tests to reach `alpha/d`; with d = 30 that requires
  **> 599 shuffles before ANY dimension can pass**. Below that the corrected mask is empty
  for arithmetic reasons and reads as a scientific null. Restrict the BH family
  (`fdr_dims`) or raise the shuffle count, and never interpret an empty FDR mask without
  checking the floor first.

- **Mask monotonicity is a property of the NULL, not of correctness.** Under a dominant-dim
  null (one scalar threshold) the significance mask must be monotone in the held-out CC.
  Under a per-dim null each dimension has its own bar and shuffled correlations fall with
  rank, so a smaller high-rank CC can legitimately pass where a larger low-rank one does
  not. A guard asserting monotonicity will false-alarm on correct per-dim data (it did).

- **A noise floor must be the SAME ESTIMATOR as the thing it gates.** Comparing a full-data
  cross-lag statistic to a half-data split-half floor mixes two noise levels. At 3 canonical
  dims the split-half floor is ~78 deg while the comparison angle is ~75 deg, so "not above
  floor" meant UNMEASURABLE, not stable — reported as a null before being caught, and the
  corrected version reversed the item-3 conclusion. For the cosine analysis both terms are
  half-data fits so the lag-0 point already carries the sampling noise.

- **Check the chance model before calling a rate a finding.** Two live examples: per-animal
  sign-mixing of 92 % was BELOW its own chance level (1 - 2*0.5^k = 94 % at k = 5); and
  label persistence of "70 % vs 50 %, p = 1e-11" had a construction floor of ~60 %, because
  the label was computed on data containing the epoch being scored (analytic
  0.5 + arcsin(sqrt(n_epoch/n_total))/pi, Monte-Carlo confirmed). Score against a
  leave-epoch-out label, or against the construction floor.

- **The 12k-bin sample cap keeps only the first ~20 trials, and the first 0.4-1.8 s of
  each.** A trial is ~2900 bins (~30 s of running). For a late learner that window is
  entirely pre-learning, and running speed rises **+12.0 cm/s naive->expert on those
  trial-onset bins** versus +6.6 cm/s over whole trials — so the cap DOUBLES a speed
  confound on any epoch contrast. The cap exists because `run_lag_curves` refits CCA 255
  times per pair; a driver that fits once (`run_cc_label_track`) should not inherit it.

- **A subspace-angle test at 3 dims is UNMEASURABLE here — check the split-half floor
  before reading any angle p-value.** The floor (two halves of the *same* data, *same*
  lag/window) is **~79°** at d=3 and **~56°** at d=1. At d=3 the comparison angle is ~75°,
  so it has nowhere to go: "not significantly above floor" means the estimate is not
  reproducible, NOT that the subspace is stable. This is the same ~80° property noted for
  the rotation-null (a ~1-D subspace, not a bug), but the consequence is stronger than
  recorded — it silently converts a power failure into an apparent null. Run cross-lag /
  cross-window subspace comparisons at **d=1**, and carry an `estimable` flag (floor
  < 70°) that is read before the p-value. Found 2026-07-29: at d=3 the FF/FB gate said
  2/8 pairs "separable"; at d=1 it says 0/8.

- **An unweighted L2 row-norm over CCA weights is partner-invariant — it is NOT a
  communication-subspace readout.** `core.cca_fit` returns `A = Vx @ diag(1/sx) @ Uc[:, :d] *
  scale`. If you take `norm(A[i, :])` across **all** `d` retained dims and `d = rank(X) ≤ rank(Y)`,
  the square-orthogonal `Uc` cancels out of every row norm exactly, so the partner area
  disappears from the arithmetic: the same neuron gets the same "contribution" whether its partner
  is perfectly coupled or pure noise (max diff 4×10⁻¹⁵ on synthetic data; median r = 0.981 across
  five real partners). Any Gini/membership metric built on it measures the population's own
  whitened-PCA loading geometry. Use `membership.subspace_contribution_connection` (CC-weighted,
  `gini_*_conn`) or restrict to significant dims (`gini_*_sig`) when you mean *connection*.
  Found 2026-07-28 — it silently underpinned the §3.0 "participation broadens" headline.

- **Capping samples by a contiguous/stride slice breaks by-trial CV and lag adjacency.**
  These sessions run to ~370k engaged 10 ms bins with ~2700 bins/trial, so `idx[:MAX_SAMPLES]`
  (first-N contiguous) spans only 3–4 trials → fails any ≥5-fold by-trial CV gate (the
  `run_trajectory_bins` 0-rows bug, 2026-06-17). But an *even-stride* subsample is the
  opposite trap: it destroys within-trial bin adjacency, so segment-aware lag pairing
  (`lagged._segment_lagged_pairs`) silently mis-pairs. Correct cap = keep a CONSECUTIVE
  within-trial block of WHOLE trials until ~MAX_SAMPLES (`run_lag_curves._capped_index`;
  `run_trajectory_bins` uses even-stride only because it does no lagging). Choose the cap by
  what the downstream step needs: trial count (CV) vs bin adjacency (lag).

- **CSV files use CRLF line endings.** The `csv` module writes `\r\n` by default even
  with `open(..., newline="")`. So `learning_changes_*.csv`, `landmark_prune_*.csv`, and
  `sweep_landmark_summary.csv` have a trailing `\r` on every field-10 value. A naive
  `awk -F, '$10=="True"'` matches **nothing** (it compares against `"True\r"`). Strip it
  first: `awk -F, '{gsub(/\r/,"")} ...'`, or read with pandas (handles it). Fixed going
  forward by passing `lineterminator="\n"` to the CSV writers (2026-06-05); already-written
  CSVs keep their CRLF until regenerated.

- **MATLAB `string`-type fields are unreadable by h5py.** `TF*_export.mat` region labels /
  learning points stored as MATLAB `string` cannot be read directly. Worked around with the
  one-off `scripts/export_cca_labels.m`, which writes `cca_labels.json` (schema
  `tom_cca_labels_v1`) that `dataio.py` reads as a companion file.

- **Judge landmark-config overfitting by `frac_cc_ge_099_*`, not `max_cc`.** A single
  saturated canonical dim pushes `max_cc` to 0.999 in otherwise-healthy configs. See
  `STATE.md` §4 and `src/tom_cca/prune_table.py`.

- **A single trial cannot re-fit pCCA/KCCA.** One trial ≈ 580 running 50 ms bins; a
  30-component fit needs ~50×30 ≈ 1500 samples. For per-trial resolution, PROJECT the
  trial onto a subspace fit on many other trials (`early_trials.reference_fit`), don't
  re-fit; fit-only metrics (weight-Gini, angles, #sig, KCCA) need ≥5-trial blocks. See
  `early_trials.py` and report §2.10c.

- **Don't nest `nohup … &` inside a `run_in_background` Bash call.** It double-backgrounds:
  the harness marks the *wrapper* complete (after the 2 s echo/sleep) while the real Python
  run keeps going untracked, so you never get a true completion signal. Launch the bare
  `python …` command with `run_in_background` (no `&`/nohup), or poll with `pgrep`.

- **The ~80° top-3 split-half "noise floor" is NOT a bug — it means the subspace is ~1-D.**
  `subspace_window.split_half_x` is the *max* principal angle over the top-3 canonical dims.
  When only CC1 carries a stable direction (dims 2–3 have cc≈0, i.e. no real structure), the
  max angle is dominated by the random noise dims → near-orthogonal. Verified: 3 genuinely
  shared dims → top-3 floor ~9°; 1 shared dim → CC1 floor ~7° but top-3 floor ~85°
  (`test_subspace_window.py::test_split_half_floor_tracks_true_dimensionality`). **Use the
  CC1-only floor (`split_half_x_cc1`) to judge dominant-direction stability**, never the
  top-3 max angle when n_sig is low. (report §3.4)

- **Running `run_trajectory` truncates its window CSV mid-run (`w` mode).** The driver opens
  `trajectory_w15_bin{25,50}{,_fsincl}.csv` in `"w"` and re-writes after each animal, so while a stage
  is in flight the file holds only the animals done so far. Regenerating trajectory figures against it
  mid-run silently builds them from 1–2 animals (caught 2026-06-12: the 23:56 fsexcl figures). **Only
  run `figs_report.py <traj-stem>` when the window CSV shows 16 animals**; to plot full data while a
  re-run is in flight, restore from git (`git show HEAD:cca/results/<stem>.csv > <stem>_safe.csv`) and
  point figs at the protected stem. The per-dim `_dims.csv` truncates the same way.

- **Too many `run_in_background` tasks → the harness evicts (kills) one.** Launching two batch
  runs as background tasks *plus* a watcher (3 tracked tasks) got a batch killed ~14 min in
  (2026-06-12) — looked like a crash but was eviction (RAM was 87 % free, not OOM). For long
  unattended compute, launch it **fully detached** so it isn't a harness task:
  `(nohup bash run.sh >/dev/null 2>&1 </dev/null &)` → reparents to launchd (PPID 1, new session),
  survives. Keep at most ~1 `run_in_background` watcher alongside.

- **`_load_temporal_streams` cache key must include every cfg field that changes the binned
  output.** It was keyed on `(animal_id, temporal_bin_ms)` only; adding `gaussian_sd_ms` (spike
  smoothing) without adding it to the key meant a second load with a different σ returned the
  STALE cached array — a smoothing validation silently showed "no effect" (identical CC) until the
  key was fixed (2026-06-13). Any new cfg field that alters `spikes_50ms` (smoothing, rebin mode,
  unit selection) must be added to the cache key in `dataio._load_temporal_streams`.

- **A "ringing" detector with no calibrated null measures nothing.** `fixed_subspace.side_peak`
  fires on **~100 % of pure-noise curves** at the `ratio > 0.5` gate — on a curve whose central
  peak is near the noise floor almost any bin qualifies as a secondary maximum. So the FRACTION
  of curves showing a side peak is uninformative and must never be quoted as evidence of
  rhythmicity (it was, on 2026-07-29). What carries information is WHERE the side peaks land:
  a real oscillation concentrates them at its period, noise spreads them across the grid. Use
  `side_peak_null` + `band_occupancy` and compare distributions (real median 140 ms, IQR
  120–148 vs noise median 170 ms, IQR 90–270; KS p = 2.7×10⁻¹¹). Found 2026-08-03.

- **Weights fitted in one residual space and applied in another silently break a "frozen"
  projection.** `run_fixed_subspace.py` fitted the subspace on data residualised with
  *fit-trial* confound coefficients, then projected data residualised with *all-trial*
  coefficients — so "the identical weights for every epoch" was not an identical transform.
  The residualisation is part of the map and must travel with it. This one manufactured the
  only significant result in the arm (CA1-RSC integration-window narrowing, p = 0.029), which
  vanished on fixing it. Found 2026-08-03 by adversarial verification, not by testing.

- **Gating a per-area statistic on the POOLED mean of two areas silently readmits the thing
  you gated out.** `analyze_lag_subspaces.stability` averaged the X and Y `angle − floor`
  terms per animal and tested estimability on the mean of the two floors, so an area with no
  usable subspace estimate still contributed. 20/71 animal-pairs had BOTH d=1 floors over the
  70° gate while the pooled test called 48/71 estimable. Direction of the bias is knowable:
  corr(floor, angle − floor) = ρ −0.57 (p ~ 1e-242), so unmeasurable areas contribute large
  negative deltas and pull toward the null. **Gate per area, then aggregate.** This one
  converted a real effect (CA1-DG cross-lag rotation) into a reported null. Found 2026-08-03.

- **A split-half floor is NOT a monotone measure of how well an area is estimated.** A
  near-rank-1 population has degenerate residual PC directions, which INFLATES its d=1
  split-half angle: a clean 4-unit area scores 34.7° while a noisy 30-unit one scores 11.3°
  (`test_floor_returns_x_area_first_not_y`). So a floor-based "estimability" gate excludes
  some well-conditioned areas and keeps some poor ones, and any count derived from such a
  gate is an analyst choice. Report the whole gate sweep and claim only what survives all of
  it. Found 2026-08-03.

- **A metric that is NaN under a condition silently conditions its own slope on that
  condition.** `gini_*_sig` (contribution Gini over significant dims only) is NaN whenever
  `n_sig` = 0 — **56 % of CA1-RSC windows**. Every per-animal slope and LMM on it is therefore
  computed on the subset of windows that HAD a significant dim, and the animal count can silently
  drop (CA1-RSC FS-excl: 8 animals for `gini_x`, 7 for `gini_x_sig`). If `n_sig` itself drifts
  over the session, that drift is inside the `_sig` slope. **Read `_sig` only beside the `n_sig`
  slope, and report both animal counts.** Found 2026-08-28 in the connection-specific Gini
  re-test.

- **An LMM over windows will star at n = 4 animals where the animals-as-n test mathematically
  cannot.** The two-sided Wilcoxon floor at n = 4 is p = 0.125, so CA1-SUB / CA3-DG / RSC-SUB can
  never reach p < 0.05 at the honest unit — yet their LMM population slopes go to p = 6×10⁻⁶ by
  pooling windows within those same 4 animals. The LMM's random slope does not buy back animals
  that were never there. **At n = 4, treat any LMM p as descriptive only.** Found 2026-08-28.

