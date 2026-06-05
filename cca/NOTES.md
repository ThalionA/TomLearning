# Tom-learning CCA -- development notes

Running log. The *spec* lives in `UNDERSTANDING.md` (decisions + edit log);
the *current state, canonical configs, and findings* live in `STATE.md`
(read that first); this file is the chronological work log.

---

## 2026-06-05 -- consolidation pass: STATE.md + verdict

Wrote `STATE.md` reconciling the spatial and landmark arms. Settled the
learning verdict from the existing `learning_changes_*.csv`: the one robust,
reproducible effect is **CA3-DG coupling strengthening at the expert stage**
(committed `landmark50_res_samp15`, pooled scope: expert>naive Δ+0.21 p=0.006
FDR-pass; expert>intermediate Δ+0.15 p=0.004 FDR-pass; reproducible in 18/40
non-overfit landmark configs). Nothing survives the per-landmark (48-test)
FDR family in any config. Excluded the 4 overfit configs
(`landmark50_res_{fix30,var75,var85,var95}`). Found + fixed a CRLF gotcha in
the learning-changes CSVs (see `GOTCHAS.md`). Spatial arm has no paired
learning test -- flagged as an open decision in `STATE.md` §5.

---

## 2026-05-24 -- scaffold: port the striatum pipeline to Tom's cohort

Theo asked to apply the new CCA analysis we developed for striatum to the
Tom-learning dataset, spatial-only, on the pairs Tom previously defined,
with FS excluded, sweeping the other parameters. New folder in the main
TomLearning repo; the data files are not in the container, so this round is
**scaffold + tests only** -- the pipeline is ready to run once the
``HC_V1_data/`` directory is mounted.

### What was done

* New package ``tom_cca`` (`cca/src/tom_cca/`, 13 modules) -- a direct port of
  ``striatum_cca`` from `StriatumACC/Striatum project/cca/`. Pure-numerical
  modules (``core, lagged, surrogate, membership, subspace, partial,
  crosspair, analysis, stage3``) are byte-identical; ``config, dataio,
  pipeline, sweep`` are Tom-specific.
* ``dataio.py`` reads Tom's per-animal exports (one ``TF*_export.mat`` per
  animal, vs the cohort-in-one-file striatum format), with
  ``analysis_spatial.firing.cued.freq`` as the spatial firing source.
  Learning points come from ``animal_behaviour.mat``
  (``period_experienced(:, 1)`` indexed by ``animal_id``); per-trial detection
  is the fallback when a cohort file is absent.
* ``config.py`` defines Tom's six areas (CA1/V1/DG/CA3/RSC/SUB) and the eight
  hard-coded pairs from Tom's MATLAB scripts (CA1-V1, CA1-DG, CA1-CA3,
  CA1-RSC, CA1-SUB, V1-RSC, RSC-SUB, CA3-DG). FS exclusion is restricted to
  V1/RSC/CA1/CA3 (Tom's convention -- DG/SUB have no FS flags), via
  ``config.FS_AREAS``.
* ``sweep.py`` drops the ``AXIS_BINS`` and ``AXIS_FS`` axes from the striatum
  sweep -- bin width is fixed at 200 (Tom's native preprocessing) and FS is
  fixed-excluded (the user's instruction). The remaining grid is:
  CCA type x z-scoring x min-units x LP criterion x PC-count rule
  = 2 x 2 x 3 x 2 x 11 = **264 configs**.
* All 14 striatum scripts ported -- ``run_stage2/3/committed/partial``,
  ``plot_stage2/3/partial/common_units/ifi_fs/parcoords/subspace_similarity``,
  ``compare_nulls``, ``committed_ifi``, ``summarise_sweep``. The pair-grid
  scripts changed from a (2 x 5 = 10 pair) layout to (2 x 4 = 8 pair); the
  ``axes[::5]/[5:]`` indexing changed to ``axes[::4]/[4:]``. The
  ``summarise_sweep`` parameter columns dropped ``bin`` and ``fs``; the
  ``committed_ifi`` committed-config tag was updated to
  ``res_z1_mu06_lp7_samp15`` (the Tom committed tag, without the bin and FS
  components). The ``--temporal`` branch was removed from ``summarise_sweep``
  and ``run_stage2/3`` since Tom is spatial-only.
* Test suite ported and extended -- **109 tests, all passing**. The
  pure-numerical tests (72) are byte-identical to the striatum tests with
  ``striatum_cca`` -> ``tom_cca``; ``test_pipeline.py``, ``test_dataio.py``,
  ``test_analysis.py``, ``test_stage3.py`` were rewritten for Tom's
  ``Animal`` constructor (no ``neurontypes`` or ``change_point`` -- a
  per-unit ``fs_mask`` instead) and the per-animal mat-file layout.
  ``test_dataio.py`` includes a .mat round-trip test that writes a tiny
  synthetic ``TF*_export.mat`` + ``animal_behaviour.mat`` and loads them
  back.
* End-to-end smoke test on synthetic Tom-shape data (two animals, 100
  trials each, 200 bins, CA1 + V1) successfully ran ``prepare_pair`` and
  ``analyse_pair`` -- produced PC-score tensors of shape
  ``(n_trials, n_bins, k) = (10, 200, 10)``, surrogate p-values, IFI per
  dimension, peak lags. The held-out CCs were ~0 (the synthetic data has
  no shared structure, as expected).
* ruff clean.

### Decisions made without Theo (flag if you disagree)

* **manual_nonlearners empty for Tom.** No analogue of the striatum
  animal-8 LP-detection artefact is known. Set
  ``cfg.manual_nonlearners`` after the first cohort run if any animal's
  detected LP is implausible.
* **Bin-width axis dropped from the sweep** (per the user's "Tom's native
  200 bins, no rebinning"). The bin-count check inside the pipeline is still
  in place -- it auto-detects the bin count from the data file, so a 50-bin
  fixture in tests works -- so this is the sweep-grid restriction, not a
  hard ceiling.
* **Spatial source field = ``freq``** (raw Hz), not ``freq_z`` (which is
  pre-z-scored). Reason: the sweep includes a ``z-scoring on/off`` axis, and
  applying our own z-scoring to ``freq_z`` would be a no-op for the
  ``z=on`` arm and a re-normalisation for the ``z=off`` arm -- both
  uninterpretable. Override via ``cfg.spatial_field = "freq_z"`` if a
  MATLAB-parity comparison is needed.
* **No `RESULTS.md` yet.** Will be written after the first cohort run, the
  way `Striatum project/cca/RESULTS.md` was. The current writeup is
  ``UNDERSTANDING.md`` -- the design + decisions; this file is the work log.

---

## 2026-05-24 -- commit pass 1: z, min_units, LP, learners-only, IFI window

Theo's three commits (single message): "commit to z-scoring. commit to 5 min
units. commit to whatever learning point tom is giving" -- followed by "Only
work with the learners. Another parameter is the window for IFI."

### What changed

* ``config.Config`` -- ``min_units`` 6 -> 5; ``zscore_units`` stays True; the
  LP-detection params (``lp_z_threshold``, ``lp_window``,
  ``lp_min_consecutive``) are deleted.
* ``dataio`` -- ``find_learning_point`` deleted; ``classify_cohort`` simplified
  to return a single dict of LEARNER ``CohortEntry`` (no yoked non-learners;
  no fallback detection). LP comes only from ``animal_behaviour.mat``.
  ``Animal.zscored_lick_errors`` removed (no longer needed).
  ``_read_behaviour_file`` no longer extracts the optional per-trial z trace.
* ``sweep`` -- axes pared back to **CCA type x PC-count rule x IFI lag
  window**:
    - AXIS_Z, AXIS_MIN_UNITS, AXIS_LP all removed
    - AXIS_LAG = (5, 10, 20) bins added -- +/-12.5 / 25 / 50 cm scan
  Grid is now **66 configs** (was 264). Tag is ``{cca}_{krule}_lag{NN}``.
* Runner scripts -- ``run_committed``, ``run_stage2``, ``run_stage3``,
  ``run_partial`` all updated for the new ``classify_cohort`` return type and
  filter ``animal.animal_id not in entries`` to skip non-learners.
* ``summarise_sweep`` -- ``PARAM_COLS`` reduced to ``["tag", "cca", "k_rule",
  "max_lag", "pair"]``; ``cfg_params`` drops the constant z / min_units /
  lp_consec columns and gains ``max_lag``.
* ``committed_ifi`` -- ``COMMITTED_TAG`` updated to
  ``res_samp15_lag{max_lag_bins:02d}`` (derived from ``config.DEFAULT``).
* ``test_dataio`` -- LP-detection tests deleted; cohort-classification tests
  rewritten for the learners-only contract; new
  ``test_committed_defaults`` pins z=True, min_units=5, FS-excluded, and
  asserts the LP-detection knobs are gone.
* ``test_pipeline`` -- ``synthetic_animal`` builder drops the
  ``zscored_lick_errors`` field.
* All 107 tests pass (down from 109 after removing the four obsolete
  LP-detection tests + adding the new committed-defaults check). ruff clean.

### IFI window note

The IFI is computed for every sub-window of the lag scan in one fit (via
``lagged.ifi_by_window``), so a config with ``max_lag_bins=10`` already
carries IFI at all of ``w = 1..10``. Sweeping ``max_lag_bins`` therefore
varies the *maximum* scan range (and the lag-curve resolution / peak-lag
detection range), not the IFI at any fixed sub-window. The per-pair x
per-config IFI summary in ``summarise_sweep`` reports the IFI at every
sub-window of each config's scan, so the user sees both axes.

### Pending / blockers

* The data directory ``TomLearning/HC_V1_data/`` is empty (only a
  ``.DS_Store``). The cohort run is pending data being mounted in the
  container.
* Once data is mounted, the natural sequence is:
  1. ``python scripts/run_committed.py --stage 2 --data-dir <path>`` -- get
     the committed-config Stage 2 result first to sanity-check the
     pipeline on real data;
  2. ``python scripts/run_committed.py --stage 3 --data-dir <path>``;
  3. ``python scripts/plot_stage2.py`` and ``plot_stage3.py``;
  4. ``python scripts/run_stage2.py --data-dir <path>`` to run the full
     264-config sweep (resumable; supports ``--max-seconds`` chunking);
  5. ``python scripts/run_stage3.py`` and ``summarise_sweep.py`` to build
     the per-pair robustness tables and grids.
* After the first run: check whether any animal's LP looks
  implausible; if so, add it to ``cfg.manual_nonlearners`` as the
  striatum project did for animal 8 (see UNDERSTANDING.md).
* Consider revisiting the committed config after the per-pair x
  per-config sweep grids are available -- the striatum project iterated
  twice on this (round 8 -> round 10).

---

## 2026-05-24 -- first real-data run attempt: MATLAB `string` blocker + fix

The 16 ``TF*_export.mat`` files + ``animal_behaviour.mat`` are now in the
container. Tried to run the full 66-config sweep; hit a data-contract blocker
on the very first animal load.

### What was found

* The ported test suite is healthy -- **112 tests pass** (107 prior + 5 new).
  Note: pytest's default temp cleanup hits the Cowork sandbox's no-delete
  restriction on mounted folders and dies with a ``RecursionError``. Run the
  suite with ``--basetemp`` pointed outside the mount, e.g.
  ``PYTHONPATH=src python -m pytest tests/ -q --basetemp=/tmp/tcca``.
* **Blocker:** ``units.regions_label`` (and ``units.region``,
  ``units.region_details``, and ``animal_behaviour.mat``'s ``animal_id``) are
  stored as MATLAB **``string``** arrays in every real file. h5py cannot
  decode the ``string`` type (nor can the ``mat73`` library -- it errors
  ``MATLAB type not supported: string``). The scaffold's ``dataio.py`` assumed
  ``regions_label`` was a cellstr because the synthetic test fixtures write it
  that way -- the scaffold was built with no real data in the container.
  Without region labels there is no idx-column -> area-name map, so no pair
  can be built and nothing runs.
* ``analysis_spatial.firing.cued.freq`` (200 bins x 78 trials x 130 units),
  ``units.idx`` (logical, 130 x 5) and ``units.idx_fs`` all load fine -- the
  break is purely the ``string``-typed label fields.
* ``animal_behaviour.mat``: ``period_experienced`` is 16 animals x **10
  columns** (the loader's old comment guessed 2). The LP is column 1 per
  Tom's MATLAB pipeline -- but this is worth confirming against the data.

### Fix (chosen with Theo: source-side re-export)

* New ``scripts/export_cca_labels.m`` -- a one-off MATLAB script Theo runs.
  It reads each file's ``string`` fields in MATLAB (where ``string`` is
  native), and writes ``HC_V1_data/cca_labels.json`` (schema
  ``tom_cca_labels_v1``): per-animal region labels in idx-column order, plus
  per-animal learning point (+ all 10 ``period_experienced`` columns for the
  record). ``lp_column`` is a script argument, default 1.
* ``dataio.py`` -- new ``_read_companion`` reads ``cca_labels.json``.
  ``load_animal`` gains a ``region_labels`` override; ``load_animals`` and
  ``_read_behaviour_file`` use the companion when it is present and fall back
  to reading the .mat fields directly when it is absent (keeps the cellstr
  test fixtures working). ``config.COMPANION_FILE = "cca_labels.json"``.
* 5 new ``test_dataio.py`` tests cover the companion path (label override,
  behaviour preference, missing-animal error, end-to-end classify). ruff clean.

### Pending

* **Theo must run ``scripts/export_cca_labels.m`` in MATLAB** (the container
  has no MATLAB) to produce ``HC_V1_data/cca_labels.json``. The sweep is
  blocked until that file exists.
* Then resume: smoke-test one animal -> ``run_stage2.py`` (66-config sweep,
  resumable) -> ``run_stage3.py`` -> ``summarise_sweep.py``.
* Confirm ``period_experienced`` column 1 is the learning point.

---

## 2026-05-24 -- companion file generated; pipeline verified; sweep handed off

* Theo ran ``export_cca_labels.m``; ``HC_V1_data/cca_labels.json`` now exists
  (schema ``tom_cca_labels_v1``, 16 animals).
* **``period_experienced`` resolved.** Its 10 columns are consecutive trial
  indices (LP, LP+1, ..., LP+9), so column 1 *is* the learning point --
  ``lp_column=1`` confirmed correct. 12 of 16 animals have a recorded LP;
  TF070/071/098/100 are NaN -> dropped as non-learners (learners-only spec).
  Cohort for all runs = **12 learner animals**.
* **End-to-end smoke test passed on real data.** TF036 loaded via the
  companion path (78 trials x 200 bins x 130 units; CA1 39, V1 29, CA3 12,
  DG 2 units), classified learner lp=28, ``fit_pair`` on CA1-V1 succeeded
  (k=26; CA1 32 / V1 26 units post-FS). Held-out CC1 ~0.01-0.09 vs in-sample
  ~0.27-0.36 -- CV behaving correctly, no red flags.
* **Sweep throughput measured:** ~0.25 (animal x pair x config) fits/sec on
  4 cores; full 66-config Stage 2 ~= 3300 fits ~= 4-6 h of compute.
* **Sweep handed to Theo for native execution.** The Cowork sandbox caps each
  shell command at 45 s and kills background processes between calls, so a
  multi-hour batch job is impractical here. The pipeline is fixed and proven;
  Theo runs ``run_stage2.py`` / ``run_stage3.py`` / ``summarise_sweep.py``
  natively (resumable, one command each). ``results/`` was left clean; the
  partial pkl from the throughput test is parked in
  ``results/.sandbox_scratch/`` (a native run starts fresh).

---

## 2026-05-25 -- sweep complete; per-config plotting

* Theo ran the full sweep natively: 66 ``stage2_*`` + 66 ``stage3_*`` pkls,
  all ``.done``; ``summarise_sweep.py`` outputs in ``figures/``.
* ``plot_stage2.py`` / ``plot_stage3.py`` gained an additive ``--tag`` arg:
  with ``--tag <sweep_tag>`` they plot ``results/stage{2,3}_<tag>.pkl`` and
  write into ``figures/<tag>/``. Without ``--tag`` behaviour is unchanged
  (committed / partial run). Default-path code untouched.
* Generated the full Stage 2 + Stage 3 figure set for the four var-85 configs
  -- ``{res,sig}_var85_lag{05,10}`` -- into ``figures/<tag>/`` (54 PNGs).
* New ``scripts/plot_sweep_explore.py`` -- config-centric exploration of the
  sweep (reads ``sweep_summary_spatial.csv``; writes ``figures/sweep_explore/``,
  6 PNGs). Key structural finding: the 66 configs collapse to **16 distinct
  headline results** -- the IFI lag-scan window (5/10/20) leaves every headline
  metric unchanged, and ``samp15/25/40`` + ``fix30`` share one effective k.
  CCA type is the dominant axis (signal ~2x CC, almost all significant effects);
  signal x variance-k is the live corner; the committed config sits low
  (~18th percentile on effect size).
* New ``scripts/plot_ifi_explore.py`` -- IFI deep-dive from the per-epoch data
  in the lag20 Stage 2 pkls (dim-0 IFI = the pipeline's headline ``.ifi``).
  4 PNGs in ``figures/ifi_explore/``: IFI across learning per pair, naive->
  expert IFI change, IFI vs lag-window (naive vs expert), peak lead/lag in cm.
  Descriptive only (configs are not independent samples).
