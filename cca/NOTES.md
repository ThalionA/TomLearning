# Tom-learning CCA -- development notes

Running log. The *spec* lives in `UNDERSTANDING.md` (decisions + edit log);
this file is the chronological work log.

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
