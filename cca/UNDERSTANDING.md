# Tom-learning CCA -- Shared Understanding

Status: **scaffolded; commit pass 1 done (z, min_units, LP, learners-only,
IFI window); awaiting first cohort run**
Last updated: 2026-05-24

---

## 1. Scientific question

How does communication between hippocampal / cortical areas change across
spatial-learning epochs in the Tom cohort? The four sub-questions mirror the
striatum project (see `StriatumACC/Striatum project/cca/UNDERSTANDING.md`):

1. Do inter-areal connections strengthen or weaken with learning?
2. Does the direction of information flow change?
3. Which units take part in the communication subspaces, and are the
   member sets shared across area pairs?
4. How do the communication subspaces reorient across learning epochs?

Method: residual canonical correlation analysis (CCA) on spatially-binned
Neuropixels activity, fit per (animal x area-pair x epoch). Implementation is a
**direct port of the striatum CCA pipeline** -- the pure-numerical modules
(core, lagged, surrogate, membership, subspace, partial, crosspair) are
byte-identical; only the data-IO and the sweep grid are Tom-specific.

---

## 2. Data reality (binding constraints)

Per Tom's MATLAB pipeline (`HC_V1_Code/HC_V1_temporal.m`,
`HC_V1_Code/legacy/CCA_HC_V1_spatial_v2.m`):

* one ``TF<id>_export.mat`` per animal in ``HC_V1_data/``;
* each file holds:
  - ``units.idx`` (n_areas x n_units) -- area membership mask
  - ``units.regions_label`` -- area name per row
  - ``units.idx_fs`` -- fast-spiking flag per unit
  - ``analysis_spatial.firing.cued.freq`` (n_units, n_trials, n_bins) --
    spatial firing rate (Hz); ``freq_z`` is the same field pre-z-scored
* cohort behaviour in ``animal_behaviour.mat`` -- ``period_experienced(:, 1)``
  is the learning point per animal (matched by ``animal_id``).

Spatial geometry: 200 bins x 2.5 cm = 500 cm corridor.

Areas: ``CA1, V1, DG, CA3, RSC, SUB``. FS flags are only set in
``V1, RSC, CA1, CA3`` (Tom's convention -- DG/SUB are skipped).

**Pairs to analyse** (8, Tom's hardcoded list, byte-identical to
`HC_V1_temporal.m` and the spatial CCA scripts):

  CA1-V1, CA1-DG, CA1-CA3, CA1-RSC, CA1-SUB, V1-RSC, RSC-SUB, CA3-DG

Per-epoch sample budget: 10 trials x 200 bins = **2000 samples**.

---

## 3. Resolved design decisions

These map onto the striatum decisions (D1-D12 in
`Striatum project/cca/UNDERSTANDING.md`); divergences are flagged.

**D1 -- Pair scope: learners only.** All eight Tom pairs, learner animals
only. A learner is an animal with a learning point recorded in
``animal_behaviour.mat`` (``period_experienced(:, 1)`` indexed by
``animal_id``). Animals without a recorded LP are dropped from the cohort
entirely -- no yoked non-learner branch. (User commit pass 1: "Only work
with the learners.")

**D2 -- "Communication" definition: residuals primary.** Subtract each unit's
per-bin trial-mean (within epoch) -> CCA on the trial-to-trial residual
fluctuations. The signal CCA variant is on the sweep grid.

**D3 -- Estimator: PCA -> CCA primary; partial CCA add-on.** Plain canonical
correlation on PCA-reduced residuals for all 8 pairs. Partial CCA conditions
each pair on every other recorded area at once (neuron-level partialling
before the per-epoch PCA -- ``pipeline.prepare_pair_partial``).

**D4 -- Dimensionality.** PCA per area, ``k = floor(n_samples / samples_per_pc)``
by default (samples-per-PC swept 15/25/40, plus fixed-k and variance-explained
modes), capped at the smaller area's unit count, by ``k_cap = 30``, and by the
smallest per-epoch numerical rank. Symmetric within a pair, fixed across the 3
epochs. PCA fitted **per epoch** (not on a shared basis) -- striatum's edit log
v3 reasoning applies here too: a shared basis can leave a component carrying
near-zero variance inside an epoch, making that epoch's CCA ill-conditioned.

**D5 -- Spatial bins: keep 200.** User instruction -- Tom's native 200-bin /
2.5-cm preprocessing, no rebinning. Bin width is **fixed**, not a sweep axis.

**D6 -- Directionality: lagged refit, +/-10-bin scan.** 10 bins x 2.5 cm =
+/-25 cm, matching the striatum pipeline's spatial window in physical units.
Readouts: per-dimension Information Flow Index (IFI) and per-dimension peak
lag, both reported per epoch.

**D7 -- Surrogate null: circshift primary.** Per-trial circular shift of the
bin axis by >= ``circshift_min_bins = 15`` bins (Gonzalez et al.). The
trial-permutation null (H&H) is available via ``--null-type trials``. 200
surrogates, per-dimension held-out-CC test.

**D8 -- Cross-validation: 5-fold over whole trials.** Whole-trial folds prevent
within-trial bin-autocorrelation leakage. Held-out CC reported as primary.

**D9 -- Unit membership: both scores, cross-checked.** Structure coefficients
+ canonical weights, top-quartile members, Gini sparsity, cross-pair Jaccard
overlap, cross-epoch stability. Same machinery as striatum.

**D10 -- Subspace rotation: principal angles, all three transitions.** Three
epochs (naive/intermediate/expert), three transitions (naive->intermediate,
intermediate->expert, naive->expert), within-epoch split-half angle as the
noise floor. ``d_sub = 1`` -- dominant canonical direction only (the
split-half stability floor at d_sub > 1 was already near-orthogonal in the
striatum work, and Tom's per-epoch budget is similar).

**D11 -- Compute: NumPy + CPU multiprocessing.**

**D12 -- Delivery.** Same staging as striatum: design doc -> full pipeline
(TDD, 109 synthetic-ground-truth tests pass) -> staged check-ins -> writeup.

### User-pinned defaults (NOT on the sweep grid)

User commit pass 1 (2026-05-24) -- the following are fixed by instruction:

* **FS exclusion: FIXED ON.** ``units.idx_fs`` applied only to
  ``V1/RSC/CA1/CA3`` (Tom's MATLAB convention).
* **Bin width: FIXED at 200 bins** (Tom's native 2.5 cm).
* **z-scoring: FIXED ON.** Each unit divided by its std over the entire
  engaged period, before epoch slicing.
* **min_units: FIXED at 5** units per area.
* **Learning point: from Tom's ``animal_behaviour.mat`` only.** No
  Python-side LP detection, no ``lp_min_consecutive`` sweep.
* **Cohort: learners only.** No yoked non-learners.

### Self-made decisions

* **Epochs:** naive = trials 1..10; intermediate = (lp-9)..lp;
  expert = (lp+1)..(lp+10).
* **manual_nonlearners:** empty for Tom -- no animal-8 analogue is known.
  Set ``cfg.manual_nonlearners`` if a known LP-detection artefact appears.
* **Spatial firing source:** ``analysis_spatial.firing.cued.freq`` (raw Hz).
  ``freq_z`` (pre-z-scored) can be selected via
  ``cfg.spatial_field = "freq_z"`` for parity with Tom's MATLAB scripts.

---

## 4. Sweep grid

After commit pass 1, ``tom_cca.sweep.build_sweep("spatial")``
(see `src/tom_cca/sweep.py`):

| axis           | values |
|----------------|--------|
| CCA type       | residual / signal |
| PC-count rule  | samples 15/25/40, fixed 3/5/10/20/30, variance 75/85/95 % |
| IFI lag window | 5 / 10 / 20 bins  (+/-12.5 / 25 / 50 cm) |

Total: 2 x 11 x 3 = **66 configs**. The tag format is
``{cca}_{krule}_lag{NN}`` (e.g. ``res_samp15_lag10`` is the committed default).

The IFI is reported at every *sub*-window of its scan, so every config
implicitly carries IFI at all of ``1..max_lag_bins`` -- the sweep axis
controls the maximum scan range and the lag-curve resolution.

Bin width (fixed 200), FS (fixed-excluded), z-scoring (fixed on), min-units
(fixed 5) and LP criterion (use Tom's recorded LP) are not on the grid, per
the user's commits.

---

## 5. Won't-Do (explicitly out of scope)

* **Temporal CCA** -- spatial bins only (user's instruction). The striatum
  pipeline's temporal arm is not ported.
* **Fast-spiking cells** -- fixed-excluded (user's instruction).
* **Bin-width rebinning** -- fixed at Tom's native 200 bins (user's
  instruction).
* **Reduced-rank regression** -- CCA only.
* **GPU backend** -- not used (matrices small).
* **Dark / ITI periods** -- corridor traversal only.

---

## 6. Open items

1. **Data files** are not in the container -- the cohort run is pending. The
   data lives in ``HC_V1_data/`` per Tom's MATLAB scripts; ``run_stage2.py``
   accepts ``--data-dir`` to override the default location.
2. **Per-pair learner counts** are unknown until the first cohort run.
3. **Committed config for Tom** (``config.DEFAULT``) after commit pass 1:
   residual CCA, z-scoring on, ``min_units = 5``, ``samples_per_pc = 15``,
   ``max_lag_bins = 10``, ``n_shuffles = 200``, learners only, LP from
   ``animal_behaviour.mat``. May need revising in light of the cohort-level
   results.

---

## 7. Pipeline architecture

Mirrors the striatum pipeline; see `src/tom_cca/`:

* ``config.py``    -- Tom areas, pairs, paths, ``Config`` dataclass
* ``core.py``      -- residualise, z-score, PCA, CCA, CV (identical)
* ``dataio.py``    -- ``TF*_export.mat`` loader, cohort/LP/epoch logic
* ``pipeline.py``  -- ``prepare_pair``, ``prepare_pair_partial``, ``fit_pair``
* ``lagged.py``    -- lag-slice CCA, IFI by window (identical)
* ``surrogate.py`` -- circshift / trial-permutation nulls (identical)
* ``analysis.py``  -- Stage 2 driver per (animal, pair)
* ``stage3.py``    -- Stage 3 driver: membership, Gini, principal angles
* ``subspace.py``  -- principal angles, split-half noise floor (identical)
* ``membership.py``-- structure coefficients, weights, Gini, Jaccard (identical)
* ``partial.py``   -- LS partial regression in (trial, bin) flat space (identical)
* ``crosspair.py`` -- within-area cross-partner subspace similarity (identical)
* ``sweep.py``     -- Tom 264-config grid

Scripts (``scripts/``):

* ``run_stage2.py``       -- sweep driver, all 264 configs, resumable
* ``run_stage3.py``       -- sweep driver (Stage 3, null-independent)
* ``run_committed.py``    -- single committed config, --stage {2,3}, --partial,
                              --null-type, --include-fs, --trials-per-epoch
* ``run_partial.py``      -- partial CCA on all 8 pairs (regress out every
                              other recorded area)
* ``compare_nulls.py``    -- circshift vs trial-perm null comparison
* ``committed_ifi.py``    -- per-pair IFI by lag window at the committed config
* ``summarise_sweep.py``  -- per-pair x per-config xlsx + grids
* ``plot_stage2.py``      -- CC strength, subspace dim, lag curves, IFI per win
* ``plot_stage3.py``      -- principal angles, Gini, membership overlap
* ``plot_common_units.py``-- member vs non-member spatial activity profiles
* ``plot_partial.py``     -- plain vs partial CC1 per pair x epoch
* ``plot_subspace_similarity.py`` -- within-area cross-partner heatmaps
* ``plot_ifi_fs.py``      -- IFI(w10) FS-excluded vs FS-included
* ``plot_parcoords.py``   -- parallel-coordinates per-pair across sweep

---

## 8. Data layout for runs

The pipeline expects:

    HC_V1_data/
      TF01_export.mat
      TF02_export.mat
      ...
      animal_behaviour.mat

To run with a custom location (e.g. an NFS mount in the container):

    python scripts/run_committed.py --stage 2 --data-dir /path/to/HC_V1_data
    python scripts/run_committed.py --stage 3 --data-dir /path/to/HC_V1_data
    python scripts/run_stage2.py --data-dir /path/to/HC_V1_data
    python scripts/run_stage3.py --data-dir /path/to/HC_V1_data

Results land in ``cca/results/``; figures in ``cca/figures/``.
