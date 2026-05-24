# Tom-learning CCA pipeline

Spatial CCA pipeline for inter-areal communication across spatial learning in
the Tom (hippocampal / cortical) cohort. Direct Python port of the striatum
CCA pipeline (`StriatumACC/Striatum project/cca/`); same analyses, same code
backbone, Tom's data layout + pair list.

Eight area pairs (the list hardcoded in
`HC_V1_Code/HC_V1_temporal.m` and the legacy spatial-CCA scripts):

  CA1-V1, CA1-DG, CA1-CA3, CA1-RSC, CA1-SUB, V1-RSC, RSC-SUB, CA3-DG

Spatial-only; 200 bins x 2.5 cm. Fast-spiking units excluded (Tom's
convention: V1/RSC/CA1/CA3 only). Sweep covers CCA type, z-scoring,
min-units, LP criterion, PC-count rule -- **264 configs**.

## Data layout

The pipeline reads:

    HC_V1_data/
      TF<id>_export.mat       # per-animal -- units.{idx,regions_label,idx_fs}
                              #               + analysis_spatial.firing.cued.freq
      animal_behaviour.mat    # cohort -- period_experienced(:, 1) + animal_id

See `UNDERSTANDING.md` for the full data contract.

## Quick start

The package is a self-contained src-layout Python package. Tests use pytest;
the package itself depends only on numpy, scipy, h5py, matplotlib and
openpyxl (for the sweep summary xlsx).

    cd cca/
    PYTHONPATH=src python -m pytest tests/ -q        # 109 tests, ~10 s

To run on real data (mount or symlink `HC_V1_data/` first):

    # Committed config, Stage 2 + Stage 3 (fast)
    PYTHONPATH=src python scripts/run_committed.py --stage 2
    PYTHONPATH=src python scripts/run_committed.py --stage 3
    PYTHONPATH=src python scripts/plot_stage2.py
    PYTHONPATH=src python scripts/plot_stage3.py

    # Full sweep (slow -- resumable)
    PYTHONPATH=src python scripts/run_stage2.py --max-seconds 1800
    PYTHONPATH=src python scripts/run_stage3.py
    PYTHONPATH=src python scripts/summarise_sweep.py

Override the data directory:

    PYTHONPATH=src python scripts/run_committed.py --stage 2 \
        --data-dir /path/to/HC_V1_data

Results go to `cca/results/`, figures to `cca/figures/`.

## Layout

    cca/
      README.md
      UNDERSTANDING.md      # design + resolved decisions
      NOTES.md              # development work log
      conftest.py
      src/tom_cca/          # the package
        analysis.py         # Stage 2 driver per (animal, pair)
        config.py           # paths, areas, pairs, FS-area convention, Config
        core.py             # residualise, z-score, PCA, CCA, CV
        crosspair.py        # within-area cross-partner subspace similarity
        dataio.py           # TF*_export.mat loader, cohort/LP/epoch logic
        lagged.py           # lag-slice CCA, IFI, per-window IFI
        membership.py       # structure coeffs, weights, Gini, Jaccard
        partial.py          # LS partial regression in (trial, bin) flat space
        pipeline.py         # prepare_pair, prepare_pair_partial, fit_pair
        stage3.py           # membership + principal angles driver
        subspace.py         # principal angles, split-half noise floor
        surrogate.py        # circshift / trial-permutation nulls
        sweep.py            # Tom sweep grid (264 configs)
      scripts/              # see UNDERSTANDING.md sec. 7
      tests/                # pytest suite (109 tests)
      results/, figures/    # outputs (gitignored or empty initially)
