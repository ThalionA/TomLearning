# Tom-learning CCA pipeline

> **Continuing the temporal/lag analyses?** Read **`HANDOFF.md`** first — pipeline
> map (script → CSV → figure), column dictionary, which numbers are trustworthy, and the
> traps. Then `PROJECT_LOG.md` (state of play) and `STATE.md` (findings).


> **Start with [`STATE.md`](STATE.md)** — it reconciles the two analysis arms
> (spatial + landmark), names the canonical configs, and states the current findings.
> This README is quick-start only.

CCA pipeline for inter-areal communication across spatial learning in
the Tom (hippocampal / cortical) cohort. Direct Python port of the striatum
CCA pipeline (`StriatumACC/Striatum project/cca/`); same analyses, same code
backbone, Tom's data layout + pair list. Two arms now exist: a **spatial**
full-corridor sweep (66 configs) and a **landmark** event-locked sweep (44
configs); see `STATE.md` and `UNDERSTANDING_temporal.md`.

Eight area pairs (the list hardcoded in
`HC_V1_Code/HC_V1_temporal.m` and the legacy spatial-CCA scripts):

  CA1-V1, CA1-DG, CA1-CA3, CA1-RSC, CA1-SUB, V1-RSC, RSC-SUB, CA3-DG

200 bins x 2.5 cm (spatial arm). Fast-spiking units excluded (Tom's
convention: V1/RSC/CA1/CA3 only). After commit pass 1 (z-scoring fixed on,
min_units=5, Tom's LP, learners only), the spatial sweep is **66 configs**
(`{res,sig}` x 11 PC-count rules x 3 lag windows); the landmark arm is a
further **44 configs**. See `STATE.md`.

## Data layout

The pipeline reads:

    HC_V1_data/
      TF<id>_export.mat       # per-animal -- units.{idx,regions_label,idx_fs}
                              #               + analysis_spatial.firing.cued.freq
      animal_behaviour.mat    # cohort -- period_experienced(:, 1) + animal_id

See `UNDERSTANDING.md` for the full data contract.

## Quick start

The package is a self-contained src-layout Python package. Tests use pytest;
the package itself depends only on numpy, scipy, h5py, matplotlib, pandas and
openpyxl (for the sweep summary xlsx).

    cd cca/
    PYTHONPATH=src python -m pytest tests/ -q        # 504 tests, ~45 s

### Temporal arm (PRIMARY) — the method that carries the write-up

Every script imports `scripts/_common.py` first, which puts `src/` on the path, so
`python scripts/<x>.py` works without `PYTHONPATH`. Defaults are the canonical
temporal config (`config.TEMPORAL`: 10 ms bins, σ 2.5 ms, K 30, 5 folds, 200 shuffles,
±25 lags); `--include-fs` writes the co-primary FS-included file.

    python scripts/run_lag_curves.py      [--include-fs]   # refit-per-lag held-out curves (uncapped; ~3 h / FS; or scripts/run_lag_curves_uncapped_batch.sh)
    python scripts/run_cc_label_track.py  [--include-fs]   # frozen axes, FF/FB labels, per-epoch curves
    python scripts/run_fixed_subspace.py  [--include-fs]   # frozen epoch-balanced subspace (historical cap: 600 bins/trial)
    python scripts/run_lag_subspaces.py   [--epochs] [--include-fs]   # refit-per-lag subspaces (historical cap 12000/600)
    python scripts/run_lag_cosine.py      [--include-fs]   # CC self-cosine across lags   (historical cap 12000/600)
    python scripts/analyze_*.py                            # CSV -> stats CSV + md tables (both FS in one run)
    python scripts/figs_*.py [fsincl]                      # -> cca/figures/HCV1_*.{svg,png} (+ vault mirror)

`HANDOFF.md` maps script → CSV → figure and lists the column dictionary; the method
layer is `src/tom_cca/{preprocess,lagpairs,core,lagged,lag_subspace,fixed_subspace,
perdim_ifi,cc_aggregate,paired_stats}.py` (see `src/tom_cca/__init__.py`).

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
        sweep.py            # Tom sweep grid (66 configs post commit-pass-1)
      scripts/              # see UNDERSTANDING.md sec. 7
      tests/                # pytest suite
      results/, figures/    # outputs (gitignored or empty initially)
