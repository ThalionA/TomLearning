# NOTES — temporal CCA sanity audit + simulation iteration

Working file for the priority-ordered fix list across
`cca_animal_sanity_check.m` and `cca_simulation.m`.

Status legend: `[ ]` pending · `[~]` in progress · `[x]` done.

## P0 — correctness (without these, plots are uninterpretable)

- [x] A1. Plot real **and** shuffle on every CC1 figure in the sanity script.
       DONE 2026-05-03. Added shuffle to per-trial method comparison plot,
       trial-averaged spatial profile, plus a new pooled-distribution
       figure (07).
- [x] A2. IFI: signed-lag definition + symmetric `abs(neg+pos) > eps` guard
       (both temporal and spatial methods). DONE 2026-05-03.
- [x] A3. Held-out CC1 per trial: pool all valid windows' samples, fit on a
       random 70%, evaluate canonical correlation on the held-out 30%.
       One honest scalar per trial. DONE 2026-05-03. New helper
       `compute_heldout_cc1` + figure 08 plus the shuffle-null variant.

## P1 — high information gain

- [x] A6. Real-vs-shuffle CC1 distribution overlay (histograms / violin).
       DONE 2026-05-03 as figure 07.
- [x] A7. Window-size sweep: CC1 (real) and CC1 (shuffle) vs
       `half_window_bins ∈ {5, 10, 15, 20}`.
       DONE 2026-05-03 as figure 09_window_sweep, with errorbars.
- [x] A4. Contiguous-valid-block-length histogram across all trials.
       DONE 2026-05-03 as figure 00_block_lengths (bins + ms axes).
- [x] C3. Synthetic-control test at top of sanity script: known coupling
       CC1=0.5, recovered before any real-data processing.
       DONE 2026-05-03. New helper `synthetic_cca_check` covers both the
       coupled and the independent-null cases; aborts the script if either
       fails. Toggle via `run_synthetic_control` flag.
- [x] S1. Kernel CCA included in every sim sweep (not only sweep 2).
       DONE 2026-05-03. New `top_cc(X, Y, kappa)` helper returns top-1
       linear and RBF-KCCA in one call; sweeps 1, 3, 4, 5, 6 now plot
       both, with shuffle nulls where applicable. Sigma via median
       heuristic on stacked (X, Y).
- [x] S2. Sim sample-size sweep with actual generative model + shuffle
       baseline (companion to the existing pure-noise plot).
       DONE 2026-05-03 as Sweep 6b (sim_6b_samplesize_generative.png).

## P2 — correctness/robustness

- [x] A5. Per-region unit count, mean firing rate (post-rebin), fraction of
       zero bins per unit. DONE 2026-05-03 (figure 00b_region_stats).
- [ ] S3. Multi-seed averaging with mean ± SEM in sim sweeps.
       DEFERRED — invasive restructure of every sweep. Surface decision.
- [x] S4. Asymptotic CC1 reference line in sim sweeps. DONE 2026-05-03
       (sweep 6b: large-N empirical estimate as horizontal yline).
- [x] S5. Disentangle smoothing from sample-count in sim sweep 3.
       DONE 2026-05-03 by sampling a fixed n_fixed=500 bins regardless of w.
- [x] S6. Sim sweep 4 lag axis in real units. DONE 2026-05-03 alongside S1.
- [x] S7. Nonlinearity sweep covers identity / tanh / squared / sin+cos.
       DONE 2026-05-03 with shuffle null bars per family.

## P3 — hygiene

- [x] A8. Name and document the hardcoded "first 5 s cropped" parameter.
       DONE — `trial_warmup_ms = 5000`.
- [x] A9. Replace `lsline` with explicit `polyfit` + NaN filtering.
       DONE — `robust_fit_line` helper, six call sites swapped.
- [x] A10/C2. Save numerical data alongside figures. DONE for both scripts
       (`sanity_data.mat`, `sim_data.mat`).
- [x] A11. Progress prints inside the sanity script trial loop. DONE
       (`fprintf` every 10 trials).
- [x] A12. Remove redundant double-centering in `get_pca_data`. DONE.
- [x] C1. Factor out hardcoded `/Users/theoamvr/...` paths. DONE — both
       scripts now resolve via `mfilename('fullpath')`.
- [x] C4. `rng` seeds for reproducibility. DONE — `rng(42)` at the top of
       both scripts; the synthetic-control helpers seed locally.
- [x] S8. Kernel CCA docstring + validation against a known-coupling case.
       DONE — `check_top_cc` runs three constructed cases (identical /
       nonlinear / independent) at the top of the sim script and aborts
       on failure.

## Deferred (need user input)

- **S3 multi-seed averaging.** Each sweep currently runs a single
  realization. Robust ± SEM lines would require wrapping every sweep in
  `for seed = 1:K` and aggregating. Easy in principle, ~20–50× longer
  runtime. Worth doing once the qualitative shape of each sweep is
  agreed.
- **S3 implication: most sweep curves should not be over-interpreted
  yet.** Single-seed CC1 wobbles are large at small sample sizes.

## Conventions

- Per CLAUDE.md: every plot needs axis labels with units, title, legend.
- Save the numerical data behind every plot as `<name>_data.mat` so the
  figure can be regenerated without rerunning the analysis.
- New plots go to `HC_V1_figures/CCA_Animal5/` (sanity) or
  `HC_V1_figures/CCA_Simulations/` (sim).
- Set `rng(42)` before any shuffle / synthetic generation.
