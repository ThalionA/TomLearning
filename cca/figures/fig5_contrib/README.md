# Fig-5 contribution thread — final figures

Copied by `scripts/collect_fig5_figures.py` from the latest renders in `cca/figures/` (regenerate with the `figs_contrib_*.py` / `figs_trial12*.py` scripts, then re-run the collector). `_fsexcl_` and `_fsincl_` are co-primary; `.svg` + `.png` pairs.

- **HCV1_contribrel_forest** — Main result: unit contribution (contrib_conn) vs ±2-trial spatial reliability; rate-partialled survivors CA1-RSC RSC, CA1-V1 V1, V1-RSC RSC.
- **HCV1_contribrel_controls** — Controls: partner-invariant contribution shows no link; rate predicts reliability at rho ~ +0.4-0.5 (the confound the partial removes).
- **HCV1_contribrel_scatter** — Per-unit scatters, median animal of each robust cell.
- **HCV1_contribtomrel_forest** — Replication on Tom's precomputed moving-window reliability (z vs his shuffle null): same three cells plus CA1-DG DG (both FS).
- **HCV1_contribtune_forest** — Dissociation: contribution vs tuning (z vs shuffle null) — raw links broad, NO rate-partialled survivor in both FS.
- **HCV1_contribpool_forest** — Pooled over ALL partners: reliability link survives rate-partialling for RSC and V1 only.
- **HCV1_contribpool_tune_forest** — Pooled vs tuning: rate-carried; only RSC survives weakly.
- **HCV1_contriboverlap** — Overlap: same units serve several subspaces (rho +0.44-0.71) below the intrinsic-geometry ceiling; partner-specific residual positive; top-quartile Jaccard above chance.
- **HCV1_contriboverlap_ca1matrix** — CA1 partner x partner similarity: residual overlap strongest among hippocampal partners, weakest toward cortex.
- **HCV1_contribrel_epochs** — Per-pair reliability link across epochs (descriptive trajectories).
- **HCV1_contribepochs_pooled_traj** — Pooled links across naive (first 10 trials) / intermediate / expert.
- **HCV1_contribepochs_delta_pooled** — Epoch contrasts, pooled rate-partialled: stable (below-chance stars).
- **HCV1_contribepochs_delta_pooled_raw** — Epoch contrasts, pooled RAW: equally stable.
- **HCV1_contribepochs_delta_overlap** — Epoch contrasts, overlap metrics: stable.
- **HCV1_contribepochs_delta_pairs** — Epoch contrasts, per-pair rate-partialled reliability link: stable.
- **HCV1_contribepochs_delta_pairs_raw** — Epoch contrasts, per-pair RAW: stable.
- **HCV1_trial12_deltas** — Trial 1 vs 2 (frozen subspace): strength and direction null, behaviour (bins, speed) robustly different.
- **HCV1_trial12_control** — The 1->2 step vs the ordinal 3..10 adjacent-step band; V1-RSC is the (uncorrected) candidate.
- **HCV1_trial12_units_deltas** — Per-unit carrying, trial 1 vs 2: strength delta, profile stability vs adjacent band, convergence to the trained membership.
- **HCV1_trial12_v1rsc** — V1-RSC deep-dive: trial 1's whole CC1 lag curve is depressed (level effect, all lags); direction delta FS/arm-fragile; all uncorrected at n=9.
- **HCV1_connlearning** — Naive-epoch communication vs learning point: no survivor; CA1-DG strength is the only repeatable lean (rho -0.69, p=0.058 both FS, leverage-sensitive).
