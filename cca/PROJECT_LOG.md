# PROJECT LOG — Tom-learning CCA

**Purpose.** Durable, cross-session memory. Read this first (with `STATE.md`) at the start of any
session; append a dated entry after any meaningful work. Newest entry on top. Keep entries terse
but self-contained — a future session must be able to resume from here alone.

**Doc map.** **`HANDOFF.md` = how the data/results/figures are organised — read it first if
you are picking up the temporal or lag analyses.** `STATE.md` = current findings + canonical configs (the verdict). `OPPORTUNITIES.md` =
the two-paper reframe + plan. `GOTCHAS.md` = bugs not to reintroduce. `NOTES.md` = older dev log.
`UNDERSTANDING.md` / `UNDERSTANDING_temporal.md` = original specs. **This file = the running
narrative + state of play.**

---

## CURRENT STATE (2026-08-30) — contributing units ARE spatially special (cortical side,
## rate-partialled, FS-robust); trial-1-vs-2 answered 08-20; ➜ READ `HANDOFF.md` FIRST

**Nothing is running.** Newest work is the 2026-08-20 entry directly below: `run_trial12.py` /
`analyze_trial12.py` compare trials 1 and 2 through a frozen subspace fit on ordinals 11+.
Headline: **0 BH survivors in 240 science tests per FS condition**, but the CIs are too wide to
call it a null for *direction*; the strength null IS informative for CA1-CA3 / CA3-DG / CA1-RSC /
V1-RSC. What is robust is that **trial 1 is longer, slower and more fragmented than trial 2**,
and bin-matching fixes only the first of those three. One candidate kept: **V1-RSC**, whose 1→2
strength step is the largest of 9 adjacent-trial steps in both FS conditions and is not speed-driven
(does not survive BH; n = 9). Three traps recorded in GOTCHAS, including that the `sig` gate passes
~100 % of dims, so "all significant CCs" means "all CCs". The 2026-08-17 entry below is unchanged
and still current for the lag-curve arm.

**Since (2026-08-28):** (a) **FINDING MOVES — the participation-broadening headline was
re-tested on the connection-specific Ginis and only half of it survives: CA1-DG is DEAD, CA1-RSC
survives but on the *performance* axis, and a new FS-robust CA1-V1 effect appears on the V1 side.**
See the entry directly below; `STATE.md` §3.0 finding 2 is rewritten. (b) Method layer only:
`subspace_window` now exports the CC1 lag curve (`lags`, `lag_cc1`) — back-ported from the
StriatumACC `striatum_tcca` port, which also took `core.cca_fit`'s covariance route, `lagpairs`,
the per-dim held-out lag curve + `perdim_significance`, and `paired_stats.paired_t` / `welch_t`
from here.

---

## (2026-08-30) — figures for the contribution × reliability finding (and a standing rule)

**Nothing is running.** `scripts/figs_contrib_reliability.py` (new) renders the 2026-08-29
finding, both FS: `HCV1_contribrel_forest_*` (16-cell forest, raw vs rate-partialled, Wilcoxon
stars), `_controls_*` (intrinsic-contrib control + rate–reliability confound), `_scatter_*`
(per-unit scatters, median animal per robust cell), `_epochs_*` (rate-partialled ρ across
epochs). One placement bug caught by looking at the render: `figstyle.star`'s 5-pt upward offset
lands between rows on an inverted categorical axis — stars are drawn va=center directly there.

**Standing rule (Theo, exasperated repeat — logged in ~/.claude/MISTAKES.md as `re-explained`):**
an analysis in this project is DONE only when its `figs_*.py` exists, the render has been looked
at, and the figure ships with a legend + method summary. Tables in `results/*_tables.md` are the
appendix, never the deliverable.

---

## (2026-08-29) — reliability × subspace membership: contributing units are spatially special on the CORTICAL side

**Nothing is running.** Tom's Fig. 5 second half ("link neuronal contribution to the subspace
with spatial representational changes") — the reliability leg, spec'd by Theo as: reliability =
mean trial-to-trial Pearson correlation of the spatially-binned map within **±2 trials**, edges
clipped.

**New (TDD, 11 tests → 546 total):** `src/tom_cca/spatial_reliability.py` —
`trial_map_reliability` (pairwise-complete over finite bins; silent trials NaN and excluded, not
zero; edges clipped) + `epoch_mean_reliability` (temporal trial ids are **1-based** — measured on
TF073, ids 1..214 for 214 spatial rows — the function subtracts 1 and the test pins it).
Computed from `Animal.spatial_fr` (`freq`) directly, NOT the export's precomputed
`analysis_spatial/reliability` (provenance unverifiable).

**Regenerated:** `run_epochs.py --bin-ms 10 --smooth-ms 2.5` both FS (~50 min each, 2 procs) —
the June epoch CSVs predated `contrib_conn` (per-neuron, 2026-07-28) and the 2026-08-19
preprocess refactor. June versions archived to `results/_archive/*_2026-06-17.csv`. New:
148/153 metric rows, 12 817/14 445 weight rows (June: 12 818/14 446 — one CA1-CA3 edge row
fewer; not chased).

**New driver:** `scripts/analyze_contrib_reliability.py` — joins per-neuron `contrib_conn` to
epoch-matched reliability (same trials the epoch's CCA was fit on; join via
`dataio.select_units` position→raw-index, FS condition enforced to match the CSV). Per
(animal, pair, epoch, area) Spearman across units; Fisher-z over epochs within animal;
animals-as-n t + Wilcoxon per (pair, area). Controls: (1) log-rate rank-partialled out of both
(rate–reliability rho ≈ +0.4–0.5 in EVERY area — half the raw link is rate); (2) the
partner-invariant `contrib` as specificity control.

**Finding (STATE §3.0 finding 7b; tables `results/contrib_reliability_tables.md`):** raw
`contrib_conn`–reliability correlation is positive nearly everywhere (strongest on the
second-named/cortical side). After rate-partialling, FS-robust: **CA1-RSC RSC side**
(+0.256/+0.255, W p = 0.016 both FS), **CA1-V1 V1 side** (+0.258/+0.170, W p = 0.002/0.037),
**V1-RSC RSC side** (+0.312/+0.411, W p = 0.031 both; V1 side 6/6 positive, t starred both).
The CA1-side link never survives the rate partial; the intrinsic-`contrib` control is n.s.
nearly everywhere → the association is communication-specific. n = 4 pairs are on the Wilcoxon
floor (descriptive only). Epoch-resolved correlations are roughly stable; CA1-CA3 CA3-side
grows naive→expert (−0.03 → +0.42 FS-excl, descriptive).

**Caveat for any write-up:** reliability and the CCA input are the same spikes; log-rate
partialling handles the SNR confound only as far as log-rate proxies SNR. A weight-shuffle null
(reliability of top-contribution units vs contribution-shuffled draws) is the natural
strengthening if this goes in a figure.

**Next:** tuning-score and COM-shift legs of the same join (loader for
`analysis_spatial/tuning_score/score` + shuffles still unwritten); figure panel (scatter +
per-pair rho forest) if Tom wants it in Fig. 5.

---

## (2026-08-28) — participation-broadening re-tested on the connection-specific Ginis: CA1-DG dies, CA1-RSC survives on *performance*, CA1-V1 appears

**Nothing is running.** This closes the question left open on 2026-07-28 ("whether
participation-broadening survives `gini_*_conn` / `gini_*_sig` is OPEN"). The trajectory CSVs
that carry those columns had been complete since 2026-07-28 (both FS, 16 animals, 1927 / 1958
window-rows) — `analyze_trajectory.py` simply never had the column names in its metric lists, so
the re-test had never been run. **Code change: four metric names added to `SLOPE_METRICS`,
`LEVEL_METRICS` and `param_metrics`.** No refits; same CSV, same rows, same windows for all five
metrics, so this is a **like-for-like** comparison of definitions, not a new run.

Full tables: **`results/gini_conn_retest_tables.md`** (LMM p per pair × axis × metric × FS).

**1. CA1-DG is dead.** `gini_x` trial_frac LMM p = 0.0088 (FS-excl) / 0.0281 (FS-incl) →
`gini_x_conn` p = 0.272 / 0.203, `gini_x_sig` p = 0.202 / 0.063. Same on `performance`
(0.0039 / 0.024 → 0.293 / 0.200). Both sides, both FS, both connection-specific definitions:
n.s. **Drop CA1-DG from the broadening claim.**

**2. CA1-RSC survives — but the axis moves from time to performance.** On `trial_frac`
(time-on-task) the effect that was p = 1.4e-05 under the area-intrinsic metric drops to
p = 0.081 / 0.223 under `_conn`; only `_sig` survives, weakly (0.019 / 0.038). On `performance`
it is robust under **every** definition, **both sides of the connection**, **both FS**:
`gini_x_conn` 6.3e-04 / 1.0e-04, `gini_x_sig` 1.9e-03 / 1.9e-04, `gini_y_conn` 0.022 / 1.3e-03,
`gini_y_sig` 1.4e-03 / 2.7e-04. This inverts the old reading: the *area-intrinsic* metric looked
like a clock (trial_frac), the *connection-specific* one looks like a performance covariate.

**3. New and previously invisible: CA1-V1, on the V1 side.** `gini_y_conn` trial_frac LMM
0.0269 (FS-excl) / 0.0056 (FS-incl); performance 0.0293 / 0.0152. `gini_x` (CA1 side) is n.s.
throughout. The old battery only ever tested `_x`, so a V1-side effect could not have been seen.
Relevant to Fig. 5, where V1 is the focus area.

**4. Learning attribution is STILL not established — the axis move does not rescue it.**
Non-learners de-sparsify on `performance` with a *larger* median slope than learners
(CA1-RSC `gini_x_conn`: learners med −0.106, non-learners med −0.287, both FS). n = 3–4 so the
p-values are at the Wilcoxon floor and uninterpretable, and there is an unresolved confound: a
non-learner's `performance` range is compressed, which inflates a slope taken w.r.t. it. Read
this as "not evidence for learning-specificity", not as "evidence against".

**5. Two traps in this table (→ GOTCHAS).** (a) `gini_*_sig` is NaN whenever `n_sig` = 0
(56 % of CA1-RSC windows), so its slopes are conditioned on windows that had a significant dim
and its animal count can drop (CA1-RSC FS-excl 8 → 7). Never read `_sig` without the `n_sig`
slope beside it. (b) The n = 4 pairs (CA1-SUB, CA3-DG, RSC-SUB) cannot star at the honest unit —
the two-sided Wilcoxon floor at n = 4 is p = 0.125 — yet their LMM p-values go to 6e-06 by
pooling windows within those 4 animals. CA1-SUB `gini_x_conn` p = 7.6e-04 is that artefact, not
a finding.

**Next:** Fig. 5's "contributions to the communication subspace with experience" bullet is now
answerable — CA1-RSC (both sides, performance axis) and CA1-V1 (V1 side) — with the
learning-vs-experience caveat attached. The Tom-facing gap that remains is the second half of
Fig. 5: `dataio` still does not load `analysis_spatial/tuning_score/score`,
`analysis_spatial/reliability/units/reliability_moving_window`, or
`analysis_spatial/reliability/corr/trial`, which are what "reliability / tuned population /
shifting responses" need. (⚠ `OPPORTUNITIES.md` §2 lists these under `units/` — wrong path,
verified against `TF073_export.mat`; they live under `analysis_spatial/`.)

---

## (2026-08-28) — two-way sync with the StriatumACC temporal port; `subspace_window` exports the lag curve

> **CORRECTION (same day, appended).** This entry's claim that the striatum-side A/B "confirms the
> covariance-route `cca_fit` on a second dataset (3.1e-13, 2.7x)" rests on ONE cell and is too
> strong. Re-measured on the whole striatum epoch grid, holding cell types fixed: **`cca_fit`'s
> route makes no difference under plain CCA (0/153 cells) but moves 7/152 cells under partial CCA**,
> worst cc1 0.40. Partialling out other areas destroys rank, and the directions the covariance route
> cuts sit at relative singular values 1e-9 to 4e-9 — the SVD route was fitting canonical dimensions
> on double-precision debris (a14's ACC: 21 units, Z=13, only 14 real directions, SVD fitted 17).
> The covariance route is still the right rule; the equivalence claim is not "identical numbers", it
> is **"identical unless you partial"**. **Relevant here:** `run_lag_subspaces.py` is the only
> TomLearning driver that partials out a third area, and `d3ffec4`'s equivalence check was run on
> the *unpartialled* lag-curve path — so the 2026-08-17 "agree to 5e-15" result does not cover it.
> Worth re-checking `lag_subspaces` before its numbers are cited.

**Nothing is running. No finding changes — this is the method layer only.**

`StriatumACC/Striatum project/tcca` (`striatum_tcca`) is a port of this package frozen on
2026-07-28; 37 commits of drift since. A module-by-module diff put 8 of its 19 numeric modules
byte-identical to ours, and the shared numerics agree to ~3e-15 (both packages imported side by
side: `cca_fit` 2.8e-15, held-out CC1 lag curve 3.3e-15, lag pairing bit-identical). Four things
were resolved.

**Taken from here into `striatum_tcca`** (their `NOTES.md` 2026-08-28 has the detail): the
covariance-route `cca_fit` (`d3ffec4`), `lagpairs.py`, `ifi_sides` +
`heldout_lag_curve_flat_perdim` + `perdim_significance`, and `paired_stats.paired_t` / `welch_t`.
Their A/B on real striatum data (animal 1, DMS-DLS, 3 epochs, full `WindowSubspace`) puts the
covariance route within **3.1e-13** of the frozen SVD route and **2.7× faster** — an independent
confirmation of our 2026-08-17 `_cca_fit_svd` pinning on a second dataset.

**Taken from `striatum_tcca` into here — the only code change in this repo:**
`subspace_window.WindowSubspace` gains two fields, `lags` and `lag_cc1` (the in-sample CC1 lag
curve and its lag axis, already computed inside `window_subspace` and previously discarded).
Drivers can now export the curve per cell and recompute IFI at any integration window offline via
`lagged.ifi_by_window`, instead of being locked to the `max_lag` the run used. Purely additive:
one construction site, all-keyword, no existing field or number touched. Their driver test came
with it (`test_lag_curve_exported_with_expected_shape`).

Also added here: 9 tests for `paired_t` / `welch_t`, which had shipped untested since 2026-08-19
(`tests/test_paired_stats.py`, identical file in both repos). **535 tests.**

**Deliberately NOT unified.** `config.py` / `dataio.py` / drivers are the dataset boundary and
stay separate. `membership.subspace_contribution_connection` and the `gini_*_conn` / `gini_*_sig`
fields were not ported — worth knowing that the striatum project's "corrected Gini" is
`gini_pearson_x/y` (the CCA-independent Pearson control, which we share), so our two
*connection-specific* corrected definitions have never been run on that dataset. A shared
installable core package is the obvious end state but should wait for the audit's P1/P3 passes
(`audit/REPORT_2026-08-19.md` §5), which will reshape exactly the modules it would hold.

---

## (2026-08-17) — lag curves re-run WHOLE-SESSION; asks 1/3 + item 1 re-answered

**Nothing is running.** Two things happened on 2026-08-17 on top of the 08-15 entry below:

### 1. `run_lag_curves` was silently answering session questions on the first ~20 trials — fixed

Theo (rightly) objected that the IFI-windows figures used the 12k-bin cap = the **first ≤600
bins of the first ~20 trials**. That cap was the driver's default since 2026-07-28; HANDOFF §2
had flagged it as "wrong for epochs" and I had left it in place for the session-level asks —
a mistake. Now: **`run_lag_curves.py` default is UNCAPPED** (all ~370 k running bins/animal;
`--max-samples/--max-per-trial` keep the old cap for smoke tests, `--animals` for subsets).
Re-run both FS via `scripts/run_lag_curves_uncapped_batch.sh` (6 processes, 3 BLAS threads
each, ~3 h; `merge_lag_curve_parts.py` stitches) — `results/lag_curves_bin10{,_fsincl}.csv`
are now whole-session (76 347 / 83 283 rows, 71 / 73 animal-pairs, 16 animals).
**Every pre-2026-08-17 number derived from `lag_curves_bin10*` (item 1, asks 1/3,
`figs_lag_curves`) is superseded.**

**Made tractable by a core change (`d3ffec4`): `core.cca_fit` now uses the k×k covariance
route** (eigen-whitening) instead of a thin SVD of the n×k data matrix — 97 ms vs 1.66 s per
fit at 300 k rows, ×1255 fits per pair (curve + 200-shuffle null). Pinned to the old
implementation (`core._cca_fit_svd`, kept for tests only) in 4 tests to 1e-12; on a real
animal-pair through the full lag machinery the held-out curves agree to 5e-15 with identical
p-values and FDR masks; 410 tests pass. Rank rule documented (`_COV_EIG_FLOOR = 1e-10`
relative eigenvalue ≈ 1e-5 relative singular value; the Gram matrix cannot resolve below that
and PCA scores never get near it). Every arm uses `cca_fit`, so all results reproduce to
floating precision — but note it in any methods text.

### 2. Whole-session answers (both FS; `cc_ifi_signs_tables.md`, `HCV1_cc_ifi_windows_*`)

**Ask 3 — overall IFI across all significant CCs at ±50 ms vs 0** (656 / 701 significant
dims, was 238 / 307; every animal now contributes; magnitudes 3–5× smaller than the capped
values, which were noise-inflated):

| pair | FS-excl (animals) | FS-incl (animals) | CCs-as-n | BH across pairs |
|---|---|---|---|---|
| **CA1→RSC** | +0.053 p=0.005 (n=12) | +0.068 p=8×10⁻⁵ | 6×10⁻⁶ / 3×10⁻¹² | **survives, every look** |
| **DG→CA1** | −0.053 p=0.014 (n=11) | −0.044 p=0.020 | 3×10⁻⁵ / 2×10⁻⁴ | weighted + CCs |
| **DG→CA3** | −0.061 p=0.027 (n=5) | −0.082 p=0.013 | 0.002 / 2×10⁻⁵ | weighted (excl) + CCs |
| **V1→RSC** | +0.038 p=0.024 (n=9) | +0.030 p=0.050 | 6×10⁻⁷ / 9×10⁻⁴ | weighted + CCs |
| CA1-CA3, CA1-V1, CA1-SUB, RSC-SUB | n.s. | n.s. | n.s. | — |

CC₁ on the same data now reproduces `bin10_tables.md` §B to 3 dp (CA1-RSC +0.081, CA1-SUB
−0.087, V1-RSC +0.021) — the earlier r ≈ 0.4 disagreement was entirely the cap. Reading:
CA1→RSC and V1→RSC feed-forward into RSC; DG leads both CA1 and CA3 (the anatomical
direction); CA1-SUB's CC₁-only SUB→CA1 flow (§B) is **not** carried by the all-CC average
(−0.025 p=0.40 / +0.008 p=0.87).

**Item 1 re-read on whole-session data:** sign-mixing is now **LESS than chance** — 64/71
animal-pairs mixed vs 70.2 expected (99 %) at ~9 sig CCs each, binomial p < 0.001 (FS-incl
69/73 vs 72.6, p = 0.001). An animal's significant CCs share a direction more often than
independent signs would — the opposite of "different CCs, different signs" and consistent
with ask 3. CA1-RSC CC₁ (+0.081, p = 6×10⁻⁴) survives per-pair FDR (CA1-CA3 CC11 FS-incl too);
no pair has two reliably-signed dims of opposite sign. The two verdict strings in the md
were hard-coded and are now data-driven — a recurring hazard, see GOTCHAS.

**Ask 2 also gained CCs-as-n** (Theo's ask): `cc_crosscorr_epochs_stats_*` and the figure
titles carry both units. IFI naive→expert stays null on both; the baseline offset is
p < 10⁻⁵ on CCs; peak−baseline gains ~5 % nominal "hits" on CCs (Δ ≈ +0.004) that are 4–6×
smaller than the offset. Also verified by hand: CA1-V1 FS-incl has a LARGER p on CCs-as-n
(0.45 vs 0.12) because two animals contribute one very negative CC each while animal 61's
seven CCs span −0.39…+0.73 — CCs-as-n is only more powerful when the CCs agree.

**Tests: 410 pass.** Figures regenerated: `HCV1_cc_ifi_windows_*` (whole-session, both
units in titles), `HCV1_perdim_ifi_*`, `HCV1_lagcc_*`, `HCV1_cc_crosscorr_epochs_*`.

**Open:** MEETING/STATE/HANDOFF numbers for item 1 and §3.0 finding 3 updated below; the
per-epoch speed export for the offset finding is still the next thing.
**2026-08-19: code audit (`audit/REPORT_2026-08-19.md`) → the DEDUP passes are done and verified like-for-like (4 commits, 504 tests; entry below). Drivers now share `preprocess.py`/`_common.py`/`config.TEMPORAL`. Open: archive decision, P4 semantic fixes, METHOD.md.**

---

## ✓ DONE (2026-08-20) — TRIAL 1 vs TRIAL 2 through the temporal method: no communication
## difference survives correction; the difference that IS robust is behavioural

**Nothing is running.** Theo's ask: apply the temporal CCA to trials 1 and 2 only and compare
them directly — strength, direction (IFI at CC₁, over all significant CCs, and across
integration windows), both FS. Gini was dropped by decision (the weight-Gini needs a ≥5-trial
refit; the per-trial participation proxy was declined).

### Method (`scripts/run_trial12.py` → `analyze_trial12.py`, `commit 249edab`)

A single trial (~500–2900 running bins) cannot refit a 30-dim CCA, so the subspace is fit ONCE
per (animal, pair) on running-trial ordinals **11+** (`early_trials.reference_fit`) and
ordinals 1..10 are PROJECTED through it. Trials 1–2 are leak-free **and** ordinals 3..10 are
equally out-of-fit, which is what makes the adjacent-trial control like-for-like. Per-trial lag
curves come from `fixed_subspace.variate_lag_curve`; the driver exports curves only and every
IFI is derived downstream by `perdim_ifi.windows_table`, so the degenerate-zero convention
cannot leak into a mean. Three arms (`matched` column): 0 raw, 1 ordinals 1/2 cut to
min(n₁,n₂) leading bins, 2 all ordinals cut to a common min.

**Cohort: 16 animals, 71 (FS-excl) / 73 (FS-incl) animal-pairs — the canonical temporal-arm
cohort.** Verification before any reading: per-trial `r0` reproduces `early_trials_projected`'s
independently-computed `cc1` across **all 720 (animal, pair, ordinal) cells to 5×10⁻⁵**
(the rounding floor), bin counts identical; no fit or significance failures. Runtime ~25 min/FS.

### Result — a failure to detect, NOT a demonstrated null

**0 BH-surviving tests of 240 science tests in each FS condition**, and **0 BH survivors** in
the window sweep under its proper nested-window test. But this must not be written up as "no
change": the 95 % CIs are far too wide to exclude a project-scale effect. Median CI half-width
on the headline IFI delta is **0.28 (FS-excl) / 0.15 (FS-incl) IFI units against a largest-ever
session-level effect of 0.061** — i.e. the direction result excludes almost nothing.
Power (MDES, 80 %, paired *t* at the real n): the **strength** null is genuinely informative for
**CA1-CA3, CA3-DG, CA1-RSC, V1-RSC** (a first-trial deficit of 8–28 % of the pair's own coupling
would have been caught) and marginal for RSC-SUB; the **direction** null is uninformative
everywhere except V1-RSC. 4/8 pairs have n ≤ 7 and Wilcoxon cannot reach p < 0.05 at CA3-DG's n = 5.

### The robust difference is behavioural, and matching only half-fixes it

| trial 1 − trial 2 | raw | matched=1 | matched=2 |
|---|---|---|---|
| `n_bins` | **+650** (7/8 pairs BH) | **0.000** (0/8) | **0.000** (0/8) |
| `vel_mean` | **−3.30 cm/s** (8/8 BH) | −3.51 (8/8) | −3.04 (7/8) |
| `n_gaps` | **+6.8** (8/8 BH) | +6.2 (8/8) | +5.2 (7/8) |

Trial 1 is longer (+6.5 s, p = 5.5×10⁻⁶), slower and more stop-start. Bin-matching removes the
duration difference **exactly** and the other two **not at all** — it keeps LEADING bins, so it
compares ~the first 80 % of trial 1 against all of trial 2, over-weighting the trial-onset phase
that HANDOFF §6 already flags for a doubled speed confound. **No arm of this analysis is
behaviour-free.** → GOTCHAS.

### The one candidate worth keeping: V1-RSC

Held-out CC₁ strength is lower in trial 1, **negative in all 6 (FS × arm) cells**, and the
per-ordinal profile is a **step, not noise**: mean r₀ by ordinal (arm 2, FS-excl)
`0.128 0.181 0.183 0.198 0.192 0.188 0.176 0.206 0.224 0.186`. The 1→2 step is **the largest of
the 9 adjacent-ordinal steps in BOTH FS conditions (rank 1/9), 3.4× / 3.8× the mean of the other
eight** (−0.053 p = 0.066 FS-excl; −0.059 p = 0.035 FS-incl, n = 9). **Not a speed artefact**:
corr(Δr₀, Δspeed) = +0.19 / +0.16 (p = 0.63 / 0.68). It does **not** survive BH and n = 9 —
a candidate for a targeted test, not a finding. CA1-DG `ifi_cc1` is the FS-included analogue
(p_perm = 0.008–0.035 at 60–100 ms windows) but is absent FS-excluded, so it is FS-fragile.

### ⚠ Three things that invalidate naive readings of these tables

1. **The `sig` gate is VACUOUS**: 100.00 % (FS-incl) / 99.56 % (FS-excl) of dims pass, ~99 % at
   the p-floor 1/201. "Average over all *significant* CCs" IS "average over all 10 CCs"; the
   `ifi_cc1` vs `ifi_allsig` contrast is CC₁ vs a 10-dim mean, nothing more.
2. **The adjacent-trial control is not an independent falsification arm**: it correlates
   r = 0.85 with the primary statistic (it is the same delta rescaled per animal), and it tests
   mean(z) vs 0 — a *directional* test — while the interpretation ("|z| ≈ 1 is ordinary jitter")
   is a *magnitude* claim that is never actually tested. Its exchangeability assumption does hold
   (no drift across ordinals 3..10: mean within-animal r = 0.005, p = 0.92).
3. **The two co-primary FS conditions disagree on direction**: median cross-FS correlation of the
   per-pair delta is 0.83 overall but only **0.40–0.51 for `ifi_cc1`** (4–6 of 8 pairs even agree
   in sign). Strength agrees; direction does not.

### Adversarial verification found all of this — it was not visible in the first pass

A 7-agent workflow was run against the finished result **before** it was written up, and
**refuted all five headline claims**. It caught two number-changing bugs (`cc_peak` meaning two
different things while both fed the same weight; `drop_degenerate` dropping dims per-ordinal so
the delta compared different dim-supports — 85/213 headline cells, and it flipped the CA1-V1
sign), the untested window sweep, the raw-arm-only covariate tests, and a frozen-projection
degeneracy (animal 70, |r| ≈ 1 at every rank, now gated). My own interim reading of the control
table was **wrong in the opposite direction** and is logged in `~/.claude/MISTAKES.md`
(2026-08-20, S, over-claimed — read one column of one pair and called the control "decisive").

**Files:** `results/trial12_{trials,delta_tests,sweep,sweep_test,control,completeness,degenerate}_bin10{,_fsincl}.csv`,
`trial12_fs_agreement_bin10.csv`, `trial12_tables.md`. Curves CSVs (~50 MB/condition) are
gitignored and regenerable. **525 tests pass.**

**Next, in order:** (1) a targeted V1-RSC test — it is the only pair with both power and a
coherent signal; pre-register it, since it is currently a post-hoc pick out of 8; (2) a
speed-matched arm (match on velocity, not bin count) — no current arm is behaviour-free;
(3) decide whether `n_sig` is worth keeping at all given the vacuous gate; (4) `figs_trial12.py`
is NOT written (the only analysis in the repo with no figure).

---

## ✓ DONE (2026-08-19, later) — DEDUP of the temporal method, zero behaviour change (4 commits)

**Nothing is running. No result changed** — every pass was verified like-for-like before
commit. Theo: "start with the massive dedup"; `duplication` logged as S in `~/.claude/MISTAKES.md`
(new tag, user-flagged critical; rule candidate: grep before writing a helper in a script —
exists in src → import, exists elsewhere → hoist with a test first).

| commit | what | verification |
|---|---|---|
| `93028ae` src | `lagpairs.py` = THE within-trial lag pairer (was 9 copies; 5 live callers repointed: `lagged._segment_lagged_pairs`, `lag_subspace.segment_lag`, `fixed_subspace.variate_lag_curve/trial_lag_moments`); `core.pca_fit_flat/pca_project/pca_scores` (was 7 copies incl. 2 byte-identical driver copies); ONE `paired_stats.fdr_bh` (was 4; `lagged.fdr_bh`, `fixed_subspace._fdr_bh` are aliases); `lagged.ifi_sides` (kills 2 inline re-derivations of the clipped means) | `tests/test_dedup_equivalence.py`: 81 tests pin the shared functions to VERBATIM copies of every old implementation (1e-12); 414 old tests unchanged |
| `6462813` drivers | `src/tom_cca/preprocess.py`: `RunningSession` + `load_running_session` (the 13-copy prep block: stream load with a printed skip REASON, running mask, ≥5-unit areas, cap, `pair(ax, ay) → X, Y, Z`), `cap_running_bins` (4 cap copies → 1; reproduces all three old selection rules), `residual_pca_scores`, `epoch_of_trial`; `config.TEMPORAL` (10 ms · σ 2.5 · K 30 · 5 folds · 200 shuf · ±25 · fdr_dims 10 · label_w 5 — ONE object), `PAIRS_TEMPORAL`/`PAIR_NAMES`, `FIGURE_MIRROR_DIR` (env `CCA_FIG_MIRROR`), `FIG_PREFIX`; `dataio.load_learning_points()` public; `scripts/_common.py` (sys.path once, `RES/FIGS`, `PAIRS`, `FS_CONDITIONS`, `temporal_parser`, `results_path`); 6 drivers rewritten (`run_lag_curves/ifi_windows/fixed_subspace/cc_label_track/lag_cosine/lag_subspaces`, ~40 % shorter, every flag shared) | one-animal (36) smoke harness, old vs new: **7 CSVs byte-identical** (incl. `run_lag_subspaces --epochs`); 9 more equivalence tests (caps, session, epoch map) |
| `67683a4` analyze | 8 live `analyze_*` take `PAIRS/RES/EPOCHS/BIN_MS/LABEL_W/EXPECTED_LEARNERS` + the FS loop from `_common`/`config` | all 8 re-run → `git status results/` clean (every tracked output byte-identical) |
| `374c0dd` figs | `figstyle.save(fig, "HCV1_x")` → `config.FIGURES_DIR` primary + vault mirror (a path stem still works for legacy scripts); 11 live `figs_*` drop `ATT/RES/PAIRS/sys.path` | 48/48 live PNGs (both FS) byte-identical before/after; vault mirror byte-identical |

**Behaviour-preserving choices to know about:** `run_lag_subspaces`/`run_lag_cosine` keep
their HISTORICAL cap (12 000 / 600 → first ~20 trials) and `run_fixed_subspace` its 600-per-trial
cap as parser DEFAULTS — now printed at start and exposed as `--max-samples/--max-per-trial`
(`0 0` uncaps). Uncapping them is a separate, numbers-changing decision (HANDOFF §2).
`run_ifi_windows` default bin is now 10 (was 25; STATE §6 always passed 10 explicitly).

**Docs fixed:** HANDOFF §1 (defaults/shared prep/caps) and §2.4 (partial-out is train-only
ONLY in `run_lag_subspaces` — the old sentence overstated it); README quick-start (temporal
arm first, 504 tests); `src/tom_cca/__init__.py` rewritten as the method map (20 modules
exported, temporal first). Audit report status block updated (`audit/REPORT_2026-08-19.md`).

**Still open from the audit (Theo's picks):** §1 archive of the spatial/landmark/kCCA arms +
the original Arm-A/B pipeline to `legacy/` (biggest remaining simplification; `bin10_tables`
chain needs a call); 1.7 `subspace_window` → `perdim_significance` swap (like-for-like first);
2.8 one circshift null honouring `temporal_circshift_min_bins`; 2.9 route the 3 hand-rolled
per-animal aggregations through `cc_aggregate` (n may change by one animal — report table);
4.1 rename `cc_heldout` → `r_frozen`; 4.3 config sidecar per results CSV; `pyproject.toml`;
4.7 `METHOD.md`.

---

## 📋 AUDIT (2026-08-19) — temporal CCA method: consolidate / simplify / document — REPORT ONLY, nothing moved

**Nothing is running. No code or results changed.** Theo asked for an audit of the temporal
CCA method with the goal of a modular, robust codebase portable to another dataset.
Deliverable: **`audit/REPORT_2026-08-19.md`** (+ four `_appendix_*` files with the raw
file:line sweeps). Baseline measured: 31 src modules / 6 954 LOC, 83 scripts / 16 828 LOC,
**414 tests pass in 42 s**.

Headlines (details + file:line in the report):
- **~⅔ of the tree is off the temporal path** — the spatial / landmark / kCCA arms and the
  *original* Arm-A/B pipeline (`lagged_temporal`, `lagged_landmark`, `segments`,
  `pipeline.prepare_pair_temporal`, `analysis.analyse_pair_temporal`, `run_temporal.py`) are
  imported by **none** of the five live drivers. Candidate `legacy/` move (~3 000 src LOC, ~45 scripts).
- **Drift risk:** within-trial lag pairing implemented **9×** (guards differ; `subspace_window`
  crosses trials); BH-FDR **4×** (identical); PCA→scores **7×** (two driver copies byte-identical);
  sample cap **4×** with **two different selection rules** (`run_fixed_subspace._cap_bins` ≠ the
  12 k rule); running-bin prep block copied in 13 drivers; IFI 0.0-vs-NaN split between
  `lagged.information_flow_index` and `perdim_ifi.curve_ifi`; per-animal aggregation bypasses
  `cc_aggregate`'s guards in 3 analyze scripts.
- **Docs ≠ code:** HANDOFF §2.4 ("train-only partial-out where the readout is cross-validated")
  is true only for `run_lag_subspaces`; `run_lag_curves`/`run_lag_cosine` residualise + PCA in-sample
  (mild, unsupervised w.r.t. X–Y — doc fix first). `config.temporal_bin_ms=50 "Locked"` vs drivers at 10.
  `cfg.temporal_circshift_min_bins` is read by no live null (global roll, min shift 1 — **not** a
  result-changing bug at n≈3e5; config/doc drift). **`cc_label_track_*.csv` column `cc_heldout` is the
  IN-SAMPLE frozen-fit r** (`run_cc_label_track.py:255` writes `frozen_perm_null(...).r_obs`) — carried,
  not tested on; rename candidate.
- **Portability:** `src/` is clean of area-name literals outside `config.py`; the coupling lives in
  the scripts — `PAIRS` redeclared 51×, vault path hard-coded in 24 figs scripts (7 absolute
  `/Users/theoamvr/…`), `sys.path.insert` in ~80 scripts (no `pyproject`), the canonical temporal
  config (10 ms · σ 2.5 · k 30 · 5 folds · 200 shuf · ±250 ms · fdr_dims 10 · label_w 5) exists nowhere
  as one object.
- **Proposed shape** (report §5): eight temporal modules (`config, dataio, preprocess*, lagpairs*,
  core, nulls*, ifi*, stats`) + `lag_subspace`/`fixed_subspace`, everything else under `legacy/`;
  `scripts/_common.py`; `pyproject.toml`; a 2-page `METHOD.md` (step → function → config → test).
  Six passes (report §6), P0 = equivalence tests first, P1+P2 zero-risk remove most of the clutter.

**Next:** Theo picks candidates by number; apply one section at a time, tests + commit per section;
rejections with a load-bearing reason go to `audit/CONVENTIONS.md`.

---

## ✓ DONE (2026-08-15) — the three 2026-08-07 meeting asks answered

**Nothing is running.** Meeting 2026-08-07 (Nathalie + Tom; handwritten note transcribed
to the vault, `Projects/Hippocampus-V1/Meetings/2026-08-07-Hippocampus-V1-Meeting.md`)
left three asks. All three are pure re-reductions of CSVs already on disk — no refits.
**⚠ Superseded on 2026-08-17 for asks 1/3: the numbers below are on the first-20-trials
cap; see the entry above.**

| ask | built | verdict |
|---|---|---|
| 1. Add the **average IFI across all CCs** | `analyze_cc_ifi_signs.overall_by_animal` → `cc_ifi_overall_bin10*.csv`; black line on `HCV1_cc_ifi_windows_*` | Sits between the ± groups, tracks whichever sign dominates the pair |
| 3. **IFI at ±50 ms — different from 0 overall?** | `overall_direction` → `cc_ifi_overall_test_bin10*.csv`; p in the panel titles; `cc_ifi_signs_tables.md` §OVERALL | **Weakly, and only for CA1→RSC and V1→RSC.** Of each pair's 4 looks (FS × weighting): CA1→RSC 2/4 nominal, all positive (+0.226 p=0.041 FS-excl unweighted; +0.157 p=0.050 FS-incl); V1→RSC 2/4, all positive (+0.120 p=0.019 FS-incl). DG→CA1 (p=0.040 FS-excl unweighted) is 1/4 looks — **not established**. **Nothing survives BH across the 8 pairs**, any look. ⚠ Data = the 12k-bin-capped curves = **first ~20 trials, ≤6 s each**; n = animals with ≥1 significant CC (CA1-RSC 8/12, CA1-V1 7/13 FS-excl) |
| 2. **Cross-correlograms naive vs expert**, also split ±/FF-FB | new `analyze_cc_crosscorr_epochs.py` + `figs_cc_crosscorr_epochs.py` → `HCV1_cc_crosscorr_epochs_{,fffb_}*` | **IFI naive→expert null 0/8, both FS.** Curve *height* is lower in naive at EVERY lag — see the trap below: a baseline offset, not learning |

**Decision taken with Theo (2026-08-15):** the "all CCs" average = **significant CCs only**
(the same population as the sign-grouped lines), unweighted per-animal mean primary,
CC-strength-weighted alongside. Every group line on the IFI-windows figure is now
**per-animal-first** (was over (animal, dim) rows — dims-as-n bands); group means barely
moved, bands widened. New shared helper `src/tom_cca/cc_aggregate.py`.

### ⚠ NEW TRAP — curve height vs the naive epoch on the ALL-TRIALS frozen fit is a baseline offset

Peak r "rises" expert − naive in **8/8 pairs, both FS** (5 at p < 0.05) on `cc_label_track`.
Too clean → checked, then adversarially verified. What it actually is:
- **~60 % of the naive deficit is a lag-INDEPENDENT offset** — the whole curve, baseline at
  |lag| ≥ 200 ms included, sits lower in naive for the cortical pairs (CA1-RSC, CA1-V1,
  V1-RSC, RSC-SUB). `curve_metrics` now splits height into `far_r` (baseline) and
  `peak_minus_far` (coupling-specific): expert − naive on peak-minus-baseline is **null in
  every pair with n > 4** (CA1-RSC p = 0.08, CA1-V1 0.49, CA1-DG 0.86, V1-RSC 0.20), while
  the baseline carries the effect (CA1-V1 Δ+0.019 p = 0.002, CA1-RSC Δ+0.014 p = 0.03). Only
  CA3-DG (n = 4) is peak-specific. A flat offset on frozen axes = slow co-modulation; **running
  speed (+6.6 cm/s naive→expert) is the obvious candidate — still unmeasured (HANDOFF §6.1)**.
- Cohort peak r is naive 0.052 → intermediate 0.067 → expert 0.067 (FS-incl 0.064 → 0.078 →
  0.079); the naive→intermediate rise is **LP-independent** (per-animal Δ vs LP, ρ ≈ −0.01);
  and the **held-out per-epoch refit arms show no naive deficit at all** (`epoch_metrics_bin10`
  CC₁ 0.149/0.139/0.145, int − naive p = 0.84; `lag_subspaces_bin10_epochs` 0.149/0.170/0.147,
  p = 0.97). That is the decisive evidence: coupling strength does not change; what changes
  is how well one whole-session set of axes fits the session's opening.
- "Intermediate already equals expert" (56 %/47 % of animal-pairs) is consistent but NOT
  decisive on its own — they are adjacent 10-trial blocks; and the balanced-trial fit
  (`fixed_subspace_stats`, Δpeak r mixed-sign n.s.) is CC₁-only on 600 bins/trial and a
  30-trial fit set, so corroboration rather than a like-for-like control.
- IFI is insensitive to proportional scaling but **not** to an additive offset (it shrinks
  |IFI|); harmless here only because the contrast is null.
**Applies equally to item 2's `cc_label_track_stats` per-label peak-r contrasts
(`p_FF_peak_r` etc.).** → `GOTCHAS.md`. **Open item:** correlate per-animal Δspeed with the
per-animal baseline offset (needs the speed trace per epoch, not yet exported).

**Registered prior vs outcome (ask 3, FS-excl):** CA1-RSC positive-significant ✓;
V1-RSC positive-significant ✗ (p = 0.16; FS-incl 0.019); CA1-SUB negative p ≈ .01–.1 ✗
(−0.041, p = 0.39); "the other 5 n.s." ✗ — DG→CA1 −0.162 p = 0.040 (1 of 4 looks). Sign
agrees with §B's CC₁ test in 7/8 pairs FS-excl (CA1-V1 flips) but only **6/8 FS-incl**
(CA1-CA3, RSC-SUB flip), and 4 of the FS-excl agreements are with §B means ≤ 0.023 — weak
evidence. ⚠ **§B is a different sample** (whole-session `run_ifi_windows`, ~370 k bins) from
these curves (12k-bin cap = first ~20 trials); per animal-pair the two CC₁ IFIs correlate
only r ≈ 0.4. The like-for-like comparator is **CC₁ of the same table**
(`cc_ifi_cc1_test_bin10*.csv`): CA1-RSC +0.149 p = 0.24 (FS-excl) / +0.224 p = 0.019
(FS-incl) — so the all-CC mean is 1.5× / 0.7× CC₁, **not** "2–3× because weak CCs" (an
earlier draft of this entry said that; wrong).

**Adversarial verification (workflow: 3 refuters recomputing from raw CSVs with their own
code + a judge; 614 k tokens):** every number in `cc_ifi_overall*`, `cc_ifi_overall_test*`,
`cc_crosscorr_epochs*` and `_stats_*` reproduced to float precision, per-animal-first
everywhere, sign convention right, BH-none right — no CSV needed regenerating. Findings were
about claims and guards, all fixed in `9dc50e5`: (1) §B comparator / "2–3×" wrong and the
12k-bin cap unstated → same-table CC₁, cap and n/N now in the md; (2) the IFI-windows figure's
"groups stay apart at other windows" was a **nesting artefact** — cumulative windows always
contain the ±50 ms labelling lags; on the disjoint lags 60–250 ms the +/− gap collapses to ~0
in every pair → reworded; (3) the naive deficit is mostly a baseline offset → `curve_metrics`;
(4) "not learning" was right for weaker reasons than stated → held-out refit arms cited;
(5)/(6) 32 looks in the OVERALL table, hits move with FS × weighting; n ≤ 4 rows no longer
bolded; SEM bands only for n ≥ 3; (7) 7/8 sign agreement is FS-excl only; (8) `cc_aggregate`
now raises if `sig_only` and no `sig` column, and on duplicate (pair, animal) rows.

**Tests: 403 pass** (`cd cca && PYTHONPATH=src python -m pytest -q`, ~50 s).

**Open threads unchanged:** speed confound (HANDOFF §6.1 — still the question to put to
Nathalie/Tom); `MEETING_2026-08-07.md` stale; `fig_rotation_floor` 3-dim angle; Gini
partner-invariance re-test; two direction discrepancies vs §3.0.

---

## ✓ DONE (2026-08-07) — Theo's four follow-ups done

**New doc: `HANDOFF.md`** — the pipeline map (script → CSV → figure), the column
dictionary, which numbers are trustworthy, and the traps. A cold session should read that
before this file. Nothing is running.

**Theo redirected the meeting items** at the start of this session: the earlier pass had
answered adjacent questions rather than the ones asked. What he actually wanted:

| ask | what was built | verdict |
|---|---|---|
| 1. Are there CCs with **different signs** of IFI? Compute IFI at different lags per CC | `analyze_cc_ifi_signs.py` over the refit-per-lag curves | **Null** — sign-mixing is at or *below* chance |
| 2. Label CCs FF/FB **once on all trials**, then track those CCs across epochs | `run_cc_label_track.py` (frozen axes, no cap) | **Labels persist (62–66 %)**, but FF/FB do not diverge with learning |
| — | Cosine similarity of significant CCs to themselves across lags | `run_lag_cosine.py` | **Canonical vectors barely identifiable** |
| 4. Integration window vs IFI curves, naive vs experienced | `analyze_ifi_windows_epochs.py` (frozen axes) | **Null** |

### The one thing that changes how everything else is read

`run_lag_cosine.py`: split-half |cos| of a canonical vector **at a fixed lag** is
**CC1 0.59, CC2 0.39, CC3 0.26, CC4 0.22** against a **0.146** random-vector baseline in
30-D. **Only CC1 has a real identity.** The best-matching dimension at another lag is a
*different* dimension **82 %** of the time, and CC1 loses half its identity by ±250 ms
(0.556 → 0.300, p = 1×10⁻⁴; FS-incl 0.665 → 0.362). This is the measured version of the
bug that hit twice, and it is why every epoch contrast now uses **frozen axes**.

### Item 1 — do different CCs have different IFI signs? NULL

Sign-mixing among significant CCs: **39/49 cells (80 %) vs 40.3 expected (82 %),
p = 0.579** FS-excl; 46/59 (78 %) vs 50.8 (86 %), p = 0.087 FS-incl — at or *below*
chance. No dimension survives FDR within its pair. Two guards were needed to get there:
gate on significance (ungated, the extreme IFIs sit on the weakest CCs, Spearman
−0.25/−0.33) and use the right chance model (P(mixed) = 1 − 2·0.5^k, 94 % at k = 5 — the
first "92 % mixed" figure was *below* chance). Strength-weighting the group means *shrinks*
the FF/FB gap by 6–20 %, a third independent line of evidence for the null.

### Item 2 — label once, track across learning: labels REAL, no divergence

470/494 labelled CCs, 12 learners. **Persistence 62 % (p = 7×10⁻⁴) / 66 % (p = 5×10⁻⁶)**
against a *genuine* 50 % (leave-epoch-out label), animals-as-n, 10/12 then 12/12 animals
above chance. **Not** a theta-ringing artefact: intra-hippocampal minus other is +0.046
(p = 0.42) then −0.032 (p = 0.54), sign flips; highest persistence is RSC-SUB 0.80 and
CA1-SUB 0.76, lowest CA1-DG 0.54. But the interaction ΔFF−ΔFB is **null**: 2/24 tests
FS-excl (1.2 expected, 0 survive BH) and **0/24** FS-incl, the two hits on different pairs.

### Item 4 — integration window vs IFI, naive vs expert: NULL

CC1 expert-vs-naive at ±50 ms: 1/8 pairs in each condition, on **different pairs**
(V1-RSC then CA1-DG). 1/8 pairs have any of 25 windows at p < 0.05 against ~1.2 expected.
FF- and FB-labelled groups likewise null. Curves are largely flat beyond ~50 ms.

### Bugs found and fixed this session (all mine, two by adversarial verification)

1. **`sig` came from a different fit than `cc`** in `run_lag_curves` — `window_subspace`'s
   mask attached by bare index. 19 % of "significant" dims had a **negative** held-out CC;
   the largest CC in the dataset was flagged non-significant. Fixed `d38a833`
   (`lagged.perdim_significance`).
2. **Then reintroduced in a new guise** in `run_cc_label_track` — per-fold refit statistics
   averaged BY RANK, attached to the frozen fit's dims. Signature: `cc_heldout` rose with
   rank in **38/38** cells, negative in 27/38. Fixed `1bae90e`
   (`fixed_subspace.frozen_perm_null`, rank-for-rank on the frozen fit's own correlations,
   ~20× cheaper via `QR(P@Y) = P@Qy`).
3. **The sample cap was harmful, not just unnecessary** — it kept the first 0.4–1.8 s of
   each ~30 s trial, where the naive→expert speed difference is **double** (+12.0 vs
   +6.6 cm/s). `run_cc_label_track` now uncapped; made tractable by
   `trial_lag_moments`/`curve_from_moments` (per-trial sufficient statistics, additive over
   trials → every epoch and leave-one-out curve by masking). 26× more data at ~40 s/animal.
4. **Per-dim null replaced the dominant-dim one** at Theo's call (`0887755`): the old one
   was conservative twice over (max-statistic *and* held-out observed vs in-sample null),
   leaving 0.7 significant dims/cell. Now 3.4 (FDR, leading 10).
5. Minor: variable shadowing (`out = ~mask` clobbering the CSV path, caught by smoke test);
   figures written only to the vault, now mirrored to `cca/figures/` (`figstyle.save`).

### ⚠ OPEN — the one that could undercut the epoch results

**Running speed rises +6.6 cm/s naive→expert** (11/12 animals, p = 0.005 on kept bins).
Nothing controls for velocity; speed sets theta frequency and every epoch contrast here is
a *timing* measure. Uncapping halved the confound, it did not remove it. **Theo has been
asked twice whether this is a known feature of the cohort — still unanswered, and it
should be settled before any epoch result is written up.**

Other open threads are listed in `HANDOFF.md` §6 (rotation-floor figure still plots the
unmeasurable 3-dim angle; two direction discrepancies vs §3.0 finding 3; the Gini
partner-invariance re-test; `MEETING_2026-08-07.md` predates this session).

**Tests: 373 pass** (`cd cca && PYTHONPATH=src python -m pytest -q`, ~50 s).

---

## ✓ DONE (2026-08-03, evening) — item 3 REVERSES: CA1-DG does rotate with lag

**A second verification pass (6 agents) found that the morning's fix was incomplete, and
correcting it flips item 3 from a null to the only positive result in the seven items.**

**The residual bug: fixed the SUBTRACTION, not the AGGREGATION.** `stability()` averaged the
X and Y `angle − floor` terms within an animal and gated estimability on the **pooled mean**
of the two floors — so an area with no usable subspace estimate was averaged in rather than
excluded. 20/71 animal-pairs (FS-excl) have **both** d=1 floors above 70°, yet the pooled
gate called 48/71 estimable.

**Not a conservative error.** Across 2840 area-lags, **corr(floor, angle − floor) = Spearman
ρ = −0.57 (p ≈ 10⁻²⁴²)** — a high floor mechanically yields a negative delta, because the
floor is built from half-data fits while the comparison uses the full window. Unmeasurable
areas therefore drag every pair toward "at floor", i.e. toward the null that was reported.

**ITEM 3, corrected: CA1-DG's communication subspace genuinely rotates with lag.** 8 of 20
lags FS-excluded and 12 of 20 FS-included at the 70° gate, Δ = 20–35°, p(Bonf) to 0.0003,
**in both FS conditions**. CA1-CA3 additionally survives FS-included.

**⚠ The lag COUNT is an analyst choice and is NOT interpretable.** The gate has no principled
value and **both directions are biased** — a tight gate selects low-floor areas (toward
rotation), no gate includes unmeasurable ones (toward the null). `gate_sensitivity()` now
sweeps it: rotating lags go 15 → 14 → 10 → 9 → 1 → 1 across gates 50/60/70/80/90/none.
**Only a pair surviving EVERY setting is claimable — that is CA1-DG (both FS) and CA1-CA3
(FS-incl).**

**⚠ And the gate's premise is itself shaky.** A high floor does not reliably mean "badly
estimated": a near-rank-1 area has degenerate residual PC directions that INFLATE its d=1
split-half angle. Demonstrated in `test_floor_returns_x_area_first_not_y` — a clean 4-unit
population scores 34.7° against a noisy 30-unit one at 11.3°.

**Test coverage added where the bugs actually lived.** Both 2026-08-03 bugs were in the
analysis/driver wiring, which had **zero** tests while `lag_subspace`'s 24 unit tests passed
throughout. New `tests/test_analyze_lag_subspaces.py` pins per-area floor pairing, exclusion
of unestimable areas, animal drop-out, gate-sweep global restoration and the Bonferroni
family; `split_half_floor`'s return order is now pinned by argument swap. **350 tests pass.**

---

## ✓ DONE (2026-08-03, earlier) — verification found 2 real bugs; one positive result was an artefact

**A 14-agent adversarial verification pass over the seven meeting items found two genuine
bugs in code written on 2026-07-28/29, and both mattered.** Fixed, all six affected drivers
re-run, commit `a49fe16`.

**BUG 1 — shared noise floor at the primary dimensionality.** `run_lag_subspaces.py`
exported only the **X-area** d=1 floor as `floor_cc1`, and `analyze_lag_subspaces.stability`
subtracted it from **both** areas' angles. This is exactly the flaw `split_half_floor`
returns two values to avoid — fixed at d=3 during the same session and then reintroduced in
the d=1 path that later became primary. Materially wrong: mean |X−Y| floor difference
**12.8°** (max 74.8°), misplacing the 70° estimability gate in **9/71** cells. Now exports
`floor_x_cc1` / `floor_y_cc1`.

**BUG 2 — fit/projection space mismatch in the fixed-subspace arm.** `fit_fixed`
residualised the third areas with coefficients from the balanced fit trials; the driver then
projected data residualised with **all-trial** coefficients. So the frozen weights were
applied in a different residual space from the one they were fitted in, and "identical
weights across epochs" was not an identical transform. The residualisation now happens once
with fit-trial coefficients and travels with the frozen map.

**⚠ CONSEQUENCE — the one candidate finding was an artefact of bug 2 and is GONE.** The
CA1-RSC integration-window narrowing (logged 2026-07-29 as "worth powering up", then
corrected to "weaker than logged") does not survive the fix. **Items 5/6/7 are now 0/32
significant in BOTH FS conditions.** Verdicts elsewhere held, but per-cell values moved a
lot — `peak_r` changed in 144/150 cells (old-vs-new r = 0.82), `ifi` 144/150 (r = 0.75),
`peak_lag_ms` **r = 0.35**. **Do not quote any pre-2026-08-03 per-cell number.**

**⚠ The theta result was argued from the wrong statistic.** The claim rested on the
*fraction* of curves that ring (33–92 %). `side_peak` fires on **~100 % of pure-noise
curves**, so that fraction is worthless. New `fixed_subspace.side_peak_null` /
`band_occupancy` supply a calibrated null, and the finding survives on the right evidence:
real side-peak offsets sit at a **median 140 ms (7.1 Hz, IQR 120–148)** against a noise
median of 170 ms (IQR 90–270), **KS p = 2.7×10⁻¹¹**, with **85 % of real offsets in the
5–12 Hz band vs 37 % of noise**.

**⚠ "Session" mode is the FIRST ~20 TRIALS, not the session.** `_capped_index` stops at
12 000 bins / 600 bins-per-trial = 20 trials. For an animal with a learning point at trial
69 that window is entirely pre-learning, so items 2 and 3 describe *early behaviour*. The
driver now prints the retained trial range. Documented, not fixed — widening it is a
deliberate decision, not a side effect.

**Corrected counts.** Items 5/6/7 are **32/32 null FS-excluded** (previously logged 31/32).

**Post-fix answers — the qualitative conclusions all held:**
- Item 2: FF/FB still **not separable, 0/8 pairs, both FS**.
- Item 3: **no interpretable rotation** across ±250 ms. Nominally 7/8 (FS-excl) and 6/8
  (FS-incl) rather than 8/8, but the exceptions are noise: CA1-DG's angle-vs-lag profile is
  non-monotonic with several 15–19° values that do not cross, and ~0.4 pairs are expected to
  hit by chance across 8 pairs. **The more informative tail is the other one** — 3 pairs
  (FS-excl) sit significantly *below* their floor, measured confirmation that the half-data
  floor is inflated relative to the full-window comparison it gates.
- Item 4: **0/8, min p = 0.068** FS-excluded, unchanged by the fix.

**Verifier observations worth raising at the meeting, beyond the bugs:** the `sig` flag's
effective per-dimension false-positive rate is **~14 %**, not 5 %; `floor_cc1` exceeds the
70° gate in **27/71** individual animal-pairs even where pair means pass; and two random
directions in the shared 30-dim PCA basis average **81.5°**, not 90°, so that is the real
ceiling for these angle comparisons.

**Standing lesson.** Two bugs, both mine, both found by adversarial verification rather than
by me, and one of them manufactured the only positive result in the set. What is robust here
is the nulls: every item null before the fix was null after, across a change that moved most
per-cell values.

---

## ✓ DONE (2026-07-29) — ALL SEVEN MEETING ITEMS ANSWERED (numbers superseded 2026-08-03)

**Every item from the 2026-07-28 meeting now has an answer, both FS conditions.** Six of
seven are null; the substantive output of the session is three methodological traps found
and one positive control. Next meeting 14:00, 2026-08-07.

| item | question | answer |
|---|---|---|
| 1 | Do different CCs have different IFI? | **No** — no rank dependence |
| 2 | Separate the subspace into FF/FB | **Not supported** — 0/8 pairs separable |
| 3 | How stable are the CCs across lags? | **No detectable rotation** in ±250 ms, 8/8 |
| 4 | FF/FB evolution with learning | **Null** — 0/8, 64 tests, min p = 0.056 |
| 5 | One CC lagged across time | Built (fixed subspace); epoch contrast null |
| 6 | Lagged curves + integration windows, naive vs exp | **Null**; widths need the theta caveat |
| 7 | Fixed subspace, all trials, naive vs exp | **Null** — 31/32 tests, both FS |

**⚠ Three traps found, all now in `GOTCHAS.md`. Each would have produced a confident
wrong answer, and the second nearly did.**
1. The `sig` flag is cleared by floor-level dims at ANY canonical rank (the circular-shift
   bar is the dominant-dim null, tiny in a weakly-coupled cell). Never select on `sig`
   alone.
2. Subspace-angle tests at 3 dims are **unmeasurable** at this N (split-half floor ~78°
   vs a comparison angle of ~75°). Run them at d=1 and gate on `estimable`. This flipped
   the FF/FB verdict from 2/8 "separable" to 0/8.
3. Half-max width on a **ringing** lag curve is a theta half-period, not an integration
   window.

**⚠ The split-half floor is conservative by construction** — it is measured from two
half-data fits while the comparison it gates uses full windows, so it sits above the
comparison's true noise level. Visible in the report's own §G, where cross-window rotations
come out significantly *below* floor. Consequence: every "at floor" verdict here excludes a
*large* effect, not a modest one.

**Cross-check on an existing claim (good news).** §3.0 finding 4's rotation-null is tested
at CC₁ in `bin10_tables.md` §G (CC₁ floors 36–69°), so it survives trap 2. But §G's top-3
column and `figs_report.fig_rotation_floor` use the unmeasurable 3-dim angle — **the figure
should be switched to CC₁; the verdict does not change.**

**Two method discrepancies to reconcile before anything is written up.**
- **CA1-SUB direction.** CA1→SUB here (5/7 animals, both FS) vs SUB→CA1 in §3.0 finding 3.
- **CA1→RSC robustness.** Here p = 0.015 FS-included but p = 0.764 FS-excluded; §3.0 has it
  robust at p = 3.9×10⁻⁴. Different methods (subspace refit at a fixed lag vs IFI integrated
  from one session-level fit) — not directly comparable, but not ignorable either.

**Lint note:** 12 pre-existing F-errors remain in older scripts (`analyze_ifi.py`,
`figs_paired.py`, `learning_changes.py`, `run_transition.py`, `lagged_landmark.py` …). None
are in this session's files; left alone to avoid bundling unrelated changes.

---

## ✓ DONE (2026-07-29, earlier) — meeting items 1, 5, 6, 7; 2/3/4 then running

**Answered so far — three nulls and one positive control.**

**Item 1 — do different CCs have different IFI? NO.** `src/tom_cca/perdim_ifi.py` +
`analyze_perdim_ifi.py` (commit `8dafedb`), a pure re-analysis of `lag_curves_bin10*.csv`
— no refits. No pair shows a CC₁-vs-significant-tail IFI difference FS-excluded; the one
FS-included hit (RSC-SUB, Δ=−0.056, p=0.039) is n=3 and is 1 of 16 tests. Sign agreement
across rank runs 0.25–0.86, nothing significant. **Direction does not depend on canonical
rank.**

**⚠ The `sig` flag does not mean what it looks like.** Falling out of item 1: significance
is scattered across all 30 canonical ranks, not concentrated at the top. Held-out CC does
fall with rank (0.170 → 0.045 → 0.024, Spearman −0.66), but only ~44 % of significant dims
sit at rank ≤ 5, ~10–15 % sit beyond rank 20 at the CC floor, and the significant set is a
contiguous leading run in only ~half of cells. Cause: the circular-shift bar is the
*dominant-dim* null, which is itself tiny in a weakly-coupled cell, so floor-level dims
clear it. **Gate any per-dim label on CC magnitude, never on `sig` alone.** This is why
FF/FB is defined below by fit lag rather than by per-dim IFI sign.

**Items 5 + 7 — fixed subspace, and item 6 — integration windows. Epoch contrast NULL.**
`src/tom_cca/fixed_subspace.py` + `run_fixed_subspace.py` (commits `aaed8be`, `46389de`,
`6ec4735`). Subspace identified once per (animal, pair) on trials **balanced across
epochs**, frozen, then each epoch projected through identical weights — so this tests the
*activity*, with the refit confound removed from every previous epoch contrast. 12
learners, 153 curves per FS condition. Expert vs naive is null on peak r, IFI, peak lag
and half-max width in 31/32 tests FS-excluded and 31/32 FS-included.

**POSITIVE CONTROL — the frozen component rings at theta.** Every pair shows a secondary
lag-curve peak at 113–186 ms (**5.4–8.9 Hz**) in 33–92 % of curves. Recovered without
being told about it; good evidence the fixed-subspace projection tracks real dynamics.

**⚠ This changes how the item-6 integration windows must be read.** Raw half-max widths
look like a clean hierarchy — intra-hippocampal CA1-DG 22 / CA3-DG 26 / CA1-CA3 29 ms vs
cortical CA1-V1 183 / RSC-SUB 220 / V1-RSC 256 ms — but on a *ringing* curve the half-max
width is the **half-period of the rhythm, not an integration window** (the region is
truncated by the first theta trough). Only **RSC-SUB (42 % ringing) and V1-RSC (33 %)**
have widths readable as integration windows. `fixed_subspace.side_peak` detects this and
`fixed_subspace_tables.md` prints the ringing fraction beside every width, so the number
cannot be quoted without its caveat.

**One candidate, and it is WEAKER than first written (2026-08-03 correction).** CA1-RSC's
integration window narrows naive→expert: Δ = −79 ms (p = 0.255) FS-excluded, Δ = −205 ms
(p = 0.029) FS-included. Same sign in both, significant in one, out of 64 tests.
- First check (passed): a low peak would lower the half-max threshold and inflate width,
  but corr(peak, width) = **+0.78 / +0.80** — the wide naive curves have *higher* peaks.
- **Second check (fails): `width_ms` is doubly CENSORED and the result rests on it.**
  ~10 % of values sit at the 0 ms floor (only the peak bin clears half-max) and ~10 % at
  the 500 ms ceiling (the curve never drops below half-max) — only ~78 % are interior.
  Dropping the single ceiling-censored animal takes CA1-RSC from **p = 0.029 (n = 8) to
  p = 0.068 (n = 7)**. **Do not quote the significant version alone.** Now enforced in
  code: `analyze_fixed_subspace.censoring_report` / `width_robustness` print the warning
  automatically and write a censoring table into `fixed_subspace_tables.md`.

**Items 2 + 3 — FS-EXCLUDED ANSWERED (16 animals, commit `705213b`).**

*Item 3, how stable are the CCs across time lags?* The **CC₁** subspace never separates
from its own split-half floor anywhere in ±250 ms — 8/8 pairs, all p(Bonf) = 1, Δ at
±50 ms spanning −9.2° to +10.2°. **State it with the floor: the floor is ~56°, so this
excludes a *large* rotation, not a modest one.**

*Item 2, separate the subspace into FF/FB?* **The data does not support the split.** The
FF-vs-FB subspace angle never clears the floor in any pair (−9.1° to +15.3°, all n.s.);
connection-specific Gini does not differ in any pair; strength differs in 1/8 (RSC-SUB
+0.012, p = 0.045, n = 7 — 1 of 8, would not survive correction). At ±50 ms these are
**one subspace read at two delays**, not two directions of flow.

**⚠ The 3-dim version of this test is UNMEASURABLE — see the new GOTCHA.** Split-half
floor ~77° at d=3 vs ~56° at d=1. At d=3 the lagged angle (~75°) has no headroom, so
"not above floor" would report a power failure as stability. The FF/FB gate said 2/8
pairs "separable" at d=3 and **0/8** at d=1. Always read the `estimable` column first.

**Sanity check worth keeping.** The group-mean direction recovers hippocampal anatomy in
**3/3** intra-hippocampal pairs — CA3→CA1, DG→CA1, DG→CA3. But per-animal sign
consistency is at chance (4/7, 6/11, 3/5), so the mean is carried by magnitude in a few
animals: suggestive, not established. The animal-consistent directions are instead
**RSC→SUB (6/7)** and **V1→RSC (7/9)**; the latter agrees with §3.0 finding 3.

**⚠ DISCREPANCY to reconcile.** CA1-SUB reads **CA1→SUB** here (5/7 animals), whereas
STATE.md §3.0 finding 3 lists **SUB→CA1** among the FS-robust tight-lag flows. Different
method (subspace refit at a fixed lag vs IFI integrated from one session-level fit), so
not directly comparable — but it should be reconciled rather than ignored.

**RUNNING (items 2, 3, 4):** `run_lag_subspaces.py`, run 1 of 4 done 00:37; remaining
(session × both FS,
then `--epochs` × both FS), ~2 h each. Also still running: the gini3 FS-included
trajectory from 2026-07-28 (PID 20601, 7.7 h elapsed, animal 75 / 1540 rows — verified
*computing*, not wedged: 45 s CPU per 5 s wall). **Do not clobber either.**

**New machinery (commit `aaed8be`, 41 new tests, 335 pass).** `lag_subspace.py` fits CCA
at a fixed segment-aware lag and exports neuron-space **weights**, which nothing did
before — that is what makes cross-lag subspace comparison possible at all. Three design
points worth not re-deriving: (1) `segment_lag` aligns the confound to each area's *own*
timepoints, since partialling a third area out of Y at X's times would inject the very lag
being measured; (2) `lag_sweep` hoists both the confound regression and the PCA out of the
lag loop — leak-free (lagging only selects rows), ~20× cheaper, and every lag shares ONE
neuron→PC basis so a cross-lag angle reflects CCA rotation rather than PCA jitter;
(3) `split_half_floor` returns a **separate floor per area** because the floor rises with
fewer units and these pairs are lopsided (CA1 vs SUB) — a shared floor would misplace the
threshold item 3's verdict turns on.

---

## ✓ DONE (2026-07-28) — Gini partner-invariance found; meeting sets a new agenda

**⚠ The headline metric was measuring the wrong thing.** `gini_x` / `gini_y` — the basis of the
§3.0 finding 2 "participation broadens" headline — is **provably partner-invariant**. In
`core.cca_fit`, `A = Vx @ diag(1/sx) @ Uc[:, :d] * scale`; `membership.subspace_contribution`
takes the *unweighted* L2 row-norm over all `d` retained dims, so when `d = rank(X) ≤ rank(Y)`
the square-orthogonal `Uc` cancels exactly out of every row norm and the partner area drops out
of the arithmetic. Verified three ways: (i) a new test pins the identity; (ii) synthetic X against
a cc₁=1.0 partner vs a pure-noise partner gives bit-identical contributions (max diff 4×10⁻¹⁵);
(iii) in the shipped data CA1's `gini_x` correlates at **median r = 0.981** across five different
partner areas (r = 1.000 for CA1-CA3 vs CA1-SUB), while the CCA-free Pearson control on the same
cells sits at r = 0.392. So the shipped metric is an **area-intrinsic** readout of a population's
own whitened-PCA loading geometry — legitimate, but *not* participation in a communication
subspace. Commits `9421092` (+ TDD) and `28f0c44`.

**Two connection-specific definitions added** (`membership.subspace_contribution_connection`),
computed in the same pass by `run_trajectory` / `run_epochs` / `run_transition`:
- `gini_*_conn` — contributions weighted by each dim's held-out canonical correlation (clipped at
  0). Always defined; reduces exactly to the area-intrinsic Gini when the canonical spectrum is
  flat (Gini is scale-invariant). Caveat: conflates concentration with how many dims survive
  clipping.
- `gini_*_sig` — restricted to significant dims only. Cleaner, free of that confound, but **NaN
  when `n_sig` = 0**, which is 56 % of CA1-RSC windows — precisely where the area-intrinsic
  version misleads most.
Both drivers now route through `membership.*` instead of an inline `np.linalg.norm`, so the two
definitions live in one place. Epoch/transition also export `contrib_conn` per neuron, which makes
the across-partners aggregate in `figs_area_gini.py` meaningful for the first time (it was
averaging near-identical copies of one partner-invariant vector).

**RUNNING (do not clobber):** `run_trajectory.py … --include-fs --out trajectory_gini3_bin10`
(PID 20601, launched 16:41, still going at 23:40 — ~1192 of an expected ~1900 window-rows).
FS-excluded finished 16:41 (`trajectory_gini3_bin10.csv` 1927 rows +
`trajectory_gini3_bin10_dims.csv` 37 602 dim-rows). Log: `results/gini3.log`.

**⚠ §3.0 finding 2 is NOT yet re-tested.** The Gini↓ trajectory result stands only for the
area-intrinsic definition. Whether participation-broadening survives under `gini_*_conn` /
`gini_*_sig` is **open** and is the first thing to check when the FS-incl run lands. Do not cite
the broadening headline as a communication-subspace result until then.

**Meeting — Nathalie + Tom, 2026-07-28 (next meeting 14:00, 2026-08-07).** Seven asks, mapped to
the repo:
1. *Do different CCs have different IFI?* — data already on disk
   (`results/lag_curves_bin10{,_fsincl}.csv`, dims 1–30 × lags ±25 bins, both FS). Pure
   re-analysis.
2. *Separate subspaces into FF/FB, track evolution with learning* — **nothing exists**; no FF/FB
   machinery anywhere in `src/`. Sits on top of (1): sign of per-dim IFI → FF/FB tag.
3. *How stable/similar are the CCs across time lags?* — `lagged.heldout_lag_curve_flat_perdim`
   returns CCs only, not weights; needs a weight export + principal angles (reuse
   `subspace_stats`).
4. *Separate lagged curves by FF/FB, and naive vs exp* — follows (1)+(2).
5. *Go back to one CC, lagged across time* — with (7), a fixed-subspace method change.
6. *Plot lagged curves + integration windows, naive vs exp* — drivers exist but are
   session-pooled; need an epoch restriction.
7. *Fixed subspaces identified across all trials, naive vs exp* — new driver.

**Stats flag raised at the meeting-mapping stage:** tagging dims FF/FB by the sign of their IFI is
a selection on a noisy statistic — per-dim held-out CC for the weak pairs is 0.02–0.09, where the
IFI sign is near a coin flip. Gate the tag on (a) dim significance vs the circular-shift null and
(b) an IFI magnitude floor from the same surrogate, or the FF/FB split is noise-splitting followed
by a circular test.

---

## ✓ DONE (2026-07-27) — doc reconciliation after a 6-week gap; nothing running

**Nothing is in flight.** Last analysis output 2026-06-18 16:11; last commit `7b6f341`
(2026-06-18). No live processes. The three entries below that were flagged ⏳ IN PROGRESS /
RUNNING have been **retitled ✓ DONE** — their batches did finish; the headers were stale, not
the work.

**Verified this session (evidence, not assumption):**
- **Tests: 265 pass** (48 s, `PYTHONPATH=src python -m pytest -q`).
- **All batches completed.** `missing_smoothed.log` died at `[3/6]`, but steps 5–6 (FS-excl
  epoch + transition weights) were relaunched and finished (`fsexcl_weights.log` DONE 19:11;
  `epoch_weights_bin10.csv` 18:58, `transition_weights_bin10.csv` 19:11). Lag-curves both FS
  and 5-trial-bin trajectories both FS all completed. **Do not re-launch these.**
- **Report figures are sound:** all **62** `![[…]]` embeds across the vault report + appendix
  resolve; the main report is **37 bin10 embeds + graphical abstract**, i.e. no 25 ms leakage
  into the smoothed body.

**`STATE.md` REWRITTEN** — it was 6 weeks stale *and wrong*: it still asserted the Arm A
temporal analysis was "designed but **not run** — no `temporal_arm_*` pkls exist", while §3
carried only landmark/spatial findings. Since STATE.md declares itself the tie-breaker, a cold
session was being actively misled. Now restructured: three arms (§2), **temporal arm promoted
to PRIMARY with its own verdict block (§3.0)** — hierarchy table, strength-null, Gini↓
participation-broadening + the learning-vs-experience caveat, CA1→RSC flow existence
(p=3.9×10⁻⁴), rotation-null, early-trials, KCCA-linearity, units-of-analysis policy — landmark
/spatial demoted to §3.1 "superseded, retained for the record" with the pseudoreplication
caveat kept; §4 gains the temporal canonical config; §5 open decisions rewritten
(learning-vs-time-on-task is now the binding question); §6 gains the temporal reproduce block
(flags verified against argparse — it is `--include-fs`, not `--fs-incl`).

**OUTSTANDING (nothing blocking):**
1. **Report number reconciliation** (Theo's domain, vault prose) — the editor's note at the top
   of the report still stands: §3.1–3.6 prose carries pre-smoothing numbers; **§3.9 tables are
   authoritative**. Two concrete discrepancies confirmed against `bin10_tables.md` §A: the
   graphical-abstract/§4 CA3-DG **CC₁ 0.36** is stale (smoothed = **0.315** FS-excl / **0.409**
   FS-incl), and "**~2.9 sig dims**" is stale (smoothed = **5.74** FS-excl / **6.14** FS-incl).
   Also §6 Provenance names branch `cca-consolidation` (we are on `main`) and cites the old
   `trajectory_windows.csv` rather than the `_bin10` set.
2. **Orphaned FS-excl bin10 early-trials.** `early_trials_blocks_bin10.csv` +
   `early_trials_projected_bin10.csv` exist (Jun 17 18:20) but `analyze_early_trials.py` never
   ran on them (no `early_trials_summary_bin10.csv`) and no `HCV1_early_*_bin10` figures exist.
   The expensive block-refit is already paid for; analysis + figs are minutes. FS-incl bin10
   early-trials stays deliberately skipped (pathologically slow, low value).
3. Vault `Hippocampus-V1-Tasks.md` last triaged 2026-05-07 — "Refine temporal CCA analysis"
   still open, plus Granger/LFP, synthetic-data CCA recovery, hyperparameter re-sweep.

**Git hygiene decided this session:** results CSVs/logs and `AGENTS.md` committed; the two
~3 MB `lag_curves_bin10*.csv` and the 27 MB `share_tom_2026-06/` deliverable **gitignored**
(never commit data — they regenerate from `run_lag_curves.py` / the figure scripts).

---

## ✓ DONE (2026-06-17 cont.) — FULL REPORT REGEN: smoothed-10 ms primary, both FS co-primary

**✓ REPORT REGENERATED (structure + tables + figures + framing).** Main report restructured
(`scripts/_regen_report.py`, one-shot, removed): **§3 = the smoothed-10 ms results (primary,
both FS)** — promoted from the old §8; **25 ms raw + KCCA → new file
`… — 25 ms Appendix.md`** (§A/§B). Regenerated `bin10_tables.md` (smoothed + paired-t) spliced
into §3.9. §3.7 figure block expanded to the full both-FS set + the new analyses (lag-CC curves,
paired panels, per-area Gini, weight CDF, 5-trial-bin trajectories). Status box, §3 intro, §2.13/§2.16
stats policy rewritten to **paired-t test of record** (Wilcoxon dropped); §5→appendix cross-refs fixed.
38 embedded figs all resolve. Backup: `…Report.md.prebin10regen.bak`. **t-test switch** applied across
figs_paired/area_gini/rotation_cc1/early_trials/units/report/trajectory_dims/stats_tables +
analyze_bin10_full + analyze_trajectory_dims; `paired_stats.paired_t`/`welch_t` added (drop-in).
**✓ COMMITTED `6f8aca2`** (paired-t switch + R2 lag curves + B/D + bin10/both-FS pipeline + GOTCHA;
21 code/doc files, no data). **Deck regenerated** (`make_smoothed_deck.py` → 113 slides, 102 figs).
**FS-excl per-area-Gini done** (weights re-run after the early-trials batch *stalled* at step 3/6 for
3 h with no worker — killed it; relaunched only the fast FS-excl epoch+transition weights).
**Tom deliverable package** at `cca/share_tom_2026-06/` (41 svg + 28 csv + README data-dictionary +
METHODS.md/.docx [pandoc, 59 native eqns] + STATS.md) — UNCOMMITTED (deliverable, has data).
**REMAINING:** (1) sentence-level **number reconciliation** of §3.1–3.6 prose + §4 cited p-values vs
the regenerated §3.9 tables (flagged in an editor's note; qualitative narrative holds — synthesis
pass, Theo's domain); (2) **FS-incl early-trials (10 ms)** deliberately skipped (pathologically slow
block-refit, low value — 25 ms early-trials is in the appendix). Report prose UNCOMMITTED (vault).


**Decision (Theo):** regenerate the ENTIRE report on smoothed-10 ms, *consistently and only*; move the
25 ms results + KCCA (§5) to a separate appendix file; FS-excl AND FS-incl **co-primary** in the body.
Plus (R2) add held-out **lagged-CC-curve** figures (the cross-correlation profile behind the IFI),
n=animals (CC1) and n=subspaces (sig dims), mean ± shaded SEM, per pair.
**Why:** §8 (10 ms) was built on the RAW run, but bin10 CSVs were later overwritten with SMOOTHED →
§8 numbers went stale & contradicted §2.2 (e.g. CA3-DG CC1 0.271→0.315); and the new paired/area-Gini/
CDF figures were never in the report.

**BUILT this session (code; commit pending at wrap):**
- R2: `lagged.heldout_lag_curve_flat_perdim` (per-dim held-out lag curve; CC1 wrapper kept) + 2 TDD
  tests (16/16 pass). `run_lag_curves.py` (held-out per-dim CC vs lag ±250 ms; n_sig via
  window_subspace; capped 12 k samples via CONSECUTIVE within-trial blocks — see new GOTCHA) +
  `figs_lag_curves.py` (8-pair grid, both n-units). VALIDATED (CC1 peaks at 0 ms).
- `figs_trajectory_bins.py` (B/D: cc1/n_sig/ifi/gini_x/mincc across the 6 bins, learners vs non).
  VALIDATED on B/D fsincl.
- Script edits for bin10 + both-FS: figs_units (epoch_*_bin10), figs_early_trials (bin10 in/out),
  run_early_trials (--smooth-ms + bin tag), figs_area_gini (FS-param), figs_rotation_cc1 +
  figs_stats_tables (both-FS naming, incl. directionality fsincl). Lint F-clean (uvx ruff --select F).
- B/D 0-rows bug FIXED (even-stride cap preserved whole trials for by-trial CV; GOTCHA logged).

**DATA STATE (smoothed 10 ms):** DONE = trajectory/ifi_windows/epoch/transition/dims (pre-existing
smoothed) · B/D fsincl (`trajectory_bins_bin10_fsincl`, 438 rows) · lag-curves FS-excl
(`lag_curves_bin10`, 76 k rows). RUNNING (detached) = `missing_smoothed.log` batch [lag-curves fsincl
→ early-trials ×2 FS → FS-excl epoch+transition WEIGHTS ×2] + B/D FS-excl (`trajbins_excl.log`).

**FIGURES regenerated smoothed, both-FS (batch-independent):** figs_report (levels/slopes/gini_traj/
ifi_traj/direction/rotation_floor/gini_control/learnervsnon), trajectory_dims (pooled), ifi_windows,
stats_tables (slopes+directionality), rotation_cc1 · PLUS lag-curves fsexcl, B/D fsincl (5 metrics).

**NEXT (post-batch):** (1) regen batch-dependent figs both-FS — figs_paired, figs_units, figs_area_gini,
figs_early_trials, figs_lag_curves fsincl, figs_trajectory_bins fsexcl; (2) re-run `analyze_bin10_full.py`
→ smoothed `bin10_tables.md` (extend w/ mincc, per-area Gini, lag-curve direction summary); (3) REWRITE
report — §2 framing → 10 ms primary; §3 rebuilt on smoothed both-FS (promote old §8 structure; fold in
early-trials/units/rotation/lag-curves/paired/area-Gini/B-D); §4 synthesis numbers; MOVE 25 ms + KCCA to
`… — 25 ms Appendix.md`; (4) deck + commit + log. Report prose stays UNCOMMITTED (vault, for review).

## ✓ DONE (2026-06-17) — Tom's meeting plots (smoothed 10 ms, FS-included)

Tom's request (Slack): paired plots on the smoothed data for a meeting next week. **Built & in the
vault (`figs_paired.py`, committed):** 12 paired-panel figures (per pair, per-animal connected dots,
mean±SEM, paired t+Wilcoxon) — naïve-vs-trained (learners) for cc1/n_sig/ifi/gini_x/mincc, and
uncued-vs-cued for the same (learners) + cc1/mincc (all animals). `mincc` = mean held-out CC over the
**min # of sig dims** shared by the two conditions.
**Tom's Qs answered:** (1) Gini(area X) = Gini coeff of the dominant-dim canonical weight vector
(membership.gini, N/(N-1) corr); (2) d(metric)/d(trial_frac) = per-animal least-squares slope vs
session-fraction; (3) page-15 V1→CA1-when-cued **does NOT hold at ±50 ms smoothed** (CA1-V1 Δ=-0.001
p=0.96; the V1→CA1 lead was a wide-window 25 ms feature). Cued trends are toward CA1 *leading*.
**Tom's 3 directives:** (1) **Gini per AREA across all partners** (not per pair) → aggregate each
area's per-neuron `contrib` across its partner pairs (partner-normalised), then Gini; (2) **save
weights** → done (contrib export added to run_epochs/run_transition, committed; **re-run IN FLIGHT**
smoothed FS-incl → `epoch_weights_bin10_fsincl`/`transition_weights_bin10_fsincl`); (3) **5-trial
bins** [1-5,6-10,…,26-30] for the trial-resolved plots.
**Methods §2 rewritten thorough (while runs go; UNCOMMITTED vault prose).** Report §2 expanded
from 11 → 17 subsections covering the *current* pipeline: Gaussian σ=2.5 ms smoothing + 10 ms
binning + z-score (§2.2), masks/epochs (2.3), partial CCA (2.4), PCA-k-per-analysis (2.5), CCA (2.6),
leak-free held-out CV (2.7), n_sig/MI/**mincc** (2.8), IFI + segment-aware held-out lag curve +
**±50 ms headline / ±250 ms window sweep** + per-dim IFI (2.9), contributions + **Gini per-pair AND
per-area-across-partners** + weight CDF (2.10), rotation (2.11), trajectory + **5-trial bins** (2.12),
epoch + **paired** contrasts (2.13), early-trials (2.14), transition (2.15), **units-of-analysis
(animals-as-n vs dims-as-n) + stats + pooled-vs-learners** (2.16), KCCA pointer (2.17). All stale
§2.x cross-refs in §3/§4/§5/§8 + status box remapped to the new numbering.

**DONE since:** contrib re-run finished (epoch_weights/transition_weights_bin10_fsincl). Built
(committed): **per-area Gini across partners** (`figs_area_gini.py`) + **weight CDF (Fig-4b Lorenz)**
for BOTH naïve-vs-trained and uncued-vs-cued (supersedes per-pair gini panel). **5-trial-bin driver**
`run_trajectory_bins.py` (metrics + dims + contrib, learners + non, fixed [1-5..26-30] bins) written
+ **launched smoothed FS-incl (detached, ~1-2 h)** → `trajectory_bins_bin10_fsincl{,_dims,_weights}`.
Methods §2 rewritten thorough (17 subsections; vault, uncommitted).
**STILL TO BUILD:** (b) B/D metric-vs-trial-bin plots (cc1/n_sig/ifi/per-area-gini/mincc over the 6
bins, learners vs non) — needs the 5-trial-bin run (in flight); (a) FS-incl/learners summary-schema
variant. Then assemble Tom's deck/figure set.

## ✓ DONE (2026-06-13) — Gaussian-smoothed (Buzsáki) 10 ms re-run, both FS

**Why.** Align preprocessing with the cited paper. **Gonzalez & Buzsáki 2026 Methods** ("Preprocessing
single-unit spiking", read from Zotero `IQYYYVBB`): *convolve each unit's spike train with a Gaussian
**2.5 ms s.d.** → smoothed spiking activity → z-score (subtract mean rate / divide by s.d.); bins
tested 0.8/8/40 ms* (our 10 ms ≈ their 8 ms primary). **NB:** the `_gf` data fields are **gap-filled**
behaviour/spatial, NOT a smoothed temporal stream — there is no pre-smoothed spike stream, so we
compute it. (The `HC_V1_Code/` MATLAB is NOT Tom's — disregarded.)

**Implemented (committed 69dfe44).** `rebin_spikes(gaussian_sd_ms=)` convolves the 1 ms train per
unit (`gaussian_filter1d`, ±4σ halo, mass-conserving) before binning; `cfg.gaussian_sd_ms` +
`--smooth-ms` on all 4 drivers; cache key now includes σ (it omitted it → first validation hit a
stale cache). σ=0 byte-identical to raw (38 dataio tests pass). **Validation** (animal 28, CA1–V1,
10 ms): σ=2.5 raises held-out CC₁ **0.121→0.151** (denoises sparse bins) while **n_sig & IFI hold**
(tight kernel preserves the ±50 ms directionality) — the intended effect.

**Launched (detached, both FS, `--smooth-ms 2.5`).** Batches `/tmp/run10sm_{fsexcl,fsincl}.sh`, logs
`results/run10sm_{fsexcl,fsincl}.log`; ifi-sweep → trajectory → epochs → transition. **Overwrites the
`*_bin10*` outputs** — 10 ms now means *Buzsáki-smoothed* (prior raw-count 10 ms is commit f4de18c).
~12–16 h (sweep is the uncapped slow stage, as before). Expect CCs up (less magnitude-compression),
directionality + de-sparsification + nulls preserved.

**NEXT when both DONE.** Re-run `analyze_bin10_full.py` (tables) + all bin10 figures + re-splice §8 +
deck + abstract on the smoothed data; update §8 prose ("10 ms now Buzsáki σ=2.5 ms smoothed; CCs no
longer magnitude-compressed"). Detached run survives turn boundaries; `run_in_background` watchers do
not — poll the logs on check-in.

## ✓ DONE (2026-06-13) — 10 ms full-suite re-run, both FS → report §8 + bin10 figures + deck

**RESULT (committed f4de18c; report §8 has the complete pair-by-pair tables, every metric ×
comparison × pair × both FS).** The 25 ms dissociation **reproduces** at 10 ms: hierarchy holds
(CCs lower — sparser bins — but same ordering); **coupling strength flat** (no slope, no
naive→expert); **no reorientation** (CC1 rotation ≤ floor); **de-sparsification robust, if anything
stronger** (Gini↓ CA1-RSC p<0.001 + most pairs). **Directionality sharpens** — three FS-robust
*tight-lag* flows (curves peak ≤±50 ms, Table C): **CA1→RSC** (p=0.001), **V1→RSC** (p=0.01–0.02 —
only a *tendency* at 25 ms; newly resolved), **SUB→CA1** (p=0.006). Flows **coherent across canonical
dims** (Table D) except **CA3–DG flips** (dominant CA3→DG, sub-dominant DG→CA3 — bidirectional).
Directional **change** with learning is **underpowered, not absent**: null by slope AND naive→expert,
but V1–RSC keeps its naive→expert decrease (Δ=−0.027, n=6, p=0.16) — same direction/size as 25 ms.
`analyze_bin10_full.py` regenerates the tables; figs are bin-tagged (`*_bin10*`); deck = 147 slides;
graphical abstract updated (3 solid flow arrows). **25 ms stays the magnitude reference** (10 ms CCs
lower). Outstanding: transition FS-incl re-ran clean (59 rows) after the time-based cap fix.

**Why.** Directionality picture needs finer temporal resolution. Decision (user): re-run the
**full suite at 10 ms, both FS**; **headline IFI integrates over ±50 ms**, but the IFI window
sweep keeps **curves out to ±250 ms**.

**Params / mapping (10 ms bins).** Headline IFI = `--max-lag 5` (±50 ms) on
trajectory/epochs/transition; sweep = `--max-lag 25` (±250 ms curves, headline read at window 5).
`information_flow_index` integrates CC1 over the full ±max_lag, so max_lag *is* the integration
window. Feasibility OK: raw data is **1 ms** `binned_spikes` (the `_50ms` names are legacy), ~8 Hz,
population activity in every 10 ms bin → CCA estimable. **Caveat:** absolute CCs drop at 10 ms
(sparser) — the strength/Gini de-sparsification headline stays at **25 ms**; 10 ms is a
finer-timescale robustness + sharper directionality view.

**Driver changes (committed 46cba64).** `--max-lag` added to run_trajectory/run_epochs/run_transition;
run_transition also gained `--bin-ms`/`--tag`/`--include-fs` (was hardcoded 25 ms). All defaults
unchanged. **Outputs are bin-tagged** (`trajectory_w15_bin10*`, `epoch_metrics_bin10*`,
`transition_*_bin10*`, `ifi_windows_bin10*`) so 10 ms coexists with the committed 25 ms results.

**Launched (DETACHED).** Two concurrent FS batches (`/tmp/run10_{fsexcl,fsincl}.sh`, logs
`results/run10_{fsexcl,fsincl}.log`), each sequential: ifi-sweep → trajectory → epochs → transition.
~12 h wall-clock (64 GB RAM ample, 16 cores, OMP capped 4/proc). **Smoke test PASSED** — 10 ms
gives valid finite CCA/IFI. Completion watcher **b2u68781m**.
- **SLOW — `run_ifi_windows` & `run_epochs` have NO sample cap** (trajectory/transition cap at
  MAX_SAMPLES=6000). At 10 ms the sweep refits CCA at 51 lags × 5 folds over the **whole ~800k-bin
  session** → ~70 min/animal on large sessions; 9 h in, only 11/16 on stage 1. **Decision (user):
  LET IT RIDE** — no cap, no restart (keep the 9 h, full-session lag curve). Revised ETA ~10–12 h
  more (done ~midday 06-13). **WATCH:** epochs (stage 3) is also uncapped — it uses max_lag=5 (11
  lags, lighter) on per-epoch data (smaller than the full session), so likely OK, but if stage 3
  stalls, surface the cap option again (don't cap unilaterally — user chose no methodology change).
  Capping for a *future* run would be ~120k whole-trial-preserving bins (ample for a stable curve).
- **GOTCHA (cost an early restart):** the first launch used `run_in_background:true` for both batches
  **plus** a smoke watcher = 3 harness-tracked tasks; the harness **evicted** a batch (≈2-task limit)
  ~14 min in (NOT OOM — 64 GB, 87 % free). Fix: launch long compute **fully detached**
  (`(nohup bash X.sh >/dev/null 2>&1 </dev/null &)` → PPID 1, new session, harness can't evict) and
  keep ≤1 run_in_background watcher. Current runs are detached; verify with `ps -o ppid` = 1.

**NEXT when both complete.** (1) `analyze_ifi_windows` — fix the **headline at ±50 ms (window 5)**,
plot curves to ±250 ms; (2) re-point figs/analysis at the `*_bin10*` outputs and regenerate figures
+ stats + deck + the directionality inventory (stable vs changing) at 10 ms; (3) compare to 25 ms.

**Directionality inventory at 25 ms (reference for the 10 ms comparison).** *Stable/established
flows:* CA1→RSC (robust, p<0.001), SUB→CA1 (p=0.009), RSC→SUB (p=0.047). *Changing with learning:*
CA1–DG IFI↑ (trajectory slope +0.020, pooled-16 **p=0.01** — the clearest), CA1–V1 IFI↑ (FS-fragile),
V1–RSC IFI↓ (underpowered-but-real, lag falls, 5/6 animals). The graphical abstract was being
revised to carry this (stable + changing) when the 10 ms decision was made.

## ✓ DONE (2026-06-12) — pooled (all-animals) replication + clearer stats + deck overhaul

Per user (deck feedback): drop spatial sweeps, KCCA→verdict-only, **replicate every figure
pooled (all 16 animals, learners + non)**, and make **stats way clearer** (both on-figure AND
summary tables). Committed **4c1b7ed**.
- **Pooled variants** for all *poolable* figures via a `pool` flag (`_cohort`/`_ptag` in
  `figs_report.py`; also `figs_trajectory_dims`, `figs_rotation_cc1`): levels, slopes-heatmap,
  gini/ifi trajectories, direction, rotation-floor, gini-control, dims-forest, transition.
  **Epochs/units stay learners-only BY DESIGN** — naive/intermediate/expert are defined relative
  to the learning point, which non-learners lack (so `run_epochs` is intrinsically learner-only;
  the `lp_rel` axis likewise). Early-trials were already all-animal. Told the user this.
- **Clearer stats:** non-occluding 2-line colour-coded panel titles (p with */**/*** tier ·
  signed slope+arrow · n; green=sig; test named once in the suptitle); asterisk-tier bar tops;
  adaptive-contrast heatmap cells; rotation-floor gained its missing paired signed-rank.
  `figs_stats_tables.py` = colour-coded summary-table FIGURES (learning slopes learners-vs-pooled;
  held-out IFI window sweep) → the "summary slide" view, globbed into a new deck section.
- **Result that fell out:** pooling **strengthens** the de-sparsification (CA1-RSC ✶✶✶ p<0.001;
  CA1-V1/CA1-SUB/V1-RSC cross into significance) — direct support for *experience, not
  learning-specific*. CA1-DG IFI↑ is the one directionality effect robust to pooling.
- **Deck:** `HCV1_all_figures.pptx` now **97 slides / 86 figs** (sweeps gone, KCCA = 2 vs-linear
  slides, Statistics-at-a-glance up front, pooled+learners throughout).
- **GOTCHA:** `figs_report` names figures by FS-tag only, **not by bin** — bin25/bin50 collide and
  overwrite. So the bar/trajectory figures stay 25 ms (report primary); bin50 robustness rides on
  the dims-forest + stats tables, which carry the bin in their filename.
- **Loose end:** bin50 FS-incl stuck at 15/16 (last animal slow/hung); its dims-forest is at 15
  animals. Not committed (in-flight). The driver (PID 44335) is still alive on that animal.

## CURRENT STATE (2026-06-12, early hours) — post-review; trajectory dims-as-n re-run IN FLIGHT

**Where things are.** The 18-comment report review (06-10/06-11) is largely landed in the vault
report (`Projects/Hippocampus-V1/…Learning Report.md`). Git policy changed to **commit straight to
`main`, no branches** (CLAUDE.md updated). Recent commits carry: W=15 trajectory re-run at **both
bins (25 & 50 ms)**; **parametric stats** (paired-t + random-slope LMM alongside Wilcoxon, §2.10);
the **held-out segment-aware IFI lag-window sweep** → upgrades CA1→RSC directionality to a clean
**p=4×10⁻⁴** (§3.2); **CC1-only rotation** → no-reorientation holds at the dominant dim (§3.4, with
a synthetic regression test proving the orthogonal top-3 split-half floor is the 1-D-subspace
signature, not a bug); **KCCA upgrade** (30 shuffles, ±8 lags) → still *largely linear* (§5);
**transition dims-as-n** (§3.5, `transition_dims.csv`, numbers re-verified this session); the
**early-trials battery** (§3.8); and a **graphical abstract** (`attachments/HCV1_graphical_abstract.svg/png`).

**RUNNING (do not disturb).** One sequential driver (wrapper **PID 44335** → `results/traj_dims.log`)
re-runs the W=15 trajectory **per-dimension** for the dims-as-n view of §3.2/§3.3, 4 stages in order:
`trajectory_w15_bin25` (FS-excl) → `…bin25 --include-fs` → `…bin50` → `…bin50 --include-fs`. Each
writes BOTH the window CSV `trajectory_w15_bin{25,50}{,_fsincl}.csv` **and** the per-dim
`…_dims.csv` (one row per canonical dim). ~5 min/animal, ~5–6 h total. Completion watcher **bxdq9o98x**
armed (prints final per-stage animal counts on driver exit).

**⚠ GOTCHA caught this session.** The running driver opens the window CSVs in `w` mode → while a
stage runs, `trajectory_w15_bin25.csv` is **truncated** to the animals done so far. The 23:56 fsexcl
trajectory figures had been built from ~2 animals. Fixed by regenerating from the **full HEAD copy**
(`git show HEAD:…trajectory_w15_bin25.csv`) via a protected stem. **Rule: only regenerate trajectory
figures when the window CSV shows 16 animals.** (Added to GOTCHAS.)

**Report §3.2 finished this session (UNCOMMITTED prose in the vault).** The §3.2 header already
promised the CA1→RSC held-out flow but the *body* still told the old "in-sample, weak, planned"
story — now wired in: a headline paragraph (CA1→RSC $+0.079$ at $\pm50$ ms, $t_{11}=5.0$,
$p=3.9\times10^{-4}$; nested-window Bonferroni $p\approx5\times10^{-3}$; SUB→CA1 & RSC→SUB exploratory),
the `HCV1_ifi_windows_bin25_fsexcl.png` embed, fixed stale caveats (§2.7 "in progress"→done; §3.2
"held-out planned"→done), and a §4 synthesis row + headline edit distinguishing *existence* of the
CA1→RSC flow (robust) from its *change* with learning (weak). Numbers re-verified vs
`ifi_windows_bin25.csv` this session. **Vault report edits are on disk, not committed — for user review.**

**STAGED & validated (committed cff6a3e).** `scripts/analyze_trajectory_dims.py` (shared
`compute_table()` → animals-as-n signed-rank | dims-as-n cluster-robust LMM | dims-as-n naive OLS,
per pair × {ifi,cc} × {trial_frac,lp_rel,performance}) and `scripts/figs_trajectory_dims.py` (three-unit
forest plot, filled = p<0.05). Partial-data preview already shows the intended read: CA1–V1 IFI rises
over trial_frac in **5/5** animals (signed-rank floored at p=0.062) but the LMM resolves it
(p=0.0039) and OLS over-inflates — yet it is **null vs lp_rel** (p=0.98) → time-on-task, not
learning-locked.

**✓ PRIMARY DONE (stage 1, bin25 FS-excl, committed 91ad829).** 16 animals. Findings:
**CA1→V1 IFI rises with experience** is the one directionality slope supported at the honest unit
(animals $p=0.012$–$0.098$ across all 3 axes, cluster-robust LMM $p=0.021$); **held-out CC slope is
null at both animals and LMM** for every pair/axis (only naive OLS inflates, to $p\sim10^{-20}$) →
strength-null survives the powerful dims unit. Written into report **§3.2** (forest figure
`HCV1_CCA_fsexcl_trajdims_bin25.png` + paragraph) and **§3.7** (CC-null-under-LMM); §4 synthesis
directionality row + re-run status box updated. fsexcl window figures regenerated from the FINAL
re-run output. **Bug fixed:** driver writes `..._dims_fsincl.csv` (`_dims` before suffix); analysis
`load()` was looking for `..._fsincl_dims.csv` → would have missed every robustness stage. Fixed.

**NEXT when bxdq9o98x fires (stages 2–4 done — robustness only).**
1. `analyze_trajectory_dims.py fsincl` / `bin50` / `bin50 fsincl` + `figs_trajectory_dims.py` same args → confirm the CA1→V1 + strength-null picture holds; add a one-line robustness note to §3.2/§3.7 (don't expect conclusion change).
2. Regenerate fsincl window figures from final data (`figs_report.py trajectory_w15_bin25` once fsincl CSV shows 16 animals); bin50 figures as robustness if wanted.
3. Commit the fsincl/bin50 `_dims*.csv`; update this log + STATE.md.
**Deferred (low priority):** figure numbers + panel letters (#9, do once figure set is final); vanilla-CCA (#2) / no-PCA (#3b) runs.
**Deferred (low priority / CPU-contended):** figure numbers + panel letters (#9); vanilla-CCA (#2,
`--no-partial`) and no-PCA (#3b, `--k`) robustness runs (driver flags ready — start only after the
trajectory queue frees up).

## ✓ DONE (2026-06-07) — NONLINEAR (kernel) CCA full suite + report §5

Complete **kernel-CCA analogue of the whole pipeline** built, run (all 16 animals, FS inc/exc,
learners/non via session+epoch scopes, cued/uncued), analysed, figured, and written into report **§5
"Nonlinear (kernel) CCA"** (4 figures: vs_linear/levels/epochs/transition × FS). 241 tests pass.
**VERDICT: the communication subspace is LARGELY LINEAR.** Kernel CCA edges linear by only a small,
consistent margin — median Δ(KCCA−linear) held-out CC = +0.016/+0.014 (FS-excl/incl), KCCA>linear in
66%/58% of cells — significant per-pair only in **V1-RSC** (FS-incl p=0.02) and **CA1-V1** (FS-excl
p=0.033); for the strongest subspaces (CA3-DG ~0.34/0.45, CA1-CA3) **KCCA≈linear**. Strength flat
naive→expert (mirrors linear); structure-coeff **Gini de-sparsifies naive→expert in CA1-V1** (FS-incl
LMM 3e-5) — the §3.3 de-sparsification echoes in the NONLINEAR membership (leading pair shifts
CA1-RSC→CA1-V1). Cued: CA1-V1 cc↑ (p=0.047), V1-RSC Gini↓ (p=0.031). No nonlinear effect overturns a
linear conclusion → modest cortically-localised nonlinear add-on, not missed dominant structure.
**Compute choices (in report):** fixed ridge reg=10 (optimism-free vs linear), n cap 900, 10 shuffles,
±4 lags, structure-coeff membership; trajectory deferred. NOTHING committed yet.
- **GOTCHA fixed this build:** KCCA transition cued cell returned NaN held-out CC at CAP=900 because
  the sample-matched cued block spanned <N_FOLDS whole trials → switched BOTH conditions to
  position-based CV folds (run_kcca_transition.py). KCCA `kcca_fit` sped up ~10× via Cholesky+ARPACK
  top-d (non-symmetric dense eig was the bottleneck; 10-min test → 25 s).
- **Review workflow done (16 findings, 11 confirmed):** fixed _heldout to pair KCCA+linear per fold
  (superiority Δ over SAME folds — verified output IDENTICAL on animal 52, no re-run needed);
  KCCAResult.r docstring (it's the ridge-penalised value, under-estimates corr; used only to scale
  beta). Report §5 fixed: per-pair superiority "neither FS-robust"; surrogate-clearance 86/82% FS-excl
  vs 100% FS-incl; "Uncued→cued" label; single-split surrogate/IFI disclosed. IFI single-split + epoch
  LMM-vs-Wilcoxon kept as documented caveats. 254 tests pass.

### (prior in-flight, now done) NONLINEAR build details
- **New code (TDD):** `src/tom_cca/kernel_cca.py` (regularised RBF KCCA — Hardoon ridge form;
  Cholesky+ARPACK top-d for speed; `kcca_fit/score/variates`, `median_gamma`, `kcca_lagged_heldout`
  for direction); `src/tom_cca/kcca_window.py` (full per-cell suite: held-out KCCA **and** linear CC
  on same folds = superiority; shift-surrogate n_sig; held-out lagged IFI; structure-coefficient
  Gini); `membership.variate_structure_coefficients` (method-agnostic membership — KCCA has no neuron
  weights). Drivers/analyze/figs as above.
- **Compute choices (state in report):** fixed ridge reg=10 (prototype sweet spot; fair vs linear, no
  selection optimism); n capped 900/cell; 10 shuffles; lags ±4; PCA→K then KCCA. Sliding-window
  trajectory DEFERRED (≈10× cost; epochs give the learning axis).
- **Prototype verdict (prototype_kcca.py, pre-build):** KCCA in-sample SATURATES at our n (gap ~0.40,
  0/36 held-out ≥0.99 → held-out mandatory); clears its surrogate in 19/36 (intra-hippocampal pairs);
  KCCA≈linear at reg=1 (Δ=−0.02), modest edge at reg=10 in some pairs → **leaning: subspace largely
  linear**, full suite quantifies it.

## ✓ DONE (2026-06-07) — Pearson control + leak-free CV (verified)

- **Pearson control** (`membership.pearson_coupling_scores` → `gini_pearson_x/y`): the CCA-free Gini
  de-sparsifies in the **same (negative) direction** as the weight-Gini (CA1-RSC trial_frac med −0.053
  vs −0.132) but **attenuated / mostly n.s.** (W 0.11, LMM 0.056) — corroborates directionally, not
  decisively (unit-count noise floor; slope-only). Headline = NOT a pure CCA-weight artifact.
- **Leak-free held-out CC** (`window_subspace(...,Z=...)` per-fold residualise+PCA+CCA): re-ran
  trajectory FS-excl/incl + transition + epochs FS-excl/incl. **`gini_x` byte-identical** (trajectory
  max|Δ|=0.0 on 636 windows; epochs 0/148) → de-sparsification headline + IFI direction UNCHANGED;
  only cc1/n_sig/mi_sig shifted (more conservative, still null). All 22+ figures regenerated.
- **Code-review fixes:** lmm_slope small-N guard; wilcoxon zero-handling; legacy `partial_cca_cv`
  caveat; dead-field removal; misc. NOTHING committed yet.

## ✓ RE-RUN COMPLETE (shuffle 20→100) — 2026-06-06

All analyses now at `config.SURROGATE_SHUFFLES = 100`. Trajectory (FS-excl 636 rows, FS-incl 647),
transition (58 rows), and epochs were all re-run/regenerated @100; all 22 vault figures refreshed
(`attachments/HCV1_CCA_*`, 2026-06-06 21:54). **No conclusion changed** — confirmed identical:
CA1-RSC Gini↓ LMM β=-0.159 p=4.25e-05; CA1→DG IFI↑ uncued→cued d=+0.019 p=0.016 (learners); rotation
at/below split-half floor. As expected: `gini`, `cc1`, `ifi`, `optimal_lag`, `rot/sh` are
shuffle-independent; only `n_sig`/`mi_sig`/dims-pools could move, and they remain null.
**Sig threshold decision (this session):** kept the 95th-pct circular-shift threshold (Han uses
mean+3sd ≈99.9th pct); divergence documented in report §2.6. Nothing committed yet.

## CURRENT STATE (2026-06-06)

**Branch:** `main` — the `cca-consolidation` arc was merged (commit `5f75022`) and pushed; working
tree clean, up to date with `origin/main`. All linear + nonlinear CCA work is committed. (Updated
2026-06-10; the prose below this line dates to 2026-06-06 — see the ✓ DONE blocks above and the newest
LOG ENTRY for what followed.)
**Tests:** 254 passing (`cd cca && PYTHONPATH=src python -m pytest -q`).

**Where the project is, in one paragraph.** The original epoch/landmark sweep was statistically
dead-on-arrival for the *learning* question: ~2,000 samples/fit with `k≈n/15` overfit (held-out
CC→0.999 in high-k configs), and the cross-animal contrast is intrinsically underpowered (1 session
/animal, 12 learners + 4 non-learners, 4–10 animals/pair). We diagnosed this, read the two
reference papers our pipeline descends from, and **reframed** onto their regime: fit pCCA on
**continuous running data** (25 ms, ~100k+ samples/session, PCA→30, samples/k ≫ 50 — no overfitting,
validated), **control the third area** (pCCA primary), and ask the questions the data can support —
**within-animal learning trajectory** and **structure** (dimensionality, direction, membership,
rotation) rather than a noisy 3-epoch magnitude contrast.

**What is built & committed (this arc):**
- `src/tom_cca/paired_stats.py` — Wilcoxon + BH-FDR (shared).
- `src/tom_cca/mixed_effects.py` — random-slope LMM + per-animal collapse (honest learning tests).
- `src/tom_cca/subspace_stats.epoch_subspace_stats` — spatial (lag-0) analogue of the landmark cell stats.
- `scripts/learning_changes_spatial.py`, `scripts/learning_changes_mixed.py` — honest learning tests.
- `scripts/prototype_continuous_pcca.py` — **validated the continuous regime** (median held-out CC
  0.18, 0/36 saturated; CA3-DG strongest up to 0.56; CA1-RSC weak 0.04–0.17).
- `src/tom_cca/trajectory.py` — sliding windows, CV pCCA held-out CC, linear slope (Frame A core).
- `src/tom_cca/subspace_window.py` — **FULL per-window readout**: held-out CC per dim, n_sig
  (circular-shift surrogate), mi_sig, IFI, optimal lag, Gini, canonical weights + member masks.
- `scripts/run_trajectory.py` — Frame A driver: full suite per window over trials, 3 learning axes
  (trial-fraction / performance / LP-relative), cross-window rotation + Jaccard, learner flag →
  writes rich `results/trajectory_windows.csv`.
- `scripts/analyze_trajectory.py` — fast slopes/sign-tests/levels from that CSV (no re-fit).
- `scripts/run_transition.py` — uncued→cued (world 4→3) full-suite comparison + cross-condition
  subspace rotation angle, sample-matched, learner-split.

**Report (external):** `ResearchVault/Projects/Hippocampus-V1/Hippocampus-V1 Communication-Subspace
Learning Report.md` — full methods (LaTeX) + 4 embedded figures (`attachments/HCV1_CCA_*.png`,
generated by `scripts/figs_report.py`), linked from the vault project hub; vault `log.md` updated.
Epistemic state: Contested.

**Running / pending right now:** nothing running. Both full-suite analyses DONE.
- Transition (uncued→cued, 13 animals, 58 rows, `results/transition_uncued_cued.csv`): NO abrupt
  strength jump (cc1/mi_sig/n_sig deltas n.s.). **CA1→DG directionality increases uncued→cued**
  (d_ifi +0.019, p=0.016 learners) — same direction as the trajectory's CA1→DG-IFI-rises-with-
  learning. Subspace rotation uncued→cued is ~75–88° for most pairs BUT that ≈ the trajectory's
  window-to-window rotation (so ~80° is the noise floor, NOT reorientation); **CA3-DG is the
  exception (41–51°) = genuinely stable across the transition**, consistent with it being the
  strongest/richest subspace. CAVEAT: rotation needs a within-condition split-half floor to claim
  reorientation (next step); uncued phase short (K=15, n=3–8/pair).
- **Speed/correctness fixes applied (2026-06-06):** window_subspace was ~11 s/window (significance
  + 21-lag scan on ~41k samples). Now cap each window to a contiguous ~6 k-bin block (grown to span
  ≥ N_FOLDS+1 trials so the CV is valid) → ~2.3 s/window. Fixed `n_sig` overcounting: significance
  now compares the *held-out* per-dim CC to the *dominant*-dim circular-shift threshold (was an
  in-sample per-dim test → 23 sig dims; now ~3). Driver writes CSV incrementally + unbuffered.
  CAVEAT: lags are in running-bin units (non-running bins removed), so IFI/optimal-lag are
  approximate (≈±250 ms); a segment-aware lag is a future refinement.

**Findings so far (honest):**
- The pooled landmark "CA3-DG strengthens" headline was **pseudoreplication** (n = animals ×
  landmarks). Honest per-animal / mixed-effects tests are ~null for magnitude. See `STATE.md` §3.
- Continuous-regime pCCA levels: intra-hippocampal (CA3-DG cc1≈0.36, n_sig≈2.9; CA1-CA3 0.31) >
  hippocampal-cortical (CA1-RSC 0.09, CA1-V1 0.12). CA3-DG strongest/richest subspace.
- **FINAL VERDICT after the full interrogation (the one robust signal, and what it is NOT):**
  - **The ONE robust effect: the communication subspace DE-SPARSIFIES over the session** (Gini↓ =
    more neurons recruited), CA1-RSC & CA1-DG, **early (naive→intermediate), plateaus post-LP**.
    All tests agree: Wilcoxon/t/LMM; trajectory LMM p=4e-5; epoch naive→int LMM p~1e-5; FS-invariant.
  - **But it is most parsimoniously EXPERIENCE / time-on-task, NOT learning-specific** —
    non-learners de-sparsify comparably; the trial_frac×learner interaction is n.s. for every pair
    (p=0.26-0.97); the LP-plateau is suggestive of a learning component but unproven (n=4 non-learners,
    pre-LP windows n=2). See 2026-06-06 learning-vs-time entry.
  - **Everything else is a NULL or an artifact:** coupling STRENGTH flat (cc1/mi_sig/n_sig; flat even
    under dims-as-n → genuine, not power); DIRECTION null (IFI *and* optimal lag null at animals-as-n;
    the dims-as-n "hits" CA1-RSC/V1-RSC are pseudoreplication — e.g. V1-RSC nai→exp lag Δ=0 p=1.0 by
    animal vs p=0.007 by dims); ROTATION at the split-half noise floor (no reorientation).
  - **dims-as-n (Buzsáki unit) lesson:** inflates N ~5-15× and manufactures significant strength/
    direction results that vanish at the animal level — included for comparison, not inference.
  - Caveats: small N (4-10/pair); no cross-pair MC correction (rely on cross-axis/cross-test
    consistency); IFI from in-sample lag curves (optimal-lag agrees); lags running-bin-approximate.
- Frame A v1 (CC1 only): within-animal trajectories clean (|r| up to 0.94) but slope sign
  animal-specific → magnitude null. Resolved by reading structure instead of magnitude.

**Data facts (Tom cohort):**
- 16 animals, **1 session each**, ~100–320 cued trials. Learners (have LP): 28,34,36,41,52,61,63,
  66,68,73,75,77. Non-learners (no LP): 70,71,98,100.
- Worlds: 1 = darkness/ITI, **3 = cued tunnel (task)**, **4 = uncued corridor** (pre-task, 3–10
  trials; present in 13/16 — not 28/34/36).
- Export has MORE than we load: `units/depth`, `units/isi/histogram`, full `waveforms`,
  cued/uncued firing — **deferred by user decision; FS-vs-regular (`idx_fs`) only for now.**
- Continuous loading reuses the temporal-arm path (`dataio.area_activity_50ms`, `_load_temporal_streams`).
- Per-trial performance = `analysis_behaviour/lick_ratio` (1-based trial idx); LP on `entries[id].lp`.

**Key conventions adopted (from the papers):** PCA→fixed ~30 comps then pCCA; ≥50 samples/variable;
control the third area; residual = subtract condition-triggered mean; **session/animal as unit**;
subsample to match neuron counts for cross-pair comparisons; surrogate everything; drop saturated
(CC≥0.99) windows.

**NEXT (when trajectory lands):**
1. `analyze_trajectory.py` → report levels (IFI direction, n_sig, Gini, rotation) + slopes vs the
   3 axes, learners vs non-learners. 2. Run `run_transition.py`, report. 3. Commit the full-suite
   drivers once validated on real data. 4. Consider: subsample-to-match-N control; split-half angle
   noise floor for rotation; (deferred) load depth/ISI for membership×properties.

---

## LOG ENTRIES (newest first)

### 2026-06-12 (cont.) — Batch DONE: IFI-window result, KCCA upgrade resolved, dims-as-n (transition), spacious figures
- **Batch finished** (all 8: trajectory W15 ×4 + KCCA upgrade ×4). Cores free.
- **IFI-window sweep (held-out, segment-aware) — CA1→RSC clean directionality:** IFI +0.079 at ±50 ms,
  t=5.0 **p=4e-4** (animals-as-n; Bonferroni/12 windows; a-priori CA1-leads-RSC). Upgrades §3.2 from "weak"
  to a clean cortico-hippocampal flow. CA1-SUB ±25 ms (SUB→CA1, p=0.009), RSC-SUB ±125 ms. **Tight windows
  win** (median ~±75 ms). Report §3.2 + header + status box; spacious figure.
- **KCCA upgrade (30 shuffles, ±8 lags) — largely-linear HOLDS:** +0.016/66% (bin25, identical to old);
  +0.001/51% (bin50). `kcca_metrics.csv` refreshed to the upgraded run; §5 ⚠ resolved.
- **Transition dims-as-n DONE (cued/uncued):** `run_transition` exports `transition_dims.csv` + both-units
  tables. §3.5: animal-level effects (CA1-V1 ΔCC 7/7, CA1-DG ΔIFI) are **dominant-dim**; dims-as-n adds
  **inflation** (CA1-CA3/DG/SUB ΔCC sig only pooled, animal-level null/opposite). Honest power-check framing.
- **Trajectory per-dim export added** (`run_trajectory` → `trajectory_w15_*_dims.csv`); **W15 re-run RUNNING
  overnight** (bin25 both FS first, then bin50) for trajectory dims-as-n + backfills bin25 CC1-rotation.
  **When it lands:** regenerate figs from the new bin25 (cc1+dims), add the trajectory both-units view to §3.2/§3.3.
- **Figures now breathe (#figs):** `figstyle.grid()` 4×2 spacious + bigger fonts/padding; applied to
  early_trials, ifi_windows, `_grid_traj` (gini/ifi_traj). Regenerated spacious. **Pending:** user confirm on
  the look (sent an example); figs_units + the 1×5 levels still to convert.

### 2026-06-12 — Meeting prep: W15 trajectory re-analysis (both bins) + CC1-rotation + parametric stats
- **Batch:** trajectory W15 DONE for bin25 (both FS) + bin50 FS-excl; bin50 FS-incl still running; **KCCA
  upgrade NOT yet started** (4 runs pending, ~hrs each).
- **De-sparsification ROBUST + SHARPER at W15/both bins:** CA1-RSC trial_frac **LMM 3.3e-7 (bin25)** /
  2.4e-5 (bin50) — stronger than W30's 4e-5; CA1-DG holds (LMM 0.005/0.013); **CA1-V1 weakens to
  borderline** (W=0.084) → dropped from the firm headline, kept suggestive. Two-bin agreement kills a
  bin-width artefact. Report §3.3 + §2.10 updated.
- **CC1-only rotation (#14) sharpens §3.4:** CC1 split-half floor **16–62°** (vs top-3 ~80°); cross-window
  CC1 rotation AT/BELOW floor for every pair — **sig BELOW** for CA1-RSC (34 vs 61°, p=0.008), CA1-DG
  (24/43, 0.008), CA1-V1 (40/62, 0.006), V1-RSC (0.031); CA3-DG at floor (11/16). → no reorientation
  confirmed at the *dominant* dim, not just a top-3 noise artefact. From bin50 (carries the cc1 cols).
- **Parametric stats (#7) DONE for early-trials:** `analyze_early_trials` now reports paired-t (primary)
  + Wilcoxon + LMM-trend. Honest: CA1-RSC IFI trial-4 t=0.054 (vs W=0.034, fragile); LMM trends null.
- **Figures:** regenerated trajectory figs from W15 bin25 (`figs_report.py <stem>` arg); all carry the
  figstyle (svg+png ≤1600px, despined, prominent red stars).
- **Report:** added a prominent **re-run status box** (final vs pending) up top for the meeting.
- **PENDING (placeholders):** KCCA upgrade (§5 figs still old 10-shuf/±4), IFI-window sweep (§2.7/§3.2),
  bin25 cc1-rotation re-run, bin50-FS-incl trajectory (running). **Do NOT commit `trajectory_w15_bin50_fsincl.csv`
  yet — still being written.**

### 2026-06-11 — Review response (18 comments) — Batch 1 correctness fixes
User reviewed the report + analysis (18 comments). Decisions: learners stay PRIMARY, pooled-all-16
added as a secondary power check (not a rewrite); trajectory window 30→**15 (step 5)**; KCCA upgrade
to **30 shuffles, ±8 lags** (moderate). Ran the `stats-rigor` checklist. **Batch 1 done (verified +
fixed; report edits in the vault, code-comment fixes here):**
- **Bin width = 25 ms for the MAIN pipeline (trajectory/epochs/transition/KCCA), 50 ms only for
  early-trials §3.8.** SUBTLE: `config.DEFAULT.temporal_bin_ms=50`, BUT `run_trajectory.py`/`run_kcca.py`
  set their own `--bin-ms` DEFAULT=25 → committed §2-5 results are 25 ms; only `run_early_trials.py`
  (inherits config.DEFAULT) ran at 50 ms. ⚠️ My first read this session wrongly concluded "50 ms
  everywhere" and the report was briefly edited to 50 ms — **REVERTED to 25 ms** with a standardisation
  note (§2.1). **Now re-running EVERYTHING at BOTH 25 & 50 ms** (`run_review_batch.sh`) to make them
  consistent and get the 25 ms lag-resolution benefit the directionality readouts want.
- **BATCH LAUNCHED (detached, ~hrs):** `run_review_batch.sh` → trajectory **window=15 step=5** + KCCA
  **upgraded (30 shuffles, ±8 lags)**, each at 25 & 50 ms × FS-excl/incl, to distinct files
  (`trajectory_w15_bin{25,50}{_fsincl}.csv`, `kcca_up_bin{25,50}{_fsincl}.csv`); log `results/review_batch.log`.
  New `lagged.heldout_lag_curve_flat` (held-out, segment-aware lag curve — fixes #12 in-sample-lag + #6
  concatenation) + `ifi_by_window` give the "which integration window is cleanest" sweep (+3 tests; 262 pass).
- **PROGRESS (later 2026-06-11, no-compute items done while batch runs):**
  - **Conceptual report write-ups DONE** (vault §2.5 PCA-in-CV leak-free, §2.6 MI=−½Σlog(1−ρ²) derivation,
    §2.7 concatenation soundness + segment-aware fix, §2.10 test-choice justification + n=6 floor, §3.1
    coupling-hierarchy quantification: intra>cortico within-animal Δ+0.178, 9/10, t p=0.004).
  - **IFI-window suite BUILT + smoke-validated** (`run_ifi_windows.py`/`analyze_ifi_windows.py`/
    `figs_ifi_windows.py`): per (animal,pair) session held-out segment-aware lag curve → `ifi_by_window`
    → per-pair "cleanest window" = max across-animal |mean/SEM|. Smoke (animal 28): CA1-DG peak +1 bin.
    **Run full at both bins AFTER the batch frees cores** (it refits CCA per lag×fold on the whole session).
  - **W15 headline check (bin25 FS-excl done):** CA1-RSC Gini slope −0.121 (p=0.008) vs W30 −0.132 (p=0.008)
    — de-sparsification **robust to window size**; W15 gives ~25 windows/animal vs 8.
  - **CC1-only rotation (#14) DONE in code** — `subspace_window` returns the CC1 split-half floor (reuses the
    same half-fits, no extra cost); `run_trajectory` writes `rot_x_cc1`/`sh_x_cc1`; `analyze_trajectory` prints
    both top-3 and CC1 rotation-vs-floor. **Motivation confirmed:** top-3 floor ~85° (noisy → the §3.4 "rotation
    = noise" result) but CC1 floor ~6° — a CC1-only rotation is far more sensitive, may revise §3.4. The
    in-progress **bin50** trajectory runs pick up the cc1 columns; **bin25 (done/running) lacks them → needs a
    re-run** for cc1 at 25 ms (or use bin50; rotation-vs-floor is bin-robust).
  - **Pooled-16 (#16) DONE** (added to report §3.3): pooling all 16 sharpens the de-sparsification (CA1-RSC
    p=0.008→0.001, CA1-V1 0.049→0.021, V1-RSC emerges) — confirms experience-driven.
  - **Figure clarity overhaul DONE (code, evergreen):** new `scripts/figstyle.py` — `apply()` sets a clean
    house style (constrained layout → no label collisions; de-spined, no gridlines → declutter) and `save()`
    writes the **`.svg`+`.png` pair with PNG capped ≤1600 px** (the repo standard — **0 SVGs existed before;
    now 55**, all PNGs ≤1600 px) + a prominent offset/bold/red `star()`. Applied across
    figs_report/kcca/units/early_trials/ifi_windows + analyze_epochs; regenerated from committed data.
    Re-applies automatically when the figures regenerate from the batch.
  - **STILL TODO when batch lands:** run IFI-window full (both bins); re-analyse trajectory(W15)+KCCA both bins
    + refresh figs (incl. the new CC1-rotation §3.4); bin25 cc1-rotation re-run; Batch-4 figure overhaul
    (numbering/panels, extra axes, heatmaps mean/parametric, dims-as-n §3.6/3.7); vanilla-CCA (#2), no-PCA (#3b).
- **V1-RSC is NOT pseudoreplication (user right).** IFI naive→expert +0.030→+0.013 (animals, n=6) vs
  +0.027→+0.014 (dims) — SAME shape/magnitude, 5/6 animals down; signed-rank p=0.16 is the n=6 FLOOR
  (0.031), not a null; the report's "Δ=0 p=1.0" was a median-of-integer-lag artefact (mean lag falls
  +1.36→+0.23). Reframed §3.2 + synthesis row + headline: directionality is **weak & underpowered**,
  not a manufactured artefact. CA1-RSC intermediate IFI remains a genuine false-positive (animal t=0.44).
- **Transition has a real effect the report missed (user right).** CA1-V1 ΔCC **7/7 animals +, +0.117,
  signed-rank p=0.016** (report said strength deltas n.s.); CA1-DG ΔIFI 7/7, t=0.012. Reframed §3.5;
  noted BH-FDR-marginal across the 8×7 grid → pre-specified-pair effects; LMM planned.
- **Confirmed PCA is leak-free in CV** (subspace_window fits residualise+PCA+CCA train-only per fold).
- **REMAINING:** Batch 1 conceptual (MI-eqn derivation, concat soundness, coupling-hierarchy quant,
  why-signed-rank); Batch 2 pooled-16 secondary; Batch 3 re-runs (window=15, held-out lag, vanilla CCA,
  no-PCA, CC1-only rotation, optional 25 ms); Batch 4 figures (numbering+panels, extra axes, dims-as-n in
  §3.6/3.7, heatmaps mean/parametric); Batch 5 KCCA upgrade.

### 2026-06-10 — Very-early-trial (pre-10) analysis — new module + driver + report §3.8
- **Question (user):** do any metrics (CC/KCCA/IFI/Gini/angles) change in the *first ~10 trials* then
  plateau? Tested trial 1 vs 4/7/10. **Key constraint:** a single trial ≈ 580 running bins → too few
  to re-fit a 30-comp pCCA (~1500 needed). So a **hybrid** (user-approved): (A) **projected** per-trial
  readout — fit subspace once on later trials, project trials 1/4/7/10 leak-free (CC, IFI, participation-
  Gini); (B) **block** refit on cumulative first-5/7/10 vs late for the fit-only metrics (held-out CC,
  n_sig, MI, weight-Gini, angle-vs-late, KCCA).
- **New code (TDD):** `src/tom_cca/early_trials.py` (`reference_fit`, `trial_cc/ifi/participation_gini`,
  `early_trial_metrics`) + `tests/test_early_trials.py` (5 ground-truth tests: held-out coupling
  detection, sparse>dense participation Gini, known lead/lag sign, leak-free fit). Driver
  `scripts/run_early_trials.py` (projected + block + KCCA, `--no-kcca`/`--include-fs`); analysis
  `scripts/analyze_early_trials.py` (1-vs-4/7/10 & block-vs-late Wilcoxon + plateau flag + angle-vs-
  split-half-floor); figures `scripts/figs_early_trials.py`. **259 tests pass.**
- **FINDINGS (fsexcl, 16 animals):**
  1. **No early strength jump** — per-trial CC & block held-out CC flat for every pair.
  2. **De-sparsification is SLOW / post-trial-10** — block weight-Gini significantly *higher* (sparser)
     in first-5/7/10 than late for CA1-DG/CA1-V1/RSC-SUB/V1-RSC(Y)/CA1-SUB, but the margin is ~constant
     across the early blocks (doesn't shrink) → within the first 10 trials it stays uniformly sparse;
     broadening happens after ~trial 10. Refines §3.6 ("naive→intermediate" spans past trial 10).
  3. **One fast (trial 1→4 then plateau) effect, cortical: CA1-RSC** — projected IFI Δ+0.254 p=0.034,
     participation-Gini(Y) Δ+0.013 p=0.027, both plateau by trial 7.
  4. **V1-RSC early subspace reorientation** — early block rotated from late ABOVE the split-half floor
     (angle−floor +7.7/+19.3/+14.6° at t1-5/7/10, all p≤0.012); every other pair at floor (cf §3.4).
- **Report:** added **§3.8** (3 figs: proj_cc1, block_gini_x, proj_ifi), methods **§2.10c**, a §4
  synthesis row, §6 [DONE] bullet, §7 provenance. Bumped report `last_enriched`. Vault `log.md` entry.
- **FS-incl DONE (16 animals):** all `early_trials_*_fsincl.csv` written, figs regenerated, §3.8
  FS-robustness paragraph + verdict + synthesis row added. **FS picture:** de-sparsification (block
  weight-Gini sparser early than late: CA1-RSC/DG/V1/SUB, RSC-SUB, V1-RSC) and flat strength are
  **FS-ROBUST**; the two FAST cortical effects are **NOT** FS-robust — CA1-RSC trial-1→4 IFI vanishes
  (p=0.91 vs 0.034), V1-RSC reorientation attenuates (only t1-10 survives FS-in). §3.8 now states the
  fast effects are FS-excluded-only/fragile. **GOTCHA:** don't wrap `python … &`/nohup inside a
  `run_in_background` Bash call — it double-backgrounds and the harness "completes" on the wrapper, not
  the real run (poll instead). Repo uncommitted (on `main`, per new no-branch policy in CLAUDE.md).

### 2026-06-10 — KCCA figures overhauled (per-animal dots + stats) + Obsidian rendering fix
- **`scripts/figs_kcca.py` rewritten** so all four KCCA figures overlay **per-animal dots** on the
  mean±SEM bars and carry a **significance annotation**, matching the linear `figs_report.py` house
  style (ported `_bar_points_sem`; added grouped `_grouped_bars_points`):
  - **vs_linear** — grouped KCCA-vs-linear bars + dots + paired-Wilcoxon superiority star (CA1-V1 *).
  - **levels** — CC₁ / IFI / Gini with dots; IFI carries a sign-test-vs-0 star (null at animal level).
  - **epochs** — naive/int/expert grouped bars + dots + **expert−naive Wilcoxon** star (CA1-V1 Gini *,
    W=0.027, matches §5.2); **stacked 3×1 vertically** so it stays legible at Obsidian note width.
  - **transition** — Δ(cued−uncued) bars + dots + signed-rank-vs-0 star (CA1-V1 ΔCC *, V1-RSC ΔGini *).
  - Regenerated all 8 PNGs (fsexcl+fsincl) into the vault attachments; stars cross-checked vs
    `analyze_kcca.py` p-values.
- **Obsidian rendering fix:** §5.2 was a bulleted list interrupted by `![[…]]` embeds (plus two
  stacked embeds with no blank line) → garbled breaks. Converted §5.2 to plain paragraphs (bold
  lead-ins, §3 style) with every embed blank-line-isolated; separated the §3.3 control embed from the
  following text. Verified all 6 edited embeds are blank-line isolated.
- 254 tests green (no `src/` change); ruff not installed in this env (script follows existing style).

### 2026-06-10 — Vault report fixes (KCCA figures embedded; Pearson control written up)
- **Vault report only — no code/results/figures changed.** Fixed two gaps in `ResearchVault/Projects/
  Hippocampus-V1/Hippocampus-V1 Communication-Subspace Learning Report.md`:
  1. **§5 KCCA figures were bare `[[…]]` wikilinks** (buried in parentheticals) → they rendered as text
     chips, not images. Converted all five to own-line `![[…]]` embeds (vs_linear, fsexcl/fsincl
     levels, epochs, transition).
  2. **Pearson (CCA-free) control was missing from §3.3** despite being built/committed/figured.
     Added `HCV1_CCA_fsexcl_gini_control.png` + a paragraph: the CCA-independent coupling-Gini
     de-sparsifies in the same (−) direction (CA1-RSC trial_frac med −0.053 vs weight-Gini −0.13) but
     attenuated/borderline (W 0.11, t 0.13, LMM β=−0.056 p=0.056) → headline is NOT a pure CCA-weight
     artifact. Numbers from `analyze_trajectory.py` (`gini_pearson_x`).
- Bumped report `last_enriched` → 2026-06-10; added a vault `log.md` audit entry.
- Repo unchanged: working tree clean, 254 tests green.
- **Still open (flagged, not actioned — user's call):** `git rm` the stray tracked
  `test_delete_check.txt` (0-byte, repo root); tick the now-answered "Kernel CCA? what kernel? what
  scale?" item in the vault `Hippocampus-V1-Tasks.md` (answered by §5: RBF, median-heuristic
  bandwidth, ridge=10); the §6 deferred items (neuron-count matching, depth/ISI membership,
  learning-vs-time with more non-learners).

### 2026-06-06 — Methods-note audit (vault CCA report); threshold decision; 20→100 re-run launched
- Read the new vault note `Methods/CCA in Systems Neuroscience — Research Report.md` (full-text audit
  of our exact lineage: Han&Helmchen = template, Gonzalez/Buzsáki = pCCA antecedent) against our code.
  **Verdict: pipeline already aligned on ~everything** (PCA→30, n≫p continuous, pCCA third-area control,
  residual-on-condition-mean, circular-shift null, held-out CV + split-half floor, session-as-unit,
  IFI/optimal-lag, low-dim 1–3 sig answer). **No reframe warranted** — the note vindicates the path.
- **One concrete discrepancy found & decided:** our sig threshold is the **95th pct** of the dominant
  circular-shifted CC (`subspace_window._significance`, `np.quantile(null_top, 1-alpha)`), whereas
  Han uses **mean+3sd** (≈99.9th pct) — much stricter. Likely why we count ~3 sig dims vs Han's 1–2.
  **DECISION (user): keep 95th pct, document divergence.** Only moves the already-null
  n_sig/mi_sig/dims-as-n metrics; Gini headline + all threshold-independent conclusions untouched.
  Documented in report §2.6.
- Three other items are caveat/optional only, NO reframe: (2) optional Pearson-correlation control for
  the Gini↓ headline (Han step 6, not done — cheap robustness add if desired); (3) temporal-autocorr/
  effective-N caveat — our circular-shift is the *recommended* mitigation, just state it explicitly;
  (4) weight-interpretation caveat — frame Gini as a distributional claim, don't over-read membership.
  Kornblith p≥n / CKA degeneracy does NOT bite us (firmly n≫p) — reassurance only.
- **Re-run COMPLETE (background, parallel, all exit 0):** run_trajectory FS-excl (636 rows, 21:50) +
  FS-incl (647, 21:53) + run_transition (58 rows, 21:09), all @ SURROGATE_SHUFFLES=100. Regenerated
  all 22 vault figures (figs_report + figs_units, 21:54). **CONFIRMED: no conclusion changed** —
  CA1-RSC Gini↓ LMM p=4.25e-05 (identical), CA1-RSC IFI p=0.0078, CA1→DG IFI↑ uncued→cued p=0.016
  (learners), rotation ≤ split-half floor. n_sig/mi_sig remain null. See the RE-RUN COMPLETE block at top.

### 2026-06-06 — Shuffle count centralised (config.SURROGATE_SHUFFLES=100); epochs re-run; re-run pending
- All drivers now read `config.SURROGATE_SHUFFLES` (=100) — no more per-driver drift. Epoch analyses
  re-run @100 (epoch_metrics/epoch_dims/figs regenerated); **Gini result unchanged** (shuffle-independent).
  dims-as-n V1-RSC "hits" persist even capped to ≤3 dims/animal (animal-level still null).
- **TRAJECTORY + TRANSITION CSVs still reflect 20 shuffles** → re-run pending (only n_sig/mi_sig/dims-pools
  affected; Gini/CC/IFI/lag/rotation + the headline are shuffle-independent). See the "RE-RUN PENDING"
  block at top for commands. Documented as to-be-continued.

### 2026-06-06 — Significance shuffles 20→100; dims-as-n capped to 3/animal; §3.6 clarified
- Significance of canonical dims = circular-shift surrogate (roll one area, refit pCCA, record
  DOMINANT shuffled CC), threshold = 95th pct, held-out CC must exceed it. **Bumped N_SHUFFLES 20→100**
  in run_epochs (re-run) — 20 was too few for a stable 95th pct. Documented in report §2.6.
- **dims-as-n now capped to top-3 sig dims/animal/epoch** (by held-out CC) in figs_units + analyze_ifi
  + analyze_dims_as_n, to curb the per-animal imbalance.
- §3.6 (epoch timing) is Gini = POPULATION metric → animals-as-n only (no per-dim analogue). The
  per-dim metrics from the same epoch fits (CC, IFI, lag) ARE shown both-units (capped) in §3.2/§3.7
  (figs_units). Clarified in report.

### 2026-06-06 — Animals-as-n vs dims-as-n comparison figures
- `figs_units.py`: per per-dim metric (held-out CC / IFI / optimal lag) a 2-row figure — ANIMALS-as-n
  (top) vs DIMS-as-n (bottom) × 8 pairs across epochs, points + mean±SEM, vs-0 stars, naive→expert p
  per panel. Embedded in report §3.2 (ifi/lag) + §3.7 (cc). Population metrics (Gini/n_sig/rotation)
  have no per-dim analogue → animals-as-n only (epochs figure). Visualises the pseudoreplication.

### 2026-06-06 — Learning vs time-on-task: de-sparsification is EXPERIENCE-driven (learning unproven)
- `analyze_learning_vs_time.py`: (1) learner vs non-learner Gini-vs-trial_frac slopes — comparable
  (CA1-RSC -0.13 vs -0.10; between-group n.s.); (2) interaction LMM gini~trial_frac*learner+(trial_frac|animal)
  — trial_frac×learner n.s. for ALL pairs (p=0.26-0.97); (3) post-LP plateau (learners) — Gini slope
  ~0 post-LP (p>0.8), doesn't keep dropping with time, but pre-LP windows too few (n=2) to clinch.
- **Verdict: the Gini↓ de-sparsification is robust but most parsimoniously EXPERIENCE/time-on-task;
  a learning-specific component (LP-locked plateau) is suggestive but NOT statistically established**
  (non-learners drop comparably; n=4 non-learners underpowers the interaction). Report §3.3/synthesis/
  headline tempered: "learning de-sparsifies" → "experience de-sparsifies; learning-specificity unproven".

### 2026-06-06 — Optimal-lag battery: confirms directionality not robust
- Generalised `analyze_ifi.py` to per-dim OPTIMAL LAG (`analyze_ifi.py lag`; +kruskal try/except).
  Agrees with IFI: animals-as-n null for every pair; dims-as-n flags V1-RSC (naive lag>0 t=0.046,
  nai-vs-exp rank-sum p=0.007, KW p=0.011). Starkest artifact: V1-RSC nai→exp lag Δ=0 (p=1.0) by
  animal vs p=0.007 by dims. Both directional readouts converge → directionality not robust;
  dims-as-n manufactures it. Report §3.2 updated.

### 2026-06-06 — IFI directionality battery (animals & dims as n): directionality NOT robust
- Added per-dim IFI/optimal-lag to `subspace_window` (+test) and `analyze_ifi.py`: per pair, IFI-vs-0
  per epoch (t + Wilcoxon), naive-vs-expert (paired t), RM-ANOVA + Friedman + Holm post-hoc —
  ANIMALS-as-n; and dims-as-n (t/Wilcoxon vs0, rank-sum, Kruskal-Wallis).
- **Result: animals-as-n IFI is NULL for every pair** (no IFI!=0 per epoch, no naive→exp change,
  RM-ANOVA n.s.; only fragile n=4 single-test hits RSC-SUB/CA1-SUB). **dims-as-n manufactures
  "significant" directionality** (CA1-RSC int IFI>0 t=5e-4; V1-RSC naive IFI>0 + naive-vs-exp + KW)
  that VANISHES at the animal level → pseudoreplication false-positives (great illustration of the
  dims-as-n hazard). Earlier trajectory IFI trends (CA1→DG/CA1→RSC) downgraded to suggestive.
- CAVEAT: epoch windows short (~10 trials) + in-sample lag curves → IFI is the noisiest readout;
  per-dim OPTIMAL LAG is more robust (battery-testable if directionality pursued).
- **Robust learning signal stands: Gini↓ (participation), early (naive→int), LMM p~1e-5.** Strength
  flat (even dims-as-n). Rotation = noise. Direction = weak/artifact. Report §3.2 + synthesis +
  headline updated.

### 2026-06-06 — Dimensions-as-n check (Buzsáki unit): strength-null is genuine
- Added per-dim export (`subspace_window.sig_mask`, `run_epochs` -> `epoch_dims.csv` 2311 dim-rows /
  224 sig) + `analyze_dims_as_n.py` (pool sig canonical dims across animals, rank-sum across epochs).
  **Result: per-dimension CC flat across learning for every pair (all U-p>0.05; V1-RSC borderline
  ~0.06).** Since dims-as-n inflates N ~5-15x, this confirms the coupling-STRENGTH null is real, not
  power-limited -> learning signal is participation (Gini) + direction (IFI), not magnitude. Report
  §3.7 + methods + synthesis updated. Caveat (nested unit) stated; not our inferential test.

### 2026-06-06 — Parametric stats + epoch analysis
- Added `mixed_effects.lmm_slope` (random-slope LMM population slope; pools all windows; +2 tests)
  and a parametric section to `analyze_trajectory.py` (per-animal-slope **t-test** + **LMM**
  beside Wilcoxon). Result: **CA1-RSC Gini↓ is decisive** — Wilcoxon 0.008 / t 0.003 / **LMM
  4e-5** (n=8, all agree); CA1-DG solid (all ~0.02–0.04); CA1-V1 borderline (~0.05–0.06).
  CA3-DG/CA1-SUB significant under LMM only → over-confident at n=4 (random slope unidentifiable;
  t-test disagrees) → suggestive only. **Rule: where t-test and LMM disagree at small n, trust the
  t-test.**
- **Epoch analysis DONE** (`results/epoch_metrics.csv`, 148 rows; `analyze_epochs.py`, fig
  `HCV1_CCA_fsexcl_epochs.png`). **TIMING: de-sparsification is EARLY (naive→intermediate), then
  plateaus.** CA1-RSC: int−nai Δ=-0.054 (W=0.039/t=0.016/**LMM=7e-5**), exp−nai (W=0.008/t=0.005/
  **LMM=8.5e-6**), exp−int **n.s.** CA1-DG same pattern (int−nai LMM=1.9e-5; exp−int n.s.). CA1-V1
  same sign, borderline. n=4 pairs: paired-t hits but LMM unidentifiable → suggestive. Report §3.6
  + methods §2.10/2.10b updated. FS-incl epoch run launched for completeness.

### 2026-06-06 — Split-half + FS + figures DONE; rotation=noise; FS-robust
- Both FS runs complete (`trajectory_windows.csv` FS-excl 636 rows w/ sh_x; `_fsincl.csv` 647).
  Regenerated full figure set both conditions (points+SEM bars, mean±SEM bands+faint lines,
  all-relationship slope heatmaps, rotation-vs-floor, learners-vs-non) → vault attachments. Updated
  the vault report (§3.1–3.5, synthesis table, methods, next-steps) + vault `log.md`.
- **Rotation = NOISE:** cross-window rotation ≤ split-half floor for every pair (all p>0.05, both
  FS) → no reorientation; retracted "CA3-DG stable backbone" (its low rotation tracks its low floor).
- **FS-robust:** FS-incl reproduces CA1-RSC Gini↓ (3/3 axes p=0.008/0.039/0.008) + rotation=noise;
  CA1-DG/V1 Gini slightly attenuated. Surviving signals: Gini↓ + CA1→DG IFI↑; strength flat.
- NEXT: neuron-count-matched subsampling; learning-vs-time-on-task (more non-learners / pre-post-LP);
  (deferred) depth/ISI membership×properties.

### 2026-06-06 — Split-half noise floor + FS toggle + figure overhaul (in progress)
- Added within-window split-half principal-angle floor to `subspace_window` (sh_x/sh_y; +2 tests).
  **Result: cross-window rotation does NOT exceed the floor for any pair** (Δ(rot−floor) negative or
  ~0, all p>0.05; CA1-RSC/CA1-CA3 even trend rotation<floor). So the ~80° window-to-window rotation
  is **estimation noise, not reorientation** — kills any "subspace rotates with learning" reading;
  real signals remain Gini↓ and IFI. (Well-controlled null.)
- `run_trajectory --include-fs` (FS-included variant) + FS-excluded re-run with split-half: FS-excl
  DONE (`trajectory_windows.csv`, now has sh_x; levels/IFI reproduce prior run exactly), FS-incl
  RUNNING (`trajectory_windows_fsincl.csv`). `figs_report.py` overhauled: per-animal points+SEM bars,
  mean±SEM bands + faint per-animal lines, ALL-relationship slope heatmaps, rotation-vs-floor,
  learners-vs-non, FS-excl/incl. `analyze_trajectory.py`: CSV arg + rotation-vs-floor table.
- TODO when FS-incl lands: regenerate both figure sets, update vault report (rotation=noise; FS
  comparison; new plot styles), commit results.

### 2026-06-06 — Vault report + figures
- Wrote the Hippocampus-V1 communication-subspace report in ResearchVault (full LaTeX methods,
  4 wiki-linked figures via `scripts/figs_report.py`), linked from the project hub, vault `log.md`
  updated. NEXT (agreed): rotation split-half noise floor; neuron-count matching; learner-vs-non
  contrast; (deferred) load depth/ISI for membership×properties.

### 2026-06-06 — Full-suite Frame A landed; structural learning effect
- Ran full-suite trajectory (16 animals, 636 windows) after speed/n_sig fixes. KEY RESULT: learning
  reshapes subspace STRUCTURE not magnitude — Gini↓ (more units recruited) in CA1-RSC/CA1-DG/CA1-V1
  (CA1-RSC robust across all 3 axes incl. LP-relative), CA1→DG IFI↑ with learning; CA1→RSC flow
  significant at baseline. Strength (cc1/mi_sig) does not track learning. Details + caveats above.
- Launched `run_transition.py` (uncued→cued full suite).

### 2026-06-06 — Full subspace suite + all-animals + transition
- Built `subspace_window` (n_sig, mi_sig, IFI, optimal lag, Gini, weights/members) — TDD, committed.
- Rewired `run_trajectory.py` to emit the full suite per window (+ rotation + Jaccard) over ALL 16
  animals (learner-tagged), 3 learning axes; split expensive fit (rich CSV) from fast analysis
  (`analyze_trajectory.py`). Rewired `run_transition.py` to full suite + cross-condition rotation.
- User directives this turn: run on all animals first (then learner/non split); check uncued→cued
  transition; **don't look at CC1 only — all sig subspaces, IFI + optimal lag, angles, membership**;
  write this project log + wire CLAUDE.md to it.

### 2026-06-05 — Reframe from two reference papers; continuous regime; Frame A
- Read Gonzalez/Buzsáki (subspace) + Han/Helmchen (top-down predictions) in full. Diagnosed our
  overfitting (sample regime) and underpower (N=sessions). Wrote `OPPORTUNITIES.md`.
- Validated continuous-regime pCCA (`prototype_continuous_pcca.py`). Built Frame A v1
  (`trajectory.py`, `run_trajectory.py`): clean within-animal trajectories, sign-heterogeneous.

### 2026-06-05 — Honest learning verdict; pseudoreplication caught
- Built spatial paired test + mixed-effects/per-animal-collapse tests. Found the pooled landmark
  test pseudoreplicates (n = animals × landmarks); honest tests are ~null. `STATE.md` §3 caveat.
- Explored the n-ladder (animals → animal×landmark → significant dims): p depends on chosen unit.

### 2026-06-05 — Initial consolidation
- Reconciled the two arms (spatial 66-config sweep vs landmark 44-config); wrote `STATE.md`,
  `GOTCHAS.md` (CRLF in CSVs), committed loose landmark outputs, quarantined-by-documentation the 4
  overfit configs (`landmark50_res_{fix30,var75,var85,var95}`).
