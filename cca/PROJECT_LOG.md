# PROJECT LOG — Tom-learning CCA

**Purpose.** Durable, cross-session memory. Read this first (with `STATE.md`) at the start of any
session; append a dated entry after any meaningful work. Newest entry on top. Keep entries terse
but self-contained — a future session must be able to resume from here alone.

**Doc map.** `STATE.md` = current findings + canonical configs (the verdict). `OPPORTUNITIES.md` =
the two-paper reframe + plan. `GOTCHAS.md` = bugs not to reintroduce. `NOTES.md` = older dev log.
`UNDERSTANDING.md` / `UNDERSTANDING_temporal.md` = original specs. **This file = the running
narrative + state of play.**

---

## CURRENT STATE (2026-06-06)

**Branch:** `cca-consolidation` (NOT merged to main; nothing pushed). ~10 commits this arc.
**Tests:** 217 passing (`cd cca && PYTHONPATH=src python -m pytest -q`).

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

**Running / pending right now:**
- `run_trajectory.py` FULL-suite over all 16 animals — IN PROGRESS (background, ~25 min after
  speed fixes), log `results/.sandbox_scratch/run_trajectory_full.log`; writes
  `results/trajectory_windows.csv` incrementally per animal. When done: run `analyze_trajectory.py`.
- `run_transition.py` — to run AFTER trajectory (avoid two 1.6 GB/animal loads at once).
- **Speed/correctness fixes applied (2026-06-06):** window_subspace was ~11 s/window (significance
  + 21-lag scan on ~41k samples). Now cap each window to a contiguous ~6 k-bin block (grown to span
  ≥ N_FOLDS+1 trials so the CV is valid) → ~2.3 s/window. Fixed `n_sig` overcounting: significance
  now compares the *held-out* per-dim CC to the *dominant*-dim circular-shift threshold (was an
  in-sample per-dim test → 23 sig dims; now ~3). Driver writes CSV incrementally + unbuffered.
  CAVEAT: lags are in running-bin units (non-running bins removed), so IFI/optimal-lag are
  approximate (≈±250 ms); a segment-aware lag is a future refinement.

**Findings so far (honest):**
- The pooled landmark "CA3-DG strengthens" headline was **pseudoreplication** (n = animals ×
  landmarks). Under honest per-animal / mixed-effects tests, no pair survives well; CA3-DG is
  directionally consistent but n=4 (signed-rank floor 0.125). See `STATE.md` §3.
- Continuous-regime pCCA: intra-hippocampal coupling (CA3-DG, CA1-DG, CA1-CA3) > hippocampal-
  cortical (CA1-RSC, CA1-V1). CA3-DG strongest.
- Frame A v1 (CC1, trial-fraction axis): within-animal trajectories are often **clean** (|r| up to
  0.94) but the **slope sign is animal-specific** → across-animal null. Hypothesis: trial-fraction
  ≠ learning stage under varying learning rates → use performance / LP-relative axes (now built).

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
