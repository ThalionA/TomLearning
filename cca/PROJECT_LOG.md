# PROJECT LOG — Tom-learning CCA

**Purpose.** Durable, cross-session memory. Read this first (with `STATE.md`) at the start of any
session; append a dated entry after any meaningful work. Newest entry on top. Keep entries terse
but self-contained — a future session must be able to resume from here alone.

**Doc map.** `STATE.md` = current findings + canonical configs (the verdict). `OPPORTUNITIES.md` =
the two-paper reframe + plan. `GOTCHAS.md` = bugs not to reintroduce. `NOTES.md` = older dev log.
`UNDERSTANDING.md` / `UNDERSTANDING_temporal.md` = original specs. **This file = the running
narrative + state of play.**

---

## ⏳ RUNNING (2026-06-12) — 10 ms full-suite re-run, both FS (directionality focus)

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
