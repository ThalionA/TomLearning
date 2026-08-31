# STATE — Tom-learning CCA project

**Last updated:** 2026-08-31 (7c/7d/7e added: membership reordering; V1-RSC gain-dip candidate; learning-speed null). This is the entry point. It reconciles the three analysis
arms, names the canonical configs, and states the findings honestly. `README.md` is
quick-start; `UNDERSTANDING.md` / `UNDERSTANDING_temporal.md` are the design specs;
`PROJECT_LOG.md` is the state of play; `NOTES.md` is the older chronological log. When they
disagree with this file, this file wins.

> **Reorientation (2026-06-13 → 06-18).** The **temporal (running-state) arm is now the
> primary analysis** and the basis of the write-up. The landmark and spatial arms below are
> retained as the earlier event-locked / full-corridor passes; their headline (CA3-DG
> strengthening) did **not** survive honest per-animal testing and is **not** the project's
> current claim. Read §3.0 before §3.1.

---

## 1. The question

How does communication between hippocampal/cortical area pairs change across spatial
learning (naive → intermediate → expert) in the Tom cohort (12 learners + 4 non-learners,
one session each)? Sub-questions: (1) does coupling *strength* change? (2) does *direction*
(who leads) change? (3) *which units* carry it? (4) does the *subspace* reorient? CCA is fit
per area pair on PCA-reduced spike-rate tensors; significance is against a circular-shift
null; canonical correlations are read on held-out trials (leak-free CV).

Eight pairs: CA1-V1, CA1-DG, CA1-CA3, CA1-RSC, CA1-SUB, V1-RSC, RSC-SUB, CA3-DG.

---

## 2. The three arms

| | **Temporal arm** (Arm A, running-state) | **Landmark arm** (Arm B) | **Spatial arm** |
|---|---|---|---|
| Question | Strength/direction/membership/orientation over running bins, resolved along the session | Trial-trial communication *around landmark events* | Strength/direction over the *whole corridor*, per epoch |
| Binning | **Gaussian σ=2.5 ms → 10 ms bins**, z-scored on the engaged reference | 25 ms or 50 ms windows, ±200 ms around landmarks 1–6 | 200 bins × 2.5 cm |
| Fit unit | per (animal, pair, **window / epoch / trial-bin**) | per (animal, pair, **landmark**, epoch) | per (animal, pair, **epoch**) |
| Data | `results/*_bin10{,_fsincl}.csv` | `results/landmark{25,50}_*.pkl` (44) | `results/stage2_*_lag*.pkl`, `stage3_*` (132) |
| Summary | **`results/bin10_tables.md`** (§A–H) | `results/sweep_landmark_summary.csv` | `figures/sweep_summary_spatial.{csv,xlsx}` |
| Test of record | **paired / one-sample *t*** (animals-as-n); LMM where nested | paired per-animal Wilcoxon + FDR | paired per-animal Wilcoxon + FDR |
| Status | **PRIMARY — carries the write-up** | Superseded; honest verdict null | Exploratory / corroborative |

**How they relate.** The spatial arm averages over the full traversal and finds *weak*
coupling that barely moves with learning. The landmark arm locks to behavioural events and
concentrates the signal, but its apparent learning effect was pseudoreplicated (§3.2) and is
null at the honest unit. The temporal arm keeps the running-state regime of the reference
papers (Gonzalez & Buzsáki 2026; Han & Helmchen 2024), adds membership/rotation readouts, and
is where the one robust effect lives. **All three agree that coupling *strength* does not
change with learning** — the landmark arm's CA3-DG "strengthening" was the pooled,
anti-conservative figure, and the temporal arm's strength-null is genuine rather than
power-limited (flat even under the inflated dims-as-n unit).

---

## 3. Findings

### 3.0 Temporal arm — the current verdict (PRIMARY, 2026-06-17)

Epistemic state: **Contested / exploratory.** N is small (12 learners + 4 non-learners, one
session each, 4–10 animals per pair); significance is a per-pair family with **no cross-pair
correction**; credibility rests on consistency across independent learning axes and across
two analyses (trajectory + transition), not on isolated p-values. Full pair-by-pair numbers
live in `results/bin10_tables.md` (§A–H) and the vault report's §3.9 tables.

**Coupling hierarchy (held-out CC₁, per-animal mean, 16 animals, FS-excluded):**

| CA3-DG | CA1-CA3 | CA1-DG | V1-RSC | RSC-SUB | CA1-V1 | CA1-RSC | CA1-SUB |
|---|---|---|---|---|---|---|---|
| 0.315 | 0.273 | 0.217 | 0.169 | 0.127 | 0.094 | 0.061 | 0.053 |

Same rank order FS-included, uniformly higher (CA3-DG 0.409). Intra-hippocampal ≫
hippocampal–cortical. `n_sig` ≈ 3.3–5.7 dims (FS-excl).

**What changes with learning — and what does not:**

1. **Strength: NULL, and genuinely so.** All trajectory strength slopes n.s.; transition
   n.s.; flat even under dims-as-n (the most anti-conservative unit). This is a real null,
   not a power failure.
2. **Membership: participation broadens for CA1-RSC and CA1-V1 — CA1-DG does NOT survive
   (re-tested 2026-08-28).** The original headline (Weight-Gini ↓ over the session; CA1-RSC LMM
   trajectory p = 4×10⁻⁵, CA1-DG p = 1.9×10⁻⁵; early then plateauing) was measured with a
   **partner-invariant** metric — see the caveat below. On the connection-specific replacements,
   like-for-like on the same windows (`results/gini_conn_retest_tables.md`):
   - **CA1-DG: dead.** n.s. under `_conn` and `_sig`, both axes, both FS. Drop it.
   - **CA1-RSC: survives, but on the `performance` axis, not `trial_frac`.** Robust under every
     definition, **both sides of the connection**, both FS (`gini_x_conn` p = 6.3×10⁻⁴ / 1.0×10⁻⁴;
     `gini_y_sig` 1.4×10⁻³ / 2.7×10⁻⁴). On `trial_frac` only `_sig` survives (p = 0.019 / 0.038).
   - **CA1-V1: new, on the V1 side.** `gini_y_conn` trial_frac p = 0.027 / 0.0056; performance
     0.029 / 0.015. Invisible to the old `_x`-only battery. CA1 side n.s.
   - ⚠ `gini_*_sig` is NaN when `n_sig` = 0 (56 % of CA1-RSC windows) — its slopes are
     conditioned on windows with a significant dim. ⚠ n = 4 pairs (CA1-SUB, CA3-DG, RSC-SUB)
     cannot star at the honest unit (Wilcoxon floor p = 0.125); their small LMM p-values come
     from pooling windows within 4 animals.
   **⚠⚠ METRIC CAVEAT (2026-07-28) — the Gini that produced this result is partner-invariant.**
   `membership.subspace_contribution` takes an unweighted L2 row-norm over all retained dims, so
   the square-orthogonal `Uc` cancels and the partner area drops out exactly (verified
   analytically, on synthetic data to 4×10⁻¹⁵, and empirically: CA1 `gini_x` correlates at median
   r = 0.981 across five partners). It is an **area-intrinsic** readout of a population's own
   whitened-PCA loading geometry, **not** participation in a communication subspace. Two
   connection-specific replacements (`gini_*_conn`, `gini_*_sig`) exist and **the re-test is now
   done** (2026-08-28) — its outcome is the bullet list above, and it is what should be cited as
   the communication-subspace result. The area-intrinsic numbers stay here only as the historical
   headline. See `PROJECT_LOG.md` (2026-07-28 and 2026-08-28).
   **⚠ Attribution to learning is NOT established** either, and the 2026-08-28 move to the
   `performance` axis does **not** rescue it — non-learners de-sparsify with a *larger* median
   slope than learners (CA1-RSC `gini_x_conn`: −0.287 vs −0.106, both FS), at n = 3–4 where the
   p-values sit on the Wilcoxon floor, and with a non-learner `performance` range compressed
   enough to inflate any slope taken w.r.t. it. The `trial_frac × learner` interaction is n.s. for
   every pair. Most parsimonious reading remains **experience / time-on-task**.
3. **Direction: a flow *exists*; its *change* is underpowered.** Held-out segment-aware IFI
   window sweep, animals-as-n: **CA1→RSC +0.079 at ±50 ms, t₁₁ = 5.0, p = 3.9×10⁻⁴** (survives
   Bonferroni across nested windows) — the Gonzalez & Buzsáki direction. **V1→RSC** and
   **SUB→CA1** are the other FS-robust tight-lag flows; CA3-DG is bidirectional. But this is
   session-pooled: it speaks to the *existence* of a flow, not its change. Directional *change*
   with learning holds its sign but is weak and underpowered (n = 4–6/pair); the CA1→V1 IFI
   rise is the one slope supported at the honest unit and it is **FS-fragile** (null with FS in).
   *Added 2026-08-15, re-done WHOLE-SESSION 2026-08-17 (meeting ask 3):* the same question on
   the **average over all significant CCs** (`cc_ifi_overall_test_bin10*.csv`, uncapped
   `run_lag_curves`, ~9 sig CCs per animal-pair): **CA1→RSC +0.053 p = 0.005 (FS-excl) /
   +0.068 p = 8×10⁻⁵ (FS-incl), surviving BH across the 8 pairs on every look**; and,
   consistent across FS and units but weaker, **DG→CA1** (−0.053 p = 0.014 / −0.044 p = 0.020),
   **DG→CA3** (−0.061 p = 0.027 / −0.082 p = 0.013) and **V1→RSC** (+0.038 p = 0.024 / +0.030
   p = 0.050); CCs-as-n p ≤ 10⁻³ for all four. CC₁ of the same table reproduces the numbers
   above to 3 dp. The CC₁-only SUB→CA1 flow is *not* carried by the all-CC average (n.s.).
   An animal's significant CCs share a sign more than chance (item 1 re-read: 64/71 mixed vs
   70.2 expected, p < 0.001). Naive→expert on the all-CC IFI: 0/8, both FS and both units
   (ask 2). ⚠ The 08-15 numbers on the first-20-trials cap are superseded.
4. **Orientation: NULL — and it also does not change with LAG.** Cross-window rotation is at
   or below the split-half noise floor for every pair, both FS conditions (all p > 0.05). No
   reorientation. **This claim is safe because §G of `bin10_tables.md` tests it at CC₁**
   (CC₁ floors 36–69°), not only at top-3. ⚠ The top-3 column of §G, and the
   `fig_rotation_floor` figure in `figs_report.py`, use the 3-dim angle whose floor is
   ~78–82° — *unmeasurable* (see `GOTCHAS.md`, 2026-07-29). **The figure should be switched
   to CC₁; the verdict does not change.**
   Added 2026-07-29 (`run_lag_subspaces.py`): the CC₁ subspace likewise never separates from
   its floor across **±250 ms of lag**, 8/8 pairs, both FS — the subspace is the same one at
   every delay in that range. Caveat in both cases: the floor is measured from half-data
   fits while the comparison uses full windows, so the floor sits above the comparison's true
   noise level and the test is conservative — it excludes a *large* rotation, not a modest one.

> **⚠ 2026-08-03 — two bugs found by adversarial verification; findings 4 and 7 below are
> the CORRECTED versions.** (a) `run_lag_subspaces` exported only the X-area d=1 floor and
> subtracted it from both areas' angles — the shared-floor flaw `split_half_floor` exists to
> avoid (mean |X−Y| floor difference 12.8°, misplacing the estimability gate in 9/71 cells).
> (b) The fixed-subspace arm fitted weights in one residual space and applied them in another,
> so "identical weights across epochs" was not an identical transform. **Consequence: the
> CA1-RSC integration-window narrowing was an ARTEFACT of (b) and is gone — items 5/6/7 are
> now 0/32 in BOTH FS conditions.** Verdicts elsewhere survived, but per-cell values moved a
> lot (`peak_lag_ms` old-vs-new correlation only 0.35), so no pre-2026-08-03 per-cell number
> should be quoted. Also corrected: the epoch contrast is **32/32 null FS-excluded**, not
> 31/32.

7. **Feedforward vs feedback: NOT separable subspaces (2026-07-29, re-verified 2026-08-03).**
   Fitting the subspace at
   +50 ms (X leads) and at −50 ms (Y leads) gives two subspaces whose principal angle never
   clears the noise floor — **0/8 pairs, both FS**. Connection-specific Gini does not differ
   between them in any pair, and the FF/FB asymmetry does not change with learning (0/8 pairs,
   both FS; 64 tests, min p = 0.056). **At ±50 ms this is one subspace read at two delays, not
   two directions of flow** — so "split the subspace into FF/FB" is not supported by this data.
   *But the strength through that one subspace is time-asymmetric*: the sign of cc₁(+50) −
   cc₁(−50) agrees across both FS conditions in **8/8** pairs and recovers hippocampal anatomy
   (CA3→CA1, DG→CA1, DG→CA3). Treat as a consistency check, not a replication — FS-included
   is FS-excluded plus fast-spiking units, so the two share most of their data. Per-animal
   sign consistency is only 4/7, 6/11, 3/5 for those three pairs; the animal-consistent flows
   are **RSC→SUB (6/7)** and **V1→RSC (7/9)**.
   ⚠ **CA1-SUB reads CA1→SUB here (5/7 animals), contradicting finding 3's SUB→CA1.** Methods
   differ (subspace refit at a fixed lag vs IFI integrated from one session-level fit) — to be
   reconciled, not ignored.
5. **Very-early trials: no first-trials jump.** Trials 1/4/7/10 and first-5/7/10 blocks —
   strength flat (FS-robust); the de-sparsification is gradual / largely post-trial-10. The
   fast early-then-plateau effects are *cortical* and FS-fragile (CA1-RSC IFI, participation-
   Gini; V1-RSC early rotation above floor).
6. **Nonlinearity: largely absent.** Kernel CCA edges linear in 58–66 % of cells but the
   median KCCA − linear gap is only +0.015, and CA3-DG (the strongest pair) is ≈ linear. The
   subspace is largely linear.
7c. **Trial 1 → 2: the channel is stable but its membership reorders (NEW 2026-08-31).**
   Frozen-subspace per-unit carrying (unit × partner-CC1 correlation, ordinals 1..10 leak-free,
   common-bin arm): carrying strength 1v2 null and the 1→2 step is not special vs adjacent
   steps — but the per-unit profile's sim(1,2) sits BELOW the adjacent-step band in 14/16
   cells; global per-animal test (n = 16): FS-excl W p = 3.1×10⁻⁵, FS-incl W p = 0.021.
   Uncorrelated with the trial-1 behaviour deltas. Convergence-to-trained-profile: null.
   Per-cell localisation underpowered (0 BH). `results/trial12_units_tables.md`;
   figure `HCV1_trial12_units_deltas_*`.
7d. **First-trial V1-RSC gain dip (candidate, NOT a claim — 2026-08-31).** Trial 1's
   frozen-CC1 lag curve is depressed at all lags in both FS (level effect; trial 2 already at
   the trials-3..10 band); paired deltas sign-consistent in every arm × FS spec but none
   W < 0.05 at n = 9; IFI delta FS/arm-fragile. Second-cohort question.
   `results/trial12_v1rsc_tables.md`.
7e. **Naive communication does not predict learning speed (exploratory NULL, 2026-08-31).**
   Naive-epoch cc1/IFI vs LP, n = 12 learners: no BH survivor, global coupling n.s.; one
   repeatable lean (CA1-DG strength ρ = −0.69, p = 0.058 both FS, leverage-sensitive).
7b. **Contributing units ARE spatially special — on the cortical side (NEW 2026-08-29).**
   Per-unit connection-specific contribution (`contrib_conn`, epoch pCCA) correlates with
   spatial reliability (mean ±2-trial map correlation, epoch-matched trials). Animals-as-n,
   Fisher-z over epochs; full tables `results/contrib_reliability_tables.md`. Rate predicts
   reliability at rho ≈ +0.4–0.5 everywhere, so the claim is made on the **rate-partialled**
   correlation; FS-robust survivors: **CA1-RSC RSC side** (+0.26 both FS, W p = 0.016),
   **CA1-V1 V1 side** (+0.26/+0.17, W p = 0.002/0.037), **V1-RSC RSC side** (+0.31/+0.41,
   W p = 0.031). The CA1-side link never survives the rate partial, and the area-intrinsic
   `contrib` control is n.s. nearly everywhere — communication-specific, not "loud units".
   Epoch-resolved stability is now TESTED (2026-08-30): 3/318 (FS-excl) and 7/318 (FS-incl)
   epoch contrasts star — below chance; only cross-FS repeat is pooled-DG reliability
   naive→intermediate (weak). Raw (un-partialled) metrics are equally stable (594 tests/FS,
   10 starred, chance ≈ 30). The links are standing properties, not learning products.
   **Tom's precomputed moving-window reliability replicates the finding** (agreement with
   ours ρ ≈ +0.67…+0.89 per cell): same three rate-partialled cells plus **CA1-DG DG**
   (both FS) — RSC/V1 metric-robust, DG-side definition-sensitive.
   **Extensions (2026-08-30):** (a) **tuning dissociates from reliability** — contribution vs
   spatial tuning (z vs the export's shuffle null) is broadly positive raw but has NO FS-robust
   rate-partialled survivor (FS-excl leaves CA1-RSC RSC / CA1-DG DG; both die FS-incl): what is
   special about contributing units is trial-to-trial reliability, not tuning strength.
   (b) **Pooled over all partners**, the rate-partialled link survives for RSC and V1 only
   (both FS) — the cortical-side story at whole-area level. (c) **The same units serve several
   subspaces, with partner-specific differentiation on top:** cross-partner contribution ρ
   +0.44…+0.71 vs an intrinsic-geometry ceiling of +0.82…+1.0; the geometry-partialled residual
   stays positive (+0.17…+0.50) and top-quartile member sets beat independent-draw Jaccard by
   +0.15…+0.35 (CA1/RSC/V1 starred, both FS). Tables `results/contrib_reliability_tables.md`
   (addendum); figures `HCV1_contribtune_forest_*`, `HCV1_contribpool_forest_*`,
   `HCV1_contriboverlap_*`.

**Headline (Contested — and now under review).** Over the task the hippocampal–cortical
communication subspace **broadens** — recruiting more neurons (Gini↓, CA1-RSC/CA1-DG, LMM
p ~ 10⁻⁵, early then plateau) — rather than changing coupling magnitude, direction, or
orientation. Whether the broadening is *learning* or *experience* is one open question; **whether
it is a property of the *connection* at all is the other** (2026-07-28: the Gini behind it is
partner-invariant — see finding 2's metric caveat). Both must resolve before this is a claim.

**Units of analysis (stats policy).** Animals-as-n is the inferential unit. Dims-as-n (the
Gonzalez & Buzsáki convention) is reported **as a power check only** — dims are nested in cells
nested in animals, so it is the most pseudoreplicated unit available. Where pooled OLS is
significant but the animal level and cluster-robust LMM are not, that gap *is* the artefact.

### 3.1 Landmark & spatial arms (earlier passes — superseded, retained for the record)

> **⚠ Pseudoreplication caveat (2026-06-05) — read before trusting the p-values below.**
> The `per_pair_pooled` test in `learning_changes.py:174` pools deltas across all 6 landmarks
> (`n` = animals × landmarks) and treats them as independent — but the 6 landmarks within an
> animal are repeated measures, not independent samples. This inflates significance: several
> "hits" have a **median delta of ≈0.000** yet `p<0.01` purely from the inflated `n` (e.g.
> CA1-RSC in `samp25`). When landmarks are **collapsed to one value per animal** (the honest
> unit, `n` = #animals = 4–10), **no pair survives within-pair correction in the committed
> config — including CA3-DG** (n=4, Wilcoxon p-floor 0.125). So the pooled headline is
> anti-conservative; the honest per-animal verdict is **null**. The real limit is N and
> pseudoreplication, not the FDR family.

**Honest verdict (mixed-effects / per-animal tests, 2026-06-05).** Two non-pseudoreplicated
tests — per-animal collapse and a random-slope LMM (`scripts/learning_changes_mixed.py`,
`src/tom_cca/mixed_effects.py`), within-pair FDR, no cross-pair correction:
- **CA1-RSC** — most reproducible: a directionality (`ifi_weighted`) shift peaking at the
  *intermediate* stage (non-monotonic), ~20 configs across both methods; plus a strength
  (`mi_sig`) increase in `landmark50_res_samp25`. n=8.
- **CA1-DG** — modest monotonic strength increase (expert>naive, +0.05) surviving the LMM in
  `landmark50_res_samp40`. n=8.
- **CA3-DG** — the pooled "headline" — **does not survive** any honest test (n=4 only;
  signed-rank floor p=0.125), though all 4 animals move in the same (+) direction.
- V1-RSC `mi_sig` hits are confined to overfit/high-k 25 ms configs and are non-monotonic —
  discount.

**Pooled (anti-conservative) figures, for the record.** Canonical config
`landmark50_res_samp15`, pooled across landmarks: expert vs naive median Δmi = +0.21,
p = 0.0059 (FDR-pass); expert vs intermediate +0.15, p = 0.0043 (FDR-pass); intermediate vs
naive ≈ 0, p = 0.65. CA3-DG expert-vs-naive survives FDR in 18/40 non-overfit landmark
configs. **These are the pseudoreplicated numbers — do not cite them as the verdict.**

**Spatial corroboration (2026-06-05).** Across all 66 spatial configs there is **1** FDR
survivor total (CA1-V1/`ifi_weighted`) — consistent with chance, i.e. **no robust spatial
learning effect**. CA3-DG expert−naive `mi_sig` is positive in 63/66 configs (committed
config Δ=+0.18) — direction-consistent but never significant (n=4, signed-rank p-floor 0.125).

**Scope caveat.** `per_landmark` scope (FDR across 8 pairs × 6 landmarks = 48 tests): nothing
survives in any config. `per_pair_pooled` (FDR across 8 pairs) is the appropriate family.

---

## 4. Canonical configs & exclusions

- **Temporal (PRIMARY — use this):** Gaussian σ = 2.5 ms smoothing on 1 ms trains → **10 ms
  bins** → z-score over the engaged reference (cued ∧ v ≥ 2 cm s⁻¹); running, in-trial bins
  only; partial CCA against all other recorded areas; PCA **k = 30** (trajectory, IFI window
  sweep), **k = 20** (epochs), **k = 15** (transition); leak-free held-out CV; **IFI headline
  window ±50 ms** (curves to ±250 ms). **FS-excluded and FS-included are co-primary** — report
  both. Outputs tagged `_bin10` (FS-excl) / `_bin10_fsincl`.
- **Landmark (earlier arm):** `landmark50_res_samp15` — healthy (frac_cc_ge_099 = 0.0, expert
  CC ≈ 0.52, ~0.8–2 sig dims/pair).
- **Spatial (exploratory default):** `res_samp15_lag10`.
- **Excluded — overfit (do not read for any verdict):** `landmark50_res_fix30`,
  `landmark50_res_var75`, `landmark50_res_var85`, `landmark50_res_var95`. These saturate
  held-out CC → 0.999 in 5/8 pairs (over-parameterised) and spuriously flag nearly every pair
  as a "survivor." Judge severity by `frac_cc_ge_099_*`, **not** `max_cc`. Confirmed in
  `results/landmark_prune_summary.csv` (`n_overfit_pairs = 5`).
- **Borderline:** `landmark50_res_fix20` (frac_cc_ge_099 ≤ 0.19, flag = False, but ~4 sig
  dims) — not quarantined; interpret with caution.

---

## 5. Open decisions

1. ~~**Spatial-arm learning gap.**~~ **Resolved 2026-06-05** — no robust spatial learning
   effect (1/66 configs has a lone FDR hit); the arm is exploratory/corroborative.
2. ~~**Arm A (running-state) not run.**~~ **Resolved 2026-06-13/17** — built, run at 10 ms
   smoothed for both FS conditions, and promoted to the primary analysis (§3.0).
3. **Learning vs time-on-task — the binding open question.** The Gini↓ broadening appears in
   non-learners too (n=3) and the `learner × time` interaction is n.s. Needs a cleaner
   contrast: more non-learners, or a within-animal pre/post-LP split. Until then the effect is
   attributable to *experience*, not *learning*.
4. **Publication scope.** The headline is now **participation broadening**, not CA3-DG
   strengthening. Decide whether that was pre-registered or is a post-hoc discovery (affects
   how the 8-pair family is framed), and whether the write-up leads with the *existence* of
   the CA1→RSC flow (robust) or the *change* results (mostly null).
5. **Membership × cell properties (deferred).** The export carries `units/depth`, `units/isi`
   (burstiness) and waveforms — load these to ask *which* cells join the subspace as it
   de-sparsifies (deep vs superficial, bursty vs regular), mirroring Gonzalez & Buzsáki 2026.
6. **Neuron-count matching.** Subsample to matched unit counts before any cross-pair
   comparison of dimensionality or strength (the hierarchy in §3.0 is unmatched).

---

## 6. Reproduce

No orchestrator (no Snakemake/Make) — scripts run manually. Tests: `cd cca && PYTHONPATH=src
python -m pytest -q` (265 pass, 2026-07-27). Browse landmark/spatial figures via
`figures/index.html`; temporal figures render into the ResearchVault attachments.

```
# Temporal arm (PRIMARY) — smoothed 10 ms, both FS. Add --include-fs for the FS-included run.
PYTHONPATH=src python scripts/run_trajectory.py      --bin-ms 10 --smooth-ms 2.5
PYTHONPATH=src python scripts/run_epochs.py          --bin-ms 10 --smooth-ms 2.5
PYTHONPATH=src python scripts/run_transition.py      --bin-ms 10 --smooth-ms 2.5
PYTHONPATH=src python scripts/run_ifi_windows.py     --bin-ms 10 --smooth-ms 2.5
PYTHONPATH=src python scripts/run_lag_curves.py      --bin-ms 10 --smooth-ms 2.5
PYTHONPATH=src python scripts/run_trajectory_bins.py --bin-ms 10 --smooth-ms 2.5
PYTHONPATH=src python scripts/analyze_bin10_full.py            # -> results/bin10_tables.md
PYTHONPATH=src python scripts/figs_report.py                   # + figs_paired / figs_units /
                                                               #   figs_area_gini / figs_lag_curves /
                                                               #   figs_trajectory_bins / figs_stats_tables

# Landmark arm — regenerate one config's figures + learning test
PYTHONPATH=src python scripts/learning_changes.py --tag landmark50_res_samp15
PYTHONPATH=src python scripts/summarise_landmark_sweep.py
PYTHONPATH=src python scripts/build_prune_table.py                # overfit table
PYTHONPATH=src python scripts/build_landmark_index.py             # index.html

# Spatial arm
PYTHONPATH=src python scripts/run_stage2.py --max-seconds 1800    # resumable
PYTHONPATH=src python scripts/run_stage3.py
PYTHONPATH=src python scripts/summarise_sweep.py
```
