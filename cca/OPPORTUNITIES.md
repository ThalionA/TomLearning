# OPPORTUNITIES — reframing the Tom CCA analysis from two reference papers

**Written 2026-06-05.** Synthesis of the methods in two papers our pipeline descends from,
mapped onto Tom's data and design, with a concrete reframe. Read alongside `STATE.md`
(current state) and `UNDERSTANDING.md` (original spec).

Reference papers:
- **Gonzalez … Buzsáki (2026)** — *Subspace communication in the hippocampal–retrosplenial
  axis* (Nature). pCCA communication subspaces in DG–CA3–CA2–CA1–RSC; membership, rotation,
  reactivation. This is the direct scientific antecedent (CA1↔RSC, hippocampal–cortical).
- **Han & Helmchen (2024)** — *Behavior-relevant top-down cross-modal predictions in mouse
  neocortex* (Nat Neurosci). CCA between S1↔PPC; residual CCA + the **information-flow index
  (IFI)** our pipeline already uses comes from here.

Both descend from the Semedo communication-subspace line; our pipeline is a port of that
machinery. So we are already in the right family — the issue is *regime* and *question*, not
the core method.

---

## 1. The binding constraint (why we kept hitting a wall)

**Tom's design: one recording session per animal, ~12–16 animals, ~100–140 trials each.**
Learning happens *within* that single session; animals learn at different rates, which is why
activity is aligned into naive / intermediate / expert epochs.

The reference papers get their statistical power from **many sessions** (Han: 40–118 sessions;
Gonzalez: 10–15) and from **within-session interleaved conditions** (matched/mismatch, two
mazes) — *not* from animal count. We have ~16 sessions total and a longitudinal (not
interleaved) learning contrast. **No statistical method makes a 3-epoch, cross-animal CC-
magnitude contrast well-powered at n=4–16.** That is the root cause of everything in `STATE.md`
§3 — the pooled test only looked significant via pseudoreplication; honest per-animal / LMM
tests are null because the design is low-N for that question. It is the design, not the test.

**Implication:** keep the questions and rigor of the papers, but reframe "learning" to extract
*within-animal* power, and lean on *structure* (which needs no learning contrast).

---

## 2. The data is much richer than we load

Our `dataio` currently pulls only `units.idx`, `units.idx_fs`, `regions_label`, and the
spatial firing tensor. The export (`TF*_export.mat`) actually contains (verified by inspection):

| Field | Shape | Use (paper analysis it enables) |
|---|---|---|
| `binned_spikes` | (7.9M, 206) | 1 ms continuous spikes — the fine-timescale regime both papers use |
| `units/depth` | (206, 1) | laminar depth — Gonzalez deep-vs-superficial CA1 membership |
| `units/isi/histogram`, `isi/mean/std` | (1000, 206) | burst index — member vs non-member burstiness |
| `units/waveform_*`, `waveforms` | … | proper cell-type classification (beyond binary FS) |
| `analysis_spatial/firing/cued`, `/uncued` | — | interleaved within-session condition contrast |
| `analysis_behaviour/masks/*` | (7.9M,1) | movement / tunnel_cued / quiescence / darkness states |
| `data_behaviour/*` | (7.9M,1) | position, velocity, licks, valve (reward), landmark, trial idx |
| `params_main/cca` | — | **Tom already ran a CCA** (his `cca_pairs`, warp, folds) |
| `units/tuning_score`, `reliability` | (206,137) | per-unit spatial tuning / reliability |

**Scope note (2026-06-05):** for now we do **not** load depth / ISI / waveforms; FS-vs-regular
(`idx_fs`) is the only cell-type split. Those fields are catalogued here as a known, deferred
opportunity — the member-characterization analyses (rate / burst / depth) become available the
moment we choose to load them.

---

## 3. The shared methodological template (both papers)

They converge on the same recipe, and it directly fixes our failures:

1. **PCA → a fixed, modest number of components (≈30), then (p)CCA on the components.** Not
   `k = n_samples / 15` capped at 30. Han states the rule of thumb explicitly: **~50 samples
   per variable** for a stable CCA. Our `samples_per_pc = 15` on ~2,000 epoch-samples is exactly
   what produced the `var95 / fix30` held-out-CC → 0.999 saturation (`STATE.md` §4). **Fixing the
   ratio to ≥50 kills the overfitting.**
2. **Fit on a large, continuous span of data**, not 10-trial epoch slices. Whole session, or a
   sliding window. This is where the samples come from.
3. **pCCA controlling for the third area** (all other simultaneously recorded areas as Z) as the
   *primary* method — isolates direct CA1↔RSC from shared DG/CA3 drive. Gonzalez Ext Fig 2d:
   pCCA correlations are substantially below CCA, i.e. the control matters.
4. **Residual = subtract the condition/stimulus-triggered average** before CCA (trial-trial
   covariation). = our residual CCA. ✓ IFI directionality (lagged CCA) = Han's exact method. ✓
5. **Unit of analysis = session**; Wilcoxon signed-rank (paired) / rank-sum (distributions);
   **subsample to match neuron counts** because subspace dimensionality and CC *scale with N*
   (Gonzalez Ext Fig 1g, R = 0.39–0.69); surrogate everywhere (circular-shift for significance,
   weight-shuffle for spatial info, trial-shuffle for #dims; threshold = mean + 3 SD or 95th pct).

---

## 4. The reframe — all three frames (per direction chosen 2026-06-05)

### Frame A — Learning as a continuous within-animal trajectory *(lead)*
Don't bin into 3 epochs. Fit the subspace on a **sliding window over the ~140 trials** (or
regress subspace metrics — dominant held-out CC, #significant dims, IFI, Gini — on trial number
or behavioural performance) so each animal yields a **slope**. Combine across animals by
sign-consistency / a mixed model with animal random effect. Uses within-animal trial count for
power, drops the arbitrary epoch binning and the landmark pseudoreplication, and directly asks
"how does the circuit reconfigure *as* the animal learns."

### Frame B — Keep 3 epochs, fix the regime
Retain naive/intermediate/expert for continuity with the existing sweep, but fit on **continuous
data** with `samples_per_pc ≥ 50` and **pCCA primary**. Honest, comparable to prior runs, but
still n ≤ 16 across animals for the contrast — so treat as confirmatory-secondary, reported with
the continuous trajectory (Frame A) as the powered version.

### Frame C — Structure-first (engagement), learning second
Lead with readouts that need **no** learning contrast and are robust per session:
- subspace **dimensionality** (# significant pCCA dims) and **Gini** sparsity per pair/region;
- **membership × intrinsic properties** — firing rate and FS-vs-regular now; depth / burst
  deferred (§2);
- **within-session interleaved contrasts**: cued vs uncued, fixed vs random reward — true paired
  contrasts per animal (the Tom analogue of matched/mismatch).
Overlay the continuous-learning trajectory (Frame A) on top of these structural metrics.

---

## 5. What we keep / change / add

**Keep (already correct):** residual CCA, IFI/lagged directionality, circular-shift null,
membership/Gini (`membership.py`), principal angles + split-half noise floor (`subspace.py`),
pCCA implementation (`partial.py`), crosspair similarity (`crosspair.py`).

**Change:**
- **Sample regime:** `samples_per_pc` 15 → ≥50, or fix PCA at ~30 comps with a hard
  samples-≥-50×k guard (`core.choose_k`). Drop / quarantine the overfit high-k configs.
- **Fit span:** fit on continuous running data (reuse `dataio.area_activity_50ms`,
  `temporal_segments`) instead of 10-trial epoch tensors.
- **pCCA primary:** make `partial.partial_cca_cv` the default path (Z = all other recorded
  areas), not a side-branch.
- **Unit/stats:** report per-animal (session) metrics; add neuron-count-matched subsampling for
  any cross-pair/region comparison.

**Add:**
- Frame-A continuous-trajectory driver (sliding-window subspace metrics vs trial/performance).
- Within-session condition contrasts (cued/uncued, reward regime) from the masks already in the
  export.
- *(Deferred)* depth / ISI / waveform loading → member-property characterization.

---

## 6. Implementation plan (phased; TDD per `CLAUDE.md`)

1. **Prototype (now):** verify continuous-data pCCA at the correct regime —
   `scripts/prototype_continuous_pcca.py`. Load continuous CA1/RSC/Z activity, running bins
   only, PCA→30, pCCA, whole-trial CV; report n_samples, samples/k, held-out CC. **Success =
   tens of thousands of samples, samples/k ≫ 50, and a stable non-saturated CC** (no 0.999),
   demonstrating the regime fix. *(No depth/ISI/waveforms.)*
2. **Regime fix in core:** add the samples-≥-50×k guard + a `pca_fixed` k-mode; test on synthetic
   undersampled data that the guard prevents CC saturation.
3. **Frame-C structure metrics** on the continuous fit: dimensionality, Gini, membership×{rate,
   FS}, with neuron-count-matched subsampling + surrogates. Session as unit.
4. **Frame-A trajectory driver:** sliding-window subspace metrics vs trial/performance; per-animal
   slope; combine across animals.
5. **Frame-B epoch contrast** on the new regime, reported as secondary.
6. **Within-session condition contrasts** (cued/uncued, reward regime).

Each step: tests first, session-as-unit stats, surrogate controls, honest N reporting.

---

## 7. Honest framing for any write-up
The powered, defensible claims will be **structural** (dimensionality, sparsity, membership,
rotation vs noise floor) and **within-animal trajectory** (does the metric track learning across
trials). A cross-animal "expert > naive CC magnitude" claim is intrinsically limited to n ≤ 16
and should be reported as such — not as the headline. This matches how both reference papers
actually make their robust claims.
