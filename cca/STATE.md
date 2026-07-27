# STATE — Tom-learning CCA project

**Last updated:** 2026-07-27. This is the entry point. It reconciles the three analysis
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
2. **Membership: participation broadens (the one robust signal).** Weight-Gini ↓ over the
   session — CA1-RSC (LMM trajectory p = 4×10⁻⁵; epoch naive→intermediate LMM p = 7×10⁻⁵,
   expert−intermediate n.s.) and CA1-DG (LMM p = 1.9×10⁻⁵). Early (by intermediate), then
   plateaus. Parametric, FS-invariant, and **invisible to CC₁**.
   **⚠ Attribution to learning is NOT established** — non-learners de-sparsify comparably and
   the `trial_frac × learner` interaction is n.s. for every pair. Most parsimonious reading is
   **experience / time-on-task**, with an LP-locked component suggestive but unproven at this N.
3. **Direction: a flow *exists*; its *change* is underpowered.** Held-out segment-aware IFI
   window sweep, animals-as-n: **CA1→RSC +0.079 at ±50 ms, t₁₁ = 5.0, p = 3.9×10⁻⁴** (survives
   Bonferroni across nested windows) — the Gonzalez & Buzsáki direction. **V1→RSC** and
   **SUB→CA1** are the other FS-robust tight-lag flows; CA3-DG is bidirectional. But this is
   session-pooled: it speaks to the *existence* of a flow, not its change. Directional *change*
   with learning holds its sign but is weak and underpowered (n = 4–6/pair); the CA1→V1 IFI
   rise is the one slope supported at the honest unit and it is **FS-fragile** (null with FS in).
4. **Orientation: NULL.** Cross-window rotation is at or below the split-half noise floor for
   every pair, both FS conditions (all p > 0.05). No reorientation. (See `GOTCHAS.md` — the
   ~80° floor is a ~1-D subspace, not a bug.)
5. **Very-early trials: no first-trials jump.** Trials 1/4/7/10 and first-5/7/10 blocks —
   strength flat (FS-robust); the de-sparsification is gradual / largely post-trial-10. The
   fast early-then-plateau effects are *cortical* and FS-fragile (CA1-RSC IFI, participation-
   Gini; V1-RSC early rotation above floor).
6. **Nonlinearity: largely absent.** Kernel CCA edges linear in 58–66 % of cells but the
   median KCCA − linear gap is only +0.015, and CA3-DG (the strongest pair) is ≈ linear. The
   subspace is largely linear.

**Headline (Contested).** Over the task the hippocampal–cortical communication subspace
**broadens** — recruiting more neurons (Gini↓, CA1-RSC/CA1-DG, LMM p ~ 10⁻⁵, early then
plateau) — rather than changing coupling magnitude, direction, or orientation. Whether the
broadening is *learning* or *experience* is the open question.

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
