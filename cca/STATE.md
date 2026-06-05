# STATE — Tom-learning CCA project

**Last updated:** 2026-06-05. This is the entry point. It reconciles the two analysis
arms, names the canonical configs, and states the findings honestly. `README.md` is
quick-start; `UNDERSTANDING.md` / `UNDERSTANDING_temporal.md` are the design specs;
`NOTES.md` is the chronological log. When they disagree with this file, this file wins.

---

## 1. The question

How does communication between hippocampal/cortical area pairs change across spatial
learning (naive → intermediate → expert) in the Tom cohort (12 learner animals)?
Sub-questions: (1) does coupling *strength* change? (2) does *direction* (who leads) change?
(3) *which units* carry it? (4) does the *subspace* reorient? CCA is fit per area pair on
PCA-reduced spike-rate tensors; significance is against a circular-shift null; canonical
correlations are read on held-out trials (5–10-fold CV).

Eight pairs: CA1-V1, CA1-DG, CA1-CA3, CA1-RSC, CA1-SUB, V1-RSC, RSC-SUB, CA3-DG.

---

## 2. The two arms (both fully run, both live)

| | **Spatial arm** | **Landmark arm** (Arm B) |
|---|---|---|
| Question | Strength/direction over the *whole corridor*, per epoch | Trial-trial communication *around landmark events* |
| Binning | 200 bins × 2.5 cm | 25 ms or 50 ms windows, ±200 ms around landmarks 1–6 |
| Fit unit | per (animal, pair, **epoch**) | per (animal, pair, **landmark**, epoch) |
| Configs | 66 (`{res,sig}` × 11 k-rules × 3 lag windows) | 44 (`{res,sig}` × 11 k-rules × {25,50} ms bins) |
| Result pkls | `results/stage2_*_lag*.pkl`, `results/stage3_*_lag*.pkl` (132) | `results/landmark{25,50}_*.pkl` (44) |
| Summary | `figures/sweep_summary_spatial.{csv,xlsx}` | `results/sweep_landmark_summary.csv` |
| Learning test | per-config `d_cc` + uncorrected `p_naive_vs_expert` only | **paired per-animal Wilcoxon + FDR** (`results/learning_changes_*.csv`) |
| Canonical config | `res_samp15_lag10` | **`landmark50_res_samp15`** |
| Status | Exploratory — see §4 gap | **Carried through to a rigorous verdict** |

**How they relate:** the spatial arm averages over the full traversal and finds *weak*
coupling (held-out CC ≈ 0.08–0.25) that barely moves with learning. The landmark arm locks
to behavioural events, *concentrates* the signal (held-out CC ≈ 0.52 expert), and is the
only arm with a paired, multiple-comparison-corrected learning test. **The landmark arm is
the canonical answer to the learning question; the spatial arm corroborates directionally.**

The Arm A (running-state) temporal analysis in `UNDERSTANDING_temporal.md` is designed but
**not run** — no `temporal_arm_*` pkls exist.

---

## 3. Findings

**Headline:** the one robust, reproducible learning effect is **CA3-DG communication
strengthening with expertise**.

- **Canonical landmark config (`landmark50_res_samp15`), pooled across landmarks:**
  - expert vs naive: median Δmi = **+0.21**, p = 0.0059, **FDR-pass** ✓
  - expert vs intermediate: median Δmi = **+0.15**, p = 0.0043, **FDR-pass** ✓
  - intermediate vs naive: Δmi ≈ 0, p = 0.65, n.s.
  - → monotonic, with the jump occurring **at the expert stage**.
- **Reproducibility:** CA3-DG expert-vs-naive survives FDR in **18 of 40** non-overfit
  landmark configs — far ahead of any other effect (next is V1-RSC intermediate-naive at
  7/40, which does *not* reproduce in the expert contrast). Source: pooled-scope rows of
  `results/learning_changes_*.csv` (CRLF-safe scan; see `GOTCHAS.md`).
- **Spatial corroboration:** in `figures/sweep_summary_spatial.csv`, CA3-DG has the largest
  positive `d_cc` (+0.053, cc_expert 0.248 — the highest spatial CC of any pair), consistent
  in sign though not significant.
- **Weaker / not robust:** CA1-CA3 (expert-naive, 2/40), CA1-RSC, CA1-DG, and scattered
  intermediate-naive hits (V1-RSC, RSC-SUB) appear in a minority of configs and do not
  reproduce across contrasts. Treat as suggestive, not established.

**Scope caveat (important):** significance depends on the test family.
- `per_landmark` scope (FDR across 8 pairs × 6 landmarks = 48 tests): **nothing survives in
  any config** — underpowered (n = 4–10 animals per cell against 48 comparisons).
- `per_pair_pooled` scope (FDR across 8 pairs): the survivors above. This is the appropriate
  family — pool landmarks within a pair, then correct across pairs.

So the honest one-liner: **CA3-DG coupling strengthens at the expert stage (pooled-significant,
FDR-robust, reproducible across configs); no per-landmark-resolved effect is detectable; all
other pairs are at best suggestive.**

---

## 4. Canonical configs & exclusions

- **Landmark (use this):** `landmark50_res_samp15` — healthy (frac_cc_ge_099 = 0.0, expert
  CC ≈ 0.52, ~0.8–2 sig dims/pair).
- **Spatial (exploratory default):** `res_samp15_lag10`.
- **Excluded — overfit (do not read for the verdict):** `landmark50_res_fix30`,
  `landmark50_res_var75`, `landmark50_res_var85`, `landmark50_res_var95`. These saturate
  held-out CC → 0.999 in 5/8 pairs (over-parameterised: ~30 canonical dims on too-few paired
  samples) and spuriously flag nearly every pair as a "survivor." Judge severity by
  `frac_cc_ge_099_*`, **not** `max_cc` (one saturated dim pushes max to 0.999 in otherwise-fine
  configs). Confirmed in `results/landmark_prune_summary.csv` (`n_overfit_pairs = 5`).
- **Borderline:** `landmark50_res_fix20` (frac_cc_ge_099 ≤ 0.19, flag = False, but ~4 sig
  dims) — not quarantined; interpret with caution.

---

## 5. Open decisions

1. **Spatial-arm learning gap.** The spatial arm has 132 result pkls and a summary, but its
   `p_naive_vs_expert` is uncorrected and there is no paired per-animal test analogous to
   `scripts/learning_changes.py`. Decide: (a) build a spatial analogue of the paired+FDR test
   to get a comparable verdict, or (b) formally declare the spatial arm exploratory/superseded
   by the landmark arm. *(Recommendation: (b) unless a corridor-wide learning claim is needed.)*
2. **Arm A (running-state).** Designed, not run. Run it or drop it from scope.
3. **Publication scope.** If CA3-DG is the headline result, decide whether the analysis was
   pre-registered to that pair or is a post-hoc discovery (affects how the 8-pair correction
   is framed).

---

## 6. Reproduce

No orchestrator (no Snakemake/Make) — scripts run manually. Tests: `cd cca && PYTHONPATH=src
python -m pytest -q`. Browse all figures via `figures/index.html`.

```
# Landmark arm (canonical) — regenerate one config's figures + learning test
PYTHONPATH=src python scripts/regen_all_landmark_figs.sh          # or per-config plot_landmark.py
PYTHONPATH=src python scripts/learning_changes.py --tag landmark50_res_samp15
PYTHONPATH=src python scripts/summarise_landmark_sweep.py
PYTHONPATH=src python scripts/build_prune_table.py                # overfit table
PYTHONPATH=src python scripts/build_landmark_index.py             # index.html

# Spatial arm
PYTHONPATH=src python scripts/run_stage2.py --max-seconds 1800    # resumable
PYTHONPATH=src python scripts/run_stage3.py
PYTHONPATH=src python scripts/summarise_sweep.py
```
