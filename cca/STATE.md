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
| Learning test | **paired Wilcoxon + FDR** (`learning_changes_spatial_*.csv`) — null; CA3-DG dir-consistent | **paired per-animal Wilcoxon + FDR** (`results/learning_changes_*.csv`) |
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

> **⚠ Pseudoreplication caveat (2026-06-05) — read before trusting the p-values below.**
> The `per_pair_pooled` test in `learning_changes.py:174` pools deltas across all 6 landmarks
> (`n` = animals × landmarks) and treats them as independent — but the 6 landmarks within an
> animal are repeated measures, not independent samples. This inflates significance: several
> "hits" have a **median delta of ≈0.000** yet `p<0.01` purely from the inflated `n` (e.g.
> CA1-RSC in `samp25`). When landmarks are **collapsed to one value per animal** (the honest
> unit, `n` = #animals = 4–10), **no pair survives within-pair correction in the committed
> config — including CA3-DG** (n=4, Wilcoxon p-floor 0.125). So the headline numbers below are
> the pooled (anti-conservative) figures; the honest per-animal verdict is **null**. CA3-DG
> remains the most consistent *direction* (positive Δ in nearly all configs and 63/66 spatial
> configs) but is **not significant under a proper per-animal / mixed-effects test**. The real
> limit is N and pseudoreplication, not the FDR family. A mixed-effects model (animal random
> effect, landmark as repeated measure) is the correct way to use the landmark-level data.

**Honest verdict (mixed-effects / per-animal tests, 2026-06-05).** Two non-pseudoreplicated
tests — per-animal collapse and a random-slope LMM (`scripts/learning_changes_mixed.py`,
`src/tom_cca/mixed_effects.py`), within-pair FDR, no cross-pair correction — reshape the
headline. Effects appear in the **well-powered pairs (n≈8–10 animals)**, not the pooled favourite:
- **CA1-RSC** — most reproducible: a directionality (`ifi_weighted`) shift that peaks at the
  *intermediate* stage (non-monotonic), in ~20 configs across both methods; plus a strength
  (`mi_sig`) increase in `landmark50_res_samp25` (both methods). n=8.
- **CA1-DG** — modest monotonic strength increase (expert>naive, +0.05) surviving the LMM in
  `landmark50_res_samp40` (both methods). n=8.
- **CA3-DG** — the pooled "headline" — **does not survive** any honest test (only n=4 animals;
  signed-rank floor p=0.125), though all 4 animals move in the same (+) direction.
- The large V1-RSC `mi_sig` hits are confined to overfit/high-k 25 ms configs and are
  non-monotonic — discount.
Bottom line: with honest statistics the learning effects are **modest, pair-specific, and
config-dependent**; the N (4–10 animals/pair) is the binding constraint. Most defensible single
result: **CA1-DG strength ↑ in `landmark50_res_samp40`** (residual, conservative k, survives the
random-slope LMM).

**Headline (pooled scope — anti-conservative, see caveat above):** the most consistent effect is
**CA3-DG communication strengthening with expertise** — direction-robust, but its significance is
pseudoreplicated and it does **not** survive the honest tests above.

- **Canonical landmark config (`landmark50_res_samp15`), pooled across landmarks:**
  - expert vs naive: median Δmi = **+0.21**, p = 0.0059, **FDR-pass** ✓
  - expert vs intermediate: median Δmi = **+0.15**, p = 0.0043, **FDR-pass** ✓
  - intermediate vs naive: Δmi ≈ 0, p = 0.65, n.s.
  - → monotonic, with the jump occurring **at the expert stage**.
- **Reproducibility:** CA3-DG expert-vs-naive survives FDR in **18 of 40** non-overfit
  landmark configs — far ahead of any other effect (next is V1-RSC intermediate-naive at
  7/40, which does *not* reproduce in the expert contrast). Source: pooled-scope rows of
  `results/learning_changes_*.csv` (CRLF-safe scan; see `GOTCHAS.md`).
- **Spatial corroboration (paired test, 2026-06-05):** the spatial arm now has its own paired
  per-animal Wilcoxon + FDR test (`scripts/learning_changes_spatial.py`, mirroring the landmark
  pooled scope). Across all 66 spatial configs there is **1** FDR survivor total
  (CA1-V1/ifi_weighted in a single config) — consistent with chance, i.e. **no robust spatial
  learning effect**. But **CA3-DG expert−naive `mi_sig` is positive in 63/66 configs** (committed
  config: Δ=+0.18) — strongly direction-consistent with the landmark finding. It never reaches
  significance because only **n=4** animals have CA3-DG, and the two-sided signed-rank p-floor at
  n=4 is 0.125. So the full-corridor arm is both diluted and underpowered; the event-locked
  landmark arm is where the effect is detectable.
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

1. ~~**Spatial-arm learning gap.**~~ **Resolved 2026-06-05** — built `learning_changes_spatial.py`
   (paired per-animal Wilcoxon + FDR). Verdict: no robust spatial learning effect (1/66 configs
   has a lone FDR hit), CA3-DG direction-consistent but underpowered (n=4). The spatial arm is
   confirmed exploratory/corroborative; the landmark arm is canonical.
2. **Arm A (running-state).** Committed config `temp50_sig_samp15` run on 2026-06-05 (see §3 once
   analysed). Full 22-config sweep still optional.
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
