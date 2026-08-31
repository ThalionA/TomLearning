# V1-RSC trial-1-vs-2 deep-dive + naive communication vs learning speed (2026-08-31)

Scripts `figs_trial12_v1rsc.py` (reuses analyze_trial12's exact chain) and
`analyze_conn_learning.py` / `figs_conn_learning.py`. Full prints:
`trial12_v1rsc_sensitivity.txt`, `conn_learning_fs*.txt`.

## V1-RSC (figure `HCV1_trial12_v1rsc_*`)

- **The robust picture is a LEVEL effect:** trial 1's frozen-CC1 lag curve is depressed at
  ALL lags (~60-70 % of the trials-3..10 band level) in BOTH FS; trial 2 already sits on the
  band. Not a peak-shape/direction reorganisation.
- **The paired tests stay weak at n = 9:** Δr0(1−2) is negative in all 6 arm × FS
  specifications (−0.015…−0.023) but never W p < 0.05 (best 0.055).
- **The direction (IFI) delta is FS- and arm-fragile:** FS-excl arms 1/2 starred (W 0.039 /
  0.020), FS-incl all n.s. Side decomposition: FS-incl the V1→RSC (positive-lag) mass drops
  (t 0.038, W 0.039); FS-excl both sides drop, n.s.
- **Reading:** on the very first traversal V1-RSC communication runs at reduced gain, with the
  V1-lead component possibly absent (FS-excl IFI ≈ 0 on trial 1, ≈ +0.05 from trial 2 on) —
  sign-consistent in every specification, significant in none robustly. n = 9; a
  second-cohort question, not a claim.

## Naive communication vs learning speed (figure `HCV1_connlearning_*`)

Exploratory; predictors = naive-epoch (first 10 trials) cc1 / IFI per pair, outcome = learning
point (n = 12 learners, LP 28-107; per pair 5-10). **No cell survives BH; global naive
coupling (within-pair z of cc1, animal mean) ρ = −0.24/−0.06, n.s.; IFI predicts nothing.**
Only repeatable lean: **CA1-DG strength ρ = −0.69, p = 0.058 in both FS** (stronger naive
CA1↔DG coupling → earlier LP) — leverage-sensitive (without animal 75: ρ = −0.57/−0.61,
p = 0.18/0.15). Confound check: naive running amount vs LP ρ = +0.39, n.s. Verdict: at n = 12
this is a null with one direction-consistent hippocampal lean; log for the second cohort.
