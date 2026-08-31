# Trial 1 vs 2, per-unit carrying — summary (2026-08-31)

Driver `run_trial12_units.py`: carry_r(unit, ordinal) = correlation of the unit's
residualised activity with the PARTNER area's frozen CC1 variate over the ordinal's running
bins (subspace fit on ordinals 11+, ordinals 1..10 leak-free; same design as run_trial12).
Analyzer `analyze_trial12_units.py`; per-cell tests in `trial12_units_tests_bin10*.csv`,
prints in `trial12_units_fs*.txt`; figure `HCV1_trial12_units_deltas_*`. Primary arm
matched=1 (all ordinals cut to a common bin count).

**1. Carrying strength (Δ mean |r|, trial 1 − 2): NULL.** All 16 cells n.s., both FS; the
1→2 step is not special vs the adjacent-step band either. GLOBAL per-animal: n.s.

**2. Membership REORDERS 1→2 (the finding).** sim(1,2) of the per-unit carrying profile sits
BELOW the adjacent-step similarity band in 14/16 cells (both FS). No single cell survives BH;
the honest global test (one value per animal, mean over its cells, n = 16) is significant in
both FS on the primary arm: **FS-excl med −0.146, W p = 3.1×10⁻⁵ (t 1.8×10⁻⁴); FS-incl med
−0.104, W p = 0.021 (t 0.079)**. Raw arm: FS-excl W p = 0.0042, FS-incl n.s. — the matched arm
is primary (removes the correlation-vs-n confound). Behaviour control: the per-animal excess
is uncorrelated with Δspeed (ρ = −0.07/+0.09) and Δbins (ρ = +0.35/+0.11, n.s.).

**3. Convergence to the trained profile (Δ sim-to-reference, 1 − 2): null** per cell and
globally (med −0.05/−0.03, W p = 0.30/0.23) — trial 1's profile is not measurably further
from the eventual membership than trial 2's at this n.

**Reading.** Between the very first and second traversal, the CHANNEL is unchanged (pair-level
strength, direction — 2026-08-20 — and per-unit carrying strength all null), but WHICH units
express it shifts more than between any later adjacent trials. Localisation to specific pairs
is underpowered. Caveat: profile similarity is estimated from ~200–2900 bins/trial; the
common-bin arm equalises that within animal, and the global test is animals-as-n.
