# Item 1 — do different CCs have different IFI (different signs)?

Lag curves are the existing **CCA-refit-at-every-lag** held-out, segment-aware per-dimension curves (`lag_curves_bin10*.csv`, dims 1–30 × lags ±250 ms). This re-reduces them per dimension over integration windows |lag| ≤ w.

> Positive IFI ⇒ the **first-named area leads**. A dimension whose held-out > CC is negative at every lag clips to zero on both sides and returns IFI = 0 > for *no coupling* rather than *balanced*; those are flagged `degenerate` > and excluded from the sign counts.

> **`sig` consistency check (FS-excluded):** OK — 0/656 flagged dims with non-positive lag-0 CC.

#### FS-excluded — sign-mixing vs chance (significant CCs only)

Each animal contributes only the CCs that beat the circular-shift null. Chance is computed **per animal** from its own count *k* of significant CCs: P(mixed) = 1 − 2·0.5^k.

- animal-pairs with ≥2 significant CCs: **71**
- observed mixed: **64 (90%)**
- expected by chance: **70.2 (99%)**
- binomial test: **p = 1.52e-05**

**Verdict: LESS sign-mixing than independent signs would give — an animal's significant CCs share a direction more often than chance (consistent with a non-zero overall IFI), the opposite of 'different CCs, different signs'.**

#### FS-excluded — OVERALL IFI at ±50 ms (all significant CCs pooled, unweighted per-animal mean) vs 0

Per animal, the mean IFI over every CC that beats the per-dim null (not split by sign); then a one-sample *t* across animals, per pair. Positive ⇒ the first-named area leads. `w·` columns = the same with each CC weighted by its held-out peak. **CCs-as-n** = every significant CC as one sample pooled over animals (power check; CCs are nested in animals).

> **Data:** the refit-per-lag curves are UNCAPPED since 2026-08-17 — every running bin of every trial (whole session, ~370 k bins/animal). Until then `run_lag_curves` kept only the first ~20 trials (12k-bin cap); numbers before that date are not comparable. **n** = animals with ≥ 1 significant CC (of `N` recorded for the pair). **CC₁ (same data)** = dim 1 of this same table, every animal — the like-for-like comparator; it reproduces `bin10_tables.md` §B. Bold = p < 0.05 on this single look; the pair has 4 looks in all (FS-excl/incl × unweighted/weighted) — see the cross-FS summary below.

| pair | n / N | mean IFI | SEM | t | p | BH (8 pairs) | weighted mean | p (weighted) | CCs-as-n: n, mean, p | CC₁ (same data): mean, p, n |
|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 12 / 12 | **+0.053** | 0.015 | +3.51 | **0.00491** | yes | +0.063 | 0.000119 | 110, +0.050, **6.44e-06** | +0.081, p=0.000572, n=12 |
| CA1-CA3 | 7 / 7 | -0.027 | 0.019 | -1.44 | 0.199 | no | -0.011 | 0.576 | 60, -0.024, 0.138 | -0.002, p=0.929, n=7 |
| CA1-DG | 11 / 11 | **-0.053** | 0.018 | -2.98 | **0.0138** | no | -0.045 | 0.016 | 101, -0.053, **2.8e-05** | -0.023, p=0.296, n=11 |
| CA1-V1 | 13 / 13 | -0.010 | 0.008 | -1.24 | 0.237 | no | -0.006 | 0.192 | 123, -0.007, 0.397 | +0.002, p=0.769, n=13 |
| CA3-DG | 5 / 5 | **-0.061** | 0.018 | -3.42 | **0.0269** | no | -0.043 | 0.0239 | 46, -0.062, **0.00249** | -0.001, p=0.962, n=5 |
| CA1-SUB | 7 / 7 | -0.025 | 0.028 | -0.91 | 0.396 | no | -0.038 | 0.0999 | 63, -0.029, 0.0535 | -0.087, p=0.00854, n=7 |
| RSC-SUB | 7 / 7 | -0.005 | 0.015 | -0.30 | 0.772 | no | -0.006 | 0.691 | 66, -0.003, 0.76 | +0.020, p=0.402, n=7 |
| V1-RSC | 9 / 9 | **+0.038** | 0.014 | +2.79 | **0.0235** | no | +0.035 | 0.0234 | 87, +0.038, **6.49e-07** | +0.021, p=0.0404, n=9 |

### FS-excluded — do CCs within a pair disagree in sign?

Leading 5 canonical dims, integration window ±50 ms. `mixed` = that animal has at least one positive-IFI CC **and** at least one negative-IFI CC in the same pair. Positive IFI ⇒ the first-named area leads.

> ⚠ **Per-animal mixing is NOT evidence on its own.** With 5 dimensions and a coin-flip sign, at least one of each is expected 94% of the time. Read the reliable-direction table below instead — mixing only means something if the individual directions reproduce across animals.

| pair | animals | mixed-sign animals | mean CCs leading | mean CCs lagging | mean frac. positive |
|---|---|---|---|---|---|
| CA1-RSC | 12 | **7/12** (58%) | 4.1 | 0.9 | 0.82 |
| CA1-CA3 | 7 | **6/7** (86%) | 2.7 | 2.3 | 0.54 |
| CA1-DG | 11 | **10/11** (91%) | 1.8 | 3.2 | 0.36 |
| CA1-V1 | 13 | **13/13** (100%) | 2.3 | 2.7 | 0.46 |
| CA3-DG | 5 | **5/5** (100%) | 1.8 | 3.2 | 0.36 |
| CA1-SUB | 7 | **6/7** (86%) | 2.1 | 2.9 | 0.43 |
| RSC-SUB | 7 | **6/7** (86%) | 2.0 | 3.0 | 0.40 |
| V1-RSC | 9 | **5/9** (56%) | 4.1 | 0.9 | 0.82 |

**Mixed-sign rate vs integration window** (leading 5 dims, pooled over pairs):

| window (ms) | 10 | 20 | 50 | 100 | 150 | 200 | 250 |
|---|---|---|---|---|---|---|---|
| mixed-sign animal-pairs | 90% | 85% | 82% | 80% | 82% | 85% | 82% |

#### FS-excluded — which CCs have a RELIABLE direction across animals?

One-sample *t* of each dimension's IFI across animals (±50 ms window), BH-FDR across dims within a pair. A pair genuinely contains CCs of opposite direction only if **two dims survive with opposite signs**.

| pair | dims with a reliable direction (FDR) | opposite signs? |
|---|---|---|
| CA1-RSC | CC1 +0.081 (p=0.000572) | no |
| CA1-CA3 | — none | no |
| CA1-DG | — none | no |
| CA1-V1 | — none | no |
| CA3-DG | — none | no |
| CA1-SUB | — none | no |
| RSC-SUB | — none | no |
| V1-RSC | — none | no |

#### FS-excluded — strongest per-CC directions across the window sweep

⚠ **Selected minima over ~150 tests per pair (25 windows × 6 dims) — descriptive pointers, not inference.** 1 dimension(s) survive the fixed-window FDR test above.

| pair | CC | window | n | mean IFI | p (uncorrected) |
|---|---|---|---|---|---|
| CA1-RSC | CC1 | ±70 ms | 12 | +0.075 | 0.0003 |
| CA1-RSC | CC1 | ±210 ms | 12 | +0.033 | 0.0003 |
| CA1-RSC | CC1 | ±200 ms | 12 | +0.034 | 0.0003 |
| CA1-RSC | CC1 | ±60 ms | 12 | +0.080 | 0.0003 |
| CA1-RSC | CC1 | ±120 ms | 12 | +0.039 | 0.0003 |
| CA1-RSC | CC1 | ±220 ms | 12 | +0.032 | 0.0003 |
| CA1-RSC | CC1 | ±130 ms | 12 | +0.037 | 0.0004 |
| CA1-RSC | CC1 | ±190 ms | 12 | +0.035 | 0.0004 |

**Pairs where two CCs of OPPOSITE sign are both nominally significant at the same window:**

- none

> **`sig` consistency check (FS-included):** OK — 0/701 flagged dims with non-positive lag-0 CC.

#### FS-included — sign-mixing vs chance (significant CCs only)

Each animal contributes only the CCs that beat the circular-shift null. Chance is computed **per animal** from its own count *k* of significant CCs: P(mixed) = 1 − 2·0.5^k.

- animal-pairs with ≥2 significant CCs: **73**
- observed mixed: **69 (95%)**
- expected by chance: **72.6 (99%)**
- binomial test: **p = 0.000716**

**Verdict: LESS sign-mixing than independent signs would give — an animal's significant CCs share a direction more often than chance (consistent with a non-zero overall IFI), the opposite of 'different CCs, different signs'.**

#### FS-included — OVERALL IFI at ±50 ms (all significant CCs pooled, unweighted per-animal mean) vs 0

Per animal, the mean IFI over every CC that beats the per-dim null (not split by sign); then a one-sample *t* across animals, per pair. Positive ⇒ the first-named area leads. `w·` columns = the same with each CC weighted by its held-out peak. **CCs-as-n** = every significant CC as one sample pooled over animals (power check; CCs are nested in animals).

> **Data:** the refit-per-lag curves are UNCAPPED since 2026-08-17 — every running bin of every trial (whole session, ~370 k bins/animal). Until then `run_lag_curves` kept only the first ~20 trials (12k-bin cap); numbers before that date are not comparable. **n** = animals with ≥ 1 significant CC (of `N` recorded for the pair). **CC₁ (same data)** = dim 1 of this same table, every animal — the like-for-like comparator; it reproduces `bin10_tables.md` §B. Bold = p < 0.05 on this single look; the pair has 4 looks in all (FS-excl/incl × unweighted/weighted) — see the cross-FS summary below.

| pair | n / N | mean IFI | SEM | t | p | BH (8 pairs) | weighted mean | p (weighted) | CCs-as-n: n, mean, p | CC₁ (same data): mean, p, n |
|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 12 / 12 | **+0.068** | 0.011 | +6.08 | **7.96e-05** | yes | +0.082 | 4.14e-05 | 118, +0.067, **2.52e-12** | +0.140, p=0.000613, n=12 |
| CA1-CA3 | 9 / 9 | +0.001 | 0.013 | +0.04 | 0.97 | no | +0.001 | 0.923 | 84, -0.001, 0.945 | +0.012, p=0.11, n=9 |
| CA1-DG | 11 / 11 | **-0.044** | 0.016 | -2.78 | **0.0195** | no | -0.035 | 0.00896 | 104, -0.048, **0.000248** | -0.012, p=0.545, n=11 |
| CA1-V1 | 13 / 13 | -0.010 | 0.009 | -1.08 | 0.3 | no | -0.006 | 0.193 | 125, -0.007, 0.316 | -0.006, p=0.186, n=13 |
| CA3-DG | 5 / 5 | **-0.082** | 0.019 | -4.29 | **0.0128** | no | -0.042 | 0.0567 | 50, -0.082, **1.66e-05** | -0.027, p=0.4, n=5 |
| CA1-SUB | 7 / 7 | +0.008 | 0.048 | +0.17 | 0.87 | no | -0.012 | 0.732 | 65, -0.010, 0.617 | -0.052, p=0.0601, n=7 |
| RSC-SUB | 7 / 7 | -0.013 | 0.013 | -1.03 | 0.344 | no | -0.011 | 0.438 | 65, -0.012, 0.291 | +0.001, p=0.964, n=7 |
| V1-RSC | 9 / 9 | +0.030 | 0.013 | +2.30 | 0.0502 | no | +0.032 | 0.0171 | 90, +0.030, **0.000867** | +0.024, p=0.00958, n=9 |

### FS-included — do CCs within a pair disagree in sign?

Leading 5 canonical dims, integration window ±50 ms. `mixed` = that animal has at least one positive-IFI CC **and** at least one negative-IFI CC in the same pair. Positive IFI ⇒ the first-named area leads.

> ⚠ **Per-animal mixing is NOT evidence on its own.** With 5 dimensions and a coin-flip sign, at least one of each is expected 94% of the time. Read the reliable-direction table below instead — mixing only means something if the individual directions reproduce across animals.

| pair | animals | mixed-sign animals | mean CCs leading | mean CCs lagging | mean frac. positive |
|---|---|---|---|---|---|
| CA1-RSC | 12 | **8/12** (67%) | 4.2 | 0.8 | 0.85 |
| CA1-CA3 | 9 | **9/9** (100%) | 2.8 | 2.2 | 0.56 |
| CA1-DG | 11 | **7/11** (64%) | 1.0 | 4.0 | 0.20 |
| CA1-V1 | 13 | **12/13** (92%) | 2.2 | 2.8 | 0.43 |
| CA3-DG | 5 | **5/5** (100%) | 1.6 | 3.4 | 0.32 |
| CA1-SUB | 7 | **6/7** (86%) | 2.0 | 3.0 | 0.40 |
| RSC-SUB | 7 | **6/7** (86%) | 2.3 | 2.7 | 0.46 |
| V1-RSC | 9 | **4/9** (44%) | 4.6 | 0.4 | 0.91 |

**Mixed-sign rate vs integration window** (leading 5 dims, pooled over pairs):

| window (ms) | 10 | 20 | 50 | 100 | 150 | 200 | 250 |
|---|---|---|---|---|---|---|---|
| mixed-sign animal-pairs | 88% | 84% | 78% | 78% | 70% | 78% | 81% |

#### FS-included — which CCs have a RELIABLE direction across animals?

One-sample *t* of each dimension's IFI across animals (±50 ms window), BH-FDR across dims within a pair. A pair genuinely contains CCs of opposite direction only if **two dims survive with opposite signs**.

| pair | dims with a reliable direction (FDR) | opposite signs? |
|---|---|---|
| CA1-RSC | CC1 +0.140 (p=0.000613) | no |
| CA1-CA3 | CC11 -0.132 (p=0.00167) | no |
| CA1-DG | — none | no |
| CA1-V1 | — none | no |
| CA3-DG | — none | no |
| CA1-SUB | — none | no |
| RSC-SUB | — none | no |
| V1-RSC | — none | no |

#### FS-included — strongest per-CC directions across the window sweep

⚠ **Selected minima over ~150 tests per pair (25 windows × 6 dims) — descriptive pointers, not inference.** 2 dimension(s) survive the fixed-window FDR test above.

| pair | CC | window | n | mean IFI | p (uncorrected) |
|---|---|---|---|---|---|
| V1-RSC | CC1 | ±190 ms | 9 | +0.015 | 0.0001 |
| V1-RSC | CC1 | ±180 ms | 9 | +0.015 | 0.0002 |
| V1-RSC | CC1 | ±200 ms | 9 | +0.015 | 0.0002 |
| V1-RSC | CC1 | ±170 ms | 9 | +0.015 | 0.0002 |
| V1-RSC | CC1 | ±160 ms | 9 | +0.016 | 0.0003 |
| V1-RSC | CC1 | ±210 ms | 9 | +0.015 | 0.0003 |
| CA1-RSC | CC3 | ±120 ms | 12 | +0.038 | 0.0004 |
| V1-RSC | CC1 | ±150 ms | 9 | +0.016 | 0.0004 |

**Pairs where two CCs of OPPOSITE sign are both nominally significant at the same window:**

- **CA1-V1** at ±20 ms: CC3 +0.020 (p=0.030), CC4 -0.024 (p=0.049)
- **RSC-SUB** at ±170 ms: CC4 -0.042 (p=0.034), CC6 +0.046 (p=0.044)
- **RSC-SUB** at ±180 ms: CC4 -0.043 (p=0.029), CC6 +0.049 (p=0.045)
- **RSC-SUB** at ±190 ms: CC4 -0.042 (p=0.025), CC6 +0.048 (p=0.037)
- **RSC-SUB** at ±200 ms: CC4 -0.043 (p=0.024), CC6 +0.046 (p=0.024)
- **RSC-SUB** at ±210 ms: CC4 -0.046 (p=0.023), CC6 +0.035 (p=0.006)
- **RSC-SUB** at ±220 ms: CC4 -0.047 (p=0.020), CC6 +0.027 (p=0.022)
- **RSC-SUB** at ±230 ms: CC4 -0.047 (p=0.020), CC6 +0.030 (p=0.012)
- **RSC-SUB** at ±240 ms: CC4 -0.046 (p=0.025), CC6 +0.039 (p=0.006)
- **RSC-SUB** at ±250 ms: CC4 -0.044 (p=0.043), CC6 +0.043 (p=0.008)

#### OVERALL IFI at ±50 ms — the four looks per pair (FS-excluded / FS-included × unweighted / weighted)

Nominal p < 0.05 counts; a direction is only worth stating when the sign is the same in all four and more than one look is nominal.

| pair | looks p<0.05 (of 4) | signs (excl-u, excl-w, incl-u, incl-w) | means |
|---|---|---|---|
| CA1-RSC | 4/4 | + + + + | +0.053 / +0.063 / +0.068 / +0.082 |
| CA1-CA3 | 0/4 | − − + + | -0.027 / -0.011 / +0.001 / +0.001 |
| CA1-DG | 4/4 | − − − − | -0.053 / -0.045 / -0.044 / -0.035 |
| CA1-V1 | 0/4 | − − − − | -0.010 / -0.006 / -0.010 / -0.006 |
| CA3-DG | 3/4 | − − − − | -0.061 / -0.043 / -0.082 / -0.042 |
| CA1-SUB | 0/4 | − − + − | -0.025 / -0.038 / +0.008 / -0.012 |
| RSC-SUB | 0/4 | − − − − | -0.005 / -0.006 / -0.013 / -0.011 |
| V1-RSC | 3/4 | + + + + | +0.038 / +0.035 / +0.030 / +0.032 |
