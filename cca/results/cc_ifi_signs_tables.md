# Item 1 — do different CCs have different IFI (different signs)?

Lag curves are the existing **CCA-refit-at-every-lag** held-out, segment-aware per-dimension curves (`lag_curves_bin10*.csv`, dims 1–30 × lags ±250 ms). This re-reduces them per dimension over integration windows |lag| ≤ w.

> Positive IFI ⇒ the **first-named area leads**. A dimension whose held-out > CC is negative at every lag clips to zero on both sides and returns IFI = 0 > for *no coupling* rather than *balanced*; those are flagged `degenerate` > and excluded from the sign counts.

> **`sig` consistency check (FS-excluded):** OK — 0/238 flagged dims with non-positive lag-0 CC.

#### FS-excluded — sign-mixing vs chance (significant CCs only)

Each animal contributes only the CCs that beat the circular-shift null. Chance is computed **per animal** from its own count *k* of significant CCs: P(mixed) = 1 − 2·0.5^k.

- animal-pairs with ≥2 significant CCs: **49**
- observed mixed: **39 (80%)**
- expected by chance: **40.3 (82%)**
- binomial test: **p = 0.579**

**Verdict: sign-mixing is indistinguishable from a coin flip.**

#### FS-excluded — OVERALL IFI at ±50 ms (all significant CCs pooled, unweighted per-animal mean) vs 0

Per animal, the mean IFI over every CC that beats the per-dim null (not split by sign); then a one-sample *t* across animals, per pair. Positive ⇒ the first-named area leads. `w·` columns = the same with each CC weighted by its held-out peak. Compare with `bin10_tables.md` §B, which is the same question for CC₁ alone.

| pair | n | mean IFI | SEM | t | p | BH (8 pairs) | weighted mean | p (weighted) |
|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | **+0.226** | 0.090 | +2.50 | **0.0409** | no | +0.218 | 0.0503 |
| CA1-CA3 | 6 | -0.095 | 0.043 | -2.20 | 0.0788 | no | -0.091 | 0.0408 |
| CA1-DG | 9 | **-0.162** | 0.066 | -2.46 | **0.0395** | no | -0.143 | 0.0674 |
| CA1-V1 | 7 | -0.111 | 0.050 | -2.22 | 0.0685 | no | -0.100 | 0.0375 |
| CA3-DG | 5 | -0.053 | 0.056 | -0.96 | 0.391 | no | -0.033 | 0.578 |
| CA1-SUB | 6 | -0.041 | 0.044 | -0.94 | 0.391 | no | -0.075 | 0.142 |
| RSC-SUB | 7 | +0.003 | 0.060 | +0.05 | 0.964 | no | +0.007 | 0.91 |
| V1-RSC | 8 | +0.089 | 0.057 | +1.56 | 0.163 | no | +0.116 | 0.0589 |

### FS-excluded — do CCs within a pair disagree in sign?

Leading 5 canonical dims, integration window ±50 ms. `mixed` = that animal has at least one positive-IFI CC **and** at least one negative-IFI CC in the same pair. Positive IFI ⇒ the first-named area leads.

> ⚠ **Per-animal mixing is NOT evidence on its own.** With 5 dimensions and a coin-flip sign, at least one of each is expected 94% of the time. Read the reliable-direction table below instead — mixing only means something if the individual directions reproduce across animals.

| pair | animals | mixed-sign animals | mean CCs leading | mean CCs lagging | mean frac. positive |
|---|---|---|---|---|---|
| CA1-RSC | 12 | **11/12** (92%) | 3.3 | 1.7 | 0.67 |
| CA1-CA3 | 7 | **7/7** (100%) | 1.9 | 3.1 | 0.37 |
| CA1-DG | 11 | **10/11** (91%) | 2.0 | 3.0 | 0.40 |
| CA1-V1 | 13 | **12/13** (92%) | 1.9 | 3.1 | 0.38 |
| CA3-DG | 5 | **3/5** (60%) | 1.4 | 3.6 | 0.28 |
| CA1-SUB | 7 | **7/7** (100%) | 1.9 | 3.1 | 0.37 |
| RSC-SUB | 7 | **7/7** (100%) | 3.0 | 2.0 | 0.60 |
| V1-RSC | 9 | **8/9** (89%) | 3.6 | 1.4 | 0.71 |

**Mixed-sign rate vs integration window** (leading 5 dims, pooled over pairs):

| window (ms) | 10 | 20 | 50 | 100 | 150 | 200 | 250 |
|---|---|---|---|---|---|---|---|
| mixed-sign animal-pairs | 87% | 90% | 92% | 90% | 86% | 80% | 79% |

#### FS-excluded — which CCs have a RELIABLE direction across animals?

One-sample *t* of each dimension's IFI across animals (±50 ms window), BH-FDR across dims within a pair. A pair genuinely contains CCs of opposite direction only if **two dims survive with opposite signs**.

| pair | dims with a reliable direction (FDR) | opposite signs? |
|---|---|---|
| CA1-RSC | — none | no |
| CA1-CA3 | — none | no |
| CA1-DG | — none | no |
| CA1-V1 | — none | no |
| CA3-DG | — none | no |
| CA1-SUB | — none | no |
| RSC-SUB | — none | no |
| V1-RSC | — none | no |

#### FS-excluded — strongest per-CC directions across the window sweep

⚠ **Selected minima over ~150 tests per pair (25 windows × 6 dims) — descriptive pointers, not inference.** Nothing here survives the fixed-window FDR test above.

| pair | CC | window | n | mean IFI | p (uncorrected) |
|---|---|---|---|---|---|
| CA3-DG | CC3 | ±110 ms | 5 | -0.152 | 0.0008 |
| CA1-SUB | CC6 | ±240 ms | 7 | -0.201 | 0.0011 |
| CA3-DG | CC3 | ±100 ms | 5 | -0.182 | 0.0014 |
| CA1-DG | CC5 | ±10 ms | 11 | -0.557 | 0.0015 |
| CA3-DG | CC3 | ±90 ms | 5 | -0.186 | 0.0023 |
| CA1-SUB | CC6 | ±250 ms | 7 | -0.158 | 0.0025 |
| CA3-DG | CC3 | ±120 ms | 5 | -0.159 | 0.0041 |
| CA1-SUB | CC6 | ±210 ms | 7 | -0.203 | 0.0046 |

**Pairs where two CCs of OPPOSITE sign are both nominally significant at the same window:**

- **CA3-DG** at ±230 ms: CC2 +0.126 (p=0.042), CC3 -0.109 (p=0.025)

> **`sig` consistency check (FS-included):** OK — 0/307 flagged dims with non-positive lag-0 CC.

#### FS-included — sign-mixing vs chance (significant CCs only)

Each animal contributes only the CCs that beat the circular-shift null. Chance is computed **per animal** from its own count *k* of significant CCs: P(mixed) = 1 − 2·0.5^k.

- animal-pairs with ≥2 significant CCs: **59**
- observed mixed: **46 (78%)**
- expected by chance: **50.8 (86%)**
- binomial test: **p = 0.087**

**Verdict: sign-mixing is indistinguishable from a coin flip.**

#### FS-included — OVERALL IFI at ±50 ms (all significant CCs pooled, unweighted per-animal mean) vs 0

Per animal, the mean IFI over every CC that beats the per-dim null (not split by sign); then a one-sample *t* across animals, per pair. Positive ⇒ the first-named area leads. `w·` columns = the same with each CC weighted by its held-out peak. Compare with `bin10_tables.md` §B, which is the same question for CC₁ alone.

| pair | n | mean IFI | SEM | t | p | BH (8 pairs) | weighted mean | p (weighted) |
|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 12 | +0.157 | 0.072 | +2.20 | 0.0501 | no | +0.169 | 0.0332 |
| CA1-CA3 | 9 | -0.047 | 0.045 | -1.06 | 0.321 | no | -0.040 | 0.185 |
| CA1-DG | 11 | -0.119 | 0.068 | -1.75 | 0.11 | no | -0.140 | 0.0365 |
| CA1-V1 | 8 | -0.201 | 0.114 | -1.76 | 0.121 | no | -0.197 | 0.132 |
| CA3-DG | 5 | -0.058 | 0.068 | -0.86 | 0.44 | no | -0.043 | 0.443 |
| CA1-SUB | 6 | -0.051 | 0.032 | -1.60 | 0.171 | no | -0.112 | 0.0266 |
| RSC-SUB | 7 | -0.019 | 0.045 | -0.41 | 0.693 | no | +0.009 | 0.88 |
| V1-RSC | 9 | **+0.120** | 0.041 | +2.94 | **0.0188** | no | +0.108 | 0.0273 |

### FS-included — do CCs within a pair disagree in sign?

Leading 5 canonical dims, integration window ±50 ms. `mixed` = that animal has at least one positive-IFI CC **and** at least one negative-IFI CC in the same pair. Positive IFI ⇒ the first-named area leads.

> ⚠ **Per-animal mixing is NOT evidence on its own.** With 5 dimensions and a coin-flip sign, at least one of each is expected 94% of the time. Read the reliable-direction table below instead — mixing only means something if the individual directions reproduce across animals.

| pair | animals | mixed-sign animals | mean CCs leading | mean CCs lagging | mean frac. positive |
|---|---|---|---|---|---|
| CA1-RSC | 12 | **10/12** (83%) | 3.2 | 1.8 | 0.65 |
| CA1-CA3 | 9 | **9/9** (100%) | 2.2 | 2.8 | 0.44 |
| CA1-DG | 11 | **9/11** (82%) | 2.0 | 3.0 | 0.40 |
| CA1-V1 | 13 | **11/13** (85%) | 1.7 | 3.3 | 0.34 |
| CA3-DG | 5 | **2/5** (40%) | 2.2 | 2.8 | 0.44 |
| CA1-SUB | 7 | **6/7** (86%) | 1.7 | 3.3 | 0.34 |
| RSC-SUB | 7 | **7/7** (100%) | 2.1 | 2.9 | 0.43 |
| V1-RSC | 9 | **6/9** (67%) | 3.7 | 1.3 | 0.73 |

**Mixed-sign rate vs integration window** (leading 5 dims, pooled over pairs):

| window (ms) | 10 | 20 | 50 | 100 | 150 | 200 | 250 |
|---|---|---|---|---|---|---|---|
| mixed-sign animal-pairs | 89% | 88% | 82% | 86% | 79% | 81% | 79% |

#### FS-included — which CCs have a RELIABLE direction across animals?

One-sample *t* of each dimension's IFI across animals (±50 ms window), BH-FDR across dims within a pair. A pair genuinely contains CCs of opposite direction only if **two dims survive with opposite signs**.

| pair | dims with a reliable direction (FDR) | opposite signs? |
|---|---|---|
| CA1-RSC | — none | no |
| CA1-CA3 | — none | no |
| CA1-DG | — none | no |
| CA1-V1 | — none | no |
| CA3-DG | — none | no |
| CA1-SUB | — none | no |
| RSC-SUB | — none | no |
| V1-RSC | — none | no |

#### FS-included — strongest per-CC directions across the window sweep

⚠ **Selected minima over ~150 tests per pair (25 windows × 6 dims) — descriptive pointers, not inference.** Nothing here survives the fixed-window FDR test above.

| pair | CC | window | n | mean IFI | p (uncorrected) |
|---|---|---|---|---|---|
| V1-RSC | CC1 | ±10 ms | 9 | +0.118 | 0.0009 |
| CA1-DG | CC6 | ±10 ms | 10 | -0.635 | 0.0009 |
| CA1-V1 | CC5 | ±210 ms | 13 | -0.104 | 0.0020 |
| V1-RSC | CC1 | ±20 ms | 9 | +0.155 | 0.0021 |
| CA1-V1 | CC5 | ±220 ms | 13 | -0.093 | 0.0040 |
| V1-RSC | CC1 | ±80 ms | 9 | +0.097 | 0.0055 |
| V1-RSC | CC1 | ±90 ms | 9 | +0.090 | 0.0057 |
| CA1-V1 | CC5 | ±230 ms | 13 | -0.091 | 0.0061 |

**Pairs where two CCs of OPPOSITE sign are both nominally significant at the same window:**

- none
