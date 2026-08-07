# Item 4 — integration window vs IFI, naive vs experienced

Frozen axes (one CCA per animal-pair on all trials, each epoch projected through identical weights), so the same component is compared across epochs and windows. IFI integrated over |lag| ≤ w for w = 10…250 ms.

> Refit-per-lag was rejected for this contrast: the best-matching dimension > at another lag is a different dimension 82 % of the time, and split-half > |cos| of a canonical vector is 0.59 (CC1) falling to ~0.2 by CC4 against a > 0.146 random baseline — so a refit contrast would compare different > components across the conditions being contrasted.

> Correlations are IN-SAMPLE by construction; these are contrast statistics.

### FS-excluded — CC1 IFI, expert vs naive, at ±50 ms

Frozen axes; animals-as-n. Positive IFI ⇒ the first-named area leads.

| pair | n | IFI naive | IFI expert | Δ | p |
|---|---|---|---|---|---|
| CA1-RSC | 8 | -0.002 | +0.034 | +0.036 | 0.766 |
| CA1-CA3 | 6 | +0.236 | -0.023 | -0.259 | 0.261 |
| CA1-DG | 8 | -0.105 | -0.123 | -0.018 | 0.896 |
| CA1-V1 | 10 | -0.037 | -0.002 | +0.036 | 0.275 |
| CA3-DG | 4 | -0.166 | -0.017 | +0.149 | 0.419 |
| CA1-SUB | 4 | -0.090 | +0.045 | +0.134 | 0.537 |
| RSC-SUB | 4 | -0.035 | +0.005 | +0.039 | 0.0905 |
| V1-RSC | 6 | +0.025 | +0.010 | -0.015 | 0.0163 |

**Swept across all 25 windows** (the smallest p per pair — a SELECTED minimum over 25 nested, highly correlated windows, so it is a pointer, not a test):

| pair | best window | Δ | p (uncorrected) |
|---|---|---|---|
| CA1-CA3 | ±100 ms | -0.266 | 0.213 |
| CA1-DG | ±130 ms | +0.069 | 0.149 |
| CA1-RSC | ±20 ms | +0.113 | 0.363 |
| CA1-SUB | ±160 ms | +0.034 | 0.177 |
| CA1-V1 | ±50 ms | +0.036 | 0.275 |
| CA3-DG | ±240 ms | +0.088 | 0.0808 |
| RSC-SUB | ±90 ms | +0.061 | 0.0723 |
| V1-RSC | ±50 ms | -0.015 | 0.0163 |

### FS-included — CC1 IFI, expert vs naive, at ±50 ms

Frozen axes; animals-as-n. Positive IFI ⇒ the first-named area leads.

| pair | n | IFI naive | IFI expert | Δ | p |
|---|---|---|---|---|---|
| CA1-RSC | 8 | +0.148 | +0.171 | +0.023 | 0.728 |
| CA1-CA3 | 7 | +0.305 | +0.213 | -0.092 | 0.482 |
| CA1-DG | 8 | +0.177 | +0.012 | -0.165 | 0.0414 |
| CA1-V1 | 10 | +0.017 | -0.011 | -0.028 | 0.306 |
| CA3-DG | 4 | -0.177 | -0.131 | +0.045 | 0.745 |
| CA1-SUB | 4 | -0.290 | +0.072 | +0.362 | 0.264 |
| RSC-SUB | 4 | -0.109 | -0.030 | +0.079 | 0.22 |
| V1-RSC | 6 | +0.019 | +0.015 | -0.005 | 0.728 |

**Swept across all 25 windows** (the smallest p per pair — a SELECTED minimum over 25 nested, highly correlated windows, so it is a pointer, not a test):

| pair | best window | Δ | p (uncorrected) |
|---|---|---|---|
| CA1-CA3 | ±20 ms | -0.116 | 0.384 |
| CA1-DG | ±90 ms | -0.165 | 0.0297 |
| CA1-RSC | ±240 ms | +0.047 | 0.129 |
| CA1-SUB | ±10 ms | +0.192 | 0.174 |
| CA1-V1 | ±180 ms | +0.027 | 0.0591 |
| CA3-DG | ±250 ms | +0.056 | 0.181 |
| RSC-SUB | ±110 ms | +0.081 | 0.16 |
| V1-RSC | ±10 ms | -0.013 | 0.135 |
