# Fixed-subspace lag curves by epoch — meeting items 5, 6, 7

One canonical component, identified once across balanced trials and frozen, then lagged across time within each learning epoch. Because the weights are identical across epochs, an epoch difference is a difference in the **activity**, not in the fit.

> **Not a coupling strength.** The frozen fit saw every epoch, so `peak r` is > in-sample and optimistic. It is a contrast statistic only; the leak-free > numbers live in `lag_subspaces_tables.md`.

### FS-excluded — expert vs naive through a FROZEN subspace

Subspace identified once on trials balanced across epochs, then both epochs projected through the identical weights; animals-as-n paired *t*. `peak r` is in-sample by construction — read the Δ, not the level.

| pair | n | peak r naive | peak r exp | Δ | p | IFI naive | IFI exp | Δ | p | peak lag (ms) naive | peak lag (ms) exp | Δ | p | half-max width (ms) naive | half-max width (ms) exp | Δ | p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | 0.113 | 0.0954 | -0.0175 | 0.109 | 0.0839 | 0.00851 | -0.0754 | 0.429 | 33.8 | 18.8 | -15 | 0.709 | 31.2 | 20 | -11.2 | 0.416 |
| CA1-CA3 | 6 | 0.261 | 0.288 | +0.0263 | 0.538 | -0.03 | -0.0378 | -0.00779 | 0.892 | 0 | 0 | +0 | nan | 31.7 | 28.3 | -3.33 | 0.175 |
| CA1-DG | 8 | 0.196 | 0.201 | +0.00474 | 0.654 | -0.114 | 0.0124 | +0.127 | 0.134 | -2.5 | 0 | +2.5 | 0.563 | 22.5 | 21.2 | -1.25 | 0.732 |
| CA1-V1 | 10 | 0.135 | 0.127 | -0.00779 | 0.693 | -0.123 | 0.0206 | +0.144 | 0.102 | -11 | -2 | +9 | 0.287 | 180 | 257 | +77 | 0.306 |
| CA3-DG | 4 | 0.291 | 0.334 | +0.0428 | 0.255 | -0.0561 | 0.039 | +0.095 | 0.176 | 0 | 0 | +0 | nan | 25 | 22.5 | -2.5 | 0.718 |
| CA1-SUB | 4 | 0.109 | 0.117 | +0.0085 | 0.693 | -0.0706 | -0.0598 | +0.0109 | 0.887 | -35 | 2.5 | +37.5 | 0.427 | 7.5 | 25 | +17.5 | 0.133 |
| RSC-SUB | 4 | 0.216 | 0.199 | -0.0172 | 0.633 | 0.108 | 0.0906 | -0.0172 | 0.734 | 2.5 | 0 | -2.5 | 0.391 | 182 | 248 | +65 | 0.559 |
| V1-RSC | 6 | 0.25 | 0.241 | -0.00864 | 0.795 | 0.0591 | -0.0116 | -0.0707 | 0.226 | 20 | 0 | -20 | 0.275 | 235 | 235 | +0 | 1 |

#### FS-excluded — is the lag curve decaying or RINGING?

A secondary peak above half the central peak means the curve oscillates rather than decays. On a ringing curve the half-max width is the **half-period of the rhythm, not a temporal integration window**, and the IFI compares two lobes of that rhythm rather than a lead.

| pair | curves with a side peak | side-peak lag | implied freq | mean width | read the width as |
|---|---|---|---|---|---|
| CA1-CA3 | 12/18 (67%) | 140 ms | 7.1 Hz | 29 ms | **half-period of a rhythm** |
| CA1-DG | 17/24 (71%) | 136 ms | 7.3 Hz | 23 ms | **half-period of a rhythm** |
| CA1-RSC | 19/24 (79%) | 150 ms | 6.7 Hz | 34 ms | **half-period of a rhythm** |
| CA1-SUB | 9/12 (75%) | 116 ms | 8.7 Hz | 24 ms | **half-period of a rhythm** |
| CA1-V1 | 14/30 (47%) | 132 ms | 7.6 Hz | 192 ms | mixed — mostly integration window |
| CA3-DG | 6/12 (50%) | 133 ms | 7.5 Hz | 26 ms | **half-period of a rhythm** |
| RSC-SUB | 5/12 (42%) | 140 ms | 7.1 Hz | 198 ms | mixed — mostly integration window |
| V1-RSC | 4/18 (22%) | 138 ms | 7.3 Hz | 242 ms | mixed — mostly integration window |

#### FS-excluded — how censored is the integration window?

`width_ms` is bounded at both ends. Of 150 (animal, pair, epoch) values, **20% sit at the 0 ms floor** (only the peak bin clears half-max) and **9% at the 500 ms ceiling** (the curve never drops below half-max); only 71% are interior. A paired *t* on this metric is therefore fragile, and any contrast should be re-run without the bounded animals before it is believed.

| pair | n | Δ width | p | censored animals | n uncens. | Δ uncens. | p uncens. |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | -11 ms | 0.416 | 3 | 5 | -4 ms | 0.789 |
| CA1-CA3 | 6 | -3 ms | 0.175 | 1 | 5 | -4 ms | 0.178 |
| CA1-DG | 8 | -1 ms | 0.732 | 1 | 7 | -3 ms | 0.457 |
| CA1-V1 | 10 | +77 ms | 0.306 | 7 | 3 | +130 ms | 0.526 |
| CA3-DG | 4 | -2 ms | 0.718 | 1 | 3 | +3 ms | 0.423 |
| CA1-SUB | 4 | +18 ms | 0.133 | 2 | 2 | +15 ms | 0.205 |
| RSC-SUB | 4 | +65 ms | 0.559 | 2 | 2 | -45 ms | 0.614 |
| V1-RSC | 6 | +0 ms | 1 | 3 | 3 | -80 ms | 0.117 |

### FS-included — expert vs naive through a FROZEN subspace

Subspace identified once on trials balanced across epochs, then both epochs projected through the identical weights; animals-as-n paired *t*. `peak r` is in-sample by construction — read the Δ, not the level.

| pair | n | peak r naive | peak r exp | Δ | p | IFI naive | IFI exp | Δ | p | peak lag (ms) naive | peak lag (ms) exp | Δ | p | half-max width (ms) naive | half-max width (ms) exp | Δ | p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | 0.137 | 0.113 | -0.0245 | 0.121 | -0.00246 | -0.00178 | +0.000676 | 0.988 | 0 | -1.25 | -1.25 | 0.685 | 71.2 | 32.5 | -38.8 | 0.186 |
| CA1-CA3 | 7 | 0.338 | 0.367 | +0.0298 | 0.283 | 0.0414 | 0.0137 | -0.0277 | 0.569 | 1.43 | 0 | -1.43 | 0.356 | 25.7 | 22.9 | -2.86 | 0.457 |
| CA1-DG | 8 | 0.232 | 0.217 | -0.0149 | 0.388 | -0.0448 | -0.0232 | +0.0215 | 0.704 | 0 | 0 | +0 | nan | 20 | 17.5 | -2.5 | 0.626 |
| CA1-V1 | 10 | 0.174 | 0.152 | -0.0229 | 0.217 | -0.0886 | 0.036 | +0.125 | 0.0789 | -12 | -20 | -8 | 0.7 | 205 | 206 | +1 | 0.923 |
| CA3-DG | 4 | 0.421 | 0.469 | +0.0488 | 0.0701 | -0.0313 | 0.019 | +0.0503 | 0.492 | 0 | 0 | +0 | nan | 20 | 17.5 | -2.5 | 0.718 |
| CA1-SUB | 4 | 0.13 | 0.124 | -0.00626 | 0.584 | -0.161 | -0.0471 | +0.114 | 0.115 | 0 | 0 | +0 | nan | 0 | 12.5 | +12.5 | 0.194 |
| RSC-SUB | 4 | 0.22 | 0.225 | +0.00511 | 0.882 | 0.102 | 0.0631 | -0.0391 | 0.45 | 2.5 | 0 | -2.5 | 0.391 | 155 | 285 | +130 | 0.371 |
| V1-RSC | 6 | 0.32 | 0.284 | -0.0353 | 0.273 | 0.039 | -0.0169 | -0.0559 | 0.312 | 3.33 | 0 | -3.33 | 0.363 | 255 | 168 | -86.7 | 0.0645 |

#### FS-included — is the lag curve decaying or RINGING?

A secondary peak above half the central peak means the curve oscillates rather than decays. On a ringing curve the half-max width is the **half-period of the rhythm, not a temporal integration window**, and the IFI compares two lobes of that rhythm rather than a lead.

| pair | curves with a side peak | side-peak lag | implied freq | mean width | read the width as |
|---|---|---|---|---|---|
| CA1-CA3 | 13/21 (62%) | 135 ms | 7.4 Hz | 24 ms | **half-period of a rhythm** |
| CA1-DG | 13/24 (54%) | 131 ms | 7.6 Hz | 20 ms | **half-period of a rhythm** |
| CA1-RSC | 18/24 (75%) | 143 ms | 7.0 Hz | 63 ms | **half-period of a rhythm** |
| CA1-SUB | 5/12 (42%) | 82 ms | 12.2 Hz | 11 ms | mixed — mostly integration window |
| CA1-V1 | 9/30 (30%) | 154 ms | 6.5 Hz | 188 ms | mixed — mostly integration window |
| CA3-DG | 5/12 (42%) | 136 ms | 7.4 Hz | 18 ms | mixed — mostly integration window |
| RSC-SUB | 4/12 (33%) | 72 ms | 13.8 Hz | 204 ms | mixed — mostly integration window |
| V1-RSC | 9/18 (50%) | 108 ms | 9.3 Hz | 203 ms | **half-period of a rhythm** |

#### FS-included — how censored is the integration window?

`width_ms` is bounded at both ends. Of 153 (animal, pair, epoch) values, **20% sit at the 0 ms floor** (only the peak bin clears half-max) and **10% at the 500 ms ceiling** (the curve never drops below half-max); only 70% are interior. A paired *t* on this metric is therefore fragile, and any contrast should be re-run without the bounded animals before it is believed.

| pair | n | Δ width | p | censored animals | n uncens. | Δ uncens. | p uncens. |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | -39 ms | 0.186 | 2 | 6 | -33 ms | 0.387 |
| CA1-CA3 | 7 | -3 ms | 0.457 | 1 | 6 | -3 ms | 0.465 |
| CA1-DG | 8 | -2 ms | 0.626 | 1 | 7 | -1 ms | 0.805 |
| CA1-V1 | 10 | +1 ms | 0.923 | 8 | 2 | +35 ms | 0.579 |
| CA3-DG | 4 | -2 ms | 0.718 | 1 | 3 | +3 ms | 0.423 |
| CA1-SUB | 4 | +12 ms | 0.194 | 4 | 0 | +nan ms | nan |
| RSC-SUB | 4 | +130 ms | 0.371 | 2 | 2 | +5 ms | 0.874 |
| V1-RSC | 6 | -87 ms | 0.0645 | 2 | 4 | -130 ms | 0.0432 |
