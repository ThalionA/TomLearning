# Fixed-subspace lag curves by epoch — meeting items 5, 6, 7

One canonical component, identified once across balanced trials and frozen, then lagged across time within each learning epoch. Because the weights are identical across epochs, an epoch difference is a difference in the **activity**, not in the fit.

> **Not a coupling strength.** The frozen fit saw every epoch, so `peak r` is > in-sample and optimistic. It is a contrast statistic only; the leak-free > numbers live in `lag_subspaces_tables.md`.

### FS-excluded — expert vs naive through a FROZEN subspace

Subspace identified once on trials balanced across epochs, then both epochs projected through the identical weights; animals-as-n paired *t*. `peak r` is in-sample by construction — read the Δ, not the level.

| pair | n | peak r naive | peak r exp | Δ | p | IFI naive | IFI exp | Δ | p | peak lag (ms) naive | peak lag (ms) exp | Δ | p | half-max width (ms) naive | half-max width (ms) exp | Δ | p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | 0.109 | 0.083 | -0.0263 | 0.0938 | 0.0774 | 0.0176 | -0.0597 | 0.432 | 35 | 23.8 | -11.2 | 0.793 | 124 | 45 | -78.8 | 0.255 |
| CA1-CA3 | 6 | 0.204 | 0.238 | +0.0345 | 0.392 | -0.17 | -0.123 | +0.0469 | 0.619 | -33.3 | -8.33 | +25 | 0.363 | 30 | 30 | +0 | 1 |
| CA1-DG | 8 | 0.174 | 0.18 | +0.00686 | 0.559 | -0.0979 | 0.0189 | +0.117 | 0.131 | 20 | -10 | -30 | 0.395 | 21.2 | 18.8 | -2.5 | 0.351 |
| CA1-V1 | 10 | 0.116 | 0.104 | -0.0116 | 0.525 | -0.0699 | 0.0553 | +0.125 | 0.0543 | 14 | -32 | -46 | 0.34 | 181 | 203 | +22 | 0.749 |
| CA3-DG | 4 | 0.209 | 0.218 | +0.00908 | 0.64 | -0.0806 | -0.0185 | +0.0621 | 0.416 | -10 | 45 | +55 | 0.484 | 25 | 22.5 | -2.5 | 0.638 |
| CA1-SUB | 4 | 0.0924 | 0.109 | +0.0171 | 0.303 | -0.0787 | -0.0658 | +0.0129 | 0.845 | -35 | 2.5 | +37.5 | 0.427 | 15 | 37.5 | +22.5 | 0.117 |
| RSC-SUB | 4 | 0.219 | 0.202 | -0.0164 | 0.673 | 0.0902 | 0.0917 | +0.00147 | 0.972 | 2.5 | 57.5 | +55 | 0.416 | 192 | 288 | +95 | 0.535 |
| V1-RSC | 6 | 0.244 | 0.226 | -0.0183 | 0.611 | 0.0552 | -0.0105 | -0.0658 | 0.249 | 33.3 | -1.67 | -35 | 0.115 | 252 | 252 | +0 | 1 |

#### FS-excluded — is the lag curve decaying or RINGING?

A secondary peak above half the central peak means the curve oscillates rather than decays. On a ringing curve the half-max width is the **half-period of the rhythm, not a temporal integration window**, and the IFI compares two lobes of that rhythm rather than a lead.

| pair | curves with a side peak | side-peak lag | implied freq | mean width | read the width as |
|---|---|---|---|---|---|
| CA1-CA3 | 14/18 (78%) | 134 ms | 7.5 Hz | 29 ms | **half-period of a rhythm** |
| CA1-DG | 18/24 (75%) | 170 ms | 5.9 Hz | 22 ms | **half-period of a rhythm** |
| CA1-RSC | 22/24 (92%) | 153 ms | 6.5 Hz | 79 ms | **half-period of a rhythm** |
| CA1-SUB | 11/12 (92%) | 137 ms | 7.3 Hz | 53 ms | **half-period of a rhythm** |
| CA1-V1 | 20/30 (67%) | 176 ms | 5.7 Hz | 183 ms | **half-period of a rhythm** |
| CA3-DG | 9/12 (75%) | 171 ms | 5.8 Hz | 26 ms | **half-period of a rhythm** |
| RSC-SUB | 5/12 (42%) | 132 ms | 7.6 Hz | 220 ms | mixed — mostly integration window |
| V1-RSC | 6/18 (33%) | 140 ms | 7.1 Hz | 256 ms | mixed — mostly integration window |

### FS-included — expert vs naive through a FROZEN subspace

Subspace identified once on trials balanced across epochs, then both epochs projected through the identical weights; animals-as-n paired *t*. `peak r` is in-sample by construction — read the Δ, not the level.

| pair | n | peak r naive | peak r exp | Δ | p | IFI naive | IFI exp | Δ | p | peak lag (ms) naive | peak lag (ms) exp | Δ | p | half-max width (ms) naive | half-max width (ms) exp | Δ | p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | 0.13 | 0.1 | -0.0298 | 0.0755 | 0.00701 | 0.00504 | -0.00197 | 0.959 | -16.2 | -1.25 | +15 | 0.435 | 260 | 55 | -205 | 0.0294 |
| CA1-CA3 | 7 | 0.307 | 0.336 | +0.0286 | 0.334 | 0.0193 | 0.0143 | -0.005 | 0.933 | -31.4 | -12.9 | +18.6 | 0.395 | 32.9 | 25.7 | -7.14 | 0.14 |
| CA1-DG | 8 | 0.202 | 0.188 | -0.0138 | 0.432 | 0.0332 | 0.0857 | +0.0525 | 0.327 | 16.2 | 18.8 | +2.5 | 0.351 | 17.5 | 18.8 | +1.25 | 0.763 |
| CA1-V1 | 10 | 0.136 | 0.114 | -0.0218 | 0.201 | -0.0457 | 0.0571 | +0.103 | 0.191 | -8 | -2 | +6 | 0.849 | 218 | 251 | +33 | 0.508 |
| CA3-DG | 4 | 0.251 | 0.266 | +0.0154 | 0.249 | 0.047 | 0.141 | +0.0941 | 0.295 | 57.5 | 12.5 | -45 | 0.391 | 25 | 20 | -5 | 0.391 |
| CA1-SUB | 4 | 0.115 | 0.117 | +0.00188 | 0.749 | -0.128 | -0.0308 | +0.0972 | 0.164 | -10 | 0 | +10 | 0.391 | 0 | 15 | +15 | 0.103 |
| RSC-SUB | 4 | 0.224 | 0.23 | +0.00565 | 0.888 | 0.0877 | 0.0663 | -0.0214 | 0.57 | 2.5 | 0 | -2.5 | 0.391 | 192 | 285 | +92.5 | 0.515 |
| V1-RSC | 6 | 0.325 | 0.272 | -0.0529 | 0.0894 | 0.036 | -0.0137 | -0.0497 | 0.341 | 5 | 0 | -5 | 0.203 | 267 | 202 | -65 | 0.11 |

#### FS-included — is the lag curve decaying or RINGING?

A secondary peak above half the central peak means the curve oscillates rather than decays. On a ringing curve the half-max width is the **half-period of the rhythm, not a temporal integration window**, and the IFI compares two lobes of that rhythm rather than a lead.

| pair | curves with a side peak | side-peak lag | implied freq | mean width | read the width as |
|---|---|---|---|---|---|
| CA1-CA3 | 16/21 (76%) | 140 ms | 7.1 Hz | 28 ms | **half-period of a rhythm** |
| CA1-DG | 16/24 (67%) | 133 ms | 7.5 Hz | 21 ms | **half-period of a rhythm** |
| CA1-RSC | 19/24 (79%) | 176 ms | 5.7 Hz | 138 ms | **half-period of a rhythm** |
| CA1-SUB | 7/12 (58%) | 113 ms | 8.9 Hz | 12 ms | **half-period of a rhythm** |
| CA1-V1 | 16/30 (53%) | 158 ms | 6.3 Hz | 213 ms | **half-period of a rhythm** |
| CA3-DG | 9/12 (75%) | 186 ms | 5.4 Hz | 23 ms | **half-period of a rhythm** |
| RSC-SUB | 4/12 (33%) | 132 ms | 7.5 Hz | 228 ms | mixed — mostly integration window |
| V1-RSC | 8/18 (44%) | 130 ms | 7.7 Hz | 236 ms | mixed — mostly integration window |
