# Item 2 — factorial: effect of epoch, effect of FF/FB, and their interaction

Replaces the ΔFF − ΔFB contrast, which discarded the intermediate epoch and could not separate main effects from the interaction.

### FS-excluded — 2-way RM-ANOVA per pair: epoch × FF/FB label

Animals-as-n (an animal's components in a cell are averaged first); animals missing any of the 6 cells are dropped. **The label main effect on signed IFI is circular by construction** — components are labelled by the sign of their whole-session IFI — and is marked ⚠.

**peak r**

| pair | n | F epoch | p | F label | p | F epoch×label | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 7 | 12.2 | 0.00128 | 3.11 | 0.128 | 0.822 | 0.463 |
| CA1-CA3 | 6 | 9.71 | 0.00454 | 4.36e-06 | 0.998 | 0.124 | 0.885 |
| CA1-DG | 8 | 5.55 | 0.0168 | 0.789 | 0.404 | 0.0287 | 0.972 |
| CA1-V1 | 10 | 14.4 | 0.000187 | 0.663 | 0.436 | 1.74 | 0.204 |
| CA3-DG | 4 | 3.82 | 0.0852 | 0.0998 | 0.773 | 0.0975 | 0.909 |
| CA1-SUB | 2 | — | — | — | — | — | — |
| RSC-SUB | 4 | 6.32 | 0.0333 | 0.869 | 0.42 | 2.37 | 0.174 |
| V1-RSC | 6 | 7.17 | 0.0117 | 1.62 | 0.26 | 0.204 | 0.819 |

**|IFI|**

| pair | n | F epoch | p | F label | p | F epoch×label | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 7 | 1.8 | 0.207 | 0.53 | 0.494 | 0.528 | 0.603 |
| CA1-CA3 | 6 | 6.55 | 0.0152 | 2.52 | 0.174 | 0.0382 | 0.963 |
| CA1-DG | 8 | 1.69 | 0.22 | 1.19 | 0.311 | 1.32 | 0.298 |
| CA1-V1 | 10 | 4.27 | 0.0303 | 1.05 | 0.332 | 1.17 | 0.333 |
| CA3-DG | 4 | 1.89 | 0.232 | 0.347 | 0.597 | 3.77 | 0.087 |
| CA1-SUB | 2 | — | — | — | — | — | — |
| RSC-SUB | 4 | 0.293 | 0.756 | 1.52 | 0.306 | 2.28 | 0.183 |
| V1-RSC | 6 | 0.713 | 0.513 | 3.15 | 0.136 | 0.08 | 0.924 |

**peak lag (ms)**

| pair | n | F epoch | p | F label | p | F epoch×label | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 7 | 0.628 | 0.55 | 0.887 | 0.383 | 0.454 | 0.646 |
| CA1-CA3 | 6 | 1.26 | 0.325 | 2.81 | 0.154 | 0.486 | 0.629 |
| CA1-DG | 8 | 0.177 | 0.839 | 0.285 | 0.61 | 0.937 | 0.415 |
| CA1-V1 | 10 | 1.73 | 0.205 | 0.379 | 0.553 | 0.852 | 0.443 |
| CA3-DG | 4 | 0.857 | 0.47 | 0.897 | 0.413 | 0.158 | 0.857 |
| CA1-SUB | 2 | — | — | — | — | — | — |
| RSC-SUB | 4 | 0.631 | 0.564 | 0.0376 | 0.859 | 1.96 | 0.222 |
| V1-RSC | 6 | 3.73 | 0.0617 | 1.25 | 0.314 | 4.7 | 0.0363 |

**IFI (signed)**

| pair | n | F epoch | p | F label | p | F epoch×label | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 7 | 0.0175 | 0.983 | 4.13 ⚠ | 0.0883 ⚠ | 3.52 | 0.0626 |
| CA1-CA3 | 6 | 0.769 | 0.489 | 10.3 ⚠ | 0.0238 ⚠ | 1.29 | 0.317 |
| CA1-DG | 8 | 0.614 | 0.555 | 15.5 ⚠ | 0.0056 ⚠ | 3.29 | 0.0675 |
| CA1-V1 | 10 | 0.212 | 0.811 | 3.99 ⚠ | 0.0769 ⚠ | 1.57 | 0.236 |
| CA3-DG | 4 | 0.599 | 0.579 | 6.56 ⚠ | 0.0831 ⚠ | 1.89 | 0.23 |
| CA1-SUB | 2 | — | — | — ⚠ | — ⚠ | — | — |
| RSC-SUB | 4 | 1.03 | 0.413 | 63.1 ⚠ | 0.00416 ⚠ | 2.59 | 0.155 |
| V1-RSC | 6 | 1.49 | 0.271 | 5.73 ⚠ | 0.0621 ⚠ | 1.13 | 0.362 |


#### Pooled across pairs — LMM `value ~ epoch * label + pair`, random intercept per animal

| DV | term | χ² | p | animals | obs |
|---|---|---|---|---|---|
| peak_r | epoch | 8.69 | 0.013 | 12 | 291 |
| peak_r | label | 0.0362 | 0.849 | 12 | 291 |
| peak_r | epoch:label | 0.0403 | 0.98 | 12 | 291 |
| abs_ifi | epoch | 13.6 | 0.00112 | 12 | 291 |
| abs_ifi | label | 2.16 | 0.141 | 12 | 291 |
| abs_ifi | epoch:label | 0.965 | 0.617 | 12 | 291 |
| peak_lag_ms | epoch | 0.833 | 0.659 | 12 | 291 |
| peak_lag_ms | label | 4.14 | 0.042 | 12 | 291 |
| peak_lag_ms | epoch:label | 2.3 | 0.316 | 12 | 291 |
| ifi | epoch | 0.535 | 0.765 | 12 | 291 |
| ifi | label ⚠ | 26.3 | 2.91e-07 | 12 | 291 |
| ifi | epoch:label | 0.485 | 0.785 | 12 | 291 |

### FS-included — 2-way RM-ANOVA per pair: epoch × FF/FB label

Animals-as-n (an animal's components in a cell are averaged first); animals missing any of the 6 cells are dropped. **The label main effect on signed IFI is circular by construction** — components are labelled by the sign of their whole-session IFI — and is marked ⚠.

**peak r**

| pair | n | F epoch | p | F label | p | F epoch×label | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 7 | 6.37 | 0.013 | 0.181 | 0.685 | 0.766 | 0.486 |
| CA1-CA3 | 7 | 9.82 | 0.00297 | 0.0343 | 0.859 | 0.57 | 0.58 |
| CA1-DG | 8 | 2.17 | 0.151 | 0.0106 | 0.921 | 0.0323 | 0.968 |
| CA1-V1 | 10 | 33 | 9.43e-07 | 0.217 | 0.653 | 0.58 | 0.57 |
| CA3-DG | 4 | 9.97 | 0.0124 | 0.863 | 0.421 | 1.11 | 0.389 |
| CA1-SUB | 3 | 2.28 | 0.218 | 5.25 | 0.149 | 1.51 | 0.325 |
| RSC-SUB | 3 | 1.38 | 0.351 | 0.179 | 0.713 | 1.22 | 0.386 |
| V1-RSC | 6 | 8.47 | 0.00704 | 3.52 | 0.12 | 0.93 | 0.426 |

**|IFI|**

| pair | n | F epoch | p | F label | p | F epoch×label | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 7 | 3.18 | 0.078 | 0.971 | 0.362 | 2.08 | 0.168 |
| CA1-CA3 | 7 | 8.34 | 0.00536 | 2.09 | 0.198 | 0.247 | 0.785 |
| CA1-DG | 8 | 1.68 | 0.222 | 0.0492 | 0.831 | 0.662 | 0.531 |
| CA1-V1 | 10 | 5.22 | 0.0163 | 0.00636 | 0.938 | 1.62 | 0.226 |
| CA3-DG | 4 | 3.87 | 0.0833 | 13.6 | 0.0345 | 2.94 | 0.129 |
| CA1-SUB | 3 | 0.416 | 0.685 | 2.1 | 0.284 | 1.58 | 0.313 |
| RSC-SUB | 3 | 2.28 | 0.219 | 0.526 | 0.544 | 0.851 | 0.492 |
| V1-RSC | 6 | 0.00425 | 0.996 | 0.443 | 0.535 | 0.434 | 0.66 |

**peak lag (ms)**

| pair | n | F epoch | p | F label | p | F epoch×label | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 7 | 2.08 | 0.167 | 4.64 | 0.0748 | 0.898 | 0.433 |
| CA1-CA3 | 7 | 0.122 | 0.886 | 3.72 | 0.102 | 0.682 | 0.524 |
| CA1-DG | 8 | 0.253 | 0.78 | 0.00165 | 0.969 | 0.701 | 0.513 |
| CA1-V1 | 10 | 1.07 | 0.364 | 0.505 | 0.496 | 1.19 | 0.327 |
| CA3-DG | 4 | 2.52 | 0.16 | 0.0221 | 0.891 | 2.65 | 0.15 |
| CA1-SUB | 3 | 0.0603 | 0.942 | 3.22 | 0.214 | 3.65 | 0.125 |
| RSC-SUB | 3 | 0.251 | 0.789 | 0.0356 | 0.868 | 0.906 | 0.474 |
| V1-RSC | 6 | 5.1 | 0.0298 | 0.128 | 0.735 | 0.732 | 0.505 |

**IFI (signed)**

| pair | n | F epoch | p | F label | p | F epoch×label | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | 7 | 6.08 | 0.015 | 15.3 ⚠ | 0.00792 ⚠ | 0.895 | 0.434 |
| CA1-CA3 | 7 | 0.222 | 0.804 | 16.9 ⚠ | 0.00631 ⚠ | 0.858 | 0.449 |
| CA1-DG | 8 | 1.06 | 0.373 | 2.12 ⚠ | 0.189 ⚠ | 0.149 | 0.863 |
| CA1-V1 | 10 | 1.5 | 0.25 | 6.07 ⚠ | 0.0359 ⚠ | 0.32 | 0.73 |
| CA3-DG | 4 | 1.19 | 0.366 | 54.3 ⚠ | 0.00517 ⚠ | 0.143 | 0.87 |
| CA1-SUB | 3 | 13.4 | 0.0169 | 8.26 ⚠ | 0.103 ⚠ | 5.57 | 0.0698 |
| RSC-SUB | 3 | 1.9 | 0.263 | 5.39 ⚠ | 0.146 ⚠ | 0.392 | 0.699 |
| V1-RSC | 6 | 0.294 | 0.752 | 15.7 ⚠ | 0.0107 ⚠ | 0.501 | 0.62 |


#### Pooled across pairs — LMM `value ~ epoch * label + pair`, random intercept per animal

| DV | term | χ² | p | animals | obs |
|---|---|---|---|---|---|
| peak_r | epoch | 5.75 | 0.0565 | 12 | 297 |
| peak_r | label | 0.28 | 0.597 | 12 | 297 |
| peak_r | epoch:label | 0.0546 | 0.973 | 12 | 297 |
| abs_ifi | epoch | 13.1 | 0.0014 | 12 | 297 |
| abs_ifi | label | 0.892 | 0.345 | 12 | 297 |
| abs_ifi | epoch:label | 1.24 | 0.539 | 12 | 297 |
| peak_lag_ms | epoch | 2.35 | 0.309 | 12 | 297 |
| peak_lag_ms | label | 4.37 | 0.0365 | 12 | 297 |
| peak_lag_ms | epoch:label | 2.56 | 0.278 | 12 | 297 |
| ifi | epoch | 2.42 | 0.298 | 12 | 297 |
| ifi | label ⚠ | 33.1 | 8.81e-09 | 12 | 297 |
| ifi | epoch:label | 1.16 | 0.561 | 12 | 297 |
