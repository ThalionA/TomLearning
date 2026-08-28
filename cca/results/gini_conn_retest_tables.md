# Connection-specific contribution Gini — re-test of the participation-broadening headline

Generated 2026-08-28 by `scripts/analyze_trajectory.py trajectory_gini3_bin10{,_fsincl}.csv`.
Same CSV, same rows, same windows for all five metrics — a like-for-like re-test, not a re-run.

- `gini_x` / `gini_y` — the shipped **area-intrinsic** (provably partner-invariant) Gini.
- `gini_*_conn` — contributions weighted by each dim's held-out canonical correlation.
- `gini_*_sig` — significant dims only. **NaN when `n_sig` = 0**, so its slopes are conditioned
  on windows that had a significant dim (and its animal count can be lower — see `n`).

Cells = LMM population-slope p (random slope over animals, all windows). `n` = animals for
`gini_x` / for `gini_*_sig`. Bold = p < 0.05, **uncorrected, per-pair family, no cross-pair
correction** (project policy, STATE.md §3.0). All slopes are negative (de-sparsifying) unless
the sign is given.

⚠ **n = 4 pairs (CA1-SUB, CA3-DG, RSC-SUB) cannot be read at the honest unit:** the two-sided
Wilcoxon floor at n = 4 is p = 0.125, so their animals-as-n test can never star, and the small
LMM p-values there come from pooling windows within 4 animals.

## axis = `trial_frac`

| pair | FS | n | `gini_x` | `gini_x_conn` | `gini_x_sig` | `gini_y_conn` | `gini_y_sig` |
|---|---|---|---|---|---|---|---|
| CA1-RSC | FS-excl | 8/7 | **1.42e-05** | 0.0806 | **0.0188** | 0.709 | **0.00265** |
| CA1-RSC | FS-incl | 8 | **1.16e-05** | 0.223 | **0.0375** | 0.905 | **0.0284** |
| CA1-DG | FS-excl | 8 | **0.00879** | 0.272 | 0.202 | 0.943 | 0.294 |
| CA1-DG | FS-incl | 8 | **0.0281** | 0.203 | 0.063 | 0.903 | 0.764 |
| CA1-CA3 | FS-excl | 6 | 0.159 | 0.102 | 0.0718 | 0.141 | 0.0986 |
| CA1-CA3 | FS-incl | 7 | 0.141 | 0.0784 | 0.0833 | 0.0749 | 0.0702 |
| CA1-V1 | FS-excl | 10/9 | 0.181 | 0.416 | 0.0993 | **0.0269** | **0.00605** |
| CA1-V1 | FS-incl | 10/9 | 0.351 | 0.335 | **0.0023** | **0.00561** | 0.0969 |
| V1-RSC | FS-excl | 6 | 0.0738 | 0.147 | 0.0607 | 0.483 | 0.0866 |
| V1-RSC | FS-incl | 6 | 0.119 | 0.383 | 0.259 | 0.803 | 0.989 |
| CA1-SUB | FS-excl | 4 | **6.02e-06** | **0.000756** | **0.0306** | 0.968 | 0.606 |
| CA1-SUB | FS-incl | 4 | **3.62e-06** | **0.00346** | **0.0257** | 0.584 | 0.614 |
| CA3-DG | FS-excl | 4 | 0.252 | 0.221 | 0.363 | 0.262 | 0.309 |
| CA3-DG | FS-incl | 4 | 0.0632 | 0.115 | **0.0286** | 0.681 | 0.917 |
| RSC-SUB | FS-excl | 4 | 0.697 | 0.704 | 0.369 | 0.45 | 0.917 |
| RSC-SUB | FS-incl | 4 | 0.696 | 0.641 | 0.11 | 0.499 | 0.996 |

## axis = `performance`

| pair | FS | n | `gini_x` | `gini_x_conn` | `gini_x_sig` | `gini_y_conn` | `gini_y_sig` |
|---|---|---|---|---|---|---|---|
| CA1-RSC | FS-excl | 8/7 | **0.00186** | **0.000634** | **0.00194** | **0.0224** | **0.00135** |
| CA1-RSC | FS-incl | 8 | **0.000612** | **0.000102** | **0.000194** | **0.00126** | **0.00027** |
| CA1-DG | FS-excl | 8 | **0.00389** | 0.293 | 0.0853 | 0.67 | 0.413 |
| CA1-DG | FS-incl | 8 | **0.024** | 0.2 | 0.0515 | 0.612 | 0.394 |
| CA1-CA3 | FS-excl | 6 | 0.136 | 0.13 | 0.191 | 0.0835 | 0.057 |
| CA1-CA3 | FS-incl | 7 | **0.0289** | 0.168 | 0.287 | 0.0538 | 0.0725 |
| CA1-V1 | FS-excl | 10/9 | 0.174 | 0.294 | 0.299 | **0.0293** | 0.0744 |
| CA1-V1 | FS-incl | 10/9 | 0.407 | 0.265 | 0.0733 | **0.0152** | 0.205 |
| V1-RSC | FS-excl | 6 | 0.651 | **0.0399** | 0.873 | **0.0159** | **0.0441** |
| V1-RSC | FS-incl | 6 | 0.717 | 0.589 | 0.633 | **0.00631** | 0.0593 |
| CA1-SUB | FS-excl | 4 | **0.00194** | **0.0249** | 0.296 | 0.465 | 0.356 |
| CA1-SUB | FS-incl | 4 | **0.00198** | **0.00836** | 0.0807 | 0.789 | 0.383 |
| CA3-DG | FS-excl | 4 | **2.86e-05** | 0.151 | 0.124 | 0.714 | 0.72 |
| CA3-DG | FS-incl | 4 | **0.0303** | 0.08 | 0.142 | 0.789 | 0.991 |
| RSC-SUB | FS-excl | 4 | **0.0457** | 0.327 | nan | 0.312 | 0.697 |
| RSC-SUB | FS-incl | 4 | 0.0612 | **0.000838** | 0.0588 | 0.531 | 0.542 |

## axis = `lp_rel`

| pair | FS | n | `gini_x` | `gini_x_conn` | `gini_x_sig` | `gini_y_conn` | `gini_y_sig` |
|---|---|---|---|---|---|---|---|
| CA1-RSC | FS-excl | 8/7 | 0.972 | 0.444 | 0.191 | 0.819 | **0.0129** |
| CA1-RSC | FS-incl | 8 | 0.595 | 0.481 | **0.0238** | 0.974 | 0.289 |
| CA1-DG | FS-excl | 8 | **0.000256** | 0.262 | 0.172 | 0.889 | 0.34 |
| CA1-DG | FS-incl | 8 | 0.0827 | 0.704 | 0.796 | 0.874 | 0.858 |
| CA1-CA3 | FS-excl | 6 | 0.954 | 0.961 | 0.973 | 0.208 | 0.158 |
| CA1-CA3 | FS-incl | 7 | 0.954 | 0.0751 | 0.973 | 0.959 | 0.158 |
| CA1-V1 | FS-excl | 10/9 | 0.498 | 0.995 | 0.66 | **0.0108** | 0.407 |
| CA1-V1 | FS-incl | 10/9 | 0.696 | 0.797 | **0.00291** | **0.00181** | 0.0802 |
| V1-RSC | FS-excl | 6 | **0.0172** | **0.0189** | **0.0147** | 0.981 | 0.737 |
| V1-RSC | FS-incl | 6 | 0.417 | 0.737 | 0.925 | 0.471 | nan |
| CA1-SUB | FS-excl | 4 | 0.796 | 0.281 | **0.0262** | 0.995 | 0.626 |
| CA1-SUB | FS-incl | 4 | 0.355 | 0.318 | 0.108 | 0.998 | 0.799 |
| CA3-DG | FS-excl | 4 | 0.679 | 0.981 | 0.713 | 0.359 | 0.245 |
| CA3-DG | FS-incl | 4 | **0.00979** | 0.0528 | 0.924 | 0.797 | 0.875 |
| RSC-SUB | FS-excl | 4 | 0.963 | 0.849 | 0.22 | 0.612 | 0.958 |
| RSC-SUB | FS-incl | 4 | 0.982 | 0.842 | 0.909 | 0.994 | 0.989 |

