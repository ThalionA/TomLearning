# 2026-08-07 ask 2 — cross-correlograms of naive and expert, whole and by FF/FB label

Frozen-axes lag curves (`cc_label_track_bin10*.csv`), averaged over each animal's significant CCs first, then across animals. `r` is in-sample for the whole session: a contrast across epochs, not a coupling strength.

### FS-excluded — cross-correlograms by epoch, all significant CCs (frozen axes, per-animal-first)

Label column for the FF/FB split: `label`. Paired *t* across animals on each animal's mean over its significant CCs of the per-CC reduction already in `cc_label_track_epoch_*` (IFI at ±50 ms; peak r). Positive IFI ⇒ first-named area leads. Per-pair families; BH across the 8 pairs is a sensitivity column.

> ⚠ **Curve-height (peak r) contrasts against naive are not a learning readout on this all-trials fit** — see the script docstring: naive is uniquely low, ~60 % of that is a lag-independent OFFSET (`far_r`, the |lag| ≥ 200 ms baseline, moves with it — slow co-modulation, speed the obvious candidate), the rise is LP-independent, and intermediate (pre-LP) already equals expert. `peak_minus_far` is the coupling-specific statistic; `intermediate-naive` sits next to `expert-naive`. IFI (shape) is the naive-vs-expert statistic and is null.

> Rows with n ≤ 4 animals (3 df) are descriptive and are not bolded.

| pair | n | metric | contrast | a | b | Δ (b − a) | t | p | BH |
|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | ifi | expert-naive | +0.013 | +0.018 | +0.004 | +0.09 | 0.93 | no |
| CA1-CA3 | 6 | ifi | expert-naive | +0.014 | +0.064 | +0.050 | +1.13 | 0.312 | no |
| CA1-DG | 8 | ifi | expert-naive | +0.015 | +0.030 | +0.014 | +0.25 | 0.809 | no |
| CA1-V1 | 10 | ifi | expert-naive | -0.016 | -0.003 | +0.014 | +1.10 | 0.298 | no |
| CA3-DG | 4 | ifi | expert-naive | -0.145 | -0.142 | +0.003 | +0.03 | 0.98 | no |
| CA1-SUB | 4 | ifi | expert-naive | -0.283 | -0.068 | +0.215 | +1.66 | 0.196 | no |
| RSC-SUB | 4 | ifi | expert-naive | +0.018 | -0.051 | -0.068 | -0.94 | 0.418 | no |
| V1-RSC | 6 | ifi | expert-naive | +0.057 | +0.070 | +0.013 | +0.49 | 0.647 | no |
| CA1-RSC | 8 | peak_r | expert-naive | +0.030 | +0.048 | **+0.018** | +3.76 | **0.00706** | yes |
| CA1-CA3 | 6 | peak_r | expert-naive | +0.058 | +0.067 | +0.009 | +2.30 | 0.0696 | no |
| CA1-DG | 8 | peak_r | expert-naive | +0.049 | +0.055 | +0.006 | +2.22 | 0.0619 | no |
| CA1-V1 | 10 | peak_r | expert-naive | +0.037 | +0.057 | **+0.020** | +4.29 | **0.00203** | yes |
| CA3-DG | 4 | peak_r | expert-naive | +0.108 | +0.122 | +0.014 | +5.25 | 0.0135 | yes |
| CA1-SUB | 4 | peak_r | expert-naive | +0.030 | +0.040 | +0.010 | +1.71 | 0.185 | no |
| RSC-SUB | 4 | peak_r | expert-naive | +0.053 | +0.069 | +0.016 | +2.69 | 0.0743 | no |
| V1-RSC | 6 | peak_r | expert-naive | +0.079 | +0.108 | +0.029 | +2.51 | 0.0538 | no |
| CA1-RSC | 8 | peak_r | intermediate-naive | +0.030 | +0.047 | **+0.017** | +3.37 | **0.0119** | yes |
| CA1-CA3 | 6 | peak_r | intermediate-naive | +0.058 | +0.066 | +0.009 | +2.31 | 0.0691 | no |
| CA1-DG | 8 | peak_r | intermediate-naive | +0.049 | +0.056 | **+0.007** | +2.83 | **0.0255** | no |
| CA1-V1 | 10 | peak_r | intermediate-naive | +0.037 | +0.053 | **+0.016** | +4.35 | **0.00186** | yes |
| CA3-DG | 4 | peak_r | intermediate-naive | +0.108 | +0.124 | +0.016 | +3.49 | 0.0398 | no |
| CA1-SUB | 4 | peak_r | intermediate-naive | +0.030 | +0.043 | +0.013 | +1.78 | 0.173 | no |
| RSC-SUB | 4 | peak_r | intermediate-naive | +0.053 | +0.071 | +0.018 | +2.51 | 0.0873 | no |
| V1-RSC | 6 | peak_r | intermediate-naive | +0.079 | +0.106 | +0.027 | +2.48 | 0.0557 | no |
| CA1-RSC | 8 | far_r | expert-naive | +0.010 | +0.025 | **+0.014** | +2.75 | **0.0283** | no |
| CA1-CA3 | 6 | far_r | expert-naive | -0.004 | +0.001 | +0.005 | +1.22 | 0.277 | no |
| CA1-DG | 8 | far_r | expert-naive | -0.001 | +0.004 | +0.005 | +1.98 | 0.0879 | no |
| CA1-V1 | 10 | far_r | expert-naive | +0.022 | +0.041 | **+0.019** | +4.43 | **0.00165** | yes |
| CA3-DG | 4 | far_r | expert-naive | +0.003 | +0.007 | +0.004 | +1.91 | 0.152 | no |
| CA1-SUB | 4 | far_r | expert-naive | +0.005 | +0.012 | +0.007 | +1.51 | 0.228 | no |
| RSC-SUB | 4 | far_r | expert-naive | +0.025 | +0.036 | +0.011 | +1.50 | 0.23 | no |
| V1-RSC | 6 | far_r | expert-naive | +0.049 | +0.074 | +0.025 | +2.17 | 0.0821 | no |
| CA1-RSC | 8 | peak_minus_far | expert-naive | +0.014 | +0.016 | +0.002 | +2.05 | 0.0799 | no |
| CA1-CA3 | 6 | peak_minus_far | expert-naive | +0.055 | +0.060 | +0.005 | +0.92 | 0.399 | no |
| CA1-DG | 8 | peak_minus_far | expert-naive | +0.046 | +0.047 | +0.000 | +0.19 | 0.856 | no |
| CA1-V1 | 10 | peak_minus_far | expert-naive | +0.009 | +0.009 | +0.001 | +0.72 | 0.489 | no |
| CA3-DG | 4 | peak_minus_far | expert-naive | +0.102 | +0.112 | +0.010 | +5.69 | 0.0108 | no |
| CA1-SUB | 4 | peak_minus_far | expert-naive | +0.021 | +0.024 | +0.003 | +0.85 | 0.457 | no |
| RSC-SUB | 4 | peak_minus_far | expert-naive | +0.024 | +0.028 | +0.004 | +1.38 | 0.261 | no |
| V1-RSC | 6 | peak_minus_far | expert-naive | +0.025 | +0.028 | +0.002 | +1.46 | 0.204 | no |
| CA1-RSC | 8 | peak_minus_far | intermediate-naive | +0.014 | +0.017 | **+0.003** | +2.48 | **0.0419** | no |
| CA1-CA3 | 6 | peak_minus_far | intermediate-naive | +0.055 | +0.062 | +0.007 | +1.40 | 0.221 | no |
| CA1-DG | 8 | peak_minus_far | intermediate-naive | +0.046 | +0.050 | +0.004 | +2.16 | 0.0677 | no |
| CA1-V1 | 10 | peak_minus_far | intermediate-naive | +0.009 | +0.010 | +0.002 | +1.05 | 0.322 | no |
| CA3-DG | 4 | peak_minus_far | intermediate-naive | +0.102 | +0.114 | +0.012 | +3.25 | 0.0473 | no |
| CA1-SUB | 4 | peak_minus_far | intermediate-naive | +0.021 | +0.027 | +0.006 | +2.03 | 0.135 | no |
| RSC-SUB | 4 | peak_minus_far | intermediate-naive | +0.024 | +0.027 | +0.003 | +1.21 | 0.314 | no |
| V1-RSC | 6 | peak_minus_far | intermediate-naive | +0.025 | +0.027 | +0.002 | +0.97 | 0.375 | no |

**Animals per (pair, epoch) — group = all / FF / FB:**

| pair | naive | intermediate | expert |
|---|---|---|---|
| CA1-RSC | 8/8/7 | 8/8/7 | 8/8/7 |
| CA1-CA3 | 6/6/6 | 6/6/6 | 6/6/6 |
| CA1-DG | 8/8/8 | 8/8/8 | 8/8/8 |
| CA1-V1 | 10/10/10 | 10/10/10 | 10/10/10 |
| CA3-DG | 4/4/4 | 4/4/4 | 4/4/4 |
| CA1-SUB | 4/2/4 | 4/2/4 | 4/2/4 |
| RSC-SUB | 4/4/4 | 4/4/4 | 4/4/4 |
| V1-RSC | 6/6/6 | 6/6/6 | 6/6/6 |

### FS-included — cross-correlograms by epoch, all significant CCs (frozen axes, per-animal-first)

Label column for the FF/FB split: `label`. Paired *t* across animals on each animal's mean over its significant CCs of the per-CC reduction already in `cc_label_track_epoch_*` (IFI at ±50 ms; peak r). Positive IFI ⇒ first-named area leads. Per-pair families; BH across the 8 pairs is a sensitivity column.

> ⚠ **Curve-height (peak r) contrasts against naive are not a learning readout on this all-trials fit** — see the script docstring: naive is uniquely low, ~60 % of that is a lag-independent OFFSET (`far_r`, the |lag| ≥ 200 ms baseline, moves with it — slow co-modulation, speed the obvious candidate), the rise is LP-independent, and intermediate (pre-LP) already equals expert. `peak_minus_far` is the coupling-specific statistic; `intermediate-naive` sits next to `expert-naive`. IFI (shape) is the naive-vs-expert statistic and is null.

> Rows with n ≤ 4 animals (3 df) are descriptive and are not bolded.

| pair | n | metric | contrast | a | b | Δ (b − a) | t | p | BH |
|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | ifi | expert-naive | -0.012 | +0.088 | +0.101 | +1.85 | 0.107 | no |
| CA1-CA3 | 7 | ifi | expert-naive | +0.024 | +0.029 | +0.005 | +0.09 | 0.929 | no |
| CA1-DG | 8 | ifi | expert-naive | -0.003 | +0.030 | +0.033 | +0.50 | 0.63 | no |
| CA1-V1 | 10 | ifi | expert-naive | -0.000 | -0.027 | -0.026 | -1.44 | 0.185 | no |
| CA3-DG | 4 | ifi | expert-naive | -0.154 | -0.129 | +0.026 | +0.55 | 0.618 | no |
| CA1-SUB | 4 | ifi | expert-naive | -0.256 | +0.014 | +0.270 | +1.35 | 0.27 | no |
| RSC-SUB | 4 | ifi | expert-naive | -0.038 | -0.062 | -0.024 | -1.17 | 0.326 | no |
| V1-RSC | 6 | ifi | expert-naive | +0.041 | +0.026 | -0.015 | -0.73 | 0.498 | no |
| CA1-RSC | 8 | peak_r | expert-naive | +0.038 | +0.055 | **+0.017** | +3.04 | **0.0189** | no |
| CA1-CA3 | 7 | peak_r | expert-naive | +0.083 | +0.095 | **+0.012** | +2.80 | **0.0313** | no |
| CA1-DG | 8 | peak_r | expert-naive | +0.058 | +0.061 | +0.003 | +0.91 | 0.395 | no |
| CA1-V1 | 10 | peak_r | expert-naive | +0.041 | +0.063 | **+0.022** | +5.03 | **0.000705** | yes |
| CA3-DG | 4 | peak_r | expert-naive | +0.130 | +0.146 | +0.016 | +7.39 | 0.00512 | yes |
| CA1-SUB | 4 | peak_r | expert-naive | +0.033 | +0.045 | +0.012 | +2.18 | 0.118 | no |
| RSC-SUB | 4 | peak_r | expert-naive | +0.064 | +0.079 | +0.015 | +2.44 | 0.0929 | no |
| V1-RSC | 6 | peak_r | expert-naive | +0.095 | +0.123 | **+0.029** | +2.84 | **0.0363** | no |
| CA1-RSC | 8 | peak_r | intermediate-naive | +0.038 | +0.055 | **+0.016** | +2.66 | **0.0322** | no |
| CA1-CA3 | 7 | peak_r | intermediate-naive | +0.083 | +0.094 | **+0.011** | +2.66 | **0.0376** | no |
| CA1-DG | 8 | peak_r | intermediate-naive | +0.058 | +0.063 | +0.005 | +1.66 | 0.141 | no |
| CA1-V1 | 10 | peak_r | intermediate-naive | +0.041 | +0.058 | **+0.017** | +4.93 | **0.000818** | yes |
| CA3-DG | 4 | peak_r | intermediate-naive | +0.130 | +0.148 | +0.017 | +4.25 | 0.0239 | no |
| CA1-SUB | 4 | peak_r | intermediate-naive | +0.033 | +0.047 | +0.014 | +2.20 | 0.115 | no |
| RSC-SUB | 4 | peak_r | intermediate-naive | +0.064 | +0.080 | +0.016 | +2.27 | 0.108 | no |
| V1-RSC | 6 | peak_r | intermediate-naive | +0.095 | +0.120 | **+0.026** | +2.91 | **0.0334** | no |
| CA1-RSC | 8 | far_r | expert-naive | +0.013 | +0.027 | **+0.014** | +2.42 | **0.0459** | no |
| CA1-CA3 | 7 | far_r | expert-naive | -0.011 | -0.002 | **+0.009** | +2.65 | **0.0382** | no |
| CA1-DG | 8 | far_r | expert-naive | -0.001 | +0.005 | **+0.006** | +2.44 | **0.045** | no |
| CA1-V1 | 10 | far_r | expert-naive | +0.024 | +0.046 | **+0.021** | +5.56 | **0.000351** | yes |
| CA3-DG | 4 | far_r | expert-naive | +0.001 | +0.003 | +0.001 | +0.81 | 0.476 | no |
| CA1-SUB | 4 | far_r | expert-naive | +0.004 | +0.011 | +0.008 | +1.87 | 0.158 | no |
| RSC-SUB | 4 | far_r | expert-naive | +0.027 | +0.038 | +0.011 | +1.36 | 0.266 | no |
| V1-RSC | 6 | far_r | expert-naive | +0.060 | +0.085 | +0.025 | +2.43 | 0.0591 | no |
| CA1-RSC | 8 | peak_minus_far | expert-naive | +0.019 | +0.022 | +0.003 | +1.70 | 0.133 | no |
| CA1-CA3 | 7 | peak_minus_far | expert-naive | +0.090 | +0.094 | +0.004 | +0.98 | 0.364 | no |
| CA1-DG | 8 | peak_minus_far | expert-naive | +0.056 | +0.053 | -0.002 | -0.97 | 0.364 | no |
| CA1-V1 | 10 | peak_minus_far | expert-naive | +0.011 | +0.011 | +0.001 | +0.52 | 0.613 | no |
| CA3-DG | 4 | peak_minus_far | expert-naive | +0.127 | +0.141 | +0.014 | +8.42 | 0.00351 | yes |
| CA1-SUB | 4 | peak_minus_far | expert-naive | +0.026 | +0.029 | +0.003 | +1.12 | 0.345 | no |
| RSC-SUB | 4 | peak_minus_far | expert-naive | +0.033 | +0.036 | +0.003 | +1.13 | 0.341 | no |
| V1-RSC | 6 | peak_minus_far | expert-naive | +0.031 | +0.033 | +0.002 | +2.34 | 0.0663 | no |
| CA1-RSC | 8 | peak_minus_far | intermediate-naive | +0.019 | +0.023 | +0.003 | +1.41 | 0.2 | no |
| CA1-CA3 | 7 | peak_minus_far | intermediate-naive | +0.090 | +0.096 | +0.006 | +1.68 | 0.144 | no |
| CA1-DG | 8 | peak_minus_far | intermediate-naive | +0.056 | +0.057 | +0.001 | +0.49 | 0.636 | no |
| CA1-V1 | 10 | peak_minus_far | intermediate-naive | +0.011 | +0.012 | +0.001 | +0.82 | 0.432 | no |
| CA3-DG | 4 | peak_minus_far | intermediate-naive | +0.127 | +0.145 | +0.018 | +8.24 | 0.00375 | yes |
| CA1-SUB | 4 | peak_minus_far | intermediate-naive | +0.026 | +0.030 | +0.005 | +1.82 | 0.166 | no |
| RSC-SUB | 4 | peak_minus_far | intermediate-naive | +0.033 | +0.034 | +0.000 | +0.30 | 0.787 | no |
| V1-RSC | 6 | peak_minus_far | intermediate-naive | +0.031 | +0.032 | +0.001 | +1.14 | 0.306 | no |

**Animals per (pair, epoch) — group = all / FF / FB:**

| pair | naive | intermediate | expert |
|---|---|---|---|
| CA1-RSC | 8/8/7 | 8/8/7 | 8/8/7 |
| CA1-CA3 | 7/7/7 | 7/7/7 | 7/7/7 |
| CA1-DG | 8/8/8 | 8/8/8 | 8/8/8 |
| CA1-V1 | 10/10/10 | 10/10/10 | 10/10/10 |
| CA3-DG | 4/4/4 | 4/4/4 | 4/4/4 |
| CA1-SUB | 4/3/4 | 4/3/4 | 4/3/4 |
| RSC-SUB | 4/3/4 | 4/3/4 | 4/3/4 |
| V1-RSC | 6/6/6 | 6/6/6 | 6/6/6 |
