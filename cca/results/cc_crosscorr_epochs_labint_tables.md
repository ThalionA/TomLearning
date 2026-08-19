# 2026-08-07 ask 2 — cross-correlograms of naive and expert, whole and by FF/FB label

Frozen-axes lag curves (`cc_label_track_bin10*.csv`), averaged over each animal's significant CCs first, then across animals. `r` is in-sample for the whole session: a contrast across epochs, not a coupling strength.

### FS-excluded — cross-correlograms by epoch, all significant CCs (frozen axes, per-animal-first)

Label column for the FF/FB split: `label_int` — sign of the INTERMEDIATE epoch's own IFI at ±50 ms — 10 trials disjoint from both plotted epochs (Tom's suggestion; noisy, few trials). Paired *t* across animals on each animal's mean over its significant CCs of the per-CC reduction already in `cc_label_track_epoch_*` (IFI at ±50 ms; peak r). Positive IFI ⇒ first-named area leads. Per-pair families; BH across the 8 pairs is a sensitivity column.

> ⚠ **Curve-height (peak r) contrasts against naive are not a learning readout on this all-trials fit** — see the script docstring: naive is uniquely low, ~60 % of that is a lag-independent OFFSET (`far_r`, the |lag| ≥ 200 ms baseline, moves with it — slow co-modulation, speed the obvious candidate), the rise is LP-independent, and intermediate (pre-LP) already equals expert. `peak_minus_far` is the coupling-specific statistic; `intermediate-naive` sits next to `expert-naive`. IFI (shape) is the naive-vs-expert statistic and is null.

> Rows with n ≤ 4 animals (3 df) are descriptive and are not bolded.

> **Two units per row.** *animals-as-n* (inferential unit): each animal's significant CCs are averaged first, then a paired *t* across animals. *CCs-as-n* (power check, field convention): every significant (animal, dim) is one paired sample, pooled across animals — CCs are nested in animals, so this over-counts and is never the inferential statement.

| pair | metric | contrast | a → b (animals) | Δ | animals: n, t, p, BH | CCs: n, Δ, t, p, BH |
|---|---|---|---|---|---|---|
| CA1-RSC | ifi | expert-naive | +0.013 → +0.018 | +0.004 | 8, +0.09, 0.93, no | 74, +0.003, +0.07, 0.945, no |
| CA1-CA3 | ifi | expert-naive | +0.014 → +0.064 | +0.050 | 6, +1.13, 0.312, no | 56, +0.059, +0.81, 0.423, no |
| CA1-DG | ifi | expert-naive | +0.015 → +0.030 | +0.014 | 8, +0.25, 0.809, no | 75, +0.003, +0.06, 0.954, no |
| CA1-V1 | ifi | expert-naive | -0.016 → -0.003 | +0.014 | 10, +1.10, 0.298, no | 97, +0.017, +0.61, 0.546, no |
| CA3-DG | ifi | expert-naive | -0.145 → -0.142 | +0.003 | 4, +0.03, 0.98, no | 40, +0.003, +0.04, 0.97, no |
| CA1-SUB | ifi | expert-naive | -0.283 → -0.068 | +0.215 | 4, +1.66, 0.196, no | 35, +0.162, +2.52, **0.0167**, no |
| RSC-SUB | ifi | expert-naive | +0.018 → -0.051 | -0.068 | 4, -0.94, 0.418, no | 36, -0.044, -0.94, 0.355, no |
| V1-RSC | ifi | expert-naive | +0.057 → +0.070 | +0.013 | 6, +0.49, 0.647, no | 57, +0.014, +0.48, 0.634, no |
| CA1-RSC | peak_r | expert-naive | +0.030 → +0.048 | **+0.018** | 8, +3.76, **0.00706**, yes | 74, +0.018, +6.85, **2.01e-09**, yes |
| CA1-CA3 | peak_r | expert-naive | +0.058 → +0.067 | +0.009 | 6, +2.30, 0.0696, no | 56, +0.009, +2.04, **0.0462**, yes |
| CA1-DG | peak_r | expert-naive | +0.049 → +0.055 | +0.006 | 8, +2.22, 0.0619, no | 75, +0.006, +2.76, **0.00738**, yes |
| CA1-V1 | peak_r | expert-naive | +0.037 → +0.057 | **+0.020** | 10, +4.29, **0.00203**, yes | 97, +0.020, +5.96, **4.14e-08**, yes |
| CA3-DG | peak_r | expert-naive | +0.108 → +0.122 | +0.014 | 4, +5.25, 0.0135, yes | 40, +0.014, +3.86, **0.000415**, yes |
| CA1-SUB | peak_r | expert-naive | +0.030 → +0.040 | +0.010 | 4, +1.71, 0.185, no | 35, +0.012, +3.06, **0.00433**, yes |
| RSC-SUB | peak_r | expert-naive | +0.053 → +0.069 | +0.016 | 4, +2.69, 0.0743, no | 36, +0.017, +3.68, **0.000785**, yes |
| V1-RSC | peak_r | expert-naive | +0.079 → +0.108 | +0.029 | 6, +2.51, 0.0538, no | 57, +0.031, +5.82, **2.94e-07**, yes |
| CA1-RSC | peak_r | intermediate-naive | +0.030 → +0.047 | **+0.017** | 8, +3.37, **0.0119**, yes | 74, +0.017, +6.51, **8.2e-09**, yes |
| CA1-CA3 | peak_r | intermediate-naive | +0.058 → +0.066 | +0.009 | 6, +2.31, 0.0691, no | 56, +0.009, +2.04, **0.046**, yes |
| CA1-DG | peak_r | intermediate-naive | +0.049 → +0.056 | **+0.007** | 8, +2.83, **0.0255**, no | 75, +0.007, +3.38, **0.00114**, yes |
| CA1-V1 | peak_r | intermediate-naive | +0.037 → +0.053 | **+0.016** | 10, +4.35, **0.00186**, yes | 97, +0.016, +5.33, **6.56e-07**, yes |
| CA3-DG | peak_r | intermediate-naive | +0.108 → +0.124 | +0.016 | 4, +3.49, 0.0398, no | 40, +0.016, +4.18, **0.000159**, yes |
| CA1-SUB | peak_r | intermediate-naive | +0.030 → +0.043 | +0.013 | 4, +1.78, 0.173, no | 35, +0.016, +3.65, **0.000879**, yes |
| RSC-SUB | peak_r | intermediate-naive | +0.053 → +0.071 | +0.018 | 4, +2.51, 0.0873, no | 36, +0.018, +3.83, **0.000511**, yes |
| V1-RSC | peak_r | intermediate-naive | +0.079 → +0.106 | +0.027 | 6, +2.48, 0.0557, no | 57, +0.029, +5.67, **5.12e-07**, yes |
| CA1-RSC | far_r | expert-naive | +0.010 → +0.025 | **+0.014** | 8, +2.75, **0.0283**, no | 74, +0.015, +5.42, **7.28e-07**, yes |
| CA1-CA3 | far_r | expert-naive | -0.004 → +0.001 | +0.005 | 6, +1.22, 0.277, no | 56, +0.006, +2.36, **0.0219**, yes |
| CA1-DG | far_r | expert-naive | -0.001 → +0.004 | +0.005 | 8, +1.98, 0.0879, no | 75, +0.005, +3.35, **0.00126**, yes |
| CA1-V1 | far_r | expert-naive | +0.022 → +0.041 | **+0.019** | 10, +4.43, **0.00165**, yes | 97, +0.019, +5.70, **1.34e-07**, yes |
| CA3-DG | far_r | expert-naive | +0.003 → +0.007 | +0.004 | 4, +1.91, 0.152, no | 40, +0.004, +2.14, **0.0384**, yes |
| CA1-SUB | far_r | expert-naive | +0.005 → +0.012 | +0.007 | 4, +1.51, 0.228, no | 35, +0.008, +2.67, **0.0115**, yes |
| RSC-SUB | far_r | expert-naive | +0.025 → +0.036 | +0.011 | 4, +1.50, 0.23, no | 36, +0.012, +2.94, **0.00578**, yes |
| V1-RSC | far_r | expert-naive | +0.049 → +0.074 | +0.025 | 6, +2.17, 0.0821, no | 57, +0.027, +5.02, **5.54e-06**, yes |
| CA1-RSC | peak_minus_far | expert-naive | +0.014 → +0.016 | +0.002 | 8, +2.05, 0.0799, no | 74, +0.004, +3.54, **0.000704**, yes |
| CA1-CA3 | peak_minus_far | expert-naive | +0.055 → +0.060 | +0.005 | 6, +0.92, 0.399, no | 56, +0.004, +0.84, 0.404, no |
| CA1-DG | peak_minus_far | expert-naive | +0.046 → +0.047 | +0.000 | 8, +0.19, 0.856, no | 75, +0.001, +0.29, 0.775, no |
| CA1-V1 | peak_minus_far | expert-naive | +0.009 → +0.009 | +0.001 | 10, +0.72, 0.489, no | 97, +0.001, +2.14, **0.0347**, no |
| CA3-DG | peak_minus_far | expert-naive | +0.102 → +0.112 | +0.010 | 4, +5.69, 0.0108, no | 40, +0.010, +2.70, **0.0102**, yes |
| CA1-SUB | peak_minus_far | expert-naive | +0.021 → +0.024 | +0.003 | 4, +0.85, 0.457, no | 35, +0.004, +2.16, **0.0379**, no |
| RSC-SUB | peak_minus_far | expert-naive | +0.024 → +0.028 | +0.004 | 4, +1.38, 0.261, no | 36, +0.005, +3.05, **0.00433**, yes |
| V1-RSC | peak_minus_far | expert-naive | +0.025 → +0.028 | +0.002 | 6, +1.46, 0.204, no | 57, +0.004, +3.35, **0.00146**, yes |
| CA1-RSC | peak_minus_far | intermediate-naive | +0.014 → +0.017 | **+0.003** | 8, +2.48, **0.0419**, no | 74, +0.004, +4.03, **0.000134**, yes |
| CA1-CA3 | peak_minus_far | intermediate-naive | +0.055 → +0.062 | +0.007 | 6, +1.40, 0.221, no | 56, +0.006, +1.48, 0.144, no |
| CA1-DG | peak_minus_far | intermediate-naive | +0.046 → +0.050 | +0.004 | 8, +2.16, 0.0677, no | 75, +0.004, +1.77, 0.0816, no |
| CA1-V1 | peak_minus_far | intermediate-naive | +0.009 → +0.010 | +0.002 | 10, +1.05, 0.322, no | 97, +0.002, +2.06, **0.0419**, no |
| CA3-DG | peak_minus_far | intermediate-naive | +0.102 → +0.114 | +0.012 | 4, +3.25, 0.0473, no | 40, +0.012, +3.62, **0.00083**, yes |
| CA1-SUB | peak_minus_far | intermediate-naive | +0.021 → +0.027 | +0.006 | 4, +2.03, 0.135, no | 35, +0.006, +3.10, **0.00391**, yes |
| RSC-SUB | peak_minus_far | intermediate-naive | +0.024 → +0.027 | +0.003 | 4, +1.21, 0.314, no | 36, +0.004, +2.95, **0.0056**, yes |
| V1-RSC | peak_minus_far | intermediate-naive | +0.025 → +0.027 | +0.002 | 6, +0.97, 0.375, no | 57, +0.003, +2.45, **0.0175**, yes |

**Animals per (pair, epoch) — group = all / FF / FB:**

| pair | naive | intermediate | expert |
|---|---|---|---|
| CA1-RSC | 8/8/8 | 8/8/8 | 8/8/8 |
| CA1-CA3 | 6/6/6 | 6/6/6 | 6/6/6 |
| CA1-DG | 8/8/8 | 8/8/8 | 8/8/8 |
| CA1-V1 | 10/10/10 | 10/10/10 | 10/10/10 |
| CA3-DG | 4/4/4 | 4/4/4 | 4/4/4 |
| CA1-SUB | 4/4/4 | 4/4/4 | 4/4/4 |
| RSC-SUB | 4/4/4 | 4/4/4 | 4/4/4 |
| V1-RSC | 6/6/6 | 6/6/6 | 6/6/6 |

**Per-label expert − naive (label = `label_int`), animals-as-n:**

| pair | label | metric | n | naive → expert | Δ | t | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | FF | ifi | 8 | -0.003 → +0.029 | +0.032 | +0.55 | 0.598 |
| CA1-CA3 | FF | ifi | 6 | +0.002 → +0.168 | +0.165 | +1.16 | 0.297 |
| CA1-DG | FF | ifi | 8 | +0.198 → +0.205 | +0.007 | +0.12 | 0.908 |
| CA1-V1 | FF | ifi | 10 | -0.026 → -0.012 | +0.015 | +0.87 | 0.408 |
| CA3-DG | FF | ifi | 4 | -0.106 → +0.042 | +0.148 | +0.62 | 0.579 |
| CA1-SUB | FF | ifi | 4 | -0.152 → +0.066 | +0.217 | +1.33 | 0.276 |
| RSC-SUB | FF | ifi | 4 | +0.097 → -0.048 | -0.145 | -1.21 | 0.312 |
| V1-RSC | FF | ifi | 6 | +0.075 → +0.078 | +0.002 | +0.07 | 0.946 |
| CA1-RSC | FB | ifi | 8 | +0.032 → -0.025 | -0.056 | -1.03 | 0.337 |
| CA1-CA3 | FB | ifi | 6 | -0.079 → -0.062 | +0.017 | +0.16 | 0.876 |
| CA1-DG | FB | ifi | 8 | -0.171 → -0.095 | +0.075 | +0.73 | 0.488 |
| CA1-V1 | FB | ifi | 10 | +0.000 → +0.005 | +0.004 | +0.24 | 0.818 |
| CA3-DG | FB | ifi | 4 | -0.190 → -0.229 | -0.038 | -0.72 | 0.522 |
| CA1-SUB | FB | ifi | 4 | -0.383 → -0.181 | +0.202 | +1.50 | 0.231 |
| RSC-SUB | FB | ifi | 4 | -0.044 → -0.052 | -0.008 | -0.34 | 0.758 |
| V1-RSC | FB | ifi | 6 | -0.029 → +0.071 | +0.100 | +1.95 | 0.109 |
| CA1-RSC | FF | peak_r | 8 | +0.032 → +0.052 | **+0.020** | +3.79 | **0.00678** |
| CA1-CA3 | FF | peak_r | 6 | +0.064 → +0.079 | +0.015 | +1.98 | 0.104 |
| CA1-DG | FF | peak_r | 8 | +0.044 → +0.048 | +0.004 | +1.78 | 0.119 |
| CA1-V1 | FF | peak_r | 10 | +0.042 → +0.061 | **+0.019** | +2.89 | **0.0178** |
| CA3-DG | FF | peak_r | 4 | +0.115 → +0.129 | +0.014 | +2.00 | 0.139 |
| CA1-SUB | FF | peak_r | 4 | +0.030 → +0.032 | +0.002 | +0.33 | 0.766 |
| RSC-SUB | FF | peak_r | 4 | +0.051 → +0.053 | +0.002 | +0.48 | 0.665 |
| V1-RSC | FF | peak_r | 6 | +0.083 → +0.110 | +0.027 | +2.00 | 0.101 |
| CA1-RSC | FB | peak_r | 8 | +0.025 → +0.040 | **+0.015** | +3.15 | **0.0162** |
| CA1-CA3 | FB | peak_r | 6 | +0.069 → +0.070 | +0.001 | +0.35 | 0.74 |
| CA1-DG | FB | peak_r | 8 | +0.051 → +0.061 | **+0.010** | +2.91 | **0.0227** |
| CA1-V1 | FB | peak_r | 10 | +0.036 → +0.055 | **+0.020** | +4.11 | **0.00263** |
| CA3-DG | FB | peak_r | 4 | +0.098 → +0.112 | +0.014 | +22.50 | 0.000192 |
| CA1-SUB | FB | peak_r | 4 | +0.030 → +0.042 | +0.013 | +1.64 | 0.199 |
| RSC-SUB | FB | peak_r | 4 | +0.055 → +0.076 | +0.021 | +2.70 | 0.0737 |
| V1-RSC | FB | peak_r | 6 | +0.058 → +0.086 | **+0.028** | +4.14 | **0.00897** |

### FS-included — cross-correlograms by epoch, all significant CCs (frozen axes, per-animal-first)

Label column for the FF/FB split: `label_int` — sign of the INTERMEDIATE epoch's own IFI at ±50 ms — 10 trials disjoint from both plotted epochs (Tom's suggestion; noisy, few trials). Paired *t* across animals on each animal's mean over its significant CCs of the per-CC reduction already in `cc_label_track_epoch_*` (IFI at ±50 ms; peak r). Positive IFI ⇒ first-named area leads. Per-pair families; BH across the 8 pairs is a sensitivity column.

> ⚠ **Curve-height (peak r) contrasts against naive are not a learning readout on this all-trials fit** — see the script docstring: naive is uniquely low, ~60 % of that is a lag-independent OFFSET (`far_r`, the |lag| ≥ 200 ms baseline, moves with it — slow co-modulation, speed the obvious candidate), the rise is LP-independent, and intermediate (pre-LP) already equals expert. `peak_minus_far` is the coupling-specific statistic; `intermediate-naive` sits next to `expert-naive`. IFI (shape) is the naive-vs-expert statistic and is null.

> Rows with n ≤ 4 animals (3 df) are descriptive and are not bolded.

> **Two units per row.** *animals-as-n* (inferential unit): each animal's significant CCs are averaged first, then a paired *t* across animals. *CCs-as-n* (power check, field convention): every significant (animal, dim) is one paired sample, pooled across animals — CCs are nested in animals, so this over-counts and is never the inferential statement.

| pair | metric | contrast | a → b (animals) | Δ | animals: n, t, p, BH | CCs: n, Δ, t, p, BH |
|---|---|---|---|---|---|---|
| CA1-RSC | ifi | expert-naive | -0.012 → +0.088 | +0.101 | 8, +1.85, 0.107, no | 79, +0.099, +2.19, **0.0316**, no |
| CA1-CA3 | ifi | expert-naive | +0.024 → +0.029 | +0.005 | 7, +0.09, 0.929, no | 70, +0.005, +0.08, 0.933, no |
| CA1-DG | ifi | expert-naive | -0.003 → +0.030 | +0.033 | 8, +0.50, 0.63, no | 75, +0.037, +0.73, 0.465, no |
| CA1-V1 | ifi | expert-naive | -0.000 → -0.027 | -0.026 | 10, -1.44, 0.185, no | 98, -0.024, -0.84, 0.404, no |
| CA3-DG | ifi | expert-naive | -0.154 → -0.129 | +0.026 | 4, +0.55, 0.618, no | 40, +0.026, +0.39, 0.702, no |
| CA1-SUB | ifi | expert-naive | -0.256 → +0.014 | +0.270 | 4, +1.35, 0.27, no | 36, +0.213, +2.48, **0.0182**, no |
| RSC-SUB | ifi | expert-naive | -0.038 → -0.062 | -0.024 | 4, -1.17, 0.326, no | 36, -0.020, -0.62, 0.54, no |
| V1-RSC | ifi | expert-naive | +0.041 → +0.026 | -0.015 | 6, -0.73, 0.498, no | 60, -0.015, -1.14, 0.26, no |
| CA1-RSC | peak_r | expert-naive | +0.038 → +0.055 | **+0.017** | 8, +3.04, **0.0189**, no | 79, +0.017, +5.47, **5.19e-07**, yes |
| CA1-CA3 | peak_r | expert-naive | +0.083 → +0.095 | **+0.012** | 7, +2.80, **0.0313**, no | 70, +0.012, +2.93, **0.00465**, yes |
| CA1-DG | peak_r | expert-naive | +0.058 → +0.061 | +0.003 | 8, +0.91, 0.395, no | 75, +0.002, +0.92, 0.362, no |
| CA1-V1 | peak_r | expert-naive | +0.041 → +0.063 | **+0.022** | 10, +5.03, **0.000705**, yes | 98, +0.022, +5.99, **3.48e-08**, yes |
| CA3-DG | peak_r | expert-naive | +0.130 → +0.146 | +0.016 | 4, +7.39, 0.00512, yes | 40, +0.016, +5.62, **1.76e-06**, yes |
| CA1-SUB | peak_r | expert-naive | +0.033 → +0.045 | +0.012 | 4, +2.18, 0.118, no | 36, +0.013, +3.47, **0.0014**, yes |
| RSC-SUB | peak_r | expert-naive | +0.064 → +0.079 | +0.015 | 4, +2.44, 0.0929, no | 36, +0.016, +3.60, **0.000981**, yes |
| V1-RSC | peak_r | expert-naive | +0.095 → +0.123 | **+0.029** | 6, +2.84, **0.0363**, no | 60, +0.029, +6.08, **9.59e-08**, yes |
| CA1-RSC | peak_r | intermediate-naive | +0.038 → +0.055 | **+0.016** | 8, +2.66, **0.0322**, no | 79, +0.017, +5.50, **4.64e-07**, yes |
| CA1-CA3 | peak_r | intermediate-naive | +0.083 → +0.094 | **+0.011** | 7, +2.66, **0.0376**, no | 70, +0.011, +2.78, **0.00703**, yes |
| CA1-DG | peak_r | intermediate-naive | +0.058 → +0.063 | +0.005 | 8, +1.66, 0.141, no | 75, +0.004, +1.90, 0.0607, no |
| CA1-V1 | peak_r | intermediate-naive | +0.041 → +0.058 | **+0.017** | 10, +4.93, **0.000818**, yes | 98, +0.017, +5.48, **3.32e-07**, yes |
| CA3-DG | peak_r | intermediate-naive | +0.130 → +0.148 | +0.017 | 4, +4.25, 0.0239, no | 40, +0.017, +5.48, **2.71e-06**, yes |
| CA1-SUB | peak_r | intermediate-naive | +0.033 → +0.047 | +0.014 | 4, +2.20, 0.115, no | 36, +0.016, +3.96, **0.000353**, yes |
| RSC-SUB | peak_r | intermediate-naive | +0.064 → +0.080 | +0.016 | 4, +2.27, 0.108, no | 36, +0.016, +3.59, **0.00102**, yes |
| V1-RSC | peak_r | intermediate-naive | +0.095 → +0.120 | **+0.026** | 6, +2.91, **0.0334**, no | 60, +0.026, +5.61, **5.61e-07**, yes |
| CA1-RSC | far_r | expert-naive | +0.013 → +0.027 | **+0.014** | 8, +2.42, **0.0459**, no | 79, +0.014, +4.46, **2.76e-05**, yes |
| CA1-CA3 | far_r | expert-naive | -0.011 → -0.002 | **+0.009** | 7, +2.65, **0.0382**, no | 70, +0.009, +4.09, **0.000114**, yes |
| CA1-DG | far_r | expert-naive | -0.001 → +0.005 | **+0.006** | 8, +2.44, **0.045**, no | 75, +0.006, +3.73, **0.000371**, yes |
| CA1-V1 | far_r | expert-naive | +0.024 → +0.046 | **+0.021** | 10, +5.56, **0.000351**, yes | 98, +0.021, +5.81, **8.03e-08**, yes |
| CA3-DG | far_r | expert-naive | +0.001 → +0.003 | +0.001 | 4, +0.81, 0.476, no | 40, +0.001, +0.56, 0.58, no |
| CA1-SUB | far_r | expert-naive | +0.004 → +0.011 | +0.008 | 4, +1.87, 0.158, no | 36, +0.009, +2.95, **0.00567**, yes |
| RSC-SUB | far_r | expert-naive | +0.027 → +0.038 | +0.011 | 4, +1.36, 0.266, no | 36, +0.012, +2.64, **0.0123**, yes |
| V1-RSC | far_r | expert-naive | +0.060 → +0.085 | +0.025 | 6, +2.43, 0.0591, no | 60, +0.025, +5.38, **1.33e-06**, yes |
| CA1-RSC | peak_minus_far | expert-naive | +0.019 → +0.022 | +0.003 | 8, +1.70, 0.133, no | 79, +0.003, +2.06, **0.0426**, no |
| CA1-CA3 | peak_minus_far | expert-naive | +0.090 → +0.094 | +0.004 | 7, +0.98, 0.364, no | 70, +0.003, +0.71, 0.477, no |
| CA1-DG | peak_minus_far | expert-naive | +0.056 → +0.053 | -0.002 | 8, -0.97, 0.364, no | 75, -0.003, -1.22, 0.226, no |
| CA1-V1 | peak_minus_far | expert-naive | +0.011 → +0.011 | +0.001 | 10, +0.52, 0.613, no | 98, +0.001, +1.56, 0.121, no |
| CA3-DG | peak_minus_far | expert-naive | +0.127 → +0.141 | +0.014 | 4, +8.42, 0.00351, yes | 40, +0.014, +4.33, **9.97e-05**, yes |
| CA1-SUB | peak_minus_far | expert-naive | +0.026 → +0.029 | +0.003 | 4, +1.12, 0.345, no | 36, +0.004, +2.22, **0.0333**, no |
| RSC-SUB | peak_minus_far | expert-naive | +0.033 → +0.036 | +0.003 | 4, +1.13, 0.341, no | 36, +0.003, +2.13, **0.0407**, no |
| V1-RSC | peak_minus_far | expert-naive | +0.031 → +0.033 | +0.002 | 6, +2.34, 0.0663, no | 60, +0.003, +3.10, **0.00295**, yes |
| CA1-RSC | peak_minus_far | intermediate-naive | +0.019 → +0.023 | +0.003 | 8, +1.41, 0.2, no | 79, +0.003, +2.33, **0.0226**, no |
| CA1-CA3 | peak_minus_far | intermediate-naive | +0.090 → +0.096 | +0.006 | 7, +1.68, 0.144, no | 70, +0.005, +1.22, 0.226, no |
| CA1-DG | peak_minus_far | intermediate-naive | +0.056 → +0.057 | +0.001 | 8, +0.49, 0.636, no | 75, +0.000, +0.16, 0.87, no |
| CA1-V1 | peak_minus_far | intermediate-naive | +0.011 → +0.012 | +0.001 | 10, +0.82, 0.432, no | 98, +0.000, +0.49, 0.622, no |
| CA3-DG | peak_minus_far | intermediate-naive | +0.127 → +0.145 | +0.018 | 4, +8.24, 0.00375, yes | 40, +0.017, +4.53, **5.41e-05**, yes |
| CA1-SUB | peak_minus_far | intermediate-naive | +0.026 → +0.030 | +0.005 | 4, +1.82, 0.166, no | 36, +0.005, +2.43, **0.0205**, no |
| RSC-SUB | peak_minus_far | intermediate-naive | +0.033 → +0.034 | +0.000 | 4, +0.30, 0.787, no | 36, +0.002, +1.48, 0.148, no |
| V1-RSC | peak_minus_far | intermediate-naive | +0.031 → +0.032 | +0.001 | 6, +1.14, 0.306, no | 60, +0.002, +1.71, 0.0924, no |

**Animals per (pair, epoch) — group = all / FF / FB:**

| pair | naive | intermediate | expert |
|---|---|---|---|
| CA1-RSC | 8/8/7 | 8/8/7 | 8/8/7 |
| CA1-CA3 | 7/7/7 | 7/7/7 | 7/7/7 |
| CA1-DG | 8/8/8 | 8/8/8 | 8/8/8 |
| CA1-V1 | 10/10/10 | 10/10/10 | 10/10/10 |
| CA3-DG | 4/4/4 | 4/4/4 | 4/4/4 |
| CA1-SUB | 4/4/4 | 4/4/4 | 4/4/4 |
| RSC-SUB | 4/4/4 | 4/4/4 | 4/4/4 |
| V1-RSC | 6/6/6 | 6/6/6 | 6/6/6 |

**Per-label expert − naive (label = `label_int`), animals-as-n:**

| pair | label | metric | n | naive → expert | Δ | t | p |
|---|---|---|---|---|---|---|---|
| CA1-RSC | FF | ifi | 8 | +0.031 → +0.119 | +0.089 | +1.27 | 0.246 |
| CA1-CA3 | FF | ifi | 7 | +0.008 → +0.213 | +0.205 | +1.60 | 0.16 |
| CA1-DG | FF | ifi | 8 | +0.121 → +0.166 | +0.045 | +0.39 | 0.707 |
| CA1-V1 | FF | ifi | 10 | -0.034 → -0.012 | +0.022 | +0.74 | 0.476 |
| CA3-DG | FF | ifi | 4 | +0.012 → -0.019 | -0.031 | -0.67 | 0.549 |
| CA1-SUB | FF | ifi | 4 | -0.139 → +0.396 | +0.536 | +2.24 | 0.111 |
| RSC-SUB | FF | ifi | 4 | +0.055 → +0.005 | -0.050 | -1.96 | 0.145 |
| V1-RSC | FF | ifi | 6 | +0.059 → +0.040 | -0.019 | -0.74 | 0.495 |
| CA1-RSC | FB | ifi | 7 | -0.112 → +0.036 | +0.148 | +2.22 | 0.0679 |
| CA1-CA3 | FB | ifi | 7 | -0.014 → -0.166 | -0.152 | -1.13 | 0.302 |
| CA1-DG | FB | ifi | 8 | -0.075 → -0.070 | +0.005 | +0.07 | 0.949 |
| CA1-V1 | FB | ifi | 10 | +0.029 → -0.032 | -0.062 | -1.52 | 0.164 |
| CA3-DG | FB | ifi | 4 | -0.260 → -0.163 | +0.097 | +0.92 | 0.426 |
| CA1-SUB | FB | ifi | 4 | -0.278 → -0.227 | +0.051 | +0.52 | 0.639 |
| RSC-SUB | FB | ifi | 4 | -0.073 → -0.090 | -0.017 | -0.47 | 0.672 |
| V1-RSC | FB | ifi | 6 | -0.033 → -0.024 | +0.009 | +0.52 | 0.624 |
| CA1-RSC | FF | peak_r | 8 | +0.041 → +0.057 | **+0.016** | +2.78 | **0.0273** |
| CA1-CA3 | FF | peak_r | 7 | +0.090 → +0.107 | +0.017 | +1.61 | 0.159 |
| CA1-DG | FF | peak_r | 8 | +0.059 → +0.058 | -0.001 | -0.14 | 0.895 |
| CA1-V1 | FF | peak_r | 10 | +0.043 → +0.066 | **+0.023** | +3.00 | **0.0149** |
| CA3-DG | FF | peak_r | 4 | +0.158 → +0.175 | +0.017 | +4.82 | 0.017 |
| CA1-SUB | FF | peak_r | 4 | +0.028 → +0.029 | +0.001 | +0.55 | 0.62 |
| RSC-SUB | FF | peak_r | 4 | +0.053 → +0.053 | +0.000 | +0.01 | 0.992 |
| V1-RSC | FF | peak_r | 6 | +0.101 → +0.131 | **+0.030** | +3.13 | **0.0261** |
| CA1-RSC | FB | peak_r | 7 | +0.035 → +0.048 | +0.013 | +1.86 | 0.113 |
| CA1-CA3 | FB | peak_r | 7 | +0.080 → +0.086 | +0.006 | +1.34 | 0.23 |
| CA1-DG | FB | peak_r | 8 | +0.071 → +0.073 | +0.003 | +0.57 | 0.587 |
| CA1-V1 | FB | peak_r | 10 | +0.039 → +0.059 | **+0.020** | +6.21 | **0.000157** |
| CA3-DG | FB | peak_r | 4 | +0.109 → +0.125 | +0.016 | +6.28 | 0.00814 |
| CA1-SUB | FB | peak_r | 4 | +0.033 → +0.048 | +0.015 | +1.89 | 0.155 |
| RSC-SUB | FB | peak_r | 4 | +0.068 → +0.083 | +0.015 | +1.53 | 0.222 |
| V1-RSC | FB | peak_r | 6 | +0.055 → +0.083 | +0.027 | +2.13 | 0.0866 |
