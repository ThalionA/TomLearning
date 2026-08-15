# 2026-08-07 ask 2 — cross-correlograms of naive and expert, whole and by FF/FB label

Frozen-axes lag curves (`cc_label_track_bin10*.csv`), averaged over each animal's significant CCs first, then across animals. `r` is in-sample for the whole session: a contrast across epochs, not a coupling strength.

### FS-excluded — cross-correlograms by epoch, all significant CCs (frozen axes, per-animal-first)

Label column for the FF/FB split: `label`. Paired *t* across animals on each animal's mean over its significant CCs of the per-CC reduction already in `cc_label_track_epoch_*` (IFI at ±50 ms; peak r). Positive IFI ⇒ first-named area leads. Per-pair families; BH across the 8 pairs is a sensitivity column.

> ⚠ **Read peak-r (level) contrasts against naive with the caveat in the script docstring:** on this all-trials fit the naive epoch is uniquely low and intermediate (pre-LP) already equals expert — `intermediate-naive` is printed next to `expert-naive` so that is visible. The balanced-trial frozen fit shows no peak-r rise. IFI (a shape ratio) is the trustworthy naive-vs-expert statistic here.

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
| CA3-DG | 4 | peak_r | expert-naive | +0.108 | +0.122 | **+0.014** | +5.25 | **0.0135** | yes |
| CA1-SUB | 4 | peak_r | expert-naive | +0.030 | +0.040 | +0.010 | +1.71 | 0.185 | no |
| RSC-SUB | 4 | peak_r | expert-naive | +0.053 | +0.069 | +0.016 | +2.69 | 0.0743 | no |
| V1-RSC | 6 | peak_r | expert-naive | +0.079 | +0.108 | +0.029 | +2.51 | 0.0538 | no |
| CA1-RSC | 8 | peak_r | intermediate-naive | +0.030 | +0.047 | **+0.017** | +3.37 | **0.0119** | yes |
| CA1-CA3 | 6 | peak_r | intermediate-naive | +0.058 | +0.066 | +0.009 | +2.31 | 0.0691 | no |
| CA1-DG | 8 | peak_r | intermediate-naive | +0.049 | +0.056 | **+0.007** | +2.83 | **0.0255** | no |
| CA1-V1 | 10 | peak_r | intermediate-naive | +0.037 | +0.053 | **+0.016** | +4.35 | **0.00186** | yes |
| CA3-DG | 4 | peak_r | intermediate-naive | +0.108 | +0.124 | **+0.016** | +3.49 | **0.0398** | no |
| CA1-SUB | 4 | peak_r | intermediate-naive | +0.030 | +0.043 | +0.013 | +1.78 | 0.173 | no |
| RSC-SUB | 4 | peak_r | intermediate-naive | +0.053 | +0.071 | +0.018 | +2.51 | 0.0873 | no |
| V1-RSC | 6 | peak_r | intermediate-naive | +0.079 | +0.106 | +0.027 | +2.48 | 0.0557 | no |

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

> ⚠ **Read peak-r (level) contrasts against naive with the caveat in the script docstring:** on this all-trials fit the naive epoch is uniquely low and intermediate (pre-LP) already equals expert — `intermediate-naive` is printed next to `expert-naive` so that is visible. The balanced-trial frozen fit shows no peak-r rise. IFI (a shape ratio) is the trustworthy naive-vs-expert statistic here.

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
| CA3-DG | 4 | peak_r | expert-naive | +0.130 | +0.146 | **+0.016** | +7.39 | **0.00512** | yes |
| CA1-SUB | 4 | peak_r | expert-naive | +0.033 | +0.045 | +0.012 | +2.18 | 0.118 | no |
| RSC-SUB | 4 | peak_r | expert-naive | +0.064 | +0.079 | +0.015 | +2.44 | 0.0929 | no |
| V1-RSC | 6 | peak_r | expert-naive | +0.095 | +0.123 | **+0.029** | +2.84 | **0.0363** | no |
| CA1-RSC | 8 | peak_r | intermediate-naive | +0.038 | +0.055 | **+0.016** | +2.66 | **0.0322** | no |
| CA1-CA3 | 7 | peak_r | intermediate-naive | +0.083 | +0.094 | **+0.011** | +2.66 | **0.0376** | no |
| CA1-DG | 8 | peak_r | intermediate-naive | +0.058 | +0.063 | +0.005 | +1.66 | 0.141 | no |
| CA1-V1 | 10 | peak_r | intermediate-naive | +0.041 | +0.058 | **+0.017** | +4.93 | **0.000818** | yes |
| CA3-DG | 4 | peak_r | intermediate-naive | +0.130 | +0.148 | **+0.017** | +4.25 | **0.0239** | no |
| CA1-SUB | 4 | peak_r | intermediate-naive | +0.033 | +0.047 | +0.014 | +2.20 | 0.115 | no |
| RSC-SUB | 4 | peak_r | intermediate-naive | +0.064 | +0.080 | +0.016 | +2.27 | 0.108 | no |
| V1-RSC | 6 | peak_r | intermediate-naive | +0.095 | +0.120 | **+0.026** | +2.91 | **0.0334** | no |

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
