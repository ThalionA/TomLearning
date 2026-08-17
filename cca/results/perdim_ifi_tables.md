# Per-dimension IFI — does directionality depend on canonical rank?

Meeting 2026-07-28 item 1. Re-analysis of the held-out segment-aware per-dim
lag curves (`lag_curves_bin10*.csv`); no refits. **Positive IFI = the
first-named area leads.** Animals-as-n throughout; 8 pairs, **no cross-pair
correction** (per-pair family, STATE.md §3.0 policy).

`sign-agree` = per-animal fraction of significant tail dims whose IFI shares
CC₁'s sign; 0.5 = chance. LMM = `ifi ~ dim + (1+dim|animal)` over significant
dims only.

> **Caveat.** CCA is refit at every lag, so `dim` is canonical *rank*, not a
> tracked component. These are statements about rank, not about a persistent
> subspace — see meeting item 3.

### FS-excluded

| pair | n | IFI CC₁ | IFI sig-tail | Δ (CC₁−tail) | t | p | sign-agree | p vs 0.5 | LMM slope/dim | LMM p |
|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 12 | +0.026 | +0.024 | +0.001 | 0.10 | 0.924 | 0.62 | 0.0857 | +0.0011 | 0.598 |
| CA1-CA3 | 7 | -0.013 | -0.029 | +0.015 | 0.77 | 0.469 | 0.60 | 0.425 | -0.0017 | 0.786 |
| CA1-DG | 11 | -0.005 | -0.028 | +0.023 | 1.66 | 0.129 | 0.49 | 0.857 | -0.0060 | 0.0219 |
| CA1-V1 | 13 | -0.018 | -0.034 | +0.016 | 1.24 | 0.238 | 0.51 | 0.86 | -0.0030 | 0.368 |
| CA3-DG | 5 | +0.021 | -0.025 | +0.046 | 3.03 | 0.0389 | 0.54 | 0.569 | -0.0073 | 0.0613 |
| CA1-SUB | 7 | +0.017 | +0.015 | +0.002 | 0.08 | 0.938 | 0.67 | 0.147 | +0.0017 | 0.615 |
| RSC-SUB | 7 | +0.015 | -0.008 | +0.023 | 2.08 | 0.0824 | 0.43 | 0.364 | -0.0028 | 0.36 |
| V1-RSC | 9 | +0.014 | +0.039 | -0.025 | -1.51 | 0.17 | 0.60 | 0.3 | +0.0019 | 0.508 |

### FS-included

| pair | n | IFI CC₁ | IFI sig-tail | Δ (CC₁−tail) | t | p | sign-agree | p vs 0.5 | LMM slope/dim | LMM p |
|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 12 | +0.059 | +0.032 | +0.027 | 1.57 | 0.146 | 0.71 | 0.000562 | +0.0015 | 0.535 |
| CA1-CA3 | 9 | +0.004 | -0.013 | +0.016 | 1.22 | 0.257 | 0.60 | 0.142 | -0.0032 | 0.296 |
| CA1-DG | 11 | -0.005 | -0.024 | +0.019 | 1.26 | 0.236 | 0.54 | 0.649 | -0.0026 | 0.267 |
| CA1-V1 | 13 | -0.026 | -0.026 | +0.001 | 0.09 | 0.93 | 0.59 | 0.236 | -0.0015 | 0.55 |
| CA3-DG | 5 | +0.011 | -0.024 | +0.035 | 3.44 | 0.0262 | 0.42 | 0.488 | -0.0085 | 0.0145 |
| CA1-SUB | 7 | +0.012 | +0.032 | -0.020 | -0.69 | 0.517 | 0.54 | 0.657 | +0.0040 | 0.257 |
| RSC-SUB | 7 | +0.003 | -0.004 | +0.007 | 0.47 | 0.655 | 0.46 | 0.564 | +0.0001 | 0.982 |
| V1-RSC | 9 | +0.015 | +0.034 | -0.020 | -1.40 | 0.199 | 0.79 | 0.00169 | -0.0002 | 0.937 |


## Is canonical rank meaningful out of sample?

**FS-excluded.** 656 significant dims. Held-out CC collapses with rank (rank 1 = 0.170, rank 5 = 0.050, rank 20 = 0.009; Spearman rank vs CC = -0.846), so rank *is* meaningful on average. But only 54% of significant dims sit at rank ≤ 5 and 0% sit beyond rank 20, where the mean held-out CC is at the floor; the significant set is a contiguous leading run in only 90% of 71 cells. **The significant tail is therefore dominated by floor-level dims** — gate on CC magnitude, not the `sig` flag alone, before assigning any per-dim label (e.g. FF/FB).

**FS-included.** 701 significant dims. Held-out CC collapses with rank (rank 1 = 0.209, rank 5 = 0.056, rank 20 = 0.009; Spearman rank vs CC = -0.865), so rank *is* meaningful on average. But only 52% of significant dims sit at rank ≤ 5 and 0% sit beyond rank 20, where the mean held-out CC is at the floor; the significant set is a contiguous leading run in only 95% of 73 cells. **The significant tail is therefore dominated by floor-level dims** — gate on CC magnitude, not the `sig` flag alone, before assigning any per-dim label (e.g. FF/FB).
