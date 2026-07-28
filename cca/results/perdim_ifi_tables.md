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
| CA1-RSC | 4 | +0.160 | -0.033 | +0.193 | 1.22 | 0.309 | 0.49 | 0.941 | -0.0109 | 0.0843 |
| CA1-CA3 | 4 | -0.084 | -0.135 | +0.050 | 0.59 | 0.594 | 0.79 | 0.112 | -0.0111 | 0.427 |
| CA1-DG | 9 | -0.062 | -0.070 | +0.008 | 0.13 | 0.898 | 0.68 | 0.114 | +nan | nan |
| CA1-V1 | 8 | -0.186 | -0.085 | -0.101 | -1.84 | 0.108 | 0.57 | 0.643 | +0.0001 | 0.991 |
| CA3-DG | 5 | +0.026 | +0.041 | -0.016 | -0.79 | 0.475 | 0.76 | 0.0815 | +0.0020 | 0.882 |
| CA1-SUB | 4 | +0.013 | +0.051 | -0.037 | -0.58 | 0.603 | 0.63 | 0.0922 | +0.0061 | 0.257 |
| RSC-SUB | 1 | -0.028 | +0.066 | -0.095 | nan | nan | 0.42 | nan | +0.0100 | 0.848 |
| V1-RSC | 8 | +0.097 | +0.063 | +0.033 | 0.65 | 0.537 | 0.52 | 0.877 | -0.0019 | 0.769 |

### FS-included

| pair | n | IFI CC₁ | IFI sig-tail | Δ (CC₁−tail) | t | p | sign-agree | p vs 0.5 | LMM slope/dim | LMM p |
|---|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 5 | +0.158 | +0.018 | +0.141 | 1.80 | 0.146 | 0.81 | 0.0626 | -0.0080 | 0.378 |
| CA1-CA3 | 9 | -0.064 | -0.041 | -0.024 | -0.55 | 0.594 | 0.61 | 0.429 | -0.0042 | 0.371 |
| CA1-DG | 8 | -0.072 | -0.083 | +0.010 | 0.19 | 0.855 | 0.64 | 0.25 | +0.0112 | 0.601 |
| CA1-V1 | 8 | -0.114 | -0.104 | -0.010 | -0.15 | 0.885 | 0.67 | 0.209 | +0.0012 | 0.845 |
| CA3-DG | 5 | +0.021 | -0.018 | +0.039 | 1.86 | 0.136 | 0.47 | 0.864 | -0.0036 | 0.656 |
| CA1-SUB | 4 | -0.105 | +0.080 | -0.185 | -1.70 | 0.187 | 0.25 | 0.215 | +nan | nan |
| RSC-SUB | 3 | -0.034 | +0.021 | -0.056 | -4.93 | 0.0388 | 0.86 | 0.118 | +nan | nan |
| V1-RSC | 8 | +0.079 | +0.104 | -0.024 | -0.67 | 0.525 | 0.76 | 0.117 | -0.0037 | 0.692 |


## Is canonical rank meaningful out of sample?

**FS-excluded.** 243 significant dims. Held-out CC collapses with rank (rank 1 = 0.170, rank 5 = 0.045, rank 20 = 0.024; Spearman rank vs CC = -0.657), so rank *is* meaningful on average. But only 44% of significant dims sit at rank ≤ 5 and 9% sit beyond rank 20, where the mean held-out CC is at the floor; the significant set is a contiguous leading run in only 46% of 56 cells. **The significant tail is therefore dominated by floor-level dims** — gate on CC magnitude, not the `sig` flag alone, before assigning any per-dim label (e.g. FF/FB).

**FS-included.** 282 significant dims. Held-out CC collapses with rank (rank 1 = 0.215, rank 5 = 0.047, rank 20 = 0.026; Spearman rank vs CC = -0.670), so rank *is* meaningful on average. But only 45% of significant dims sit at rank ≤ 5 and 15% sit beyond rank 20, where the mean held-out CC is at the floor; the significant set is a contiguous leading run in only 51% of 63 cells. **The significant tail is therefore dominated by floor-level dims** — gate on CC magnitude, not the `sig` flag alone, before assigning any per-dim label (e.g. FF/FB).
