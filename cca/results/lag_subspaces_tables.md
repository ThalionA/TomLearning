# Lagged communication subspaces — meeting items 2, 3, 4

Animals-as-n throughout; 8 pairs, no cross-pair correction (per-pair family, STATE.md §3.0 policy). Angles are the LARGEST principal angle over 3 canonical dims, so a subspace that matches on its dominant axis but diverges elsewhere is not scored as stable.

### FS-excluded — item 3: subspace stability across lag (CC₁ only)

`angle − floor` is the principal angle between the lag-0 and lagged subspace minus that pair's own split-half floor, averaged over the two areas; animals-as-n, Bonferroni across |lag| within a pair.

| pair | estimable? | stability width | mean angle @ ±50 ms | floor | Δ | p (Bonf) |
|---|---|---|---|---|---|---|
| CA1-CA3 | yes | ≥ 250 ms (censored) | 54.6° | 48.7° | +5.8° | 1 |
| CA1-DG | yes | ≥ 250 ms (censored) | 59.8° | 55.6° | +4.3° | 1 |
| CA1-RSC | yes | ≥ 250 ms (censored) | 66.5° | 63.6° | +2.8° | 1 |
| CA1-SUB | yes | ≥ 250 ms (censored) | 72.2° | 62.0° | +10.2° | 1 |
| CA1-V1 | yes | ≥ 250 ms (censored) | 55.4° | 64.6° | -9.2° | 1 |
| CA3-DG | yes | ≥ 250 ms (censored) | 47.5° | 50.3° | -2.7° | 1 |
| RSC-SUB | yes | ≥ 250 ms (censored) | 59.3° | 63.6° | -4.3° | 1 |
| V1-RSC | yes | ≥ 250 ms (censored) | 43.3° | 42.6° | +0.7° | 1 |


### FS-excluded — item 3: subspace stability across lag (3 canonical dims)

`angle − floor` is the principal angle between the lag-0 and lagged subspace minus that pair's own split-half floor, averaged over the two areas; animals-as-n, Bonferroni across |lag| within a pair.

> **⚠ This table is a power check, not a result.** The split-half floor at 3 dims is ~78°, i.e. two halves of the *same* data at the *same* lag are nearly orthogonal — the 3-dim subspace is not estimable at this N. A lagged angle that fails to exceed that floor means UNMEASURABLE, not stable. Read the CC₁ table instead.

| pair | estimable? | stability width | mean angle @ ±50 ms | floor | Δ | p (Bonf) |
|---|---|---|---|---|---|---|
| CA1-CA3 | **NO** | n/a — not estimable | 78.0° | 80.5° | -2.6° | 1 |
| CA1-DG | **NO** | n/a — not estimable | 78.6° | 77.7° | +0.9° | 1 |
| CA1-RSC | **NO** | n/a — not estimable | 81.1° | 81.9° | -0.8° | 1 |
| CA1-SUB | **NO** | n/a — not estimable | 81.0° | 83.6° | -2.7° | 1 |
| CA1-V1 | **NO** | n/a — not estimable | 73.2° | 82.4° | -9.2° | 0.205 |
| CA3-DG | yes | ≥ 250 ms (censored) | 68.3° | 61.3° | +7.0° | 0.549 |
| RSC-SUB | **NO** | n/a — not estimable | 74.0° | 77.9° | -3.8° | 1 |
| V1-RSC | **NO** | n/a — not estimable | 70.4° | 71.9° | -1.5° | 1 |


### FS-excluded — items 2/4: feedforward (+50 ms) vs feedback (−50 ms)

Positive `Δcc₁` = the first-named area leading is more strongly coupled. `Δgini_conn` uses the CONNECTION-SPECIFIC Gini. `FF/FB angle − floor` is the gate: at or below 0 the two are one subspace read at two delays.

| pair | n | cc₁ FF | cc₁ FB | Δcc₁ | p | Δgini_conn | p | FF/FB angle − floor | p |
|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 12 | 0.032 | 0.029 | +0.003 | 0.764 | +0.000 | 0.987 | +0.7° | 0.937 |
| CA1-CA3 | 7 | 0.089 | 0.123 | -0.034 | 0.212 | -0.003 | 0.924 | +15.3° | 0.137 |
| CA1-DG | 11 | 0.112 | 0.123 | -0.011 | 0.471 | -0.013 | 0.16 | -7.1° | 0.459 |
| CA1-V1 | 13 | 0.034 | 0.042 | -0.008 | 0.153 | +0.014 | 0.448 | -8.7° | 0.241 |
| CA3-DG | 5 | 0.113 | 0.121 | -0.007 | 0.729 | -0.016 | 0.635 | +9.5° | 0.538 |
| CA1-SUB | 7 | 0.029 | 0.018 | +0.011 | 0.193 | +0.003 | 0.879 | +7.4° | 0.497 |
| RSC-SUB | 7 | 0.065 | 0.053 | +0.012 | 0.0447 | +0.011 | 0.3 | -9.1° | 0.245 |
| V1-RSC | 9 | 0.141 | 0.122 | +0.019 | 0.193 | +0.008 | 0.246 | -4.3° | 0.707 |

