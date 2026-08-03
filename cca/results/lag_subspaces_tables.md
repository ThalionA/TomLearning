# Lagged communication subspaces — meeting items 2, 3, 4

Animals-as-n throughout; 8 pairs, no cross-pair correction (per-pair family, STATE.md §3.0 policy). Angles are the LARGEST principal angle over 3 canonical dims, so a subspace that matches on its dominant axis but diverges elsewhere is not scored as stable.

### FS-excluded — item 3: subspace stability across lag (CC₁ only)

`angle − floor` is the principal angle between the lag-0 and lagged subspace minus that pair's own split-half floor, averaged over the two areas; animals-as-n, Bonferroni across |lag| within a pair.

| pair | estimable? | stability width | mean angle @ ±50 ms | floor | Δ | p (Bonf) |
|---|---|---|---|---|---|---|
| CA1-CA3 | yes | 30 ms | 47.8° | 33.8° | +14.1° | 0.248 |
| CA1-DG | yes | 20 ms | 57.2° | 35.1° | +22.4° | 0.0011 |
| CA1-RSC | yes | ≥ 250 ms (censored) | 64.4° | 52.3° | +12.1° | 0.263 |
| CA1-SUB | yes | 150 ms | 71.7° | 42.9° | +27.0° | 0.16 |
| CA1-V1 | yes | ≥ 250 ms (censored) | 48.3° | 48.2° | +3.1° | 1 |
| CA3-DG | yes | ≥ 250 ms (censored) | 48.0° | 37.9° | +12.1° | 1 |
| RSC-SUB | yes | ≥ 250 ms (censored) | 52.0° | 51.0° | +1.0° | 1 |
| V1-RSC | yes | ≥ 250 ms (censored) | 45.6° | 39.3° | +7.5° | 1 |


#### FS-excluded — item 3 depends on the estimability threshold

The gate is an analyst choice with no principled value, and **both directions are biased**: across area-lags corr(floor, angle−floor) is Spearman ρ = −0.57, so no gate drags toward the null (unmeasurable areas contribute large negative deltas) while a tight gate selects low-floor areas and drags toward rotation.

| gate | rotating lags | pairs |
|---|---|---|
| 50° | 15 | CA1-CA3, CA1-DG, CA1-RSC, CA3-DG, RSC-SUB |
| 60° | 14 | CA1-CA3, CA1-DG, CA1-SUB, CA3-DG |
| 70° | 10 | CA1-CA3, CA1-DG, CA1-SUB |
| 80° | 9 | CA1-DG |
| 90° | 1 | CA1-DG |
| no gate | 1 | CA1-DG |

**Robust to every gate, including none: CA1-DG.** That is the only claim this test supports; the lag COUNT is not interpretable.


### FS-excluded — item 3: subspace stability across lag (3 canonical dims)

`angle − floor` is the principal angle between the lag-0 and lagged subspace minus that pair's own split-half floor, averaged over the two areas; animals-as-n, Bonferroni across |lag| within a pair.

> **⚠ This table is a power check, not a result.** The split-half floor at 3 dims is ~78°, i.e. two halves of the *same* data at the *same* lag are nearly orthogonal — the 3-dim subspace is not estimable at this N. A lagged angle that fails to exceed that floor means UNMEASURABLE, not stable. Read the CC₁ table instead.

| pair | estimable? | stability width | mean angle @ ±50 ms | floor | Δ | p (Bonf) |
|---|---|---|---|---|---|---|
| CA1-CA3 | yes | ≥ 250 ms (censored) | 85.9° | 67.7° | +17.9° | 0.24 |
| CA1-DG | yes | ≥ 250 ms (censored) | 84.4° | 57.3° | +24.3° | 0.626 |
| CA1-RSC | yes | ≥ 250 ms (censored) | 84.7° | 57.5° | +27.2° | 1 |
| CA1-SUB | **NO** | n/a — not estimable | nan° | nan° | +nan° | nan |
| CA1-V1 | **NO** | n/a — not estimable | nan° | nan° | +nan° | nan |
| CA3-DG | yes | ≥ 250 ms (censored) | 60.2° | 49.7° | +12.0° | 0.21 |
| RSC-SUB | yes | ≥ 250 ms (censored) | 62.0° | 49.5° | +19.4° | 1 |
| V1-RSC | yes | ≥ 250 ms (censored) | 67.1° | 57.9° | +7.0° | 1 |


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


### FS-excluded — item 2/4: does the FF/FB picture change with learning?

Expert − naive, animals-as-n paired *t*, learners only. FF = +50 ms (first area leads), FB = −50 ms. `asym` = cc₁(FF) − cc₁(FB), the directional asymmetry.

> Read against the session-level gate: if FF and FB are not separable > subspaces, a change in `asym` is a change *within one* subspace, not a > shift between two streams.

| pair | n | Δ cc₁ FF | p | Δ cc₁ FB | p | Δ asym | p | Δ FF/FB angle | p |
|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | -0.011 | 0.399 | +0.008 | 0.337 | -0.020 | 0.352 | -3.1° | 0.605 |
| CA1-CA3 | 6 | +0.000 | 0.996 | +0.023 | 0.392 | -0.023 | 0.236 | +0.6° | 0.939 |
| CA1-DG | 8 | +0.022 | 0.104 | +0.004 | 0.741 | +0.017 | 0.387 | +1.2° | 0.861 |
| CA1-V1 | 10 | -0.007 | 0.655 | -0.016 | 0.295 | +0.009 | 0.715 | -5.6° | 0.482 |
| CA3-DG | 4 | +0.018 | 0.469 | -0.022 | 0.391 | +0.040 | 0.1 | +8.4° | 0.608 |
| CA1-SUB | 4 | +0.001 | 0.889 | +0.010 | 0.352 | -0.009 | 0.609 | +3.2° | 0.79 |
| RSC-SUB | 4 | -0.004 | 0.863 | -0.008 | 0.706 | +0.004 | 0.835 | +15.8° | 0.122 |
| V1-RSC | 6 | -0.027 | 0.226 | -0.021 | 0.068 | -0.006 | 0.775 | -18.9° | 0.274 |

### FS-included — item 3: subspace stability across lag (CC₁ only)

`angle − floor` is the principal angle between the lag-0 and lagged subspace minus that pair's own split-half floor, averaged over the two areas; animals-as-n, Bonferroni across |lag| within a pair.

| pair | estimable? | stability width | mean angle @ ±50 ms | floor | Δ | p (Bonf) |
|---|---|---|---|---|---|---|
| CA1-CA3 | yes | 200 ms | 51.6° | 30.8° | +20.8° | 0.0709 |
| CA1-DG | yes | 20 ms | 56.7° | 35.0° | +20.2° | 0.00755 |
| CA1-RSC | yes | 70 ms | 62.0° | 52.8° | +9.3° | 0.272 |
| CA1-SUB | yes | ≥ 250 ms (censored) | 74.6° | 48.6° | +23.0° | 0.829 |
| CA1-V1 | yes | ≥ 250 ms (censored) | 43.8° | 47.8° | -3.3° | 1 |
| CA3-DG | yes | ≥ 250 ms (censored) | 51.9° | 33.1° | +18.8° | 0.872 |
| RSC-SUB | yes | ≥ 250 ms (censored) | 56.3° | 51.1° | +3.6° | 1 |
| V1-RSC | yes | ≥ 250 ms (censored) | 44.6° | 34.9° | +9.7° | 1 |


#### FS-included — item 3 depends on the estimability threshold

The gate is an analyst choice with no principled value, and **both directions are biased**: across area-lags corr(floor, angle−floor) is Spearman ρ = −0.57, so no gate drags toward the null (unmeasurable areas contribute large negative deltas) while a tight gate selects low-floor areas and drags toward rotation.

| gate | rotating lags | pairs |
|---|---|---|
| 50° | 20 | CA1-CA3, CA1-DG, CA1-RSC, CA3-DG |
| 60° | 20 | CA1-CA3, CA1-DG, CA1-RSC, CA1-SUB |
| 70° | 16 | CA1-CA3, CA1-DG, CA1-RSC |
| 80° | 10 | CA1-CA3, CA1-DG |
| 90° | 4 | CA1-CA3, CA1-DG |
| no gate | 4 | CA1-CA3, CA1-DG |

**Robust to every gate, including none: CA1-CA3, CA1-DG.** That is the only claim this test supports; the lag COUNT is not interpretable.


### FS-included — item 3: subspace stability across lag (3 canonical dims)

`angle − floor` is the principal angle between the lag-0 and lagged subspace minus that pair's own split-half floor, averaged over the two areas; animals-as-n, Bonferroni across |lag| within a pair.

> **⚠ This table is a power check, not a result.** The split-half floor at 3 dims is ~78°, i.e. two halves of the *same* data at the *same* lag are nearly orthogonal — the 3-dim subspace is not estimable at this N. A lagged angle that fails to exceed that floor means UNMEASURABLE, not stable. Read the CC₁ table instead.

| pair | estimable? | stability width | mean angle @ ±50 ms | floor | Δ | p (Bonf) |
|---|---|---|---|---|---|---|
| CA1-CA3 | yes | ≥ 250 ms (censored) | 68.8° | 57.4° | +7.7° | 1 |
| CA1-DG | yes | ≥ 250 ms (censored) | 78.8° | 49.5° | +27.0° | 0.647 |
| CA1-RSC | yes | ≥ 250 ms (censored) | 80.9° | 59.2° | +21.7° | 0.495 |
| CA1-SUB | **NO** | n/a — not estimable | nan° | nan° | +nan° | nan |
| CA1-V1 | **NO** | n/a — not estimable | nan° | nan° | +nan° | nan |
| CA3-DG | yes | ≥ 250 ms (censored) | 67.0° | 45.7° | +18.1° | 1 |
| RSC-SUB | yes | ≥ 250 ms (censored) | 53.5° | 62.8° | -11.4° | 1 |
| V1-RSC | yes | 200 ms | 57.9° | 54.6° | +2.0° | 1 |


### FS-included — items 2/4: feedforward (+50 ms) vs feedback (−50 ms)

Positive `Δcc₁` = the first-named area leading is more strongly coupled. `Δgini_conn` uses the CONNECTION-SPECIFIC Gini. `FF/FB angle − floor` is the gate: at or below 0 the two are one subspace read at two delays.

| pair | n | cc₁ FF | cc₁ FB | Δcc₁ | p | Δgini_conn | p | FF/FB angle − floor | p |
|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 12 | 0.072 | 0.049 | +0.023 | 0.0154 | -0.002 | 0.836 | +9.1° | 0.075 |
| CA1-CA3 | 9 | 0.182 | 0.215 | -0.033 | 0.143 | +0.022 | 0.359 | +2.5° | 0.848 |
| CA1-DG | 11 | 0.096 | 0.107 | -0.012 | 0.426 | -0.012 | 0.406 | -7.6° | 0.414 |
| CA1-V1 | 13 | 0.041 | 0.046 | -0.005 | 0.374 | +0.031 | 0.121 | -11.6° | 0.196 |
| CA3-DG | 5 | 0.134 | 0.157 | -0.022 | 0.149 | -0.032 | 0.188 | +13.2° | 0.143 |
| CA1-SUB | 7 | 0.028 | 0.022 | +0.005 | 0.588 | +0.010 | 0.515 | +3.9° | 0.664 |
| RSC-SUB | 7 | 0.075 | 0.062 | +0.013 | 0.116 | +0.012 | 0.266 | -11.3° | 0.0622 |
| V1-RSC | 9 | 0.188 | 0.158 | +0.030 | 0.0554 | +0.008 | 0.195 | -5.0° | 0.693 |


### FS-included — item 2/4: does the FF/FB picture change with learning?

Expert − naive, animals-as-n paired *t*, learners only. FF = +50 ms (first area leads), FB = −50 ms. `asym` = cc₁(FF) − cc₁(FB), the directional asymmetry.

> Read against the session-level gate: if FF and FB are not separable > subspaces, a change in `asym` is a change *within one* subspace, not a > shift between two streams.

| pair | n | Δ cc₁ FF | p | Δ cc₁ FB | p | Δ asym | p | Δ FF/FB angle | p |
|---|---|---|---|---|---|---|---|---|---|
| CA1-RSC | 8 | -0.012 | 0.375 | +0.001 | 0.918 | -0.013 | 0.422 | +1.7° | 0.814 |
| CA1-CA3 | 7 | +0.020 | 0.556 | +0.031 | 0.29 | -0.012 | 0.365 | +4.9° | 0.656 |
| CA1-DG | 8 | +0.017 | 0.243 | -0.004 | 0.644 | +0.021 | 0.235 | +4.2° | 0.142 |
| CA1-V1 | 10 | +0.002 | 0.874 | -0.015 | 0.39 | +0.017 | 0.358 | +5.3° | 0.403 |
| CA3-DG | 4 | +0.053 | 0.0745 | -0.005 | 0.787 | +0.058 | 0.0682 | -2.4° | 0.8 |
| CA1-SUB | 4 | +0.004 | 0.619 | +0.005 | 0.835 | -0.001 | 0.967 | +6.7° | 0.0881 |
| RSC-SUB | 4 | -0.008 | 0.625 | -0.001 | 0.956 | -0.007 | 0.778 | +6.6° | 0.639 |
| V1-RSC | 6 | -0.064 | 0.0988 | -0.053 | 0.0556 | -0.011 | 0.607 | -12.9° | 0.455 |
