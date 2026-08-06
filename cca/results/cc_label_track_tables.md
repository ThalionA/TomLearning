# Item 2 — FF/FB CCs labelled once, tracked across learning

Canonical dimensions labelled by the sign of their IFI on the **whole-session fit**, then followed through naive/intermediate/expert through **frozen axes** — so a change is a change in the activity, not in the fit, and the label cannot be re-drawn by noise each epoch.

> Correlations are IN-SAMPLE by construction (the frozen fit saw every > epoch). They are contrast statistics, not coupling strengths.

### FS-excluded — do FF- and FB-labelled CCs change differently with learning?

223 FF and 247 FB canonical dimensions (significant CCs only), labelled once on the whole-session fit and followed through the epochs on frozen axes. Animals-as-n: an animal's CCs of one label are averaged before testing. `Δ` is expert − naive; `interaction` is ΔFF − ΔFB, which is the question item 2 asks.

| pair | metric | n | Δ FF | p | Δ FB | p | interaction | p |
|---|---|---|---|---|---|---|---|---|
| CA1-RSC | peak r | 7 | +0.0182 | 0.0112 | +0.0164 | 0.0121 | +0.00287 | 0.603 |
| CA1-RSC | IFI | 7 | +0.0351 | 0.586 | -0.0291 | 0.727 | +0.0771 | 0.461 |
| CA1-RSC | peak lag (ms) | 7 | +10.5 | 0.594 | -26.5 | 0.498 | +40.3 | 0.411 |
| CA1-CA3 | peak r | 6 | +0.0135 | 0.137 | +0.00995 | 0.0741 | +0.00359 | 0.738 |
| CA1-CA3 | IFI | 6 | +0.174 | 0.0765 | -0.0945 | 0.436 | +0.269 | 0.181 |
| CA1-CA3 | peak lag (ms) | 6 | +32.4 | 0.322 | +7.78 | 0.663 | +24.6 | 0.491 |
| CA1-DG | peak r | 8 | +0.00588 | 0.105 | +0.00657 | 0.129 | -0.000689 | 0.865 |
| CA1-DG | IFI | 8 | -0.0503 | 0.545 | +0.101 | 0.238 | -0.152 | 0.277 |
| CA1-DG | peak lag (ms) | 8 | +20.6 | 0.281 | -6.51 | 0.645 | +27.1 | 0.29 |
| CA1-V1 | peak r | 10 | +0.0172 | 0.00303 | +0.0239 | 0.00229 | -0.00677 | 0.177 |
| CA1-V1 | IFI | 10 | +0.0454 | 0.0284 | -0.0203 | 0.256 | +0.0657 | 0.0191 |
| CA1-V1 | peak lag (ms) | 10 | +48.2 | 0.0981 | +46.1 | 0.399 | +2.12 | 0.972 |
| CA3-DG | peak r | 4 | +0.0166 | 0.346 | +0.0135 | 0.00033 | +0.0031 | 0.852 |
| CA3-DG | IFI | 4 | +0.185 | 0.355 | -0.044 | 0.617 | +0.229 | 0.124 |
| CA3-DG | peak lag (ms) | 4 | +4.17 | 0.912 | -30.3 | 0.239 | +34.5 | 0.542 |
| CA1-SUB | peak r | 2 | +0.0138 | 0.485 | +0.00874 | 0.268 | +0.0124 | 0.394 |
| CA1-SUB | IFI | 2 | +0.414 | 0.284 | +0.204 | 0.189 | +0.0436 | 0.305 |
| CA1-SUB | peak lag (ms) | 2 | -1.67 | 0.5 | -46.2 | 0.444 | +80.2 | 0.616 |
| RSC-SUB | peak r | 4 | +0.00642 | 0.15 | +0.0185 | 0.0674 | -0.012 | 0.252 |
| RSC-SUB | IFI | 4 | -0.233 | 0.164 | +0.0425 | 0.554 | -0.276 | 0.0458 |
| RSC-SUB | peak lag (ms) | 4 | -53 | 0.29 | +6.88 | 0.667 | -59.9 | 0.245 |
| V1-RSC | peak r | 6 | +0.0286 | 0.0757 | +0.031 | 0.0261 | -0.00236 | 0.731 |
| V1-RSC | IFI | 6 | -0.000929 | 0.939 | +0.0668 | 0.308 | -0.0678 | 0.264 |
| V1-RSC | peak lag (ms) | 6 | -14.7 | 0.191 | -53.3 | 0.0459 | +38.7 | 0.0627 |

**Label persistence (leave-epoch-out).** Across 470 labelled CCs in 12 animals, an epoch's own IFI sign matches the label computed WITHOUT that epoch 62% of the time (chance 50 %, p = 0.000735; animals-as-n). Excluding the scored epoch is what makes 50 % the right baseline — against the whole-session label the construction floor is ~60 %.


### FS-included — do FF- and FB-labelled CCs change differently with learning?

234 FF and 260 FB canonical dimensions (significant CCs only), labelled once on the whole-session fit and followed through the epochs on frozen axes. Animals-as-n: an animal's CCs of one label are averaged before testing. `Δ` is expert − naive; `interaction` is ΔFF − ΔFB, which is the question item 2 asks.

| pair | metric | n | Δ FF | p | Δ FB | p | interaction | p |
|---|---|---|---|---|---|---|---|---|
| CA1-RSC | peak r | 7 | +0.0166 | 0.0272 | +0.0198 | 0.0863 | -0.0069 | 0.485 |
| CA1-RSC | IFI | 7 | +0.0928 | 0.25 | +0.366 | 0.0538 | -0.231 | 0.296 |
| CA1-RSC | peak lag (ms) | 7 | -0.299 | 0.989 | -0.816 | 0.986 | +4.62 | 0.941 |
| CA1-CA3 | peak r | 7 | +0.0157 | 0.118 | +0.00641 | 0.4 | +0.00932 | 0.546 |
| CA1-CA3 | IFI | 7 | +0.129 | 0.336 | -0.0645 | 0.559 | +0.193 | 0.342 |
| CA1-CA3 | peak lag (ms) | 7 | -4.16 | 0.831 | +13.6 | 0.59 | -17.8 | 0.61 |
| CA1-DG | peak r | 8 | +0.00259 | 0.452 | +0.00242 | 0.478 | +0.000172 | 0.962 |
| CA1-DG | IFI | 8 | +0.0472 | 0.562 | +0.0705 | 0.45 | -0.0233 | 0.804 |
| CA1-DG | peak lag (ms) | 8 | +1.46 | 0.908 | -19.4 | 0.369 | +20.8 | 0.284 |
| CA1-V1 | peak r | 10 | +0.0185 | 0.000385 | +0.0254 | 0.00249 | -0.00682 | 0.375 |
| CA1-V1 | IFI | 10 | -0.0323 | 0.6 | -0.0611 | 0.19 | +0.0289 | 0.762 |
| CA1-V1 | peak lag (ms) | 10 | +50.8 | 0.161 | -0.0218 | 0.998 | +50.8 | 0.155 |
| CA3-DG | peak r | 4 | +0.0119 | 0.0258 | +0.0173 | 0.00421 | -0.00534 | 0.0637 |
| CA3-DG | IFI | 4 | +0.053 | 0.457 | +0.0234 | 0.74 | +0.0296 | 0.797 |
| CA3-DG | peak lag (ms) | 4 | -54.2 | 0.185 | -1.79 | 0.611 | -52.4 | 0.174 |
| CA1-SUB | peak r | 3 | +0.0139 | 0.334 | +0.00969 | 0.27 | +0.00402 | 0.594 |
| CA1-SUB | IFI | 3 | +0.82 | 0.0098 | +0.19 | 0.329 | +0.514 | 0.0717 |
| CA1-SUB | peak lag (ms) | 3 | +59.3 | 0.46 | -58.6 | 0.36 | +137 | 0.208 |
| RSC-SUB | peak r | 3 | -8.06e-05 | 0.992 | +0.0168 | 0.0627 | -0.0115 | 0.258 |
| RSC-SUB | IFI | 3 | -0.0613 | 0.276 | -0.0138 | 0.595 | -0.0477 | 0.125 |
| RSC-SUB | peak lag (ms) | 3 | +86.7 | 0.541 | -20.3 | 0.434 | +117 | 0.438 |
| V1-RSC | peak r | 6 | +0.028 | 0.0378 | +0.0309 | 0.0344 | -0.00294 | 0.559 |
| V1-RSC | IFI | 6 | -0.00817 | 0.73 | -0.0244 | 0.328 | +0.0162 | 0.473 |
| V1-RSC | peak lag (ms) | 6 | -40.9 | 0.0103 | -45.5 | 0.0859 | +4.59 | 0.794 |

**Label persistence (leave-epoch-out).** Across 494 labelled CCs in 12 animals, an epoch's own IFI sign matches the label computed WITHOUT that epoch 66% of the time (chance 50 %, p = 4.78e-06; animals-as-n). Excluding the scored epoch is what makes 50 % the right baseline — against the whole-session label the construction floor is ~60 %.

