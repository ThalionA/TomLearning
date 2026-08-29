# Contribution × spatial reliability — are communicating units spatially special?

Generated 2026-08-29 by `scripts/analyze_contrib_reliability.py` (full prints:
`results/contrib_reliability_fs{excl,incl}.txt`; per-fit CSVs `contrib_reliability_bin10*.csv`,
per-unit join `contrib_reliability_units_bin10*.csv`).

**Method.** Per (animal, pair, epoch, area): Spearman across units of the connection-specific
contribution `contrib_conn` (run_epochs, one pCCA per animal × pair × epoch, bin10 σ2.5, K=20)
vs spatial reliability = mean trial-to-trial Pearson r of the unit's spatially-binned map within
±2 trials (`spatial_reliability.trial_map_reliability`), averaged over the SAME trials the
epoch's CCA was fit on. Aggregation: Fisher-z mean over epochs within animal → animals-as-n
one-sample t + Wilcoxon vs 0, per (pair, area). Learners only (epochs need an LP).

**Controls.** (1) `ratepart` = log mean rate rank-partialled out of both variables — rate
predicts reliability at rho ≈ +0.4–0.5 in every area, so roughly half the raw link is rate.
(2) The area-intrinsic `contrib` (partner-invariant) is n.s. nearly everywhere — the link is
specific to the connection-weighted contribution, not 'loud units load on everything'.

**Caveats.** Per-pair family, no cross-pair correction (project policy); 16 (pair, area) cells
per metric per FS. n = 4 pairs (CA3-DG, CA1-SUB, RSC-SUB) sit on the Wilcoxon floor (0.125) —
descriptive only. At n = 6 the floor is 0.0312, i.e. a starred W there means 6/6 same sign.
Reliability and CCA input share the same spikes: rate-partialling handles the SNR confound only
to the extent log-rate captures SNR.

**FS-robust survivors after rate-partialling** (starred in both FS conditions):
**CA1-RSC RSC side** (+0.256/+0.255, W p=0.016 both), **CA1-V1 V1 side** (+0.258/+0.170,
W p=0.002/0.037), **V1-RSC RSC side** (+0.312/+0.411, W p=0.031 both; V1 side positive 6/6,
t starred both FS). The pattern: the CORTICAL/target-side units' reliability tracks their
contribution; the CA1-side link never survives the rate partial.

## FS-excluded — contrib_conn vs reliability (raw)

```
  CA1-RSC   CA1  n=8  mean_rho= +0.221 med= +0.269 up=6/8 | t p=0.0169* | W p=0.0391*
  CA1-RSC   RSC  n=8  mean_rho= +0.379 med= +0.346 up=8/8 | t p=0.000389* | W p=0.00781*
  CA1-CA3   CA1  n=6  mean_rho= +0.055 med= +0.225 up=4/6 | t p=0.724  | W p=1 
  CA1-CA3   CA3  n=6  mean_rho= +0.216 med= +0.227 up=5/6 | t p=0.0334* | W p=0.0625 
  CA1-DG    CA1  n=8  mean_rho= +0.030 med= -0.004 up=4/8 | t p=0.697  | W p=0.844 
  CA1-DG    DG   n=8  mean_rho= +0.300 med= +0.367 up=7/8 | t p=0.0208* | W p=0.0547 
  CA1-V1    CA1  n=10 mean_rho= +0.192 med= +0.295 up=7/10 | t p=0.103  | W p=0.105 
  CA1-V1    V1   n=10 mean_rho= +0.353 med= +0.375 up=10/10 | t p=0.000108* | W p=0.00195*
  CA3-DG    CA3  n=4  mean_rho= +0.203 med= +0.256 up=3/4 | t p=0.124  | W p=0.25 
  CA3-DG    DG   n=4  mean_rho= +0.351 med= +0.312 up=4/4 | t p=0.00727* | W p=0.125 
  CA1-SUB   CA1  n=4  mean_rho= +0.065 med= +0.003 up=3/4 | t p=0.79  | W p=0.625 
  CA1-SUB   SUB  n=4  mean_rho= +0.404 med= +0.395 up=4/4 | t p=0.0811  | W p=0.125 
  RSC-SUB   RSC  n=4  mean_rho= +0.411 med= +0.389 up=4/4 | t p=0.00298* | W p=0.125 
  RSC-SUB   SUB  n=4  mean_rho= +0.397 med= +0.402 up=4/4 | t p=0.00901* | W p=0.125 
  V1-RSC    V1   n=6  mean_rho= +0.235 med= +0.270 up=5/6 | t p=0.0142* | W p=0.0625 
  V1-RSC    RSC  n=6  mean_rho= +0.362 med= +0.404 up=5/6 | t p=0.0154* | W p=0.0625 
```

## FS-excluded — contrib_conn vs reliability (rate-partialled)

```
  CA1-RSC   CA1  n=8  mean_rho= +0.137 med= +0.138 up=6/8 | t p=0.156  | W p=0.148 
  CA1-RSC   RSC  n=8  mean_rho= +0.256 med= +0.258 up=7/8 | t p=0.0111* | W p=0.0156*
  CA1-CA3   CA1  n=6  mean_rho= -0.156 med= -0.142 up=3/6 | t p=0.253  | W p=0.438 
  CA1-CA3   CA3  n=6  mean_rho= -0.030 med= +0.013 up=4/6 | t p=0.764  | W p=1 
  CA1-DG    CA1  n=8  mean_rho= -0.180 med= -0.208 up=1/8 | t p=0.0176* | W p=0.0234*
  CA1-DG    DG   n=8  mean_rho= +0.096 med= +0.094 up=6/8 | t p=0.105  | W p=0.148 
  CA1-V1    CA1  n=10 mean_rho= +0.115 med= +0.203 up=7/10 | t p=0.228  | W p=0.232 
  CA1-V1    V1   n=10 mean_rho= +0.258 med= +0.263 up=10/10 | t p=0.000202* | W p=0.00195*
  CA3-DG    CA3  n=4  mean_rho= -0.007 med= -0.027 up=1/4 | t p=0.953  | W p=0.875 
  CA3-DG    DG   n=4  mean_rho= +0.040 med= -0.002 up=2/4 | t p=0.675  | W p=1 
  CA1-SUB   CA1  n=4  mean_rho= -0.005 med= -0.077 up=2/4 | t p=0.978  | W p=1 
  CA1-SUB   SUB  n=4  mean_rho= +0.298 med= +0.350 up=3/4 | t p=0.268  | W p=0.375 
  RSC-SUB   RSC  n=4  mean_rho= +0.263 med= +0.281 up=4/4 | t p=0.00355* | W p=0.125 
  RSC-SUB   SUB  n=4  mean_rho= +0.157 med= +0.362 up=3/4 | t p=0.595  | W p=0.875 
  V1-RSC    V1   n=6  mean_rho= +0.152 med= +0.187 up=5/6 | t p=0.0149* | W p=0.0625 
  V1-RSC    RSC  n=6  mean_rho= +0.312 med= +0.313 up=6/6 | t p=0.00292* | W p=0.0312*
```

## FS-excluded — area-intrinsic contrib control (raw)

```
  CA1-RSC   CA1  n=8  mean_rho= -0.107 med= -0.277 up=3/8 | t p=0.403  | W p=0.461 
  CA1-RSC   RSC  n=8  mean_rho= +0.190 med= +0.241 up=6/8 | t p=0.104  | W p=0.148 
  CA1-CA3   CA1  n=6  mean_rho= -0.042 med= +0.003 up=3/6 | t p=0.722  | W p=0.844 
  CA1-CA3   CA3  n=6  mean_rho= -0.165 med= -0.134 up=2/6 | t p=0.248  | W p=0.312 
  CA1-DG    CA1  n=8  mean_rho= +0.019 med= +0.109 up=5/8 | t p=0.887  | W p=1 
  CA1-DG    DG   n=8  mean_rho= +0.029 med= +0.039 up=4/8 | t p=0.681  | W p=0.742 
  CA1-V1    CA1  n=10 mean_rho= -0.071 med= -0.066 up=5/10 | t p=0.471  | W p=0.432 
  CA1-V1    V1   n=10 mean_rho= +0.156 med= +0.167 up=7/10 | t p=0.0411* | W p=0.0645 
  CA3-DG    CA3  n=4  mean_rho= -0.118 med= -0.063 up=2/4 | t p=0.574  | W p=0.875 
  CA3-DG    DG   n=4  mean_rho= +0.063 med= +0.018 up=2/4 | t p=0.511  | W p=0.875 
  CA1-SUB   CA1  n=4  mean_rho= -0.213 med= -0.333 up=1/4 | t p=0.373  | W p=0.625 
  CA1-SUB   SUB  n=4  mean_rho= +0.144 med= +0.190 up=3/4 | t p=0.129  | W p=0.25 
  RSC-SUB   RSC  n=4  mean_rho= +0.411 med= +0.406 up=4/4 | t p=0.00235* | W p=0.125 
  RSC-SUB   SUB  n=4  mean_rho= +0.140 med= +0.176 up=3/4 | t p=0.141  | W p=0.25 
  V1-RSC    V1   n=6  mean_rho= -0.040 med= +0.073 up=4/6 | t p=0.778  | W p=1 
  V1-RSC    RSC  n=6  mean_rho= +0.236 med= +0.289 up=5/6 | t p=0.136  | W p=0.219 
```

## FS-excluded — rate vs reliability (the confound's size)

```
  CA1-RSC   CA1  n=8  mean_rho= +0.387 med= +0.459 up=7/8 | t p=0.000959* | W p=0.0156*
  CA1-RSC   RSC  n=8  mean_rho= +0.406 med= +0.483 up=7/8 | t p=0.00555* | W p=0.0156*
  CA1-CA3   CA1  n=6  mean_rho= +0.366 med= +0.439 up=5/6 | t p=0.0123* | W p=0.0625 
  CA1-CA3   CA3  n=6  mean_rho= +0.453 med= +0.392 up=6/6 | t p=0.00267* | W p=0.0312*
  CA1-DG    CA1  n=8  mean_rho= +0.412 med= +0.471 up=7/8 | t p=0.00053* | W p=0.0156*
  CA1-DG    DG   n=8  mean_rho= +0.569 med= +0.599 up=8/8 | t p=7.5e-05* | W p=0.00781*
  CA1-V1    CA1  n=10 mean_rho= +0.495 med= +0.506 up=10/10 | t p=1.74e-07* | W p=0.00195*
  CA1-V1    V1   n=10 mean_rho= +0.476 med= +0.493 up=10/10 | t p=2.47e-05* | W p=0.00195*
  CA3-DG    CA3  n=4  mean_rho= +0.436 med= +0.392 up=4/4 | t p=0.0297* | W p=0.125 
  CA3-DG    DG   n=4  mean_rho= +0.462 med= +0.501 up=4/4 | t p=0.0135* | W p=0.125 
  CA1-SUB   CA1  n=4  mean_rho= +0.481 med= +0.506 up=4/4 | t p=0.00141* | W p=0.125 
  CA1-SUB   SUB  n=4  mean_rho= +0.401 med= +0.277 up=4/4 | t p=0.0924  | W p=0.125 
  RSC-SUB   RSC  n=4  mean_rho= +0.498 med= +0.483 up=4/4 | t p=0.00825* | W p=0.125 
  RSC-SUB   SUB  n=4  mean_rho= +0.401 med= +0.277 up=4/4 | t p=0.0924  | W p=0.125 
  V1-RSC    V1   n=6  mean_rho= +0.420 med= +0.447 up=6/6 | t p=0.00345* | W p=0.0312*
  V1-RSC    RSC  n=6  mean_rho= +0.485 med= +0.483 up=6/6 | t p=0.000644* | W p=0.0312*
```

## FS-included — contrib_conn vs reliability (raw)

```
  CA1-RSC   CA1  n=8  mean_rho= +0.250 med= +0.308 up=6/8 | t p=0.0257* | W p=0.0391*
  CA1-RSC   RSC  n=8  mean_rho= +0.358 med= +0.431 up=8/8 | t p=0.000377* | W p=0.00781*
  CA1-CA3   CA1  n=7  mean_rho= +0.188 med= +0.223 up=6/7 | t p=0.122  | W p=0.297 
  CA1-CA3   CA3  n=7  mean_rho= +0.253 med= +0.303 up=6/7 | t p=0.00693* | W p=0.0312*
  CA1-DG    CA1  n=8  mean_rho= +0.220 med= +0.257 up=7/8 | t p=0.0259* | W p=0.0781 
  CA1-DG    DG   n=8  mean_rho= +0.354 med= +0.429 up=7/8 | t p=0.00536* | W p=0.0156*
  CA1-V1    CA1  n=10 mean_rho= +0.320 med= +0.397 up=9/10 | t p=0.0174* | W p=0.0645 
  CA1-V1    V1   n=10 mean_rho= +0.266 med= +0.362 up=9/10 | t p=0.0137* | W p=0.0195*
  CA3-DG    CA3  n=4  mean_rho= +0.139 med= +0.171 up=3/4 | t p=0.2  | W p=0.25 
  CA3-DG    DG   n=4  mean_rho= +0.406 med= +0.402 up=4/4 | t p=0.000601* | W p=0.125 
  CA1-SUB   CA1  n=4  mean_rho= +0.334 med= +0.527 up=3/4 | t p=0.259  | W p=0.25 
  CA1-SUB   SUB  n=4  mean_rho= +0.472 med= +0.410 up=4/4 | t p=0.0149* | W p=0.125 
  RSC-SUB   RSC  n=4  mean_rho= +0.457 med= +0.440 up=4/4 | t p=0.00135* | W p=0.125 
  RSC-SUB   SUB  n=4  mean_rho= +0.433 med= +0.421 up=4/4 | t p=0.00175* | W p=0.125 
  V1-RSC    V1   n=6  mean_rho= +0.263 med= +0.299 up=6/6 | t p=0.00517* | W p=0.0312*
  V1-RSC    RSC  n=6  mean_rho= +0.526 med= +0.536 up=6/6 | t p=1.38e-05* | W p=0.0312*
```

## FS-included — contrib_conn vs reliability (rate-partialled)

```
  CA1-RSC   CA1  n=8  mean_rho= +0.089 med= +0.137 up=5/8 | t p=0.361  | W p=0.312 
  CA1-RSC   RSC  n=8  mean_rho= +0.255 med= +0.290 up=7/8 | t p=0.00241* | W p=0.0156*
  CA1-CA3   CA1  n=7  mean_rho= -0.155 med= -0.162 up=3/7 | t p=0.146  | W p=0.297 
  CA1-CA3   CA3  n=7  mean_rho= +0.044 med= +0.022 up=4/7 | t p=0.612  | W p=0.688 
  CA1-DG    CA1  n=8  mean_rho= -0.063 med= -0.044 up=3/8 | t p=0.327  | W p=0.547 
  CA1-DG    DG   n=8  mean_rho= +0.132 med= +0.164 up=6/8 | t p=0.0516  | W p=0.0781 
  CA1-V1    CA1  n=10 mean_rho= +0.176 med= +0.221 up=8/10 | t p=0.106  | W p=0.16 
  CA1-V1    V1   n=10 mean_rho= +0.170 med= +0.259 up=8/10 | t p=0.0305* | W p=0.0371*
  CA3-DG    CA3  n=4  mean_rho= -0.099 med= -0.051 up=2/4 | t p=0.363  | W p=0.625 
  CA3-DG    DG   n=4  mean_rho= +0.049 med= -0.029 up=2/4 | t p=0.648  | W p=1 
  CA1-SUB   CA1  n=4  mean_rho= +0.197 med= +0.238 up=3/4 | t p=0.284  | W p=0.375 
  CA1-SUB   SUB  n=4  mean_rho= +0.352 med= +0.345 up=3/4 | t p=0.105  | W p=0.25 
  RSC-SUB   RSC  n=4  mean_rho= +0.290 med= +0.293 up=4/4 | t p=0.0033* | W p=0.125 
  RSC-SUB   SUB  n=4  mean_rho= +0.279 med= +0.399 up=3/4 | t p=0.173  | W p=0.25 
  V1-RSC    V1   n=6  mean_rho= +0.174 med= +0.213 up=6/6 | t p=0.00682* | W p=0.0312*
  V1-RSC    RSC  n=6  mean_rho= +0.411 med= +0.395 up=6/6 | t p=8.38e-05* | W p=0.0312*
```

## FS-included — area-intrinsic contrib control (raw)

```
  CA1-RSC   CA1  n=8  mean_rho= -0.055 med= +0.097 up=4/8 | t p=0.673  | W p=0.641 
  CA1-RSC   RSC  n=8  mean_rho= +0.175 med= +0.238 up=6/8 | t p=0.15  | W p=0.25 
  CA1-CA3   CA1  n=7  mean_rho= +0.122 med= +0.239 up=5/7 | t p=0.292  | W p=0.375 
  CA1-CA3   CA3  n=7  mean_rho= -0.059 med= -0.085 up=3/7 | t p=0.667  | W p=0.812 
  CA1-DG    CA1  n=8  mean_rho= +0.126 med= +0.217 up=7/8 | t p=0.239  | W p=0.195 
  CA1-DG    DG   n=8  mean_rho= +0.030 med= +0.036 up=4/8 | t p=0.661  | W p=0.641 
  CA1-V1    CA1  n=10 mean_rho= +0.103 med= +0.205 up=8/10 | t p=0.219  | W p=0.131 
  CA1-V1    V1   n=10 mean_rho= +0.123 med= +0.152 up=7/10 | t p=0.0926  | W p=0.105 
  CA3-DG    CA3  n=4  mean_rho= -0.110 med= +0.020 up=2/4 | t p=0.639  | W p=1 
  CA3-DG    DG   n=4  mean_rho= +0.045 med= +0.024 up=2/4 | t p=0.635  | W p=0.875 
  CA1-SUB   CA1  n=4  mean_rho= +0.055 med= +0.095 up=2/4 | t p=0.711  | W p=0.875 
  CA1-SUB   SUB  n=4  mean_rho= +0.144 med= +0.183 up=3/4 | t p=0.135  | W p=0.25 
  RSC-SUB   RSC  n=4  mean_rho= +0.435 med= +0.436 up=4/4 | t p=0.00138* | W p=0.125 
  RSC-SUB   SUB  n=4  mean_rho= +0.142 med= +0.179 up=3/4 | t p=0.138  | W p=0.25 
  V1-RSC    V1   n=6  mean_rho= -0.073 med= +0.012 up=3/6 | t p=0.586  | W p=1 
  V1-RSC    RSC  n=6  mean_rho= +0.248 med= +0.297 up=5/6 | t p=0.109  | W p=0.156 
```

## FS-included — rate vs reliability (the confound's size)

```
  CA1-RSC   CA1  n=8  mean_rho= +0.439 med= +0.482 up=8/8 | t p=0.000359* | W p=0.00781*
  CA1-RSC   RSC  n=8  mean_rho= +0.424 med= +0.503 up=7/8 | t p=0.00187* | W p=0.0156*
  CA1-CA3   CA1  n=7  mean_rho= +0.417 med= +0.477 up=7/7 | t p=0.00278* | W p=0.0156*
  CA1-CA3   CA3  n=7  mean_rho= +0.396 med= +0.372 up=7/7 | t p=0.00121* | W p=0.0156*
  CA1-DG    CA1  n=8  mean_rho= +0.419 med= +0.482 up=8/8 | t p=0.000737* | W p=0.00781*
  CA1-DG    DG   n=8  mean_rho= +0.569 med= +0.599 up=8/8 | t p=7.5e-05* | W p=0.00781*
  CA1-V1    CA1  n=10 mean_rho= +0.500 med= +0.506 up=10/10 | t p=4.83e-07* | W p=0.00195*
  CA1-V1    V1   n=10 mean_rho= +0.494 med= +0.522 up=10/10 | t p=7.13e-06* | W p=0.00195*
  CA3-DG    CA3  n=4  mean_rho= +0.377 med= +0.423 up=4/4 | t p=0.0244* | W p=0.125 
  CA3-DG    DG   n=4  mean_rho= +0.462 med= +0.501 up=4/4 | t p=0.0135* | W p=0.125 
  CA1-SUB   CA1  n=4  mean_rho= +0.522 med= +0.519 up=4/4 | t p=0.00153* | W p=0.125 
  CA1-SUB   SUB  n=4  mean_rho= +0.401 med= +0.277 up=4/4 | t p=0.0924  | W p=0.125 
  RSC-SUB   RSC  n=4  mean_rho= +0.533 med= +0.522 up=4/4 | t p=0.00208* | W p=0.125 
  RSC-SUB   SUB  n=4  mean_rho= +0.401 med= +0.277 up=4/4 | t p=0.0924  | W p=0.125 
  V1-RSC    V1   n=6  mean_rho= +0.454 med= +0.461 up=6/6 | t p=0.000474* | W p=0.0312*
  V1-RSC    RSC  n=6  mean_rho= +0.500 med= +0.508 up=6/6 | t p=0.000144* | W p=0.0312*
```

