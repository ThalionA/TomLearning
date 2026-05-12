# Figure legends — Temporal vs Spatial CCA comparison

Self-contained legends for each comparison figure. Each legend states
panel layout, what the axes represent, the underlying unit of
observation, error/statistic conventions, and `n`. Companion to
`METHODS.md`.

## Figure 1 — Combined paired scatter (real values)

**File**: `HC_V1_figures/CCA_Compare/TemporalVsSpatial_combined.svg`.

Trial-by-trial agreement between temporal CCA v4 and spatial CCA v2,
all eight region pairs co-plotted. Each datum is one (animal, pair,
trial id) cell intersected between the two pipelines. Two metrics
(CC1, IFI) × three epochs (early, pre, post) = six panels arranged in a
2 × 3 grid. Within each panel, dots are colour-coded by region pair;
the dashed line is the identity *y = x*. Spearman ρ in each panel is
computed across all paired points in that panel (i.e. all pairs, all
animals, all trials in the relevant epoch); panel-level *n* is given
in the title. Pairs that did not produce a per-trial vector in either
pipeline (typically due to insufficient units in one region for that
animal) are absent from the panel. Note that CC1 is bounded to [0, 1]
and IFI to [-1, 1] by construction (after the 2026-05-08 IFI fix).

## Figure 2 — Per-pair paired scatter (real values)

**File**: `HC_V1_figures/CCA_Compare/TemporalVsSpatial_per_pair.svg`.

Per-pair detail of Figure 1. Eight rows (one per region pair) × six
columns (CC1 in early/pre/post, then IFI in early/pre/post). Each panel
plots the trial-by-trial paired values for that pair × epoch × metric,
split by group: blue = learner, red = non-learner. The dashed identity
line and per-panel Spearman ρ are computed over all (animal, trial) data
points in the panel, ignoring group. Group labels appear via the colour
key in the leftmost panels. *n* per panel in the panel title is the
number of finite paired points contributing to ρ. Sparse panels (e.g.
pairs containing CA3 or DG, where some animals have <5 units in those
regions) reflect dropout in one or both pipelines, not zero coupling.

## Figure 3 — Combined paired scatter (shuffle-corrected)

**File**: `HC_V1_figures/CCA_Compare/TemporalVsSpatial_excess_combined.svg`.

Same layout and conventions as Figure 1, but axes plot **excess** =
real - shuffle for both temporal and spatial. The shuffle is the
per-trial circular row permutation (temporal CC1, IFI) or the per-trial
spatial-bin shift (spatial CC1, precession_idx) described in METHODS
§3.5 and §4.4. Subtracting the shuffle removes upward bias from
small-sample CCA inflation, so a non-zero excess on either axis
reflects coupling beyond chance. The methodologically interesting
pattern is **whether the temporal-vs-spatial relationship visible in
Figure 1 survives shuffle correction**: if it does, the asymmetry is
biological; if it collapses to a tight diagonal cluster near the
origin, the asymmetry was driven by per-method bias. Dotted axes
through the origin help locate quadrant identity.

## Figure 4 — Per-pair paired scatter (shuffle-corrected)

**File**: `HC_V1_figures/CCA_Compare/TemporalVsSpatial_excess_per_pair.svg`.

Per-pair detail of Figure 3, conventions matching Figure 2 (eight rows
× six columns, learner/non-learner split). Particularly useful for
spotting pair-specific dissociations: e.g. an RSC-containing pair with
spatial excess CC1 substantially above zero but temporal excess CC1
near zero would indicate genuine spatial co-tuning without
moment-to-moment co-fluctuation, distinct from the case where both
methods show shuffle-level coupling.

## Figure 5 — Marginal distributions per epoch

**Files**: `HC_V1_figures/CCA_Compare/TemporalVsSpatial_marginals_early.svg`,
`_marginals_pre.svg`, `_marginals_post.svg` (one per epoch).

Per-pair marginal distributions of each metric, split by group. Each
figure has eight rows (region pairs) × four columns (temporal CC1,
spatial CC1, temporal IFI, spatial IFI). Within a panel, points are
jittered along x at column 1 (learner) and column 2 (non-learner);
the horizontal bar marks the group median; column counts in the
x-tick labels report the number of finite trial-level data points per
group. Two p-values are reported in each panel: `p_KS` from a two-
sample Kolmogorov-Smirnov test against the null of identical
distributions, and `p_RS` from a Wilcoxon rank-sum test against the
null of equal medians. Significance asterisks (`*`, `**`, `***` for
p < 0.05, 0.01, 0.001) reflect the more conservative of the two; no
multiple-comparison correction is applied at the figure level. These
panels test for **group differences within a method** — a question
distinct from "do the methods agree", which is what Figures 1–4
address.

## Figure 6 — Asymmetry summary across pairs

**File**: `HC_V1_figures/CCA_Compare/TemporalVsSpatial_asymmetry_summary.svg`.

Four-panel summary quantifying the divergence between the two methods,
one bar per epoch, grouped by region pair (eight pair groups on the
x-axis).

(A) Mean trial-level *temporal CC1 − spatial CC1*. Bars above zero
indicate temporal-biased coupling; below zero, spatial-biased.
Computed across all (animal, trial) finite paired points per pair ×
epoch.

(B) Mean trial-level *temporal excess CC1 − spatial excess CC1*
(i.e. shuffle-corrected version of (A)). Direct test of whether the
asymmetry in (A) survives bias correction. If panels (A) and (B)
disagree in sign or magnitude on some pair, that pair's asymmetry is
methodological in origin.

(C) Sign-agreement rate of trial-level IFI between the two methods.
Per (animal, trial), the IFI signs from the two methods are compared;
the bar reports the fraction of trials where signs agree. Dashed line
at 0.5 marks chance. Trials with IFI exactly zero in either method
are excluded.

(D) Spearman ρ of trial-level |IFI| between the two methods. Tests
whether the methods agree on which trials show stronger lag
asymmetry, regardless of direction. ρ > 0 indicates magnitude
agreement; ρ near 0 indicates the two methods select independent
trials as "high-asymmetry".

The numerical content of each panel is also saved to
`HC_V1_data/compare_asymmetry.mat`.
