# HANDOFF — how the temporal/lag analyses are organised

**For a cold agent or a future session.** `PROJECT_LOG.md` tells you *what happened* and
`STATE.md` *what is believed*; this file tells you **where things live, which script made
them, and which numbers you are allowed to trust**. Read this before running anything.

Last updated 2026-08-07, after the seven 2026-07-28 meeting items plus the four follow-ups
Theo asked for (per-CC IFI signs, FF/FB label-and-track, cosine-across-lags, integration
window vs IFI by epoch).

---

## 0. The one-minute orientation

- Everything below is the **temporal arm** (running-state, 10 ms bins). The landmark and
  spatial arms are older and separate — see `STATE.md` §2.
- **Two families of analysis, and confusing them is the main hazard:**
  - **REFIT-PER-LAG** — CCA refit at every lag. `run_lag_curves`, `run_lag_subspaces`,
    `run_lag_cosine`. A canonical dimension here is a **rank**, not a tracked component.
  - **FROZEN AXES** — one CCA fit, then project subsets through identical weights.
    `run_fixed_subspace`, `run_cc_label_track`. A dimension here **is** a stable
    component, which is why every epoch contrast uses this family.
- **Why it matters:** measured on this data, the best-matching dimension at another lag
  is a *different* dimension **82 %** of the time. Any analysis that matches dimensions
  across two fits by bare index is comparing different things. This has caused two real
  bugs (`d38a833`, `1bae90e`).

---

## 1. Pipeline map — script → output → figure

Every driver takes `--bin-ms 10 --smooth-ms 2.5` and writes both FS conditions
(`_fsincl` suffix = fast-spiking units included). **FS-excluded and FS-included are
co-primary — report both.**

| # | Question | Driver | Analysis | Figures |
|---|---|---|---|---|
| 1 | Do different CCs have different IFI (different signs)? | `run_lag_curves.py` → `lag_curves_bin10*.csv` | `analyze_cc_ifi_signs.py` → `cc_ifi_windows_*`, `cc_ifi_signs_*`, `cc_ifi_direction*`, `cc_ifi_mixing_*`, `cc_ifi_signs_tables.md` | `figs_cc_ifi_signs.py` → `HCV1_cc_ifi_windows_*`, `HCV1_cc_ifi_signmap_*` |
| 2 | Separate subspaces into FF/FB, track with learning | `run_cc_label_track.py` → `cc_label_track_bin10*.csv` | `analyze_cc_label_track.py` → `cc_label_track_epoch_*`, `_stats_*`, `_tables.md` | `figs_cc_label_track.py` → `HCV1_cc_label_track_*`, `HCV1_cc_label_persistence_*` |
| 3 | How stable are the CCs across lags? | `run_lag_subspaces.py` → `lag_subspaces_bin10*.csv` | `analyze_lag_subspaces.py` → `lag_subspaces_stability_*`, `_fffb_*`, `_gate_sensitivity_*`, `_evolution_*` | `figs_lag_subspaces.py`, `figs_lag_subspaces_epochs.py` |
| — | Cosine of each CC to itself across lags | `run_lag_cosine.py` → `lag_cosine_bin10*.csv` | (inline) | `figs_lag_cosine.py` → `HCV1_lag_cosine_*`, `HCV1_lag_cosine_swap_*` |
| 4 | Integration window vs IFI, naive vs exp | *(reuses item 2's CSV)* | `analyze_ifi_windows_epochs.py` → `ifi_windows_epochs_*` | `figs_ifi_windows_epochs.py` |
| 5/6/7 | One CC lagged across time; windows; fixed subspace naive vs exp | `run_fixed_subspace.py` → `fixed_subspace_bin10*.csv` | `analyze_fixed_subspace.py` → `fixed_subspace_epoch_*`, `_stats_*` | `figs_fixed_subspace.py`, `figs_integration_windows.py` |

**Written write-up of items 1–7:** `results/MEETING_2026-08-07.md` (methods + numbers +
caveats, one section per item). ⚠ It predates the four follow-ups and the significance
rework below; treat its per-CC significance statements as superseded.

**Figures** are written to `~/Documents/ResearchVault/attachments/` **and mirrored to
`cca/figures/`** by `figstyle.save` (one render, copied — byte-identical). Both are
gitignored; regenerate from the scripts.

---

## 2. Shared preprocessing (identical everywhere unless stated)

1. 1 ms spike trains → Gaussian **σ = 2.5 ms** → **10 ms** bins.
2. Engaged running bins only: in-trial and **velocity ≥ 2 cm/s**.
3. An area needs **≥ 5 units**. Eight pairs (see any driver's `PAIRS`).
4. Third-area control: all *other* recorded areas regressed out of both X and Y
   (`partial.partial_out_cv`, train-only coefficients where the readout is cross-validated).
5. PCA to **k = 30** per area.
6. CV: **5 folds split by WHOLE TRIALS**.
7. **Sign convention: positive lag ⇒ the FIRST-named area leads.** `lag_slice` pairs
   `X[t]` with `Y[t+lag]`.
8. Inferential unit: **animals-as-n**, paired *t*. 8 pairs = a per-pair family, no
   cross-pair correction (project policy, `STATE.md` §3.0).

### ⚠ The sample cap differs by driver and it matters

- `run_lag_curves`, `run_lag_subspaces`, `run_lag_cosine` cap at **12 000 bins** by taking
  consecutive within-trial blocks. Because a trial is ~2 900 bins (~30 s of running),
  that stops inside the **first ~20 trials** — for a late learner that window is entirely
  pre-learning. Fine for questions about lag geometry; **wrong for anything comparing
  epochs**.
- `run_cc_label_track` has **no cap** — every running bin of every trial. It fits CCA
  once, so the cap bought nothing and actively hurt: the old cap kept only the first
  0.4–1.8 s of each trial, and **running speed rises +12.0 cm/s naive→expert on those
  trial-onset bins versus +6.6 cm/s over whole trials**, doubling a confound on the very
  contrast being made.

---

## 3. Column dictionary for the files you will actually open

**`lag_curves_bin10*.csv`** — refit-per-lag, held-out, segment-aware, one row per
(animal, pair, lag, dim). 76 k rows.
`cc` held-out canonical correlation at that lag · `sig` per-dim significance **at lag 0**
· `sig_uncorr`, `p_perdim` the uncorrected mask and empirical p · `n_sig`.

**`cc_label_track_bin10*.csv`** — frozen axes, one row per (animal, pair, dim, epoch, lag).
`r` correlation of the frozen variates (**IN-SAMPLE** — a contrast statistic, never a
coupling strength) · `label` FF/FB from the whole-session fit · `label_loo`, `ifi_loo` the
same recomputed **excluding the epoch being scored** (use this one for persistence) ·
`sig`, `p_perdim`, `cc_heldout` from `frozen_perm_null` · `n_bins`, `n_fit_trials`.

**`lag_cosine_bin10*.csv`** — one row per (animal, pair, dim, lag).
`cos_same_*`, `cos_best_*`, `best_dim_*` full-data cosine vs the lag-0 fit ·
**`cos_split_*`, `split_best_*` the matched split-half version — use these**, since the
full-data ones start at 1.0 by construction at lag 0.

**`lag_subspaces_bin10*.csv`** — refit-per-lag subspaces. `angle_x/y` (3-dim) and
`angle_x_cc1/angle_y_cc1` (CC1) principal angles vs lag 0, with **per-area** floors
`floor_x/floor_y` and `floor_x_cc1/floor_y_cc1`. **Use the CC1 columns**; the 3-dim floor
is ~78°, i.e. unmeasurable.

---

## 4. Statistical machinery specific to this project

| Function | Where | What it is for |
|---|---|---|
| `lagged.information_flow_index` | `src/tom_cca/lagged.py` | IFI = (X-leads − Y-leads)/(sum), CC clipped at 0, lag 0 excluded |
| `perdim_ifi.ifi_windows_by_dim` | `perdim_ifi.py` | IFI per dim over integration windows; `degenerate` flags "0 because both sides clipped", ≠ "balanced" |
| `lagged.perdim_significance` | `lagged.py` | Per-dim **held-out** circular-shift null for the refit-per-lag arm |
| `fixed_subspace.frozen_perm_null` | `fixed_subspace.py` | Permutation test on a **frozen** fit's own canonical correlations, rank for rank |
| `fixed_subspace.trial_lag_moments` / `curve_from_moments` | `fixed_subspace.py` | Per-(trial,lag) sufficient statistics; any subset's curve by masking. Makes uncapped runs tractable |
| `lag_subspace.split_half_floor` | `lag_subspace.py` | Per-area noise floor for subspace angles |

**Two arithmetic traps that will bite you again if you forget them:**

1. **Permutation p-floor.** An empirical p cannot go below `1/(n_shuffles+1)`, while BH
   needs `alpha/d`. With d = 30 that needs **> 599 shuffles before anything can pass** —
   otherwise the corrected mask is empty for arithmetic reasons and reads as a scientific
   null. Hence `fdr_dims` restricts the BH family to the leading 10 dims and drivers
   default to 200 shuffles.
2. **Mask monotonicity depends on the null.** Under a *dominant-dim* null (one scalar
   threshold) the mask must be monotone in cc. Under a *per-dim* null each dimension has
   its own bar and shuffled correlations fall with rank, so a smaller high-rank CC can
   legitimately pass where a larger low-rank one does not. A guard that checks
   monotonicity will false-alarm on correct per-dim data.

---

## 5. What is believed, and how strongly

**Positive results (2):**
- **FF/FB labels persist across epochs**: 62 % (FS-excl, p = 7×10⁻⁴) / 66 % (FS-incl,
  p = 5×10⁻⁶) against a *genuine* 50 % (leave-epoch-out label), animals-as-n, 10/12 then
  12/12 animals above chance. **Not** driven by theta-ringing pairs (intra − other is
  +0.05 then −0.03, both n.s., sign flips).
- **CA1-DG rotates with lag** (item 3, `11c0278`): survives every estimability gate
  including none, both FS, Δ 20–35°, p(Bonf) to 0.0003.

**Nulls (everything else):** per-CC IFI sign is at chance (item 1); FF/FB not separable as
subspaces, 0/8 (item 2 session-level); FF/FB do not diverge with learning (item 2
interaction, 0 survive FDR); no epoch effect through a frozen subspace, 0/32 (items 5/7);
integration window vs IFI unchanged with learning (item 4).

**Measured methodological facts you should not re-derive:**
- Split-half |cos| of a canonical vector at a fixed lag: **CC1 0.59, CC2 0.39, CC3 0.26,
  CC4 0.22** against a **0.146** random-vector baseline in 30-D. Only CC1 has real
  identity.
- Rank-swap rate across lags: **82 %**, both FS.
- CC1 loses about half its identity by ±250 ms (0.556 → 0.300, p = 1×10⁻⁴).

---

## 6. Open threads, in the order I would pick them up

1. **The running-speed confound is unresolved.** Speed rises **+6.6 cm/s naive→expert**
   over whole trials (11/12 animals, p = 0.005 on the kept bins). Nothing in the pipeline
   controls for velocity, and every epoch contrast here is a *timing* measure while speed
   sets theta frequency. Uncapping halved it; it did not remove it. **Ask Theo whether
   this is a known feature of the cohort before building a speed-matched analysis.**
2. **`figs_report.fig_rotation_floor` still plots the unmeasurable 3-dim angle.** The
   rotation-null verdict is safe (§G of `bin10_tables.md` tests it at CC1) but the figure
   should switch.
3. **Two direction discrepancies vs `STATE.md` §3.0 finding 3**, unreconciled: CA1-SUB
   reads CA1→SUB here (5/7 animals) vs SUB→CA1 there; and CA1→RSC is FS-fragile under the
   subspace-refit method (p = 0.015 FS-incl, 0.764 FS-excl) where §3.0 has it at
   p = 3.9×10⁻⁴. Different methods, so not directly comparable — but not ignorable.
4. **Gini partner-invariance re-test** (`STATE.md` §3.0 finding 2 is held open on it). The
   `trajectory_gini3_*` outputs exist with `gini_*_conn` / `gini_*_sig` populated; nobody
   has re-run the broadening result on them.
5. **`MEETING_2026-08-07.md` needs a refresh** — it predates the per-dim significance
   rework and the four follow-ups.

---

## 7. How to work here without repeating our mistakes

- **Never gate one fit's dimensions with another fit's statistics.** Compute significance
  in the same driver, on the same fit, from the same scores. Two separate bugs came from
  exactly this.
- **Every "at floor" verdict needs a floor you have actually measured**, and the floor
  must be the *same estimator* as the thing it gates (half-data floor vs full-data
  comparison is a mismatch — it made item 3's 3-dim test unmeasurable and inverted its
  conclusion once).
- **Check a chance model before calling a rate a finding.** Sign-mixing "92 %" was *below*
  its own chance level; label persistence "70 % vs 50 %, p = 1e-11" had a construction
  floor of ~60 %.
- **A result that looks too clean is a bug until checked** — `cc_heldout` rising with rank
  in 38/38 cells was the tell that dimensions were mis-attributed.
- **Smoke-test a driver on a few shuffles before launching the long run.** A variable
  shadowing bug (`out = ~mask` clobbering the CSV path) surfaced only at the final flush,
  after 191 s of real work.
- Run the full suite before committing: `cd cca && PYTHONPATH=src python -m pytest -q`
  (**373 tests**, ~50 s). Lint with `uvx ruff check --select F`.
