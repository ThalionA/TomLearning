"""Trial 1 vs trial 2 — analysis of the frozen-projection per-trial readout.

Reads results/trial12_curves_bin10{,_fsincl}.csv + trial12_trials_bin10{,_fsincl}.csv
(from run_trial12.py). Derives per-dim IFIs at every integration window through
`perdim_ifi.windows_table` (degenerate flags intact), collapses dims per animal
(`cc_aggregate.per_animal_mean`: CC1 and all-significant-CC tables), forms per-animal
trial-1-minus-2 deltas, and tests them animals-as-n (`one_sample_by_pair`: paired t
primary, Wilcoxon cross-check, BH across pairs as sensitivity) at the pre-registered
headline window (±LABEL_W bins = ±50 ms). The full window sweep is DESCRIPTIVE
(cumulative windows are nested — GOTCHAS.md); the disjoint 60-250 ms band is the one
independent wide-lag column.

Control: the 1-vs-2 delta is calibrated against ordinary adjacent-trial variability
on the common-matched arm — z = (m1 - m2) / (sqrt(2) * SD(m over ordinals 3..10))
per (animal, pair), plus the assumption-light |d12| - mean|adjacent delta|. Signed
adjacent deltas centre on 0 under any hypothesis, so a mean-excess "control" would
be vacuous: SCALE is the question. A trend across ordinals 3..10 inflates the SD,
making the z conservative.

Writes: trial12_delta_tests_bin10{,_fsincl}.csv    headline tests, all metrics/arms
        trial12_sweep_bin10{,_fsincl}.csv          per-pair delta vs window
        trial12_completeness_bin10{,_fsincl}.csv
        trial12_tables.md                          both FS, verdict strings generated

Usage: PYTHONPATH=src python scripts/analyze_trial12.py
"""

from __future__ import annotations

import zlib

import numpy as np
import pandas as pd

from _common import FS_CONDITIONS, RES, TEMPORAL, config
from tom_cca import cc_aggregate, lagged, paired_stats, perdim_ifi

PAIRS = list(config.PAIR_NAMES)
BIN_MS = TEMPORAL.bin_ms
HEADLINE_W = TEMPORAL.label_w_bins          # 5 bins = ±50 ms
WIDE_MIN = HEADLINE_W + 1                   # disjoint band: 60..250 ms
CTRL_ORDS = list(range(3, 11))              # the out-of-fit adjacent-trial band
MIN_CTRL = 6                                # finite control ordinals required
MIN_N = 4                                   # animals per pair for a test to print
UNDERPOWERED_N = 6                          # n <= this is labelled underpowered


def per_trial_windows(curves: pd.DataFrame) -> pd.DataFrame:
    """One row per (animal, pair, ordinal, matched, dim, window): IFI + degenerate."""
    return perdim_ifi.windows_table(
        curves, value="r",
        group_keys=("animal", "learner", "pair", "ordinal", "matched"),
        bin_ms=BIN_MS)


def wide_band(curves: pd.DataFrame) -> pd.DataFrame:
    """Per (animal, pair, ordinal, matched, dim): IFI on the DISJOINT |lag| >= WIDE_MIN
    band — the one wide-lag statement independent of the (nested) headline window.

    ``cc_peak`` is the per-dim max over ALL lags, matching
    :func:`perdim_ifi.windows_table` exactly. It must: both frames feed the same
    ``per_animal_mean(weight="cc_peak")``, and the weight means "does this dim carry
    any coupling at all" — a band-local max is a different quantity under the same
    name, and using one silently changed both wide-band results (2026-08-20 audit).
    """
    df = curves.copy()
    df["r"] = pd.to_numeric(df["r"], errors="coerce")
    keys = ["animal", "learner", "pair", "ordinal", "matched", "dim"]
    peak_all = df.groupby(keys)["r"].max()               # ALL lags, as windows_table
    band = df[np.abs(df["lag_bins"]) >= WIDE_MIN]
    rows = []
    for key, g in band.groupby(keys, sort=True):
        lg, c = g["lag_bins"].to_numpy(float), g["r"].to_numpy(float)
        m = np.isfinite(c)
        if not (np.any(lg[m] > 0) and np.any(lg[m] < 0)):
            continue
        pm, nm = lagged.ifi_sides(lg[m], c[m])
        rows.append({**dict(zip(keys, key)),
                     "ifi": lagged.information_flow_index(lg[m], c[m]),
                     "degenerate": int((pm + nm) <= 0),
                     "cc_peak": float(peak_all.get(key, np.nan)),
                     "sig": int(g["sig"].max())})
    return pd.DataFrame(rows)


def dim_matched(tab: pd.DataFrame, value: str = "ifi", ordinals=(1, 2)) -> pd.DataFrame:
    """Keep only dims that are usable in BOTH ordinals of the contrast.

    ``per_animal_mean(drop_degenerate=True)`` drops degenerate dims PER ORDINAL, so a
    delta would subtract a mean over trial 2's dim-set from a mean over trial 1's —
    different supports, which is exactly the like-for-like failure the project's
    non-negotiables name. Measured on this data (2026-08-20 audit): 85/213 headline
    (animal, pair, arm) cells had mismatched dim-sets, and matching them FLIPS the
    sign of the CA1-V1 all-CC IFI delta in all three arms. Apply this BEFORE
    ``per_animal_mean(..., drop_degenerate=False)``.
    """
    sub = tab[tab["ordinal"].isin(ordinals)].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    ok = np.isfinite(sub[value])
    if "degenerate" in sub.columns:
        ok &= pd.to_numeric(sub["degenerate"], errors="coerce") == 0
    if "sig" in sub.columns:
        ok &= pd.to_numeric(sub["sig"], errors="coerce") == 1
    sub = sub[ok]
    keys = [c for c in ("animal", "pair", "matched", "window_bins") if c in sub.columns]
    n_ord = sub.groupby([*keys, "dim"])["ordinal"].nunique()
    keep = n_ord[n_ord == len(ordinals)].index
    return sub.set_index([*keys, "dim"]).loc[
        sub.set_index([*keys, "dim"]).index.isin(keep)].reset_index()


TAIL_DIMS = (6, 10)          # dims whose held-out r must be near zero
TAIL_MAX = 0.5               # mean |r0| over TAIL_DIMS above this = degenerate trial


def degenerate_trials(trials: pd.DataFrame) -> pd.DataFrame:
    """(animal, pair, ordinal, matched) cells whose frozen projection has collapsed.

    A genuine held-out canonical correlation DECAYS with rank: CC1 carries the
    coupling and the tail sits near zero (project-wide, CC1 is 0.05-0.4). A trial
    reading |r| ~ 1 at rank 8 is a rank/collinearity degeneracy of the projection on
    that trial's bins, not coupling. Found 2026-08-20 in animal 70, ordinals 1 and 9,
    pairs CA1-DG / CA1-SUB / CA1-V1 (CC1 -0.941 with CC8 -0.996 on the same trial);
    those cells alone inflated the delta SD ~10x and made the strength null formally
    uninformative for exactly those three pairs.
    """
    t = trials.copy()
    t["r0"] = pd.to_numeric(t["r0"], errors="coerce")
    tail = t[(t["dim"] >= TAIL_DIMS[0]) & (t["dim"] <= TAIL_DIMS[1])]
    g = (tail.groupby(["animal", "pair", "ordinal", "matched"])["r0"]
         .apply(lambda s: s.abs().mean()).rename("tail_abs_r").reset_index())
    g["degenerate_trial"] = (g["tail_abs_r"] > TAIL_MAX).astype(int)
    return g


def drop_degenerate_trials(df: pd.DataFrame, flags: pd.DataFrame) -> pd.DataFrame:
    """Remove the flagged (animal, pair, ordinal, matched) cells from any per-dim or
    per-animal frame that carries those keys."""
    keys = [c for c in ("animal", "pair", "ordinal", "matched") if c in df.columns]
    if not flags.empty and set(keys) <= set(flags.columns):
        bad = flags[flags["degenerate_trial"] == 1][keys]
        if not bad.empty:
            merged = df.merge(bad.assign(_bad=1), on=keys, how="left")
            return merged[merged["_bad"].isna()].drop(columns="_bad")
    return df


def cc1_per_animal(win: pd.DataFrame, value: str = "ifi") -> pd.DataFrame:
    """CC1 rows as a per-animal frame (dim 1 is one value per animal already)."""
    sub = win[(win["dim"] == 1) & (win["degenerate"] == 0)].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    return sub[np.isfinite(sub[value])]


def delta12(pa: pd.DataFrame, value: str) -> pd.DataFrame:
    """Per-animal trial-1-minus-2 delta (and the |·| delta) — animals present in
    BOTH ordinals only; the completeness table counts the drops."""
    sub = pa[pa["ordinal"].isin([1, 2])].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    idx = [c for c in ("pair", "animal", "matched", "window_bins") if c in sub.columns]
    wide = sub.pivot_table(index=idx, columns="ordinal", values=value, aggfunc="first")
    if 1 not in wide.columns or 2 not in wide.columns:
        return pd.DataFrame(columns=idx + ["delta", "delta_abs"])
    wide = wide.dropna(subset=[1, 2])
    out = wide.reset_index()
    out["delta"] = out[1] - out[2]
    out["delta_abs"] = out[1].abs() - out[2].abs()
    return out.drop(columns=[1, 2])


def control_stats(pa: pd.DataFrame, value: str) -> pd.DataFrame:
    """Per (animal, pair) on the common-matched arm: the standardised 1-vs-2 delta
    z = d12 / (sqrt(2)·SD over ordinals 3..10), and |d12| - mean|adjacent delta|."""
    sub = pa[pa["matched"] == 2].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    rows = []
    for (pair, animal), g in sub.groupby(["pair", "animal"]):
        m = {int(o): float(v) for o, v in zip(g["ordinal"], g[value])
             if np.isfinite(v)}
        ctrl = [m[o] for o in CTRL_ORDS if o in m]
        if 1 not in m or 2 not in m or len(ctrl) < MIN_CTRL:
            continue
        sd = float(np.std(ctrl, ddof=1))
        d12 = m[1] - m[2]
        adj = [m[o + 1] - m[o] for o in CTRL_ORDS[:-1] if o in m and o + 1 in m]
        rows.append({"pair": pair, "animal": animal, "n_ctrl": len(ctrl),
                     "z": d12 / (np.sqrt(2.0) * sd) if sd > 0 else np.nan,
                     "abs_excess": (abs(d12) - float(np.mean(np.abs(adj)))
                                    if adj else np.nan)})
    return pd.DataFrame(rows)


def _test(per_animal: pd.DataFrame, value: str, metric: str,
          matched) -> pd.DataFrame:
    out = cc_aggregate.one_sample_by_pair(per_animal, value, PAIRS, min_n=MIN_N,
                                          wilcoxon=True)
    out.insert(0, "metric", metric)
    out.insert(1, "matched", matched)
    return out


def headline_tests(win, wideb, allsig, allsig_wide, trials, trials_allsig
                   ) -> pd.DataFrame:
    """Every (metric, matched-arm) test at the headline window, one frame."""
    tests = []
    for matched in (0, 1, 2):
        # -- strength through the frozen axes, held out per trial
        t1 = trials[(trials["dim"] == 1) & (trials["matched"] == matched)]
        tests.append(_test(delta12(t1, "r0"), "delta", f"r0_cc1", matched))
        ts = trials_allsig[trials_allsig["matched"] == matched]
        for col, lab in (("mean", "r0_allsig"), ("wmean", "r0_allsig_w")):
            if col in ts.columns:
                tests.append(_test(delta12(ts, col), "delta", lab, matched))
        # -- IFI at the headline window
        w1 = cc1_per_animal(win[(win["window_bins"] == HEADLINE_W)
                                & (win["matched"] == matched)])
        d1 = delta12(w1, "ifi")
        tests.append(_test(d1, "delta", "ifi_cc1", matched))
        tests.append(_test(d1, "delta_abs", "absifi_cc1", matched))
        aw = allsig[(allsig["window_bins"] == HEADLINE_W)
                    & (allsig["matched"] == matched)]
        for col, lab in (("mean", "ifi_allsig"), ("wmean", "ifi_allsig_w")):
            da = delta12(aw, col)
            tests.append(_test(da, "delta", lab, matched))
            if col == "mean":
                tests.append(_test(da, "delta_abs", "absifi_allsig", matched))
        # -- disjoint wide band
        wb1 = cc1_per_animal(wideb[wideb["matched"] == matched])
        tests.append(_test(delta12(wb1, "ifi"), "delta", "ifi_wide_cc1", matched))
        ab = allsig_wide[allsig_wide["matched"] == matched]
        if "wmean" in ab.columns:
            tests.append(_test(delta12(ab, "wmean"), "delta", "ifi_wide_allsig_w",
                               matched))
    # -- the confounds, on EVERY arm. Matching equalises bin count exactly, but it
    # does NOT equalise speed or fragmentation: it keeps the LEADING bins, so it
    # compares roughly the first 80 % of trial 1 against all of trial 2 — different
    # phases of a trial, and the phase it over-weights is the trial-onset phase that
    # HANDOFF.md §6 documents as carrying a DOUBLED speed confound. Testing these on
    # the raw arm only (as this script did until the 2026-08-20 audit) hid that.
    for matched in (0, 1, 2):
        t1a = trials[(trials["dim"] == 1) & (trials["matched"] == matched)]
        for col in ("n_bins", "vel_mean", "n_gaps"):
            tests.append(_test(delta12(t1a, col), "delta", f"d_{col}", matched))
    return pd.concat([t for t in tests if not t.empty], ignore_index=True)


def direction_exists(win: pd.DataFrame, allsig: pd.DataFrame) -> pd.DataFrame:
    """Per ordinal 1 and 2 (pairwise-matched arm): is there a direction AT ALL?
    One-sample t of the per-animal headline IFI against 0."""
    rows = []
    for o in (1, 2):
        w1 = cc1_per_animal(win[(win["window_bins"] == HEADLINE_W)
                                & (win["matched"] == 1) & (win["ordinal"] == o)])
        t = cc_aggregate.one_sample_by_pair(w1, "ifi", PAIRS, min_n=MIN_N,
                                            wilcoxon=True)
        t.insert(0, "metric", "ifi_cc1"); t.insert(1, "ordinal", o)
        rows.append(t)
        aw = allsig[(allsig["window_bins"] == HEADLINE_W) & (allsig["matched"] == 1)
                    & (allsig["ordinal"] == o)]
        if "wmean" in aw.columns:
            t = cc_aggregate.one_sample_by_pair(aw, "wmean", PAIRS, min_n=MIN_N,
                                                wilcoxon=True)
            t.insert(0, "metric", "ifi_allsig_w"); t.insert(1, "ordinal", o)
            rows.append(t)
    return pd.concat([r for r in rows if not r.empty], ignore_index=True)


def completeness(trials: pd.DataFrame) -> pd.DataFrame:
    """Per (pair, matched): animals with a finite CC1 r0 at each ordinal 1/2, and
    surviving 1-vs-2 deltas — the missingness is non-random (short first trials)."""
    t1 = trials[trials["dim"] == 1].copy()
    t1["r0"] = pd.to_numeric(t1["r0"], errors="coerce")
    rows = []
    for (pair, matched), g in t1.groupby(["pair", "matched"]):
        fin = {o: g[(g["ordinal"] == o) & np.isfinite(g["r0"])]["animal"].nunique()
               for o in (1, 2)}
        n_animals = g["animal"].nunique()
        d = delta12(g, "r0")
        rows.append({"pair": pair, "matched": matched, "n_animals": n_animals,
                     "finite_t1": fin[1], "finite_t2": fin[2],
                     "n_deltas": int(len(d))})
    return pd.DataFrame(rows)


def sweep(win: pd.DataFrame, allsig: pd.DataFrame) -> pd.DataFrame:
    """Descriptive: per (pair, matched, window) mean ± sem of the 1-vs-2 delta."""
    rows = []
    for label, pa, col in (("ifi_cc1", cc1_per_animal(win), "ifi"),
                           ("ifi_allsig_w", allsig, "wmean")):
        if col not in pa.columns:
            continue
        d = delta12(pa, col)
        if d.empty:
            continue
        g = d.groupby(["pair", "matched", "window_bins"])["delta"]
        agg = g.agg(["mean", "sem", "count"]).reset_index()
        agg.insert(0, "metric", label)
        rows.append(agg)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def sweep_test(win: pd.DataFrame, allsig: pd.DataFrame, n_perm: int = 10_000,
               seed: int = 0) -> pd.DataFrame:
    """TEST the window sweep, honouring the nesting (the user's 4th ask).

    Integration windows are CUMULATIVE (|lag| <= w), so the 25 per-window tests are
    heavily dependent and a plain BH over them would be anti-conservative — which is
    why the sweep shipped descriptive-only on 2026-08-20. The defensible statistic is
    **max |t| over the windows** with a sign-flip permutation null: flipping the sign
    of each animal's whole delta-vs-window PROFILE preserves the nesting exactly and
    needs no independence assumption. One p per (metric, arm, pair); BH across the 8
    pairs within each (metric, arm), as everywhere else.

    Columns: ``metric, matched, pair, n, best_window_bins, best_window_ms, mean_at_best,
    max_abs_t, p_perm, bh_pass``.
    """
    rows = []
    for label, pa, col in (("ifi_cc1", cc1_per_animal(win), "ifi"),
                           ("ifi_allsig_w", allsig, "wmean")):
        if col not in pa.columns:
            continue
        d = delta12(pa, col)
        if d.empty:
            continue
        for (pair, arm), g in d.groupby(["pair", "matched"]):
            piv = g.pivot_table(index="animal", columns="window_bins", values="delta")
            piv = piv.dropna(axis=0, how="any")
            if piv.shape[0] < MIN_N or piv.shape[1] < 2:
                continue
            M = piv.to_numpy(float)                      # (n_animals, n_windows)
            n = M.shape[0]

            def _absT(A):
                sd = A.std(axis=0, ddof=1)
                sd = np.where(sd > 0, sd, np.nan)
                return np.abs(A.mean(axis=0) / (sd / np.sqrt(n)))

            t_obs = _absT(M)
            if not np.any(np.isfinite(t_obs)):
                continue
            j = int(np.nanargmax(t_obs))
            obs = float(t_obs[j])
            # Seed PER GROUP, not once for the whole sweep: a single stream is
            # consumed in group order, so dropping any upstream row (e.g. a
            # degenerate trial) silently changes every later group's p. Two
            # borderline results moved across 0.05 that way before this fix.
            # (zlib.crc32, not hash(): Python randomises string hashing per process,
            # so hash() would make the p-values irreproducible across runs.)
            rng = np.random.default_rng(
                zlib.crc32(f"{seed}|{label}|{pair}|{arm}".encode()))
            null = np.empty(n_perm)
            for i in range(n_perm):
                s = rng.choice([-1.0, 1.0], size=(n, 1))
                null[i] = np.nanmax(_absT(M * s))
            p = float((1 + np.sum(null >= obs)) / (1 + n_perm))
            w = int(piv.columns[j])
            rows.append({"metric": label, "matched": int(arm), "pair": pair, "n": n,
                         "best_window_bins": w, "best_window_ms": w * BIN_MS,
                         "mean_at_best": float(M[:, j].mean()),
                         "max_abs_t": obs, "p_perm": p})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    keep = []
    for _, g in out.groupby(["metric", "matched"]):
        g = g.copy()
        g["bh_pass"] = paired_stats.fdr_bh(g["p_perm"].to_numpy(float))
        keep.append(g)
    return pd.concat(keep, ignore_index=True)


def fs_agreement(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    """Do the two co-primary FS conditions agree? Per (metric, arm): correlation of
    the per-pair delta means and sign-agreement count. Nothing else checks this, and
    the two conditions share animals and most units, so agreement is weak evidence
    while DISagreement is strong evidence of fragility."""
    m = a.merge(b, on=["metric", "matched", "pair"], suffixes=("_x", "_i"))
    rows = []
    for (metric, arm), g in m.groupby(["metric", "matched"]):
        x, y = g["mean_x"].to_numpy(float), g["mean_i"].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 3:
            continue
        rows.append({"metric": metric, "matched": int(arm), "n_pairs": int(ok.sum()),
                     "r": float(np.corrcoef(x[ok], y[ok])[0, 1]),
                     "same_sign": int(np.sum(np.sign(x[ok]) == np.sign(y[ok])))})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, digits: int = 4) -> str:
    return df.round(digits).to_markdown(index=False)


def main():
    md = ["# Trial 1 vs trial 2 — frozen-subspace per-trial readout", "",
          "Generated by `analyze_trial12.py`; every number below is derived from "
          "`trial12_*_bin10*.csv` at run time (no hand-written verdicts).", "",
          "> Subspace fit on running-trial ordinals 11+ (`early_trials.reference_fit`"
          "), ordinals 1..10 projected leak-free. Positive IFI ⇒ the first-named "
          "area leads. `arm`: 0 = raw bins, 1 = ordinals 1/2 cut to min(n1,n2) "
          "leading bins, 2 = all ordinals cut to a common min (the control arm). "
          f"Headline window ±{HEADLINE_W * BIN_MS} ms (pre-registered); the wide "
          f"band {WIDE_MIN * BIN_MS}-250 ms is disjoint from it. Pairs with "
          f"n <= {UNDERPOWERED_N} animals are UNDERPOWERED for a one-trial effect, "
          "not evidence of absence. Wilcoxon cannot reach p < 0.05 below n = 6.", ""]
    per_pair_means: dict[str, pd.DataFrame] = {}
    for suf, fs in FS_CONDITIONS:
        cpath = RES / f"trial12_curves_bin{BIN_MS}{suf}.csv"
        tpath = RES / f"trial12_trials_bin{BIN_MS}{suf}.csv"
        if not cpath.exists() or not tpath.exists():
            print(f"skip {fs}: {cpath.name} / {tpath.name} not found"); continue
        curves, trials = pd.read_csv(cpath), pd.read_csv(tpath)
        print(f"\n[{fs}] {cpath.name}: {len(curves)} rows, "
              f"{curves['animal'].nunique()} animals, "
              f"{curves.groupby(['animal', 'pair']).ngroups} animal-pairs")

        # Degenerate frozen projections (|r| ~ 1 at high rank) are an artefact, not
        # coupling — drop the affected trial cells from BOTH the curve and trial sides
        # before anything is averaged, and report what went.
        flags = degenerate_trials(trials)
        flags.to_csv(RES / f"trial12_degenerate_bin{BIN_MS}{suf}.csv", index=False,
                     lineterminator="\n")
        n_bad = int(flags["degenerate_trial"].sum())
        if n_bad:
            bad = flags[flags["degenerate_trial"] == 1]
            who = (bad.groupby(["animal", "pair"])["ordinal"]
                   .apply(lambda s: sorted(set(s))).to_dict())
            print(f"[{fs}] DEGENERATE frozen projections dropped: {n_bad} "
                  f"(animal, pair, ordinal, arm) cells — {who}")
        curves = drop_degenerate_trials(curves, flags)
        trials = drop_degenerate_trials(trials, flags)

        win = per_trial_windows(curves)
        wideb = wide_band(curves)
        # DIM-MATCHED before collapsing: the 1-vs-2 delta must average the SAME dims
        # in both trials (2026-08-20 audit — 85/213 headline cells had mismatched
        # dim-sets, which flipped the CA1-V1 sign). dim_matched applies the sig and
        # degenerate gates itself, so per_animal_mean must not re-apply them.
        allsig = cc_aggregate.per_animal_mean(
            dim_matched(win, "ifi"), value="ifi",
            by=["ordinal", "matched", "window_bins"], weight="cc_peak",
            sig_only=False, drop_degenerate=False)
        allsig_wide = cc_aggregate.per_animal_mean(
            dim_matched(wideb, "ifi"), value="ifi", by=["ordinal", "matched"],
            weight="cc_peak", sig_only=False, drop_degenerate=False)
        sig_rate = pd.to_numeric(win["sig"], errors="coerce").mean()
        print(f"[{fs}] sig gate pass rate: {sig_rate:.4f} — "
              f"{'VACUOUS: `all significant CCs` is `all CCs`' if sig_rate > 0.95 else 'selective'}")
        trials_allsig = cc_aggregate.per_animal_mean(
            trials, value="r0", by=["ordinal", "matched"], weight="r_frozen")

        comp = completeness(trials)
        comp.to_csv(RES / f"trial12_completeness_bin{BIN_MS}{suf}.csv", index=False,
                    lineterminator="\n")
        tests = headline_tests(win, wideb, allsig, allsig_wide, trials, trials_allsig)
        tests.to_csv(RES / f"trial12_delta_tests_bin{BIN_MS}{suf}.csv", index=False,
                     lineterminator="\n")
        sw = sweep(win, allsig)
        sw.to_csv(RES / f"trial12_sweep_bin{BIN_MS}{suf}.csv", index=False,
                  lineterminator="\n")
        swt = sweep_test(win, allsig)
        swt.to_csv(RES / f"trial12_sweep_test_bin{BIN_MS}{suf}.csv", index=False,
                   lineterminator="\n")
        if not swt.empty:
            print(f"[{fs}] window sweep (max|t|, sign-flip null): "
                  f"{int((swt['p_perm'] < 0.05).sum())} nominal, "
                  f"{int(swt['bh_pass'].sum())} BH-surviving of {len(swt)}")
        exists = direction_exists(win, allsig)
        per_pair_means[suf] = tests[["metric", "matched", "pair", "mean"]].copy()

        # controls on the two headline metrics
        ctrl_frames = []
        w1 = cc1_per_animal(win[win["window_bins"] == HEADLINE_W])
        for label, pa, col in (("ifi_cc1", w1, "ifi"),
                               ("ifi_allsig_w",
                                allsig[allsig["window_bins"] == HEADLINE_W], "wmean"),
                               ("r0_cc1", trials[trials["dim"] == 1], "r0")):
            cs = control_stats(pa, col)
            if cs.empty:
                continue
            for v in ("z", "abs_excess"):
                t = cc_aggregate.one_sample_by_pair(cs, v, PAIRS, min_n=MIN_N,
                                                    wilcoxon=True)
                if not t.empty:
                    t.insert(0, "metric", label); t.insert(1, "stat", v)
                    ctrl_frames.append(t)
        ctrl = (pd.concat(ctrl_frames, ignore_index=True)
                if ctrl_frames else pd.DataFrame())
        if not ctrl.empty:
            ctrl.to_csv(RES / f"trial12_control_bin{BIN_MS}{suf}.csv", index=False,
                        lineterminator="\n")

        n_hits = int((pd.to_numeric(tests["p"], errors="coerce") < 0.05).sum())
        print(f"[{fs}] {len(tests)} (metric, arm, pair) tests; "
              f"{n_hits} at p < 0.05 before any correction; "
              f"BH-across-pairs survivors: {int(tests['bh_pass'].sum())}")

        md += [f"## {fs}", "",
               f"### Completeness (CC1 r0; missingness points at short trials)", "",
               _md(comp), "",
               "### The confound: Δ(trial 1 − trial 2) in bins, speed, gaps "
               "(raw arm)", "",
               _md(tests[tests["metric"].str.startswith("d_")]), "",
               "### Strength (per-trial held-out r0 through frozen axes)", "",
               _md(tests[tests["metric"].str.startswith("r0")]), "",
               f"### Direction: IFI at ±{HEADLINE_W * BIN_MS} ms and the "
               f"{WIDE_MIN * BIN_MS}-250 ms band", "",
               _md(tests[tests["metric"].str.startswith(("ifi", "absifi"))]), "",
               "### Is there a direction at all, per trial? (matched arm 1)", "",
               _md(exists), "",
               "### Control: 1-vs-2 delta vs adjacent-trial variability "
               "(ordinals 3..10, common-matched arm)", "",
               "`z` = d12 / (√2·SD of the metric over ordinals 3..10) — |z| ≈ 1 is "
               "ordinary adjacent-trial jitter; `abs_excess` = |d12| − mean"
               "|adjacent delta|. A trend across 3..10 inflates the SD (conservative).",
               "", _md(ctrl) if not ctrl.empty else "_(insufficient control data)_",
               "",
               "### Integration-window sweep, TESTED (max |t| over windows, "
               "sign-flip permutation null)", "",
               "Cumulative windows are nested, so a per-window BH would be "
               "anti-conservative. Flipping the sign of each animal's whole "
               "delta-vs-window profile preserves the nesting exactly.", "",
               _md(swt) if not swt.empty else "_(insufficient data)_", ""]
    if len(per_pair_means) == 2:
        agree = fs_agreement(per_pair_means[""], per_pair_means["_fsincl"])
        agree.to_csv(RES / f"trial12_fs_agreement_bin{BIN_MS}.csv", index=False,
                     lineterminator="\n")
        md += ["## Do the two co-primary FS conditions agree?", "",
               "Per (metric, arm): correlation and sign-agreement of the per-pair "
               "delta means across FS conditions. They share animals and most units, "
               "so agreement is WEAK evidence; disagreement is strong evidence of "
               "FS-fragility.", "", _md(agree), ""]
        n_low = int((agree["r"] < 0.5).sum())
        print(f"\nFS agreement: {len(agree)} (metric, arm) cells, "
              f"{n_low} with r < 0.5 between conditions")
    (RES / "trial12_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / 'trial12_tables.md'}")


if __name__ == "__main__":
    main()
