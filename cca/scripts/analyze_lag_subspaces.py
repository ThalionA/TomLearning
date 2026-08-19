"""Lagged-subspace analysis — meeting 2026-07-28 items 2, 3 and 4.

Reads results/lag_subspaces_bin10{,_fsincl}.csv (written by run_lag_subspaces.py).

ITEM 3 — how stable are the CCs across time lags?
    Per pair and lag, the principal angle between the lag-0 subspace and the lagged one
    is compared to that pair's own split-half floor, animals-as-n (paired t on
    angle - floor). Bonferroni across the 10 non-zero |lag| values within a pair — the
    same nested-window correction the report already applies to the IFI sweep. The
    headline readout is the STABILITY WIDTH: the largest |lag| at which the subspace is
    still indistinguishable from its own noise floor.

ITEMS 2/4 — feedforward vs feedback.
    FF = the fit at +TAU (X leads), FB = the fit at -TAU. Three contrasts, all
    animals-as-n paired t:
      * cc1(FF) - cc1(FB)          is one direction more strongly coupled?
      * gini_conn(FF) - gini_conn(FB)  does one direction recruit more neurons?
        (connection-specific Gini — the partner-invariant one is a GOTCHA)
      * angle(FF, FB) vs floor     are they even separable subspaces, or one subspace
                                   read at two delays?
    The third is the gating question: if FF and FB sit at the floor, splitting the
    subspace into FF/FB components is not supported by this data and the first two
    contrasts describe one subspace, not two.

Writes results/lag_subspaces_tables.md + per-pair stats CSVs.

Usage: PYTHONPATH=src python scripts/analyze_lag_subspaces.py
"""
from __future__ import annotations


import numpy as np
import pandas as pd


from _common import FS_CONDITIONS, RES, TEMPORAL, config
from tom_cca import paired_stats
PAIRS = list(config.PAIR_NAMES)
TAU_BINS = TEMPORAL.label_w_bins   # must match run_lag_subspaces.TAU_BINS


def _paired(deltas):
    d = np.asarray([v for v in deltas if np.isfinite(v)], dtype=float)
    if d.size < 2:
        return len(d), np.nan, np.nan, np.nan
    n, _, t, p = paired_stats.paired_t(d)
    return n, float(np.mean(d)), t, p


NOT_ESTIMABLE_DEG = 70.0     # a split-half floor above this = no usable subspace estimate


def stability(df: pd.DataFrame, dims: int = 1) -> pd.DataFrame:
    """Item 3: angle-vs-floor at each lag, per pair, at ``dims`` canonical dimensions.

    **Read the floor before the p-value.** The floor is the angle between two halves of
    the SAME data at the SAME lag — i.e. how reproducible the subspace estimate is at
    all. At d=3 it comes out at ~78 deg, essentially orthogonal, so the lagged angle
    (~75 deg) has nowhere to go and "not different from floor" means UNMEASURABLE, not
    stable. At d=1 the floor is ~53 deg and there is real headroom. CC1 is therefore the
    primary readout and d=3 is reported only to document the power failure — hence the
    `estimable` column, which must be checked before any verdict is read off `p_bonf`.

    **The floor is CONSERVATIVE, and by a known amount.** It is measured from two
    half-data fits, whereas the lag-0-vs-lagged comparison uses the full window — fewer
    samples per fit means a noisier subspace and a larger angle, so the floor sits above
    the true noise level of the comparison it gates. The same asymmetry is visible in the
    report's existing reorientation table (§G of `bin10_tables.md`), where cross-window
    rotations come out significantly *below* their floors. Consequence here: "at floor"
    excludes a LARGE rotation with lag, but the test can miss a modest one.
    """
    ax, ay, fx, fy = (("angle_x_cc1", "angle_y_cc1", "floor_x_cc1", "floor_y_cc1")
                      if dims == 1 else
                      ("angle_x", "angle_y", "floor_x", "floor_y"))
    rows = []
    for pair in PAIRS:
        sub = df[df["pair"] == pair]
        if sub.empty:
            continue
        lags = sorted(x for x in sub["lag_ms"].unique() if x != 0)
        n_tests = len({abs(x) for x in lags})
        for lag in lags:
            at = sub[sub["lag_ms"] == lag]
            per_animal, kept_angles, kept_floors, n_areas = [], [], [], 0
            for _, g in at.groupby("animal"):
                # GATE PER AREA, then average — not average, then gate. An area whose OWN
                # floor exceeds the threshold has no usable subspace estimate, and its
                # (angle - floor) term is large and negative, so averaging it in drags the
                # animal toward "at floor". That biases TOWARD the null this test reports,
                # i.e. it is not a conservative error. The animal drops out entirely when
                # neither area is estimable.
                terms = []
                for a_col, f_col in ((ax, fx), (ay, fy)):
                    a_v, f_v = float(g[a_col].iloc[0]), float(g[f_col].iloc[0])
                    if np.isfinite(f_v) and f_v < NOT_ESTIMABLE_DEG and np.isfinite(a_v):
                        terms.append(a_v - f_v)
                        kept_angles.append(a_v); kept_floors.append(f_v)
                if terms:
                    per_animal.append(float(np.mean(terms)))
                    n_areas += len(terms)
            n, mean, t, p = _paired(per_animal)
            rows.append({
                "pair": pair, "dims": dims, "lag_ms": lag, "n_animals": n,
                "n_areas_used": n_areas,
                "angle_minus_floor": mean, "t": t, "p": p,
                "p_bonf": min(1.0, p * n_tests) if np.isfinite(p) else np.nan,
                "mean_angle": float(np.mean(kept_angles)) if kept_angles else np.nan,
                "mean_floor": float(np.mean(kept_floors)) if kept_floors else np.nan,
                "estimable": bool(n_areas > 0),
            })
    return pd.DataFrame(rows)


GATE_SWEEP = [50.0, 60.0, 70.0, 80.0, 90.0, float("inf")]


def gate_sensitivity(df: pd.DataFrame, dims: int = 1) -> pd.DataFrame:
    """How much of item 3's answer is the analyst's choice of estimability threshold?

    **Both directions of this choice are biased, and by a measurable amount.** Across
    area-lags, corr(floor, angle - floor) is Spearman rho = -0.57 (p ~ 1e-242): a high
    floor mechanically produces a negative delta, because the floor is built from
    half-data fits while the comparison uses the full window. So

      * NO gate  -> unmeasurable areas (floor near 90 deg) contribute strongly negative
                    deltas and drag every pair toward "at floor": biased toward the NULL.
      * TIGHT gate -> selects low-floor areas, whose deltas are positive by the same
                    correlation: biased toward finding ROTATION.

    There is no threshold-free answer, so the honest output is the whole sweep. A pair
    that appears at EVERY gate, including no gate, is the only kind of claim this test
    supports.
    """
    out = []
    global NOT_ESTIMABLE_DEG
    original = NOT_ESTIMABLE_DEG
    try:
        for thr in GATE_SWEEP:
            NOT_ESTIMABLE_DEG = thr
            s = stability(df, dims=dims)
            up = s[(s["p_bonf"] < 0.05) & (s["angle_minus_floor"] > 0)]
            out.append({
                "gate_deg": thr, "n_rotating_lags": int(len(up)),
                "pairs": ", ".join(sorted(up["pair"].unique())) or "—",
                "n_pairs": int(up["pair"].nunique()),
            })
    finally:
        NOT_ESTIMABLE_DEG = original
    return pd.DataFrame(out)


def _md_gate(sens: pd.DataFrame, fs: str) -> str:
    robust = None
    sets = [set(r["pairs"].split(", ")) - {"—"} for _, r in sens.iterrows()]
    if sets:
        inter = set.intersection(*sets) if all(sets) else set()
        robust = ", ".join(sorted(inter)) if inter else "none"
    lines = [f"#### {fs} — item 3 depends on the estimability threshold", "",
             "The gate is an analyst choice with no principled value, and **both "
             "directions are biased**: across area-lags corr(floor, angle−floor) is "
             "Spearman ρ = −0.57, so no gate drags toward the null (unmeasurable areas "
             "contribute large negative deltas) while a tight gate selects low-floor "
             "areas and drags toward rotation.", "",
             "| gate | rotating lags | pairs |", "|---|---|---|"]
    for _, r in sens.iterrows():
        g = "no gate" if not np.isfinite(r["gate_deg"]) else f"{r['gate_deg']:.0f}°"
        lines.append(f"| {g} | {r['n_rotating_lags']:.0f} | {r['pairs']} |")
    lines += ["", f"**Robust to every gate, including none: {robust}.** That is the only "
              "claim this test supports; the lag COUNT is not interpretable."]
    return "\n".join(lines) + "\n"


def stability_width(stab: pd.DataFrame) -> pd.DataFrame:
    """Largest |lag| (ms) at which the subspace is still AT its own noise floor.

    Scanned outward from 0: the width ends at the first |lag| whose Bonferroni-corrected
    paired t says the subspace has moved. Reported as `>= max_lag` when no lag in the
    swept range separates from the floor.
    """
    rows = []
    for pair, g in stab.groupby("pair"):
        by_abs = (g.assign(abs_lag=g["lag_ms"].abs())
                   .groupby("abs_lag")
                   .agg(p_bonf=("p_bonf", "min"),
                        angle_minus_floor=("angle_minus_floor", "mean"))
                   .sort_index())
        width, censored = None, False
        for abs_lag, r in by_abs.iterrows():
            if np.isfinite(r["p_bonf"]) and r["p_bonf"] < 0.05 and \
                    r["angle_minus_floor"] > 0:
                break
            width = abs_lag
        else:
            censored = True
        rows.append({"pair": pair, "stability_width_ms": width,
                     "censored_at_max_lag": censored,
                     "max_lag_swept_ms": float(by_abs.index.max())})
    return pd.DataFrame(rows)


def ff_fb(df: pd.DataFrame) -> pd.DataFrame:
    """Items 2/4: feedforward (+TAU) vs feedback (-TAU)."""
    tau_ms = TAU_BINS * int(df["bin_ms"].iloc[0])
    rows = []
    for pair in PAIRS:
        sub = df[df["pair"] == pair]
        if sub.empty:
            continue
        cc_d, gini_d, ang_d, ff_cc, fb_cc = [], [], [], [], []
        for _, g in sub.groupby("animal"):
            ff = g[g["lag_ms"] == tau_ms]
            fb = g[g["lag_ms"] == -tau_ms]
            if ff.empty or fb.empty:
                continue
            ff, fb = ff.iloc[0], fb.iloc[0]
            cc_d.append(float(ff["cc1"]) - float(fb["cc1"]))
            ff_cc.append(float(ff["cc1"])); fb_cc.append(float(fb["cc1"]))
            gini_d.append(np.nanmean([float(ff["gini_x_conn"]),
                                      float(ff["gini_y_conn"])]) -
                          np.nanmean([float(fb["gini_x_conn"]),
                                      float(fb["gini_y_conn"])]))
            # CC1, not 3 dims: the 3-dim split-half floor is ~78 deg (the subspace is
            # not estimable there), so an FF/FB angle measured at 3 dims has no
            # headroom and the comparison would be vacuous. Both terms are x-area.
            ang_d.append(float(ff["angle_ff_fb_cc1"]) - float(ff["floor_x_cc1"]))
        n_cc, m_cc, t_cc, p_cc = _paired(cc_d)
        n_g, m_g, t_g, p_g = _paired(gini_d)
        n_a, m_a, t_a, p_a = _paired(ang_d)
        rows.append({
            "pair": pair, "n_animals": n_cc,
            "cc1_ff": float(np.mean(ff_cc)) if ff_cc else np.nan,
            "cc1_fb": float(np.mean(fb_cc)) if fb_cc else np.nan,
            "d_cc1": m_cc, "t_cc1": t_cc, "p_cc1": p_cc,
            "d_gini_conn": m_g, "t_gini": t_g, "p_gini": p_g,
            "n_angle": n_a, "ff_fb_angle_minus_floor": m_a, "t_angle": t_a,
            "p_angle": p_a,
        })
    return pd.DataFrame(rows)


def ff_fb_evolution(df: pd.DataFrame) -> pd.DataFrame:
    """Item 2's second half + item 4: does the FF/FB picture CHANGE with learning?

    Operates on the ``--epochs`` output. Per pair, expert vs naive paired t on four
    quantities: the feedforward and feedback CC1 separately, their difference (the
    directional asymmetry), and the FF/FB subspace angle. Animals-as-n, learners only.

    A change in `d_cc1` is the interesting one — it is the directional asymmetry moving
    with learning, which is what "check their evolution with learning" asks for. But read
    it against the session-level gate first: if FF and FB are not separable subspaces to
    begin with, a change in their difference is a change within ONE subspace.
    """
    tau_ms = TAU_BINS * int(df["bin_ms"].iloc[0])
    rows = []
    for pair in PAIRS:
        sub = df[df["pair"] == pair]
        if sub.empty:
            continue
        acc = {k: [] for k in ("ff", "fb", "asym", "angle")}
        for _, g in sub.groupby("animal"):
            vals = {}
            for epoch in ("naive", "expert"):
                e = g[g["epoch"] == epoch]
                ff = e[e["lag_ms"] == tau_ms]
                fb = e[e["lag_ms"] == -tau_ms]
                if ff.empty or fb.empty:
                    break
                vals[epoch] = {
                    "ff": float(ff["cc1"].iloc[0]), "fb": float(fb["cc1"].iloc[0]),
                    "asym": float(ff["cc1"].iloc[0]) - float(fb["cc1"].iloc[0]),
                    "angle": float(ff["angle_ff_fb_cc1"].iloc[0]),
                }
            if len(vals) == 2:
                for k in acc:
                    acc[k].append(vals["expert"][k] - vals["naive"][k])
        row = {"pair": pair}
        for k, label in [("ff", "ff"), ("fb", "fb"), ("asym", "asym"),
                         ("angle", "angle")]:
            n, m, t, p = _paired(acc[k])
            row |= {f"n_{label}": n, f"d_{label}": m, f"t_{label}": t, f"p_{label}": p}
        rows.append(row)
    return pd.DataFrame(rows)


def _md_evolution(tab: pd.DataFrame, fs: str, tau_ms: int) -> str:
    lines = [f"### {fs} — item 2/4: does the FF/FB picture change with learning?", "",
             f"Expert − naive, animals-as-n paired *t*, learners only. FF = +{tau_ms} ms "
             f"(first area leads), FB = −{tau_ms} ms. `asym` = cc₁(FF) − cc₁(FB), the "
             "directional asymmetry.", "",
             "> Read against the session-level gate: if FF and FB are not separable "
             "> subspaces, a change in `asym` is a change *within one* subspace, not a "
             "> shift between two streams.", "",
             "| pair | n | Δ cc₁ FF | p | Δ cc₁ FB | p | Δ asym | p | Δ FF/FB angle | p |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in tab.iterrows():
        lines.append(
            f"| {r['pair']} | {r['n_asym']:.0f} | {r['d_ff']:+.3f} | {r['p_ff']:.3g} | "
            f"{r['d_fb']:+.3f} | {r['p_fb']:.3g} | {r['d_asym']:+.3f} | "
            f"{r['p_asym']:.3g} | {r['d_angle']:+.1f}° | {r['p_angle']:.3g} |")
    return "\n".join(lines) + "\n"


def _md_stability(stab, width, fs, dims):
    what = "CC₁ only" if dims == 1 else "3 canonical dims"
    lines = [f"### {fs} — item 3: subspace stability across lag ({what})", "",
             "`angle − floor` is the principal angle between the lag-0 and lagged "
             "subspace minus that pair's own split-half floor, averaged over the two "
             "areas; animals-as-n, Bonferroni across |lag| within a pair.", ""]
    if dims != 1:
        lines += ["> **⚠ This table is a power check, not a result.** The split-half "
                  "floor at 3 dims is ~78°, i.e. two halves of the *same* data at the "
                  "*same* lag are nearly orthogonal — the 3-dim subspace is not "
                  "estimable at this N. A lagged angle that fails to exceed that floor "
                  "means UNMEASURABLE, not stable. Read the CC₁ table instead.", ""]
    lines += ["| pair | estimable? | stability width | mean angle @ ±50 ms | floor | "
              "Δ | p (Bonf) |", "|---|---|---|---|---|---|---|"]
    for _, w in width.iterrows():
        g = stab[(stab["pair"] == w["pair"]) & (stab["lag_ms"].abs() == 50)]
        ang = g["mean_angle"].mean() if len(g) else np.nan
        flr = g["mean_floor"].mean() if len(g) else np.nan
        d = g["angle_minus_floor"].mean() if len(g) else np.nan
        p = g["p_bonf"].min() if len(g) else np.nan
        est = bool(g["estimable"].all()) if len(g) else False
        wid = ("≥ %.0f ms (censored)" % w["max_lag_swept_ms"]
               if w["censored_at_max_lag"] else
               ("%.0f ms" % w["stability_width_ms"]
                if pd.notna(w["stability_width_ms"]) else "< smallest lag"))
        if not est:
            wid = "n/a — not estimable"
        lines.append(f"| {w['pair']} | {'yes' if est else '**NO**'} | {wid} | "
                     f"{ang:.1f}° | {flr:.1f}° | {d:+.1f}° | {p:.3g} |")
    return "\n".join(lines) + "\n"


def _md_fffb(tab, fs, tau_ms):
    lines = [f"### {fs} — items 2/4: feedforward (+{tau_ms} ms) vs feedback "
             f"(−{tau_ms} ms)", "",
             "Positive `Δcc₁` = the first-named area leading is more strongly coupled. "
             "`Δgini_conn` uses the CONNECTION-SPECIFIC Gini. `FF/FB angle − floor` is "
             "the gate: at or below 0 the two are one subspace read at two delays.", "",
             "| pair | n | cc₁ FF | cc₁ FB | Δcc₁ | p | Δgini_conn | p | "
             "FF/FB angle − floor | p |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in tab.iterrows():
        lines.append(
            f"| {r['pair']} | {r['n_animals']:.0f} | {r['cc1_ff']:.3f} | "
            f"{r['cc1_fb']:.3f} | {r['d_cc1']:+.3f} | {r['p_cc1']:.3g} | "
            f"{r['d_gini_conn']:+.3f} | {r['p_gini']:.3g} | "
            f"{r['ff_fb_angle_minus_floor']:+.1f}° | {r['p_angle']:.3g} |")
    return "\n".join(lines) + "\n"


def main():
    md = ["# Lagged communication subspaces — meeting items 2, 3, 4", "",
          "Animals-as-n throughout; 8 pairs, no cross-pair correction (per-pair family, "
          "STATE.md §3.0 policy). Angles are the LARGEST principal angle over 3 "
          "canonical dims, so a subspace that matches on its dominant axis but diverges "
          "elsewhere is not scored as stable.", ""]
    for suf, fs in FS_CONDITIONS:
        src = RES / f"lag_subspaces_bin10{suf}.csv"
        if not src.exists():
            print(f"skip {fs}: {src.name} not found"); continue
        df = pd.read_csv(src)
        for c in ("cc1", "cc_mean3", "angle_x", "angle_y", "angle_x_cc1",
                  "angle_y_cc1", "floor_x", "floor_y", "floor_x_cc1",
                  "floor_y_cc1", "gini_x_conn", "gini_y_conn", "angle_ff_fb_x",
                  "angle_ff_fb_y", "angle_ff_fb_cc1"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        stab1 = stability(df, dims=1)
        stab3 = stability(df, dims=3)
        w1, w3 = stability_width(stab1), stability_width(stab3)
        tab = ff_fb(df)
        pd.concat([stab1, stab3]).to_csv(
            RES / f"lag_subspaces_stability_bin10{suf}.csv", index=False,
            lineterminator="\n")
        tab.to_csv(RES / f"lag_subspaces_fffb_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        tau_ms = TAU_BINS * int(df["bin_ms"].iloc[0])
        sens = gate_sensitivity(df, dims=1)
        sens.to_csv(RES / f"lag_subspaces_gate_sensitivity_bin10{suf}.csv", index=False,
                    lineterminator="\n")
        md += [_md_stability(stab1, w1, fs, 1), "", _md_gate(sens, fs), "",
               _md_stability(stab3, w3, fs, 3), "", _md_fffb(tab, fs, tau_ms), ""]
        print(f"\n{fs}: {df['animal'].nunique()} animals, {len(df)} rows")
        print(f"  split-half floor: CC₁ {stab1['mean_floor'].mean():.1f}° | "
              f"3-dim {stab3['mean_floor'].mean():.1f}°"
              f"  (>{NOT_ESTIMABLE_DEG:.0f}° = subspace not estimable)")
        n_est = stab1.groupby('pair')['estimable'].all().sum()
        print(f"  pairs with an estimable CC₁ subspace: {n_est}/{len(PAIRS)}")
        sets = [set(r["pairs"].split(", ")) - {"—"}
                for _, r in gate_sensitivity(df, dims=1).iterrows()]
        rob = set.intersection(*sets) if sets and all(sets) else set()
        print(f"  rotating at EVERY estimability gate (the only robust claim): "
              f"{', '.join(sorted(rob)) if rob else 'none'}")
        sep = tab[tab["p_angle"] < 0.05]
        print(f"  FF/FB separable from floor: "
              f"{len(sep)}/{len(tab)} pairs" +
              (" — " + ", ".join(sep['pair']) if len(sep) else ""))
        hits = tab[tab["p_cc1"] < 0.05]
        print(f"  FF vs FB strength differs: {len(hits)}/{len(tab)} pairs" +
              (" — " + ", ".join(hits['pair']) if len(hits) else ""))
        est1 = stab1.groupby("pair")["estimable"].all()
        cens = w1[w1["censored_at_max_lag"] & w1["pair"].map(est1).fillna(False)]
        print(f"  CC₁ subspace at floor across the WHOLE swept range: "
              f"{len(cens)}/{int(est1.sum())} estimable pairs")

        # items 2/4 — the FF/FB evolution, if the --epochs run has landed
        ep_src = RES / f"lag_subspaces_bin10_epochs{suf}.csv"
        needed = {"cc1", "angle_ff_fb_cc1", "floor_x_cc1", "floor_y_cc1",
                  "gini_x_conn", "gini_y_conn"}
        if ep_src.exists() and needed <= set(pd.read_csv(ep_src, nrows=0).columns):
            ep = pd.read_csv(ep_src)
            for c in needed:
                ep[c] = pd.to_numeric(ep[c], errors="coerce")
            evo = ff_fb_evolution(ep)
            evo.to_csv(RES / f"lag_subspaces_evolution_bin10{suf}.csv", index=False,
                       lineterminator="\n")
            md.append(_md_evolution(evo, fs, tau_ms))
            hits = evo[evo["p_asym"] < 0.05]
            print(f"  FF/FB asymmetry changes with learning: {len(hits)}/{len(evo)} "
                  f"pairs" + (" — " + ", ".join(hits["pair"]) if len(hits) else ""))
        elif ep_src.exists():
            print(f"  (skipping {ep_src.name}: written before the per-area d=1 floor "
                  "fix — re-run run_lag_subspaces.py --epochs)")
        else:
            print(f"  (epoch run pending — {ep_src.name} not found)")
    (RES / "lag_subspaces_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / 'lag_subspaces_tables.md'}")


if __name__ == "__main__":
    main()
