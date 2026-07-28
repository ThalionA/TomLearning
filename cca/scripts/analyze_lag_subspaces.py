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

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import paired_stats  # noqa: E402

RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
TAU_BINS = 5


def _paired(deltas):
    d = np.asarray([v for v in deltas if np.isfinite(v)], dtype=float)
    if d.size < 2:
        return len(d), np.nan, np.nan, np.nan
    n, _, t, p = paired_stats.paired_t(d)
    return n, float(np.mean(d)), t, p


def stability(df: pd.DataFrame) -> pd.DataFrame:
    """Item 3: angle-vs-floor at each lag, per pair."""
    rows = []
    for pair in PAIRS:
        sub = df[df["pair"] == pair]
        if sub.empty:
            continue
        lags = sorted(x for x in sub["lag_ms"].unique() if x != 0)
        n_tests = len({abs(x) for x in lags})
        for lag in lags:
            at = sub[sub["lag_ms"] == lag]
            # each area contributes its own angle against its own floor
            deltas = np.concatenate([
                (at["angle_x"] - at["floor_x"]).to_numpy(float),
                (at["angle_y"] - at["floor_y"]).to_numpy(float)])
            per_animal = {}
            for animal, g in at.groupby("animal"):
                per_animal[animal] = np.nanmean([
                    float(g["angle_x"].iloc[0]) - float(g["floor_x"].iloc[0]),
                    float(g["angle_y"].iloc[0]) - float(g["floor_y"].iloc[0])])
            n, mean, t, p = _paired(list(per_animal.values()))
            rows.append({
                "pair": pair, "lag_ms": lag, "n_animals": n,
                "angle_minus_floor": mean, "t": t, "p": p,
                "p_bonf": min(1.0, p * n_tests) if np.isfinite(p) else np.nan,
                "mean_angle": float(np.nanmean([at["angle_x"].mean(),
                                                at["angle_y"].mean()])),
                "mean_floor": float(np.nanmean([at["floor_x"].mean(),
                                                at["floor_y"].mean()])),
                "n_obs": int(np.sum(np.isfinite(deltas))),
            })
    return pd.DataFrame(rows)


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
            ang_d.append(np.nanmean([float(ff["angle_ff_fb_x"]),
                                     float(ff["angle_ff_fb_y"])]) -
                         np.nanmean([float(ff["floor_x"]), float(ff["floor_y"])]))
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


def _md_stability(stab, width, fs):
    lines = [f"### {fs} — item 3: subspace stability across lag", "",
             "`angle − floor` is the principal angle between the lag-0 and lagged "
             "subspace minus that pair's own split-half floor, averaged over the two "
             "areas; animals-as-n, Bonferroni across |lag| within a pair.", "",
             "| pair | stability width | mean angle @ ±50 ms | floor | Δ | p (Bonf) |",
             "|---|---|---|---|---|---|"]
    for _, w in width.iterrows():
        g = stab[(stab["pair"] == w["pair"]) & (stab["lag_ms"].abs() == 50)]
        ang = g["mean_angle"].mean() if len(g) else np.nan
        flr = g["mean_floor"].mean() if len(g) else np.nan
        d = g["angle_minus_floor"].mean() if len(g) else np.nan
        p = g["p_bonf"].min() if len(g) else np.nan
        wid = ("≥ %.0f ms (censored)" % w["max_lag_swept_ms"]
               if w["censored_at_max_lag"] else
               ("%.0f ms" % w["stability_width_ms"]
                if pd.notna(w["stability_width_ms"]) else "< smallest lag"))
        lines.append(f"| {w['pair']} | {wid} | {ang:.1f}° | {flr:.1f}° | "
                     f"{d:+.1f}° | {p:.3g} |")
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
    for suf, fs in [("", "FS-excluded"), ("_fsincl", "FS-included")]:
        src = RES / f"lag_subspaces_bin10{suf}.csv"
        if not src.exists():
            print(f"skip {fs}: {src.name} not found"); continue
        df = pd.read_csv(src)
        for c in ("cc1", "cc_mean3", "angle_x", "angle_y", "floor_x", "floor_y",
                  "gini_x_conn", "gini_y_conn", "angle_ff_fb_x", "angle_ff_fb_y"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        stab = stability(df)
        width = stability_width(stab)
        tab = ff_fb(df)
        stab.to_csv(RES / f"lag_subspaces_stability_bin10{suf}.csv", index=False,
                    lineterminator="\n")
        tab.to_csv(RES / f"lag_subspaces_fffb_bin10{suf}.csv", index=False,
                   lineterminator="\n")
        tau_ms = TAU_BINS * int(df["bin_ms"].iloc[0])
        md += [_md_stability(stab, width, fs), "", _md_fffb(tab, fs, tau_ms), ""]
        print(f"\n{fs}: {df['animal'].nunique()} animals, {len(df)} rows")
        sep = tab[tab["p_angle"] < 0.05]
        print(f"  FF/FB separable from floor: "
              f"{len(sep)}/{len(tab)} pairs" +
              (" — " + ", ".join(sep['pair']) if len(sep) else ""))
        hits = tab[tab["p_cc1"] < 0.05]
        print(f"  FF vs FB strength differs: {len(hits)}/{len(tab)} pairs" +
              (" — " + ", ".join(hits['pair']) if len(hits) else ""))
        cens = width[width["censored_at_max_lag"]]
        print(f"  subspace at floor across the WHOLE swept range: "
              f"{len(cens)}/{len(width)} pairs")
    (RES / "lag_subspaces_tables.md").write_text("\n".join(md))
    print(f"\nwrote {RES / 'lag_subspaces_tables.md'}")


if __name__ == "__main__":
    main()
