"""Item 1 figures — IFI of each CC across integration windows, and the sign structure.

Reads results/cc_ifi_windows_bin10{,_fsincl}.csv (analyze_cc_ifi_signs.py).

Figure 1  HCV1_cc_ifi_windows_<fs>_bin10
          "Calculate IFI at different lags for each CC and look at them." Per pair,
          IFI vs integration window (|lag| <= w): the significant CCs split into a
          +IFI and a -IFI group at +/-50 ms, PLUS (2026-08-07 ask 1) the ungrouped
          mean over ALL significant CCs, with its one-sample test vs 0 at +/-50 ms
          (ask 3) in the panel title. Every line is a per-ANIMAL-first mean +/- SEM
          across animals (cc_aggregate) — an animal with five CCs counts once.

Figure 2  HCV1_cc_ifi_signmap_<fs>_bin10
          The same thing per animal, at the headline +/-50 ms window: a diverging
          heatmap of IFI, rows = animals, columns = canonical dimension. Red = the
          first-named area leads, blue = it lags. A row with both colours is an animal
          whose CCs disagree in direction.

Positive IFI => the FIRST-named area leads.

Usage: PYTHONPATH=src python scripts/figs_cc_ifi_signs.py [fsincl]
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402
from tom_cca import cc_aggregate  # noqa: E402

figstyle.apply()
ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"
RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
HEADLINE_W = 5         # bins = +/-50 ms
REF_W = 5              # window (bins) at which CCs are labelled +IFI / -IFI
BIN_MS = 10
MAP_DIMS = 8


def _group_curves(g: pd.DataFrame):
    """Per-animal-first mean, SEM and CC-strength-weighted mean of IFI vs window.

    ``g`` is one row per (animal, dim, window). Each animal's CCs are collapsed to one
    value per window first (cc_aggregate.per_animal_mean), then averaged across
    animals — so the SEM band is animals-as-n, the same unit as the test in the
    title. Returns ``(windows_ms, mean, sem, wmean, n_animals)``.
    """
    pa = cc_aggregate.per_animal_mean(g, value="ifi", by=["window_ms"],
                                      weight="cc_peak", sig_only=False,
                                      drop_degenerate=False)
    piv = pa.pivot_table(index="animal", columns="window_ms", values="mean")
    wins = np.array(sorted(piv.columns), dtype=float)
    M = piv[sorted(piv.columns)].to_numpy(float)
    n = np.sum(np.isfinite(M), axis=0)
    mean = np.nanmean(M, axis=0)
    sd = np.nanstd(M, axis=0, ddof=1)
    sem = np.where(n > 1, sd / np.sqrt(np.maximum(n, 1)), np.nan)
    W = pa.pivot_table(index="animal", columns="window_ms", values="wmean")
    wmean = np.nanmean(W[sorted(W.columns)].to_numpy(float), axis=0)
    return wins, mean, sem, wmean, int(piv.shape[0])


def fig_windows(df: pd.DataFrame, fs: str, overall_test: pd.DataFrame | None = None):
    """+IFI, −IFI and ALL significant CCs, per-animal-first mean +/- SEM vs window.

    Each significant CC is labelled by the SIGN of its IFI at the reference window
    (+/-REF_W bins), then both groups — and the ungrouped union (2026-08-07 ask 1) —
    are plotted across all windows. The title carries the ungrouped mean's one-sample
    test vs 0 at +/-50 ms (ask 3; `analyze_cc_ifi_signs.overall_direction`).

    ⚠ The +/− split is CIRCULAR at the reference window itself — those CCs were
    selected for being above/below zero there, so a gap at the dotted line is
    guaranteed and means nothing. ⚠ And because the windows are CUMULATIVE
    (|lag| <= w), every wider window still CONTAINS the labelling lags, which carry
    the largest CCs — so the groups staying apart at w > 50 ms is mostly the same
    circularity leaking outward, not independent evidence (verified 2026-08-15: on
    the DISJOINT lags 60–250 ms alone the +/− gap collapses to ~0 in every pair).
    The all-CC line is not sign-circular: it never looked at the sign — though its
    population is still 'CCs significant at lag 0', i.e. the better-coupled cells.
    """
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        # ONLY CCs that beat the circular-shift null. Ungated, CC3-CC8 clear it in just
        # 14-17% of cells with held-out CC at the floor, so their IFI sign is noise.
        sub = df[(df["pair"] == pair) & (df["degenerate"] == 0) & (df["sig"] == 1)]
        if sub.empty:
            ax.set_title(f"{pair}\n(no significant CCs)", fontsize=9)
            ax.axis("off"); continue
        # label each (animal, CC) by its sign at the reference window; the ALL group
        # takes every finite CC there (an exact 0 has no sign but is still a CC)
        ref = sub[sub["window_bins"] == REF_W].set_index(["animal", "dim"])["ifi"]
        ref = ref[np.isfinite(ref)]
        if ref.empty:
            ax.set_title(f"{pair}\n(no usable CCs)", fontsize=9)
            ax.axis("off"); continue
        keyed = sub.set_index(["animal", "dim"])
        x_lead = pair.split("-")[0]
        groups = [(1, "#c0392b", f"+IFI ({x_lead} leads)"),
                  (-1, "#2c6fbb", "−IFI"),
                  (0, "#222222", "ALL sig. CCs")]
        for sign, colour, lab in groups:
            members = ref.index if sign == 0 else ref.index[np.sign(ref.to_numpy()) == sign]
            if len(members) < 2:
                continue
            g = keyed.loc[members].reset_index()
            wins, mean, sem, wmean, n_an = _group_curves(g)
            # Strength-weighted mean alongside the plain one. |IFI| falls with coupling
            # strength (Spearman -0.25 / -0.33, p < 1e-4), i.e. the extreme IFIs sit on
            # the WEAKEST CCs — the signature of a noisy estimator, since a poorly
            # estimated lag curve yields a large |IFI| by chance. Weighting by CC
            # strength therefore SHRINKS the +/− gap (6-20%) rather than sharpening it.
            lw_dash, lw_solid = (1.6, 3.0) if sign == 0 else (1.2, 2.4)
            ax.plot(wins, mean, "--", color=colour, lw=lw_dash, alpha=0.8, zorder=3,
                    label=f"{lab} · {len(members)} CCs, {n_an} an. (unweighted)")
            ax.fill_between(wins, mean - sem, mean + sem, color=colour, alpha=0.13,
                            zorder=2)
            ax.plot(wins, wmean, "-", color=colour, lw=lw_solid, zorder=4,
                    label=f"{lab} · CC-strength weighted")
        ax.axhline(0.0, color="k", lw=1.0, ls="--", zorder=1)
        ax.axvline(REF_W * BIN_MS, color="#666", lw=1.0, ls=":", zorder=1)
        title = f"{pair}  ({sub['animal'].nunique()} animals)"
        if overall_test is not None:
            row = overall_test[overall_test["pair"] == pair]
            if not row.empty:
                r = row.iloc[0]
                star = "*" if r["p"] < 0.05 else ""
                title += (f"\nall-CC IFI @±{int(r['window_ms'])}: {r['mean']:+.3f}, "
                          f"p={r['p']:.3f}{star} ({int(r['n'])} animals)")
                if "ccs_n" in row.columns and np.isfinite(r["ccs_n"]):
                    cs = "*" if r["ccs_p"] < 0.05 else ""
                    title += f"\nCCs as n: p={r['ccs_p']:.3g}{cs} ({int(r['ccs_n'])} CCs)"
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("integration window ±w (ms)", fontsize=8)
        ax.set_ylabel(f"IFI   +: {x_lead} leads", fontsize=8)
        ax.legend(fontsize=5.2, loc="best", framealpha=0.6)
    fig.suptitle(f"IFI vs integration window, significant CCs — {fs} | red/blue = CCs "
                 f"grouped by IFI sign at the dotted line (±{REF_W * BIN_MS} ms): circular "
                 f"there AND leaking outward, since every wider window still contains those "
                 f"lags; black = ALL significant CCs (not sign-circular)\n"
                 f"solid = CC-strength weighted, dashed = unweighted ± SEM, per-animal-first "
                 f"| title: one-sample t of the all-CC mean vs 0 at ±{REF_W * BIN_MS} ms, "
                 f"animals-as-n (n = animals with ≥1 sig. CC), no cross-pair correction | "
                 f"data = whole session (uncapped)", fontsize=10)
    figstyle.save(fig, ATT / f"HCV1_cc_ifi_windows_{fs}_bin10.png")
    print(f"wrote HCV1_cc_ifi_windows_{fs}_bin10.png")


def fig_signmap(df: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    sub_all = df[(df["window_bins"] == HEADLINE_W) & (df["dim"] <= MAP_DIMS) &
                 (df["sig"] == 1)]
    vmax = float(np.nanpercentile(np.abs(sub_all["ifi"].to_numpy(float)), 95)) or 1.0
    im = None
    for ax, pair in zip(axes, PAIRS):
        sub = sub_all[sub_all["pair"] == pair]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        piv = sub.pivot_table(index="animal", columns="dim", values="ifi")
        piv = piv.reindex(columns=range(1, MAP_DIMS + 1))
        M = piv.to_numpy(float)
        # NaN = that CC did not beat the null for that animal. In RdBu_r a NaN and a
        # genuine IFI ~ 0 would both render near-white, so non-significant cells get an
        # explicit grey instead — "no reliable coupling" must not look like "balanced".
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#dcdcdc")
        im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-vmax, vmax=vmax,
                       aspect="auto")
        # mark animals whose CCs disagree in sign
        n_mixed = 0
        for r in range(M.shape[0]):
            row = M[r][np.isfinite(M[r]) & (M[r] != 0)]
            if row.size and (row > 0).any() and (row < 0).any():
                n_mixed += 1
                ax.text(-0.9, r, "◀", va="center", ha="center", fontsize=7,
                        color="#c0392b")
        x_lead = pair.split("-")[0]
        ax.set_title(f"{pair}   {n_mixed}/{M.shape[0]} animals mixed", fontsize=9)
        ax.set_xticks(range(MAP_DIMS))
        ax.set_xticklabels([f"CC{d}" for d in range(1, MAP_DIMS + 1)], fontsize=6)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels([f"{a} ({int(np.isfinite(M[i]).sum())})"
                            for i, a in enumerate(piv.index)], fontsize=6)
        ax.set_xlabel(f"red: {x_lead} leads", fontsize=7)
    if im is not None:
        cb = fig.colorbar(im, ax=list(axes), fraction=0.015, pad=0.01)
        cb.set_label("IFI (red = first area leads)", fontsize=8)
    fig.suptitle(f"Sign of each CC's IFI per animal, ±{HEADLINE_W * 10} ms window — "
                 f"{fs} | ONLY CCs beating the null (grey = not significant) | "
                 "◀ = that animal's significant CCs disagree in direction", fontsize=12)
    figstyle.save(fig, ATT / f"HCV1_cc_ifi_signmap_{fs}_bin10.png")
    print(f"wrote HCV1_cc_ifi_signmap_{fs}_bin10.png")


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    src = RES / f"cc_ifi_windows_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/analyze_cc_ifi_signs.py first")
    df = pd.read_csv(src)
    df["ifi"] = pd.to_numeric(df["ifi"], errors="coerce")
    test_csv = RES / f"cc_ifi_overall_test_bin10{suf}.csv"
    overall_test = pd.read_csv(test_csv) if test_csv.exists() else None
    if overall_test is None:
        print(f"  ({test_csv.name} missing — panel titles will not carry the ask-3 test)")
    fig_windows(df, fs, overall_test)
    fig_signmap(df, fs)


if __name__ == "__main__":
    main()
