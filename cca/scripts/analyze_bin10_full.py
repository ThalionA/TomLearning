"""Comprehensive pair-by-pair 10 ms results — every metric × every comparison × both FS,
emitted as Markdown tables for the report (exact numbers, no shortcuts).

Sections:
  A. Levels (hierarchy): CC1, n_sig, MI_sig, Gini, mean IFI(±50ms) per pair.
  B. Directionality — existence at ±50 ms: who leads, IFI, t, p(t), n.
  C. Directionality — integration curves to ±250 ms: IFI by window.
  D. Directionality — per canonical dimension (lower CCs): IFI for dim0..3.
  E. Change vs learning — slopes: CC1, Gini, IFI vs trial_frac / performance / lp_rel.
  F. Epoch contrasts: naive→expert and naive→intermediate for CC1, Gini, IFI, n_sig.
  G. Reorientation: CC1 and top-3 cross-window rotation vs split-half floor.
  H. Transition (uncued→cued): Δ for CC1, Gini, IFI, n_sig.

Writes results/bin10_tables.md.  Run: PYTHONPATH=src python scripts/analyze_bin10_full.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figs_report as F  # noqa: E402
from tom_cca import config, paired_stats  # noqa: E402

PAIRS = F.PAIRS
R = config.RESULTS_DIR
OUT = R / "bin10_tables.md"
FS = [("FS-excluded", ""), ("FS-included", "_fsincl")]
WINS = [10, 30, 50, 100, 150, 200, 250]
out_lines: list[str] = []


def w(s=""):
    out_lines.append(s)


def _num(s):
    return pd.to_numeric(s, errors="coerce")


def ptt(v):
    """Paired t-test of per-animal deltas vs 0 → (median, p, n)."""
    v = [x for x in v if np.isfinite(x)]
    if len(v) < 3:
        return np.nan, np.nan, len(v)
    _, med, _, p = paired_stats.paired_t(v)
    return med, p, len(v)


def star(p):
    return "**" if (np.isfinite(p) and p < 0.05) else ""


def lead(pair, ifi):
    a, b = pair.split("-")
    return f"{a}→{b}" if ifi > 0 else f"{b}→{a}"


# ---------------------------------------------------------------- A. levels
def sec_levels():
    w("## A. Levels — the coupling hierarchy (10 ms)\n")
    for fslab, suf in FS:
        tr = pd.read_csv(R / f"trajectory_w15_bin10{suf}.csv")
        win = pd.read_csv(R / f"ifi_windows_bin10{suf}.csv")
        w50 = win[win.window_ms == 50]
        w(f"**{fslab}** (per-animal mean, pooled over 16 animals; IFI at ±50 ms):\n")
        w("| pair | CC₁ | n_sig | MI_sig | Gini_x | mean IFI (±50 ms) |")
        w("|---|---|---|---|---|---|")
        # order by CC1 desc
        rows = []
        for p in PAIRS:
            cc = np.nanmean(F.per_animal_level(tr, p, "cc1"))
            if not np.isfinite(cc):
                continue
            ns = np.nanmean(F.per_animal_level(tr, p, "n_sig"))
            mi = np.nanmean(F.per_animal_level(tr, p, "mi_sig"))
            gi = np.nanmean(F.per_animal_level(tr, p, "gini_x"))
            ifi = _num(w50[w50.pair == p]["ifi"]).mean()
            rows.append((cc, p, ns, mi, gi, ifi))
        for cc, p, ns, mi, gi, ifi in sorted(rows, reverse=True):
            w(f"| {p} | {cc:.3f} | {ns:.2f} | {mi:.3f} | {gi:.3f} | {ifi:+.3f} |")
        w()


# ------------------------------------------------ B. directionality existence
def sec_existence():
    w("## B. Directionality — existence at ±50 ms (held-out IFI window sweep, 10 ms)\n")
    w("IFI > 0 ⇒ first-named area leads. *t* = one-sample (paired) t vs 0 across animals.\n")
    for fslab, suf in FS:
        win = pd.read_csv(R / f"ifi_windows_bin10{suf}.csv")
        w50 = win[win.window_ms == 50]
        w(f"**{fslab}** ({win.animal.nunique()} animals):\n")
        w("| pair | flow | mean IFI | t | p (t) | n |")
        w("|---|---|---|---|---|---|")
        for p in PAIRS:
            v = _num(w50[w50.pair == p]["ifi"]).dropna().to_numpy()
            if v.size < 3:
                continue
            t, pt = st.ttest_1samp(v, 0)
            w(f"| {p} | {lead(p, v.mean())} | {v.mean():+.3f} | {t:+.2f} | "
              f"{pt:.3f}{star(pt)} | {v.size} |")
        w()


# ------------------------------------------------ C. directionality curves
def sec_curves():
    w("## C. Directionality — integration curves to ±250 ms (10 ms)\n")
    w("Mean IFI across animals at each integration window (`**` = t-test p<0.05). "
      "Shows whether each flow is tight-lag (peaks ≤±50 ms) or diffuse.\n")
    for fslab, suf in FS:
        win = pd.read_csv(R / f"ifi_windows_bin10{suf}.csv")
        w(f"**{fslab}**:\n")
        w("| pair | " + " | ".join(f"±{x} ms" for x in WINS) + " |")
        w("|---|" + "---|" * len(WINS))
        for p in PAIRS:
            g = win[win.pair == p]
            cells = []
            ok = False
            for x in WINS:
                v = _num(g[g.window_ms == x]["ifi"]).dropna().to_numpy()
                if v.size >= 3 and v.std() > 0:
                    t, pv = st.ttest_1samp(v, 0)
                    cells.append(f"{v.mean():+.3f}{'**' if pv < 0.05 else ''}")
                    ok = True
                else:
                    cells.append("—")
            if ok:
                w(f"| {p} | " + " | ".join(cells) + " |")
        w()


# ------------------------------------------------ D. per-dim IFI (lower CCs)
def sec_perdim():
    w("## D. Directionality — by canonical dimension (lower CCs, 10 ms, ±50 ms)\n")
    w("dim 0 = dominant CC₁ (the headline). `**` = pooled-dim t-test p<0.05 "
      "(dims-as-n, anti-conservative). Mean CC per dim in the last column.\n")
    for fslab, suf in FS:
        d = pd.read_csv(R / f"trajectory_w15_bin10_dims{suf}.csv")
        w(f"**{fslab}**:\n")
        w("| pair | dim0 (CC₁) | dim1 | dim2 | dim3 | mean CC dim0/1/2 |")
        w("|---|---|---|---|---|---|")
        for p in PAIRS:
            g = d[d.pair == p]
            cells, ccs = [], []
            ok = False
            for dim in [0, 1, 2, 3]:
                v = _num(g[g.dim == dim]["ifi"]).dropna().to_numpy()
                ccs.append(_num(g[g.dim == dim]["cc"]).dropna().mean())
                if v.size >= 20:
                    t, pv = st.ttest_1samp(v, 0)
                    cells.append(f"{v.mean():+.3f}{'**' if pv < 0.05 else ''}")
                    ok = True
                else:
                    cells.append("—")
            if ok:
                w(f"| {p} | {cells[0]} | {cells[1]} | {cells[2]} | {cells[3]} | "
                  f"{ccs[0]:.2f}/{ccs[1]:.2f}/{ccs[2]:.2f} |")
        w()


# ------------------------------------------------ E. change-slopes
def sec_slopes():
    w("## E. Change vs learning — per-animal slopes (10 ms, paired t on per-animal slopes)\n")
    w("Slope of each metric vs the learning axis; `**` = paired t p<0.05.\n")
    for metric, mlab in [("cc1", "CC₁ (strength)"), ("gini_x", "Gini (sparsity)"),
                         ("ifi", "IFI (direction, ±50 ms)")]:
        w(f"### {mlab}\n")
        for fslab, suf in FS:
            tr = pd.read_csv(R / f"trajectory_w15_bin10{suf}.csv")
            w(f"**{fslab}**:\n")
            w("| pair | vs trial-fraction | vs performance | vs LP-relative | n |")
            w("|---|---|---|---|---|")
            for p in PAIRS:
                cells = []
                n = 0
                for axis in ["trial_frac", "performance", "lp_rel"]:
                    s, pv, nn = F._slope_quant(tr, metric, axis, p)
                    n = nn
                    cells.append(f"{s:+.4f} (p={pv:.3f}){star(pv)}" if np.isfinite(s) else "—")
                if any(c != "—" for c in cells):
                    w(f"| {p} | {cells[0]} | {cells[1]} | {cells[2]} | {n} |")
            w()


# ------------------------------------------------ F. epoch contrasts
def sec_epochs():
    w("## F. Epoch contrasts — naive→expert and naive→intermediate (10 ms)\n")
    w("Per-animal paired Δ; *t* = paired t vs 0. Learners only "
      "(epochs are defined relative to the learning point).\n")
    for metric, mlab in [("cc1", "CC₁"), ("gini_x", "Gini"), ("ifi", "IFI"), ("n_sig", "n_sig")]:
        w(f"### {mlab}\n")
        for fslab, suf in FS:
            ep = pd.read_csv(R / f"epoch_metrics_bin10{suf}.csv")
            w(f"**{fslab}**:\n")
            w("| pair | naive | int | expert | Δ(exp−nai) t-p | Δ(int−nai) t-p | n |")
            w("|---|---|---|---|---|---|---|")
            for p in PAIRS:
                g = ep[ep.pair == p]
                per = {}
                for an, gg in g.groupby("animal"):
                    per[an] = {e: _num(gg[gg.epoch == e][metric]).mean()
                               for e in ["naive", "intermediate", "expert"]}
                nv = [v["naive"] for v in per.values() if np.isfinite(v["naive"])]
                it = [v["intermediate"] for v in per.values() if np.isfinite(v["intermediate"])]
                xp = [v["expert"] for v in per.values() if np.isfinite(v["expert"])]
                # paired exp-nai
                dxe = [v["expert"] - v["naive"] for v in per.values()
                       if np.isfinite(v["expert"]) and np.isfinite(v["naive"])]
                dxi = [v["intermediate"] - v["naive"] for v in per.values()
                       if np.isfinite(v["intermediate"]) and np.isfinite(v["naive"])]
                if len(dxe) < 3:
                    continue
                te = st.ttest_1samp(dxe, 0)[1] if len(dxe) >= 2 else np.nan
                ti = st.ttest_1samp(dxi, 0)[1] if len(dxi) >= 2 else np.nan
                w(f"| {p} | {np.mean(nv):+.3f} | {np.mean(it) if it else float('nan'):+.3f} | "
                  f"{np.mean(xp):+.3f} | {te:.3f}{star(te)} | "
                  f"{ti:.3f}{star(ti)} | {len(dxe)} |")
            w()


# ------------------------------------------------ G. reorientation
def sec_rotation():
    w("## G. Reorientation — cross-window rotation vs split-half floor (10 ms)\n")
    w("Rotation > floor ⇒ genuine reorientation. `**` = rotation significantly "
      "*below* floor (paired t), the no-reorientation signature.\n")
    for fslab, suf in FS:
        tr = pd.read_csv(R / f"trajectory_w15_bin10{suf}.csv")
        w(f"**{fslab}** (degrees):\n")
        w("| pair | CC₁ rot | CC₁ floor | p(rot<floor) | top-3 rot | top-3 floor | p |")
        w("|---|---|---|---|---|---|---|")
        for p in PAIRS:
            g = tr[tr.pair == p]
            res = {}
            for rc, sc, key in [("rot_x_cc1", "sh_x_cc1", "cc1"), ("rot_x", "sh_x", "t3")]:
                rs, fs, dz = [], [], []
                for an, gg in g.groupby("animal"):
                    r = _num(gg[rc]).mean(); h = _num(gg[sc]).mean()
                    if np.isfinite(r):
                        rs.append(r)
                    if np.isfinite(h):
                        fs.append(h)
                    if np.isfinite(r) and np.isfinite(h):
                        dz.append(r - h)
                _, pv, _ = ptt(dz)
                res[key] = (np.nanmean(rs) if rs else np.nan,
                            np.nanmean(fs) if fs else np.nan, pv)
            if not np.isfinite(res["cc1"][0]):
                continue
            c, t3 = res["cc1"], res["t3"]
            w(f"| {p} | {c[0]:.1f} | {c[1]:.1f} | {c[2]:.3f}{star(c[2])} | "
              f"{t3[0]:.1f} | {t3[1]:.1f} | {t3[2]:.3f}{star(t3[2])} |")
        w()


# ------------------------------------------------ H. transition
def sec_transition():
    w("## H. Transition — uncued→cued (10 ms, Δ = cued − uncued)\n")
    w("Per-animal paired Δ, paired t vs 0. `**` = p<0.05.\n")
    for fslab, suf in FS:
        path = R / f"transition_uncued_cued_bin10{suf}.csv"
        if not path.is_file() or path.stat().st_size < 50:
            w(f"**{fslab}**: (no rows yet — re-run pending)\n")
            continue
        t = pd.read_csv(path)
        if len(t) == 0:
            w(f"**{fslab}**: (0 rows)\n")
            continue
        w(f"**{fslab}** ({t.animal.nunique()} animals):\n")
        w("| pair | ΔCC₁ | ΔGini | ΔIFI | Δn_sig | n |")
        w("|---|---|---|---|---|---|")
        for p in PAIRS:
            g = t[t.pair == p]
            if len(g) < 3:
                continue
            cells = []
            for m in ["d_cc1", "d_gini_x", "d_ifi", "d_n_sig"]:
                med, pv, n = ptt(_num(g[m]).dropna().tolist())
                cells.append(f"{med:+.3f} (p={pv:.3f}){star(pv)}" if np.isfinite(med) else "—")
            w(f"| {p} | {cells[0]} | {cells[1]} | {cells[2]} | {cells[3]} | {g.animal.nunique()} |")
        w()


def main():
    w("# 10 ms re-run — comprehensive pair-by-pair results\n")
    w("All values are held-out cross-validated. Bin = 10 ms; headline IFI integrated "
      "over **±50 ms** (curves shown to ±250 ms, §C). pCCA controls all other recorded "
      "areas. `**` flags p<0.05 (uncorrected; per-pair family).\n")
    for fn in (sec_levels, sec_existence, sec_curves, sec_perdim,
               sec_slopes, sec_epochs, sec_rotation, sec_transition):
        fn()
    OUT.write_text("\n".join(out_lines) + "\n")
    print(f"wrote {OUT}  ({len(out_lines)} lines)")


if __name__ == "__main__":
    main()
