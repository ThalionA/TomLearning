"""Item 2 figures — FF/FB CCs labelled once, followed across learning.

Reads results/cc_label_track_epoch_bin10{,_fsincl}.csv (analyze_cc_label_track.py).

Figure 1  HCV1_cc_label_track_<fs>_bin10
          Per pair, the IFI of the FF-labelled and FB-labelled CCs across
          naive/intermediate/expert, mean +/- SEM across animals. Each title carries the
          three terms of the 2-way RM-ANOVA (epoch, FF/FB, interaction) rather than the
          old dFF - dFB contrast, which discarded the intermediate epoch and could not
          separate a main effect from an interaction. The label term on SIGNED IFI is
          circular -- components are labelled BY the sign of their whole-session IFI --
          and is marked as such.

Figure 3  HCV1_cc_label_anova_<fs>_bin10
          The ANOVA itself: one panel per dependent variable, rows = pair, columns = the
          three terms, cells coloured by -log10(p) and annotated with p. The circular
          label column for signed IFI is hatched out. Reading it: the EPOCH column
          carries the effect (peak r in 6/7 pairs), the interaction column is empty.

Figure 2  HCV1_cc_label_persistence_<fs>_bin10
          Does a CC keep its direction? Per animal, the fraction of epochs whose own IFI
          sign matches the LEAVE-EPOCH-OUT label. Chance is a genuine 0.5 only because
          the label excludes the epoch being scored — against the whole-session label the
          construction floor is ~0.60.

Usage: PYTHONPATH=src python scripts/figs_cc_label_track.py [fsincl]
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle  # noqa: E402

figstyle.apply()
ATT = Path.home() / "Documents" / "ResearchVault" / "attachments"
RES = Path(__file__).resolve().parents[1] / "results"
PAIRS = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG", "CA1-SUB",
         "RSC-SUB", "V1-RSC"]
EPOCHS = ["naive", "intermediate", "expert"]
EX = {"naive": 0, "intermediate": 1, "expert": 2}
C_FF, C_FB = "#c0392b", "#2c6fbb"


def fig_track(tab: pd.DataFrame, anova: pd.DataFrame, fs: str):
    fig, axes = figstyle.grid(len(PAIRS), ncols=4)
    for ax, pair in zip(axes, PAIRS):
        sub = tab[tab["pair"] == pair]
        if sub.empty:
            ax.set_title(f"{pair}\n(no data)", fontsize=9); ax.axis("off"); continue
        x_lead = pair.split("-")[0]
        for label, colour in [("FF", C_FF), ("FB", C_FB)]:
            g = sub[sub["label"] == label]
            if g.empty:
                continue
            # one value per animal per epoch first — an animal contributes many CCs
            per = g.pivot_table(index="animal", columns="epoch", values="ifi",
                                aggfunc="mean")
            xs, ms, es = [], [], []
            for e in EPOCHS:
                if e not in per.columns:
                    continue
                v = per[e].to_numpy(float); v = v[np.isfinite(v)]
                if v.size == 0:
                    continue
                xs.append(EX[e]); ms.append(v.mean())
                es.append(v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else np.nan)
            if xs:
                ax.errorbar(xs, ms, yerr=es, fmt="-o", color=colour, lw=2.0, ms=5,
                            capsize=3, zorder=3,
                            label=f"{label} · {g['animal'].nunique()} animals")
        # 2-way RM-ANOVA terms for the plotted DV (signed IFI), not the old
        # difference-of-differences: that discarded the intermediate epoch and could
        # not separate a main effect from an interaction.
        row = anova[(anova["pair"] == pair) & (anova["dv"] == "ifi")]
        def _p(term):
            if row.empty or f"p_{term}" not in row:
                return np.nan
            return float(row[f"p_{term}"].iloc[0])
        pe, pl, pi = _p("epoch"), _p("label"), _p("epoch:label")
        ax.axhline(0.0, color="k", lw=0.9, ls="--", zorder=1)
        ax.set_xlim(-0.3, 2.3); ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["naive", "int", "exp"], fontsize=8)
        # the LABEL term is circular for signed IFI — components are labelled BY the
        # sign of their whole-session IFI — so it is shown struck through, not as a result
        ax.set_title(f"{pair}\nepoch p={pe:.2g} · label p={pl:.2g}(circ) · "
                     f"int p={pi:.2g}", fontsize=8)
        ax.set_ylabel(f"IFI   +: {x_lead} leads", fontsize=8)
        ax.legend(fontsize=6.5, loc="best", framealpha=0.6)
    fig.suptitle(f"FF/FB CCs labelled once, followed across learning — {fs} | frozen "
                 f"axes | titles carry the 2-way RM-ANOVA (epoch × label); the label "
                 f"term on signed IFI is circular by construction", fontsize=10)
    figstyle.save(fig, ATT / f"HCV1_cc_label_track_{fs}_bin10.png")
    print(f"wrote HCV1_cc_label_track_{fs}_bin10.png")


def fig_persistence(tab: pd.DataFrame, fs: str):
    per_cc = []
    for (animal, pair, dim), g in tab.groupby(["animal", "pair", "dim"]):
        gg = g[g["label_loo"].isin(["FF", "FB"])]
        if gg.empty:
            continue
        want = np.where(gg["label_loo"].to_numpy() == "FF", 1, -1)
        v = gg["ifi"].to_numpy(float)
        ok = np.isfinite(v) & (v != 0)
        if ok.any():
            per_cc.append({"animal": animal, "pair": pair,
                           "match": float(np.mean(np.sign(v[ok]) == want[ok]))})
    if not per_cc:
        print("  no persistence data"); return
    df = pd.DataFrame(per_cc)
    by_animal = df.groupby("animal")["match"].mean().sort_values()
    by_pair = df.groupby("pair")["match"].mean()

    fig, axes = figstyle.grid(2, ncols=2, panel=(6.4, 4.0))
    ax = axes[0]
    ax.barh(range(len(by_animal)), by_animal.to_numpy(), color="#7f8c8d")
    ax.axvline(0.5, color="#c0392b", lw=1.6, ls="--",
               label="chance (leave-epoch-out label)")
    ax.set_yticks(range(len(by_animal)))
    ax.set_yticklabels(by_animal.index, fontsize=7)
    ax.set_xlabel("fraction of epochs matching the CC's label", fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_title(f"per animal (n={len(by_animal)}, mean {by_animal.mean():.2f})",
                 fontsize=10)
    ax.legend(fontsize=7, loc="lower right", framealpha=0.6)

    ax = axes[1]
    order = by_pair.reindex([p for p in PAIRS if p in by_pair.index])
    ax.barh(range(len(order)), order.to_numpy(), color="#95a5a6")
    ax.axvline(0.5, color="#c0392b", lw=1.6, ls="--")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order.index, fontsize=8)
    ax.set_xlabel("fraction of epochs matching the CC's label", fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_title("per pair — is it carried by the theta-ringing\n"
                 "intra-hippocampal pairs, or general?", fontsize=10)
    fig.suptitle(f"Does a CC keep the direction it was labelled with? — {fs} | label "
                 f"recomputed WITHOUT the epoch being scored, so 0.5 is genuine chance",
                 fontsize=11)
    figstyle.save(fig, ATT / f"HCV1_cc_label_persistence_{fs}_bin10.png")
    print(f"wrote HCV1_cc_label_persistence_{fs}_bin10.png")


def fig_anova(anova: pd.DataFrame, fs: str):
    """All three ANOVA terms, per DV and pair, as -log10(p).

    One panel per dependent variable; rows = area pair, columns = the three terms.
    Cell text is the p-value; the dashed contour marks p = 0.05. The LABEL column for
    signed IFI is hatched and greyed because that term is circular by construction —
    components are labelled BY the sign of their whole-session IFI — so it carries no
    information and must not be read as a result.
    """
    from matplotlib.patches import Rectangle
    dvs = [("peak_r", "peak r", False), ("abs_ifi", "|IFI|", False),
           ("peak_lag_ms", "peak lag", False), ("ifi", "IFI (signed)", True)]
    terms = [("epoch", "epoch"), ("label", "FF/FB"), ("epoch:label", "interaction")]
    fig, axes = figstyle.grid(len(dvs), ncols=4, panel=(4.6, 4.4))
    im = None
    for ax, (dv, name, circ) in zip(axes, dvs):
        sub = anova[anova["dv"] == dv]
        pairs = [p for p in PAIRS if p in set(sub["pair"])]
        if not pairs:
            ax.set_title(f"{name}\n(no data)", fontsize=9); ax.axis("off"); continue
        M = np.full((len(pairs), len(terms)), np.nan)
        for i, pr in enumerate(pairs):
            row = sub[sub["pair"] == pr]
            for j, (term, _) in enumerate(terms):
                col = f"p_{term}"
                if not row.empty and col in row:
                    v = row[col].iloc[0]
                    if np.isfinite(v):
                        M[i, j] = -np.log10(max(v, 1e-12))
        im = ax.imshow(np.ma.masked_invalid(M), cmap="YlOrRd", vmin=0, vmax=4,
                       aspect="auto")
        for i in range(len(pairs)):
            for j in range(len(terms)):
                if not np.isfinite(M[i, j]):
                    continue
                pv = 10 ** (-M[i, j])
                ax.text(j, i, f"{pv:.2g}", ha="center", va="center", fontsize=6,
                        color="white" if M[i, j] > 2.2 else "black",
                        fontweight="bold" if pv < 0.05 else "normal")
        if circ:                                   # grey out the circular column
            ax.add_patch(Rectangle((0.5, -0.5), 1, len(pairs), facecolor="#dcdcdc",
                                   alpha=0.72, hatch="///", edgecolor="#888",
                                   lw=0.8, zorder=5))
            ax.text(1, len(pairs) / 2 - 0.5, "circular", rotation=90, ha="center",
                    va="center", fontsize=7, color="#444", zorder=6)
        ax.set_xticks(range(len(terms)))
        ax.set_xticklabels([t[1] for t in terms], fontsize=8)
        ax.set_yticks(range(len(pairs)))
        ax.set_yticklabels(pairs, fontsize=7)
        ax.set_title(name, fontsize=10)
    if im is not None:
        cb = fig.colorbar(im, ax=list(axes), fraction=0.014, pad=0.01)
        cb.set_label("−log₁₀(p)   (1.3 = 0.05)", fontsize=8)
    fig.suptitle(f"2-way RM-ANOVA per pair: epoch × FF/FB label — {fs} | bold = "
                 f"p < 0.05 | the epoch column is where the effect lives; the "
                 f"interaction is empty", fontsize=10.5)
    figstyle.save(fig, ATT / f"HCV1_cc_label_anova_{fs}_bin10.png")
    print(f"wrote HCV1_cc_label_anova_{fs}_bin10.png")


def main():
    fsincl = "fsincl" in sys.argv[1:]
    suf = "_fsincl" if fsincl else ""
    fs = "fsincl" if fsincl else "fsexcl"
    src = RES / f"cc_label_track_epoch_bin10{suf}.csv"
    st = RES / f"cc_label_anova_bin10{suf}.csv"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/analyze_cc_label_track.py first")
    tab = pd.read_csv(src)
    anova = pd.read_csv(st) if st.exists() else pd.DataFrame()
    fig_track(tab, anova, fs)
    fig_persistence(tab, fs)
    if not anova.empty:
        fig_anova(anova, fs)


if __name__ == "__main__":
    main()
