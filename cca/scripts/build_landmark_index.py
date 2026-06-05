"""Build a single browsable HTML contact sheet for the landmark-arm figures.

Scans ``figures/`` for the sweep-summary figures and the per-config figure sets
produced by ``plot_landmark.py`` and ``learning_changes.py``, and writes
``figures/index.html`` -- a self-contained page (relative <img> paths, no
external assets) with a table of contents, the cross-config sweep up top, then
one collapsible section per config with captioned figures.

Only canonical figure filenames are embedded; stale leftovers from older runs
(e.g. ``cc_strength.png``) are ignored. Operates purely on existing pngs.

Run:  python scripts/build_landmark_index.py
"""

from __future__ import annotations

import html
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config  # noqa: E402

FIG_DIR = Path(config.FIGURES_DIR)
SWEEP_DIR = FIG_DIR / "landmark_sweep"
COMMITTED = "landmark50_res_samp15"

# --- captions -------------------------------------------------------------

SWEEP_FIGS = [
    ("sweep_strength_naive.png", "Communication strength (naive)",
     "Median held-out peak canonical correlation across significant subspace "
     "dimensions, naive epoch. Rows = config, cols = HC↔cortex pair."),
    ("sweep_strength_intermediate.png", "Communication strength (intermediate)",
     "Median held-out peak CC across significant dims, intermediate epoch."),
    ("sweep_strength_expert.png", "Communication strength (expert)",
     "Median held-out peak CC across significant dims, expert epoch."),
    ("sweep_n_sig_naive.png", "Subspace dimensionality (naive)",
     "Mean number of significant canonical dimensions per fit cell, naive epoch."),
    ("sweep_n_sig_intermediate.png", "Subspace dimensionality (intermediate)",
     "Mean number of significant dims per cell, intermediate epoch."),
    ("sweep_n_sig_expert.png", "Subspace dimensionality (expert)",
     "Mean number of significant dims per cell, expert epoch."),
    ("sweep_frac_sig_cells_naive.png", "Fraction of significant cells (naive)",
     "Fraction of fit cells with ≥1 significant canonical dimension, naive epoch."),
    ("sweep_ifi_expert.png", "Information-flow index (expert)",
     "Median information-flow index (IFI) across significant dims over the full "
     "lag window, expert epoch. Sign indicates net HC↔cortex direction."),
    ("sweep_lr_change_cc.png", "Learning change: Δ CC",
     "Median Δ peak CC (expert − naive) across configs × pairs — robustness "
     "of the strength learning effect."),
    ("sweep_lr_change_n_sig.png", "Learning change: Δ n_sig",
     "Mean Δ number of significant dims (expert − naive) across configs × pairs."),
    ("sweep_ifi_window_per_pair_expert.png", "IFI vs lag-window (expert)",
     "Per-pair median IFI as a function of lag-window size, expert epoch. "
     "Rows = config, cols = lag window."),
    ("sweep_ifi_window_lr_change.png", "IFI vs lag-window: learning change",
     "Per-pair Δ median IFI (expert − naive) as a function of lag-window size."),
    ("sweep_pvalues_summary.png", "Headline-test p-values",
     "Per-pair p-values (−log₁₀ p) for the headline tests: CC vs 0 (naive, "
     "expert), CC naive-vs-expert; same trio for IFI."),
]

# Per-config: (filename-relative-to-config-dir-template, title, caption).
# {tag} is substituted; NN is the IFI lag-window index.
PERCFG_FIGS = [
    ("coverage_{tag}.png", "Coverage (diagnostic)",
     "Per-(epoch, landmark) mean number of landmark crossings per pair. "
     "Sample-size sanity check before interpreting fits."),
    ("stage2_comm_strength_{tag}.png", "Communication strength",
     "Held-out peak CC across significant dims, box+jitter per epoch per pair. "
     "* = Wilcoxon vs 0 (p<0.05); ANOVA + Tukey across epochs in panel titles."),
    ("stage2_subspace_dim_{tag}.png", "Subspace dimensionality",
     "Number of significant canonical dims per cell, per epoch per pair. "
     "Paired and unpaired naive-vs-expert tests."),
    ("stage2_lag_curves_{tag}.png", "Lag curves",
     "Held-out CC vs lag (ms), mean ± SEM across significant dims; one line per "
     "epoch, per pair. Peak lag indicates HC↔cortex lead/lag."),
    ("stage2_per_landmark_{tag}.png", "Per-landmark detail",
     "Per-pair × per-landmark heatmap of median held-out CC across significant "
     "dims, per epoch — exposes landmark-specific structure."),
    ("stage2_ifi_win01_{tag}.png", "IFI (|lag| ≤ 1 bin)",
     "IFI distribution across significant dims per epoch, restricted to |lag| ≤ 1 "
     "bin. * = Wilcoxon vs 0; ANOVA + Tukey."),
    ("stage2_ifi_win02_{tag}.png", "IFI (|lag| ≤ 2 bins)",
     "IFI across significant dims per epoch, |lag| ≤ 2 bins."),
    ("stage2_ifi_win03_{tag}.png", "IFI (|lag| ≤ 3 bins)",
     "IFI across significant dims per epoch, |lag| ≤ 3 bins."),
    ("stage2_ifi_win04_{tag}.png", "IFI (|lag| ≤ 4 bins)",
     "IFI across significant dims per epoch, |lag| ≤ 4 bins."),
    ("stage2_ifi_win05_{tag}.png", "IFI (|lag| ≤ 5 bins)",
     "IFI across significant dims per epoch, |lag| ≤ 5 bins."),
    ("learning_changes_perlandmark_mi.png", "Learning change: Δ MI per landmark",
     "Per-(pair, landmark) median Δ MI_sig (expert − naive) with FDR significance "
     "markers."),
    ("learning_changes_perlandmark_ifi.png", "Learning change: Δ IFI per landmark",
     "Per-(pair, landmark) median Δ CC-weighted IFI (expert − naive) with FDR "
     "markers."),
    ("learning_changes_pooled.png", "Learning change: pooled per pair",
     "Per-pair learning-change bars (expert − naive), pooled over landmarks."),
]

KRULE_DESC = {
    "fix03": "fixed k=3", "fix05": "fixed k=5", "fix10": "fixed k=10",
    "fix20": "fixed k=20", "fix30": "fixed k=30",
    "samp15": "sampled, 15 units", "samp25": "sampled, 25 units",
    "samp40": "sampled, 40 units",
    "var75": "variance 75%", "var85": "variance 85%", "var95": "variance 95%",
}


def parse_tag(tag: str):
    """landmark{bin}_{cca}_{krule} -> (bin_ms:int, cca:str, krule:str)."""
    m = re.match(r"^landmark(\d+)_(res|sig)_(\w+)$", tag)
    if not m:
        return (9999, "zzz", tag)
    return (int(m.group(1)), m.group(2), m.group(3))


def config_dirs():
    """All per-config figure dirs, committed default first, then sorted."""
    dirs = [
        p.name for p in FIG_DIR.iterdir()
        if p.is_dir() and re.match(r"^landmark(25|50)_(res|sig)_\w+$", p.name)
    ]
    dirs.sort(key=lambda t: (parse_tag(t)[0], parse_tag(t)[1], t))
    if COMMITTED in dirs:
        dirs.remove(COMMITTED)
        dirs.insert(0, COMMITTED)
    return dirs


def fig_block(rel_src: str, title: str, caption: str) -> str:
    s = html.escape(rel_src)
    return (
        '<figure>\n'
        f'  <a href="{s}" target="_blank"><img loading="lazy" src="{s}" '
        f'alt="{html.escape(title)}"></a>\n'
        f'  <figcaption><b>{html.escape(title)}</b><br>{html.escape(caption)}'
        '</figcaption>\n</figure>\n'
    )


def build():
    if not FIG_DIR.exists():
        sys.exit(f"figures dir not found: {FIG_DIR}")

    parts: list[str] = []
    parts.append(HEAD)

    # --- sweep section ----------------------------------------------------
    parts.append('<section id="sweep"><h2>Cross-config sweep summary</h2>')
    parts.append(
        '<p class="lede">Robustness of the landmark-arm findings across all '
        'preprocessing configs (25/50 ms bins × residual/signal CCA × '
        'subspace-selection rule). Statistical unit = significant canonical '
        'subspace dimension, pooled across animals and landmarks within each '
        'HC↔cortex pair.</p>')
    parts.append('<div class="grid">')
    n_sweep = 0
    for fname, title, cap in SWEEP_FIGS:
        if (SWEEP_DIR / fname).exists():
            parts.append(fig_block(f"landmark_sweep/{fname}", title, cap))
            n_sweep += 1
    parts.append('</div></section>')

    # --- per-config sections ---------------------------------------------
    dirs = config_dirs()
    toc = ['<nav id="toc"><b>Configs:</b> ']
    sections = []
    n_cfg_figs = 0
    for tag in dirs:
        bin_ms, cca, krule = parse_tag(tag)
        cca_full = "residual" if cca == "res" else "signal"
        krule_full = KRULE_DESC.get(krule, krule)
        badge = " ★ committed default" if tag == COMMITTED else ""
        toc.append(f'<a href="#{html.escape(tag)}">{html.escape(tag)}</a> ')

        body = [
            f'<section id="{html.escape(tag)}" class="cfg">',
            f'<h2>{html.escape(tag)}<span class="badge">{html.escape(badge)}'
            '</span></h2>',
            f'<p class="meta">{bin_ms} ms bins &middot; {cca_full} CCA &middot; '
            f'subspace rule: {html.escape(krule_full)}</p>',
            '<div class="grid">',
        ]
        cfg_dir = FIG_DIR / tag
        for tmpl, title, cap in PERCFG_FIGS:
            fname = tmpl.format(tag=tag)
            if (cfg_dir / fname).exists():
                body.append(fig_block(f"{tag}/{fname}", title, cap))
                n_cfg_figs += 1
        body.append('</div><p class="top"><a href="#top">↑ top</a></p>'
                     '</section>')
        sections.append("\n".join(body))
    toc.append('</nav>')

    parts.append("".join(toc))
    parts.extend(sections)
    parts.append(FOOT)

    out = FIG_DIR / "index.html"
    out.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {out}")
    print(f"  sweep figures embedded: {n_sweep}/{len(SWEEP_FIGS)}")
    print(f"  configs: {len(dirs)}; per-config figures embedded: {n_cfg_figs}")


HEAD = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Landmark arm — figure index</title>
<style>
 :root { --bg:#0f1115; --card:#171a21; --fg:#e6e6e6; --mut:#9aa4b2;
         --acc:#6ea8fe; --bord:#262b36; }
 * { box-sizing:border-box; }
 body { margin:0; background:var(--bg); color:var(--fg);
        font:15px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif; }
 header { padding:24px 28px; border-bottom:1px solid var(--bord); }
 h1 { margin:0 0 6px; font-size:22px; }
 h2 { font-size:18px; border-bottom:1px solid var(--bord); padding-bottom:6px; }
 .sub { color:var(--mut); max-width:70ch; }
 section { padding:18px 28px; }
 #toc { padding:14px 28px; border-bottom:1px solid var(--bord);
        line-height:2.1; }
 #toc a { color:var(--acc); text-decoration:none; margin-right:4px;
          font-size:12.5px; white-space:nowrap; }
 #toc a:hover { text-decoration:underline; }
 .lede, .meta { color:var(--mut); max-width:80ch; }
 .meta { margin-top:-4px; font-size:13px; }
 .badge { color:#ffd166; font-weight:600; font-size:13px; margin-left:8px; }
 .grid { display:grid; gap:18px;
         grid-template-columns:repeat(auto-fill,minmax(420px,1fr)); }
 figure { margin:0; background:var(--card); border:1px solid var(--bord);
          border-radius:10px; overflow:hidden; }
 figure img { width:100%; display:block; background:#fff; }
 figcaption { padding:10px 12px; color:var(--mut); font-size:12.5px; }
 figcaption b { color:var(--fg); }
 .cfg { border-top:1px solid var(--bord); }
 .top { margin:10px 0 0; }
 .top a, a.up { color:var(--acc); text-decoration:none; font-size:12.5px; }
</style></head><body><a id="top"></a>
<header>
 <h1>Landmark arm — figure index</h1>
 <p class="sub">Per-landmark HC↔cortex communication-subspace analysis
 (Arm B): 500 ms windows centred on visual-landmark entries, residual CCA,
 10-fold leave-one-trial-out CV, fit per (animal, pair, epoch, landmark).
 Cross-config sweep first, then one section per preprocessing config
 (committed default starred).</p>
</header>
"""

FOOT = "</body></html>"


if __name__ == "__main__":
    build()
