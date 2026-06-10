"""Config-pruning tables for the landmark arm.

Turns the per-(config, pair, epoch) pooled sig-dim values produced by the sweep
aggregation into a flat, sortable table for deciding which preprocessing configs
to keep or discard. Pure functions only -- no I/O, no pickle, no matplotlib --
so the numeric logic is unit-testable against synthetic ground truth.

The statistical unit throughout is the **significant canonical subspace
dimension**, pooled across animals and landmarks within a (config, pair) cell --
the same convention as ``summarise_landmark_sweep.py``.

Key signal for pruning: a *held-out* peak canonical correlation that saturates
near 1.0 is an overfitting artifact (too many canonical dims relative to paired
samples), not real communication. ``cc_diag`` exposes that via ``max_cc`` and
``frac_cc_ge_099`` so over-parameterised configs (large fixed-k / high
variance-threshold rules) can be thrown out on evidence.

Input bucket shape (one per (config, pair, epoch)), as built by
``summarise_landmark_sweep._per_pair_stats``::

    {
      "cc":            [held-out peak CC, one per significant dim],
      "ifi":           [full-window IFI, one per significant dim],
      "ifi_win":       {window_int -> [IFI values across sig dims]},
      "n_sig_cells":   [per-cell count of sig dims],
      "n_cells_total":     int,
      "n_cells_with_sig":  int,
    }
"""

from __future__ import annotations

import re
from typing import Optional

import numpy as np
from scipy import stats

TAG_RE = re.compile(r"^landmark(\d+)_(res|sig)_(\w+)$")

# Epoch order used for column suffixes and the naive->expert contrast.
EPOCHS = ("naive", "intermediate", "expert")

# A (config, pair) cell is flagged as an overfitting suspect if, in ANY epoch,
# the median held-out CC across sig dims reaches this, or at least this fraction
# of sig-dim CCs are essentially perfect.
OVERFIT_MEDIAN_CC = 0.99
OVERFIT_FRAC_GE_099 = 0.5


def parse_tag(tag: str) -> Optional[tuple[int, str, str]]:
    """``landmark{bin}_{cca}_{krule}`` -> ``(bin_ms, cca, krule)`` or None."""
    m = TAG_RE.match(tag)
    if not m:
        return None
    return int(m.group(1)), m.group(2), m.group(3)


def wilcoxon_vs0(v) -> float:
    """Two-sided Wilcoxon signed-rank p-value vs 0; NaN if underpowered.

    Matches ``summarise_landmark_sweep._wilcoxon_vs0`` (>=6 finite, non-all-zero).
    """
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 6 or not np.any(v != 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(v).pvalue)
    except ValueError:
        return float("nan")


def mwu(a, b) -> float:
    """Two-sided Mann-Whitney U p-value; NaN if either group < 3 finite."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return float("nan")
    try:
        return float(stats.mannwhitneyu(a, b).pvalue)
    except ValueError:
        return float("nan")


def cc_diag(cc_list) -> dict:
    """Overfitting diagnostics on a pool of held-out peak-CC values.

    Returns median, max, the fraction >= 0.99 (the overfitting tell), and n.
    Empty / all-non-finite input yields NaNs and n=0.
    """
    cc = np.asarray([v for v in np.asarray(cc_list, float) if np.isfinite(v)],
                    float)
    if cc.size == 0:
        return {"median": float("nan"), "max": float("nan"),
                "frac_ge_099": float("nan"), "n": 0}
    return {"median": float(np.median(cc)), "max": float(cc.max()),
            "frac_ge_099": float(np.mean(cc >= 0.99)), "n": int(cc.size)}


def _max_window_ifi(bucket) -> list:
    """IFI values at the largest available lag window for a bucket."""
    iw = bucket.get("ifi_win") or {}
    if not iw:
        return []
    return iw[max(iw)]


def headline_pvalues(bucket_naive, bucket_expert) -> dict:
    """The six headline tests for one (config, pair).

    CC and IFI(at max |lag| window) each get: naive vs 0, expert vs 0 (Wilcoxon
    signed-rank), and naive vs expert (Mann-Whitney on the pooled sig-dim
    values). Mirrors the tests drawn in ``sweep_pvalues_summary.png`` exactly.
    Missing epoch buckets yield NaN p-values.
    """
    out = {"cc_p_naive_vs0": float("nan"), "cc_p_expert_vs0": float("nan"),
           "cc_p_naive_vs_expert": float("nan"),
           "ifi_p_naive_vs0": float("nan"), "ifi_p_expert_vs0": float("nan"),
           "ifi_p_naive_vs_expert": float("nan")}
    if bucket_naive is None or bucket_expert is None:
        return out
    out["cc_p_naive_vs0"] = wilcoxon_vs0(bucket_naive["cc"])
    out["cc_p_expert_vs0"] = wilcoxon_vs0(bucket_expert["cc"])
    out["cc_p_naive_vs_expert"] = mwu(bucket_naive["cc"], bucket_expert["cc"])

    iw_n, iw_e = bucket_naive.get("ifi_win") or {}, bucket_expert.get("ifi_win") or {}
    max_w_n = max(iw_n) if iw_n else None
    max_w_e = max(iw_e) if iw_e else None
    if max_w_n is not None:
        out["ifi_p_naive_vs0"] = wilcoxon_vs0(iw_n[max_w_n])
    if max_w_e is not None:
        out["ifi_p_expert_vs0"] = wilcoxon_vs0(iw_e[max_w_e])
    if max_w_n is not None and max_w_n == max_w_e:
        out["ifi_p_naive_vs_expert"] = mwu(iw_n[max_w_n], iw_e[max_w_e])
    return out


def _epoch_block(bucket) -> dict:
    """Per-epoch metric block for one (config, pair, epoch) bucket."""
    if bucket is None:
        bucket = {}
    n_total = int(bucket.get("n_cells_total", 0))
    n_sig_cells = np.asarray(bucket.get("n_sig_cells", []), float)
    ifi = np.asarray([v for v in np.asarray(bucket.get("ifi", []), float)
                      if np.isfinite(v)], float)
    d = cc_diag(bucket.get("cc", []))
    return {
        "n_cells": n_total,
        "frac_sig": (bucket.get("n_cells_with_sig", 0) / n_total
                     if n_total > 0 else float("nan")),
        "mean_n_sig": float(np.mean(n_sig_cells)) if n_sig_cells.size
        else float("nan"),
        "median_cc": d["median"],
        "max_cc": d["max"],
        "frac_cc_ge_099": d["frac_ge_099"],
        "median_ifi": float(np.median(ifi)) if ifi.size else float("nan"),
        "n_sig_dims": d["n"],
    }


# Column order for the per-(config, pair) detailed table.
def prune_fieldnames() -> list[str]:
    base = ["tag", "bin_ms", "cca", "krule", "pair"]
    per_epoch = ["n_cells", "frac_sig", "mean_n_sig", "median_cc", "max_cc",
                 "frac_cc_ge_099", "median_ifi", "n_sig_dims"]
    epoch_cols = [f"{m}_{e}" for e in EPOCHS for m in per_epoch]
    deltas = ["delta_cc_expert_naive", "delta_n_sig_expert_naive"]
    pvals = ["cc_p_naive_vs0", "cc_p_expert_vs0", "cc_p_naive_vs_expert",
             "ifi_p_naive_vs0", "ifi_p_expert_vs0", "ifi_p_naive_vs_expert"]
    return base + epoch_cols + deltas + pvals + ["overfit_flag"]


def build_prune_rows(per_cpe: dict) -> list[dict]:
    """One row per (config, pair); wide over epochs, with p-values + flags.

    ``per_cpe`` maps ``(tag, pair, epoch) -> bucket``. ``pair`` may be a tuple
    ``(area_x, area_y)`` or a string; it is rendered as ``"x-y"``.
    """
    cps = sorted({(tag, pair) for (tag, pair, _e) in per_cpe})
    rows: list[dict] = []
    for tag, pair in cps:
        pt = parse_tag(tag)
        if pt is None:
            continue
        bin_ms, cca, krule = pt
        pair_str = f"{pair[0]}-{pair[1]}" if isinstance(pair, tuple) else str(pair)
        row = {"tag": tag, "bin_ms": bin_ms, "cca": cca, "krule": krule,
               "pair": pair_str}
        blocks = {}
        for e in EPOCHS:
            blk = _epoch_block(per_cpe.get((tag, pair, e)))
            blocks[e] = blk
            for m, v in blk.items():
                row[f"{m}_{e}"] = v
        row["delta_cc_expert_naive"] = (blocks["expert"]["median_cc"]
                                        - blocks["naive"]["median_cc"])
        row["delta_n_sig_expert_naive"] = (blocks["expert"]["mean_n_sig"]
                                           - blocks["naive"]["mean_n_sig"])
        row.update(headline_pvalues(per_cpe.get((tag, pair, "naive")),
                                    per_cpe.get((tag, pair, "expert"))))
        row["overfit_flag"] = any(
            (np.isfinite(blocks[e]["median_cc"])
             and blocks[e]["median_cc"] >= OVERFIT_MEDIAN_CC)
            or (np.isfinite(blocks[e]["frac_cc_ge_099"])
                and blocks[e]["frac_cc_ge_099"] >= OVERFIT_FRAC_GE_099)
            for e in EPOCHS)
        rows.append(row)
    return rows


def rollup_fieldnames() -> list[str]:
    return ["tag", "bin_ms", "cca", "krule", "n_pairs",
            "max_cc_any", "n_overfit_pairs", "frac_overfit_pairs",
            "median_cc_expert_med", "mean_n_sig_expert_med",
            "n_pairs_cc_expert_sig", "n_pairs_cc_learn_up", "alpha"]


def rollup_rows(prune_rows: list[dict], alpha: float = 0.05) -> list[dict]:
    """Collapse the per-(config, pair) table to one row per config.

    Surfaces the headline pruning signals: how many pairs overfit, the worst
    held-out CC anywhere, typical expert strength/dimensionality, and how many
    pairs show a significant expert CC and a significant naive->expert increase.
    """
    by_tag: dict[str, list[dict]] = {}
    for r in prune_rows:
        by_tag.setdefault(r["tag"], []).append(r)
    out: list[dict] = []
    for tag in sorted(by_tag):
        rs = by_tag[tag]
        pt = parse_tag(tag)
        bin_ms, cca, krule = pt if pt else (float("nan"), "", "")
        max_cc_vals = [r[f"max_cc_{e}"] for r in rs for e in EPOCHS
                       if np.isfinite(r[f"max_cc_{e}"])]
        med_cc_exp = [r["median_cc_expert"] for r in rs
                      if np.isfinite(r["median_cc_expert"])]
        n_sig_exp = [r["mean_n_sig_expert"] for r in rs
                     if np.isfinite(r["mean_n_sig_expert"])]
        n_overfit = sum(1 for r in rs if r["overfit_flag"])
        n_cc_sig = sum(1 for r in rs
                       if np.isfinite(r["cc_p_expert_vs0"])
                       and r["cc_p_expert_vs0"] < alpha)
        n_learn_up = sum(1 for r in rs
                         if np.isfinite(r["cc_p_naive_vs_expert"])
                         and r["cc_p_naive_vs_expert"] < alpha
                         and np.isfinite(r["delta_cc_expert_naive"])
                         and r["delta_cc_expert_naive"] > 0)
        out.append({
            "tag": tag, "bin_ms": bin_ms, "cca": cca, "krule": krule,
            "n_pairs": len(rs),
            "max_cc_any": max(max_cc_vals) if max_cc_vals else float("nan"),
            "n_overfit_pairs": n_overfit,
            "frac_overfit_pairs": n_overfit / len(rs) if rs else float("nan"),
            "median_cc_expert_med": float(np.median(med_cc_exp))
            if med_cc_exp else float("nan"),
            "mean_n_sig_expert_med": float(np.median(n_sig_exp))
            if n_sig_exp else float("nan"),
            "n_pairs_cc_expert_sig": n_cc_sig,
            "n_pairs_cc_learn_up": n_learn_up,
            "alpha": alpha,
        })
    return out
