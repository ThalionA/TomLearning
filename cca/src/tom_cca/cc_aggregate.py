"""Average across an animal's canonical components FIRST, then across animals.

The per-CC tables (`cc_ifi_windows_*`, `cc_label_track_*`) have one row per
(animal, pair, dim, ...). Any statistic that treats those rows as independent is
dims-as-n — dims are nested in animals, so it is the most pseudoreplicated unit
available (STATE.md §3.0). These two helpers make the honest unit the easy one:

* :func:`per_animal_mean` — collapse an animal's CCs to one value per
  (animal, pair, *by), optionally weighted by each CC's held-out strength;
* :func:`one_sample_by_pair` — the across-animal one-sample *t* vs 0, per pair.

They are deliberately dumb pandas so both the analysis scripts and the figure
scripts call the same code and cannot drift apart.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from tom_cca import paired_stats


def per_animal_mean(tab: pd.DataFrame, value: str, by: Sequence[str],
                    weight: str | None = None, sig_only: bool = True,
                    drop_degenerate: bool = True) -> pd.DataFrame:
    """One row per (animal, pair, *by): the mean of ``value`` over that animal's CCs.

    ``sig_only`` keeps rows with ``sig == 1`` and **raises** if the column is missing
    — a table that was already gated upstream must say so with ``sig_only=False``,
    never silently pass every rank through. ``drop_degenerate`` drops rows with
    ``degenerate == 1`` (present only in the IFI tables — ignored if the column is
    absent). Non-finite ``value`` rows are dropped; a cell left with no finite value
    is omitted rather than returned as NaN.

    Columns: ``animal, pair, *by, mean, n_ccs`` and, when ``weight`` is given,
    ``wmean`` — the mean weighted by ``max(weight, 0)`` (a negative held-out CC is
    "no coupling", not negative evidence; a NaN weight also counts as 0, so that CC
    still counts in ``n_ccs`` but not in ``wmean``). ``wmean`` is NaN when every
    weight in the cell clips to zero.
    """
    keys = ["animal", "pair", *by]
    sub = tab.copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    if sig_only:
        if "sig" not in sub.columns:
            raise KeyError("per_animal_mean(sig_only=True) needs a `sig` column; pass "
                           "sig_only=False for a table that was gated upstream")
        sub = sub[pd.to_numeric(sub["sig"], errors="coerce") == 1]
    if drop_degenerate and "degenerate" in sub.columns:
        sub = sub[pd.to_numeric(sub["degenerate"], errors="coerce") == 0]
    sub = sub[np.isfinite(sub[value])]
    if sub.empty:
        cols = keys + ["mean", "n_ccs"] + (["wmean"] if weight else [])
        return pd.DataFrame(columns=cols)
    if weight is not None:
        w = pd.to_numeric(sub[weight], errors="coerce").fillna(0.0).clip(lower=0.0)
        sub = sub.assign(_w=w.to_numpy(), _wv=(w * sub[value]).to_numpy())
    g = sub.groupby(keys, sort=True)
    out = g[value].agg(mean="mean", n_ccs="size").reset_index()
    if weight is not None:
        ws = g["_w"].sum()
        wv = g["_wv"].sum()
        wmean = (wv / ws.where(ws > 0)).rename("wmean").reset_index()
        out = out.merge(wmean, on=keys, how="left")
    out["n_ccs"] = out["n_ccs"].astype(int)
    return out


def one_sample_by_pair(per_animal: pd.DataFrame, value: str, pairs: Sequence[str],
                       min_n: int = 3) -> pd.DataFrame:
    """Per pair: one-sample *t* of the per-animal ``value`` against 0.

    ``per_animal`` must be one row per animal within each pair (the output of
    :func:`per_animal_mean` at a fixed window / epoch); a duplicate (pair, animal)
    **raises**, because letting it through would silently make windows-as-n or
    epochs-as-n out of an animals-as-n test. Pairs with fewer than ``min_n`` finite
    animals are skipped. Rows come back in the order of ``pairs``.

    Columns: ``pair, n, mean, sem, t, p, bh_pass`` — ``bh_pass`` is Benjamini-
    Hochberg across the pairs returned. Project policy treats the 8 pairs as
    separate families with no cross-pair correction, so ``bh_pass`` is a
    sensitivity column, not the inferential statement.
    """
    dup = per_animal.duplicated(subset=["pair", "animal"])
    if dup.any():
        raise ValueError(f"one_sample_by_pair: {int(dup.sum())} duplicate (pair, animal) "
                         "rows — collapse to one value per animal first")
    rows = []
    for pair in pairs:
        v = pd.to_numeric(per_animal.loc[per_animal["pair"] == pair, value],
                          errors="coerce").to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size < min_n:
            continue
        n, _, t, p = paired_stats.paired_t(v)
        sem = float(np.std(v, ddof=1) / np.sqrt(v.size)) if v.size > 1 else np.nan
        rows.append({"pair": pair, "n": int(n), "mean": float(np.mean(v)),
                     "sem": sem, "t": t, "p": p})
    out = pd.DataFrame(rows, columns=["pair", "n", "mean", "sem", "t", "p"])
    if out.empty:
        out["bh_pass"] = pd.Series(dtype=bool)
        return out
    out["bh_pass"] = paired_stats.fdr_bh(out["p"].to_numpy(float))
    return out
