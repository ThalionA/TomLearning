"""Merge the per-animal-group parts written by run_lag_curves_uncapped_batch.sh into
results/lag_curves_bin10{,_fsincl}.csv — animal order ascending, within-animal row
order as the driver wrote it (pair -> lag -> dim). Refuses to merge if any animal
appears in two parts or a part is empty; prints the row count and animal list so the
result can be checked against the previous (capped) file's coverage.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

RES = Path(__file__).resolve().parents[1] / "results"


def merge(suffix: str) -> int:
    parts = sorted(RES.glob(f"lag_curves_bin10{suffix}_part*.csv"))
    if not parts:
        print(f"no parts for '{suffix or 'fsexcl'}'"); return 1
    frames = []
    for p in parts:
        d = pd.read_csv(p)
        if d.empty:
            print(f"EMPTY part {p.name} — refusing to merge"); return 1
        frames.append(d)
    all_ = pd.concat(frames, ignore_index=True)
    per_part = [set(f["animal"].unique()) for f in frames]
    for i in range(len(per_part)):
        for j in range(i + 1, len(per_part)):
            if per_part[i] & per_part[j]:
                print(f"animal in two parts: {per_part[i] & per_part[j]} — refusing"); return 1
    all_ = all_.sort_values("animal", kind="mergesort").reset_index(drop=True)
    out = RES / f"lag_curves_bin10{suffix}.csv"
    all_.to_csv(out, index=False, lineterminator="\n")
    print(f"{out.name}: {len(all_)} rows, animals {sorted(all_['animal'].unique())}, "
          f"{all_.groupby(['animal', 'pair']).ngroups} animal-pairs")
    return 0


if __name__ == "__main__":
    rc = merge("") | merge("_fsincl")
    sys.exit(rc)
