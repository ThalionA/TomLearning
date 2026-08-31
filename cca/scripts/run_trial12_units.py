"""Per-UNIT carrying of the frozen communication subspace, trial 1 vs 2.

Companion to run_trial12.py (same gates, same reference design: subspace fit
once per (animal, pair) on running-trial ordinals > N_HELD, ordinals 1..N_HELD
read out leak-free). Where run_trial12 exports pair-level lag curves, this
driver asks WHICH UNITS carry the communication on each early trial:

  carry_r(unit, ordinal) = Pearson r between the unit's residualised activity
  (fit.Xr / fit.Yr — third area partialled with train-only coefficients) and
  the PARTNER area's frozen CC1 variate, over that ordinal's running bins
  (membership.variate_structure_coefficients — the method-agnostic membership).

Ordinal 0 = the REFERENCE profile (same correlation over the fit rows,
ordinals > N_HELD). No permutation nulls here — inference lives downstream in
analyze_trial12_units.py (deltas + adjacent-trial control band, animals-as-n).

MEASURED CONFOUND carried over from run_trial12: trial 1 has more running bins
than later trials, and correlation estimates sharpen with n, so every score is
exported in two arms (`matched`):
  0  raw     — all of the ordinal's running bins (ordinal 0 = all fit rows)
  1  common  — every ordinal 1..N_HELD cut to the min bin count over the gated
               ordinals (leading CONSECUTIVE bins; ordinal 0 stays raw)

Writes results/trial12_units_bin10{,_fsincl}.csv
Usage: python scripts/run_trial12_units.py [--include-fs] [--animals 36]
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (PAIRS, animals_filter, cfg_from_args, results_path,  # noqa: E402
                     temporal_parser)
from run_trial12 import K, MIN_TRIALS, N_HELD  # noqa: E402
from tom_cca import config, dataio, early_trials, membership, preprocess  # noqa: E402

MIN_BINS = 200                    # same per-trial gate as run_trial12's lag=0 floor
FIELDS = ["animal", "learner", "pair", "area", "side", "unit", "raw_unit",
          "ordinal", "trial_id", "matched", "carry_r", "n_bins", "n_fit_trials"]


def carry_profile(res_units, partner_variate, rows):
    """(n_units,) correlation of each unit with the partner's CC1 on ``rows``."""
    S = membership.variate_structure_coefficients(
        res_units[rows], partner_variate[rows, :1])
    return S[:, 0]


def main():
    ap = temporal_parser("per-unit frozen-subspace carrying, ordinals 1..10",
                         cap=False, shuffles=False, max_lag=False)
    args = ap.parse_args()
    cfg = cfg_from_args(args)
    out = results_path("trial12_units", fs=args.include_fs, bin_ms=args.bin_ms)
    print(f"TRIAL12 UNITS | bin={args.bin_ms}ms K={K} | fit on ordinals >{N_HELD} "
          f"| gate >={MIN_BINS} bins/trial | FS={'incl' if args.include_fs else 'excl'}")

    animals = dataio.load_animals(config.DATA_DIR)
    entries = dataio.classify_cohort(animals, cfg,
                                     behaviour_lookup=dataio.load_learning_points())
    only = animals_filter(args.animals)
    rows_out = []

    def _flush():
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
            w.writeheader(); w.writerows(rows_out)

    for a in animals:
        if only is not None and int(a.animal_id) not in only:
            continue
        sess = preprocess.load_running_session(a, cfg, entries, min_trials=MIN_TRIALS)
        if isinstance(sess, str):
            print(f"  animal {a.animal_id}: SKIP — {sess}"); continue
        uniq = sess.trials
        train_ids = uniq[N_HELD:]
        train_mask = np.isin(sess.trial_ids, train_ids)
        sel = {area: dataio.select_units(a, area, cfg)
               for area in config.AREAS if area in sess.areas}

        for ax, ay in PAIRS:
            pair = sess.pair(ax, ay)
            if pair is None:
                continue
            X, Y, Z = pair
            try:
                fit = early_trials.reference_fit(X, Y, train_mask, Z=Z, k=K)
            except Exception as e:                              # noqa: BLE001
                print(f"  fit fail {a.animal_id} {ax}-{ay}: {e}"); continue
            u_all, v_all = early_trials.variates(fit, slice(None))

            # gated early ordinals and their leading-consecutive row blocks
            blocks = {}
            for o in range(1, N_HELD + 1):
                rows = np.flatnonzero(sess.trial_ids == uniq[o - 1])
                if rows.size >= MIN_BINS:
                    blocks[o] = rows
            if not blocks:
                continue
            common = min(r.size for r in blocks.values())

            base = {"animal": a.animal_id, "learner": int(sess.learner),
                    "pair": f"{ax}-{ay}", "n_fit_trials": int(train_ids.size)}
            # side x carries via the PARTNER's variate v, and vice versa
            for area, side, res, partner_var in ((ax, "x", fit.Xr, v_all),
                                                 (ay, "y", fit.Yr, u_all)):
                raw_idx = sel[area]
                jobs = [(0, -1, 0, np.flatnonzero(train_mask))]
                for o, rows in blocks.items():
                    tid = int(uniq[o - 1])
                    jobs.append((o, tid, 0, rows))
                    jobs.append((o, tid, 1, rows[:common]))
                for o, tid, matched, rows in jobs:
                    prof = carry_profile(res, partner_var, rows)
                    for ui, r in enumerate(prof):
                        rows_out.append({
                            **base, "area": area, "side": side, "unit": ui,
                            "raw_unit": int(raw_idx[ui]), "ordinal": o,
                            "trial_id": tid, "matched": matched,
                            "carry_r": round(float(r), 5),
                            "n_bins": int(rows.size)})
        _flush()
        print(f"  animal {a.animal_id}: done ({len(rows_out)} rows)")
    _flush()
    print(f"-> {out}")


if __name__ == "__main__":
    main()
