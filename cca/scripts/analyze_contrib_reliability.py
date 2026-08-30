"""Are the units that carry the communication subspace spatially special?

Joins the per-neuron connection-specific contribution (`contrib_conn` from
run_epochs.py, one CCA per animal × pair × epoch) with per-unit spatial
reliability = mean trial-to-trial Pearson correlation of the spatially-binned
map within ±2 trials (spatial_reliability.trial_map_reliability), averaged over
the SAME trials the epoch's CCA was fit on.

Join key: the weights CSV stores `unit` as the position within the area's kept
units; `dataio.select_units(animal, area, cfg)` re-derives the raw export
index deterministically, so the FS condition of the CSV and of the cfg MUST
match (enforced via --include-fs).

Statistics (animals-as-n, never pooled units — units are nested in animals):
- per (animal, pair, epoch, area): Spearman(contrib, reliability) across units,
  raw AND partialled for log mean rate (rank-OLS residuals), for BOTH
  contrib_conn and the area-intrinsic contrib (if reliability correlates with
  the intrinsic contribution just as strongly, the link is "loud units", not
  communication-specific).
- per (pair, area): Fisher-z mean over epochs within animal -> one value per
  animal -> one-sample t + Wilcoxon vs 0.

Writes results/contrib_reliability_bin10{,_fsincl}.csv   (per-fit correlations)
       results/contrib_reliability_units_bin10{,_fsincl}.csv (per-unit join)
       prints the animals-as-n tables (redirect to keep).

Usage: python scripts/analyze_contrib_reliability.py [--include-fs]
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tom_cca import config, dataio, paired_stats, spatial_reliability  # noqa: E402

EPOCHS = ["naive", "intermediate", "expert"]
PAIR_ORDER = ["CA1-RSC", "CA1-CA3", "CA1-DG", "CA1-V1", "CA3-DG",
              "CA1-SUB", "RSC-SUB", "V1-RSC"]
MIN_UNITS_CORR = 5          # a Spearman over fewer units is noise
HALF_WINDOW = 2             # ±2 trials (spec, 2026-08-29)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--include-fs", action="store_true")
    p.add_argument("--tag", default="_bin10")
    return p.parse_args()


def partial_spearman(x, y, z):
    """Spearman r of x vs y with z rank-partialled out of both (rank-OLS)."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    zc = np.column_stack([np.ones_like(rz), rz])
    res_x = rx - zc @ np.linalg.lstsq(zc, rx, rcond=None)[0]
    res_y = ry - zc @ np.linalg.lstsq(zc, ry, rcond=None)[0]
    if res_x.std() == 0 or res_y.std() == 0:
        return float("nan")
    return float(np.corrcoef(res_x, res_y)[0, 1])


def epoch_trial_sets(animal, cfg, lp):
    """Replicate run_epochs' epoch -> 1-based trial-id sets exactly."""
    streams = dataio._load_temporal_streams(animal, cfg)
    run = ((~np.isnan(streams.trial_idx_50ms))
           & (streams.vel_50ms >= cfg.velocity_thresh_cm_s))
    trial_ids = streams.trial_idx_50ms[run]
    uniq = np.unique(trial_ids)
    ew = dataio.epoch_windows(int(lp), uniq.size, cfg)
    if ew is None:
        return None
    out = {}
    for epoch in EPOCHS:
        pos = ew[epoch]
        pos = pos[(pos >= 0) & (pos < uniq.size)]
        out[epoch] = uniq[pos].astype(int)
    return out


def main():
    args = parse_args()
    suffix = "_fsincl" if args.include_fs else ""
    cfg = dataclasses.replace(config.DEFAULT, temporal_bin_ms=10,
                              gaussian_sd_ms=2.5,
                              exclude_fast_spiking=not args.include_fs)
    wpath = config.RESULTS_DIR / f"epoch_weights{args.tag}{suffix}.csv"
    wdf = pd.read_csv(wpath)
    if "contrib_conn" not in wdf.columns:
        sys.exit(f"{wpath.name} has no contrib_conn column — stale pre-2026-07-28 "
                 f"export; re-run scripts/run_epochs.py first.")
    animals = {a.animal_id: a for a in dataio.load_animals(config.DATA_DIR)}
    behaviour = dataio._read_behaviour_file(config.DATA_DIR / "animal_behaviour.mat")
    entries = dataio.classify_cohort(list(animals.values()), cfg,
                                     behaviour_lookup=behaviour)
    print(f"{wpath.name}: {len(wdf)} weight-rows, "
          f"{wdf['animal'].nunique()} animals | reliability = ±{HALF_WINDOW}-trial "
          f"map correlation on spatial_fr (freq), epoch-matched trials\n")

    unit_rows, corr_rows = [], []
    for aid, adf in wdf.groupby("animal"):
        a = animals[aid]
        rel = spatial_reliability.trial_map_reliability(
            a.spatial_fr, half_window=HALF_WINDOW)
        trial_sets = epoch_trial_sets(a, cfg, int(entries[aid].lp))
        if trial_sets is None:
            print(f"  animal {aid}: epoch_windows None (skip)")
            continue
        # per-epoch per-unit mean reliability, tuning z and mean rate over the
        # SAME trials. Tuning z is (n_units, n_trials) — transposed (see
        # dataio.load_tuning_scores) — so hand epoch_mean_reliability its .T.
        tuning_z = dataio.load_tuning_scores(a.streams_path)["z"]
        # Tom's precomputed moving-window reliability, z vs HIS shuffle null;
        # (trials, units) orientation, so no transpose (unlike tuning_z)
        tomrel_z = dataio.load_reliability_tom(a.streams_path)["z"]
        rel_by_epoch = {e: spatial_reliability.epoch_mean_reliability(rel, t)
                        for e, t in trial_sets.items()}
        tune_by_epoch = {e: spatial_reliability.epoch_mean_reliability(tuning_z.T, t)
                         for e, t in trial_sets.items()}
        tomrel_by_epoch = {e: spatial_reliability.epoch_mean_reliability(tomrel_z, t)
                           for e, t in trial_sets.items()}
        rate_by_epoch = {}
        for e, t in trial_sets.items():
            with np.errstate(invalid="ignore"):
                rate_by_epoch[e] = np.nanmean(a.spatial_fr[t - 1], axis=(0, 1))
        sel = {area: dataio.select_units(a, area, cfg)
               for area in config.AREAS if area in a.area_masks}
        for (pair, epoch, area), g in adf.groupby(["pair", "epoch", "area"]):
            raw_idx = sel[area][g["unit"].to_numpy()]
            r_ep = rel_by_epoch[epoch][raw_idx]
            t_ep = tune_by_epoch[epoch][raw_idx]
            tr_ep = tomrel_by_epoch[epoch][raw_idx]
            rate = rate_by_epoch[epoch][raw_idx]
            for u, ridx in enumerate(raw_idx):
                unit_rows.append({
                    "animal": aid, "pair": pair, "epoch": epoch, "area": area,
                    "unit": int(g["unit"].iloc[u]), "raw_unit": int(ridx),
                    "contrib": g["contrib"].iloc[u],
                    "contrib_conn": g["contrib_conn"].iloc[u],
                    "reliability": round(float(r_ep[u]), 6),
                    "tuning_z": round(float(t_ep[u]), 6),
                    "reliability_tom_z": round(float(tr_ep[u]), 6),
                    "mean_rate": round(float(rate[u]), 6)})
            ok = (np.isfinite(r_ep) & np.isfinite(rate)
                  & np.isfinite(g["contrib_conn"].to_numpy()))
            if ok.sum() < MIN_UNITS_CORR:
                continue
            row = {"animal": aid, "pair": pair, "epoch": epoch, "area": area,
                   "n_units": int(ok.sum())}
            lograte = np.log10(rate[ok] + 1e-3)
            for metric in ("contrib_conn", "contrib"):
                c = g[metric].to_numpy()[ok]
                row[f"rho_{metric}"] = round(
                    float(stats.spearmanr(c, r_ep[ok]).statistic), 4)
                row[f"rho_{metric}_ratepart"] = round(
                    partial_spearman(c, r_ep[ok], lograte), 4)
            row["rho_rate_rel"] = round(
                float(stats.spearmanr(lograte, r_ep[ok]).statistic), 4)
            okt = (np.isfinite(t_ep) & np.isfinite(rate)
                   & np.isfinite(g["contrib_conn"].to_numpy()))
            if okt.sum() >= MIN_UNITS_CORR:
                lograte_t = np.log10(rate[okt] + 1e-3)
                for metric in ("contrib_conn", "contrib"):
                    c = g[metric].to_numpy()[okt]
                    row[f"rho_{metric}_tune"] = round(
                        float(stats.spearmanr(c, t_ep[okt]).statistic), 4)
                    row[f"rho_{metric}_tune_ratepart"] = round(
                        partial_spearman(c, t_ep[okt], lograte_t), 4)
                row["rho_rate_tune"] = round(
                    float(stats.spearmanr(lograte_t, t_ep[okt]).statistic), 4)
            okr = (np.isfinite(tr_ep) & np.isfinite(rate)
                   & np.isfinite(g["contrib_conn"].to_numpy()))
            if okr.sum() >= MIN_UNITS_CORR:
                lograte_r = np.log10(rate[okr] + 1e-3)
                c = g["contrib_conn"].to_numpy()[okr]
                row["rho_contrib_conn_tomrel"] = round(
                    float(stats.spearmanr(c, tr_ep[okr]).statistic), 4)
                row["rho_contrib_conn_tomrel_ratepart"] = round(
                    partial_spearman(c, tr_ep[okr], lograte_r), 4)
                row["rho_rate_tomrel"] = round(
                    float(stats.spearmanr(lograte_r, tr_ep[okr]).statistic), 4)
                both = okr & np.isfinite(r_ep)
                if both.sum() >= MIN_UNITS_CORR:
                    row["rho_rel_tomrel"] = round(float(stats.spearmanr(
                        r_ep[both], tr_ep[both]).statistic), 4)
            corr_rows.append(row)

    udf = pd.DataFrame(unit_rows)
    cdf = pd.DataFrame(corr_rows)
    u_out = config.RESULTS_DIR / f"contrib_reliability_units{args.tag}{suffix}.csv"
    c_out = config.RESULTS_DIR / f"contrib_reliability{args.tag}{suffix}.csv"
    udf.to_csv(u_out, index=False)
    cdf.to_csv(c_out, index=False)
    print(f"wrote {u_out.name} ({len(udf)} unit-rows) + {c_out.name} "
          f"({len(cdf)} fit-rows)\n")

    # ---- animals-as-n: Fisher-z mean over epochs within animal ----
    def animals_as_n(col):
        print("=" * 78)
        print(f"{col}: per-(pair, area) one-sample test of per-animal Fisher-z "
              f"mean rho vs 0")
        print("=" * 78)
        for pair in PAIR_ORDER:
            for area in pair.split("-"):
                g = cdf[(cdf["pair"] == pair) & (cdf["area"] == area)]
                per_an = []
                for _, sub in g.groupby("animal"):
                    z = np.arctanh(np.clip(sub[col].astype(float), -0.999, 0.999))
                    z = z[np.isfinite(z)]
                    if z.size:
                        per_an.append(float(np.tanh(np.mean(z))))
                if len(per_an) < 3:
                    continue
                arr = np.asarray(per_an)
                t_p = float(stats.ttest_1samp(arr, 0.0).pvalue)
                _, med, _, w_p = paired_stats.wilcoxon_signed(per_an)
                mark = lambda p: "*" if (np.isfinite(p) and p < 0.05) else " "
                print(f"  {pair:9s} {area:4s} n={len(arr):<2d} "
                      f"mean_rho={arr.mean():>+7.3f} med={med:>+7.3f} "
                      f"up={int((arr > 0).sum())}/{len(arr)} "
                      f"| t p={t_p:.3g}{mark(t_p)} | W p={w_p:.3g}{mark(w_p)}")
        print()

    for col in ("rho_contrib_conn", "rho_contrib_conn_ratepart",
                "rho_contrib", "rho_contrib_ratepart", "rho_rate_rel",
                "rho_contrib_conn_tune", "rho_contrib_conn_tune_ratepart",
                "rho_contrib_tune", "rho_contrib_tune_ratepart", "rho_rate_tune",
                "rho_contrib_conn_tomrel", "rho_contrib_conn_tomrel_ratepart",
                "rho_rate_tomrel", "rho_rel_tomrel"):
        animals_as_n(col)

    # ---- per-epoch means (descriptive): does the link move with experience? ----
    print("=" * 78)
    print("rho_contrib_conn by epoch (per-animal value, then mean±sem across animals)")
    print("=" * 78)
    for pair in PAIR_ORDER:
        for area in pair.split("-"):
            g = cdf[(cdf["pair"] == pair) & (cdf["area"] == area)]
            if g.empty:
                continue
            cells = []
            for epoch in EPOCHS:
                vals = g[g["epoch"] == epoch].groupby("animal")[
                    "rho_contrib_conn"].mean().to_numpy()
                vals = vals[np.isfinite(vals)]
                cells.append(f"{np.mean(vals):>+6.3f}±{stats.sem(vals):.3f}"
                             f"(n={vals.size})" if vals.size >= 3 else f"{'-':>16}")
            print(f"  {pair:9s} {area:4s} " + "  ".join(cells))


if __name__ == "__main__":
    main()
