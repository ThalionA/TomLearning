"""Run the committed temporal-arm config across the cohort -- Tom.

Two arms (pick via ``--arm``):
* ``runstate`` (Arm A) -- full-traversal, signal-only, segment-fragmented.
  Per-segment circshift null. Output: ``results/temp50_<tag>.pkl``.
* ``landmark`` (Arm B) -- per-landmark 500 ms windows, residual or signal.
  Per-window circshift null. Output: ``results/landmark50_<tag>.pkl``.

Single-config driver (config.DEFAULT plus the CLI overrides) -- the cohort-wide
sweep harness is the follow-up (run_temporal_sweep.py).

Resumable in the same way as run_committed.py: re-running picks up where the
previous run left off, ``--max-seconds`` chunks long runs, ``--fresh`` discards
prior state.

See ``cca/UNDERSTANDING_temporal.md`` for the design contract.

Run:  python scripts/run_temporal.py --arm runstate --data-dir <path>
      python scripts/run_temporal.py --arm landmark --data-dir <path>
"""

from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import dataclasses  # noqa: E402
import functools  # noqa: E402
import multiprocessing as mp  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tom_cca import analysis, config, dataio, pipeline, sweep  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--arm", choices=("runstate", "landmark"), required=True,
                   help="runstate = Arm A; landmark = Arm B")
    p.add_argument("--cca-type", choices=("residual", "signal"), default=None,
                   help="Arm B only. Defaults: residual for landmark, signal "
                        "for runstate (signal forced for Arm A).")
    p.add_argument("--bin-ms", type=int, default=50,
                   help="Single-config bin width (committed run). Use --sweep "
                        "to span 25 and 50 ms together.")
    p.add_argument("--shuffles", type=int, default=200)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--max-seconds", type=float, default=0.0,
                   help="0 = no time cap")
    p.add_argument("--fresh", action="store_true")
    p.add_argument("--include-fs", action="store_true")
    p.add_argument("--data-dir", type=str, default=None)
    p.add_argument("--trials-per-epoch", type=int, default=0,
                   help="0 = use config default")
    p.add_argument("--sweep", action="store_true",
                   help="Iterate the full bin_ms x cca x k_rule grid "
                        "(landmark: 44 configs; runstate: 22 configs). "
                        "Each config writes its own pkl; resumable per config.")
    p.add_argument("--tags", type=str, default=None,
                   help="Comma-separated tag filter for --sweep "
                        "(e.g. landmark50_res_samp15,landmark25_res_samp15)")
    return p.parse_args()


def _key(o):
    return (o.animal_id, o.area_x, o.area_y)


def _save(path: Path, done: dict, cfg) -> None:
    tmp = path.with_suffix(".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump({"results": list(done.values()), "cfg": cfg}, fh)
    os.replace(tmp, path)


def _build_cfg(args) -> config.Config:
    repl: dict = {
        "n_shuffles": int(args.shuffles),
        "temporal_bin_ms": int(args.bin_ms),
    }
    if args.arm == "runstate":
        # Arm A: signal-only forced (ragged trials, no PSTH).
        repl["bin_mode"] = "temporal_runstate"
        repl["subtract_trial_mean"] = False
        if args.cca_type == "residual":
            print("WARNING: --cca-type residual is invalid for Arm A "
                  "(ragged trials); forcing signal.", file=sys.stderr)
    else:
        # Arm B: residual is the committed default; signal available.
        repl["bin_mode"] = "landmark"
        if args.cca_type == "signal":
            repl["subtract_trial_mean"] = False
    if args.include_fs:
        repl["exclude_fast_spiking"] = False
    if args.trials_per_epoch:
        repl["trials_per_epoch"] = int(args.trials_per_epoch)
    return dataclasses.replace(config.DEFAULT, **repl)


def _tag(cfg: config.Config, args) -> str:
    cca = "res" if cfg.subtract_trial_mean else "sig"
    fs = "_fsincl" if args.include_fs else ""
    prefix = "temp" if args.arm == "runstate" else "landmark"
    return (f"{prefix}{int(cfg.temporal_bin_ms):02d}_{cca}"
            f"_samp{cfg.samples_per_pc:02d}{fs}")


def _worker_runstate(arg, cfg: config.Config):
    animal, entry, area_x, area_y = arg
    prep = pipeline.prepare_pair_temporal(animal, area_x, area_y, entry, cfg)
    if isinstance(prep, pipeline.SkippedPair):
        return _key(prep), prep
    res = analysis.analyse_pair_temporal(prep, cfg)
    return _key(res), res


def _worker_landmark(arg, cfg: config.Config):
    animal, entry, area_x, area_y = arg
    prep = pipeline.prepare_pair_landmark(animal, area_x, area_y, entry, cfg)
    if isinstance(prep, pipeline.SkippedPair):
        return _key(prep), prep
    res = analysis.analyse_pair_landmark(prep, cfg)
    return _key(res), res


def _run_one_config(tag: str, cfg: config.Config, learners, entries, args):
    """Fit one config across the cohort, writing/resuming results/<tag>.pkl."""
    out = config.RESULTS_DIR / f"{tag}.pkl"
    done_marker = config.RESULTS_DIR / f"{tag}.done"
    out.parent.mkdir(parents=True, exist_ok=True)

    todo = []
    for a in learners:
        for area_x, area_y in config.PAIRS:
            todo.append((a, entries[a.animal_id], area_x, area_y))

    if not args.fresh and out.is_file():
        with open(out, "rb") as fh:
            prior = pickle.load(fh)
        done = {_key(r): r for r in prior["results"]}
        todo = [j for j in todo if (j[0].animal_id, j[2], j[3]) not in done]
        print(f"[{tag}] resuming; {len(done)} prior, {len(todo)} todo")
    else:
        done = {}
        if args.fresh and out.is_file():
            out.unlink()

    worker = _worker_runstate if args.arm == "runstate" else _worker_landmark
    f_worker = functools.partial(worker, cfg=cfg)
    start = time.time()

    if args.jobs == 1:
        for j in todo:
            key, res = f_worker(j)
            done[key] = res
            if args.max_seconds and time.time() - start > args.max_seconds:
                break
    else:
        with mp.Pool(args.jobs) as pool:
            for key, res in pool.imap_unordered(f_worker, todo):
                done[key] = res
                if args.max_seconds and time.time() - start > args.max_seconds:
                    pool.terminate()
                    break
                if len(done) % 10 == 0:
                    _save(out, done, cfg)
                    print(f"  [{tag}] {len(done)} done, "
                          f"{time.time()-start:.0f}s elapsed")

    _save(out, done, cfg)
    done_marker.write_text(f"{len(done)} results\n{time.ctime()}\n")
    print(f"[{tag}] wrote {out} ({len(done)} pair results, "
          f"{time.time()-start:.0f}s)")
    return out


def _sweep_configs(args) -> list[tuple[str, config.Config]]:
    """Resolve the sweep grid (filtered by --tags / --cca-type / --bin-ms)."""
    name = ("temporal_runstate" if args.arm == "runstate"
            else "temporal_landmark")
    grid = sweep.build_sweep(name)
    # Apply n_shuffles + FS + trials_per_epoch overrides on every cfg.
    repl_common: dict = {"n_shuffles": int(args.shuffles)}
    if args.include_fs:
        repl_common["exclude_fast_spiking"] = False
    if args.trials_per_epoch:
        repl_common["trials_per_epoch"] = int(args.trials_per_epoch)
    grid = [(tag, dataclasses.replace(cfg, **repl_common)) for tag, cfg in grid]
    if args.tags:
        wanted = {t.strip() for t in args.tags.split(",") if t.strip()}
        grid = [(t, c) for t, c in grid if t in wanted]
        if not grid:
            raise SystemExit(f"--tags filter matched no configs: {wanted}")
    return grid


def main():
    args = parse_args()

    data_dir = Path(args.data_dir) if args.data_dir else config.DATA_DIR
    base_cfg = _build_cfg(args)            # used for cohort classification
    animals = dataio.load_animals(data_dir)
    behaviour = dataio._read_behaviour_file(data_dir / "animal_behaviour.mat")
    entries = dataio.classify_cohort(animals, base_cfg, behaviour_lookup=behaviour)
    learners = [a for a in animals if a.animal_id in entries]
    print(f"cohort: {len(learners)} learner animals; "
          f"{len(config.PAIRS)} pairs/animal")

    if args.sweep:
        grid = _sweep_configs(args)
        print(f"sweep: {len(grid)} configs to run "
              f"({args.shuffles} shuffles each)")
        for tag, cfg in grid:
            _run_one_config(tag, cfg, learners, entries, args)
            if args.max_seconds:
                # Per-config budget; if we hit it on this config, stop.
                pass
    else:
        cfg = base_cfg
        tag = _tag(cfg, args)
        _run_one_config(tag, cfg, learners, entries, args)


if __name__ == "__main__":
    main()
