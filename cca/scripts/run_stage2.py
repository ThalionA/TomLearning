"""Stage 2 driver for the Tom-learning CCA sweep.

Runs lagged CCA + held-out-CC surrogate significance for every config in the
spatial sweep (see ``tom_cca.sweep``):

  * residual / signal CCA   x  z-scoring on/off
  * min-units 4 / 6 / 10    x  LP criterion 7 / 8
  * PC-count rule           x  AXIS_KRULE  (samples/fixed/variance, 11 modes)

= 264 configs. FS is fixed-excluded, bin width fixed at 200 (per the user's
brief; the bin-width and FS sweep axes from the striatum pipeline are removed).

Each config saves results/stage2_<tag>.pkl, resumable per config; completed
pairs are saved incrementally (atomically). ``--max-seconds`` stops cleanly
after a time budget so a long sweep can be driven in short chunks -- just
re-run to resume.

Run:  python scripts/run_stage2.py [--shuffles N] [--jobs N]
                                    [--max-seconds S] [--fresh]
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

EPOCHS = config.EPOCH_NAMES


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sweep", choices=("spatial",), default="spatial")
    p.add_argument("--shuffles", type=int, default=200)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--max-seconds", type=float, default=0.0,
                   help="stop cleanly after this many seconds (0 = no limit)")
    p.add_argument("--fresh", action="store_true", help="ignore saved results")
    p.add_argument("--data-dir", type=str, default=None,
                   help="override config.DATA_DIR (where TF*_export.mat lives)")
    return p.parse_args()


def _key(obj):
    return (obj.animal_id, obj.area_x, obj.area_y)


def _save(path, done, skipped, cfg):
    """Atomic pickle write -- a kill mid-write cannot corrupt the resume file."""
    tmp = path.with_suffix(".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump({"results": list(done.values()), "skipped": skipped,
                     "cfg": cfg}, fh)
    os.replace(tmp, path)


def main() -> None:
    args = parse_args()
    deadline = time.time() + args.max_seconds if args.max_seconds else None
    grid = sweep.build_sweep(args.sweep)
    print(f"Sweep '{args.sweep}': {len(grid)} configs; "
          f"{args.shuffles} surrogates; jobs={args.jobs}.\n")

    animals = dataio.load_animals(args.data_dir,
                                  spatial_field=config.DEFAULT.spatial_field)
    behaviour = dataio._read_behaviour_file(
        (Path(args.data_dir) if args.data_dir else config.DATA_DIR)
        / config.LEARNING_FILE
    )

    n_done = 0
    for tag, cfg in grid:
        marker = config.RESULTS_DIR / f"stage2_{tag}.done"
        if marker.exists() and not args.fresh:        # O(1) resume skip
            n_done += 1
            continue
        cfg = dataclasses.replace(cfg, n_shuffles=args.shuffles)
        entries = dataio.classify_cohort(animals, cfg, behaviour)
        if run_config(animals, entries, cfg, tag, args, deadline):
            n_done += 1
        if deadline is not None and time.time() > deadline:
            break

    done = n_done == len(grid)
    print(f"\n{n_done}/{len(grid)} configs complete -- "
          + ("SWEEP COMPLETE." if done else "re-run to resume."))


def run_config(animals, entries, cfg, tag: str, args, deadline) -> bool:
    """Prepare, analyse and save one configuration. Returns True if complete."""
    print(f"=== config '{tag}' ({pipeline.config_label(cfg)}) ===")
    prepared, skipped = [], []
    for animal in animals:
        if animal.animal_id not in entries:        # learners only
            continue
        entry = entries[animal.animal_id]
        for area_x, area_y in config.PAIRS:
            result = pipeline.prepare_pair(animal, area_x, area_y, entry, cfg)
            (prepared if isinstance(result, pipeline.PreparedPair)
             else skipped).append(result)

    out_path = config.RESULTS_DIR / f"stage2_{tag}.pkl"
    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    done = {}
    if out_path.exists() and not args.fresh:
        with open(out_path, "rb") as fh:
            done = {_key(r): r for r in pickle.load(fh)["results"]}
    todo = [p for p in prepared if _key(p) not in done]
    print(f"  {len(prepared)} pairs prepared, {len(done)} done, "
          f"{len(todo)} to do.")

    if todo:
        worker = functools.partial(analysis.analyse_pair, cfg=cfg)
        t0 = time.time()
        paused = False
        with mp.Pool(args.jobs) as pool:
            for i, res in enumerate(pool.imap_unordered(worker, todo), start=1):
                done[_key(res)] = res
                if i % 6 == 0:
                    _save(out_path, done, skipped, cfg)
                    print(f"    {len(done)}/{len(prepared)} "
                          f"({time.time() - t0:.0f}s)")
                if deadline is not None and time.time() > deadline:
                    paused = True
                    break
        _save(out_path, done, skipped, cfg)
        print(f"  {'paused' if paused else 'saved'} {out_path.name} "
              f"({time.time() - t0:.0f}s)")

    complete = len(done) == len(prepared)
    if complete:
        out_path.with_suffix(".done").touch()          # O(1) resume marker
    print(f"  config '{tag}' {'complete' if complete else 'INCOMPLETE'}.\n")
    return complete


if __name__ == "__main__":
    main()
